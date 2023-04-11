#pragma once

#include <fstream>
#include <onnx.pb.h>
#include <vector>
#include <gravity/model.h>
#include <gravity/func.h>

typedef std::map<std::string, gravity::param<float>> Initializers;
typedef std::map<std::string, gravity::func<float>> HiddenStates;

// Define the base Layer class
class Layer {
public:
    Layer(const onnx::NodeProto& node, const Initializers& global_initializers) {
        for (const auto& input : node.input()) {
            if (global_initializers.find(input) != global_initializers.end()) {
                this->initializers[input] = global_initializers.at(input);
            }
            this->inputs.push_back(input);
        }
        for (const auto& output : node.output()) {
            this->outputs.push_back(output);
        }
        this->name = node.name();
    }

    virtual ~Layer() = default;

    virtual void forward(HiddenStates& hidden_states) = 0;
    virtual void print() = 0;

    void print_io() {
        std::cout << "| inputs: " << std::endl;
        for (const auto& input : this->inputs) {
            std::cout << "|   " << input << std::endl;
        }
        std::cout << "| outputs: " << std::endl;
        for (const auto& output : this->outputs) {
            std::cout << "|   " << output << std::endl;
        }
    }

    const onnx::AttributeProto* find_attribute(const onnx::NodeProto& node, const std::string& name) {
        for (const auto& attr : node.attribute()) {
            if (attr.name() == name) {
                return &attr;
            }
        }
        return nullptr;
    }

    gravity::func<float> get_parameter(const std::string& name, const HiddenStates& hidden_states) {
        std::cout << "  Getting parameter: " << name << std::endl;
        // Check if its an initializer
        if (this->initializers.find(name) != this->initializers.end()) {
            std::cout << "    Found initializer: " << name << std::endl;
            auto ret = this->initializers.at(name);
            ret.print();
            ret.print_symbolic();
            return this->initializers.at(name);
        } else if (hidden_states.find(name) != hidden_states.end()) {
            std::cout << "    Found hidden state: " << name << std::endl;
            auto ret = hidden_states.at(name);
            ret.print();
            ret.print_symbolic();
            return ret;
        } else {
            throw std::runtime_error("Could not find parameter: " + name);
        }
    }

    // Name of the layer
    std::string name;

    std::vector<std::string> inputs;
    std::vector<std::string> outputs;

    // Map of input names to initializers (generally weights, biases, etc.)
    Initializers initializers;
};

class GEMM : public Layer {
public:
    /*
        Y = alpha * A’ * B’ + beta * C
    */
    GEMM(const onnx::NodeProto& node, Initializers& global_initializers): Layer(node, global_initializers) {
        const auto* alpha_attr = find_attribute(node, "alpha");
        // print
        if (alpha_attr) {
            this->alpha = alpha_attr->f();
        }
        const auto* beta_attr = find_attribute(node, "beta");
        if (beta_attr) {
            this->beta = beta_attr->f();
        }
        const auto *transA_attr = find_attribute(node, "transA");
        if (transA_attr) {
            this->transA = transA_attr->i();
        }
        const auto *transB_attr = find_attribute(node, "transB");
        if (transB_attr) {
            this->transB = transB_attr->i();
        }

        if (this->inputs.size() == 3) {
            this->has_optional_C = true;
        }
    }

    void forward(HiddenStates& hidden_states) override {
        
        auto A = this->get_parameter(this->inputs.at(0), hidden_states);
        auto B = this->get_parameter(this->inputs.at(1), hidden_states);
        A._is_transposed = transA;
        B._is_transposed = transB;
        auto Y = this->alpha * A * B;
        Y.eval_all();
        if (this->has_optional_C) {
            auto C = this->get_parameter(this->inputs[2], hidden_states);
            Y += this->beta * C;
        }
        hidden_states[this->outputs[0]] = Y;
    }

    void print() override {
        std::cout << "---------------------------------" << std::endl;
        std::cout << "| GEMM: " << this->name << std::endl;
        std::cout << "---------------------------------" << std::endl;
        std::cout << "| alpha: "  << this->alpha << std::endl;
        std::cout << "| beta: "   << this->beta << std::endl;
        std::cout << "| transA: " << this->transA << std::endl;
        std::cout << "| transB: " << this->transB << std::endl;
        this->print_io();
        std::cout << "---------------------------------" << std::endl;
    }

    float alpha = 1.0;
    float beta = 1.0;
    bool transA = false;
    bool transB = false;
    bool has_optional_C = false;
};

class ReLU : public Layer {
public:
    ReLU(const onnx::NodeProto& node, const Initializers& global_initializers): Layer(node, global_initializers) {}
    void forward(HiddenStates& hidden_states) override {
        auto x = this->get_parameter(this->inputs[0], hidden_states);
        gravity::param<float> zero("zero");
        zero.set_size(x.get_dim());
        auto output = gravity::max(x, zero);
        hidden_states[this->outputs[0]] = output;
    }

    void print() override {
        std::cout << "---------------------------------" << std::endl;
        std::cout << "| ReLU: " << this->name << std::endl;
        std::cout << "---------------------------------" << std::endl;
        this->print_io();
        std::cout << "---------------------------------" << std::endl;
    }
};

gravity::param<float> parse_tensor(const onnx::TensorProto& tensor) {
    std::list<gravity::indices> dims;
    std::vector<size_t> shape;
    gravity::indices ids;
    for (int i = 0; i < tensor.dims_size(); i++) {
        dims.push_back(gravity::range(0,tensor.dims(i)-1));
        shape.push_back(tensor.dims(i));
    }
    ids = gravity::indices(dims);
    gravity::param<float> weight(tensor.name());
    weight.in(ids);

    std::cout << "Parsing weight: [";
    for (int i = 0; i < shape.size(); i++) {
        std::cout << shape[i];
        if (i != shape.size() - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "] | Name: " << tensor.name() << std::endl;

    std::vector<float> data_;
    // Sometimes the data is stored in raw_data, sometimes in float_data
    // Unclear why this is the case
    if (tensor.raw_data().size() > 0) {
        const void* raw_data = tensor.raw_data().data();
        const size_t num_bytes = tensor.raw_data().size();
        const size_t num_floats = num_bytes / sizeof(float);
        data_ = std::vector<float>(num_floats);
        std::memcpy(data_.data(), raw_data, num_bytes);
    } else {
        const size_t num_elements = tensor.float_data_size();
        data_ = std::vector<float>(num_elements);
        for (size_t i = 0; i < num_elements; ++i) {
            data_.at(i) = tensor.float_data(i);
        }
    }
    for (int i = 0; i < data_.size(); i++) {
        weight.set_val(i, data_.at(i));
    }
    return weight;
}
