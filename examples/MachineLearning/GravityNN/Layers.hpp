#pragma once

#include <fstream>
#include <onnx.pb.h>
#include <vector>
#include <gravity/model.h>
#include <gravity/func.h>

typedef std::map<std::string, onnx::TensorProto> Initializers;
typedef std::map<std::string, gravity::func<float>> HiddenStates;

typedef enum { _gemm, _relu, _conv, _matmul, _add} OType; /* Operator Type */

vector<size_t> get_input_dim(const ::google::protobuf::RepeatedPtrField< ::onnx::ValueInfoProto > &info)
{
    vector<size_t> dims;
  for (auto input_data: info)
  {
    auto shape = input_data.type().tensor_type().shape();
    std::cout << "  " << input_data.name() << ":";
    std::cout << "[";
    if (shape.dim_size() != 0)
    {
      int size = shape.dim_size();
      for (int i = 0; i < size; ++i)
      {
        dims.push_back(shape.dim(i).dim_value());
      }
    }
  }
    return dims;
}

class Tensor {
public:
    Tensor() {}

    Tensor(std::string name, const Initializers& global_initializers) {
        this->name = name;

        if (global_initializers.count(name) == 0) {
            this->is_initializer = false;
            return;
        }

        auto& tensor = global_initializers.at(name);

        this->is_initializer = true;
        this->shape = std::vector<size_t>(tensor.dims().begin(), tensor.dims().end());

        if (tensor.raw_data().empty()) {
            throw std::runtime_error("Tensor " + tensor.name() + " has no data");
        }

        const void* raw_data = tensor.raw_data().data();
        const size_t num_bytes = tensor.raw_data().size();
        const size_t num_floats = num_bytes / sizeof(float);

        data.resize(num_floats);
        std::memcpy(data.data(), raw_data, num_bytes);
    }

    gravity::func<float> get(HiddenStates& hidden_states) const {
        if (!is_initializer) {
            return hidden_states.at(name);
        }

        std::list<gravity::indices> dims;
        for (auto dim : shape) {
            dims.push_back(gravity::range(0, dim - 1));
        }
        gravity::indices ids(dims);

        gravity::param<float> weight(name);
        weight.in(ids);

        size_t i = 0;
        for (auto val : data) {
            weight.set_val(i++, val);
        }

        return weight;
    }

    bool is_initializer;
    std::string name;
    std::vector<float> data;
    std::vector<size_t> shape;
    
};

// Define the base Layer class
class Layer {
public:
    Layer(const onnx::NodeProto& node, const Initializers& global_initializers) {
        for (const auto& input : node.input()) {
            this->inputs.push_back(input);
        }
        for (const auto& output : node.output()) {
            this->outputs.push_back(output);
        }
        this->name = node.name();
    }

    virtual ~Layer() = default;

    virtual void forward(HiddenStates& hidden_states) = 0;
    virtual void print() const = 0;

    void print_io() const {
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

    // Name of the layer
    std::string name;
    OType operator_type;
    bool is_activation_func = false;
    std::vector<std::string> inputs;
    std::vector<std::string> outputs;
    std::vector<size_t> var_dims;
};

/* Matrix multiply results from A * B */
class MatMul : public Layer {
public:
    /*
        Y = A * B
    */
    MatMul(const onnx::NodeProto& node, const Initializers& global_initializers): Layer(node, global_initializers) {
        operator_type = _matmul;
        
        
        this->A = Tensor(this->inputs.at(0), global_initializers);
        this->B = Tensor(this->inputs.at(1), global_initializers);
        if(this->B.shape.size()==2){
            this->var_dims = this->B.shape;
        }
    }

    void forward(HiddenStates& hidden_states) override {
        auto fA = this->A.get(hidden_states);
        auto fB = this->B.get(hidden_states);
        auto Y = fA * fB;
        Y.eval_all();
        hidden_states[this->outputs[0]] = Y;
    }

    void print() const override{
        std::cout << "---------------------------------" << std::endl;
        std::cout << "| MatMul: " << this->name << std::endl;
        std::cout << "---------------------------------" << std::endl;
        this->print_io();
        std::cout << "---------------------------------" << std::endl;
    }

    Tensor A;
    Tensor B;

};

/* Performs element-wise binary addition */
class Add : public Layer {
public:
    /*
        Y = A + B
    */
    Add(const onnx::NodeProto& node, const Initializers& global_initializers): Layer(node, global_initializers) {
        operator_type = _add;
        this->A = Tensor(this->inputs.at(0), global_initializers);
        this->B = Tensor(this->inputs.at(1), global_initializers);
        this->var_dims = this->B.shape;
    }

    void forward(HiddenStates& hidden_states) override {
        auto fA = this->A.get(hidden_states);
        auto fB = this->B.get(hidden_states);
        auto Y = fA + fB;
        Y.eval_all();
        hidden_states[this->outputs[0]] = Y;
    }

    void print() const override{
        std::cout << "---------------------------------" << std::endl;
        std::cout << "| Add: " << this->name << std::endl;
        std::cout << "---------------------------------" << std::endl;
        this->print_io();
        std::cout << "---------------------------------" << std::endl;
    }

    Tensor A;
    Tensor B;

};

class GEMM : public Layer {
public:
    /*
        Y = alpha * A’ * B’ + beta * C
    */
    GEMM(const onnx::NodeProto& node, const Initializers& global_initializers): Layer(node, global_initializers) {
        operator_type = _gemm;
        const auto* alpha_attr = find_attribute(node, "alpha");
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
        
        this->A = Tensor(this->inputs.at(0), global_initializers);
        this->B = Tensor(this->inputs.at(1), global_initializers);
        if(this->B.shape.size()==2){
            this->var_dims = this->B.shape;
            if(this->transB){
                this->var_dims[0] = this->B.shape[1];
                this->var_dims[1] = this->B.shape[0];
            }
        }
        if (this->inputs.size() == 3) {
            this->C = Tensor(this->inputs.at(2), global_initializers);
            this->has_optional_C = true;
        }

        
    }

    void forward(HiddenStates& hidden_states) override {
        auto fA = this->A.get(hidden_states);
        auto fB = this->B.get(hidden_states);
        fA._is_transposed = transA;
        fB._is_transposed = transB;
        auto Y = this->alpha * fA * fB;
        Y.eval_all();
        if (this->has_optional_C) {
            auto fC = this->C.get(hidden_states);
            Y += this->beta * fC;
        }
        hidden_states[this->outputs[0]] = Y;
    }

    void print() const override{
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

    Tensor A;
    Tensor B;
    Tensor C;

    float alpha = 1.0;
    float beta = 1.0;
    bool transA = false;
    bool transB = false;
    bool has_optional_C = false;
};


class ReLU : public Layer {
public:
    ReLU(const onnx::NodeProto& node, const Initializers& global_initializers): Layer(node, global_initializers) {
        operator_type = _relu;
        this->X = Tensor(this->inputs.at(0), global_initializers);
        this->is_activation_func = true;
    }
    void forward(HiddenStates& hidden_states) override {
        auto fX = this->X.get(hidden_states);

        gravity::param<float> zero("zero");
        zero.set_size(fX.get_dim());

        auto output = gravity::max(fX, zero);
        hidden_states[this->outputs[0]] = output;
    }

    void print() const override{
        std::cout << "---------------------------------" << std::endl;
        std::cout << "| ReLU: " << this->name << std::endl;
        std::cout << "---------------------------------" << std::endl;
        this->print_io();
        std::cout << "---------------------------------" << std::endl;
    }

    Tensor X;
};
