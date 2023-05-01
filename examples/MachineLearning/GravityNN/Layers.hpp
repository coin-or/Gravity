#pragma once

#include <fstream>
#include <onnx.pb.h>
#include <vector>
#include <gravity/model.h>
#include <gravity/func.h>

typedef std::map<std::string, onnx::TensorProto> Initializers;
typedef std::map<std::string, onnx::ValueInfoProto> TensorShapes;
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

    Tensor(std::string name, const Initializers& global_initializers, const TensorShapes& tensor_shapes) {
        this->name = name;

        this->is_initializer = global_initializers.count(name) != 0;

        this->_get_metadata(tensor_shapes, global_initializers);

        if (!this->is_initializer) {
            return;
        }

        auto& tensor = global_initializers.at(name);
        if (!tensor.raw_data().empty()) {
            // throw std::runtime_error("Tensor " + tensor.name() + " has no data");
            const void* raw_data = tensor.raw_data().data();
            const size_t num_bytes = tensor.raw_data().size();
            const size_t num_floats = num_bytes / sizeof(float);

            data.resize(num_floats);
            std::memcpy(data.data(), raw_data, num_bytes);
        } else {
//            std::cout << "Reading float data" << std::endl;
            this->data.clear();
            for (auto val : tensor.float_data()) {
                data.push_back(val);
            }
        }
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

    void _transpose() {
        if (this->shape.size() != 2) {
            throw std::runtime_error("Cannot transpose tensor with shape " + std::to_string(this->shape.size()));
        }
        if (!this->is_initializer) {
            throw std::runtime_error("Cannot transpose non-initializer tensor");
        }

        // Tranpose data if needed, it is stored row-major
        std::vector<float> temp_data = this->data;
        for (size_t i = 0; i < this->shape[0]; i++) {
            for (size_t j = 0; j < this->shape[1]; j++) {
                this->data[j * this->shape[0] + i] = temp_data[i * this->shape[1] + j];
            }
        }

        // Swap shape
        std::swap(this->shape[0], this->shape[1]);
    }

    void _get_metadata(const TensorShapes& tensor_shapes, const Initializers& global_initializers) {
        if (tensor_shapes.count(name) > 0) {
            auto& vinfo = tensor_shapes.at(name);
            this->shape.clear();
            for (auto dim : vinfo.type().tensor_type().shape().dim()) {
                this->shape.push_back(dim.dim_value());
            }
        } else if (global_initializers.count(name) > 0) {
            auto& tensor = global_initializers.at(name);
            this->shape = std::vector<size_t>(tensor.dims().begin(), tensor.dims().end());
        } else {
            throw std::runtime_error("Tensor " + name + " is not an initializer and has no value info. Please run shape inference on model with *NO* dynamic axes.");
        }
    }

    std::string name;
    std::vector<float> data;
    std::vector<size_t> shape;

    bool is_initializer;
};

// Define the base Layer class
class Layer {
public:
    Layer(const onnx::NodeProto& node, const Initializers& global_initializers, const TensorShapes& shapes) {
        for (const auto& input : node.input()) {
            this->inputs.push_back(Tensor(input, global_initializers, shapes));
        }
        for (const auto& output : node.output()) {
            this->outputs.push_back(Tensor(output, global_initializers, shapes));
        }
        this->name = node.name();

        // Try to find upper and lower bounds for each output
        for (auto out: this->outputs)
        {
            auto lower_name = out.name + "_lower";
            if (global_initializers.count(lower_name) != 0) {
                auto lower = Tensor(lower_name, global_initializers, shapes);
                this->lowers.push_back(lower);
            }
            auto upper_name = out.name + "_upper";
            if (global_initializers.count(lower_name) != 0) {
                auto upper = Tensor(upper_name, global_initializers, shapes);
                this->uppers.push_back(upper);
            }
        }
    }

    virtual ~Layer() = default;
    virtual void print() const = 0;

    void print_io() const {
        std::cout << "| inputs: " << std::endl;
        for (const auto& input : this->inputs) {
            std::cout << "|   " << input.name << std::endl;
        }
        std::cout << "| outputs: " << std::endl;
        for (const auto& output : this->outputs) {
            std::cout << "|   " << output.name << std::endl;
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
    std::vector<Tensor> inputs;
    std::vector<Tensor> outputs;
    std::vector<size_t> var_dims;

    std::vector<Tensor> lowers;
    std::vector<Tensor> uppers;
};

/* Matrix multiply results from A * B */
class MatMul : public Layer {
public:
    /*
        Y = A * B
    */
    MatMul(const onnx::NodeProto& node, const Initializers& global_initializers, const TensorShapes& shapes): Layer(node, global_initializers, shapes) {
        operator_type = _matmul;
        
        this->A = this->inputs.at(0);
        this->B = this->inputs.at(1);
        if(this->B.shape.size()==2){
            this->var_dims = this->B.shape;
        }
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
    Add(const onnx::NodeProto& node, const Initializers& global_initializers, const TensorShapes& shapes): Layer(node, global_initializers, shapes) {
        operator_type = _add;
        this->A = this->inputs.at(0);
        this->B = this->inputs.at(1);
        this->var_dims = this->B.shape;
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
    GEMM(const onnx::NodeProto& node, const Initializers& global_initializers, const TensorShapes& shapes): Layer(node, global_initializers, shapes) {
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
        
        this->A = this->inputs.at(0);
        this->B = this->inputs.at(1);
        if(this->B.shape.size()==2){
            if(this->transB){
                this->B._transpose();
            }
        }
        if (this->inputs.size() == 3) {
            this->C = this->inputs.at(2);
            this->has_optional_C = true;
        }

        this->var_dims = this->B.shape;
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

class Relu : public Layer {
public:
    Relu(const onnx::NodeProto& node, const Initializers& global_initializers, const TensorShapes& shapes): Layer(node, global_initializers, shapes) {
        operator_type = _relu;
        this->X = this->inputs.at(0);
        this->is_activation_func = true;

        this->var_dims = this->X.shape;
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
