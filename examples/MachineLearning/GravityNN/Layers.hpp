#pragma once

#include <fstream>
#include <onnx.pb.h>
#include <vector>
#include <gravity/model.h>
#include <gravity/func.h>


typedef enum { _gemm, _relu, _conv, _matmul, _add} OType; /* Operator Type */

class Tensor {
public:
    Tensor() {}

    Tensor(const onnx::TensorProto& tensor) {
        this->name = tensor.name();
        this->is_initializer = true;

        this->shape = std::vector<size_t>(tensor.dims().begin(), tensor.dims().end());
        this->numel = std::accumulate(this->shape.begin(), this->shape.end(), 1, std::multiplies<size_t>());
        if (!tensor.raw_data().empty()) {
            const void* raw_data = tensor.raw_data().data();
            data.resize(this->numel);
            std::memcpy(data.data(), raw_data, this->numel * sizeof(float));
        } else if (!tensor.float_data().empty()) {
            this->data = std::vector<float>(tensor.float_data().begin(), tensor.float_data().end());
        } else {
            throw std::runtime_error("Tensor " + tensor.name() + " has data in neither raw_data nor float_data.");
        }
    }

    Tensor(const onnx::ValueInfoProto& vinfo) {
        this->name = vinfo.name();
        this->is_initializer = false;

        this->shape.clear();
        for (auto dim : vinfo.type().tensor_type().shape().dim()) {
            this->shape.push_back(dim.dim_value());
        }
        this->numel = std::accumulate(this->shape.begin(), this->shape.end(), 1, std::multiplies<size_t>());
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

    float operator()(size_t i) const {
        if (!this->is_initializer) {
            throw std::runtime_error("Reading from non-initializer tensor. Perhaps you're assuming this tensor is a weight when it's actually an output of a previous layer?");
        }
        return this->data.at(i);
    }

    std::string name;
    std::vector<size_t> shape;
    size_t numel;
    bool is_initializer;

private:
    std::vector<float> data;
};
typedef std::map<std::string, Tensor> Tensors;

// Define the base Layer class
class Layer {
public:
    Layer(const onnx::NodeProto& node, const Tensors& tensors) {
        for (const auto& input : node.input()) {
            this->inputs.push_back(tensors.at(input));
        }
        for (const auto& output : node.output()) {
            this->outputs.push_back(tensors.at(output));
        }
        this->name = node.name();

        // Try to find upper and lower bounds for each output
        for (auto out: this->outputs)
        {
            auto lower_name = out.name + "_lower";
            if (tensors.count(lower_name) != 0) {
                this->lowers.push_back(tensors.at(lower_name));
            }
            auto upper_name = out.name + "_upper";
            if (tensors.count(lower_name) != 0) {
                this->uppers.push_back(tensors.at(upper_name));
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
    MatMul(const onnx::NodeProto& node, const Tensors& tensors): Layer(node, tensors) {
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
    Add(const onnx::NodeProto& node, const Tensors& tensors): Layer(node, tensors) {
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
    GEMM(const onnx::NodeProto& node, const Tensors& tensors): Layer(node, tensors) {
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
    Relu(const onnx::NodeProto& node, const Tensors& tensors): Layer(node, tensors) {
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
