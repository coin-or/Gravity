#pragma once

#include <fstream>
#include <onnx.pb.h>
#include <vector>
#include <gravity/model.h>
#include <gravity/func.h>
#include "Tensor.hpp"

typedef enum { _gemm, _relu } OType; /* Operator Type */

class Layer {
public:
    Layer(const onnx::NodeProto& node, const Tensors& tensors) {
        this->_node_proto = node;

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

    const onnx::AttributeProto* find_attribute(const std::string& name) const {
        for (const auto& attr : this->_node_proto.attribute()) {
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
    bool is_pre_activation = false;

    std::vector<Tensor> inputs;
    std::vector<Tensor> outputs;

    std::vector<Tensor> lowers;
    std::vector<Tensor> uppers;

    onnx::NodeProto _node_proto;
};

class GEMM : public Layer {
public:
    /*
        Y = alpha * A’ * B’ + beta * C
    */
    GEMM(const onnx::NodeProto& node, const Tensors& tensors): Layer(node, tensors) {
        operator_type = _gemm;

        const auto* alpha_attr = find_attribute("alpha");
        if (alpha_attr) {
            this->alpha = alpha_attr->f();
        }
        const auto* beta_attr = find_attribute("beta");
        if (beta_attr) {
            this->beta = beta_attr->f();
        }
        const auto *transA_attr = find_attribute("transA");
        if (transA_attr) {
            this->transA = transA_attr->i();
        }
        const auto *transB_attr = find_attribute("transB");
        if (transB_attr) {
            this->transB = transB_attr->i();
        }
        
        this->A = this->inputs.at(0);
        this->B = this->inputs.at(1);
        if(this->transB){
            this->B._transpose();
//            this->inputs[1]._transpose();
        }
        if (this->inputs.size() == 3) {
            this->C = this->inputs.at(2);
            this->has_optional_C = true;
        }
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
