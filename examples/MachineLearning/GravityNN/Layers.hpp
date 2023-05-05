#pragma once

#include <fstream>
#include <onnx.pb.h>
#include <vector>
#include <gravity/model.h>
#include <gravity/func.h>
#include "Tensor.hpp"

typedef enum { _gemm, _relu, _conv } OType; /* Operator Type */

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

    virtual void add_parameters(std::vector<gravity::param<>*> params) const = 0;

    const onnx::AttributeProto* find_attribute(const std::string& name, const onnx::NodeProto& node) const {
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
    bool is_pre_activation = false;

    std::vector<Tensor> inputs;
    std::vector<Tensor> outputs;

    std::vector<Tensor> lowers;
    std::vector<Tensor> uppers;
};

class GEMM : public Layer {
public:
    /*
        Y = alpha * A’ * B’ + beta * C
    */
    GEMM(const onnx::NodeProto& node, const Tensors& tensors): Layer(node, tensors) {
        operator_type = _gemm;

        if (const auto* alpha_attr = find_attribute("alpha", node)) {
            this->alpha = alpha_attr->f();
        }
        if (const auto* beta_attr = find_attribute("beta", node)) {
            this->beta = beta_attr->f();
        }
        if (const auto *transA_attr = find_attribute("transA", node)) {
            this->transA = transA_attr->i();
        }
        if (const auto *transB_attr = find_attribute("transB", node)) {
            this->transB = transB_attr->i();
        }
        
        this->A = &this->inputs.at(0);
        this->B = &this->inputs.at(1);
        this->Y = &this->outputs.at(0);
        if(this->transB){
            this->B->_transpose();
        }
        if (this->inputs.size() == 3) {
            this->C = &this->inputs.at(2);
            this->has_optional_C = true;
        }

        this->in_dim = this->A->shape[1];
        this->out_dim = this->B->shape[1];
    }

    void add_parameters(std::vector<gravity::param<>*> params) const override {
        auto B = params[0];
        auto C = params[1];
        this->B->add_params(B, this->name);
        if (this->has_optional_C) {
            this->C->add_params(C, this->name);
        }
    }

    Tensor *A, *B, *C; // Inputs
    Tensor *Y; // Output

    size_t in_dim, out_dim;

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
        this->X = &this->inputs.at(0);
        this->is_activation_func = true;
    }

    void add_parameters(std::vector<gravity::param<>*> params) const override {}

    Tensor* X; // Input
};

class Conv : public Layer {
public:
    Conv(const onnx::NodeProto& node, const Tensors& tensors): Layer(node, tensors) {
        operator_type = _conv;
        // -2 because we don't count the batch and channel dimensions
        size_t num_spatial_dims = this->inputs.at(0).shape.size() - 2;
        if (num_spatial_dims != 2) {
            throw std::runtime_error("Conv: Only 2D convolutions are supported");
        }

        this->dilations = std::vector<size_t>(num_spatial_dims, 1);
        this->pads = std::vector<size_t>(num_spatial_dims*2, 0);
        this->strides = std::vector<size_t>(num_spatial_dims, 1);

        if (const auto* auto_pad_attr = find_attribute("auto_pad", node)) {
            this->auto_pad = auto_pad_attr->s();
            if (this->auto_pad != "NOTSET") {
                throw std::runtime_error("Conv: Only auto_pad=NOTSET is supported");
            }
        }
        
        if (const auto* group_attr = find_attribute("group", node)) {
            this->group = group_attr->i();
            if (this->group != 1) {
                throw std::runtime_error("Conv: Only group=1 is supported");
            }
        }

        if (const auto* dilations_attr = find_attribute("dilations", node)) {
            this->dilations = std::vector<size_t>(dilations_attr->ints().begin(), dilations_attr->ints().end());
        }
        
        if (const auto* kernel_shape_attr = find_attribute("kernel_shape", node)) {
            this->kernel_shape = std::vector<size_t>(kernel_shape_attr->ints().begin(), kernel_shape_attr->ints().end());
            if (this->kernel_shape.size() != 2) {
                throw std::runtime_error("Conv: Only 2D kernel_shape is supported");
            }
        } else {
            throw std::runtime_error("Conv: kernel_shape attribute is required for us. If you see this error, go annoy Haydn.");
        }

        if (const auto* pads_attr = find_attribute("pads", node)) {
            this->pads = std::vector<size_t>(pads_attr->ints().begin(), pads_attr->ints().end());
        }

        if (const auto* strides_attr = find_attribute("strides", node)) {
            this->strides = std::vector<size_t>(strides_attr->ints().begin(), strides_attr->ints().end());
        }

        this->X = &this->inputs.at(0);
        this->W = &this->inputs.at(1);
        this->Y = &this->outputs.at(0);

        if (this->inputs.size() == 3) {
            this->has_optional_B = true;
            this->B = &this->inputs.at(2);
        }

        this->out_c = this->outputs[0].shape[1];
        this->out_h = this->outputs[0].shape[2];
        this->out_w = this->outputs[0].shape[3];

        this->inp_c = this->inputs[0].shape[1];
        this->inp_h = this->inputs[0].shape[2];
        this->inp_w = this->inputs[0].shape[3];

        this->kern_c = this->inputs[1].shape[1];
        this->kern_h = this->inputs[1].shape[2];
        this->kern_w = this->inputs[1].shape[3];
    }

    void add_parameters(std::vector<gravity::param<>*> params) const override {
        auto W = params[0];
        auto B = params[1];
        this->W->add_params(W, this->name);

        // We have to do this because the bias is applied per "pixel"
        // and does not have the same number of elements as the output
        if (this->has_optional_B) {
            for (auto j = 0; j < this->outputs.at(0).numel; j++) {
                B->add_val(this->name+","+to_string(j), this->B->operator()(j % this->out_c));
            }
        } else {
            for (auto j = 0; j < this->outputs.at(0).numel; j++) {
                B->add_val(this->name+","+to_string(j), 0.0);
            }
        }
    }

    std::string auto_pad = "NOTSET";
    size_t group = 1;

    std::vector<size_t> dilations;
    std::vector<size_t> kernel_shape;
    std::vector<size_t> pads;
    std::vector<size_t> strides;

    Tensor *W, *B; // Weight and bias
    Tensor *X, *Y; // Input and output

    bool has_optional_B = false;

    size_t out_c, out_h, out_w;
    size_t inp_c, inp_h, inp_w;
    size_t kern_c, kern_h, kern_w;
};