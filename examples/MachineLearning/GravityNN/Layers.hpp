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

    virtual void add_parameters(std::vector<gravity::param<>*> params) = 0;

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
        }
        if (this->inputs.size() == 3) {
            this->C = this->inputs.at(2);
            this->has_optional_C = true;
        }
    }

    void add_parameters(std::vector<gravity::param<>*> params) override {
        auto B = params[0];
        auto C = params[1];
        this->B.add_params(B, this->name);
        if (this->has_optional_C) {
            this->C.add_params(C, this->name);
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

    void add_parameters(std::vector<gravity::param<>*> params) override {}

    void print() const override{
        std::cout << "---------------------------------" << std::endl;
        std::cout << "| ReLU: " << this->name << std::endl;
        std::cout << "---------------------------------" << std::endl;
        this->print_io();
        std::cout << "---------------------------------" << std::endl;
    }

    Tensor X;
};

class Conv : public Layer {
public:
    Conv(const onnx::NodeProto& node, const Tensors& tensors): Layer(node, tensors) {
        operator_type = _conv;
        // -2 because we don't count the batch and channel dimensions
        size_t num_spatial_dims = this->inputs.at(0).shape.size() - 2;

        this->dilations = std::vector<size_t>(num_spatial_dims, 1);
        this->pads = std::vector<size_t>(num_spatial_dims*2, 0);
        this->strides = std::vector<size_t>(num_spatial_dims, 1);

        if (const auto* auto_pad_attr = find_attribute("auto_pad")) {
            this->auto_pad = auto_pad_attr->s();
            if (this->auto_pad != "NOTSET") {
                throw std::runtime_error("Conv: Only auto_pad=NOTSET is supported");
            }
        }
        
        if (const auto* group_attr = find_attribute("group")) {
            this->group = group_attr->i();
            if (this->group != 1) {
                throw std::runtime_error("Conv: Only group=1 is supported");
            }
        }

        if (const auto* dilations_attr = find_attribute("dilations")) {
            this->dilations = std::vector<size_t>(dilations_attr->ints().begin(), dilations_attr->ints().end());
            if (this->dilations.size() != 2) {
                throw std::runtime_error("Conv: Only 2D dilations is supported");
            }
            if (this->dilations[0] != 1 || this->dilations[1] != 1) {
                throw std::runtime_error("Conv: Only dilations=1 is supported");
            }
        }
        
        if (const auto* kernel_shape_attr = find_attribute("kernel_shape")) {
            this->kernel_shape = std::vector<size_t>(kernel_shape_attr->ints().begin(), kernel_shape_attr->ints().end());
            if (this->kernel_shape.size() != 2) {
                throw std::runtime_error("Conv: Only 2D kernel_shape is supported");
            }
        } else {
            throw std::runtime_error("Conv: kernel_shape attribute is required for us. If you see this error, go annoy Haydn.");
        }

        if (const auto* pads_attr = find_attribute("pads")) {
            this->pads = std::vector<size_t>(pads_attr->ints().begin(), pads_attr->ints().end());
            if (this->pads.size() != 4) {
                throw std::runtime_error("Conv: Only 4D pads is supported");
            }
        }

        if (const auto* strides_attr = find_attribute("strides")) {
            this->strides = std::vector<size_t>(strides_attr->ints().begin(), strides_attr->ints().end());
            if (this->strides.size() != 2) {
                throw std::runtime_error("Conv: Only 2D strides is supported");
            }
        }

        this->W = this->inputs.at(1);

        if (this->inputs.size() == 3) {
            this->has_optional_B = true;
            this->B = this->inputs.at(2);
        }
    }

    void print() const override{
        std::cout << "---------------------------------" << std::endl;
        std::cout << "| Conv: " << this->name << std::endl;
        std::cout << "---------------------------------" << std::endl;
        std::cout << "| auto_pad: "  << this->auto_pad << std::endl;
        std::cout << "| group: "   << this->group << std::endl;
        std::cout << "| dilations: " << print_vector(this->dilations) << std::endl;
        std::cout << "| kernel_shape: " << print_vector(this->kernel_shape) << std::endl;
        std::cout << "| pads: " << print_vector(this->pads) << std::endl;
        std::cout << "| strides: " << print_vector(this->strides) << std::endl;
        this->print_io();
        std::cout << "---------------------------------" << std::endl;
    }

    void add_parameters(std::vector<gravity::param<>*> params) override {
        auto W = params[0];
        auto B = params[1];
        this->W.add_params(W, this->name);

        if (this->has_optional_B) {
            this->B.add_params(B, this->name);
        }
    }

    std::string auto_pad = "NOTSET";
    size_t group = 1;

    std::vector<size_t> dilations;
    std::vector<size_t> kernel_shape;
    std::vector<size_t> pads;
    std::vector<size_t> strides;

    Tensor W;
    Tensor B;

    bool has_optional_B = false;

    /*
    auto_pad : string (default is NOTSET)
        auto_pad must be either NOTSET, SAME_UPPER, SAME_LOWER or VALID. Where default value is NOTSET, 
        which means explicit padding is used. SAME_UPPER or SAME_LOWER mean pad the input so that
        `output_shape[i] = ceil(input_shape[i] / strides[i])` for each axis `i`. The padding is split 
        between the two sides equally or almost equally (depending on whether it is even or odd). 
        In case the padding is an odd number, the extra padding is added at the end for SAME_UPPER 
        and at the beginning for SAME_LOWER.
    dilations : list of ints
        dilation value along each spatial axis of the filter. If not present, the dilation defaults is 
        1 along each spatial axis.
    group : int (default is 1)
        number of groups input channels and output channels are divided into.
    kernel_shape : list of ints
        The shape of the convolution kernel. If not present, should be inferred from input W.
    pads : list of ints
        Padding for the beginning and ending along each spatial axis, it can take any value greater than
        or equal to 0. The value represent the number of pixels added to the beginning and end part of 
        the corresponding axis. `pads` format should be as follow [x1_begin, x2_begin...x1_end, x2_end,...],
        where xi_begin the number of pixels added at the beginning of axis `i` and xi_end, the number of 
        pixels added at the end of axis `i`. This attribute cannot be used simultaneously with auto_pad 
        attribute. If not present, the padding defaults to 0 along start and end of each spatial axis.
    strides : list of ints
        Stride along each spatial axis. If not present, the stride defaults is 1 along each spatial axis.
    */
};