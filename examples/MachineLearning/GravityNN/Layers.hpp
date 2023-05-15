#pragma once

#include <fstream>
#include <onnx.pb.h>
#include <vector>
#include <gravity/Node.h>
#include "Tensor.hpp"

using namespace gravity;

typedef enum { _gemm, _relu, _conv, _input, _noop, _add, _sub, _cos, _sin, _neg, _pow, _mul, _div} OType; /* Operator Type */

class IndexSet {
public:
    IndexSet(std::vector<string> names) {
        for (const auto& name : names) {
            this->_indices[name] = indices(name);
        }
    }

    indices& operator[](const std::string& name) {
        return this->_indices.at(name);
    }

    std::map<std::string, indices> _indices;
    size_t row_id = 0;
};

class Layer: public Node {
public:
    Layer(const onnx::NodeProto& node, Tensors& tensors): Node(node.name()) {
        for (const auto& input : node.input()) {
            if (!tensors.at(input).is_initializer) {
                this->input_names.push_back(input);
            }
            this->inputs.push_back(&tensors.at(input));
        }
        for (const auto& output : node.output()) {
            this->output_names.push_back(output);
            this->outputs.push_back(&tensors.at(output));
        }

        this->name = node.name();

        this->_load_bounds(tensors);
        this->_load_forward(tensors);
    }

    // Input layer constructor
    Layer(const onnx::ValueInfoProto& input_node, Tensors& tensors): Node(input_node.name()) {
        this->name = input_node.name();
        this->operator_type = _input;
        this->output_names.push_back(input_node.name());
        this->outputs.push_back(&tensors.at(input_node.name()));
        this->_load_bounds(tensors);
        this->_load_forward(tensors);
    }

    void _load_bounds(Tensors& tensors) {
        for (auto out: this->output_names)
        {
            auto lower_name = out + "_lower";
            if (tensors.count(lower_name) != 0) {
                this->lowers.push_back(&tensors.at(lower_name));
            }
            auto upper_name = out + "_upper";
            if (tensors.count(lower_name) != 0) {
                this->uppers.push_back(&tensors.at(upper_name));
            }
        }
    }

    void _load_forward(Tensors& tensors) {
        for (auto out: this->output_names)
        {
            auto forward_name = out + "_forward";
            if (tensors.count(forward_name) != 0) {
                this->forward_values.push_back(&tensors.at(forward_name));
            }
        }
    }

    virtual ~Layer() = default;

    virtual void add_parameters(std::vector<gravity::param<>*> params) const = 0;
    virtual void build_constraint(IndexSet& inds, std::vector<gravity::param<>*> params)= 0;

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

    std::vector<Tensor*> lowers;
    std::vector<Tensor*> uppers;
    std::vector<Tensor*> forward_values;

    // Names for inputs from other layers, NOT including initializers
    std::vector<std::string> input_names;
    std::vector<std::string> output_names;

    std::vector<Tensor*> outputs;
    std::vector<Tensor*> inputs;
};

class GEMM : public Layer {
public:
    /*
        Y = alpha * A’ * B’ + beta * C
    */
    GEMM(const onnx::NodeProto& node, Tensors& tensors): Layer(node, tensors) {
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

        this->A = &tensors.at(node.input(0));
        this->B = &tensors.at(node.input(1));
        this->Y = &tensors.at(node.output(0));

        if (node.input_size() == 3) {
            this->C = &tensors.at(node.input(2));
        }
    }

    void add_parameters(std::vector<gravity::param<>*> params) const override {
        auto B = params[0];
        auto C = params[1];

        if (this->transB) {
            Tensor tb = Tensor::transpose(*this->B);
            tb.add_params(B);
        } else {
            this->B->add_params(B);
        }

        if (this->C) {
            this->C->add_params(C);
        }
    }

    void build_constraint(IndexSet& inds, std::vector<gravity::param<>*> params) override {
        this->add_parameters(params);

        Tensor tb = *this->B;
        if (this->transB) {
            tb = Tensor::transpose(*this->B);
        }

        Tensor ta = *this->A;
        if (this->transA) {
            ta = Tensor::transpose(*this->A);
        }

        // Ensure dimensions match
        if (ta.shape[1] != tb.shape[0]) {
            throw std::runtime_error("GEMM: A and B inner dimensions do not match");
        }
        if ((ta.shape[0] != this->Y->shape[0]) || (tb.shape[1] != this->Y->shape[1])) {
            throw std::runtime_error("GEMM: A and B outer dimensions do not match Y");
        }

        for (size_t out_row = 0; out_row < this->Y->shape[0]; out_row++) {
            for (size_t out_col = 0; out_col < this->Y->shape[1]; out_col++) {
                inds["Constr"].add(this->Y->strkey(out_row, out_col));
                inds["Out"].add_ref(this->Y->strkey(out_row, out_col));

                for (size_t i = 0; i < ta.shape[1]; i++) {
                    inds["In"].add_in_row(inds.row_id, ta.strkey(out_row, i));
                    inds["B"].add_in_row(inds.row_id,   tb.strkey(i, out_col));
                }

                // Add bias
                if (this->C) {
                    inds["C"].add_in_row(inds.row_id, this->C->strkey(out_col));
                } else {
                    inds["C"].add_empty_row();
                }
                inds.row_id++;
            }
        }
    }

    Tensor *A, *B, *C = nullptr; // Inputs
    Tensor *Y; // Output

    float alpha = 1.0;
    float beta = 1.0;
    bool transA = false;
    bool transB = false;
};

class Relu : public Layer {
public:
    Relu(const onnx::NodeProto& node, Tensors& tensors): Layer(node, tensors) {
        operator_type = _relu;
        this->X = &tensors.at(node.input(0));
        this->Y = &tensors.at(node.output(0));
    }

    void add_parameters(std::vector<gravity::param<>*> params) const override {}

    void build_constraint(IndexSet& inds, std::vector<gravity::param<>*> params) override {
        for(auto j = 0; j < this->X->numel;j++){
            inds["Constr"].add(this->Y->strkey(j));
            inds["In"].add_ref(this->X->strkey(j));
            inds["Out"].add_ref(this->Y->strkey(j));
        }
    }

    Tensor* X; // Input
    Tensor* Y; // Output
};

class Conv : public Layer {
public:
    Conv(const onnx::NodeProto& node, Tensors& tensors): Layer(node, tensors) {
        operator_type = _conv;
        this->X = &tensors.at(node.input(0));
        this->W = &tensors.at(node.input(1));
        if (node.input_size() == 3) {
            this->B = &tensors.at(node.input(2));
        }

        this->Y = &tensors.at(node.output(0));

        // -2 because we don't count the batch and channel dimensions
        size_t num_spatial_dims = this->X->shape.size() - 2;
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

        this->out_c = this->Y->shape[1];
        this->out_h = this->Y->shape[2];
        this->out_w = this->Y->shape[3];

        this->inp_c = this->X->shape[1];
        this->inp_h = this->X->shape[2];
        this->inp_w = this->X->shape[3];

        this->kern_c = this->W->shape[1];
        this->kern_h = this->W->shape[2];
        this->kern_w = this->W->shape[3];
    }

    void build_constraint(IndexSet& inds, std::vector<gravity::param<>*> params) override {
        this->add_parameters(params);

        // Output indexing
        for (auto j = 0; j < this->Y->numel; j++) {
            inds["Constr"].add(this->name + "," + to_string(j));
        }

        for (int ob = 0; ob < this->Y->shape[0]; ob++) {
            for (int oh = 0; oh < this->out_h; oh++) {
                for (int ow = 0; ow < this->out_w; ow++) {
                    for (int oc = 0; oc < this->out_c; oc++) {
                        inds["Out"].add_ref(this->Y->strkey(ob, oc, oh, ow));
                        for (int kh = 0; kh < this->kern_h; kh++) {
                            for (int kw = 0; kw < this->kern_w; kw++) {
                                for (int kc = 0; kc < this->kern_c; kc++) {
                                    int h_ind = (this->strides[0]*oh + this->dilations[0]*kh - this->pads[0]);
                                    int w_ind = (this->strides[1]*ow + this->dilations[1]*kw - this->pads[3]);
                                    if ((h_ind < this->inp_h) && (h_ind >= 0) && (w_ind < this->inp_w) && (w_ind >= 0)) {
                                        inds["W"].add_in_row(inds.row_id, this->W->strkey(oc, kc, kh, kw));
                                        inds["In"].add_in_row(inds.row_id, this->X->strkey(ob, kc, h_ind, w_ind));
                                    }
                                }
                            }
                        }
                        // Add bias
                        if (this->B) {
                            inds["B"].add_in_row(inds.row_id, this->B->strkey(oc));
                        } else {
                            inds["B"].add_empty_row();
                        }
                        inds.row_id++;
                    }
                }
            }
        }

    }

    void add_parameters(std::vector<gravity::param<>*> params) const override {
        auto W = params[0];
        auto B = params[1];
        this->W->add_params(W);
        if (this->B) {
            this->B->add_params(B);
        }
    }

    std::string auto_pad = "NOTSET";
    size_t group = 1;

    std::vector<size_t> dilations;
    std::vector<size_t> kernel_shape;
    std::vector<size_t> pads;
    std::vector<size_t> strides;

    Tensor *W, *B = nullptr; // Weight and bias
    Tensor *X, *Y; // Input and output

    size_t out_c, out_h, out_w;
    size_t inp_c, inp_h, inp_w;
    size_t kern_c, kern_h, kern_w;
};

class Input : public Layer {
public:
    Input(const onnx::ValueInfoProto& node, Tensors& tensors): Layer(node, tensors) {
        operator_type = _input;
    }

    void add_parameters(std::vector<gravity::param<>*> params) const override {}
    void build_constraint(IndexSet& inds, std::vector<gravity::param<>*> params) override {}
};

class NoOp : public Layer {
public:
    NoOp(const onnx::NodeProto& node, Tensors& tensors): Layer(node, tensors) {
        operator_type = _noop;
        this->X = &tensors[node.input(0)];
        this->Y = &tensors[node.output(0)];
    }

    void add_parameters(std::vector<gravity::param<>*> params) const override {}

    void build_constraint(IndexSet& inds, std::vector<gravity::param<>*> params) override {
        std::cout << "NoOp: " << this->name << std::endl;
        for(auto j = 0; j < this->X->numel;j++){
            inds["Constr"].add(this->Y->strkey(j));
            inds["In"].add_ref(this->X->strkey(j));
            inds["Out"].add_ref(this->Y->strkey(j));
        }
    }

    Tensor *X, *Y; // Input and output
};

class Split : public NoOp {
public:
    Split(const onnx::NodeProto& node, Tensors& tensors): NoOp(node, tensors) {
        if (const auto* axis_attr = find_attribute("axis", node)) {
            int ax = axis_attr->i();
            if (ax < 0) {
                ax += this->X->shape.size();
            }
            this->axis = ax;
        }

        if (const auto* split_attr = find_attribute("split", node)) {
            this->split = std::vector<size_t>(split_attr->ints().begin(), split_attr->ints().end());
        } else {
            throw std::runtime_error("Split: split attribute not found. This is optional but I don't know how to handle it yet.");
        }
    }

    void build_constraint(IndexSet& inds, std::vector<gravity::param<>*> params) override {
        size_t cur_axis_idx = 0;
        for (auto out: this->outputs) {
            for (size_t out_idx = 0; out_idx < out->numel; out_idx++) {
                auto inp_vec_idx = out->unflatten_index(out_idx);
                inp_vec_idx[this->axis] += cur_axis_idx;
                auto inp_idx = this->X->flatten_index(inp_vec_idx);

                inds["Constr"].add(out->strkey(out_idx));
                inds["In"].add_ref(this->X->strkey(inp_idx));
                inds["Out"].add_ref(out->strkey(out_idx));
            }

            cur_axis_idx += out->shape[this->axis];
        }
    }

    size_t axis;
    std::vector<size_t> split;
};

class Concat : public NoOp {
public:
    Concat(const onnx::NodeProto& node, Tensors& tensors): NoOp(node, tensors) {
        if (const auto* axis_attr = find_attribute("axis", node)) {
            int ax = axis_attr->i();
            if (ax < 0) {
                ax += this->X->shape.size();
            }
            this->axis = ax;
        } else {
            throw std::runtime_error("Concat: axis attribute not found.");
        }
    }

    void build_constraint(IndexSet& inds, std::vector<gravity::param<>*> params) override {
        size_t cur_axis_idx = 0;
        for (auto inp: this->inputs) {
            for (size_t in_idx = 0; in_idx < inp->numel; in_idx++) {
                auto out_vec_idx = inp->unflatten_index(in_idx);
                out_vec_idx[this->axis] += cur_axis_idx;
                auto out_idx = this->Y->flatten_index(out_vec_idx);

                inds["Constr"].add(this->Y->strkey(out_idx));
                inds["In"].add_ref(inp->strkey(in_idx));
                inds["Out"].add_ref(this->Y->strkey(out_idx));
            }

            cur_axis_idx += inp->shape[this->axis];
        }
    }

    size_t axis;
};

class Add : public Layer {
public:
    Add(const onnx::NodeProto& node, Tensors& tensors): Layer(node, tensors) {
        operator_type = _add;
        this->A = &tensors[node.input(0)];
        this->B = &tensors[node.input(1)];
        this->Y = &tensors[node.output(0)];

        if ((this->A->is_initializer == true) || (this->B->is_initializer == true)) {
            throw std::runtime_error("Add: initializer not supported.");
        }
    }

    void add_parameters(std::vector<gravity::param<>*> params) const override {}

    void build_constraint(IndexSet& inds, std::vector<gravity::param<>*> params) override {
        for(auto j = 0; j < this->A->numel;j++){
            inds["Constr"].add(this->Y->strkey(j));
            inds["A"].add_ref(this->A->strkey(j));
            inds["B"].add_ref(this->B->strkey(j));
            inds["Out"].add_ref(this->Y->strkey(j));
        }
    }

    Tensor *A, *B, *Y; // Input and output
};

class Sub : public Layer {
public:
    Sub(const onnx::NodeProto& node, Tensors& tensors): Layer(node, tensors) {
        operator_type = _sub;
        this->A = &tensors[node.input(0)];
        this->B = &tensors[node.input(1)];
        this->Y = &tensors[node.output(0)];

        if ((this->A->is_initializer == true) || (this->B->is_initializer == true)) {
            throw std::runtime_error("Sub: initializer not supported.");
        }
    }

    void add_parameters(std::vector<gravity::param<>*> params) const override {}

    void build_constraint(IndexSet& inds, std::vector<gravity::param<>*> params) override {
        for(auto j = 0; j < this->A->numel;j++){
            inds["Constr"].add(this->Y->strkey(j));
            inds["A"].add_ref(this->A->strkey(j));
            inds["B"].add_ref(this->B->strkey(j));
            inds["Out"].add_ref(this->Y->strkey(j));
        }
    }

    Tensor *A, *B, *Y; // Input and output
};

class Cos : public Layer {
public:
    Cos(const onnx::NodeProto& node, Tensors& tensors): Layer(node, tensors) {
        operator_type = _cos;
        this->X = &tensors[node.input(0)];
        this->Y = &tensors[node.output(0)];
    }

    void add_parameters(std::vector<gravity::param<>*> params) const override {}

    void build_constraint(IndexSet& inds, std::vector<gravity::param<>*> params) override {
        for(auto j = 0; j < this->X->numel;j++){
            inds["Constr"].add(this->Y->strkey(j));
            inds["In"].add_ref(this->X->strkey(j));
            inds["Out"].add_ref(this->Y->strkey(j));
        }
    }

    Tensor *X, *Y; // Input and output
};

class Sin : public Layer {
public:
    Sin(const onnx::NodeProto& node, Tensors& tensors): Layer(node, tensors) {
        operator_type = _sin;
        this->X = &tensors[node.input(0)];
        this->Y = &tensors[node.output(0)];
    }

    void add_parameters(std::vector<gravity::param<>*> params) const override {}

    void build_constraint(IndexSet& inds, std::vector<gravity::param<>*> params) override {
        for(auto j = 0; j < this->X->numel;j++){
            inds["Constr"].add(this->Y->strkey(j));
            inds["In"].add_ref(this->X->strkey(j));
            inds["Out"].add_ref(this->Y->strkey(j));
        }
    }

    Tensor *X, *Y; // Input and output
};

class Neg : public Layer {
public:
    Neg(const onnx::NodeProto& node, Tensors& tensors): Layer(node, tensors) {
        operator_type = _neg;
        this->X = &tensors[node.input(0)];
        this->Y = &tensors[node.output(0)];
    }

    void add_parameters(std::vector<gravity::param<>*> params) const override {}

    void build_constraint(IndexSet& inds, std::vector<gravity::param<>*> params) override {
        for(auto j = 0; j < this->X->numel;j++){
            inds["Constr"].add(this->Y->strkey(j));
            inds["In"].add_ref(this->X->strkey(j));
            inds["Out"].add_ref(this->Y->strkey(j));
        }
    }

    Tensor *X, *Y; // Input and output
};

class Pow : public Layer {
public:
    Pow(const onnx::NodeProto& node, Tensors& tensors): Layer(node, tensors) {
        operator_type = _pow;
        this->X = &tensors[node.input(0)];
        this->exp = &tensors[node.input(1)];
        this->Y = &tensors[node.output(0)];

        if (this->exp->numel != 1) {
            throw std::runtime_error("Pow: exponent must be a scalar.");
        }

        if (this->exp->operator()(0) != 2.0) {
            throw std::runtime_error("Pow: exponent must be 2.");
        }
    }

    void add_parameters(std::vector<gravity::param<>*> params) const override {}

    void build_constraint(IndexSet& inds, std::vector<gravity::param<>*> params) override {
        for(auto j = 0; j < this->X->numel;j++){
            inds["Constr"].add(this->Y->strkey(j));
            inds["In"].add_ref(this->X->strkey(j));
            inds["Out"].add_ref(this->Y->strkey(j));
        }
    }

    Tensor *X, *exp, *Y; // Input and output
};

class Mul : public Layer {
public:
    Mul(const onnx::NodeProto& node, Tensors& tensors): Layer(node, tensors) {
        operator_type = _mul;
        this->A = &tensors[node.input(0)];
        this->B = &tensors[node.input(1)];
        this->Y = &tensors[node.output(0)];

        if ((this->A->is_initializer == true) || (this->B->is_initializer == true)) {
            throw std::runtime_error("Add: initializer not supported.");
        }
    }

    void add_parameters(std::vector<gravity::param<>*> params) const override {}

    void build_constraint(IndexSet& inds, std::vector<gravity::param<>*> params) override {
        for(auto j = 0; j < this->A->numel;j++){
            inds["Constr"].add(this->Y->strkey(j));
            inds["A"].add_ref(this->A->strkey(j));
            inds["B"].add_ref(this->B->strkey(j));
            inds["Out"].add_ref(this->Y->strkey(j));
        }
    }

    Tensor *A, *B, *Y; // Input and output
};

class Transpose : public NoOp {
public:
    Transpose(const onnx::NodeProto& node, Tensors& tensors): NoOp(node, tensors) {
        this->X = &tensors[node.input(0)];
        this->Y = &tensors[node.output(0)];

        if (const auto* attr = this->find_attribute("perm", node)) {
            this->perm = std::vector<size_t>(attr->ints().begin(), attr->ints().end());
        } else {
            throw std::runtime_error("Transpose: perm attribute not found.");
        }
    }

    void build_constraint(IndexSet& inds, std::vector<gravity::param<>*> params) override {
        // Copy this->X
        Tensor trx = *this->X;
        trx.shape = apply_permutation(this->X->shape, this->perm);

        for(auto j = 0; j < this->X->numel;j++){
            auto xunflat = apply_permutation(trx.unflatten_index(j), this->perm);
            auto xindex = this->X->flatten_index(xunflat);

            inds["Constr"].add(this->Y->strkey(j));
            inds["In"].add_ref(this->X->strkey(xindex));
            inds["Out"].add_ref(this->Y->strkey(j));
        }
    }

    Tensor *X, *Y; // Input and output
    std::vector<size_t> perm;
};

class Slice: public NoOp {
public:
    Slice(const onnx::NodeProto& node, Tensors& tensors): NoOp(node, tensors) {
        // Required inputs
        this->X = &tensors[node.input(0)];
        this->Y = &tensors[node.output(0)];

        this->starts = tensors[node.input(1)].get_int_data();
        this->ends = tensors[node.input(2)].get_int_data();

        if (node.input_size() == 5) {
            this->axes  = tensors[node.input(3)].get_int_data();
            this->steps = tensors[node.input(4)].get_int_data();
        } else {
            throw std::runtime_error("Slice: optional inputs not supported.");
        }

        /*
        All negative elements of axes are made non-negatve by adding r to them, where r =rank(input).
        All negative values in starts[i] and ends[i] have dims[axes[i]] added to them, where dims are the dimensions of input.
        */

        for (auto& v : this->axes) {
            if (v < 0) {
                v += this->X->ndims;
            }
        }

        for (size_t i = 0; i < this->starts.size(); i++) {
            if (this->starts[i] < 0) {
                this->starts[i] += this->X->shape.at(this->axes.at(i)) + 1;
            }
            if (this->ends[i] < 0) {
                this->ends[i] += this->X->shape.at(this->axes.at(i)) + 1;
            }
            this->ends[i] = std::min(this->ends[i], (int64_t)this->X->shape.at(this->axes.at(i)));
            this->ends[i] = std::max(this->ends[i], (int64_t)0);
        }

        std::cout << "Slice: " << this->starts.size() << " " << this->ends.size() << " " << this->axes.size() << " " << this->steps.size() << std::endl;
    }

    void build_constraint(IndexSet& inds, std::vector<gravity::param<>*> params) override {
        /*
            Slice uses the starts, ends, axes and steps inputs to select a sub-tensor of its input data tensor.
            An effective start[i], end[i], and step[i] must be computed for each i in [0, ... r-1] where r = rank(input) as follows:
            If axes are omitted, they are set to [0, ..., r-1]. If steps are omitted, they are set to [1, ..., 1] of length len(starts)
            The effective values are initialized as start[i] = 0, end[i] = dims[i] where dims are the dimensions of input and step[i] = 1.
            All negative elements of axes are made non-negatve by adding r to them, where r =rank(input).
            All negative values in starts[i] and ends[i] have dims[axes[i]] added to them, where dims are the dimensions of input. Then start[axes[i]] is the adjusted starts[i] is clamped into the range [0, dims[axes[i]]] for positive stepping and [0, dims[axes[i]]-1] for negative stepping.
            The clamping for the adjusted ends[i] depends on the sign of steps[i] and must accommodate copying 0 through dims[axes[i]] elements, so for positive stepping end[axes[i]] is clamped to [0, dims[axes[i]]], while for negative stepping it is clamped to [-1, dims[axes[i]]-1].
            Finally, step[axes[i]] = steps[i].
            For slicing to the end of a dimension with unknown size, it is recommended to pass in INT_MAX when slicing forward and ‘INT_MIN’ when slicing backward.
        */

        std::map<int64_t, int64_t> eff_start;
        std::map<int64_t, int64_t> eff_end;
        std::map<int64_t, int64_t> eff_step;
        for (auto i = 0; i < this->X->ndims; i++) {
            eff_start[i] = 0;
            eff_end[i] = this->X->shape[i];
            eff_step[i] = 1;
        }
        for (auto i = 0; i < this->starts.size(); i++) {
            eff_start[this->axes[i]] = this->starts[i];
            eff_end[this->axes[i]] = this->ends[i];
            eff_step[this->axes[i]] = this->steps[i];
        }

        if (this->X->ndims == 2) {
            size_t outr = 0;
            for (auto r = eff_start[0]; r < eff_end[0]; r += eff_step[0]) {
                size_t outc = 0;
                for (auto c = eff_start[1]; c < eff_end[1]; c += eff_step[1]) {
                    inds["Constr"].add(this->Y->strkey(outr, outc));
                    inds["In"].add_ref(this->X->strkey(r, c));
                    inds["Out"].add_ref(this->Y->strkey(outr, outc));
                    outc += 1;
                }
                outr += 1;
            }
        } else {
            throw std::runtime_error("Slice: ndims > 2 not supported.");
        }
    }

    Tensor *X, *Y;
    std::vector<int64_t> starts;
    std::vector<int64_t> ends;
    std::vector<int64_t> axes;
    std::vector<int64_t> steps;
};

class Div : public Layer {
public:
    Div(const onnx::NodeProto& node, Tensors& tensors): Layer(node, tensors) {
        operator_type = _div;
        this->A = &tensors[node.input(0)];
        this->B = &tensors[node.input(1)];
        this->Y = &tensors[node.output(0)];

        if ((this->A->is_initializer == true) || (this->B->is_initializer == true)) {
            throw std::runtime_error("Add: initializer not supported.");
        }
    }

    void add_parameters(std::vector<gravity::param<>*> params) const override {}

    void build_constraint(IndexSet& inds, std::vector<gravity::param<>*> params) override {
        for(auto j = 0; j < this->A->numel;j++){
            inds["Constr"].add(this->Y->strkey(j));
            inds["A"].add_ref(this->A->strkey(j));
            inds["B"].add_ref(this->B->strkey(j));
            inds["Out"].add_ref(this->Y->strkey(j));
        }
    }

    Tensor *A, *B, *Y; // Input and output
};

class Gather : public NoOp {
public:
    Gather(const onnx::NodeProto& node, Tensors& tensors): NoOp(node, tensors) {
        this->X = &tensors[node.input(0)];
        this->Y = &tensors[node.output(0)];
        this->indices = &tensors[node.input(1)];

        if (const auto* axis = find_attribute("axis", node)) {
            int64_t i = axis->i();
            if (i < 0) {
                i += this->X->ndims;
            }
            this->axis = i;
        }

    }

    void build_constraint(IndexSet& inds, std::vector<gravity::param<>*> params) override {
        /*
        Given data tensor of rank r >= 1, and indices tensor of rank q, gather entries of the axis dimension
        of data (by default outer-most one as axis=0) indexed by indices, and concatenates them in an output
        tensor of rank q + (r - 1).

        axis = 0:
            Let k = indices[i_{0}, …, i_{q-1}]  (Rank k = q)
            output[i_{0}, …, i_{q-1}, j_{0}, …, j_{r-2}] = input[k , j_{0}, …, j_{r-2}]
        axis = 1:
            Let k = indices[i_{0}, …, i_{q-1}]
            output[j_{0}, i_{0}, …, i_{q-1}, j_{1}, …, j_{r-2}] = input[j_{0}, k, j_{1}, …, j_{r-2}]
        */

        size_t r = this->X->ndims;
        size_t q = this->indices->ndims;
        // We may have to calculate q manually because we 
        // do not handle rank-0 tensors correctly
        if (this->X->ndims != this->Y->ndims) {
            q = this->Y->ndims - (r - 1);
            std::cout << "Contraction!" << std::endl;
        }

        for (auto o = 0; o < this->Y->numel; o++) {
            auto unflat = this->Y->unflatten_index(o);
            std::vector<size_t> i(unflat.begin() + this->axis, unflat.begin() + this->axis + q);
            std::vector<size_t> j;
            for (auto _i = 0; _i < unflat.size(); _i++) {
                if ((_i < this->axis) || (_i >= this->axis + q)) {
                    j.push_back(unflat[_i]);
                }
            }

            if (i.size() == 0) {
                i.push_back(0);
            }
            size_t k = this->indices->get_int_data().at(this->indices->flatten_index(i));
            // insert k at position this->axis in j
            j.insert(j.begin() + this->axis, k);

            auto in_index = this->X->flatten_index(j);
            auto out_index = this->Y->flatten_index(unflat);
            inds["Constr"].add(this->Y->strkey(out_index));
            inds["In"].add_ref(this->X->strkey(in_index));
            inds["Out"].add_ref(this->Y->strkey(out_index));
        }
    }

    Tensor *X, *Y; // Input and output
    Tensor* indices;
    size_t axis = 0;
};

class Clip : public NoOp {
public:
    Clip(const onnx::NodeProto& node, Tensors& tensors): NoOp(node, tensors) {
        this->X = &tensors.at(node.input(0));
        this->Y = &tensors.at(node.output(0));

        if (this->inputs.size() > 1) {
            this->min = tensors.at(node.input(1))(0);
        }
        if (this->inputs.size() > 2) {
            this->max = tensors.at(node.input(2))(0);
        }
    }

    void build_constraint(IndexSet& inds, std::vector<gravity::param<>*> params) override {
        // this->add_parameters(params);
        for(auto j = 0; j < this->X->numel;j++){
            inds["Constr"].add(this->Y->strkey(j));
            inds["In"].add_ref(this->X->strkey(j));
            inds["Out"].add_ref(this->Y->strkey(j));
        }
    }

    Tensor* X; // Input
    Tensor* Y; // Output

    float min = std::numeric_limits<float>::lowest();
    float max = std::numeric_limits<float>::max();
};