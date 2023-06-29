#pragma once

#include <network/Layers/LayerBase.hpp>
#include <utils/shape_iterator.hpp>

using namespace gravity;

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
            if (!this->C->is_initializer) {
                throw std::runtime_error("GEMM: C must be an initializer");
            }
        }

        if (this->A->is_initializer) {
            throw std::runtime_error("GEMM: A cannot be an initializer");
        }

        if (!this->B->is_initializer) {
            throw std::runtime_error("GEMM: B must be an initializer");
        }
    }

    std::vector<std::vector<std::string>> get_indices() const override {
        return {{"In", "Out"}, {"B", "C"}};
    }

    void add_parameters(gravity::param<>& w) const override {
        if (this->transB) {
            Tensor tb = Tensor::transpose(*this->B);
            tb.add_params(w);
        } else {
            this->B->add_params(w);
        }

        if (this->C) {
            this->C->add_params(w);
        }
    }

    void index_constraint(IndexSet& inds) override {
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
                    inds["B"].add_in_row(inds.row_id,  tb.strkey(i, out_col));
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

    void add_constraints(gravity::Model<>& NN, IndexSet& inds, gravity::param<>& w, gravity::var<>& x, gravity::var<int>& y) override {
        Constraint<> Gemm(this->lname() + "_Gemm");
        Gemm = x.in(inds["Out"]) - (x.in(inds["In"])*w.in(inds["B"]) + w.in(inds["C"]));
        NN.add(Gemm.in(inds["Constr"]) == 0);
    }

    Tensor *A, *B, *C = nullptr; // Inputs
    Tensor *Y; // Output

    double alpha = 1.0;
    double beta = 1.0;
    bool transA = false;
    bool transB = false;
};

class MatMul : public Layer {
public:
    /*
        Y = A * B
    */
    MatMul(const onnx::NodeProto& node, Tensors& tensors): Layer(node, tensors) {
        operator_type = _matmul;

        this->A = &tensors.at(node.input(0));
        this->B = &tensors.at(node.input(1));
        this->Y = &tensors.at(node.output(0));

        this->rdim = this->Y->shape.size() - 2;
        this->cdim = this->Y->shape.size() - 1;

        if (this->A->is_initializer && this->B->is_initializer) {
            throw std::runtime_error("MatMul: A and B are both initializers. Why wasn't this optimized out?");
        }

        // Either B must be two dimensional or must have the exact same shape as A
        // This is simply because we dont support broadcasting over dims of B yet,
        // Only A
        if ((this->B->ndims != 2)) {
            // Leading dimensions of B must match A
            for (size_t i = 0; i < this->A->shape.size() - 2; i++) {
                if (this->A->shape.at(i) != this->B->shape.at(i)) {
                    throw std::runtime_error("MatMul: B must be two dimensional or have the same shape as A");
                }
            }
        }
    }

    std::vector<std::vector<std::string>> get_indices() const override {
        return {{"pIn", "pOut", "hLeft", "hRight", "hOut"}, {"Param"}};
    }

    void add_parameters(gravity::param<>& w) const override {
        if (this->A->is_initializer) {
            this->A->add_params(w);
        } 

        if (this->B->is_initializer) {
            this->B->add_params(w);
        }
    }

    void _index_constraint(std::vector<size_t> lmat_broad_ind, indices& lmat, indices& rmat, indices& out, indices& constr, size_t& row_id) {
        for (size_t out_row = 0; out_row < this->Y->shape.at(this->rdim); out_row++) {
            for (size_t out_col = 0; out_col < this->Y->shape.at(this->cdim); out_col++) {
                auto y_ind = concat(lmat_broad_ind, {out_row, out_col});

                constr.add_in_row(row_id, this->Y->strkey(y_ind));
                out.add_in_row(row_id, this->Y->strkey(y_ind));

                for (size_t i = 0; i < this->A->shape.at(this->cdim); i++) {
                    auto a_ind = concat(lmat_broad_ind, {out_row, i});
                    auto b_ind = std::vector<size_t>{i, out_col};
                    if (this->B->ndims != 2) {
                        b_ind = concat(lmat_broad_ind, b_ind);
                    }
                
                    lmat.add_in_row(row_id, this->A->strkey(a_ind));
                    rmat.add_in_row(row_id, this->B->strkey(b_ind));
                }
                row_id++;
            }
        }
    }

    void index_constraint(IndexSet& inds) override {
        std::vector<size_t> lmat_broad_ind;
        for (size_t i = 0; i < this->A->shape.size() - 2; i++) {
            lmat_broad_ind.push_back(this->A->shape.at(i));
        }

        for (auto broad_ind: ShapeIter(lmat_broad_ind)) {
            if (this->A->is_initializer) {
                this->_index_constraint(broad_ind, inds["Param"], inds["pIn"], inds["pOut"], inds["Constr"], inds.row_id);
            } else if (this->B->is_initializer) {                                                                        
                this->_index_constraint(broad_ind, inds["pIn"], inds["Param"], inds["pOut"], inds["Constr"], inds.row_id);
            } else {
                this->_index_constraint(broad_ind, inds["hLeft"], inds["hRight"], inds["hOut"], inds["ConstrB"], inds.row_id2);
            }
        }
    }

    void add_constraints(gravity::Model<>& NN, IndexSet& inds, gravity::param<>& w, gravity::var<>& x, gravity::var<int>& y) override {
        Constraint<> MatMul_P(this->lname() + "_MatMul_Param");
        MatMul_P = x.in(inds["pOut"]) - (w.in(inds["Param"])*x.in(inds["pIn"]));
        NN.add(MatMul_P.in(inds["Constr"]) == 0);

        Constraint<> MatMul_H(this->lname() + "_MatMul_Hidden");
        MatMul_H = x.in(inds["hOut"]) - (x.in(inds["hLeft"])*x.in(inds["hRight"]));
        NN.add(MatMul_H.in(inds["ConstrB"]) == 0);
    }

    Tensor *A, *B; // Inputs
    Tensor *Y; // Output

    size_t rdim, cdim;
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
        this->conv_dim = this->X->shape.size() - 2;
        if (this->conv_dim != 2 && this->conv_dim != 1) {
            throw std::runtime_error("Conv: Only 1D/2D convolutions are supported");
        }

        this->dilations = std::vector<size_t>(this->conv_dim, 1);
        this->pads = std::vector<size_t>(this->conv_dim*2, 0);
        this->strides = std::vector<size_t>(this->conv_dim, 1);

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
        this->inp_c = this->X->shape[1];
        this->kern_c = this->W->shape[1];
        if (this->conv_dim == 2) {
            this->out_h = this->Y->shape[2];
            this->out_w = this->Y->shape[3];

            this->inp_h = this->X->shape[2];
            this->inp_w = this->X->shape[3];

            this->kern_h = this->W->shape[2];
            this->kern_w = this->W->shape[3];
        } else {
            this->out_w = this->Y->shape[2];
            this->inp_w = this->X->shape[2];
            this->kern_w = this->W->shape[2];
        }

        if (this->X->is_initializer) {
            throw std::runtime_error("Conv: X cannot be an initializer");
        }

        if (!this->W->is_initializer) {
            throw std::runtime_error("Conv: W must be an initializer");
        }

        if ((this->B != nullptr) && (!this->B->is_initializer)) {
            throw std::runtime_error("Conv: B must be an initializer");
        }
    }

    std::vector<std::vector<std::string>> get_indices() const override {
        return {{"In", "Out"}, {"W", "B"}};
    }

    void index_constraint(IndexSet& inds) override {
        if (this->conv_dim == 2) {
            this->index_constraint_2d(inds);
        } else {
            this->index_constraint_1d(inds);
        }
    }

    void index_constraint_1d(IndexSet& inds) {
        // Output indexing
        for (auto j = 0; j < this->Y->numel; j++) {
            inds["Constr"].add(this->name + "," + to_string(j));
        }

        for (int ob = 0; ob < this->Y->shape[0]; ob++) {
            for (int ow = 0; ow < this->out_w; ow++) {
                for (int oc = 0; oc < this->out_c; oc++) {
                    inds["Out"].add_ref(this->Y->strkey(ob, oc, ow));
                    for (int kw = 0; kw < this->kern_w; kw++) {
                        for (int kc = 0; kc < this->kern_c; kc++) {
                            int w_ind = (this->strides[0]*ow + this->dilations[0]*kw - this->pads[0]);
                            if ((w_ind < this->inp_w) && (w_ind >= 0)) {
                                inds["W"].add_in_row(inds.row_id, this->W->strkey(oc, kc, kw));
                                inds["In"].add_in_row(inds.row_id, this->X->strkey(ob, kc, w_ind));
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

    void index_constraint_2d(IndexSet& inds) {
        auto& oInd = inds["Out"];
        auto& wInd = inds["W"];
        auto& iInd = inds["In"];
        auto& bInd = inds["B"];
        auto& cInd = inds["Constr"];

        for (int ob = 0; ob < this->Y->shape[0]; ob++) {
            for (int oh = 0; oh < this->out_h; oh++) {
                for (int ow = 0; ow < this->out_w; ow++) {
                    for (int oc = 0; oc < this->out_c; oc++) {
                        auto okey = this->Y->strkey(ob, oc, oh, ow);
                        oInd.add_ref(okey);
                        cInd.add(okey);
                        for (int kh = 0; kh < this->kern_h; kh++) {
                            for (int kw = 0; kw < this->kern_w; kw++) {
                                for (int kc = 0; kc < this->kern_c; kc++) {
                                    int h_ind = (this->strides[0]*oh + this->dilations[0]*kh - this->pads[0]);
                                    int w_ind = (this->strides[1]*ow + this->dilations[1]*kw - this->pads[3]);
                                    if ((h_ind < this->inp_h) && (h_ind >= 0) && (w_ind < this->inp_w) && (w_ind >= 0)) {
                                        wInd.add_in_row(inds.row_id, this->W->strkey(oc, kc, kh, kw));
                                        iInd.add_in_row(inds.row_id, this->X->strkey(ob, kc, h_ind, w_ind));
                                    }
                                }
                            }
                        }
                        // Add bias
                        if (this->B) {
                            bInd.add_in_row(inds.row_id, this->B->strkey(oc));
                        } else {
                            bInd.add_empty_row();
                        }
                        inds.row_id++;
                    }
                }
            }
        }
    }

    void add_parameters(gravity::param<>& w) const override {
        this->W->add_params(w);
        if (this->B != nullptr) {
            this->B->add_params(w);
        }
    }

    void add_constraints(gravity::Model<>& NN, IndexSet& inds, gravity::param<>& w, gravity::var<>& x, gravity::var<int>& y) override {
        Constraint<> Conv_(this->lname() + "_Conv");
        Conv_ = x.in(inds["Out"]) - (x.in(inds["In"])*w.in(inds["W"]) + w.in(inds["B"]));
        NN.add(Conv_.in(inds["Constr"]) == 0);
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

    // 1D/2D
    size_t conv_dim;
};

class BatchNorm: public Layer {
public:
    BatchNorm(const onnx::NodeProto& node, Tensors& tensors): Layer(node, tensors) {
        operator_type = _batchnorm;
        this->X = &tensors[node.input(0)];

        this->scale = &tensors[node.input(1)];
        this->B = &tensors[node.input(2)];
        this->mean = &tensors[node.input(3)];
        this->var = &tensors[node.input(4)];

        this->Y = &tensors[node.output(0)];

        // find epsilon
        if (const auto* attr = this->find_attribute("epsilon", node)) {
            this->epsilon = attr->f();
        }
    }

    std::vector<std::vector<std::string>> get_indices() const override {
        return {{"Out", "In"}, {"scale", "B", "mean", "var"}};
    }

    void index_hidden_states(indices& hidden_states, indices& y_ids) override {
        for (auto output: this->outputs) {
            for (auto i = 0; i < output->numel; i++) {
                hidden_states.add(output->strkey(i));
            }
        }
    }

    void index_constraint(IndexSet& inds) override {
        for (auto outer_ind: ShapeIter(this->Y->shape)) {
            size_t inner_ind = outer_ind.at(1);

            inds["Constr"].add(this->Y->strkey(outer_ind));
            inds["In"].add_ref(this->X->strkey(outer_ind));
            inds["Out"].add_ref(this->Y->strkey(outer_ind));

            inds["scale"].add_ref(this->scale->strkey(inner_ind));
            inds["B"].add_ref(this->B->strkey(inner_ind));
            inds["mean"].add_ref(this->mean->strkey(inner_ind));
            inds["var"].add_ref(this->var->strkey(inner_ind));
        }
    }

    void add_parameters(gravity::param<>& w) const override {
        this->scale->add_params(w);
        this->B->add_params(w);
        this->mean->add_params(w);
        this->var->add_params(w);
    }

    void add_constraints(gravity::Model<>& NN, IndexSet& inds, gravity::param<>& w, gravity::var<>& x, gravity::var<int>& y) override {
        auto Y = x.in(inds["Out"]);
        auto X = x.in(inds["In"]);
        auto var = w.in(inds["var"]);
        auto mean = w.in(inds["mean"]);
        auto B = w.in(inds["B"]);
        auto scale = w.in(inds["scale"]);

        Constraint<> BatchNorm(this->lname() + "_BatchNorm");
        BatchNorm = (Y-B)*(1.0/scale)*gravity::sqrt(var + this->epsilon) + mean - X;
        NN.add(BatchNorm.in(inds["Constr"]) == 0);
    }

    Tensor *X, *Y; // Input and output
    Tensor *scale, *B, *mean, *var;
    double epsilon = 1e-05;
};