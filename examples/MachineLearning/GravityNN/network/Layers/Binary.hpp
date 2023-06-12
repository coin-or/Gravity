#pragma once

#include <network/Layers/LayerBase.hpp>

using namespace gravity;

std::vector<bool> get_matching_dims(const Tensor* A, const Tensor* B) {
    auto a_dims = A->shape;
    auto b_dims = B->shape;
    // If one is smaller than the other, prepend with 1s
    if (a_dims.size() < b_dims.size()) {
        a_dims.insert(a_dims.begin(), b_dims.size() - a_dims.size(), 1);
    } else if (b_dims.size() < a_dims.size()) {
        b_dims.insert(b_dims.begin(), a_dims.size() - b_dims.size(), 1);
    }
    // Check that the dimensions match
    if (a_dims.size() != b_dims.size()) {
        throw std::runtime_error("get_matching_dims: dimensions do not match.");
    }

    std::vector<bool> matching_dims;
    for (auto i = 0; i < a_dims.size(); i++) {
        matching_dims.push_back(a_dims[i] == b_dims[i]);
    }

    return matching_dims;
}

class UnaryBase: public Layer {
public:
    UnaryBase(const onnx::NodeProto& node, Tensors& tensors, OType op): Layer(node, tensors) {
        this->operator_type = op;

        this->A = &tensors[node.input(0)];
        this->B = &tensors[node.input(1)];
        this->Y = &tensors[node.output(0)];

        if ((this->A->is_initializer == true) && (this->B->is_initializer == true)) {
            throw std::runtime_error("Mul: both args being initializer not supported.");
        }

        // Make sure that B is the initializer
        if (this->A->is_initializer == true) {
            std::swap(this->A, this->B);
        }

        this->broadcast_dims = get_matching_dims(this->A, this->B);
    }

    void add_parameters(gravity::param<>& w) const override {
        if (this->B->is_initializer) {
            this->B->add_params(w);
        }
    }

    std::vector<std::vector<std::string>> get_indices() const override {
        return {{"hOut", "pOut", "A", "B", "pA"}, {"UnaryVar"}};
    }

    // Returns A index and B index taking into account broadcasting
    std::pair<size_t, size_t> broadcast_index(size_t j) {
        // If A has fewer elements than B, then we need to broadcast A
        size_t A_ind = j;
        size_t B_ind = j;
        if (this->A->numel < this->B->numel) {
            auto B_ind_unflat = this->B->unflatten_index(j);
            auto A_ind_unflat = this->B->unflatten_index(j);
            for (auto i = 0; i < this->broadcast_dims.size(); i++) {
                if (this->broadcast_dims[i] == false) {
                    A_ind_unflat[i] = 0;
                }
            }
            A_ind = this->A->flatten_index(A_ind_unflat);
            B_ind = this->B->flatten_index(B_ind_unflat);
        } else {
            auto B_ind_unflat = this->A->unflatten_index(j);
            auto A_ind_unflat = this->A->unflatten_index(j);
            for (auto i = 0; i < this->broadcast_dims.size(); i++) {
                if (this->broadcast_dims[i] == false) {
                    B_ind_unflat[i] = 0;
                }
            }
            A_ind = this->A->flatten_index(A_ind_unflat);
            B_ind = this->B->flatten_index(B_ind_unflat);
        }

        return {A_ind, B_ind};
    }

    void index_constraint(IndexSet& inds) override {
        for(size_t j = 0; j < this->Y->numel;j++){
            auto broad_inds = broadcast_index(j);
            auto A_ind = broad_inds.first;
            auto B_ind = broad_inds.second;
            if (this->B->is_initializer) {
                inds["ConstrB"].add(this->Y->strkey(j));
                inds["pA"].add_ref(this->A->strkey(A_ind));
                inds["UnaryVar"].add_ref(this->B->strkey(B_ind));
                inds["pOut"].add_ref(this->Y->strkey(j));
            } else {
                inds["Constr"].add(this->Y->strkey(j));
                inds["A"].add_ref(this->A->strkey(A_ind));
                inds["B"].add_ref(this->B->strkey(B_ind));
                inds["hOut"].add_ref(this->Y->strkey(j));
            }
        }
    }

    Tensor *A, *B, *Y; // Input and output
    std::vector<bool> broadcast_dims;
};

class Add: public UnaryBase {
public:
    Add(const onnx::NodeProto& node, Tensors& tensors): UnaryBase(node, tensors, _add) {}

    void add_constraints(gravity::Model<>& NN, IndexSet& inds, gravity::param<>& w, gravity::var<>& x, gravity::var<int>& y) override {
        // Constraint for addition of hidden states
        Constraint<> Add_("Add");
        Add_ = x.in(inds["hOut"]) - (x.in(inds["A"]) + x.in(inds["B"]));
        NN.add(Add_.in(inds["Constr"]) == 0);

        // Constraint where B is a parameter
        Constraint<> Add_Param_("Add_Param");
        Add_Param_ = x.in(inds["pOut"]) - (x.in(inds["pA"]) + w.in(inds["UnaryVar"]));
        NN.add(Add_Param_.in(inds["ConstrB"]) == 0);
    }
};

class Sub: public UnaryBase {
public:
    Sub(const onnx::NodeProto& node, Tensors& tensors): UnaryBase(node, tensors, _sub) {}

    void add_constraints(gravity::Model<>& NN, IndexSet& inds, gravity::param<>& w, gravity::var<>& x, gravity::var<int>& y) override {
        // Constraint for addition of hidden states
        Constraint<> Sub_("Sub");
        Sub_ = x.in(inds["hOut"]) - (x.in(inds["A"]) - x.in(inds["B"]));
        NN.add(Sub_.in(inds["Constr"]) == 0);

        // Constraint where B is a parameter
        Constraint<> Sub_Param_("Sub_Param");
        Sub_Param_ = x.in(inds["pOut"]) - (x.in(inds["pA"]) - w.in(inds["UnaryVar"]));
        NN.add(Sub_Param_.in(inds["ConstrB"]) == 0);
    }
};

class Mul: public UnaryBase {
public:
    Mul(const onnx::NodeProto& node, Tensors& tensors): UnaryBase(node, tensors, _mul) {}

    void add_constraints(gravity::Model<>& NN, IndexSet& inds, gravity::param<>& w, gravity::var<>& x, gravity::var<int>& y) override {
        // Constraint for mul of hidden states
        Constraint<> Mul_("Mul");
        Mul_ = x.in(inds["hOut"]) - (x.in(inds["A"]) * x.in(inds["B"]));
        NN.add(Mul_.in(inds["Constr"]) == 0);

        // Constraint where B is a parameter
        Constraint<> Mul_Param_("Mul_Param");
        Mul_Param_ = x.in(inds["pOut"]) - (x.in(inds["pA"]) * w.in(inds["UnaryVar"]));
        NN.add(Mul_Param_.in(inds["ConstrB"]) == 0);
    }
};

class Div: public UnaryBase {
public:
    Div(const onnx::NodeProto& node, Tensors& tensors): UnaryBase(node, tensors, _div) {}

    void add_constraints(gravity::Model<>& NN, IndexSet& inds, gravity::param<>& w, gravity::var<>& x, gravity::var<int>& y) override {
        Constraint<> Div_("Div");
        Div_ = x.in(inds["Out"])*x.in(inds["B"]) - x.in(inds["A"]);
        NN.add(Div_.in(inds["Constr"]) == 0);
    }
};
