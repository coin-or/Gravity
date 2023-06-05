#pragma once

#include <network/Layers/LayerBase.hpp>

using namespace gravity;

class Add: public Layer {
public:
    Add(const onnx::NodeProto& node, Tensors& tensors): Layer(node, tensors) {
        operator_type = _add;
        this->A = &tensors[node.input(0)];
        this->B = &tensors[node.input(1)];
        this->Y = &tensors[node.output(0)];

        if ((this->A->is_initializer == true) && (this->B->is_initializer == true)) {
            throw std::runtime_error("Add: both args being initializer not supported.");
        }

        // Make sure that B is the initializer
        if (this->A->is_initializer == true) {
            std::swap(this->A, this->B);
        }

        if (this->B->numel != this->A->numel) {
            throw std::runtime_error("Add: initializer must have same size as input.");
        }
    }

    void add_parameters(gravity::param<>& w) const override {
        if (this->B->is_initializer) {
            this->B->add_params(w);
        }
    }

    std::vector<std::vector<std::string>> get_indices() const override {
        return {{"hOut", "pOut", "A", "B", "pA"}, {"AddVar"}};
    }

    void build_constraint(IndexSet& inds) override {
        for(auto j = 0; j < this->A->numel;j++){

            if (this->B->is_initializer) {
                inds["ConstrB"].add(this->Y->strkey(j));
                inds["pA"].add_ref(this->A->strkey(j));
                inds["AddVar"].add_ref(this->B->strkey(j));
                inds["pOut"].add_ref(this->Y->strkey(j));
            } else {
                inds["Constr"].add(this->Y->strkey(j));
                inds["A"].add_ref(this->A->strkey(j));
                inds["B"].add_ref(this->B->strkey(j));
                inds["hOut"].add_ref(this->Y->strkey(j));
            }
        }
    }

    void add_constraints(gravity::Model<>& NN, IndexSet& inds, gravity::param<>& w, gravity::var<>& x, gravity::var<int>& y) override {
        // Constraint for addition of hidden states
        Constraint<> Add_("Add");
        Add_ = x.in(inds["hOut"]) - (x.in(inds["A"]) + x.in(inds["B"]));
        NN.add(Add_.in(inds["Constr"]) == 0);

        // Constraint where B is a parameter
        Constraint<> Add_Param_("Add_Param");
        Add_Param_ = x.in(inds["pOut"]) - (x.in(inds["pA"]) + w.in(inds["AddVar"]));
        NN.add(Add_Param_.in(inds["ConstrB"]) == 0);
    }

    Tensor *A, *B, *Y; // Input and output
};

class Sub: public Layer {
public:
    Sub(const onnx::NodeProto& node, Tensors& tensors): Layer(node, tensors) {
        operator_type = _sub;
        this->A = &tensors[node.input(0)];
        this->B = &tensors[node.input(1)];
        this->Y = &tensors[node.output(0)];

        if ((this->A->is_initializer == true) && (this->B->is_initializer == true)) {
            throw std::runtime_error("Sub: both args being initializer not supported.");
        }

        // Make sure that B is the initializer
        if (this->A->is_initializer == true) {
            std::swap(this->A, this->B);
        }

        if (this->B->numel != this->A->numel) {
            throw std::runtime_error("Sub: initializer must have same size as input.");
        }
    }

    void add_parameters(gravity::param<>& w) const override {
        if (this->B->is_initializer) {
            this->B->add_params(w);
        }
    }

    std::vector<std::vector<std::string>> get_indices() const override {
        return {{"hOut", "pOut", "A", "B", "pA"}, {"SubVar"}};
    }

    void build_constraint(IndexSet& inds) override {
        for(auto j = 0; j < this->A->numel;j++){

            if (this->B->is_initializer) {
                inds["ConstrB"].add(this->Y->strkey(j));
                inds["pA"].add_ref(this->A->strkey(j));
                inds["SubVar"].add_ref(this->B->strkey(j));
                inds["pOut"].add_ref(this->Y->strkey(j));
            } else {
                inds["Constr"].add(this->Y->strkey(j));
                inds["A"].add_ref(this->A->strkey(j));
                inds["B"].add_ref(this->B->strkey(j));
                inds["hOut"].add_ref(this->Y->strkey(j));
            }
        }
    }

    void add_constraints(gravity::Model<>& NN, IndexSet& inds, gravity::param<>& w, gravity::var<>& x, gravity::var<int>& y) override {
        // Constraint for sub of hidden states
        Constraint<> Sub_("Sub");
        Sub_ = x.in(inds["hOut"]) - (x.in(inds["A"]) - x.in(inds["B"]));
        NN.add(Sub_.in(inds["Constr"]) == 0);

        // Constraint where B is a parameter
        Constraint<> Sub_Param_("Sub_Param");
        Sub_Param_ = x.in(inds["pOut"]) - (x.in(inds["pA"]) - w.in(inds["SubVar"]));
        NN.add(Sub_Param_.in(inds["ConstrB"]) == 0);
    }

    Tensor *A, *B, *Y; // Input and output
};

class Mul: public Layer {
public:
    Mul(const onnx::NodeProto& node, Tensors& tensors): Layer(node, tensors) {
        operator_type = _mul;
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

        if (this->B->numel != this->A->numel) {
            throw std::runtime_error("Mul: initializer must have same size as input.");
        }
    }

    void add_parameters(gravity::param<>& w) const override {
        if (this->B->is_initializer) {
            this->B->add_params(w);
        }
    }

    std::vector<std::vector<std::string>> get_indices() const override {
        return {{"hOut", "pOut", "A", "B", "pA"}, {"MulVar"}};
    }

    void build_constraint(IndexSet& inds) override {
        for(auto j = 0; j < this->A->numel;j++){

            if (this->B->is_initializer) {
                inds["ConstrB"].add(this->Y->strkey(j));
                inds["pA"].add_ref(this->A->strkey(j));
                inds["MulVar"].add_ref(this->B->strkey(j));
                inds["pOut"].add_ref(this->Y->strkey(j));
            } else {
                inds["Constr"].add(this->Y->strkey(j));
                inds["A"].add_ref(this->A->strkey(j));
                inds["B"].add_ref(this->B->strkey(j));
                inds["hOut"].add_ref(this->Y->strkey(j));
            }
        }
    }

    void add_constraints(gravity::Model<>& NN, IndexSet& inds, gravity::param<>& w, gravity::var<>& x, gravity::var<int>& y) override {
        // Constraint for mul of hidden states
        Constraint<> Mul_("Mul");
        Mul_ = x.in(inds["hOut"]) - (x.in(inds["A"]) * x.in(inds["B"]));
        NN.add(Mul_.in(inds["Constr"]) == 0);

        // Constraint where B is a parameter
        Constraint<> Mul_Param_("Mul_Param");
        Mul_Param_ = x.in(inds["pOut"]) - (x.in(inds["pA"]) * w.in(inds["MulVar"]));
        NN.add(Mul_Param_.in(inds["ConstrB"]) == 0);
    }

    Tensor *A, *B, *Y; // Input and output
};

class Div : public Layer {
public:
    Div(const onnx::NodeProto& node, Tensors& tensors): Layer(node, tensors) {
        operator_type = _div;
        this->A = &tensors[node.input(0)];
        this->B = &tensors[node.input(1)];
        this->Y = &tensors[node.output(0)];

        if ((this->A->is_initializer == true) || (this->B->is_initializer == true)) {
            throw std::runtime_error("Div: initializer not supported.");
        }
    }

    std::vector<std::vector<std::string>> get_indices() const override {
        return {{"Out", "A", "B"}, {}};
    }

    void build_constraint(IndexSet& inds) override {
        for(auto j = 0; j < this->A->numel;j++){
            inds["Constr"].add(this->Y->strkey(j));
            inds["A"].add_ref(this->A->strkey(j));
            inds["B"].add_ref(this->B->strkey(j));
            inds["Out"].add_ref(this->Y->strkey(j));
        }
    }

    void add_constraints(gravity::Model<>& NN, IndexSet& inds, gravity::param<>& w, gravity::var<>& x, gravity::var<int>& y) override {
        Constraint<> Div_("Div");
        Div_ = x.in(inds["Out"])*x.in(inds["B"]) - x.in(inds["A"]);
        NN.add(Div_.in(inds["Constr"]) == 0);
    }

    Tensor *A, *B, *Y; // Input and output
};
