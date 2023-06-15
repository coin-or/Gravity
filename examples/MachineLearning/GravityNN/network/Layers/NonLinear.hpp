#pragma once

#include <network/Layers/LayerBase.hpp>

using namespace gravity;

class Cos : public Layer {
public:
    Cos(const onnx::NodeProto& node, Tensors& tensors): Layer(node, tensors) {
        operator_type = _cos;
        this->X = &tensors[node.input(0)];
        this->Y = &tensors[node.output(0)];

        this->range_lower = -1.0;
        this->range_upper =  1.0;
    }

    std::vector<std::vector<std::string>> get_indices() const override {
        return {{"Out", "In"}, {}};
    }

    void index_constraint(IndexSet& inds) const override {
        for(auto j = 0; j < this->X->numel;j++){
            inds["Constr"].add(this->Y->strkey(j));
            inds["In"].add_ref(this->X->strkey(j));
            inds["Out"].add_ref(this->Y->strkey(j));
        }
    }

    void add_constraints(gravity::Model<>& NN, IndexSet& inds, gravity::param<>& w, gravity::var<>& x, gravity::var<int>& y) const override {
        Constraint<> Cos_("Cos");
        Cos_ = x.in(inds["Out"]) - cos(x.in(inds["In"]));
        NN.add(Cos_.in(inds["Constr"]) == 0);
    }

    Tensor *X, *Y; // Input and output
};

class Sin : public Layer {
public:
    Sin(const onnx::NodeProto& node, Tensors& tensors): Layer(node, tensors) {
        operator_type = _sin;
        this->X = &tensors[node.input(0)];
        this->Y = &tensors[node.output(0)];

        this->range_lower = -1.0;
        this->range_upper =  1.0;
    }

    std::vector<std::vector<std::string>> get_indices() const override {
        return {{"Out", "In"}, {}};
    }

    void index_constraint(IndexSet& inds) const override {
        for(auto j = 0; j < this->X->numel;j++){
            inds["Constr"].add(this->Y->strkey(j));
            inds["In"].add_ref(this->X->strkey(j));
            inds["Out"].add_ref(this->Y->strkey(j));
        }
    }

    void add_constraints(gravity::Model<>& NN, IndexSet& inds, gravity::param<>& w, gravity::var<>& x, gravity::var<int>& y) const override {
        Constraint<> Sin_("Sin");
        Sin_ = x.in(inds["Out"]) - sin(x.in(inds["In"]));
        NN.add(Sin_.in(inds["Constr"]) == 0);
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

        this->range_lower = 0.0;
    }

    std::vector<std::vector<std::string>> get_indices() const override {
        return {{"Out", "In"}, {}};
    }

    void index_constraint(IndexSet& inds) const override {
        for(auto j = 0; j < this->X->numel;j++){
            inds["Constr"].add(this->Y->strkey(j));
            inds["In"].add_ref(this->X->strkey(j));
            inds["Out"].add_ref(this->Y->strkey(j));
        }
    }

    void add_constraints(gravity::Model<>& NN, IndexSet& inds, gravity::param<>& w, gravity::var<>& x, gravity::var<int>& y) const override {
        Constraint<> Pow_("Pow");
        Pow_ = x.in(inds["Out"]) - pow(x.in(inds["In"]), 2.0);
        NN.add(Pow_.in(inds["Constr"]) == 0);
    }


    Tensor *X, *exp, *Y; // Input and output
};

class Exp: public Layer {
public:
    Exp(const onnx::NodeProto& node, Tensors& tensors): Layer(node, tensors) {
        operator_type = _exp;
        this->X = &tensors[node.input(0)];
        this->Y = &tensors[node.output(0)];

        this->range_lower = 0.0;
    }

    std::vector<std::vector<std::string>> get_indices() const override {
        return {{"Out", "In"}, {}};
    }

    void index_constraint(IndexSet& inds) const override {
        for(auto j = 0; j < this->X->numel;j++){
            inds["Constr"].add(this->Y->strkey(j));
            inds["In"].add_ref(this->X->strkey(j));
            inds["Out"].add_ref(this->Y->strkey(j));
        }
    }

    void add_constraints(gravity::Model<>& NN, IndexSet& inds, gravity::param<>& w, gravity::var<>& x, gravity::var<int>& y) const override {
        Constraint<> Exp_("Exp");
        Exp_ = x.in(inds["Out"]) - gravity::exp(x.in(inds["In"]));
        NN.add(Exp_.in(inds["Constr"]) == 0);
    }

    Tensor *X, *Y; // Input and output
};

class Sigmoid: public Layer {
public:
    Sigmoid(const onnx::NodeProto& node, Tensors& tensors): Layer(node, tensors) {
        operator_type = _sigmoid;
        this->X = &tensors[node.input(0)];
        this->Y = &tensors[node.output(0)];

        this->range_lower = 0.0;
        this->range_upper = 1.0;
    }

    std::vector<std::vector<std::string>> get_indices() const override {
        return {{"Out", "In", "ExpAux", "NegAux"}, {}};
    }

    void index_hidden_states(indices& hidden_states, indices& y_ids) const override {
        for (auto output: this->outputs) {
            for (auto i = 0; i < output->numel; i++) {
                hidden_states.add(output->strkey(i));
                hidden_states.add(output->strkey(i) + "_exp_aux");
                hidden_states.add(output->strkey(i) + "_neg_aux");
            }
        }
    }

    void index_constraint(IndexSet& inds) const override {
        for(auto j = 0; j < this->X->numel;j++){
            inds["Constr"].add(this->Y->strkey(j));
            inds["In"].add_ref(this->X->strkey(j));
            inds["ExpAux"].add_ref(this->Y->strkey(j) + "_exp_aux");
            inds["NegAux"].add_ref(this->Y->strkey(j) + "_neg_aux");
            inds["Out"].add_ref(this->Y->strkey(j));
        }
    }

    void add_constraints(gravity::Model<>& NN, IndexSet& inds, gravity::param<>& w, gravity::var<>& x, gravity::var<int>& y) const override {
        // Negate x
        Constraint<> Neg_("Sigmoid_Neg_Aux");
        Neg_ = x.in(inds["NegAux"]) + x.in(inds["In"]);
        NN.add(Neg_.in(inds["Constr"]) == 0);

        // Exp constraint
        Constraint<> Exp_("Sigmoid_Aux");
        Exp_ = x.in(inds["ExpAux"]) - gravity::exp(x.in(inds["NegAux"]));
        NN.add(Exp_.in(inds["Constr"]) == 0);

        // Rest of sigmoid
        Constraint<> Sigmoid_("Sigmoid");
        Sigmoid_ = x.in(inds["Out"])*(1.0 + x.in(inds["ExpAux"])) - 1.0;
        NN.add(Sigmoid_.in(inds["Constr"]) == 0);
    }

    Tensor *X, *Y; // Input and output
};

class Softmax: public Layer {
public:
    Softmax(const onnx::NodeProto& node, Tensors& tensors): Layer(node, tensors) {
        operator_type = _softmax;
        this->X = &tensors[node.input(0)];
        this->Y = &tensors[node.output(0)];

        if (const auto* axis = find_attribute("axis", node)) {
            int64_t i = axis->i();
            if (i < 0) {
                i += this->X->ndims;
            }
            this->axis = i;
        } else {
            throw std::runtime_error("Softmax: axis attribute not found");
        }

        this->range_lower = 0.0;
        this->range_upper = 1.0;
    }

    std::vector<std::vector<std::string>> get_indices() const override {
        return {{"Out", "In", "ExpAux", "ExpSum", "SumAux", "SumProd"}, {}};
    }

    void index_hidden_states(indices& hidden_states, indices& y_ids) const override {
        auto outer_shape = this->X->shape;
        outer_shape.at(this->axis) = 1;
        for (auto outer_ind: ShapeIter(outer_shape)) {
            hidden_states.add(this->Y->strkey(outer_ind) + "_sum_aux");

            for (size_t axis_ind = 0; axis_ind < this->X->shape.at(this->axis); axis_ind++) {
                auto inner_ind = outer_ind;
                inner_ind.at(this->axis) = axis_ind;

                hidden_states.add(this->Y->strkey(inner_ind));
                hidden_states.add(this->Y->strkey(inner_ind) + "_exp_aux");
            }
        }
    }

    void index_constraint(IndexSet& inds) const override {
        auto outer_shape = this->X->shape;
        outer_shape.at(this->axis) = 1;
        for (auto outer_ind: ShapeIter(outer_shape)) {
            inds["ConstrB"].add(this->Y->strkey(outer_ind) + "_sum");
            inds["SumAux"].add_ref(this->Y->strkey(outer_ind) + "_sum_aux");

            for (size_t axis_ind = 0; axis_ind < this->X->shape.at(this->axis); axis_ind++) {
                auto inner_ind = outer_ind;
                inner_ind.at(this->axis) = axis_ind;

                inds["Constr"].add(this->Y->strkey(inner_ind));
                inds["In"].add_ref(this->X->strkey(inner_ind));
                inds["ExpAux"].add_ref(this->Y->strkey(inner_ind) + "_exp_aux");
                inds["Out"].add_ref(this->Y->strkey(inner_ind));

                inds["ExpSum"].add_in_row(inds.row_id, this->Y->strkey(inner_ind) + "_exp_aux");
                inds["SumProd"].add_ref(this->Y->strkey(outer_ind) + "_sum_aux");
            }
            inds.row_id++;
        }
    }

    void add_constraints(gravity::Model<>& NN, IndexSet& inds, gravity::param<>& w, gravity::var<>& x, gravity::var<int>& y) const override {
        // Exp constraint
        Constraint<> Exp_("Softmax_Exp");
        Exp_ = x.in(inds["ExpAux"]) - gravity::exp(x.in(inds["In"]));
        NN.add(Exp_.in(inds["Constr"]) == 0);

        // Sum constraint
        Constraint<> Sum_("Softmax_Sum");
        Sum_ = x.in(inds["SumAux"]) - x.in(inds["ExpSum"]);
        Sum_.print();
        NN.add(Sum_.in(inds["ConstrB"]) == 0);

        // Out constraint
        Constraint<> Out_("Softmax_Out");
        Out_ = x.in(inds["Out"])*x.in(inds["SumProd"]) - x.in(inds["ExpAux"]);
        NN.add(Out_.in(inds["Constr"]) == 0);
    }

    Tensor *X, *Y; // Input and output
    size_t axis;
};