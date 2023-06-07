#pragma once

#include <network/Layers/LayerBase.hpp>

using namespace gravity;

class Cos : public Layer {
public:
    Cos(const onnx::NodeProto& node, Tensors& tensors): Layer(node, tensors) {
        operator_type = _cos;
        this->X = &tensors[node.input(0)];
        this->Y = &tensors[node.output(0)];
    }

    std::vector<std::vector<std::string>> get_indices() const override {
        return {{"Out", "In"}, {}};
    }

    void index_constraint(IndexSet& inds) override {
        for(auto j = 0; j < this->X->numel;j++){
            inds["Constr"].add(this->Y->strkey(j));
            inds["In"].add_ref(this->X->strkey(j));
            inds["Out"].add_ref(this->Y->strkey(j));
        }
    }

    void add_constraints(gravity::Model<>& NN, IndexSet& inds, gravity::param<>& w, gravity::var<>& x, gravity::var<int>& y) override {
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
    }

    std::vector<std::vector<std::string>> get_indices() const override {
        return {{"Out", "In"}, {}};
    }

    void index_constraint(IndexSet& inds) override {
        for(auto j = 0; j < this->X->numel;j++){
            inds["Constr"].add(this->Y->strkey(j));
            inds["In"].add_ref(this->X->strkey(j));
            inds["Out"].add_ref(this->Y->strkey(j));
        }
    }

    void add_constraints(gravity::Model<>& NN, IndexSet& inds, gravity::param<>& w, gravity::var<>& x, gravity::var<int>& y) override {
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
    }

    std::vector<std::vector<std::string>> get_indices() const override {
        return {{"Out", "In"}, {}};
    }

    void index_constraint(IndexSet& inds) override {
        for(auto j = 0; j < this->X->numel;j++){
            inds["Constr"].add(this->Y->strkey(j));
            inds["In"].add_ref(this->X->strkey(j));
            inds["Out"].add_ref(this->Y->strkey(j));
        }
    }

    void add_constraints(gravity::Model<>& NN, IndexSet& inds, gravity::param<>& w, gravity::var<>& x, gravity::var<int>& y) override {
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
    }

    std::vector<std::vector<std::string>> get_indices() const override {
        return {{"Out", "In"}, {}};
    }

    void index_constraint(IndexSet& inds) override {
        for(auto j = 0; j < this->X->numel;j++){
            inds["Constr"].add(this->Y->strkey(j));
            inds["In"].add_ref(this->X->strkey(j));
            inds["Out"].add_ref(this->Y->strkey(j));
        }
    }

    void add_constraints(gravity::Model<>& NN, IndexSet& inds, gravity::param<>& w, gravity::var<>& x, gravity::var<int>& y) override {
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
    }

    std::vector<std::vector<std::string>> get_indices() const override {
        return {{"Out", "In", "ExpAux"}, {}};
    }

    void index_hidden_states(indices& hidden_states, indices& y_ids) override {
        for (auto output: this->outputs) {
            for (auto i = 0; i < output->numel; i++) {
                hidden_states.add(output->strkey(i));
                hidden_states.add(output->strkey(i) + "_exp_aux");
            }
        }
    }

    void index_constraint(IndexSet& inds) override {
        for(auto j = 0; j < this->X->numel;j++){
            inds["Constr"].add(this->Y->strkey(j));
            inds["In"].add_ref(this->X->strkey(j));
            inds["ExpAux"].add_ref(this->Y->strkey(j) + "_exp_aux");
            inds["Out"].add_ref(this->Y->strkey(j));
        }
    }

    void add_constraints(gravity::Model<>& NN, IndexSet& inds, gravity::param<>& w, gravity::var<>& x, gravity::var<int>& y) override {
        // Exp constraint
        // Constraint<> Exp_("Sigmoid_Aux");
        // Exp_ = x.in(inds["ExpAux"]) - gravity::exp(x.in(inds["In"]));
        // NN.add(Exp_.in(inds["Constr"]) == 0);

        // y >= x+1, y>= exp(-1)*x + 2*exp(-1)
        Constraint<> ExpLin1_("ExpLin1");
        ExpLin1_ = x.in(inds["ExpAux"]) - (x.in(inds["In"]) + 1.0);
        NN.add(ExpLin1_.in(inds["Constr"]) >= 0);

        Constraint<> ExpLin2_("ExpLin2");
        ExpLin2_ = x.in(inds["ExpAux"]) - (exp(-1.0)*x.in(inds["In"]) + 2.0*exp(-1.0));
        NN.add(ExpLin2_.in(inds["Constr"]) >= 0);

        // Rest of sigmoid
        Constraint<> Sigmoid_("Sigmoid");
        Sigmoid_ = x.in(inds["Out"])*(x.in(inds["ExpAux"]) + 1.0) - x.in(inds["ExpAux"]);
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

        if ((this->axis != 1) || (this->X->ndims != 2)) {
            throw std::runtime_error("Softmax: axis != 1 or ndims != 2 not supported yet");
        }

        if (this->X->shape[0] != 1) {
            throw std::runtime_error("Softmax: batch size > 1 not supported yet");
        }
    }

    std::vector<std::vector<std::string>> get_indices() const override {
        return {{"Out", "In", "ExpAux", "ExpSum", "SumAux", "SumProd"}, {}};
    }

    void index_hidden_states(indices& hidden_states, indices& y_ids) override {
        for (auto output: this->outputs) {
            for (auto i = 0; i < output->numel; i++) {
                hidden_states.add(output->strkey(i));
                hidden_states.add(output->strkey(i) + "_exp_aux");
            }
        }
        hidden_states.add(this->Y->name + "_sum_aux");
    }

    void index_constraint(IndexSet& inds) override {
        for(auto j = 0; j < this->X->numel;j++){
            inds["Constr"].add(this->Y->strkey(j));
            inds["In"].add_ref(this->X->strkey(j));
            inds["ExpAux"].add_ref(this->Y->strkey(j) + "_exp_aux");
            inds["Out"].add_ref(this->Y->strkey(j));

            inds["ExpSum"].add_in_row(inds.row_id, this->Y->strkey(j) + "_exp_aux");
            inds["SumProd"].add_ref(this->Y->name + "_sum_aux");
        }
        inds["ConstrB"].add(this->Y->name + "_sum");
        inds["SumAux"].add_ref(this->Y->name + "_sum_aux");
        inds.row_id++;
    }

    void add_constraints(gravity::Model<>& NN, IndexSet& inds, gravity::param<>& w, gravity::var<>& x, gravity::var<int>& y) override {
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