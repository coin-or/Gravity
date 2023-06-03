#pragma once

#include <network/Layers/LayerBase.hpp>

using namespace gravity;

class Clip : public Layer {
public:
    Clip(const onnx::NodeProto& node, Tensors& tensors): Layer(node, tensors) {
        operator_type = _clip;
        this->X = &tensors.at(node.input(0));
        this->Y = &tensors.at(node.output(0));

        if (node.input_size() > 1) {
            this->min = tensors.at(node.input(1))(0);
        }
        if (node.input_size() > 2) {
            this->max = tensors.at(node.input(2))(0);
        }
    }

    std::vector<std::vector<std::string>> get_indices() const override {
        return {{"In", "Out"}, {"MinP", "MaxP"}, {"MinB", "MaxB", "Eq"}};
    }

    void index_hidden_states(indices& hidden_states, indices& y_ids) override {
        for (auto output: this->outputs) {
            for (auto i = 0; i < output->numel; i++) {
                hidden_states.add(output->strkey(i));
                y_ids.add(output->strkey(i) + "_min");
                y_ids.add(output->strkey(i) + "_max");
                y_ids.add(output->strkey(i) + "_eq");
            }
        }
    }

    void add_parameters(gravity::param<>& w) const override {
        for (auto i = 0; i < this->Y->numel; i++) {
            w.add_val(this->Y->strkey(i) + "_min", this->min);
            w.add_val(this->Y->strkey(i) + "_max", this->max);
        }
    }

    void build_constraint(IndexSet& inds) override {
        for(auto j = 0; j < this->X->numel; j++){
            inds["Constr"].add(this->Y->strkey(j));

            inds["In"].add_ref(this->X->strkey(j));
            inds["Out"].add_ref(this->Y->strkey(j));

            inds["MinP"].add_ref(this->Y->strkey(j) + "_min");
            inds["MaxP"].add_ref(this->Y->strkey(j) + "_max");

            inds["MinB"].add_ref(this->Y->strkey(j) + "_min");
            inds["MaxB"].add_ref(this->Y->strkey(j) + "_max");
            inds["Eq"].add_ref(this->Y->strkey(j)  + "_eq");
        }
    }

    void add_constraints(gravity::Model<>& NN, IndexSet& inds, gravity::param<>& w, gravity::var<>& x, gravity::var<int>& y) override {
        Constraint<> ClipU("Upper");
        ClipU = x.in(inds["Out"]) - w.in(inds["MaxP"]);
        NN.add_on_off(ClipU.in(inds["Constr"]) == 0, y.in(inds["MaxB"]), true);

        Constraint<> ClipL("Lower");
        ClipL = x.in(inds["Out"]) - w.in(inds["MinP"]);
        NN.add_on_off(ClipL.in(inds["Constr"]) == 0, y.in(inds["MinB"]), true);

        Constraint<> ClipEQ("Equal");
        ClipEQ = x.in(inds["Out"]) - x.in(inds["In"]);
        NN.add_on_off(ClipEQ.in(inds["Constr"]) == 0, y.in(inds["Eq"]), true);

        Constraint<> ClipBin("Bin");
        ClipBin = y.in(inds["MaxB"]) + y.in(inds["MinB"]) + y.in(inds["Eq"]);
        NN.add(ClipBin.in(inds["Constr"]) == 1);

        Constraint<> ClipInU("InU");
        ClipInU = x.in(inds["In"]);
        NN.add_on_off(ClipInU.in(inds["Constr"]) <= w.in(inds["MaxP"]), y.in(inds["Eq"]), true);

        Constraint<> ClipInL("InL");
        ClipInL = x.in(inds["In"]);
        NN.add_on_off(ClipInL.in(inds["Constr"]) >= w.in(inds["MinP"]), y.in(inds["Eq"]), true);

        Constraint<> ClipInMax("InMax");
        ClipInMax = x.in(inds["In"]);
        NN.add_on_off(ClipInMax.in(inds["Constr"]) >= w.in(inds["MaxP"]), y.in(inds["MaxB"]), true);

        Constraint<> ClipInMin("InMin");
        ClipInMin = x.in(inds["In"]);
        NN.add_on_off(ClipInMin.in(inds["Constr"]) <= w.in(inds["MinP"]), y.in(inds["MinB"]), true);
    }

    Tensor* X; // Input
    Tensor* Y; // Output

    float min = std::numeric_limits<float>::lowest();
    float max = std::numeric_limits<float>::max();
};

class Relu : public Layer {
public:
    Relu(const onnx::NodeProto& node, Tensors& tensors): Layer(node, tensors) {
        operator_type = _relu;
        this->X = &tensors.at(node.input(0));
        this->Y = &tensors.at(node.output(0));
    }

    std::vector<std::vector<std::string>> get_indices() const override {
        return {{"In", "Out"}, {}};
    }

    void index_hidden_states(indices& hidden_states, indices& y_ids) override {
        for (auto output: this->outputs) {
            for (auto i = 0; i < output->numel; i++) {
                hidden_states.add(output->strkey(i));
                y_ids.add(output->strkey(i));
            }
        }
    }

    void build_constraint(IndexSet& inds) override {
        for(auto j = 0; j < this->X->numel;j++){
            inds["Constr"].add(this->Y->strkey(j));
            inds["In"].add_ref(this->X->strkey(j));
            inds["Out"].add_ref(this->Y->strkey(j));
        }
    }

    void add_constraints(gravity::Model<>& NN, IndexSet& inds, gravity::param<>& w, gravity::var<>& x, gravity::var<int>& y) override {
        /* Constraints */
        Constraint<> ReLU("ReLU");
        ReLU = x.in(inds["Out"]) - x.in(inds["In"]);
        NN.add(ReLU.in(inds["Constr"]) >= 0);

        Constraint<> ReLU_on("ReLU_on");
        ReLU_on = x.in(inds["Out"]) - x.in(inds["In"]);
        NN.add_on_off(ReLU_on.in(inds["Constr"]) <= 0, y.in(inds["Constr"]), true);

        Constraint<> ReLU_y_off("ReLU_y_off");
        ReLU_y_off = x.in(inds["Out"]);
        NN.add_on_off(ReLU_y_off.in(inds["Constr"]) <= 0, y.in(inds["Constr"]), false);
    }

    Tensor* X; // Input
    Tensor* Y; // Output
};

class Input : public Layer {
public:
    Input(const onnx::ValueInfoProto& node, Tensors& tensors): Layer(node, tensors) {
        operator_type = _input;
    }

    void build_constraint(IndexSet& inds) override {}
    std::vector<std::vector<std::string>> get_indices() const override {
        return {{}, {}};
    }

    void add_constraints(gravity::Model<>& NN, IndexSet& inds, gravity::param<>& w, gravity::var<>& x, gravity::var<int>& y) override {}
};

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

class Neg : public Layer {
public:
    Neg(const onnx::NodeProto& node, Tensors& tensors): Layer(node, tensors) {
        operator_type = _neg;
        this->X = &tensors[node.input(0)];
        this->Y = &tensors[node.output(0)];
    }

    std::vector<std::vector<std::string>> get_indices() const override {
        return {{"Out", "In"}, {}};
    }

    void build_constraint(IndexSet& inds) override {
        for(auto j = 0; j < this->X->numel;j++){
            inds["Constr"].add(this->Y->strkey(j));
            inds["In"].add_ref(this->X->strkey(j));
            inds["Out"].add_ref(this->Y->strkey(j));
        }
    }

    void add_constraints(gravity::Model<>& NN, IndexSet& inds, gravity::param<>& w, gravity::var<>& x, gravity::var<int>& y) override {
        Constraint<> Neg_("Neg");
        Neg_ = x.in(inds["Out"]) + x.in(inds["In"]);
        NN.add(Neg_.in(inds["Constr"]) == 0);
    }

    Tensor *X, *Y; // Input and output
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

        // sqrt the variance
        auto data = this->var->get_data();
        for (auto& v: data) {
            v = sqrt(v);
        }
        this->var->_set_data(data);

        if (this->X->ndims > 2) {
            throw std::runtime_error("BatchNorm: ndims > 2 not supported yet");
        }

        if (this->X->shape[0] != 1) {
            throw std::runtime_error("BatchNorm: batch size > 1 not supported yet");
        }
    }

    std::vector<std::vector<std::string>> get_indices() const override {
        return {{"Out", "In", "AuxA", "AuxB", "AuxC"}, {"scale", "B", "mean", "var"}};
    }

    void index_hidden_states(indices& hidden_states, indices& y_ids) override {
        for (auto output: this->outputs) {
            for (auto i = 0; i < output->numel; i++) {
                hidden_states.add(output->strkey(i));
                hidden_states.add(output->strkey(i) + "_aux_a");
                hidden_states.add(output->strkey(i) + "_aux_b");
                hidden_states.add(output->strkey(i) + "_aux_c");
            }
        }
    }

    void build_constraint(IndexSet& inds) override {
        for(auto j = 0; j < this->X->numel;j++){
            inds["Constr"].add(this->Y->strkey(j));
            inds["In"].add_ref(this->X->strkey(j));
            inds["AuxA"].add_ref(this->Y->strkey(j) + "_aux_a");
            inds["AuxB"].add_ref(this->Y->strkey(j) + "_aux_b");
            inds["AuxC"].add_ref(this->Y->strkey(j) + "_aux_c");
            inds["Out"].add_ref(this->Y->strkey(j));

            inds["scale"].add_ref(this->scale->strkey(j));
            inds["B"].add_ref(this->B->strkey(j));
            inds["mean"].add_ref(this->mean->strkey(j));
            inds["var"].add_ref(this->var->strkey(j));
        }
    }

    void add_parameters(gravity::param<>& w) const override {
        this->scale->add_params(w);
        this->B->add_params(w);
        this->mean->add_params(w);
        this->var->add_params(w);
    }

    void add_constraints(gravity::Model<>& NN, IndexSet& inds, gravity::param<>& w, gravity::var<>& x, gravity::var<int>& y) override {
        Constraint<> Aux_A("BN_AuxA");
        Aux_A = x.in(inds["AuxA"]) - (x.in(inds["In"]) - w.in(inds["mean"]));
        NN.add(Aux_A.in(inds["Constr"]) == 0);

        Constraint<> Aux_B("BN_AuxB");
        Aux_B = x.in(inds["AuxB"])*w.in(inds["var"]) - x.in(inds["AuxA"]);
        NN.add(Aux_B.in(inds["Constr"]) == 0);

        Constraint<> Aux_C("BN_AuxC");
        Aux_C = x.in(inds["AuxC"]) - x.in(inds["AuxB"])*w.in(inds["scale"]);
        NN.add(Aux_C.in(inds["Constr"]) == 0);

        Constraint<> Aux_D("BN_Out");
        Aux_D = x.in(inds["Out"]) - (x.in(inds["AuxC"]) + w.in(inds["B"]));
        NN.add(Aux_D.in(inds["Constr"]) == 0);
    }

    Tensor *X, *Y; // Input and output
    Tensor *scale, *B, *mean, *var;
};