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

        this->range_lower = this->min;
        this->range_upper = this->max;
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

    void index_constraint(IndexSet& inds) override {
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

        this->range_lower = 0.0;
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

    void index_constraint(IndexSet& inds) override {
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

    void index_constraint(IndexSet& inds) override {}
    std::vector<std::vector<std::string>> get_indices() const override {
        return {{}, {}};
    }

    void add_constraints(gravity::Model<>& NN, IndexSet& inds, gravity::param<>& w, gravity::var<>& x, gravity::var<int>& y) override {}
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

    void index_constraint(IndexSet& inds) override {
        for(auto j = 0; j < this->X->numel; j++){
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

class ReduceSum : public Layer {
public:
    ReduceSum(const onnx::NodeProto& node, Tensors& tensors): Layer(node, tensors) {
        operator_type = _reduce_sum;
        this->X = &tensors[node.input(0)];
        this->Y = &tensors[node.output(0)];
        if (const auto* axis_attr = find_attribute("axes", node)) {
            this->axes = std::vector<int>(axis_attr->ints().begin(), axis_attr->ints().end());
        }

        for (auto& ax: this->axes) {
            if (ax < 0) {
                ax += this->X->ndims;
            }
        }

        if (this->axes.size() == 0) {
            throw std::runtime_error("ReduceSum: Reduction over the full tensor is not supported");
        }
        if ((this->axes.size() != 1) || (this->axes.at(0) != 1)) {
            throw std::runtime_error("ReduceSum: Reduction over axis != 1 is not supported");
        }
    }

    std::vector<std::vector<std::string>> get_indices() const override {
        return {{"Out", "In"}, {}};
    }

    void index_constraint(IndexSet& inds) override {
        inds["Constr"].add(this->Y->strkey(0));
        inds["Out"].add_ref(this->Y->strkey(0));
        for(auto j = 0; j < this->X->numel; j++){
            inds["In"].add_in_row(inds.row_id, this->X->strkey(j));
        }
        inds.row_id++;
    }

    void add_constraints(gravity::Model<>& NN, IndexSet& inds, gravity::param<>& w, gravity::var<>& x, gravity::var<int>& y) override {
        Constraint<> RDSum("ReduceSum");
        RDSum = x.in(inds["Out"]) - x.in(inds["In"]);
        NN.add(RDSum.in(inds["Constr"]) == 0);
    }

    Tensor *X, *Y; // Input and output
    std::vector<int> axes;
};