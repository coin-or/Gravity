#pragma once

#include <network/Layers/LayerBase.hpp>
#include <utils/shape_iterator.hpp>

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

    void index_hidden_states(indices& hidden_states, indices& y_ids) const override {
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

    void index_constraint(IndexSet& inds) const override {
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

    void add_constraints(gravity::Model<>& NN, IndexSet& inds, gravity::param<>& w, gravity::var<>& x, gravity::var<int>& y) const override {
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

    void index_hidden_states(indices& hidden_states, indices& y_ids) const override {
        for (auto output: this->outputs) {
            for (auto i = 0; i < output->numel; i++) {
                hidden_states.add(output->strkey(i));
                y_ids.add(output->strkey(i));
            }
        }
    }

    void index_constraint(IndexSet& inds) const override {
        for(auto j = 0; j < this->X->numel;j++){
            inds["Constr"].add(this->Y->strkey(j));
            inds["In"].add_ref(this->X->strkey(j));
            inds["Out"].add_ref(this->Y->strkey(j));
        }
    }

    void add_constraints(gravity::Model<>& NN, IndexSet& inds, gravity::param<>& w, gravity::var<>& x, gravity::var<int>& y) const override {
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

    void index_constraint(IndexSet& inds) const override {}
    std::vector<std::vector<std::string>> get_indices() const override {
        return {{}, {}};
    }

    void add_constraints(gravity::Model<>& NN, IndexSet& inds, gravity::param<>& w, gravity::var<>& x, gravity::var<int>& y) const override {}

    void set_bounds(gravity::param<>& x_lb, gravity::param<>& x_ub) const override {
        for (auto o: this->outputs) {
            for(auto j = 0; j < o->numel; j++){
                auto key  = o->strkey(j);
                auto lb = std::max(
                    o->lb.at(j),
                    this->range_lower
                );
                auto ub = std::min(
                    o->ub.at(j),
                    this->range_upper
                );

                x_lb.set_val(key, std::max(lb, -100.0f));
                x_ub.set_val(key, std::min(ub,  100.0f));
            }
        }
    }
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

    void index_constraint(IndexSet& inds) const override {
        for(auto j = 0; j < this->X->numel; j++){
            inds["Constr"].add(this->Y->strkey(j));
            inds["In"].add_ref(this->X->strkey(j));
            inds["Out"].add_ref(this->Y->strkey(j));
        }
    }

    void add_constraints(gravity::Model<>& NN, IndexSet& inds, gravity::param<>& w, gravity::var<>& x, gravity::var<int>& y) const override {
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
            auto tmp_axes = std::vector<int>(axis_attr->ints().begin(), axis_attr->ints().end());
            for (auto ax: tmp_axes) {
                if (ax < 0) {
                    ax += this->X->ndims;
                }
                this->axes.push_back(ax);
            }
        }

        if (const auto* keepdims_attr = find_attribute("keepdims", node)) {
            this->keepdims = keepdims_attr->i();
        }
    }

    std::vector<std::vector<std::string>> get_indices() const override {
        return {{"Out", "In"}, {}};
    }

    void index_constraint(IndexSet& inds) const override {
        std::vector<size_t> reduce_shape;
        for (auto& ax: this->axes) {
            reduce_shape.push_back(this->X->shape.at(ax));
        }

        for (auto index: ShapeIter(this->Y->shape)) {
            inds["Constr"].add(this->Y->strkey(index));
            inds["Out"].add_ref(this->Y->strkey(index));
            for (auto subind: ShapeIter(reduce_shape)) {
                auto xind = index;
                for (int ax_id = 0; ax_id < this->axes.size(); ax_id++) {
                    if (keepdims == 0) {
                        xind.insert(xind.begin() + this->axes.at(ax_id), subind.at(ax_id));
                    } else {
                        xind.at(this->axes.at(ax_id)) = subind.at(ax_id);
                    }
                }
                inds["In"].add_in_row(inds.row_id, this->X->strkey(xind));
            }
            inds.row_id++;
        }
    }

    void add_constraints(gravity::Model<>& NN, IndexSet& inds, gravity::param<>& w, gravity::var<>& x, gravity::var<int>& y) const override {
        Constraint<> RDSum("ReduceSum");
        RDSum = x.in(inds["Out"]) - x.in(inds["In"]);
        NN.add(RDSum.in(inds["Constr"]) == 0);
    }

    Tensor *X, *Y; // Input and output
    std::vector<size_t> axes;
    int keepdims = 1;
};

class ReduceMean : public Layer {
public:
    ReduceMean(const onnx::NodeProto& node, Tensors& tensors): Layer(node, tensors) {
        operator_type = _reduce_mean;
        this->X = &tensors[node.input(0)];
        this->Y = &tensors[node.output(0)];
        if (const auto* axis_attr = find_attribute("axes", node)) {
            auto tmp_axes = std::vector<int>(axis_attr->ints().begin(), axis_attr->ints().end());
            for (auto ax: tmp_axes) {
                if (ax < 0) {
                    ax += this->X->ndims;
                }
                this->axes.push_back(ax);
            }
        }

        if (const auto* keepdims_attr = find_attribute("keepdims", node)) {
            this->keepdims = keepdims_attr->i();
        }
    }

    std::vector<std::vector<std::string>> get_indices() const override {
        return {{"Out", "In"}, {}};
    }

    void index_constraint(IndexSet& inds) const override {
        std::vector<size_t> reduce_shape;
        for (auto& ax: this->axes) {
            reduce_shape.push_back(this->X->shape.at(ax));
        }

        for (auto index: ShapeIter(this->Y->shape)) {
            inds["Constr"].add(this->Y->strkey(index));
            inds["Out"].add_ref(this->Y->strkey(index));
            for (auto subind: ShapeIter(reduce_shape)) {
                auto xind = index;
                for (int ax_id = 0; ax_id < this->axes.size(); ax_id++) {
                    if (keepdims == 0) {
                        xind.insert(xind.begin() + this->axes.at(ax_id), subind.at(ax_id));
                    } else {
                        xind.at(this->axes.at(ax_id)) = subind.at(ax_id);
                    }
                }
                inds["In"].add_in_row(inds.row_id, this->X->strkey(xind));
            }
            inds.row_id++;
        }
    }

    void add_constraints(gravity::Model<>& NN, IndexSet& inds, gravity::param<>& w, gravity::var<>& x, gravity::var<int>& y) const override {
        size_t reduce_size = 1;
        for (auto& ax: this->axes) {
            reduce_size *= this->X->shape.at(ax);
        }

        Constraint<> RDSum("ReduceSum");
        RDSum = x.in(inds["Out"]) - (x.in(inds["In"]) * (1.0 / reduce_size));
        NN.add(RDSum.in(inds["Constr"]) == 0);
    }

    Tensor *X, *Y; // Input and output
    std::vector<size_t> axes;
    int keepdims = 1;
};