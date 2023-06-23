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

    double min = HMIN;
    double max = HMAX;
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

        // Using == 0 instead of <= 0 on these constraints seems better
        Constraint<> ReLU_on("ReLU_on");
        ReLU_on = x.in(inds["Out"]) - x.in(inds["In"]);
        NN.add_on_off(ReLU_on.in(inds["Constr"]) == 0, y.in(inds["Constr"]), true);

        Constraint<> ReLU_y_off("ReLU_y_off");
        ReLU_y_off = x.in(inds["Out"]);
        NN.add_on_off(ReLU_y_off.in(inds["Constr"]) == 0, y.in(inds["Constr"]), false);
    }

    Tensor* X; // Input
    Tensor* Y; // Output
};

class LeakyRelu : public Layer {
public:
    LeakyRelu(const onnx::NodeProto& node, Tensors& tensors): Layer(node, tensors) {
        operator_type = _leaky_relu;
        this->X = &tensors.at(node.input(0));
        this->Y = &tensors.at(node.output(0));

        if (const auto* alpha_attr = find_attribute("alpha", node)) {
            this->alpha = alpha_attr->f();
        }

        if (this->alpha > 1.0 || this->alpha < 0.0) {
            throw std::runtime_error("LeakyRelu: alpha must be between 0 and 1 in our formulation.");
        }
    }

    std::vector<std::vector<std::string>> get_indices() const override {
        return {{"In", "Out"}, {"alpha"}};
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
        for(auto j = 0; j < this->X->numel; j++){
            inds["Constr"].add(this->Y->strkey(j));
            inds["In"].add_ref(this->X->strkey(j));
            inds["Out"].add_ref(this->Y->strkey(j));
            inds["alpha"].add_ref(this->Y->strkey(j) + "_alpha");
        }
    }

    void add_parameters(gravity::param<>& w) const override {
        for (auto i = 0; i < this->Y->numel; i++) {
            w.add_val(this->Y->strkey(i) + "_alpha", this->alpha);
        }
    }

    void add_constraints(gravity::Model<>& NN, IndexSet& inds, gravity::param<>& w, gravity::var<>& x, gravity::var<int>& y) override {
        /* Constraints */
        Constraint<> LeakyReLU("LeakyReLU");
        LeakyReLU = x.in(inds["Out"]) - w.in(inds["alpha"])*x.in(inds["In"]);
        NN.add(LeakyReLU.in(inds["Constr"]) >= 0);

        Constraint<> LeakyReLU2("LeakyReLU2");
        LeakyReLU2 = x.in(inds["Out"]) - x.in(inds["In"]);
        NN.add(LeakyReLU2.in(inds["Constr"]) >= 0);

        Constraint<> LeakyReLU_on("LeakyReLU_on");
        LeakyReLU_on = x.in(inds["Out"]) - x.in(inds["In"]);
        NN.add_on_off(LeakyReLU_on.in(inds["Constr"]) == 0, y.in(inds["Constr"]), true);

        Constraint<> LeakyReLU_off("LeakyReLU_off");
        LeakyReLU_off = x.in(inds["Out"]) - w.in(inds["alpha"]) * x.in(inds["In"]);
        NN.add_on_off(LeakyReLU_off.in(inds["Constr"]) == 0, y.in(inds["Constr"]), false);
    }

    Tensor* X; // Input
    Tensor* Y; // Output

    double alpha = 0.01;
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

    void index_constraint(IndexSet& inds) override {
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

    void add_constraints(gravity::Model<>& NN, IndexSet& inds, gravity::param<>& w, gravity::var<>& x, gravity::var<int>& y) override {
        Constraint<> RDSum("ReduceSum");
        RDSum = x.in(inds["Out"]) - x.in(inds["In"]);
        NN.add(RDSum.in(inds["Constr"]) == 0);
    }

    Tensor *X, *Y; // Input and output
    std::vector<size_t> axes;
    int keepdims = 1;
};

class AveragePool : public Layer {
public:
    AveragePool(const onnx::NodeProto& node, Tensors& tensors): Layer(node, tensors) {
        operator_type = _average_pool;
        this->X = &tensors[node.input(0)];
        this->Y = &tensors[node.output(0)];

        if (const auto* auto_pad_attr = find_attribute("auto_pad", node)) {
            this->auto_pad = auto_pad_attr->s();
        }
        if (this->auto_pad != "NOTSET") {
            throw std::runtime_error("AveragePool: auto_pad " + this->auto_pad + " not implemented");
        }

        if (const auto* ceil_mode_attr = find_attribute("ceil_mode", node)) {
            this->ceil_mode = ceil_mode_attr->i();
        }

        if (const auto* count_include_pad_attr = find_attribute("count_include_pad", node)) {
            this->count_include_pad = count_include_pad_attr->i();
        }

        if (const auto* kernel_shape_attr = find_attribute("kernel_shape", node)) {
            this->kernel_shape = std::vector<size_t>(kernel_shape_attr->ints().begin(), kernel_shape_attr->ints().end());
        } else {
            throw std::runtime_error("AveragePool: kernel_shape not found");
        }

        if (const auto* pads_attr = find_attribute("pads", node)) {
            this->pads = std::vector<int>(pads_attr->ints().begin(), pads_attr->ints().end());
        } else {
            this->pads = std::vector<int>(this->kernel_shape.size() * 2, 0);
        }
        if (!std::all_of(this->pads.begin(), this->pads.end(), [](int i){return i == 0;})) {
            throw std::runtime_error("AveragePool: padding != 0 not supported");
        }

        if (const auto* strides_attr = find_attribute("strides", node)) {
            this->strides = std::vector<int>(strides_attr->ints().begin(), strides_attr->ints().end());
        } else {
            this->strides = std::vector<int>(this->kernel_shape.size(), 1);
        }
    }

    void add_parameters(gravity::param<>& w) const override {
        double scale = 1.0 / ((double)vecprod(this->kernel_shape));
        for (auto i = 0; i < this->Y->numel; i++) {
            w.add_val(this->Y->strkey(i) + "_scale", scale);
        }
    }


    std::vector<std::vector<std::string>> get_indices() const override {
        return {{"Out", "In"}, {"scale"}};
    }

    void index_constraint(IndexSet& inds) override {
        for (auto out_idx: ShapeIter(this->Y->shape)) {
            inds["Constr"].add(this->Y->strkey(out_idx));
            inds["Out"].add_ref(this->Y->strkey(out_idx));
            inds["scale"].add_ref(this->Y->strkey(out_idx) + "_scale");

            std::vector<size_t> tl;
            for (int ax = 2; ax < this->Y->ndims; ax++) {
                tl.push_back(out_idx.at(ax) * this->strides.at(ax-2));
            }

            for (auto window_idx: ShapeIter(this->kernel_shape)) {
                std::vector<size_t> x_feature = vecsum(tl, window_idx);
                std::vector<size_t> xindex = concat({out_idx.at(0), out_idx.at(1)}, x_feature);
                inds["In"].add_in_row(inds.row_id, this->X->strkey(xindex));
            }
            inds.row_id++;
        }
    }

    void add_constraints(gravity::Model<>& NN, IndexSet& inds, gravity::param<>& w, gravity::var<>& x, gravity::var<int>& y) override {
        Constraint<> AveragePool_("AveragePool");
        AveragePool_ = x.in(inds["Out"]) - (x.in(inds["In"]) * w.in(inds["scale"]));
        NN.add(AveragePool_.in(inds["Constr"]) == 0);
    }

    Tensor *X, *Y; // Input and output

    std::string auto_pad = "NOTSET";
    int ceil_mode = 0;
    int count_include_pad = 0;
    std::vector<size_t> kernel_shape;
    std::vector<int> pads;
    std::vector<int> strides;
};