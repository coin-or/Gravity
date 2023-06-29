#pragma once

#include <network/Layers/LayerBase.hpp>

using namespace gravity;

// These are operators that do not need to generate constraints / add params
std::set<OType> PURE_NOOPS = {_flatten, _squeeze, _reshape, _split, _transpose, _gather};
// These are operators that do
std::set<OType> IMPURE_NOOPS = {_concat};

class NoOp : public Layer {
public:
    NoOp(const onnx::NodeProto& node, Tensors& tensors): Layer(node, tensors) {
        this->folded = true;

        // operator_type = _noop;
        this->X = &tensors[node.input(0)];
        this->Y = &tensors[node.output(0)];
    }

    std::vector<std::vector<std::string>> get_indices() const override {
        return {{"In", "Out"}, {}};
    }

    void set_bounds(gravity::param<>& x_lb, gravity::param<>& x_ub) override {}
    void index_hidden_states(indices& hidden_states, indices& y_ids) override {}
    void add_constraints(gravity::Model<>& NN, IndexSet& inds, gravity::param<>& w, gravity::var<>& x, gravity::var<int>& y) override {}
    void index_constraint(IndexSet& inds) override {}

    virtual void remap_indices(IndexContainer& inds) {
        for(auto j = 0; j < this->X->numel;j++){
            inds.add_remap(this->X->strkey(j), this->Y->strkey(j));
        }
    }

    Tensor *X, *Y; // Input and output
};

class Flatten: public NoOp {
public:
    Flatten(const onnx::NodeProto& node, Tensors& tensors): NoOp(node, tensors) {
        this->operator_type = _flatten;
    }
};


class Reshape: public NoOp {
public:
    Reshape(const onnx::NodeProto& node, Tensors& tensors): NoOp(node, tensors) {
        this->operator_type = _reshape;
    }
};

class Squeeze: public NoOp {
public:
    Squeeze(const onnx::NodeProto& node, Tensors& tensors): NoOp(node, tensors) {
        this->operator_type = _squeeze;
    }
};

class Split : public NoOp {
public:
    Split(const onnx::NodeProto& node, Tensors& tensors): NoOp(node, tensors) {
        this->operator_type = _split;

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
            /*
                Split a tensor into a list of tensors, along the specified ‘axis’. Lengths of the parts can be specified using argument ‘split’.
                Otherwise, the tensor is split to equal sized parts.
            */
            size_t split_size = this->X->shape[this->axis] / this->outputs.size();
            for (size_t i = 0; i < this->outputs.size(); i++) {
                this->split.push_back(split_size);
            }
        }
    }

    void remap_indices(IndexContainer& inds) override {
        size_t cur_axis_idx = 0;
        for (auto out: this->outputs) {
            for (auto y_ind: ShapeIter(out->shape)) {
                auto x_ind = y_ind;
                x_ind.at(this->axis) += cur_axis_idx;
                inds.add_remap(this->X->strkey(x_ind), out->strkey(y_ind));
            }
            cur_axis_idx += out->shape.at(this->axis);
        }
    }

    size_t axis;
    std::vector<size_t> split;
};

class Concat : public NoOp {
public:
    Concat(const onnx::NodeProto& node, Tensors& tensors): NoOp(node, tensors) {
        this->operator_type = _concat;
        this->Y = &tensors[node.output(0)];

        if (const auto* axis_attr = find_attribute("axis", node)) {
            int ax = axis_attr->i();
            if (ax < 0) {
                ax += this->inputs.at(0)->shape.size();
            }
            this->axis = ax;
        } else {
            throw std::runtime_error("Concat: axis attribute not found.");
        }
    }

    std::vector<std::vector<std::string>> get_indices() const override {
        return {{"pOut"}, {"pIn"}};
    }

    void remap_indices(IndexContainer& inds) override {
        size_t cur_axis_idx = 0;
        for (auto inp: this->inputs) {
            if (inp->is_initializer) {
                cur_axis_idx += inp->shape[this->axis];
                continue;
            }

            for (auto index: ShapeIter(inp->shape)) {
                auto yindex = index;
                yindex.at(this->axis) += cur_axis_idx;
                inds.add_remap(inp->strkey(index), this->Y->strkey(yindex));
            }
            cur_axis_idx += inp->shape[this->axis];
        }
    }

    void index_hidden_states(indices& hidden_states, indices& y_ids) override {
        size_t cur_axis_idx = 0;
        for (auto inp: this->inputs) {
            if (!inp->is_initializer) {
                cur_axis_idx += inp->shape[this->axis];
                continue;
            }
            for (auto index: ShapeIter(inp->shape)) {
                auto yindex = index;
                yindex.at(this->axis) += cur_axis_idx;
                hidden_states.add(this->Y->strkey(yindex));
            }
            cur_axis_idx += inp->shape[this->axis];
        }

    }

    void index_constraint(IndexSet& inds) override {
        size_t cur_axis_idx = 0;
        for (auto inp: this->inputs) {
            if (!inp->is_initializer) {
                cur_axis_idx += inp->shape[this->axis];
                continue;
            }
            for (auto index: ShapeIter(inp->shape)) {
                auto yindex = index;
                yindex.at(this->axis) += cur_axis_idx;

                inds["Constr"].add(this->Y->strkey(yindex));
                inds["pIn"].add_ref(inp->strkey(index));
                inds["pOut"].add_ref(this->Y->strkey(yindex));
            }
            cur_axis_idx += inp->shape[this->axis];
        }
    }

    void add_parameters(gravity::param<>& w) const override {
        for (auto inp: this->inputs) {
            if (inp->is_initializer) {
                inp->add_params(w);
            }
        }
    }

    void add_constraints(gravity::Model<>& NN, IndexSet& inds, gravity::param<>& w, gravity::var<>& x, gravity::var<int>& y) override {
        Constraint<> Concat_P(this->lname() + "_Concat_Param");
        Concat_P = x.in(inds["pOut"]) - w.in(inds["pIn"]);
        NN.add(Concat_P.in(inds["Constr"]) == 0);
    }

    size_t axis;
    Tensor* Y;
};

class Transpose : public NoOp {
public:
    Transpose(const onnx::NodeProto& node, Tensors& tensors): NoOp(node, tensors) {
        this->operator_type = _transpose;
        this->X = &tensors[node.input(0)];
        this->Y = &tensors[node.output(0)];

        if (const auto* attr = this->find_attribute("perm", node)) {
            this->perm = std::vector<size_t>(attr->ints().begin(), attr->ints().end());
        } else {
            throw std::runtime_error("Transpose: perm attribute not found.");
        }
    }

    void remap_indices(IndexContainer& inds) override {
        for (auto xindex: ShapeIter(this->X->shape)) {
            auto yindex = apply_permutation(xindex, this->perm);
            inds.add_remap(this->X->strkey(xindex), this->Y->strkey(yindex));
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
    }

    void remap_indices(IndexContainer& inds) override {
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

        // Loop over every point, unflatten index, and make sure it is in the slice
        size_t out_idx = 0;
        for (auto i = 0; i < this->X->numel; i++) {
            auto xunflat = this->X->unflatten_index(i);
            // Begin checking if this point is in the slice
            std::vector<bool> in_slice(this->X->ndims, false);
            for (auto dim = 0; dim < eff_start.size(); dim++) {
                // Check if xuflat[dim] is in range(start, end, step)
                for (auto tmp = eff_start[dim]; tmp < eff_end[dim]; tmp += eff_step[dim]) {
                    if (tmp == xunflat[dim]) {
                        in_slice[dim] = true;
                        break;
                    }
                }
                if (!in_slice[dim]) {
                    break;
                }
            }
            // If all dimensions are in the slice, add this point to the constraint
            if (std::all_of(in_slice.begin(), in_slice.end(), [](bool v) { return v; })) {
                inds.add_remap(this->X->strkey(i), this->Y->strkey(out_idx));
                out_idx += 1;
            }
        }
    }

    Tensor *X, *Y;
    std::vector<int64_t> starts;
    std::vector<int64_t> ends;
    std::vector<int64_t> axes;
    std::vector<int64_t> steps;
};

class Gather : public NoOp {
public:
    Gather(const onnx::NodeProto& node, Tensors& tensors): NoOp(node, tensors) {
        this->operator_type = _gather;

        this->X = &tensors[node.input(0)];
        this->Y = &tensors[node.output(0)];

        if (const auto* axis = find_attribute("axis", node)) {
            int64_t i = axis->i();
            if (i < 0) {
                i += this->X->ndims;
            }
            this->axis = i;
        }

        // Pull out the indices
        auto& ind_ten = tensors[node.input(1)];
        for (auto v : ind_ten.get_int_data()) {
            if (v < 0) {
                v += this->X->shape[this->axis];
            }
            this->indices.push_back(v);
        }
    }

    void remap_indices(IndexContainer& inds) override {
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
        size_t ind_ptr = 0;
        size_t out_numel = 0;
        while (ind_ptr < this->indices.size()) {
            for (auto i = 0; i < this->X->numel; i++) {
                auto xunflat = this->X->unflatten_index(i);
                if (xunflat.at(this->axis) == this->indices.at(ind_ptr)) {
                    inds.add_remap(this->X->strkey(i), this->Y->strkey(out_numel));
                    out_numel += 1;
                }
            }
            ind_ptr += 1;
        }
    }

    Tensor *X, *Y; // Input and output
    size_t axis = 0;

    std::vector<size_t> indices;
};