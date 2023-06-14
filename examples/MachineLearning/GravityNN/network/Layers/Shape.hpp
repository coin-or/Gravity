#pragma once

#include <network/Layers/LayerBase.hpp>

using namespace gravity;

class NoOp : public Layer {
public:
    NoOp(const onnx::NodeProto& node, Tensors& tensors): Layer(node, tensors) {
        operator_type = _noop;
        this->X = &tensors[node.input(0)];
        this->Y = &tensors[node.output(0)];
    }

    std::vector<std::vector<std::string>> get_indices() const override {
        return {{"In", "Out"}, {}};
    }

    void index_constraint(IndexSet& inds) override {
        for(auto j = 0; j < this->X->numel;j++){
            inds["Constr"].add(this->Y->strkey(j));
            inds["In"].add_ref(this->X->strkey(j));
            inds["Out"].add_ref(this->Y->strkey(j));
        }
    }

    void add_constraints(gravity::Model<>& NN, IndexSet& inds, gravity::param<>& w, gravity::var<>& x, gravity::var<int>& y) override {
        Constraint<> NoOp("NoOp");
        NoOp = x.in(inds["Out"]) - x.in(inds["In"]);
        NN.add(NoOp.in(inds["Constr"]) == 0);
    }

    Tensor *X, *Y; // Input and output
};

class Split : public NoOp {
public:
    Split(const onnx::NodeProto& node, Tensors& tensors): NoOp(node, tensors) {
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

    void index_constraint(IndexSet& inds) override {
        size_t cur_axis_idx = 0;
        for (auto out: this->outputs) {
            for (size_t out_idx = 0; out_idx < out->numel; out_idx++) {
                auto inp_vec_idx = out->unflatten_index(out_idx);
                inp_vec_idx[this->axis] += cur_axis_idx;
                auto inp_idx = this->X->flatten_index(inp_vec_idx);

                inds["Constr"].add(out->strkey(out_idx));
                inds["In"].add_ref(this->X->strkey(inp_idx));
                inds["Out"].add_ref(out->strkey(out_idx));
            }

            cur_axis_idx += out->shape[this->axis];
        }
    }

    size_t axis;
    std::vector<size_t> split;
};

class Concat : public Layer {
public:
    Concat(const onnx::NodeProto& node, Tensors& tensors): Layer(node, tensors) {
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
        return {{"hIn", "hOut", "pOut"}, {"pIn"}};
    }

    void index_constraint(IndexSet& inds) override {
        size_t cur_axis_idx = 0;
        for (auto inp: this->inputs) {
            for (auto index: ShapeIter(inp->shape)) {
                auto yindex = index;
                yindex.at(this->axis) += cur_axis_idx;

                if (inp->is_initializer) {
                    inds["ConstrB"].add(this->Y->strkey(yindex));
                    inds["pIn"].add_ref(inp->strkey(index));
                    inds["pOut"].add_ref(this->Y->strkey(yindex));
                } else {
                    inds["Constr"].add(this->Y->strkey(yindex));
                    inds["hIn"].add_ref(inp->strkey(index));
                    inds["hOut"].add_ref(this->Y->strkey(yindex));
                }
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
        Constraint<> Concat_H("Concat_Hidden");
        Concat_H = x.in(inds["hOut"]) - x.in(inds["hIn"]);
        NN.add(Concat_H.in(inds["Constr"]) == 0);

        Constraint<> Concat_P("Concat_Param");
        Concat_P = x.in(inds["pOut"]) - w.in(inds["pIn"]);
        NN.add(Concat_P.in(inds["ConstrB"]) == 0);
    }

    size_t axis;
    Tensor* Y;
};

class Transpose : public NoOp {
public:
    Transpose(const onnx::NodeProto& node, Tensors& tensors): NoOp(node, tensors) {
        this->X = &tensors[node.input(0)];
        this->Y = &tensors[node.output(0)];

        if (const auto* attr = this->find_attribute("perm", node)) {
            this->perm = std::vector<size_t>(attr->ints().begin(), attr->ints().end());
        } else {
            throw std::runtime_error("Transpose: perm attribute not found.");
        }
    }

    void index_constraint(IndexSet& inds) override {
        for (auto yindex: ShapeIter(this->Y->shape)) {
            auto xindex = apply_permutation(yindex, this->perm);

            inds["Constr"].add(this->Y->strkey(yindex));
            inds["In"].add_ref(this->X->strkey(xindex));
            inds["Out"].add_ref(this->Y->strkey(yindex));
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

    void index_constraint(IndexSet& inds) override {
        /*
            Slice uses the starts, ends, axes and steps inputs to select a sub-tensor of its input data tensor.
            An effective start[i], end[i], and step[i] must be computed for each i in [0, ... r-1] where r = rank(input) as follows:
            If axes are omitted, they are set to [0, ..., r-1]. If steps are omitted, they are set to [1, ..., 1] of length len(starts)
            The effective values are initialized as start[i] = 0, end[i] = dims[i] where dims are the dimensions of input and step[i] = 1.
            All negative elements of axes are made non-negatve by adding r to them, where r =rank(input).
            All negative values in starts[i] and ends[i] have dims[axes[i]] added to them, where dims are the dimensions of input. Then start[axes[i]] is the adjusted starts[i] is clamped into the range [0, dims[axes[i]]] for positive stepping and [0, dims[axes[i]]-1] for negative stepping.
            The clamping for the adjusted ends[i] depends on the sign of steps[i] and must accommodate copying 0 through dims[axes[i]] elements, so for positive stepping end[axes[i]] is clamped to [0, dims[axes[i]]], while for negative stepping it is clamped to [-1, dims[axes[i]]-1].
            Finally, step[axes[i]] = steps[i].
            For slicing to the end of a dimension with unknown size, it is recommended to pass in INT_MAX when slicing forward and ‘INT_MIN’ when slicing backward.
        */

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
                inds["Constr"].add(this->Y->strkey(out_idx));
                inds["In"].add_ref(this->X->strkey(i));
                inds["Out"].add_ref(this->Y->strkey(out_idx));
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

    void index_constraint(IndexSet& inds) override {
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
                    inds["Constr"].add(this->Y->strkey(out_numel));
                    inds["In"].add_ref(this->X->strkey(i));
                    inds["Out"].add_ref(this->Y->strkey(out_numel));
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