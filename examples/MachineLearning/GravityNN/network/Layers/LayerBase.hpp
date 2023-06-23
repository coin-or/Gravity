#pragma once

#include <fstream>
#include <onnx.pb.h>
#include <vector>
#include <network/Tensor.hpp>
#include <network/IndexContainer.hpp>
#include <network/types.hpp>
#include <gravity/constraint.h>
#include <gravity/model.h>

using namespace gravity;

class Layer {
public:
    Layer(const onnx::NodeProto& node, Tensors& tensors) {
        this->name = node.name();
        this->opname = node.op_type();

        for (const auto& input : node.input()) {
            this->inputs.push_back(&tensors.at(input));
        }
        for (const auto& output : node.output()) {
            auto& out_ten = tensors.at(output);
            this->outputs.push_back(&out_ten);

            if (out_ten.lb.size() == 0) {
                throw std::runtime_error("Output tensor " + out_ten.name + " has no lower bound");
            }
            if (out_ten.ub.size() == 0) {
                throw std::runtime_error("Output tensor " + out_ten.name + " has no upper bound");
            }
        }

        this->name = node.name();
    }

    // Input layer constructor
    Layer(const onnx::ValueInfoProto& input_node, Tensors& tensors) {
        this->name = input_node.name();
        this->operator_type = _input;
        this->outputs.push_back(&tensors.at(input_node.name()));
    }

    virtual ~Layer() = default;

    virtual void index_constraint(IndexSet& inds) = 0;
    virtual std::vector<std::vector<std::string>> get_indices() const = 0;
    virtual void add_constraints(gravity::Model<>& NN, IndexSet& inds, gravity::param<>& w, gravity::var<>& x, gravity::var<int>& y) = 0;

    // Override this if you have to introduce auxiliary variables or node
    // uses int vars
    // i.e. node uses exponential or node is a relu
    virtual void index_hidden_states(indices& hidden_states, indices& y_ids) {
        for (auto output: this->outputs) {
            for (auto i = 0; i < output->numel; i++) {
                hidden_states.add(output->strkey(i));
            }
        }
    }

    // Override this if the layer has parameters
    virtual void add_parameters(gravity::param<>& w) const {}

    // Override this if the layer has aux vars that can be bounded
    // i.e. sigmoid produces an aux exp which has a lower bound of 0.0
    virtual void set_bounds(gravity::param<>& x_lb, gravity::param<>& x_ub) {
        for (size_t i = 0; i < this->outputs.size(); i++) {
            auto o = this->outputs.at(i);
            for(auto j = 0; j < o->numel; j++){
                auto key = o->strkey(j);
                auto lb  = this->get_lb(i, j);
                auto ub  = this->get_ub(i, j);

                x_lb.add_val(key, lb);
                x_ub.add_val(key, ub);
            }
        }
    }

    double get_lb(size_t output_index, size_t neuron_index) const {
        return std::max(this->outputs.at(output_index)->lb.at(neuron_index) - 1e-6, this->range_lower);
    }

    double get_ub(size_t output_index, size_t neuron_index) const {
        return std::min(this->outputs.at(output_index)->ub.at(neuron_index) + 1e-6, this->range_upper);
    }

    bool is_bounded(size_t output_index) {
        for (auto i = 0; i < this->outputs.at(output_index)->numel; i++) {
            if (this->get_lb(output_index, i) == HMIN || this->get_ub(output_index, i) == HMAX) {
                return false;
            }
        }
        return true;
    }

    const onnx::AttributeProto* find_attribute(const std::string& name, const onnx::NodeProto& node) const {
        for (const auto& attr : node.attribute()) {
            if (attr.name() == name) {
                return &attr;
            }
        }
        return nullptr;
    }


    // Name of the layer
    std::string name;
    OType operator_type;

    std::vector<Tensor*> outputs;
    std::vector<Tensor*> inputs;

    double range_lower = HMIN;
    double range_upper = HMAX;

    std::string opname;
};

class UnsupportedLayer: public Layer {
public:
    UnsupportedLayer(const onnx::NodeProto& node, Tensors& tensors) : Layer(node, tensors) {
        this->operator_type = _unsupported;
        std::cout << "WARNING: Unsupported layer: " << node.op_type() << " | " << node.name() << std::endl;
    }

    std::vector<std::vector<std::string>> get_indices() const override {
        throw std::runtime_error("Cannot get indices of unsupported layer " + this->name);
    }

    void index_constraint(IndexSet& inds) override {
        throw std::runtime_error("Cannot index unsupported layer " + this->name);

    }

    void add_constraints(gravity::Model<>& NN, IndexSet& inds, gravity::param<>& w, gravity::var<>& x, gravity::var<int>& y) override {
        throw std::runtime_error("Cannot add constraints to unsupported layer " + this->name);
    }
};

class Input : public Layer {
public:
    Input(const onnx::ValueInfoProto& node, Tensors& tensors): Layer(node, tensors) {
        operator_type = _input;
        this->opname = "Input";

        // If the input layer is unbounded, set the range to [-1, 1]
        for (auto o: this->outputs) {
            for (size_t i = 0; i < o->numel; i++) {
                if ((o->lb.at(i) == HMIN) || (o->ub.at(i) == HMAX)) {
                    this->range_lower = -1.0;
                    this->range_upper =  1.0;
                }
            }
        }

        if (this->range_lower == -1.0 || this->range_upper == 1.0) {
            std::cout << "WARNING: Input layer " << this->name << " is unbounded. Setting range to [-1, 1]" << std::endl;
        }
    }

    void index_constraint(IndexSet& inds) override {}
    std::vector<std::vector<std::string>> get_indices() const override {
        return {{}, {}};
    }

    void add_constraints(gravity::Model<>& NN, IndexSet& inds, gravity::param<>& w, gravity::var<>& x, gravity::var<int>& y) override {}
};
