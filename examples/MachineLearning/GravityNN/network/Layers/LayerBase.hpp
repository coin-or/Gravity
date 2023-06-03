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
        for (const auto& input : node.input()) {
            this->inputs.push_back(&tensors.at(input));
        }
        for (const auto& output : node.output()) {
            this->output_names.push_back(output);
            this->outputs.push_back(&tensors.at(output));
        }

        this->name = node.name();
    }

    // Input layer constructor
    Layer(const onnx::ValueInfoProto& input_node, Tensors& tensors) {
        this->name = input_node.name();
        this->operator_type = _input;
        this->output_names.push_back(input_node.name());
        this->outputs.push_back(&tensors.at(input_node.name()));
    }

    virtual ~Layer() = default;

    virtual void build_constraint(IndexSet& inds) = 0;
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

    std::vector<std::string> output_names;
    std::vector<Tensor*> outputs;
    std::vector<Tensor*> inputs;
};