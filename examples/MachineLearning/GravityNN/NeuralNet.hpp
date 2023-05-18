#pragma once

#include <string>
#include <onnx.pb.h>
#include <vector>
#include "Layers.hpp"

std::set<std::string> noops = {"Flatten", "Reshape", "Squeeze"};

class NeuralNet {
public:
    NeuralNet(const std::string& onnx_path) {
        onnx::GraphProto graph = _open_file(onnx_path);
        if (graph.input_size() > 1) {
            throw std::runtime_error("Network has more than one input. Not supported.");
        }

        this->tensors = get_tensors(graph);

        this->input_numel = 0;
        for (const auto& input : graph.input()) {
            Layer* inp_layer = new Input(input, tensors);
            this->layers.push_back(inp_layer);
            this->input_numel += inp_layer->outputs[0]->numel;
        }

        this->indices = IndexContainer();

        for (const auto& node : graph.node()) {
            if (node.op_type() == "Constant") {
                // We've stuffed all constants into the tensors map
                continue;
            }

            Layer* node_ptr;
            if ((node.op_type() == "Gemm") || (node.op_type() == "MatMul")) {
                node_ptr = new GEMM(node, this->tensors);
            } else if (node.op_type() == "Relu") {
                node_ptr = new Relu(node, this->tensors);
            } else if (node.op_type() == "Conv") {
                node_ptr = new Conv(node, this->tensors);
            } else if (noops.count(node.op_type())) {
                node_ptr = new NoOp(node, this->tensors);
            } else if (node.op_type() == "Split") {
                node_ptr = new Split(node, this->tensors);
            } else if (node.op_type() == "Concat") {
                node_ptr = new Concat(node, this->tensors);
            } else if (node.op_type() == "Add") {
                node_ptr = new Add(node, this->tensors);
            } else if (node.op_type() == "Sub") {
                node_ptr = new Sub(node, this->tensors);
            } else if (node.op_type() == "Cos") {
                node_ptr = new Cos(node, this->tensors);
            } else if (node.op_type() == "Sin") {
                node_ptr = new Sin(node, this->tensors);
            } else if (node.op_type() == "Neg") {
                node_ptr = new Neg(node, this->tensors);
            } else if (node.op_type() == "Pow") {
                node_ptr = new Pow(node, this->tensors);
            } else if (node.op_type() == "Mul") {
                node_ptr = new Mul(node, this->tensors);
            } else if (node.op_type() == "Transpose") {
                node_ptr = new Transpose(node, this->tensors);
            } else if (node.op_type() == "Slice") {
                node_ptr = new Slice(node, this->tensors);
            } else if (node.op_type() == "Div") {
                node_ptr = new Div(node, this->tensors);
            } else if (node.op_type() == "Gather") {
                node_ptr = new Gather(node, this->tensors);
            } else {
                throw std::runtime_error("Unsupported operator " + node.op_type());
            }

            this->layers.push_back(node_ptr);
        }
    }

    /*
        Adds all required indices to the model.
        This includes hidden states, binaries, weight indices, etc.
    */
    void build_indexing() {
        // Add index sets for each layer
        for (auto l: this->layers) {
            auto optype = l->operator_type;
            this->indices.add(l->operator_type, l->get_indices());
        }

        // First, index all hidden states
        for (auto l: this->layers) {
            l->index_hidden_states(this->indices.hidden_states, this->indices.y_ids);
        }

        // Index parameters
        for (auto l: this->layers) {
            l->add_parameters(this->indices.w);
        }
    }

    // Builds constraints for each layer
    void build_constraints() {
        for (auto l: this->layers) {
            l->build_constraint(this->indices(l->operator_type));
        }
    }

    void set_bounds(gravity::param<>& x_lb, gravity::param<>& x_ub) {
        for (auto l: this->layers) {
            if (l->operator_type != _input) {
                continue;
            }
            for (size_t o = 0; o < l->outputs.size(); o++) {
                for(auto j = 0; j < l->outputs[o]->numel; j++){
                    x_lb.set_val(l->outputs[o]->strkey(j), -1.0);
                    x_ub.set_val(l->outputs[o]->strkey(j),  1.0);
                }
            }
        }
        for (auto l: this->layers) {
            // Enforce lb of 0 for relu
            if (l->operator_type == _relu) {
                for(auto j = 0; j < l->outputs[0]->numel;j++){
                    x_lb.set_val(l->outputs[0]->strkey(j), 0.0);
                }
            }

            // Set provided bounds, skip if not provided
            if (l->lowers.size() == 0) {
                continue;
            }

            for (size_t o = 0; o < l->outputs.size(); o++) {
                for(auto j = 0; j < l->lowers[o]->numel; j++){
                    x_lb.set_val(l->outputs[o]->strkey(j), (*l->lowers[o])(j));
                    x_ub.set_val(l->outputs[o]->strkey(j), (*l->uppers[o])(j));
                }
            }
        }
    }

    void initialize_state(gravity::var<>& x, gravity::var<int>& y) {
        for (auto l: this->layers) {
            // Set provided forward values, skip if not provided
            if (l->forward_values.size() == 0) {
                continue;
            }

            for (size_t o = 0; o < l->outputs.size(); o++) {
                for(auto j = 0; j < l->forward_values[o]->numel; j++){
                    x.param<double>::set_val(l->outputs[o]->strkey(j), (*l->forward_values[o])(j));
                    if (l->operator_type == _relu) {
                        y.param<int>::set_val(l->outputs[o]->strkey(j), (int)((*l->forward_values[o])(j) > 0));
                    }
                }
            }
        }
    }

    void add_constraints(gravity::Model<>& NN, gravity::var<>& x, gravity::var<int>& y, IndexContainer& indices) {
        // Add constraints. Only add constraints for each operator type once.
        std::set<OType> visited;
        for (auto l: this->layers) {
            if (visited.find(l->operator_type) != visited.end()) {
                continue;
            }
            visited.insert(l->operator_type);
            l->add_constraints(NN, indices(l->operator_type), indices.w, x, y);
        }
    }

    Tensors tensors;
    IndexContainer indices;
    size_t input_numel;

    std::vector<Layer*> layers;
};
