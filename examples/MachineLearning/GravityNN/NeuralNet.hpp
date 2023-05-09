#pragma once

#include <string>
#include <onnx.pb.h>
#include <vector>
#include "Layers.hpp"
#include <gravity/Net.h>

std::set<std::string> noops = {"Flatten"};

class NeuralNet: public Net {
public:
    NeuralNet(const std::string& onnx_path) {
        onnx::GraphProto graph = _open_file(onnx_path);
        if (graph.input_size() > 1) {
            throw std::runtime_error("Network has more than one input. Not supported.");
        }

        this->tensors = get_tensors(graph);

        for (const auto& input : graph.input()) {
            Layer* inp_layer = new Input(input, tensors);
            this->layers.push_back(inp_layer);
            this->add_node(inp_layer);
        }

        for (const auto& node : graph.node()) {
            Layer* node_ptr;
            if (node.op_type() == "Gemm") {
                node_ptr = new GEMM(node, this->tensors);
            } else if (node.op_type() == "Relu") {
                node_ptr = new Relu(node, this->tensors);
            } else if (node.op_type() == "Conv") {
                node_ptr = new Conv(node, this->tensors);
            } else if (noops.count(node.op_type())) {
                node_ptr = new NoOp(node, this->tensors);
            } else {
                throw std::runtime_error("Unsupported operator " + node.op_type());
            }

            this->layers.push_back(node_ptr);
            this->add_node(node_ptr);
        }

        this->_build_arcs();

        this->input_numel = vecprod(tensors[graph.input(0).name()].shape);
    }

    void index_hidden_states(indices& hidden_states, indices& y_ids) {
        for (auto l: this->layers) {
            for (auto output: l->outputs) {
                for (auto i = 0; i < output->numel; i++) {
                    hidden_states.add(output->strkey(i));
                    if (l->operator_type == _relu) {
                        y_ids.add(l->name + "," + to_string(i));
                    }
                }
            }
        }
    }

    void set_bounds(gravity::param<>& x_lb, gravity::param<>& x_ub) {
        for (auto l: this->layers) {
            // Enforce lb of 0 for relu
            if (l->operator_type == _relu) {
                auto relu = dynamic_cast<Relu*>(l);
                for(auto j = 0; j < relu->Y->numel;j++){
                    x_lb.set_val(relu->Y->name + "," + to_string(j), 0.0);
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
                        y.param<int>::set_val(l->name + "," + to_string(j), (int)((*l->forward_values[o])(j) > 0));
                    }
                }
            }
        }
    }

    void _build_arcs() {
        // build output name to layer map
        std::map<std::string, Layer*> name_to_layer;
        for (auto layer : this->layers) {
            for (auto output : layer->output_names) {
                name_to_layer.insert({output, layer});
            }
        }

        // build arcs
        for (auto dest : this->layers) {
            for (auto input : dest->input_names) {
                auto src = name_to_layer[input];
                std::cout << "Input: " << dest->name << std::endl;
                Arc* arc = new Arc(src, dest);
                this->add_arc(arc);
                arc->connect();
            }
        }
        this->print();
    }

    Tensors tensors;
    size_t input_numel;

    std::vector<Layer*> layers;
};
