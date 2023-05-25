#pragma once

#include <string>
#include <onnx.pb.h>
#include <vector>
#include <network/Layers.hpp>

std::set<std::string> noops = {"Flatten", "Reshape", "Squeeze"};

class NeuralNet {
public:
    NeuralNet(const std::string& onnx_path, std::string final_node = "") {
        onnx::GraphProto graph = _open_file(onnx_path);
        this->layer_names = subgraph_extraction(graph, final_node);

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

            // If node is not in requested subgraph, skip it
            if (this->layer_names.count(node.name()) == 0) {
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
            } else if (node.op_type() == "Clip") {
                node_ptr = new Clip(node, this->tensors);
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
            for (auto o: l->outputs) {
                for(auto j = 0; j < o->numel; j++){
                    x_lb.set_val(o->strkey(j), o->lb.at(j));
                    x_ub.set_val(o->strkey(j), o->ub.at(j));

                    if (l->operator_type == _clip) {
                        auto clip = static_cast<Clip*>(l);
                        x_lb.set_val(o->strkey(j), clip->min);
                        x_ub.set_val(o->strkey(j), clip->max);
                    }

                    if (l->operator_type == _relu) {
                        x_lb.set_val(o->strkey(j), 0.0);
                    }

                    if (l->operator_type == _input) {
                        x_lb.set_val(o->strkey(j), std::max(-10.0f, o->lb.at(j)));
                        x_ub.set_val(o->strkey(j), std::min( 10.0f, o->ub.at(j)));
                    }
                }
            }
        }
    }

    void initialize_state(gravity::var<>& x, gravity::var<int>& y) {
        for (auto l: this->layers) {
            for (auto o: l->outputs) {
                for(auto j = 0; j < o->numel; j++){
                    auto fv = o->forward.at(j);
                    auto key = o->strkey(j);

                    x.param<double>::set_val(key, fv);
                    if (l->operator_type == _relu) {
                        y.param<int>::set_val(key, (int)(fv > 0));
                    }

                    if (l->operator_type == _clip) {
                        auto clip = static_cast<Clip*>(l);
                        y.param<int>::set_val(key + "_min", 0);
                        y.param<int>::set_val(key + "_max", 0);
                        y.param<int>::set_val(key + "_eq",  0);
                        // 3 cases: x < min, min <= x <= max, x > max
                        if (fv < clip->min) {
                            y.param<int>::set_val(key + "_min", 1);
                        } else if (fv > clip->max) {
                            y.param<int>::set_val(key + "_max", 1);
                        } else {
                            y.param<int>::set_val(key + "_eq",  1);
                        }
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
    std::set<std::string> layer_names;
};
