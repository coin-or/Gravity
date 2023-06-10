#pragma once

#include <string>
#include <onnx.pb.h>
#include <vector>
#include <network/Layers/LayerBase.hpp>
#include <network/Layers/Linear.hpp>
#include <network/Layers/Binary.hpp>
#include <network/Layers/Unary.hpp>
#include <network/Layers/NonLinear.hpp>
#include <network/Layers/Shape.hpp>

std::set<std::string> noops = {"Flatten", "Reshape", "Squeeze"};

class NeuralNet {
public:
    NeuralNet(const std::string& onnx_path, std::string start_node = "",std::string final_node = "") {
        onnx::GraphProto graph = _open_file(onnx_path);
        this->layer_names = subgraph_extraction(graph, start_node, final_node);

        if (graph.input_size() > 1) {
            throw std::runtime_error("Network has more than one input. Not supported.");
        }

        this->tensors = get_tensors(graph);
        this->build_layers(graph);

        this->indices = IndexContainer();

        // Initialize MIP variables
        this->NN = Model<>(graph.name());

        this->x_lb = param<>("x_lb");
        this->x_ub = param<>("x_ub");

        this->x = var<>("x");
        this->y = var<int>("y", 0, 1);

        if (this->tensors.count("obj_spec_matrix") != 0) {
            this->obj_spec = &this->tensors.at("obj_spec_matrix");
        }

        if (this->tensors.count("obj_spec_values") != 0) {
            this->obj_val = &this->tensors.at("obj_spec_values");
        }

    }

    void build_layers(const onnx::GraphProto& graph) {
        this->input_numel = 0;
        for (const auto& input : graph.input()) {
            Layer* inp_layer = new Input(input, tensors);
            this->layers.push_back(inp_layer);
            this->input_numel += inp_layer->outputs[0]->numel;
        }

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
            if (node.op_type() == "Gemm") {
                node_ptr = new GEMM(node, this->tensors);
            } else if (node.op_type() == "MatMul") {
                node_ptr = new MatMul(node, this->tensors);
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
            } else if (node.op_type() == "Exp") {
                node_ptr = new Exp(node, this->tensors);
            } else if (node.op_type() == "Sigmoid") {
                node_ptr = new Sigmoid(node, this->tensors);
            } else if (node.op_type() == "BatchNormalization") {
                node_ptr = new BatchNorm(node, this->tensors);
            } else if (node.op_type() == "Softmax") {
                node_ptr = new Softmax(node, this->tensors);
            } else if (node.op_type() == "ReduceSum") {
                node_ptr = new ReduceSum(node, this->tensors);
            } else {
                throw std::runtime_error("Unsupported operator " + node.op_type());
            }

            this->layers.push_back(node_ptr);
        }
    }

    // Set obj_index to -1 if you want to use a custom objective, otherwise the index of the objective
    Model<>& build_model(int obj_index) {
        this->build_indexing();
        this->index_constraints();

        this->x_lb.in(this->indices.hidden_states);
        this->x_ub.in(this->indices.hidden_states);
        this->x_lb = std::numeric_limits<double>::lowest();
        this->x_ub = std::numeric_limits<double>::max();

        this->set_bounds(x_lb, x_ub);

        this->x.add_bounds(this->x_lb, this->x_ub);
        this->x.in(this->indices.hidden_states);
        this->y.in(this->indices.y_ids);
        this->initialize_state(x, y);

        this->NN.add(this->x);
        this->NN.add(this->y);
        this->add_constraints();

        this->set_objective(obj_index);

        return this->NN;
    }

    void set_objective(int obj_index) {
        if (obj_index < 0) {
            return;
        }

        if (this->obj_spec == nullptr || this->obj_val == nullptr) {
            throw std::runtime_error("Objective specification and value must be described in the ONNX model for this function.");
        }

        auto& spec = *this->obj_spec;
        auto& val = *this->obj_val;

        gravity::func<> obj = 0;
        for (size_t i = 0; i < spec.shape[1]; i++) {
            auto coeff = spec(spec.flatten_index({(size_t)obj_index, i}));
            std::string key = this->layers.back()->outputs[0]->strkey(i);
            obj += coeff * this->x(key);
        }
        // add the constant
        obj += val(obj_index);

        std::cout << "Objective: " << obj_index << std::endl;
        obj.print();

        NN.min(obj);
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
    void index_constraints() {
        for (auto l: this->layers) {
            l->index_constraint(this->indices(l->operator_type));
        }
    }

    void set_bounds(gravity::param<>& x_lb, gravity::param<>& x_ub) {
        for (auto l: this->layers) {
            l->set_bounds(x_lb, x_ub);
        }
    }

    void set_aux_bounds(const std::vector<Bound>& aux_bounds) {
        // use newbounds
        for (auto& v: aux_bounds) {
            if (this->layer_names.count(v.layer_name) == 0) {
                continue;
            }

            auto tensor_name = v.neuron_name.substr(0, v.neuron_name.find_last_of(","));
            size_t neuron_idx = std::stoi(v.neuron_name.substr(v.neuron_name.find_last_of(",") + 1));
            auto& ten = this->tensors.at(tensor_name);
            if (v.side == Side::LOWER) {
                ten.lb.at(neuron_idx) = std::max(v.value, ten.lb.at(neuron_idx));
            } else {
                ten.ub.at(neuron_idx) = std::min(v.value, ten.ub.at(neuron_idx));
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
                    } else if (l->operator_type == _clip) {
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

    void add_constraints() {
        // Add constraints. Only add constraints for each operator type once.
        std::set<OType> visited;
        for (auto l: this->layers) {
            if (visited.find(l->operator_type) != visited.end()) {
                continue;
            }
            visited.insert(l->operator_type);
            l->add_constraints(this->NN, this->indices(l->operator_type), this->indices.w, this->x, this->y);
        }
    }

    Tensors tensors;
    IndexContainer indices;
    size_t input_numel;

    Tensor* obj_spec = nullptr;
    Tensor* obj_val = nullptr;

    std::vector<Layer*> layers;
    std::set<std::string> layer_names;

    Model<> NN;
    param<> x_lb, x_ub;
    var<> x;
    var<int> y;
};
