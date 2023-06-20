#pragma once

#include <vector>
#include <queue>
#include <stack>
#include <onnx.pb.h>
#include <network/NeuralNet.hpp>
#include <network/types.hpp>

class Graph {
public:
    Graph(std::vector<Layer*> layers): layers(layers) {
        this->numVertices = this->layers.size();
    }

    const Layer* get_layer(std::string name) {
        for (auto& layer: this->layers) {
            if (layer->name == name) {
                return layer;
            }
        }
        throw std::runtime_error("Layer " + name + " not found.");
    }

    std::vector<Layer*> get_inputs(std::string layer_name) {
        auto layer = this->get_layer(layer_name);
        std::set<std::string> inputs;
        for (auto& inp: layer->inputs) {
            inputs.insert(inp->name);
        }

        std::vector<Layer*> input_layers;
        for (auto& layer: this->layers) {
            for (auto output: layer->outputs) {
                if (inputs.find(output->name) != inputs.end()) {
                    input_layers.push_back(layer);
                }
            }
        }
        return input_layers;
    }

    std::pair<Layer*, size_t> get_producer(Tensor* tensor) {
        for (auto& layer: this->layers) {
            for (size_t i = 0; i < layer->outputs.size(); i++) {
                if (layer->outputs.at(i) == tensor) {
                    return {layer, i};
                }
            }
        }
        throw std::runtime_error("Layer producing tensor " + tensor->name + " not found.");
    }

    std::vector<std::string> shortest_path(std::string startVertex, std::string endVertex) {
        std::map<std::string, std::string> predecessors;
        std::queue<std::string> queue;
        std::set<std::string> visited;

        queue.push(endVertex);
        visited.insert(endVertex);

        while (!queue.empty()) {
            auto current = queue.front();
            queue.pop();

            if (current == startVertex) {
                std::vector<std::string> path;
                while (current != endVertex) {
                    path.push_back(current);
                    current = predecessors[current];
                }
                path.push_back(endVertex);
                return path;
            }

            auto inputs = get_inputs(current);
            for (const auto& inp_layer: inputs) {
                if (visited.find(inp_layer->name) == visited.end()) {
                    queue.push(inp_layer->name);
                    visited.insert(inp_layer->name);
                    predecessors[inp_layer->name] = current;
                }
            }
        }

        throw std::runtime_error("No path found from " + startVertex + " to " + endVertex + ".");
    }

    void recursive_add_unbounded(std::vector<std::string>& added, const Layer* cur_layer) {
        // loop over inputs
        for (auto& inp: cur_layer->inputs) {
            if (inp->is_initializer) {
                continue;
            }

            auto tup = this->get_producer(inp);
            auto producer = tup.first;
            auto index = tup.second;

            // Is producer already in the subgraph?
            if (std::find(added.begin(), added.end(), producer->name) != added.end()) {
                continue;
            }

            // If this child has unbounded inputs, add it to the subgraph
            if (producer->is_bounded(index) == false) {
                added.push_back(producer->name);
                recursive_add_unbounded(added, producer);
            }
        }
    }

    int numVertices;
    std::vector<Layer*> layers;
};

/*
    * Extract a subgraph from an ONNX graph.
    * BFS from the final node to a start node for the shortest path between these nodes, including immediate parents of the intermediate and final nodes.
    * If start_node and final_node is empty, return all layers.
*/
std::vector<Layer*> subgraph_extraction(std::vector<Layer*> layers, std::string start_node, std::string final_node) {
    // If start_node and final_node is empty, return all layers
    if (start_node.empty() && final_node.empty()) {
        return layers;
    }

    // If start_node is empty, set it to the input
    if (start_node.empty()) {
        for (auto& layer: layers) {
            if (layer->operator_type == _input) {
                start_node = layer->name;
                break;
            }
        }
    }
    if (start_node.empty()) {
        throw std::runtime_error("No input layer found.");
    }

    // If final_node is empty, set it to the last node in the graph
    if (final_node.empty()) {
       final_node = layers.back()->name;
    }

    auto graph = Graph(layers);
    std::vector<std::string> path = graph.shortest_path(start_node, final_node);
    std::vector<std::string> full_graph = {start_node};
    // Add first level children to path
    for (size_t i = 1; i < path.size(); i++) {
        auto layer = graph.get_layer(path[i]);
        full_graph.push_back(layer->name);
        for (auto inp: graph.get_inputs(layer->name)) {
            if (std::find(path.begin(), path.end(), inp->name) == path.end()) {
                full_graph.push_back(inp->name);
            }
        }
    }

    // Add unbounded layers
    for (auto& layer_name: full_graph) {
        auto layer = graph.get_layer(layer_name);
        graph.recursive_add_unbounded(full_graph, layer);
    }

    // Remove duplicates
    full_graph.erase(std::unique(full_graph.begin(), full_graph.end()), full_graph.end());

    // Reorder according to original order in layers
    std::vector<Layer*> subgraph;
    for (auto& layer: layers) {
        if (std::find(full_graph.begin(), full_graph.end(), layer->name) != full_graph.end()) {
            subgraph.push_back(layer);
        }
    }

    return subgraph;
}