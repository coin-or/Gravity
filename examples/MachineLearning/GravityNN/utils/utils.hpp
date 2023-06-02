#pragma once

#include <sstream>
#include <vector>
#include <numeric>
#include <queue>
#include <thread>

template <typename T>
std::string print_vector(const std::vector<T>& v) {
    std::stringstream ss;
    ss << "[";
    for (size_t i = 0; i < v.size(); ++i) {
        ss << v[i];
        if (i != v.size() - 1) {
            ss << ", ";
        }
    }
    ss << "]";
    return ss.str();
}

template <typename T>
std::string vec_to_index(const std::vector<T>& v) {
    std::string out;
    out.reserve(v.size()*2);
    for (auto& val : v) {
        out += "," + std::to_string(val);
    }
    return out;
}

template <typename T>
T vecprod (const std::vector<T>& v) {
    return std::accumulate(v.begin(), v.end(), 1, std::multiplies<T>());
}


onnx::GraphProto _open_file(const std::string& onnx_path) {
    // Check if file exists
    std::ifstream f(onnx_path.c_str());
    if (!f.good()) {
        throw std::runtime_error("File " + onnx_path + " does not exist.");
    }

    std::ifstream onnx_file(onnx_path, std::ios::binary);
    std::string onnx_str((std::istreambuf_iterator<char>(onnx_file)), std::istreambuf_iterator<char>());

    onnx::ModelProto model;
    bool isSuccess = model.ParseFromString(onnx_str);
    if (!isSuccess) {
        throw std::runtime_error("Failed to parse onnx file.");
    }
    onnx::GraphProto graph = model.graph();

    return graph;
}

template<typename T>
std::vector<T> apply_permutation(const std::vector<T>& v, const std::vector<T>& indices) {
    std::vector<T> v2;
    v2.reserve(v.size());
    for (size_t i = 0; i < v.size(); i++) {
        v2.push_back(v[indices[i]]);
    }
    return v2;
}

/*
    * Extract a subgraph from an ONNX graph.
    * The subgraph is defined by the path from the final_node to the input nodes.
    * If final_node is empty, return all layers.
*/
std::set<std::string> subgraph_extraction(onnx::GraphProto& graph, std::string final_node) {
    std::vector<std::string> all_layers;
    for (auto& node: graph.node()) {
        all_layers.push_back(node.name());
    }

    if (final_node.empty()) {
        return std::set<std::string>(all_layers.begin(), all_layers.end());
    }

    std::map<std::string, std::string> output_to_layer;
    std::map<std::string, int> layer_to_index;

    for (auto i = 0; i < graph.node_size(); ++i) {
        auto node = graph.node(i);
        layer_to_index[node.name()] = i;
        for (const auto& output : node.output()) {
            output_to_layer[output] = node.name();
        }
    }

    std::queue<std::string> queue;
    queue.push(final_node);

    std::set<std::string> subgraph;
    std::vector<std::string> order;
    while (!queue.empty()) {
        std::string node = queue.front();
        queue.pop();
        subgraph.insert(node);
        order.push_back(node);

        auto layer = graph.node(layer_to_index[node]);
        for (const auto& input : layer.input()) {
            if (output_to_layer.find(input) != output_to_layer.end()) {
                queue.push(output_to_layer[input]);
            }
        }
    }

    // Reverse order
    std::reverse(order.begin(), order.end());

    // Print all layers
    std::cout << "All layers: ";
    std::cout << print_vector(all_layers) << std::endl;
    std::cout << "Subgraph: ";
    std::cout << print_vector(order) << std::endl;

    // get unused layers
    std::vector<std::string> unused_layers;
    std::set_difference(all_layers.begin(), all_layers.end(), subgraph.begin(), subgraph.end(),
                        std::inserter(unused_layers, unused_layers.begin()));

    std::cout << "Unused layers: ";
    std::cout << print_vector(unused_layers) << std::endl;

    return subgraph;
}

bool ends_with(std::string const& value, std::string const & ending)
{
    if (ending.size() > value.size()) {
        return false;
    }
    return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

size_t get_num_threads() {
    return std::thread::hardware_concurrency();
}

std::string ftostr(float v) {
    if (v == std::numeric_limits<float>::lowest()) {
        return "-∞";
    }
    if (v == std::numeric_limits<float>::max()) {
        return "∞";
    }
    return std::to_string(v);
}