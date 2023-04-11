#include <iostream>
#include <fstream>
#include <onnx.pb.h>
#include <vector>
#include "Layers.hpp"
#include <gravity/model.h>
#include <gravity/func.h>

void build_graph(const onnx::GraphProto& graph) {
    std::vector<Layer*> layers;
    Initializers initializers;

    // Parse initializers
    std::cout << "Initializers: " << std::endl;
    for (const auto& initializer : graph.initializer()) {
        auto param = parse_tensor(initializer);
        initializers[param._name] = param;
    }

    for (const auto& node : graph.node()) {
        if (node.op_type() == "Gemm") {
            layers.push_back(new GEMM(node, initializers));
        } else if (node.op_type() == "Relu") {
            layers.push_back(new ReLU(node, initializers));
        }
    }

    for (const auto& layer : layers) {
        layer->print();
    }

    gravity::param<float> input("input");
    input = {1, 2, 3, 4};
    HiddenStates hidden_states;
    hidden_states[layers[0]->inputs[0]] = input;

    for (const auto& layer : layers) {
        layer->forward(hidden_states);
    }

    auto output = hidden_states[(*layers.end())->outputs[0]];
    output.print();
}

int main() {
    std::fstream input("/mnt/trail_test/haydnj/vnn/mip_maker/simple.onnx", std::ios::in | std::ios::binary);
    onnx::ModelProto model;
    bool isSuccess = model.ParseFromIstream(&input);
    onnx::GraphProto graph = model.graph();

    // Print the graph.
    build_graph(graph);
    return 0;
}
