#pragma once
#include <string>
#include <onnx.pb.h>
#include <vector>
#include "Layers.hpp"

class NeuralNet {
public:
    NeuralNet(const std::string& onnx_path) {
        onnx::GraphProto graph = _open_file(onnx_path);
        if (graph.input_size() > 1) {
            throw std::runtime_error("Network has more than one input. Not supported.");
        }

        // Tensors with data
        for (const auto& initializer : graph.initializer()) {
            this->tensors[initializer.name()] = Tensor(initializer);
        }

        // Tensor with shape/metadata only
        for (const auto& vinfo : graph.value_info()) {
            this->tensors[vinfo.name()] = Tensor(vinfo);
        }
        for (const auto& input : graph.input()) {
            this->tensors[input.name()] = Tensor(input);
        }
        for (const auto& output : graph.output()) {
            this->tensors[output.name()] = Tensor(output);
        }

        for (const auto& node : graph.node()) {
            if (node.op_type() == "Gemm") {
                this->layers.push_back(new GEMM(node, tensors));
            } else if (node.op_type() == "Relu") {
                this->layers.push_back(new Relu(node, tensors));
            } else if (node.op_type() == "Conv") {
                this->layers.push_back(new Conv(node, tensors));
            } else {
                throw std::runtime_error("Unsupported operator " + node.op_type());
            }
        }

        for (auto i = 0; i < layers.size()-1; i++) {
            this->layers[i]->is_pre_activation = this->layers[i+1]->is_activation_func;
        }

        this->input_numel = vecprod(tensors[graph.input(0).name()].shape);
    }

    ~NeuralNet() {
        for (auto layer : layers) {
            delete layer;
        }
    }

    std::vector<Layer*> layers;
    Tensors tensors;
    size_t input_numel;
};