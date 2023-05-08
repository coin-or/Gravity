#pragma once

#include <fstream>
#include <onnx.pb.h>
#include <vector>
#include "utils.hpp"
#include <gravity/param.h>

class Tensor {
public:
    Tensor() {}

    Tensor(const onnx::TensorProto& tensor) {
        this->name = tensor.name();
        this->is_initializer = true;

        this->_set_shape(std::vector<size_t>(tensor.dims().begin(), tensor.dims().end()));
        if (!tensor.raw_data().empty()) {
            const void* raw_data = tensor.raw_data().data();
            data.resize(this->numel);
            std::memcpy(data.data(), raw_data, this->numel * sizeof(float));
        } else if (!tensor.float_data().empty()) {
            this->data = std::vector<float>(tensor.float_data().begin(), tensor.float_data().end());
        } else {
            throw std::runtime_error("Tensor " + tensor.name() + " has data in neither raw_data nor float_data.");
        }
    }

    Tensor(const onnx::ValueInfoProto& vinfo) {
        this->name = vinfo.name();
        this->is_initializer = false;

        std::vector<size_t> shape;
        for (auto dim : vinfo.type().tensor_type().shape().dim()) {
            shape.push_back(dim.dim_value());
        }
        this->_set_shape(shape);
    }

    void _transpose() {
        if (this->shape.size() != 2) {
            throw std::runtime_error("Cannot transpose tensor with shape " + std::to_string(this->shape.size()));
        }
        if (!this->is_initializer) {
            throw std::runtime_error("Cannot transpose non-initializer tensor");
        }

        // Tranpose data if needed, it is stored row-major
        std::vector<float> temp_data = this->data;
        for (size_t i = 0; i < this->shape[0]; i++) {
            for (size_t j = 0; j < this->shape[1]; j++) {
                this->data[j * this->shape[0] + i] = temp_data[i * this->shape[1] + j];
            }
        }

        // Swap shape
        std::swap(this->shape[0], this->shape[1]);
    }

    float operator()(size_t i) const {
        if (!this->is_initializer) {
            throw std::runtime_error("Reading from non-initializer tensor. Perhaps you're assuming this tensor is a weight when it's actually an output of a previous layer?");
        }
        return this->data.at(i);
    }

    void add_params(gravity::param<>* p) const {
        if (!this->is_initializer) {
            throw std::runtime_error("Reading from non-initializer tensor. Perhaps you're assuming this tensor is a weight when it's actually an output of a previous layer?");
        }

        for (size_t i = 0; i < this->numel; i++) {
            p->add_val(this->strkey(i), this->data[i]);
        }
    }

    std::string strkey(size_t idx) const {
        size_t flat = idx;
        this->_boundcheck(flat);
        return this->name + "," + std::to_string(flat);
    }

    std::string strkey(size_t i, size_t j) const {
        size_t flat = i * this->shape[1] + j;
        this->_boundcheck(flat);
        return this->name + "," + std::to_string(flat);
    }
    
    std::string strkey(size_t i, size_t j, size_t k) const {
        size_t flat = i * this->shape[1] * this->shape[2] + j * this->shape[2] + k;
        this->_boundcheck(flat);
        return this->name + "," + std::to_string(flat);
    }

    std::string strkey(size_t i, size_t j, size_t k, size_t l) const {
        size_t flat = i * this->shape[1] * this->shape[2] * this->shape[3] + j * this->shape[2] * this->shape[3] + k * this->shape[3] + l;
        this->_boundcheck(flat);
        return this->name + "," + std::to_string(flat);
    }

    void _set_shape(const std::vector<size_t>& shape) {
        this->shape = shape;
        this->numel = vecprod(this->shape);
        this->ndims = this->shape.size();
    }

    void _boundcheck(size_t flat_idx) const {
        if (flat_idx >= this->numel) {
            throw std::runtime_error("Index " + std::to_string(flat_idx) + " out of bounds for tensor of size " + std::to_string(this->numel) + ".");
        }
    }

    std::string name;
    bool is_initializer;

    std::vector<size_t> shape;
    size_t numel;
    size_t ndims;

private:
    std::vector<float> data;
};

typedef std::map<std::string, Tensor> Tensors;

Tensors get_tensors(onnx::GraphProto& graph) {
    Tensors tensors;

    // Tensor with shape/metadata only
    for (const auto& vinfo : graph.value_info()) {
        tensors[vinfo.name()] = Tensor(vinfo);
    }

    // Tensors with data
    for (const auto& initializer : graph.initializer()) {
        tensors[initializer.name()] = Tensor(initializer);
    }

    // Output tensors
    for (const auto& output : graph.output()) {
        tensors[output.name()] = Tensor(output);
    }

    // Input tensors
    for (const auto& input : graph.input()) {
        tensors[input.name()] = Tensor(input);
    }

    return tensors;
}