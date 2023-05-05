#pragma once

#include <fstream>
#include <onnx.pb.h>
#include <vector>
#include "utils.hpp"

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

    void add_params(gravity::param<>* p, std::string& name) const {
        if (!this->is_initializer) {
            throw std::runtime_error("Reading from non-initializer tensor. Perhaps you're assuming this tensor is a weight when it's actually an output of a previous layer?");
        }
        for (size_t i = 0; i < this->numel; i++) {
            std::string key = name + this->unflatten(i);
            std::cout << "Adding param " << key << std::endl;
            p->add_val(key, this->data[i]);
        }
    }

    std::string unflatten(size_t idx) const {
        if (this->ndims == 1) {
            return "," + std::to_string(idx);
        }

        std::vector<size_t> indices = std::vector<size_t>(this->ndims, 0);

        size_t d = this->numel;
        size_t r = 0;
        for(size_t i = 0; i <= ndims - 2; i ++) {
            d /= shape[i];
            indices[i] = idx / d;

            r = idx % d;
            idx = r;
        }
        indices[this->ndims - 1] = r;

        return vec_to_index(indices);
    }

    size_t flatten(size_t i) {
        if (this->ndims != 1) {
            throw std::runtime_error("Cannot flatten indices of size " + std::to_string(this->ndims) + " into tensor of size 1");
        }
        return i;
    }

    size_t flatten(size_t i, size_t j) {
        if (this->ndims != 2) {
            throw std::runtime_error("Cannot flatten indices of size " + std::to_string(this->ndims) + " into tensor of size 2");
        }
        return i * this->shape[1] + j;
    }

    size_t flatten(size_t i, size_t j, size_t k) {
        if (this->ndims != 3) {
            throw std::runtime_error("Cannot flatten indices of size " + std::to_string(this->ndims) + " into tensor of size 3");
        }
        return i * this->shape[1] * this->shape[2] + j * this->shape[2] + k;
    }

    size_t flatten(size_t i, size_t j, size_t k, size_t l) {
        if (this->ndims != 4) {
            throw std::runtime_error("Cannot flatten indices of size " + std::to_string(this->ndims) + " into tensor of size 4");
        }
        return i * this->shape[1] * this->shape[2] * this->shape[3] + j * this->shape[2] * this->shape[3] + k * this->shape[3] + l;
    }

    void _set_shape(const std::vector<size_t>& shape) {
        this->shape = shape;
        this->numel = vecprod(this->shape);
        this->ndims = this->shape.size();
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