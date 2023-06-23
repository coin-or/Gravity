#pragma once

#include <fstream>
#include <onnx.pb.h>
#include <vector>
#include <utils/utils.hpp>
#include <gravity/param.h>
#include <cstring>

template <typename T>
std::vector<T> read_data(const onnx::TensorProto& tensor, const ::google::protobuf::RepeatedField<T>& field);

class Tensor {
public:
    Tensor() {}

    Tensor(const onnx::TensorProto& tensor) {
        this->name = tensor.name();
        this->is_initializer = true;

        this->_set_shape(std::vector<size_t>(tensor.dims().begin(), tensor.dims().end()));

        switch (tensor.data_type()) {
            case onnx::TensorProto::FLOAT: {
                auto data = read_data<float>(tensor, tensor.float_data());
                this->data.clear();
                for (auto d : data) {
                    this->data.push_back(d);
                }
                break;
            }
            case onnx::TensorProto::INT64:
                this->int_data = read_data<int64_t>(tensor, tensor.int64_data());
                break;
            default:
                throw std::runtime_error("Tensor " + this->name + " has unsupported type " + std::to_string(tensor.data_type()));
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

    static Tensor transpose(const Tensor& t) {
        if (t.shape.size() != 2) {
            throw std::runtime_error("Cannot transpose tensor with shape " + std::to_string(t.shape.size()));
        }

        Tensor newt = t;

        if (t.is_initializer) {
            // Tranpose data if needed, it is stored row-major
            for (size_t i = 0; i < t.shape[0]; i++) {
                for (size_t j = 0; j < t.shape[1]; j++) {
                    newt.data[j * t.shape[0] + i] = t.data[i * t.shape[1] + j];
                }
            }
        }

        // Swap shape
        std::swap(newt.shape[0], newt.shape[1]);

        return newt;
    }

    double operator()(size_t i) const {
        if (!this->is_initializer) {
            throw std::runtime_error("Reading from non-initializer tensor. Perhaps you're assuming this tensor is a weight when it's actually an output of a previous layer?");
        }
        return this->data.at(i);
    }

    std::vector<int64_t> get_int_data() {
        return this->int_data;
    }

    void add_params(gravity::param<>& p) {
        if (!this->is_initializer) {
            throw std::runtime_error("Reading from non-initializer tensor. Perhaps you're assuming this tensor is a weight when it's actually an output of a previous layer?");
        }

        if (this->added_to_vals) {
            return;
        }
        this->added_to_vals = true;

        for (size_t i = 0; i < this->numel; i++) {
            p.add_val(this->strkey(i), this->data.at(i));
        }
    }

    std::string strkey(std::vector<size_t> indices) const {
        size_t flat = this->flatten_index(indices);
        return Tensor::_strkey(this->name, flat, this->numel);
    }

    std::string strkey(size_t idx) const {
        size_t flat = idx;
        return Tensor::_strkey(this->name, flat, this->numel);
    }

    std::string strkey(size_t i, size_t j) const {
        if (this->ndims != 2) {
            throw std::runtime_error("Cannot use 2D indexing on tensor with " + std::to_string(this->ndims) + " dimensions");
        }
        size_t flat = i * this->shape[1] + j;
        return Tensor::_strkey(this->name, flat, this->numel);
    }

    std::string strkey(size_t i, size_t j, size_t k) const {
        if (this->ndims != 3) {
            throw std::runtime_error("Cannot use 3D indexing on tensor with " + std::to_string(this->ndims) + " dimensions");
        }
        size_t flat = i * this->shape[1] * this->shape[2] + j * this->shape[2] + k;
        return Tensor::_strkey(this->name, flat, this->numel);
    }

    std::string strkey(size_t i, size_t j, size_t k, size_t l) const {
        if (this->ndims != 4) {
            throw std::runtime_error("Cannot use 4D indexing on tensor with " + std::to_string(this->ndims) + " dimensions");
        }
        size_t flat = i * this->shape[1] * this->shape[2] * this->shape[3] + j * this->shape[2] * this->shape[3] + k * this->shape[3] + l;
        return Tensor::_strkey(this->name, flat, this->numel);
    }

    size_t flatten_index(const std::vector<size_t>& indices) const {
        return Tensor::flatten_index(indices, this->shape);
    }

    static std::string _strkey(const std::string& name, size_t idx, size_t numel) {
        if (idx >= numel) {
            throw std::runtime_error("Index " + std::to_string(idx) + " out of bounds for tensor of size " + std::to_string(numel) + ".");
        }
        std::string str;
        str.reserve(name.size() + 1 + 5);
        str.append(name);
        str.push_back(',');
        str.append(std::to_string(idx));
        return str;
    }

    static size_t flatten_index(const std::vector<size_t>& indices, const std::vector<size_t>& shape) {
        if (indices.size() != shape.size()) {
            throw std::runtime_error("Cannot flatten index with " + std::to_string(indices.size()) + " indices for a tensor with " + std::to_string(shape.size()) + " dimensions");
        }

        size_t index = 0;
        size_t stride = 1;
        for (int64_t i = shape.size()-1; i>=0; --i) {
            index += indices[i] * stride;
            stride *= shape[i];
        }
        return index;
    }

    std::vector<size_t> unflatten_index(size_t index) {
        return Tensor::unflatten_index(index, this->shape);
    }

    static std::vector<size_t> unflatten_index(size_t index, const std::vector<size_t>& shape) {
        if (shape.empty()) {
            return {};
        }

        std::vector<size_t> result;
        std::vector<size_t> revshape = shape;
        std::reverse(revshape.begin(), revshape.end());
        for (const auto& size : revshape) {
            result.push_back(index % size);
            index = index / size;
        }
        std::reverse(result.begin(), result.end());
        return result;
    }

    void _set_shape(const std::vector<size_t>& shape) {
        this->shape = shape;
        if (this->shape.size() == 0) {
            this->shape.push_back(1);
        }
        this->numel = vecprod(this->shape);
        this->ndims = std::max(this->shape.size(), (size_t)1);

        this->lb = std::vector<double>(this->numel, HMIN);
        this->ub = std::vector<double>(this->numel, HMAX);
        this->forward = std::vector<double>(this->numel, 0);
    }

    std::vector<double> get_data() {
        return this->data;
    }

    std::string name;
    bool is_initializer;

    std::vector<size_t> shape;
    size_t numel;
    size_t ndims;

    // LB/UB/Forward
    std::vector<double> lb;
    std::vector<double> ub;
    std::vector<double> forward;

private:
    std::vector<double> data;
    std::vector<int64_t> int_data;
    bool added_to_vals = false;
};

typedef std::map<std::string, Tensor> Tensors;

Tensors get_tensors(const onnx::GraphProto& graph) {
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

    // Find constants
    for (const auto& node : graph.node()) {
        if (node.op_type() != "Constant") {
            continue;
        }

        std::string name = node.output(0);
        if (tensors.find(name) == tensors.end()) {
            throw std::runtime_error("Constant " + name + " not found in graph");
        }

        Tensor& tensor = tensors[name];
        if (tensor.is_initializer) {
            throw std::runtime_error("Constant " + name + " is already initialized");
        }
        Tensor datat = Tensor(node.attribute(0).t());
        datat.name = name;
        tensors[name] = datat;
    }

    // Search for tensors that end in _lower, _upper, or _forward
    std::set<std::string> todel;
    for (auto& tup: tensors) {
        auto name = tup.first;
        auto& tensor = tup.second;
        
        if (ends_with(name, "_lower")) {
            auto base_name = name.substr(0, name.size() - 6);
            if (tensors.find(base_name) == tensors.end()) {
                throw std::runtime_error("Tensor " + base_name + " not found in graph");
            }
            Tensor& base_tensor = tensors.at(base_name);
            base_tensor.lb = tensor.get_data();
            todel.insert(name);
        } else if (ends_with(name, "_upper")) {
            auto base_name = name.substr(0, name.size() - 6);
            if (tensors.find(base_name) == tensors.end()) {
                throw std::runtime_error("Tensor " + base_name + " not found in graph");
            }
            Tensor& base_tensor = tensors.at(base_name);
            base_tensor.ub = tensor.get_data();
            todel.insert(name);
        } else if (ends_with(name, "_forward")) {
            auto base_name = name.substr(0, name.size() - 8);
            if (tensors.find(base_name) == tensors.end()) {
                throw std::runtime_error("Tensor " + base_name + " not found in graph");
            }
            Tensor& base_tensor = tensors.at(base_name);
            base_tensor.forward = tensor.get_data();
            todel.insert(name);
        }
    }

    // Delete _lower, _upper, and _forward tensors
    for (const auto& name : todel) {
        tensors.erase(name);
    }

    return tensors;
}

template <typename T>
std::vector<T> read_data(const onnx::TensorProto& tensor, const ::google::protobuf::RepeatedField<T>& field) {
    if (field.size() > 0) {
        return {field.begin(), field.end()};
    } else if (!tensor.raw_data().empty()) {
        size_t elem_size = sizeof(T);
        size_t num_elems = tensor.raw_data().size() / elem_size;

        if (tensor.raw_data().size() % elem_size != 0) {
            throw std::runtime_error("Size of raw data is not a multiple of element size.");
        }

        std::vector<T> vec(num_elems);
        std::memcpy(vec.data(), tensor.raw_data().data(), tensor.raw_data().size());
        return vec;
    } else {
        throw std::runtime_error("Empty tensor.");
    }
}
