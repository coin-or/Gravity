#pragma once

#include <sstream>
#include <vector>
#include <numeric>

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

    std::fstream input(onnx_path, std::ios::in | std::ios::binary);
    onnx::ModelProto model;
    bool isSuccess = model.ParseFromIstream(&input);
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