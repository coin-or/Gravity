#pragma once

#include <sstream>
#include <vector>
#include <numeric>
#include <queue>
#include <thread>

// We are using float limits as we cannot perform
// arithmetic ops on double limits (see extended_add I believe)
constexpr double HMAX = std::numeric_limits<float>::max(); // Max value for hidden states
constexpr double HMIN = std::numeric_limits<float>::lowest(); // Min value for hidden states

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
std::vector<T> concat(std::vector<T> a, std::vector<T> b) {
    std::vector<T> out;
    out.reserve(a.size() + b.size());
    out.insert(out.end(), a.begin(), a.end());
    out.insert(out.end(), b.begin(), b.end());
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
    std::vector<T> perm;
    for (auto i : indices) {
        perm.push_back(v.at(i));
    }

    return perm;
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

std::string ftostr(double v) {
    if (v == HMIN) {
        return "-∞";
    }
    if (v == HMAX) {
        return "∞";
    }
    return std::to_string(v);
}
