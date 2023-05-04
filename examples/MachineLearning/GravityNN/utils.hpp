#pragma once

#include <sstream>
#include <vector>

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