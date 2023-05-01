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