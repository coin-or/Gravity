#pragma once

#include <unordered_map>
#include <vector>
#include <gravity/types.h>
#include <gravity/var.h>
#include <gravity/param.h>
#include <network/types.hpp>

using namespace gravity;

class IndexSet {
public:
    IndexSet() = default;

    IndexSet(std::vector<std::vector<std::string>> names, indices& hidden_states, indices& w_ids, indices& y_ids) {
        this->_indices["Constr"] = indices("Constr");
        this->_indices["ConstrB"] = indices("ConstrB");

        for (const auto& name : names.at(0)) {
            this->_indices[name] = indices(name);
            this->_indices[name] = hidden_states;
        }

        for (const auto& name : names.at(1)) {
            this->_indices[name] = indices(name);
            this->_indices[name] = w_ids;
        }

        for (const auto& name : names.at(2)) {
            this->_indices[name] = indices(name);
            this->_indices[name] = y_ids;
        }
    }

    indices& operator[](const std::string& name) {
        return this->_indices.at(name);
    }

    std::unordered_map<std::string, indices> _indices;
    size_t row_id = 0;
    size_t row_id2 = 0;
};

class IndexContainer {
public:
    IndexContainer() {
        // Global indices
        this->hidden_states = indices("hidden_states");
        this->y_ids = indices("y_ids");

        // Params
        this->w = param<>("weights");
        this->w_ids = indices("w_ids");
        this->w.in(this->w_ids);
    }
    
    void add(OType op, std::vector<std::vector<std::string>> names) {
        if (this->_indices.find(op) != this->_indices.end()) {
            return;
        }

        names.resize(3);
        this->_indices[op] = IndexSet(
            names,
            this->hidden_states,
            this->w_ids,
            this->y_ids
        );
    }

    IndexSet& operator()(OType op) {
        return this->_indices.at(op);
    }

    // Weights and indices for weights
    param<> w;
    indices w_ids;

    // Indices for hidden states and binaries (y)
    indices hidden_states, y_ids;

    // Map from operator type to index set
    std::unordered_map<OType, IndexSet> _indices;
};