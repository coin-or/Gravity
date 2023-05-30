#pragma once

#include <map>
#include <vector>
#include <gravity/types.h>
#include <gravity/var.h>
#include <gravity/param.h>
#include <network/types.hpp>

using namespace gravity;

class IndexSet {
public:
    IndexSet() {}

    IndexSet(std::vector<string> state_names, std::vector<std::string> w_names, std::vector<std::string> y_names, indices& hidden_states, indices& w_ids, indices& y_ids) {
        this->_indices["Constr"] = indices("Constr");
        this->_indices["ConstrB"] = indices("ConstrB");

        for (const auto& name : state_names) {
            this->_indices[name] = indices(name);
            this->_indices[name] = hidden_states;
        }

        for (const auto& name : w_names) {
            this->_indices[name] = indices(name);
            this->_indices[name] = w_ids;
        }

        for (const auto& name : y_names) {
            this->_indices[name] = indices(name);
            this->_indices[name] = y_ids;
        }
    }

    indices& operator[](const std::string& name) {
        if (this->_indices.count(name) == 0) {
            throw std::runtime_error("Unknown index set: " + name);
        }
        return this->_indices.at(name);
    }

    std::map<std::string, indices> _indices;
    size_t row_id = 0;
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
        if (this->_indices.count(op) > 0) {
            return;
        }

        if (names.size() == 2) {
            this->_indices[op] = IndexSet(names[0], names[1], {}, this->hidden_states, this->w_ids, this->y_ids);
        } else {
            this->_indices[op] = IndexSet(names[0], names[1], names[2], this->hidden_states, this->w_ids, this->y_ids);
        }
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
    std::map<OType, IndexSet> _indices;
};