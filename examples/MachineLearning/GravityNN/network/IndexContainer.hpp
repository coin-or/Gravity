#pragma once

#include <unordered_map>
#include <vector>
#include <gravity/types.h>
#include <gravity/var.h>
#include <gravity/param.h>
#include <network/types.hpp>

using namespace gravity;

typedef std::unordered_map<std::string, std::string> NoOpRemap;

class SetWrapper {
public:
    SetWrapper() = default;

    SetWrapper(const std::string& name, std::shared_ptr<NoOpRemap> noop_remap): noop_remap(noop_remap) {
        this->_indices = indices(name);
    }

    SetWrapper& operator=(const indices& other) {
        this->_indices = other;
        return *this;
    }

    void add(const std::string& name) {
        this->_indices.add(this->_remap_key(name));
    }

    void add_ref(const std::string& name) {
        this->_indices.add_ref(this->_remap_key(name));
    }

    void add_in_row(size_t row_nb, const std::string& name) {
        this->_indices.add_in_row(row_nb, this->_remap_key(name));
    }

    void add_empty_row() {
        this->_indices.add_empty_row();
    }

    // Implicit conversion to gravity::indices
    operator indices&() {
        return this->_indices;
    }

    std::string _remap_key(const std::string& key) {
        if (this->noop_remap->find(key) != this->noop_remap->end()) {
            return this->noop_remap->at(key);
        }

        return key;
    }

    indices _indices;
    std::shared_ptr<NoOpRemap> noop_remap;
};

class IndexSet {
public:
    IndexSet() = default;

    IndexSet(
        std::vector<std::vector<std::string>> names,
        indices& hidden_states,
        indices& w_ids, indices& y_ids,
        std::shared_ptr<NoOpRemap> noop_remap
    ) {
        this->noop_remap = noop_remap;

        this->_indices["Constr"] = SetWrapper("Constr", noop_remap);
        this->_indices["ConstrB"] = SetWrapper("ConstrB", noop_remap);

        for (const auto& name : names.at(0)) {
            this->_indices[name] = SetWrapper(name, noop_remap);
            this->_indices[name] = hidden_states.subset();
        }

        for (const auto& name : names.at(1)) {
            this->_indices[name] = SetWrapper(name, noop_remap);
            this->_indices[name] = w_ids.subset();
        }

        for (const auto& name : names.at(2)) {
            this->_indices[name] = SetWrapper(name, noop_remap);
            this->_indices[name] = y_ids.subset();
        }
    }

    SetWrapper& operator[](const std::string& name) {
        return this->_indices.at(name);
    }

    std::unordered_map<std::string, SetWrapper> _indices;
    std::shared_ptr<NoOpRemap> noop_remap;

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

        this->noop_remap = std::make_shared<NoOpRemap>();
    }
    
    void add(std::string layer_name, std::vector<std::vector<std::string>> names) {
        if (this->_indices.find(layer_name) != this->_indices.end()) {
            return;
        }

        names.resize(3);
        this->_indices[layer_name] = IndexSet(
            names,
            this->hidden_states,
            this->w_ids,
            this->y_ids,
            this->noop_remap
        );
    }

    IndexSet& operator()(std::string layer_name) {
        return this->_indices.at(layer_name);
    }

    void add_remap(std::string true_key, std::string aux_key) {
        // Check if true_key is already in the map
        // If it is, we want aux_key to point to the value
        // that true_key points to
        if (this->noop_remap->find(true_key) != this->noop_remap->end()) {
            true_key = this->noop_remap->at(true_key);
        }
        this->noop_remap->insert({aux_key, true_key});
    }

    // Weights and indices for weights
    param<> w;
    indices w_ids;

    // Indices for hidden states and binaries (y)
    indices hidden_states, y_ids;

    // Map from operator type to index set
    std::unordered_map<std::string, IndexSet> _indices;

    /*
        Some operators like Reshape do nothing. When we reference a variable from the output of reshape
        we actually want to reference the input to reshape
    */
    std::shared_ptr<NoOpRemap> noop_remap;
};
