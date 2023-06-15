#pragma once

#include <vector>
#include <gravity/solver.h>
#include <network/NeuralNet.hpp>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

std::vector<Bound> tighten_layer(NeuralNet& nn, std::string layer_name, std::string fname) {
    Layer* l = nullptr;
    for (auto layer: nn.layers) {
        if (layer->name == layer_name) {
            l = layer;
            break;
        }
    }
    if (l == nullptr) {
        throw std::runtime_error("Layer not found");
    }

    std::vector<Bound> local_bounds;
    for (auto o: l->outputs) {
        for (auto i = 0; i < o->numel; i++) {
            float lb  = o->lb.at(i);
            float ub  = o->ub.at(i);
            auto name = o->strkey(i);

            // If both LB and UB are on the same side of 0, we can skip this neuron
            if ((lb < 0 && ub < 0) || (lb > 0 && ub > 0)) {
                continue;
            }

            local_bounds.push_back(Bound(l->name, name, lb, LOWER));
            local_bounds.push_back(Bound(l->name, name, ub, UPPER));
        }
    }

    int bak, new_;
    fflush(stdout);
    bak = dup(1);
    new_ = open("/dev/null", O_WRONLY);
    dup2(new_, 1);
    close(new_);

    #pragma omp parallel for
    for (auto& neuron: local_bounds) {
        auto new_bound = bound_neuron(fname, neuron, global_bounds);
        auto prev_bound = neuron.value;
        if (neuron.side == LOWER) {
            neuron.value = std::max(neuron.value, new_bound);
        } else {
            neuron.value = std::min(neuron.value, new_bound);
        }
    }
    fflush(stdout);
    dup2(bak, 1);
    close(bak);

    return local_bounds;
}