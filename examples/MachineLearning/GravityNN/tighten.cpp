#include <iostream>
#include <fstream>
#include <onnx.pb.h>
#include <vector>
#include <gravity/solver.h>
#include <network/NeuralNet.hpp>

using namespace gravity;

float bound_neuron(std::string fname, Bound neuron, const std::vector<Bound>& global_bounds) {
    NeuralNet nn(fname, neuron.layer_name);

    nn.build_indexing();
    nn.build_constraints();

    Model<> NN("NN_"+fname.substr(fname.find_last_of("/")));
    param<> x_lb("x_lb"), x_ub("x_ub");
    x_lb.in(nn.indices.hidden_states);
    x_ub.in(nn.indices.hidden_states);
    x_lb = std::numeric_limits<double>::lowest();
    x_ub = std::numeric_limits<double>::max();

    nn.set_bounds(x_lb, x_ub);

    var<> x("x", x_lb, x_ub);
    var<int> y("y", 0, 1);

    // use newbounds
    for (auto v: global_bounds) {
        if (nn.layer_names.count(v.layer_name) == 0) {
            continue;
        }

        if (v.side == LOWER) {
            x_lb.set_val(v.neuron_name, v.value);
        } else {
            x_ub.set_val(v.neuron_name, v.value);
        }
    }

    NN.add(x.in(nn.indices.hidden_states));
    NN.add(y.in(nn.indices.y_ids));
    nn.initialize_state(x, y);
    nn.add_constraints(NN, x, y, nn.indices);

    float mult = (neuron.side == LOWER) ? -1.0 : 1.0;
    NN.max(x(neuron.neuron_name) * mult);

    solver<> S(NN,gurobi);
    auto grb_prog = (GurobiProgram*)(S._prog.get());
    auto grb_mod = grb_prog->grb_mod;
    grb_mod->set(GRB_IntParam_Threads, 192);
    // grb_mod->set(GRB_IntParam_NonConvex,2);
    grb_mod->set(GRB_IntParam_OutputFlag, 0);
    grb_mod->set(GRB_IntParam_MIPFocus,3);
    grb_mod->set(GRB_DoubleParam_BestBdStop, -1e-4);
    grb_mod->set(GRB_DoubleParam_BestObjStop, 1e-4);

    int retval = S.run(1e-4, 10.0);

    if (retval == 3) {
        throw std::runtime_error("Infeasible");
    }

    // Return relative objective value because we're dumping once we hit negative
    // Otherwise objval will be the current initialization if we didn't
    // find a feasible solution
    return mult * NN._rel_obj_val;
}

int main(int argc, char * argv[]) {
    string fname = string(prj_dir)+"/data_sets/VNN/tll_bound.onnx";
    if(argc >= 2) {
        fname = argv[1];
    }

    NeuralNet nn(fname);
    std::vector<Layer*> layers_to_optimize;

    bool found_relu = false;
    // Only optimize a layer once we have found a ReLU
    for (auto i = 0; i < nn.layers.size()-1; i++) {
        if (nn.layers[i]->operator_type == _relu) {
            found_relu = true;
        }

        if (!found_relu) {
            continue;
        }

        if (nn.layers[i+1]->operator_type != _relu) {
            continue;
        }

        layers_to_optimize.push_back(nn.layers[i]);
    }

    std::cout << "Layers to optimize: " << std::endl;
    for (auto l: layers_to_optimize) {
        std::cout << l->name << std::endl;
    }

    std::vector<Bound> global_bounds;
    for (auto l: layers_to_optimize) {
        for (auto i = 0; i < l->outputs[0]->numel; i++) {
            auto lb = (*l->lowers[0])(i);
            auto ub = (*l->uppers[0])(i);
            auto name = l->outputs[0]->strkey(i);

            // If both LB and UB are on the same side of 0, we can skip this neuron
            if ((lb < 0 && ub < 0) || (lb > 0 && ub > 0)) {
                continue;
            }

            // Optimize the smaller bound
            if (std::abs(lb) < std::abs(ub)) {
                global_bounds.push_back(Bound(l->name, name, lb, LOWER));
            } else {
                global_bounds.push_back(Bound(l->name, name, ub, UPPER));
            }
        }

    }

    // sort neurons by abs(value)
    std::sort(global_bounds.begin(), global_bounds.end(), [](const Bound& a, const Bound& b) {
        return std::abs(a.value) < std::abs(b.value);
    });

    size_t fixed_count = 0;
    size_t fail_count = 0;
    for (auto& neuron: global_bounds) {
        auto new_bound = bound_neuron(fname, neuron, global_bounds);
        auto prev_bound = neuron.value;
        if (neuron.side == LOWER) {
            neuron.value = std::max(neuron.value, new_bound);
        } else {
            neuron.value = std::min(neuron.value, new_bound);
        }

        // if the sign of the bound changed, we can fix this neuron
        if (std::signbit(prev_bound) != std::signbit(neuron.value)) {
            fixed_count++;
        } else {
            fail_count++;
        }
        std::cout << "#################################################" << std::endl;
        std::cout << "Neuron: " << neuron.neuron_name << std::endl;
        std::cout << prev_bound << " -> " << neuron.value << std::endl;
        std::cout << "Fixed: " << fixed_count << "/" << global_bounds.size() << std::endl;
        std::cout << "Failed: " << fail_count << "/" << global_bounds.size() << std::endl;
        if (neuron.side == LOWER) {
            std::cout << "Lower bound" << std::endl;
        } else {
            std::cout << "Upper bound" << std::endl;
        }
        std::cout << "#################################################" << std::endl;
    }

    // print out the final bounds
    std::cout << "Final bounds: " << std::endl;
    for (auto& neuron: global_bounds) {
        if (neuron.side == LOWER) {
            std::cout << neuron.neuron_name << ", " << neuron.value << ", LOWER" << std::endl;
        } else {
            std::cout << neuron.neuron_name << ", " << neuron.value << ", UPPER" << std::endl;
        }
    }

    return 0;
}