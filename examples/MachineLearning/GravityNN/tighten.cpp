#include <iostream>
#include <fstream>
#include <onnx.pb.h>
#include <vector>
#include <gravity/solver.h>
#include <network/NeuralNet.hpp>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

using namespace gravity;

void final_run(std::string fname, const std::vector<Bound>& global_bounds) {
    std::cout << "A" << std::endl;
    NeuralNet nn(fname);
    std::cout << "B" << std::endl;

    std::cout << "C" << std::endl;
    nn.build_indexing();
    std::cout << "D" << std::endl;
    nn.build_constraints();
    std::cout << "E" << std::endl;

    Model<> NN("NN_"+fname.substr(fname.find_last_of("/")));
    std::cout << "F" << std::endl;
    param<> x_lb("x_lb"), x_ub("x_ub");
    std::cout << "G" << std::endl;
    x_lb.in(nn.indices.hidden_states);
    std::cout << "H" << std::endl;
    x_ub.in(nn.indices.hidden_states);
    std::cout << "I" << std::endl;
    x_lb = std::numeric_limits<double>::lowest();
    std::cout << "J" << std::endl;
    x_ub = std::numeric_limits<double>::max();
    std::cout << "K" << std::endl;

    nn.set_bounds(x_lb, x_ub);
    std::cout << "L" << std::endl;
    // use newbounds
    for (auto& v: global_bounds) {
        if (nn.layer_names.count(v.layer_name) == 0) {
            continue;
        }

        if (v.side == LOWER) {
            x_lb.set_val(v.neuron_name, v.value);
        } else {
            x_ub.set_val(v.neuron_name, v.value);
        }
    }
    std::cout << "M" << std::endl;

    var<> x("x", x_lb, x_ub);
    std::cout << "N" << std::endl;
    var<int> y("y", 0, 1);
    std::cout << "O" << std::endl;
    nn.initialize_state(x, y);
    std::cout << "P" << std::endl;

    NN.add(x.in(nn.indices.hidden_states));
    std::cout << "Q" << std::endl;
    NN.add(y.in(nn.indices.y_ids));
    std::cout << "R" << std::endl;
    nn.add_constraints(NN, x, y, nn.indices);
    std::cout << "S" << std::endl;

    NN.max(x(nn.layers.back()->outputs[0]->strkey(0)));
    std::cout << "T" << std::endl;

    solver<> S(NN,gurobi);
    auto grb_prog = (GurobiProgram*)(S._prog.get());
    auto grb_mod = grb_prog->grb_mod;
    grb_mod->set(GRB_IntParam_Threads, 64);
    grb_mod->set(GRB_IntParam_OutputFlag, 1);
    std::cout << "U" << std::endl;

    int retval = S.run();
    std::cout << "V" << std::endl;

    if (retval == 3) {
        throw std::runtime_error("Infeasible");
    }
}

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
    // use newbounds
    for (auto& v: global_bounds) {
        if (nn.layer_names.count(v.layer_name) == 0) {
            continue;
        }

        if (v.side == LOWER) {
            x_lb.set_val(v.neuron_name, v.value);
        } else {
            x_ub.set_val(v.neuron_name, v.value);
        }
    }

    var<> x("x", x_lb, x_ub);
    var<int> y("y", 0, 1);

    NN.add(x.in(nn.indices.hidden_states));
    NN.add(y.in(nn.indices.y_ids));
    nn.initialize_state(x, y);
    nn.add_constraints(NN, x, y, nn.indices);

    float mult = (neuron.side == LOWER) ? -1.0 : 1.0;
    NN.max(x(neuron.neuron_name) * mult);

    solver<> S(NN,gurobi);
    auto grb_prog = (GurobiProgram*)(S._prog.get());
    auto grb_mod = grb_prog->grb_mod;
    grb_mod->set(GRB_IntParam_Threads, 1);
    // grb_mod->set(GRB_IntParam_NonConvex,2);
    grb_mod->set(GRB_IntParam_OutputFlag, 0);
    grb_mod->set(GRB_IntParam_MIPFocus, 3);
    // grb_mod->set(GRB_DoubleParam_BestBdStop, -1e-4);
    // grb_mod->set(GRB_DoubleParam_BestObjStop, 1e-4);

    int retval = S.run(1e-4, 5.0);

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

        // if (!found_relu) {
            // continue;
        // }

        if (
            (nn.layers[i+1]->operator_type != _relu) &&
            (nn.layers[i+1]->operator_type != _clip)
        ) {
            continue;
        }

        layers_to_optimize.push_back(nn.layers[i]);
    }

    std::vector<Bound> global_bounds;
    for (auto l: layers_to_optimize) {
        std::cout << "################################################" << std::endl;
        std::cout << "Optimizing layer: " << l->name << std::endl;
        std::vector<Bound> local_bounds;
        for (auto i = 0; i < l->outputs[0]->numel; i++) {
            float lb = std::numeric_limits<float>::lowest();
            float ub = std::numeric_limits<float>::max();
            if (l->lowers.size() > 0) {
                lb = (*l->lowers[0])(i);
                ub = (*l->uppers[0])(i);
            }

            auto name = l->outputs[0]->strkey(i);

            // If both LB and UB are on the same side of 0, we can skip this neuron
            if ((lb < 0 && ub < 0) || (lb > 0 && ub > 0)) {
                continue;
            }

            local_bounds.push_back(Bound(l->name, name, lb, LOWER));
            local_bounds.push_back(Bound(l->name, name, ub, UPPER));
        }

        std::cout << "Number of neurons to optimize: " << local_bounds.size() << std::endl;
        int bak, new_;
        fflush(stdout);
        bak = dup(1);
        new_ = open("/dev/null", O_WRONLY);
        dup2(new_, 1);
        close(new_);

        #pragma omp parallel for num_threads(64)
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

        // print out the final bounds
        for (auto& neuron: local_bounds) {
            if (neuron.side == LOWER) {
                std::cout << neuron.neuron_name << ", " << neuron.value << ", LOWER" << std::endl;
            } else {
                std::cout << neuron.neuron_name << ", " << neuron.value << ", UPPER" << std::endl;
            }
        }
        global_bounds.insert(global_bounds.end(), local_bounds.begin(), local_bounds.end());
    }

    std::cout << "Starting final run" << std::endl;
    final_run(fname, global_bounds);

    return 0;
}