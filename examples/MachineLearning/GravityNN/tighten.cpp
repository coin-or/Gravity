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
    NeuralNet nn(fname);
    nn.set_aux_bounds(global_bounds);
    Model<>& NN = nn.build_model();

    NN.max(nn.x(nn.layers.back()->outputs[0]->strkey(0)));

    solver<> S(NN,gurobi);
    auto grb_prog = (GurobiProgram*)(S._prog.get());
    auto grb_mod = grb_prog->grb_mod;
    grb_mod->set(GRB_IntParam_Threads, get_num_threads() / 2);
    grb_mod->set(GRB_IntParam_OutputFlag, 1);

    int retval = S.run();
}

float bound_neuron(std::string fname, Bound neuron, const std::vector<Bound>& global_bounds) {
    NeuralNet nn(fname, neuron.layer_name);
    nn.set_aux_bounds(global_bounds);

    Model<>& NN = nn.build_model();

    float mult = (neuron.side == LOWER) ? -1.0 : 1.0;
    NN.max(nn.x(neuron.neuron_name) * mult);

    solver<> S(NN,gurobi);
    auto grb_prog = (GurobiProgram*)(S._prog.get());
    auto grb_mod = grb_prog->grb_mod;
    grb_mod->set(GRB_IntParam_Threads, 1);
    // grb_mod->set(GRB_IntParam_NonConvex,2);
    grb_mod->set(GRB_IntParam_OutputFlag, 1);
    grb_mod->set(GRB_IntParam_MIPFocus, 3);
    // grb_mod->set(GRB_DoubleParam_BestBdStop, -1e-4);
    // grb_mod->set(GRB_DoubleParam_BestObjStop, 1e-4);

    int retval = S.run(1e-4, 2.0);

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

    for (auto i = 0; i < nn.layers.size()-1; i++) {
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

        std::cout << "Number of neurons to optimize: " << local_bounds.size() << std::endl;
        int bak, new_;
        fflush(stdout);
        bak = dup(1);
        new_ = open("/dev/null", O_WRONLY);
        dup2(new_, 1);
        close(new_);

        auto start_time = std::chrono::high_resolution_clock::now();
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
        auto end_time = std::chrono::high_resolution_clock::now();
        fflush(stdout);
        dup2(bak, 1);
        close(bak);

        // print out the final bounds
        for (int i = 0; i < local_bounds.size()-1; i+=2) {
            auto lb = local_bounds[i];
            auto ub = local_bounds[i+1];
            std::cout << lb.neuron_name << ": ";
            std::cout << "[" << lb.old_value << ", " << ub.old_value << "] -> ";
            std::cout << "[" << lb.value << ", " << ub.value << "]" << std::endl;
        }
        auto time_taken = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
        std::cout << "Time: " << time_taken << "ms" << std::endl;
        global_bounds.insert(global_bounds.end(), local_bounds.begin(), local_bounds.end());
    }

    std::cout << "Starting final run" << std::endl;
    final_run(fname, global_bounds);

    return 0;
}