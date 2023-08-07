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

void final_run(std::string fname, const std::vector<Bound>& global_bounds, size_t obj_idx) {
    NeuralNet nn(fname);
    nn.set_aux_bounds(global_bounds);

    Model<>& NN = nn.build_model(obj_idx, "", "");

    solver<> S(NN,gurobi);
    auto grb_prog = (GurobiProgram*)(S._prog.get());
    auto grb_mod = grb_prog->grb_mod;
    grb_mod->set(GRB_IntParam_Threads, thread::hardware_concurrency() / 2);
    grb_mod->set(GRB_IntParam_OutputFlag, 1);
    grb_mod->set(GRB_IntParam_NonConvex, 2);

    int retval = S.run();
}

double bound_neuron(std::string fname, std::string start_node, Bound neuron, const std::vector<Bound>& global_bounds) {
    NeuralNet nn(fname);
    nn.set_aux_bounds(global_bounds);

    // Passing -1 means we will write a custom objective rather
    // than use one in the model
    Model<>& NN = nn.build_model(-1, start_node, neuron.layer_name);

    double mult = (neuron.side == LOWER) ? -1.0 : 1.0;
    NN.max(mult*nn.x(neuron.neuron_name));

    solver<> S(NN,gurobi);
    auto grb_prog = (GurobiProgram*)(S._prog.get());
    auto grb_mod = grb_prog->grb_mod;
    grb_mod->set(GRB_IntParam_Threads, 1);
    grb_mod->set(GRB_IntParam_NonConvex,2);
    grb_mod->set(GRB_IntParam_OutputFlag, 1);
    grb_mod->set(GRB_IntParam_MIPFocus, 3);
    // grb_mod->set(GRB_DoubleParam_BestBdStop, -1e-4);
    // grb_mod->set(GRB_DoubleParam_BestObjStop, 1e-4);

    // int retval = S.run(1e-4, 5.0);
    int retval = S.run(1e-6, 120.0);
    if (retval == -1) {
        // throw std::runtime_error("Infeasible");
        return mult * HMAX;
    }

    // Return relative objective value because we're dumping once we hit negative
    // Otherwise objval will be the current initialization if we didn't
    // find a feasible solution
    return mult * NN._rel_obj_val;
}

void run_MPI_bound_neuron(std::vector<Bound>& local_bounds, const std::string& fname, const std::string& start_node, const std::vector<Bound>& global_bounds) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

//    std::cout << "Rank " << rank << " of " << size << " started." << std::endl;

    // Determine the range of neurons this rank should handle.
    size_t neurons_per_rank = local_bounds.size() / size;
    size_t start_idx = rank * neurons_per_rank;
    size_t end_idx = (rank == size - 1) ? local_bounds.size() : start_idx + neurons_per_rank;

//    std::cout << "Rank " << rank << " is handling neurons from " << start_idx << " to " << end_idx << std::endl;

    // Each rank computes bounds for its subset of neurons.
    for (size_t i = start_idx; i < end_idx; i++) {
        Bound& neuron = local_bounds[i];
        auto new_bound = bound_neuron(fname, start_node, neuron, global_bounds);
        if (neuron.side == LOWER) {
            neuron.value = std::max(neuron.value, new_bound);
        } else {
            neuron.value = std::min(neuron.value, new_bound);
        }

//        std::cout << "Rank " << rank << " updated bound for neuron " << i << " with value " << neuron.value << std::endl;
    }

    // Gather the results at rank 0.
    if (rank == 0) {
        for (int src_rank = 1; src_rank < size; src_rank++) {
            size_t src_start_idx = src_rank * neurons_per_rank;
            size_t src_end_idx = (src_rank == size - 1) ? local_bounds.size() : src_start_idx + neurons_per_rank;
            MPI_Recv(&local_bounds[src_start_idx], (src_end_idx - src_start_idx) * sizeof(Bound), MPI_BYTE, src_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//            std::cout << "Rank 0 received data from rank " << src_rank << std::endl;
        }
    } else {
        MPI_Send(&local_bounds[start_idx], (end_idx - start_idx) * sizeof(Bound), MPI_BYTE, 0, 0, MPI_COMM_WORLD);
//        std::cout << "Rank " << rank << " sent data to rank 0." << std::endl;
    }

    // If necessary, rank 0 can then broadcast the complete results back to all other ranks.
    MPI_Bcast(&local_bounds[0], local_bounds.size() * sizeof(Bound), MPI_BYTE, 0, MPI_COMM_WORLD);
//    std::cout << "Rank " << rank << " completed broadcasting." << std::endl;
}

int main(int argc, char * argv[]) {
    string fname = string(prj_dir)+"/data_sets/VNN/tll_new_old.onnx";
    if(argc >= 2) {
        fname = argv[1];
    }

    NeuralNet nn(fname);
    std::vector<Layer*> layers_to_optimize;

    for (auto i = 1; i < nn._all_layers.size() - 1; i++) {
        if (
            (nn._all_layers[i+1]->operator_type != _relu) &&
            (nn._all_layers[i+1]->operator_type != _clip)
        ) {
            continue;
        }

        layers_to_optimize.push_back(nn._all_layers[i]);
    }
    
    // include last layer
    layers_to_optimize.push_back(nn._all_layers[nn._all_layers.size()-1]);

    std::cout << "Optimizing layers:" << std::endl;
    for (auto l: layers_to_optimize) {
        std::cout << l->lname() << std::endl;
    }

    std::vector<Bound> global_bounds;

    int rolling_horizon = 2;
    if (rolling_horizon > layers_to_optimize.size()){
        rolling_horizon = layers_to_optimize.size();
    }
    
    auto start_time_all = std::chrono::high_resolution_clock::now();
    
    MPI_Init(nullptr, nullptr);

    for (auto lidx = 0; lidx < layers_to_optimize.size(); lidx++) {
        auto l = layers_to_optimize[lidx];
        std::cout << "################################################" << std::endl;
        std::cout << "Optimizing layer: " << l->lname() << std::endl;
        std::cout << "Layer " << lidx+1 << "/" << layers_to_optimize.size() << std::endl;
        std::vector<Bound> local_bounds;
        for (auto o: l->outputs) {
            for (auto i = 0; i < o->numel; i++) {
                double lb  = o->lb.at(i);
                double ub  = o->ub.at(i);
                auto name = o->strkey(i);

                // If both LB and UB are on the same side of 0, we can skip this neuron
                if ((lb < 0 && ub < 0) || (lb > 0 && ub > 0)) {
                    continue;
                }

                local_bounds.push_back(Bound(l->lname(), name, lb, LOWER));
                local_bounds.push_back(Bound(l->lname(), name, ub, UPPER));
            }
        }

        std::cout << "Number of neurons to optimize: " << local_bounds.size()/2 << std::endl;
        
        // skip the rest of this loop iteration if local_bounds is of size 0
        if (local_bounds.size() == 0) {
            continue;
        }
        
        int bak, new_;
        fflush(stdout);
        bak = dup(1);
        new_ = open("/dev/null", O_WRONLY);
        dup2(new_, 1);
        close(new_);

        auto start_time = std::chrono::high_resolution_clock::now();

        std::string start_node = "";
        if (lidx > rolling_horizon - 1) {
            start_node = layers_to_optimize[lidx - (rolling_horizon - 1)]->name;
        }

        auto num_threads = std::thread::hardware_concurrency();
        
        // Use the run_MPI_bound_neuron function for parallel neuron bound calculation
        run_MPI_bound_neuron(local_bounds, fname, start_node, global_bounds);
        auto end_time = std::chrono::high_resolution_clock::now();
        fflush(stdout);
        dup2(bak, 1);
        close(bak);

        // print out the final bounds
        for (int i = 0; i < local_bounds.size()-1; i+=2) {
            auto lb = local_bounds[i];
            auto ub = local_bounds[i+1];
            std::cout << lb.neuron_name << ": ";
            std::cout << "[" << ftostr(lb.old_value) << ", " << ftostr(ub.old_value) << "] -> ";
            std::cout << "[" << ftostr(lb.value) << ", " << ftostr(ub.value) << "]";
            std::cout << std::endl;
        }
        auto time_taken = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
        std::cout << "Time: " << time_taken << "ms" << std::endl;
        global_bounds.insert(global_bounds.end(), local_bounds.begin(), local_bounds.end());
    }
    MPI_Finalize();

    std::cout << "Starting final runs" << std::endl;

    for (size_t obj_idx = 0; obj_idx < nn.obj_spec->shape[0]; obj_idx++) {
        std::cout << "########################################" << std::endl;
        final_run(fname, global_bounds, obj_idx);
    }
    
    auto end_time_all = std::chrono::high_resolution_clock::now();
    auto time_taken_all = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_all - start_time_all).count();
    std::cout << "Time in total: " << time_taken_all << "ms" << std::endl;

    return 0;
}
