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
#include <cstdlib>

using namespace gravity;

void final_run(std::string fname, const std::vector<Bound>& global_bounds, size_t obj_idx) {
    // char hostname[128];
    // gethostname(hostname, sizeof(hostname));
    // std::cout << "Running final_run on host: " << hostname << " with process ID: " << getpid() << std::endl;

    NeuralNet nn(fname);
    nn.set_aux_bounds(global_bounds);

    Model<>& NN = nn.build_model(obj_idx, "", "");

    solver<> S(NN,gurobi);
    auto grb_prog = (GurobiProgram*)(S._prog.get());
    auto grb_mod = grb_prog->grb_mod;
    grb_mod->set(GRB_IntParam_Threads, thread::hardware_concurrency() / 2);
    grb_mod->set(GRB_IntParam_OutputFlag, 1);
    grb_mod->set(GRB_IntParam_NonConvex, 2);
    
//    grb_mod->set(GRB_DoubleParam_Heuristics, 0);
//    grb_mod->set(GRB_IntParam_CutPasses, 5);
//    grb_mod->set(GRB_IntParam_VarBranch, 0);
//
//    double tolerance = grb_mod->get(GRB_DoubleParam_MIPGap);
//
//    // Stop the program if the lower bound exceeds 0
//    grb_mod->set(GRB_DoubleParam_BestBdStop, tolerance);
//
//    // Stop the program if a feasible solution with objective less than 0 is found
//    grb_mod->set(GRB_DoubleParam_BestObjStop, -tolerance);

    int retval = S.run();
}

double bound_neuron(std::string fname, std::string start_node, Bound neuron, const std::vector<Bound>& global_bounds) {
    // char hostname[128];
    // gethostname(hostname, sizeof(hostname));
    // std::cout << "Running bound_neuron on host: " << hostname << " with process ID: " << getpid() << " bounding " << neuron.side << " of neuron " << neuron.neuron_name << std::endl;

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

constexpr int MAX_LAYER_NAME_LEN = 100;
constexpr int MAX_NEURON_NAME_LEN = 100;

MPI_Datatype MPI_PROXY_BOUND_TYPE;
const int nitems = 5;
int blocklengths[5] = {MAX_LAYER_NAME_LEN, MAX_LAYER_NAME_LEN, 1, 1, 1};
MPI_Datatype types[5] = {MPI_CHAR, MPI_CHAR, MPI_DOUBLE, MPI_DOUBLE, MPI_INT};
MPI_Aint offsets[5];

class ProxyBound {
public:
    char layer_name[MAX_LAYER_NAME_LEN];
    char neuron_name[MAX_NEURON_NAME_LEN];
    double value;
    double old_value;
    Side side;

    ProxyBound(const Bound& b) {
        strncpy(layer_name, b.layer_name.c_str(), MAX_LAYER_NAME_LEN);
        layer_name[MAX_LAYER_NAME_LEN - 1] = '\0';

        strncpy(neuron_name, b.neuron_name.c_str(), MAX_NEURON_NAME_LEN);
        neuron_name[MAX_NEURON_NAME_LEN - 1] = '\0';

        value = b.value;
        old_value = b.old_value;
        side = b.side;
    }

    Bound toBound() const {
        return Bound(std::string(layer_name), std::string(neuron_name), value, old_value, side);
    }
};

void run_MPI_bound_neuron(std::vector<Bound>& local_bounds, const std::string& fname, const std::string& start_node, const std::vector<Bound>& global_bounds) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Adjust how many neurons each rank handles
    size_t total_bounds = local_bounds.size();
    size_t neurons_per_rank = total_bounds / size;
    size_t extra_neurons = total_bounds % size;

    size_t start_idx = rank * neurons_per_rank + std::min(static_cast<size_t>(rank), extra_neurons);
    size_t end_idx = start_idx + neurons_per_rank + (rank < extra_neurons);

    // Compute bounds for local neurons
    #pragma omp parallel for
    for (size_t i = start_idx; i < end_idx; i++) {
        Bound& neuron = local_bounds[i];
        auto new_bound = bound_neuron(fname, start_node, neuron, global_bounds);
        if (neuron.side == LOWER) {
            neuron.value = std::max(neuron.value, new_bound);
        } else {
            neuron.value = std::min(neuron.value, new_bound);
        }
    }

    // Convert local_bounds to ProxyBound for MPI communication
    std::vector<ProxyBound> proxyLocalBounds;
    for (const auto& bound : local_bounds) {
        proxyLocalBounds.emplace_back(bound);
    }

    // Gather the results at rank 0
    if (rank == 0) {
        for (int src_rank = 1; src_rank < size; src_rank++) {
            size_t src_start_idx = src_rank * neurons_per_rank + std::min(static_cast<size_t>(src_rank), extra_neurons);
            size_t src_end_idx = src_start_idx + neurons_per_rank + (src_rank < extra_neurons);

            MPI_Recv(&proxyLocalBounds[src_start_idx], src_end_idx - src_start_idx, MPI_PROXY_BOUND_TYPE, src_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    } else {
        MPI_Send(&proxyLocalBounds[start_idx], end_idx - start_idx, MPI_PROXY_BOUND_TYPE, 0, 0, MPI_COMM_WORLD);
    }

    // Broadcast the updated bounds back to all ranks
    MPI_Bcast(&proxyLocalBounds[0], proxyLocalBounds.size(), MPI_PROXY_BOUND_TYPE, 0, MPI_COMM_WORLD);

    // Convert back the received ProxyBounds to Bounds
    for (size_t i = 0; i < local_bounds.size(); ++i) {
        local_bounds[i] = proxyLocalBounds[i].toBound();
    }
}

int main(int argc, char * argv[]) {
    string fname = string(prj_dir)+"/data_sets/VNN/tll_new_old.onnx";
    if(argc >= 2) {
        fname = argv[1];
    }

    MPI_Init(nullptr, nullptr);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

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
    if (rank == 0) {
        std::cout << "Optimizing layers:" << std::endl;
        for (auto l: layers_to_optimize) {
            std::cout << l->lname() << std::endl;
        }
    }

    std::vector<Bound> global_bounds;

    int rolling_horizon = 2;
    if (rolling_horizon > layers_to_optimize.size()){
        rolling_horizon = layers_to_optimize.size();
    }

    for (auto lidx = 0; lidx < layers_to_optimize.size(); lidx++) {
        auto l = layers_to_optimize[lidx];
        if (rank == 0) {
            std::cout << "################################################" << std::endl;
            std::cout << "Optimizing layer: " << l->lname() << std::endl;
            std::cout << "Layer " << lidx+1 << "/" << layers_to_optimize.size() << std::endl;
        }
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

        if (rank == 0) {
            std::cout << "Number of neurons to optimize: " << local_bounds.size()/2 << std::endl;
        }

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

        // Use the run_MPI_bound_neuron function for parallel neuron bound calculation
        offsets[0] = offsetof(ProxyBound, layer_name);
        offsets[1] = offsetof(ProxyBound, neuron_name);
        offsets[2] = offsetof(ProxyBound, value);
        offsets[3] = offsetof(ProxyBound, old_value);
        offsets[4] = offsetof(ProxyBound, side);

        MPI_Type_create_struct(nitems, blocklengths, offsets, types, &MPI_PROXY_BOUND_TYPE);
        MPI_Type_commit(&MPI_PROXY_BOUND_TYPE);

        run_MPI_bound_neuron(local_bounds, fname, start_node, global_bounds);
        auto end_time = std::chrono::high_resolution_clock::now();
        fflush(stdout);
        dup2(bak, 1);
        close(bak);

        // print out the final bounds
        if (rank == 0) {
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
        }
        global_bounds.insert(global_bounds.end(), local_bounds.begin(), local_bounds.end());
    }
    MPI_Type_free(&MPI_PROXY_BOUND_TYPE);

    if (rank == 0) {
        std::cout << "Starting final runs" << std::endl;
    }

    size_t total_objs = nn.obj_spec->shape[0];
    size_t objs_per_rank = total_objs / size;
    size_t extra_objs = total_objs % size;

    size_t start_obj_idx = rank * objs_per_rank + std::min(static_cast<size_t>(rank), extra_objs);
    size_t end_obj_idx = start_obj_idx + objs_per_rank + (rank < extra_objs);

    for (size_t obj_idx = start_obj_idx; obj_idx < end_obj_idx; obj_idx++) {
        std::cout << "########################################" << std::endl;
        final_run(fname, global_bounds, obj_idx);
    }

    // FIXME: gather final_run results

    MPI_Finalize();
    return 0;
}
