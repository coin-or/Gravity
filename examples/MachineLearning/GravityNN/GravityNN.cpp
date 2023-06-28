#include <iostream>
#include <fstream>
#include <onnx.pb.h>
#include <vector>
#include <network/NeuralNet.hpp>
#include <gravity/solver.h>

using namespace gravity;

int main(int argc, char * argv[]){
    string fname = string(prj_dir)+"/data_sets/VNN/tll_new_old.onnx";
    if(argc >= 2) {
        fname = argv[1];
    }
    int idx = 0;
    if (argc >= 3) {
        idx = atoi(argv[2]);
    }

    // Empty string means we build the entire network, otherwise we build up to the specified node
    std::string start_node = "";
    std::string final_node = "";
    NeuralNet nn(fname);

    Model<>& NN = nn.build_model(idx, start_node, final_node);

    if (idx < 0) {
        auto tensor = nn.subgraph.back()->outputs.at(0);
        gravity::func<> expr = 0.0;
        for (auto index: ShapeIter(tensor->shape)) {
            expr += nn.x(tensor->strkey(index));
        }
        NN.max(expr);
    }

    // NN.write();

    solver<> S(NN,gurobi);
    auto grb_prog = (GurobiProgram*)(S._prog.get());
    auto grb_mod = grb_prog->grb_mod;
    grb_mod->set(GRB_IntParam_Threads, get_num_threads() / 2);
    grb_mod->set(GRB_IntParam_NonConvex,2);
    grb_mod->set(GRB_IntParam_Presolve,2);
    // grb_mod->set(GRB_DoubleParam_BestBdStop, -1e-4);

    S.run();

    auto sol = std::vector<double>();
    NN.print_solution();
    NN.get_solution(sol);
    sol.resize(nn.input_numel);

    std::cout << "Solution: " << print_vector(sol) << std::endl;
    std::cout << "Obj. value: " << std::setprecision(8) << NN.get_obj_val() << std::endl;

    return 0;
}
