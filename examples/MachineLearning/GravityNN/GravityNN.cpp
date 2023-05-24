#include <iostream>
#include <fstream>
#include <onnx.pb.h>
#include <vector>
#include <network/NeuralNet.hpp>
#include <gravity/solver.h>

using namespace gravity;

int main(int argc, char * argv[]){
    string fname = string(prj_dir)+"/data_sets/VNN/tll_bound.onnx";
    if(argc >= 2) {
        fname = argv[1];
    }

    // Empty string means we build the entire network, otherwise we build up to the specified node
    std::string final_node = "";
    NeuralNet nn(fname, final_node);

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

    NN.add(x.in(nn.indices.hidden_states));
    NN.add(y.in(nn.indices.y_ids));
    nn.initialize_state(x, y);
    nn.add_constraints(NN, x, y, nn.indices);

    ///* Objective function */
    //NN.min(
    //    x(nn.layers.back()->outputs[0]->strkey(4))
    //   -x(nn.layers.back()->outputs[0]->strkey(9))
    //);
    NN.max(
        x(nn.layers.back()->outputs[0]->strkey(0))
        // x("Gemm3/0,152")
    );

    NN.print();
    NN.write();

    solver<> S(NN,gurobi);
    auto grb_prog = (GurobiProgram*)(S._prog.get());
    auto grb_mod = grb_prog->grb_mod;
    grb_mod->set(GRB_IntParam_Threads, 8);
    grb_mod->set(GRB_IntParam_NonConvex,2);
    // grb_mod->set(GRB_IntParam_MIPFocus,3);
    // grb_mod->set(GRB_DoubleParam_BestBdStop, -1e-4);
    // grb_mod->set(GRB_DoubleParam_BestObjStop, 1e-4);

    S.run();

    NN.print_solution();

    auto sol = std::vector<double>();
    NN.get_solution(sol);
    sol.resize(nn.input_numel);

    std::cout << "Solution: " << print_vector(sol) << std::endl;
    std::cout << "Obj. value: " << std::setprecision(8) << NN.get_obj_val() << std::endl;

    return 0;
}
