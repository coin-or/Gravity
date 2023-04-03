#include <string>
#include <gravity/solver.h>
#include "MPSParser.hpp"

int main(int argc, char** argv) {
    string fname = string(prj_dir)+"/data_sets/VNN/mip_lay-16_0.mps";
    Model<> model("new_mip");
    model.readMPS(fname);
//    model = model_from_file(fname);
    solver<> MIP_Solver(model, gurobi);
    MIP_Solver.run();
//    model.print();
    model.print_symbolic();
    model.print_solution();
    model.print_int_solution();

    // ReLU upper bound: ReLUX_Y
    // ReLU Bool: aReLUX_Y

    // for (auto v : model._vars_name) {
        // std::cout << v.first << std::endl;
    // }
    
    
	return 0;
}
