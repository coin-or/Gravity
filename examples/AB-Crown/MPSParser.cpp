#include <string>
#include <gurobi_c++.h>
#include <gravity/solver.h>
#include "MPSParser.hpp"

int main(int argc, char** argv) {
    string fname = string(prj_dir)+"/data_sets/VNN/new_mip.mps";
    Model<> model = model_from_file(fname);
    solver<> MIP_Solver(model, ipopt);
    MIP_Solver.run();

    // ReLU upper bound: ReLUX_Y
    // ReLU Bool: aReLUX_Y

    // for (auto v : model._vars_name) {
        // std::cout << v.first << std::endl;
    // }
    //model.print();
	return 0;
}
