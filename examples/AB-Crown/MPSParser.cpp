#include <string>
#include <gurobi_c++.h>
#include <gravity/solver.h>
#include "MPSParser.hpp"

int main(int argc, char** argv) {
    Model<> model = model_from_file("/mnt/trail_test/haydnj/vnn/Gravity/new_mip.mps");
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