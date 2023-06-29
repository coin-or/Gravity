#include <iostream>
#include <fstream>
#include <onnx.pb.h>
#include <vector>
#include <network/NeuralNet.hpp>
#include <gravity/solver.h>
#include <CLI/CLI.hpp>

using namespace gravity;

class Config {
public:
    std::string fname;
    std::string start_node = "";
    std::string final_node = "";
    bool w_gurobi = false;
    bool w_gravity = false;
    int obj_idx;

    Config(int argc, char * argv[]): app("GravityNN", "GravityNN") {
        this->app.add_option("-f,--file", this->fname, "ONNX file path")
            ->required()
            ->check(CLI::ExistingFile);

        this->app.add_option("-s,--start",     this->start_node, "Start node name");
        this->app.add_option("-e,--end",       this->final_node, "Final node name");
        this->app.add_option("-i,--index",     this->obj_idx,    "Objective index (-1 will sum output of last layer)")
            ->required()
            ->check(CLI::Number);

        this->app.add_flag("-g,--gurobi",      this->w_gurobi,   "Write Gurobi model");
        this->app.add_flag("-v,--gravity",     this->w_gravity,  "Write Gravity model");

        try {
            this->app.parse(argc, argv);
        } catch(const CLI::ParseError &e) {
            std::exit(this->app.exit(e));
        }
    }

    CLI::App app;
};

int main(int argc, char * argv[]){
    Config config(argc, argv);
    NeuralNet nn(config.fname);

    Model<>& NN = nn.build_model(config.obj_idx, config.start_node, config.final_node);


    solver<> S(NN,gurobi);
    auto grb_prog = (GurobiProgram*)(S._prog.get());
    auto grb_mod = grb_prog->grb_mod;
    grb_mod->set(GRB_IntParam_Threads, get_num_threads() / 2);
    grb_mod->set(GRB_IntParam_NonConvex,2);
    grb_mod->set(GRB_IntParam_Presolve,2);
    // grb_mod->set(GRB_DoubleParam_BestBdStop, -1e-4);

    if (config.w_gravity) {
        NN.write();
    }
    if (config.w_gurobi) {
        grb_mod->write("gurobiprint.lp");
    }

    S.run();

    auto sol = std::vector<double>();
    NN.print_solution();
    NN.get_solution(sol);
    sol.resize(nn.input_numel);

    std::cout << "Solution: " << print_vector(sol) << std::endl;
    std::cout << "Obj. value: " << std::setprecision(8) << NN.get_obj_val() << std::endl;

    return 0;
}
