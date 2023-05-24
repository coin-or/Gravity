#include <iostream>
#include <fstream>
#include <onnx.pb.h>
#include <vector>
#include <network/NeuralNet.hpp>
#include <gravity/solver.h>

using namespace gravity;

std::vector<std::pair<std::string, double>> get_solution(GRBModel* model, GurobiProgram* prog) {
    size_t vid, vid_inst;
    GRBVar gvar;
    param_* v;
    std::vector<std::pair<std::string, double>> sol;
    for (auto& v_p: prog->_model->_vars)
    {
        v = v_p.second.get();
        auto idx = v->get_id();
        auto dim = v->_dim[0];
        for (auto i = 0; i < dim; i++) {
            auto vid = idx + v->get_id_inst(i);
            gvar = prog->_grb_vars.at(vid);

            auto val = gvar.get(GRB_DoubleAttr_ScenNX);
            auto name = v->get_name(i);
            sol.push_back({name, val});
        }
    }
    return sol;
}

template<typename T = double>
std::pair<GRBLinExpr, std::vector<GRBVar>> build_objective(const func<T>& obj, GurobiProgram* grb_prog) {
    GRBLinExpr lterm;
    GRBVar gvar1;
    lterm = 0;

    std::vector<GRBVar> vars;
    for (auto& it1: obj.get_lterms()) {
        gvar1 = grb_prog->_grb_vars[it1.second._p->get_id() + it1.second._p->get_id_inst()];
        lterm += 1*gvar1;
        vars.push_back(gvar1);
    }
    return {lterm, vars};
}

int main(int argc, char * argv[]){
    string fname = string(prj_dir)+"/data_sets/VNN/tll_bound.onnx";
    if(argc >= 2) {
        fname = argv[1];
    }

    // Empty string means we build the entire network, otherwise we build up to the specified node
    std::string final_node = "/Gemm_2/Gemm";
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
    NN.write();

    solver<> S(NN,gurobi);
    auto grb_prog = (GurobiProgram*)(S._prog.get());
    auto grb_mod = grb_prog->grb_mod;
    grb_mod->set(GRB_IntParam_Threads, 1);

    grb_prog->_model->fill_in_maps();
    grb_prog->_model->compute_funcs();
    grb_prog->fill_in_grb_vmap();
    grb_prog->create_grb_constraints();
    grb_mod->set(GRB_IntParam_OutputFlag, 1);
    grb_mod->set(GRB_IntParam_MIPFocus,3);
    // grb_mod->set(GRB_DoubleParam_BestBdStop, -1e-4);
    // grb_mod->set(GRB_DoubleParam_BestObjStop, 1e-4);

    gravity::func<> obj = x("/Gemm_2/Gemm_output_0,73");
    size_t nvals = 2;

    auto out = build_objective(obj, grb_prog);
    auto grb_obj = out.first;
    auto vars = out.second;

    grb_mod->set(GRB_IntAttr_NumScenarios, nvals);
    grb_mod->setObjective(grb_obj, GRB_MAXIMIZE);

    // Scen 0 (ub)
    grb_mod->set(GRB_IntParam_ScenarioNumber, 0);
    vars[0].set(GRB_DoubleAttr_ScenNObj, 1.0);
    // Scen 1 (lb)
    grb_mod->set(GRB_IntParam_ScenarioNumber, 1);
    vars[0].set(GRB_DoubleAttr_ScenNObj, -1.0);

    grb_mod->set(GRB_DoubleParam_TimeLimit, 10.0);
    grb_mod->update();
    grb_mod->write("multi.lp");
    grb_mod->optimize();

    // Print out obj vals for each scenario
    for (int i = 0; i < nvals; i++) {
        if(grb_mod->get(GRB_IntAttr_SolCount) == 0) {
            continue;
        }

        grb_mod->set(GRB_IntParam_ScenarioNumber, i);
        // grb_prog->update_solution();
        double scenarioObjVal   = grb_mod->get(GRB_DoubleAttr_ScenNObjVal);
        double scenarioObjBound = grb_mod->get(GRB_DoubleAttr_ScenNObjBound);
        std::cout << "######################################################" << std::endl;
        std::cout << "Scen " << i << ": " << scenarioObjVal << " (" << scenarioObjBound << ")" << std::endl;
        std::cout << "Solution: " << std::endl;
        // auto soln = get_solution(grb_mod, grb_prog);
        // for (auto& s: soln) {
            // std::cout << s.first << ": " << s.second << std::endl;
        // }
    }

    return 0;
}
