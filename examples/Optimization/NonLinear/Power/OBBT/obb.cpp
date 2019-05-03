//
// Created by kbestuzheva on 12/11/17.
//
#include <stdio.h>
#include <iostream>
#include <string>
#include <stdio.h>
#include <cstring>
#include <fstream>
#include "../PowerNet.h"
#include <gravity/solver.h>
#include <stdio.h>
#include <stdlib.h>
#include <optionParser.hpp>

using namespace std;
using namespace gravity;


/* main */
int main (int argc, char * argv[]) {
    int output = 0;
    bool relax = false, sdp_cuts = true, soc=true;
    bool loss_from = false, llnc=true;
    size_t num_bags = 0;
    string num_bags_s = "100";
    string solver_str = "ipopt";
    string sdp_cuts_s = "yes";
    string loss_from_s = "yes";
    string lazy_s = "no";
    bool lazy_bool = false;
    SolverType solv_type = ipopt;
    double tol = 1e-6;
    string mehrotra = "no";
    
    //string fname = string(prj_dir)+"/data_sets/Power/nesta_case5_pjm.m";
    
    string fname = string(prj_dir)+"/data_sets/Power/nesta_case5_pjm.m";
    
    // create a OptionParser with options
    op::OptionParser opt;
    opt.add_option("h", "help",
                   "shows option help"); // no default value means boolean options, which default value is false
    opt.add_option("f", "file", "Input file name", fname);
    opt.add_option("s", "solver", "Solvers: ipopt/cplex/gurobi, default = ipopt", solver_str);
    opt.add_option("b", "numbags", "Number of bags per iteration", num_bags_s);
    opt.add_option("l", "losses", "add loss constraints", loss_from_s);
    opt.add_option("lz", "lazy", "Generate 3d SDP cuts in a lazy fashion, default = no", lazy_s);
    // parse the options and verify that all went well. If not, errors and help will be shown
    bool correct_parsing = opt.parse_options(argc, argv);
    
    if (!correct_parsing) {
        return EXIT_FAILURE;
    }
    
    fname = opt["f"];
    bool has_help = op::str2bool(opt["h"]);
    if (has_help) {
        opt.show_help();
        exit(0);
    }
    solver_str = opt["s"];
    if (solver_str.compare("gurobi")==0) {
        solv_type = gurobi;
    }
    else if(solver_str.compare("cplex")==0) {
        solv_type = cplex;
    }else if(solver_str.compare("Mosek")==0) {
        solv_type = mosek;
    }
    lazy_s = opt["lz"];
    if (lazy_s.compare("no")==0) {
        lazy_bool = false;
    }
    
    loss_from_s = opt["l"];
    if (loss_from_s.compare("no")==0) {
        loss_from = false;
    }
    else {
        loss_from = true;
    }
    
    num_bags = atoi(opt["b"].c_str());
    
    cout << "\nnum bags = " << num_bags << endl;
    
    PowerNet grid;
    grid.readgrid(fname);
    grid.get_tree_decomp_bags(false,true);
    
    double solver_time_start = get_wall_time();
    vector<shared_ptr<Model<>>> batch_models;
    DebugOn("Machine has " << std::thread::hardware_concurrency() << " threads." << endl);
    int nb_threads = std::thread::hardware_concurrency();
    
    auto SDP= build_SDPOPF(grid, loss_from);
    /** Build model */
   auto var_map=SDP->_vars_name;
   auto count_map=0;
    auto count_var=0;
        for(auto it:var_map)
        {
            count_map++;
        std::string vname=it.first;
        var<> v=SDP->get_var<double>(vname);
        auto v_keys=v.get_keys();
            count_var=0;
        
        for(auto key: *v_keys)
        {
            count_var++;
            auto SDP1= build_SDPOPF(grid, loss_from);
            var<> v1=SDP1->get_var<double>(vname);
        SDP1->min(v1(key));
      //  SDP1->print();
        batch_models.push_back(SDP1);
                if (batch_models.size()==nb_threads || (count_map==var_map.size() && count_var==v.get_dim()))
                {
                    double batch_time_start = get_wall_time();
                    run_parallel(batch_models,ipopt,1e-6,nb_threads,"ma27");
                    double batch_time_end = get_wall_time();
                    auto batch_time = batch_time_end - batch_time_start;
                    DebugOn("Done running batch models, solve time = " << to_string(batch_time) << endl);
                    batch_models.clear();
                
               /// auto SDP= build_SDPOPF(grid, loss_from);
                }
            }
                    }
            

    return 0;
    
}


