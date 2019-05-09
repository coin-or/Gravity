//
// Created by kbestuzheva on 12/11/17.
//
#include <stdio.h>
#include <iostream>
#include <string>
#include <stdio.h>
#include <cstring>
#include <fstream>
#include <utility>
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
    
    double start_tol=1E-3, interval_tol=1E-3;
    
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
    
    
    
    
    
    DebugOn("Machine has " << thread::hardware_concurrency() << " threads." << endl);
    /*Code Assumes number of threads is even and greater than or equal to 2*/
    int nb_threads = thread::hardware_concurrency();
    double upper_bound = grid.solve_acopf();
    auto SDP= build_SDPOPF(grid, loss_from, upper_bound);
    solver<> SDPLB(SDP,solv_type);
    SDPLB.run(output = 5, tol = 1e-6);
    double lower_bound=SDP->get_obj_val();
    
    double solver_time_start = get_wall_time();
    vector<shared_ptr<Model<>>> batch_models;
    map<pair<string, string>, int> fixed_point;
    string vname;
    string mname;
    string key;
    string dir_array[2]={"LB", "UB"};
    bool terminate=false;
    const double upp_low_tol=1e-3, range_tol=0.02;
    if (upper_bound-lower_bound>=upp_low_tol)
    {
        for(auto &it:SDP->_vars_name)
        {
            string vname=it.first;
            var<> v=SDP->get_var<double>(vname);
            auto v_keys=v.get_keys();
            for(auto &key: *v_keys)
            {
                auto p=make_pair(vname, key);
                fixed_point[p]=false;
            }
            
        }
        //        if(SDP->_obj->is_linear())
        //        {
        //            Constraint<> obj_lower("obj_lower");
        //            obj_lower = *(SDP->_obj);
        //            SDP->add(obj_lower>=lower_bound);
        //        }
        
        while(!terminate)
        {
            for(auto &it:SDP->_vars_name)
            {
                vname=it.first;
                var<> v=SDP->get_var<double>(vname);
                for(auto &key: *(v.get_keys()))
                {
                    auto p=make_pair(vname, key);
                    if(fixed_point[p]==false)
                    {
                        for(auto &dir: dir_array)
                        {
                            auto modelk = make_shared<Model<>>(*SDP);
                            // auto modelk=build_SDPOPF(grid, loss_from, upper_bound);
                            mname=vname+"."+key+"."+dir;
                            modelk->set_name(mname);
                            if(dir=="LB")
                            {
                                modelk->min(v(key));
                            }
                            else
                            {
                                modelk->max(v(key));
                            }
                            batch_models.push_back(modelk);
                            if (batch_models.size()==nb_threads || (it==*SDP->_vars_name.end() && key==v.get_keys()->back() && dir=="UB"))
                            {
                                double batch_time_start = get_wall_time();
                                run_parallel(batch_models,ipopt,1e-6,nb_threads,"ma27");
                                double batch_time_end = get_wall_time();
                                auto batch_time = batch_time_end - batch_time_start;
                                DebugOn("Done running batch models, solve time = " << to_string(batch_time) << endl);
                                
                            }
                        }
                    }
                }
            }
            terminate = true;
        }
    }
    return 0;
    
}


