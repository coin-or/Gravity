//
// Created by kbestuzheva on 12/11/17.
//
#include <stdio.h>
#include <iostream>
#include <string>
#include <stdio.h>
#include <cstring>
#include <fstream>
#include "../PoolNet.h"
#include <gravity/solver.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef USE_OPT_PARSER
#include <optionParser.hpp>
#endif

using namespace std;
using namespace gravity;


/* main */
int main (int argc, char * argv[]) {
    
    PoolNet poolnet;
    
    string fname=string(prj_dir)+"/data_sets/Pooling/Adhya1_gms.txt";
   // string fname=argv[1];

    poolnet.readgrid(fname);
    SolverType solv_type = ipopt;
    
    //do bounds on x,y,z using preprocessign in paper!
    //This is p-q-formulaiton of pooling problem!
    
    
    //    auto SPP=build_pool_qform(poolnet);
    
//    auto SPP=build_pool_pqform(poolnet, solv_type);
    
    //   auto SPP=build_pool_qform(poolnet);
    
    
    
//    SPP.print();
    Model<> SPP;
    
//    string NL_file = string(prj_dir)+"/data_sets/NL/pooling_adhya1pq.nl";
    string NL_file = string(prj_dir)+"/data_sets/NL/pooling_digabel16.nl";
    int status = SPP.readNL(NL_file);
    
//    auto g = SPP.get_interaction_graph();
//    g.get_tree_decomp_bags();
//    if(g._tree){
//        DebugOn("Interaction graph is a tree!" << endl);
//    }
//    else {
//        auto bags_3d=g.decompose_bags_3d();
//    }
//    SPP.restructure();
    SPP.project();
    bool LB = true, UB = true;
    int nb = 1e4;
    SPP.add_bound_RLTs(LB,UB,nb);
//    SPP.print();
    auto start=get_wall_time();
    solver<> SPP_solv(SPP,solv_type);
    SPP_solv.run(5, 1e-6);
    auto end=get_wall_time();
    auto comp_time=end-start;
//    SPP.print();
    auto Rel = SPP.relax(3);
//    Rel->print();
    solver<> Rel_solv(Rel,solv_type);
    Rel_solv.run(5, 1e-6);
    Rel->print_constraints_stats(1e-6);
    double upper_bound = SPP.get_obj_val(), lower_bound = Rel->get_obj_val();
    auto gap = 100*(upper_bound - lower_bound)/std::abs(upper_bound);
    DebugOn("UB = " << upper_bound << endl);
    DebugOn("LB = " << lower_bound << endl);
    DebugOn("Root node gap = " << gap << "%" << endl);
    
//    exit(1);
    double max_time = 7200,ub_solver_tol=1e-8, lb_solver_tol=1e-8, range_tol=1e-5, opt_rel_tol=1e-3, opt_abs_tol=1e6;
    unsigned max_iter=3000, nb_threads = thread::hardware_concurrency();
    SolverType ub_solver_type = ipopt, lb_solver_type = ipopt;
//        LB->print();
//        LB->initialize_zero();
//        solver<> sLB(LB,gurobi);
//        auto solve_status = sLB.run(5,1e-6);
//        CHECK(solve_status==0);
//        LB->scale_coefs(1e3);
//        LB->scale_vars(1e2);
    auto res = SPP.run_obbt(Rel, max_time, max_iter, opt_rel_tol, opt_abs_tol, nb_threads-2, ub_solver_type, lb_solver_type, ub_solver_tol, lb_solver_tol, range_tol);
//    DebugOn("obj"<<SPP._obj->eval()<<endl);
//     string result_name=string(prj_dir)+"results_pool.txt";
//    ofstream fout(result_name.c_str(), ios::app);
//    fout<<poolnet._name<<"\t"<<std::fixed<<std::setprecision(5)<<SPP._obj->eval()<<"\t"<<std::setprecision(5)<<comp_time<<"\t"<<SPP._status<<endl;
//    fout.close();
//    SPP.print();
    if(false){
        auto fk_old=SPP.get_obj_val();
        while (SPP._status==0) {
            
            auto sumyk=poolnet.sumyk;
            
            auto sumyk_val=sumyk.eval();
            
            auto con=SPP.get_constraint("sumy_con");
            con->uneval();
            auto con_val=con->eval();
            
            sumyk.set_val(0, (con_val+0.1+sumyk_val));
//            auto vx= SPP.get_var<double>("x");
//            vx.initialize_all(2.0);
            
            
            SPP_solv.run(5, 1e-6);
            if(SPP._status==0){
                auto fk_new=SPP.get_obj_val();
                if(fk_new>=fk_old)
                {
                    DebugOn("Reached same objective value again");
                    break;
                }
                fk_old=fk_new;
            }
        }
    }
    
    //    SPP.print_solution();
    //    SPP.print_constraints_stats(1e-7);
    
    
}
