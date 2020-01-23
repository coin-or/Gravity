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
    
    auto SPP=build_pool_pqform(poolnet, solv_type);
    
    //   auto SPP=build_pool_qform(poolnet);
    
    
    
    SPP->print();
    
//    auto g = SPP->get_interaction_graph();
//    g.get_tree_decomp_bags();
//    if(g._tree){
//        DebugOn("Interaction graph is a tree!" << endl);
//    }
//    else {
//        auto bags_3d=g.decompose_bags_3d();
//    }
    auto start=get_wall_time();
    solver<> SPP_solv(SPP,solv_type);
    SPP_solv.run(5, 1e-6);
    auto end=get_wall_time();
    auto comp_time=end-start;
    SPP->project();
//    DebugOn("obj"<<SPP->_obj->eval()<<endl);
//     string result_name=string(prj_dir)+"results_pool.txt";
//    ofstream fout(result_name.c_str(), ios::app);
//    fout<<poolnet._name<<"\t"<<std::fixed<<std::setprecision(5)<<SPP->_obj->eval()<<"\t"<<std::setprecision(5)<<comp_time<<"\t"<<SPP->_status<<endl;
//    fout.close();
//    SPP->print();
    if(false){
        auto fk_old=SPP->get_obj_val();
        while (SPP->_status==0) {
            
            auto sumyk=poolnet.sumyk;
            
            auto sumyk_val=sumyk.eval();
            
            auto con=SPP->get_constraint("sumy_con");
            con->uneval();
            auto con_val=con->eval();
            
            sumyk.set_val(0, (con_val+0.1+sumyk_val));
//            auto vx= SPP->get_var<double>("x");
//            vx.initialize_all(2.0);
            
            
            SPP_solv.run(5, 1e-6);
            if(SPP->_status==0){
                auto fk_new=SPP->get_obj_val();
                if(fk_new>=fk_old)
                {
                    DebugOn("Reached same objective value again");
                    break;
                }
                fk_old=fk_new;
            }
        }
    }
    
    //    SPP->print_solution();
    //    SPP->print_constraints_stats(1e-7);
    
    
}
