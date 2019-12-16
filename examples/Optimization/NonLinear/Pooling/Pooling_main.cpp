//
// Created by kbestuzheva on 12/11/17.
//
#include <stdio.h>
#include <iostream>
#include <string>
#include <stdio.h>
#include <cstring>
#include <fstream>
#include "PoolNet.h"
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
     //string fname=string(prj_dir)+"/data_sets/Pooling/Adhya1_gms.txt";
    string fname="/Users/smitha/Desktop/Pooling_instances/Adhya3_gms.txt";
    poolnet.readgrid(fname);
    SolverType solv_type = ipopt;
    
    //do bounds on x,y,z using preprocessign in paper!
    //This is p-q-formulaiton of pooling problem!
    
//    auto SPP=build_pool_qform(poolnet);
     auto SPP=build_pool_pform(poolnet, solv_type);
    

    
    SPP->print();
    auto g = SPP->get_interaction_graph();
    g.print();
    if(g._tree){
        DebugOn("Interaction graph is a tree!" << endl);
    }
    else {
        g.get_tree_decomp_bags();
        auto bags_3d=g.decompose_bags_3d();
        
    }

    solver<> SPP_solv(SPP,solv_type);
        SPP_solv.run(5, 1e-7);
            auto fk_old=SPP->get_obj_val();
    while (SPP->_status==0) {

        auto sumyk=poolnet.sumyk;
        auto sumyk_val=sumyk.eval(0);
    
        auto con=SPP->get_constraint("sumy_con");
            con->uneval();
        auto con_val=con->eval(0);

            sumyk.set_val(0, (con_val+0.1+sumyk_val));
        //SPP->reset_constrs();
        SPP->print();
    SPP_solv.run(5, 1e-7);
      SPP->print_solution();
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

    


    SPP->print_solution();
    SPP->print_constraints_stats(1e-7);
    
    
}
