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
    poolnet.readgrid();
    SolverType solv_type = gurobi;
    
    //do bounds on x,y,z using preprocessign in paper!
    //This is p-q-formulaiton of pooling problem!
    
  //  auto SPP=build_pool_qform(poolnet);
     auto SPP=build_pool_pform(poolnet);
    
   
    
  SPP->print();

    solver<> SPP_solv(SPP,solv_type);
    SPP_solv.run(5, 1e-7);

    SPP->print_solution();
    SPP->print_constraints_stats(1e-7);
    
    
    
}
