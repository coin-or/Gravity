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
    int output = 0;
    PoolNet pool;
    pool.readgrid();

//do bounds on x,y,z using preprocessign in paper!
    //This is p-formulaiton of pooling problem, yet to explore, p-q and q!

    
    
    Model<double> SPP("Std-Pooling-Prob-P");
    indices Inputs=pool.Inputs;
    indices Pools=pool.Pools;
    indices Outputs=pool.Outputs;
    indices Nodes=pool.nodes;
    
    indices inputs_pools=pool.inputs_pools;
    indices pools_outputs=pool.pools_outputs;
    indices inputs_outputs=pool.inputs_outputs;
    
    
    auto out_arcs_ip=pool.out_arcs_per_node(inputs_pools);
    auto out_arcs_io=pool.out_arcs_per_node(inputs_outputs);
    auto in_arcs_ip=pool.in_arcs_per_node(inputs_pools);
    
    auto x_min=pool.x_min.in(inputs_pools);
    auto x_max=pool.x_max.in(inputs_pools);
    auto inqual=pool.inqual.in(inputs_pools);
    
    auto y_min=pool.y_min.in(pools_outputs);
    auto y_max=pool.y_max.in(pools_outputs);
    
    auto z_min=pool.z_min.in(inputs_outputs);
    auto z_max=pool.z_max.in(inputs_outputs);
    
    auto cost=pool.cost.in(Inputs);
    auto avail_min=pool.avail_min.in(Inputs);
    auto avail_max=pool.avail_max.in(Inputs);
    
    auto rev=pool.rev.in(Outputs);
    auto dem_min=pool.dem_min.in(Outputs);
    auto dem_max=pool.dem_max.in(Outputs);
    auto outqual_min=pool.outqual_min.in(Outputs);
    auto outqual_max=pool.outqual_max.in(Outputs);
    
    auto pool_cap=pool.pool_cap.in(Pools);
    
    
    var<> x("x", x_min, x_max);
    
    var<> y("y", y_min, y_max), z("z", z_min, z_max);
    SPP.add(x.in(inputs_pools));
    SPP.add(y.in(pools_outputs));
    SPP.add(z.in(inputs_outputs));
    
     //auto avail = pool.avail.in(inputs_pools);
    
//    Constraint<> mass_balance("mass_balance");
//    mass_balance=

    Constraint<> avail_ub("avail_ub");
    avail_ub=sum(x, out_arcs_ip)+sum(z, out_arcs_io)-avail_max;
    SPP.add(avail_ub.in(Inputs)<=0);
    
    Constraint<> avail_lb("avail_lb");
    avail_lb=sum(x, out_arcs_ip)+sum(z, out_arcs_io)-avail_min;
    SPP.add(avail_lb.in(Inputs)>=0);
    
    Constraint<> pool_capacity("pool_capacity");
    pool_capacity=sum(x, in_arcs_ip)-pool_cap;
    SPP.add(pool_capacity.in(Nodes)<=0);
  //  SPP.add(pool_capacity.in(*(Pools._keys))<=0);

    SPP.print();
    
}
