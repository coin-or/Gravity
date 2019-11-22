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
    PoolNet poolnet;
    poolnet.readgrid();

//do bounds on x,y,z using preprocessign in paper!
    //This is p-formulaiton of pooling problem, yet to explore, p-q and q!

    
    
    Model<> SPP("Std-Pooling-Prob-P");
    indices Inputs=poolnet.Inputs;
    indices Pools=poolnet.Pools;
    indices Outputs=poolnet.Outputs;
    //indices Nodes=pool.nodes;
    
    indices inputs_pools=poolnet.inputs_pools;
    indices pools_outputs=poolnet.pools_outputs;
    indices inputs_outputs=poolnet.inputs_outputs;
    
    auto out_arcs_to_pool_per_input = poolnet.out_arcs_to_pool_per_input();
    auto out_arcs_to_output_per_input = poolnet.out_arcs_to_output_per_input();
    auto in_arcs_per_pool = poolnet.in_arcs_per_pool();
    auto out_arcs_per_pool = poolnet.out_arcs_per_pool();
    auto in_arcs_from_pool_per_output = poolnet.in_arcs_from_pool_per_output();
    auto in_arcs_from_input_per_output = poolnet.in_arcs_from_input_per_output();

    auto x_min=poolnet.x_min.in(inputs_pools);
    auto x_max=poolnet.x_max.in(inputs_pools);
    
    
    auto y_min=poolnet.y_min.in(pools_outputs);
    auto y_max=poolnet.y_max.in(pools_outputs);
    
    auto z_min=poolnet.z_min.in(inputs_outputs);
    auto z_max=poolnet.z_max.in(inputs_outputs);
    
    auto cost=poolnet.cost.in(Inputs);
    auto avail_min=poolnet.avail_min.in(Inputs);
    auto avail_max=poolnet.avail_max.in(Inputs);
    auto p_in=poolnet.inqual.in(Inputs);
    
    auto rev=poolnet.rev.in(Outputs);
    auto dem_min=poolnet.dem_min.in(Outputs);
    auto dem_max=poolnet.dem_max.in(Outputs);
    auto p_out_min=poolnet.outqual_min.in(Outputs);
    auto p_out_max=poolnet.outqual_max.in(Outputs);
    
    auto pool_cap=poolnet.pool_cap.in(Pools);
    
    
    var<> x("x", x_min, x_max);
    
    var<> y("y", y_min, y_max), z("z", z_min, z_max);
    var<> p_pool("p_pool", 0, 1);
    SPP.add(x.in(inputs_pools));
    SPP.add(y.in(pools_outputs));
    SPP.add(z.in(inputs_outputs));
    SPP.add(p_pool.in(Pools));
    
    Constraint<> avail_lb("avail_lb");
    avail_lb=sum(x, out_arcs_to_pool_per_input)+sum(z, out_arcs_to_output_per_input)-avail_min;
    SPP.add(avail_lb.in(Inputs)>=0);
    

    Constraint<> avail_ub("avail_ub");
    avail_ub=sum(x, out_arcs_to_pool_per_input)+sum(z, out_arcs_to_output_per_input)-avail_max;
    SPP.add(avail_ub.in(Inputs)<=0);
    
 
    Constraint<> pool_capacity("pool_capacity");
    pool_capacity=sum(x, in_arcs_per_pool)-pool_cap;
    SPP.add(pool_capacity.in(Pools)<=0);
  
    
    Constraint<> demand_lb("demand_lb");
    demand_lb=sum(y, in_arcs_from_pool_per_output)+sum(z,in_arcs_from_input_per_output)-dem_min;
    SPP.add(demand_lb.in(Outputs)>=0);
    
    Constraint<> demand_ub("demand_ub");
    demand_ub=sum(y, in_arcs_from_pool_per_output)+sum(z,in_arcs_from_input_per_output)-dem_max;
    SPP.add(demand_ub.in(Outputs)<=0);
    
    Constraint<> mass_balance("mass_balance");
    mass_balance=sum(x, in_arcs_per_pool)-sum(y, out_arcs_per_pool);
    SPP.add(mass_balance.in(Pools)==0);
    
        int row_id = 0;
        indices pool_matrix = indices("pool_matrix");
        for (const string& pool_key:*Pools._keys) {
            for (auto i = 0; i<Outputs.size(); i++) {
                pool_matrix.add_in_row(row_id, pool_key);
            }
            row_id++;
        }

        Constraint<> quality_balance("quality_balance");//TODO debug transpose version
        quality_balance=x.in(in_arcs_per_pool)*p_in - p_pool.in(pool_matrix)*y.in(out_arcs_per_pool);// - p_pool* sum(y, out_arcs_per_pool)
        SPP.add(quality_balance.in(Pools)==0);
    
    row_id = 0;
    indices outpool_matrix = indices("outpool_matrix");
    for (const string& out_key:*Outputs._keys) {
        for (auto i = 0; i<Pools.size(); i++) {
            outpool_matrix.add_in_row(row_id, out_key);
        }
        row_id++;
    }
    
    row_id = 0;
    indices outinput_matrix = indices("outinput_matrix");
    for (const string& out_key:*Outputs._keys) {
        for (auto i = 0; i<Inputs.size(); i++) {
            outinput_matrix.add_in_row(row_id, out_key);
        }
        row_id++;
    }
    
//    Constraint<> product_quality("product_quality");//TODO debug transpose version
//    product_quality=y.in(in_arcs_from_pool_per_output)*p_pool+z.in(in_arcs_from_input_per_output)*p_in - p_out.in(outpool_matrix)*y.in(in_arcs_from_pool_per_output)-p_out.in(outinput_matrix)*z.in(in_arcs_from_input_per_output);
//    SPP.add(product_quality.in(Outputs)==0);
//
//    Constraint<> product_quality_min("product_quality_min");
//    product_quality_min=p_out-p_out_min;
//    SPP.add(product_quality_min.in(Outputs)>=0);
    
    
    
    Constraint<> product_quality_lb("product_quality_lb");//TODO debug transpose version and propagate matrix indexing to function
    product_quality_lb=y.in(in_arcs_from_pool_per_output)*p_pool+z.in(in_arcs_from_input_per_output)*(p_in-p_out_min.in(outinput_matrix)).in(outinput_matrix) - p_out_min.in(outpool_matrix)*y.in(in_arcs_from_pool_per_output);//-p_out_min.in(outinput_matrix)*z.in(in_arcs_from_input_per_output);
    SPP.add(product_quality_lb.in(Outputs)>=0);
    
    Constraint<> product_quality_ub("product_quality_ub");//TODO debug transpose version and propagate matrix indexing to function
    product_quality_ub=y.in(in_arcs_from_pool_per_output)*p_pool+z.in(in_arcs_from_input_per_output)*(p_in-p_out_max.in(outinput_matrix)).in(outinput_matrix) - p_out_max.in(outpool_matrix)*y.in(in_arcs_from_pool_per_output);//-p_out_min.in(outinput_matrix)*z.in(in_arcs_from_input_per_output);
    SPP.add(product_quality_ub.in(Outputs)<=0);
    



    SPP.print();
    
}
