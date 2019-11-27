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
    indices Attr=poolnet.Attr;
    
    //indices Nodes=pool.nodes;
    
    indices inputs_pools=poolnet.inputs_pools;
    indices pools_outputs=poolnet.pools_outputs;
    indices inputs_outputs=poolnet.inputs_outputs;
    indices inputs_attr=poolnet.inputs_attr;
    indices outputs_attr=poolnet.outputs_attr;
    indices pool_attr = indices(Pools,Attr);
    
    auto out_arcs_to_pool_per_input = poolnet.out_arcs_to_pool_per_input();
    auto out_arcs_to_output_per_input = poolnet.out_arcs_to_output_per_input();
    auto in_arcs_per_pool = poolnet.in_arcs_per_pool();
    auto in_arcs_per_pool_attr = poolnet.in_arcs_per_pool_attr();
    auto in_arcs_attr_per_pool = poolnet.in_arcs_attr_per_pool();
    auto out_arcs_per_pool = poolnet.out_arcs_per_pool();
    auto out_arcs_per_pool_attr = poolnet.out_arcs_per_pool_attr();
    auto in_arcs_from_pool_per_output = poolnet.in_arcs_from_pool_per_output();
    auto in_arcs_from_pool_per_output_attr=poolnet.in_arcs_from_pool_per_output_attr();
    auto in_arcs_from_input_per_output = poolnet.in_arcs_from_input_per_output();
    auto in_arcs_from_input_per_output_attr = poolnet.in_arcs_from_input_per_output_attr();
    
    auto x_min=poolnet.x_min.in(inputs_pools);
    auto x_max=poolnet.x_max.in(inputs_pools);
    
    
    auto y_min=poolnet.y_min.in(pools_outputs);
    auto y_max=poolnet.y_max.in(pools_outputs);
    
    auto z_min=poolnet.z_min.in(inputs_outputs);
    auto z_max=poolnet.z_max.in(inputs_outputs);
    
    // auto cost=poolnet.cost.in(Inputs);
    auto avail_min=poolnet.avail_min.in(Inputs);
    auto avail_max=poolnet.avail_max.in(Inputs);
    auto p_in=poolnet.inqual.in(inputs_attr);
    
    auto rev=poolnet.rev.in(Outputs);
    auto dem_min=poolnet.dem_min.in(Outputs);
    auto dem_max=poolnet.dem_max.in(Outputs);
    auto p_out_min=poolnet.outqual_min.in(outputs_attr);
    auto p_out_max=poolnet.outqual_max.in(outputs_attr);
    
    auto pool_cap=poolnet.pool_cap.in(Pools);
    
    auto cost_ip=poolnet.cost_ip.in(inputs_pools);
    auto cost_io=poolnet.cost_io.in(inputs_outputs);
    auto cost_po=poolnet.cost_po.in(pools_outputs);
    
    var<> x("x", x_min, x_max);
    
    var<> y("y", y_min, y_max), z("z", z_min, z_max);
    var<> p_pool("p_pool", 0, 5);
    SPP.add(x.in(inputs_pools));
    SPP.add(y.in(pools_outputs));
    SPP.add(z.in(inputs_outputs));
    SPP.add(p_pool.in(pool_attr));
    
    x.initialize_all(2.0);
    y.initialize_all(1.0);
    
    Constraint<> avail_lb("avail_lb");
    avail_lb=sum(x, out_arcs_to_pool_per_input)+sum(z, out_arcs_to_output_per_input)-avail_min;
    //    avail_lb=sum(x, out_arcs_to_pool_per_input)-avail_min;
    SPP.add(avail_lb.in(Inputs)>=0);
    SPP.print();
    
    
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
    for (const string& pool_key:*pool_attr._keys) {
        auto pos = nthOccurrence(pool_key, ",", 1);
        auto pool = pool_key.substr(0,pos);
        pool_matrix.add_empty_row();
        for (auto &out:*Outputs._keys) {
            auto po_out=pool+","+out;
            if(pools_outputs._keys_map->find(po_out) !=pools_outputs._keys_map->end()){
                pool_matrix.add_in_row(row_id, pool_key);
            }
        }
        row_id++;
    }
    
    
    Constraint<> quality_balance("quality_balance");//TODO debug transpose version
    quality_balance=p_in.in(in_arcs_attr_per_pool)*x.in(in_arcs_per_pool_attr) - p_pool.in(pool_matrix)*y.in(out_arcs_per_pool_attr);// - p_pool* sum(y, out_arcs_per_pool)
    SPP.add(quality_balance.in(pool_attr)==0);
    SPP.print();
    
    
    
    
    
    row_id = 0;
    indices pool_attr_per_output_attr_matrix = indices("pool_attr_per_output_attr_matrix");
    for (const string& outat_key:*outputs_attr._keys) {
        auto pos = nthOccurrence(outat_key, ",", 1);
        auto out=outat_key.substr(0, pos);
        auto attr = outat_key.substr(pos+1);
        pool_attr_per_output_attr_matrix.add_empty_row();
        for (auto &pool:*Pools._keys) {
            auto po_out= pool+","+out;
            if(pools_outputs._keys_map->find(po_out) !=pools_outputs._keys_map->end()){
                auto key=pool+","+attr;
                pool_attr_per_output_attr_matrix.add_in_row(row_id, key);
            }
        }
        row_id++;
    }
    
    row_id = 0;
    indices input_attr_per_output_attr_matrix = indices("input_attr_per_output_attr_matrix");
    for (const string& outat_key:*outputs_attr._keys) {
        auto pos = nthOccurrence(outat_key, ",", 1);
        auto out=outat_key.substr(0, pos);
        auto attr = outat_key.substr(pos+1);
        input_attr_per_output_attr_matrix.add_empty_row();
        for (auto &input:*Inputs._keys) {
            auto io_out= input+","+out;
            if(inputs_outputs._keys_map->find(io_out) !=inputs_outputs._keys_map->end()){
                auto key=input+","+attr;
                input_attr_per_output_attr_matrix.add_in_row(row_id, key);
            }
        }
        row_id++;
    }
    
    row_id = 0;
    indices output_attr_per_ypo_matrix = indices("output_attr_per_ypo_matrix");
    for (const string& outat_key:*outputs_attr._keys) {
        auto pos = nthOccurrence(outat_key, ",", 1);
        auto out = outat_key.substr(0,pos);
        output_attr_per_ypo_matrix.add_empty_row();
        for (auto &pool:*Pools._keys) {
            auto po_out=pool+","+out;
            if(pools_outputs._keys_map->find(po_out) !=pools_outputs._keys_map->end()){
                output_attr_per_ypo_matrix.add_in_row(row_id, outat_key);
            }
        }
        row_id++;
    }
    
    row_id = 0;
    indices output_attr_per_zio_matrix = indices("output_attr_per_zio_matrix");
    for (const string& outat_key:*outputs_attr._keys) {
        auto pos = nthOccurrence(outat_key, ",", 1);
        auto out = outat_key.substr(0,pos);
        output_attr_per_zio_matrix.add_empty_row();
        for (auto &input:*Inputs._keys) {
            auto io_out=input+","+out;
            if(inputs_outputs._keys_map->find(io_out) !=inputs_outputs._keys_map->end()){
                output_attr_per_zio_matrix.add_in_row(row_id, outat_key);
            }
        }
        row_id++;
    }
    
    
    //
    //    Constraint<> product_quality_lb("product_quality_lb");//TODO debug transpose version and propagate matrix indexing to function
    //    product_quality_lb=y.in(in_arcs_from_pool_per_output)*p_pool+z.in(in_arcs_from_input_per_output)*(p_in-p_out_min.in(outinput_matrix)).in(outinput_matrix) - p_out_min.in(outpool_matrix)*y.in(in_arcs_from_pool_per_output);//-p_out_min.in(outinput_matrix)*z.in(in_arcs_from_input_per_output);
    //    SPP.add(product_quality_lb.in(Outputs)>=0);
    //
    //    Constraint<> product_quality_ub("product_quality_ub");//TODO debug transpose version and propagate matrix indexing to function
    //    product_quality_ub=y.in(in_arcs_from_pool_per_output)*p_pool+z.in(in_arcs_from_input_per_output)*(p_in-p_out_max.in(outinput_matrix)).in(outinput_matrix) - p_out_max.in(outpool_matrix)*y.in(in_arcs_from_pool_per_output);//-p_out_min.in(outinput_matrix)*z.in(in_arcs_from_input_per_output);
    //    SPP.add(product_quality_ub.in(Outputs)<=0);
    
    
    Constraint<> product_quality_ub("product_quality_ub");//TODO debug transpose version and propagate matrix indexing to function
   
//    product_quality_ub = z.in(in_arcs_from_input_per_output_attr);
    product_quality_ub=y.in(in_arcs_from_pool_per_output_attr)*p_pool.in(pool_attr_per_output_attr_matrix)+z.in(in_arcs_from_input_per_output_attr)*p_in.in(input_attr_per_output_attr_matrix)-p_out_max.in(output_attr_per_ypo_matrix)*y.in(in_arcs_from_pool_per_output_attr)-p_out_max.in(output_attr_per_zio_matrix)*z.in(in_arcs_from_input_per_output_attr);
    SPP.add(product_quality_ub.in(outputs_attr)<=0);
    
    Constraint<> sumy("sumy");
    sumy=sum(y);
    //SPP.add(sumy>=1);
    
    auto obj= product(cost_ip, x)+product(cost_io, z)+product(cost_po, y);
    SPP.min(obj);
    
    SPP.print();
    
    solver<> SPP_solv(SPP,ipopt);
    SPP_solv.run(output = 5, 1e-7);
    
    
    
    
}
