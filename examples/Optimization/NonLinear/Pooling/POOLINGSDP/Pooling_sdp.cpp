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
    auto SPP= make_shared<Model<>>("Std-Pooling-Prob-PQ");
    
    PoolNet poolnet;
    
    string fname=string(prj_dir)+"/data_sets/Pooling/Adhya1_gms.txt";
    poolnet.readgrid(fname);
    SolverType solv_type = ipopt;
    
    indices I=poolnet.Inputs;
    indices L=poolnet.Pools;
    indices J=poolnet.Outputs;
    indices K=poolnet.Attr;
    
    //indices Nodes=pool.nodes;
    
    
    indices Tx=poolnet.inputs_pools;
    indices Ty=poolnet.pools_outputs;
    indices Tz=poolnet.inputs_outputs;
    indices I_K=indices(I,K);
    indices J_K=indices(J,K);
    indices L_K = indices(L,K);
    
    indices inputs_pools_outputs=poolnet.inputs_pools_outputs;
    auto out_arcs_to_output_per_input = poolnet.out_arcs_to_output_per_input();
    auto in_arcs_from_input_per_output=poolnet.in_arcs_from_input_per_output();
    auto in_arcs_from_input_per_output_attr=poolnet.in_arcs_from_input_per_output_attr();
    
    indices pool_x_matrix = poolnet.pool_x_matrix_fill();
    indices input_x_matrix = poolnet.input_x_matrix_fill();
    indices output_x_matrix = poolnet.output_x_matrix_fill();
    auto x_p_indices =poolnet.outattr_x_p_matrix_fill();
    indices outattr_x_matrix= x_p_indices[0];
    indices outattr_pin_matrix= x_p_indices[1];
    indices outattr_pout_matrix=x_p_indices[2];
    auto p_indices=poolnet.outattrz_p_matrix_fill();
    indices outattrz_pin_matrix = p_indices[0];
    indices outattrz_pout_matrix = p_indices[1];
    indices pool_q_matrix = poolnet.pool_q_matrix_fill();
    indices pooloutput_x_matrix = poolnet.pooloutput_x_matrix_fill();
    
    auto x_q_indices=poolnet.inputpool_x_q_S_matrix_fill();
    indices inputpool_x_matrix=x_q_indices[0];
    indices inputpool_q_matrix=x_q_indices[1];
    indices inputpool_poolcap_matrix=x_q_indices[2];
    
    auto y_q_indices=poolnet.inpoolout_y_q_matrix_fill();
    indices inpoolout_y_matrix = y_q_indices[0];
    indices inpoolout_q_matrix = y_q_indices[1];
    
    auto x_c_indices=poolnet.inpoolout_x_c_matrix_fill();
    indices inpoolout_x_matrix=x_c_indices[0];
    indices inpoolout_cip_matrix=x_c_indices[1];
    indices inpoolout_cpo_matrix=x_c_indices[2];
    
    
    
    
    // auto cost=poolnet.cost.in(Inputs);
    auto A_L=poolnet.A_L.in(I);
    auto A_U=poolnet.A_U.in(I);
    auto C=poolnet.C.in(I_K);
    
    auto D_L=poolnet.D_L.in(J);
    auto D_U=poolnet.D_U.in(J);
    auto P_U=poolnet.P_U.in(J_K);
    
    auto S=poolnet.S.in(L);
    
    auto c_tx=poolnet.c_tx.in(Tx);
    auto c_tz=poolnet.c_tz.in(Tz);
    auto c_ty=poolnet.c_ty.in(Ty);
    auto sumyk=poolnet.sumyk;
    
    auto x_min=poolnet.x_min.in(inputs_pools_outputs);
    auto x_max=poolnet.x_max.in(inputs_pools_outputs);
    
    
    auto y_min=poolnet.y_min.in(Ty);
    auto y_max=poolnet.y_max.in(Ty);
    
    auto z_min=poolnet.z_min.in(Tz);
    auto z_max=poolnet.z_max.in(Tz);
    
    var<> x("x",x_min, x_max), y("y", y_min, y_max);
    var<> q("q", 0, 1), z("z", z_min, z_max);
    var<> objvar("objvar");
    
    SPP->add(x.in(inputs_pools_outputs));
    SPP->add(q.in(Tx));
    SPP->add(z.in(Tz));
    SPP->add(y.in(Ty));
    SPP->add(objvar);
    
    //    SPP->add(sumyk);
    //    sumyk.set_lb(0);
    // SPP->initialize_zero();
    x.initialize_all(1.0);
    y.initialize_all(1.0);
    z.initialize_all(1.0);
    q.initialize_all(0.5);
    
    
    Constraint<> feed_availability("feed_availability");
    feed_availability=sum(x, input_x_matrix)+sum(z, out_arcs_to_output_per_input)-A_U;
    SPP->add(feed_availability.in(I)<=0);
    
    Constraint<> pool_capacity("pool_capacity");
    pool_capacity=sum(x, pool_x_matrix)-S;
    SPP->add(pool_capacity.in(L)<=0);
    
    Constraint<> product_demand("product_demand");
    product_demand=sum(x, output_x_matrix)+sum(z,in_arcs_from_input_per_output)-D_U;
    SPP->add(product_demand.in(J)<=0);
    
    
    Constraint<> product_quality("product_quality");//TODO debug transpose version
    product_quality=(C.in(outattr_pin_matrix)-P_U.in(outattr_pout_matrix)).in(outattr_x_matrix)*x.in(outattr_x_matrix)+(C.in(outattrz_pin_matrix)-P_U.in(outattrz_pout_matrix)).in(in_arcs_from_input_per_output_attr)*z.in(in_arcs_from_input_per_output_attr);
    SPP->add(product_quality.in(J_K)<=0);
    
    Constraint<> simplex("simplex");//TODO debug transpose version
    simplex=q.in(pool_q_matrix)-1;
    SPP->add(simplex.in(L)==0);
    
    
    //    Constraint<> PQ("PQ");//TODO debug transpose version
    //    PQ=x.in(pooloutput_x_matrix)-y;
    //    SPP->add(PQ.in(Ty)==0);
    
    
    //    Constraint<> PQ1("PQ1");//TODO debug transpose version
    //    PQ1=x.in(inputpool_x_matrix)-q.in(inputpool_q_matrix)*S.in(inputpool_poolcap_matrix);
    //    SPP->add(PQ1.in(Tx)<=0);
    
    
    
    if(solv_type==gurobi){
        Constraint<> mass_balance_le("mass_balance_le");
        mass_balance_le=x.in(inputs_pools_outputs)-q.in(inpoolout_q_matrix)*y.in(inpoolout_y_matrix);
        SPP->add(mass_balance_le.in(inputs_pools_outputs)<=0);
        Constraint<> mass_balance_ge("mass_balance_ge");
        mass_balance_ge=x.in(inputs_pools_outputs)-q.in(inpoolout_q_matrix)*y.in(inpoolout_y_matrix);
        SPP->add(mass_balance_ge.in(inputs_pools_outputs)>=0);
    }
    else{
        Constraint<> mass_balance("mass_balance");
        mass_balance=x.in(inputs_pools_outputs)-q.in(inpoolout_q_matrix)*y.in(inpoolout_y_matrix);
        SPP->add(mass_balance.in(inputs_pools_outputs)==0, true);   
    }
    
    
    //    Constraint<> sumy_con("sumy_con");
    //    sumy_con=sum(x)-sumyk;
    //    SPP->add(sumy_con.in(range(0,0))>=0);
    
    
    
    
    Constraint<> obj_eq("obj_eq");
    obj_eq=objvar-x.in(inpoolout_x_matrix)*(c_tx.in(inpoolout_cip_matrix)+c_ty.in(inpoolout_cpo_matrix)).in(inpoolout_x_matrix)+product(c_tz, z);
    SPP->add(obj_eq==0);
    
    
    SPP->min(objvar);
    SPP->print();

//
//    solver<> SPP_solv(SPP,ipopt);
//    SPP_solv.run(5, 1e-6);
//
////    poolnet.get_interaction_graph();
////
//    poolnet.get_tree_decomp_bags();
////    auto bags_3d=grid.decompose_bags_3d();
//
//
////        auto g = SPP->get_interaction_graph();
////        g.get_tree_decomp_bags();
////        if(g._tree){
////            DebugOn("Interaction graph is a tree!" << endl);
////        }
////        else {
////            auto bags_3d=g.decompose_bags_3d();
////        }

    
    auto SPP_NC=build_pool_pqform(poolnet, ipopt);
    auto g = SPP_NC->get_interaction_graph();
    g.get_tree_decomp_bags();
    
       auto bags_3d=g.decompose_bags_3d();
    DebugOn("bags \n");
    for(auto bag:bags_3d){
        DebugOn(bag.second[0]->_name<<"\t"<<bag.second[1]->_name<<"\t"<<bag.second[2]->_name<<"\n");
    }
    
    g.print();
    solver<> SPP_ub_solv(SPP_NC,ipopt);
    SPP_ub_solv.run(5, 1e-6);
    double upper_bound=SPP_NC->get_obj_val();
    DebugOn("upper_bound \t"<<upper_bound<<endl);
//    solver<> SPP_lb_solv(SPP,ipopt);
//    SPP_lb_solv.run(5, 1e-6);
    //    double lower_bound=SPP->get_obj_val();
//    DebugOn("lower_bound \t"<<lower_bound<<endl);
//    double gap=(upper_bound-lower_bound)/std::abs(upper_bound)*100;
//    DebugOn("Gap before sdp \t"<<gap<<endl);
    
    auto pairs=g.get_bus_pairs();
    DebugOn("bus pairs \n");
    for(auto k: *pairs._keys){
        DebugOn(k<<endl);
    }
    var<> W("W");
    SPP->add(W.in(pairs));
    
    auto pairs1=g.get_pairs_chord(bags_3d);
    DebugOn("bus pairs chordal \n");
    for(auto k: *pairs1[0]._keys){
        DebugOn(k<<endl);
    }

    
    //auto pairs_chord=g.get_bus_pairs_chord(bags_3d);
    
//    SPP->project();

//    SPP->print();
    

    
}
