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
    fname="/Users/smitha/Desktop/Pooling_instances/Adhya3.gms";

    if(argc==2){
        fname=argv[1];
    }

    poolnet.readgrid(fname);
    SolverType solv_type = ipopt;
    
    auto SPP_NC=build_pool_pqform(poolnet, ipopt);
    auto g = SPP_NC->get_interaction_graph();
    g.pool_get_tree_decomp_bags();
    
    auto bags_3d=g.pool_decompose_bags_3d();

    
    DebugOn("bags \n");
    for(auto bag:bags_3d){
        DebugOn(bag.second[0]->_name<<"\t"<<bag.second[1]->_name<<"\t"<<bag.second[2]->_name<<"\n");
    }
    
    g.print();
    solver<> SPP_ub_solv(SPP_NC,ipopt);
    SPP_ub_solv.run(5, 1e-6);
    double upper_bound=SPP_NC->get_obj_val();
    pair<bool,double> ub = {true,upper_bound};
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
    
    
    auto res=g.get_pairs_chord(bags_3d);
    auto pairs_chordal=res[0];
    DebugOn("bus pairs chordal \n");
    for(auto k: *pairs_chordal._keys){
        DebugOn(k<<endl);
    }
    
    auto qq=res[1];
    auto yy=res[2];
    auto pairs_chordal_from = res[3];
    auto pairs_chordal_to = res[4];
    
    auto vec_param=fill_wijbounds(poolnet, res);
    auto Wij_min=vec_param[0].in(pairs_chordal);
    auto Wij_max=vec_param[1].in(pairs_chordal);
    
    auto SPP= make_shared<Model<>>("Std-Pooling-Prob-PQ");
    
    
    
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
    
    auto txty_indices=poolnet.TxplusTy_fill();
    indices TxplusTy=txty_indices[0];
    indices q_diag=txty_indices[1];
    indices y_diag=txty_indices[2];

    
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
    auto q_min=poolnet.q_min.in(Tx);
    auto q_max=poolnet.q_max.in(Tx);
    
    auto y_min=poolnet.y_min.in(Ty);
    auto y_max=poolnet.y_max.in(Ty);
    
    auto z_min=poolnet.z_min.in(Tz);
    auto z_max=poolnet.z_max.in(Tz);
    
    auto Wii_min=poolnet.Wii_min.in(TxplusTy);
    auto Wii_max=poolnet.Wii_max.in(TxplusTy);
    

    auto objvar_min = param<>("objvar_min");
    objvar_min.in(R(1));
    objvar_min = -1e6;
    auto objvar_max = param<>("objvar_max");
    objvar_max.in(R(1));
    objvar_max = upper_bound;
    var<> x("x",x_min, x_max), y("y", y_min, y_max);
    var<> q("q", q_min, q_max), z("z", z_min, z_max);
    var<> objvar("objvar",objvar_min, objvar_max);
//    auto Wij_min = param<>("Wij_min");
//    Wij_min.in(pairs_chordal);
//    Wij_min = 0;
//    auto Wij_max = param<>("Wij_max");
//    Wij_max.in(pairs_chordal);
//    Wij_max = 1e2;
//    auto Wii_min = param<>("Wii_min");
//    Wii_min.in(TxplusTy);
//    Wii_min = 0;
//    auto Wii_max = param<>("Wii_max");
//    Wii_max.in(TxplusTy);
//    Wii_max = 1e3;
//
    var<> Wij("Wij", Wij_min, Wij_max);
    var<> Wii("Wii", Wii_min, Wii_max);
    
    SPP->add(x.in(inputs_pools_outputs));
    SPP->add(q.in(Tx));
    SPP->add(z.in(Tz));
    SPP->add(y.in(Ty));
    SPP->add(objvar.in(R(1)));
    SPP->add(Wij.in(pairs_chordal));
    SPP->add(Wii.in(TxplusTy));
    
    //    SPP->add(sumyk);
    //    sumyk.set_lb(0);
    // SPP->initialize_zero();
    //    x.initialize_all(1.0);
    //    y.initialize_all(1.0);
    //    z.initialize_all(1.0);
    //    q.initialize_all(0.5);
    //    Wii.initialize_all(1.0);
    //    Wij.initialize_all(1.0);
    
    Constraint<> feed_availability("feed_availability");
    feed_availability=sum(x, input_x_matrix)+sum(z, out_arcs_to_output_per_input)-A_U;
    SPP->add(feed_availability.in(I)<=0);
    
    Constraint<> pool_capacity("pool_capacity");
    pool_capacity=sum(x, pool_x_matrix)-S;
    SPP->add(pool_capacity.in(L)<=0);
    
    Constraint<> product_demand("product_demand");
    product_demand=sum(x, output_x_matrix)+sum(z,in_arcs_from_input_per_output)-D_U;
    SPP->add(product_demand.in(J)<=0);
    
    
    Constraint<> product_quality("product_quality");    product_quality=(C.in(outattr_pin_matrix)-P_U.in(outattr_pout_matrix)).in(outattr_x_matrix)*x.in(outattr_x_matrix)+(C.in(outattrz_pin_matrix)-P_U.in(outattrz_pout_matrix)).in(in_arcs_from_input_per_output_attr)*z.in(in_arcs_from_input_per_output_attr);
    SPP->add(product_quality.in(J_K)<=0);
    
    Constraint<> simplex("simplex");//TODO debug transpose version
    simplex=q.in(pool_q_matrix)-1;
    SPP->add(simplex.in(L)==0);
    
    
        Constraint<> PQ("PQ");//TODO debug transpose version
        PQ=x.in(pooloutput_x_matrix)-y;
        SPP->add(PQ.in(Ty)==0);
    
    
        Constraint<> PQ1("PQ1");//TODO debug transpose version
        PQ1=x.in(inputpool_x_matrix)-q.in(inputpool_q_matrix)*S.in(inputpool_poolcap_matrix);
        SPP->add(PQ1.in(Tx)<=0);
    
    
    
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
    
    indices inpoolout_W_matrix = indices("inpoolout_W_matrix");
    int row_id=0;
    
    for (const string& inpoout_key:*inputs_pools_outputs._keys) {
        auto pos = nthOccurrence(inpoout_key, ",", 2);
        auto qid=inpoout_key.substr(0, pos);
        auto pos1 = nthOccurrence(inpoout_key, ",", 1);
        auto yid=inpoout_key.substr(pos1+1);
        inpoolout_W_matrix.add("q["+qid+"],y["+yid+"]");
    }
    
    
    Constraint<> x_Wij("x_Wij");
    x_Wij=x-Wij.in(inpoolout_W_matrix);
    SPP->add(x_Wij.in(inputs_pools_outputs)==0);
    
    indices qq_Wa_matrix = indices("qq_Wa_matrix");
    indices qq_Wb_matrix = indices("qq_Wb_matrix");
    

    auto q_from = Tx;
    auto q_to = Tx;
    for (const string& key:*qq._keys) {

        auto key_pos=(pairs_chordal._keys_map)->at(key);
        auto from_ref=pairs_chordal_from._ids->at(0)[key_pos];
        auto from=pairs_chordal_from._keys->at(from_ref);
        
        auto pos=from.find_first_of("[");
        auto pos1=from.find_first_of("]");
        auto a_key=from.substr(pos+1,pos1-pos-1);
        q_from.add_ref(a_key);
        
        auto to_ref=pairs_chordal_to._ids->at(0)[key_pos];
        auto to=pairs_chordal_from._keys->at(to_ref);
        
        auto pos2=to.find_first_of("[");
        auto pos3=to.find_first_of("]");
        auto b_key=to.substr(pos2+1,pos3-pos2-1);
        q_to.add_ref(b_key);

    }

    Constraint<> q_W("q_W");
    q_W=Wij.in(qq)-q.in(q_from)*q.in(q_to);
    SPP->add(q_W.in(qq)==0);
    
    q_W.print();
    
    
    

    
    //    Constraint<> sumy_con("sumy_con");
    //    sumy_con=sum(x)-sumyk;
    //    SPP->add(sumy_con.in(range(0,0))>=0);
    
    Constraint<> q2("q2");
    q2 = Wii.in(q_diag) - pow(q,2);
    SPP->add(q2.in(q_diag)==0, true, "on/off", false);
    
    Constraint<> y2("y2");
    y2 = Wii.in(y_diag) - pow(y,2);
    SPP->add(y2.in(y_diag)==0,true, "on/off", false);
    
    
    Constraint<> SOC("SOC");
    SOC = pow(Wij, 2) - Wii.in(pairs_chordal_from)*Wii.in(pairs_chordal_to);
    SPP->add(SOC.in(pairs_chordal) == 0, true);
    
//    Constraint<> SOC_nc("SOC_nc");
//    SOC_nc = pow(Wij, 2) - Wii.in(pairs_chordal_from)*Wii.in(pairs_chordal_to);
//    SPP->add(SOC_nc.in(pairs_chordal) >= 0, true);

    
    
    Constraint<> obj_eq("obj_eq");
    obj_eq = objvar - (c_tx.in(inpoolout_cip_matrix)+c_ty.in(inpoolout_cpo_matrix)).tr()*x.in(inpoolout_x_matrix).in(inpoolout_x_matrix)-product(c_tz, z);
    //    obj_eq=objvar-x.in(inpoolout_x_matrix)*(c_tx.in(inpoolout_cip_matrix)+c_ty.in(inpoolout_cpo_matrix)).in(inpoolout_x_matrix)+product(c_tz, z);
    SPP->add(obj_eq==0);
    
    
    SPP->min(objvar);
    
    auto bag_size = bags_3d.size();
    Constraint<> SDP3("SDP_3D");
    //      Constraint<> SDPD("SDPD");
    if(!poolnet._tree)
    {
        DebugOn("\nNum of bags = " << bag_size << endl);
        DebugOn("Adding 3d determinant polynomial cuts\n");
        auto Wij_ = Wij.pairs_in_bags(bags_3d, 3);
        auto Wii_ = Wii.in_bags(bags_3d, 3);
        auto nb_bags3 = Wij_[0]._indices->size();
        
        
        SDP3 = 2 * Wij_[0] * Wij_[1] * Wij_[2];
        SDP3 -= pow(Wij_[0], 2) * Wii_[2];
        SDP3 -= pow(Wij_[1], 2) * Wii_[0];
        SDP3 -= pow(Wij_[2], 2) * Wii_[1];
        SDP3 += Wii_[0] * Wii_[1] * Wii_[2];
        
       // SPP->add(SDP3.in(range(0, bag_size-1)) >= 0);
        
        DebugOn("Number of 3d determinant cuts = " << SDP3.get_nb_instances() << endl);
        
//        Constraint<> Rank_type2a("RankType2a");
//        Rank_type2a=Wij_[0]*Wij_[1]-Wii_[1]*Wij_[2];
//        SPP->add(Rank_type2a.in(range(1,nb_bags3))==0, true);
//
//        Constraint<> Rank_type2b("RankType2b");
//        Rank_type2b=Wij_[2]*(Wij_[1])-Wii_[2]*Wij_[0];
//        SPP->add(Rank_type2b.in(range(1,nb_bags3))==0, true);
//
//        Constraint<> Rank_type2c("RankType2c");
//        Rank_type2c=Wij_[2]*(Wij_[0])-Wii_[0]*Wij_[1];
//        SPP->add(Rank_type2c.in(range(1,nb_bags3))==0, true);
    }
    
    
    
    
    
    
    
    SPP->print();
    
    
    solver<> SPP_solv(SPP,solv_type);
    SPP_solv.run(5, 1e-6);
    
    auto gapl = 100*(upper_bound - SPP->get_obj_val())/std::abs(upper_bound);
    
    DebugOn("Initial Gap  = " << gapl << "%" << endl);
    double max_time = 600, precision = 1e-6;
    bool nonlin = true;
    int max_iter = 50;
    auto status = SPP->run_obbt(max_time, max_iter, ub, precision, SPP_NC, SPP, nonlin);
    status = SPP->run_obbt(max_time, max_iter, ub, precision, SPP_NC, SPP, nonlin);
    SPP->print();
//    SPP->print_solution();
    
    
    
    
    
    
    
    //    SPP->project();
    
    //    SPP->print();
    
    
    
}
