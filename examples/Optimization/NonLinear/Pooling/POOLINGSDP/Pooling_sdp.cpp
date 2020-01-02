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
    
    //string fname=string(prj_dir)+"/data_sets/Pooling/Adhya1_gms.txt";
   // string fname=argv[1];
    string fname="/Users/smitha/Desktop/Pooling_instances/Adhya1.gms";
    
    poolnet.readgrid(fname);
    SolverType solv_type = ipopt;
    
        
        auto SPP= make_shared<Model<>>("Std-Pooling-Prob-P");
        
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
        
        indices inputs_pools_outputs=poolnet.inputs_pools_outputs();
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
        
        SPP->add(x.in(inputs_pools_outputs));
        SPP->add(q.in(Tx));
        SPP->add(z.in(Tz));
        SPP->add(y.in(Ty));
        
        //    SPP->add(sumyk);
        //    sumyk.set_lb(0);
        // SPP->initialize_zero();
        x.initialize_all(1.0);
        q.initialize_all(0.5);
        
        int row_id = 0;
        indices input_x_matrix = indices("input_x_matrix");
        for (const string& input_key:*I._keys) {
            
            input_x_matrix.add_empty_row();
            for (auto &ipo:*inputs_pools_outputs._keys) {
                
                auto pos = nthOccurrence(ipo, ",", 1);
                auto input1 = ipo.substr(0,pos);
                
                if(input_key==input1){
                    input_x_matrix.add_in_row(row_id, ipo);
                }
            }
            row_id++;
        }
        
        
        
        
        Constraint<> feed_availability("feed_availability");
        feed_availability=sum(x, input_x_matrix)+sum(z, out_arcs_to_output_per_input)-A_U;
        SPP->add(feed_availability.in(I)<=0);
        
        row_id = 0;
        indices pool_x_matrix = indices("pool_x_matrix");
        for (const string& pool_key:*L._keys) {
            
            pool_x_matrix.add_empty_row();
            for (auto &ipo:*inputs_pools_outputs._keys) {
                
                auto pos = nthOccurrence(ipo, ",", 1);
                auto pos1 = nthOccurrence(ipo, ",", 2);
                auto pool1 = ipo.substr(pos+1, pos1-(pos+1));
                
                if(pool_key==pool1){
                    pool_x_matrix.add_in_row(row_id, ipo);
                }
            }
            row_id++;
        }
        
        Constraint<> pool_capacity("pool_capacity");
        pool_capacity=sum(x, pool_x_matrix)-S;
        SPP->add(pool_capacity.in(L)<=0);
        
        
        
        row_id = 0;
        indices output_x_matrix = indices("pool_x_matrix");
        for (const string& out_key:*J._keys) {
            
            output_x_matrix.add_empty_row();
            for (auto &ipo:*inputs_pools_outputs._keys) {
                
                auto pos = nthOccurrence(ipo, ",", 2);
                auto pos1 = nthOccurrence(ipo, ",", 3);
                auto out1 = ipo.substr(pos+1, pos1-(pos+1));
                
                if(out_key==out1){
                    output_x_matrix.add_in_row(row_id, ipo);
                }
            }
            row_id++;
        }
        
        
        Constraint<> product_demand("product_demand");
        product_demand=sum(x, output_x_matrix)+sum(z,in_arcs_from_input_per_output)-D_U;
        SPP->add(product_demand.in(J)<=0);
        
        
        
        row_id = 0;
        indices outattr_x_matrix = indices("outattr_x_matrix");
        for (const string& outat_key:*J_K._keys) {
            
            outattr_x_matrix.add_empty_row();
            auto out_key = outat_key.substr(0, outat_key.find_first_of(","));
            for (auto &ipo:*inputs_pools_outputs._keys) {
                
                auto pos = nthOccurrence(ipo, ",", 2);
                auto pos1 = nthOccurrence(ipo, ",", 3);
                auto out1 = ipo.substr(pos+1, pos1-(pos+1));
                
                if(out_key==out1){
                    outattr_x_matrix.add_in_row(row_id, ipo);
                }
            }
            row_id++;
        }
        
        row_id = 0;
        indices outattr_pin_matrix = indices("outattr_pin_matrix");
        for (const string& outat_key:*J_K._keys) {
            
            outattr_pin_matrix.add_empty_row();
            auto out_key = outat_key.substr(0, outat_key.find_first_of(","));
            auto at_key = outat_key.substr(outat_key.find_first_of(",")+1);
            for (auto &ipo:*inputs_pools_outputs._keys) {
                auto pos1 = nthOccurrence(ipo, ",", 1);
                auto pos2 = nthOccurrence(ipo, ",", 2);
                auto pos3 = nthOccurrence(ipo, ",", 3);
                auto out1 = ipo.substr(pos2+1, pos3-(pos2+1));
                auto in=ipo.substr(0, pos1);
                if(out_key==out1){
                    outattr_pin_matrix.add_in_row(row_id, in+","+at_key);
                }
            }
            row_id++;
        }
        
        row_id = 0;
        indices outattr_pout_matrix = indices("outattr_pout_matrix");
        for (const string& outat_key:*J_K._keys) {
            
            outattr_pout_matrix.add_empty_row();
            auto out_key = outat_key.substr(0, outat_key.find_first_of(","));
            auto at_key = outat_key.substr(outat_key.find_first_of(",")+1);
            for (auto &ipo:*inputs_pools_outputs._keys) {
                auto pos1 = nthOccurrence(ipo, ",", 1);
                auto pos2 = nthOccurrence(ipo, ",", 2);
                auto pos3 = nthOccurrence(ipo, ",", 3);
                auto out1 = ipo.substr(pos2+1, pos3-(pos2+1));
                if(out_key==out1){
                    outattr_pout_matrix.add_in_row(row_id, outat_key);
                }
            }
            row_id++;
        }
        
        row_id = 0;
        indices outattrz_pin_matrix = indices("outattrz_pin_matrix");
        for (const string& outat_key:*J_K._keys) {
            
            outattrz_pin_matrix.add_empty_row();
            auto out_key = outat_key.substr(0, outat_key.find_first_of(","));
            auto at_key = outat_key.substr(outat_key.find_first_of(",")+1);
            for (auto &ipo:*Tz._keys) {
                auto pos1 = nthOccurrence(ipo, ",", 1);
                auto out1 = ipo.substr(pos1+1);
                auto in=ipo.substr(0, pos1);
                if(out_key==out1){
                    outattrz_pin_matrix.add_in_row(row_id, in+","+at_key);
                }
            }
            row_id++;
        }
        
        row_id = 0;
        indices outattrz_pout_matrix = indices("outattrz_pout_matrix");
        for (const string& outat_key:*J_K._keys) {
            
            outattrz_pout_matrix.add_empty_row();
            auto out_key = outat_key.substr(0, outat_key.find_first_of(","));
            auto at_key = outat_key.substr(outat_key.find_first_of(",")+1);
            for (auto &ipo:*Tz._keys) {
                auto pos1 = nthOccurrence(ipo, ",", 1);
                auto out1=ipo.substr(pos1+1);
                if(out_key==out1){
                    outattrz_pout_matrix.add_in_row(row_id, outat_key);
                }
            }
            row_id++;
        }
        
        Constraint<> product_quality("product_quality");//TODO debug transpose version
        product_quality=(C.in(outattr_pin_matrix)-P_U.in(outattr_pout_matrix)).in(outattr_x_matrix)*x.in(outattr_x_matrix)+(C.in(outattrz_pin_matrix)-P_U.in(outattrz_pout_matrix)).in(in_arcs_from_input_per_output_attr)*z.in(in_arcs_from_input_per_output_attr);
        SPP->add(product_quality.in(J_K)<=0);
        
        
        row_id = 0;
        indices pool_q_matrix = indices("pool_q_matrix");
        for (const string& pool_key:*L._keys) {
            
            pool_q_matrix.add_empty_row();
            for (auto &ip:*Tx._keys) {
                
                auto pos = nthOccurrence(ip, ",", 1);
                auto pool1 = ip.substr(pos+1);
                
                if(pool_key==pool1){
                    pool_q_matrix.add_in_row(row_id, ip);
                }
            }
            row_id++;
        }
        
        Constraint<> simplex("simplex");//TODO debug transpose version
        simplex=q.in(pool_q_matrix)-1;
        SPP->add(simplex.in(L)==0);
        
        
        row_id = 0;
        indices pooloutput_x_matrix = indices("pooloutput_x_matrix");
        for (const string& poolout_key:*Ty._keys) {
            
            pooloutput_x_matrix.add_empty_row();
            for (auto &ipo:*inputs_pools_outputs._keys) {
                
                auto pos = nthOccurrence(ipo, ",", 1);
                auto poolout1 = ipo.substr(pos+1);
                
                if(poolout_key==poolout1){
                    pooloutput_x_matrix.add_in_row(row_id, ipo);
                }
            }
            row_id++;
        }
        
        
        //    Constraint<> PQ("PQ");//TODO debug transpose version
        //    PQ=x.in(pooloutput_x_matrix)-y;
        //    SPP->add(PQ.in(Ty)==0);
        
        
        row_id = 0;
        indices inputpool_x_matrix = indices("inputpool_x_matrix");
        for (const string& inputpool_key:*Tx._keys) {
            
            inputpool_x_matrix.add_empty_row();
            for (auto &ipo:*inputs_pools_outputs._keys) {
                
                auto pos = nthOccurrence(ipo, ",", 2);
                auto inputpool1 = ipo.substr(0, pos);
                
                if(inputpool_key==inputpool1){
                    inputpool_x_matrix.add_in_row(row_id, ipo);
                }
            }
            row_id++;
        }
        row_id = 0;
        indices inputpool_q_matrix = indices("inputpool_q_matrix");
        indices inputpool_poolcap_matrix = indices("inputpool_poolcap_matrix");
        for (const string& inputpool_key:*Tx._keys) {
            
            inputpool_q_matrix.add_empty_row();
            inputpool_poolcap_matrix.add_empty_row();
            inputpool_q_matrix.add_in_row(row_id, inputpool_key);
            auto pos=nthOccurrence(inputpool_key, ",", 1);
            auto pool_key=inputpool_key.substr(pos+1);
            inputpool_poolcap_matrix.add_in_row(row_id, pool_key);
            row_id++;
        }
        
        
        
        
        
        
        //    Constraint<> PQ1("PQ1");//TODO debug transpose version
        //    PQ1=x.in(inputpool_x_matrix)-q.in(inputpool_q_matrix)*S.in(inputpool_poolcap_matrix);
        //    SPP->add(PQ1.in(Tx)<=0);
        
        
        row_id = 0;
        indices inpoolout_y_matrix = indices("inpoolout_y_matrix");
        indices inpoolout_q_matrix = indices("inpoolout_q_matrix");
        for (const string& inpoout_key:*inputs_pools_outputs._keys) {
            auto pos = nthOccurrence(inpoout_key, ",", 1);
            auto poout=inpoout_key.substr(pos+1);
            auto pos1 = nthOccurrence(inpoout_key, ",", 2);
            auto inpo=inpoout_key.substr(0,pos1);
            inpoolout_y_matrix.add_in_row(row_id, poout);
            inpoolout_q_matrix.add_in_row(row_id, inpo);
            row_id++;
        }
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
//            SPP->add(mass_balance.in(inputs_pools_outputs)==0, true);
            
        }
        
        Constraint<> sumy_con("sumy_con");
        sumy_con=sum(x)-sumyk;
        SPP->add(sumy_con.in(range(0,0))>=0);
        
        row_id = 0;
        indices inpoolout_x_matrix = indices("inpoolout_x_matrix");
        indices inpoolout_cip_matrix = indices("inpoolout_cip_matrix");
        indices inpoolout_cpo_matrix = indices("inpoolout_cpo_matrix");
        for (const string& inpoout_key:*inputs_pools_outputs._keys) {
            inpoolout_x_matrix.add_empty_row();
            inpoolout_x_matrix.add_in_row(row_id, inpoout_key);
            auto pos = nthOccurrence(inpoout_key, ",", 1);
            auto poout=inpoout_key.substr(pos+1);
            auto pos1 = nthOccurrence(inpoout_key, ",", 2);
            auto inpo=inpoout_key.substr(0,pos1);
            inpoolout_cip_matrix.add_in_row(row_id,inpo);
            inpoolout_cpo_matrix.add_in_row(row_id,poout);
            
            
        }
        row_id++;
        
        auto obj=x.in(inpoolout_x_matrix)*(c_tx.in(inpoolout_cip_matrix)+c_ty.in(inpoolout_cpo_matrix)).in(inpoolout_x_matrix)+product(c_tz, z);
        
        SPP->min(obj);
        SPP->print();
        
        
    
    auto start=get_wall_time();
    solver<> SPP_solv(SPP,solv_type);
    SPP_solv.run(5, 1e-6);
    auto end=get_wall_time();
    auto comp_time=end-start;

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
