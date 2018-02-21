//  Decompose.cpp
//
//  Gravity
//
//  Created by Guanglei Wang on 6/9/17.
//
//
#include "Partition.hpp"
bool node_id_compare(const Node* n1, const Node* n2) {
    return n1->_id < n2->_id;
}

bool arc_id_compare(const Arc* n1, const Arc* n2) {
    return n1->_id < n2->_id;
}

bool bus_pair_compare(const index_pair* n1, const index_pair* n2) {
    return n1->_name < n2->_name;
}

struct net_param {
    param<Real> c0, c1, c2; /**< Generation costs */
    param<Real> tan_th_min, tan_th_max;
    param<Real> g_ff, g_ft, g_tt, g_tf, b_ff, b_ft, b_tf, b_tt;
    param<Real> S_max;
};
double getdual_relax(PowerNet& grid, const unsigned T, const Partition& P, const unsigned nbparts,
                     const param<Real>& rate_ramp, const param<Real>& rate_switch,
                     const param<Real>& min_up, const param<Real>& min_down,
                     const param<Real>& cost_up, const param<Real>& cost_down,
                     vector<var<Real>>& Pg, vector<var<Real>>& Qg,
                     vector<var<bool>>& Start_up, vector<var<bool>>& Shut_down, vector<var<bool>>& On_off,
                     vector<var<Real>>& Xii, vector<var<Real>>& R_Xij, vector<var<Real>>& Im_Xij,
                     param<Real>& R_lambda, param<Real>& Im_lambda, param<Real>& lambda)
{
    // Grid Parameters
    Model ACUC("ACUC Model");
    for (int c = 0; c < nbparts; c++) {
        ACUC.add_var(Pg[c].in(P.bag_gens[c], T));
        ACUC.add_var(Qg[c].in(P.bag_gens[c], T));
        //ACUC.add_var(Start_up[c].in(P.bag_gens[c], T));
        //ACUC.add_var(Shut_down[c].in(P.bag_gens[c], T));
        //ACUC.add_var(On_off[c].in(P.bag_gens[c], T));
        ACUC.add_var(Xii[c].in(P.bag_bus_union_out[c], T));
        ACUC.add_var(R_Xij[c].in(P.bag_bus_pairs_union[c], T));
        cout << "R_Xij dim: " << R_Xij[c].in(P.bag_bus_pairs_union[c], T).get_dim() << endl;
        ACUC.add_var(Im_Xij[c].in(P.bag_bus_pairs_union[c], T));
        cout << "Im_Xij dim: " << Im_Xij[c].in(P.bag_bus_pairs_union[c], T).get_dim() << endl;
    }
    //power flow vars are treated as auxiliary vars.
    vector<var<Real>> Pf_from;
    vector<var<Real>> Qf_from;
    vector<var<Real>> Pf_to;
    vector<var<Real>> Qf_to;
    for (int c = 0; c < nbparts; c++) {
        if (P.bag_arcs_union[c].size() > 0) {
            var<Real> bag_Pf_from("Pf_from" +to_string(c), grid.S_max.in(P.bag_arcs_union_out[c], T));
            var<Real> bag_Qf_from("Qf_from"+to_string(c), grid.S_max.in(P.bag_arcs_union_out[c], T));
            var<Real> bag_Pf_to("Pf_to"+to_string(c), grid.S_max.in(P.bag_arcs_union_out[c], T));
            var<Real> bag_Qf_to("Qf_to"+to_string(c), grid.S_max.in(P.bag_arcs_union_out[c], T));
            Pf_from.push_back(bag_Pf_from);
            Pf_to.push_back(bag_Pf_to);
            Qf_from.push_back(bag_Qf_from);
            Qf_to.push_back(bag_Qf_to);
            ACUC.add_var(Pf_from[c].in(P.bag_arcs_union_out[c], T));
            ACUC.add_var(Qf_from[c].in(P.bag_arcs_union_out[c], T));
            ACUC.add_var(Pf_to[c].in(P.bag_arcs_union_in[c], T));
            ACUC.add_var(Qf_to[c].in(P.bag_arcs_union_in[c], T));
        }
        else {
            var<Real> empty("empty");
            empty.set_size(0);
            Pf_from.push_back(empty);
            Pf_to.push_back(empty);
            Qf_from.push_back(empty);
            Qf_to.push_back(empty);
        }
    }
    /* Construct the objective function*/
    func_ obj;
    for (int c = 0; c < nbparts; c++) {
        for (auto g: P.bag_gens[c]) {
            if (g->_active) {
                string name = g->_name+",0";
                obj += grid.c1(name)*Pg[c](name)+ grid.c2(name)*Pg[c](name)*Pg[c](name) +grid.c0(name)*On_off[c](name);
            }
            for (int t = 1; t < T; t++) {
                if (g->_active) {
                    string name = g->_name + ","+ to_string(t);
                    obj += grid.c1(name)*Pg[c](name)+ grid.c2(name)*Pg[c](name)*Pg[c](name) + grid.c0(name)*On_off[c](name);
                    obj += cost_up.getvalue()*Start_up[c](name)+ cost_down.getvalue()*Shut_down[c](name);
                }
            }
        }
    }
    ACUC.set_objective(min(obj));
    for (int c = 0; c < nbparts; c++) {
        if (P.bag_bus_pairs_union_directed[c].size() > 0){
            Constraint SOC("SOC_" + to_string(c));
            SOC =  power(R_Xij[c], 2)+ power(Im_Xij[c], 2) - Xii[c].from()*Xii[c].to() ;
            ACUC.add_constraint(SOC.in(P.bag_bus_pairs_union_directed[c], T) <= 0);
        }
    }
    //KCL
 for (int c = 0; c < nbparts; c++) {
        Constraint KCL_P("KCL_P"+ to_string(c));
        Constraint KCL_Q("KCL_Q"+ to_string(c));
        KCL_P = sum(Pf_from[c].out_arcs()) + sum(Pf_to[c].in_arcs())+ grid.pl -sum(Pg[c].in_gens())+ grid.gs*Xii[c];
        ACUC.add_constraint(KCL_P.in(P.bag_bus[c], T) == 0);
        KCL_P.in(P.bag_bus[c], T).print();
        cout << "# constraints: " << KCL_P.in(P.bag_bus[c], T).get_nb_instances()<< endl;
        KCL_Q  = sum(Qf_from[c].out_arcs()) + sum(Qf_to[c].in_arcs())+ grid.ql -sum(Qg[c].in_gens())-grid.bs*Xii[c];
        ACUC.add_constraint(KCL_Q.in(P.bag_bus[c], T) == 0);
        cout << "# constraints: " << KCL_Q.in(P.bag_bus[c], T).get_nb_instances()<< endl;

        
        Constraint Flow_P_From("Flow_P_From_" + to_string(c));
        Flow_P_From = Pf_from[c]- (grid.g_ff*Xii[c].from()+ grid.g_ft*R_Xij[c].in_pairs() + grid.b_ft*Im_Xij[c].in_pairs());
        ACUC.add_constraint(Flow_P_From.in(P.bag_arcs_union_out[c], T) == 0);
        cout << "# constraints: " << Flow_P_From.in(P.bag_arcs_union_out[c], T).get_nb_instances()<< endl;

        Constraint Flow_P_To("Flow_P_To" + to_string(c));
        Flow_P_To = Pf_to[c] - (grid.g_tt*Xii[c].to() + grid.g_tf*R_Xij[c].in_pairs() - grid.b_tf*Im_Xij[c].in_pairs());
        ACUC.add_constraint(Flow_P_To.in(P.bag_arcs_union_in[c], T) == 0);
        cout << "# constraints: " << Flow_P_To.in(P.bag_arcs_union_in[c], T).get_nb_instances()<< endl;

        Constraint Flow_Q_From("Flow_Q_From" + to_string(c));
        Flow_Q_From = Qf_from[c] - (grid.g_ft*Im_Xij[c].in_pairs() - grid.b_ff*Xii[c].from() - grid.b_ft*R_Xij[c].in_pairs());
        ACUC.add_constraint(Flow_Q_From.in(P.bag_arcs_union_out[c], T) == 0);

        Constraint Flow_Q_To("Flow_Q_To" + to_string(c));
        Flow_Q_To = Qf_to[c] + (grid.b_tt*Xii[c].to() + grid.b_tf*R_Xij[c].in_pairs() + grid.g_tf*Im_Xij[c].in_pairs());
        ACUC.add_constraint(Flow_Q_To.in(P.bag_arcs_union_in[c], T) == 0);

////        Constraint Thermal_Limit_from("Thermal_Limit_from" + to_string(c));
////        Thermal_Limit_from = power(Pf_from[c], 2) + power(Qf_from[c], 2);
////        Thermal_Limit_from <= power(grid.S_max,2);
////        ACUC.add_constraint(Thermal_Limit_from.in(P.bag_arcs_union_out[c], T));
////
////        Constraint Thermal_Limit_to("Thermal_Limit_to" + to_string(c));
////        Thermal_Limit_to = power(Pf_to[c], 2) + power(Qf_to[c], 2);
////        Thermal_Limit_to <= power(grid.S_max,2);
////        ACUC.add_constraint(Thermal_Limit_to.in(P.bag_arcs_union_in[c], T));
    }
//    ///* Phase Angle Bounds constraints */
    //for (int c = 0; c < nbparts; c++) {
    //    if (P.bag_bus_pairs_union_directed[c].size() > 0){
    //        Constraint PAD_UB("PAD_UB_"+to_string(c));
    //        PAD_UB = Im_Xij[c]- grid.tan_th_max*R_Xij[c];
    //        ACUC.add_constraint(PAD_UB.in(P.bag_bus_pairs_union_directed[c], T) <= 0);

    //        Constraint PAD_LB("PAD_LB_"+to_string(c));
    //        PAD_LB = Im_Xij[c]- grid.tan_th_min*R_Xij[c];
    //        ACUC.add_constraint(PAD_LB.in(P.bag_bus_pairs_union_directed[c], T) >= 0);
    //    }
    //}
    // COMMITMENT CONSTRAINTS
    // Inter-temporal constraints 3a, 3d
//    for (int c = 0; c < nbparts; c++) {
//        if (P.bag_gens[c].size() > 0) {
//            for (int t = 1; t < T; t++) {
//                Constraint MC1("MC1_"+ to_string(c)+ ","+to_string(t));
//                Constraint MC2("MC2_"+ to_string(c)+ ","+to_string(t));
//                MC1 = On_off[c].in_at(P.bag_gens[c], t)- On_off[c].in_at(P.bag_gens[c], t-1)-  Start_up[c].in_at(P.bag_gens[c], t);
//                MC2 = On_off[c].in_at(P.bag_gens[c], t-1) - On_off[c].in_at(P.bag_gens[c], t) - Shut_down[c].in_at(P.bag_gens[c], t);
//                ACUC.add_constraint(MC1 <= 0);
//                ACUC.add_constraint(MC2 <= 0);
//            }
//        }
//    }
//    // Min-up constraints  4a
//    for (int c = 0; c < nbparts; c++) {
//        if (P.bag_gens[c].size() > 0) {
//            for (int t = 1; t < T; t++) {
//                Constraint Min_up1("Min_up1_"+ to_string(c) + "_"+to_string(t));
//                Min_up1 = On_off[c].in_at(P.bag_gens[c], t) - On_off[c].in_at(P.bag_gens[c], t-1) - Start_up[c].in_at(P.bag_gens[c], t) + Shut_down[c].in_at(P.bag_gens[c], t);
//                ACUC.add_constraint(Min_up1 == 0);
//            }
//        }
//    }
//    // 4b
//    for (int c = 0; c < nbparts; c++) {
//        if (P.bag_gens[c].size() > 0) {
//            for (int t = min_up.getvalue(); t < T; t++) {
//                Constraint Min_Up("Min_Up_constraint" + to_string(c) + "_"+ to_string(t));
//                for (int l = t-min_up.getvalue()+1; l < t+1; l++) {
//                    Min_Up   += Start_up[c].in_at(P.bag_gens[c], l);
//                }
//                Min_Up -= On_off[c].in_at(P.bag_gens[c], t);
//                ACUC.add_constraint(Min_Up <= 0);
//            }
//        }
//    }
//    // 4c
//    for (int c = 0; c < nbparts; c++) {
//        if (P.bag_gens[c].size() > 0) {
//            for (int t = min_down.getvalue(); t < T; t++) {
//                Constraint Min_Down("Min_Down_constraint_" + to_string(c) + "_"+ to_string(t));
//                for (int l = t-min_down.getvalue()+1; l < t +1; l++) {
//                    Min_Down   += Shut_down[c].in_at(P.bag_gens[c], l);
//                }
//                Min_Down -= 1 - On_off[c].in_at(P.bag_gens[c], t);
//                ACUC.add_constraint(Min_Down <= 0);
//            }
//        }
//    }
//
//   ////Ramp Rate
//    for (int c = 0; c < nbparts; c++) {
//        if (P.bag_gens[c].size() > 0) {
//            Constraint Production_P_LB("Production_P_LB_"+ to_string(c));
//            Constraint Production_P_UB("Production_P_UB_"+ to_string(c));
//            Constraint Production_Q_LB("Production_Q_LB_"+ to_string(c));
//            Constraint Production_Q_UB("Production_Q_UB_"+ to_string(c));
//            // 5A
//            Production_P_UB = Pg[c].in(P.bag_gens[c], T) - grid.pg_max.in(P.bag_gens[c], T)*On_off[c].in(P.bag_gens[c],T);
//            Production_P_LB = Pg[c].in(P.bag_gens[c], T) - grid.pg_min.in(P.bag_gens[c], T)*On_off[c].in(P.bag_gens[c],T);
//            ACUC.add_constraint(Production_P_UB <=0);
//            ACUC.add_constraint(Production_P_LB >= 0);
//
//            Production_Q_UB = Qg[c].in(P.bag_gens[c], T) - grid.qg_max.in(P.bag_gens[c], T)*On_off[c].in(P.bag_gens[c],T);
//            Production_Q_LB = Qg[c].in(P.bag_gens[c], T) - grid.qg_min.in(P.bag_gens[c], T)*On_off[c].in(P.bag_gens[c],T);
//            ACUC.add_constraint(Production_Q_UB <= 0);
//            ACUC.add_constraint(Production_Q_LB >= 0);
//        }
//    }
////    // 5C
//    for (int c = 0; c < nbparts; c++) {
//        if (P.bag_gens[c].size() > 0) {
//            for (int t = 1; t < T; t++) {
//                Constraint Ramp_up("Ramp_up_constraint_" +  to_string(c) + "_" + to_string(t));
//                Constraint Ramp_down("Ramp_down_constraint_" + to_string(c) + "_" +to_string(t));
//                Ramp_up =  Pg[c].in_at(P.bag_gens[c], t);
//                Ramp_up -= Pg[c].in_at(P.bag_gens[c], t-1);
//                Ramp_up -= rate_ramp*On_off[c].in_at(P.bag_gens[c], t-1);
//                Ramp_up -= rate_switch*(1 - On_off[c].in_at(P.bag_gens[c], t));
//
//                Ramp_down =  Pg[c].in_at(P.bag_gens[c], t-1);
//                Ramp_down -= Pg[c].in_at(P.bag_gens[c], t);
//                Ramp_down -= rate_ramp*On_off[c].in_at(P.bag_gens[c], t);
//                Ramp_down -= rate_switch*(1 - On_off[c].in_at(P.bag_gens[c], t-1));
//
//                ACUC.add_constraint(Ramp_up <= 0);
//                ACUC.add_constraint(Ramp_down <= 0);
//            }
//        }
//    }
////    // set the initial state of generators.
//    for (int c = 0; c < nbparts; c++) {
//        for(auto g: P.bag_gens[c]) {
//            Constraint gen_initial("gen_" + g->_name + to_string(c) + ",0");
//            gen_initial +=  On_off[c]( g->_name + ",0");
//            ACUC.add_constraint(gen_initial == 1);
//        }
//    }
    // Linking constraints
    for (const auto& a: P.G_part.arcs) {
//        Constraint Link_R("link_R_"+a->_name);
//        Link_R = R_Xij[a->_src->_id].in_pairs() - R_Xij[a->_dest->_id].in_pairs();
//        ACUC.add_constraint(Link_R.in(a->_intersection_clique, T) ==0);

//        Constraint Link_Im("link_Im_"+a->_name);
//        Link_Im = Im_Xij[a->_src->_id].in_pairs() - Im_Xij[a->_dest->_id].in_pairs();
//        ACUC.add_constraint(Link_Im.in(a->_intersection_clique, T)==0);
          for (const auto& pair: a->_intersection_clique){
              Constraint Link_Im("link_Im_"+a->_name + "," + pair->_name);
              Constraint Link_R("link_R_"+a->_name + "," + pair->_name);
              Link_Im = Im_Xij[a->_src->_id](pair->_name+","+to_string(0)) -Im_Xij[a->_dest->_id](pair->_name+","+to_string(0));
              Link_R = R_Xij[a->_src->_id](pair->_name+","+to_string(0)) -R_Xij[a->_dest->_id](pair->_name+","+to_string(0));
              ACUC.add_constraint(Link_Im ==0);
              ACUC.add_constraint(Link_R ==0);
              Constraint Link_Xii("link_Xii_" + a->_name + "," + pair->_name);
              Link_Xii = Xii[a->_src->_id](pair->_dest->_name+","+to_string(0)) - Xii[a->_dest->_id](pair->_dest->_name+","+to_string(0));
              ACUC.add_constraint(Link_Xii==0);
          }
//        Constraint Link_Xii("link_Xii_" + a->_name);
//        Link_Xii = Xii[a->_src->_id].to() - Xii[a->_dest->_id].to();
//        ACUC.add_constraint(Link_Xii.in(a->_intersection_clique, T)==0);
    }
    //for (const auto& a: P.G_part.arcs) {
    //   // Constraint Link_R("link_R_"+a->_name);
    //   // Link_R = R_Xij[a->_src->_id].in_pairs(a->_intersection_clique, 1) - R_Xij[a->_dest->_id].in_pairs(a->_intersection_clique, 1);
    //    //ACUC.add_constraint(Link_R ==0);
    //    //Link_Im = Im_Xij[a->_src->_id].in_pairs()-Im_Xij[a->_dest->_id].in_pairs();
    //    //Link_Im.print(true);
    //    for (const auto& pair: a->_intersection_clique){
    //        cout << "pair_name: " << pair->_name << endl;
    //        Constraint Link_Im("link_Im_"+a->_name + "," + pair->_name);
    //        Constraint Link_R("link_R_"+a->_name + "," + pair->_name);
    //        Link_Im = Im_Xij[a->_src->_id](pair->_name+","+to_string(0)) -Im_Xij[a->_dest->_id](pair->_name+","+to_string(0));
    //        Link_R = R_Xij[a->_src->_id](pair->_name+","+to_string(0)) -R_Xij[a->_dest->_id](pair->_name+","+to_string(0));
    //        Link_Im.print(true);
    //        ACUC.add_constraint(Link_Im ==0);
    //        ACUC.add_constraint(Link_R ==0);
    //    }
    //    
    //    Constraint Link_Xii("link_Xii_" + a->_name);
    //    Link_Xii = Xii[a->_src->_id].to() - Xii[a->_dest->_id].to();
    //    ACUC.add_constraint(Link_Xii.in(a->_intersection_clique, 1)== 0);
    //}
    /* Solver selection */
    solver cpx_acuc(ACUC, ipopt);
    auto solver_time_start = get_wall_time();
    auto cpu_time_start = get_wall_time();

    bool relax =true;
    int output = 1;
    double tol = 1e-6;
    cpx_acuc.run(output, relax, tol);
    cout << "the continuous relaxation bound is: " << ACUC._obj_val << endl;
    // return dual vars.
//    for (const auto& a: P.G_part.arcs) {
//        auto cons = *ACUC.get_constraint("link_R_"+a->_name);
//        for (auto i = 0; i <cons._dual.size(); i++){
//            cout << "dual vars: " << cons._dual[i] << endl;
//        }
//    }
    return ACUC._obj_val;
}
/** INITIALISE SUBPROBLEM MODEL */
double subproblem(PowerNet& grid,  unsigned T, const Partition& P, unsigned c,
                  const param<Real>& rate_ramp,const param<Real>& rate_switch,
                  const param<Real>& min_up, const param<Real>& min_down,
                  const param<Real>& cost_up, const param<Real>& cost_down,
                  var<Real>& Pg, var<Real>& Qg, var<bool>& Start_up, var<bool>& Shut_down, var<bool>& On_off,
                  var<Real>& Xii, var<Real>& R_Xij, var<Real>& Im_Xij,
                  param<Real>& R_lambda_sep, param<Real>& Im_lambda_sep,param<Real>& lambda_sep,
                  param<Real>& Xii_log, param<Real>& R_Xij_log, param<Real>& Im_Xij_log)
{
    Model Subr("Subr");
    Subr.add_var(Pg.in(P.bag_gens[c], T));
    Subr.add_var(Qg.in(P.bag_gens[c], T));
    Subr.add_var(On_off.in(P.bag_gens[c], T));
    Subr.add_var(Start_up.in(P.bag_gens[c], T));
    Subr.add_var(Shut_down.in(P.bag_gens[c], T));
    Subr.add_var(Xii.in(P.bag_bus[c], T));
    Subr.add_var(R_Xij.in(P.bag_bus_pairs_union[c], T));
    Subr.add_var(Im_Xij.in(P.bag_bus_pairs_union[c], T));
    
    //power flow
    var<Real> Pf_from("Pf_from", grid.S_max);
    var<Real> Qf_from("Qf_from", grid.S_max);
    var<Real> Pf_to("Pf_to", grid.S_max);
    var<Real> Qf_to("Qf_to", grid.S_max);
    if (P.bag_arcs_union[c].size() > 0) {
        Subr.add_var(Pf_from.in(P.bag_arcs_union[c], T));
        Subr.add_var(Qf_from.in(P.bag_arcs_union[c], T));
        Subr.add_var(Pf_to.in(P.bag_arcs_union[c], T));
        Subr.add_var(Qf_to.in(P.bag_arcs_union[c], T));
    }
    /* Construct the objective function*/
    func_ obj;
    for (auto g:P.bag_gens[c]) {
        if (g->_active) {
            string name = g->_name + ",0";
            obj += grid.c1(name)*Pg(name)+ grid.c2(name)*Pg(name)*Pg(name) + grid.c0(name)*On_off(name);

            for (int t = 1; t < T; t++) {
                string name2 = g->_name + ","+ to_string(t);
                obj += grid.c1(name)*Pg(name) + grid.c2(name)*Pg(name)*Pg(name) + grid.c0(name)*On_off(name);
                obj += cost_up.getvalue()*Start_up(name)+ cost_down.getvalue()*Shut_down(name);
            }
        }
    }
    //RELAXED TERMS
    auto bag = P.G_part.get_node(to_string(c));
    for (const auto& a: bag->get_out()) {
        for (int t = 1; t < T; t++) {
            for (auto pair: a->_intersection_clique) {
                string name = pair->_name + ","+ to_string(t);
                string dest = pair->_dest->_name + ","+ to_string(t);
                obj += R_lambda_sep(name)*R_Xij(name);
                obj += Im_lambda_sep(name)*Im_Xij(name);
                obj += lambda_sep(name)*Xii(dest);
            }
        }
    }
    for (const auto& a: bag->get_in()) {
        for (int t = 1; t < T; t++) {
            for (auto pair: a->_intersection_clique) {
                string name = pair->_name + ","+ to_string(t);
                string dest = pair->_dest->_name + ","+ to_string(t);
                obj -= R_lambda_sep(name)*R_Xij(name);
                obj -= Im_lambda_sep(name)*Im_Xij(name);
                obj -= lambda_sep(name)*Xii(dest);
            }
        }
    }
    Subr.min(obj);
    if (P.bag_bus_pairs_union_directed[c].size() > 0) {
        Constraint SOC("SOC_" + to_string(c));
        SOC =  power(R_Xij.in(P.bag_bus_pairs_union_directed[c], T), 2)+ power(Im_Xij.in(P.bag_bus_pairs_union_directed[c], T), 2)
               - Xii.from(P.bag_bus_pairs_union_directed[c], T)*Xii.to(P.bag_bus_pairs_union_directed[c], T) ;
        Subr.add_constraint(SOC <= 0);
    }
    Constraint KCL_P("KCL_P");
    Constraint KCL_Q("KCL_Q");
    KCL_P = sum(Pf_from.out_arcs()) + sum(Pf_to.in_arcs())+grid.pl-sum(Pg.in_gens()) + grid.gs*Xii;
    Subr.add_constraint(KCL_P.in(P.bag_bus[c], T) == 0);
    KCL_Q = sum(Qf_from.out_arcs()) + sum(Qf_to.in_arcs()) + grid.ql - sum(Qg.in_gens()) - grid.bs*Xii;
    Subr.add_constraint(KCL_Q.in(P.bag_bus[c], T) == 0);

    Constraint Flow_P_From("Flow_P_From");
    Flow_P_From = Pf_from -(grid.g_ff*Xii.from() + grid.g_ft*R_Xij.in_pairs() + grid.b_ft*Im_Xij.in_pairs());
    Subr.add_constraint(Flow_P_From.in(P.bag_arcs_union_out[c], T) == 0);

    Constraint Flow_P_To("Flow_P_To");
    Flow_P_To = Pf_to - (grid.g_tt*Xii.to() + grid.g_tf*R_Xij.in_pairs() - grid.b_tf*Im_Xij.in_pairs());
    Subr.add_constraint(Flow_P_To.in(P.bag_arcs_union[c], T) == 0);

    Constraint Flow_Q_From("Flow_Q_From");
    Flow_Q_From = Qf_from - (grid.g_ft*Im_Xij.in_pairs() - grid.b_ff*Xii.from() - grid.b_ft*R_Xij.in_pairs());
    Subr.add_constraint(Flow_Q_From.in(P.bag_arcs_union[c], T) == 0);

    Constraint Flow_Q_To("Flow_Q_To");
    Flow_Q_To = Qf_to + (grid.b_tt*Xii.to() + grid.b_tf*R_Xij.in_pairs() + grid.g_tf*Im_Xij.in_pairs());
    Subr.add_constraint(Flow_Q_To.in(P.bag_arcs_union[c], T) == 0);

    Constraint Thermal_Limit_from("Thermal_Limit_from");
    Thermal_Limit_from = power(Pf_from, 2) + power(Qf_from, 2);
    Thermal_Limit_from <= power(grid.S_max,2);
    Subr.add_constraint(Thermal_Limit_from.in(P.bag_arcs_union[c], T));

    Constraint Thermal_Limit_to("Thermal_Limit_to");
    Thermal_Limit_to = power(Pf_to, 2) + power(Qf_to, 2);
    Thermal_Limit_to <= power(grid.S_max, 2);
    Subr.add_constraint(Thermal_Limit_to.in(P.bag_arcs_union[c], T));
    /* Phase Angle Bounds constraints */
    vector<index_pair*> pairs = P.bag_bus_pairs_union_directed[c];
    Constraint PAD_UB("PAD_UB_"+to_string(c));
    PAD_UB = Im_Xij.in(pairs, T)- grid.tan_th_max.in(pairs, T)*R_Xij.in(pairs, T);
    Subr.add_constraint(PAD_UB <= 0);

    Constraint PAD_LB("PAD_LB_"+to_string(c));
    PAD_LB = Im_Xij.in(pairs, T)- grid.tan_th_min.in(pairs, T)*R_Xij.in(pairs, T);
    Subr.add_constraint(PAD_LB >= 0);

    // COMMITMENT CONSTRAINTS
    // Inter-temporal constraints 3a, 3d
    if (P.bag_gens[c].size() > 0) {
        for (int t = 1; t < T; t++) {
            Constraint MC1("MC1_"+to_string(t));
            Constraint MC2("MC2_"+to_string(t));
            MC1 = On_off.in_at(P.bag_gens[c], t)- On_off.in_at(P.bag_gens[c], t-1)-  Start_up.in_at(P.bag_gens[c], t);
            MC2 = On_off.in_at(P.bag_gens[c], t-1) - On_off.in_at(P.bag_gens[c], t) - Shut_down.in_at(P.bag_gens[c], t);
            Subr.add_constraint(MC1 <= 0);
            Subr.add_constraint(MC2 <= 0);
        }
    }
    // Min-up constraints  4a
    if (P.bag_gens[c].size() > 0) {
        for (int t = 1; t < T; t++) {
            Constraint Min_up1("Min_up1_"+ to_string(c) + "_"+to_string(t));
            Min_up1 = On_off.in_at(P.bag_gens[c], t) - On_off.in_at(P.bag_gens[c], t-1) - Start_up.in_at(P.bag_gens[c], t) + Shut_down.in_at(P.bag_gens[c], t);
            Subr.add_constraint(Min_up1 == 0);
        }
    }
    // 4b
    if (P.bag_gens[c].size() > 0) {
        for (int t = min_up.getvalue(); t < T; t++) {
            Constraint Min_Up("Min_Up_constraint" + to_string(c) + "_"+ to_string(t));
            for (int l = t-min_up.getvalue()+1; l < t+1; l++) {
                Min_Up   += Start_up.in_at(P.bag_gens[c], l);
            }
            Min_Up -= On_off.in_at(P.bag_gens[c], t);
            Subr.add_constraint(Min_Up <= 0);
        }
    }
    // 4c
    if (P.bag_gens[c].size() > 0) {
        for (int t = min_down.getvalue(); t < T; t++) {
            Constraint Min_Down("Min_Down_constraint_"+ to_string(t));
            for (int l = t-min_down.getvalue()+1; l < t +1; l++) {
                Min_Down   += Shut_down.in_at(P.bag_gens[c], l);
            }
            Min_Down -= 1 - On_off.in_at(P.bag_gens[c], t);
            Subr.add_constraint(Min_Down <= 0);
        }
    }

    ////Ramp Rate
    if (P.bag_gens[c].size() > 0) {
        Constraint Production_P_LB("Production_P_LB");
        Constraint Production_P_UB("Production_P_UB");
        Constraint Production_Q_LB("Production_Q_LB");
        Constraint Production_Q_UB("Production_Q_UB");
        // 5A
        Production_P_UB = Pg.in(P.bag_gens[c], T) - grid.pg_max.in(P.bag_gens[c], T)*On_off.in(P.bag_gens[c],T);
        Production_P_LB = Pg.in(P.bag_gens[c], T) - grid.pg_min.in(P.bag_gens[c], T)*On_off.in(P.bag_gens[c],T);
        Subr.add_constraint(Production_P_UB <=0);
        Subr.add_constraint(Production_P_LB >= 0);

        Production_Q_UB = Qg.in(P.bag_gens[c], T) - grid.qg_max.in(P.bag_gens[c], T)*On_off.in(P.bag_gens[c],T);
        Production_Q_LB = Qg.in(P.bag_gens[c], T) - grid.qg_min.in(P.bag_gens[c], T)*On_off.in(P.bag_gens[c],T);
        Subr.add_constraint(Production_Q_UB <= 0);
        Subr.add_constraint(Production_Q_LB >= 0);
    }
    // 5C
    if (P.bag_gens[c].size() > 0) {
        for (int t = 1; t < T; t++) {
            Constraint Ramp_up("Ramp_up_constraint_" + to_string(t));
            Constraint Ramp_down("Ramp_down_constraint_" +to_string(t));
            Ramp_up =  Pg.in_at(P.bag_gens[c], t);
            Ramp_up -= Pg.in_at(P.bag_gens[c], t-1);
            Ramp_up -= rate_ramp*On_off.in_at(P.bag_gens[c], t-1);
            Ramp_up -= rate_switch*(1 - On_off.in_at(P.bag_gens[c], t));

            Ramp_down =  Pg.in_at(P.bag_gens[c], t-1);
            Ramp_down -= Pg.in_at(P.bag_gens[c], t);
            Ramp_down -= rate_ramp*On_off.in_at(P.bag_gens[c], t);
            Ramp_down -= rate_switch*(1 - On_off.in_at(P.bag_gens[c], t-1));

            Subr.add_constraint(Ramp_up <= 0);
            Subr.add_constraint(Ramp_down <= 0);
        }
    }
    // set the initial state of generators.
    for(auto g: P.bag_gens[c]) {
        Constraint gen_initial("gen_" + g->_name + to_string(c) + ",0");
        gen_initial +=  On_off( g->_name + ",0");
        Subr.add_constraint(gen_initial == 1);
    }
    /* solve it! */
    solver solve_Subr(Subr, cplex);
    solve_Subr.run();
    // COLLECT THE LINKED VARIABLES
    //std::string name = Xii.in(P.bag_bus[c],T).get_name();
    //Xii_log= (*(var<Real>*) Subr.get_var(name));
    //name = R_Xij.in(P.bag_bus_pairs_union[c], T).get_name();
    //R_Xij_log = (*(var<Real>*) Subr.get_var(name));
    //name = Im_Xij.in(P.bag_bus_pairs_union[c], T).get_name();
    //Im_Xij_log = (*(var<Real>*) Subr.get_var(name));
    return Subr._obj_val;
}

int main (int argc, const char * argv[])
{
    // Decompose
    const char* fname;
    double l = 0;
    if (argc >= 2) {
        fname = argv[1];
        l = atof(argv[2]);
    }
    else {
        //fname = "../../data_sets/Power/nesta_case5_pjm.m";
        //fname = "../../data_sets/Power/nesta_case30_ieee.m";
        //fname = "../../data_sets/Power/nesta_case6_c.m";
        fname = "../../data_sets/Power/nesta_case5_pjm.m";
        //fname = "../../data_sets/Power/nesta_case3_lmbd.m";
        //fname = "../../data_sets/Power/nesta_case300_ieee.m";
        //fname = "../../data_sets/Power/nesta_case1354_pegase.m";
        //fname = "../../data_sets/Power/nesta_case14_ieee.m";
        //fname = "../../data_sets/Power/nesta_case57_ieee.m";
        l = 1;
    }
    PowerNet grid;
    grid.readgrid(fname);
    int nbparts = 2;
    //GRAPH PARTITION
    auto bus_pairs = grid.get_bus_pairs();
    Partition P;
    P.get_ncut(grid, nbparts);
    // Schedule Parameters
    unsigned T = 1;
    param<Real> rate_ramp("rate_ramp");
    param<Real> rate_switch("rate_switch");
    param<Real> min_up("min_up");
    param<Real> min_down("min_down");
    param<Real> cost_up("cost_up");
    param<Real> cost_down("cost_down");
    for (auto g: grid.gens) {
        rate_ramp(g->_name) = max(grid.pg_min(g->_name).getvalue(), 0.25*grid.pg_max(g->_name).getvalue());
        rate_switch(g->_name) = max(grid.pg_min(g->_name).getvalue(), 0.25*grid.pg_max(g->_name).getvalue());
    }
    min_up = 2;
    min_down = 2;
    cost_up = 50;
    cost_down = 30;
    grid.time_expand(T);
    //rate_ramp.time_expand(T);
    //rate_switch.time_expand(T);
    /** Variables */
    vector<var<Real>> R_Xij;
    vector<var<Real>> Im_Xij;
    vector<var<Real>> Xii;
    vector<var<Real>> Pg;
    vector<var<Real>> Qg;
    // Commitment variables
    vector<var<bool>>  On_off;
    vector<var<bool>>  Start_up;
    vector<var<bool>>  Shut_down;
    for (int c = 0; c < nbparts; c++) {
        var<Real>  bag_Xii("Xii_"+ to_string(c), grid.w_min, grid.w_max);
        //var<Real>  bag_Xii("Xii_"+ to_string(c), grid.w_min.in(P.bag_bus_union_out[c], T), grid.w_max.in(P.bag_bus_union_out[c], T));
        //bag_Xii.set_size(T*P.bag_bus_union_out[c].size());
        bag_Xii.initialize_all(1.001);
        Xii.push_back(bag_Xii);

        var<Real>  bag_R_Xij("R_Xij_"+ to_string(c), grid.wr_min, grid.wr_max);
        var<Real>  bag_Im_Xij("Im_Xij_"+ to_string(c), grid.wi_min, grid.wi_max);
        //var<Real>  bag_R_Xij("R_Xij_"+ to_string(c), grid.wr_min.in(P.bag_bus_pairs_union[c], T), grid.wr_max.in(P.bag_bus_pairs_union[c], T));
        //var<Real>  bag_Im_Xij("Im_Xij_"+ to_string(c), grid.wi_min.in(P.bag_bus_pairs_union[c], T), grid.wi_max.in(P.bag_bus_pairs_union[c], T));
        //bag_R_Xij.set_size(T*P.bag_bus_pairs_union[c].size());
        //bag_Im_Xij.set_size(T*P.bag_bus_pairs_union[c].size());
        R_Xij.push_back(bag_R_Xij);
        bag_R_Xij.initialize_all(1.0);
        Im_Xij.push_back(bag_Im_Xij);
        if (P.bag_gens[c].size() > 0) {
            var<Real>  bag_Pg("Pg_" + to_string(c), grid.pg_min, grid.pg_max);
            var<Real>  bag_Qg("Qg_" + to_string(c), grid.qg_min, grid.qg_max);
            var<bool>  bag_Onoff("On_off_" + to_string(c));
            var<bool>  bag_Up("Start_up_" + to_string(c));
            var<bool>  bag_Down("Shut_down_" + to_string(c));
            //bag_Pg.set_size(T*P.bag_gens[c].size());
            //bag_Qg.set_size(T*P.bag_gens[c].size());
            Pg.push_back(bag_Pg);
            Qg.push_back(bag_Qg);
            On_off.push_back(bag_Onoff);
            Start_up.push_back(bag_Up);
            Shut_down.push_back(bag_Down);
        }
        else {
            var<Real> empty("empty");
            var<bool> empty2("empty");
            empty.set_size(0);
            Pg.push_back(empty);
            Qg.push_back(empty);
            On_off.push_back(empty2);
            Start_up.push_back(empty2);
            Shut_down.push_back(empty2);
        }
    }
/////////////////// DEFINE LAGRANGE MULTIPLIERS  ////////////////////////////////
    param<Real> R_lambda("R_lambda");
    param<Real> Im_lambda("Im_lambda");
    param<Real> lambda("lambda");
    R_lambda.initialize_all(0);
    Im_lambda.initialize_all(0);
    lambda.initialize_all(0);
//////////////////  CONVERGENCE INFORMATION /////////////////////////
    param<Real> R_lambda_log("R_lambda_log");
    param<Real> Im_lambda_log("Im_lambda_log");
    param<Real> lambda_log("lambda_log");
    R_lambda_log.in(P.inter_pairs, T);
    Im_lambda_log.in(P.inter_pairs, T);
    lambda_log.in(P.inter_pairs, T);

    vector<param<Real>> Xii_log;
    vector<param<Real>> R_Xij_log;
    vector<param<Real>> Im_Xij_log;
    for (int c =0; c < nbparts; c++) {
        param<Real> Xii_C_log("Xii_C_log_" + to_string(c));
        Xii_log.push_back(Xii_C_log);

        param<Real> Im_Xij_C_log("Im_Xij_C_log" + to_string(c));
        param<Real> R_Xij_C_log("R_Xij_C_log" +  to_string(c));
        R_Xij_log.push_back(R_Xij_C_log);
        Im_Xij_log.push_back(Im_Xij_C_log);
    }
/////////////////////////////////// INITIALIZATION ///////////////////////////////////////////
    auto wall0 = get_wall_time();
    auto cpu0  = get_cpu_time();
    double lb_cts = getdual_relax(grid, T, P, nbparts, rate_ramp, rate_switch, min_up, min_down, cost_up, cost_down, Pg, Qg,
                                  Start_up, Shut_down, On_off, Xii, R_Xij,  Im_Xij, R_lambda, Im_lambda, lambda);
//    double dual = 0.0;
//    double value_dual[nbparts];
//    for (int c = 0; c < nbparts; c++) {
//        value_dual[c] = subproblem(grid, T, P, c, rate_ramp, rate_switch, min_up, min_down, cost_up, cost_down, Pg[c], Qg[c],
//                                   Start_up[c], Shut_down[c], On_off[c], Xii[c], R_Xij[c], Im_Xij[c],
//                                   R_lambda, Im_lambda, lambda, Xii_log[c], R_Xij_log[c], Im_Xij_log[c]);
//        dual += value_dual[c];
//    }
    return 0;
}
