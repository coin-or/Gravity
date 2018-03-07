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
double getdual_relax(PowerNet& grid, const unsigned T,
                     const param<Real>& rate_ramp, const param<Real>& rate_switch,
                     const param<Real>& min_up, const param<Real>& min_down,
                     const param<Real>& cost_up, const param<Real>& cost_down,
                     var<Real>& Pg, var<Real>& Qg, var<bool>& Start_up, var<bool>& Shut_down, var<bool>& On_off,
                     var<Real>& Xii, var<Real>& R_Xij, var<Real>& Im_Xij,
                     param<Real>& R_lambda, param<Real>& Im_lambda, param<Real>& lambda)
{
    // Grid Parameters
    const auto bus_pairs = grid.get_bus_pairs();
    Model ACUC("ACUC Model");
    ACUC.add_var(Pg.in(grid.gens,T));
    ACUC.add_var(Qg.in(grid.gens,T));
    ACUC.add_var(On_off.in(grid.gens, T));
    ACUC.add_var(Start_up.in(grid.gens, T));
    ACUC.add_var(Shut_down.in(grid.gens, T));
    ACUC.add_var(Xii.in(grid.nodes, T));
    ACUC.add_var(R_Xij.in(bus_pairs, T));
    ACUC.add_var(Im_Xij.in(bus_pairs, T));
    //power flow vars are treated as auxiliary vars.
    var<Real> Pf_from("Pf_from", grid.S_max);
    var<Real> Qf_from("Qf_from", grid.S_max);
    var<Real> Pf_to("Pf_to", grid.S_max);
    var<Real> Qf_to("Qf_to", grid.S_max);
    ACUC.add_var(Pf_from.in(grid.arcs, T));
    ACUC.add_var(Qf_from.in(grid.arcs, T));
    ACUC.add_var(Pf_to.in(grid.arcs, T));
    ACUC.add_var(Qf_to.in(grid.arcs, T));
    /* Construct the objective function*/
    func_ obj;
    for (int t = 0; t < T; t++) {
        for (auto g:grid.gens) {
            if (g->_active) {
                if (t == 0) {
                    string name = g->_name + ",0";
                    obj += grid.c1(name)*Pg(name)+ grid.c2(name)*Pg(name)*Pg(name) +grid.c0(name)*On_off(name);
                }
                else {
                    if (g->_active) {
                        string name = g->_name + ","+ to_string(t);
                        obj += grid.c1(name)*Pg(name)+ grid.c2(name)*Pg(name)*Pg(name) + grid.c0(name)*On_off(name);
                        obj += cost_up.getvalue()*Start_up(name)+ cost_down.getvalue()*Shut_down(name);
                    }
                }
            }
        }
    }
    ACUC.set_objective(min(obj));
    for (int t= 0; t < T; t++){
        Constraint SOC("SOC_" + to_string(t));
        SOC =  power(R_Wij, 2) + power(Im_Wij, 2) - Wii.from()*Wii.to() ;
        ACUC.add_constraint(SOC.in_at(bus_pairs, t) <= 0);
    }
    //KCL
    for (int c = 0; c < nbparts; c++) {
        Constraint KCL_P("KCL_P"+ to_string(c));
        Constraint KCL_Q("KCL_Q"+ to_string(c));
        KCL_P = sum(Pf_from[c].out_arcs()) + sum(Pf_to[c].in_arcs())+ grid.pl -sum(Pg[c].in_gens())+ grid.gs*Xii[c];
        ACUC.add_constraint(KCL_P.in(P.bag_bus[c], T) == 0);
        KCL_Q  = sum(Qf_from[c].out_arcs()) + sum(Qf_to[c].in_arcs())+ grid.ql -sum(Qg[c].in_gens())-grid.bs*Xii[c];
        ACUC.add_constraint(KCL_Q.in(P.bag_bus[c], T) == 0);

        Constraint Flow_P_From("Flow_P_From_" + to_string(c));
        Flow_P_From = Pf_from[c]- (grid.g_ff*Xii[c].from()+ grid.g_ft*R_Xij[c].in_pairs() + grid.b_ft*Im_Xij[c].in_pairs());
        ACUC.add_constraint(Flow_P_From.in(P.bag_arcs_union_out[c], T) == 0);

        Constraint Flow_P_To("Flow_P_To" + to_string(c));
        Flow_P_To = Pf_to[c] - (grid.g_tt*Xii[c].to() + grid.g_tf*R_Xij[c].in_pairs() - grid.b_tf*Im_Xij[c].in_pairs());
        ACUC.add_constraint(Flow_P_To.in(P.bag_arcs_union_in[c], T) == 0);

        Constraint Flow_Q_From("Flow_Q_From" + to_string(c));
        Flow_Q_From = Qf_from[c] - (grid.g_ft*Im_Xij[c].in_pairs() - grid.b_ff*Xii[c].from() - grid.b_ft*R_Xij[c].in_pairs());
        ACUC.add_constraint(Flow_Q_From.in(P.bag_arcs_union_out[c], T) == 0);

        Constraint Flow_Q_To("Flow_Q_To" + to_string(c));
        Flow_Q_To = Qf_to[c] + (grid.b_tt*Xii[c].to() + grid.b_tf*R_Xij[c].in_pairs() + grid.g_tf*Im_Xij[c].in_pairs());
        ACUC.add_constraint(Flow_Q_To.in(P.bag_arcs_union_in[c], T) == 0);

        Constraint Thermal_Limit_from("Thermal_Limit_from" + to_string(c));
        Thermal_Limit_from = power(Pf_from[c], 2) + power(Qf_from[c], 2);
        Thermal_Limit_from <= power(grid.S_max,2);
        ACUC.add_constraint(Thermal_Limit_from.in(P.bag_arcs_union_out[c], T));

        Constraint Thermal_Limit_to("Thermal_Limit_to" + to_string(c));
        Thermal_Limit_to = power(Pf_to[c], 2) + power(Qf_to[c], 2);
        Thermal_Limit_to <= power(grid.S_max,2);
        ACUC.add_constraint(Thermal_Limit_to.in(P.bag_arcs_union_in[c], T));
    }
//    ///* Phase Angle Bounds constraints */
    for (int c = 0; c < nbparts; c++) {
        if (P.bag_bus_pairs_union_directed[c].size() > 0) {
            Constraint PAD_UB("PAD_UB_"+to_string(c));
            PAD_UB = Im_Xij[c]- grid.tan_th_max*R_Xij[c];
            ACUC.add_constraint(PAD_UB.in(P.bag_bus_pairs_union_directed[c], T) <= 0);

            Constraint PAD_LB("PAD_LB_"+to_string(c));
            PAD_LB = Im_Xij[c]- grid.tan_th_min*R_Xij[c];
            ACUC.add_constraint(PAD_LB.in(P.bag_bus_pairs_union_directed[c], T) >= 0);
        }
    }
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
//    ////Ramp Rate
//    for (int c = 0; c < nbparts; c++) {
//        if (P.bag_gens[c].size() > 0) {
//            Constraint Production_P_LB("Production_P_LB_"+ to_string(c));
//            Constraint Production_P_UB("Production_P_UB_"+ to_string(c));
//            Constraint Production_Q_LB("Production_Q_LB_"+ to_string(c));
//            Constraint Production_Q_UB("Production_Q_UB_"+ to_string(c));
//            // 5A
//            Production_P_UB = Pg[c] - grid.pg_max*On_off[c];
//            Production_P_LB = Pg[c] - grid.pg_min*On_off[c];
//            ACUC.add_constraint(Production_P_UB.in(P.bag_gens[c], T) <=0);
//            ACUC.add_constraint(Production_P_LB.in(P.bag_gens[c], T) >= 0);
//
//            Production_Q_UB = Qg[c] - grid.qg_max*On_off[c];
//            Production_Q_LB = Qg[c] - grid.qg_min*On_off[c];
//            ACUC.add_constraint(Production_Q_UB.in(P.bag_gens[c], T) <= 0);
//            ACUC.add_constraint(Production_Q_LB.in(P.bag_gens[c], T) >= 0);
//        }
//    }
//    // 5C
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
//    // set the initial state of generators.
//    for (int c = 0; c < nbparts; c++) {
//        for(auto g: P.bag_gens[c]) {
//            Constraint gen_initial("gen_" + g->_name + to_string(c) + ",0");
//            gen_initial +=  On_off[c]( g->_name + ",0");
//            ACUC.add_constraint(gen_initial == 1);
//        }
//    }
    // Linking constraints
    for (const auto& a: P.G_part.arcs) {
        Constraint Link_R("link_R_"+a->_name);
        Link_R = R_Xij[a->_src->_id].in_pairs() - R_Xij[a->_dest->_id].in_pairs();
        ACUC.add_constraint(Link_R.in(a->_intersection_clique, T) ==0);

        Constraint Link_Im("link_Im_"+a->_name);
        Link_Im = Im_Xij[a->_src->_id].in_pairs() - Im_Xij[a->_dest->_id].in_pairs();
        ACUC.add_constraint(Link_Im.in(a->_intersection_clique, T)==0);
        Constraint Link_Xii("link_Xii_" + a->_name);
        Link_Xii = Xii[a->_src->_id].to() - Xii[a->_dest->_id].to();
        ACUC.add_constraint(Link_Xii.in(a->_intersection_clique, T)==0);
    }
    /* Solver selection */
    solver cpx_acuc(ACUC, cplex);
    //solver cpx_acuc(ACUC, ipopt);

    bool relax =true;
    int output = 1;
    double tol = 1e-6;
    cpx_acuc.run(output, relax, tol);
    cout << "the continuous relaxation bound is: " << ACUC._obj_val << endl;

    for (const auto& a: P.G_part.arcs) {
        auto consR = ACUC.get_constraint("link_R_"+a->_name);
        auto consIm = ACUC.get_constraint("link_Im_"+a->_name);
        auto cons = ACUC.get_constraint("link_Xii_"+a->_name);
        int i = 0;
        for (unsigned t = 0; t < T; t++) {
            for (auto& line: a->_intersection_clique) {
                string name =line->_name+","+to_string(t);
                R_lambda(name) = -consR->_dual[i];
                Im_lambda(name) = -consIm->_dual[i];
                lambda(name) = -cons->_dual[i];
                cout << "dual of  R_lambda_" << name << " " << consR->_dual[i] << endl;
                cout << "dual of  Im_lambda_" << name << " " << consIm->_dual[i] << endl;
                cout << "dual of  lambda_" << name << " " << cons->_dual[i] << endl;
                i++;
            }
        }
    }
    return ACUC._obj_val;
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
    //GRAPH PARTITION
    auto bus_pairs = grid.get_bus_pairs();
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
    rate_ramp.time_expand(T);
    rate_switch.time_expand(T);
    /** Variables */
    // POWER GENERATION
    var<Real> Pg("Pg", grid->pg_min, grid->pg_max);
    var<Real> Qg ("Qg", grid->qg_min, grid->qg_max);

    //Lifted variables.
    var<Real>  R_Wij("R_Wij", grid->wr_min, grid->wr_max); // real part of Wij
    var<Real>  Im_Wij("Im_Wij", grid->wi_min, grid->wi_max); // imaginary part of Wij.
    var<Real>  Wii("Wii", grid->w_min, grid->w_max);
    R_Wij.initialize_all(1.0);
    Wii.initialize_all(1.001);

    // Commitment variables
    var<bool>  On_off("On_off");
    var<bool>  Start_up("Start_up");
    var<bool>  Shut_down("Shut_down");
    ACUC.add_var(On_off.in(grid->gens, T));
    ACUC.add_var(Start_up.in(grid->gens, T));
    ACUC.add_var(Shut_down.in(grid->gens, T));
/////////////////// DEFINE LAGRANGE MULTIPLIERS  ////////////////////////////////
    param<Real> R_lambda("R_lambda");
    param<Real> Im_lambda("Im_lambda");
    param<Real> lambda("lambda");
    R_lambda.set_size(T*P.inter_pairs.size());
    Im_lambda.set_size(T*P.inter_pairs.size());
    lambda.set_size(T*P.inter_pairs.size());

    R_lambda.initialize_all(0);
    Im_lambda.initialize_all(0);
    lambda.initialize_all(0);
    double lb_cts = getdual_relax(grid, T, rate_ramp, rate_switch, min_up, min_down, cost_up, cost_down, Pg, Qg,
                                  Start_up, Shut_down, On_off, Xii, R_Xij,  Im_Xij, R_lambda, Im_lambda, lambda);

    cout << "The initial Lower bound of the ACUC problem is: " << LB  << endl;
    return 0;
}
