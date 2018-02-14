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
/** INITIALISE SUBPROBLEM MODEL */
//double subproblem(PowerNet& grid, const Partition& P, unsigned c, var<Real>& Pg, var<Real>& Qg,
//                  var<Real>& Xii, var<Real>& R_Xij, var<Real>& Im_Xij,
//                  param<Real>& R_lambda_sep, param<Real>& Im_lambda_sep, param<Real>& lambda_sep,
//                  param<Real>& Xii_log, param<Real>& R_Xij_log, param<Real>& Im_Xij_log)
//{
//    DebugOff("Solving subproblem associated with maximal clique "<< c << endl);
//    Model Subr("Subr");
//    Subr.add_var(Pg);
//    Subr.add_var(Qg);
//    Subr.add_var(Xii);
//    Subr.add_var(R_Xij);
//    Subr.add_var(Im_Xij);
//    /* Construct the objective function*/
//    func_ obj;
//    for (const auto &g: P.bag_gens[c]) {
//        if (g->_active) {
//            obj += grid.c1(g->_name)*Pg(g->_name)+ grid.c2(g->_name)*Pg(g->_name)*Pg(g->_name)+grid.c0(g->_name);
//        }
//    }
//    //RELAXED TERMS
//    auto bag = P.G_part.get_node(to_string(c));
//    for (const auto &a: bag->get_out()) {
//        for (const auto &pair: a->_intersection_clique) {
//            obj += R_lambda_sep(pair->_name)*R_Xij(pair->_name);
//            obj += Im_lambda_sep(pair->_name)*Im_Xij(pair->_name);
//            obj += lambda_sep(pair->_name)*Xii(pair->_dest->_name);
//            DebugOn("+Xii[" << c <<"]("<< pair->_dest->_name <<")" << endl);
//        }
//    }
//    for (auto a: bag->get_in()) {
//        for (auto pair: a->_intersection_clique) {
//            obj -= R_lambda_sep(pair->_name)*R_Xij(pair->_name);
//            obj -= Im_lambda_sep(pair->_name)*Im_Xij(pair->_name);
//            obj -= lambda_sep(pair->_name)*Xii(pair->_dest->_name);
//            DebugOn("-Xii[" << c <<"]("<< pair->_dest->_name <<")" << endl);
//        }
//    }
//
//    Subr.set_objective(min(obj));
//    if (P.bag_bus_pairs_union_directed[c].size() > 0) {
//        Constraint SOC("SOC_" + to_string(c));
//        SOC =  power(R_Xij, 2)+ power(Im_Xij, 2) - Xii.from()*Xii.to() ;
//        Subr.add_constraint(SOC.in(P.bag_bus_pairs_union_directed[c]) <= 0);
//    }
//    //* KCL */
//    for (const auto& bus: P.bag_bus[c]) {
//        Constraint KCL_P("KCL_P"+bus->_name);
//        Constraint KCL_Q("KCL_Q"+bus->_name);
//
//        for (const auto& a: bus->get_out()) {
//            KCL_P += grid.g_ff(a->_name)*Xii(a->_src->_name)
//                     +grid.g_ft(a->_name)*R_Xij(a->_src->_name+","+a->_dest->_name)
//                     +grid.b_ft(a->_name)*Im_Xij(a->_src->_name+","+a->_dest->_name);
//
//            KCL_Q += -1*grid.b_ff(a->_name)*Xii(a->_src->_name)
//                     - grid.b_ft(a->_name)*R_Xij(a->_src->_name+","+a->_dest->_name)
//                     +grid.g_ft(a->_name)*Im_Xij(a->_src->_name+","+a->_dest->_name);
//        }
//
//        for (auto a: bus->get_in()) {
//            KCL_P  += grid.g_tt(a->_name)*Xii(a->_dest->_name)
//                      +grid.g_tf(a->_name)*R_Xij(a->_src->_name+","+a->_dest->_name)
//                      -grid.b_tf(a->_name)*Im_Xij(a->_src->_name+","+a->_dest->_name);
//
//            KCL_Q  -= grid.b_tt(a->_name)*Xii(a->_dest->_name)
//                      + grid.b_tf(a->_name)*R_Xij(a->_src->_name+","+a->_dest->_name)
//                      + grid.g_tf(a->_name)*Im_Xij(a->_src->_name+","+a->_dest->_name);
//        }
//        if(bus->_has_gen) {
//            KCL_P += bus->pl()- sum(Pg.in(bus->_gen));
//            KCL_Q += bus->ql()- sum(Qg.in(bus->_gen));
//        }
//        else {
//            KCL_P += bus->pl();
//            KCL_Q += bus->ql();
//        }
//
//        /* Shunts */
//        KCL_P += bus->gs()*(Xii(bus->_name));
//        KCL_Q -= bus->bs()*(Xii(bus->_name));
//
//        Subr.add_constraint(KCL_P == 0);
//        Subr.add_constraint(KCL_Q == 0);
//    }
//    /* Phase Angle Bounds constraints */
//    vector<index_pair*> pairs = P.bag_bus_pairs_union_directed[c];
//    Constraint PAD_UB("PAD_UB_"+to_string(c));
//    PAD_UB = Im_Xij- grid.tan_th_max*R_Xij;
//    Subr.add_constraint(PAD_UB.in(pairs) <= 0);
//
//    Constraint PAD_LB("PAD_LB_"+to_string(c));
//    PAD_LB = Im_Xij- grid.tan_th_min*R_Xij;
//    Subr.add_constraint(PAD_LB.in(pairs) >= 0);
//
//    /* Thermal Limit Constraints */
//    if (P.bag_arcs_union_out[c].size() > 0) {
//        /* Thermal Limit Constraints */
//        Constraint Thermal_Limit_from("Thermal_Limit_from"+to_string(c));
//        Constraint Aux1("Aux1_"+to_string(c));
//        Constraint Aux2("Aux2_"+to_string(c));
//        var<Real> t1("t1_"+to_string(c));
//        var<Real> t2("t2_"+to_string(c));
//        Subr.add_var(t1^(P.bag_arcs_union_out[c].size()));
//        Subr.add_var(t2^(P.bag_arcs_union_out[c].size()));
//        Aux1 += t1 - grid.g_ff*Xii.from()- grid.g_ft*R_Xij.in_pairs() -grid.b_ft*Im_Xij.in_pairs();
//        Aux2 += t2 - grid.b_ff*Xii.from() - grid.b_ft*R_Xij.in_pairs() + grid.g_ft*Im_Xij.in_pairs();
//        Subr.add_constraint(Aux1.in(P.bag_arcs_union_out[c]) = 0);
//        Subr.add_constraint(Aux2.in(P.bag_arcs_union_out[c]) = 0);
//
//        Thermal_Limit_from += power(t1, 2)+ power(t2,2) - power(grid.S_max, 2);
//        Subr.add_constraint(Thermal_Limit_from.in(P.bag_arcs_union_out[c]) <= 0);
//    }
//    if (P.bag_arcs_union_in[c].size() > 0) {
//        Constraint Thermal_Limit_to("Thermal_Limit_to" + to_string(c));
//        Constraint Aux3("Aux3_"+to_string(c));
//        Constraint Aux4("Aux4_"+to_string(c));
//        var<Real> t3("t3_"+to_string(c));
//        var<Real> t4("t4_"+to_string(c));
//        Subr.add_var(t3^(P.bag_arcs_union_in[c].size()));
//        Subr.add_var(t4^(P.bag_arcs_union_in[c].size()));
//        Aux3 += t3 - grid.g_tt*Xii.to() - grid.g_tf*R_Xij.in_pairs() - grid.b_tf*Im_Xij.in_pairs();
//        Aux4 += t4 - grid.b_tt*Xii.to() -  grid.b_tf*R_Xij.in_pairs() -  grid.g_tf*Im_Xij.in_pairs();
//        Subr.add_constraint(Aux3.in(P.bag_arcs_union_in[c]) = 0);
//        Subr.add_constraint(Aux4.in(P.bag_arcs_union_in[c]) = 0);
//        Thermal_Limit_to += power(t3, 2) + power(t4, 2)  - power(grid.S_max, 2);
//        Subr.add_constraint(Thermal_Limit_to.in(P.bag_arcs_union_in[c]) <= 0);
//    }
//    /* solve it! */
//    solver solve_Subr(Subr, cplex);
//    //solver solve_Subr(Subr, ipopt);
//    solve_Subr.run();
//
//    // COLLECT THE LINKED VARIABLES
//    //cout << "values: " << Xii_log.getvalue() << endl;
//    auto val = (*(var<Real>*) Subr.get_var("Xii_"+ to_string(c))); //.in(bag_bus_out);
//    for (const auto& b: P.bag_bus[c]) {
//        Xii_log(b->_name)  = val(b->_name).getvalue();
//    }
//    for (const auto& b: P.bag_bus_out[c]) {
//        Xii_log(b->_name)  = val(b->_name).getvalue();
//    }
//    auto R_val = (*(var<Real>*) Subr.get_var("R_Xij_"+to_string(c))).in(P.bag_bus_pairs_neighbour[c]);
//    auto Im_val = (*(var<Real>*) Subr.get_var("Im_Xij_"+to_string(c))).in(P.bag_bus_pairs_neighbour[c]);
//    for(const auto& b: P.bag_bus_pairs_neighbour[c]) {
//        R_Xij_log(b->_name) = R_val(b->_name).getvalue();
//        Im_Xij_log(b->_name) = Im_val(b->_name).getvalue();
//    }
//    return Subr._obj_val;
//}

void reform_inout (PowerNet& grid, unsigned nbparts) {
    Partition P;
    P.get_ncut(grid, nbparts);
    /** Build model */
    Model CLT("A hierarchicial Model");
    /** Variables */
    vector<var<Real>> R_Xij;
    vector<var<Real>> Im_Xij;
    vector<var<Real>> Xii;
    vector<var<Real>> Pg;
    vector<var<Real>> Qg;
    for (int c = 0; c < nbparts; c++) {
        var<Real>  bag_Xii("Xii_"+ to_string(c), grid.w_min, grid.w_max);
        bag_Xii.set_size(P.bag_bus_union[c].size());
        bag_Xii.initialize_all(1.001);
        Xii.push_back(bag_Xii);

        var<Real>  bag_R_Xij("R_Xij_"+ to_string(c), grid.wr_min, grid.wr_max);
        var<Real>  bag_Im_Xij("Im_Xij_"+ to_string(c), grid.wi_min, grid.wi_max);
        CLT.add_var(bag_R_Xij.in(P.bag_bus_pairs_union[c]));
        CLT.add_var(bag_Im_Xij.in(P.bag_bus_pairs_union[c]));
        bag_R_Xij.initialize_all(1.0);
        R_Xij.push_back(bag_R_Xij);
        Im_Xij.push_back(bag_Im_Xij);
        if (P.bag_gens[c].size() > 0) {
            var<Real>  bag_Pg("Pg_" + to_string(c), grid.pg_min, grid.pg_max);
            var<Real>  bag_Qg("Qg_" + to_string(c), grid.qg_min, grid.qg_max);
            CLT.add_var(bag_Pg.in(P.bag_gens[c]));
            CLT.add_var(bag_Qg.in(P.bag_gens[c]));
            Pg.push_back(bag_Pg);
            Qg.push_back(bag_Qg);
        }
        else {
            var<Real> empty("empty");
            empty.set_size(0);
            Pg.push_back(empty);
            Qg.push_back(empty);
        }
    }
    func_ obj;
    for (int c = 0; c < nbparts; c++) {
        for (auto g: P.bag_gens[c]) {
            if (g->_active) {
                obj += grid.c1(g->_name)*Pg[c](g->_name)+ grid.c2(g->_name)*Pg[c](g->_name)*Pg[c](g->_name)+grid.c0(g->_name);
            }
        }
    }
    CLT.set_objective(min(obj));
    for (int c = 0; c < nbparts; c++) {
        if (P.bag_bus_pairs_union_directed[c].size() > 0) {
            Constraint SOC("SOC_" + to_string(c));
            SOC =  power(R_Xij[c], 2)+ power(Im_Xij[c], 2) - Xii[c].from()*Xii[c].to() ;
            CLT.add_constraint(SOC.in(P.bag_bus_pairs_union_directed[c]) <= 0);
        }
    }
//     //KCL
//    for (int c = 0; c < nbparts; c++) {
//        for (const auto& bus: P.bag_bus[c]) {
//            Constraint KCL_P("KCL_P"+bus->_name);
//            Constraint KCL_Q("KCL_Q"+bus->_name);
//            for (const auto& a: bus->get_out()) {
//                KCL_P += grid.g_ff(a->_name)*Xii[c](a->_src->_name)
//                         +grid.g_ft(a->_name)*R_Xij[c](a->_src->_name+","+a->_dest->_name)
//                         +grid.b_ft(a->_name)*Im_Xij[c](a->_src->_name+","+a->_dest->_name);
//
//                KCL_Q += -1*grid.b_ff(a->_name)*Xii[c](a->_src->_name)
//                         - grid.b_ft(a->_name)*R_Xij[c](a->_src->_name+","+a->_dest->_name)
//                         +grid.g_ft(a->_name)*Im_Xij[c](a->_src->_name+","+a->_dest->_name);
//            }
//            for (auto a: bus->get_in()) {
//                KCL_P  += grid.g_tt(a->_name)*Xii[c](a->_dest->_name)
//                          +grid.g_tf(a->_name)*R_Xij[c](a->_src->_name+","+a->_dest->_name)
//                          -grid.b_tf(a->_name)*Im_Xij[c](a->_src->_name+","+a->_dest->_name);
//
//                KCL_Q  -= grid.b_tt(a->_name)*Xii[c](a->_dest->_name)
//                          + grid.b_tf(a->_name)*R_Xij[c](a->_src->_name+","+a->_dest->_name)
//                          + grid.g_tf(a->_name)*Im_Xij[c](a->_src->_name+","+a->_dest->_name);
//            }
//            if(bus->_has_gen) {
//                KCL_P += bus->pl()- sum(Pg[c].in(bus->_gen));
//                KCL_Q += bus->ql()- sum(Qg[c].in(bus->_gen));
//            }
//            else {
//                KCL_P += bus->pl();
//                KCL_Q += bus->ql();
//            }
//            /* Shunts */
//            KCL_P += bus->gs()*(Xii[c](bus->_name));
//            KCL_Q -= bus->bs()*(Xii[c](bus->_name));
//            CLT.add_constraint(KCL_P == 0);
//            CLT.add_constraint(KCL_Q == 0);
//        }
//    }
    //Phase Angle Difference Constraints
//    for (int c = 0; c < nbparts; c++) {
//        vector<index_pair*> pairs = P.bag_bus_pairs_union_directed[c];
//        Constraint PAD_UB("PAD_UB_"+to_string(c));
//        PAD_UB = Im_Xij[c]- grid.tan_th_max*R_Xij[c];
//        CLT.add_constraint(PAD_UB.in(pairs) <= 0);
//
//        Constraint PAD_LB("PAD_LB_"+to_string(c));
//        PAD_LB = Im_Xij[c]- grid.tan_th_min*R_Xij[c];
//        CLT.add_constraint(PAD_LB.in(pairs) >= 0);
//    }
//    //THERMAL LIMIT CONSTRAINTS.
//    for (int c= 0; c < nbparts; c++) {
//        if (P.bag_arcs_union_out[c].size() > 0) {
//            /* Thermal Limit Constraints */
//            Constraint Thermal_Limit_from("Thermal_Limit_from"+to_string(c));
//            Constraint Aux1("Aux1_"+to_string(c));
//            Constraint Aux2("Aux2_"+to_string(c));
//            var<Real> t1("t1_"+to_string(c));
//            var<Real> t2("t2_"+to_string(c));
//            CLT.add_var(t1.in(P.bag_arcs_union_out[c]));
//            CLT.add_var(t2.in(P.bag_arcs_union_out[c]));
//            Aux1 += t1 - grid.g_ff*Xii[c].from()- grid.g_ft*R_Xij[c].in_pairs() -grid.b_ft*Im_Xij[c].in_pairs();
//            Aux2 += t2 - grid.b_ff*Xii[c].from() - grid.b_ft*R_Xij[c].in_pairs() + grid.g_ft*Im_Xij[c].in_pairs();
//            CLT.add_constraint(Aux1.in(P.bag_arcs_union_out[c]) == 0);
//            CLT.add_constraint(Aux2.in(P.bag_arcs_union_out[c]) == 0);
//
//            Thermal_Limit_from += power(t1, 2)+ power(t2,2) - power(grid.S_max, 2);
//            CLT.add_constraint(Thermal_Limit_from.in(P.bag_arcs_union_out[c]) <= 0);
//        }
//        if (P.bag_arcs_union_in[c].size() > 0) {
//            Constraint Thermal_Limit_to("Thermal_Limit_to" + to_string(c));
//            Constraint Aux3("Aux3_"+to_string(c));
//            Constraint Aux4("Aux4_"+to_string(c));
//            var<Real> t3("t3_"+to_string(c));
//            var<Real> t4("t4_"+to_string(c));
//            CLT.add_var(t3.in(P.bag_arcs_union_in[c]));
//            CLT.add_var(t4.in(P.bag_arcs_union_in[c]));
//            Aux3 += t3 - grid.g_tt*Xii[c].to() - grid.g_tf*R_Xij[c].in_pairs() - grid.b_tf*Im_Xij[c].in_pairs();
//            Aux4 += t4 - grid.b_tt*Xii[c].to() - grid.b_tf*R_Xij[c].in_pairs() - grid.g_tf*Im_Xij[c].in_pairs();
//            CLT.add_constraint(Aux3.in(P.bag_arcs_union_in[c]) == 0);
//            CLT.add_constraint(Aux4.in(P.bag_arcs_union_in[c]) == 0);
//            Thermal_Limit_to += power(t3, 2) + power(t4, 2)  - power(grid.S_max, 2);
//            CLT.add_constraint(Thermal_Limit_to.in(P.bag_arcs_union_in[c]) <= 0);
//        }
//    }
//    //Linking Constraints.
//    for (const auto& a: P.G_part.arcs) {
//        Constraint Link_R("link_R_"+a->_name);
//        Link_R = R_Xij[a->_src->_id].in_pairs() - R_Xij[a->_dest->_id].in_pairs();
//        CLT.add_constraint(Link_R.in(a->_intersection_clique) == 0);
//
//        Constraint Link_Im("link_Im_"+a->_name);
//        Link_Im = Im_Xij[a->_src->_id].in_pairs() - Im_Xij[a->_dest->_id].in_pairs();
//        CLT.add_constraint(Link_Im.in(a->_intersection_clique)== 0);
//
//        Constraint Link_Xii("link_Xii_" + a->_name);
//        Link_Xii = Xii[a->_src->_id].to() - Xii[a->_dest->_id].to();
//        CLT.add_constraint(Link_Xii.in((a->_intersection_clique))== 0);
//    }
    /* solve it! */
    solver solve_CLT(CLT, cplex);
    solve_CLT.run();
}

//int inout(PowerNet& grid, unsigned nbparts, unsigned iter_limit) {
//    Partition P;
//    P.get_ncut(grid, nbparts);
//    /** Variables */
//    vector<var<Real>> R_Xij;
//    vector<var<Real>> Im_Xij;
//    vector<var<Real>> Xii;
//    vector<var<Real>> Pg;
//    vector<var<Real>> Qg;
//    for (int c = 0; c < nbparts; c++) {
//        var<Real>  bag_Xii("Xii_"+ to_string(c), grid.w_min.in(P.bag_bus_union[c]), grid.w_max.in(P.bag_bus_union[c]));
//        bag_Xii^(P.bag_bus_union[c].size());
//        bag_Xii.initialize_all(1.001);
//        Xii.push_back(bag_Xii);
//
//        var<Real>  bag_R_Xij("R_Xij_"+ to_string(c), grid.wr_min.in(P.bag_bus_pairs_union[c]), grid.wr_max.in(P.bag_bus_pairs_union[c]));
//        var<Real>  bag_Im_Xij("Im_Xij_"+ to_string(c), grid.wi_min.in(P.bag_bus_pairs_union[c]), grid.wi_max.in(P.bag_bus_pairs_union[c]));
//        bag_R_Xij^(P.bag_bus_pairs_union[c].size());
//        bag_Im_Xij^(P.bag_bus_pairs_union[c].size());
//        bag_R_Xij.initialize_all(1.0);
//        R_Xij.push_back(bag_R_Xij);
//        Im_Xij.push_back(bag_Im_Xij);
//        if (P.bag_gens[c].size() > 0) {
//            var<Real>  bag_Pg("Pg_" + to_string(c), grid.pg_min.in(P.bag_gens[c]), grid.pg_max.in(P.bag_gens[c]));
//            var<Real>  bag_Qg("Qg_" + to_string(c), grid.qg_min.in(P.bag_gens[c]), grid.qg_max.in(P.bag_gens[c]));
//            bag_Pg^(P.bag_gens[c].size());
//            bag_Qg^(P.bag_gens[c].size());
//            Pg.push_back(bag_Pg);
//            Qg.push_back(bag_Qg);
//        }
//        else {
//            var<Real> empty("empty");
//            empty.set_size(0);
//            Pg.push_back(empty);
//            Qg.push_back(empty);
//        }
//    }
///////////////////// DEFINE LAGRANGE MULTIPLIERS  ////////////////////////////////
//    param<Real> R_lambda_in("R_lambda_in");
//    param<Real> Im_lambda_in("Im_lambda_in");
//    param<Real> lambda_in("lambda_in");
//    R_lambda_in^(P.inter_pairs.size());
//    Im_lambda_in^(P.inter_pairs.size());
//    lambda_in^(P.inter_pairs.size());
//
//    param<Real> R_lambda_out("R_lambda_out");
//    param<Real> Im_lambda_out("Im_lambda_out");
//    param<Real> lambda_out("lambda_out");
//    R_lambda_out^(P.inter_pairs.size());
//    Im_lambda_out^(P.inter_pairs.size());
//    lambda_out^(P.inter_pairs.size());
//
//    param<Real> R_lambda_sep("R_lambda_sep");
//    param<Real> Im_lambda_sep("Im_lambda_sep");
//    param<Real> lambda_sep("lambda_sep");
//    R_lambda_sep^(P.inter_pairs.size());
//    Im_lambda_sep^(P.inter_pairs.size());
//    lambda_sep^(P.inter_pairs.size());
//
//    R_lambda_in.initialize_all(0);
//    Im_lambda_in.initialize_all(0);
//    lambda_in.initialize_all(0);
//    R_lambda_out.initialize_all(0);
//    Im_lambda_out.initialize_all(0);
//    lambda_out.initialize_all(0);
//    R_lambda_sep.initialize_all(0);
//    Im_lambda_sep.initialize_all(0);
//    lambda_sep.initialize_all(0);
/////////////////////////////////////// Master Problem ///////////////////////////////////
//    Model Master("Master");
//    /** param **/
//    param<Real> gamma_in("gamma_C_in");
//    param<Real> gamma_out("gamma_C_out");
//    param<Real> gamma_sep("gamma_C_sep");
//    gamma_in^nbparts;
//    gamma_out^nbparts;
//    gamma_sep^nbparts;
//    gamma_in.initialize_all(0.0);
//
//    /** Variables  */
//    var<Real> gamma_C("gamma_C");
//    Master.add_var(gamma_C^nbparts);
//
//    var<Real> R_lambda_var("R_lambda_var");
//    var<Real> Im_lambda_var("Im_lambda_var");
//    var<Real> lambda_var("lambda_var");
//    Master.add_var(R_lambda_var^(P.inter_pairs.size()));
//    Master.add_var(Im_lambda_var^(P.inter_pairs.size()));
//    Master.add_var(lambda_var^(P.inter_pairs.size()));
//
//    /////////** OBJ*//////////////
//    func_ master_obj;
//    master_obj += sum(gamma_C);
//    Master.set_objective(max(master_obj));
//    double bound = 100000;
//
//    Constraint UB("Upper_bound");
//    UB = sum(gamma_C) - bound;
//    Constraint LOB("Lower_bound");
//    LOB = sum(gamma_C);
//    Master.add_constraint(UB <= 0);
//    Master.add_constraint(LOB >= 0);
////////////////////  CONVERGENCE INFORMATION /////////////////////////
//    double alpha = 0.9;
//    double LBlog[iter_limit];
//    double UBlog[iter_limit];
//
//    double LDlog[iter_limit];
//    // LOG OF SOLUTIONS: Log here means the previous primal and dual solution.
//    param<Real> R_lambda_log("R_lambda_log");
//    param<Real> Im_lambda_log("Im_lambda_log");
//    param<Real> lambda_log("lambda_log");
//    R_lambda_log^(P.inter_pairs.size());
//    Im_lambda_log^(P.inter_pairs.size());
//    lambda_log^(P.inter_pairs.size());
//
//    vector<param<Real>> Xii_log;
//    vector<param<Real>> R_Xij_log;
//    vector<param<Real>> Im_Xij_log;
//    for (int c =0; c < nbparts; c++) {
//        param<Real> Xii_C_log("Xii_C_log_" + to_string(c));
//        Xii_C_log^(P.bag_bus_union[c].size());
//        Xii_log.push_back(Xii_C_log);
//
//        param<Real> Im_Xij_C_log("Im_Xij_C_log" + to_string(c));
//        param<Real> R_Xij_C_log("R_Xij_C_log" +  to_string(c));
//        R_Xij_C_log^(P.bag_bus_pairs_neighbour[c].size());
//        Im_Xij_C_log^(P.bag_bus_pairs_neighbour[c].size());
//        R_Xij_log.push_back(R_Xij_C_log);
//        Im_Xij_log.push_back(Im_Xij_C_log);
//    }
///////////////////////////////////// INITIALIZATION ///////////////////////////////////////////
//    auto wall0 = get_wall_time();
//    auto cpu0  = get_cpu_time();
//    double dual = 0.0;
//    double value_dual[nbparts];
//    for (int c = 0; c < nbparts; c++) {
//        value_dual[c] = subproblem(grid, P, c, Pg[c], Qg[c], Xii[c], R_Xij[c], Im_Xij[c],
//                                   R_lambda_sep, Im_lambda_sep, lambda_sep, Xii_log[c], R_Xij_log[c], Im_Xij_log[c]);
//        dual += value_dual[c];
//        gamma_in(c) = value_dual[c];
//    }
//    cout << "Initialization_value:   " << dual <<endl;
//    LBlog[0] = std::max(0.0, dual);
//    
/////////////////////// APPEND MORE CONSTRAINTS TO MAIN //////////////////////////////////
//    if (iter_limit > 0) {
//        for (const auto& bag: P.G_part.nodes) {
//            unsigned c = bag->_id;
//            Constraint Concavity("Iter_0_Concavity_" + to_string(c));
//            Concavity += gamma_C(c) - value_dual[c];
//            for (auto a: bag->get_out()) {
//                for (auto pair: a->_intersection_clique) {
//                    Concavity -=(R_lambda_var(pair->_name) - R_lambda_out(pair->_name).getvalue())*R_Xij_log[c](pair->_name).getvalue();
//                    Concavity -= (Im_lambda_var(pair->_name)- Im_lambda_out(pair->_name).getvalue())*Im_Xij_log[c](pair->_name).getvalue();
//                    Concavity -=  (lambda_var(pair->_name)- lambda_out(pair->_name).getvalue())*Xii_log[c](pair->_dest->_name).getvalue();
//                    DebugOn("lambda_out["<< c << "]("<<  pair->_name << ") = " << lambda_out(pair->_name).getvalue()<< endl);
//                    DebugOn("Xii_log[" << c << "]("<< pair->_dest->_name << ")= " << Xii_log[c](pair->_dest->_name).getvalue() << endl);
//                }
//            }
//
//            for (auto a: bag->get_in()) {
//                for (auto pair: a->_intersection_clique) {
//                    Concavity += (R_lambda_var(pair->_name) - R_lambda_out(pair->_name).getvalue())*R_Xij_log[c](pair->_name).getvalue();
//                    Concavity += (Im_lambda_var(pair->_name)- Im_lambda_out(pair->_name).getvalue())*Im_Xij_log[c](pair->_name).getvalue();
//                    Concavity += (lambda_var(pair->_name)- lambda_out(pair->_name).getvalue())*Xii_log[c](pair->_dest->_name).getvalue();
//                    DebugOn("lambda_out["<< c << "]("<<  pair->_name << ") = " << lambda_out(pair->_name).getvalue()<< endl);
//                    DebugOn("Xii_log[" << c << "]("<< pair->_dest->_name << ")= " << Xii_log[c](pair->_dest->_name).getvalue() << endl);
//                }
//            }
//            Master.add_constraint(Concavity <= 0);
//        }
//    }
//    solver solve_Master(Master, cplex);
//    solve_Master.run();
//    // Initialise the out point.
//    auto val = (*(var<Real>*) Master.get_var("gamma_C"));
//    for (const auto& bag: P.G_part.nodes) {
//        gamma_out(bag->_id) = val(bag->_id).getvalue();
//    }
//    auto R_val = (*(var<Real>*) Master.get_var("R_lambda_var"));
//    auto Im_val  = (*(var<Real>*)Master.get_var("Im_lambda_var"));
//    auto Val = (*(var<Real>*) Master.get_var("lambda_var"));
//    for (const auto& b: P.inter_pairs) {
//        R_lambda_out(b->_name) = R_val(b->_name).getvalue();
//        Im_lambda_out(b->_name) = Im_val(b->_name).getvalue();
//        lambda_out(b->_name) = Val(b->_name).getvalue();
//    }
//    cout << "value: " << Master._obj_val  <<endl;
//    UBlog[0] = Master._obj_val;
//////////////////////////////// BEGIN LAGRANGE ITERATIONS HERE /////////////////////////////////////
//    DebugOn("<<<<<<<<<<< Lagrange decomposition algorithm >>>>>>>>>"<< endl);
//    double epsilon = 0.005;
//    int itcount = 1;
//    int cout_inpoints = 0;
//    while ((UBlog[itcount-1] -LBlog[itcount-1] > epsilon*std::abs(LBlog[itcount-1])) && itcount < iter_limit) {
//        ////// CONSTRUCT SEPARATION POINTS
//        for (int c = 0; c < nbparts; c++) {
//            gamma_sep(c) = alpha*gamma_out(c).getvalue() + (1 - alpha)*gamma_in(c).getvalue();
//        }
//
//        for (const auto& a: P.inter_pairs) {
//            R_lambda_sep(a->_name) = alpha*R_lambda_out(a->_name).getvalue()
//                                     + (1 - alpha)*R_lambda_in(a->_name).getvalue();
//            Im_lambda_sep(a->_name) = alpha*Im_lambda_out(a->_name).getvalue()
//                                      + (1 - alpha)*Im_lambda_in(a->_name).getvalue();
//            lambda_sep(a->_name) = alpha*lambda_out(a->_name).getvalue()
//                                   + (1 - alpha)*lambda_in(a->_name).getvalue();
//        }
//        /////////////////SLOVE SUBPROBLEM/////////////////////
//        dual = 0;
//        for (int c = 0; c < nbparts; c++) {
//            value_dual[c] = subproblem(grid, P, c, Pg[c], Qg[c], Xii[c], R_Xij[c], Im_Xij[c],
//                                       R_lambda_sep, Im_lambda_sep, lambda_sep, Xii_log[c], R_Xij_log[c], Im_Xij_log[c]);
//            dual += value_dual[c] ;
//        }
//        DebugOn("dual: " << dual << endl);
//// UPDATE POINTS of Kelly using in-out algorithm (Ben-Ameur and Neto)
//        if (dual- sum(gamma_sep).eval()< 0) {
//            for (auto bag: P.G_part.nodes) {
//                unsigned c = bag->_id;
//                Constraint Concavity("Iter_" + to_string(itcount) + "_Concavity_" + to_string(c));
//                Concavity += gamma_C(c) - value_dual[c];
//                for (auto a: bag->get_out()) {
//                    for (auto pair: a->_intersection_clique) {
//                        Concavity -=(R_lambda_var(pair->_name) - R_lambda_sep(pair->_name).getvalue())*R_Xij_log[c](pair->_name).getvalue();
//                        Concavity -= (Im_lambda_var(pair->_name)- Im_lambda_sep(pair->_name).getvalue())*Im_Xij_log[c](pair->_name).getvalue();
//                        Concavity -=  (lambda_var(pair->_name)- lambda_sep(pair->_name).getvalue())*Xii_log[c](pair->_dest->_name).getvalue();
//                    }
//                }
//
//                for (auto a: bag->get_in()) {
//                    for (auto pair: a->_intersection_clique) {
//                        Concavity += (R_lambda_var(pair->_name) - R_lambda_sep(pair->_name).getvalue())*R_Xij_log[c](pair->_name).getvalue();
//                        Concavity += (Im_lambda_var(pair->_name)- Im_lambda_sep(pair->_name).getvalue())*Im_Xij_log[c](pair->_name).getvalue();
//                        Concavity += (lambda_var(pair->_name)- lambda_sep(pair->_name).getvalue())*Xii_log[c](pair->_dest->_name).getvalue();
//                    }
//                }
//                Master.add_constraint(Concavity <= 0);
//            }
//            if (dual > sum(gamma_in).eval()) {
//                cout_inpoints += 1;
//                for (int c = 0; c < nbparts; c++)
//                    gamma_in(c) = value_dual[c];
//                for (auto pair: P.inter_pairs) {
//                    R_lambda_in(pair->_name)  = R_lambda_sep(pair->_name).getvalue();
//                    Im_lambda_in(pair->_name) = Im_lambda_sep(pair->_name).getvalue();
//                    lambda_in(pair->_name) = lambda_sep(pair->_name).getvalue();
//                }
//            }
//            solver solve_master(Master, cplex);
//            solve_master.run();
//            DebugOn("master problem value: " << Master._obj_val << endl);
//            auto R_val = (*(var<Real>*) Master.get_var("R_lambda_var"));
//            auto Im_val  = (*(var<Real>*)Master.get_var("Im_lambda_var"));
//            auto Val = (*(var<Real>*) Master.get_var("lambda_var"));
//            for (auto b: P.inter_pairs) {
//                R_lambda_out(b->_name) = R_val(b->_name).getvalue();
//                Im_lambda_out(b->_name) = Im_val(b->_name).getvalue();
//                lambda_out(b->_name) = Val(b->_name).getvalue();
//            }
//            auto val = (*(var<Real>*) Master.get_var("gamma_C"));
//            for (const auto& bag: P.G_part.nodes) {
//                gamma_out(bag->_id) = val(bag->_id);
//            }
//        }
//        else {
//            cout_inpoints += 1;
//            for (int c = 0; c < nbparts; c++)
//                gamma_in(c) = value_dual[c];
//            for (const auto& pair: P.inter_pairs) {
//                R_lambda_in(pair->_name)  = R_lambda_sep(pair->_name).getvalue();
//                Im_lambda_in(pair->_name) = Im_lambda_sep(pair->_name).getvalue();
//                lambda_in(pair->_name) = lambda_sep(pair->_name).getvalue();
//            }
//        }
//        LDlog[itcount] = dual;
//        LBlog[itcount] = std::max(dual, LBlog[itcount-1]);
//        if (itcount > 1)
//            UBlog[itcount] = std::min(Master._obj_val, UBlog[itcount-1]);
//        else
//            UBlog[itcount] = Master._obj_val;
//        itcount +=1;
//    }
//    auto wall1 = get_wall_time();
//    auto cpu1 = get_cpu_time();
//    DebugOn("CPU time: " << cpu1-cpu0 << endl);
//    DebugOn("Wall time: " << wall1- wall0<< endl);
//
//    cout<< setw(15) << left <<"ITERATION" << setw(15) << "Current" << setw(15) << "LB" << setw(15)  << "UB" << endl;
//    for(int i = 0; i < itcount; i++) {
//        cout<< setw(15) << left <<i << setw(15) << LDlog[i]<<setw(15)<< LBlog[i] << setw(15) << UBlog[i] << endl;
//    }
//    cout << "updating in points " << cout_inpoints << "times. " << endl;
//    return 0;
//}

int main (int argc, const char * argv[])
{
    // Decompose
    const char* fname;
    double l = 0;
    unsigned iter_limit = 50;

    if (argc >= 2) {
        fname = argv[1];
        l = atof(argv[2]);
        iter_limit = (int) atof(argv[3]);
    }
    else {
        //fname = "../../data_sets/Power/nesta_case5_pjm.m";
        //fname = "../../data_sets/Power/nesta_case30_ieee.m";
        //fname = "../../data_sets/Power/nesta_case6_c.m";
        fname = "../../data_sets/Power/nesta_case5_pjm.m";
        //fname = "../../data_sets/Power/nesta_case3_lmbd.m";
        //fname = "../../data_sets/Power/nesta_case300_ieee.m";
        //fname = "../../data_sets/Power/nesta_case1354_pegase.m";
        //fname = "../../data_sets/Power/nesta_case118_ieee.m";
        //fname = "../../data_sets/Power/nesta_case57_ieee.m";
        //fname = "../../data_sets/Power/nesta_case14_ieee.m";
        //fname = "../../data_sets/Power/nesta_case57_ieee.m";
        l = 1;
    }
    PowerNet grid;
    grid.readgrid(fname);
    cout << "////////////////////////////////////////" << endl;
    reform_inout(grid, 1);
    //inout(grid, 2, 40);
    return 0;
}
