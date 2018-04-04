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
                     vector<var<Real>>& Pg, vector<var<Real>>& Qg, vector<var<bool>>& Start_up, vector<var<bool>>& Shut_down, vector<var<bool>>& On_off,
                     vector<var<Real>>& Xii, vector<var<Real>>& R_Xij, vector<var<Real>>& Im_Xij,
                     param<Real>& lambda_up, param<Real>& lambda_down,
                     param<Real>& zeta_up, param<Real>& zeta_down)
{
    // Grid Parameters
    const auto bus_pairs = grid.get_bus_pairs();
    Model ACUC("ACUC Model");
    for (int t = 0; t < T; t++) {
        ACUC.add_var(Pg[t].in_at(grid.gens,t));
        ACUC.add_var(Qg[t].in_at(grid.gens,t));
        ACUC.add_var(On_off[t].in_at(grid.gens, t));
        ACUC.add_var(Start_up[t].in_at(grid.gens, t));
        ACUC.add_var(Shut_down[t].in_at(grid.gens, t));
        ACUC.add_var(Xii[t].in_at(grid.nodes, t));
        ACUC.add_var(R_Xij[t].in_at(bus_pairs, t));
        ACUC.add_var(Im_Xij[t].in_at(bus_pairs, t));
    }
    //power flow vars are treated as auxiliary vars.
    var<Real> Pf_from("Pf_from", grid.S_max.in(grid.arcs, T));
    var<Real> Qf_from("Qf_from", grid.S_max.in(grid.arcs, T));
    var<Real> Pf_to("Pf_to", grid.S_max.in(grid.arcs, T));
    var<Real> Qf_to("Qf_to", grid.S_max.in(grid.arcs, T));
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
                    obj += grid.c1(name)*Pg[t](name)+ grid.c2(name)*Pg[t](name)*Pg[t](name) +grid.c0(name)*On_off[t](name);
                }
                else {
                    string name = g->_name + ","+ to_string(t);
                    obj += grid.c1(name)*Pg[t](name)+ grid.c2(name)*Pg[t](name)*Pg[t](name) + grid.c0(name)*On_off[t](name);
                    obj += cost_up.getvalue()*Start_up[t](name)+ cost_down.getvalue()*Shut_down[t](name);
                }
            }
        }
    }
    ACUC.set_objective(min(obj));
    for (int t= 0; t < T; t++) {
        Constraint SOC("SOC_" + to_string(t));
        SOC =  power(R_Xij[t], 2) + power(Im_Xij[t], 2) - Xii[t].from()*Xii[t].to() ;
        ACUC.add_constraint(SOC.in_at(bus_pairs, t) <= 0);
    }
    //KCL
    for (int t= 0; t < T; t++) {
        Constraint KCL_P("KCL_P_"+ to_string(t));
        Constraint KCL_Q("KCL_Q_"+ to_string(t));
        KCL_P =  sum(Pf_from.out_arcs()) + sum(Pf_to.in_arcs()) +grid.pl- sum(Pg[t].in_gens())+ grid.gs*Xii[t];
        ACUC.add_constraint(KCL_P.in_at(grid.nodes,t) == 0);

        KCL_Q  = sum(Qf_from.out_arcs()) + sum(Qf_to.in_arcs())+ grid.ql -sum(Qg[t].in_gens())-grid.bs*Xii[t];
        ACUC.add_constraint(KCL_Q.in_at(grid.nodes, t) == 0);

        Constraint Flow_P_From("Flow_P_From" + to_string(t));
        Flow_P_From = Pf_from- (grid.g_ff*Xii[t].from()+ grid.g_ft*R_Xij[t].in_pairs() + grid.b_ft*Im_Xij[t].in_pairs());
        ACUC.add_constraint(Flow_P_From.in_at(grid.arcs, t) == 0);

        Constraint Flow_P_To("Flow_P_To" + to_string(t));
        Flow_P_To = Pf_to - (grid.g_tt*Xii[t].to()+ grid.g_tf*R_Xij[t].in_pairs()- grid.b_tf*Im_Xij[t].in_pairs());
        ACUC.add_constraint(Flow_P_To.in_at(grid.arcs, t) == 0);

        Constraint Flow_Q_From("Flow_Q_From" + to_string(t));
        Flow_Q_From = Qf_from-(grid.g_ft*Im_Xij[t].in_pairs ()- grid.b_ff*Xii[t].from()- grid.b_ft*R_Xij[t].in_pairs());
        ACUC.add_constraint(Flow_Q_From.in_at(grid.arcs, t) == 0);

        Constraint Flow_Q_To("Flow_Q_To" + to_string(t));
        Flow_Q_To = Qf_to + (grid.b_tt*Xii[t].to()+ grid.b_tf*R_Xij[t].in_pairs() + grid.g_tf*Im_Xij[t].in_pairs());
        ACUC.add_constraint(Flow_Q_To.in_at(grid.arcs, t) == 0);

        Constraint Thermal_Limit_from("Thermal_Limit_from" + to_string(t));
        Thermal_Limit_from = power(Pf_from, 2) + power(Qf_from, 2);
        Thermal_Limit_from <= power(grid.S_max,2);
        ACUC.add_constraint(Thermal_Limit_from.in_at(grid.arcs, t));

        Constraint Thermal_Limit_to("Thermal_Limit_to" + to_string(t));
        Thermal_Limit_to = power(Pf_to, 2) + power(Qf_to, 2);
        Thermal_Limit_to <= power(grid.S_max, 2);
        ACUC.add_constraint(Thermal_Limit_to.in_at(grid.arcs, t));
    }
    ///* Phase Angle Bounds constraints */
    for (int t= 0; t < T; t++) {
        Constraint PAD_UB("PAD_UB_"+to_string(t));
        PAD_UB = Im_Xij[t]- grid.tan_th_max*R_Xij[t];
        ACUC.add_constraint(PAD_UB.in_at(bus_pairs, t) <= 0);

        Constraint PAD_LB("PAD_LB_"+to_string(t));
        PAD_LB = Im_Xij[t]- grid.tan_th_min*R_Xij[t];
        ACUC.add_constraint(PAD_LB.in_at(bus_pairs, t) >= 0);
    }
    // COMMITMENT CONSTRAINTS
    // Inter-temporal constraints 3a, 3d
    for (int t = 1; t < T; t++) {
        Constraint MC1("MC1_" + to_string(t));
        Constraint MC2("MC2_" + to_string(t));
        MC1 = On_off[t].in_at(grid.gens, t) -  On_off[t-1].in_at(grid.gens, t-1) -Start_up[t].in_at(grid.gens, t);
        MC2 = On_off[t-1].in_at(grid.gens, t-1) -  On_off[t].in_at(grid.gens, t) -Shut_down[t].in_at(grid.gens, t);
        ACUC.add_constraint(MC1 <= 0);
        ACUC.add_constraint(MC2 <= 0);
    }
    // Min-up constraints  4a
    for (int t = 1; t < T; t++) {
        Constraint Min_up1("Min_up1_"+to_string(t));
        Min_up1 = On_off[t].in_at(grid.gens, t) - On_off[t-1].in_at(grid.gens, t-1) - Start_up[t].in_at(grid.gens, t) + Shut_down[t].in_at(grid.gens, t);
        ACUC.add_constraint(Min_up1 == 0);
    }
    // 4b
    for (int t = min_up.getvalue(); t < T; t++) {
        Constraint Min_Up("Min_Up_constraint_"+ to_string(t));
        for (int l = t-min_up.getvalue()+1; l < t+1; l++) {
            Min_Up   += Start_up[t].in_at(grid.gens, l);
        }
        Min_Up -= On_off[t].in_at(grid.gens, t);
        ACUC.add_constraint(Min_Up <= 0);
    }
    // 4c
    for (int t = min_down.getvalue(); t < T; t++) {
        Constraint Min_Down("Min_Down_constraint_" + to_string(t));
        for (int l = t-min_down.getvalue()+1; l < t +1; l++) {
            Min_Down   += Shut_down[t].in_at(grid.gens, l);
        }
        Min_Down -= 1 - On_off[t].in_at(grid.gens, t);
        ACUC.add_constraint(Min_Down <= 0);
    }
    // Ramp Rate
    for (int t = 0; t < T; t++) {
        Constraint Production_P_LB("Production_P_LB_"+ to_string(t));
        Constraint Production_P_UB("Production_P_UB_"+ to_string(t));
        Constraint Production_Q_LB("Production_Q_LB_"+ to_string(t));
        Constraint Production_Q_UB("Production_Q_UB_"+ to_string(t));
        // 5A
        Production_P_UB = Pg[t]- grid.pg_max*On_off[t];
        Production_P_LB = Pg[t]- grid.pg_min*On_off[t];
        ACUC.add_constraint(Production_P_UB.in_at(grid.gens, t)<=0);
        ACUC.add_constraint(Production_P_LB.in_at(grid.gens, t)>= 0);

        Production_Q_UB = Qg[t] - grid.qg_max*On_off[t];
        Production_Q_LB = Qg[t] - grid.qg_min*On_off[t];
        ACUC.add_constraint(Production_Q_UB.in_at(grid.gens, t) <= 0);
        ACUC.add_constraint(Production_Q_LB.in_at(grid.gens, t) >= 0);
    }
    // 5C
    for (int t = 1; t < T; t++) {
        Constraint Ramp_up("Ramp_up_constraint_"  + to_string(t));
        Constraint Ramp_down("Ramp_down_constraint_"+ to_string(t));
        Ramp_up =  Pg[t].in_at(grid.gens, t);
        Ramp_up -= Pg[t-1].in_at(grid.gens, t-1);
        Ramp_up -= rate_ramp*On_off[t-1].in_at(grid.gens, t-1);
        Ramp_up -= rate_switch*(1 - On_off[t].in_at(grid.gens, t));

        Ramp_down =  Pg[t-1].in_at(grid.gens, t-1);
        Ramp_down -= Pg[t].in_at(grid.gens, t);
        Ramp_down -= rate_ramp*On_off[t].in_at(grid.gens, t);
        Ramp_down -= rate_switch*(1 - On_off[t-1].in_at(grid.gens, t-1));

        ACUC.add_constraint(Ramp_up <= 0);
        ACUC.add_constraint(Ramp_down <= 0);
    }
    // set the initial state of generators.
    Constraint gen_initial("initial_state");
    gen_initial +=  On_off[0].in_at(grid.gens, 0);
    ACUC.add_constraint(gen_initial == 1);
    /* Solver selection */
    solver cpx_acuc(ACUC, cplex);
    //solver cpx_acuc(ACUC, ipopt);
    bool relax =true;
    int output = 1;
    double tol = 1e-6;
    cpx_acuc.run(output, relax, tol);
    cout << "the continuous relaxation bound is: " << ACUC._obj_val << endl;
    for (int t = 1; t < T; t++) {
        auto MC1 = ACUC.get_constraint("MC1_" + to_string(t));
        auto MC2 = ACUC.get_constraint("MC2_" + to_string(t));
        auto Ramp_up = ACUC.get_constraint("Ramp_up_constraint_"  + to_string(t));
        auto Ramp_down = ACUC.get_constraint("Ramp_down_constraint_"  + to_string(t));
        int i = 0;
        for (auto& g: grid.gens) {
            string name = g->_name + "," + to_string(t);
            lambda_up(name) = abs(MC1->_dual.at(i));
            lambda_down(name) = abs(MC2->_dual.at(i));
            zeta_up(name) = abs(Ramp_up->_dual.at(i));
            zeta_down(name) = abs(Ramp_down->_dual.at(i));
            DebugOn("dual of  lambda_up " << name << " " << abs(MC1->_dual[i]) << endl);
            DebugOn("dual of  lambda_down " << name << " " << abs(MC2->_dual[i]) << endl);
            DebugOn("dual of  zeta_up " << name << " " << abs(Ramp_up->_dual[i]) << endl);
            DebugOn("dual of  zeta_down " << name << " " << abs(Ramp_down->_dual[i]) << endl);
            ++i;
        }
    }
    return ACUC._obj_val;
}

double subproblem(PowerNet& grid,  unsigned t, unsigned T, param<Real>& rate_ramp, param<Real>& rate_switch,
                  const param<Real>& min_up, const param<Real>& min_down,
                  const param<Real>& cost_up, const param<Real>& cost_down,
                  var<Real>& Pg, var<Real>& Qg, var<bool>& Start_up, var<bool>& Shut_down, var<bool>& On_off,
                  var<Real>& Xii, var<Real>& R_Xij, var<Real>& Im_Xij,
                  param<Real>& lambda_up, param<Real>& lambda_down, param<Real>& zeta_up, param<Real>& zeta_down)
{
    //Grid Parameters
    const auto bus_pairs = grid.get_bus_pairs();
    Model Sub("Sub"+to_string(t));
    Sub.add_var(Pg.in_at(grid.gens,t));
    Pg.in_at(grid.gens, t).print(true);
    cout << "dim: " << Pg.in_at(grid.gens, t).get_dim() << endl;;
    Sub.add_var(Qg.in_at(grid.gens,t));
    Sub.add_var(On_off.in_at(grid.gens, t));
    Sub.add_var(Start_up.in_at(grid.gens, t));
    Sub.add_var(Shut_down.in_at(grid.gens, t));
    Sub.add_var(Xii.in_at(grid.nodes, t));
    Sub.add_var(R_Xij.in_at(bus_pairs, t));
    Sub.add_var(Im_Xij.in_at(bus_pairs, t));
    //power flow
    var<Real> Pf_from("Pf_from", grid.S_max.in_at(grid.arcs, t));
    var<Real> Qf_from("Qf_from", grid.S_max.in_at(grid.arcs, t));
    var<Real> Pf_to("Pf_to", grid.S_max.in_at(grid.arcs, t));
    var<Real> Qf_to("Qf_to", grid.S_max.in_at(grid.arcs, t));
    Sub.add_var(Pf_from.in_at(grid.arcs, t));
    Sub.add_var(Qf_from.in_at(grid.arcs, t));
    Sub.add_var(Pf_to.in_at(grid.arcs, t));
    Sub.add_var(Qf_to.in_at(grid.arcs, t));
    /* Construct the objective function*/
    func_ obj;
    for (auto g:grid.gens) {
        if (g->_active) {
            if (t == 0) {
                string name = g->_name + ",0";
                string name1 = g->_name + ",1";
                obj += (grid.c1(name) + zeta_down(name1) - zeta_up(name1))*Pg(name) + grid.c2(name)*Pg(name)*Pg(name);
                obj +=(grid.c0(name) + lambda_down(name1) - lambda_down(name1) + zeta_down(name1)*rate_ramp(name1)
                      - zeta_up(name1)*rate_switch(name1))*On_off(name);
            }
            else if (t == T-1){
                string name = g->_name + ","+ to_string(t);
                obj += (grid.c1(name) + zeta_up(name)- zeta_down(name))*Pg(name)
                + grid.c2(name)*Pg(name)*Pg(name);
                
                obj += (grid.c0(name)+lambda_up(name) -lambda_down(name)
                        - zeta_down(name)*rate_ramp(name) + zeta_down(name)*rate_switch(name))*On_off(name);
                
                obj += (cost_up.getvalue()-lambda_up(name))*Start_up(name);
                obj += (cost_down.getvalue()-lambda_down(name))*Shut_down(name);
            }
            else{
                string name = g->_name + ","+ to_string(t);
                string name1 = g->_name + ","+ to_string(t+1);
                obj += (grid.c1(name) + zeta_up(name)+zeta_down(name1)- zeta_down(name)-zeta_up(name1))*Pg(name)
                + grid.c2(name)*Pg(name)*Pg(name);
                
                obj += (grid.c0(name)+lambda_up(name) -lambda_up(name1) + lambda_down(name1) -lambda_down(name)
                    - zeta_down(name)*rate_ramp(name) - zeta_up(name1)*rate_ramp(name1)
                    + zeta_up(name1)*rate_switch(name1)+ zeta_down(name)*rate_switch(name))*On_off(name);
                
                obj += (cost_up.getvalue()-lambda_up(name))*Start_up(name);
                obj += (cost_down.getvalue()-lambda_down(name))*Shut_down(name);
            }
        }
    }
    Sub.set_objective(min(obj));
    Constraint SOC("SOC_" + to_string(t));
    SOC =  power(R_Xij, 2) + power(Im_Xij, 2) - Xii.from()*Xii.to() ;
    Sub.add_constraint(SOC.in_at(bus_pairs, t) <= 0);

    Constraint KCL_P("KCL_P_"+ to_string(t));
    Constraint KCL_Q("KCL_Q_"+ to_string(t));
    KCL_P =  sum(Pf_from.out_arcs()) + sum(Pf_to.in_arcs()) +grid.pl- sum(Pg.in_gens())+ grid.gs*Xii;
    Sub.add_constraint(KCL_P.in_at(grid.nodes,t) == 0);

    KCL_Q  = sum(Qf_from.out_arcs()) + sum(Qf_to.in_arcs())+ grid.ql -sum(Qg.in_gens())-grid.bs*Xii;
    Sub.add_constraint(KCL_Q.in_at(grid.nodes, t) == 0);

    Constraint Flow_P_From("Flow_P_From" + to_string(t));
    Flow_P_From = Pf_from- (grid.g_ff*Xii.from()+ grid.g_ft*R_Xij.in_pairs() + grid.b_ft*Im_Xij.in_pairs());
    Sub.add_constraint(Flow_P_From.in_at(grid.arcs, t) == 0);

    Constraint Flow_P_To("Flow_P_To" + to_string(t));
    Flow_P_To = Pf_to - (grid.g_tt*Xii.to()+ grid.g_tf*R_Xij.in_pairs()- grid.b_tf*Im_Xij.in_pairs());
    Sub.add_constraint(Flow_P_To.in_at(grid.arcs, t) == 0);

    Constraint Flow_Q_From("Flow_Q_From" + to_string(t));
    Flow_Q_From = Qf_from-(grid.g_ft*Im_Xij.in_pairs ()- grid.b_ff*Xii.from()- grid.b_ft*R_Xij.in_pairs());
    Sub.add_constraint(Flow_Q_From.in_at(grid.arcs, t) == 0);

    Constraint Flow_Q_To("Flow_Q_To" + to_string(t));
    Flow_Q_To = Qf_to + (grid.b_tt*Xii.to()+ grid.b_tf*R_Xij.in_pairs() + grid.g_tf*Im_Xij.in_pairs());
    Sub.add_constraint(Flow_Q_To.in_at(grid.arcs, t) == 0);

    Constraint Thermal_Limit_from("Thermal_Limit_from" + to_string(t));
    Thermal_Limit_from = power(Pf_from, 2) + power(Qf_from, 2);
    Thermal_Limit_from <= power(grid.S_max,2);
    Sub.add_constraint(Thermal_Limit_from.in_at(grid.arcs, t));

    Constraint Thermal_Limit_to("Thermal_Limit_to" + to_string(t));
    Thermal_Limit_to = power(Pf_to, 2) + power(Qf_to, 2);
    Thermal_Limit_to <= power(grid.S_max, 2);
    Sub.add_constraint(Thermal_Limit_to.in_at(grid.arcs, t));

    Constraint PAD_UB("PAD_UB_"+to_string(t));
    PAD_UB = Im_Xij- grid.tan_th_max*R_Xij;
    Sub.add_constraint(PAD_UB.in_at(bus_pairs, t) <= 0);

    Constraint PAD_LB("PAD_LB_"+to_string(t));
    PAD_LB = Im_Xij- grid.tan_th_min*R_Xij;
    Sub.add_constraint(PAD_LB.in_at(bus_pairs, t) >= 0);

    Constraint Production_P_LB("Production_P_LB_"+ to_string(t));
    Constraint Production_P_UB("Production_P_UB_"+ to_string(t));
    Constraint Production_Q_LB("Production_Q_LB_"+ to_string(t));
    Constraint Production_Q_UB("Production_Q_UB_"+ to_string(t));
    Production_P_UB = Pg- grid.pg_max*On_off;
    Production_P_LB = Pg- grid.pg_min*On_off;
    Sub.add_constraint(Production_P_UB.in_at(grid.gens, t)<=0);
    Sub.add_constraint(Production_P_LB.in_at(grid.gens, t)>= 0);

    Production_Q_UB = Qg - grid.qg_max*On_off;
    Production_Q_LB = Qg - grid.qg_min*On_off;
    Sub.add_constraint(Production_Q_UB.in_at(grid.gens, t) <= 0);
    Sub.add_constraint(Production_Q_LB.in_at(grid.gens, t) >= 0);

    if (t  == 0) {
        Constraint gen_initial("initial_state");
        gen_initial +=  On_off.in_at(grid.gens, 0);
        Sub.add_constraint(gen_initial == 1);
    }
    /* Solver selection */
    solver cpx_acuc(Sub, cplex);
    bool relax = false;
    int output = 1;
    double tol = 1e-6;
    cpx_acuc.run(output, relax, tol);
    return Sub._obj_val;
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
        //fname = "../../data_sets/Power/nesta_case5_pjm.m";
        //fname = "../../data_sets/Power/nesta_case3_lmbd.m";
        //fname = "../../data_sets/Power/nesta_case300_ieee.m";
        //fname = "../../data_sets/Power/nesta_case1354_pegase.m";
        fname = "../../data_sets/Power/nesta_case14_ieee.m";
        //fname = "../../data_sets/Power/nesta_case57_ieee.m";
        l = 1;
    }
    PowerNet grid;
    grid.readgrid(fname);
    //GRAPH PARTITION
    auto bus_pairs = grid.get_bus_pairs();
    // Schedule Parameters
    unsigned T = 2;
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
    //var<Real> Pg("Pg", grid.pg_min.in(grid.gens, T), grid.pg_max.in(grid.gens, T));
    //var<Real> Qg ("Qg", grid.qg_min.in(grid.gens, T), grid.qg_max.in(grid.gens, T));
    vector<var<Real>> Pg;
    vector<var<Real>> Qg;
    for (int t = 0; t < T; t++) {
        var<Real> pgt("Pg" + to_string(t), grid.pg_min.in_at(grid.gens, t), grid.pg_max.in_at(grid.gens, t));
        var<Real> qgt("Qg" + to_string(t), grid.pg_min.in_at(grid.gens, t), grid.pg_max.in_at(grid.gens, t));
        Pg.push_back(pgt);
        Qg.push_back(qgt);
    }

    //Lifted variables.
    vector<var<Real>> R_Xij;
    vector<var<Real>> Im_Xij;
    vector<var<Real>> Xii;
    for (int t = 0; t < T; t++) {
        var<Real>  R_Xijt("R_Wij" + to_string(t), grid.wr_min, grid.wr_max); // real part of Wij
        var<Real>  Im_Xijt("Im_Wij" + to_string(t), grid.wi_min, grid.wi_max); // imaginary part of Wij.
        var<Real>  Xiit("Wii" + to_string(t), grid.w_min, grid.w_max);
        R_Xijt.initialize_all(1.0);
        Xiit.initialize_all(1.001);
        R_Xij.push_back(R_Xijt);
        Im_Xij.push_back(Im_Xijt);
        Xii.push_back(Xiit);
    }
    // Commitment variables
    vector<var<bool>> On_off;
    vector<var<bool>> Start_up;
    vector<var<bool>> Shut_down;
    for (int t = 0; t < T; t++) {
        var<bool>  On_offt("On_off" + to_string(t));
        var<bool>  Start_upt("Start_up" + to_string(t));
        var<bool>  Shut_downt("Shut_down" + to_string(t));
        On_off.push_back(On_offt);
        Start_up.push_back(Start_upt);
        Shut_down.push_back(Shut_downt);
    }
/////////////////// DEFINE LAGRANGE MULTIPLIERS  ////////////////////////////////
    param<Real> lambda_up("lambda_up");
    param<Real> lambda_down("lambda_down");
    param<Real> zeta_up("zeta_up");
    param<Real> zeta_down("zeta_down");
    lambda_up.in(grid.gens, T);
    lambda_down.in(grid.gens, T);
    zeta_up.in(grid.gens, T);
    zeta_down.in(grid.gens, T);

    lambda_up.initialize_all(0);
    lambda_down.initialize_all(0);
    zeta_up.initialize_all(0);
    zeta_down.initialize_all(0);
    double lb_cts = getdual_relax(grid, T, rate_ramp, rate_switch, min_up, min_down, cost_up, cost_down, Pg, Qg,
                                  Start_up, Shut_down, On_off, Xii, R_Xij,  Im_Xij, lambda_up, lambda_down, zeta_up,
                                  zeta_down);
    //Improve the lower bound using MISCOP.
    std::vector<double> Subs;
    Subs.resize(T);
    double LB = 0;
    for(int t = 0; t < T; t++) {
        Subs[t]= subproblem(grid, t, T, rate_ramp, rate_switch, min_up, min_down, cost_up, cost_down, Pg[t], Qg[t],
                            Start_up[t], Shut_down[t], On_off[t], Xii[t], R_Xij[t], Im_Xij[t], lambda_up, lambda_down, zeta_up, zeta_down);
        LB += Subs[t];
    }
    for (int t = 1; t < T; t++){
        for (auto& g: grid.gens) {
            string name = g->_name + "," + to_string(t);
            LB -= zeta_down(name).getvalue()*rate_switch(name).getvalue();
        }
    }
    for (int t = 0; t < T-1; t++){
        for (auto& g: grid.gens) {
            string name = g->_name + "," + to_string(t);
            LB -= zeta_up(name).getvalue()*rate_switch(name).getvalue();
        }
    }
    cout << "The initial Lower bound of the ACUC problem is: " << LB << endl;
    // now we need to solve it faster
    return 0;
}
