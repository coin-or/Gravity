
//
//  Global.cpp
//  projects
//
//  Created by Guanglei on 5/3/18.
//
//
#include "Global.hpp"
#include <queue>
#include <algorithm>
Global::Global() {
    rate_ramp.set_name("rate_ramp");
    rate_switch.set_name("rate_switch");
    min_up.set_name("min_up");
    min_down.set_name("min_down");
    cost_up.set_name("cost_up");
    cost_down.set_name("cost_down");
    min_up = 2;
    min_down = 2;
    cost_up = 50;
    cost_down = 30;
};
Global::~Global() {
    delete grid;
};

Global::Global(PowerNet* net, int parts, int T) {
    grid = net;
    P_ = new Partition();
    Num_parts = parts;
    if (Num_parts >1 && Num_parts < grid->nodes.size()) {
        P_->get_ncut(*grid, Num_parts);
    }
    else {
        throw std::invalid_argument("please input a valid partition size (>=1 and <= graph size)!");
    }

    R_lambda_.set_name("R_lambda");
    Im_lambda_.set_name("Im_lambda");
    lambda_.set_name("lambda");
    R_lambda_.initialize_all(0);
    Im_lambda_.initialize_all(0);
    lambda_.initialize_all(0);

    Num_time  =  T;
    rate_ramp.set_name("rate_ramp");
    rate_switch.set_name("rate_switch");
    min_up.set_name("min_up");
    min_down.set_name("min_down");
    cost_up.set_name("cost_up");
    cost_down.set_name("cost_down");
    On_off_initial.set_name("On_off_initial");
    Pg_initial.set_name("Pg_initial");
    //parameters
    min_up = 2;
    min_down = 2;
    cost_up = 50;
    cost_down = 30;
    for (auto g: grid->gens) {
        rate_ramp(g->_name) = max(grid->pg_min(g->_name).getvalue(), 1*grid->pg_max(g->_name).getvalue());
        rate_switch(g->_name) = max(grid->pg_min(g->_name).getvalue(), 1*grid->pg_max(g->_name).getvalue());
    }
    min_up = 2;
    min_down = 2;
    cost_up = 50;
    cost_down = 30;
    grid->time_expand(T);
    rate_ramp.time_expand(T);
    rate_switch.time_expand(T);
    // set initial state or read the initial state
    for (auto g: grid->gens) {
        Pg_initial(g->_name) = 0;
        On_off_initial(g->_name) = 0;
    }

    auto bus_pairs = grid->get_bus_pairs();
    for (int t = 0; t < T; t++) {
        var<Real> pgt("Pg_" + to_string(t), grid->pg_min, grid->pg_max);
        //var<Real> pgt("Pg" + to_string(t), grid->pg_min.in_at(grid->gens, t), grid->pg_max.in_at(grid->gens, t));
        var<Real> qgt("Qg" + to_string(t), grid->qg_min.in_at(grid->gens, t), grid->qg_max.in_at(grid->gens, t));
        var<Real> pg2t("Pg2_" + to_string(t), non_neg_); // new var introduced for the perspective formulation.
        Pg.push_back(pgt);
        Pg2.push_back(pg2t);
        Qg.push_back(qgt);
    }
    //Lifted variables.
    for (int t = 0; t < T; t++) {
        var<Real>  R_Xijt("R_Wij" + to_string(t), grid->wr_min.in_at(bus_pairs, t), grid->wr_max.in_at(bus_pairs, t)); // real part of Wij
        var<Real>  Im_Xijt("Im_Wij" + to_string(t), grid->wi_min.in_at(bus_pairs, t), grid->wi_max.in_at(bus_pairs, t));
        var<Real>  Xiit("Wii" + to_string(t), grid->w_min.in_at(grid->nodes,t), grid->w_max.in_at(grid->nodes,t));
        R_Xijt.initialize_all(1.0);
        Xiit.initialize_all(1.001);
        R_Xij.push_back(R_Xijt);
        Im_Xij.push_back(Im_Xijt);
        Xii.push_back(Xiit);
    }
    // Commitment variables
    for (int t = 0; t < T; t++) {
        var<bool>  On_offt("On_off_" + to_string(t-1));
        var<bool>  Start_upt("Start_up_" + to_string(t));
        var<bool>  Shut_downt("Shut_down_" + to_string(t));
        On_off.push_back(On_offt);
        Start_up.push_back(Start_upt);
        Shut_down.push_back(Shut_downt);
    }
    var<bool>  On_offT("On_off_" + to_string(T-1));
    On_off.push_back(On_offT);
    On_off_sol_.resize(T);
    Pg_sol_.resize(T);
    Start_up_sol_.resize(T);
    Shut_down_sol_.resize(T);
    // multipliers.
    lambda_up.set_name("lambda_up");
    lambda_down.set_name("lambda_down");
    zeta_up.set_name("zeta_up");
    zeta_down.set_name("zeta_down");
    mu.set_name("mu");
    mu_up.set_name("mu_up");
    mu_down.set_name("mu_down");

    lambda_up.in(grid->gens, T);
    lambda_down.in(grid->gens, T);
    zeta_up.in(grid->gens, T);
    zeta_down.in(grid->gens, T);

    lambda_down.initialize_all(0);
    zeta_up.initialize_all(0);
    zeta_down.initialize_all(0);
    mu.in(grid->gens, T);
    mu_up.in(grid->gens, T);
    mu_down.in(grid->gens, T);
    Sub_.resize(T);
}

double Global::getdual_relax_time_(bool include) {
    include_min_updown_ = include;
    const auto bus_pairs = grid->get_bus_pairs();
    Model ACUC("ACUC Model");
    for (int t = 0; t < Num_time; t++) {
        ACUC.add_var(Pg[t].in_at(grid->gens, t));
        ACUC.add_var(Qg[t].in_at(grid->gens,t));
        ACUC.add_var(On_off[t].in_at(grid->gens, t-1));
        ACUC.add_var(Start_up[t].in_at(grid->gens, t));
        ACUC.add_var(Shut_down[t].in_at(grid->gens, t));
        ACUC.add_var(Xii[t].in_at(grid->nodes,  t));
        ACUC.add_var(R_Xij[t].in_at(bus_pairs, t));
        ACUC.add_var(Im_Xij[t].in_at(bus_pairs,t));
        Xii[t].initialize_all(1.001);
        R_Xij[t].initialize_all(1.0);
    }
    ACUC.add_var(On_off[Num_time].in_at(grid->gens, Num_time-1));
    //power flow vars are treated as auxiliary vars.
    var<Real> Pf_from("Pf_from", grid->S_max.in(grid->arcs, Num_time));
    var<Real> Qf_from("Qf_from", grid->S_max.in(grid->arcs, Num_time));
    var<Real> Pf_to("Pf_to", grid->S_max.in(grid->arcs, Num_time));
    var<Real> Qf_to("Qf_to", grid->S_max.in(grid->arcs, Num_time));
    ACUC.add_var(Pf_from.in(grid->arcs, Num_time));
    ACUC.add_var(Qf_from.in(grid->arcs, Num_time));
    ACUC.add_var(Pf_to.in(grid->arcs, Num_time));
    ACUC.add_var(Qf_to.in(grid->arcs, Num_time));
    /* Construct the objective function*/
    func_ obj;
    for (int t = 0; t < Num_time; t++) {
        for (auto g:grid->gens) {
            if (g->_active) {
                string name = g->_name + ","+ to_string(t);
                obj += grid->c1(name)*Pg[t](name)+ grid->c2(name)*Pg[t](name)*Pg[t](name) + grid->c0(name)*On_off[t+1](name);
                obj += cost_up.getvalue()*Start_up[t](name)+ cost_down.getvalue()*Shut_down[t](name);
            }
        }
    }
    ACUC.set_objective(min(obj));
    for (int t= 0; t < Num_time; t++) {
        Constraint SOC("SOC_" + to_string(t));
        SOC =  power(R_Xij[t], 2) + power(Im_Xij[t], 2) - Xii[t].from()*Xii[t].to() ;
        ACUC.add_constraint(SOC.in_at(bus_pairs, t) <= 0);
    }
    //KCL
    for (int t= 0; t < Num_time; t++) {
        Constraint KCL_P("KCL_P_"+ to_string(t));
        Constraint KCL_Q("KCL_Q_"+ to_string(t));
        KCL_P =  sum(Pf_from.out_arcs()) + sum(Pf_to.in_arcs()) +grid->pl- sum(Pg[t].in_gens())+ grid->gs*Xii[t];
        ACUC.add_constraint(KCL_P.in_at(grid->nodes,t) == 0);

        KCL_Q  = sum(Qf_from.out_arcs()) + sum(Qf_to.in_arcs())+ grid->ql -sum(Qg[t].in_gens())-grid->bs*Xii[t];
        ACUC.add_constraint(KCL_Q.in_at(grid->nodes, t) == 0);

        Constraint Flow_P_From("Flow_P_From" + to_string(t));
        Flow_P_From = Pf_from- (grid->g_ff*Xii[t].from()+ grid->g_ft*R_Xij[t].in_pairs() + grid->b_ft*Im_Xij[t].in_pairs());
        ACUC.add_constraint(Flow_P_From.in_at(grid->arcs, t) == 0);

        Constraint Flow_P_To("Flow_P_To" + to_string(t));
        Flow_P_To = Pf_to - (grid->g_tt*Xii[t].to()+ grid->g_tf*R_Xij[t].in_pairs()- grid->b_tf*Im_Xij[t].in_pairs());
        ACUC.add_constraint(Flow_P_To.in_at(grid->arcs, t) == 0);

        Constraint Flow_Q_From("Flow_Q_From" + to_string(t));
        Flow_Q_From = Qf_from-(grid->g_ft*Im_Xij[t].in_pairs ()- grid->b_ff*Xii[t].from()- grid->b_ft*R_Xij[t].in_pairs());
        ACUC.add_constraint(Flow_Q_From.in_at(grid->arcs, t) == 0);

        Constraint Flow_Q_To("Flow_Q_To" + to_string(t));
        Flow_Q_To = Qf_to + (grid->b_tt*Xii[t].to()+ grid->b_tf*R_Xij[t].in_pairs() + grid->g_tf*Im_Xij[t].in_pairs());
        ACUC.add_constraint(Flow_Q_To.in_at(grid->arcs, t) == 0);

        Constraint Thermal_Limit_from("Thermal_Limit_from" + to_string(t));
        Thermal_Limit_from = power(Pf_from, 2) + power(Qf_from, 2);
        Thermal_Limit_from <= power(grid->S_max,2);
        ACUC.add_constraint(Thermal_Limit_from.in_at(grid->arcs, t));

        Constraint Thermal_Limit_to("Thermal_Limit_to" + to_string(t));
        Thermal_Limit_to = power(Pf_to, 2) + power(Qf_to, 2);
        Thermal_Limit_to <= power(grid->S_max, 2);
        ACUC.add_constraint(Thermal_Limit_to.in_at(grid->arcs, t));
    }
    ///* Phase Angle Bounds constraints */
    for (int t= 0; t < Num_time; t++) {
        Constraint PAD_UB("PAD_UB_"+to_string(t));
        PAD_UB = Im_Xij[t]- grid->tan_th_max*R_Xij[t];
        ACUC.add_constraint(PAD_UB.in_at(bus_pairs, t) <= 0);

        Constraint PAD_LB("PAD_LB_"+to_string(t));
        PAD_LB = Im_Xij[t]- grid->tan_th_min*R_Xij[t];
        ACUC.add_constraint(PAD_LB.in_at(bus_pairs, t) >= 0);
    }
    // COMMITMENT CONSTRAINTS
    for (int t = 0; t < Num_time; t++) {
        Constraint Production_P_LB("Production_P_LB_"+ to_string(t));
        Constraint Production_P_UB("Production_P_UB_"+ to_string(t));
        Constraint Production_Q_LB("Production_Q_LB_"+ to_string(t));
        Constraint Production_Q_UB("Production_Q_UB_"+ to_string(t));
        Production_P_UB = Pg[t]- grid->pg_max*On_off[t+1];
        Production_P_LB = Pg[t]- grid->pg_min*On_off[t+1];
        ACUC.add_constraint(Production_P_UB.in_at(grid->gens, t)<=0);
        ACUC.add_constraint(Production_P_LB.in_at(grid->gens, t)>= 0);

        Production_Q_UB = Qg[t] - grid->qg_max*On_off[t+1];
        Production_Q_LB = Qg[t] - grid->qg_min*On_off[t+1];
        ACUC.add_constraint(Production_Q_UB.in_at(grid->gens, t) <= 0);
        ACUC.add_constraint(Production_Q_LB.in_at(grid->gens, t) >= 0);
    }

    for (int t = 0; t < Num_time; t++) {
        for (auto& g: grid->gens) {
            if (g->_active) {
                Constraint MC1("Inter_temporal_MC1_" + to_string(t)+ ","+ g->_name);
                Constraint MC2("Inter_temporal_MC2_" + to_string(t)+ ","+ g->_name);
                string name = g->_name +"," + to_string(t);
                string name1 = g->_name +"," + to_string(t-1);
                MC1 = On_off[t+1](name) -  On_off[t](name1) -Start_up[t](name);
                MC2 = On_off[t](name1) -  On_off[t+1](name) -Shut_down[t](name);
                ACUC.add_constraint(MC1 <= 0);
                ACUC.add_constraint(MC2 <= 0);
            }
        }
//        Constraint MC1("Inter_temporal_MC1_" + to_string(t));
//        Constraint MC2("Inter_temporal_MC2_" + to_string(t));
//        MC1 = On_off[t+1].in_at(grid.gens, t) - On_off[t].in_at(grid.gens, t-1) - Start_up[t].in_at(grid.gens, t);
//        MC2 = On_off[t].in_at(grid.gens, t-1) - On_off[t+1].in_at(grid.gens, t) - Shut_down[t].in_at(grid.gens, t);
//        ACUC.add_constraint(MC1 <= 0);
//        ACUC.add_constraint(MC2 <= 0);
    }
    for (int t = 0; t < Num_time; t++) {
        for (auto& g: grid->gens) {
            if (g->_active) {
                Constraint OnOffStartupShutdown("OnOffStartupShutdown_"+ to_string(t) + ","+ g->_name);
                string name = g->_name +"," + to_string(t);
                string name1 = g->_name +"," + to_string(t-1);
                OnOffStartupShutdown = On_off[t+1](name) - On_off[t](name1)
                                       - Start_up[t](name) + Shut_down[t](name);
                ACUC.add_constraint(OnOffStartupShutdown == 0);
            }
        }
    }
    //for (auto& g: grid->gens) {
    //    if (g->_active) {
//        Constraint OnOffStartupShutdown("OnOffStartupShutdown_"+ to_string(t));
//        OnOffStartupShutdown = On_off[t+1].in_at(grid.gens, t) - On_off[t].in_at(grid.gens, t-1)
//                               - Start_up[t].in_at(grid.gens, t) + Shut_down[t].in_at(grid.gens, t);
//        ACUC.add_constraint(OnOffStartupShutdown == 0);
    //    }
    //}

    if (include_min_updown_) {
        // Min-up constraints  4b
        for (int t = min_up.getvalue()-1; t < Num_time; t++) {
            Constraint Min_Up("Min_Up_constraint_"+ to_string(t));
            for (int l = t-min_up.getvalue()+1; l < t+1; l++) {
                Min_Up   += Start_up[l].in_at(grid->gens, l);
            }
            Min_Up -= On_off[t+1].in_at(grid->gens, t);
            ACUC.add_constraint(Min_Up <= 0);
        }
        // 4c
        for (int t = min_down.getvalue()-1; t < Num_time; t++) {
            Constraint Min_Down("Min_Down_constraint_" + to_string(t));
            for (int l = t-min_down.getvalue()+1; l < t +1; l++) {
                Min_Down   += Shut_down[l].in_at(grid->gens, l);
            }
            Min_Down -= 1 - On_off[t+1].in_at(grid->gens, t);
            ACUC.add_constraint(Min_Down <= 0);
        }
    }
    for (int t = 1; t < Num_time; t++) {
        for (auto& g: grid->gens) {
            if (g->_active) {
                Constraint Ramp_up("Ramp_up_constraint_"  + to_string(t) + "," + g->_name);
                Constraint Ramp_down("Ramp_down_constraint_"+ to_string(t) + "," + g->_name);
                string name = g->_name +"," + to_string(t);
                string name1 = g->_name +"," + to_string(t-1);
                Ramp_up =  Pg[t](name) - Pg[t-1](name1) -  rate_ramp.getvalue()*On_off[t](name1) - rate_switch.getvalue()*(1 - On_off[t](name1));
                Ramp_down =  Pg[t-1](name1) - Pg[t](name) - rate_ramp.getvalue()*On_off[t+1](name)- rate_switch.getvalue()*(1 - On_off[t+1](name));
                ACUC.add_constraint(Ramp_up <= 0);
                ACUC.add_constraint(Ramp_down <= 0);
            }
        }
    }
    // t =0, we have ramp up constraint.
    for (auto& g: grid->gens) {
        Constraint Ramp_up("Ramp_up_constraint_0," + g->_name);
        string name = g->_name +",0";
        Ramp_up =  Pg[0](name) - Pg_initial(g->_name)
                   - rate_ramp(name)*On_off_initial(g->_name);
        Ramp_up -= rate_switch(name)*(1 - On_off_initial(g->_name));
        ACUC.add_constraint(Ramp_up <= 0);

        Constraint Ramp_down("Ramp_down_constraint_0," + g->_name);
        Ramp_down =   -1*Pg[0](name) + Pg_initial.in(g->_name);
        Ramp_down -= rate_ramp(name)*On_off[1](name);
        Ramp_down -= rate_switch(name)*(1 - On_off[1](name));
        ACUC.add_constraint(Ramp_down <= 0);
    }
    // set the initial state of generators.
    Constraint gen_initial("initial_state");
    gen_initial +=  On_off[0].in_at(grid->gens, -1) - On_off_initial.in(grid->nodes);
    ACUC.add_constraint(gen_initial == 0);
    /* Solver selection */
    bool relax = true;
    int output = 1;
    solver cpx_acuc(ACUC, cplex);
    double tol = 10e-6;
    cpx_acuc.run(output, relax, tol);
    cout << "the continuous relaxation bound is: " << ACUC._obj_val << endl;
    for (int t = 1; t < Num_time; t++) {
        for (auto& g: grid->gens) {
            if (g->_active) {
                auto MC1 = ACUC.get_constraint("Inter_temporal_MC1_" + to_string(t)+","+ g->_name);
                auto MC2 = ACUC.get_constraint("Inter_temporal_MC2_" + to_string(t)+","+ g->_name);
                string name = g->_name + "," + to_string(t);
                lambda_up(name) = abs(MC1->_dual.at(0));
                lambda_down(name) = abs(MC2->_dual.at(0));
                DebugOn("dual of  lambda_up_" << name << " " << abs(MC1->_dual[0]) << endl);
                DebugOn("dual of  lambda_down_" << name << " " << abs(MC2->_dual[0]) << endl);
            }
        }
    }
    // we do not relax first step constraint
    for (int t = 1; t < Num_time; t++) {
        for (auto& g: grid->gens) {
            if (g->_active) {
                auto Ramp_up = ACUC.get_constraint("Ramp_up_constraint_" + to_string(t) + "," + g->_name);
                string name = g->_name +","+ to_string(t);
                zeta_up(name) = abs(Ramp_up->_dual.at(0));

                auto Ramp_down = ACUC.get_constraint("Ramp_down_constraint_"  + to_string(t)+"," + g->_name);
                zeta_down(name) = abs(Ramp_down->_dual.at(0));

                DebugOff("dual of  zeta_up_" << name << " " << abs(Ramp_up->_dual[0]) << endl);
                DebugOff("dual of  zeta_down_" << name << " " << abs(Ramp_down->_dual[0]) << endl);
            }
        }
    }
    // we do not relax first step constraint
    for (int t= 1; t < Num_time; t++) {
        for (auto& g: grid->gens) {
            if (g->_active) {
                auto OOSS  = ACUC.get_constraint("OnOffStartupShutdown_"+ to_string(t) + ","+ g->_name);
                string name = g->_name + "," + to_string(t);
                cout << "val: " << to_string(OOSS->eval(0)) << endl;
                DebugOff("mu: " << -OOSS->_dual.at(0)  << endl);
                mu(name) = -OOSS->_dual.at(0); // when positive and when negative?
            }
        }
    }

    if (include_min_updown_) {
        if (min_up.getvalue() >1) {
            for (int t = min_up.getvalue()-1; t < Num_time; t++) {
                auto Min_up = ACUC.get_constraint("Min_Up_constraint_"+ to_string(t));
                auto Min_down = ACUC.get_constraint("Min_Down_constraint_"+ to_string(t));
                int l=0;
                for (auto& g: grid->gens) {
                    if (g->_active) {
                        string name = g->_name + "," + to_string(t);
                        mu_up(name) = abs(Min_up->_dual.at(l));
                        mu_down(name) = abs(Min_down->_dual.at(l));
                        l += 1;
                    }
                }
            }
        }
    }
    return ACUC._obj_val;
}

double Global::getdual_relax_spatial() {
    R_Xij.clear();
    Im_Xij.clear();
    Xii.clear();
    Pg.clear();
    Pg2.clear();
    Qg.clear();
    On_off.clear();
    Start_up.clear();
    Shut_down.clear();

    for (int c = 0; c < Num_parts; c++) {
        var<Real>  bag_Xii("Xii_"+ to_string(c), grid->w_min.in(P_->bag_bus_union_out[c], Num_time),
                           grid->w_max.in(P_->bag_bus_union_out[c], Num_time));
        Xii.push_back(bag_Xii);
        var<Real>  bag_R_Xij("R_Xij_"+ to_string(c), grid->wr_min.in(P_->bag_bus_pairs_union[c], Num_time),
                             grid->wr_max.in(P_->bag_bus_pairs_union[c], Num_time));
        var<Real>  bag_Im_Xij("Im_Xij_"+ to_string(c), grid->wi_min.in(P_->bag_bus_pairs_union[c], Num_time),
                              grid->wi_max.in(P_->bag_bus_pairs_union[c], Num_time));
        R_Xij.push_back(bag_R_Xij);
        Im_Xij.push_back(bag_Im_Xij);

        var<Real>  bag_Pg("Pg_" + to_string(c), grid->pg_min.in(P_->bag_gens[c], Num_time), grid->pg_max.in(P_->bag_gens[c],Num_time));
        var<Real>  bag_Qg("Qg_" + to_string(c), grid->qg_min.in(P_->bag_gens[c], Num_time), grid->qg_max.in(P_->bag_gens[c], Num_time));
        var<Real> bag_Pg2("Pg2_" + to_string(c), non_neg_); // new var introduced for the perspective formulation.
        Pg.push_back(bag_Pg);
        Pg2.push_back(bag_Pg2.in(P_->bag_gens[c], Num_time));
        Qg.push_back(bag_Qg);

        var<bool>  bag_Onoff("On_off_" + to_string(c));
        var<bool>  bag_Up("Start_up_" + to_string(c));
        var<bool>  bag_Down("Shut_down_" + to_string(c));

        On_off.push_back(bag_Onoff);
        Start_up.push_back(bag_Up);
        Shut_down.push_back(bag_Down);
    }

    Model ACUC("ACUC Model");
    for (int c = 0; c < Num_parts; c++) {
        ACUC.add_var(Pg[c].in(P_->bag_gens[c], Num_time));
        ACUC.add_var(Pg2[c].in(P_->bag_gens[c], Num_time));
        ACUC.add_var(Qg[c].in(P_->bag_gens[c], Num_time));
        ACUC.add_var(Start_up[c].in(P_->bag_gens[c], Num_time));
        ACUC.add_var(Shut_down[c].in(P_->bag_gens[c], Num_time));
        ACUC.add_var(On_off[c].in(P_->bag_gens[c], -1,  Num_time));
        ACUC.add_var(Xii[c].in(P_->bag_bus_union_out[c], Num_time));
        Xii[c].initialize_all(1.001);
        ACUC.add_var(R_Xij[c].in(P_->bag_bus_pairs_union[c], Num_time));
        R_Xij[c].initialize_all(1.0);
        ACUC.add_var(Im_Xij[c].in(P_->bag_bus_pairs_union[c], Num_time));
    }
    //power flow vars are treated as auxiliary vars.
    vector<var<Real>> Pf_from;
    vector<var<Real>> Qf_from;
    vector<var<Real>> Pf_to;
    vector<var<Real>> Qf_to;
    for (int c = 0; c < Num_parts; c++) {
        if (P_->bag_arcs_union[c].size() > 0) {
            var<Real> bag_Pf_from("Pf_from"+to_string(c), grid->S_max.in(P_->bag_arcs_union_out[c], Num_time));
            var<Real> bag_Qf_from("Qf_from"+to_string(c), grid->S_max.in(P_->bag_arcs_union_out[c], Num_time));
            var<Real> bag_Pf_to("Pf_to"+to_string(c), grid->S_max.in(P_->bag_arcs_union_in[c], Num_time));
            var<Real> bag_Qf_to("Qf_to"+to_string(c), grid->S_max.in(P_->bag_arcs_union_in[c], Num_time));
            Pf_from.push_back(bag_Pf_from);
            Qf_from.push_back(bag_Qf_from);
            Pf_to.push_back(bag_Pf_to);
            Qf_to.push_back(bag_Qf_to);
            ACUC.add_var(bag_Pf_from.in(P_->bag_arcs_union_out[c], Num_time));
            ACUC.add_var(bag_Qf_from.in(P_->bag_arcs_union_out[c], Num_time));
            ACUC.add_var(bag_Pf_to.in(P_->bag_arcs_union_in[c], Num_time));
            ACUC.add_var(bag_Qf_to.in(P_->bag_arcs_union_in[c], Num_time));
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
    for (int c = 0; c < Num_parts; c++) {
        for (auto g: P_->bag_gens[c]) {
            for (int t = 0; t < Num_time; t++) {
                if (g->_active) {
                    string name = g->_name + ","+ to_string(t);
                    obj += grid->c1(name)*Pg[c](name)+ grid->c2(name)*Pg[c](name)*Pg[c](name) +grid->c0(name)*On_off[c](name);
                    //obj += grid.c1(name)*Pg[c](name)+ grid.c2(name)*Pg2[c](name) +grid.c0(name)*On_off[c](name);
                    obj += cost_up*Start_up[c](name)+ cost_down*Shut_down[c](name);
                }
            }
        }
    }
    ACUC.set_objective(min(obj));
    for (int c = 0; c < Num_parts; c++) {
        if (P_->bag_bus_pairs_union_directed[c].size() > 0) {
            Constraint SOC("SOC_" + to_string(c));
            SOC =  power(R_Xij[c], 2)+ power(Im_Xij[c], 2) - Xii[c].from()*Xii[c].to() ;
            ACUC.add_constraint(SOC.in(P_->bag_bus_pairs_union_directed[c], Num_time) <= 0);
        }
    }
    /* perspective formulation of Pg^2 */
    //for (int c = 0; c < nbparts; c++) {
    //    if (P_->bag_gens[c].size() >0){
    //        Constraint Perspective_OnOff_Pg2("Perspective_OnOff_Pg2_" + to_string(c));
    //        Perspective_OnOff_Pg2 = power(Pg[c], 2) - Pg2[c]*On_off[c];
    //        ACUC.add(Perspective_OnOff_Pg2.in(P_->bag_gens[c], Num_time) <= 0);
    //    }
    //}
    //KCL
    for (int c = 0; c < Num_parts; c++) {
        Constraint KCL_P("KCL_P"+ to_string(c));
        Constraint KCL_Q("KCL_Q"+ to_string(c));
        KCL_P = sum(Pf_from[c].out_arcs()) + sum(Pf_to[c].in_arcs()) -sum(Pg[c].in_gens())+ grid->gs*Xii[c] + grid->pl;
        ACUC.add_constraint(KCL_P.in(P_->bag_bus[c], Num_time) == 0);
        KCL_Q  = sum(Qf_from[c].out_arcs()) + sum(Qf_to[c].in_arcs())+ grid->ql -sum(Qg[c].in_gens())-grid->bs*Xii[c];
        ACUC.add_constraint(KCL_Q.in(P_->bag_bus[c], Num_time) == 0);

        Constraint Flow_P_From("Flow_P_From_" + to_string(c));
        Flow_P_From = Pf_from[c]- (grid->g_ff*Xii[c].from()+ grid->g_ft*R_Xij[c].in_pairs() + grid->b_ft*Im_Xij[c].in_pairs());
        ACUC.add_constraint(Flow_P_From.in(P_->bag_arcs_union_out[c], Num_time) == 0);

        Constraint Flow_P_Num_timeo("Flow_P_Num_timeo" + to_string(c));
        Flow_P_Num_timeo = Pf_to[c] - (grid->g_tt*Xii[c].to() + grid->g_tf*R_Xij[c].in_pairs() - grid->b_tf*Im_Xij[c].in_pairs());
        ACUC.add_constraint(Flow_P_Num_timeo.in(P_->bag_arcs_union_in[c], Num_time) == 0);

        Constraint Flow_Q_From("Flow_Q_From" + to_string(c));
        Flow_Q_From = Qf_from[c] - (grid->g_ft*Im_Xij[c].in_pairs() - grid->b_ff*Xii[c].from() - grid->b_ft*R_Xij[c].in_pairs());
        ACUC.add_constraint(Flow_Q_From.in(P_->bag_arcs_union_out[c], Num_time) == 0);

        Constraint Flow_Q_Num_timeo("Flow_Q_Num_timeo" + to_string(c));
        Flow_Q_Num_timeo = Qf_to[c] + (grid->b_tt*Xii[c].to() + grid->b_tf*R_Xij[c].in_pairs() + grid->g_tf*Im_Xij[c].in_pairs());
        ACUC.add_constraint(Flow_Q_Num_timeo.in(P_->bag_arcs_union_in[c], Num_time) == 0);

        Constraint Num_timehermal_Limit_from("Num_timehermal_Limit_from" + to_string(c));
        Num_timehermal_Limit_from = power(Pf_from[c], 2) + power(Qf_from[c], 2);
        Num_timehermal_Limit_from <= power(grid->S_max,2);
        ACUC.add_constraint(Num_timehermal_Limit_from.in(P_->bag_arcs_union_out[c], Num_time));

        Constraint Num_timehermal_Limit_to("Num_timehermal_Limit_to" + to_string(c));
        Num_timehermal_Limit_to = power(Pf_to[c], 2) + power(Qf_to[c], 2);
        Num_timehermal_Limit_to <= power(grid->S_max,2);
        ACUC.add_constraint(Num_timehermal_Limit_to.in(P_->bag_arcs_union_in[c], Num_time));
    }
    ///* Phase Angle Bounds constraints */
    for (int c = 0; c < Num_parts; c++) {
        if (P_->bag_bus_pairs_union_directed[c].size() > 0) {
            Constraint PAD_UB("PAD_UB_"+to_string(c));
            PAD_UB = Im_Xij[c]- grid->tan_th_max*R_Xij[c];
            ACUC.add_constraint(PAD_UB.in(P_->bag_bus_pairs_union_directed[c], Num_time) <= 0);

            Constraint PAD_LB("PAD_LB_"+to_string(c));
            PAD_LB = Im_Xij[c]- grid->tan_th_min*R_Xij[c];
            ACUC.add_constraint(PAD_LB.in(P_->bag_bus_pairs_union_directed[c], Num_time) >= 0);
        }
    }
    // COMMINum_timeMENNum_time CONSNum_timeRAINNum_timeS
    // Inter-temporal constraints 3a, 3d
    for (int c = 0; c < Num_parts; c++) {
        if (P_->bag_gens[c].size() > 0) {
            for (int t = 0; t < Num_time; t++) {
                Constraint MC1("MC1_"+ to_string(c)+ ","+to_string(t));
                Constraint MC2("MC2_"+ to_string(c)+ ","+to_string(t));
                MC1 = On_off[c].in_at(P_->bag_gens[c], t)- On_off[c].in_at(P_->bag_gens[c], t-1)-  Start_up[c].in_at(P_->bag_gens[c], t);
                MC2 = On_off[c].in_at(P_->bag_gens[c], t-1) - On_off[c].in_at(P_->bag_gens[c], t) - Shut_down[c].in_at(P_->bag_gens[c], t);
                ACUC.add_constraint(MC1 <= 0);
                ACUC.add_constraint(MC2 <= 0);
            }
        }
    }
    for (int c = 0; c < Num_parts; c++) {
        if (P_->bag_gens[c].size() > 0) {
            for (int t = 0; t < Num_time; t++) {
                Constraint OnOffStartupShutdown("OnOffStartupShutdown_"+ to_string(t));
                OnOffStartupShutdown = On_off[c].in_at(P_->bag_gens[c], t) - On_off[c].in_at(P_->bag_gens[c], t-1)
                                       - Start_up[c].in_at(P_->bag_gens[c], t) + Shut_down[c].in_at(P_->bag_gens[c], t);
                ACUC.add_constraint(OnOffStartupShutdown == 0);
            }
        }
    }
    // 4b
    for (int c = 0; c < Num_parts; c++) {
        if (P_->bag_gens[c].size() > 0) {
            for (int t = min_up.getvalue()-1; t < Num_time; t++) {
                Constraint Min_Up("Min_Up_constraint" + to_string(c) + "_"+ to_string(t));
                for (int l = t-min_up.getvalue()+1; l < t+1; l++) {
                    Min_Up   += Start_up[c].in_at(P_->bag_gens[c], l);
                }
                Min_Up -= On_off[c].in_at(P_->bag_gens[c], t);
                ACUC.add_constraint(Min_Up <= 0);
            }
        }
    }
    // 4c
    for (int c = 0; c < Num_parts; c++) {
        if (P_->bag_gens[c].size() > 0) {
            for (int t = min_down.getvalue()-1; t < Num_time; t++) {
                Constraint Min_Down("Min_Down_constraint_" + to_string(c) + "_"+ to_string(t));
                for (int l = t-min_down.getvalue()+1; l < t +1; l++) {
                    Min_Down   += Shut_down[c].in_at(P_->bag_gens[c], l);
                }
                Min_Down -= 1 - On_off[c].in_at(P_->bag_gens[c], t);
                ACUC.add_constraint(Min_Down <= 0);
            }
        }
    }
    for (int c = 0; c < Num_parts; c++) {
        if (P_->bag_gens[c].size() > 0) {
            Constraint Production_P_LB("Production_P_LB_"+ to_string(c));
            Constraint Production_P_UB("Production_P_UB_"+ to_string(c));
            Constraint Production_Q_LB("Production_Q_LB_"+ to_string(c));
            Constraint Production_Q_UB("Production_Q_UB_"+ to_string(c));
            // 5A
            Production_P_UB = Pg[c] - grid->pg_max*On_off[c];
            Production_P_LB = Pg[c] - grid->pg_min*On_off[c];
            ACUC.add_constraint(Production_P_UB.in(P_->bag_gens[c], Num_time) <=0);
            ACUC.add_constraint(Production_P_LB.in(P_->bag_gens[c], Num_time) >= 0);

            Production_Q_UB = Qg[c] - grid->qg_max*On_off[c];
            Production_Q_LB = Qg[c] - grid->qg_min*On_off[c];
            ACUC.add_constraint(Production_Q_UB.in(P_->bag_gens[c], Num_time) <= 0);
            ACUC.add_constraint(Production_Q_LB.in(P_->bag_gens[c], Num_time) >= 0);
        }
    }
    // 5C
    for (int c = 0; c < Num_parts; c++) {
        if (P_->bag_gens[c].size() > 0) {
            for (int t = 1; t < Num_time; t++) {
                Constraint Ramp_up("Ramp_up_constraint_" +  to_string(c) + "_" + to_string(t));
                Constraint Ramp_down("Ramp_down_constraint_" + to_string(c) + "_" +to_string(t));
                Ramp_up =  Pg[c].in_at(P_->bag_gens[c], t);
                Ramp_up -= Pg[c].in_at(P_->bag_gens[c], t-1);
                Ramp_up -= rate_ramp.in_at(P_->bag_gens[c], t)*On_off[c].in_at(P_->bag_gens[c], t-1);
                Ramp_up -= rate_switch.in_at(P_->bag_gens[c], t)*(1 - On_off[c].in_at(P_->bag_gens[c], t));

                Ramp_down =  Pg[c].in_at(P_->bag_gens[c], t-1);
                Ramp_down -= Pg[c].in_at(P_->bag_gens[c], t);
                Ramp_down -= rate_ramp.in_at(P_->bag_gens[c], t)*On_off[c].in_at(P_->bag_gens[c], t);
                Ramp_down -= rate_switch.in_at(P_->bag_gens[c], t)*(1 - On_off[c].in_at(P_->bag_gens[c], t-1));

                ACUC.add_constraint(Ramp_up <= 0);
                ACUC.add_constraint(Ramp_down <= 0);
            }
            Constraint Ramp_up("Ramp_up_constraint0" + to_string(c));
            Ramp_up =  Pg[c].in_at(P_->bag_gens[c], 0) - Pg_initial.in(P_->bag_gens[c])
                       - rate_ramp.in_at(P_->bag_gens[c], 0)*On_off_initial.in(P_->bag_gens[c]);
            Ramp_up -= rate_switch.in_at(P_->bag_gens[c], 0)*(1 - On_off_initial.in(P_->bag_gens[c]));
            ACUC.add_constraint(Ramp_up <= 0);

            Constraint Ramp_down("Ramp_down_constraint0");
            Ramp_down =   -1*Pg[c].in_at(P_->bag_gens[c],0) + Pg_initial.in(P_->bag_gens[c]);
            Ramp_down -= rate_ramp.in_at(P_->bag_gens[c], 0)*On_off[c].in_at(P_->bag_gens[c], 0);
            Ramp_down -= rate_switch.in_at(P_->bag_gens[c], 0)*(1 - On_off[c].in_at(P_->bag_gens[c], 0));
            ACUC.add_constraint(Ramp_down <= 0);
        }
    }
    // set the initial state of generators.
    for (int c = 0; c < Num_parts; c++) {
        if (P_->bag_gens[c].size() >0) {
            Constraint gen_initial("gen_initial_"+to_string(c));
            gen_initial +=  On_off[c].in_at(P_->bag_gens[c], -1) - On_off_initial.in(P_->bag_gens[c]);
            ACUC.add_constraint(gen_initial == 0);
        }
    }
    // Linking constraints
    for (const auto& a: P_->G_part.arcs) {
        Constraint Link_R("link_R_"+a->_name);
        Link_R = R_Xij[a->_src->_id].in_pairs() - R_Xij[a->_dest->_id].in_pairs();
        ACUC.add_constraint(Link_R.in(a->_intersection_clique, Num_time) ==0);

        Constraint Link_Im("link_Im_"+a->_name);
        Link_Im = Im_Xij[a->_src->_id].in_pairs() - Im_Xij[a->_dest->_id].in_pairs();
        ACUC.add_constraint(Link_Im.in(a->_intersection_clique, Num_time)==0);

        Constraint Link_Xii("link_Xii_" + a->_name);
        Link_Xii = Xii[a->_src->_id].to() - Xii[a->_dest->_id].to();
        ACUC.add_constraint(Link_Xii.in(a->_intersection_clique, Num_time)==0);
    }
    /* Solver selection */
    solver cpx_acuc(ACUC, cplex);
    //solver cpx_acuc(ACUC, ipopt);
    bool relax =true;
    int output = 1;
    double tol = 1e-6;
    cpx_acuc.run(output, relax, tol);
    cout << "the continuous relaxation bound is: " << ACUC._obj_val << endl;
    for (const auto& a: P_->G_part.arcs) {
        auto consR = ACUC.get_constraint("link_R_"+a->_name);
        auto consIm = ACUC.get_constraint("link_Im_"+a->_name);
        auto cons = ACUC.get_constraint("link_Xii_"+a->_name);
        for (unsigned t = 0; t < Num_time; t++) {
            for (auto& line: a->_intersection_clique) {
                string name =line->_name+","+to_string(t);
                auto cR = (*consR)(name);
                auto cIm = (*consIm)(name);
                auto  c= (*cons)(name);
                R_lambda_(name) = -cR._dual.at(0);
                Im_lambda_(name)= -cIm._dual.at(0);
                lambda_(name) = -c._dual.at(0);
                DebugOff("R_lambda_" << -cR._dual.at(0) << endl);
                DebugOff("Im_lambda_" << -cIm._dual.at(0)<< endl);
                DebugOff("lambda_" << -c._dual.at(0)<< endl);
            }
        }
    }
    return ACUC._obj_val;
}


void Global::add_var_Sub_time(Model& Sub, int t) {
    const auto bus_pairs = grid->get_bus_pairs();
    Sub.add_var(Pg[t]);
    //Sub.add_var(Pg[t].in_at(grid->gens,t));
    Sub.add_var(Pg2[t].in_at(grid->gens,t));
    Sub.add_var(Qg[t].in_at(grid->gens,t));
    //Sub.add_var(On_off[t+1].in_at(grid->gens, t));// On_off has ranges from -1 to 1.
    Sub.add_var(On_off[t+1]);// On_off has ranges from -1 to 1.
    Sub.add_var(Start_up[t]);
    //Sub.add_var(Start_up[t].in_at(grid->gens, t));
    Sub.add_var(Shut_down[t]);
    //Sub.add_var(Shut_down[t].in_at(grid->gens, t));
    Sub.add_var(Xii[t].in_at(grid->nodes, t));
    Sub.add_var(R_Xij[t].in_at(bus_pairs, t));
    Sub.add_var(Im_Xij[t].in_at(bus_pairs, t));
    Xii[t].initialize_all(1.001);
    R_Xij[t].initialize_all(1.0);
    //power flow
    Pf_from = var<Real>("Pf_from", grid->S_max.in_at(grid->arcs, t));
    Qf_from = var<Real>("Qf_from", grid->S_max.in_at(grid->arcs, t));
    Pf_to = var<Real>("Pf_to", grid->S_max.in_at(grid->arcs, t));
    Qf_to = var<Real>("Qf_to", grid->S_max.in_at(grid->arcs, t));
    Sub.add_var(Pf_from.in_at(grid->arcs, t));
    Sub.add_var(Qf_from.in_at(grid->arcs, t));
    Sub.add_var(Pf_to.in_at(grid->arcs, t));
    Sub.add_var(Qf_to.in_at(grid->arcs, t));
}

void Global::add_obj_Sub_time(gravity::Model& Sub, int t) {
    // great
    func_ obj;
    if (t == 0) {
        for (auto g:grid->gens) {
            if (g->_active) {
                string name = g->_name + ",0";
                string name1 = g->_name + ",1";
                obj += (grid->c1(name) + zeta_down(name1) - zeta_up(name1))*Pg[t](name) + grid->c2(name)*Pg2[t](name);
                obj +=(grid->c0(name) + lambda_down(name1) - lambda_up(name1) + zeta_up(name1)*rate_switch(name1)
                       - zeta_up(name1)*rate_ramp(name1)-mu(name1))*On_off[t+1](name);
                obj += cost_up.getvalue()*Start_up[t](name) + cost_down.getvalue()*Shut_down[t](name);
                if (include_min_updown_) {
                    string name2 = g->_name +"," + to_string(min_up.getvalue()-1);
                    if (min_up.getvalue() >1) {
                        obj += mu_up(name2)*Start_up[t](name);
                    }
                    if(min_down.getvalue()>1) {
                        obj += mu_down(name2)*Shut_down[t](name);
                    }
                }
            }
        }
    }
    else if(t == Num_time-1) {
        for (auto g:grid->gens) {
            if (g->_active) {
                string name = g->_name + ","+ to_string(t);
                obj += (grid->c1(name) + zeta_up(name)- zeta_down(name))*Pg[t](name);
                obj += (grid->c0(name)+lambda_up(name) -lambda_down(name)
                        - zeta_down(name)*rate_ramp(name) + zeta_down(name)*rate_switch(name)
                        +mu(name))*On_off[t+1](name);
                if (include_min_updown_) {
                    obj += (mu_down(name)- mu_up(name))*On_off[t+1](name);
                }

                obj += (cost_up-lambda_up(name) - mu(name))*Start_up[t](name);
                obj += (cost_down-lambda_down(name)+mu(name))*Shut_down[t](name);
                if (include_min_updown_) {
                    if (min_up.getvalue() > 1) {
                        obj += mu_up(name)*Start_up[t](name);
                    }
                    if(min_down.getvalue()> 1) {
                        obj += mu_down(name)*Shut_down[t](name);
                    }
                }
            }
        }
    }
    else {
        for (auto g:grid->gens) {
            if (g->_active) {
                string name = g->_name + ","+ to_string(t);
                string name1 = g->_name + ","+ to_string(t+1);
                obj += (grid->c1(name) + zeta_up(name)+zeta_down(name1)- zeta_down(name)-zeta_up(name1))*Pg[t](name)
                       + grid->c2(name)*Pg[t](name)*Pg[t](name);

                obj += (grid->c0(name)+lambda_up(name) -lambda_up(name1) + lambda_down(name1) -lambda_down(name)
                        - zeta_down(name)*rate_ramp(name) - zeta_up(name1)*rate_ramp(name1)
                        + zeta_up(name1)*rate_switch(name1)+ zeta_down(name)*rate_switch(name)
                        + mu(name) - mu(name1))*On_off[t+1](name);

                if (include_min_updown_) {
                    if (t >= min_up.getvalue() -1) {
                        obj -= mu_up(name)*On_off[t+1](name);
                    }
                    if (t >= min_down.getvalue() -1) {
                        obj += mu_down(name)*On_off[t+1](name);
                    }
                }
                obj += (cost_up.getvalue()-lambda_up(name)-mu(name))*Start_up[t](name);
                obj += (cost_down.getvalue()-lambda_down(name)+mu(name))*Shut_down[t](name);
                if (include_min_updown_) {
                    if (min_up.getvalue() >1) {
                        int start = std::max(min_up.getvalue()-1, t);
                        int end = std::min(min_up.getvalue()+t, Num_time); // both increase 1.
                        for (int l = start; l < end; l++) {
                            string name2 = g->_name+","+to_string(l);
                            obj +=mu_up(name2)*Start_up[t](name);
                        }
                    }
                    if(min_down.getvalue()>1) {
                        int start = std::max(min_down.getvalue()-1, t);
                        int end = std::min(min_down.getvalue()+t, Num_time);
                        for (int l = start; l < end; l++) {
                            string name2 = g->_name+","+to_string(l);
                            obj +=mu_down(name2)*Shut_down[t](name);
                        }
                    }
                }
            }
        }
    }
    Sub.set_objective(min(obj));
}

void Global::add_obj_Sub_upper_time(gravity::Model& Sub, int t) {
    // great
    func_ obj;
    for (auto g:grid->gens) {
        if (g->_active) {
            string name = g->_name + ","+ to_string(t);
            obj += grid->c1(name)*Pg[t](name)+ grid->c2(name)*Pg[t](name)*Pg[t](name) + grid->c0(name)*On_off[t+1](name);
            obj += cost_up.getvalue()*Start_up[t](name)+ cost_down.getvalue()*Shut_down[t](name);
        }
    }
    Sub.set_objective(min(obj));
}

void Global::add_perspective_OnOff_Sub_time(Model& Sub, int t) {
    /* Construct the objective function*/
    Constraint Perspective_OnOff_Pg2("Perspective_OnOff_Pg2_");
    Perspective_OnOff_Pg2 = power(Pg[t], 2) - Pg2[t]*On_off[t+1];
    Sub.add(Perspective_OnOff_Pg2.in_at(grid->gens, t) <= 0);
}
void Global::add_SOCP_Sub_time(Model& Sub, int t) {
    const auto bus_pairs = grid->get_bus_pairs();
    Constraint SOC("SOC_" + to_string(t));
    SOC =  power(R_Xij[t], 2) + power(Im_Xij[t], 2) - Xii[t].from()*Xii[t].to() ;
    Sub.add_constraint(SOC.in_at(bus_pairs, t) <= 0);
}
void Global::add_KCL_Sub_time(Model& Sub, int t) {
    Constraint KCL_P("KCL_P_"+ to_string(t));
    Constraint KCL_Q("KCL_Q_"+ to_string(t));
    KCL_P =  sum(Pf_from.out_arcs()) + sum(Pf_to.in_arcs()) +grid->pl- sum(Pg[t].in_gens())+ grid->gs*Xii[t];
    Sub.add_constraint(KCL_P.in_at(grid->nodes,t) == 0);

    KCL_Q  = sum(Qf_from.out_arcs()) + sum(Qf_to.in_arcs())+ grid->ql -sum(Qg[t].in_gens())-grid->bs*Xii[t];
    Sub.add_constraint(KCL_Q.in_at(grid->nodes, t) == 0);

    Constraint Flow_P_From("Flow_P_From" + to_string(t));
    Flow_P_From = Pf_from- (grid->g_ff*Xii[t].from()+ grid->g_ft*R_Xij[t].in_pairs() + grid->b_ft*Im_Xij[t].in_pairs());
    Sub.add_constraint(Flow_P_From.in_at(grid->arcs, t) == 0);

    Constraint Flow_P_To("Flow_P_To" + to_string(t));
    Flow_P_To = Pf_to - (grid->g_tt*Xii[t].to()+ grid->g_tf*R_Xij[t].in_pairs()- grid->b_tf*Im_Xij[t].in_pairs());
    Sub.add_constraint(Flow_P_To.in_at(grid->arcs, t) == 0);

    Constraint Flow_Q_From("Flow_Q_From" + to_string(t));
    Flow_Q_From = Qf_from-(grid->g_ft*Im_Xij[t].in_pairs ()- grid->b_ff*Xii[t].from()- grid->b_ft*R_Xij[t].in_pairs());
    Sub.add_constraint(Flow_Q_From.in_at(grid->arcs, t) == 0);

    Constraint Flow_Q_To("Flow_Q_To" + to_string(t));
    Flow_Q_To = Qf_to + (grid->b_tt*Xii[t].to()+ grid->b_tf*R_Xij[t].in_pairs() + grid->g_tf*Im_Xij[t].in_pairs());
    Sub.add_constraint(Flow_Q_To.in_at(grid->arcs, t) == 0);

}
void Global::add_thermal_Sub_time(Model& Sub, int t) {
    const auto bus_pairs = grid->get_bus_pairs();
    Constraint Thermal_Limit_from("Thermal_Limit_from" + to_string(t));
    Thermal_Limit_from = power(Pf_from, 2) + power(Qf_from, 2);
    Thermal_Limit_from <= power(grid->S_max,2);
    Sub.add_constraint(Thermal_Limit_from.in_at(grid->arcs, t));

    Constraint Thermal_Limit_to("Thermal_Limit_to" + to_string(t));
    Thermal_Limit_to = power(Pf_to, 2) + power(Qf_to, 2);
    Thermal_Limit_to <= power(grid->S_max, 2);
    Sub.add_constraint(Thermal_Limit_to.in_at(grid->arcs, t));

    Constraint PAD_UB("PAD_UB_"+to_string(t));
    PAD_UB = Im_Xij[t]- grid->tan_th_max*R_Xij[t];
    Sub.add_constraint(PAD_UB.in_at(bus_pairs, t) <= 0);

    Constraint PAD_LB("PAD_LB_"+to_string(t));
    PAD_LB = Im_Xij[t]- grid->tan_th_min*R_Xij[t];
    Sub.add_constraint(PAD_LB.in_at(bus_pairs, t) >= 0);

    Constraint Production_P_LB("Production_P_LB_"+ to_string(t));
    Constraint Production_P_UB("Production_P_UB_"+ to_string(t));
    Constraint Production_Q_LB("Production_Q_LB_"+ to_string(t));
    Constraint Production_Q_UB("Production_Q_UB_"+ to_string(t));
    Production_P_UB = Pg[t]- grid->pg_max*On_off[t+1];
    Production_P_LB = Pg[t]- grid->pg_min*On_off[t+1];
    Sub.add_constraint(Production_P_UB.in_at(grid->gens, t)<=0);
    Sub.add_constraint(Production_P_LB.in_at(grid->gens, t)>= 0);

    Production_Q_UB = Qg[t] - grid->qg_max*On_off[t+1];
    Production_Q_LB = Qg[t] - grid->qg_min*On_off[t+1];
    Sub.add_constraint(Production_Q_UB.in_at(grid->gens, t) <= 0);
    Sub.add_constraint(Production_Q_LB.in_at(grid->gens, t) >= 0);
}

void Global::add_MC_upper_Sub_time(Model& Sub, int t) {
    Constraint MC_upper1("MC_upper1_constraint_"+ to_string(t));
    param<bool> On_off_val("On_off_val");
    //On_off_val = On_off[t+1];
    MC_upper1  = Start_up[t] - On_off[t+1];
    Sub.add_constraint(MC_upper1.in_at(grid->gens, t)<=0);

    Constraint MC_upper2("MC_upper2_constraint_"+ to_string(t));
    MC_upper2  = Shut_down[t] -1 + On_off[t+1];
    Sub.add_constraint(MC_upper1.in_at(grid->gens, t)<=0);
}

void  Global::add_MC_intertemporal_Sub_upper_time(Model& ACUC, int t) {
    for (auto& g: grid->gens) {
        if (g->_active) {
            Constraint MC1("Inter_temporal_MC1_" + to_string(t)+ ","+ g->_name);
            Constraint MC2("Inter_temporal_MC2_" + to_string(t)+ ","+ g->_name);
            if (t >0) {
                string name = g->_name +"," + to_string(t);
                string name1 = g->_name +"," + to_string(t-1);
                MC1 = On_off[t+1](name) - On_off_sol_[t-1](name1) -Start_up[t](name);
                MC2 = -1*On_off[t+1](name) -Shut_down[t](name) + On_off_sol_[t-1](name1) ;
            }
            else {
                string name = g->_name +",0" ;
                MC1 = On_off[t+1](name) -  On_off_initial(g->_name) -Start_up[t](name);
                MC2 = -1*On_off[t+1](name)+ On_off_initial(g->_name)  -Shut_down[t](name);

            }
            ACUC.add_constraint(MC1 <= 0);
            ACUC.add_constraint(MC2 <= 0);
        }
    }
}
void Global::add_OnOff_Sub_upper_time(Model& ACUC, int t) {
    for (auto& g: grid->gens) {
        if (g->_active) {
            Constraint OnOffStartupShutdown("OnOffStartupShutdown_"+ to_string(t) + ","+ g->_name);
            string name = g->_name +"," + to_string(t);
            if (t >0) {
                string name1 = g->_name +"," + to_string(t-1);
                OnOffStartupShutdown = On_off[t+1](name) - On_off_sol_[t-1](name1)
                                       - Start_up[t](name) + Shut_down[t](name);
            }
            else {
                OnOffStartupShutdown = On_off[t+1](name) - On_off_initial(g->_name)
                                       - Start_up[t](name) + Shut_down[t](name);
            }
            ACUC.add_constraint(OnOffStartupShutdown == 0);
        }
    }
}

void Global::add_Ramp_Sub_upper_time(Model& ACUC, int t) {
    if (t > 0) {
        for (auto& g: grid->gens) {
            if (g->_active) {
                Constraint Ramp_up("Ramp_up_constraint_"  + to_string(t) + "," + g->_name);
                Constraint Ramp_down("Ramp_down_constraint_"+ to_string(t) + "," + g->_name);
                string name = g->_name +"," + to_string(t);
                string name1 = g->_name +"," + to_string(t-1);
                Ramp_up =  Pg[t](name) - Pg_sol_[t-1](name1) -  rate_ramp.getvalue()*On_off_sol_[t-1](name1)
                           - rate_switch.getvalue()*(1 - On_off_sol_[t](name1));
                Ramp_down =   -1* Pg[t](name) + Pg_sol_[t-1](name1) - rate_ramp.getvalue()*On_off[t+1](name)- rate_switch.getvalue()*(1 - On_off[t+1](name));
                ACUC.add_constraint(Ramp_up <= 0);
                ACUC.add_constraint(Ramp_down <= 0);
            }
        }
    }
    else {
        for (auto& g: grid->gens) {
            if (g->_active) {
                Constraint Ramp_up("Ramp_up_constraint_0," + g->_name);
                string name = g->_name +",0";
                Ramp_up =  Pg[0](name) - Pg_initial(g->_name)
                           - rate_ramp(name)*On_off_initial(g->_name);
                Ramp_up -= rate_switch(name)*(1 - On_off_initial(g->_name));
                ACUC.add_constraint(Ramp_up <= 0);

                Constraint Ramp_down("Ramp_down_constraint_0," + g->_name);
                Ramp_down =   -1*Pg[0](name) + Pg_initial.in(g->_name);
                Ramp_down -= rate_ramp(name)*On_off[1](name);
                Ramp_down -= rate_switch(name)*(1 - On_off[1](name));
                ACUC.add_constraint(Ramp_down <= 0);
            }
        }
    }
}

void Global::add_minupdown_Sub_upper_time(Model& ACUC, int t) {
    if ( t >= min_up.getvalue()-1 && t < Num_time) {
        Constraint Min_Up("Min_Up_constraint_"+ to_string(t));
        Min_Up -= On_off[t+1].in_at(grid->gens, t);
        for (int l = t-min_up.getvalue()+1; l < t; l++) {
            Min_Up   += Start_up_sol_[l].in_at(grid->gens, l);
        }
        Min_Up   += Start_up[t].in_at(grid->gens, t);
        ACUC.add_constraint(Min_Up <= 0);
    }
    if ( t >= min_up.getvalue()-1 && t < Num_time) {
        Constraint Min_Down("Min_Down_constraint_" + to_string(t));
        Min_Down -= 1 - On_off[t+1].in_at(grid->gens, t);
        for (int l = t-min_down.getvalue()+1; l < t; l++) {
            Min_Down   += Shut_down_sol_[l].in_at(grid->gens, l);
        }
        Min_Down   += Shut_down[t].in_at(grid->gens, t);
        ACUC.add_constraint(Min_Down <= 0);
    }
}


double Global::Subproblem_time_(int t) {
    //Grid Parameters
    const auto bus_pairs = grid->get_bus_pairs();
    Model Sub("Sub_" + to_string(t));
    add_var_Sub_time(Sub, t);
    add_obj_Sub_time(Sub, t);
    add_perspective_OnOff_Sub_time(Sub, t);
    add_SOCP_Sub_time(Sub, t);
    add_KCL_Sub_time(Sub, t);
    add_thermal_Sub_time(Sub, t);
    add_MC_upper_Sub_time(Sub, t);
    if (t == 0) {
        for (auto& g: grid->gens) {
            if (g->_active) {
                Constraint MC1("Inter_temporal_MC1_0,"+ g->_name);
                Constraint MC2("Inter_temporal_MC2_0,"+ g->_name);
                string name = g->_name +",0" ;
                MC1 = On_off[t+1](name) -  On_off_initial(g->_name) -Start_up[t](name);
                MC2 = -1*On_off[t+1](name)+ On_off_initial(g->_name)  -Shut_down[t](name);
                Sub.add_constraint(MC1 <= 0);
                Sub.add_constraint(MC2 <= 0);
            }
        }

        for (auto& g: grid->gens) {
            Constraint Ramp_up("Ramp_up_constraint_0," + g->_name);
            string name = g->_name +",0";
            Ramp_up =  Pg[t](name) - Pg_initial(g->_name)
                       - rate_ramp(name)*On_off_initial(g->_name);
            Ramp_up -= rate_switch(name)*(1 - On_off_initial(g->_name));
            Sub.add_constraint(Ramp_up <= 0);

            Constraint Ramp_down("Ramp_down_constraint_0," + g->_name);
            Ramp_down =   -1*Pg[t](name) + Pg_initial.in(g->_name);
            Ramp_down -= rate_ramp(name)*On_off[t+1](name);
            Ramp_down -= rate_switch(name)*(1 - On_off[t+1](name));
            Sub.add_constraint(Ramp_down <= 0);
        }
        for (auto& g: grid->gens) {
            Constraint OnOffStartupShutdown("OnOffStartupShutdown_"+ to_string(t) + ","+ g->_name);
            string name = g->_name +"," + to_string(t);
            OnOffStartupShutdown = On_off[t+1](name) - On_off_initial(g->_name)
                                   - Start_up[t](name) + Shut_down[t](name);
            Sub.add_constraint(OnOffStartupShutdown == 0);
        }
    }

    /* Solver selection */
    solver cpx_acuc(Sub, cplex);
    bool relax = false;
    int output = 1;
    double tol = 1e-6;
    cpx_acuc.run(output, relax, tol);
    if (t==0) {
        // fill solution
        On_off_sol_[t].set_name(" On_off_sol_"+to_string(t));
        On_off_sol_[t] = *(param<bool>*)(Sub.get_var("On_off_"+to_string(t)));
        //On_off_sol_[t].print(true);
        Start_up_sol_[t] = *(param<bool>*)(Sub.get_var("Start_up_"+to_string(t)));
        //Start_up_sol_[t].print(true);
        Shut_down_sol_[t] = *(param<bool>*)(Sub.get_var("Shut_down_"+to_string(t)));
        //Shut_down_sol_[t].print(true);
        Pg_sol_[t] = *(param<Real>*)(Sub.get_var("Pg_"+to_string(t)));
    }

    return Sub._obj_val;
}

double Global::Subproblem_upper_time_(int t) {
    const auto bus_pairs = grid->get_bus_pairs();
    Model Sub("Sub_" + to_string(t));
    add_var_Sub_time(Sub, t);
    add_obj_Sub_time(Sub, t);
    add_obj_Sub_upper_time(Sub, t);
    //add_perspective_OnOff_Sub_time(Sub, t);
    add_SOCP_Sub_time(Sub, t);
    add_KCL_Sub_time(Sub, t);
    add_thermal_Sub_time(Sub, t);
    add_MC_upper_Sub_time(Sub, t);
    // add interesting constraints
    add_MC_intertemporal_Sub_upper_time(Sub, t);
    add_OnOff_Sub_upper_time(Sub, t);
    add_Ramp_Sub_upper_time(Sub, t);
    add_minupdown_Sub_upper_time(Sub, t);

    double ub = 0;
    /* Solver selection */
    solver cpx_acuc(Sub, cplex);
    bool relax = false;
    int output = 1;
    double tol = 1e-6;
    cpx_acuc.run(output, relax, tol);
    // fill solution
    On_off_sol_[t].set_name(" On_off_sol_"+to_string(t));
    On_off_sol_[t] = *(param<bool>*)(Sub.get_var("On_off_"+to_string(t)));
    //On_off_sol_[t].print(true);
    Start_up_sol_[t] = *(param<bool>*)(Sub.get_var("Start_up_"+to_string(t)));
    //Start_up_sol_[t].print(true);
    Shut_down_sol_[t] = *(param<bool>*)(Sub.get_var("Shut_down_"+to_string(t)));
    //Shut_down_sol_[t].print(true);
    Pg_sol_[t] = *(param<Real>*)(Sub.get_var("Pg_"+to_string(t)));
    //Pg_sol_[t].print(true);
    return Sub._obj_val;
}


double Global::LR_bound_time_(bool included_min_up_down) {
    include_min_updown_ = included_min_up_down;
    double LB = 0;
    for(int t = 0; t < Num_time; t++) {
        Sub_[t]= Subproblem_time_(t);
        LB += Sub_[t];
    }
    for (int t = 1; t < Num_time; t++) {
        for (auto& g: grid->gens) {
            string name = g->_name + "," + to_string(t);
            LB -= zeta_down(name).getvalue()*rate_switch(name).getvalue();
            LB -= zeta_up(name).getvalue()*rate_switch(name).getvalue();
        }
    }
    if (include_min_updown_) {
        if (min_down.getvalue() - 1.0 > 0) {
            for (int t = min_down.getvalue()-1;  t < Num_time; t++ ) {
                for (auto& g: grid->gens) {
                    string name = g->_name + ","+to_string(t);
                    LB -= mu_down(name).getvalue();
                }
            }
        }
    }
    return LB;
}

double Global::Upper_bound_sequence_(bool included_min_up_down) {
    include_min_updown_ = included_min_up_down;
    On_off_sol_.resize(Num_time);
    Start_up_sol_.resize(Num_time);
    Shut_down_sol_.resize(Num_time);
    double UB = 0;
    for(int t = 0; t < Num_time; t++) {
        if (t >=1) {
            Sub_[t]= Subproblem_upper_time_(t);
        }
        else {
            func_ obj;
            for (auto g:grid->gens) {
                if (g->_active) {
                    string name = g->_name + ","+ to_string(t);
                    obj += (grid->c1(name)*Pg_sol_[t](name)+ grid->c2(name)*Pg_sol_[t](name)*Pg_sol_[t](name) + grid->c0(name)*On_off_sol_[t](name));
                    obj += cost_up.getvalue()*Start_up_sol_[t](name)+ cost_down.getvalue()*Shut_down_sol_[t](name);
                }
            }
            Sub_[t] = poly_eval(&obj);
        }
        UB += Sub_[t];
    }
    return UB;
}

double Global::Subproblem_spatial_(int l) {
    assert(Pg.size() == Num_parts);
    Model Subr("Subr");
    Subr.add_var(Pg[l].in(P_->bag_gens[l], Num_time));
    Subr.add_var(Pg2[l].in(P_->bag_gens[l], Num_time));
    Subr.add_var(Qg[l].in(P_->bag_gens[l], Num_time));
    Subr.add_var(On_off[l].in(P_->bag_gens[l], -1, Num_time));
    Subr.add_var(Start_up[l].in(P_->bag_gens[l], Num_time));
    Subr.add_var(Shut_down[l].in(P_->bag_gens[l], Num_time));
    Subr.add_var(Xii[l].in(P_->bag_bus_union_out[l], Num_time));
    Xii[l].initialize_all(1.001);
    Subr.add_var(R_Xij[l].in(P_->bag_bus_pairs_union[l], Num_time));
    R_Xij[l].initialize_all(1.0);
    Subr.add_var(Im_Xij[l].in(P_->bag_bus_pairs_union[l], Num_time));
    //power flow
    var<Real> Pf_from("Pf_from", grid->S_max.in(P_->bag_arcs_union_out[l], Num_time));
    var<Real> Qf_from("Qf_from", grid->S_max.in(P_->bag_arcs_union_out[l], Num_time));
    var<Real> Pf_to("Pf_to", grid->S_max.in(P_->bag_arcs_union_in[l], Num_time));
    var<Real> Qf_to("Qf_to", grid->S_max.in(P_->bag_arcs_union_in[l], Num_time));
    if (P_->bag_arcs_union[l].size() > 0) {
        Subr.add_var(Pf_from.in(P_->bag_arcs_union_out[l], Num_time));
        Subr.add_var(Qf_from.in(P_->bag_arcs_union_out[l], Num_time));
        Subr.add_var(Pf_to.in(P_->bag_arcs_union_in[l], Num_time));
        Subr.add_var(Qf_to.in(P_->bag_arcs_union_in[l], Num_time));
    }
    /* Construct the objective function*/
    func_ obj;
    for (auto g:P_->bag_gens[l]) {
        if (g->_active) {
            for (int t = 0; t < Num_time; t++) {
                string name = g->_name + ","+ to_string(t);
                obj += grid->c1(name)*Pg[l](name)+ grid->c2(name)*Pg2[l](name) +grid->c0(name)*On_off[l](name);
                obj += cost_up*Start_up[l](name)+ cost_down*Shut_down[l](name);
            }
        }
    }
    const auto& bag = P_->G_part.get_node(to_string(l));
    for (const auto& a: bag->get_out()) {
        for (int t = 0; t < Num_time; t++) {
            for (auto& pair: a->_intersection_clique) {
                string name = pair->_name + ","+ to_string(t);
                string dest = pair->_dest->_name + ","+ to_string(t);
                obj += R_lambda_(name)*R_Xij[l](name) +  Im_lambda_(name)*Im_Xij[l](name) + lambda_(name)*Xii[l](dest);
            }
        }
    }

    for (const auto& a: bag->get_in()) {
        for (unsigned t = 0; t < Num_time; t++) {
            for (auto pair: a->_intersection_clique) {
                string name = pair->_name + ","+ to_string(t);
                string dest = pair->_dest->_name + ","+ to_string(t);
                obj -= R_lambda_(name)*R_Xij[l](name) +  Im_lambda_(name)*Im_Xij[l](name) + lambda_(name)*Xii[l](dest);
            }
        }
    }
    Subr.min(obj);
    /* perspective formulation of Pg[l]^2 */
    if (P_->bag_gens[l].size() >0) {
        Constraint Perspective_OnOff_Pg2("Perspective_OnOff_Pg2_");
        Perspective_OnOff_Pg2 = power(Pg[l], 2) - Pg2[l]*On_off[l];
        Subr.add(Perspective_OnOff_Pg2.in(P_->bag_gens[l], Num_time) <= 0);
    }

    if (P_->bag_bus_pairs_union_directed[l].size() > 0) {
        Constraint SOC("SOC_" + to_string(l));
        SOC =  power(R_Xij[l], 2)+ power(Im_Xij[l], 2)-Xii[l].from()*Xii[l].to() ;
        Subr.add_constraint(SOC.in(P_->bag_bus_pairs_union_directed[l], Num_time) <= 0);
    }
    Constraint KCL_P("KCL_P");
    Constraint KCL_Q("KCL_Q");
    KCL_P = sum(Pf_from.out_arcs()) + sum(Pf_to.in_arcs()) + grid->pl-sum(Pg[l].in_gens()) + grid->gs*Xii[l];
    Subr.add_constraint(KCL_P.in(P_->bag_bus[l], Num_time) == 0);
    KCL_Q = sum(Qf_from.out_arcs()) + sum(Qf_to.in_arcs()) + grid->ql - sum(Qg[l].in_gens()) - grid->bs*Xii[l];
    Subr.add_constraint(KCL_Q.in(P_->bag_bus[l], Num_time) == 0);

    Constraint Flow_P_From("Flow_P_From");
    Flow_P_From = Pf_from -(grid->g_ff*Xii[l].from() + grid->g_ft*R_Xij[l].in_pairs() + grid->b_ft*Im_Xij[l].in_pairs());
    Subr.add_constraint(Flow_P_From.in(P_->bag_arcs_union_out[l], Num_time) == 0);

    Constraint Flow_P_Num_timeo("Flow_P_Num_time");
    Flow_P_Num_timeo = Pf_to - (grid->g_tt*Xii[l].to() + grid->g_tf*R_Xij[l].in_pairs() - grid->b_tf*Im_Xij[l].in_pairs());
    Subr.add_constraint(Flow_P_Num_timeo.in(P_->bag_arcs_union_in[l], Num_time) == 0);

    Constraint Flow_Q_From("Flow_Q_From");
    Flow_Q_From = Qf_from - (grid->g_ft*Im_Xij[l].in_pairs() - grid->b_ff*Xii[l].from() - grid->b_ft*R_Xij[l].in_pairs());
    Subr.add_constraint(Flow_Q_From.in(P_->bag_arcs_union_out[l], Num_time) == 0);

    Constraint Flow_Q_Num_timeo("Flow_Q_Num_time");
    Flow_Q_Num_timeo = Qf_to + (grid->b_tt*Xii[l].to() + grid->b_tf*R_Xij[l].in_pairs() + grid->g_tf*Im_Xij[l].in_pairs());
    Subr.add_constraint(Flow_Q_Num_timeo.in(P_->bag_arcs_union_in[l], Num_time) == 0);

    Constraint Num_timehermal_Limit_from("Num_timehermal_Limit_from");
    Num_timehermal_Limit_from = power(Pf_from, 2) + power(Qf_from, 2);
    Num_timehermal_Limit_from <= power(grid->S_max, 2);
    Subr.add_constraint(Num_timehermal_Limit_from.in(P_->bag_arcs_union_out[l], Num_time));

    Constraint Num_timehermal_Limit_to("Num_timehermal_Limit_to");
    Num_timehermal_Limit_to = power(Pf_to, 2) + power(Qf_to, 2);
    Num_timehermal_Limit_to <= power(grid->S_max, 2);
    Subr.add_constraint(Num_timehermal_Limit_to.in(P_->bag_arcs_union_in[l], Num_time));
    /* Phase Angle Bounds constraints */
    Constraint PAD_UB("PAD_UB_"+to_string(l));
    PAD_UB = Im_Xij[l]- grid->tan_th_max*R_Xij[l];
    Subr.add_constraint(PAD_UB.in(P_->bag_bus_pairs_union_directed[l], Num_time) <= 0);

    Constraint PAD_LB("PAD_LB_"+to_string(l));
    PAD_LB = Im_Xij[l] - grid->tan_th_min*R_Xij[l];
    Subr.add_constraint(PAD_LB.in(P_->bag_bus_pairs_union_directed[l], Num_time) >= 0);
    //COMMINum_timeMENNum_time CONSNum_timeRAINNum_timeS
    // Inter-temporal constraints 3a, 3d
    if (P_->bag_gens[l].size() > 0) {
        for (int t = 0; t < Num_time; t++) {
            Constraint MC1("Inter_temporal_MC1_"+to_string(t));
            Constraint MC2("Inter_temporal_MC2_"+to_string(t));
            MC1 = On_off[l].in_at(P_->bag_gens[l], t)- On_off[l].in_at(P_->bag_gens[l], t-1)-  Start_up[l].in_at(P_->bag_gens[l], t);
            MC2 = On_off[l].in_at(P_->bag_gens[l], t-1) - On_off[l].in_at(P_->bag_gens[l], t) - Shut_down[l].in_at(P_->bag_gens[l], t);
            Subr.add_constraint(MC1 <= 0);
            Subr.add_constraint(MC2 <= 0);
        }
    }
    if (P_->bag_gens[l].size() > 0) {
        for (int t = 0; t < Num_time; t++) {
            Constraint OnOffStartupShutdown("OnOffStartupShutdown_"+ to_string(t));
            OnOffStartupShutdown = On_off[l].in_at(P_->bag_gens[l], t) - On_off[l].in_at(P_->bag_gens[l], t-1)
                                   - Start_up[l].in_at(P_->bag_gens[l], t) + Shut_down[l].in_at(P_->bag_gens[l], t);
            Subr.add_constraint(OnOffStartupShutdown == 0);
        }
    }

    if (P_->bag_gens[l].size() > 0) {
        for (int t = min_up.getvalue()-1; t < Num_time; t++) {
            Constraint Min_Up("Min_Up_constraint" + to_string(l) + "_"+ to_string(t));
            for (int ll = t-min_up.getvalue()+1; ll < t+1; ll++) {
                Min_Up   += Start_up[l].in_at(P_->bag_gens[l], ll);
            }
            Min_Up -= On_off[l].in_at(P_->bag_gens[l], t);
            Subr.add_constraint(Min_Up <= 0);
        }
    }
    // 4c
    if (P_->bag_gens[l].size() > 0) {
        for (int t = min_down.getvalue()-1; t < Num_time; t++) {
            Constraint Min_Down("Min_Down_constraint_"+ to_string(t));
            for (int ll = t-min_down.getvalue()+1; ll < t +1; ll++) {
                Min_Down   += Shut_down[l].in_at(P_->bag_gens[l], ll);
            }
            Min_Down -= 1 - On_off[l].in_at(P_->bag_gens[l], t);
            Subr.add_constraint(Min_Down <= 0);
        }
    }

    ////Ramp Rate
    if (P_->bag_gens[l].size() > 0) {
        Constraint Production_P_LB("Production_P_LB");
        Constraint Production_P_UB("Production_P_UB");
        Constraint Production_Q_LB("Production_Q_LB");
        Constraint Production_Q_UB("Production_Q_UB");
        // 5A
        Production_P_UB = Pg[l] - grid->pg_max*On_off[l];
        Production_P_LB = Pg[l] - grid->pg_min*On_off[l];
        Subr.add_constraint(Production_P_UB.in(P_->bag_gens[l], Num_time) <=0);
        Subr.add_constraint(Production_P_LB.in(P_->bag_gens[l], Num_time) >= 0);

        Production_Q_UB = Qg[l] - grid->qg_max*On_off[l];
        Production_Q_LB = Qg[l] - grid->qg_min*On_off[l];
        Subr.add_constraint(Production_Q_UB.in(P_->bag_gens[l], Num_time) <= 0);
        Subr.add_constraint(Production_Q_LB.in(P_->bag_gens[l], Num_time) >= 0);
    }
    // 5C
    if (P_->bag_gens[l].size() > 0) {
        for (int t = 1; t < Num_time; t++) {
            Constraint Ramp_up("Ramp_up_constraint_" + to_string(t));
            Constraint Ramp_down("Ramp_down_constraint_" +to_string(t));
            Ramp_up =  Pg[l].in_at(P_->bag_gens[l], t);
            Ramp_up -= Pg[l].in_at(P_->bag_gens[l], t-1);
            Ramp_up -= rate_ramp*On_off[l].in_at(P_->bag_gens[l], t-1);
            Ramp_up -= rate_switch*(1 - On_off[l].in_at(P_->bag_gens[l], t));

            Ramp_down =  Pg[l].in_at(P_->bag_gens[l], t-1);
            Ramp_down -= Pg[l].in_at(P_->bag_gens[l], t);
            Ramp_down -= rate_ramp*On_off[l].in_at(P_->bag_gens[l], t);
            Ramp_down -= rate_switch*(1 - On_off[l].in_at(P_->bag_gens[l], t-1));

            Subr.add_constraint(Ramp_up <= 0);
            Subr.add_constraint(Ramp_down <= 0);
        }

        Constraint Ramp_up("Ramp_up_constraint0" + to_string(l));
        Ramp_up =  Pg[l].in_at(P_->bag_gens[l], 0) - Pg_initial.in(P_->bag_gens[l]) - rate_ramp.in_at(P_->bag_gens[l], 0)*On_off_initial.in(P_->bag_gens[l]);
        Ramp_up -= rate_switch.in_at(P_->bag_gens[l], 0)*(1 - On_off_initial.in(P_->bag_gens[l]));
        Subr.add_constraint(Ramp_up <= 0);

        Constraint Ramp_down("Ramp_down_constraint0");
        Ramp_down =   -1*Pg[l].in_at(P_->bag_gens[l],0) + Pg_initial.in(P_->bag_gens[l]);
        Ramp_down -= rate_ramp.in_at(P_->bag_gens[l], 0)*On_off[l].in_at(P_->bag_gens[l], 0);
        Ramp_down -= rate_switch.in_at(P_->bag_gens[l], 0)*(1 - On_off[l].in_at(P_->bag_gens[l], 0));
        Subr.add_constraint(Ramp_down <= 0);
    }
    //// set the initial state of generators.
    if (P_->bag_gens[l].size() >0) {
        Constraint gen_initial("gen_initial_"+to_string(l));
        gen_initial +=  On_off[l].in_at(P_->bag_gens[l], -1) - On_off_initial.in(P_->bag_gens[l]);
        Subr.add_constraint(gen_initial == 0);
    }
    /* solve it! */
    solver solve_Subr(Subr, cplex);
    bool relax = true;
    int output = 1;
    double tol = 1e-6;
    solve_Subr.run(output, relax, tol);
    // LINKED VARIABLES
    // std::string name = Xii[l].in(P_->bag_bus[l],Num_time).get_name();
//    Xii_log= (*(var<Real>*) Subr.get_var(name));
//    name = R_Xij_.in(P_->bag_bus_pairs_union[l], Num_time).get_name();
//    R_Xij_log = (*(var<Real>*) Subr.get_var(name));
//    name = Im_Xij.in(P_->bag_bus_pairs_union[l], Num_time).get_name();
//    Im_Xij_log = (*(var<Real>*) Subr.get_var(name));
    return Subr._obj_val;
}

double Global::LR_bound_spatial_() {
    double LB = 0;
    Sub_.clear();
    Sub_.resize(Num_parts);
    for(int l = 0; l < Num_parts; l++) {
        Sub_[l]= Subproblem_spatial_(l);
        LB += Sub_[l];
    }
    return LB;
}
