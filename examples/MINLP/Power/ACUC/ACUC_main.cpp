//
//  ACUC.cpp
//  Gravity
//
//  Created by Guanglei Wang on 6/9/17.
//
//
#include <stdio.h>
#include <iostream>
#include <string>
#include <stdio.h>
#include <cstring>
#include <fstream>
#include "../PowerNet.h"
#include <gravity/solver.h>
#include <stdio.h>
#include <stdlib.h>
using namespace std;
using namespace gravity;

int main (int argc, const char * argv[])
{
    const char* fname;
    if (argc >= 2) {
        fname = argv[1];
    }
    else {
        // fname = "../../data_sets/Power/nesta_case3_lmbd.m";
        // fname = "/Users/hh/Dropbox/Work/Dev/nesta-0.7.0/opf/nesta_case2383wp_mp.m";
        // fname = "../../data_sets/Power/nesta_case3_lmbd.m";
        // fname = "../../data_sets/Power/nesta_case2383wp_mp.m";
         fname = "../../data_sets/Power/nesta_case5_pjm.m";
        //fname = "../../data_sets/Power/nesta_case300_ieee.m";
    }
    // ACUC
    PowerNet* grid = new PowerNet();
    grid->readgrid(fname);

    // Grid Parameters
    auto bus_pairs = grid->get_bus_pairs();
    auto nb_bus_pairs = bus_pairs.size();
    auto nb_gen = grid->get_nb_active_gens();
    auto nb_lines = grid->get_nb_active_arcs();
    auto nb_buses = grid->get_nb_active_nodes();

    // Schedule
    unsigned T = 1;
    param<Real> rate_ramp("rate_ramp");
    param<Real> rate_switch("rate_switch");
    param<Real> min_up("min_up");
    param<Real> min_down("min_down");
    param<Real> cost_up("cost_up");
    param<Real> cost_down("cost_down");
    for (auto g: grid->gens) {
        rate_ramp(g->_name) = max(grid->pg_min(g->_name).getvalue(), 0.25*grid->pg_max(g->_name).getvalue());
        rate_switch(g->_name) = max(grid->pg_min(g->_name).getvalue(), 0.25*grid->pg_max(g->_name).getvalue());
    }
    min_up = 1;
    min_down = 1;
    cost_up = 50;
    cost_down = 30;
    
    grid->time_expand(T);
    rate_ramp.time_expand(T);
    rate_switch.time_expand(T);

    /** build model */
    Model ACUC("ACUC Model");

    /** Variables */
    // POWER GENERATION
    var<Real> Pg("Pg", grid->pg_min.in(grid->gens, T), grid->pg_max.in(grid->gens, T));
    var<Real> Qg ("Qg", grid->qg_min.in(grid->gens, T), grid->qg_max.in(grid->gens, T));
    ACUC.add_var(Pg^(T*nb_gen));
    ACUC.add_var(Qg^(T*nb_gen));

     //power flow
    var<Real> Pf_from("Pf_from", grid->S_max.in(grid->arcs, T));
    var<Real> Qf_from("Qf_from", grid->S_max.in(grid->arcs, T));
    var<Real> Pf_to("Pf_to", grid->S_max.in(grid->arcs, T));
    var<Real> Qf_to("Qf_to", grid->S_max.in(grid->arcs, T));
    ACUC.add_var(Pf_from^(T*nb_lines));
    ACUC.add_var(Qf_from^(T*nb_lines));
    ACUC.add_var(Pf_to^(T*nb_lines));
    ACUC.add_var(Qf_to^(T*nb_lines));

     //Lifted variables.
    var<Real>  R_Wij("R_Wij", grid->wr_min.in(bus_pairs, T), grid->wr_max.in(bus_pairs, T)); // real part of Wij
    var<Real>  Im_Wij("Im_Wij", grid->wi_min.in(bus_pairs, T), grid->wi_max.in(bus_pairs, T)); // imaginary part of Wij.
    var<Real>  Wii("Wii", grid->w_min.in(grid->nodes, T), grid->w_max.in(grid->nodes, T));
    R_Wij.print(true);
    Im_Wij.print(true);
    
    ACUC.add_var(Wii^(T*nb_buses));
    ACUC.add_var(R_Wij^(T*nb_bus_pairs));
    ACUC.add_var(Im_Wij^(T*nb_bus_pairs));
    R_Wij.initialize_all(1.0);
    Wii.initialize_all(1.001);

    // Commitment variables
    var<Real>  On_off("On_off");
    var<Real>  Start_up("Start_up", 0, 1);
    var<Real>  Shut_down("Shut_down", 0, 1);
    ACUC.add_var(On_off^(T*nb_gen));
    ACUC.add_var(Start_up^(T*nb_gen));
    ACUC.add_var(Shut_down^(T*nb_gen));

    /* Construct the objective function*/
    func_ obj;
    for (auto g:grid->gens) {
        if (g->_active) {
            obj += grid->c1.in(g,T)*Pg.in(g,T)+ grid->c2.in(g,T)*Pg.in(g,T)*Pg.in(g,T) + grid->c0.in(g,T);
            obj += cost_up.getvalue()*Start_up.in(g, T)+ cost_down.getvalue()*Shut_down.in(g, T);
        }
    }
    ACUC.set_objective(min(obj));
    
    /** Define constraints */
    /* SOCP constraints */
    Constraint SOC("SOC");
    SOC =  power(R_Wij.in(bus_pairs, T), 2) + power(Im_Wij.in(bus_pairs, T), 2) - Wii.from(bus_pairs, T)*Wii.to(bus_pairs, T) ;
    ACUC.add_constraint(SOC <= 0);
    
    //KCL
    for (int t = 0; t < T; t++)
        for (auto b: grid->nodes) {
            if (!b->_active) {
                continue;
            }
            Bus* bus = (Bus*) b;
            Constraint KCL_P("KCL_P"+bus->_name+ "time_" + to_string(t));
            Constraint KCL_Q("KCL_Q"+bus->_name+ "time_" + to_string(t));

            /* Power Conservation */
            KCL_P  = sum(Pf_from.in_at(b->get_out(), t)) + sum(Pf_to.in_at(b->get_in(), t)) + grid->pl(bus->_name, to_string(t))- sum(Pg.in_at(bus->_gen, t));
            KCL_Q  = sum(Qf_from.in_at(b->get_out(), t)) + sum(Qf_to.in_at(b->get_in(), t)) + grid->ql(bus->_name, to_string(t))- sum(Qg.in_at(bus->_gen, t));

            /* Shunts */
            KCL_P +=  grid->gs(bus->_name, to_string(t))*Wii(bus->_name, to_string(t));
            KCL_Q -=  grid->bs(bus->_name, to_string(t))*Wii(bus->_name, to_string(t));

            ACUC.add_constraint(KCL_P == 0);
            ACUC.add_constraint(KCL_Q == 0);
        }

    ////AC Power Flow.
    Constraint Flow_P_From("Flow_P_From");
    Flow_P_From += Pf_from.in(grid->arcs, T);
    Flow_P_From -= grid->g_ff.in(grid->arcs, T)*Wii.from(grid->arcs, T);
    Flow_P_From -= grid->g_ft.in(grid->arcs, T)*R_Wij.in_pairs(grid->arcs, T);
    Flow_P_From -= grid->b_ft.in(grid->arcs, T)*Im_Wij.in_pairs(grid->arcs, T);
    ACUC.add_constraint(Flow_P_From == 0);
    

    Constraint Flow_P_To("Flow_P_To");
    Flow_P_To += Pf_to.in(grid->arcs, T);
    Flow_P_To -= grid->g_tt.in(grid->arcs, T)*Wii.to(grid->arcs, T);
    Flow_P_To -= grid->g_tf.in(grid->arcs, T)*R_Wij.in_pairs(grid->arcs, T);
    Flow_P_To += grid->b_tf.in(grid->arcs, T)*Im_Wij.in_pairs(grid->arcs, T);
    ACUC.add_constraint(Flow_P_To == 0);

    Constraint Flow_Q_From("Flow_Q_From");
    Flow_Q_From += Qf_from.in(grid->arcs, T);
    Flow_Q_From += grid->b_ff.in(grid->arcs, T)*Wii.from(grid->arcs, T);
    Flow_Q_From += grid->b_ft.in(grid->arcs, T)*R_Wij.in_pairs(grid->arcs, T);
    Flow_Q_From += grid->g_ft.in(grid->arcs, T)*Im_Wij.in_pairs(grid->arcs, T);
    ACUC.add_constraint(Flow_Q_From == 0);

    Constraint Flow_Q_To("Flow_Q_To");
    Flow_Q_To += Qf_to.in(grid->arcs, T);
    Flow_Q_To += grid->b_tt.in(grid->arcs, T)*Wii.to(grid->arcs, T);
    Flow_Q_To += grid->b_tf.in(grid->arcs, T)*R_Wij.in_pairs(grid->arcs, T);
    Flow_Q_To -= grid->g_tf.in(grid->arcs, T)*Im_Wij.in_pairs(grid->arcs, T);
    ACUC.add_constraint(Flow_Q_To == 0);
//
    /* Phase Angle Bounds constraints */
    Constraint PAD_UB("PAD_UB");
    PAD_UB = Im_Wij.in(bus_pairs, T);
    PAD_UB -= (grid->tan_th_max).in(bus_pairs, T)*R_Wij.in(bus_pairs, T);
    ACUC.add_constraint(PAD_UB <= 0);

    Constraint PAD_LB("PAD_LB");
    PAD_LB =  Im_Wij.in(bus_pairs, T);
    PAD_LB -= grid->tan_th_min.in(bus_pairs, T)*R_Wij.in(bus_pairs, T);
    ACUC.add_constraint(PAD_LB >= 0);

    ///* Thermal Limit Constraints */
    Constraint Thermal_Limit_from("Thermal_Limit_from");
    Thermal_Limit_from += power(Pf_from.in(grid->arcs, T),  2) + power(Qf_from.in(grid->arcs, T), 2);
    Thermal_Limit_from -= power(grid->S_max.in(grid->arcs, T), 2);
    ACUC.add_constraint(Thermal_Limit_from <= 0);

    Constraint Thermal_Limit_to("Thermal_Limit_to");
    Thermal_Limit_to += power(Pf_to.in(grid->arcs, T), 2) + power(Qf_to.in(grid->arcs, T), 2);
    Thermal_Limit_to -= power(grid->S_max.in(grid->arcs, T),2);
    ACUC.add_constraint(Thermal_Limit_to <= 0);

    // COMMITMENT CONSTRAINTS
    // Inter-temporal constraints
    for (int t = 1; t < T; t++) {
        Constraint MC1("MC1_"+ to_string(t));
        Constraint MC2("MC2_"+ to_string(t));
        MC1 = On_off.in_at(grid->gens, t)- On_off.in_at(grid->gens, t-1)-  Start_up.in_at(grid->gens, t);
        MC2 = On_off.in_at(grid->gens, t-1) - On_off.in_at(grid->gens, t) - Shut_down.in_at(grid->gens, t);
        ACUC.add_constraint(MC1 <= 0);
        ACUC.add_constraint(MC2 <= 0);
    }

    //// Min-up constraints
    for (int t = 1; t < T; t++) {
        Constraint Min_up1("Min_up1_"+ to_string(t));
        Min_up1 = On_off.in_at(grid->gens, t) - On_off.in_at(grid->gens, t-1) - Start_up.in_at(grid->gens, t) + Shut_down.in_at(grid->gens, t);
        ACUC.add_constraint(Min_up1 == 0);
    }

    for (int t = min_up.getvalue(); t < T; t++) {
        Constraint Min_Up("Min_Up_constraint" + to_string(t));
        for (int l = t-min_up.getvalue()+1; l < t+1; l++) {
            Min_Up   += Start_up.in_at(grid->gens, l);
        }
        Min_Up -= On_off.in_at(grid->gens, t);
        ACUC.add_constraint(Min_Up <= 0);
    }

    for (int t = min_down.getvalue(); t < T; t++) {
        Constraint Min_Down("Min_Down_constraint" + to_string(t));
        for (int l = t-min_down.getvalue()+1; l < t +1; l++) {
            Min_Down   += Shut_down.in_at(grid->gens, l);
        }
        Min_Down -= 1 - On_off.in_at(grid->gens, t);
        ACUC.add_constraint(Min_Down <= 0);
    }

    ////Ramp rate
    Constraint Production_P_LB("Production_P_LB");
    Constraint Production_P_UB("Production_P_UB");
    Constraint Production_Q_LB("Production_Q_LB");
    Constraint Production_Q_UB("Production_Q_UB");

    Production_P_UB = Pg.in(grid->gens, T) - grid->pg_max.in(grid->gens, T)*On_off.in(grid->gens,T);
    Production_P_LB = Pg.in(grid->gens, T) - grid->pg_min.in(grid->gens, T)*On_off.in(grid->gens,T);
    ACUC.add_constraint(Production_P_UB <=0);
    ACUC.add_constraint(Production_P_LB >= 0);

    grid->qg_max.print(true);
    grid->qg_min.print(true);

    Production_Q_UB = Qg.in(grid->gens, T) - grid->qg_max.in(grid->gens, T)*On_off.in(grid->gens,T);
    Production_Q_LB = Qg.in(grid->gens, T) - grid->qg_min.in(grid->gens, T)*On_off.in(grid->gens,T);
    ACUC.add_constraint(Production_Q_UB <= 0);
    ACUC.add_constraint(Production_Q_LB >= 0);

    for (int t = 1; t < T; t++) {
        Constraint Ramp_up("Ramp_up_constraint" + to_string(t));
        Constraint Ramp_down("Ramp_down_constraint" + to_string(t));

        Ramp_up = Pg.in_at(grid->gens, t);
        Ramp_up -= Pg.in_at(grid->gens, t-1);
        Ramp_up -= rate_ramp*On_off.in_at(grid->gens, t-1);
        Ramp_up -= rate_switch*(1 - On_off.in_at(grid->gens, t));

        Ramp_down = Pg.in_at(grid->gens, t-1);
        Ramp_down -= Pg.in_at(grid->gens, t);
        Ramp_down -= rate_ramp*On_off.in_at(grid->gens, t);
        Ramp_down -= rate_switch*(1 - On_off.in_at(grid->gens, t-1));

        ACUC.add_constraint(Ramp_up <= 0);
        ACUC.add_constraint(Ramp_down <= 0);

    }

    /* Resolve it! */
    //solver OPF(ACUC,ipopt);
    solver OPF(ACUC, cplex);
    OPF.run();

    return 0;
}
