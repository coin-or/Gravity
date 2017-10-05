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
#include <gravity/model.h>
#include <gravity/solver.h>
#include <stdio.h>
#include <stdlib.h>
using namespace std;
using namespace gravity;

int main (int argc, const char * argv[])
{
    // ACUC
    PowerNet* grid = new PowerNet();
    const char* fname;
    fname = "../../data_sets/Power/nesta_case5_pjm.m";
    grid->readgrid(fname);
    
    // Grid Parameters
    int nb_gen = grid->gens.size();
    int nb_lines = grid->arcs.size();
    int nb_buses = grid->nodes.size();
    
    // Schedule periods
    unsigned T = 2;
    grid->c0.time_expand(T);
    grid->c1.time_expand(T);
    grid->c2.time_expand(T);
    grid->S_max.time_expand(T);
    grid->tan_th_min.time_expand(T);
    grid->tan_th_max.time_expand(T);
    grid->g_tt.time_expand(T);
    grid->g_ff.time_expand(T);
    grid->g_ft.time_expand(T);
    grid->g_tf.time_expand(T);
    grid->b_tt.time_expand(T);
    grid->b_ff.time_expand(T);
    grid->b_ft.time_expand(T);
    grid->b_tf.time_expand(T);
    grid->pg_min.time_expand(T);
    grid->pg_max.time_expand(T);
    grid->qg_min.time_expand(T);
    grid->qg_max.time_expand(T);
    grid->w_min.time_expand(T);
    grid->w_max.time_expand(T);

    /** build model */
    Model ACUC("ACUC Model");

    /** Variables */
    // power generation
    //var<Real> Pg("Pg", grid->pg_min, grid->pg_max);
    var<Real> Pg("Pg", grid->pg_min.in(grid->gens, T), grid->pg_max.in(grid->gens, T));
    var<Real> Qg ("Qg", grid->qg_min.in(grid->gens, T), grid->qg_max.in(grid->gens, T));
    ACUC.add_var(Pg^(T*nb_gen));
    ACUC.add_var(Qg^(T*nb_gen));

    // power flow
    var<double> Pf_from("Pf_from");
    var<double> Qf_from("Qf_from");
    var<double> Pf_to("Pf_to");
    var<double> Qf_to("Qf_to");
    ACUC.add_var(Pf_from^(T*nb_lines));
    ACUC.add_var(Qf_from^(T*nb_lines));
    ACUC.add_var(Pf_to^(T*nb_lines));
    ACUC.add_var(Qf_to^(T*nb_lines));

    // Lifted variables.
    var<Real>  R_Wij("R_Wij");   // real part of Wij
    var<Real>  Im_Wij("Im_Wij"); // imaginary part of Wij.
    var<Real>  Wii("Wii", grid->w_min, grid->w_max);
    ACUC.add_var(Wii^(T*nb_buses));
    ACUC.add_var(R_Wij^(T*nb_lines));
    ACUC.add_var(Im_Wij^(T*nb_lines));

    // Commitment variables
    var<Real>  On_off("On_off", 0, 1);
    var<Real>  Start_up("Start_up", 0, 1);
    var<Real>  Shut_down("Shut_down", 0, 1);
    ACUC.add_var(On_off^(T*nb_gen));
    ACUC.add_var(Start_up^(T*nb_gen));
    ACUC.add_var(Shut_down^(T*nb_gen));


    /* Construct the objective function*/
    func_ obj;
    obj  = sum(grid->c0.in(grid->gens, T));
    obj += sum(grid->c1.in(grid->gens, T), Pg.in(grid->gens, T));
    obj += sum(grid->c2.in(grid->gens, T), power(Pg.in(grid->gens, T), 2));
    ACUC.set_objective(min(obj));

    grid->c1.in(grid->gens, T).print(true);

    /** Define constraints */
    /* SOCP constraints */
    Constraint SOC("SOC");
    SOC =  power(R_Wij.in(grid->arcs, T), 2) + power(Im_Wij.in(grid->arcs, T), 2) - Wii.from(grid->arcs, T)*Wii.to(grid->arcs, T) ;
    ACUC.add_constraint(SOC <= 0);
    //KCL
    for (int t = 0; t < T; t++)
        for (auto b: grid->nodes) {
            Bus* bus = (Bus*) b;
            Constraint KCL_P("KCL_P"+bus->_name+ "time_" + to_string(t));
            Constraint KCL_Q("KCL_Q"+bus->_name+ "time_" + to_string(t));

            /* Power Conservation */
            KCL_P  = sum(Pf_from.in_at(b->get_out(), t)) + sum(Pf_to.in_at(b->get_in(), t)) + bus->pl()- sum(Pg.in_at(bus->_gen, t));
            KCL_Q  = sum(Qf_from.in_at(b->get_out(), t)) + sum(Qf_to.in_at(b->get_in(), t)) + bus->ql()- sum(Qg.in_at(bus->_gen, t));

            /* Shunts */
            //KCL_P +=  bus->gs()*Wii(bus->_name);
            //KCL_Q -=  bus->bs()*Wii(bus->_name);

            ACUC.add_constraint(KCL_P = 0);
            ACUC.add_constraint(KCL_Q = 0);
        }
    //AC Power Flow.
    Constraint Flow_P_From("Flow_P_From");
    Flow_P_From += Pf_from.in(grid->arcs, T);
    Flow_P_From -= grid->g_ff.in(grid->arcs, T)*Wii.from(grid->arcs, T);
    Flow_P_From -= grid->g_ft.in(grid->arcs, T)*R_Wij.in(grid->arcs, T);
    Flow_P_From -= grid->b_ft.in(grid->arcs, T)*Im_Wij.in(grid->arcs, T);
    Flow_P_From = 0;
    ACUC.add_constraint(Flow_P_From);

    Constraint Flow_P_To("Flow_P_To");
    Flow_P_To += Pf_to.in(grid->arcs, T);
    Flow_P_To -= grid->g_tt.in(grid->arcs, T)*Wii.to(grid->arcs, T);
    Flow_P_To -= grid->g_tf.in(grid->arcs, T)*R_Wij.in(grid->arcs, T);
    Flow_P_To += grid->b_tf.in(grid->arcs, T)*Im_Wij.in(grid->arcs, T);
    Flow_P_To = 0;
    ACUC.add_constraint(Flow_P_To);

    Constraint Flow_Q_From("Flow_Q_From");
    Flow_Q_From += Qf_from.in(grid->arcs, T);
    Flow_Q_From += grid->b_ff.in(grid->arcs, T)*Wii.from(grid->arcs, T);
    Flow_Q_From += grid->b_ft.in(grid->arcs, T)*R_Wij.in(grid->arcs, T);
    Flow_Q_From += grid->g_ft.in(grid->arcs, T)*Im_Wij.in(grid->arcs, T);
    Flow_Q_From = 0;
    ACUC.add_constraint(Flow_Q_From);

    Constraint Flow_Q_To("Flow_Q_To");
    Flow_Q_To += Qf_to.in(grid->arcs, T);
    Flow_Q_To += grid->b_tt.in(grid->arcs, T)*Wii.to(grid->arcs, T);
    Flow_Q_To += grid->b_tf.in(grid->arcs, T)*R_Wij.in(grid->arcs, T);
    Flow_Q_To -= grid->g_tf.in(grid->arcs, T)*Im_Wij.in(grid->arcs, T);
    Flow_Q_To = 0;
    ACUC.add_constraint(Flow_Q_To);

    /* Phase Angle Bounds constraints */
    Constraint PAD_UB("PAD_UB");
    PAD_UB = Im_Wij.in(grid->arcs, T);
    PAD_UB -= (grid->tan_th_max).in(grid->arcs, T)*R_Wij.in(grid->arcs, T);
    ACUC.add_constraint(PAD_UB <= 0);

    Constraint PAD_LB("PAD_LB");
    PAD_LB =  Im_Wij.in(grid->arcs, T);
    PAD_LB -= grid->tan_th_min.in(grid->arcs, T)*R_Wij.in(grid->arcs, T);
    ACUC.add_constraint(PAD_LB >= 0);

    /* Thermal Limit Constraints */
    Constraint Thermal_Limit_from("Thermal_Limit_from");
    Thermal_Limit_from += power(Pf_from.in(grid->arcs, T),  2) + power(Qf_from.in(grid->arcs, T), 2);
    Thermal_Limit_from -= power(grid->S_max.in(grid->arcs, T), 2);
    ACUC.add_constraint(Thermal_Limit_from <= 0);

    Constraint Thermal_Limit_to("Thermal_Limit_to");
    Thermal_Limit_to += power(Pf_to.in(grid->arcs, T), 2) + power(Qf_to.in(grid->arcs, T), 2);
    Thermal_Limit_to -= power(grid->S_max.in(grid->arcs, T),2);
    ACUC.add_constraint(Thermal_Limit_to <= 0);

    /* Commitment constraints */
    // Inter-temporal constraints
    for (int t = 1; t < T; t++){
        Constraint MC1("MC1_"+ to_string(t));
        Constraint MC2("MC2_"+ to_string(t));
        MC1 = On_off.in_at(grid->gens, 1)- On_off.in_at(grid->gens, 0)-  Start_up.in_at(grid->gens, 1);
        MC2 = On_off.in_at(grid->gens, t-1) - On_off.in_at(grid->gens, t) - Shut_down.in_at(grid->gens, t);
        ACUC.add_constraint(MC1 <= 0);
        ACUC.add_constraint(MC2 <= 0);
    }
    
    // Min-up constraints
    for (int t = 1; t < T; t++) {
        Constraint Min_up1("Min_up1_"+ to_string(t));
        Min_up1 = On_off.in_at(grid->gens, t) - On_off.in_at(grid->gens, t-1) - Start_up.in_at(grid->gens, t) + Shut_down.in_at(grid->gens, t);
        ACUC.add_constraint(Min_up1 = 0);
    }

    /* Resolve it! */
    solver OPF(ACUC,ipopt);
    OPF.run();
    return 0;
}
