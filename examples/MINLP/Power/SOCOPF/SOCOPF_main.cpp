//
//  SOCOPF.cpp
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
        fname = "/Users/hh/Dropbox/Work/Dev/nesta-0.7.0/opf/nesta_case2383wp_mp.m";
//        fname = "../../data_sets/Power/nesta_case3_lmbd.m";
    }
    PowerNet* grid = new PowerNet();
    grid->readgrid(fname);
//    grid->get_tree_decomp_bags(true);
    
    // Grid Parameters
    auto bus_pairs = grid->get_bus_pairs();
    auto nb_bus_pairs = bus_pairs.size();
    auto nb_gen = grid->get_nb_active_gens();
    auto nb_lines = grid->get_nb_active_arcs();
    auto nb_buses = grid->get_nb_active_nodes();
    /** build model */
    Model SOCP("SOCP Model");

    /** Variables */
    // power generation
    var<Real> Pg("Pg", grid->pg_min.in(grid->gens), grid->pg_max.in(grid->gens));
    var<Real> Qg ("Qg", grid->qg_min.in(grid->gens), grid->qg_max.in(grid->gens));
    SOCP.add_var(Pg^(nb_gen));
    SOCP.add_var(Qg^(nb_gen));
    
    // power flow
    var<Real> Pf_from("Pf_from");
    var<Real> Qf_from("Qf_from");
    var<Real> Pf_to("Pf_to");
    var<Real> Qf_to("Qf_to");
    SOCP.add_var(Pf_from^(nb_lines));
    SOCP.add_var(Qf_from^(nb_lines));
    SOCP.add_var(Pf_to^(nb_lines));
    SOCP.add_var(Qf_to^(nb_lines));
    
    // Lifted variables.
    var<Real>  R_Wij("R_Wij", grid->wr_min.in(bus_pairs), grid->wr_max.in(bus_pairs)); // real part of Wij
    var<Real>  Im_Wij("Im_Wij", grid->wi_min.in(bus_pairs), grid->wi_max.in(bus_pairs)); // imaginary part of Wij.
//    var<Real>  R_Wij("R_Wij"); // real part of Wij
//    var<Real>  Im_Wij("Im_Wij"); // imaginary part of Wij.

    var<Real>  Wii("Wii", grid->w_min.in(grid->nodes), grid->w_max.in(grid->nodes));
//    var<Real>  Wii("Wii");
    SOCP.add_var(Wii^nb_buses);
    SOCP.add_var(R_Wij^nb_bus_pairs);
    SOCP.add_var(Im_Wij^nb_bus_pairs);
    R_Wij.initialize_all(1.0);
    Wii.initialize_all(1.001);
    /** Construct the objective function*/
    func_ obj;
    for (auto g:grid->gens) {
        if (g->_active) {
            obj += grid->c1(g->_name)*Pg(g->_name) + grid->c2(g->_name)*Pg(g->_name)*Pg(g->_name) + grid->c0(g->_name);
        }
    }
    //func_ obj = sum(grid->c0) +sum(grid->c1.in(grid->gens),Pg.in(grid->gens)) + sum(grid->c2.in(grid->gens), power(Pg.in(grid->gens),2));
    SOCP.set_objective(min(obj));
    
    /** Define constraints */
    /* SOCP constraints */
    Constraint SOC("SOC");
    SOC =  power(R_Wij.in(bus_pairs), 2) + power(Im_Wij.in(bus_pairs), 2) - Wii.from(bus_pairs)*Wii.to(bus_pairs) ;
    SOCP.add_constraint(SOC <= 0);
    
    //KCL
    for (auto b: grid->nodes) {
        Bus* bus = (Bus*) b;
        Constraint KCL_P("KCL_P"+bus->_name);
        Constraint KCL_Q("KCL_Q"+bus->_name);

        /* Power Conservation */
        KCL_P  = sum(Pf_from.in(b->get_out())) + sum(Pf_to.in(b->get_in())) + bus->pl()- sum(Pg.in(bus->_gen));
        KCL_Q  = sum(Qf_from.in(b->get_out())) + sum(Qf_to.in(b->get_in())) + bus->ql()- sum(Qg.in(bus->_gen));

        /* Shunts */
        KCL_P +=  bus->gs()*(Wii(bus->_name));
        KCL_Q +=  bus->bs()*(Wii(bus->_name));
        
        SOCP.add_constraint(KCL_P = 0);
        SOCP.add_constraint(KCL_Q = 0);
    }

    //AC Power Flow
    Constraint Flow_P_From("Flow_P_From");
    Flow_P_From += Pf_from.in(grid->arcs);
    Flow_P_From -= grid->g_ff.in(grid->arcs)*Wii.from(grid->arcs);
    Flow_P_From -= grid->g_ft.in(grid->arcs)*R_Wij.in_pairs(grid->arcs);
    Flow_P_From -= grid->b_ft.in(grid->arcs)*Im_Wij.in_pairs(grid->arcs);
    SOCP.add_constraint(Flow_P_From = 0);
    
    Constraint Flow_P_To("Flow_P_To");
    Flow_P_To += Pf_to.in(grid->arcs);
    Flow_P_To -= grid->g_tt.in(grid->arcs)*Wii.to(grid->arcs);
    Flow_P_To -= grid->g_tf.in(grid->arcs)*R_Wij.in_pairs(grid->arcs);
    Flow_P_To += grid->b_tf.in(grid->arcs)*Im_Wij.in_pairs(grid->arcs);
    SOCP.add_constraint(Flow_P_To = 0);
    
    Constraint Flow_Q_From("Flow_Q_From");
    Flow_Q_From += Qf_from.in(grid->arcs);
    Flow_Q_From += grid->b_ff.in(grid->arcs)*Wii.from(grid->arcs);
    Flow_Q_From += grid->b_ft.in(grid->arcs)*R_Wij.in_pairs(grid->arcs);
    Flow_Q_From -= grid->g_ft.in(grid->arcs)*Im_Wij.in_pairs(grid->arcs);
    SOCP.add_constraint(Flow_Q_From = 0);
    
    Constraint Flow_Q_To("Flow_Q_To");
    Flow_Q_To += Qf_to.in(grid->arcs);
    Flow_Q_To += grid->b_tt.in(grid->arcs)*Wii.to(grid->arcs);
    Flow_Q_To += grid->b_tf.in(grid->arcs)*R_Wij.in_pairs(grid->arcs);
    Flow_Q_To += grid->g_tf.in(grid->arcs)*Im_Wij.in_pairs(grid->arcs);
    SOCP.add_constraint(Flow_Q_To = 0);

    /* Phase Angle Bounds constraints */
    Constraint PAD_UB("PAD_UB");
    PAD_UB = Im_Wij.in(bus_pairs);
    PAD_UB -= (grid->tan_th_max).in(bus_pairs)*R_Wij.in(bus_pairs);
    SOCP.add_constraint(PAD_UB <= 0);
    
    Constraint PAD_LB("PAD_LB");
    PAD_LB =  Im_Wij.in(bus_pairs);
    PAD_LB -= grid->tan_th_min.in(bus_pairs)*R_Wij.in(bus_pairs);
    SOCP.add_constraint(PAD_LB >= 0);
    
    /* Thermal Limit Constraints */
    Constraint Thermal_Limit_from("Thermal_Limit_from");
    Thermal_Limit_from += power(Pf_from.in(grid->arcs), 2) + power(Qf_from.in(grid->arcs), 2);
    Thermal_Limit_from -= power(grid->S_max.in(grid->arcs),2);
    SOCP.add_constraint(Thermal_Limit_from <= 0);
    
    Constraint Thermal_Limit_to("Thermal_Limit_to");
    Thermal_Limit_to += power(Pf_to.in(grid->arcs), 2) + power(Qf_to.in(grid->arcs), 2);
    Thermal_Limit_to -= power(grid->S_max.in(grid->arcs),2);
    SOCP.add_constraint(Thermal_Limit_to <= 0);
    
    solver SCOPF(SOCP,ipopt);
    SCOPF.run();
    return 0;
}
