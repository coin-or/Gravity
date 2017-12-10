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
#include <optionParser.hpp>

using namespace std;
using namespace gravity;

int main (int argc, char * argv[])
{
    int output = 0;
    bool relax = false, use_cplex = false;
    double tol = 1e-6;
    string mehrotra = "no";
//    string fname = "/users/hh/Dropbox/Work/Dev/pglib-opf/sad/pglib_opf_case6470_rte__sad.m";
    string fname = "../../data_sets/Power/nesta_case5_pjm.m";
    
    // create a OptionParser with options
    op::OptionParser opt;
    opt.add_option("h", "help", "shows option help"); // no default value means boolean options, which default value is false
    opt.add_option("f", "file", "Input file name", fname );
    
    // parse the options and verify that all went well. If not, errors and help will be shown
    bool correct_parsing = opt.parse_options(argc, argv);
    
    if(!correct_parsing){
        return EXIT_FAILURE;
    }
    
    fname = opt["f"];
    bool has_help = op::str2bool(opt["h"]);
    if(has_help) {
        opt.show_help();
        exit(0);
    }
    double total_time_start = get_wall_time();
    PowerNet* grid = new PowerNet();
    grid->readgrid(fname.c_str());
    
    // Grid Parameters
    auto bus_pairs = grid->get_bus_pairs();
    auto nb_bus_pairs = grid->get_nb_active_bus_pairs();
    auto nb_gen = grid->get_nb_active_gens();
    auto nb_lines = grid->get_nb_active_arcs();
    auto nb_buses = grid->get_nb_active_nodes();
    DebugOff("nb gens = " << nb_gen << endl);
    DebugOff("nb lines = " << nb_lines << endl);
    DebugOff("nb buses = " << nb_buses << endl);
    DebugOff("nb bus_pairs = " << nb_bus_pairs << endl);
    
    /** Build model */
    Model SOCP("SOCP Model");
    
    /** Variables */
    /* power generation variables */
    var<Real> Pg("Pg", grid->pg_min.in(grid->gens), grid->pg_max.in(grid->gens));
    var<Real> Qg ("Qg", grid->qg_min.in(grid->gens), grid->qg_max.in(grid->gens));
    SOCP.add_var(Pg^(nb_gen));
    SOCP.add_var(Qg^(nb_gen));
    
    /* power flow variables */
    var<Real> Pf_from("Pf_from", grid->S_max.in(grid->arcs));
    var<Real> Qf_from("Qf_from", grid->S_max.in(grid->arcs));
    var<Real> Pf_to("Pf_to", grid->S_max.in(grid->arcs));
    var<Real> Qf_to("Qf_to", grid->S_max.in(grid->arcs));
    SOCP.add_var(Pf_from^(nb_lines));
    SOCP.add_var(Qf_from^(nb_lines));
    SOCP.add_var(Pf_to^(nb_lines));
    SOCP.add_var(Qf_to^(nb_lines));
    
    /* Real part of Wij = ViVj */
    var<Real>  R_Wij("R_Wij", grid->wr_min.in(bus_pairs), grid->wr_max.in(bus_pairs));
    /* Imaginary part of Wij = ViVj */
    var<Real>  Im_Wij("Im_Wij", grid->wi_min.in(bus_pairs), grid->wi_max.in(bus_pairs));
    /* Magnitude of Wii = Vi^2 */
    var<Real>  Wii("Wii", grid->w_min.in(grid->nodes), grid->w_max.in(grid->nodes));
    SOCP.add_var(Wii^nb_buses);
    SOCP.add_var(R_Wij^nb_bus_pairs);
    SOCP.add_var(Im_Wij^nb_bus_pairs);
    
    /* Initialize variables */
    R_Wij.initialize_all(1.0);
    Wii.initialize_all(1.001);
    
    /**  Objective */
    func_ obj;
    for (auto g:grid->gens) {
        if (g->_active) {
            obj += grid->c1(g->_name)*Pg(g->_name) + grid->c2(g->_name)*Pg(g->_name)*Pg(g->_name) + grid->c0(g->_name);
        }
    }
    SOCP.min(obj);
    
    
    /** Constraints */
    /* Second-order cone constraints */
    Constraint SOC("SOC");
    SOC =  power(R_Wij, 2) + power(Im_Wij, 2) - Wii.from()*Wii.to() ;
    SOCP.add_constraint(SOC.in(bus_pairs) <= 0);

    /* Flow conservation */
    for (auto b: grid->nodes) {
        Bus* bus = (Bus*) b;
        Constraint KCL_P("KCL_P"+bus->_name);
        Constraint KCL_Q("KCL_Q"+bus->_name);
        KCL_P  = sum(Pf_from.in(b->get_out())) + sum(Pf_to.in(b->get_in())) + bus->pl() - sum(Pg.in(bus->_gen));
        KCL_Q  = sum(Qf_from.in(b->get_out())) + sum(Qf_to.in(b->get_in())) + bus->ql() - sum(Qg.in(bus->_gen));

        /* Shunts */
        KCL_P +=  bus->gs()*(Wii(bus->_name));
        KCL_Q -=  bus->bs()*(Wii(bus->_name));
        
        SOCP.add_constraint(KCL_P = 0);
        SOCP.add_constraint(KCL_Q = 0);
    }

    /* AC Power Flow */
    Constraint Flow_P_From("Flow_P_From");
    Flow_P_From += Pf_from;
    Flow_P_From -= grid->g_ff*Wii.from();
    Flow_P_From -= grid->g_ft*R_Wij.in_pairs();
    Flow_P_From -= grid->b_ft*Im_Wij.in_pairs();
    SOCP.add_constraint(Flow_P_From.in(grid->arcs) = 0);
    
    Constraint Flow_P_To("Flow_P_To");
    Flow_P_To += Pf_to;
    Flow_P_To -= grid->g_tt*Wii.to();
    Flow_P_To -= grid->g_tf*R_Wij.in_pairs();
    Flow_P_To += grid->b_tf*Im_Wij.in_pairs();
    SOCP.add_constraint(Flow_P_To.in(grid->arcs) = 0);
    
    Constraint Flow_Q_From("Flow_Q_From");
    Flow_Q_From += Qf_from;
    Flow_Q_From += grid->b_ff*Wii.from();
    Flow_Q_From += grid->b_ft*R_Wij.in_pairs();
    Flow_Q_From -= grid->g_ft*Im_Wij.in_pairs();
    SOCP.add_constraint(Flow_Q_From.in(grid->arcs) = 0);

    Constraint Flow_Q_To("Flow_Q_To");
    Flow_Q_To += Qf_to;
    Flow_Q_To += grid->b_tt*Wii.to();
    Flow_Q_To += grid->b_tf*R_Wij.in_pairs();
    Flow_Q_To += grid->g_tf*Im_Wij.in_pairs();
    SOCP.add_constraint(Flow_Q_To.in(grid->arcs) = 0);

    /* Phase Angle Bounds constraints */
    Constraint PAD_UB("PAD_UB");
    PAD_UB = Im_Wij;
    PAD_UB -= (grid->tan_th_max)*R_Wij;
    SOCP.add_constraint(PAD_UB.in(bus_pairs) <= 0);
    
    Constraint PAD_LB("PAD_LB");
    PAD_LB =  Im_Wij;
    PAD_LB -= grid->tan_th_min*R_Wij;
    SOCP.add_constraint(PAD_LB.in(bus_pairs) >= 0);
    
    /* Lifted Nonlinear Cuts */
    Constraint LNC1("LNC1");
//    LNC1 += (grid->v_min.from()+grid->v_max.from())*(grid->v_min.to()+grid->v_max.to())*(sin(0.5*(grid->th_max+grid->th_min))*Im_Wij + cos(0.5*(grid->th_max+grid->th_min))*R_Wij);
//    LNC1 -= grid->v_max.to()*cos(0.5*(grid->th_max-grid->th_min))*(grid->v_min.to()+grid->v_max.to())*Wii.from();
//    LNC1 -= grid->v_max.from()*cos(0.5*(grid->th_max-grid->th_min))*(grid->v_min.from()+grid->v_max.from())*Wii.to();
//    LNC1 -= grid->v_max.from()*grid->v_max.to()*cos(0.5*(grid->th_max-grid->th_min))*(grid->v_min.from()*grid->v_min.to() - grid->v_max.from()*grid->v_max.to());
    LNC1 += (grid->v_min.from()+grid->v_max.from())*(grid->v_min.to()+grid->v_max.to())*(grid->sphi*Im_Wij + grid->cphi*R_Wij);
    LNC1 -= grid->v_max.to()*grid->cos_d*(grid->v_min.to()+grid->v_max.to())*Wii.from();
    LNC1 -= grid->v_max.from()*grid->cos_d*(grid->v_min.from()+grid->v_max.from())*Wii.to();
    LNC1 -= grid->v_max.from()*grid->v_max.to()*grid->cos_d*(grid->v_min.from()*grid->v_min.to() - grid->v_max.from()*grid->v_max.to());
    SOCP.add_constraint(LNC1.in(bus_pairs) >= 0);
    
    Constraint LNC2("LNC2");
//    LNC2 += (grid->v_min.from()+grid->v_max.from())*(grid->v_min.to()+grid->v_max.to())*(sin(0.5*(grid->th_max+grid->th_min))*Im_Wij + cos(0.5*(grid->th_max+grid->th_min))*R_Wij);
//    LNC2 -= grid->v_min.to()*cos(0.5*(grid->th_max-grid->th_min))*(grid->v_min.to()+grid->v_max.to())*Wii.from();
//    LNC2 -= grid->v_min.from()*cos(0.5*(grid->th_max-grid->th_min))*(grid->v_min.from()+grid->v_max.from())*Wii.to();
//    LNC2 += grid->v_min.from()*grid->v_min.to()*cos(0.5*(grid->th_max-grid->th_min))*(grid->v_min.from()*grid->v_min.to() - grid->v_max.from()*grid->v_max.to());
    LNC2 += (grid->v_min.from()+grid->v_max.from())*(grid->v_min.to()+grid->v_max.to())*(grid->sphi*Im_Wij + grid->cphi*R_Wij);
    LNC2 -= grid->v_min.to()*grid->cos_d*(grid->v_min.to()+grid->v_max.to())*Wii.from();
    LNC2 -= grid->v_min.from()*grid->cos_d*(grid->v_min.from()+grid->v_max.from())*Wii.to();
    LNC2 += grid->v_min.from()*grid->v_min.to()*grid->cos_d*(grid->v_min.from()*grid->v_min.to() - grid->v_max.from()*grid->v_max.to());
    SOCP.add_constraint(LNC2.in(bus_pairs) >= 0);
    
    /* Thermal Limit Constraints */
    Constraint Thermal_Limit_from("Thermal_Limit_from");
    Thermal_Limit_from += power(Pf_from, 2) + power(Qf_from, 2);
    Thermal_Limit_from -= power(grid->S_max,2);
    SOCP.add_constraint(Thermal_Limit_from.in(grid->arcs) <= 0);
    
    Constraint Thermal_Limit_to("Thermal_Limit_to");
    Thermal_Limit_to += power(Pf_to, 2) + power(Qf_to, 2);
    Thermal_Limit_to -= power(grid->S_max,2);
    SOCP.add_constraint(Thermal_Limit_to.in(grid->arcs) <= 0);

    
    /* Solver selection */
//    if (use_cplex) {
//        solver SCOPF_CPX(SOCP, cplex);
//        SCOPF_CPX.run(output = 0, relax = false, tol = 1e-6);
//    }
//    else {
        solver SCOPF(SOCP,ipopt);
        double solver_time_start = get_wall_time();
        SCOPF.run(output = 0, relax = false, tol = 1e-6, "ma27", mehrotra = "no");
        double solver_time_end = get_wall_time();
        double total_time_end = get_wall_time();
        auto solve_time = solver_time_end - solver_time_start;
        auto total_time = total_time_end - total_time_start;
//    }
    string out = "DATA_OPF, " + grid->_name + ", " + to_string(nb_buses) + ", " + to_string(nb_lines) +", " + to_string(SOCP._obj_val) + ", " + to_string(-numeric_limits<double>::infinity()) + ", " + to_string(solve_time) + ", LocalOptimal, " + to_string(total_time);
    DebugOn(out <<endl);
    return 0;
}
