//
// Created by kbestuzheva on 12/11/17.
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


/* main */
int main (int argc, char * argv[]) {
    int output = 0;
    bool relax = false, sdp_cuts = true;
    size_t num_bags = 0;
    string num_bags_s = "100";
    string solver_str = "ipopt";
    string sdp_cuts_s = "yes";
    string lazy = "no";
    bool lazy_bool = true;
    SolverType solv_type = ipopt;
    double tol = 1e-6;
    string mehrotra = "no";

    string fname = string(prj_dir)+"/data_sets/Power/nesta_case5_pjm.m";

    // create a OptionParser with options
    op::OptionParser opt;
    opt.add_option("h", "help",
                   "shows option help"); // no default value means boolean options, which default value is false
    opt.add_option("f", "file", "Input file name", fname);
    opt.add_option("s", "solver", "Solvers: ipopt/cplex/gurobi, default = ipopt", solver_str);
    opt.add_option("b", "numbags", "Number of bags per iteration", num_bags_s);
    opt.add_option("lz", "lazy", "Generate 3d SDP cuts in a lazy fashion, default = no", lazy);
    // parse the options and verify that all went well. If not, errors and help will be shown
    bool correct_parsing = opt.parse_options(argc, argv);

    if (!correct_parsing) {
        return EXIT_FAILURE;
    }

    fname = opt["f"];
    bool has_help = op::str2bool(opt["h"]);
    if (has_help) {
        opt.show_help();
        exit(0);
    }
    solver_str = opt["s"];
    if (solver_str.compare("gurobi")==0) {
        solv_type = gurobi;
    }
    else if(solver_str.compare("cplex")==0) {
        solv_type = cplex;
    }else if(solver_str.compare("Mosek")==0) {
        solv_type = mosek;
    }
    lazy = opt["lz"];
    if (lazy.compare("no")==0) {
        lazy_bool = false;
    }
    
    num_bags = atoi(opt["b"].c_str());

    cout << "\nnum bags = " << num_bags << endl;

    double total_time_start = get_wall_time();
    PowerNet grid;
    grid.readgrid(fname);

    grid.get_tree_decomp_bags(false,true);
//    grid.update_net();

    // Grid Parameters
    auto bus_pairs = grid.get_bus_pairs();
    auto bus_pairs_chord = grid.get_bus_pairs_chord();
    auto nb_lines = grid.get_nb_active_arcs();
    auto nb_buses = grid.get_nb_active_nodes();
    DebugOff("nb gens = " << nb_gen << endl);
    DebugOff("nb lines = " << nb_lines << endl);
    DebugOff("nb buses = " << nb_buses << endl);
    DebugOff("nb bus_pairs = " << nb_bus_pairs_chord << endl);


    double upper_bound = grid.solve_acopf();


    /** Build model */
    Model SDP("SDP Model");

    /** Variables */
    /* power generation variables */
    var<double> Pg("Pg", grid.pg_min, grid.pg_max);
    var<double> Qg ("Qg", grid.qg_min, grid.qg_max);
    SDP.add(Pg.in(grid.gens));
    SDP.add(Qg.in(grid.gens));
    /* power flow variables */
    var<double> Pf_from("Pf_from", grid.S_max);
    var<double> Qf_from("Qf_from", grid.S_max);
    var<double> Pf_to("Pf_to", grid.S_max);
    var<double> Qf_to("Qf_to", grid.S_max);
    SDP.add(Pf_from.in(grid.arcs));
    SDP.add(Qf_from.in(grid.arcs));
    SDP.add(Pf_to.in(grid.arcs));
    SDP.add(Qf_to.in(grid.arcs));

    /* Real part of Wij = ViVj */
    var<double>  R_Wij("R_Wij", grid.wr_min, grid.wr_max);
    /* Imaginary part of Wij = ViVj */
    var<double>  Im_Wij("Im_Wij", grid.wi_min, grid.wi_max);
    /* Magnitude of Wii = Vi^2 */
    var<double>  Wii("Wii", grid.w_min, grid.w_max);
    SDP.add(Wii.in(grid.nodes));

    if (sdp_cuts) {
        SDP.add(R_Wij.in(bus_pairs_chord));
        SDP.add(Im_Wij.in(bus_pairs_chord));
    }
    else{
        SDP.add(R_Wij.in(bus_pairs));
        SDP.add(Im_Wij.in(bus_pairs));
    }
    /* Initialize variables */
    R_Wij.initialize_all(1.0);
    Wii.initialize_all(1.001);

    /**  Objective */
    auto obj = product(grid.c1, Pg) + product(grid.c2, power(Pg,2)) + sum(grid.c0);
    SDP.min(obj.in(grid.gens));

    Constraint obj_cstr("obj_ub");
    obj_cstr += SDP._obj - 1.05*upper_bound;
    SDP.add(obj_cstr <= 0);

    
    /** Constraints */
    if(grid.add_3d_nlin && sdp_cuts) {
        auto bag_size = grid._bags.size();
        DebugOn("\nNum of bags = " << bag_size << endl);
        DebugOn("Adding 3d determinant polynomial cuts\n");
        auto R_Wij_ = R_Wij.pairs_in_bags(grid._bags, 3);
        auto Im_Wij_ = Im_Wij.pairs_in_bags(grid._bags, 3);
        auto Wii_ = Wii.in_bags(grid._bags, 3);
        
        Constraint SDP3("SDP_3D");
        SDP3 = 2 * R_Wij_[0] * (R_Wij_[1] * R_Wij_[2] + Im_Wij_[1] * Im_Wij_[2]);
        SDP3 -= 2 * Im_Wij_[0] * (R_Wij_[2] * Im_Wij_[1] - Im_Wij_[2] * R_Wij_[1]);
        SDP3 -= (power(R_Wij_[0], 2) + power(Im_Wij_[0], 2)) * Wii_[2];
        SDP3 -= (power(R_Wij_[1], 2) + power(Im_Wij_[1], 2)) * Wii_[0];
        SDP3 -= (power(R_Wij_[2], 2) + power(Im_Wij_[2], 2)) * Wii_[1];
        SDP3 += Wii_[0] * Wii_[1] * Wii_[2];
        if (lazy_bool) {
            SDP.add_lazy(SDP3 >= 0);
        }
        else {
            SDP.add(SDP3 >= 0);
            DebugOn("Number of 3d determinant cuts = " << SDP3.get_nb_instances() << endl);
        }
    }
    /* Second-order cone constraints */
    Constraint SOC("SOC");
    SOC = power(R_Wij, 2) + power(Im_Wij, 2) - Wii.from()*Wii.to();
    if (sdp_cuts) {
        SDP.add(SOC.in(bus_pairs_chord) <= 0);
    }
    else{
        SDP.add(SOC.in(bus_pairs) <= 0);
    }
    

    
    /* Flow conservation */
    Constraint KCL_P("KCL_P");
    KCL_P  = sum(Pf_from.out_arcs()) + sum(Pf_to.in_arcs()) + grid.pl - sum(Pg.in_gens()) + grid.gs*Wii;
    SDP.add(KCL_P.in(grid.nodes) == 0);
    
    Constraint KCL_Q("KCL_Q");
    KCL_Q  = sum(Qf_from.out_arcs()) + sum(Qf_to.in_arcs()) + grid.ql - sum(Qg.in_gens()) - grid.bs*Wii;
    SDP.add(KCL_Q.in(grid.nodes) == 0);
    
    /* AC Power Flow */
    Constraint Flow_P_From("Flow_P_From");
    Flow_P_From = Pf_from - (grid.g_ff*Wii.from() + grid.g_ft*R_Wij.in_pairs() + grid.b_ft*Im_Wij.in_pairs());
    SDP.add(Flow_P_From.in(grid.arcs) == 0);
    
    Constraint Flow_P_To("Flow_P_To");
    Flow_P_To = Pf_to - (grid.g_tt*Wii.to() + grid.g_tf*R_Wij.in_pairs() - grid.b_tf*Im_Wij.in_pairs());
    SDP.add(Flow_P_To.in(grid.arcs) == 0);
    
    Constraint Flow_Q_From("Flow_Q_From");
    Flow_Q_From = Qf_from - (grid.g_ft*Im_Wij.in_pairs() - grid.b_ff*Wii.from() - grid.b_ft*R_Wij.in_pairs());
    SDP.add(Flow_Q_From.in(grid.arcs) == 0);
    
    Constraint Flow_Q_To("Flow_Q_To");
    Flow_Q_To = Qf_to + (grid.b_tt*Wii.to() + grid.b_tf*R_Wij.in_pairs() + grid.g_tf*Im_Wij.in_pairs());
    SDP.add(Flow_Q_To.in(grid.arcs) == 0);
    

    /* Phase Angle Bounds constraints */

    Constraint PAD_UB("PAD_UB");
    PAD_UB = Im_Wij;
    PAD_UB <= grid.tan_th_max*R_Wij;
    SDP.add_lazy(PAD_UB.in(bus_pairs));
//    SDP.add(PAD_UB.in(bus_pairs));

    Constraint PAD_LB("PAD_LB");
    PAD_LB =  Im_Wij;
    PAD_LB >= grid.tan_th_min*R_Wij;
    SDP.add_lazy(PAD_LB.in(bus_pairs));
//    SDP.add(PAD_LB.in(bus_pairs));

    /* Thermal Limit Constraints */
    Constraint Thermal_Limit_from("Thermal_Limit_from");
    Thermal_Limit_from = power(Pf_from, 2) + power(Qf_from, 2);
    Thermal_Limit_from <= power(grid.S_max,2);
    SDP.add(Thermal_Limit_from.in(grid.arcs));
    
    
    Constraint Thermal_Limit_to("Thermal_Limit_to");
    Thermal_Limit_to = power(Pf_to, 2) + power(Qf_to, 2);
    Thermal_Limit_to <= power(grid.S_max,2);
    SDP.add(Thermal_Limit_to.in(grid.arcs));
    

    /* Lifted Nonlinear Cuts */
    Constraint LNC1("LNC1");
    LNC1 = (grid.v_min.from()+grid.v_max.from())*(grid.v_min.to()+grid.v_max.to())*(sin(0.5*(grid.th_max+grid.th_min))*Im_Wij + cos(0.5*(grid.th_max+grid.th_min))*R_Wij);
    LNC1 -= grid.v_max.to()*cos(0.5*(grid.th_max-grid.th_min))*(grid.v_min.to()+grid.v_max.to())*Wii.from();
    LNC1 -= grid.v_max.from()*cos(0.5*(grid.th_max-grid.th_min))*(grid.v_min.from()+grid.v_max.from())*Wii.to();
    LNC1 -= grid.v_max.from()*grid.v_max.to()*cos(0.5*(grid.th_max-grid.th_min))*(grid.v_min.from()*grid.v_min.to() - grid.v_max.from()*grid.v_max.to());
//    SDP.add(LNC1.in(bus_pairs) >= 0);
    SDP.add_lazy(LNC1.in(bus_pairs) >= 0);


    Constraint LNC2("LNC2");
    LNC2 = (grid.v_min.from()+grid.v_max.from())*(grid.v_min.to()+grid.v_max.to())*(sin(0.5*(grid.th_max+grid.th_min))*Im_Wij + cos(0.5*(grid.th_max+grid.th_min))*R_Wij);
    LNC2 -= grid.v_min.to()*cos(0.5*(grid.th_max-grid.th_min))*(grid.v_min.to()+grid.v_max.to())*Wii.from();
    LNC2 -= grid.v_min.from()*cos(0.5*(grid.th_max-grid.th_min))*(grid.v_min.from()+grid.v_max.from())*Wii.to();
    LNC2 += grid.v_min.from()*grid.v_min.to()*cos(0.5*(grid.th_max-grid.th_min))*(grid.v_min.from()*grid.v_min.to() - grid.v_max.from()*grid.v_max.to());
//    SDP.add(LNC2.in(bus_pairs) >= 0);
    SDP.add_lazy(LNC2.in(bus_pairs) >= 0);



    total_time_start = get_wall_time();
    /* Solver selection */
    solver SDPOPF(SDP,solv_type);
    double solver_time_start = get_wall_time();
    SDPOPF.run(output = 5, relax = false, tol = 1e-6, 1e-6, "mumps", mehrotra = "no");
    SDP.print();
  double gap = 100*(upper_bound - SDP._obj_val)/upper_bound;
    double solver_time_end = get_wall_time();
    double total_time_end = get_wall_time();
    auto solve_time = solver_time_end - solver_time_start;
    auto total_time = total_time_end - total_time_start;
    string out = "\nDATA_OPF, " + grid._name + ", " + to_string(nb_buses) + ", " + to_string(nb_lines) +", " + to_string(SDP._obj_val) + ", " + to_string(-numeric_limits<double>::infinity()) + ", " + to_string(solve_time) + ", LocalOptimal, " + to_string(total_time);
    DebugOn(out <<endl);
    DebugOn("Final Gap = " << to_string(gap) << "%."<<endl);
    DebugOn("Upper bound = " << to_string(upper_bound) << "."<<endl);
    DebugOn("Lower bound = " << to_string(SDP._obj_val) << "."<<endl);
    DebugOn("\nResults: " << grid._name << " " << to_string(SDP._obj_val) << " " << to_string(total_time)<<endl);
//    DebugOn("\nTime in nfp: " << time_in_all_nfp << endl);
    return 0;
}
