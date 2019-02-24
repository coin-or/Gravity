//
//  SOCOPF.cpp
//  Gravity
//
//  Created by Hassan Hijazi on 21 January 2018.
//
//
#include <stdio.h>
#include <PowerNet.h>
#include <gravity/solver.h>
#include <optionParser.hpp>

using namespace std;
using namespace gravity;

int main (int argc, char * argv[])
{
    int output = 0;
    bool use_cplex = false, use_gurobi = false;
    double tol = 1e-6;
    double solver_time_end, total_time_end, solve_time, total_time;
    string mehrotra = "no", log_level="0";
    string fname = string(prj_dir)+"/data_sets/Power/nesta_case5_pjm.m";
    
    string path = argv[0];
    string solver_str="ipopt";
    
    /** Create a OptionParser with options */
    op::OptionParser opt;
    opt.add_option("h", "help", "shows option help"); // no default value means boolean options, which default value is false
    opt.add_option("f", "file", "Input file name", fname);
    opt.add_option("l", "log", "Log level (def. 0)", log_level );
    opt.add_option("s", "solver", "Solvers: ipopt/cplex/gurobi, default = ipopt", solver_str);
    
    /** Parse the options and verify that all went well. If not, errors and help will be shown */
    bool correct_parsing = opt.parse_options(argc, argv);
    
    if(!correct_parsing){
        return EXIT_FAILURE;
    }
    
    output = op::str2int(opt["l"]);
    
    fname = opt["f"];
    bool has_help = op::str2bool(opt["h"]);
    if(has_help) {
        opt.show_help();
        exit(0);
    }
    solver_str = opt["s"];
    if (solver_str.compare("gurobi")==0) {
        use_gurobi = true;
    }
    else if(solver_str.compare("cplex")==0) {
        use_cplex = true;
    }
    double total_time_start = get_wall_time();
    PowerNet grid;
    grid.readgrid(fname);
    
    /* Grid Stats */
    auto bus_pairs = grid.get_bus_pairs();
    auto nb_gen = grid.get_nb_active_gens();
    auto nb_lines = grid.get_nb_active_arcs();
    auto nb_buses = grid.get_nb_active_nodes();
    DebugOn("nb active gens = " << nb_gen << endl);
    DebugOn("nb active lines = " << nb_lines << endl);
    DebugOn("nb active buses = " << nb_buses << endl);
    
    /** Sets */
    auto nodes = indices(grid.nodes);
    auto arcs = indices(grid.arcs);
    auto gens = indices(grid.gens);
    auto gen_nodes = grid.gens_per_node();
    auto out_arcs = grid.out_arcs_per_node();
    auto in_arcs = grid.in_arcs_per_node();
    
    /* Grid Parameters */
    auto pg_min = grid.pg_min.in(gens);
    auto pg_max = grid.pg_max.in(gens);
    auto qg_min = grid.qg_min.in(gens);
    auto qg_max = grid.qg_max.in(gens);
    auto c1 = grid.c1.in(gens);
    auto c2 = grid.c2.in(gens);
    auto c0 = grid.c0.in(gens);
    

    
    /** MODEL DECLARATION */
    Model<> SOCP("SCOPF Model");
    /** Variables */
    /* power generation variables */
    var<> Pg("Pg", pg_min, pg_max);
    var<> Qg ("Qg", qg_min, qg_max);
    SOCP.add(Pg.in(gens),Qg.in(gens));
    
    
    /* power flow variables */
    var<> Pf_from("Pf_from", -1*grid.S_max,grid.S_max);
    var<> Qf_from("Qf_from", -1*grid.S_max,grid.S_max);
    var<> Pf_to("Pf_to", -1*grid.S_max,grid.S_max);
    var<> Qf_to("Qf_to", -1*grid.S_max,grid.S_max);
    SOCP.add(Pf_from.in(arcs), Qf_from.in(arcs), Pf_to.in(arcs), Qf_to.in(arcs));
    
    /* Real part of Wij = ViVj */
    var<>  R_Wij("R_Wij", grid.wr_min, grid.wr_max);
    /* Imaginary part of Wij = ViVj */
    var<>  Im_Wij("Im_Wij", grid.wi_min, grid.wi_max);
    /* Magnitude of Wii = Vi^2 */
    var<>  Wii("Wii", grid.w_min, grid.w_max);
    SOCP.add(Wii.in(nodes));
    SOCP.add(R_Wij.in(bus_pairs));
    SOCP.add(Im_Wij.in(bus_pairs));
    
    /* Initialize variables */
    R_Wij.initialize_all(1.0);
    Wii.initialize_all(1.001);
    
    
    
    /**  Objective */
    auto obj = c1.tr()*Pg + c2.tr()*pow(Pg.vec(),2) + sum(c0);
    SOCP.min(obj);
    
    /** Constraints */
    /* Second-order cone constraints */
    Constraint<> SOC("SOC");
    SOC = pow(R_Wij, 2) + pow(Im_Wij, 2) - Wii.from()*Wii.to();
    SOCP.add(SOC.in(bus_pairs) <= 0);
    
    /* Flow conservation */
    Constraint<> KCL_P("KCL_P");
    KCL_P  = sum(Pf_from, out_arcs) + sum(Pf_to, in_arcs) + grid.pl - sum(Pg, gen_nodes) + grid.gs*Wii;
    SOCP.add(KCL_P.in(nodes) == 0);
    
    Constraint<> KCL_Q("KCL_Q");
    KCL_Q  = sum(Qf_from, out_arcs) + sum(Qf_to, in_arcs) + grid.ql - sum(Qg, gen_nodes) - grid.bs*Wii;
    SOCP.add(KCL_Q.in(nodes) == 0);
    
    /* AC Power Flow */
    Constraint<> Flow_P_From("Flow_P_From");
    Flow_P_From = Pf_from - (grid.g_ff*Wii.from() + grid.g_ft*R_Wij.in_pairs() + grid.b_ft*Im_Wij.in_pairs());
    SOCP.add(Flow_P_From.in(arcs) == 0);
    
    Constraint<> Flow_P_To("Flow_P_To");
    Flow_P_To = Pf_to - (grid.g_tt*Wii.to() + grid.g_tf*R_Wij.in_pairs() - grid.b_tf*Im_Wij.in_pairs());
    SOCP.add(Flow_P_To.in(arcs) == 0);
    
    Constraint<> Flow_Q_From("Flow_Q_From");
    Flow_Q_From = Qf_from - (grid.g_ft*Im_Wij.in_pairs() - grid.b_ff*Wii.from() - grid.b_ft*R_Wij.in_pairs());
    SOCP.add(Flow_Q_From.in(arcs) == 0);
    
    Constraint<> Flow_Q_To("Flow_Q_To");
    Flow_Q_To = Qf_to + (grid.b_tt*Wii.to() + grid.b_tf*R_Wij.in_pairs() + grid.g_tf*Im_Wij.in_pairs());
    SOCP.add(Flow_Q_To.in(arcs) == 0);
    
    /* Phase Angle Bounds constraints */
    Constraint<> PAD_UB("PAD_UB");
    PAD_UB = Im_Wij;
    PAD_UB <= grid.tan_th_max*R_Wij;
    SOCP.add(PAD_UB.in(arcs));
    
    Constraint<> PAD_LB("PAD_LB");
    PAD_LB =  Im_Wij;
    PAD_LB >= grid.tan_th_min*R_Wij;
    SOCP.add(PAD_LB.in(arcs));
    
    /* Thermal Limit Constraints */
    Constraint<> Thermal_Limit_from("Thermal_Limit_from");
    Thermal_Limit_from = pow(Pf_from, 2) + pow(Qf_from, 2);
    Thermal_Limit_from <= pow(grid.S_max,2);
    SOCP.add(Thermal_Limit_from.in(arcs));
    
    
    Constraint<> Thermal_Limit_to("Thermal_Limit_to");
    Thermal_Limit_to = pow(Pf_to, 2) + pow(Qf_to, 2);
    Thermal_Limit_to <= pow(grid.S_max,2);
    SOCP.add(Thermal_Limit_to.in(arcs));
    
    /* Lifted Nonlinear Cuts */
    Constraint<> LNC1("LNC1");
    LNC1 += (grid.v_min.from()+grid.v_max.from())*(grid.v_min.to()+grid.v_max.to())*(grid.sphi*Im_Wij + grid.cphi*R_Wij);
    LNC1 -= grid.v_max.to()*grid.cos_d*(grid.v_min.to()+grid.v_max.to())*Wii.from();
    LNC1 -= grid.v_max.from()*grid.cos_d*(grid.v_min.from()+grid.v_max.from())*Wii.to();
    LNC1 -= grid.v_max.from()*grid.v_max.to()*grid.cos_d*(grid.v_min.from()*grid.v_min.to() - grid.v_max.from()*grid.v_max.to());
    SOCP.add(LNC1.in(bus_pairs) >= 0);
    
    Constraint<> LNC2("LNC2");
    LNC2 += (grid.v_min.from()+grid.v_max.from())*(grid.v_min.to()+grid.v_max.to())*(grid.sphi*Im_Wij + grid.cphi*R_Wij);//Moveaway from from()
    LNC2 -= grid.v_min.to()*grid.cos_d*(grid.v_min.to()+grid.v_max.to())*Wii.from();
    LNC2 -= grid.v_min.from()*grid.cos_d*(grid.v_min.from()+grid.v_max.from())*Wii.to();
    LNC2 += grid.v_min.from()*grid.v_min.to()*grid.cos_d*(grid.v_min.from()*grid.v_min.to() - grid.v_max.from()*grid.v_max.to());
    SOCP.add(LNC2.in(bus_pairs) >= 0);
    
    /* Solver selection */
    /* TODO: declare only one solver and one set of time measurment functions for all solvers. */
//    if (use_cplex) {
//        solver<> SCOPF_CPX(SOCP, cplex);
//        auto solver_time_start = get_wall_time();
//        SCOPF_CPX.run(output = 0, relax = false, tol = 1e-6);
//        solver_time_end = get_wall_time();
//        total_time_end = get_wall_time();
//        solve_time = solver_time_end - solver_time_start;
//        total_time = total_time_end - total_time_start;
//    }
//    else if (use_gurobi) {
//        solver<> SCOPF_GRB(SOCP, gurobi);
//        auto solver_time_start = get_wall_time();
//        SCOPF_GRB.run(output, relax = false, tol = 1e-6);
//        solver_time_end = get_wall_time();
//        total_time_end = get_wall_time();
//        solve_time = solver_time_end - solver_time_start;
//        total_time = total_time_end - total_time_start;
//    }
//    else {
//        SOCP.print();
        solver<> SOCOPF(SOCP,ipopt);
        auto solver_time_start = get_wall_time();
        SOCOPF.run(output, tol=1e-6);
        solver_time_end = get_wall_time();
        total_time_end = get_wall_time();
        solve_time = solver_time_end - solver_time_start;
        total_time = total_time_end - total_time_start;
//    }
//    SOCP.print_solution(12);
    /** Uncomment next line to print expanded model */
    /* SOCP.print(); */
    string out = "DATA_OPF, " + grid._name + ", " + to_string(nb_buses) + ", " + to_string(nb_lines) +", " + to_string_with_precision(SOCP.get_obj_val(),10) + ", " + to_string(-numeric_limits<double>::infinity()) + ", " + to_string(solve_time) + ", GlobalOptimal, " + to_string(total_time);
    DebugOn(out <<endl);
    return 0;
}
