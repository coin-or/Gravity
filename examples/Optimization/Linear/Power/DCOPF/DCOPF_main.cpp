//
//  DCOPF.cpp
//  Gravity
//
//  Created by Hassan Hijazi on 19 Jan 18.
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
    bool relax = false, use_cplex = false, use_gurobi = false;
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
    
    /* Grid Parameters */
    auto bus_pairs = grid.get_bus_pairs();
    auto nb_bus_pairs = grid.get_nb_active_bus_pairs();
    auto nb_gen = grid.get_nb_active_gens();
    auto nb_lines = grid.get_nb_active_arcs();
    auto nb_buses = grid.get_nb_active_nodes();
    DebugOn("nb gens = " << nb_gen << endl);
    DebugOn("nb lines = " << nb_lines << endl);
    DebugOn("nb buses = " << nb_buses << endl);
    DebugOn("nb bus_pairs = " << nb_bus_pairs << endl);
    
    /** Declare model */
    Model DCOPF("DCOPF Model");
    
    /** Variables */
    /* Power generation variables */
    var<double> Pg("Pg", grid.pg_min, grid.pg_max);
    DCOPF.add(Pg.in(grid.gens));
    
    /* Power flow variables */
    var<double> Pf("Pf", grid.S_max);
    DCOPF.add(Pf.in(grid.arcs));
    
    /* Phase angle variables */
    var<double> theta("theta");
    DCOPF.add(theta.in(grid.nodes));
    
    /**  Objective */
    auto obj = product(grid.c1, Pg) + product(grid.c2, power(Pg,2)) + sum(grid.c0);
    DCOPF.min(obj.in(grid.gens));
    
    /** Constraints */
    
    /* REF BUS */
    Constraint Ref_Bus("Ref_Bus");
    Ref_Bus = theta(grid.ref_bus);
    DCOPF.add(Ref_Bus == 0);
    
    /* Flow conservation */
    Constraint KCL_P("KCL_P");
    KCL_P  = sum(Pf.out_arcs()) - sum(Pf.in_arcs()) + grid.pl + grid.gs - sum(Pg.in_gens());
    DCOPF.add(KCL_P.in(grid.nodes) == 0);
    
    /* AC Power Flow */
    Constraint Flow_P("Flow_P");
    Flow_P = Pf + grid.b*(theta.from() - theta.to());
    DCOPF.add(Flow_P.in(grid.arcs) == 0);

    /* Phase Angle Bounds constraints */
    Constraint PAD_UB("PAD_UB");
    PAD_UB = theta.from() - theta.to();
    PAD_UB -= grid.th_max;
    DCOPF.add(PAD_UB.in(bus_pairs) <= 0);
    Constraint PAD_LB("PAD_LB");
    PAD_LB = theta.from() - theta.to();
    PAD_LB -= grid.th_min;
    DCOPF.add(PAD_LB.in(bus_pairs) >= 0);
    DCOPF.print();
    /* Solver selection */
    if (use_cplex) {
        solver DCOPF_CPX(DCOPF, cplex);
        auto solver_time_start = get_wall_time();
        DCOPF_CPX.run(output = 0, relax = false, tol = 1e-6);
        solver_time_end = get_wall_time();
        total_time_end = get_wall_time();
        solve_time = solver_time_end - solver_time_start;
        total_time = total_time_end - total_time_start;
    }
    else if (use_gurobi) {
        solver DCOPF_GRB(DCOPF, gurobi);
        auto solver_time_start = get_wall_time();
        DCOPF_GRB.run(output, relax = false, tol = 1e-6);
        solver_time_end = get_wall_time();
        total_time_end = get_wall_time();
        solve_time = solver_time_end - solver_time_start;
        total_time = total_time_end - total_time_start;
    }
    else {
        solver DCOPF_IPT(DCOPF,ipopt);
        auto solver_time_start = get_wall_time();
        DCOPF_IPT.run(output, relax = false, tol = 1e-6, 0.01, "ma27", mehrotra = "no");
        solver_time_end = get_wall_time();
        total_time_end = get_wall_time();
        solve_time = solver_time_end - solver_time_start;
        total_time = total_time_end - total_time_start;
    }
    /** Uncomment next line to print expanded model */
    /* DCOPF.print(); */
    string out = "DATA_OPF, " + grid._name + ", " + to_string(nb_buses) + ", " + to_string(nb_lines) +", " + to_string(DCOPF._obj_val) + ", " + to_string(-numeric_limits<double>::infinity()) + ", " + to_string(solve_time) + ", GlobalOptimal, " + to_string(total_time);
    DebugOn(out <<endl);
    return 0;
}
