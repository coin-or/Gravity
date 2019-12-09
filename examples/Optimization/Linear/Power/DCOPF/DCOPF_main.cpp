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
#ifdef USE_OPT_PARSER
#include <optionParser.hpp>
#endif
using namespace std;
using namespace gravity;

int main (int argc, char * argv[])
{
    int output = 0;
    bool projected = false, use_cplex = false, use_gurobi = false;
    double tol = 1e-6;
    double solver_time_end, total_time_end, solve_time, total_time;
    string mehrotra = "no", log_level="0";
    string fname = string(prj_dir)+"/data_sets/Power/nesta_case5_pjm.m";
    
    string path = argv[0];
    string solver_str="ipopt";
    string proj_str="0";
    
#ifdef USE_OPT_PARSER
    /** Create a OptionParser with options */
    auto options = readOptions(argc, argv);
    options.add_option("f", "file", "Input file name", fname);
    options.add_option("p", "project", "Project the power flow variables", proj_str);
    options.add_option("s", "solver", "Solvers: ipopt/cplex/gurobi, default = ipopt", solver_str);
    
    /** Parse the options and verify that all went well. If not, errors and help will be shown */
    bool correct_parsing = options.parse_options(argc, argv);
    
    if(!correct_parsing){
        return EXIT_FAILURE;
    }
    
    output = op::str2int(options["l"]);
    
    fname = options["f"];
    bool has_help = op::str2bool(options["h"]);
    if(has_help) {
        options.show_help();
        exit(0);
    }
    solver_str = options["s"];
    if (solver_str.compare("gurobi")==0) {
        use_gurobi = true;
    }
    else if(solver_str.compare("cplex")==0) {
        use_cplex = true;
    }
    solver_str = options["p"];
    if (solver_str.compare("1")==0) {
        projected = true;
    }
#endif
    double total_time_start = get_wall_time();
    PowerNet grid;
    grid.readgrid(fname,false);
    
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
    auto pl = grid.pl.in(nodes);
    auto c1 = grid.c1.in(gens);
    auto c2 = grid.c2.in(gens);
    auto c0 = grid.c0.in(gens);
    auto gs = grid.gs.in(nodes);
    auto S_max = grid.S_max.in(arcs);
    auto b = grid.b.in(arcs);
    auto th_min = grid.th_min.in(bus_pairs);
    auto th_max = grid.th_max.in(bus_pairs);
    
    /** Declare model */
    Model<> DCOPF("DCOPF Model");
    
    if(projected) {/* Project out the power flow variables*/
        /** Variables */
        /* Power generation variables */
        var<> Pg("Pg", pg_min, pg_max);
        DCOPF.add(Pg.in(gens));
        
        
        /* Phase angle variables */
        var<> theta("𝛉");
        DCOPF.add(theta.in(nodes));
        
        /**  Objective */
        auto obj = product(c1,Pg) + product(c2,pow(Pg,2)) + sum(c0);
        DCOPF.min(obj);
        
        /** Constraints */
        
        /* REF BUS */
        Constraint<> Ref_Bus("Ref_Bus");
        Ref_Bus = theta(grid.ref_bus);
        DCOPF.add(Ref_Bus == 0);
        
        /* Flow conservation */
        Constraint<> KCL_P("KCL_P");
        KCL_P  = b.tr().in(in_arcs)*(theta.from(in_arcs)-theta.to(in_arcs)) - b.tr().in(out_arcs)*(theta.from(out_arcs)-theta.to(out_arcs)) + pl + gs - sum(Pg, gen_nodes);
        DCOPF.add(KCL_P.in(nodes) == 0);
        
        /* Phase Angle Bounds constraints */
        Constraint<> PAD_UB("PAD_UB");
        PAD_UB = theta.from(bus_pairs) - theta.to(bus_pairs);
        PAD_UB -= th_max;
        DCOPF.add(PAD_UB.in(bus_pairs) <= 0);
        Constraint<> PAD_LB("PAD_LB");
        PAD_LB = theta.from(bus_pairs) - theta.to(bus_pairs);
        PAD_LB -= th_min;
        DCOPF.add(PAD_LB.in(bus_pairs) >= 0);
        
        /* Line Limits constraints */
        Constraint<> Thermal_UB("Thermal_UB");
        Thermal_UB = b*(theta.to(arcs) - theta.from(arcs));
        Thermal_UB -= S_max;
        DCOPF.add(Thermal_UB.in(arcs) <= 0);
        Constraint<> Thermal_LB("Thermal_LB");
        Thermal_LB = b*(theta.to(arcs) - theta.from(arcs));
        Thermal_LB += S_max;
        DCOPF.add(Thermal_LB.in(arcs) >= 0);
    }
    else {
        /** Variables */
        /* Power generation variables */
        var<> Pg("Pg", pg_min, pg_max);
        DCOPF.add(Pg.in(gens));
        
        /* Power flow variables */
        var<> Pf("Pf", -1*S_max, S_max);
        DCOPF.add(Pf.in(arcs));
        
        /* Phase angle variables */
        var<> theta("𝛉");
        DCOPF.add(theta.in(nodes));
        
        /**  Objective */
        auto obj = product(c1,Pg) + product(c2,pow(Pg,2)) + sum(c0);
        DCOPF.min(obj);
        
        /** Constraints */
        
        /* REF BUS */
        Constraint<> Ref_Bus("Ref_Bus");
        Ref_Bus = theta(grid.ref_bus);
        DCOPF.add(Ref_Bus == 0);
        
        /* Flow conservation */
        Constraint<> KCL_P("KCL_P");
        KCL_P  = sum(Pf, out_arcs) - sum(Pf, in_arcs) + pl + gs - sum(Pg, gen_nodes);
        DCOPF.add(KCL_P.in(nodes) == 0);
        
        /* AC Power Flow */
        Constraint<> Flow_P("Flow_P");
        Flow_P = Pf + b*(theta.from(arcs) - theta.to(arcs));
        DCOPF.add(Flow_P.in(arcs) == 0);
        
        /* Phase Angle Bounds constraints */
        Constraint<> PAD_UB("PAD_UB");
        PAD_UB = theta.from(bus_pairs) - theta.to(bus_pairs);
        PAD_UB -= th_max;
        DCOPF.add(PAD_UB.in(bus_pairs) <= 0);
        Constraint<> PAD_LB("PAD_LB");
        PAD_LB = theta.from(bus_pairs) - theta.to(bus_pairs);
        PAD_LB -= th_min;
        DCOPF.add(PAD_LB.in(bus_pairs) >= 0);
    }
    /* Solver selection */
    if (use_cplex) {
        solver<> DCOPF_CPX(DCOPF, cplex);
        auto solver_time_start = get_wall_time();
        DCOPF_CPX.run(output = 0, tol = 1e-6);
        solver_time_end = get_wall_time();
        total_time_end = get_wall_time();
        solve_time = solver_time_end - solver_time_start;
        total_time = total_time_end - total_time_start;
    }
    //    else if (use_gurobi) {
    //        solver DCOPF_GRB(DCOPF, gurobi);
    //        auto solver_time_start = get_wall_time();
    //        DCOPF_GRB.run(output, relax = false, tol = 1e-6);
    //        solver_time_end = get_wall_time();
    //        total_time_end = get_wall_time();
    //        solve_time = solver_time_end - solver_time_start;
    //        total_time = total_time_end - total_time_start;
    //    }
    else {
        solver<> DCOPF_IPT(DCOPF,ipopt);
        auto solver_time_start = get_wall_time();
        DCOPF_IPT.run(output = 5, tol = 1e-6);
        /** Warm starting Cplex
         solver<> DCOPF_CPX(DCOPF, cplex);
         DCOPF_CPX.run(output = 5, tol = 1e-6);
         **/
        solver_time_end = get_wall_time();
        total_time_end = get_wall_time();
        solve_time = solver_time_end - solver_time_start;
        total_time = total_time_end - total_time_start;
    }
    /** Uncomment next line to print expanded model */
    /* DCOPF.print(); */
    string out = "DATA_OPF, " + grid._name + ", " + to_string(nb_buses) + ", " + to_string(nb_lines) +", " + to_string(DCOPF.get_obj_val()) + ", " + to_string(-numeric_limits<double>::infinity()) + ", " + to_string(solve_time) + ", GlobalOptimal, " + to_string(total_time);
    DebugOn(out <<endl);
    return 0;
}
