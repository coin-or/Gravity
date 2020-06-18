//
// RSCED.cpp
// Gravity
//
// Created by Avinash Madavan on 5/20/20.
//
//

#include <stdio.h>
#include <gravity/solver.h>
#include <PowerNet.h>
#ifdef USE_OPT_PARSER
#include <optionParser.hpp>
#endif

int main(int argc, char *argv[]) {
    int output = 0;
    bool projected = false, use_cplex = false, use_gurobi = false;
    double tol = 1e-6;
    double solver_time_end, total_time_end, solve_time, total_time;
    string mehrotra = "no", log_level = "0";
    string fname = string(prj_dir) + "/data_sets/Power/nesta_case5_pjm.m";
    
    string path = argv[0];
    string solver_str = "ipopt";
    string proj_str = "0";
    
    // TODO: properly include CVaR parameter
    double cvar_param = 0.;
    
#ifdef USE_OPT_PARSER
    /** Create a OptionParser with options */
    auto options = readOptions(argc, argv);
    options.add_option("f", "file", "Input file name", fname);
    options.add_option("p", "project", "Project the power flow variables", proj_str);
    options.add_option("s", "solver", "Solvers: ipopt/cplex/gurobi, default = ipopt", solver_str);
    
    /** Parse the options and verify that all went well. If not, errors and help will be shown */
    bool correct_parsing = options.parse_options(argc, argv);
    
    if (!correct_parsing) {
        return EXIT_FAILURE;
    }
    
    output = op::str2int(options["l"]);
    
    fname = options["f"];
    bool has_help = op::str2bool(options["h"]);
    if (has_help) {
        options.show_help();
        exit(0);
    }
    solver_str = options["s"];
    if (solver_str.compare("gurobi") == 0) {
        use_gurobi = true;
    } else if (solver_str.compare("cplex") == 0) {
        use_cplex = true;
    }
    solver_str = options["p"];
    if (solver_str.compare("1") == 0) {
        projected = true;
    }
#else
    if(argc==2){
        fname=argv[1];
    }
    else{
        fname=string(prj_dir)+"/data_sets/Power/nesta_case5_pjm.m";
    }
#endif
    
    double total_time_start = get_wall_time();
    PowerNet grid;
    grid.readgrid(fname, true);
    
    /* Grid Stats */
    auto node_pairs = grid.get_node_pairs();
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
    auto th_min = grid.th_min.in(node_pairs);
    auto th_max = grid.th_max.in(node_pairs);
    
    // TODO: Include value of lost load (voll)
    auto voll = c1;
    
    // TODO: DA/STE limits for line capacities S_max
    double DAL_multiplier = 1.5;
    double STE_multiplier = 1.2;
    
    // TODO: Line failure probabilities
    double failure_probability = 0.05;
    
    // TODO: Include CVaR parameter (alpha' = 1/(1 - alpha))
    double alpha_prime = 1. / (1. - cvar_param);
    
    // TODO: Include generation ramp capacity
    auto ramp_max = pg_max;
    
    /** Declare model */
    Model<> RSCED("RSCED Model");
    
    /** Variables */
    /* Power generation variables */
    var<> Pg("Pg", pg_min, pg_max);
    RSCED.add(Pg.in(gens));
    
    /* Phase angle variables */
    var<> theta("𝛉");
    RSCED.add(theta.in(nodes));
    
    /* CVaR reformulation variable */
    // y0, z0(?)
    var<> z("z");
    var<> y("y");
    RSCED.add(z);
    RSCED.add(y);
    z.add_lb_only(0);
    y.add_lb_only(0);
    
    /** Constraints */
    /* CVaR lower bound */
    Constraint<> OBJ_BOUND("OBJ_BOUND");
    OBJ_BOUND = y + z - product(c2, Pg);
    RSCED.add(OBJ_BOUND >= 0);
    
    /* REF BUS */
    Constraint<> Ref_Bus("Ref_Bus");
    Ref_Bus = theta(grid.ref_bus);
    RSCED.add(Ref_Bus == 0);
    
    /* Power balance constraint */
    Constraint<> KCL_P("KCL_P");
    KCL_P = b.tr().in(in_arcs) * (theta.from(in_arcs) - theta.to(in_arcs))
    - b.tr().in(out_arcs) * (theta.from(out_arcs) - theta.to(out_arcs)) + pl + gs - sum(Pg, gen_nodes);
    RSCED.add(KCL_P.in(nodes) == 0);
    KCL_P.print();
    
    /* Phase Angle Bounds constraints */
    Constraint<> PAD_UB("PAD_UB");
    PAD_UB = theta.from(node_pairs) - theta.to(node_pairs);
    PAD_UB -= th_max;
    RSCED.add(PAD_UB.in(node_pairs) <= 0);
    Constraint<> PAD_LB("PAD_LB");
    PAD_LB = theta.from(node_pairs) - theta.to(node_pairs);
    PAD_LB -= th_min;
    RSCED.add(PAD_LB.in(node_pairs) >= 0);
    
    /* Line Limits constraints */
    Constraint<> Thermal_UB("Thermal_UB");
    Thermal_UB = b * (theta.to(arcs) - theta.from(arcs));
    Thermal_UB -= S_max;
    RSCED.add(Thermal_UB.in(arcs) <= 0);
    Constraint<> Thermal_LB("Thermal_LB");
    Thermal_LB = b * (theta.to(arcs) - theta.from(arcs));
    Thermal_LB += S_max;
    RSCED.add(Thermal_LB.in(arcs) >= 0);
    
    /** Objective */
    auto obj = alpha_prime * (1 - nb_lines * failure_probability) * product(c1, Pg);
    
    /** Handling contingencies */
    std::vector<std::pair<std::string, std::pair<Arc *, Gen *>>> conting_lines;
    for (auto cont : grid.arcs)
        conting_lines.emplace_back("cont_"+to_string(cont->_id), std::make_pair(cont, nullptr));
    
    auto arcs_c = grid.get_conting_arcs(conting_lines);
    auto node_pairs_c = grid.get_node_pairs_cont(conting_lines);
    indices contingencies("conting");
    for (auto cont : conting_lines)
        contingencies.insert(cont.first);
    
    auto nodes_c = indices(contingencies, nodes);
    auto gens_c = indices(contingencies, gens);
    
    auto gen_nodes_c = grid.gens_per_node_cont(conting_lines, gens_c);
    auto out_arcs_c = grid.out_arcs_per_node_cont(conting_lines, arcs_c);
    auto in_arcs_c = grid.in_arcs_per_node_cont(conting_lines, arcs_c);
    
    std::cout << grid.arcs.size() << std::endl;
    grid.S_max.print_vals(3);
    auto S_max_c = grid.S_max.from_ith(1, arcs_c);
    auto b_c = grid.b.from_ith(1, arcs_c);
    
    /** Contingency variable declaration */
    /* Positive recourse generation */
    var<> dPg_p("dPg_p");
    RSCED.add(dPg_p.in(gens_c));
    
    /* Negative recourse generation */
    var<> dPg_n("dPg_n");
    RSCED.add(dPg_n.in(gens_c));
    
    /* Power flows */
    var<> Pf_c("Pf_c");
    RSCED.add(Pf_c.in(arcs_c));
    var<> Pf_recourse_c("Pf_recourse_c");
    RSCED.add(Pf_recourse_c.in(arcs_c));
    
    /* Phase angle variables */
    var<> theta_c("theta_c");
    RSCED.add(theta_c.in(nodes_c));
    var<> theta_recourse_c("theta_recourse_c");
    RSCED.add(theta_recourse_c.in(nodes_c));
    
    /* Load shed variables */
    var<> dD_c("dD_c");
    RSCED.add(dD_c.in(nodes_c));
    
    /* CVaR reformulation variable */
    var<> y_c("y_c");
    RSCED.add(y_c.in(contingencies));
    y_c.add_lb_only(0);
    
    /** Post-failure constraints */
    /* REF BUS */
    Constraint<> Ref_Bus_c("Ref_Bus_c");
    Ref_Bus_c = theta_c.in(grid.ref_bus);
    RSCED.add(Ref_Bus_c == 0);
    
    /* Power flow constraint */
    Constraint<> Flow_P_c("Flow_P_c");
    Flow_P_c = Pf_c + b_c * (theta_c.from(arcs_c) - theta_c.to(arcs_c));
    RSCED.add(Flow_P_c.in(arcs_c) == 0);
    
    /* Power balance constraint */
    Constraint<> KCL_P_c("KCL_P_c");
    KCL_P_c  = sum(Pf_c, out_arcs_c) - sum(Pf_c, in_arcs_c) + pl.from_ith(1,nodes_c) - sum(Pg.from_ith(1,gen_nodes_c));
    RSCED.add(KCL_P_c.in(nodes_c) == 0);
//    KCL_P_c.print();
    
    Constraint<> PAD_UB_c("PAD_UB_c");
    PAD_UB_c = theta_c.from_ith(0,node_pairs_c) - theta_c.in_ignore_ith(1,1,node_pairs_c);
    PAD_UB_c -= th_max.from_ith(1, node_pairs_c);
    RSCED.add(PAD_UB_c.in(node_pairs_c) <= 0);
    PAD_UB_c.print();
    /* Phase Angle Bounds constraints */
//    Constraint<> PAD_UB_c("PAD_UB_c");
//    PAD_UB_c = theta_c.from(nodes_c) - theta_c.to(nodes_c);
//    PAD_UB_c -= th_max;
//    RSCED.add(PAD_UB_c.in(nodes_c) <= 0);
//    Constraint<> PAD_LB_c("PAD_LB_c");
//    PAD_LB_c = theta_c.from(nodes_c) - theta_c.to(nodes_c);
//    PAD_LB_c -= th_min;
//    RSCED.add(PAD_LB_c.in(nodes_c) >= 0);
    
    /* Line Limits constraints */
    Constraint<> Thermal_UB_c("Thermal_UB_c");
    Thermal_UB_c = b_c * (theta_c.to(arcs_c) - theta_c.from(arcs_c));
    Thermal_UB_c -= S_max_c * DAL_multiplier;
    RSCED.add(Thermal_UB_c.in(arcs_c) <= 0);
    Constraint<> Thermal_LB_c("Thermal_LB_c");
    Thermal_LB_c = b_c * (theta_c.to(arcs_c) - theta_c.from(arcs_c));
    Thermal_LB_c += S_max_c * DAL_multiplier;
    RSCED.add(Thermal_LB_c.in(arcs_c) >= 0);
    
    /** Post-recourse constraints */
    /* REF BUS */
    Constraint<> Ref_Bus_recourse_c("Ref_Bus_recourse_c");
    Ref_Bus_recourse_c = theta_c.in(grid.ref_bus);
    RSCED.add(Ref_Bus_recourse_c == 0);
    
    //  /* CVaR lower bound */
    //  Constraint<> OBJ_BOUND_c("OBJ_BOUND_c");
    //  OBJ_BOUND_c = y_c.in(contingencies) + z - product(c2, dPg_p) - product(c2, dPg_n) - product(voll, dD_c) - product(c2, Pg);
    //  RSCED.add(OBJ_BOUND_c.in(contingencies) >= 0);
    
    /* Power flow constraint */
    Constraint<> Flow_P_recourse_c("Flow_P_recourse_c");
    Flow_P_recourse_c = Pf_recourse_c + b_c * (theta_recourse_c.from(arcs_c) - theta_recourse_c.to(arcs_c));
    RSCED.add(Flow_P_recourse_c.in(arcs_c) == 0);
    
    /* Power balance constraint */
    Constraint<> KCL_P_recourse_c("KCL_P_recourse_c");
    KCL_P_recourse_c  = Pf_recourse_c + pl + gs - sum(Pg, gen_nodes) - sum(dPg_p, gen_nodes) + sum(dPg_n, gen_nodes_c) - sum(dD_c, gen_nodes_c);
    RSCED.add(KCL_P_recourse_c.in(nodes_c) == 0);
    
    /* Phase Angle Bounds constraints */
    Constraint<> PAD_UB_recourse_c("PAD_UB_recourse_c");
    PAD_UB_recourse_c = theta_recourse_c.from(nodes_c) - theta_recourse_c.to(nodes_c);
    PAD_UB_recourse_c -= th_max;
    RSCED.add(PAD_UB_recourse_c.in(nodes_c) <= 0);
    Constraint<> PAD_LB_recourse_c("PAD_LB_recourse_c");
    PAD_LB_recourse_c = theta_recourse_c.from(nodes_c) - theta_recourse_c.to(nodes_c);
    PAD_LB_recourse_c -= th_min;
    RSCED.add(PAD_LB_recourse_c.in(nodes_c) >= 0);
    
    /* Line Limits constraints */
    Constraint<> Thermal_UB_recourse_c("Thermal_UB_recourse_c");
    Thermal_UB_recourse_c = b_c * (theta_recourse_c.to(arcs_c) - theta_recourse_c.from(arcs_c));
    Thermal_UB_recourse_c -= S_max_c * STE_multiplier;
    RSCED.add(Thermal_UB_recourse_c.in(arcs_c) <= 0);
    Constraint<> Thermal_LB_recourse_c("Thermal_LB_recourse_c");
    Thermal_LB_recourse_c = b_c * (theta_recourse_c.to(arcs_c) - theta_recourse_c.from(arcs_c));
    Thermal_LB_recourse_c += S_max_c * STE_multiplier;
    RSCED.add(Thermal_LB_recourse_c.in(arcs_c) >= 0);
    
    /* Generation capacity limits */
    Constraint<> Gen_UB_c("Gen_recourse_UB_c");
    Gen_UB_c = sum(Pg, gens) + sum(dPg_p, gens_c) - sum(dPg_n, gens_c);
    Gen_UB_c -= pg_max;
    RSCED.add(Gen_UB_c.in(gens_c) <= 0);
    Constraint<> Gen_LB_c("Gen_recourse_LB_c");
    Gen_LB_c = sum(Pg, gens) + sum(dPg_p, gens_c) - sum(dPg_n, gens_c);
    Gen_LB_c += pg_min;
    RSCED.add(Gen_LB_c.in(gens_c) >= 0);
    
    /* Recourse ramp capacity limits */
    Constraint<> Recourse_UB_c("Recourse_UB_c");
    Recourse_UB_c = sum(dPg_p, gens_c) - sum(dPg_n, gens_c);
    Recourse_UB_c -= ramp_max;
    RSCED.add(Recourse_UB_c.in(gens_c) <= 0);
    Constraint<> Recourse_LB_c("Recourse_LB_c");
    Recourse_LB_c = sum(dPg_p, gens_c) - sum(dPg_n, gens_c);
    Recourse_LB_c += ramp_max;
    RSCED.add(Recourse_LB_c.in(gens_c) >= 0);
    
    /** Objective */
    obj += alpha_prime * failure_probability * sum(y_c, contingencies);
    
    /* Set objective */
    RSCED.min(obj);
    
    RSCED.print();
    
    /** Solve */
    solver<> RSCED_SOLVER(RSCED, ipopt);
    
    if (use_gurobi)
        RSCED_SOLVER = solver<>(RSCED, gurobi);
    else if (use_cplex)
        RSCED_SOLVER = solver<>(RSCED, cplex);
    else
        RSCED_SOLVER = solver<>(RSCED, ipopt);
    
    auto solver_time_start = get_wall_time();
    RSCED_SOLVER.run(output = 0, tol = 1e-6);
    solver_time_end = get_wall_time();
    total_time_end = get_wall_time();
    solve_time = solver_time_end - solver_time_start;
    total_time = total_time_end - total_time_start;
    
    /** Uncomment next line to print expanded model */
    /* DCOPF.print(); */
    string out = "DATA_OPF, " + grid._name + ", " + to_string(nb_buses) + ", " + to_string(nb_lines) + ", "
    + to_string(RSCED.get_obj_val()) + ", " + to_string(-numeric_limits<double>::infinity()) + ", "
    + to_string(solve_time) + ", GlobalOptimal, " + to_string(total_time);
    DebugOn(out << endl);
    return 0;
}

////
////  DCOPF.cpp
////  Gravity
////
////  Created by Hassan Hijazi on 19 Jan 18.
////
////
//#include <stdio.h>
//#include <PowerNet.h>
//#include <gravity/solver.h>
//#ifdef USE_OPT_PARSER
//#include <optionParser.hpp>
//#endif
//using namespace std;
//using namespace gravity;
//
//int main (int argc, char * argv[])
//{
//    int output = 0;
//    bool projected = false, use_cplex = false, use_gurobi = false;
//    double tol = 1e-6;
//    double solver_time_end, total_time_end, solve_time, total_time;
//    string mehrotra = "no", log_level="0";
//    string fname = string(prj_dir)+"/data_sets/Power/nesta_case5_pjm.m";
//
//    string path = argv[0];
//    string solver_str="ipopt";
//    string proj_str="0";
//
//#ifdef USE_OPT_PARSER
//    /** Create a OptionParser with options */
//    auto options = readOptions(argc, argv);
//    options.add_option("f", "file", "Input file name", fname);
//    options.add_option("p", "project", "Project the power flow variables", proj_str);
//    options.add_option("s", "solver", "Solvers: ipopt/cplex/gurobi, default = ipopt", solver_str);
//
//    /** Parse the options and verify that all went well. If not, errors and help will be shown */
//    bool correct_parsing = options.parse_options(argc, argv);
//
//    if(!correct_parsing){
//        return EXIT_FAILURE;
//    }
//
//    output = op::str2int(options["l"]);
//
//    fname = options["f"];
//    bool has_help = op::str2bool(options["h"]);
//    if(has_help) {
//        options.show_help();
//        exit(0);
//    }
//    solver_str = options["s"];
//    if (solver_str.compare("gurobi")==0) {
//        use_gurobi = true;
//    }
//    else if(solver_str.compare("cplex")==0) {
//        use_cplex = true;
//    }
//    solver_str = options["p"];
//    if (solver_str.compare("1")==0) {
//        projected = true;
//    }
//#else
//    if(argc==2){
//        fname=argv[1];
//    }
//    else{
//        fname=string(prj_dir)+"/data_sets/Power/nesta_case5_pjm.m";
//    }
//#endif
//    double total_time_start = get_wall_time();
//    PowerNet grid;
//    grid.readgrid(fname,false);
//
//    /* Grid Stats */
//    auto node_pairs = grid.get_node_pairs();
//    auto nb_gen = grid.get_nb_active_gens();
//    auto nb_lines = grid.get_nb_active_arcs();
//    auto nb_buses = grid.get_nb_active_nodes();
//    DebugOn("nb active gens = " << nb_gen << endl);
//    DebugOn("nb active lines = " << nb_lines << endl);
//    DebugOn("nb active buses = " << nb_buses << endl);
//
//    /** Sets */
//    auto nodes = indices(grid.nodes);
//    auto arcs = indices(grid.arcs);
//    auto gens = indices(grid.gens);
//    auto gen_nodes = grid.gens_per_node();
//    auto out_arcs = grid.out_arcs_per_node();
//    auto in_arcs = grid.in_arcs_per_node();
//
//    /* Grid Parameters */
//    auto pg_min = grid.pg_min.in(gens);
//    auto pg_max = grid.pg_max.in(gens);
//    auto pl = grid.pl.in(nodes);
//    auto c1 = grid.c1.in(gens);
//    auto c2 = grid.c2.in(gens);
//    auto c0 = grid.c0.in(gens);
//    auto gs = grid.gs.in(nodes);
//    auto S_max = grid.S_max.in(arcs);
//    auto b = grid.b.in(arcs);
//    auto th_min = grid.th_min.in(node_pairs);
//    auto th_max = grid.th_max.in(node_pairs);
//
//    /** Declare model */
//    Model<> DCOPF("DCOPF Model");
//
//    if(projected) {/* Project out the power flow variables*/
//        /** Variables */
//        /* Power generation variables */
//        var<> Pg("Pg", pg_min, pg_max);
//        DCOPF.add(Pg.in(gens));
//
//
//        /* Phase angle variables */
//        var<> theta("𝛉");
//        DCOPF.add(theta.in(nodes));
//
//        /**  Objective */
//        auto obj = product(c1,Pg) + product(c2,pow(Pg,2)) + sum(c0);
//        DCOPF.min(obj);
//
//        /** Constraints */
//
//        /* REF BUS */
//        Constraint<> Ref_Bus("Ref_Bus");
//        Ref_Bus = theta(grid.ref_bus);
//        DCOPF.add(Ref_Bus == 0);
//
//        /* Flow conservation */
//        Constraint<> KCL_P("KCL_P");
//        KCL_P  = b.tr().in(in_arcs)*(theta.from(in_arcs)-theta.to(in_arcs)) - b.tr().in(out_arcs)*(theta.from(out_arcs)-theta.to(out_arcs)) + pl + gs - sum(Pg, gen_nodes);
//        DCOPF.add(KCL_P.in(nodes) == 0);
//
//        /* Phase Angle Bounds constraints */
//        Constraint<> PAD_UB("PAD_UB");
//        PAD_UB = theta.from(node_pairs) - theta.to(node_pairs);
//        PAD_UB -= th_max;
//        DCOPF.add(PAD_UB.in(node_pairs) <= 0);
//        Constraint<> PAD_LB("PAD_LB");
//        PAD_LB = theta.from(node_pairs) - theta.to(node_pairs);
//        PAD_LB -= th_min;
//        DCOPF.add(PAD_LB.in(node_pairs) >= 0);
//
//        /* Line Limits constraints */
//        Constraint<> Thermal_UB("Thermal_UB");
//        Thermal_UB = b*(theta.to(arcs) - theta.from(arcs));
//        Thermal_UB -= S_max;
//        DCOPF.add(Thermal_UB.in(arcs) <= 0);
//        Constraint<> Thermal_LB("Thermal_LB");
//        Thermal_LB = b*(theta.to(arcs) - theta.from(arcs));
//        Thermal_LB += S_max;
//        DCOPF.add(Thermal_LB.in(arcs) >= 0);
//    }
//    else {
//        /** Variables */
//        /* Power generation variables */
//        var<> Pg("Pg", pg_min, pg_max);
//        DCOPF.add(Pg.in(gens));
//
//        /* Power flow variables */
//        var<> Pf("Pf", -1*S_max, S_max);
//        DCOPF.add(Pf.in(arcs));
//
//        /* Phase angle variables */
//        var<> theta("𝛉");
//        DCOPF.add(theta.in(nodes));
//
//        /**  Objective */
//        auto obj = product(c1,Pg) + product(c2,pow(Pg,2)) + sum(c0);
//        DCOPF.min(obj);
//
//        /** Constraints */
//
//        /* REF BUS */
//        Constraint<> Ref_Bus("Ref_Bus");
//        Ref_Bus = theta(grid.ref_bus);
//        DCOPF.add(Ref_Bus == 0);
//
//        /* Flow conservation */
//        Constraint<> KCL_P("KCL_P");
//        KCL_P  = sum(Pf, out_arcs) - sum(Pf, in_arcs) + pl + gs - sum(Pg, gen_nodes);
//        DCOPF.add(KCL_P.in(nodes) == 0);
//
//        /* AC Power Flow */
//        Constraint<> Flow_P("Flow_P");
//        Flow_P = Pf + b*(theta.from(arcs) - theta.to(arcs));
//        DCOPF.add(Flow_P.in(arcs) == 0);
//
//        /* Phase Angle Bounds constraints */
//        Constraint<> PAD_UB("PAD_UB");
//        PAD_UB = theta.from(node_pairs) - theta.to(node_pairs);
//        PAD_UB -= th_max;
//        DCOPF.add(PAD_UB.in(node_pairs) <= 0);
//        Constraint<> PAD_LB("PAD_LB");
//        PAD_LB = theta.from(node_pairs) - theta.to(node_pairs);
//        PAD_LB -= th_min;
//        DCOPF.add(PAD_LB.in(node_pairs) >= 0);
//    }
//    /* Solver selection */
//    if (use_cplex) {
//        solver<> DCOPF_CPX(DCOPF, cplex);
//        auto solver_time_start = get_wall_time();
//        DCOPF_CPX.run(output = 0, tol = 1e-6);
//        solver_time_end = get_wall_time();
//        total_time_end = get_wall_time();
//        solve_time = solver_time_end - solver_time_start;
//        total_time = total_time_end - total_time_start;
//    }
//    //    else if (use_gurobi) {
//    //        solver DCOPF_GRB(DCOPF, gurobi);
//    //        auto solver_time_start = get_wall_time();
//    //        DCOPF_GRB.run(output, relax = false, tol = 1e-6);
//    //        solver_time_end = get_wall_time();
//    //        total_time_end = get_wall_time();
//    //        solve_time = solver_time_end - solver_time_start;
//    //        total_time = total_time_end - total_time_start;
//    //    }
//    else {
//        solver<> DCOPF_IPT(DCOPF,ipopt);
//        auto solver_time_start = get_wall_time();
//        DCOPF_IPT.run(output = 5, tol = 1e-6);
//        /** Warm starting Cplex
//         solver<> DCOPF_CPX(DCOPF, cplex);
//         DCOPF_CPX.run(output = 5, tol = 1e-6);
//         **/
//        solver_time_end = get_wall_time();
//        total_time_end = get_wall_time();
//        solve_time = solver_time_end - solver_time_start;
//        total_time = total_time_end - total_time_start;
//    }
//    /** Uncomment next line to print expanded model */
//    /* DCOPF.print(); */
//    string out = "DATA_OPF, " + grid._name + ", " + to_string(nb_buses) + ", " + to_string(nb_lines) +", " + to_string(DCOPF.get_obj_val()) + ", " + to_string(-numeric_limits<double>::infinity()) + ", " + to_string(solve_time) + ", GlobalOptimal, " + to_string(total_time);
//    DebugOn(out <<endl);
//    return 0;
//}
