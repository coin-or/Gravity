//
//  ACOPF.cpp
//  Gravity
//
//  Created by Hassan Hijazi on Dec 7 2017
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
    string fname = string(prj_dir)+"/data_sets/Power/nesta_case5_pjm.m", mtype = "ACPOL";
    int output = 0;
    bool relax = false;
    double tol = 1e-6;
    string mehrotra = "no", log_level="0";
    
    /** create a OptionParser with options */
    op::OptionParser opt;
    opt.add_option("h", "help", "shows option help"); // no default value means boolean options, which default value is false
    opt.add_option("f", "file", "Input file name (def. ../data_sets/Power/nesta_case5_pjm.m)", fname );
    opt.add_option("l", "log", "Log level (def. 0)", log_level );
    opt.add_option("m", "model", "power flow model: ACPOL/ACRECT (def. ACPOL)", mtype );
    
    /** parse the options and verify that all went well. If not, errors and help will be shown */
    bool correct_parsing = opt.parse_options(argc, argv);
    
    if(!correct_parsing){
        return EXIT_FAILURE;
    }
    
    fname = opt["f"];
    mtype = opt["m"];
    output = op::str2int(opt["l"]);
    output = 5;
    bool has_help = op::str2bool(opt["h"]);
    /** show help */
    if(has_help) {
        opt.show_help();
        exit(0);
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
    

    
    PowerModelType pmt = ACRECT;
    if(!strcmp(mtype.c_str(),"ACRECT")) pmt = ACRECT;
    bool polar = (pmt==ACPOL);
    if (polar) {
        DebugOn("Using polar model\n");
    }
    else {
        DebugOn("Using rectangular model\n");
    }
    Model<> ACOPF("AC-OPF Model");
    /** Variables */
    /* Power generation variables */
    var<> Pg("Pg", pg_min, pg_max);
    var<> Qg ("Qg", qg_min, qg_max);
    ACOPF.add(Pg.in(gens),Qg.in(gens));
    
    /* Power flow variables */
    var<> Pf_from("Pf_from", -1*grid.S_max,grid.S_max);
    var<> Qf_from("Qf_from", -1*grid.S_max,grid.S_max);
    var<> Pf_to("Pf_to", -1*grid.S_max,grid.S_max);
    var<> Qf_to("Qf_to", -1*grid.S_max,grid.S_max);
    
    ACOPF.add(Pf_from.in(arcs), Qf_from.in(arcs),Pf_to.in(arcs),Qf_to.in(arcs));
    
    
    /** Voltage related variables */
    var<> theta("theta");
    var<> v("|V|", grid.v_min, grid.v_max);
    var<> vr("vr", -1*grid.v_max,grid.v_max);
    var<> vi("vi", -1*grid.v_max,grid.v_max);
    
    if (polar) {
        ACOPF.add(v.in(nodes));
        ACOPF.add(theta.in(nodes));
        v.initialize_all(1);
    }
    else {
        ACOPF.add(vr.in(nodes));
        ACOPF.add(vi.in(nodes));
        vr.initialize_all(1.0);
    }
    
    
    /** Construct the objective function */
    /**  Objective */
    auto obj = c1.tr()*Pg + c2.tr()*pow(Pg.vec(),2) + sum(c0);
    ACOPF.min(obj);

    /** Define constraints */
    
    /* REF BUS */
    Constraint<> Ref_Bus("Ref_Bus");
    if (polar) {
        Ref_Bus = theta(grid.ref_bus);
    }
    else {
        Ref_Bus = vi(grid.ref_bus);
    }
    ACOPF.add(Ref_Bus == 0);
    
    /** KCL Flow conservation */
    Constraint<> KCL_P("KCL_P");
    Constraint<> KCL_Q("KCL_Q");
    KCL_P  = sum(Pf_from, out_arcs) + sum(Pf_to, in_arcs) + grid.pl - sum(Pg, gen_nodes);
    KCL_Q  = sum(Qf_from, out_arcs) + sum(Qf_to, in_arcs) + grid.ql - sum(Qg, gen_nodes);
    /* Shunts */
    if (polar) {
        KCL_P +=  grid.gs*pow(v,2);
        KCL_Q -=  grid.bs*pow(v,2);
    }
    else {
        KCL_P +=  grid.gs*(pow(vr,2)+pow(vi,2));
        KCL_Q -=  grid.bs*(pow(vr,2)+pow(vi,2));
    }
    ACOPF.add(KCL_P.in(nodes) == 0);
    ACOPF.add(KCL_Q.in(nodes) == 0);
    
    /** AC Power Flows */
    /** TODO write the constraints in Complex form */
    Constraint<> Flow_P_From("Flow_P_From");
    Flow_P_From += Pf_from;
    if (polar) {
        Flow_P_From -= grid.g/pow(grid.tr,2)*pow(v.from(),2);
        Flow_P_From += grid.g/grid.tr*(v.from()*v.to()*cos(theta.from() - theta.to() - grid.as));
        Flow_P_From += grid.b/grid.tr*(v.from()*v.to()*sin(theta.from() - theta.to() - grid.as));
    }
    else {
        Flow_P_From -= grid.g_ff*(pow(vr.from(), 2) + pow(vi.from(), 2));
        Flow_P_From -= grid.g_ft*(vr.from()*vr.to() + vi.from()*vi.to());
        Flow_P_From -= grid.b_ft*(vi.from()*vr.to() - vr.from()*vi.to());
    }
    ACOPF.add(Flow_P_From.in(arcs)==0);
    
    Constraint<> Flow_P_To("Flow_P_To");
    Flow_P_To += Pf_to;
    if (polar) {
        Flow_P_To -= grid.g*pow(v.to(), 2);
        Flow_P_To += grid.g/grid.tr*(v.from()*v.to()*cos(theta.to() - theta.from() + grid.as));
        Flow_P_To += grid.b/grid.tr*(v.from()*v.to()*sin(theta.to() - theta.from() + grid.as));
    }
    else {
        Flow_P_To -= grid.g_tt*(pow(vr.to(), 2) + pow(vi.to(), 2));
        Flow_P_To -= grid.g_tf*(vr.from()*vr.to() + vi.from()*vi.to());
        Flow_P_To -= grid.b_tf*(vi.to()*vr.from() - vr.to()*vi.from());
    }
    ACOPF.add(Flow_P_To.in(arcs)==0);
    
    Constraint<> Flow_Q_From("Flow_Q_From");
    Flow_Q_From += Qf_from;
    if (polar) {
        Flow_Q_From += (0.5*grid.ch+grid.b)/pow(grid.tr,2)*pow(v.from(),2);
        Flow_Q_From -= grid.b/grid.tr*(v.from()*v.to()*cos(theta.from() - theta.to() - grid.as));
        Flow_Q_From += grid.g/grid.tr*(v.from()*v.to()*sin(theta.from() - theta.to() - grid.as));
    }
    else {
        Flow_Q_From += grid.b_ff*(pow(vr.from(), 2) + pow(vi.from(), 2));
        Flow_Q_From += grid.b_ft*(vr.from()*vr.to() + vi.from()*vi.to());
        Flow_Q_From -= grid.g_ft*(vi.from()*vr.to() - vr.from()*vi.to());
    }
    ACOPF.add(Flow_Q_From.in(arcs)==0);
    
    Constraint<> Flow_Q_To("Flow_Q_To");
    Flow_Q_To += Qf_to;
    if (polar) {
        Flow_Q_To += (0.5*grid.ch+grid.b)*pow(v.to(),2);
        Flow_Q_To -= grid.b/grid.tr*(v.from()*v.to()*cos(theta.to() - theta.from() + grid.as));
        Flow_Q_To += grid.g/grid.tr*(v.from()*v.to()*sin(theta.to() - theta.from() + grid.as));
    }
    else {
        Flow_Q_To += grid.b_tt*(pow(vr.to(), 2) + pow(vi.to(), 2));
        Flow_Q_To += grid.b_tf*(vr.from()*vr.to() + vi.from()*vi.to());
        Flow_Q_To -= grid.g_tf*(vi.to()*vr.from() - vr.to()*vi.from());
    }
    ACOPF.add(Flow_Q_To.in(arcs)==0);
    
    /** AC voltage limit constraints. */
    if (!polar) {
        Constraint<> Vol_limit_UB("Vol_limit_UB");
        Vol_limit_UB = pow(vr, 2) + pow(vi, 2);
        Vol_limit_UB -= pow(grid.v_max, 2);
        ACOPF.add(Vol_limit_UB.in(nodes) <= 0);
        
        Constraint<> Vol_limit_LB("Vol_limit_LB");
        Vol_limit_LB = pow(vr, 2) + pow(vi, 2);
        Vol_limit_LB -= pow(grid.v_min,2);
        ACOPF.add(Vol_limit_LB.in(nodes) >= 0);
    }
    
    
    /* Phase Angle Bounds constraints */
    Constraint<> PAD_UB("PAD_UB");
    Constraint<> PAD_LB("PAD_LB");
    if (polar) {
        PAD_UB = theta.from() - theta.to();
        PAD_UB -= grid.th_max;
        PAD_LB = theta.from() - theta.to();
        PAD_LB -= grid.th_min;
    }
    else {
        DebugOff("Number of bus_pairs = " << bus_pairs.size() << endl);
        PAD_UB = vi.from()*vr.to() - vr.from()*vi.to();
        PAD_UB -= grid.tan_th_max*(vr.from()*vr.to() + vi.from()*vi.to());
        
        PAD_LB = vi.from()*vr.to() - vr.from()*vi.to();
        PAD_LB -= grid.tan_th_min*(vr.from()*vr.to() + vi.from()*vi.to());
    }
        ACOPF.add(PAD_UB.in(bus_pairs) <= 0);
        ACOPF.add(PAD_LB.in(bus_pairs) >= 0);
    
    
    /*  Thermal Limit Constraints */
    Constraint<> Thermal_Limit_from("Thermal_Limit_from");
    Thermal_Limit_from += pow(Pf_from, 2) + pow(Qf_from, 2);
    Thermal_Limit_from -= pow(grid.S_max, 2);
    ACOPF.add(Thermal_Limit_from.in(arcs) <= 0);
    
    Constraint<> Thermal_Limit_to("Thermal_Limit_to");
    Thermal_Limit_to += pow(Pf_to, 2) + pow(Qf_to, 2);
    Thermal_Limit_to -= pow(grid.S_max,2);
    ACOPF.add(Thermal_Limit_to.in(arcs) <= 0);
    solver<> OPF(ACOPF,ipopt);
    double solver_time_start = get_wall_time();
    OPF.run(output=5, tol = 1e-6);
    double solver_time_end = get_wall_time();
    double total_time_end = get_wall_time();
    auto solve_time = solver_time_end - solver_time_start;
    auto total_time = total_time_end - total_time_start;
    /* Uncomment lines below to print model */
    /*
    ACOPF.print_symbolic();
    ACOPF.print();
     */
    /** Terminal output */
    string out = "DATA_OPF, " + grid._name + ", " + to_string(nb_buses) + ", " + to_string(nb_lines) +", " + to_string_with_precision(ACOPF.get_obj_val(),10) + ", " + to_string(-numeric_limits<double>::infinity()) + ", " + to_string(solve_time) + ", LocalOptimal, " + to_string(total_time);
    DebugOn(out <<endl);
    return 0;
}
