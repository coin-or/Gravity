//
//  ACOPF.cpp
//  Gravity
//
//  Created by Hassan Hijazi on Dec 7 2017
//
//
#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include "../PowerNet.h"
#include <gravity/solver.h>
#include <stdlib.h>
#include <optionParser.hpp>


using namespace std;
using namespace gravity;


int main (int argc, char * argv[])
{
//    string fname = "../../data_sets/Power/nesta_case5_pjm.m", mtype = "ACRECT";
    string fname = "/Users/hlh/Dropbox/Work/Dev/pglib-opf/pglib_opf_case30_ieee.m", mtype = "ACRECT";
    int output = 0;
    bool relax = false;
    double tol = 1e-6;
    string mehrotra = "no";
    
    // create a OptionParser with options
    op::OptionParser opt;
    opt.add_option("h", "help", "shows option help"); // no default value means boolean options, which default value is false
    opt.add_option("f", "file", "Input file name", fname );
    opt.add_option("m", "model", "power flow model: ACPOL/ACRECT", mtype );
    
    // parse the options and verify that all went well. If not, errors and help will be shown
    bool correct_parsing = opt.parse_options(argc, argv);
    
    if(!correct_parsing){
        return EXIT_FAILURE;
    }
    
    fname = opt["f"];
    mtype = opt["m"];
    
    bool has_help = op::str2bool(opt["h"]);
    // show help
    if(has_help) {
        opt.show_help();
        exit(0);
    }
    
    // ACOPF
    double total_time_start = get_wall_time();
    PowerNet grid;
    grid.readgrid(fname.c_str());

    // Grid Parameters
    unsigned nb_gen = grid.get_nb_active_gens();
    unsigned nb_lines = grid.get_nb_active_arcs();
    unsigned nb_buses = grid.get_nb_active_nodes();


    DebugOn("nb gens = " << nb_gen << endl);
    DebugOn("nb lines = " << 2*nb_lines << endl);
    DebugOn("nb buses = " << nb_buses << endl);

    PowerModelType pmt = ACPOL;
    if(!strcmp(mtype.c_str(),"ACRECT")) pmt = ACRECT;
    bool polar = (pmt==ACPOL);
    if (polar) {
        DebugOn("Using polar model\n");
    }
    else {
        DebugOn("Using rectangular model\n");
    }
    Model ACOPF("AC-OPF Model");
    /** Variables */
    // power generation
    var<Real> Pg("Pg", grid.pg_min.in(grid.gens), grid.pg_max.in(grid.gens));
    var<Real> Qg ("Qg", grid.qg_min.in(grid.gens), grid.qg_max.in(grid.gens));
    ACOPF.add_var(Pg^(nb_gen));
    ACOPF.add_var(Qg^(nb_gen));

    // power flow
    var<Real> Pf_from("Pf_from", grid.S_max.in(grid.arcs));
    var<Real> Qf_from("Qf_from", grid.S_max.in(grid.arcs));
    var<Real> Pf_to("Pf_to", grid.S_max.in(grid.arcs));
    var<Real> Qf_to("Qf_to", grid.S_max.in(grid.arcs));

//    var<Real> Pf_to("Pf_to");
//    var<Real> Qf_to("Qf_to");
//    var<Real> Pf_from("Pf_from");
//    var<Real> Qf_from("Qf_from");


    ACOPF.add_var(Pf_from^(nb_lines));
    ACOPF.add_var(Qf_from^(nb_lines));
    ACOPF.add_var(Pf_to^(nb_lines));
    ACOPF.add_var(Qf_to^(nb_lines));

    // voltage related variables.
        var<Real> theta("theta");
        var<Real> v("|V|", grid.v_min.in(grid.nodes), grid.v_max.in(grid.nodes));
        var<Real> vr("vr");
        var<Real> vi("vi");
//        var<Real> vr("vr", grid.v_max.in(grid.nodes));
//        var<Real> vi("vi", grid.v_max.in(grid.nodes));
    
    if (polar) {
        ACOPF.add_var(v^(nb_buses));
        ACOPF.add_var(theta^(nb_buses));
        v.initialize_all(1);
    }
    else {
        ACOPF.add_var(vr^(nb_buses));
        ACOPF.add_var(vi^(nb_buses));
        vr.initialize_all(1.0);
    }

    /** Construct the objective function */
//    func_ obj = product(grid.c1, Pg) + product(grid.c2, power(Pg,2)) + sum(grid.c0);
//    ACOPF.min(obj.in(grid.gens));
    func_ obj;
    for (auto g:grid.gens) {
        if (g->_active) {
            obj += grid.c1(g->_name)*Pg(g->_name) + grid.c2(g->_name)*Pg(g->_name)*Pg(g->_name) + grid.c0(g->_name);
        }
    }
    ACOPF.min(obj);
    
    /** Define constraints */
    
    /* REF BUS */
    Constraint Ref_Bus("Ref_Bus");
    if (polar) {
        Ref_Bus = theta(grid.get_ref_bus());
    }
    else {
        Ref_Bus = vi(grid.get_ref_bus());
    }
    ACOPF.add_constraint(Ref_Bus = 0);
    
    
    //KCL
    for (auto b: grid.nodes) {
        if (!b->_active) {
            continue;
        }
        Bus* bus = (Bus*) b;
        Constraint KCL_P("KCL_P"+bus->_name);
        Constraint KCL_Q("KCL_Q"+bus->_name);

        /* Power Conservation */
//        KCL_P = Pf_from.in() + Pf_to.out() + grid.pl - Pg.in() + grid.gs*power(v,2);
        KCL_P  = sum(Pf_from.in(b->get_out())) + sum(Pf_to.in(b->get_in())) + bus->pl()- sum(Pg.in(bus->_gen));
        KCL_Q  = sum(Qf_from.in(b->get_out())) + sum(Qf_to.in(b->get_in())) + bus->ql()- sum(Qg.in(bus->_gen));
        
        /* Shunts */
        if (bus->gs()!=0) {
            if (polar) {
                KCL_P +=  bus->gs()*(power(v(bus->_name), 2));
            }
            else {
                KCL_P +=  bus->gs()*(power(vr(bus->_name), 2) + power(vi(bus->_name), 2));
            }
            DebugOff("Bus" << bus->_name << " : Shunt gs = " << bus->gs() << endl);
        }
        if (bus->bs()!=0) {
            if (polar) {
                KCL_Q -=  bus->bs()*(power(v(bus->_name), 2));
            }
            else {
                KCL_Q -=  bus->bs()*(power(vr(bus->_name), 2) + power(vi(bus->_name), 2));
            }
            DebugOff("Bus" << bus->_name << " : Shunt bs = " << bus->bs() << endl);
        }
        ACOPF.add_constraint(KCL_P = 0);
        ACOPF.add_constraint(KCL_Q = 0);
    }
    //AC Power Flow.

     /** TODO write the constraints in Complex form */
    Constraint Flow_P_From("Flow_P_From");
    Flow_P_From += Pf_from;
    if (polar) {
        Flow_P_From -= grid.g/power(grid.tr,2)*power(v.from(),2);
        Flow_P_From += grid.g/grid.tr*(v.from()*v.to()*cos(theta.from() - theta.to() - grid.as));
        Flow_P_From += grid.b/grid.tr*(v.from()*v.to()*sin(theta.from() - theta.to() - grid.as));
    }
    else {
        Flow_P_From -= grid.g_ff*(power(vr.from(), 2) + power(vi.from(), 2));
        Flow_P_From -= grid.g_ft*(vr.from()*vr.to() + vi.from()*vi.to());
        Flow_P_From -= grid.b_ft*(vi.from()*vr.to() - vr.from()*vi.to());
    }
    ACOPF.add_constraint(Flow_P_From.in(grid.arcs)=0);

    Constraint Flow_P_To("Flow_P_To");
    Flow_P_To += Pf_to;
    if (polar) {
        Flow_P_To -= grid.g*power(v.to(), 2);
        Flow_P_To += grid.g/grid.tr*(v.from()*v.to()*cos(theta.to() - theta.from() + grid.as));
        Flow_P_To += grid.b/grid.tr*(v.from()*v.to()*sin(theta.to() - theta.from() + grid.as));
    }
    else {
        Flow_P_To -= grid.g_tt*(power(vr.to(), 2) + power(vi.to(), 2));
        Flow_P_To -= grid.g_tf*(vr.from()*vr.to() + vi.from()*vi.to());
        Flow_P_To -= grid.b_tf*(vi.to()*vr.from() - vr.to()*vi.from());
    }
    ACOPF.add_constraint(Flow_P_To.in(grid.arcs)=0);

    Constraint Flow_Q_From("Flow_Q_From");
    Flow_Q_From += Qf_from;
    if (polar) {
        Flow_Q_From += (0.5*grid.ch+grid.b)/power(grid.tr,2)*power(v.from(),2);
        Flow_Q_From -= grid.b/grid.tr*(v.from()*v.to()*cos(theta.from() - theta.to() - grid.as));
        Flow_Q_From += grid.g/grid.tr*(v.from()*v.to()*sin(theta.from() - theta.to() - grid.as));
    }
    else {
        Flow_Q_From += grid.b_ff*(power(vr.from(), 2) + power(vi.from(), 2));
        Flow_Q_From += grid.b_ft*(vr.from()*vr.to() + vi.from()*vi.to());
        Flow_Q_From -= grid.g_ft*(vi.from()*vr.to() - vr.from()*vi.to());
    }
    ACOPF.add_constraint(Flow_Q_From.in(grid.arcs)=0);
    Constraint Flow_Q_To("Flow_Q_To");
    Flow_Q_To += Qf_to;
    if (polar) {
        Flow_Q_To += (0.5*grid.ch+grid.b)*power(v.to(),2);
        Flow_Q_To -= grid.b/grid.tr*(v.from()*v.to()*cos(theta.to() - theta.from() + grid.as));
        Flow_Q_To += grid.g/grid.tr*(v.from()*v.to()*sin(theta.to() - theta.from() + grid.as));
    }
    else {
        Flow_Q_To += grid.b_tt*(power(vr.to(), 2) + power(vi.to(), 2));
        Flow_Q_To += grid.b_tf*(vr.from()*vr.to() + vi.from()*vi.to());
        Flow_Q_To -= grid.g_tf*(vi.to()*vr.from() - vr.to()*vi.from());
    }
    ACOPF.add_constraint(Flow_Q_To.in(grid.arcs)=0);
    
    // AC voltage limit constraints.
    if (!polar) {
        Constraint Vol_limit_UB("Vol_limit_UB");
        Vol_limit_UB = power(vr, 2) + power(vi, 2);
        Vol_limit_UB -= power(grid.v_max, 2);
        ACOPF.add_constraint(Vol_limit_UB.in(grid.nodes) <= 0);

        Constraint Vol_limit_LB("Vol_limit_LB");
        Vol_limit_LB = power(vr, 2) + power(vi, 2);
        Vol_limit_LB -= power(grid.v_min,2);
        ACOPF.add_constraint(Vol_limit_LB.in(grid.nodes) >= 0);
        DebugOff(grid.v_min.to_str(true) << endl);
        DebugOff(grid.v_max.to_str(true) << endl);
    }

    
    /* Phase Angle Bounds constraints */
    Constraint PAD_UB("PAD_UB");
    Constraint PAD_LB("PAD_LB");
    auto bus_pairs = grid.get_bus_pairs();
    if (polar) {
        PAD_UB = theta.from() - theta.to();
        PAD_UB -= grid.th_max;
        PAD_LB = theta.from() - theta.to();
        PAD_LB -= grid.th_min;
        DebugOff(grid.th_min.to_str(true) << endl);
        DebugOff(grid.th_max.to_str(true) << endl);
    }
    else {        
        DebugOff("Number of bus_pairs = " << bus_pairs.size() << endl);
        PAD_UB = vi.from()*vr.to() - vr.from()*vi.to();
        PAD_UB -= grid.tan_th_max*(vr.from()*vr.to() + vi.from()*vi.to());
        
        PAD_LB = vi.from()*vr.to() - vr.from()*vi.to();
        PAD_LB -= grid.tan_th_min*(vr.from()*vr.to() + vi.from()*vi.to());
        DebugOff(grid.th_min.to_str(true) << endl);
        DebugOff(grid.th_max.to_str(true) << endl);
    }
    ACOPF.add_constraint(PAD_UB.in(bus_pairs) <= 0);
    ACOPF.add_constraint(PAD_LB.in(bus_pairs) >= 0);


//  Thermal Limit Constraints 
    Constraint Thermal_Limit_from("Thermal_Limit_from");
    Thermal_Limit_from += power(Pf_from, 2) + power(Qf_from, 2);
    Thermal_Limit_from -= power(grid.S_max, 2);
    ACOPF.add_constraint(Thermal_Limit_from.in(grid.arcs) <= 0);

    Constraint Thermal_Limit_to("Thermal_Limit_to");
    Thermal_Limit_to += power(Pf_to, 2) + power(Qf_to, 2);
    Thermal_Limit_to -= power(grid.S_max,2);
    ACOPF.add_constraint(Thermal_Limit_to.in(grid.arcs) <= 0);
    DebugOff(grid.S_max.in(grid.arcs).to_str(true) << endl);

    //solver OPF(ACOPF,cplex);
    solver OPF(ACOPF,ipopt);
    double solver_time_start = get_wall_time();
    OPF.run(output = 0, relax = false, tol = 1e-6, "ma27", mehrotra = "no");
    double solver_time_end = get_wall_time();
    double total_time_end = get_wall_time();
    auto solve_time = solver_time_end - solver_time_start;
    auto total_time = total_time_end - total_time_start;
//    ACOPF.print_expanded();
    string out = "DATA_OPF, " + grid._name + ", " + to_string(nb_buses) + ", " + to_string(nb_lines) +", " + to_string(ACOPF._obj_val) + ", " + to_string(-numeric_limits<double>::infinity()) + ", " + to_string(solve_time) + ", LocalOptimal, " + to_string(total_time);
    DebugOn(out <<endl);
    return 0;
}
