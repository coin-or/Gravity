//
//  ACOPF.cpp
//  Gravity
//
//  Created by Guanglei Wang on 6/9/17.
//
//
#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include "../PowerNet.h"
#include <gravity/solver.h>
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
//            fname = "../../data_sets/Power/nesta_case5_pjm.m";
           //fname = "../../data_sets/Power/nesta_case14_ieee.m";
           //fname = "../../data_sets/Power/nesta_case9241_pegase.m";
           //fname = "../../data_sets/Power/nesta_case2383wp_mp.m";
           //fname = "../../data_sets/Power/nesta_case1354_pegase_api.m";
           //fname = "../../data_sets/Power/nesta_case118_ieee.m";
//           fname = "/Users/hh/Dropbox/Work/Dev/pglib-opf/pglib_opf_case2383wp_k.m";
        fname = "/Users/hh/Dropbox/Work/Dev/pglib-opf/pglib_opf_case5_pjm.m";
//        fname = "/Users/hh/Dropbox/Work/Dev/pglib-opf/pglib_opf_case6495_rte.m";
//        fname = "/Users/hh/Dropbox/Work/Dev/nesta-0.7.0/opf/nesta_case3_lmbd.m";
//        fname = "/Users/hh/Dropbox/Work/Dev/nesta-0.7.0/opf/nesta_case5_pjm.m";
    }
    // ACOPF
    PowerNet grid;
//    fname = "../../data_sets/Power/nesta_case3_lmbd.m";
//    fname = "../../data_sets/Power/nesta_case14_ieee.m";
//    fname = "../../data_sets/Power/nesta_case9241_pegase.m";
//    fname = "/Users/hh/Dropbox/Work/Dev/nesta-0.7.0/opf/nesta_case3375wp_mp.m";
//    fname = "../../data_sets/Power/nesta_case300_ieee.m";
//    fname = "../../data_sets/Power/nesta_case2383wp_mp.m";
    grid.readgrid(fname);

    // Grid Parameters
    unsigned nb_gen = grid.get_nb_active_gens();
    unsigned nb_lines = grid.get_nb_active_arcs();
    unsigned nb_buses = grid.get_nb_active_nodes();


    DebugOn("nb gens = " << nb_gen << endl);
    DebugOn("nb lines = " << 2*nb_lines << endl);
    DebugOn("nb buses = " << nb_buses << endl);

    PowerModelType pmt = ACPOL;
    /** build model */
    if (argc >= 3) {
        if(!strcmp(argv[2],"ACPOL")) pmt = ACPOL;
        else if(!strcmp(argv[2],"ACRECT")) pmt = ACRECT;
        else if(!strcmp(argv[2],"QC")) pmt = QC;
        else if(!strcmp(argv[2],"QC_SDP")) pmt = QC_SDP;
        else if(!strcmp(argv[2],"OTS")) pmt = OTS;
        else if(!strcmp(argv[2],"SOCP")) pmt = SOCP;
        else if(!strcmp(argv[2],"SDP")) pmt = SDP;
        else if(!strcmp(argv[2],"DC")) pmt = DC;
        else if(!strcmp(argv[2],"QC_OTS_O")) pmt = QC_OTS_O;
        else if(!strcmp(argv[2],"QC_OTS_N")) pmt = QC_OTS_N;
        else if(!strcmp(argv[2],"QC_OTS_L")) pmt = QC_OTS_L;
        else if(!strcmp(argv[2],"SOCP_OTS")) pmt = SOCP_OTS;
        else{
            throw invalid_argument("Unknown model type.\n");            
        }
    }
    bool polar = (pmt==ACPOL);
    if (polar) {
        DebugOn("Using polar model\n");
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

    ACOPF.add_var(Pf_from^(nb_lines));
    ACOPF.add_var(Qf_from^(nb_lines));
    ACOPF.add_var(Pf_to^(nb_lines));
    ACOPF.add_var(Qf_to^(nb_lines));

    // voltage related variables.
        var<Real> theta("theta");
        var<Real> v("|V|", grid.v_min.in(grid.nodes), grid.v_max.in(grid.nodes));
        var<Real> vr("vr", grid.v_max.in(grid.nodes));
        var<Real> vi("vi", grid.v_max.in(grid.nodes));
    
    if (polar) {
        ACOPF.add_var(v^(nb_buses));
        ACOPF.add_var(theta^(nb_buses));
        v.initialize_all(1);
    }
    else {
        ACOPF.add_var(vr^(nb_buses));
        ACOPF.add_var(vi^(nb_buses));
        vr.initialize_all(1);
    }

    /** Construct the objective function */
    func_ obj;
    for (auto g:grid.gens) {
        if (g->_active) {
            obj += grid.c1(g->_name)*Pg(g->_name) + grid.c2(g->_name)*Pg(g->_name)*Pg(g->_name) + grid.c0(g->_name);
        }
    }

    ACOPF.set_objective(min(obj));
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
        KCL_P = 0;
        KCL_Q = 0;
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
        Vol_limit_UB = power(vr.in(grid.nodes), 2) + power(vi.in(grid.nodes), 2);
        Vol_limit_UB -= power(grid.v_max.in(grid.nodes), 2);
        ACOPF.add_constraint(Vol_limit_UB <= 0);

        Constraint Vol_limit_LB("Vol_limit_LB");
        Vol_limit_LB = power(vr.in(grid.nodes), 2) + power(vi.in(grid.nodes), 2);
        Vol_limit_LB -= power(grid.v_min.in(grid.nodes),2);    
        ACOPF.add_constraint(Vol_limit_LB >= 0);
    }

    
    /* Phase Angle Bounds constraints */
    Constraint PAD_UB("PAD_UB");
    Constraint PAD_LB("PAD_LB:");
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
        PAD_LB -= grid.tan_th_min*(vr.from()*vr.to() + vi.to()*vi.from());
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

    //solver OPF(ACOPF,cplex);
    solver OPF(ACOPF,ipopt);
    OPF.run();
    return 0;
}
