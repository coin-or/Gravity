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


//int main (int argc, char * argv[])
//{
//    string fname = "../data_sets/Power/nesta_case5_pjm.m", mtype = "ACPOL";
//    DebugOn("argv[0] =" << argv[0] << endl);
//    string path = argv[0];
//    int output = 0;
//    bool relax = false;
//    double tol = 1e-6;
//    string mehrotra = "no", log_level="0";
//    
//    /** create a OptionParser with options */
//    op::OptionParser opt;
//    opt.add_option("h", "help", "shows option help"); // no default value means boolean options, which default value is false
//    opt.add_option("f", "file", "Input file name (def. ../data_sets/Power/nesta_case5_pjm.m)", fname );
//    opt.add_option("l", "log", "Log level (def. 0)", log_level );
//    opt.add_option("m", "model", "power flow model: ACPOL/ACRECT (def. ACPOL)", mtype );
//    
//    /** parse the options and verify that all went well. If not, errors and help will be shown */
//    bool correct_parsing = opt.parse_options(argc, argv);
//    
//    if(!correct_parsing){
//        return EXIT_FAILURE;
//    }
//    
//    fname = opt["f"];
//    mtype = opt["m"];
//    output = op::str2int(opt["l"]);
//    bool has_help = op::str2bool(opt["h"]);
//    /** show help */
//    if(has_help) {
//        opt.show_help();
//        exit(0);
//    }
//    
//    double total_time_start = get_wall_time();
//    PowerNet grid;
//    grid.readgrid(fname.c_str());
//
//    /* Grid Parameters */
//    unsigned nb_gen = grid.get_nb_active_gens();
//    unsigned nb_lines = grid.get_nb_active_arcs();
//    unsigned nb_buses = grid.get_nb_active_nodes();
//
//
//    DebugOn("nb gens = " << nb_gen << endl);
//    DebugOn("nb lines = " << nb_lines << endl);
//    DebugOn("nb buses = " << nb_buses << endl);
//
//    PowerModelType pmt = ACPOL;
//    if(!strcmp(mtype.c_str(),"ACRECT")) pmt = ACRECT;
//    bool polar = (pmt==ACPOL);
//    if (polar) {
//        DebugOn("Using polar model\n");
//    }
//    else {
//        DebugOn("Using rectangular model\n");
//    }
//    Model ACOPF("AC-OPF Model");
//    /** Variables */
//    /* Power generation variables */
//    var<double> Pg("Pg", grid.pg_min, grid.pg_max);
//    var<double> Qg("Qg", grid.qg_min, grid.qg_max);
//
//
//    /** Voltage related variables */
//    var<double> vr("vr", grid.v_max);
//    var<double> vi("vi", grid.v_max);
//    
//    ACOPF.add(vr.in(grid.nodes));
//    ACOPF.add(vi.in(grid.nodes));
//    vr.initialize_all(1.0);
//
//    /** Construct the objective function */
//    indices gin;
//    for (auto n:grid.nodes) {
//        auto b = (Bus*)n;
//        if (b->_gen.size()>1) {
//            auto nb_gens = b->_gen.size();
//            for(unsigned i = 1; i< nb_gens; i++){
//                gin.add(b->_gen[i]->_name);
//            }
//        }
//    }
//    
//    ACOPF.add(Pg.in(gin));
//    ACOPF.add(Qg.in(gin));
//    
//    auto obj = product(grid.c1.in(gin), Pg.in(gin)) + product(grid.c2.in(gin), power(Pg.in(gin),2)) + sum(grid.c0.in(gin));
//    
//    for (auto n:grid.nodes) {
//        auto b = (Bus*)n;
//        
//            //If Gmin = G_max, transform into equation
//        for (auto a:b->get_out()) {
//            b->p_injec += grid.g_ff(a->_name)*(power(vr(a->_src->_name), 2) + power(vi(a->_src->_name), 2)) + grid.g_ft(a->_name)*(vr(a->_src->_name)*vr(a->_dest->_name) + vi(a->_src->_name)*vi(a->_dest->_name)) + grid.b_ft(a->_name)*(vi(a->_src->_name)*vr(a->_dest->_name) - vr(a->_src->_name)*vi(a->_dest->_name));//p_from
//            b->q_injec += grid.g_ft*(vi(a->_src->_name)*vr(a->_dest->_name) - vr(a->_src->_name)*vi(a->_dest->_name)) - grid.b_ff*(power(vr(a->_src->_name), 2) + power(vi(a->_src->_name), 2))- grid.b_ft*(vr(a->_src->_name)*vr(a->_dest->_name) + vi(a->_src->_name)*vi(a->_dest->_name));//q_from
//        }
//        for (auto a:b->get_in()) {
//            b->p_injec += grid.g_tt(a->_name)*(power(vr(a->_dest->_name), 2) + power(vi(a->_dest->_name), 2)) + grid.g_tf(a->_name)*(vr(a->_src->_name)*vr(a->_dest->_name) + vi(a->_src->_name)*vi(a->_dest->_name)) + grid.b_tf(a->_name)*(vi(a->_src->_name)*vr(a->_dest->_name) - vr(a->_src->_name)*vi(a->_dest->_name));//p_to
//            b->q_injec += grid.g_ft*(vi(a->_src->_name)*vr(a->_dest->_name) - vr(a->_src->_name)*vi(a->_dest->_name)) - grid.b_ff*(power(vr(a->_src->_name), 2) + power(vi(a->_src->_name), 2))- grid.b_ft*(vr(a->_src->_name)*vr(a->_dest->_name) + vi(a->_src->_name)*vi(a->_dest->_name));//q_to
//        }
//        b->p_injec += b->pl() + b->gs()*(power(vr(b->_name),2)+power(vi(b->_name),2));
//        b->q_injec += b->ql() - b->bs()*(power(vr(b->_name),2)+power(vi(b->_name),2));
//        if (b->_gen.size()==0) {
//            Constraint KCL_P("KCL_P"+b->_name);
//            Constraint KCL_Q("KCL_Q"+b->_name);
//            KCL_P  = b->p_injec;
//            KCL_Q  = b->q_injec;
//            ACOPF.add(KCL_P == 0);
//            ACOPF.add(KCL_Q == 0);
//        }
//        
//        if (b->_gen.size()==1) {
//            obj += b->_gen.back()->_cost->c2*power(b->p_injec,2);
//            obj += b->_gen.back()->_cost->c1*(b->p_injec);
//        }
//    }
//    
//    
//    //    obj.print_expanded();
//    
//    Constraint Thermal_Limit_from("Thermal_Limit_from");
//    Thermal_Limit_from += power(grid.g_ff*(power(vr.from(), 2) + power(vi.from(), 2)) + grid.g_ft*(vr.from()*vr.to() + vi.from()*vi.to()) + grid.b_ft*(vi.from()*vr.to() - vr.from()*vi.to()), 2);// p_from^2
//    Thermal_Limit_from += power(grid.g_ft*(vi.from()*vr.to() - vr.from()*vi.to()) - grid.b_ff*(power(vr.from(), 2) + power(vi.from(), 2))- grid.b_ft*(vr.from()*vr.to() + vi.from()*vi.to()), 2);//q_from^2
//    Thermal_Limit_from -= power(grid.S_max, 2);
////    ACOPF.add(Thermal_Limit_from.in(grid.arcs) <= 0);
////    Thermal_Limit_from.print_expanded();
//    
//    Constraint Thermal_Limit_to("Thermal_Limit_to");
//    Thermal_Limit_to += power(grid.g_tt*(power(vr.to(), 2) + power(vi.to(), 2)) + grid.g_tf*(vr.from()*vr.to() + vi.from()*vi.to()) + grid.b_tf*(vi.to()*vr.from() - vr.to()*vi.from()), 2);// p_to^2
//    Thermal_Limit_to += power(grid.g_tf*(vi.to()*vr.from() - vr.to()*vi.from()) - grid.b_tt*(power(vr.to(), 2) + power(vi.to(), 2)) - grid.b_tf*(vr.from()*vr.to() + vi.from()*vi.to()), 2);// q_to^2
//    Thermal_Limit_to -= power(grid.S_max,2);
////    ACOPF.add(Thermal_Limit_to.in(grid.arcs) <= 0);
////    Thermal_Limit_to.print_expanded();
//    
//    
//    
//    //    /* REF BUS */
//    Constraint Ref_Bus("Ref_Bus");
//    Ref_Bus = vi(grid.get_ref_bus());
//    ACOPF.add(Ref_Bus == 0);
//    
//    
//    Constraint Vol_limit_UB("Vol_limit_UB");
//    Vol_limit_UB = power(vr, 2) + power(vi, 2);
//    Vol_limit_UB -= power(grid.v_max, 2);
//    ACOPF.add(Vol_limit_UB.in(grid.nodes) <= 0);
//
//    Constraint Vol_limit_LB("Vol_limit_LB");
//    Vol_limit_LB = power(vr, 2) + power(vi, 2);
//    Vol_limit_LB -= power(grid.v_min,2);
//    ACOPF.add(Vol_limit_LB.in(grid.nodes) >= 0);
//    
//    
//    for (auto n:grid.nodes) {
//        auto b = (Bus*)n;
//        if (b->_gen.size()>=1) {
//            auto g = b->_gen.front();
//            if (b->_gen.size()>1) {
//                auto nb_gens = b->_gen.size();
//                double agg_pg_max = g->_pbound.max, agg_qg_max= g->_qbound.max, agg_pg_min = g->_pbound.min, agg_qg_min= g->_qbound.min;
//                for(unsigned i = 1; i< nb_gens; i++){
//                    b->p_injec -= Pg(b->_gen[i]->_name);
//                    b->q_injec -= Qg(b->_gen[i]->_name);
//                    agg_pg_max += b->_gen[i]->_pbound.max;
//                    agg_qg_max += b->_gen[i]->_qbound.max;
//                    agg_pg_min += b->_gen[i]->_pbound.min;
//                    agg_qg_min += b->_gen[i]->_qbound.min;
//                    
//                }
////                Constraint Agg_PGen_limit_UB("Agg_PGen_limit_UB"+b->_name);
////                Agg_PGen_limit_UB += b->p_injec;
////                ACOPF.add(Agg_PGen_limit_UB <= agg_pg_max);
////                
////                Constraint PGen_limit_LB("PGen_limit_LB"+b->_name);
////                PGen_limit_LB += b->p_injec;
////                ACOPF.add(PGen_limit_LB >= agg_pg_min);
////                
////                Constraint QGen_limit_UB("QGen_limit_UB"+b->_name);
////                QGen_limit_UB += b->q_injec;
////                ACOPF.add(QGen_limit_UB <= agg_qg_max);
////                
////                Constraint QGen_limit_LB("QGen_limit_LB"+b->_name);
////                QGen_limit_LB += b->q_injec;
////                ACOPF.add(QGen_limit_LB >= agg_qg_min);
//            }
//            Constraint PGen_limit_UB("PGen_limit_UB"+b->_name);
//            PGen_limit_UB += b->p_injec;
//            ACOPF.add(PGen_limit_UB <= g->_pbound.max);
//            
//            Constraint PGen_limit_LB("PGen_limit_LB"+b->_name);
//            PGen_limit_LB += b->p_injec;
//            ACOPF.add(PGen_limit_LB >= g->_pbound.min);
//            
//            Constraint QGen_limit_UB("QGen_limit_UB"+b->_name);
//            QGen_limit_UB += b->q_injec;
//            ACOPF.add(QGen_limit_UB <= g->_qbound.max);
//            
//            Constraint QGen_limit_LB("QGen_limit_LB"+b->_name);
//            QGen_limit_LB += b->q_injec;
//            ACOPF.add(QGen_limit_LB >= g->_qbound.min);
//            
//        }
//    }
//    
//    solver OPF(ACOPF,ipopt);
//    double solver_time_start = get_wall_time();
//    OPF.run(output, relax = false, tol = 1e-6, 0.01, "ma27", mehrotra = "no");
//    double solver_time_end = get_wall_time();
//    double total_time_end = get_wall_time();
//    auto solve_time = solver_time_end - solver_time_start;
//    auto total_time = total_time_end - total_time_start;
//
//    /** Terminal output */
//    string out = "DATA_OPF, " + grid._name + ", " + to_string(nb_buses) + ", " + to_string(nb_lines) +", " + to_string(ACOPF._obj_val) + ", " + to_string(-numeric_limits<double>::infinity()) + ", " + to_string(solve_time) + ", LocalOptimal, " + to_string(total_time);
//    DebugOn(out <<endl);
//    
//    
//    return 0;
//    
//    
////    ACOPF.min(obj.in(grid.gens));
////
////    /** Define constraints */
////    
////    /* REF BUS */
////    Constraint Ref_Bus("Ref_Bus");
////    if (polar) {
////        Ref_Bus = theta(grid.get_ref_bus());
////    }
////    else {
////        Ref_Bus = vi(grid.get_ref_bus());
////    }
////    ACOPF.add(Ref_Bus == 0);
////    
////    /** KCL Flow conservation */
////    Constraint KCL_P("KCL_P");
////    Constraint KCL_Q("KCL_Q");
////    KCL_P  = sum(Pf_from.out_arcs()) + sum(Pf_to.in_arcs()) + grid.pl - sum(Pg.in_gens());
////    KCL_Q  = sum(Qf_from.out_arcs()) + sum(Qf_to.in_arcs()) + grid.ql - sum(Qg.in_gens());
////    /* Shunts */
////    KCL_P +=  grid.gs*(power(vr,2)+power(vi,2));
////    KCL_Q -=  grid.bs*(power(vr,2)+power(vi,2));
////    
////    ACOPF.add(KCL_P.in(grid.nodes) == 0);
////    ACOPF.add(KCL_Q.in(grid.nodes) == 0);
////    
////    /** AC Power Flows */
////    /** TODO write the constraints in Complex form */
////    Constraint Flow_P_From("Flow_P_From");
////    Flow_P_From += Pf_from;
////    Flow_P_From -= grid.g_ff*(power(vr.from(), 2) + power(vi.from(), 2));
////    Flow_P_From -= grid.g_ft*(vr.from()*vr.to() + vi.from()*vi.to());
////    Flow_P_From -= grid.b_ft*(vi.from()*vr.to() - vr.from()*vi.to());
////    ACOPF.add(Flow_P_From.in(grid.arcs)==0);
////
////    Constraint Flow_P_To("Flow_P_To");
////    Flow_P_To += Pf_to;
////    Flow_P_To -= grid.g_tt*(power(vr.to(), 2) + power(vi.to(), 2));
////    Flow_P_To -= grid.g_tf*(vr.from()*vr.to() + vi.from()*vi.to());
////    Flow_P_To -= grid.b_tf*(vi.to()*vr.from() - vr.to()*vi.from());
////    ACOPF.add(Flow_P_To.in(grid.arcs)==0);
////
////    Constraint Flow_Q_From("Flow_Q_From");
////    Flow_Q_From += Qf_from;
////    if (polar) {
////        Flow_Q_From += (0.5*grid.ch+grid.b)/power(grid.tr,2)*power(v.from(),2);
////        Flow_Q_From -= grid.b/grid.tr*(v.from()*v.to()*cos(theta.from() - theta.to() - grid.as));
////        Flow_Q_From += grid.g/grid.tr*(v.from()*v.to()*sin(theta.from() - theta.to() - grid.as));
////    }
////    else {
////        Flow_Q_From += grid.b_ff*(power(vr.from(), 2) + power(vi.from(), 2));
////        Flow_Q_From += grid.b_ft*(vr.from()*vr.to() + vi.from()*vi.to());
////        Flow_Q_From -= grid.g_ft*(vi.from()*vr.to() - vr.from()*vi.to());
////    }
////    ACOPF.add(Flow_Q_From.in(grid.arcs)==0);
////    Constraint Flow_Q_To("Flow_Q_To");
////    Flow_Q_To += Qf_to;
////    if (polar) {
////        Flow_Q_To += (0.5*grid.ch+grid.b)*power(v.to(),2);
////        Flow_Q_To -= grid.b/grid.tr*(v.from()*v.to()*cos(theta.to() - theta.from() + grid.as));
////        Flow_Q_To += grid.g/grid.tr*(v.from()*v.to()*sin(theta.to() - theta.from() + grid.as));
////    }
////    else {
////        Flow_Q_To += grid.b_tt*(power(vr.to(), 2) + power(vi.to(), 2));
////        Flow_Q_To += grid.b_tf*(vr.from()*vr.to() + vi.from()*vi.to());
////        Flow_Q_To -= grid.g_tf*(vi.to()*vr.from() - vr.to()*vi.from());
////    }
////    ACOPF.add(Flow_Q_To.in(grid.arcs)==0);
////    
////    /** AC voltage limit constraints. */
////    if (!polar) {
////        Constraint Vol_limit_UB("Vol_limit_UB");
////        Vol_limit_UB = power(vr, 2) + power(vi, 2);
////        Vol_limit_UB -= power(grid.v_max, 2);
////        ACOPF.add(Vol_limit_UB.in(grid.nodes) <= 0);
////
////        Constraint Vol_limit_LB("Vol_limit_LB");
////        Vol_limit_LB = power(vr, 2) + power(vi, 2);
////        Vol_limit_LB -= power(grid.v_min,2);
////        ACOPF.add(Vol_limit_LB.in(grid.nodes) >= 0);
////        DebugOff(grid.v_min.to_str(true) << endl);
////        DebugOff(grid.v_max.to_str(true) << endl);
////    }
////
////    
////    /* Phase Angle Bounds constraints */
////    Constraint PAD_UB("PAD_UB");
////    Constraint PAD_LB("PAD_LB");
////    auto bus_pairs = grid.get_bus_pairs();
////    if (polar) {
////        PAD_UB = theta.from() - theta.to();
////        PAD_UB -= grid.th_max;
////        PAD_LB = theta.from() - theta.to();
////        PAD_LB -= grid.th_min;
////        DebugOff(grid.th_min.to_str(true) << endl);
////        DebugOff(grid.th_max.to_str(true) << endl);
////    }
////    else {        
////        DebugOff("Number of bus_pairs = " << bus_pairs.size() << endl);
////        PAD_UB = vi.from()*vr.to() - vr.from()*vi.to();
////        PAD_UB -= grid.tan_th_max*(vr.from()*vr.to() + vi.from()*vi.to());
////        
////        PAD_LB = vi.from()*vr.to() - vr.from()*vi.to();
////        PAD_LB -= grid.tan_th_min*(vr.from()*vr.to() + vi.from()*vi.to());
////        DebugOff(grid.th_min.to_str(true) << endl);
////        DebugOff(grid.th_max.to_str(true) << endl);
////    }
//////    ACOPF.add(PAD_UB.in(bus_pairs) <= 0);
//////    ACOPF.add(PAD_LB.in(bus_pairs) >= 0);
////
////
////    /*  Thermal Limit Constraints */
////    Constraint Thermal_Limit_from("Thermal_Limit_from");
////    Thermal_Limit_from += power(Pf_from, 2) + power(Qf_from, 2);
////    Thermal_Limit_from -= power(grid.S_max, 2);
//////    ACOPF.add(Thermal_Limit_from.in(grid.arcs) <= 0);
////
////    Constraint Thermal_Limit_to("Thermal_Limit_to");
////    Thermal_Limit_to += power(Pf_to, 2) + power(Qf_to, 2);
////    Thermal_Limit_to -= power(grid.S_max,2);
//////    ACOPF.add(Thermal_Limit_to.in(grid.arcs) <= 0);
////    DebugOff(grid.S_max.in(grid.arcs).to_str(true) << endl);
////    
////    solver OPF(ACOPF,ipopt);
////    double solver_time_start = get_wall_time();
////    OPF.run(output, relax = false, tol = 1e-6, "ma27", mehrotra = "no");
////    double solver_time_end = get_wall_time();
////    double total_time_end = get_wall_time();
////    auto solve_time = solver_time_end - solver_time_start;
////    auto total_time = total_time_end - total_time_start;
////
////    /** Terminal output */
////    string out = "DATA_OPF, " + grid._name + ", " + to_string(nb_buses) + ", " + to_string(nb_lines) +", " + to_string(ACOPF._obj_val) + ", " + to_string(-numeric_limits<double>::infinity()) + ", " + to_string(solve_time) + ", LocalOptimal, " + to_string(total_time);
////    DebugOn(out <<endl);
//    
//    return 0;
//}

int main (int argc, char * argv[])
{
    string fname = "../data_sets/Power/nesta_case5_pjm.m", mtype = "ACPOL";
    DebugOn("argv[0] =" << argv[0] << endl);
    string path = argv[0];
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
    grid.readgrid(fname.c_str());
    
    /* Grid Parameters */
    unsigned nb_gen = grid.get_nb_active_gens();
    unsigned nb_lines = grid.get_nb_active_arcs();
    unsigned nb_buses = grid.get_nb_active_nodes();
    
    
    DebugOn("nb gens = " << nb_gen << endl);
    DebugOn("nb lines = " << nb_lines << endl);
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
    /* Power generation variables */
    var<double> Pg("Pg", grid.pg_min, grid.pg_max);
    var<double> Qg ("Qg", grid.qg_min, grid.qg_max);
    ACOPF.add(Pg.in(grid.gens));
    ACOPF.add(Qg.in(grid.gens));
    
    /* Power flow variables */
    var<double> Pf_from("Pf_from", grid.S_max);
    var<double> Qf_from("Qf_from", grid.S_max);
    var<double> Pf_to("Pf_to", grid.S_max);
    var<double> Qf_to("Qf_to", grid.S_max);
    
    ACOPF.add(Pf_from.in(grid.arcs));
    ACOPF.add(Qf_from.in(grid.arcs));
    ACOPF.add(Pf_to.in(grid.arcs));
    ACOPF.add(Qf_to.in(grid.arcs));
    
    
    /** Voltage related variables */
    var<double> theta("theta");
    var<double> v("|V|", grid.v_min, grid.v_max);
//    var<double> vr("vr", grid.v_max);
//    var<double> vi("vi", grid.v_max);
    var<double> vr("vr");
    var<double> vi("vi");
    
    if (polar) {
        ACOPF.add(v.in(grid.nodes));
        ACOPF.add(theta.in(grid.nodes));
        v.initialize_all(1);
    }
    else {
        ACOPF.add(vr.in(grid.nodes));
        ACOPF.add(vi.in(grid.nodes));
        vr.initialize_all(1.0);
    }
    
    /** Construct the objective function */
    func_ obj = product(grid.c1, Pg) + product(grid.c2, power(Pg,2)) + sum(grid.c0);
    ACOPF.min(obj.in(grid.gens));
    
    /** Define constraints */
    
    /* REF BUS */
    Constraint Ref_Bus("Ref_Bus");
    if (polar) {
        Ref_Bus = theta(grid.get_ref_bus());
    }
    else {
        Ref_Bus = vi(grid.get_ref_bus());
    }
    ACOPF.add(Ref_Bus == 0);
    
    /** KCL Flow conservation */
    Constraint KCL_P("KCL_P");
    Constraint KCL_Q("KCL_Q");
    KCL_P  = sum(Pf_from.out_arcs()) + sum(Pf_to.in_arcs()) + grid.pl - sum(Pg.in_gens());
    KCL_Q  = sum(Qf_from.out_arcs()) + sum(Qf_to.in_arcs()) + grid.ql - sum(Qg.in_gens());
    /* Shunts */
    if (polar) {
        KCL_P +=  grid.gs*power(v,2);
        KCL_Q -=  grid.bs*power(v,2);
    }
    else {
        KCL_P +=  grid.gs*(power(vr,2)+power(vi,2));
        KCL_Q -=  grid.bs*(power(vr,2)+power(vi,2));
    }
    ACOPF.add(KCL_P.in(grid.nodes) == 0);
    ACOPF.add(KCL_Q.in(grid.nodes) == 0);
    
    /** AC Power Flows */
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
    ACOPF.add(Flow_P_From.in(grid.arcs)==0);
    
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
    ACOPF.add(Flow_P_To.in(grid.arcs)==0);
    
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
    ACOPF.add(Flow_Q_From.in(grid.arcs)==0);
    
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
    ACOPF.add(Flow_Q_To.in(grid.arcs)==0);
    
    /** AC voltage limit constraints. */
    if (!polar) {
        Constraint Vol_limit_UB("Vol_limit_UB");
        Vol_limit_UB = power(vr, 2) + power(vi, 2);
        Vol_limit_UB -= power(grid.v_max, 2);
        ACOPF.add(Vol_limit_UB.in(grid.nodes) <= 0);
        
        Constraint Vol_limit_LB("Vol_limit_LB");
        Vol_limit_LB = power(vr, 2) + power(vi, 2);
        Vol_limit_LB -= power(grid.v_min,2);
        ACOPF.add(Vol_limit_LB.in(grid.nodes) >= 0);
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
    }
    else {
        DebugOff("Number of bus_pairs = " << bus_pairs.size() << endl);
        PAD_UB = vi.from()*vr.to() - vr.from()*vi.to();
        PAD_UB -= grid.tan_th_max*(vr.from()*vr.to() + vi.from()*vi.to());
        
        PAD_LB = vi.from()*vr.to() - vr.from()*vi.to();
        PAD_LB -= grid.tan_th_min*(vr.from()*vr.to() + vi.from()*vi.to());
    }
        ACOPF.add_lazy(PAD_UB.in(bus_pairs) <= 0);
        ACOPF.add_lazy(PAD_LB.in(bus_pairs) >= 0);
    
    
    /*  Thermal Limit Constraints */
    Constraint Thermal_Limit_from("Thermal_Limit_from");
    Thermal_Limit_from += power(Pf_from, 2) + power(Qf_from, 2);
    Thermal_Limit_from -= power(grid.S_max, 2);
    ACOPF.add(Thermal_Limit_from.in(grid.arcs) <= 0);
    
    Constraint Thermal_Limit_to("Thermal_Limit_to");
    Thermal_Limit_to += power(Pf_to, 2) + power(Qf_to, 2);
    Thermal_Limit_to -= power(grid.S_max,2);
    ACOPF.add(Thermal_Limit_to.in(grid.arcs) <= 0);
    
    solver OPF(ACOPF,ipopt);
    double solver_time_start = get_wall_time();
    OPF.run(output, relax = false, tol = 1e-6, 0.02, "ma27", mehrotra = "no");
    double solver_time_end = get_wall_time();
    double total_time_end = get_wall_time();
    auto solve_time = solver_time_end - solver_time_start;
    auto total_time = total_time_end - total_time_start;
    
//    ACOPF.print_expanded();
    /** Terminal output */
    string out = "DATA_OPF, " + grid._name + ", " + to_string(nb_buses) + ", " + to_string(nb_lines) +", " + to_string(ACOPF._obj_val) + ", " + to_string(-numeric_limits<double>::infinity()) + ", " + to_string(solve_time) + ", LocalOptimal, " + to_string(total_time);
    DebugOn(out <<endl);
    
    return 0;
}
