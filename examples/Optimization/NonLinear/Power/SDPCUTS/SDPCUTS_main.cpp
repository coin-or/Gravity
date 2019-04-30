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
    bool relax = false, sdp_cuts = true, soc=true;
    bool loss = true, llnc=false;
    size_t num_bags = 0;
    string num_bags_s = "100";
    string solver_str = "ipopt";
    string sdp_cuts_s = "yes";
    string lazy = "no";
    bool lazy_bool = false;
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
    
    /* Grid Stats */
    auto nb_gen = grid.get_nb_active_gens();
    auto nb_lines = grid.get_nb_active_arcs();
    auto nb_buses = grid.get_nb_active_nodes();
    DebugOn("nb active gens = " << nb_gen << endl);
    DebugOn("nb active lines = " << nb_lines << endl);
    DebugOn("nb active buses = " << nb_buses << endl);
    
    /** Sets */
    auto bus_pairs = grid.get_bus_pairs();
    auto bus_pairs_chord = grid.get_bus_pairs_chord();
    if (!sdp_cuts) {
        bus_pairs_chord = bus_pairs;
    }
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
    auto pl = grid.pl.in(nodes);
    auto ql = grid.ql.in(nodes);
    auto gs = grid.gs.in(nodes);
    auto bs = grid.bs.in(nodes);
    auto b = grid.b.in(arcs);
    auto g = grid.g.in(arcs);
    auto as = grid.as.in(arcs);
    auto ch = grid.ch.in(arcs);
    auto tr = grid.tr.in(arcs);
    auto th_min = grid.th_min.in(bus_pairs);
    auto th_max = grid.th_max.in(bus_pairs);
    auto g_ft = grid.g_ft.in(arcs);
    auto g_ff = grid.g_ff.in(arcs);
    auto g_tt = grid.g_tt.in(arcs);
    auto g_tf = grid.g_tf.in(arcs);
    auto b_ft = grid.b_ft.in(arcs);
    auto b_ff = grid.b_ff.in(arcs);
    auto b_tf = grid.b_tf.in(arcs);
    auto b_tt = grid.b_tt.in(arcs);
    auto S_max = grid.S_max.in(arcs);
    auto w_max = grid.w_max.in(nodes);
    auto w_min = grid.w_min.in(nodes);
    auto tan_th_min = grid.tan_th_min.in(bus_pairs);
    auto tan_th_max = grid.tan_th_max.in(bus_pairs);
    auto wr_min = grid.wr_min.in(bus_pairs_chord);
    auto wr_max = grid.wr_max.in(bus_pairs_chord);
    auto wi_min = grid.wi_min.in(bus_pairs_chord);
    auto wi_max = grid.wi_max.in(bus_pairs_chord);
    auto cc=grid.cc.in(arcs);
    auto dd=grid.dd.in(arcs);
    auto bt=grid.bt.in(arcs);
    auto cht=grid.cht.in(arcs);
    auto ch_half=grid.ch_half.in(arcs);
    auto cht_half=grid.cht_half.in(arcs);
    /*auto Y=grid.Y.in(arcs);
     auto Ych=grid.Ych.in(arcs);*/
    //    auto V_sq_max=1.21;
    auto V_sq_min=0.81;
    
    double upper_bound = grid.solve_acopf();
    
    
    /** Build model */
    Model<> SDP("SDP Model");
    
    /** Variables */
    /* Power generation variables */
    var<> Pg("Pg", pg_min, pg_max);
    var<> Qg ("Qg", qg_min, qg_max);
    SDP.add(Pg.in(gens),Qg.in(gens));
    
    /* Power flow variables */
    var<> Pf_from("Pf_from", -1.*S_max,S_max);
    var<> Qf_from("Qf_from", -1.*S_max,S_max);
    var<> Pf_to("Pf_to", -1.*S_max,S_max);
    var<> Qf_to("Qf_to", -1.*S_max,S_max);
    var<> lij("lij", pos_);
    if(loss){
        SDP.add(lij.in(arcs));
    }
    SDP.add(Pf_from.in(arcs), Qf_from.in(arcs),Pf_to.in(arcs),Qf_to.in(arcs));
    
    /*
     //Real part of Wij = ViVj
     var<>  R_Wij("R_Wij", -1.21, 1.21);
     //Imaginary part of Wij = ViVj
     var<>  Im_Wij("Im_Wij", -1.21, 1.21);
     //Magnitude of Wii = Vi^2
     */
    
    /* Real part of Wij = ViVj */
    var<>  R_Wij("R_Wij", wr_min, wr_max);
    /* Imaginary part of Wij = ViVj */
    var<>  Im_Wij("Im_Wij", wi_min, wi_max);
    /* Magnitude of Wii = Vi^2 */
    
    var<>  Wii("Wii", w_min, w_max);
    SDP.add(Wii.in(nodes));
    SDP.add(R_Wij.in(bus_pairs_chord));
    SDP.add(Im_Wij.in(bus_pairs_chord));
    
    
    /* Initialize variables */
    R_Wij.initialize_all(1.0);
    Wii.initialize_all(1.001);
    
    /**  Objective */
    auto obj = product(c1,Pg) + product(c2,pow(Pg,2)) + sum(c0);
    SDP.min(obj);
    
    //    Constraint obj_cstr("obj_ub");
    //    obj_cstr += SDP._obj - 1.05*upper_bound;
    //    SDP.add(obj_cstr <= 0);
    
    
    /** Constraints */
    if(grid.add_3d_nlin && sdp_cuts) {
        auto bag_size = grid._bags.size();
        DebugOn("\nNum of bags = " << bag_size << endl);
        DebugOn("Adding 3d determinant polynomial cuts\n");
        auto R_Wij_ = R_Wij.pairs_in_bags(grid._bags, 3);
        auto Im_Wij_ = Im_Wij.pairs_in_bags(grid._bags, 3);
        auto Wii_ = Wii.in_bags(grid._bags, 3);
        
        Constraint<> SDP3("SDP_3D");
        SDP3 = 2 * R_Wij_[0] * (R_Wij_[1] * R_Wij_[2] + Im_Wij_[1] * Im_Wij_[2]);
        SDP3 -= 2 * Im_Wij_[0] * (R_Wij_[2] * Im_Wij_[1] - Im_Wij_[2] * R_Wij_[1]);
        SDP3 -= (pow(R_Wij_[0], 2) + pow(Im_Wij_[0], 2)) * Wii_[2];
        SDP3 -= (pow(R_Wij_[1], 2) + pow(Im_Wij_[1], 2)) * Wii_[0];
        SDP3 -= (pow(R_Wij_[2], 2) + pow(Im_Wij_[2], 2)) * Wii_[1];
        SDP3 += Wii_[0] * Wii_[1] * Wii_[2];
        if (lazy_bool) {
            SDP.add_lazy(SDP3 >= 0);
        }
        else {
            SDP.add(SDP3 >= 0);
            DebugOn("Number of 3d determinant cuts = " << SDP3.get_nb_instances() << endl);
        }
    }
    
    /** Constraints */
    /* Second-order cone constraints */
    if(soc)
    {
        Constraint<> SOC("SOC");
        SOC = pow(R_Wij, 2) + pow(Im_Wij, 2) - Wii.from(bus_pairs_chord)*Wii.to(bus_pairs_chord);
        SDP.add(SOC.in(bus_pairs_chord) <= 0);
    }
    
    /* Flow conservation */
    Constraint<> KCL_P("KCL_P");
    KCL_P  = sum(Pf_from, out_arcs) + sum(Pf_to, in_arcs) + pl - sum(Pg, gen_nodes) + gs*Wii;
    SDP.add(KCL_P.in(nodes) == 0);
    
    Constraint<> KCL_Q("KCL_Q");
    KCL_Q  = sum(Qf_from, out_arcs) + sum(Qf_to, in_arcs) + ql - sum(Qg, gen_nodes) - bs*Wii;
    SDP.add(KCL_Q.in(nodes) == 0);
    
    /* AC Power Flow */
    Constraint<> Flow_P_From("Flow_P_From");
    Flow_P_From = Pf_from - (g_ff*Wii.from(arcs) + g_ft*R_Wij + b_ft*Im_Wij);
    SDP.add(Flow_P_From.in(arcs) == 0);
    
    Constraint<> Flow_P_To("Flow_P_To");
    Flow_P_To = Pf_to - (g_tt*Wii.to(arcs) + g_tf*R_Wij - b_tf*Im_Wij);
    SDP.add(Flow_P_To.in(arcs) == 0);
    
    Constraint<> Flow_Q_From("Flow_Q_From");
    Flow_Q_From = Qf_from - (g_ft*Im_Wij - b_ff*Wii.from(arcs) - b_ft*R_Wij);
    SDP.add(Flow_Q_From.in(arcs) == 0);
    
    Constraint<> Flow_Q_To("Flow_Q_To");
    Flow_Q_To = Qf_to + b_tt*Wii.to(arcs) + b_tf*R_Wij + g_tf*Im_Wij;
    SDP.add(Flow_Q_To.in(arcs) == 0);
    
    /* Phase Angle Bounds constraints */
    Constraint<> PAD_UB("PAD_UB");
    PAD_UB = Im_Wij.in(bus_pairs);
    PAD_UB <= tan_th_max*R_Wij.in(bus_pairs);
    SDP.add_lazy(PAD_UB.in(bus_pairs));
    
    Constraint<> PAD_LB("PAD_LB");
    PAD_LB =  Im_Wij.in(bus_pairs);
    PAD_LB >= tan_th_min*R_Wij.in(bus_pairs);
    SDP.add_lazy(PAD_LB.in(bus_pairs));
    
    /* Thermal Limit Constraints */
    Constraint<> Thermal_Limit_from("Thermal_Limit_from");
    Thermal_Limit_from = pow(Pf_from, 2) + pow(Qf_from, 2);
    Thermal_Limit_from <= pow(S_max,2);
    SDP.add(Thermal_Limit_from.in(arcs));
    
    
    Constraint<> Thermal_Limit_to("Thermal_Limit_to");
    Thermal_Limit_to = pow(Pf_to, 2) + pow(Qf_to, 2);
    Thermal_Limit_to <= pow(S_max,2);
    SDP.add(Thermal_Limit_to.in(arcs));
    
    
    param<Cpx> T("T"), Y("Y"), Yt("Yt"), Ych("Ych"), Ycht("Ycht");
    var<Cpx> S("S"), L("L"), W("W");
    T.real_imag(cc.in(arcs), dd.in(arcs));
    Yt.real_imag(g_tt.in(arcs), bt.in(arcs));
    Y.real_imag(g.in(arcs), b.in(arcs));
    Ych.set_imag(ch_half.in(arcs));
    Ycht.set_imag(cht_half.in(arcs));
    
    /* ba.set_real(b.in(arcs)*dd.in(arcs));
     ga.set_real(g.in(arcs));
     cha.set_real(ch.in(arcs));*/
    
    L.set_real(lij.in(arcs));
    W.real_imag(R_Wij.in(arcs), Im_Wij.in(arcs));
    
    Constraint<Cpx> C1("C1");
    //        C1 = (Yt+Ycht)*(conj(Y)+conj(Ych))*Wii.from(arcs)-(conj(Yt)+conj(Ycht))*Y*T*conj(W)-conj(Yt)*(Y+Ych)*conj(T)*W+Y*conj(Y)*Wii.to(arcs);
    //        SDP.add(C1==L);
    
    
    C1 = (Yt+Ycht)*(conj(Y)+conj(Ych))*Wii.from(arcs)-(conj(Yt)+conj(Ycht))*Y*conj(W)-conj(Y)*(Yt+Ycht)*W+Y*conj(Y)*Wii.to(arcs);
    SDP.add_real(C1==L);
    //         C1 = (pow(ga,2)+(ba+cha)*(ba+cha))*Wii.from(arcs)-(conj(Yt)+conj(Ycht))*Y*conj(W)-conj(Yt)*(Y+Ych)*W+(pow(ga,2)+pow(ba,2))*Wii.to(arcs);
    //        SDP.add(C1.in(arcs)==L.in(arcs));
    
    Constraint<> Loss_from("Loss_from");
    Loss_from = (pow(Pf_from, 2) + pow(Qf_from, 2))*pow(tr,2)-w_max.from(arcs)*lij;
    SDP.add(Loss_from.in(arcs) <= 0);
    
    Constraint<> Loss_to("Loss_to");
    Loss_to = (pow(Pf_to, 2) + pow(Qf_to, 2))*pow(tr,2)-w_max.to(arcs)*lij;
    SDP.add(Loss_to.in(arcs) <= 0);
    
    Constraint<> Loss_U("Loss_U");
    Loss_U = w_min.from(arcs)*lij - pow(tr,2)*pow(S_max,2);
    SDP.add(Loss_U.in(arcs) <= 0);
    
    
    
    Constraint<> Loss_C("Loss_C");
    
    Loss_C = lij - pow(tr,2)*((pow(g_ff,2)+pow(b_ff,2))*Wii.from(arcs) + (pow(g_ft,2)+pow(b_ft,2))*Wii.to(arcs)
                              +(g_ff*g_ft+b_ff*b_ft)*2*R_Wij +(g_ff*b_ft-b_ff*g_ft)*2*Im_Wij);
    
    //SDP.add(Loss_C.in(arcs) == 0);
    
    
    Constraint<> Loss_C1("Loss_C1");
    
    Loss_C1 = lij - ((pow(g,2)+pow(b+ch*0.5,2))*Wii.from(arcs) + (pow(g,2)+pow(b,2))*Wii.to(arcs)
                     -2*((g*g+b*b+b*ch*0.5)*R_Wij +g*ch*0.5*Im_Wij));
    
    //SDP.add(Loss_C1.in(arcs) == 0);
    
    
    
    
    
    if (llnc)
    {
        /* Lifted Nonlinear Cuts */
        Constraint<> LNC1("LNC1");
        LNC1 += (grid.v_min.from(bus_pairs)+grid.v_max.from(bus_pairs))*(grid.v_min.to(bus_pairs)+grid.v_max.to(bus_pairs))*(grid.sphi*Im_Wij + grid.cphi*R_Wij);
        LNC1 -= grid.v_max.to(bus_pairs)*grid.cos_d*(grid.v_min.to(bus_pairs)+grid.v_max.to(bus_pairs))*Wii.from(bus_pairs);
        LNC1 -= grid.v_max.from(bus_pairs)*grid.cos_d*(grid.v_min.from(bus_pairs)+grid.v_max.from(bus_pairs))*Wii.to(bus_pairs);
        LNC1 -= grid.v_max.from(bus_pairs)*grid.v_max.to(bus_pairs)*grid.cos_d*(grid.v_min.from(bus_pairs)*grid.v_min.to(bus_pairs) - grid.v_max.from(bus_pairs)*grid.v_max.to(bus_pairs));
        SDP.add_lazy(LNC1.in(bus_pairs) >= 0);
        
        Constraint<> LNC2("LNC2");
        LNC2 += (grid.v_min.from(bus_pairs)+grid.v_max.from(bus_pairs))*(grid.v_min.to(bus_pairs)+grid.v_max.to(bus_pairs))*(grid.sphi*Im_Wij + grid.cphi*R_Wij);
        LNC2 -= grid.v_min.to(bus_pairs)*grid.cos_d*(grid.v_min.to(bus_pairs)+grid.v_max.to(bus_pairs))*Wii.from(bus_pairs);
        LNC2 -= grid.v_min.from(bus_pairs)*grid.cos_d*(grid.v_min.from(bus_pairs)+grid.v_max.from(bus_pairs))*Wii.to(bus_pairs);
        LNC2 += grid.v_min.from(bus_pairs)*grid.v_min.to(bus_pairs)*grid.cos_d*(grid.v_min.from(bus_pairs)*grid.v_min.to(bus_pairs) - grid.v_max.from(bus_pairs)*grid.v_max.to(bus_pairs));
        SDP.add_lazy(LNC2.in(bus_pairs) >= 0);
    }
    SDP.print();
    
    
    
    total_time_start = get_wall_time();
    /* Solver selection */
    solver<> SDPOPF(SDP,solv_type);
    double solver_time_start = get_wall_time();
    SDPOPF.run(output = 5, tol = 1e-6);
    //    SDP.print();
    double gap = 100*(upper_bound - SDP.get_obj_val())/upper_bound;
    double solver_time_end = get_wall_time();
    double total_time_end = get_wall_time();
    auto solve_time = solver_time_end - solver_time_start;
    auto total_time = total_time_end - total_time_start;
    string out = "\nDATA_OPF, " + grid._name + ", " + to_string(nb_buses) + ", " + to_string(nb_lines) +", " + to_string(SDP.get_obj_val()) + ", " + to_string(-numeric_limits<double>::infinity()) + ", " + to_string(solve_time) + ", LocalOptimal, " + to_string(total_time);
    DebugOn(out <<endl);
    DebugOn("Final Gap = " << to_string(gap) << "%."<<endl);
    DebugOn("Upper bound = " << to_string(upper_bound) << "."<<endl);
    DebugOn("Lower bound = " << to_string(SDP.get_obj_val()) << "."<<endl);
    DebugOn("\nResults: " << grid._name << " " << to_string(SDP.get_obj_val()) << " " << to_string(total_time)<<endl);
    //    DebugOn("\nTime in nfp: " << time_in_all_nfp << endl);
    
    
    for(auto it:*(R_Wij.get_keys()))
        DebugOn("R_Wij "<<R_Wij.eval(it)<<endl);
    
    for(auto it:*lij.get_vals())
        DebugOn("LIj "<<it<<endl);
    
    
    /*Loss_C.print();
     Loss_C.print_symbolic();*/
    
    
    
    g_ff.print();
    g.print();
    b_ff.print();
    auto c_E=b+ch*0.5;
    c_E.print();
    
    auto a=pow(tr,2)*(pow(g_ff,2)+pow(b_ff,2));
    a.print();
    
    auto b_E=(pow(g,2)+pow(b+ch*0.5,2));
    b_E.print();
    
    tr.print();
    
    //C1.print();
    
    w_min.print();
    w_max.print();
    
    
    for(size_t it=0; it<Loss_to.get_nb_instances(); it++)
        DebugOn("C "<<Loss_to.eval(it)<<endl);
    
    
    
    return 0;
    
}

