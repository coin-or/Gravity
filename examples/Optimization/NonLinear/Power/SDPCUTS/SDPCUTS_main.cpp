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
    bool loss_from = false, llnc=true;
    size_t num_bags = 0;
    string num_bags_s = "100";
    string solver_str = "ipopt";
    string sdp_cuts_s = "yes";
    string loss_from_s = "yes";
    string lazy_s = "no";
    bool lazy_bool = false;
    SolverType solv_type = ipopt;
    double tol = 1e-6;
    string mehrotra = "no";
    
    //string fname = string(prj_dir)+"/data_sets/Power/nesta_case5_pjm.m";
    
    string fname = string(prj_dir)+"/data_sets/Power/nesta_case5_pjm.m";
    
    // create a OptionParser with options
    op::OptionParser opt;
    opt.add_option("h", "help",
                   "shows option help"); // no default value means boolean options, which default value is false
    opt.add_option("f", "file", "Input file name", fname);
    opt.add_option("s", "solver", "Solvers: ipopt/cplex/gurobi, default = ipopt", solver_str);
    opt.add_option("b", "numbags", "Number of bags per iteration", num_bags_s);
    opt.add_option("l", "losses", "add loss constraints", loss_from_s);
    opt.add_option("lz", "lazy", "Generate 3d SDP cuts in a lazy fashion, default = no", lazy_s);
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
    lazy_s = opt["lz"];
    if (lazy_s.compare("no")==0) {
        lazy_bool = false;
    }
    
    loss_from_s = opt["l"];
    if (loss_from_s.compare("no")==0) {
        loss_from = false;
    }
    else {
        loss_from = true;
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
    auto ch_half=grid.ch_half.in(arcs);


    
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
    var<> lij("lij");
    if(loss_from){
        SDP.add(lij.in(arcs));
    }
    SDP.add(Pf_from.in(arcs), Qf_from.in(arcs),Pf_to.in(arcs),Qf_to.in(arcs));
    
    
      //Real part of Wij = ViVj
    var<>  R_Wij("R_Wij", wr_min, wr_max);
     //Imaginary part of Wij = ViVj
    var<>  Im_Wij("Im_Wij", wi_min, wi_max);
     //Magnitude of Wii = Vi^2
    

    
    var<>  Wii("Wii", w_min, w_max);
    SDP.add(Wii.in(nodes));
    SDP.add(R_Wij.in(bus_pairs_chord));
    SDP.add(Im_Wij.in(bus_pairs_chord));
    
    
    // Initialize variables
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
    Flow_P_From = Pf_from - (g_ff*Wii.from(arcs) + g_ft*R_Wij.in(arcs) + b_ft*Im_Wij.in(arcs));
    SDP.add(Flow_P_From.in(arcs) == 0);
    
    Constraint<> Flow_P_To("Flow_P_To");
    Flow_P_To = Pf_to - (g_tt*Wii.to(arcs) + g_tf*R_Wij.in(arcs) - b_tf*Im_Wij.in(arcs));
    SDP.add(Flow_P_To.in(arcs) == 0);
    
    Constraint<> Flow_Q_From("Flow_Q_From");
    Flow_Q_From = Qf_from - (g_ft*Im_Wij.in(arcs) - b_ff*Wii.from(arcs) - b_ft*R_Wij.in(arcs));
    SDP.add(Flow_Q_From.in(arcs) == 0);
    
    Constraint<> Flow_Q_To("Flow_Q_To");
    Flow_Q_To = Qf_to + b_tt*Wii.to(arcs) + b_tf*R_Wij.in(arcs) + g_tf*Im_Wij.in(arcs);
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
    
    if(loss_from){
        param<Cpx> T("T"), Y("Y"), Ych("Ych");
        var<Cpx> L_from("L_from"), W("W");
        T.real_imag(cc.in(arcs), dd.in(arcs));
        Y.real_imag(g.in(arcs), b.in(arcs));
        Ych.set_imag(ch_half.in(arcs));


        L_from.set_real(lij.in(arcs));
        W.real_imag(R_Wij.in(arcs), Im_Wij.in(arcs));
  
        Constraint<Cpx> I_from("I_from");
        I_from=(Y+Ych)*(conj(Y)+conj(Ych))*Wii.from(arcs)-T*Y*(conj(Y)+conj(Ych))*conj(W)-conj(T)*conj(Y)*(Y+Ych)*W+pow(tr,2)*Y*conj(Y)*Wii.to(arcs);
        SDP.add_real(I_from.in(arcs)==pow(tr,2)*L_from.in(arcs));

    
        Constraint<> I_from_L("I_from_L");
        I_from_L = (pow(Pf_from, 2) + pow(Qf_from, 2))*pow(tr,2)-w_max.from(arcs)*lij;
        SDP.add(I_from_L.in(arcs) <= 0);
    

        
        Constraint<> I_from_U("I_from_U");
        I_from_U = w_min.from(arcs)*lij - pow(tr,2)*pow(S_max,2);
        SDP.add(I_from_U.in(arcs) <= 0);
    
    }
    
        
      
    
    
    
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
    
    return 0;

}

