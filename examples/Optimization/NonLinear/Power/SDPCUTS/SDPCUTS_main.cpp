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
    bool sdp_cuts = true;
    bool current_from = true, llnc=true, current_to=true, loss=false, loss_bounds=true;
    size_t num_bags = 0;
    string num_bags_s = "100";
    string solver_str = "ipopt";
    string sdp_cuts_s = "yes";
    string current_from_s = "yes";
    string orig_s = "yes";
    string current_to_s="yes";
    string lazy_s = "yes";
    bool lazy_bool = true;
    SolverType solv_type = ipopt;
    double tol = 1e-6;
    string mehrotra = "no";
    double ilb, iub;
    
    //string fname = string(prj_dir)+"/data_sets/Power/nesta_case5_pjm.m";
    
    string fname = string(prj_dir)+"/data_sets/Power/nesta_case5_pjm.m";
    
    // create a OptionParser with options
    op::OptionParser opt;
    opt.add_option("h", "help",
                   "shows option help"); // no default value means boolean options, which default value is false
    opt.add_option("f", "file", "Input file name", fname);
    opt.add_option("s", "solver", "Solvers: ipopt/cplex/gurobi, default = ipopt", solver_str);
    opt.add_option("b", "numbags", "Number of bags per iteration", num_bags_s);
    opt.add_option("If", "current_from", "add from current constraints", current_from_s);
    opt.add_option("o", "original", "add original variables and linking constraints", orig_s);
    opt.add_option("It", "current_to", "add to current constraints", current_to_s);
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
    
    current_from_s = opt["If"];
    if (current_from_s.compare("no")==0) {
        current_from = false;
    }
    else {
        current_from = true;
    }
    bool add_original = true;
    
    orig_s = opt["o"];
    if (orig_s.compare("no")==0) {
        add_original = false;
    }
    else {
        add_original = true;
    }
    
    
    
    current_to_s = opt["It"];
    if (current_to_s.compare("no")==0) {
        current_to = false;
    }
    else {
        current_to = true;
    }
    
    current_from=true;
    current_to=true;
//    loss=true;
    
    num_bags = atoi(opt["b"].c_str());
    
    cout << "\nnum bags = " << num_bags << endl;
    
    double total_time_start = get_wall_time();
    PowerNet grid;
    grid.readgrid(fname);
    grid.update_ref_bus();
    
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
    auto v_max = grid.v_max.in(nodes);
    auto v_min = grid.v_min.in(nodes);
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
    
    
    
    double upper_bound = grid.solve_acopf(ACRECT);
    
    
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
    var<> lji("lji", pos_);
    if(current_from){
        SDP.add(lij.in(arcs));
    }
    
    if(current_from){
        SDP.add(lji.in(arcs));
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
    
    
    var<>  R_Vi("R_Vi", -1*v_max, v_max);
    var<>  Im_Vi("Im_Vi", -1*v_max, v_max);
    
        //var<> sWf_sWt("sWf_sWt", w_min, w_max), Wf_Wt("Wf_Wt", w_min*w_min, w_max*w_max);
    func<> theta_L1;
    theta_L1.allocate_mem();
    theta_L1=(acos(R_Wij.get_lb().in(bus_pairs)/sqrt(Wii.get_ub().from(bus_pairs)*Wii.get_ub().to(bus_pairs))))*(-1);
    func<> inter_1;
    inter_1.allocate_mem();
    inter_1=min(Im_Wij.get_lb().in(bus_pairs)/sqrt(Wii.get_ub().from(bus_pairs)*Wii.get_ub().to(bus_pairs)), Im_Wij.get_lb().in(bus_pairs)/sqrt(Wii.get_lb().from(bus_pairs)*Wii.get_lb().to(bus_pairs)));
    func<> theta_L2;
    theta_L2.allocate_mem();
    theta_L2= asin(inter_1.in(bus_pairs));
    func<> inter_2;
    inter_2.allocate_mem();
    inter_2=max(Im_Wij.get_ub().in(bus_pairs)/sqrt(Wii.get_ub().from(bus_pairs)*Wii.get_ub().to(bus_pairs)), Im_Wij.get_ub().in(bus_pairs)/sqrt(Wii.get_lb().from(bus_pairs)*Wii.get_lb().to(bus_pairs)));
    func<> theta_U2;
    theta_U2.allocate_mem();
    theta_U2= asin(inter_2.in(bus_pairs));
    
    func<> inter_U;
    inter_U.allocate_mem();
    inter_U=min(theta_L1.in(bus_pairs)*(-1.0), theta_U2.in(bus_pairs));
    func<> theta_U;
    theta_U.allocate_mem();
    theta_U=min(inter_U.in(bus_pairs), th_max.in(bus_pairs));
    
    func<> inter_L;
    inter_L.allocate_mem();
    inter_L=max(theta_L1.in(bus_pairs), theta_L2.in(bus_pairs));
    func<> theta_L;
    theta_L.allocate_mem();
    theta_L=max(inter_L.in(bus_pairs), th_min.in(bus_pairs));
    func<>  cos_L;
    cos_L.allocate_mem();
         cos_L= (min(cos(theta_L.in(bus_pairs)), cos(theta_U.in(bus_pairs)))).in(bus_pairs);
    cos_L.print();
    theta_L.print();
    theta_U.print();
   

    
    if(add_original){
        SDP.add(R_Vi.in(nodes),Im_Vi.in(nodes));
        R_Vi.initialize_all(1);
    }
    
    if(add_original){
        Im_Vi.set_lb((grid.ref_bus),0);
        Im_Vi.set_ub((grid.ref_bus),0);
        
        R_Vi.set_lb((grid.ref_bus),v_min(grid.ref_bus).eval());
        R_Vi.set_ub((grid.ref_bus),v_max(grid.ref_bus).eval());
        
       auto ref_bus_pairs_from=grid.get_ref_bus_pairs_from();
       auto ref_bus_pairs_to=grid.get_ref_bus_pairs_to();
        
        for (auto &p: *ref_bus_pairs_from._keys)
        {
            auto ngb = p.substr(0,p.find_first_of(","));
            R_Vi.set_lb(ngb, v_min.eval(ngb)*cos_L.eval(p));
            R_Vi.set_ub(ngb, v_max.eval(ngb));
            ilb=gravity::min(v_min.eval(ngb)*sin(theta_L.eval(p)), v_max.eval(ngb)*sin(theta_L.eval(p)));
            iub=gravity::max(v_min.eval(ngb)*sin(theta_U.eval(p)), v_max.eval(ngb)*sin(theta_U.eval(p)));
            Im_Vi.set_lb(ngb, ilb);
            Im_Vi.set_ub(ngb, iub);
        }
        for (auto &p: *ref_bus_pairs_to._keys)
        {
            auto ngb = p.substr(p.find_first_of(",")+1);
            R_Vi.set_lb(ngb, v_min.eval(ngb)*cos_L.eval(p));
            R_Vi.set_ub(ngb, v_max.eval(ngb));
            ilb=gravity::min(v_min.eval(ngb)*sin(theta_L.eval(p)), v_max.eval(ngb)*sin(theta_L.eval(p)));
            iub=gravity::max(v_min.eval(ngb)*sin(theta_U.eval(p)), v_max.eval(ngb)*sin(theta_U.eval(p)));
            Im_Vi.set_lb(ngb, ilb);
            Im_Vi.set_ub(ngb, iub);
        }

        bool convexify = true;
        var<Cpx> Vi("Vi"), Vj("Vj"), Wij("Wij");
        Vi.real_imag(R_Vi.from(bus_pairs_chord), Im_Vi.from(bus_pairs_chord));
        Vj.real_imag(R_Vi.to(bus_pairs_chord), Im_Vi.to(bus_pairs_chord));
        Wij.real_imag(R_Wij.in(bus_pairs_chord), Im_Wij.in(bus_pairs_chord));
        Constraint<Cpx> Linking_Wij("Linking_Wij");
        Linking_Wij = Wij - Vi*conj(Vj);
        SDP.add(Linking_Wij.in(bus_pairs_chord)==0, convexify);
        Vi.real_imag(R_Vi.in(nodes), Im_Vi.in(nodes));
        Constraint<Cpx> Linking_Wi("Linking_Wi");
        Linking_Wi = Wii - Vi*conj(Vi);
        SDP.add(Linking_Wi.in(nodes)==0, convexify);
        
//        Constraint<Cpx> Vol_limit_LB("Vol_limit_LB");
//        Vol_limit_LB = Vi*conj(Vi);
//        SDP.add(Vol_limit_LB.in(nodes)>=pow(v_min,2), convexify);
//        SDP.print();
        auto Im_L = SDP.get_var<double>("Lift_Im_Vi.in(Nodes)_Im_Vi.in(Nodes).in(Nodes)");
        auto R_L = SDP.get_var<double>("Lift_R_Vi.in(Nodes)_R_Vi.in(Nodes).in(Nodes)");
        Constraint<> Vol_limit_LB("Vol_limit_LB");
        Vol_limit_LB = Im_L + R_L - pow(v_min.in(nodes),2);
        SDP.add(Vol_limit_LB.in(nodes)>=0);
        
        Constraint<Cpx> Rank_type1("RankType1");
        Rank_type1 += Wij*conj(Wij) - Wii.from(bus_pairs_chord)*Wii.to(bus_pairs_chord);
        SDP.add(Rank_type1.in(bus_pairs_chord)>=0, convexify);
        //                SDP.print();
        //        SDP.print_symbolic();
        
        
//
//        var<> RWij_ij("RWij_ij"), RWii_jj("RWii_jj");
//        SDP.add(RWij_ij.in(bus_pairs));
//        SDP.add(RWii_jj.in(bus_pairs));
//
//
//        Constraint<Cpx> Linking_Wij_ij("Linking_Wij_ij");
//        Linking_Wij_ij = RWij_ij -Wij*conj(Wij);
//        SDP.add(Linking_Wij_ij.in(bus_pairs)==0, convexify);
//
//        Constraint<Cpx> Linking_Wii_jj("Linking_Wii_jj");
//        Linking_Wii_jj = RWii_jj -Wii.from(bus_pairs)*Wii.to(bus_pairs);
//        SDP.add(Linking_Wii_jj.in(bus_pairs)==0, convexify);
//
//        Constraint<> Rank1a("Rank1a");
//        Rank1a=RWij_ij-RWii_jj;
//        SDP.add(Rank1a.in(bus_pairs)==0);
        
//         auto ref_bus_pairs_ijkl=grid.get_pairsof_bus_pairs_ijkl();
//
//        
//         for (auto &p: *ref_bus_pairs_ijkl._keys)
//         {
//             auto p1 = p.substr(0,p.find_first_of("|"));
//             auto p2 = p.substr(p.find_first_of("|")+1);
//             auto n1=p1.substr(0, p1.find_first_of(","));
//             auto n2=p1.substr(p1.find_first_of(",")+1);
//             auto n3=p2.substr(0, p2.find_first_of(","));
//             auto n4=p2.substr(p2.find_first_of(",")+1);
//             auto p3=n1+","+n4;
//             auto p4=n2+","+n3;
//
//
//         }
        
        
//        var<Cpx> Wij_kl("Wij_kl"), Wii_jj("Wii_jj");
//
//        Constraint<Cpx> Linking_Wij_ij("Linking_Wij_ij");
//        Linking_Wij_ij = RWij_ij -Wij*conj(Wij);
//        SDP.add(Linking_Wij_ij.in(bus_pairs)==0, convexify);
//
//        Constraint<Cpx> Linking_Wii_jj("Linking_Wii_jj");
//        Linking_Wii_jj = RWii_jj -Wii.from(bus_pairs)*Wii.to(bus_pairs);
//        SDP.add(Linking_Wii_jj.in(bus_pairs)==0, convexify);
//
//        Constraint<> Rank1b("Rank1b");
//        Rank1a=RWij_ij-RWii_jj;
//        SDP.add(Rank1a.in(bus_pairs)==0);
        
    }
 
 
   
        
    
    
    /**  Objective */
    auto obj = (product(c1,Pg) + product(c2,pow(Pg,2)) + sum(c0));
    SDP.min(obj);
    
    
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
    Constraint<> SOC("SOC");
    SOC = pow(R_Wij, 2) + pow(Im_Wij, 2) - Wii.from(bus_pairs_chord)*Wii.to(bus_pairs_chord);
    SDP.add(SOC.in(bus_pairs_chord) <= 0);
    
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
    SDP.add(PAD_UB.in(bus_pairs));
    
    Constraint<> PAD_LB("PAD_LB");
    PAD_LB =  Im_Wij.in(bus_pairs);
    PAD_LB >= tan_th_min*R_Wij.in(bus_pairs);
    SDP.add(PAD_LB.in(bus_pairs));
    
    /* Thermal Limit Constraints */
    Constraint<> Thermal_Limit_from("Thermal_Limit_from");
    Thermal_Limit_from = pow(Pf_from, 2) + pow(Qf_from, 2);
    Thermal_Limit_from <= pow(S_max,2);
    SDP.add(Thermal_Limit_from.in(arcs));
    
    
    Constraint<> Thermal_Limit_to("Thermal_Limit_to");
    Thermal_Limit_to = pow(Pf_to, 2) + pow(Qf_to, 2);
    Thermal_Limit_to <= pow(S_max,2);
    SDP.add(Thermal_Limit_to.in(arcs));
    
    if(current_from){
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
        I_from_L = (pow(Pf_from, 2) + pow(Qf_from, 2))*pow(tr,2)-Wii.get_ub().from(arcs)*lij;
        SDP.add(I_from_L.in(arcs) <= 0);
        
        Constraint<> I_from_U("I_from_U");
        I_from_U = Wii.get_lb().from(arcs)*lij - pow(tr,2)*(max(pow(Pf_from.get_ub(), 2),pow(Pf_from.get_lb(), 2))+max(pow(Qf_from.get_ub(),2),pow(Qf_from.get_lb(),2)));
        SDP.add(I_from_U.in(arcs) <= 0);
        
        Constraint<> I_from_U1("I_from_U1");
        I_from_U1 = Wii.get_lb().from(arcs)*lij - pow(tr,2)*pow(S_max,2);
         SDP.add(I_from_U1.in(arcs) <= 0);
    }
    if(current_to){
        
        param<Cpx> T("T"), Y("Y"), Ych("Ych");
        var<Cpx> L_to("L_to"), W("W");
        T.real_imag(cc.in(arcs), dd.in(arcs));
        Y.real_imag(g.in(arcs), b.in(arcs));
        Ych.set_imag(ch_half.in(arcs));
        
        
        L_to.set_real(lji.in(arcs));
        W.real_imag(R_Wij.in(arcs), Im_Wij.in(arcs));
        
        Constraint<Cpx> I_to("I_to");
        I_to=pow(tr,2)*(Y+Ych)*(conj(Y)+conj(Ych))*Wii.to(arcs)-T*conj(Y)*(Y+Ych)*conj(W)-conj(T)*Y*(conj(Y)+conj(Ych))*W+Y*conj(Y)*Wii.from(arcs);
        SDP.add_real(I_to.in(arcs)==pow(tr,2)*L_to.in(arcs));
        
        
        //        func<> lji=lij.in(arcs)-(pow(ch_half.in(arcs), 2)+2*grid.b.in(arcs)*ch_half.in(arcs))/pow(tr.in(arcs),2)*Wii.from(arcs)+(pow(ch_half.in(arcs), 2)+2*grid.b.in(arcs)*ch_half.in(arcs))*Wii.to(arcs)
        //        +(4*ch_half.in(arcs)*g.in(arcs))/pow(tr.in(arcs),2)*(dd.in(arcs)*R_Wij.in(arcs)-cc.in(arcs)*Im_Wij.in(arcs));
        
        Constraint<> I_to_L("I_to_L");
        I_to_L = (pow(Pf_to, 2) + pow(Qf_to, 2))-Wii.get_ub().to(arcs)*lji;
        SDP.add(I_to_L.in(arcs) <= 0);
        
        
        Constraint<> I_to_U("I_to_U");
        I_to_U = Wii.get_lb().to(arcs)*lji - (max(pow(Pf_to.get_ub(), 2),pow(Pf_to.get_lb(), 2))+max(pow(Qf_to.get_ub(),2),pow(Qf_to.get_lb(),2)));
        SDP.add(I_to_U.in(arcs) <= 0);
        
        Constraint<> I_to_U1("I_to_U1");
        I_to_U1 = Wii.get_lb().to(arcs)*lji - pow(S_max,2);
        SDP.add(I_to_U1.in(arcs) <= 0);
        
        //        func<> lji=lij.in(arcs)-(pow(ch_half.in(arcs), 2)+2*grid.b.in(arcs)*ch_half.in(arcs))/pow(tr.in(arcs),2)*Wii.from(arcs)+(pow(ch_half.in(arcs), 2)+2*grid.b.in(arcs)*ch_half.in(arcs))*Wii.to(arcs)
        //        +(4*ch_half.in(arcs)*g.in(arcs))/pow(tr.in(arcs),2)*(dd.in(arcs)*R_Wij.in(arcs)-cc.in(arcs)*Im_Wij.in(arcs));
        
        
    }
    
    
    

    
    
    Constraint<> Im_U("Im_U");
    Im_U=pow(Im_Wij, 2) - max(pow(sin(theta_L.in(bus_pairs)),2),pow(sin(theta_U.in(bus_pairs)),2))*Wii.from(bus_pairs)*Wii.to(bus_pairs);
    SDP.add(Im_U.in(bus_pairs)<=0);
    
    Constraint<> Im_L("Im_L");
    Im_L=min(sin(theta_L.in(bus_pairs))*sqrt(Wii.get_lb().from(bus_pairs)*Wii.get_lb().to(bus_pairs)),sin(theta_L.in(bus_pairs))*
             sqrt(Wii.get_ub().from(bus_pairs)*Wii.get_ub().to(bus_pairs))) -Im_Wij;
    SDP.add(Im_L.in(bus_pairs)<=0);
    
    
    Constraint<> R_L("R_L");
    R_L=min(cos(theta_L.in(bus_pairs)),cos(theta_U.in(bus_pairs)))*sqrt(Wii.get_lb().from(bus_pairs)*
                                                                        Wii.get_lb().to(bus_pairs))-R_Wij;
    SDP.add(R_L.in(bus_pairs)<=0);

    
    if(loss){
        func<> cosl=min(cos(theta_L.in(arcs)-as.in(arcs)),cos(theta_U.in(arcs)-as.in(arcs)));
    
        Constraint<> p_U("p_U");
        p_U=(Pf_from+Pf_to)*pow(tr,2)-g*max(pow(sqrt(Wii.get_ub().from(arcs))-tr*sqrt(Wii.get_lb().to(arcs)), 2),pow(sqrt(Wii.get_lb().from(arcs))-tr*sqrt(Wii.get_ub().to(arcs)), 2))-g*(1-cosl)*(Wii.from(arcs)+pow(tr,2)*Wii.to(arcs));
        SDP.add(p_U.in(arcs)<=0);
        
        Constraint<> p_L("p_L");
        p_L=(Pf_from+Pf_to)*pow(tr,2)-g*cosl*min(pow(sqrt(Wii.get_ub().from(arcs))-tr*sqrt(Wii.get_lb().to(arcs)), 2),pow(sqrt(Wii.get_lb().from(arcs))-tr*sqrt(Wii.get_ub().to(arcs)), 2));
        SDP.add(p_L.in(arcs)>=0);
        
        Constraint<> q_U("q_U");
        q_U=(Qf_from+Qf_to)*pow(tr,2)+b*max(pow(sqrt(Wii.get_ub().from(arcs))-tr*sqrt(Wii.get_lb().to(arcs)), 2),pow(sqrt(Wii.get_lb().from(arcs))-tr*sqrt(Wii.get_ub().to(arcs)), 2))+(b*(1-cosl)+ch_half)*(Wii.from(arcs)+pow(tr,2)*Wii.to(arcs));
        SDP.add(q_U.in(arcs)<=0);
        
        Constraint<> q_L("q_L");
        q_L=(Qf_from+Qf_to)*pow(tr,2)+b*cosl*min(pow(sqrt(Wii.get_ub().from(arcs))-tr*sqrt(Wii.get_lb().to(arcs)), 2),pow(sqrt(Wii.get_lb().from(arcs))-tr*sqrt(Wii.get_ub().to(arcs)), 2))+ch_half*(Wii.from(arcs)+pow(tr,2)*Wii.to(arcs));
        SDP.add(q_L.in(arcs)>=0);
        
        
        
    }
    loss_bounds=false;
    if(false){
        func<> cosl=min(cos(theta_L.in(arcs)-as.in(arcs)),cos(theta_U.in(arcs)-as.in(arcs)));
        b.print();
        
        Constraint<> p_U1("p_U1");
        p_U1=Pf_from+Pf_to-g*(Wii.from(arcs)*1.0/tr*1.0/tr+Wii.to(arcs))+2*g/tr*sqrt(Wii.get_lb().to(arcs))*sqrt(Wii.get_lb().from(arcs))*cosl;
        SDP.add(p_U1.in(arcs)<=0);
        
        Constraint<> p_L1("p_L1");
        p_L1=Pf_from+Pf_to-g*(Wii.from(arcs)*1.0/tr*1.0/tr+Wii.to(arcs))+2*g/tr*sqrt(Wii.get_ub().to(arcs))*sqrt(Wii.get_ub().from(arcs));
        SDP.add(p_L1.in(arcs)>=0);
        
        Constraint<> q_U1("q_U1");
        q_U1=Qf_from+Qf_to+(b+ch_half)*(Wii.from(arcs)*1.0/tr*1.0/tr+Wii.to(arcs))-2*b/tr*sqrt(Wii.get_lb().to(arcs))*sqrt(Wii.get_lb().from(arcs))*cosl;
        SDP.add(q_U1.in(arcs)<=0);
        
        Constraint<> q_L1("q_L1");
        q_L1=Qf_from+Qf_to+(b+ch_half)*(Wii.from(arcs)*1.0/tr*1.0/tr+Wii.to(arcs))-2*b/tr*sqrt(Wii.get_ub().to(arcs))*sqrt(Wii.get_ub().from(arcs));
        SDP.add(q_L1.in(arcs)>=0);
        
        
    }
//    bool loss_bounds1=false;
//    if(false){
//        
//    
//        SDP.add(sWf_sWt.in(arcs), Wf_Wt.in(arcs));
//        
//           bool convexify = true;
//        
////var<Cpx> Wif, Wit, Wift, sW;
//        
////        Wif.set_real(Wii.from(arcs));
////        Wit.set_real(Wii.to(arcs));
////        Wift.set_real(Wf_Wt.in(arcs));
////        sW.set_real(sWf_sWt.in(arcs));
//        
//        
////        Constraint<Cpx> Linking_sqW("Linking_sqW");
////        Linking_sqW = Wift - sW*sW;
////        SDP.add(Linking_sqW.in(arcs)==0, convexify);
////
////        Constraint<Cpx> Linking_WiiWjj("Linking_WiiWjj");
////        Linking_WiiWjj = Wift - Wif*Wit;
////        SDP.add(Linking_WiiWjj.in(arcs)==0, convexify);
//        
//        Constraint<Cpx> Linking_sqW("Linking_sqW");
//        Linking_sqW = Wf_Wt - sWf_sWt.in(arcs)*sWf_sWt.in(arcs);
//        SDP.add(Linking_sqW.in(arcs)==0, convexify);
//        
//        Constraint<Cpx> Linking_WiiWjj("Linking_WiiWjj");
//        Linking_WiiWjj = Wf_Wt - Wii.from(arcs)*Wii.to(arcs);
//        SDP.add(Linking_WiiWjj.in(arcs)==0, convexify);
//        
//        func<> cosl=min(cos(theta_L.in(arcs)-as.in(arcs)),cos(theta_U.in(arcs)-as.in(arcs)));
//        b.print();
//        
//        Constraint<> p_U2("p_U2");
//        p_U2=Pf_from+Pf_to-g*(Wii.from(arcs)*1.0/tr*1.0/tr+Wii.to(arcs))+2*g/tr*sWf_sWt*cosl;
//        SDP.add(p_U2.in(arcs)<=0);
//        
//        Constraint<> p_L2("p_L2");
//        p_L2=Pf_from+Pf_to-g*(Wii.from(arcs)*1.0/tr*1.0/tr+Wii.to(arcs))+2*g/tr*sWf_sWt;
//        SDP.add(p_L2.in(arcs)>=0);
//        
//        Constraint<> q_U2("q_U2");
//        q_U2=Qf_from+Qf_to+(b+ch_half)*(Wii.from(arcs)*1.0/tr*1.0/tr+Wii.to(arcs))-2*b/tr*sWf_sWt*cosl;
//        SDP.add(q_U2.in(arcs)<=0);
//
//        Constraint<> q_L2("q_L2");
//        q_L2=Qf_from+Qf_to+(b+ch_half)*(Wii.from(arcs)*1.0/tr*1.0/tr+Wii.to(arcs))-2*b/tr*sWf_sWt;
//        SDP.add(q_L2.in(arcs)>=0);
//        
//        
//    }
//    bool loss_bounds2=false;
//    if(false){
//        
//        
//        SDP.add(sWf_sWt.in(arcs), Wf_Wt.in(arcs));
//        
//        bool convexify = true;
//        
//        Constraint<Cpx> Linking_sqW("Linking_sqW");
//        Linking_sqW = Wf_Wt - sWf_sWt.in(arcs)*sWf_sWt.in(arcs);
//        SDP.add(Linking_sqW.in(arcs)==0, convexify);
//        
//        Constraint<Cpx> Linking_WiiWjj("Linking_WiiWjj");
//        Linking_WiiWjj = Wf_Wt - Wii.from(arcs)*Wii.to(arcs);
//        SDP.add(Linking_WiiWjj.in(arcs)==0, convexify);
//        
//        func<> cosl=min(cos(theta_L.in(arcs)-as.in(arcs)),cos(theta_U.in(arcs)-as.in(arcs)));
//        b.print();
//        
//        Constraint<> p_U2("p_U2");
//        p_U2=Pf_from+Pf_to-g*(Wii.from(arcs)*1.0/tr*1.0/tr+Wii.to(arcs))+2*g/tr*sWf_sWt*cosl;
//        SDP.add(p_U2.in(arcs)<=0);
//        
//        Constraint<> p_L2("p_L2");
//        p_L2=Pf_from+Pf_to-g*(Wii.from(arcs)*1.0/tr*1.0/tr+Wii.to(arcs))+2*g/tr*sWf_sWt;
//        SDP.add(p_L2.in(arcs)>=0);
//        
//        Constraint<> q_U2("q_U2");
//        q_U2=Qf_from+Qf_to+(b+ch_half)*(Wii.from(arcs)*1.0/tr*1.0/tr+Wii.to(arcs))-2*b/tr*sWf_sWt*cosl;
//        SDP.add(q_U2.in(arcs)<=0);
//        
//        Constraint<> q_L2("q_L2");
//        q_L2=Qf_from+Qf_to+(b+ch_half)*(Wii.from(arcs)*1.0/tr*1.0/tr+Wii.to(arcs))-2*b/tr*sWf_sWt;
//        SDP.add(q_L2.in(arcs)>=0);
//        
//        
//    }
    
    if (llnc)
    {
        
        func<> phi=(theta_U.in(bus_pairs)+theta_L.in(bus_pairs))/2.0;
        DebugOn("phi");
        phi.print();
        func<> del=(theta_U.in(bus_pairs)-theta_L.in(bus_pairs))/2.0;
        DebugOn("del");
        del.print();
        
        Constraint<> LNC1("LNC1");
        LNC1 += (sqrt(Wii.get_lb().from(bus_pairs))+sqrt(Wii.get_ub().from(bus_pairs)))*(sqrt(Wii.get_lb().to(bus_pairs))+sqrt(Wii.get_ub().to(bus_pairs)))*(Im_Wij.in(bus_pairs)*sin(phi.in(bus_pairs)) + R_Wij.in(bus_pairs)*cos(phi.in(bus_pairs)));
        LNC1 -= sqrt(Wii.get_ub().to(bus_pairs))*cos(del.in(bus_pairs))*(sqrt(Wii.get_lb().to(bus_pairs))+sqrt(Wii.get_ub().to(bus_pairs)))*Wii.from(bus_pairs);
        LNC1 -= sqrt(Wii.get_ub().from(bus_pairs))*cos(del.in(bus_pairs))*(sqrt(Wii.get_lb().from(bus_pairs))+sqrt(Wii.get_ub().from(bus_pairs)))*Wii.to(bus_pairs);
        LNC1-=sqrt(Wii.get_ub().from(bus_pairs))*sqrt(Wii.get_ub().to(bus_pairs))*cos(del)*(sqrt(Wii.get_lb().from(bus_pairs))*
                                                                                            sqrt(Wii.get_lb().to(bus_pairs)) - sqrt(Wii.get_ub().from(bus_pairs))*sqrt(Wii.get_ub().to(bus_pairs)));
      //  SDP.add(LNC1.in(bus_pairs) >= 0);
        Constraint<> LNC2("LNC2");
        LNC2 += (sqrt(Wii.get_lb().from(bus_pairs))+sqrt(Wii.get_ub().from(bus_pairs)))*(sqrt(Wii.get_lb().to(bus_pairs))+sqrt(Wii.get_ub().to(bus_pairs)))*(sin(phi.in(bus_pairs))*Im_Wij.in(bus_pairs) + cos(phi.in(bus_pairs))*R_Wij.in(bus_pairs));
        LNC2 -= sqrt(Wii.get_lb().to(bus_pairs))*cos(del.in(bus_pairs))*(sqrt(Wii.get_lb().to(bus_pairs))+sqrt(Wii.get_ub().to(bus_pairs)))*Wii.from(bus_pairs);
        LNC2 -= sqrt(Wii.get_lb().from(bus_pairs))*cos(del.in(bus_pairs))*(sqrt(Wii.get_lb().from(bus_pairs))+sqrt(Wii.get_ub().from(bus_pairs)))*Wii.to(bus_pairs);
        LNC2+=sqrt(Wii.get_lb().from(bus_pairs))*sqrt(Wii.get_lb().to(bus_pairs))*cos(del.in(bus_pairs))*(sqrt(Wii.get_lb().from(bus_pairs))*
                                                                                                          sqrt(Wii.get_lb().to(bus_pairs))-sqrt(Wii.get_ub().from(bus_pairs))*sqrt(Wii.get_ub().to(bus_pairs)));
        //SDP.add(LNC2.in(bus_pairs) >= 0);
        //
        //                        Constraint<> LNC2("LNC2");
        //                        LNC2 += (sqrt(Wii.get_lb().from(bus_pairs))+sqrt(Wii.get_ub().from(bus_pairs)))*(sqrt(Wii.get_lb().to(bus_pairs))+sqrt(Wii.get_ub().to(bus_pairs)))*(grid.sphi*Im_Wij + grid.cphi*R_Wij);
        //                        LNC2 -= sqrt(Wii.get_lb().to(bus_pairs))*grid.cos_d*(sqrt(Wii.get_lb().to(bus_pairs))+sqrt(Wii.get_ub().to(bus_pairs)))*Wii.from(bus_pairs);
        //                        LNC2 -= sqrt(Wii.get_lb().from(bus_pairs))*grid.cos_d*(sqrt(Wii.get_lb().from(bus_pairs))+sqrt(Wii.get_ub().from(bus_pairs)))*Wii.to(bus_pairs);
        //                        LNC2 += sqrt(Wii.get_lb().from(bus_pairs))*sqrt(Wii.get_lb().to(bus_pairs))*grid.cos_d*(sqrt(Wii.get_lb().from(bus_pairs))*sqrt(Wii.get_lb().to(bus_pairs)) - sqrt(Wii.get_ub().from(bus_pairs))*sqrt(Wii.get_ub().to(bus_pairs)));
        //                        SDP.add(LNC2.in(bus_pairs) >= 0);
        //
        //                Constraint<> LNC1("LNC1");
        //                LNC1 += (sqrt(Wii.get_lb().from(bus_pairs))+sqrt(Wii.get_ub().from(bus_pairs)))*(sqrt(Wii.get_lb().to(bus_pairs))+sqrt(Wii.get_ub().to(bus_pairs)))*(grid.sphi*Im_Wij + grid.cphi*R_Wij);
        //                LNC1 -= sqrt(Wii.get_ub().to(bus_pairs))*grid.cos_d*(sqrt(Wii.get_lb().to(bus_pairs))+sqrt(Wii.get_ub().to(bus_pairs)))*Wii.from(bus_pairs);
        //                LNC1 -= sqrt(Wii.get_ub().from(bus_pairs))*grid.cos_d*(sqrt(Wii.get_lb().from(bus_pairs))+sqrt(Wii.get_ub().from(bus_pairs)))*Wii.to(bus_pairs);
        //                LNC1 -=sqrt(Wii.get_ub().from(bus_pairs))*sqrt(Wii.get_ub().to(bus_pairs))*grid.cos_d*(sqrt(Wii.get_lb().from(bus_pairs))*
        //                    sqrt(Wii.get_lb().to(bus_pairs)) - sqrt(Wii.get_ub().from(bus_pairs))*sqrt(Wii.get_ub().to(bus_pairs)));
        //                SDP.add(LNC1.in(bus_pairs) >= 0);
        //
        
    }
    
    
    total_time_start = get_wall_time();
    /* Solver selection */
    solver<> SDPOPF(SDP,solv_type);
    double solver_time_start = get_wall_time();
    
    //SDP.print();
    SDPOPF.run(output = 5, tol = 1e-6);
    SDP.print();
    
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
    
    
    //
    //        string result_name="/Users/smitha/Desktop/Results/SDPCUTS.txt";
    //        ofstream fout(result_name.c_str(), ios_base::app);
    //    fout<<grid._name<<"\t"<<std::fixed<<std::setprecision(5)<<gap<<"\t"<<std::setprecision(5)<<upper_bound<<"\t"<<std::setprecision(5)<<SDP.get_obj_val()<<"\t"<<std::setprecision(5)<<solve_time<<endl;
    
        SDP.print_solution();
    
    
    
    return 0;
    
}

