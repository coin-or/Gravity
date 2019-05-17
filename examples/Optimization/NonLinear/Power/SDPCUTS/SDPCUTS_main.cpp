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
    bool loss_from = true, llnc=true, loss_to=true;
    size_t num_bags = 0;
    string num_bags_s = "100";
    string solver_str = "ipopt";
    string sdp_cuts_s = "yes";
    string loss_from_s = "yes";
    string orig_s = "yes";
    string loss_to_s="yes";
    string lazy_s = "yes";
    bool lazy_bool = true;
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
    opt.add_option("lf", "losses_from", "add from loss constraints", loss_from_s);
    opt.add_option("o", "original", "add original variables and linking constraints", orig_s);
    opt.add_option("lt", "losses_to", "add to loss constraints", loss_to_s);
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
    
    loss_from_s = opt["lf"];
    if (loss_from_s.compare("no")==0) {
        loss_from = false;
    }
    else {
        loss_from = true;
    }
    bool add_original = true;
    
    orig_s = opt["o"];
    if (orig_s.compare("no")==0) {
        add_original = false;
    }
    else {
        add_original = true;
    }

    
    
    loss_to_s = opt["lt"];
    if (loss_to_s.compare("no")==0) {
        loss_to = false;
    }
    else {
        loss_to = true;
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
    auto v_max = grid.v_max.in(nodes);
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
    if(loss_from){
        SDP.add(lij.in(arcs));
    }
    var<> lji("lji", pos_);
    if(loss_to){
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

   

    if(add_original){
        SDP.add(R_Vi.in(nodes),Im_Vi.in(nodes));
        R_Vi.initialize_all(1);
    }
    
    if(add_original){

//        Constraint<> Ref_Bus("Ref_Bus");
//        Ref_Bus = Im_Vi(grid.ref_bus);
//        SDP.add(Ref_Bus == 0);
        Im_Vi.set_lb((grid.ref_bus),0);
        Im_Vi.set_ub((grid.ref_bus),0);
        
        bool convexify = true;
        var<Cpx> Vi("Vi"), Vj("Vj"), Wij("Wij"), Wi("Wi");
        Vi.real_imag(R_Vi.from(bus_pairs_chord), Im_Vi.from(bus_pairs_chord));
        Vj.real_imag(R_Vi.to(bus_pairs_chord), Im_Vi.to(bus_pairs_chord));
        Wij.real_imag(R_Wij.in(bus_pairs_chord), Im_Wij.in(bus_pairs_chord));
        Constraint<Cpx> Linking_Wij("Linking_Wij");
        Linking_Wij = Wij - Vi*conj(Vj);
        SDP.add(Linking_Wij.in(bus_pairs_chord)==0, convexify);
        Wi.set_real(Wii.in(nodes));
        Vi.real_imag(R_Vi.in(nodes), Im_Vi.in(nodes));
        Constraint<Cpx> Linking_Wi("Linking_Wi");
        Linking_Wi = Wii - Vi*conj(Vi);
        SDP.add(Linking_Wi.in(nodes)==0, convexify);
//                SDP.print();
//        SDP.print_symbolic();
    }
    
    /**  Objective */
    auto obj = product(c1,Pg) + product(c2,pow(Pg,2)) + sum(c0);
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
        I_from_U = w_min.from(arcs)*lij - pow(tr,2)*(max(pow(Pf_from.get_ub(), 2),pow(Pf_from.get_lb(), 2))+max(pow(Qf_from.get_ub(),2),pow(Qf_from.get_lb(),2)));
//        SDP.add_lazy(I_from_U.in(arcs) <= 0);
        
        Constraint<> I_from_U1("I_from_U1");
        I_from_U1 = w_min.from(arcs)*lij - pow(tr,2)*pow(S_max,2);
        SDP.add(I_from_U1.in(arcs) <= 0);
    
    }
    
    if(loss_to){
        
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
            
            
            Constraint<> I_to_L("I_to_L");
            I_to_L = (pow(Pf_to, 2) + pow(Qf_to, 2))-w_max.to(arcs)*lji;
            SDP.add_lazy(I_to_L.in(arcs) <= 0);
            
        
            Constraint<> I_to_U("I_to_U");
            I_to_U = w_min.to(arcs)*lji - (max(pow(Pf_to.get_ub(), 2),pow(Pf_to.get_lb(), 2))+max(pow(Qf_to.get_ub(),2),pow(Qf_to.get_lb(),2)));
//            SDP.add_lazy(I_to_U.in(arcs) <= 0);
        
            Constraint<> I_to_U1("I_to_U1");
            I_to_U1 = w_min.to(arcs)*lji - pow(S_max,2);
            SDP.add_lazy(I_to_U1.in(arcs) <= 0);
            
        }
      
    
    
    
    if (llnc)
    {
    /* Lifted Nonlinear Cuts */
    Constraint<> LNC1("LNC1");
    LNC1 += (grid.v_min.from(bus_pairs)+grid.v_max.from(bus_pairs))*(grid.v_min.to(bus_pairs)+grid.v_max.to(bus_pairs))*(grid.sphi*Im_Wij + grid.cphi*R_Wij);
    LNC1 -= grid.v_max.to(bus_pairs)*grid.cos_d*(grid.v_min.to(bus_pairs)+grid.v_max.to(bus_pairs))*Wii.from(bus_pairs);
    LNC1 -= grid.v_max.from(bus_pairs)*grid.cos_d*(grid.v_min.from(bus_pairs)+grid.v_max.from(bus_pairs))*Wii.to(bus_pairs);
    LNC1 -= grid.v_max.from(bus_pairs)*grid.v_max.to(bus_pairs)*grid.cos_d*(grid.v_min.from(bus_pairs)*grid.v_min.to(bus_pairs) - grid.v_max.from(bus_pairs)*grid.v_max.to(bus_pairs));
    SDP.add(LNC1.in(bus_pairs) >= 0);
    
    Constraint<> LNC2("LNC2");
    LNC2 += (grid.v_min.from(bus_pairs)+grid.v_max.from(bus_pairs))*(grid.v_min.to(bus_pairs)+grid.v_max.to(bus_pairs))*(grid.sphi*Im_Wij + grid.cphi*R_Wij);
    LNC2 -= grid.v_min.to(bus_pairs)*grid.cos_d*(grid.v_min.to(bus_pairs)+grid.v_max.to(bus_pairs))*Wii.from(bus_pairs);
    LNC2 -= grid.v_min.from(bus_pairs)*grid.cos_d*(grid.v_min.from(bus_pairs)+grid.v_max.from(bus_pairs))*Wii.to(bus_pairs);
    LNC2 += grid.v_min.from(bus_pairs)*grid.v_min.to(bus_pairs)*grid.cos_d*(grid.v_min.from(bus_pairs)*grid.v_min.to(bus_pairs) - grid.v_max.from(bus_pairs)*grid.v_max.to(bus_pairs));
    SDP.add(LNC2.in(bus_pairs) >= 0);
    }
    
    
    total_time_start = get_wall_time();
    /* Solver selection */
    solver<> SDPOPF(SDP,solv_type);
    double solver_time_start = get_wall_time();



    //SDPOPF.run(output = 5, tol = 1e-6);

    SDPOPF.run(output = 5, tol = 1e-6, "ma97");


    SDP.print();

    SDPOPF.run(output = 5, tol = 1e-6, "ma97");

//    SDP.print();
//    SDP.print_symbolic();

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

//    SDP.print_solution();


    
    return 0;

}

