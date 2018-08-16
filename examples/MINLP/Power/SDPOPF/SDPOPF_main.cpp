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
#include "Bag.h"

using namespace std;
using namespace gravity;

double time_in_create_constrs;

void nearest_point(int i, int j, vector<param<>>& w_hat, const vector<Bag>& bags){
    for (unsigned index = i; index <j; index++) {
        auto b = bags[index];
        w_hat[index] = b.nfp();
    }
}


//int main_new (int argc, char * argv[]) {
//    int output = 0;
//    bool relax = false;
//    string solver_str = "ipopt";
//    SolverType solv_type = ipopt;
//    double tol = 1e-6;
//    string mehrotra = "no";
//    string fname = "../data_sets/Power/nesta_case5_pjm.m";
//
//    // create a OptionParser with options
//    op::OptionParser opt;
//    opt.add_option("h", "help",
//                   "shows option help"); // no default value means boolean options, which default value is false
//    opt.add_option("f", "file", "Input file name", fname);
//    opt.add_option("s", "solver", "Solvers: ipopt/cplex/gurobi, default = ipopt", solver_str);
//    /* parse the options and verify that all went well. If not, errors and help will be shown */
//    bool correct_parsing = opt.parse_options(argc, argv);
//
//    if (!correct_parsing) {
//        return EXIT_FAILURE;
//    }
//
//    fname = opt["f"];
//    bool has_help = op::str2bool(opt["h"]);
//    if (has_help) {
//        opt.show_help();
//        exit(0);
//    }
//    solver_str = opt["s"];
//    if (solver_str.compare("gurobi")==0) {
//        solv_type = gurobi;
//    }
//    else if(solver_str.compare("cplex")==0) {
//        solv_type = cplex;
//    }
//
//    double total_time_start = get_wall_time();
//    PowerNet grid;
//    grid.readgrid(fname.c_str());
//
//    grid.get_tree_decomp_bags();
//    grid.update_net();
//
//    /* Grid Parameters */
//    auto bus_pairs = grid.get_bus_pairs();
//    auto bus_pairs_chord = grid.get_bus_pairs_chord();
//    auto nb_lines = grid.get_nb_active_arcs();
//    auto nb_buses = grid.get_nb_active_nodes();
//    DebugOff("nb gens = " << nb_gen << endl);
//    DebugOff("nb lines = " << nb_lines << endl);
//    DebugOff("nb buses = " << nb_buses << endl);
//    DebugOff("nb bus_pairs = " << nb_bus_pairs_chord << endl);
//
//
//    double upper_bound = grid.solve_acopf();
//
//
//    /** Build model */
//    Model SDP("SDP Model");
//
//    /** Variables */
//    /* power generation variables */
//    var<Real> Pg("Pg", grid.pg_min, grid.pg_max);
//    var<Real> Qg ("Qg", grid.qg_min, grid.qg_max);
//    SDP.add_var(Pg.in(grid.gens));
//    SDP.add_var(Qg.in(grid.gens));
//
//
//    /* power flow variables */
//    var<Real> Pf_from("Pf_from", grid.S_max);
//    var<Real> Qf_from("Qf_from", grid.S_max);
//    var<Real> Pf_to("Pf_to", grid.S_max);
//    var<Real> Qf_to("Qf_to", grid.S_max);
//    SDP.add_var(Pf_from.in(grid.arcs));
//    SDP.add_var(Qf_from.in(grid.arcs));
//    SDP.add_var(Pf_to.in(grid.arcs));
//    SDP.add_var(Qf_to.in(grid.arcs));
//
//    /* Real part of Wij = ViVj */
//    var<Real>  R_Wij("R_Wij", grid.wr_min, grid.wr_max);
//    /* Imaginary part of Wij = ViVj */
//    var<Real>  Im_Wij("Im_Wij", grid.wi_min, grid.wi_max);
//    /* Magnitude of Wii = Vi^2 */
//    var<Real>  Wii("Wii", grid.w_min, grid.w_max);
//    SDP.add_var(Wii.in(grid.nodes));
//    SDP.add_var(R_Wij.in(bus_pairs_chord));
//    SDP.add_var(Im_Wij.in(bus_pairs_chord));
//
//    /* Initialize variables */
//    R_Wij.initialize_all(1.0);
//    Wii.initialize_all(1.001);
//
//    /**  Objective */
//    auto obj = product(grid.c1, Pg) + product(grid.c2, power(Pg,2)) + sum(grid.c0);
//    SDP.min(obj.in(grid.gens));
//
//    /** Constraints */
//
//    DebugOn("Adding Cholesky cuts\n");
//    auto R_Wij_ = R_Wij.pairs_in_directed(grid, grid._bags, 3);
//    auto Im_Wij_ = Im_Wij.pairs_in_directed(grid, grid._bags, 3);
//    auto Wii_ = Wii.in(grid._bags, 3);
//    auto I_sgn = signs(grid, grid._bags);
//    Constraint Wii_LL("Wii_LL");
//    Wii_LL = Wii_[0];
//    SDP.add(Wii_LL == 0);
//
//
//    /* Second-order cone constraints */
//    Constraint SOC("SOC");
//    SOC = power(R_Wij, 2) + power(Im_Wij, 2) - Wii.from()*Wii.to();
//    SDP.add_constraint(SOC.in(bus_pairs) <= 0);
//
//    /* Flow conservation */
//    Constraint KCL_P("KCL_P");
//    KCL_P  = sum(Pf_from.out_arcs()) + sum(Pf_to.in_arcs()) + grid.pl - sum(Pg.in_gens()) + grid.gs*Wii;
//    SDP.add_constraint(KCL_P.in(grid.nodes) == 0);
//
//    Constraint KCL_Q("KCL_Q");
//    KCL_Q  = sum(Qf_from.out_arcs()) + sum(Qf_to.in_arcs()) + grid.ql - sum(Qg.in_gens()) - grid.bs*Wii;
//    SDP.add_constraint(KCL_Q.in(grid.nodes) == 0);
//
//    /* AC Power Flow */
//    Constraint Flow_P_From("Flow_P_From");
//    Flow_P_From = Pf_from - (grid.g_ff*Wii.from() + grid.g_ft*R_Wij.in_pairs() + grid.b_ft*Im_Wij.in_pairs());
//    SDP.add_constraint(Flow_P_From.in(grid.arcs) == 0);
//
//    Constraint Flow_P_To("Flow_P_To");
//    Flow_P_To = Pf_to - (grid.g_tt*Wii.to() + grid.g_tf*R_Wij.in_pairs() - grid.b_tf*Im_Wij.in_pairs());
//    SDP.add_constraint(Flow_P_To.in(grid.arcs) == 0);
//
//    Constraint Flow_Q_From("Flow_Q_From");
//    Flow_Q_From = Qf_from - (grid.g_ft*Im_Wij.in_pairs() - grid.b_ff*Wii.from() - grid.b_ft*R_Wij.in_pairs());
//    SDP.add_constraint(Flow_Q_From.in(grid.arcs) == 0);
//
//    Constraint Flow_Q_To("Flow_Q_To");
//    Flow_Q_To = Qf_to + (grid.b_tt*Wii.to() + grid.b_tf*R_Wij.in_pairs() + grid.g_tf*Im_Wij.in_pairs());
//    SDP.add_constraint(Flow_Q_To.in(grid.arcs) == 0);
//
//    /* Phase Angle Bounds constraints */
//
//    Constraint PAD_UB("PAD_UB");
//    PAD_UB = Im_Wij;
//    PAD_UB <= grid.tan_th_max*R_Wij;
//    SDP.add_lazy(PAD_UB.in(bus_pairs));
//
//    Constraint PAD_LB("PAD_LB");
//    PAD_LB =  Im_Wij;
//    PAD_LB >= grid.tan_th_min*R_Wij;
//    SDP.add_lazy(PAD_LB.in(bus_pairs));
//
//    /* Thermal Limit Constraints */
//    Constraint Thermal_Limit_from("Thermal_Limit_from");
//    Thermal_Limit_from = power(Pf_from, 2) + power(Qf_from, 2);
//    Thermal_Limit_from <= power(grid.S_max,2);
//    SDP.add_lazy(Thermal_Limit_from.in(grid.arcs));
//
//
//    Constraint Thermal_Limit_to("Thermal_Limit_to");
//    Thermal_Limit_to = power(Pf_to, 2) + power(Qf_to, 2);
//    Thermal_Limit_to <= power(grid.S_max,2);
//    SDP.add_lazy(Thermal_Limit_to.in(grid.arcs));
//
//    /* Lifted Nonlinear Cuts */
//    Constraint LNC1("LNC1");
//    LNC1 = (grid.v_min.from()+grid.v_max.from())*(grid.v_min.to()+grid.v_max.to())*(sin(0.5*(grid.th_max+grid.th_min))*Im_Wij + cos(0.5*(grid.th_max+grid.th_min))*R_Wij);
//    LNC1 -= grid.v_max.to()*cos(0.5*(grid.th_max-grid.th_min))*(grid.v_min.to()+grid.v_max.to())*Wii.from();
//    LNC1 -= grid.v_max.from()*cos(0.5*(grid.th_max-grid.th_min))*(grid.v_min.from()+grid.v_max.from())*Wii.to();
//    LNC1 -= grid.v_max.from()*grid.v_max.to()*cos(0.5*(grid.th_max-grid.th_min))*(grid.v_min.from()*grid.v_min.to() - grid.v_max.from()*grid.v_max.to());
//    SDP.add_lazy(LNC1.in(bus_pairs) >= 0);
//
//    Constraint LNC2("LNC2");
//    LNC2 = (grid.v_min.from()+grid.v_max.from())*(grid.v_min.to()+grid.v_max.to())*(sin(0.5*(grid.th_max+grid.th_min))*Im_Wij + cos(0.5*(grid.th_max+grid.th_min))*R_Wij);
//    LNC2 -= grid.v_min.to()*cos(0.5*(grid.th_max-grid.th_min))*(grid.v_min.to()+grid.v_max.to())*Wii.from();
//    LNC2 -= grid.v_min.from()*cos(0.5*(grid.th_max-grid.th_min))*(grid.v_min.from()+grid.v_max.from())*Wii.to();
//    LNC2 += grid.v_min.from()*grid.v_min.to()*cos(0.5*(grid.th_max-grid.th_min))*(grid.v_min.from()*grid.v_min.to() - grid.v_max.from()*grid.v_max.to());
//    SDP.add_lazy(LNC2.in(bus_pairs) >= 0);
//
//    vector<Bag> bags;
//    int n3 = 0;
//    int bagid = 0;
//    for(auto& b: grid._bags){
//        bags.push_back(Bag(bagid,grid,b));
//        if(b.size()==3) n3++;
//        bagid++;
//    }
//
//    DebugOn("\nNum of 3d bags = " << n3);
//
//
//    /* Solver selection */
//    solver SDPOPF(SDP,solv_type);
//    double solver_time_start = get_wall_time();
//    SDPOPF.run(output = 5, relax = false);
//    //    SDPOPF.run(output = 0, relax = false, tol = 1e-6, "ma27", mehrotra = "no");
//
//    for(auto& arc: grid.arcs){
//        ((Line*)arc)->wr = R_Wij(arc->_src->_name+","+arc->_dest->_name).eval();
//        ((Line*)arc)->wi = Im_Wij(arc->_src->_name+","+arc->_dest->_name).eval();
//        if(grid.add_3d_nlin && !arc->_free) arc->_active = true;
//    }
//    for(auto& node: grid.nodes) ((Bus*)node)->w = Wii(node->_name).eval();
//
//
//    //    param<double> w_hat;
//    string namew, namewr, namewi;
//    int numcuts = 0;
//    node_pairs bus_pairs_sdp;
//    //    double prev_opt = 0;
//    //    double fp_tol = 0.0001;
//    DebugOff("\nNumber of bags = " << bags.size());
//
////    vector<param<double>> R_star, I_star, W_star;
//    vector<param<double>> R_diff, I_diff, W_diff;
//    vector<param<double>> R_hat, I_hat, W_hat;
//
//    //    int unchanged = 0;
//    int iter = 0, hdim_cuts = 0, cuts_added = 1;
//
//    //    while(unchanged < 3) {
//    unsigned nr_threads = 1;
//    double gap = 100*(upper_bound - SDP._obj_val)/upper_bound;
//    while(cuts_added > 0 && gap > 1) {
//        DebugOn("Current Gap = " << to_string(gap) << "%."<<endl);
//        cuts_added = 0;
//        vector<param<>> w_hat_vec;
//
//        iter++;
//        //    while((SDP._obj_val - prev_opt)/SDP._obj_val > fp_tol) {
//        //        prev_opt = SDP._obj_val;
//
//        /* Update _is_psd */
//        vector<Bag> violated_bags;
//        for(auto& b: bags){
//            b.update_PSD();
//            if (!b._is_psd) {
//                violated_bags.push_back(b);
//            }
//        }
//        auto nb_bags = violated_bags.size();
//        w_hat_vec.resize(nb_bags);
//        vector<thread> threads;
//        /* Split subproblems into nr_threads parts */
//        vector<int> limits = bounds(nr_threads, nb_bags);
//        /* Launch all threads in parallel */
//        for (int i = 0; i < nr_threads; ++i) {
//            threads.push_back(thread(nearest_point, limits[i], limits[i+1], ref(w_hat_vec), ref(violated_bags)));
//        }
//        /* Join the threads with the main thread */
//        for(auto &t : threads){
//            t.join();
//        }
//        for (unsigned index = 0; index < nb_bags; index++) {
//            auto b = violated_bags[index];
//            DebugOff("\nBag number " << b._id);
//            DebugOff("\nNodes: n1 = " << b._nodes[0]->_name << ", n2 = " << b._nodes[1]->_name << ", n3 = " << b._nodes[2]->_name);
//            auto w_hat(move(w_hat_vec[index]));
////            R_star.resize(numcuts+1);
////            R_star[numcuts] = param<>("R_star"+to_string(numcuts));
////            I_star.resize(numcuts+1);
////            I_star[numcuts] = param<>("I_star"+to_string(numcuts));
////            W_star.resize(numcuts+1);
////            W_star[numcuts] = param<>("W_star"+to_string(numcuts));
//            R_diff.resize(numcuts+1);
//            R_diff[numcuts] = param<>("R_diff"+to_string(numcuts));
//            I_diff.resize(numcuts+1);
//            I_diff[numcuts] = param<>("I_diff"+to_string(numcuts));
//            W_diff.resize(numcuts+1);
//            W_diff[numcuts] = param<>("W_diff"+to_string(numcuts));
//            R_hat.resize(numcuts+1);
//            R_hat[numcuts] = param<>("R_hat"+to_string(numcuts));
//            I_hat.resize(numcuts+1);
//            I_hat[numcuts] = param<>("I_hat"+to_string(numcuts));
//            W_hat.resize(numcuts+1);
//            W_hat[numcuts] = param<>("W_hat"+to_string(numcuts));
//
//            node_pairs b_pairs("node_pairs");
//            Constraint sdpcut("sdpcut_" + to_string(numcuts));
//            Node *ni;
//            Arc *aij;
//            for (int i = 0; i < b._nodes.size(); i++) {
//                for (int j = i; j < b._nodes.size(); j++) {
//                    if (i == j) {
//                        namew = "w(" + b._nodes[i]->_name + ")";
//                        ni = b._nodes[i];
//                        W_diff[numcuts].set_val(ni->_name,((Bus *) ni)->w - w_hat(namew).eval());
//                        W_hat[numcuts].set_val(ni->_name,w_hat(namew).eval());
////                        W_star[numcuts].set_val(ni->_name,((Bus *) ni)->w);
//                    }
//                    else {
//                        aij = grid.get_arc(b._nodes[i]->_name, b._nodes[j]->_name);
//                        if(aij->_imaginary && !aij->_active) {
//                            bus_pairs_sdp._keys.push_back(new index_pair(index_(aij->_src->_name), index_(aij->_dest->_name)));
//                            aij->_active = true;
//                        }
//                        b_pairs._keys.push_back(new index_pair(index_(aij->_src->_name), index_(aij->_dest->_name)));
//
//                        namewr = "wr(" + b._nodes[i]->_name + "," + b._nodes[j]->_name + ")";
//                        namewi = "wi(" + b._nodes[i]->_name + "," + b._nodes[j]->_name + ")";
//
//                        R_diff[numcuts].set_val(aij->_src->_name + "," + aij->_dest->_name,((Line *) aij)->wr - w_hat(namewr).eval());
//                        I_diff[numcuts].set_val(aij->_src->_name + "," + aij->_dest->_name,((Line *) aij)->wi - w_hat(namewi).eval());
//                        R_hat[numcuts].set_val(aij->_src->_name + "," + aij->_dest->_name,w_hat(namewr).eval());
//                        I_hat[numcuts].set_val(aij->_src->_name + "," + aij->_dest->_name,w_hat(namewi).eval());
////                        R_star[numcuts].set_val(aij->_src->_name + "," + aij->_dest->_name,((Line *) aij)->wr);
////                        I_star[numcuts].set_val(aij->_src->_name + "," + aij->_dest->_name,((Line *) aij)->wi);
//                    }
//                }
//            }
//            if(b._nodes.size()>3) hdim_cuts++;
//            DebugOn("W_hat = \n");
//            W_hat[numcuts].print(true);
//            cout << "\n";
//            Constraint lin("lin"+to_string(numcuts));
//            lin = product(R_diff[numcuts].in(b_pairs),(R_Wij.in(b_pairs) - R_hat[numcuts].in(b_pairs)));
//            lin += product(I_diff[numcuts].in(b_pairs),(Im_Wij.in(b_pairs) - I_hat[numcuts].in(b_pairs)));
//            lin += product(W_diff[numcuts].in(b._nodes),(Wii.in(b._nodes) - W_hat[numcuts].in(b._nodes)));
//            //            lin.print_expanded();
//            SDP.add_constraint(lin <= 0);
//
//            cuts_added++; numcuts++;
//        }
//
//        if(!bus_pairs_sdp._keys.empty()) {
//            SDP.add_indices("SOC",bus_pairs_sdp);
//            bus_pairs_sdp.clear();
//        }
//
//
//        for(auto& a: grid.arcs){
//            if(a->_imaginary && !a->_active) a->_free = true;
//        }
//
//        SDPOPF.run(output = 5, relax = false, tol = 1e-6, "ma27", mehrotra = "no");
//        gap = 100*(upper_bound - SDP._obj_val)/upper_bound;
//        DebugOff("\nPrev = " << prev_opt << ", cur = " << SDP._obj_val);
//        DebugOff("\n(opt - prev)/opt = " << (SDP._obj_val - prev_opt)/SDP._obj_val << endl);
//
//        /* Update values for w_star */
//        for(auto& arc: grid.arcs){
//            ((Line*)arc)->wr = R_Wij(arc->_src->_name+","+arc->_dest->_name).eval();
//            ((Line*)arc)->wi = Im_Wij(arc->_src->_name+","+arc->_dest->_name).eval();
//        }
//        for(auto& node: grid.nodes) ((Bus*)node)->w = Wii(node->_name).eval();
//
//
//        cout << "\nNum of iterations = " << iter << ", number of cuts = " << numcuts << endl;
//        //        if((SDP._obj_val - prev_opt)/SDP._obj_val < fp_tol) unchanged++;
//    }
//    string lin = "lin";
//    int num_act_cuts = 0;
//    for (auto &cp: SDP._cons) {
//        if(cp.second->get_name().find(lin) != std::string::npos) {
//            for (unsigned inst = 0; inst < cp.second->_nb_instances; inst++) {
//                cp.second->eval(inst);
//                if(cp.second->is_active(inst)) num_act_cuts++;
//            }
//        }
//    }
//
//    cout << "\nEnd: num of iterations = " << iter << ", number of cuts = " << numcuts << ", number of higher dim cuts = " << hdim_cuts;
//    cout << "\nNumber of active cuts = " << num_act_cuts;
//
//    double solver_time_end = get_wall_time();
//    double total_time_end = get_wall_time();
//    auto solve_time = solver_time_end - solver_time_start;
//    auto total_time = total_time_end - total_time_start;
//    string out = "\nDATA_OPF, " + grid._name + ", " + to_string(nb_buses) + ", " + to_string(nb_lines) +", " + to_string(SDP._obj_val) + ", " + to_string(-numeric_limits<double>::infinity()) + ", " + to_string(solve_time) + ", LocalOptimal, " + to_string(total_time);
//    DebugOn(out <<endl);
//    DebugOn("Final Gap = " << to_string(gap) << "%."<<endl);
//    DebugOn("Upper bound = " << to_string(upper_bound) << "."<<endl);
//    DebugOn("Lower bound = " << to_string(SDP._obj_val) << "."<<endl);
//    DebugOn("\nResults: " << grid._name << " " << to_string(SDP._obj_val) << " " << to_string(total_time)<<endl);
//    return 0;
//}

const bool bag_compare(const Bag b1, const Bag b2) {
    return b1._violation < b2._violation;
}

/* main */
int main (int argc, char * argv[]) {
    int output = 0;
    bool relax = false, sdp_cuts = true;
    size_t num_bags = 0;
    string num_bags_s = "100";
    string solver_str = "ipopt";
    string sdp_cuts_s = "yes";
    SolverType solv_type = ipopt;
    double tol = 1e-6;
    string mehrotra = "no";

    string fname = "../data_sets/Power/nesta_case5_pjm.m";

    // create a OptionParser with options
    op::OptionParser opt;
    opt.add_option("h", "help",
                   "shows option help"); // no default value means boolean options, which default value is false
    opt.add_option("f", "file", "Input file name", fname);
    opt.add_option("s", "solver", "Solvers: ipopt/cplex/gurobi, default = ipopt", solver_str);
    opt.add_option("b", "numbags", "Number of bags per iteration", num_bags_s);
    opt.add_option("c", "sdpcuts", "Generate 3d SDP cuts, default = yes", sdp_cuts_s);
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
        solv_type = Mosek;
    }
    solver_str = opt["c"];
    if (solver_str.compare("no")==0) {
        sdp_cuts = false;
    }
    
    num_bags = atoi(opt["b"].c_str());

    cout << "\nnum bags = " << num_bags;

    double total_time_start = get_wall_time();
    PowerNet grid;
    grid.readgrid(fname.c_str());

    grid.get_tree_decomp_bags();
    grid.update_net();

    // Grid Parameters
    auto bus_pairs = grid.get_bus_pairs();
    auto bus_pairs_chord = grid.get_bus_pairs_chord();
    auto nb_gen = grid.get_nb_active_gens();
    auto nb_lines = grid.get_nb_active_arcs();
    auto nb_buses = grid.get_nb_active_nodes();
    DebugOff("nb gens = " << nb_gen << endl);
    DebugOff("nb lines = " << nb_lines << endl);
    DebugOff("nb buses = " << nb_buses << endl);
    DebugOff("nb bus_pairs = " << nb_bus_pairs_chord << endl);


    double upper_bound = grid.solve_acopf();


    /** Build model */
    Model SDP("SDP Model");

    /** Variables */
    /* power generation variables */
    var<Real> Pg("Pg", grid.pg_min, grid.pg_max);
    var<Real> Qg ("Qg", grid.qg_min, grid.qg_max);
    SDP.add(Pg.in(grid.gens));
    SDP.add(Qg.in(grid.gens));
    /* power flow variables */
    var<Real> Pf_from("Pf_from", grid.S_max);
    var<Real> Qf_from("Qf_from", grid.S_max);
    var<Real> Pf_to("Pf_to", grid.S_max);
    var<Real> Qf_to("Qf_to", grid.S_max);
    SDP.add(Pf_from.in(grid.arcs));
    SDP.add(Qf_from.in(grid.arcs));
    SDP.add(Pf_to.in(grid.arcs));
    SDP.add(Qf_to.in(grid.arcs));

    /* Real part of Wij = ViVj */
    var<Real>  R_Wij("R_Wij", grid.wr_min, grid.wr_max);
    /* Imaginary part of Wij = ViVj */
    var<Real>  Im_Wij("Im_Wij", grid.wi_min, grid.wi_max);
    /* Magnitude of Wii = Vi^2 */
    var<Real>  Wii("Wii", grid.w_min, grid.w_max);
    SDP.add(Wii.in(grid.nodes));

    SDP.add(R_Wij.in(bus_pairs_chord));
    SDP.add(Im_Wij.in(bus_pairs_chord));
//    SDP.add(R_Wij.in(bus_pairs));
//    SDP.add(Im_Wij.in(bus_pairs));
    /* Initialize variables */
    R_Wij.initialize_all(1.0);
    Wii.initialize_all(1.001);

    /**  Objective */
    auto obj = product(grid.c1, Pg) + product(grid.c2, power(Pg,2)) + sum(grid.c0);
//    obj.print_expanded();
    SDP.min(obj.in(grid.gens));

    /** Constraints */

//    if(grid.add_3d_nlin && sdp_cuts) {
    if (true) {
        auto bags_copy = grid._bags;
        auto bag_size = grid._bags.size();
        DebugOn("\nNum of bags = " << bag_size << endl);
//        if (bag_size>1000) {
//            bags_copy.resize(1000);
//            DebugOn("\nOnly adding the first 1000 3d cuts." << endl);
//        }
        DebugOn("Adding 3d determinant polynomial cuts\n");
        auto R_Wij_ = R_Wij.pairs_in_directed(grid, bags_copy, 3);
        auto Im_Wij_ = Im_Wij.pairs_in_directed(grid, bags_copy, 3);
        auto Wii_ = Wii.in(bags_copy, 3);
        
        Constraint SDP3("SDP_3D");
        SDP3 = 2 * R_Wij_[0] * (R_Wij_[1] * R_Wij_[2] + Im_Wij_[1] * Im_Wij_[2]);
        SDP3 += 2 * Im_Wij_[0] * (R_Wij_[1] * Im_Wij_[2] - Im_Wij_[1] * R_Wij_[2]);
        SDP3 -= (power(R_Wij_[0], 2) + power(Im_Wij_[0], 2)) * Wii_[2];
        SDP3 -= (power(R_Wij_[1], 2) + power(Im_Wij_[1], 2)) * Wii_[0];
        SDP3 -= (power(R_Wij_[2], 2) + power(Im_Wij_[2], 2)) * Wii_[1];
        SDP3 += Wii_[0] * Wii_[1] * Wii_[2];
        SDP.add(SDP3 >= 0);
        
        DebugOn("\nsdp nb inst = " << SDP3.get_nb_instances() << endl);
//        SDP3.print_expanded();
    }


    /* Second-order cone constraints */
    Constraint SOC("SOC");
    SOC = power(R_Wij, 2) + power(Im_Wij, 2) - Wii.from()*Wii.to();
    if (sdp_cuts) {
        SDP.add(SOC.in(bus_pairs_chord) <= 0);
    }
    else{
        SDP.add(SOC.in(bus_pairs) <= 0);
    }
    
    
    /* Flow conservation */
    Constraint KCL_P("KCL_P");
    KCL_P  = sum(Pf_from.out_arcs()) + sum(Pf_to.in_arcs()) + grid.pl - sum(Pg.in_gens()) + grid.gs*Wii;
    SDP.add(KCL_P.in(grid.nodes) == 0);
    
    Constraint KCL_Q("KCL_Q");
    KCL_Q  = sum(Qf_from.out_arcs()) + sum(Qf_to.in_arcs()) + grid.ql - sum(Qg.in_gens()) - grid.bs*Wii;
    SDP.add(KCL_Q.in(grid.nodes) == 0);
    
    /* AC Power Flow */
    Constraint Flow_P_From("Flow_P_From");
    Flow_P_From = Pf_from - (grid.g_ff*Wii.from() + grid.g_ft*R_Wij.in_pairs() + grid.b_ft*Im_Wij.in_pairs());
    SDP.add(Flow_P_From.in(grid.arcs) == 0);
    
    Constraint Flow_P_To("Flow_P_To");
    Flow_P_To = Pf_to - (grid.g_tt*Wii.to() + grid.g_tf*R_Wij.in_pairs() - grid.b_tf*Im_Wij.in_pairs());
    SDP.add(Flow_P_To.in(grid.arcs) == 0);
    
    Constraint Flow_Q_From("Flow_Q_From");
    Flow_Q_From = Qf_from - (grid.g_ft*Im_Wij.in_pairs() - grid.b_ff*Wii.from() - grid.b_ft*R_Wij.in_pairs());
    SDP.add(Flow_Q_From.in(grid.arcs) == 0);
    
    Constraint Flow_Q_To("Flow_Q_To");
    Flow_Q_To = Qf_to + (grid.b_tt*Wii.to() + grid.b_tf*R_Wij.in_pairs() + grid.g_tf*Im_Wij.in_pairs());
    SDP.add(Flow_Q_To.in(grid.arcs) == 0);
    

    /* Phase Angle Bounds constraints */

    Constraint PAD_UB("PAD_UB");
    PAD_UB = Im_Wij;
    PAD_UB <= grid.tan_th_max*R_Wij;
//    SDP.add_lazy(PAD_UB.in(bus_pairs));
    SDP.add_constraint(PAD_UB.in(bus_pairs));

    Constraint PAD_LB("PAD_LB");
    PAD_LB =  Im_Wij;
    PAD_LB >= grid.tan_th_min*R_Wij;
//    SDP.add_lazy(PAD_LB.in(bus_pairs));
    SDP.add_constraint(PAD_LB.in(bus_pairs));

    /* Thermal Limit Constraints */
    Constraint Thermal_Limit_from("Thermal_Limit_from");
    Thermal_Limit_from = power(Pf_from, 2) + power(Qf_from, 2);
    Thermal_Limit_from <= power(grid.S_max,2);
    SDP.add(Thermal_Limit_from.in(grid.arcs));
    
    
    Constraint Thermal_Limit_to("Thermal_Limit_to");
    Thermal_Limit_to = power(Pf_to, 2) + power(Qf_to, 2);
    Thermal_Limit_to <= power(grid.S_max,2);
    SDP.add(Thermal_Limit_to.in(grid.arcs));
    

    /* Lifted Nonlinear Cuts */
    Constraint LNC1("LNC1");
    LNC1 = (grid.v_min.from()+grid.v_max.from())*(grid.v_min.to()+grid.v_max.to())*(sin(0.5*(grid.th_max+grid.th_min))*Im_Wij + cos(0.5*(grid.th_max+grid.th_min))*R_Wij);
    LNC1 -= grid.v_max.to()*cos(0.5*(grid.th_max-grid.th_min))*(grid.v_min.to()+grid.v_max.to())*Wii.from();
    LNC1 -= grid.v_max.from()*cos(0.5*(grid.th_max-grid.th_min))*(grid.v_min.from()+grid.v_max.from())*Wii.to();
    LNC1 -= grid.v_max.from()*grid.v_max.to()*cos(0.5*(grid.th_max-grid.th_min))*(grid.v_min.from()*grid.v_min.to() - grid.v_max.from()*grid.v_max.to());
//    SDP.add(LNC1.in(bus_pairs) >= 0);
//    SDP.add_lazy(LNC1.in(bus_pairs) >= 0);

    Constraint LNC2("LNC2");
    LNC2 = (grid.v_min.from()+grid.v_max.from())*(grid.v_min.to()+grid.v_max.to())*(sin(0.5*(grid.th_max+grid.th_min))*Im_Wij + cos(0.5*(grid.th_max+grid.th_min))*R_Wij);
    LNC2 -= grid.v_min.to()*cos(0.5*(grid.th_max-grid.th_min))*(grid.v_min.to()+grid.v_max.to())*Wii.from();
    LNC2 -= grid.v_min.from()*cos(0.5*(grid.th_max-grid.th_min))*(grid.v_min.from()+grid.v_max.from())*Wii.to();
    LNC2 += grid.v_min.from()*grid.v_min.to()*cos(0.5*(grid.th_max-grid.th_min))*(grid.v_min.from()*grid.v_min.to() - grid.v_max.from()*grid.v_max.to());
//    SDP.add(LNC2.in(bus_pairs) >= 0);
//    SDP.add_lazy(LNC2.in(bus_pairs) >= 0);

    vector<Bag> bags;
    int n3 = 0;
    int bagid = 0;
    for(auto& b: grid._bags){
        bags.push_back(Bag(bagid,grid,b));
        if(b.size()==3) n3++;
        bagid++;
    }

    num_bags = bags.size()*num_bags / 100;
    double num_bags_tol = num_bags/2;

    if(bags.size() < 20) {
        num_bags = bags.size();
        num_bags_tol = num_bags;
    }

    DebugOn("\nNum of 3d bags = " << n3 << endl);

    total_time_start = get_wall_time();
    /* Solver selection */
    solver SDPOPF(SDP,solv_type);
    double solver_time_start = get_wall_time();
    SDPOPF.run(output = 5, relax = false);
//    SDP.print_expanded();
//    SDPOPF.run(output = 0, relax = false, tol = 1e-6, "ma27", mehrotra = "no");

    for(auto& arc: grid.arcs){
        ((Line*)arc)->wr = R_Wij(arc->_src->_name+","+arc->_dest->_name).eval();
        ((Line*)arc)->wi = Im_Wij(arc->_src->_name+","+arc->_dest->_name).eval();
        if(grid.add_3d_nlin && !arc->_free) arc->_active = true;
    }
    for(auto& node: grid.nodes) ((Bus*)node)->w = Wii(node->_name).eval();

    string namew, namewr, namewi;
    int numcuts = 0;
    node_pairs bus_pairs_sdp;
    DebugOff("\nNumber of bags = " << bags.size());

//    vector<param<double>> R_star, I_star, W_star;
    vector<param<double>> R_diff, I_diff, W_diff;
    vector<param<double>> R_hat, I_hat, W_hat;

    int iter = 0, hdim_cuts = 0, cuts_added = 1;
    double time_in_all_nfp = 0, time_in_nfp;

    unsigned nr_threads = 5;
    double gap = 100*(upper_bound - SDP._obj_val)/upper_bound;
    while(cuts_added > 0 && gap > 1) {
        DebugOn("Current Gap = " << to_string(gap) << "%."<<endl);
        cuts_added = 0;
        vector<param<>> w_hat_vec;

        iter++;

//        for(auto& b: bags) {
//            if(b._nodes.size()==3) b.add_lines();
//        }

        /* Update _is_psd */
        vector<Bag> violated_bags;
        for(auto& b: bags){
            b.update_PSD();
            if (!b._is_psd) {
                violated_bags.push_back(b);
            }
        }
	if(violated_bags.size() > num_bags) std::partial_sort(violated_bags.begin(), violated_bags.begin()+num_bags, violated_bags.end(), bag_compare);
        else if(violated_bags.size() < num_bags_tol) break;

	//for(auto& bv: violated_bags) {
        //    cout << "\nViolation: " << bv._violation;
        //}

        auto nb_bags = violated_bags.size();
//        size_t nb_bags = min(violated_bags.size(),size_t(num_bags));
        w_hat_vec.resize(nb_bags);
        vector<thread> threads;
        /* Split subproblems into nr_threads parts */
        vector<int> limits = bounds(nr_threads, nb_bags);
        time_in_nfp = get_wall_time();
        /* Launch all threads in parallel */
        for (int i = 0; i < nr_threads; ++i) {
            threads.push_back(thread(nearest_point, limits[i], limits[i+1], ref(w_hat_vec), ref(violated_bags)));
        }
        /* Join the threads with the main thread */
        for(auto &t : threads){
            t.join();
        }
        time_in_nfp = get_wall_time() - time_in_nfp;
        time_in_all_nfp += time_in_nfp;
        for (unsigned index = 0; index < nb_bags; index++) {
            auto b = violated_bags[index];
            DebugOff("\nBag number " << b._id);
            DebugOff("\nNodes: n1 = " << b._nodes[0]->_name << ", n2 = " << b._nodes[1]->_name << ", n3 = " << b._nodes[2]->_name);
            auto w_hat(move(w_hat_vec[index]));
//            R_star.resize(numcuts+1);
//            R_star[numcuts] = param<>("R_star"+to_string(numcuts));
//            I_star.resize(numcuts+1);
//            I_star[numcuts] = param<>("I_star"+to_string(numcuts));
//            W_star.resize(numcuts+1);
//            W_star[numcuts] = param<>("W_star"+to_string(numcuts));
            R_diff.resize(numcuts+1);
            R_diff[numcuts] = param<>("R_diff"+to_string(numcuts));
            I_diff.resize(numcuts+1);
            I_diff[numcuts] = param<>("I_diff"+to_string(numcuts));
            W_diff.resize(numcuts+1);
            W_diff[numcuts] = param<>("W_diff"+to_string(numcuts));
            R_hat.resize(numcuts+1);
            R_hat[numcuts] = param<>("R_hat"+to_string(numcuts));
            I_hat.resize(numcuts+1);
            I_hat[numcuts] = param<>("I_hat"+to_string(numcuts));
            W_hat.resize(numcuts+1);
            W_hat[numcuts] = param<>("W_hat"+to_string(numcuts));

            node_pairs b_pairs("node_pairs");
            Constraint sdpcut("sdpcut_" + to_string(numcuts));
            Node *ni;
            Arc *aij;
            for (int i = 0; i < b._nodes.size(); i++) {
                for (int j = i; j < b._nodes.size(); j++) {
                    if (i == j) {
                        namew = "w(" + b._nodes[i]->_name + ")";
                        ni = b._nodes[i];
                        W_diff[numcuts].set_val(ni->_name,((Bus *) ni)->w - w_hat(namew).eval());
                        W_hat[numcuts].set_val(ni->_name,w_hat(namew).eval());
//                        W_star[numcuts].set_val(ni->_name,((Bus *) ni)->w);
                    }
                    else {
                        aij = grid.get_arc(b._nodes[i]->_name, b._nodes[j]->_name);
                        if(aij->_imaginary && !aij->_active) {
                            bus_pairs_sdp._keys.push_back(new index_pair(index_(aij->_src->_name), index_(aij->_dest->_name)));
                            aij->_active = true;
                        }
                        b_pairs._keys.push_back(new index_pair(index_(aij->_src->_name), index_(aij->_dest->_name)));

                        namewr = "wr(" + b._nodes[i]->_name + "," + b._nodes[j]->_name + ")";
                        namewi = "wi(" + b._nodes[i]->_name + "," + b._nodes[j]->_name + ")";

                        R_diff[numcuts].set_val(aij->_src->_name + "," + aij->_dest->_name,((Line *) aij)->wr - w_hat(namewr).eval());
                        I_diff[numcuts].set_val(aij->_src->_name + "," + aij->_dest->_name,((Line *) aij)->wi - w_hat(namewi).eval());
                        R_hat[numcuts].set_val(aij->_src->_name + "," + aij->_dest->_name,w_hat(namewr).eval());
                        I_hat[numcuts].set_val(aij->_src->_name + "," + aij->_dest->_name,w_hat(namewi).eval());
//                        R_star[numcuts].set_val(aij->_src->_name + "," + aij->_dest->_name,((Line *) aij)->wr);
//                        I_star[numcuts].set_val(aij->_src->_name + "," + aij->_dest->_name,((Line *) aij)->wi);
                    }
                }
            }
            if(b._nodes.size()>3) hdim_cuts++;

//            cout << "W_HAT=   \n";
//            W_hat[numcuts].print(true);
//            cout << "WI_HAT= \n";
//            I_hat[numcuts].print(true);
            Constraint lin("lin"+to_string(numcuts));
            lin = product(R_diff[numcuts].in(b_pairs),(R_Wij.in(b_pairs) - R_hat[numcuts].in(b_pairs)));
            lin += product(I_diff[numcuts].in(b_pairs),(Im_Wij.in(b_pairs) - I_hat[numcuts].in(b_pairs)));
            lin += product(W_diff[numcuts].in(b._nodes),(Wii.in(b._nodes) - W_hat[numcuts].in(b._nodes)));
//            lin.print_expanded();
            SDP.add(lin <= 0);

            cuts_added++; numcuts++;
        }

//        if(!bus_pairs_sdp._keys.empty()) {
//            SDP.add_indices("SOC",bus_pairs_sdp);
//            bus_pairs_sdp.clear();
//        }


        for(auto& a: grid.arcs){
            if(a->_imaginary && !a->_active) a->_free = true;
        }

        SDPOPF.run(output = 0, relax = false, tol = 1e-6, "ma27", mehrotra = "no");
        gap = 100*(upper_bound - SDP._obj_val)/upper_bound;
        DebugOff("\nPrev = " << prev_opt << ", cur = " << SDP._obj_val);
        DebugOff("\n(opt - prev)/opt = " << (SDP._obj_val - prev_opt)/SDP._obj_val << endl);

        /* Update values for w_star */
        for(auto& arc: grid.arcs){
            ((Line*)arc)->wr = R_Wij(arc->_src->_name+","+arc->_dest->_name).eval();
            ((Line*)arc)->wi = Im_Wij(arc->_src->_name+","+arc->_dest->_name).eval();
        }
        for(auto& node: grid.nodes) ((Bus*)node)->w = Wii(node->_name).eval();

        if(get_wall_time() - solver_time_start > 300) {
            cout << "\nTime limit reached";
            break;
        }

        cout << "\nNum of iterations = " << iter << ", number of cuts = " << numcuts << ", " << cuts_added << endl;
//        if((SDP._obj_val - prev_opt)/SDP._obj_val < fp_tol) unchanged++;
    }

    string lin = "lin";
    int num_act_cuts = 0;
    for (auto &cp: SDP._cons) {
        if(cp.second->get_name().find(lin) != std::string::npos) {
            for (unsigned inst = 0; inst < cp.second->_nb_instances; inst++) {
                cp.second->eval(inst);
                if(cp.second->is_active(inst)) num_act_cuts++;
            }
        }
    }

    cout << "\nEnd: num of iterations = " << iter << ", number of cuts = " << numcuts << ", number of higher dim cuts = " << hdim_cuts;
    cout << "\nNumber of active cuts = " << num_act_cuts;

    double solver_time_end = get_wall_time();
    double total_time_end = get_wall_time();
    auto solve_time = solver_time_end - solver_time_start;
    auto total_time = total_time_end - total_time_start;
    string out = "\nDATA_OPF, " + grid._name + ", " + to_string(nb_buses) + ", " + to_string(nb_lines) +", " + to_string(SDP._obj_val) + ", " + to_string(-numeric_limits<double>::infinity()) + ", " + to_string(solve_time) + ", LocalOptimal, " + to_string(total_time);
    DebugOn(out <<endl);
    DebugOn("Final Gap = " << to_string(gap) << "%."<<endl);
    DebugOn("Upper bound = " << to_string(upper_bound) << "."<<endl);
    DebugOn("Lower bound = " << to_string(SDP._obj_val) << "."<<endl);
    DebugOn("\nResults: " << grid._name << " " << to_string(SDP._obj_val) << " " << to_string(total_time)<<endl);
    DebugOn("\nTime in nfp: " << time_in_all_nfp << endl);
    return 0;
}
