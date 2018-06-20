//
// Created by Ksenia on 11/05/18.
//

#include <stdio.h>
#include <iostream>
#include <string>
#include <stdio.h>
#include <cstring>
#include <fstream>
#include "PowerNet.h"
#include <gravity/solver.h>
#include <stdio.h>
#include <stdlib.h>
#include <optionParser.hpp>
#include <gravity/param.h>


int main (int argc, char * argv[]) {
    int output = 0;
    bool relax = false;
    bool decompose = true;
    string solver_str = "mosek";
    SolverType solv_type = Mosek;
    double tol = 1e-6;
    string mehrotra = "no";
//    string fname = "../data_sets/Power/nesta_case3_lmbd.m";
    string fname = "../data_sets/Power/nesta_case5_pjm.m";
//    string fname = "../nesta-0.7.0/opf/api/nesta_case24_ieee_rts__api.m";

    // create a OptionParser with options
    op::OptionParser opt;
    opt.add_option("h", "help",
                   "shows option help"); // no default value means boolean options, which default value is false
    opt.add_option("f", "file", "Input file name", fname);
    opt.add_option("s", "solver", "Solvers: ipopt/cplex/gurobi, default = ipopt", solver_str);
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
    }

    double total_time_start = get_wall_time();
    PowerNet grid;
    grid.readgrid(fname.c_str());

    grid.get_tree_decomp_bags();
    grid.update_net();

    if(!decompose) grid.fill_wbnds();

    auto bus_pairs = grid.get_bus_pairs();
    auto bus_pairs_all = grid.get_bus_pairs_all();
    auto bus_pairs_chord = grid.get_bus_pairs_chord();
//    auto nb_bus_pairs_chord = bus_pairs_chord.size();
    auto nb_gen = grid.get_nb_active_gens();
    auto nb_lines = grid.get_nb_active_arcs();
    auto nb_buses = grid.get_nb_active_nodes();

    /** Build model */
    Model SDP("SDP Model");

    /** Variables */
    /* power generation variables */
    var<Real> Pg("Pg", grid.pg_min, grid.pg_max);
    var<Real> Qg ("Qg", grid.qg_min, grid.qg_max);
    SDP.add_var(Pg.in(grid.gens));
    SDP.add_var(Qg.in(grid.gens));

    /* power flow variables */
    var<Real> Pf_from("Pf_from", grid.S_max);
    var<Real> Qf_from("Qf_from", grid.S_max);
    var<Real> Pf_to("Pf_to", grid.S_max);
    var<Real> Qf_to("Qf_to", grid.S_max);
    SDP.add_var(Pf_from.in(grid.arcs));
    SDP.add_var(Qf_from.in(grid.arcs));
    SDP.add_var(Pf_to.in(grid.arcs));
    SDP.add_var(Qf_to.in(grid.arcs));

    /* Real part of Wij = ViVj */
    var<Real>  R_Wij("R_Wij", grid.wr_min, grid.wr_max);
    /* Imaginary part of Wij = ViVj */
    var<Real>  Im_Wij("Im_Wij", grid.wi_min, grid.wi_max);
    /* Magnitude of Wii = Vi^2 */
    var<Real>  Wii("Wii", grid.w_min, grid.w_max);
    SDP.add_var(Wii.in(grid.nodes));
    if(decompose) {
        SDP.add_var(R_Wij.in(bus_pairs_chord));
        SDP.add_var(Im_Wij.in(bus_pairs_chord));
    }else{
        SDP.add_var(R_Wij.in(bus_pairs_all));
        SDP.add_var(Im_Wij.in(bus_pairs_all));
    }

    /* Initialize variables */
    R_Wij.initialize_all(1.0);
    Wii.initialize_all(1.001);

    /**  Objective */
    auto obj = product(grid.c1, Pg) + product(grid.c2, power(Pg,2)) + sum(grid.c0);
    obj.print_expanded();
    SDP.min(obj.in(grid.gens));

    /** Constraints */

    /* Second-order cone constraints */
//    Constraint SOC("SOC");
//    SOC = power(R_Wij, 2) + power(Im_Wij, 2) - Wii.from()*Wii.to();
//    SDP.add_constraint(SOC.in(bus_pairs) <= 0);
//    SOC.print_expanded();

    /* Flow conservation */
    Constraint KCL_P("KCL_P");
    KCL_P  = sum(Pf_from.out_arcs()) + sum(Pf_to.in_arcs()) + grid.pl - sum(Pg.in_gens()) + grid.gs*Wii;
    SDP.add_constraint(KCL_P.in(grid.nodes) == 0);
//    KCL_P.print_expanded();

    Constraint KCL_Q("KCL_Q");
    KCL_Q  = sum(Qf_from.out_arcs()) + sum(Qf_to.in_arcs()) + grid.ql - sum(Qg.in_gens()) - grid.bs*Wii;
    SDP.add_constraint(KCL_Q.in(grid.nodes) == 0);
//    KCL_Q.print_expanded();

    /* AC Power Flow */
    Constraint Flow_P_From("Flow_P_From");
    Flow_P_From = Pf_from - (grid.g_ff*Wii.from() + grid.g_ft*R_Wij.in_pairs() + grid.b_ft*Im_Wij.in_pairs());
    SDP.add_constraint(Flow_P_From.in(grid.arcs) == 0);
//    Flow_P_From.print_expanded();

    Constraint Flow_P_To("Flow_P_To");
    Flow_P_To = Pf_to - (grid.g_tt*Wii.to() + grid.g_tf*R_Wij.in_pairs() - grid.b_tf*Im_Wij.in_pairs());
    SDP.add_constraint(Flow_P_To.in(grid.arcs) == 0);
//    Flow_P_To.print_expanded();

    Constraint Flow_Q_From("Flow_Q_From");
    Flow_Q_From = Qf_from - (grid.g_ft*Im_Wij.in_pairs() - grid.b_ff*Wii.from() - grid.b_ft*R_Wij.in_pairs());
    SDP.add_constraint(Flow_Q_From.in(grid.arcs) == 0);
//    Flow_Q_From.print_expanded();

    Constraint Flow_Q_To("Flow_Q_To");
    Flow_Q_To = Qf_to + (grid.b_tt*Wii.to() + grid.b_tf*R_Wij.in_pairs() + grid.g_tf*Im_Wij.in_pairs());
    SDP.add_constraint(Flow_Q_To.in(grid.arcs) == 0);
//    Flow_Q_To.print_expanded();

    /* Phase Angle Bounds constraints */

    Constraint PAD_UB("PAD_UB");
    PAD_UB = Im_Wij;
    PAD_UB <= grid.tan_th_max*R_Wij;
//    SDP.add_constraint(PAD_UB.in(bus_pairs));
//    PAD_UB.print_expanded();

    Constraint PAD_LB("PAD_LB");
    PAD_LB =  Im_Wij;
    PAD_LB >= grid.tan_th_min*R_Wij;
//    SDP.add_constraint(PAD_LB.in(bus_pairs));
//    PAD_LB.print_expanded();

    /* Thermal Limit Constraints */ //TODO
    Constraint Thermal_Limit_from("Thermal_Limit_from");
    Thermal_Limit_from = power(Pf_from, 2) + power(Qf_from, 2);
    Thermal_Limit_from <= power(grid.S_max,2);
    SDP.add_constraint(Thermal_Limit_from.in(grid.arcs));
//    Thermal_Limit_from.print_expanded();

    Constraint Thermal_Limit_to("Thermal_Limit_to");
    Thermal_Limit_to = power(Pf_to, 2) + power(Qf_to, 2);
    Thermal_Limit_to <= power(grid.S_max,2);
    SDP.add_constraint(Thermal_Limit_to.in(grid.arcs));

    /* Lifted Nonlinear Cuts */
    Constraint LNC1("LNC1");
    LNC1 = (grid.v_min.from()+grid.v_max.from())*(grid.v_min.to()+grid.v_max.to())*(sin(0.5*(grid.th_max+grid.th_min))*Im_Wij + cos(0.5*(grid.th_max+grid.th_min))*R_Wij);
    LNC1 -= grid.v_max.to()*cos(0.5*(grid.th_max-grid.th_min))*(grid.v_min.to()+grid.v_max.to())*Wii.from();
    LNC1 -= grid.v_max.from()*cos(0.5*(grid.th_max-grid.th_min))*(grid.v_min.from()+grid.v_max.from())*Wii.to();
    LNC1 -= grid.v_max.from()*grid.v_max.to()*cos(0.5*(grid.th_max-grid.th_min))*(grid.v_min.from()*grid.v_min.to() - grid.v_max.from()*grid.v_max.to());
//    SDP.add_constraint(LNC1.in(bus_pairs) >= 0);
//    LNC1.print_expanded();

    Constraint LNC2("LNC2");
    LNC2 = (grid.v_min.from()+grid.v_max.from())*(grid.v_min.to()+grid.v_max.to())*(sin(0.5*(grid.th_max+grid.th_min))*Im_Wij + cos(0.5*(grid.th_max+grid.th_min))*R_Wij);
    LNC2 -= grid.v_min.to()*cos(0.5*(grid.th_max-grid.th_min))*(grid.v_min.to()+grid.v_max.to())*Wii.from();
    LNC2 -= grid.v_min.from()*cos(0.5*(grid.th_max-grid.th_min))*(grid.v_min.from()+grid.v_max.from())*Wii.to();
    LNC2 += grid.v_min.from()*grid.v_min.to()*cos(0.5*(grid.th_max-grid.th_min))*(grid.v_min.from()*grid.v_min.to() - grid.v_max.from()*grid.v_max.to());
//    SDP.add_constraint(LNC2.in(bus_pairs) >= 0);
//    LNC2.print_expanded();

    vector<var<double>> W; // store the matrix variables
    int bagid = 0;

    for (auto& a: grid.arcs) {
        vector<Node*> v;
        v.push_back(a->_src);
        v.push_back(a->_dest);
        grid._bags.push_back(v);
    }

    if(!decompose) {
        grid._bags.clear();
        vector<Node*> bag = grid.nodes;
        grid._bags.push_back(bag);
    }

    //todo: for SOCP: add lines to bags?
    double solver_time_start = get_wall_time();
    for(auto& b: grid._bags){
        int n = b.size();
        cout << "\nn = " << n;
        var<double> W_b("W"+to_string(bagid));
        W_b._psd = true;
        SDP.add_var(W_b.in(R(2*n,2*n)));
        W.push_back(W_b);

        //connect W to the w vars: if all nodes are sorted and all arcs are proper direction, no need to worry about reversed etc?

        vector<index_> wij_idxs, w_idxs, zero_idxs; //todo: need to save all indices?
        vector<index_> d_u_idxs, ul_u_idxs, ur_u_idxs, ur_l_idxs, lr_u_idxs, d_l_idxs;

        for(int i = 0; i < n; i++) {
            for(int j = i; j < n; j++) {
                if(i==j) { //diagonal
                    w_idxs.push_back(index_(b[i]->_name));
                    d_u_idxs.push_back(index_(to_string(i)+","+to_string(j)));
                    d_l_idxs.push_back(index_(to_string(i+n)+","+to_string(j+n)));
                    zero_idxs.push_back(index_(to_string(i)+","+to_string(j+n)));
                }else{ // off diagonal
                    wij_idxs.push_back(index_(b[i]->_name+","+b[j]->_name));
                    ul_u_idxs.push_back(index_(to_string(i)+","+to_string(j)));
                    ur_u_idxs.push_back(index_(to_string(i)+","+to_string(j+n)));
                    ur_l_idxs.push_back(index_(to_string(j)+","+to_string(i+n)));
                    lr_u_idxs.push_back(index_(to_string(i+n)+","+to_string(j+n)));
                }
            }
        }

        Constraint upper_diag("upper_diag"+to_string(bagid));
        upper_diag = W_b.in(d_u_idxs) - Wii.in(w_idxs);
        SDP.add_constraint(upper_diag == 0);
//        upper_diag.print_expanded();

        Constraint lower_diag("lower_diag"+to_string(bagid));
        lower_diag = W_b.in(d_l_idxs) - Wii.in(w_idxs);
        SDP.add_constraint(lower_diag == 0);
//        lower_diag.print_expanded();

        Constraint zeros("zeros"+to_string(bagid));
        zeros = W_b.in(zero_idxs);
        SDP.add_constraint(zeros == 0);
//        zeros.print_expanded();

        Constraint W_wR1("W_wR1"+to_string(bagid));
        W_wR1 = W_b.in(ul_u_idxs) - R_Wij.in(wij_idxs);
        SDP.add_constraint(W_wR1 == 0);
//        W_wR1.print_expanded();

        Constraint W_wR2("W_wR2"+to_string(bagid));
        W_wR2 = W_b.in(lr_u_idxs) - R_Wij.in(wij_idxs);
        SDP.add_constraint(W_wR2 == 0);
//        W_wR2.print_expanded();

        Constraint W_wI1("W_wI1"+to_string(bagid));
        W_wI1 = W_b.in(ur_u_idxs) - Im_Wij.in(wij_idxs);
        SDP.add_constraint(W_wI1 == 0);

        Constraint W_wI2("W_wI2"+to_string(bagid));
        W_wI2 = W_b.in(ur_l_idxs) + Im_Wij.in(wij_idxs);
        SDP.add_constraint(W_wI2 == 0);

        bagid++;
    }

    solver s(SDP,Mosek);
    s.run(0,0);

    double solver_time_end = get_wall_time();
    double total_time_end = get_wall_time();
    auto solve_time = solver_time_end - solver_time_start;
    auto total_time = total_time_end - total_time_start;
    string out = "\nDATA_OPF, " + grid._name + ", " + to_string(nb_buses) + ", " + to_string(nb_lines) +", " + to_string(SDP._obj_val) + ", " + to_string(-numeric_limits<double>::infinity()) + ", " + to_string(solve_time) + ", LocalOptimal, " + to_string(total_time);
    DebugOn(out <<endl);

    return 0;
}