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

int main (int argc, char * argv[]) {
    int output = 0;
    bool relax = false;
    string solver_str = "ipopt";
    SolverType solv_type = ipopt;
    double tol = 1e-6;
    string mehrotra = "no";
    string fname = "../data_sets/Power/nesta_case3_lmbd.m";
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
    PowerNet *grid = new PowerNet();
    grid->readgrid(fname.c_str());

    grid->get_tree_decomp_bags();
    grid->update_net();

    // Grid Parameters
    auto bus_pairs = grid->get_bus_pairs();
    auto bus_pairs_chord = grid->get_bus_pairs_chord();
    auto nb_bus_pairs_chord = bus_pairs_chord.size();
    auto nb_gen = grid->get_nb_active_gens();
    auto nb_lines = grid->get_nb_active_arcs();
    auto nb_buses = grid->get_nb_active_nodes();
    DebugOff("nb gens = " << nb_gen << endl);
    DebugOff("nb lines = " << nb_lines << endl);
    DebugOff("nb buses = " << nb_buses << endl);
    DebugOff("nb bus_pairs = " << nb_bus_pairs_chord << endl);



    /** Build model */
    Model SDP("SDP Model");

    /** Variables */
    /* power generation variables */
    var<Real> Pg("Pg", grid->pg_min.in(grid->gens), grid->pg_max.in(grid->gens));
    var<Real> Qg ("Qg", grid->qg_min.in(grid->gens), grid->qg_max.in(grid->gens));
    SDP.add_var(Pg^(nb_gen));
    SDP.add_var(Qg^(nb_gen));

    /* power flow variables */
    var<Real> Pf_from("Pf_from", grid->S_max.in(grid->arcs));
    var<Real> Qf_from("Qf_from", grid->S_max.in(grid->arcs));
    var<Real> Pf_to("Pf_to", grid->S_max.in(grid->arcs));
    var<Real> Qf_to("Qf_to", grid->S_max.in(grid->arcs));
    
    SDP.add_var(Pf_from^(nb_lines));
    SDP.add_var(Qf_from^(nb_lines));
    SDP.add_var(Pf_to^(nb_lines));
    SDP.add_var(Qf_to^(nb_lines));

    /* Real part of Wij = ViVj */
    var<Real>  R_Wij("R_Wij", grid->wr_min.in(bus_pairs_chord), grid->wr_max.in(bus_pairs_chord));
    /* Imaginary part of Wij = ViVj */
    var<Real>  Im_Wij("Im_Wij", grid->wi_min.in(bus_pairs_chord), grid->wi_max.in(bus_pairs_chord));
    /* Magnitude of Wii = Vi^2 */
    var<Real>  Wii("Wii", grid->w_min.in(grid->nodes), grid->w_max.in(grid->nodes));
    SDP.add_var(Wii^nb_buses);
    SDP.add_var(R_Wij^nb_bus_pairs_chord);
    SDP.add_var(Im_Wij^nb_bus_pairs_chord);

    for(auto& bp: bus_pairs_chord) {
        cout << "Bp: " << bp->_name;
    }

    /* Initialize variables */
    R_Wij.initialize_all(1.0);
//    Im_Wij.initialize_all(0.0);
    Wii.initialize_all(1.001);

    /**  Objective */
    auto obj = product(grid->c1, Pg) + product(grid->c2, power(Pg,2)) + sum(grid->c0);
    obj.print_expanded();
    SDP.min(obj.in(grid->gens));
    
    /** Constraints */
    /* Second-order cone constraints */
    Constraint SOC("SOC");
    SOC = power(R_Wij, 2) + power(Im_Wij, 2) - Wii.from()*Wii.to();
    SDP.add_constraint(SOC.in(bus_pairs) <= 0);
    
    /* Flow conservation */
    Constraint KCL_P("KCL_P");
    KCL_P  = sum(Pf_from.out_arcs()) + sum(Pf_to.in_arcs()) + grid->pl - sum(Pg.in_gens()) + grid->gs*Wii;
    SDP.add_constraint(KCL_P.in(grid->nodes) == 0);
    
    Constraint KCL_Q("KCL_Q");
    KCL_Q  = sum(Qf_from.out_arcs()) + sum(Qf_to.in_arcs()) + grid->ql - sum(Qg.in_gens()) - grid->bs*Wii;
    SDP.add_constraint(KCL_Q.in(grid->nodes) == 0);
    
    /* AC Power Flow */
    Constraint Flow_P_From("Flow_P_From");
    Flow_P_From = Pf_from - (grid->g_ff*Wii.from() + grid->g_ft*R_Wij.in_pairs() + grid->b_ft*Im_Wij.in_pairs());
    SDP.add_constraint(Flow_P_From.in(grid->arcs) == 0);
    
    Constraint Flow_P_To("Flow_P_To");
    Flow_P_To = Pf_to - (grid->g_tt*Wii.to() + grid->g_tf*R_Wij.in_pairs() - grid->b_tf*Im_Wij.in_pairs());
    SDP.add_constraint(Flow_P_To.in(grid->arcs) == 0);
    
    Constraint Flow_Q_From("Flow_Q_From");
    Flow_Q_From = Qf_from - (grid->g_ft*Im_Wij.in_pairs() - grid->b_ff*Wii.from() - grid->b_ft*R_Wij.in_pairs());
    SDP.add_constraint(Flow_Q_From.in(grid->arcs) == 0);
    
    Constraint Flow_Q_To("Flow_Q_To");
    Flow_Q_To = Qf_to + (grid->b_tt*Wii.to() + grid->b_tf*R_Wij.in_pairs() + grid->g_tf*Im_Wij.in_pairs());
    SDP.add_constraint(Flow_Q_To.in(grid->arcs) == 0);
    
    /* Phase Angle Bounds constraints */
    Constraint PAD_UB("PAD_UB");
    PAD_UB = Im_Wij;
    PAD_UB <= grid->tan_th_max*R_Wij;
    SDP.add_constraint(PAD_UB.in(bus_pairs));
    
    Constraint PAD_LB("PAD_LB");
    PAD_LB =  Im_Wij;
    PAD_LB >= grid->tan_th_min*R_Wij;
    SDP.add_constraint(PAD_LB.in(bus_pairs));
    
    /* Thermal Limit Constraints */
    Constraint Thermal_Limit_from("Thermal_Limit_from");
    Thermal_Limit_from = power(Pf_from, 2) + power(Qf_from, 2);
    Thermal_Limit_from <= power(grid->S_max,2);
    SDP.add_constraint(Thermal_Limit_from.in(grid->arcs));
    
    
    Constraint Thermal_Limit_to("Thermal_Limit_to");
    Thermal_Limit_to = power(Pf_to, 2) + power(Qf_to, 2);
    Thermal_Limit_to <= power(grid->S_max,2);
    SDP.add_constraint(Thermal_Limit_to.in(grid->arcs));
    
    /* Lifted Nonlinear Cuts */
//    Constraint LNC1("LNC1");
//    LNC1 = (grid->v_min.from()+grid->v_max.from())*(grid->v_min.to()+grid->v_max.to())*(sin(0.5*(grid->th_max+grid->th_min))*Im_Wij + cos(0.5*(grid->th_max+grid->th_min))*R_Wij);
//    LNC1 -= grid->v_max.to()*cos(0.5*(grid->th_max-grid->th_min))*(grid->v_min.to()+grid->v_max.to())*Wii.from();
//    LNC1 -= grid->v_max.from()*cos(0.5*(grid->th_max-grid->th_min))*(grid->v_min.from()+grid->v_max.from())*Wii.to();
//    LNC1 -= grid->v_max.from()*grid->v_max.to()*cos(0.5*(grid->th_max-grid->th_min))*(grid->v_min.from()*grid->v_min.to() - grid->v_max.from()*grid->v_max.to());
//    SDP.add_constraint(LNC1.in(bus_pairs) >= 0);
//    
//    Constraint LNC2("LNC2");
//    LNC2 = (grid->v_min.from()+grid->v_max.from())*(grid->v_min.to()+grid->v_max.to())*(sin(0.5*(grid->th_max+grid->th_min))*Im_Wij + cos(0.5*(grid->th_max+grid->th_min))*R_Wij);
//    LNC2 -= grid->v_min.to()*cos(0.5*(grid->th_max-grid->th_min))*(grid->v_min.to()+grid->v_max.to())*Wii.from();
//    LNC2 -= grid->v_min.from()*cos(0.5*(grid->th_max-grid->th_min))*(grid->v_min.from()+grid->v_max.from())*Wii.to();
//    LNC2 += grid->v_min.from()*grid->v_min.to()*cos(0.5*(grid->th_max-grid->th_min))*(grid->v_min.from()*grid->v_min.to() - grid->v_max.from()*grid->v_max.to());
//    SDP.add_constraint(LNC2.in(bus_pairs) >= 0);

    /* Solver selection */
    solver SDPOPF(SDP,solv_type);
    double solver_time_start = get_wall_time();
    SDPOPF.run(output = 0, relax = false);
//    SDPOPF.run(output = 0, relax = false, tol = 1e-6, "ma27", mehrotra = "no");

    for(auto& arc: grid->arcs){
        ((Line*)arc)->wr = R_Wij(arc->_src->_name+","+arc->_dest->_name).eval();
        ((Line*)arc)->wi = Im_Wij(arc->_src->_name+","+arc->_dest->_name).eval();
    }
    for(auto& node: grid->nodes) ((Bus*)node)->w = Wii(node->_name).eval();

    vector<Bag*> bags;
    int bagid = 0;
    for(auto& b: grid->_bags){
//        cout << "\nSize = " << b.size();
        bags.push_back(new Bag(bagid,grid,b));
        bagid++;
    }
    param<double> what;
    string namew, namewr, namewi;
    int numcuts = 0;
    node_pairs bus_pairs_sdp;
    double prev_opt = 0;
    double fp_tol = 0.0001;
    cout << "\nBags size = " << bags.size();

    vector<param<double>> R_star, I_star, W_star;
    vector<param<double>> R_diff, I_diff, W_diff;
    vector<param<double>> R_hat, I_hat, W_hat;

    
    while((SDP._obj_val - prev_opt)/SDP._obj_val > fp_tol) {
        prev_opt = SDP._obj_val;
        for (auto &b: bags) {
            if (b->is_PSD()) {
                continue;
            }
            R_star.resize(numcuts+1);
            R_star[numcuts] = param<>("R_star"+to_string(numcuts));
            I_star.resize(numcuts+1);
            I_star[numcuts] = param<>("I_star"+to_string(numcuts));
            W_star.resize(numcuts+1);
            W_star[numcuts] = param<>("W_star"+to_string(numcuts));
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

            what = b->nfp();
            node_pairs b_pairs;
//            param<> Wii_star("Wii_star"), Wii_hat("Wii_hat"), W_diff("W_Diff");
            Constraint sdpcut("sdpcut_" + to_string(numcuts));
            Node *ni;
            Arc *aij;
            double sdp_cst = 0;
            for (int i = 0; i < b->_nodes.size(); i++) {
                for (int j = i; j < b->_nodes.size(); j++) {
                    if (i == j) {
                        namew = "w(" + b->_nodes[i]->_name + ")";
                        ni = b->_nodes[i];
                        W_diff[numcuts].set_val(ni->_name,((Bus *) ni)->w - what(namew).eval());
                        W_hat[numcuts].set_val(ni->_name,what(namew).eval());
                        W_star[numcuts].set_val(ni->_name,((Bus *) ni)->w);
                    }
                    else {
                        aij = grid->get_arc(b->_nodes[i]->_name, b->_nodes[j]->_name);
                        if(aij->_imaginary && !aij->_active) {
                            bus_pairs_sdp._keys.push_back(new index_pair(index_(aij->_src->_name), index_(aij->_dest->_name)));
                            aij->_active = true;
                        }
                        b_pairs._keys.push_back(new index_pair(index_(aij->_src->_name), index_(aij->_dest->_name)));

                        namewr = "wr(" + b->_nodes[i]->_name + "," + b->_nodes[j]->_name + ")";
                        namewi = "wi(" + b->_nodes[i]->_name + "," + b->_nodes[j]->_name + ")";

                        R_diff[numcuts].set_val(aij->_src->_name + "," + aij->_dest->_name,((Line *) aij)->wr - what(namewr).eval());
                        I_diff[numcuts].set_val(aij->_src->_name + "," + aij->_dest->_name,((Line *) aij)->wi - what(namewi).eval());
                        R_hat[numcuts].set_val(aij->_src->_name + "," + aij->_dest->_name,what(namewr).eval());
                        I_hat[numcuts].set_val(aij->_src->_name + "," + aij->_dest->_name,what(namewi).eval());
                        R_star[numcuts].set_val(aij->_src->_name + "," + aij->_dest->_name,((Line *) aij)->wr);
                        I_star[numcuts].set_val(aij->_src->_name + "," + aij->_dest->_name,((Line *) aij)->wi);
                    }
                }
            }

            Constraint lin("lin"+to_string(numcuts));
            lin = product(R_diff[numcuts].in(b_pairs),(R_Wij.in(b_pairs) - R_hat[numcuts].in(b_pairs)));
            lin += product(I_diff[numcuts].in(b_pairs),(Im_Wij.in(b_pairs) - I_hat[numcuts].in(b_pairs)));
            lin += product(W_diff[numcuts].in(b->_nodes),(Wii.in(b->_nodes) - W_hat[numcuts].in(b->_nodes)));
//            lin = product(R_star.in(b_pairs._keys) - R_hat.in(b_pairs._keys),(R_Wij.in(b_pairs._keys) - R_hat.in(b_pairs._keys)));
//            lin += product(I_star.in(b_pairs._keys) - I_hat.in(b_pairs._keys),(Im_Wij.in(b_pairs._keys) - I_hat.in(b_pairs._keys)));
//            lin += product(W_star.in(b->_nodes) - W_hat.in(b->_nodes),(Wii.in(b->_nodes) - W_hat.in(b->_nodes)));
//            lin.print();
            lin.print_expanded();
            SDP.add_constraint(lin <= 0);

            numcuts++;
        }

        if(!bus_pairs_sdp._keys.empty()) {
            for(auto& bp: bus_pairs_sdp._keys) {
                cout << "Bp sdp: " << bp->_name;
            }

            Constraint SOC_im("SOC_im"+to_string(numcuts));
            SOC_im = power(R_Wij, 2) + power(Im_Wij, 2) - Wii.from() * Wii.to();
            SDP.add_constraint(SOC_im.in(bus_pairs_sdp._keys) <= 0);
            bus_pairs_sdp._keys.clear();
        }


        for(auto& a: grid->arcs){
            if(a->_imaginary && !a->_active) a->_free = true;
        }
        
        SDPOPF.run(output = 0, relax = false, tol = 1e-6, "ma27", mehrotra = "no");
        DebugOn("\n(opt - prev)/opt = " << (SDP._obj_val - prev_opt)/SDP._obj_val << endl);

        /* Update values for w_star */
        for(auto& arc: grid->arcs){
            ((Line*)arc)->wr = R_Wij(arc->_src->_name+","+arc->_dest->_name).eval();
            ((Line*)arc)->wi = Im_Wij(arc->_src->_name+","+arc->_dest->_name).eval();
        }
        for(auto& node: grid->nodes) ((Bus*)node)->w = Wii(node->_name).eval();
    }



    double solver_time_end = get_wall_time();
    double total_time_end = get_wall_time();
    auto solve_time = solver_time_end - solver_time_start;
    auto total_time = total_time_end - total_time_start;
    string out = "DATA_OPF, " + grid->_name + ", " + to_string(nb_buses) + ", " + to_string(nb_lines) +", " + to_string(SDP._obj_val) + ", " + to_string(-numeric_limits<double>::infinity()) + ", " + to_string(solve_time) + ", LocalOptimal, " + to_string(total_time);
    DebugOn(out <<endl);
    return 0;
}
