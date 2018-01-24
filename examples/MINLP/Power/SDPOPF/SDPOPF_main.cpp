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

vector<param<>> signs(Net* net, const std::vector<std::vector<Node*>>& bags) {
    vector<param<>> res;
    string key;
    res.resize(3);
    for (int i = 0; i<3; i++) {
        res[i].set_name("Im_Wij_signs_"+to_string(i));
    }
    set<vector<unsigned>> ids;
    for (auto &bag: net->_bags) {
        if (bag.size() != 3) {
            continue;
        }
        vector<unsigned> ids_bag;
        for (int i = 0; i<3; i++) {
            ids_bag.push_back(bag[i]->_id);
        }
        if(ids.count(ids_bag)==0) {
            ids.insert(ids_bag);
        } else {
            continue;
        }
        for (int i = 0; i< 2; i++) {
            if(net->has_directed_arc(bag[i], bag[i+1])) {
                key = bag[i]->_name + "," + bag[i+1]->_name;
                res[i].set_val(key,1.0);
            }
            else {
                key = bag[i+1]->_name + "," + bag[i]->_name;
                res[i].set_val(key,-1.0);
            }
        }
        /* Loop back pair */
        if(net->has_directed_arc(bag[0], bag[2])) {
            key = bag[0]->_name + "," + bag[2]->_name;
            res[2].set_val(key,1.0);
        }
        else{
            key = bag[2]->_name + "," + bag[0]->_name;
            res[2].set_val(key,-1.0);
        }
    }
    return res;
}

int main (int argc, char * argv[]) {
    int output = 0;
    bool relax = false;
    string solver_str = "ipopt";
    SolverType solv_type = ipopt;
    double tol = 1e-6;
    string mehrotra = "no";
//    string fname = "../data_sets/Power/nesta_case3_lmbd.m";
    string fname = "../data_sets/Power/nesta_case5_pjm.m";
    fname = "/Users/hlh/Dropbox/Work/Dev/power_models/data/nesta_api/nesta_case30_fsr__api.m";
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

    // Grid Parameters
    auto bus_pairs = grid.get_bus_pairs();
    auto bus_pairs_chord = grid.get_bus_pairs_chord();
//    auto nb_bus_pairs_chord = bus_pairs_chord.size();
    auto nb_gen = grid.get_nb_active_gens();
    auto nb_lines = grid.get_nb_active_arcs();
    auto nb_buses = grid.get_nb_active_nodes();
    DebugOff("nb gens = " << nb_gen << endl);
    DebugOff("nb lines = " << nb_lines << endl);
    DebugOff("nb buses = " << nb_buses << endl);
    DebugOff("nb bus_pairs = " << nb_bus_pairs_chord << endl);



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
    SDP.add_var(R_Wij.in(bus_pairs_chord));
    SDP.add_var(Im_Wij.in(bus_pairs_chord));
//    SDP.add_var(R_Wij.in(bus_pairs));
//    SDP.add_var(Im_Wij.in(bus_pairs));
    
    /* Initialize variables */
    R_Wij.initialize_all(1.0);
    Wii.initialize_all(1.001);

    /**  Objective */
    auto obj = product(grid.c1, Pg) + product(grid.c2, power(Pg,2)) + sum(grid.c0);
//    obj.print_expanded();
    SDP.min(obj.in(grid.gens));
    
    /** Constraints */
    /* Second-order cone constraints */
    Constraint SOC("SOC");
    SOC = power(R_Wij, 2) + power(Im_Wij, 2) - Wii.from()*Wii.to();
    SDP.add_constraint(SOC.in(bus_pairs) <= 0);
    
    /* Flow conservation */
    Constraint KCL_P("KCL_P");
    KCL_P  = sum(Pf_from.out_arcs()) + sum(Pf_to.in_arcs()) + grid.pl - sum(Pg.in_gens()) + grid.gs*Wii;
    SDP.add_constraint(KCL_P.in(grid.nodes) == 0);
    
    Constraint KCL_Q("KCL_Q");
    KCL_Q  = sum(Qf_from.out_arcs()) + sum(Qf_to.in_arcs()) + grid.ql - sum(Qg.in_gens()) - grid.bs*Wii;
    SDP.add_constraint(KCL_Q.in(grid.nodes) == 0);
    
    /* AC Power Flow */
    Constraint Flow_P_From("Flow_P_From");
    Flow_P_From = Pf_from - (grid.g_ff*Wii.from() + grid.g_ft*R_Wij.in_pairs() + grid.b_ft*Im_Wij.in_pairs());
    SDP.add_constraint(Flow_P_From.in(grid.arcs) == 0);
    
    Constraint Flow_P_To("Flow_P_To");
    Flow_P_To = Pf_to - (grid.g_tt*Wii.to() + grid.g_tf*R_Wij.in_pairs() - grid.b_tf*Im_Wij.in_pairs());
    SDP.add_constraint(Flow_P_To.in(grid.arcs) == 0);
    
    Constraint Flow_Q_From("Flow_Q_From");
    Flow_Q_From = Qf_from - (grid.g_ft*Im_Wij.in_pairs() - grid.b_ff*Wii.from() - grid.b_ft*R_Wij.in_pairs());
    SDP.add_constraint(Flow_Q_From.in(grid.arcs) == 0);
    
    Constraint Flow_Q_To("Flow_Q_To");
    Flow_Q_To = Qf_to + (grid.b_tt*Wii.to() + grid.b_tf*R_Wij.in_pairs() + grid.g_tf*Im_Wij.in_pairs());
    SDP.add_constraint(Flow_Q_To.in(grid.arcs) == 0);
    
    /* Phase Angle Bounds constraints */

//    Constraint PAD_UB("PAD_UB");
//    PAD_UB = Im_Wij;
//    PAD_UB <= grid.tan_th_max*R_Wij;
//    SDP.add_constraint(PAD_UB.in(bus_pairs));
//
//    Constraint PAD_LB("PAD_LB");
//    PAD_LB =  Im_Wij;
//    PAD_LB >= grid.tan_th_min*R_Wij;
//    SDP.add_constraint(PAD_LB.in(bus_pairs));
    
    /* Thermal Limit Constraints */
    Constraint Thermal_Limit_from("Thermal_Limit_from");
    Thermal_Limit_from = power(Pf_from, 2) + power(Qf_from, 2);
    Thermal_Limit_from <= power(grid.S_max,2);
    SDP.add_constraint(Thermal_Limit_from.in(grid.arcs));
    
    
    Constraint Thermal_Limit_to("Thermal_Limit_to");
    Thermal_Limit_to = power(Pf_to, 2) + power(Qf_to, 2);
    Thermal_Limit_to <= power(grid.S_max,2);
    SDP.add_constraint(Thermal_Limit_to.in(grid.arcs));
    
    /* Lifted Nonlinear Cuts */
//    Constraint LNC1("LNC1");
//    LNC1 = (grid.v_min.from()+grid.v_max.from())*(grid.v_min.to()+grid.v_max.to())*(sin(0.5*(grid.th_max+grid.th_min))*Im_Wij + cos(0.5*(grid.th_max+grid.th_min))*R_Wij);
//    LNC1 -= grid.v_max.to()*cos(0.5*(grid.th_max-grid.th_min))*(grid.v_min.to()+grid.v_max.to())*Wii.from();
//    LNC1 -= grid.v_max.from()*cos(0.5*(grid.th_max-grid.th_min))*(grid.v_min.from()+grid.v_max.from())*Wii.to();
//    LNC1 -= grid.v_max.from()*grid.v_max.to()*cos(0.5*(grid.th_max-grid.th_min))*(grid.v_min.from()*grid.v_min.to() - grid.v_max.from()*grid.v_max.to());
//    SDP.add_constraint(LNC1.in(bus_pairs) >= 0);
//    
//    Constraint LNC2("LNC2");
//    LNC2 = (grid.v_min.from()+grid.v_max.from())*(grid.v_min.to()+grid.v_max.to())*(sin(0.5*(grid.th_max+grid.th_min))*Im_Wij + cos(0.5*(grid.th_max+grid.th_min))*R_Wij);
//    LNC2 -= grid.v_min.to()*cos(0.5*(grid.th_max-grid.th_min))*(grid.v_min.to()+grid.v_max.to())*Wii.from();
//    LNC2 -= grid.v_min.from()*cos(0.5*(grid.th_max-grid.th_min))*(grid.v_min.from()+grid.v_max.from())*Wii.to();
//    LNC2 += grid.v_min.from()*grid.v_min.to()*cos(0.5*(grid.th_max-grid.th_min))*(grid.v_min.from()*grid.v_min.to() - grid.v_max.from()*grid.v_max.to());
//    SDP.add_constraint(LNC2.in(bus_pairs) >= 0);

    vector<Bag> bags;
    int n3 = 0;
    int bagid = 0;
    for(auto& b: grid._bags){
        bags.push_back(Bag(bagid,grid,b));
        if(b.size()==3) n3++;
        bagid++;
    }

    //        double SDP = wr12*(wr23*a13->wr + wi23*a13->wi) + wi12*(-wi23*a13->wr + wr23*a13->wi);
//        SDP *= 2;
//        SDP -= (wr12*wr12 + wi12*wi12)*w3 + (a13->wr*a13->wr + a13->wi*a13->wi)*w2 + (wr23*wr23 + wi23*wi23)*w1;
//        SDP += w1*w2*w3;

//    DebugOn("\nNum of 3d bags = " << n3);

//    auto R_Wij_ = R_Wij.pairs_in_directed(grid, grid._bags, 3);
//    auto Im_Wij_ = Im_Wij.pairs_in_directed(grid, grid._bags, 3);
//    auto Wii_ = Wii.in(grid._bags, 3);
//    auto I_sgn = signs(grid,grid._bags);
//    DebugOn("\n" << I_sgn[0].to_str(true) << ", " << I_sgn[0].get_nb_instances() << endl);
//    DebugOn("\n" << I_sgn[1].to_str(true) << ", " << I_sgn[1].get_nb_instances() << endl);
//    DebugOn("\n" << I_sgn[2].to_str(true) << ", " << I_sgn[2].get_nb_instances() << endl);
//    DebugOn("\n" << R_Wij_[0].to_str(true) << ", " << R_Wij_[0].get_nb_instances() << endl);
//    DebugOn("\n" << R_Wij_[1].to_str(true) << ", " << R_Wij_[1].get_nb_instances() << endl);
//    DebugOn("\n" << R_Wij_[2].to_str(true) << ", " << R_Wij_[2].get_nb_instances() << endl);
//    DebugOn("\n" << Im_Wij_[0].to_str(true) << ", " << Im_Wij_[0].get_nb_instances()  << endl);
//    DebugOn("\n" << Im_Wij_[1].to_str(true) << ", " << Im_Wij_[1].get_nb_instances()  << endl);
//    DebugOn("\n" << Im_Wij_[2].to_str(true) << ", " << Im_Wij_[2].get_nb_instances()  << endl);
//    DebugOn("\n" << Wii_[1].to_str(true) << ", " << Wii_[1].get_nb_instances()  << endl);
//    Constraint SDP3("SDP_3D");

//    SDP3 = 2*R_Wij_[0]*(R_Wij_[1]*R_Wij_[2] + I_sgn[1]*I_sgn[2]*Im_Wij_[1]*Im_Wij_[2]);
//    SDP3 += 2*I_sgn[0]*Im_Wij_[0]*(R_Wij_[1]*I_sgn[2]*Im_Wij_[2] - I_sgn[1]*Im_Wij_[1]*R_Wij_[2]);
//    SDP3 = R_Wij_[1]*R_Wij_[2];
//    DebugOn("\nsdp nb inst = " << SDP3.get_nb_instances() << endl);
//    SDP3 = 2*R_Wij_[0]*(R_Wij_[1]*R_Wij_[2] + Im_Wij_[1]*Im_Wij_[2]);
//    SDP3 += 2*Im_Wij_[0]*(R_Wij_[1]*Im_Wij_[2] - Im_Wij_[1]*R_Wij_[2]);
//    SDP3 -= (power(R_Wij_[0],2) + power(Im_Wij_[0],2))*Wii_[2];
//    SDP3 -= (power(R_Wij_[1],2) + power(Im_Wij_[1],2))*Wii_[0];
//    SDP3 -= (power(R_Wij_[2],2) + power(Im_Wij_[2],2))*Wii_[1];
//    SDP3 += Wii_[0]*Wii_[1]*Wii_[2];
//    SDP.add_constraint(SDP3 >= 0);
//    SDP3.print_expanded();
//
//    exit(0);

//    SDP3 = -2*Xij_[0]*Xij_[1]*Xij_[2];
//    SDP3 -= Xii_[0]*Xii_[1]*Xii_[2];
//    SDP3 += power(Xij_[0],2)*Xii_[2];
//    SDP3 += power(Xij_[2],2)*Xii_[1];
//    SDP3 += power(Xij_[1],2)*Xii_[0];


    /* Solver selection */
    solver SDPOPF(SDP,solv_type);
    double solver_time_start = get_wall_time();
    SDPOPF.run(output = 0, relax = false);
//    SDPOPF.run(output = 0, relax = false, tol = 1e-6, "ma27", mehrotra = "no");

    for(auto& arc: grid.arcs){
        ((Line*)arc)->wr = R_Wij(arc->_src->_name+","+arc->_dest->_name).eval();
        ((Line*)arc)->wi = Im_Wij(arc->_src->_name+","+arc->_dest->_name).eval();
    }
    for(auto& node: grid.nodes) ((Bus*)node)->w = Wii(node->_name).eval();


    param<double> what;
    string namew, namewr, namewi;
    int numcuts = 0;
    node_pairs bus_pairs_sdp;
//    double prev_opt = 0;
//    double fp_tol = 0.0001;
    DebugOff("\nNumber of bags = " << bags.size());

    vector<param<double>> R_star, I_star, W_star;
    vector<param<double>> R_diff, I_diff, W_diff;
    vector<param<double>> R_hat, I_hat, W_hat;

//    int unchanged = 0;
    int iter = 0, hdim_cuts = 0, cuts_added = 1;

//    while(unchanged < 3) {
    while(cuts_added > 0) {
        cuts_added = 0;
        iter++;
//    while((SDP._obj_val - prev_opt)/SDP._obj_val > fp_tol) {
//        prev_opt = SDP._obj_val;
        for (auto &b: bags) {
            DebugOff("\nBag number " << b._id);
            DebugOff("\nNodes: n1 = " << b._nodes[0]->_name << ", n2 = " << b._nodes[1]->_name << ", n3 = " << b._nodes[2]->_name);

            if (b.is_PSD()) {
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

            what = b.nfp();
            node_pairs b_pairs("node_pairs");
            Constraint sdpcut("sdpcut_" + to_string(numcuts));
            Node *ni;
            Arc *aij;
            for (int i = 0; i < b._nodes.size(); i++) {
                for (int j = i; j < b._nodes.size(); j++) {
                    if (i == j) {
                        namew = "w(" + b._nodes[i]->_name + ")";
                        ni = b._nodes[i];
                        W_diff[numcuts].set_val(ni->_name,((Bus *) ni)->w - what(namew).eval());
                        W_hat[numcuts].set_val(ni->_name,what(namew).eval());
                        W_star[numcuts].set_val(ni->_name,((Bus *) ni)->w);
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

                        R_diff[numcuts].set_val(aij->_src->_name + "," + aij->_dest->_name,((Line *) aij)->wr - what(namewr).eval());
                        I_diff[numcuts].set_val(aij->_src->_name + "," + aij->_dest->_name,((Line *) aij)->wi - what(namewi).eval());
                        R_hat[numcuts].set_val(aij->_src->_name + "," + aij->_dest->_name,what(namewr).eval());
                        I_hat[numcuts].set_val(aij->_src->_name + "," + aij->_dest->_name,what(namewi).eval());
                        R_star[numcuts].set_val(aij->_src->_name + "," + aij->_dest->_name,((Line *) aij)->wr);
                        I_star[numcuts].set_val(aij->_src->_name + "," + aij->_dest->_name,((Line *) aij)->wi);
                    }
                }
            }
            if(b._nodes.size()>3) hdim_cuts++;

            Constraint lin("lin"+to_string(numcuts));
            lin = product(R_diff[numcuts].in(b_pairs),(R_Wij.in(b_pairs) - R_hat[numcuts].in(b_pairs)));
            lin += product(I_diff[numcuts].in(b_pairs),(Im_Wij.in(b_pairs) - I_hat[numcuts].in(b_pairs)));
            lin += product(W_diff[numcuts].in(b._nodes),(Wii.in(b._nodes) - W_hat[numcuts].in(b._nodes)));
//            lin = product(R_star.in(b_pairs._keys) - R_hat.in(b_pairs._keys),(R_Wij.in(b_pairs._keys) - R_hat.in(b_pairs._keys)));
//            lin += product(I_star.in(b_pairs._keys) - I_hat.in(b_pairs._keys),(Im_Wij.in(b_pairs._keys) - I_hat.in(b_pairs._keys)));
//            lin += product(W_star.in(b._nodes) - W_hat.in(b._nodes),(Wii.in(b._nodes) - W_hat.in(b._nodes)));
//            lin.print();
//            lin.print_expanded();
            SDP.add_constraint(lin <= 0);

            cuts_added++; numcuts++;
        }

        if(!bus_pairs_sdp._keys.empty()) {
//            cout << "Adding SOC indices: {";
//            for(auto& bp: bus_pairs_sdp._keys) {
//                cout << bp->_name << ", ";
//            }
//            cout << "}"  << endl;

//            Constraint SOC_im("SOC_im"+to_string(numcuts));
//            SOC_im = power(R_Wij, 2) + power(Im_Wij, 2) - Wii.from() * Wii.to();
//            SDP.add_constraint(SOC_im.in(bus_pairs_sdp) <= 0);
            SDP.add_indices("SOC",bus_pairs_sdp);
//            SDP.get_constraint("SOC")->print_expanded();
            bus_pairs_sdp.clear();
        }


        for(auto& a: grid.arcs){
            if(a->_imaginary && !a->_active) a->_free = true;
        }
        
        SDPOPF.run(output = 0, relax = false, tol = 1e-6, "ma27", mehrotra = "no");
        DebugOff("\nPrev = " << prev_opt << ", cur = " << SDP._obj_val);
        DebugOff("\n(opt - prev)/opt = " << (SDP._obj_val - prev_opt)/SDP._obj_val << endl);

        /* Update values for w_star */
        for(auto& arc: grid.arcs){
            ((Line*)arc)->wr = R_Wij(arc->_src->_name+","+arc->_dest->_name).eval();
            ((Line*)arc)->wi = Im_Wij(arc->_src->_name+","+arc->_dest->_name).eval();
        }
        for(auto& node: grid.nodes) ((Bus*)node)->w = Wii(node->_name).eval();


        cout << "\nNum of iterations = " << iter << ", number of cuts = " << numcuts << endl;
//        if((SDP._obj_val - prev_opt)/SDP._obj_val < fp_tol) unchanged++;
    }
    string lin = "lin";
    int num_act_cuts = 0;
    for (auto &cp: SDP._cons) {
        if(cp.second->get_name().find(lin) != std::string::npos) {
            for (unsigned inst = 0; inst < cp.second->_nb_instances; inst++) {
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
    DebugOn("\nResults: " << grid._name << " " << to_string(SDP._obj_val) << " " << to_string(total_time)<<endl);
    return 0;
}
