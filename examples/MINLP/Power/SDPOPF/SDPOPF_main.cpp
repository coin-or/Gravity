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
    bool relax = false, use_cplex = false;
    double tol = 1e-6;
    string mehrotra = "no";
    string fname = "../data_sets/Power/nesta_case5_pjm.m";
    // create a OptionParser with options
    op::OptionParser opt;
    opt.add_option("h", "help",
                   "shows option help"); // no default value means boolean options, which default value is false
    opt.add_option("f", "file", "Input file name", fname);

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


//    Model NPP("NPP model");
//
//    int n = 3;
//
//    sdpvar<double> W("W");
//    NPP.add_var(W ^ (2*n));
//
//    double* wstar = new double[n*n]; // w* = (w1, wI12, wI13, wR12, w2, wI23, wR13, wR23, w3)
//    wstar[0] = 1.1682; wstar[1] = -0.320016; wstar[2] = 0.337146;
//    wstar[3] = 1.05257; wstar[4] = 1.03605; wstar[5] = 0.30723;
//    wstar[6] = 1.02769; wstar[7] = 0.965777; wstar[8] = 1.00139;
//
//    var<double> z("z"); //difference between w* and w (as vectors)
//    z.in_q_cone();
//    NPP.add_var(z ^ (n*n+1));
//
//    func_ obj;
//    obj = z(0);
//    NPP.set_objective(min(obj));
//
//    cout << "\n" << NPP._obj_val << endl;
//
////    min t  s.t.
////    ||z|| <= t
////    w*-w = z
////    X is PSD, L <= X <= U
//
//    /** Constraints **/
//
//    for(int i = 0; i < n; i++){
//        for(int j = 0; j < n; j++){
//            if(i==j){
//                /* w*-w = z */
//                cout << "\nW(" << i << "," << i << ") goes to " << 3*i+i;
//                Constraint svec("svec"+to_string(i)+to_string(i));
//                svec = wstar[n*i+i] - W(i,i) - z(n*i+i+1);
//                NPP.add_constraint(svec=0);
//
//                /* bounds */
//                Constraint lbound("lbound("+to_string(i)+ "," + to_string(i) +")");
//                Constraint ubound("ubound("+to_string(i)+ "," + to_string(i) +")");
//                lbound = W(i,i) - 0.81;
//                ubound = W(i,i) - 1.21;
//                NPP.add_constraint(lbound >=0);
//                NPP.add_constraint(ubound <=0);
//
//                /* zeros */
//                Constraint zero("zero"+to_string(i)+to_string(i+n));
//                zero = W(i,i+n);
//                NPP.add_constraint(zero=0);
//                Constraint zero2("zero"+to_string(i+n)+to_string(i));
//                zero2 = W(i+n,i);
//                NPP.add_constraint(zero2=0);
//
//                /* matrix structure */
//                Constraint mstruct("mstruct="+to_string(i+n)+to_string(i+n));
//                mstruct = W(i+n,i+n) - W(i,i);
//                NPP.add_constraint(mstruct=0);
//            }
//            else if(i>j){
//                /* w*-w = z */
//                cout << "\nW(" << i << "," << j << ") goes to " << 3*i+j;
//                Constraint svec("svec"+to_string(i)+to_string(j));
//                svec = wstar[n*i+j] - W(i,j) - z(n*i+j+1);
//                NPP.add_constraint(svec=0);
//
//                /* bounds */
//                Constraint lbound("lbound("+to_string(i)+ "," + to_string(j) +")");
//                Constraint ubound("ubound("+to_string(i)+ "," + to_string(j) +")");
//                lbound = W(i,j) - 0.701;
//                ubound = W(i,j) - 1.21;
//                NPP.add_constraint(lbound >=0);
//                NPP.add_constraint(ubound <=0);
//
//                /* matrix structure */
//                Constraint mstruct("mstruct>"+to_string(j)+to_string(i));
//                mstruct = W(j,i) - W(i,j);
//                NPP.add_constraint(mstruct=0);
//                Constraint mstruct2("mstruct>"+to_string(i+n)+to_string(j+n));
//                mstruct2 = W(i+n,j+n) - W(i,j);
//                NPP.add_constraint(mstruct2=0);
//                Constraint mstruct3("mstruct>"+to_string(j+n)+to_string(i+n));
//                mstruct3 = W(j+n,i+n) - W(i,j);
//                NPP.add_constraint(mstruct3=0);
//            }else{
//                /* w*-w = z */
//                cout << "\neW(" << i << "," << j+n << ") goes to " << n*i+j;
//                Constraint svec("svecI"+to_string(i)+to_string(j));
//                svec = wstar[n*i+j] - W(i,j+n) - z(n*i+j+1);
//                NPP.add_constraint(svec=0);
//
//                /* bounds */
//                Constraint lbound("lbound("+to_string(i)+ "," + to_string(j+n) +")");
//                lbound = W(i,j+n) + 0.605;
//                NPP.add_constraint(lbound >=0);
//                Constraint ubound("ubound("+to_string(i)+ "," + to_string(j+n) +")");
//                ubound = W(i,j+n) - 0.605;
//                NPP.add_constraint(ubound <=0);
//
//                /* matrix structure */
//                Constraint mstruct("mstruct<"+to_string(j)+to_string(i+n));
//                mstruct = W(j,i+n) + W(i,j+n);
//                NPP.add_constraint(mstruct=0);
//                Constraint mstruct2("mstruct<"+to_string(j+n)+to_string(i));
//                mstruct2 = W(j+n,i) - W(i,j+n);
//                NPP.add_constraint(mstruct2=0);
//                Constraint mstruct3("mstruct<"+to_string(i+n)+to_string(j));
//                mstruct3 = W(i+n,j) + W(i,j+n);
//                NPP.add_constraint(mstruct3=0);
//            }
//        }
//    }
//
//    solver s(NPP,mosek_);
//    s.run(1,0);
//    W.print(true);

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

    /* Initialize variables */
    R_Wij.initialize_all(1.0);
    Im_Wij.initialize_all(0.0);
    Wii.initialize_all(1.001);

    /**  Objective */
    func_ obj;
    for (auto g:grid->gens) {
        if (g->_active) {
            obj += grid->c1(g->_name)*Pg(g->_name) + grid->c2(g->_name)*Pg(g->_name)*Pg(g->_name) + grid->c0(g->_name);
        }
    }
    SDP.min(obj);


    /** Constraints */
    /* Second-order cone constraints */
    Constraint SOC("SOC");
    SOC =  power(R_Wij, 2) + power(Im_Wij, 2) - Wii.from()*Wii.to() ;
    SDP.add_constraint(SOC.in(bus_pairs) <= 0);

    /* Flow conservation */
    for (auto b: grid->nodes) {
        Bus* bus = (Bus*) b;
        Constraint KCL_P("KCL_P"+bus->_name);
        Constraint KCL_Q("KCL_Q"+bus->_name);
        KCL_P  = sum(Pf_from.in(b->get_out())) + sum(Pf_to.in(b->get_in())) + bus->pl() - sum(Pg.in(bus->_gen));
        KCL_Q  = sum(Qf_from.in(b->get_out())) + sum(Qf_to.in(b->get_in())) + bus->ql() - sum(Qg.in(bus->_gen));

        /* Shunts */
        KCL_P +=  bus->gs()*(Wii(bus->_name));
        KCL_Q -=  bus->bs()*(Wii(bus->_name));

        SDP.add_constraint(KCL_P = 0);
        SDP.add_constraint(KCL_Q = 0);
    }

    /* AC Power Flow */
    Constraint Flow_P_From("Flow_P_From");
    Flow_P_From += Pf_from;
    Flow_P_From -= grid->g_ff*Wii.from();
    Flow_P_From -= grid->g_ft*R_Wij.in_pairs();
    Flow_P_From -= grid->b_ft*Im_Wij.in_pairs();
    SDP.add_constraint(Flow_P_From.in(grid->arcs) = 0);

    Constraint Flow_P_To("Flow_P_To");
    Flow_P_To += Pf_to;
    Flow_P_To -= grid->g_tt*Wii.to();
    Flow_P_To -= grid->g_tf*R_Wij.in_pairs();
    Flow_P_To += grid->b_tf*Im_Wij.in_pairs();
    SDP.add_constraint(Flow_P_To.in(grid->arcs) = 0);

    Constraint Flow_Q_From("Flow_Q_From");
    Flow_Q_From += Qf_from;
    Flow_Q_From += grid->b_ff*Wii.from();
    Flow_Q_From += grid->b_ft*R_Wij.in_pairs();
    Flow_Q_From -= grid->g_ft*Im_Wij.in_pairs();
    SDP.add_constraint(Flow_Q_From.in(grid->arcs) = 0);

    Constraint Flow_Q_To("Flow_Q_To");
    Flow_Q_To += Qf_to;
    Flow_Q_To += grid->b_tt*Wii.to();
    Flow_Q_To += grid->b_tf*R_Wij.in_pairs();
    Flow_Q_To += grid->g_tf*Im_Wij.in_pairs();
    SDP.add_constraint(Flow_Q_To.in(grid->arcs) = 0);

    /* Phase Angle Bounds constraints */
    Constraint PAD_UB("PAD_UB");
    PAD_UB = Im_Wij;
    PAD_UB -= (grid->tan_th_max)*R_Wij;
    SDP.add_constraint(PAD_UB.in(bus_pairs) <= 0);

    Constraint PAD_LB("PAD_LB");
    PAD_LB =  Im_Wij;
    PAD_LB -= grid->tan_th_min*R_Wij;
    SDP.add_constraint(PAD_LB.in(bus_pairs) >= 0);

    /* Thermal Limit Constraints */
    Constraint Thermal_Limit_from("Thermal_Limit_from");
    Thermal_Limit_from += power(Pf_from, 2) + power(Qf_from, 2);
    Thermal_Limit_from -= power(grid->S_max,2);
    SDP.add_constraint(Thermal_Limit_from.in(grid->arcs) <= 0);

    Constraint Thermal_Limit_to("Thermal_Limit_to");
    Thermal_Limit_to += power(Pf_to, 2) + power(Qf_to, 2);
    Thermal_Limit_to -= power(grid->S_max,2);
    SDP.add_constraint(Thermal_Limit_to.in(grid->arcs) <= 0);


    /* Solver selection */
    solver SDPOPF(SDP,ipopt);
    double solver_time_start = get_wall_time();
    SDPOPF.run(output = 0, relax = false, tol = 1e-6, "ma27", mehrotra = "no");

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

    while((SDP._obj_val - prev_opt)/SDP._obj_val > fp_tol) {
        prev_opt = SDP._obj_val;
        for (auto &b: bags) {
//            b->add_lines();
            if (b->is_PSD()) {
                continue;
            }
            what = b->nfp();
            node_pairs b_pairs;
//            param<> R_Wij_star("R_Wij_star"), R_Wij_hat("I_Wij_hat"),
            param<double> R_diff("R_Diff"), I_diff("I_diff"), W_diff("W_diff");
            param<double> R_hat("R_hat"), I_hat("I_hat"), W_hat("W_hat");
//            param<> Wii_star("Wii_star"), Wii_hat("Wii_hat"), W_diff("W_Diff");
            // the cuts for different dimensions don't have the same form...
            Constraint sdpcut("sdpcut_" + to_string(numcuts));
            Node *ni;
            Arc *aij;
            for (int i = 0; i < b->_nodes.size(); i++) {
                for (int j = i; j < b->_nodes.size(); j++) {
                    if (i == j) {
                        namew = "w(" + b->_nodes[i]->_name + ")";
                        ni = b->_nodes[i];
                        sdpcut += (((Bus *) ni)->w - what(namew)) * (Wii(ni->_name) - what(namew));

                        W_diff.set_val(ni->_name,((Bus *) ni)->w - what(namew).eval());
                        W_hat.set_val(ni->_name,what(namew).eval());
                    } else {
                        aij = grid->get_arc(b->_nodes[i]->_name, b->_nodes[j]->_name);
                        if(aij->_imaginary && !aij->_active) {
                            bus_pairs_sdp._keys.push_back(new index_pair(index_(aij->_src->_name), index_(aij->_dest->_name)));
                            aij->_active = true;
                        }
                        b_pairs._keys.push_back(new index_pair(index_(aij->_src->_name), index_(aij->_dest->_name)));

                        namewr = "wr(" + b->_nodes[i]->_name + "," + b->_nodes[j]->_name + ")";
                        namewi = "wi(" + b->_nodes[i]->_name + "," + b->_nodes[j]->_name + ")";
                        sdpcut += (((Line *) aij)->wr - what(namewr)) *
                                  (R_Wij(aij->_src->_name + "," + aij->_dest->_name) - what(namewr));
                        sdpcut += (((Line *) aij)->wi - what(namewi)) *
                                  (Im_Wij(aij->_src->_name + "," + aij->_dest->_name) - what(namewi));

                        R_diff.set_val(aij->_src->_name + "," + aij->_dest->_name,((Line *) aij)->wr - what(namewr).eval());
                        I_diff.set_val(aij->_src->_name + "," + aij->_dest->_name,((Line *) aij)->wi - what(namewi).eval());
                        R_hat.set_val(aij->_src->_name + "," + aij->_dest->_name,what(namewr).eval());
                        I_hat.set_val(aij->_src->_name + "," + aij->_dest->_name,what(namewi).eval());
                    }
                }
            }
            sdpcut.print();
            SDP.add_constraint(sdpcut <= 0);


            Constraint lin("lin"+to_string(numcuts));
            cout << "\nbpairs size = " << b_pairs._keys.size() << endl;
//            lin = (R_diff.in(b_pairs._keys) + R_Wij.in(b_pairs._keys) - R_hat.in(b_pairs._keys));
            lin = product(R_diff.in(b_pairs._keys),(R_Wij.in(b_pairs._keys) - R_hat.in(b_pairs._keys)));
//            lin += product(I_diff.in(b_pairs._keys),(Im_Wij.in(b_pairs._keys) - I_hat.in(b_pairs._keys)));
//            lin += product(W_diff.in(b->_nodes),(Wii.in(b->_nodes) - W_hat.in(b->_nodes)));
            SDP.add_constraint(lin <= 0);
            lin.print_expanded();

            numcuts++;
        }

//        if(!bus_pairs_sdp._keys.empty()) {
//            Constraint SOC_im("SOC_im");
//            SOC_im = power(R_Wij, 2) + power(Im_Wij, 2) - Wii.from() * Wii.to();
//            SDP.add_constraint(SOC_im.in(bus_pairs_sdp._keys) <= 0);
//            bus_pairs_sdp._keys.clear();
//        }


//        param<> R_Wij_star("R_Wij_star"), R_Wij_hat("I_Wij_hat"), R_diff("R_Diff");
//        param<> Wii_star("Wii_star"), Wii_hat("Wii_hat"), W_diff("W_Diff");

        for(auto& a: grid->arcs){
            if(a->_imaginary && !a->_active) a->_free = true;
        }

        SDPOPF.run(output = 0, relax = false, tol = 1e-6, "ma27", mehrotra = "no");
        DebugOn("\n(opt - prev)/opt = " << (SDP._obj_val - prev_opt)/SDP._obj_val << endl);
    }



    double solver_time_end = get_wall_time();
    double total_time_end = get_wall_time();
    auto solve_time = solver_time_end - solver_time_start;
    auto total_time = total_time_end - total_time_start;
    string out = "DATA_OPF, " + grid->_name + ", " + to_string(nb_buses) + ", " + to_string(nb_lines) +", " + to_string(SDP._obj_val) + ", " + to_string(-numeric_limits<double>::infinity()) + ", " + to_string(solve_time) + ", LocalOptimal, " + to_string(total_time);
    DebugOn(out <<endl);
    return 0;
}
