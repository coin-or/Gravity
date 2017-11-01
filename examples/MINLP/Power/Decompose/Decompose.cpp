//
//  Decompose.cpp
//
//  Gravity
//
//  Created by Guanglei Wang on 6/9/17.
//
//
#include <stdio.h>
#include <string>
#include <stdio.h>
#include <cstring>
#include <fstream>
#include "../PowerNet.h"
#include <gravity/solver.h>
#include <stdio.h>
#include <stdlib.h>
#include "functions.h"

#ifdef USE_BOOST
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <deque>
#include <iterator>
#endif

/** INITIALISE SUBPROBLEM MODEL */
double  subproblem(PowerNet& grid, Net* chordal, unsigned c, Net* cliquetree,
                   vector<Bus*>  bag_bus,  vector<Line*> bag_arcs, vector<Gen*> bag_gens,node_pairs* bag_bus_pairs,
                   vector<Bus*>  bag_bus_disjoint,  vector<Line*> bag_arcs_disjoint,
                   vector<param<Real>>& R_lambda_sep, vector<param<Real>>& Im_lambda_sep, vector<param<Real>>& lambda_sep,
                   param<Real>& R_rho_sep, param<Real>& Im_rho_sep)
{
    cout << "Solving subproblem associated with maximal clique .........." << c << endl;

    if (bag_arcs.size() == 0) {
        return 0;
    }
    Model Subr("Subr");

    // POWER FLOW
    DebugOn("bag_arcs " << c << " has " << bag_arcs.size() << " lines." << endl);
    DebugOn("bag_bus  " << c << " has " << bag_bus.size() << " bus." << endl);
    DebugOn("bag_gens " << c << " has " << bag_gens.size() << " gens." << endl);


    // LIFTED VARIABLES..
    var<Real>  R_Wij("R_Wij", grid.wr_min.in(bag_bus_pairs->_keys), grid.wr_max.in(bag_bus_pairs->_keys));
    var<Real>  Im_Wij("Im_Wij", grid.wi_min.in(bag_bus_pairs->_keys), grid.wi_max.in(bag_bus_pairs->_keys));
    var<Real>  Wii("Wii", grid.w_min.in(bag_bus), grid.w_max.in(bag_bus));
    Subr.add_var(Wii^(bag_bus.size()));
    Subr.add_var(R_Wij^(bag_bus.size()*(bag_bus.size() - 1)/2));
    Subr.add_var(Im_Wij^(bag_bus.size()*(bag_bus.size() - 1)/2));

    R_Wij.initialize_all(1.0);
    Wii.initialize_all(1.001);

    /* Construct the objective function*/
    func_ obj;
    for (auto g: bag_gens) {
        if (g->_active) {
            Constraint Production_P_UB("Production_P_UB" + g->_name);
            Constraint Production_P_LB("Production_P_LB" + g->_name);
            Constraint Production_Q_UB("Production_Q_UB" + g->_name);
            Constraint Production_Q_LB("Production_Q_LB" + g->_name);
            auto bus = g->_bus;
            obj += grid.c1(g->_name).getvalue()*bus->pl() + grid.c0(g->_name).getvalue();
            DebugOff("the constant value is: " << bus->pl()*grid.c1(g->_name).getvalue()+grid.c0(g->_name).getvalue() << endl);
            DebugOff("the constant of this polynomial function is: " << poly_eval(obj.get_cst()) << endl);
            Production_P_UB += bus->gs()*Wii(bus->_name) + bus->pl() - grid.pg_max(g->_name).getvalue();
            Production_P_LB += bus->gs()*Wii(bus->_name) + bus->pl() - grid.pg_min(g->_name).getvalue();
            Production_Q_UB += bus->bs()*Wii(bus->_name) + bus->ql() - grid.qg_max(g->_name).getvalue();
            Production_Q_LB += -bus->bs()*Wii(bus->_name) + bus->ql()- grid.qg_min(g->_name).getvalue();
            //shunt
            obj  += grid.c1(g->_name)*bus->gs()*Wii(bus->_name);
            for (auto &a: bus->get_out()) {
                if (std::find(bag_bus_disjoint.begin(), bag_bus_disjoint.end(),a->_src) != bag_bus_disjoint.end()) {
                    obj  += grid.c1(g->_name)*grid.g_ff(a->_name)*Wii(a->_src->_name);
                    Production_P_UB += grid.g_ff(a->_name)*Wii(a->_src->_name);
                    Production_P_LB += grid.g_ff(a->_name)*Wii(a->_src->_name);
                    Production_Q_UB  += -1*grid.b_ff(a->_name)*Wii(a->_src->_name);
                    Production_Q_LB  += -1*grid.b_ff(a->_name)*Wii(a->_src->_name);
                }
                if (std::find(bag_arcs_disjoint.begin(), bag_arcs_disjoint.end(), a) != bag_arcs_disjoint.end()) {
                    obj += grid.c1(g->_name)*grid.g_ft(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name)
                           + grid.c1(g->_name)*grid.b_ft(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);

                    Production_P_UB +=  grid.g_ft(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name)
                                        + grid.b_ft(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);

                    Production_P_LB += grid.g_ft(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name)
                                       + grid.b_ft(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);

                    Production_Q_UB  +=  -1*grid.b_ft(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name)
                                         + grid.g_ft(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);

                    Production_Q_LB  +=  -1*grid.b_ft(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name)
                                         + grid.g_ft(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);
                }
            }

            for (auto &a: bus->get_in()) {
                if (std::find(bag_bus_disjoint.begin(), bag_bus_disjoint.end(),a->_dest) != bag_bus_disjoint.end()) {
                    obj  += grid.c1(g->_name)*grid.g_ff(a->_name)*Wii(a->_dest->_name);
                    Production_P_UB += grid.g_tt(a->_name)*Wii(a->_dest->_name);
                    Production_P_LB += grid.g_tt(a->_name)*Wii(a->_dest->_name);
                    Production_Q_UB -= grid.b_tt(a->_name)*Wii(a->_dest->_name);
                    Production_Q_LB -= grid.b_tt(a->_name)*Wii(a->_dest->_name);
                }
                if (std::find(bag_arcs_disjoint.begin(), bag_arcs_disjoint.end(),a) != bag_arcs_disjoint.end()) {
                    obj  +=  grid.c1(g->_name)*grid.g_tf(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name)
                             - grid.c1(g->_name)*grid.b_tf(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);
                    Production_P_UB += grid.g_tf(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name)
                                       - grid.b_tf(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);

                    Production_P_LB += grid.g_tf(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name)
                                       - grid.b_tf(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);

                    Production_Q_UB -=  grid.b_tf(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name)
                                        + grid.g_tf(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);

                    Production_Q_LB -= grid.b_tf(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name)
                                       + grid.g_tf(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);
                }
            }
            //Subr.add_constraint(Production_P_UB <= 0);
            //Subr.add_constraint(Production_P_LB >= 0);
            //Subr.add_constraint(Production_Q_UB <= 0);
            //Subr.add_constraint(Production_Q_LB >= 0);
        }
    }

//    for (auto b: bag_bus_disjoint) {
//        if (!b->_has_gen) {
////            obj += 2*R_rho_sep(b->_name)*b->pl() + Im_rho_sep(b->_name)*b->ql();
////            //shunt
////            obj  += 2*R_rho_sep(b->_name)*b->gs()*Wii(b->_name);
////            obj  -= 2*R_rho_sep(b->_name)*b->bs()*Wii(b->_name);
//        }
//    }
//
//    for (auto a: bag_arcs_disjoint) {
//        auto b1 = (Bus*) a->_src;
//        if (!b1->_has_gen) {
////            obj  += 2*R_rho_sep(b1->_name)*grid.g_ff(a->_name)*Wii(a->_src->_name)
////                    + 2*R_rho_sep(b1->_name)*grid.g_ft(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name)
////                    + 2*R_rho_sep(b1->_name)*grid.b_ft(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);
////
////            obj  += -2*R_rho_sep(b1->_name)*grid.b_ff(a->_name)*Wii(a->_src->_name)
////                    - 2*R_rho_sep(b1->_name)*grid.b_ft(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name)
////                    + 2*R_rho_sep(b1->_name)*grid.g_ft(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);
//        }
//        auto b2 = (Bus*) a->_dest;
//        if (!b2->_has_gen) {
//            obj +=  2*R_rho_sep(b2->_name)*grid.g_tt(a->_name)*Wii(a->_dest->_name)
//                    +2*R_rho_sep(b2->_name)*grid.g_tf(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name)
//                    -2*R_rho_sep(b2->_name)* grid.b_tf(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);
//
//            obj -=  2*R_rho_sep(b2->_name)*grid.b_tt(a->_name)*Wii(a->_dest->_name)
//                    +2*R_rho_sep(b2->_name)*grid.b_tf(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name)
//                    +2*R_rho_sep(b2->_name)* grid.g_tf(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);
//        }
//    }
    Node* Cr = cliquetree->get_node(to_string(c));
    Arc* arc = nullptr;
    Bus* bus = nullptr;
    for (auto a: Cr->get_out()) {
        Debug("a->_intersection.size " << a->_intersection.size() << endl);
        for (int i = 0; i < a->_intersection.size(); i ++) {
            bus = (Bus*)a->_intersection.at(i);
            obj += lambda_sep[a->_id](bus->_name)*Wii(bus->_name);
            string name;
            for (int j = i + 1; j < a->_intersection.size(); j ++) {
                arc = chordal->get_arc(bus, a->_intersection.at(j));
                name = arc->_src->_name + "," + arc->_dest->_name;
                obj += R_lambda_sep[a->_id](name)*R_Wij(name);
                obj += Im_lambda_sep[a->_id](name)*Im_Wij(name);
            }
        }
    }

    for (auto a: Cr->get_in()) {
        for (int i = 0; i < a->_intersection.size(); i ++) {
            bus = (Bus*) a->_intersection.at(i);
            obj -= lambda_sep[a->_id](bus->_name)*Wii(bus->_name);
            string name;
            for (int j = i + 1; j < a->_intersection.size(); j ++) {
                arc = chordal->get_arc(bus, a->_intersection.at(j));
                name = arc->_src->_name + "," + arc->_dest->_name;
                obj -= R_lambda_sep[a->_id](name)*R_Wij(name);
                obj -= Im_lambda_sep[a->_id](name)*Im_Wij(name);
            }
        }
    }

    Subr.set_objective(min(obj));

    /* Subr constraints */
    Constraint SOC("SOC");
    SOC =  power(R_Wij.in(bag_bus_pairs->_keys), 2) + power(Im_Wij.in(bag_bus_pairs->_keys), 2)
           - Wii.from(bag_bus_pairs->_keys)*Wii.to(bag_bus_pairs->_keys);
    Subr.add_constraint(SOC <= 0);

    //* Phase Angle Bounds constraints */
    Constraint PAD_UB("PAD_UB");
    PAD_UB = Im_Wij.in(bag_bus_pairs->_keys);
    PAD_UB -= (grid.tan_th_max).in(bag_bus_pairs->_keys)*R_Wij.in(bag_bus_pairs->_keys);
    Subr.add_constraint(PAD_UB <= 0);

    Constraint PAD_LB("PAD_LB");
    PAD_LB =  Im_Wij.in(bag_bus_pairs->_keys);
    PAD_LB -= grid.tan_th_min.in(bag_bus_pairs->_keys)*R_Wij.in(bag_bus_pairs->_keys);
    Subr.add_constraint(PAD_LB >= 0);

    /* Thermal Limit Constraints */
    Constraint Thermal_Limit_from("Thermal_Limit_from");
    Thermal_Limit_from += power(grid.g_ff.in(bag_arcs)*Wii.from(bag_arcs)
                                + grid.g_ft.in(bag_arcs)*R_Wij.in_pairs(bag_arcs)+ grid.b_ft.in(bag_arcs)*Im_Wij.in_pairs(bag_arcs), 2)
                          + power(grid.g_ft.in(bag_arcs)*Im_Wij.in_pairs(bag_arcs)-grid.b_ff.in(bag_arcs)*Wii.to(bag_arcs)
                                  - grid.b_ft.in(bag_arcs)*R_Wij.in_pairs(bag_arcs), 2);
    Thermal_Limit_from -= power(grid.S_max.in(bag_arcs),2);
    Subr.add_constraint(Thermal_Limit_from <= 0);

    Constraint Thermal_Limit_to("Thermal_Limit_to");
    Thermal_Limit_to += power(grid.g_tt.in(bag_arcs)*Wii.from(bag_arcs) + grid.g_tf.in(bag_arcs)*R_Wij.in_pairs(bag_arcs)
                              + grid.b_tf.in(bag_arcs)*Im_Wij.in_pairs(bag_arcs), 2)
                        + power(grid.b_tt.in(bag_arcs)*Wii.to(bag_arcs) + grid.b_tf.in(bag_arcs)*R_Wij.in_pairs(bag_arcs)
                                + grid.g_tf.in(bag_arcs)*Im_Wij.in_pairs(bag_arcs), 2);

    Thermal_Limit_to -= power(grid.S_max.in(bag_arcs),2);
    Subr.add_constraint(Thermal_Limit_to <= 0);


    /* solve it! */
    solver solve_Subr(Subr,cplex);
    solve_Subr.run();

    return Subr._obj_val;
}

int main (int argc, const char * argv[])
{
    // Decompose
    const char* fname;
    if (argc >= 2) {
        fname = argv[1];
    }
    else {
        //            fname = "../../data_sets/Power/nesta_case5_pjm.m";
        //fname = "../../data_sets/Power/nesta_case14_ieee.m";
        //fname = "../../data_sets/Power/nesta_case9241_pegase.m";
        //fname = "../../data_sets/Power/nesta_case2383wp_mp.m";
        //fname = "../../data_sets/Power/nesta_case1354_pegase_api.m";
        //fname = "../../data_sets/Power/nesta_case118_ieee.m";
        fname = "/Users/hh/Dropbox/Work/Dev/pglib-opf/pglib_opf_case6495_rte.m";
        //        fname = "/Users/hh/Dropbox/Work/Dev/nesta-0.7.0/opf/nesta_case3_lmbd.m";
        //        fname = "/Users/hh/Dropbox/Work/Dev/nesta-0.7.0/opf/nesta_case5_pjm.m";
    }
    // ACOPF
    PowerNet grid;
    //    fname = "../../data_sets/Power/nesta_case3_lmbd.m";
    //    fname = "../../data_sets/Power/nesta_case14_ieee.m";
    //    fname = "../../data_sets/Power/nesta_case9241_pegase.m";
    //    fname = "/Users/hh/Dropbox/Work/Dev/nesta-0.7.0/opf/nesta_case3375wp_mp.m";
    //    fname = "../../data_sets/Power/nesta_case300_ieee.m";
    //    fname = "../../data_sets/Power/nesta_case2383wp_mp.m";
    grid.readgrid(fname);
    
//    OPF_Clique_W(grid);
    scopf_W(grid, true);
    return 0;
    
    // Grid Parameters
    unsigned nb_gen = grid.get_nb_active_gens();
    unsigned nb_lines = grid.get_nb_active_arcs();
    unsigned nb_buses = grid.get_nb_active_nodes();

    DebugOn("nb gens = " << nb_gen << endl);
    DebugOn("nb lines = " << 2*nb_lines << endl);
    DebugOn("nb buses = " << nb_buses << endl);

    /** Clique tree decomposition **/
    Net* chordal = grid.get_chordal_extension();
    grid.get_clique_tree();
    const unsigned nb_cliques = grid._bags.size();

    vector<vector<Bus*>> bag_bus; // each clique contains just nodes, not buses! Fixed this by modifying the bag definition.
    vector<vector<Bus*>> bag_bus_disjoint;
    vector<vector<Gen*>> bag_gens;
    vector<vector<Line*>> bag_arcs; //bag_arcs contains the arcs of the power grid while variables associated with W are defined on chordal graph.
    vector<vector<Line*>> bag_arcs_disjoint;
    vector<gravity::node_pairs*> bag_bus_pairs; // bus_pairs in each bag.
    map<Line*, unsigned> indexij;
    map<Bus*, unsigned> indexii;

    for (int c = 0; c < nb_cliques; c++) {
        vector<Bus*> bag_B;
        vector<Bus*> bag_B_disjoint;
        vector<Gen*> bag_G;
        vector<Line*> bag_A;
        vector<Line*> bag_A_disjoint;
        node_pairs* bag_BP = new node_pairs();
        for (int i = 0; i < grid._bags[c].size(); i++) {
            Bus* b = (Bus*) grid.get_node(grid._bags[c].at(i)->_name);
            if (b !=nullptr) {
                bag_B.push_back(b);
                auto pp = indexii.insert(make_pair<>(b, c));
                if (pp.second) {
                    bag_B_disjoint.push_back(b);
                }
            }
            if (b->_has_gen) {
                bag_G.insert(bag_G.end(), b->_gen.begin(), b->_gen.end());
            }
            for (int j = i+1; j < grid._bags[c].size(); j++) {
                Line* a = (Line*)grid.get_arc(b, grid.get_node(grid._bags[c].at(j)->_name));
                if (a != nullptr) {
                    bag_A.push_back(a);
                    bag_BP->_keys.push_back(new index_pair(index_(a->_src->_name), index_(a->_dest->_name), a->_active));
                    auto pp= indexij.insert(make_pair<>(a, c));
                    if (pp.second) {
                        bag_A_disjoint.push_back(a);
                    }
                }
            }
        }
        bag_bus.push_back(bag_B);
        bag_gens.push_back(bag_G);
        bag_arcs.push_back(bag_A);
        bag_bus_pairs.push_back(bag_BP);
        bag_bus_disjoint.push_back(bag_B_disjoint);
        bag_arcs_disjoint.push_back(bag_A_disjoint);
    }

    auto cliquetree = grid.get_clique_tree();
///////////////// DEFINE LAGRANGE MULTIPLIERS  ////////////////////////////////
    vector<param<Real>> R_lambda_in;
    vector<param<Real>> Im_lambda_in;
    vector<param<Real>> lambda_in;
    param<Real> R_rho_in("R_rho_in");
    param<Real> Im_rho_in("Im_rho_in");
    R_rho_in^nb_buses;
    Im_rho_in^nb_buses;
    R_rho_in.initialize_all(0);
    Im_rho_in.initialize_all(0);


    vector<param<Real>> R_lambda_out;
    vector<param<Real>> Im_lambda_out;
    vector<param<Real>> lambda_out;
    param<Real> R_rho_out("R_rho_out");
    param<Real> Im_rho_out("Im_rho_out");
    R_rho_out^nb_buses;
    Im_rho_out^nb_buses;
    R_rho_out.initialize_all(0);
    Im_rho_out.initialize_all(0);

    vector<param<Real>> R_lambda_sep;
    vector<param<Real>> Im_lambda_sep;
    vector<param<Real>> lambda_sep;
    param<Real> R_rho_sep("R_rho_sep");
    param<Real> Im_rho_sep("Im_rho_sep");
    R_rho_sep^nb_buses;
    Im_rho_sep^nb_buses;
    R_rho_sep.initialize_all(0);
    Im_rho_sep.initialize_all(0);

    vector<param<Real>> R_lambda_grad;
    vector<param<Real>> Im_lambda_grad;
    vector<param<Real>> lambda_grad;
    param<Real> R_rho_grad("R_rho_grad");
    param<Real> Im_rho_grad("R_rho_grad");
    R_rho_grad^nb_buses;
    Im_rho_grad^nb_buses;
    R_rho_grad.initialize_all(0);
    Im_rho_grad.initialize_all(0);

    for (auto a: cliquetree->arcs) {
        auto l = a->_id;
        param<Real> R_lambda_arc_in("R_lambda_arc_in" + to_string(l));
        param<Real> Im_lambda_arc_in("Im_lambda_arc_in" + to_string(l));
        param<Real> lambda_arc_in("R_lambda_arc_in" + to_string(l));
        R_lambda_arc_in^(a->_weight*(a->_weight - 1)/2);
        Im_lambda_arc_in^(a->_weight*(a->_weight - 1)/2);
        lambda_arc_in^(a->_weight*(a->_weight - 1)/2);

        param<Real> R_lambda_arc_out("R_lambda_arc_out" + to_string(l));
        param<Real> Im_lambda_arc_out("Im_lambda_arc_out" + to_string(l));
        param<Real> lambda_arc_out("lambda_arc_out" + to_string(l));
        R_lambda_arc_out^(a->_weight*(a->_weight-1)/2);
        Im_lambda_arc_out^(a->_weight*(a->_weight-1)/2);
        lambda_arc_out^(a->_weight*(a->_weight-1)/2);

        param<Real> R_lambda_arc_sep("R_lambda_arc_sep" + to_string(l));
        param<Real> Im_lambda_arc_sep("Im_lambda_arc_sep" + to_string(l));
        param<Real> lambda_arc_sep("lambda_arc_sep" + to_string(l));
        R_lambda_arc_sep^(a->_weight*(a->_weight-1)/2);
        Im_lambda_arc_sep^(a->_weight*(a->_weight-1)/2);
        lambda_arc_sep^(a->_weight*(a->_weight-1)/2);

        param<Real> R_lambda_arc_grad("R_lambda_arc_grad" + to_string(l));
        param<Real> Im_lambda_arc_grad("Im_lambda_arc_grad" + to_string(l));
        param<Real> lambda_arc_grad("lambda_arc_grad" + to_string(l));

        for (int i = 0; i < a->_weight; i++) {
            lambda_arc_sep(i) = 0;
        }

        R_lambda_in.push_back(R_lambda_arc_in);
        Im_lambda_in.push_back(Im_lambda_arc_in);
        lambda_in.push_back(lambda_arc_in);

        R_lambda_out.push_back(R_lambda_arc_out);
        Im_lambda_out.push_back(Im_lambda_arc_out);
        lambda_out.push_back(lambda_arc_out);

        R_lambda_sep.push_back(R_lambda_arc_sep);
        Im_lambda_sep.push_back(Im_lambda_arc_sep);
        lambda_sep.push_back(lambda_arc_sep);

        R_lambda_grad.push_back(R_lambda_arc_grad);
        Im_lambda_grad.push_back(Im_lambda_arc_grad);
        lambda_grad.push_back(lambda_arc_grad);
    }

    Arc* arc = nullptr;
    Bus* bus = nullptr;
    for (auto a: cliquetree->arcs) {
        int l = a->_id;
        for (int i = 0; i < a->_intersection.size(); i ++) {
            bus = (Bus*) a->_intersection.at(i);
            for (int j = i + 1; j < a->_intersection.size(); j ++) {
                arc = chordal->get_arc(bus, a->_intersection.at(j)); // we have to use the arc name.
                R_lambda_sep[l](arc->_name)  =  0;
                R_lambda_in[l](arc->_name)  =  0;
                R_lambda_out[l](arc->_name)  =  0;
                Im_lambda_sep[l](arc->_name)  =  0;
                Im_lambda_in[l](arc->_name)  =  0;
                Im_lambda_out[l](arc->_name)  =  0;
            }

            for (auto g: bus->_gen) {
                if (g->_active) {
                    lambda_sep[l](g->_name) = 0;
                    lambda_out[l](g->_name) = 0;
                    lambda_in[l] (g->_name) = 0;
                }
            }
        }
    }

///////////////////////////////// INITIALIZATION ///////////////////////////////////////////
    vector<double> value_dual;
    double val = 0;
    double dual = 0;
    for (int c = 0; c < nb_cliques; c++) {
        val = subproblem(grid, chordal, c, cliquetree,
                         bag_bus[c],  bag_arcs[c],bag_gens[c], bag_bus_pairs[c],
                         bag_bus_disjoint[c], bag_arcs_disjoint[c], R_lambda_sep,
                         Im_lambda_sep, lambda_sep, R_rho_sep, Im_rho_sep);
        value_dual.push_back(val);
        dual +=val;
    }
    value_dual.resize(nb_cliques);
    cout << "................  Initialization value:  " << dual <<endl;
    return 0;
}
