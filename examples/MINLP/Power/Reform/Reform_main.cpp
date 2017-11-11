//
//  Reform.cpp
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

#ifdef USE_BOOST
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <deque>
#include <iterator>
#endif


void reform1(PowerNet& grid){
    Net* grid_augment = grid.clone();
    // making each node simplicial.
    
    for (auto &node: grid.nodes){
        vector<Node*> neighbours = node->get_neighbours();
        for (unsigned i = 0; i < neighbours.size(); i ++){
            auto n1 = grid_augment->get_node(neighbours.at(i)->_name);
            for ( unsigned j = i+1; j < neighbours.size(); j++){
                auto n2 = grid_augment->get_node(neighbours.at(j)->_name);
                string name = to_string((int)grid_augment->arcs.size()+1);
                Arc* arc = new Arc(name);
                arc->_id = grid_augment->arcs.size();
                arc->_src = n1;
                arc->_dest = n2;
                arc->connect();
                grid_augment->add_undirected_arc(arc);
            }
        }
    }
    grid_augment->get_tree_decomp_bags();
    auto cliquetree = grid_augment->get_clique_tree();
    const unsigned nb_cliques = grid_augment->_bags.size();
}




int main (int argc, const char * argv[])
{
    // Reform
    const char* fname;
    if (argc >= 2) {
        fname = argv[1];
    }
    else {
         //fname = "../../data_sets/Power/nesta_case5_pjm.m";
         //fname = "../../data_sets/Power/nesta_case14_ieee.m";
        // fname = "../../data_sets/Power/nesta_case30_ieee.m";
        // fname = "../../data_sets/Power/nesta_case6_c.m";
        // fname = "../../data_sets/Power/nesta_case5_pjm.m";
        // fname = "../../data_sets/Power/nesta_case3_lmbd.m";

        //fname = "../../data_sets/Power/nesta_case9241_pegase.m";
        //fname = "../../data_sets/Power/nesta_case2383wp_mp.m";
        //fname = "../../data_sets/Power/nesta_case1354_pegase_api.m";
        fname = "../../data_sets/Power/nesta_case300_ieee.m";
         // fname = "../../data_sets/Power/nesta_case118_ieee.m";
        //fname = "/Users/hh/Dropbox/Work/Dev/pglib-opf/pglib_opf_case6495_rte.m";
        //        fname = "/Users/hh/Dropbox/Work/Dev/nesta-0.7.0/opf/nesta_case3_lmbd.m";
        //        fname = "/Users/hh/Dropbox/Work/Dev/nesta-0.7.0/opf/nesta_case5_pjm.m";
    }
    // ACOPF
    PowerNet grid;
    grid.readgrid(fname);
    cout << "////////////////////////////////////////" << endl;
    
    reform1(grid);
    return 0;
    /** Clique tree decomposition **/
    //grid.chol_decompose(true);
    grid.get_chordal_extension();
    auto cliquetree = grid.get_clique_tree();
    const unsigned nb_cliques = grid._bags.size();


    vector<vector<Bus*>> bag_bus;
    vector<vector<Bus*>> bag_bus_disjoint;
    vector<vector<Gen*>> bag_gens;
    vector<vector<Gen*>> bag_gens_disjoint;
    vector<vector<Line*>> bag_arcs; //bag_arcs contains the arcs of the power grid while variables associated with W are defined on chordal graph.
    vector<vector<Line*>> bag_arcs_disjoint;
    vector<gravity::node_pairs*> bag_bus_pairs; // bus_pairs in each bag.
    vector<gravity::node_pairs*> bag_bus_pairs_disjoint; // bus_pairs in each bag.

    map<Line*, unsigned> indexij;
    map<Bus*, unsigned> indexii;

    for (int c = 0; c < nb_cliques; c++) {
        vector<Bus*> bag_B;
        vector<Bus*> bag_B_disjoint;
        vector<Gen*> bag_G;
        vector<Gen*> bag_G_disjoint;
        vector<Line*> bag_A;
        vector<Line*> bag_A_disjoint;
        node_pairs* bag_BP = new node_pairs();
        node_pairs* bag_BP_disjoint = new node_pairs();
        for (int i = 0; i < grid._bags[c].size(); i++) {
            Bus* b = (Bus*) grid.get_node(grid._bags[c].at(i)->_name);
            if (b !=nullptr) {
                DebugOn(b->_name << " ");
                bag_B.push_back(b);
                if (b->_has_gen) {
                    bag_G.insert(bag_G.end(), b->_gen.begin(), b->_gen.end());
                }
                
                auto pp = indexii.insert(make_pair<>(b, c));
               
                if (pp.second) {
                    bag_B_disjoint.push_back(b);
                    if (b->_has_gen) {
                        bag_G_disjoint.insert(bag_G_disjoint.end(), b->_gen.begin(), b->_gen.end());
                    }
                }
                
            }
            for (int j = i+1; j < grid._bags[c].size(); j++) {
                Line* a = (Line*)grid.get_arc(b, grid.get_node(grid._bags[c].at(j)->_name));
                if (a != nullptr) {
                    bag_A.push_back(a);
                    bag_BP->_keys.push_back(new index_pair(index_(a->_src->_name), index_(a->_dest->_name), a->_active));
                    auto pp= indexij.insert(make_pair<>(a, c));
                    if (pp.second) {
                        bag_A_disjoint.push_back(a);
                        bag_BP_disjoint->_keys.push_back(new index_pair(index_(a->_src->_name), index_(a->_dest->_name), a->_active));

                    }
                }
            }
        }

        DebugOn("] " << endl);
        bag_bus.push_back(bag_B);
        bag_gens.push_back(bag_G);
        bag_arcs.push_back(bag_A);
        bag_bus_pairs.push_back(bag_BP);
        bag_bus_pairs_disjoint.push_back(bag_BP_disjoint);
        bag_bus_disjoint.push_back(bag_B_disjoint);
        bag_gens_disjoint.push_back(bag_G_disjoint);
        bag_arcs_disjoint.push_back(bag_A_disjoint);
    }

    /** build model */
    Model CLT("Clique tree based Model");

    /** Variables */
    vector<var<Real>> R_Wij;
    vector<var<Real>> Im_Wij;
    vector<var<Real>> Wii;
    vector<var<Real>> Pg;
    vector<var<Real>> Qg;
    for (int c = 0; c < nb_cliques; c++) {
        var<Real>  bag_R_Wij("R_Wij_" + to_string(c), grid.wr_min.in(bag_bus_pairs[c]->_keys), grid.wr_max.in(bag_bus_pairs[c]->_keys)); // real part of Wij
        var<Real>  bag_Im_Wij("Im_Wij_" + to_string(c), grid.wi_min.in(bag_bus_pairs[c]->_keys), grid.wi_max.in(bag_bus_pairs[c]->_keys)); // imaginary part of Wij.
        var<Real>  bag_Wii("Wii_" + to_string(c), grid.w_min.in(bag_bus[c]), grid.w_max.in(bag_bus[c]));
        var<Real>  bag_Pg("Pg_" + to_string(c), grid.pg_min.in(bag_gens[c]), grid.pg_max.in(bag_gens[c]));
        var<Real>  bag_Qg ("Qg_" + to_string(c), grid.qg_min.in(bag_gens[c]), grid.qg_max.in(bag_gens[c]));
        CLT.add_var(bag_Wii^(bag_bus[c].size()));
        CLT.add_var(bag_R_Wij^(bag_bus[c].size()*(bag_bus[c].size()-1)/2));
        CLT.add_var(bag_Im_Wij^(bag_bus[c].size()*(bag_bus[c].size()-1)/2));
        CLT.add_var(bag_Pg^(bag_gens[c].size()));
        CLT.add_var(bag_Qg^(bag_gens[c].size()));

        bag_R_Wij.initialize_all(1.0);
        bag_Wii.initialize_all(1.001);
        R_Wij.push_back(bag_R_Wij);
        Im_Wij.push_back(bag_Im_Wij);
        Wii.push_back(bag_Wii);
        Pg.push_back(bag_Pg);
        Qg.push_back(bag_Qg);
    }

    /* Construct the objective function with generations bound constraints */
    func_ obj;
    for (int c = 0; c < nb_cliques; c++) {
        for (auto g:bag_gens_disjoint[c]) {
            if (g->_active) {
                obj += grid.c1(g->_name)*Pg[c](g->_name) + grid.c0(g->_name) + grid.c2(g->_name)*Pg[c](g->_name)*Pg[c](g->_name);
            }
        }
    }

    CLT.set_objective(min(obj));

    /** Define constraints */
    /* CLT constraints */
//    for (int c = 0; c < nb_cliques; c++) {
//        Constraint SOC("SOC_" + to_string(c));
//        if (bag_bus_pairs[c]->_keys.size() > 0){
//            SOC =  power(R_Wij[c].in(bag_bus_pairs[c]->_keys), 2)
//               + power(Im_Wij[c].in(bag_bus_pairs[c]->_keys), 2)
//               - Wii[c].from(bag_bus_pairs[c]->_keys)*Wii[c].to(bag_bus_pairs[c]->_keys) ;
//            CLT.add_constraint(SOC <= 0);
//        }
//    }


    for (int c = 0; c < nb_cliques; c++) {
        //* Phase Angle Bounds constraints */
        if (bag_bus_pairs_disjoint[c]->_keys.size() > 0) {
//            Constraint PAD_UB("PAD_UB" + to_string(c));
//            PAD_UB = Im_Wij[c].in(bag_bus_pairs[c]->_keys);
//            PAD_UB -= (grid.tan_th_max).in(bag_bus_pairs[c]->_keys)*R_Wij[c].in(bag_bus_pairs[c]->_keys);
//            CLT.add_constraint(PAD_UB <= 0);
//
//            Constraint PAD_LB("PAD_LB" + to_string(c));
//            PAD_LB = Im_Wij[c].in(bag_bus_pairs[c]->_keys);
//            PAD_LB -= (grid.tan_th_min).in(bag_bus_pairs[c]->_keys)*R_Wij[c].in(bag_bus_pairs[c]->_keys);
//            CLT.add_constraint(PAD_LB >= 0);

            /* Thermal Limit Constraints */
            Constraint Thermal_Limit_from("Thermal_Limit_from" + to_string(c));
            Thermal_Limit_from += power(grid.g_ff.in(bag_arcs[c])*Wii[c].from(bag_arcs[c])+
                                        grid.g_ft.in(bag_arcs[c])*R_Wij[c].in_pairs(bag_arcs[c])
                                        + grid.b_ft.in(bag_arcs[c])*Im_Wij[c].in_pairs(bag_arcs[c]), 2)
                                + power(grid.g_ft.in(bag_arcs[c])*Im_Wij[c].in_pairs(bag_arcs[c])-
                                        grid.b_ff.in(bag_arcs[c])*Wii[c].from(bag_arcs[c])
                                          - grid.b_ft.in(bag_arcs[c])*R_Wij[c].in_pairs(bag_arcs[c]), 2);
            Thermal_Limit_from -= power(grid.S_max.in(bag_arcs[c]),2);
            CLT.add_constraint(Thermal_Limit_from <= 0);

            Constraint Thermal_Limit_to("Thermal_Limit_to" + to_string(c));
            Thermal_Limit_to += power(grid.g_tt.in(bag_arcs[c])*Wii[c].to(bag_arcs[c]) + grid.g_tf.in(bag_arcs[c])*R_Wij[c].in_pairs(bag_arcs[c])
                                      + grid.b_tf.in(bag_arcs[c])*Im_Wij[c].in_pairs(bag_arcs[c]), 2)
                                + power(grid.b_tt.in(bag_arcs[c])*Wii[c].to(bag_arcs[c]) + grid.b_tf.in(bag_arcs[c])*R_Wij[c].in_pairs(bag_arcs[c])
                                        + grid.g_tf.in(bag_arcs[c])*Im_Wij[c].in_pairs(bag_arcs[c]), 2);

            Thermal_Limit_to -= power(grid.S_max.in(bag_arcs[c]), 2);
            CLT.add_constraint(Thermal_Limit_to <= 0);
        }
    }

    //// COUPLING CONSTRAINTS
    for (auto b: grid.nodes) {
        Bus* bus = (Bus*) b;
        Constraint KCL_P("KCL_P"+bus->_name);
        Constraint KCL_Q("KCL_Q"+bus->_name);
        for (int c = 0; c < nb_cliques; c++) {
            if (std::find(bag_bus_disjoint[c].begin(), bag_bus_disjoint[c].end(), bus) != bag_bus_disjoint[c].end()) {
                KCL_P += bus->pl()- sum(Pg[c].in(bus->_gen))+ bus->gs()*(Wii[c](bus->_name));
                KCL_Q += bus->ql()- sum(Qg[c].in(bus->_gen))- bus->bs()*(Wii[c](bus->_name));
            }
            for (auto &a: b->get_out()) {
                if (std::find(bag_arcs_disjoint[c].begin(), bag_arcs_disjoint[c].end(), a) != bag_arcs_disjoint[c].end()) {
                    KCL_P  += grid.g_ff(a->_name)*Wii[c](a->_src->_name)
                              +grid.g_ft(a->_name)*R_Wij[c](a->_src->_name+","+a->_dest->_name)
                              +grid.b_ft(a->_name)*Im_Wij[c](a->_src->_name+","+a->_dest->_name);

                    KCL_Q  += -1*grid.b_ff(a->_name)*Wii[c](a->_src->_name)
                              - grid.b_ft(a->_name)*R_Wij[c](a->_src->_name+","+a->_dest->_name)
                              +grid.g_ft(a->_name)*Im_Wij[c](a->_src->_name+","+a->_dest->_name);
                }
            }

            for (auto &a: b->get_in()) {
                if (std::find(bag_arcs_disjoint[c].begin(), bag_arcs_disjoint[c].end(), a) != bag_arcs_disjoint[c].end()) {
                    KCL_P  += grid.g_tt(a->_name)*Wii[c](a->_dest->_name)
                              +grid.g_tf(a->_name)*R_Wij[c](a->_src->_name+","+a->_dest->_name)-grid.b_tf(a->_name)*Im_Wij[c](a->_src->_name+","+a->_dest->_name);

                    KCL_Q  -= grid.b_tt(a->_name)*Wii[c](a->_dest->_name)
                              + grid.b_tf(a->_name)*R_Wij[c](a->_src->_name+","+a->_dest->_name)
                              + grid.g_tf(a->_name)*Im_Wij[c](a->_src->_name+","+a->_dest->_name);
                }
            }
        }
        CLT.add_constraint(KCL_P = 0);
        CLT.add_constraint(KCL_Q = 0);
    }


    for (auto a: cliquetree->arcs) {
        Constraint Link_Wii("Link_Wii_" + to_string(a->_id));
        Link_Wii += Wii[a->_src->_id].in(a->_intersection);
        Link_Wii -= Wii[a->_dest->_id].in(a->_intersection);
        CLT.add_constraint(Link_Wii = 0);

        if (a->_intersection_clique.size()>0) {
            Constraint Link_Im_Wij("Link_Im_Wij_" + to_string(a->_id));
            Link_Im_Wij += Im_Wij[a->_src->_id].in(a->_intersection_clique);
            Link_Im_Wij -= Im_Wij[a->_dest->_id].in(a->_intersection_clique);
            CLT.add_constraint(Link_Im_Wij = 0);

            Constraint Link_R_Wij("Link_R_Wij_" + to_string(a->_id));
            Link_R_Wij += R_Wij[a->_src->_id].in(a->_intersection_clique);
            Link_R_Wij -= R_Wij[a->_dest->_id].in(a->_intersection_clique);
            CLT.add_constraint(Link_R_Wij = 0);
        }
    }
    //solver SCOPF(CLT, cplex);
    solver SCOPF(CLT, ipopt);
    SCOPF.run();

    return 0;
}
