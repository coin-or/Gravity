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

bool node_id_compare(const Node* n1, const Node* n2) {
    return n1->_id < n2->_id;
}

bool arc_id_compare(const Arc* n1, const Arc* n2) {
    return n1->_id < n2->_id;
}

bool bus_pair_compare(const index_pair* n1, const index_pair* n2) {
    return n1->_name < n2->_name;
}

void LP(PowerNet& grid) {
    //Grid parameters
    auto bus_pairs = grid.get_bus_pairs();
    Net* grid_augment = grid.clone_undirected();
    DebugOn("number of edges " << grid_augment->arcs.size() << endl);

    for (auto &node: grid.nodes) {
        vector<Node*> neighbours = node->get_neighbours();
        for (int i = 0; i < neighbours.size()-1; i++) {
            auto n1 = grid_augment->get_node(neighbours.at(i)->_name);
            for (int j = i+1; j < neighbours.size(); j++) {
                auto n2 = grid_augment->get_node(neighbours.at(j)->_name);
                if (grid_augment->get_undirected_arc(n1, n2) !=nullptr) {
                    continue;
                }
                else {
                    string name = to_string((int)grid_augment->arcs.size()) + "," +n1->_name + "," + n2->_name;
                    Arc* arc = new Arc(name);
                    arc->_id = grid_augment->arcs.size();
                    arc->_src = n1;
                    arc->_dest = n2;
                    arc->connect();
                    grid_augment->add_undirected_arc(arc);
                }
            }
        }
    }
// note that grid_augment may NOT be chordal.
    grid_augment->get_tree_decomp_bags();
    auto cliquetree = grid_augment->get_clique_tree_kruskal(); // Note: it builds on grid_augment.
    const unsigned nb_cliques = grid_augment->_bags.size();

    vector<vector<Bus*>> bag_bus;
    vector<vector<Bus*>> bag_bus_disjoint;
    vector<vector<Gen*>> bag_gens;
    vector<vector<Gen*>> bag_gens_disjoint;
    vector<vector<Line*>> bag_arcs; //bag_arcs contains the arcs of the power grid while variables associated with W are defined on chordal graph.
    vector<vector<Line*>> bag_arcs_disjoint;
    vector<gravity::node_pairs*> bag_bus_pairs; // bus_pairs in each bag.
    vector<gravity::node_pairs*> bag_bus_pairs_disjoint; // bus_pairs in each bag.

    map<string, unsigned> indexij;
    map<Bus*, unsigned> indexii;

    for (int c = 0; c < nb_cliques; c++) {
        vector<Bus*> bag_B;
        vector<Gen*> bag_G;
        DebugOff("bag " << c << ": ");
        for (int i = 0; i < grid_augment->_bags[c].size(); i++) {
            Node* n = grid_augment->_bags[c].at(i);
            Bus* b = (Bus*) grid.get_node(n->_name);
            DebugOff(b->_name << " ");
            bag_B.push_back(b);
            if (b->_has_gen) {
                bag_G.insert(bag_G.end(), b->_gen.begin(), b->_gen.end());
            }
        }
        bag_bus.push_back(bag_B);
        bag_gens.push_back(bag_G);
    }

    for (int c = 0; c < nb_cliques; c++) {
        vector<Bus*> bag_B_disjoint;
        vector<Gen*> bag_G_disjoint;

        vector<Line*> bag_A;
        vector<Line*> bag_A_disjoint;

        node_pairs* bag_BP = new node_pairs();
        node_pairs* bag_BP_disjoint = new node_pairs();

        sort(bag_bus[c].begin(), bag_bus[c].end(),node_id_compare);
        for (int i = 0; i < bag_bus[c].size(); i++) {
            Bus* b = bag_bus[c].at(i);
            vector<Node*> N = b->get_neighbours();
            sort(N.begin(), N.end(),node_id_compare);

            if (indexii.find(b)==indexii.end()) {
                bool inclusion = std::includes(bag_bus[c].begin(), bag_bus[c].end(), N.begin(),N.end(),node_id_compare);
                if (inclusion) {
                    indexii.insert(make_pair<>(b,c));
                    bag_B_disjoint.push_back(b);
                    if (b->_has_gen) {
                        bag_G_disjoint.insert(bag_G_disjoint.end(), b->_gen.begin(), b->_gen.end());
                    }
                }
            }
            for (int j = i+1; j < bag_bus[c].size(); j++) {
                std::vector<Arc*> vec_arcs1 = grid.get_arcs(b, bag_bus[c].at(j));
                std::vector<Arc*> vec_arcs2 = grid.get_arcs(bag_bus[c].at(j), b);
                if (vec_arcs1.size() +vec_arcs2.size() > 0) {
                    if (vec_arcs1.size() > 0) {
                        bag_BP->_keys.push_back(new index_pair(index_(b->_name), index_(bag_bus[c].at(j)->_name), true));
                    }
                    else {
                        bag_BP->_keys.push_back(new index_pair(index_(bag_bus[c].at(j)->_name),index_(b->_name), true));
                    }

                    for (auto a: vec_arcs1)
                        bag_A.push_back((Line*)a);
                    for (auto a: vec_arcs2)
                        bag_A.push_back((Line*)a);
                    string key = b->_name + "," +bag_bus[c].at(j)->_name ;
                    string key_inv = bag_bus[c].at(j)->_name + ","+b->_name;
                    if (indexij.find(key) == indexij.end() && indexij.find(key_inv) == indexij.end()) {
                        if (vec_arcs1.size() > 0) {
                            indexij.insert(make_pair<>(key, c));
                            bag_BP_disjoint->_keys.push_back(new index_pair(index_(b->_name), index_(bag_bus[c].at(j)->_name), true));
                        }
                        else {
                            indexij.insert(make_pair<>(key_inv, c));
                            bag_BP_disjoint->_keys.push_back(new index_pair(index_(bag_bus[c].at(j)->_name),index_(b->_name), true));
                        }

                        for (auto a: vec_arcs1) {
                            bag_A_disjoint.push_back((Line*)a);
                        }

                        for (auto a: vec_arcs2) {
                            bag_A_disjoint.push_back((Line*)a);
                        }
                    }
                }
            }
        }

        sort(bag_BP->_keys.begin(), bag_BP->_keys.end(), bus_pair_compare);
        sort(bag_BP_disjoint->_keys.begin(), bag_BP_disjoint->_keys.end(), bus_pair_compare);
        bag_bus_disjoint.push_back(bag_B_disjoint);
        bag_gens_disjoint.push_back(bag_G_disjoint);
        bag_arcs.push_back(bag_A);
        bag_arcs_disjoint.push_back(bag_A_disjoint);
        bag_bus_pairs.push_back(bag_BP);
        bag_bus_pairs_disjoint.push_back(bag_BP_disjoint);
    }

    unsigned  nb_arcs = 0;
    unsigned  nb_buses = 0;
    unsigned  nb_bus_pairs = 0;
    unsigned  nb_gens = 0;

    for (int c =0 ; c < nb_cliques; c++) {
        nb_buses += bag_bus_disjoint[c].size();
        nb_arcs += bag_arcs_disjoint[c].size();
        nb_gens += bag_gens_disjoint[c].size();
        nb_bus_pairs += bag_bus_pairs_disjoint[c]->_keys.size();
    }

    DebugOn("the number of total arcs: " <<  nb_arcs << endl);
    DebugOn("the number of total bus pairs: " <<  nb_bus_pairs << endl);
    DebugOn("the number of total buses: " <<  nb_buses << endl);
    DebugOn("the number of total gens: " <<  nb_gens << endl);
    // intersection_cliques
    for (auto a: cliquetree->arcs){
        std::vector<gravity::index_pair*> v3;
        std::set_intersection(bag_bus_pairs[a->_src->_id]->_keys.begin(), bag_bus_pairs[a->_src->_id]->_keys.end(),
                              bag_bus_pairs[a->_dest->_id]->_keys.begin(), bag_bus_pairs[a->_dest->_id]->_keys.end(),
                              back_inserter(v3), bus_pair_compare);
        if (v3.size() > 0) {
            a->_intersection_clique = v3;
        }
    }

    /** build model */
    Model CLT("Clique tree based Model");

    /** Variables */
    vector<var<Real>> R_Wij;
    vector<var<Real>> Im_Wij;
    vector<var<Real>> Wii;
    var<Real> Sij("Sij");
    var<Real> Cij("Cij");
    Sij^bus_pairs.size();
    Cij^bus_pairs.size();
    vector<var<Real>> Pg;
    vector<var<Real>> Qg;
    for (int c = 0; c < nb_cliques; c++) {
        var<Real>  bag_Wii("Wii_" + to_string(c), grid.w_min.in(bag_bus[c]), grid.w_max.in(bag_bus[c]));
        CLT.add_var(bag_Wii^(bag_bus[c].size()));
        bag_Wii.initialize_all(1.001);
        Wii.push_back(bag_Wii);

        if (bag_bus_pairs[c]->_keys.size() > 0){
            var<Real>  bag_R_Wij("R_Wij_" + to_string(c), grid.wr_min.in(bag_bus_pairs[c]->_keys), grid.wr_max.in(bag_bus_pairs[c]->_keys));
            var<Real>  bag_Im_Wij("Im_Wij_" + to_string(c), grid.wi_min.in(bag_bus_pairs[c]->_keys), grid.wi_max.in(bag_bus_pairs[c]->_keys));
            CLT.add_var(bag_R_Wij^(bag_bus_pairs[c]->_keys.size()));
            CLT.add_var(bag_Im_Wij^(bag_bus_pairs[c]->_keys.size()));
            bag_R_Wij.initialize_all(1.0);
            R_Wij.push_back(bag_R_Wij);
            Im_Wij.push_back(bag_Im_Wij);
        }
        else{
            var<Real> empty("empty");
            empty.set_size(0);
            R_Wij.push_back(empty);
            Im_Wij.push_back(empty);
        }
        
        if (bag_gens_disjoint[c].size() > 0){
            var<Real>  bag_Pg("Pg_" + to_string(c), grid.pg_min.in(bag_gens_disjoint[c]), grid.pg_max.in(bag_gens_disjoint[c]));
            var<Real>  bag_Qg ("Qg_" + to_string(c), grid.qg_min.in(bag_gens_disjoint[c]), grid.qg_max.in(bag_gens_disjoint[c]));
            CLT.add_var(bag_Pg^(bag_gens_disjoint[c].size()));
            CLT.add_var(bag_Qg^(bag_gens_disjoint[c].size()));
            Pg.push_back(bag_Pg);
            Qg.push_back(bag_Qg);
        }
        else{
            var<Real> empty("empty");
            empty.set_size(0);
            Pg.push_back(empty);
            Qg.push_back(empty);
        }
    }

    /* Construct the objective function with generations bound constraints */
    func_ obj;
    for (int c = 0; c < nb_cliques; c++) {
        for (auto g:bag_gens_disjoint[c]) {
            if (g->_active) {
                obj += grid.c1(g->_name)*Pg[c](g->_name)+ grid.c2(g->_name)*Pg[c](g->_name)*Pg[c](g->_name)+grid.c0(g->_name);
            }
        }
    }
    CLT.set_objective(min(obj));

    for (int c = 0; c < nb_cliques; c++) {
        if (bag_bus_pairs_disjoint[c]->_keys.size() > 0) {
            Constraint MC("MC" + to_string(c));
            MC = R_Wij[c].in(bag_bus_pairs_disjoint[c]->_keys) - Cij.in(bag_bus_pairs_disjoint[c]->_keys)
                -Sij.in(bag_bus_pairs_disjoint[c]->_keys);
            CLT.add_constraint(MC = 0);
        }
    }
    for (int c = 0; c < nb_cliques; c++) {
        //* KCL */
        for (auto bus:  bag_bus_disjoint[c]) {
            Constraint KCL_P("KCL_P"+bus->_name);
            Constraint KCL_Q("KCL_Q"+bus->_name);

            for (auto a: bus->get_out()) {
                KCL_P += grid.g_ff(a->_name)*Wii[c](a->_src->_name)
                         +grid.g_ft(a->_name)*R_Wij[c](a->_src->_name+","+a->_dest->_name)
                         +grid.b_ft(a->_name)*Im_Wij[c](a->_src->_name+","+a->_dest->_name);

                KCL_Q += -1*grid.b_ff(a->_name)*Wii[c](a->_src->_name)
                         - grid.b_ft(a->_name)*R_Wij[c](a->_src->_name+","+a->_dest->_name)
                         +grid.g_ft(a->_name)*Im_Wij[c](a->_src->_name+","+a->_dest->_name);
            }

            for (auto a: bus->get_in()) {
                KCL_P  += grid.g_tt(a->_name)*Wii[c](a->_dest->_name)
                          +grid.g_tf(a->_name)*R_Wij[c](a->_src->_name+","+a->_dest->_name)
                          -grid.b_tf(a->_name)*Im_Wij[c](a->_src->_name+","+a->_dest->_name);

                KCL_Q  -= grid.b_tt(a->_name)*Wii[c](a->_dest->_name)
                          + grid.b_tf(a->_name)*R_Wij[c](a->_src->_name+","+a->_dest->_name)
                          + grid.g_tf(a->_name)*Im_Wij[c](a->_src->_name+","+a->_dest->_name);
            }
            if(bus->_has_gen) {
                KCL_P += bus->pl()- sum(Pg[c].in(bus->_gen));
                KCL_Q += bus->ql()- sum(Qg[c].in(bus->_gen));
            }
            else {
                KCL_P += bus->pl();
                KCL_Q += bus->ql();
            }

            /* Shunts */
            KCL_P += bus->gs()*(Wii[c](bus->_name));
            KCL_Q -= bus->bs()*(Wii[c](bus->_name));

            CLT.add_constraint(KCL_P = 0);
            CLT.add_constraint(KCL_Q = 0);
        }

        /* Phase Angle Bounds constraints */
        if (bag_bus_pairs_disjoint[c]->_keys.size() > 0) {
            Constraint PAD_UB("PAD_UB" + to_string(c));
            PAD_UB = Im_Wij[c].in(bag_bus_pairs_disjoint[c]->_keys);
            PAD_UB -= (grid.tan_th_max).in(bag_bus_pairs_disjoint[c]->_keys)*R_Wij[c].in(bag_bus_pairs_disjoint[c]->_keys);
            CLT.add_constraint(PAD_UB <= 0);

            Constraint PAD_LB("PAD_LB" + to_string(c));
            PAD_LB = Im_Wij[c].in(bag_bus_pairs_disjoint[c]->_keys);
            PAD_LB -= (grid.tan_th_min).in(bag_bus_pairs_disjoint[c]->_keys)*R_Wij[c].in(bag_bus_pairs_disjoint[c]->_keys);
            CLT.add_constraint(PAD_LB >= 0);
        }
    }

// COUPLING CONSTRAINTS
    for (auto a: cliquetree->arcs) {
        Constraint Link_Wii("Link_Wii_" + to_string(a->_id));
        Link_Wii += Wii[a->_src->_id].in(a->_intersection);
        Link_Wii -= Wii[a->_dest->_id].in(a->_intersection);
        CLT.add_constraint(Link_Wii = 0);

        if (a->_intersection_clique.size() > 0) {
            auto v3 = a->_intersection_clique;
            Constraint Link_Im_Wij("Link_Im_Wij_" + to_string(a->_id));
            Link_Im_Wij += Im_Wij[a->_src->_id].in(v3);
            Link_Im_Wij -= Im_Wij[a->_dest->_id].in(v3);
            CLT.add_constraint(Link_Im_Wij = 0);

            Constraint Link_R_Wij("Link_R_Wij_" + to_string(a->_id));
            Link_R_Wij += R_Wij[a->_src->_id].in(v3);
            Link_R_Wij -= R_Wij[a->_dest->_id].in(v3);
            CLT.add_constraint(Link_R_Wij = 0);
        }
    }
    //solver SCOPF(CLT, cplex);
    solver SCOPF(CLT, ipopt);
    SCOPF.run();
}

// standard tree decomposition
void reform(PowerNet& grid) {
    // Grid Parameters
    auto bus_pairs = grid.get_bus_pairs();
    auto nb_bus_pairs = bus_pairs.size();
    auto nb_gen = grid.get_nb_active_gens();
    auto nb_buses = grid.get_nb_active_nodes();

    /** build model */
    Model SOCP("SOCP Model");
    /** Variables */
    // lIFTED VARIABLES.
    var<Real>  R_Wij("R_Wij", grid.wr_min.in(bus_pairs), grid.wr_max.in(bus_pairs)); // real part of Wij
    var<Real>  Im_Wij("Im_Wij", grid.wi_min.in(bus_pairs), grid.wi_max.in(bus_pairs)); // imaginary part of Wij.
    var<Real>  Wii("Wii", grid.w_min.in(grid.nodes), grid.w_max.in(grid.nodes));
    SOCP.add_var(Wii^nb_buses);
    SOCP.add_var(R_Wij^nb_bus_pairs);
    SOCP.add_var(Im_Wij^nb_bus_pairs);

    R_Wij.initialize_all(1.0);
    Wii.initialize_all(1.001);

    // power generation
    var<Real> Pg("Pg", grid.pg_min.in(grid.gens), grid.pg_max.in(grid.gens));
    var<Real> Qg ("Qg", grid.qg_min.in(grid.gens), grid.qg_max.in(grid.gens));
    SOCP.add_var(Pg^(nb_gen));
    SOCP.add_var(Qg^(nb_gen));

    /** Construct the objective function*/
    func_ obj;
    for (auto g:grid.gens) {
        if (g->_active) {
            obj += grid.c1(g->_name)*Pg(g->_name) + grid.c2(g->_name)*Pg(g->_name)*Pg(g->_name) + grid.c0(g->_name);
        }
    }
    SOCP.set_objective(min(obj));

    /** Define constraints */
    /* SOCP constraints */
    Constraint SOC("SOC");
    SOC =  power(R_Wij.in(bus_pairs), 2) + power(Im_Wij.in(bus_pairs), 2) - Wii.from(bus_pairs)*Wii.to(bus_pairs) ;
    SOCP.add_constraint(SOC <= 0);

    //KCL
    for (auto b: grid.nodes) {
        Bus* bus = (Bus*) b;
        Constraint KCL_P("KCL_P"+bus->_name);
        Constraint KCL_Q("KCL_Q"+bus->_name);
        for (auto a: b->get_out()) {
            KCL_P  += grid.g_ff(a->_name)*Wii(a->_src->_name)
                      +grid.g_ft(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name)
                      +grid.b_ft(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);
            KCL_Q  += -1*grid.b_ff(a->_name)*Wii(a->_src->_name)
                      - grid.b_ft(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name)
                      +grid.g_ft(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);
        }

        for (auto a: b->get_in()) {
            KCL_P  += grid.g_tt(a->_name)*Wii(a->_dest->_name)
                      + grid.g_tf(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name)
                      - grid.b_tf(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);

            KCL_Q  -= grid.b_tt(a->_name)*Wii(a->_dest->_name)
                      + grid.b_tf(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name)
                      + grid.g_tf(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);
        }

        KCL_P += bus->pl()- sum(Pg.in(bus->_gen));;
        KCL_Q += bus->ql()- sum(Qg.in(bus->_gen));

        /* Shunts */
        KCL_P +=  bus->gs()*(Wii(bus->_name));
        KCL_Q -=  bus->bs()*(Wii(bus->_name));

        SOCP.add_constraint(KCL_P = 0);
        SOCP.add_constraint(KCL_Q = 0);
    }
    //* Phase Angle Bounds constraints */
    Constraint PAD_UB("PAD_UB");
    PAD_UB = Im_Wij.in(bus_pairs);
    PAD_UB -= (grid.tan_th_max).in(bus_pairs)*R_Wij.in(bus_pairs);
    SOCP.add_constraint(PAD_UB <= 0);

    Constraint PAD_LB("PAD_LB");
    PAD_LB =  Im_Wij.in(bus_pairs);
    PAD_LB -= grid.tan_th_min.in(bus_pairs)*R_Wij.in(bus_pairs);
    SOCP.add_constraint(PAD_LB >= 0);

    /* Thermal Limit Constraints */
    Constraint Thermal_Limit_from("Thermal_Limit_from");
    Thermal_Limit_from += power(grid.g_ff.in(grid.arcs)*Wii.from(grid.arcs) + grid.g_ft.in(grid.arcs)*R_Wij.in_pairs(grid.arcs)
                                + grid.b_ft.in(grid.arcs)*Im_Wij.in_pairs(grid.arcs), 2)
                          + power(grid.b_ff.in(grid.arcs)*Wii.from(grid.arcs)+grid.b_ft.in(grid.arcs)*R_Wij.in_pairs(grid.arcs)
                                  -grid.g_ft.in(grid.arcs)*Im_Wij.in_pairs(grid.arcs),2);

    Thermal_Limit_from -= power(grid.S_max.in(grid.arcs),2);
    SOCP.add_constraint(Thermal_Limit_from <= 0);

    Constraint Thermal_Limit_to("Thermal_Limit_to");
    Thermal_Limit_to += power(grid.g_tt.in(grid.arcs)*Wii.to(grid.arcs) + grid.g_tf.in(grid.arcs)*R_Wij.in_pairs(grid.arcs)
                              + grid.b_tf.in(grid.arcs)*Im_Wij.in_pairs(grid.arcs), 2)
                        + power(grid.b_tt.in(grid.arcs)*Wii.to(grid.arcs)+ grid.b_tf.in(grid.arcs)*R_Wij.in_pairs(grid.arcs)
                                + grid.g_tf.in(grid.arcs)*Im_Wij.in_pairs(grid.arcs), 2);

    Thermal_Limit_to -= power(grid.S_max.in(grid.arcs),2);
    SOCP.add_constraint(Thermal_Limit_to <= 0);

    //solver SCOPF(SOCP, cplex);
    solver SCOPF(SOCP, ipopt);
    SCOPF.run();
}

void reform1v(PowerNet& grid) {
    Net* grid_augment = grid.clone_undirected();
    DebugOn("number of edges " << grid_augment->arcs.size() << endl);

    for (auto &node: grid.nodes) {
        vector<Node*> neighbours = node->get_neighbours();
        for (int i=0; i < neighbours.size()-1; i++) {
            auto n1 = grid_augment->get_node(neighbours.at(i)->_name);
            for (int j = i+1; j < neighbours.size(); j++) {
                auto n2 = grid_augment->get_node(neighbours.at(j)->_name);
                if (grid_augment->get_undirected_arc(n1, n2) !=nullptr) {
                    continue;
                }
                else {
                    string name = to_string((int)grid_augment->arcs.size()) + "," +n1->_name + "," + n2->_name;
                    Arc* arc = new Arc(name);
                    arc->_id = grid_augment->arcs.size();
                    arc->_src = n1;
                    arc->_dest = n2;
                    arc->connect();
                    grid_augment->add_undirected_arc(arc);
                }
            }
        }
    }
// note that grid_augment may NOT be chordal.
    grid_augment->get_tree_decomp_bags();
    auto cliquetree = grid_augment->get_clique_tree_kruskal(); // Note: it builds on grid_augment.
    const unsigned nb_cliques = grid_augment->_bags.size();

    vector<vector<Bus*>> bag_bus;
    vector<vector<Bus*>> bag_bus_disjoint;
    vector<vector<Gen*>> bag_gens;
    vector<vector<Gen*>> bag_gens_disjoint;
    vector<vector<Line*>> bag_arcs; //bag_arcs contains the arcs of the power grid while variables associated with W are defined on chordal graph.
    vector<vector<Line*>> bag_arcs_disjoint;
    vector<gravity::node_pairs*> bag_bus_pairs; // bus_pairs in each bag.
    vector<gravity::node_pairs*> bag_bus_pairs_disjoint; // bus_pairs in each bag.

    map<string, unsigned> indexij;
    map<Bus*, unsigned> indexii;

    for (int c = 0; c < nb_cliques; c++) {
        vector<Bus*> bag_B;
        vector<Gen*> bag_G;
        DebugOff("bag " << c << ": ");
        for (int i = 0; i < grid_augment->_bags[c].size(); i++) {
            Node* n = grid_augment->_bags[c].at(i);
            Bus* b = (Bus*) grid.get_node(n->_name);
            DebugOff(b->_name << " ");
            bag_B.push_back(b);
            if (b->_has_gen) {
                bag_G.insert(bag_G.end(), b->_gen.begin(), b->_gen.end());
            }
        }
        DebugOff(" "<< endl);
        bag_bus.push_back(bag_B);
        bag_gens.push_back(bag_G);
    }


    for (int c = 0; c < nb_cliques; c++) {
        vector<Bus*> bag_B_disjoint;
        vector<Gen*> bag_G_disjoint;

        vector<Line*> bag_A;
        vector<Line*> bag_A_disjoint;

        node_pairs* bag_BP = new node_pairs();
        node_pairs* bag_BP_disjoint = new node_pairs();

        sort(bag_bus[c].begin(), bag_bus[c].end(),node_id_compare);
        for (int i = 0; i < bag_bus[c].size(); i++) {
            Bus* b = bag_bus[c].at(i);
            vector<Node*> N = b->get_neighbours();
            sort(N.begin(), N.end(),node_id_compare);

            if (indexii.find(b)==indexii.end()) {
                bool inclusion = std::includes(bag_bus[c].begin(), bag_bus[c].end(), N.begin(),N.end(),node_id_compare);
                if (inclusion) {
                    indexii.insert(make_pair<>(b,c));
                    bag_B_disjoint.push_back(b);
                    if (b->_has_gen) {
                        bag_G_disjoint.insert(bag_G_disjoint.end(), b->_gen.begin(), b->_gen.end());
                    }
                }
            }
            for (int j = i+1; j < bag_bus[c].size(); j++) {
                std::vector<Arc*> vec_arcs1 = grid.get_arcs(b, bag_bus[c].at(j));
                std::vector<Arc*> vec_arcs2 = grid.get_arcs(bag_bus[c].at(j), b);
                if (vec_arcs1.size() +vec_arcs2.size() > 0) {
                    if (vec_arcs1.size() > 0) {
                        bag_BP->_keys.push_back(new index_pair(index_(b->_name), index_(bag_bus[c].at(j)->_name), true));
                    }
                    else {
                        bag_BP->_keys.push_back(new index_pair(index_(bag_bus[c].at(j)->_name),index_(b->_name), true));
                    }

                    for (auto a: vec_arcs1)
                        bag_A.push_back((Line*)a);
                    for (auto a: vec_arcs2)
                        bag_A.push_back((Line*)a);
                    string key = b->_name + "," +bag_bus[c].at(j)->_name ;
                    string key_inv = bag_bus[c].at(j)->_name + ","+b->_name;
                    if (indexij.find(key) == indexij.end() && indexij.find(key_inv) == indexij.end()) {
                        if (vec_arcs1.size() > 0) {
                            indexij.insert(make_pair<>(key, c));
                            bag_BP_disjoint->_keys.push_back(new index_pair(index_(b->_name), index_(bag_bus[c].at(j)->_name), true));
                        }
                        else {
                            indexij.insert(make_pair<>(key_inv, c));
                            bag_BP_disjoint->_keys.push_back(new index_pair(index_(bag_bus[c].at(j)->_name),index_(b->_name), true));
                        }

                        for (auto a: vec_arcs1) {
                            bag_A_disjoint.push_back((Line*)a);
                        }

                        for (auto a: vec_arcs2) {
                            bag_A_disjoint.push_back((Line*)a);
                        }
                    }
                }
            }
        }

        sort(bag_BP->_keys.begin(), bag_BP->_keys.end(), bus_pair_compare);
        sort(bag_BP_disjoint->_keys.begin(), bag_BP_disjoint->_keys.end(), bus_pair_compare);
        bag_bus_disjoint.push_back(bag_B_disjoint);
        bag_gens_disjoint.push_back(bag_G_disjoint);
        bag_arcs.push_back(bag_A);
        bag_arcs_disjoint.push_back(bag_A_disjoint);
        bag_bus_pairs.push_back(bag_BP);
        bag_bus_pairs_disjoint.push_back(bag_BP_disjoint);
    }

    unsigned  nb_arcs = 0;
    unsigned  nb_buses = 0;
    unsigned  nb_bus_pairs = 0;
    unsigned  nb_gens = 0;

    for (int c = 0 ; c < nb_cliques; c++) {
        nb_buses += bag_bus_disjoint[c].size();
        nb_arcs += bag_arcs_disjoint[c].size();
        nb_gens += bag_gens_disjoint[c].size();
        nb_bus_pairs += bag_bus_pairs_disjoint[c]->_keys.size();
    }

    DebugOn("the number of total arcs: " <<  nb_arcs << endl);
    DebugOn("the number of total bus pairs: " <<  nb_bus_pairs << endl);
    DebugOn("the number of total buses: " <<  nb_buses << endl);
    DebugOn("the number of total gens: " <<  nb_gens << endl);

    /** build model */
    Model CLT("Clique tree based Model");

    /** Variables */
    vector<var<Real>> R_V;
    vector<var<Real>> Im_V;
    vector<var<Real>> Pg;
    vector<var<Real>> Qg;
    for (int c = 0; c < nb_cliques; c++) {
        var<Real>  bag_Pg("Pg_" + to_string(c), grid.pg_min.in(bag_gens_disjoint[c]), grid.pg_max.in(bag_gens_disjoint[c]));
        var<Real>  bag_Qg ("Qg_" + to_string(c), grid.qg_min.in(bag_gens_disjoint[c]), grid.qg_max.in(bag_gens_disjoint[c]));
        CLT.add_var(bag_Pg^(bag_gens_disjoint[c].size()));
        CLT.add_var(bag_Qg^(bag_gens_disjoint[c].size()));
        var<Real> Vr("vr", grid.v_max.in(bag_bus[c]));
        var<Real> Vi("vi", grid.v_max.in(bag_bus[c]));
        CLT.add_var(Vr^(bag_bus[c].size()));
        CLT.add_var(Vi^(bag_bus[c].size()));

        Vr.initialize_all(1.0);
        R_V.push_back(Vr);
        Im_V.push_back(Vi);
        Pg.push_back(bag_Pg);
        Qg.push_back(bag_Qg);
    }

    /* Construct the objective function with generations bound constraints */
    func_ obj;
    for (int c = 0; c < nb_cliques; c++) {
        for (auto g:bag_gens_disjoint[c]) {
            if (g->_active) {
                obj += grid.c1(g->_name)*Pg[c](g->_name)+ grid.c2(g->_name)*Pg[c](g->_name)*Pg[c](g->_name)+grid.c0(g->_name);
            }
        }
    }
    CLT.set_objective(min(obj));

    for (int c = 0; c < nb_cliques; c++) {
        //* KCL */
        for (auto bus:  bag_bus_disjoint[c]) {
            Constraint KCL_P("KCL_P"+bus->_name);
            Constraint KCL_Q("KCL_Q"+bus->_name);

            for (auto a: bus->get_out()) {
                string sn = a->_src->_name;
                string dn = a->_dest->_name;
                //assert(std::find(bag_arcs[c].begin(), bag_arcs[c].end(), (Line*) a) != bag_arcs.end());
                KCL_P += grid.g_ff(a->_name)*(power(R_V[c](sn), 2) + power(Im_V[c](sn),2))
                         +grid.g_ft(a->_name)*(R_V[c](sn)*R_V[c](dn) + Im_V[c](sn)*Im_V[c](dn))
                         +grid.b_ft(a->_name)*(Im_V[c](sn)*R_V[c](dn) - R_V[c](sn)*Im_V[c](dn));

                KCL_Q += -1*grid.b_ff(a->_name)*(power(R_V[c](sn), 2) + power(Im_V[c](sn), 2))
                         - grid.b_ft(a->_name)*(R_V[c](sn)*R_V[c](dn) + Im_V[c](sn)*Im_V[c](dn))
                         +grid.g_ft(a->_name)*(Im_V[c](sn)*R_V[c](dn) - R_V[c](sn)*Im_V[c](dn));
            }

            for (auto a: bus->get_in()) {
                string sn = a->_src->_name;
                string dn = a->_dest->_name;
               // assert(std::find(bag_arcs[c].begin(), bag_arcs[c].end(), (Line*) a) != bag_arcs.end());
                KCL_P  += grid.g_tt(a->_name)*(power(R_V[c](dn), 2) + power(Im_V[c](dn), 2))
                          +grid.g_tf(a->_name)*(R_V[c](sn)*R_V[c](dn) + Im_V[c](sn)*Im_V[c](dn))
                          -grid.b_tf(a->_name)*(Im_V[c](sn)*R_V[c](dn) -R_V[c](sn)*Im_V[c](dn));

                KCL_Q  -= grid.b_tt(a->_name)*(power(R_V[c](sn), 2) + power(Im_V[c](sn), 2))
                          + grid.b_tf(a->_name)*(R_V[c](sn)*R_V[c](dn) + Im_V[c](sn)*Im_V[c](dn))
                          + grid.g_tf(a->_name)*(Im_V[c](sn)*R_V[c](dn) - R_V[c](sn)*Im_V[c](dn));
            }

            if(bus->_has_gen) {
                KCL_P += bus->pl()- sum(Pg[c].in(bus->_gen));
                KCL_Q += bus->ql()- sum(Qg[c].in(bus->_gen));
            }
            else {
                KCL_P += bus->pl();
                KCL_Q += bus->ql();
            }

            /* Shunts */
            KCL_P +=  bus->gs()*(power(R_V[c](bus->_name), 2) + power(Im_V[c](bus->_name), 2));
            KCL_Q -=  bus->bs()*(power(R_V[c](bus->_name), 2) + power(Im_V[c](bus->_name), 2));

            CLT.add_constraint(KCL_P = 0);
            CLT.add_constraint(KCL_Q = 0);
        }

        /* Phase Angle Bounds constraints */
        if (bag_bus_pairs_disjoint[c]->_keys.size() > 0) {
            vector<index_pair*> pairs = bag_bus_pairs_disjoint[c]->_keys;
            Constraint PAD_UB("PAD_UB" + to_string(c));
            PAD_UB = Im_V[c].from(pairs)*R_V[c].to(pairs) - R_V[c].from(pairs)*Im_V[c].to(pairs);
            PAD_UB -= grid.tan_th_max.in(pairs)*(R_V[c].from(pairs)*R_V[c].to(pairs)
                                                 + Im_V[c].from(pairs)*Im_V[c].to(pairs));
            CLT.add_constraint(PAD_UB <= 0);

            Constraint PAD_LB("PAD_LB" + to_string(c));
            PAD_LB = Im_V[c].from(pairs)*R_V[c].to(pairs)- R_V[c].from(pairs)*Im_V[c].to(pairs);
            PAD_LB -= grid.tan_th_min.in(pairs)*(R_V[c].from(pairs)*R_V[c].to(pairs)
                                                 + Im_V[c].to(pairs)*Im_V[c].from(pairs));
            CLT.add_constraint(PAD_LB >= 0);
        }

        if (bag_arcs_disjoint[c].size() > 0) {
            /* Thermal Limit Constraints */
            Constraint Thermal_Limit_from("Thermal_Limit_from"+to_string(c));
            Thermal_Limit_from += power(grid.g_ff.in(bag_arcs_disjoint[c])*(power(R_V[c].from(bag_arcs_disjoint[c]), 2) + power(Im_V[c].from(bag_arcs_disjoint[c]),2))
                                        +grid.g_ft.in(bag_arcs_disjoint[c])*(R_V[c].from(bag_arcs_disjoint[c])*R_V[c].to(bag_arcs_disjoint[c])
                                                +Im_V[c].from(bag_arcs_disjoint[c])*Im_V[c].to(bag_arcs_disjoint[c]))
                                        +grid.b_ft.in(bag_arcs_disjoint[c])*(Im_V[c].from(bag_arcs_disjoint[c])*R_V[c].to(bag_arcs_disjoint[c])
                                                - R_V[c].from(bag_arcs_disjoint[c])*Im_V[c].to(bag_arcs_disjoint[c])), 2)
                                  + power(grid.b_ff.in(bag_arcs_disjoint[c])*(power(R_V[c].from(bag_arcs_disjoint[c]), 2) + power(Im_V[c].from(bag_arcs_disjoint[c]),2))
                                          +grid.b_ft.in(bag_arcs_disjoint[c])*(R_V[c].from(bag_arcs_disjoint[c])*R_V[c].to(bag_arcs_disjoint[c])
                                                  + Im_V[c].from(bag_arcs_disjoint[c])*Im_V[c].to(bag_arcs_disjoint[c]))
                                          -grid.g_ft.in(bag_arcs_disjoint[c])*(Im_V[c].from(bag_arcs_disjoint[c])*R_V[c].to(bag_arcs_disjoint[c])
                                                  - R_V[c].from(bag_arcs_disjoint[c])*Im_V[c].to(bag_arcs_disjoint[c])), 2);

            Thermal_Limit_from -= power(grid.S_max.in(bag_arcs_disjoint[c]), 2);
            CLT.add_constraint(Thermal_Limit_from <= 0);


            //Constraint Thermal_Limit_to("Thermal_Limit_to" + to_string(c));
            //Thermal_Limit_to += power(grid.g_tt.in(bag_arcs_disjoint[c])*Wii[c].to(bag_arcs_disjoint[c])
            //                          + grid.g_tf.in(bag_arcs_disjoint[c])*R_Wij[c].in_pairs(bag_arcs_disjoint[c])
            //                          + grid.b_tf.in(bag_arcs_disjoint[c])*Im_Wij[c].in_pairs(bag_arcs_disjoint[c]), 2)
            //                    + power(grid.b_tt.in(bag_arcs_disjoint[c])*Wii[c].to(bag_arcs_disjoint[c])
            //                            + grid.b_tf.in(bag_arcs_disjoint[c])*R_Wij[c].in_pairs(bag_arcs_disjoint[c])
            //                            + grid.g_tf.in(bag_arcs_disjoint[c])*Im_Wij[c].in_pairs(bag_arcs_disjoint[c]), 2);

            //Thermal_Limit_to -= power(grid.S_max.in(bag_arcs_disjoint[c]), 2);
            //CLT.add_constraint(Thermal_Limit_to <= 0);
        }
    }
// COUPLING CONSTRAINTS
    for (auto a: cliquetree->arcs) {
        Constraint Link_Im_V("Link_Im_V_" + to_string(a->_id));
        Link_Im_V += Im_V[a->_src->_id].in(a->_intersection);
        Link_Im_V -= Im_V[a->_dest->_id].in(a->_intersection);
        CLT.add_constraint(Link_Im_V = 0);

        Constraint Link_R_V("Link_R_V_" + to_string(a->_id));
        Link_R_V += R_V[a->_src->_id].in(a->_intersection);
        Link_R_V -= R_V[a->_dest->_id].in(a->_intersection);
        CLT.add_constraint(Link_R_V = 0);
    }
    //solver SCOPF(CLT, cplex);
    solver SCOPF(CLT, ipopt);
    SCOPF.run();
}


// Reform1 augmented graph
void reform1(PowerNet& grid) {
    Net* grid_augment = grid.clone_undirected();
    DebugOn("number of edges " << grid_augment->arcs.size() << endl);

    for (auto &node: grid.nodes) {
        vector<Node*> neighbours = node->get_neighbours();
        for (int i = 0; i < neighbours.size()-1; i++) {
            auto n1 = grid_augment->get_node(neighbours.at(i)->_name);
            for (int j = i+1; j < neighbours.size(); j++) {
                auto n2 = grid_augment->get_node(neighbours.at(j)->_name);
                if (grid_augment->get_undirected_arc(n1, n2) !=nullptr) {
                    continue;
                }
                else {
                    string name = to_string((int)grid_augment->arcs.size()) + "," +n1->_name + "," + n2->_name;
                    Arc* arc = new Arc(name);
                    arc->_id = grid_augment->arcs.size();
                    arc->_src = n1;
                    arc->_dest = n2;
                    arc->connect();
                    grid_augment->add_undirected_arc(arc);
                }
            }
        }
    }
// note that grid_augment may NOT be chordal.
    grid_augment->get_tree_decomp_bags();
    auto cliquetree = grid_augment->get_clique_tree_kruskal(); // Note: it builds on grid_augment.
    const unsigned nb_cliques = grid_augment->_bags.size();

    vector<vector<Bus*>> bag_bus;
    vector<vector<Bus*>> bag_bus_disjoint;
    vector<vector<Gen*>> bag_gens;
    vector<vector<Gen*>> bag_gens_disjoint;
    vector<vector<Line*>> bag_arcs; //bag_arcs contains the arcs of the power grid while variables associated with W are defined on chordal graph.
    vector<vector<Line*>> bag_arcs_disjoint;
    vector<gravity::node_pairs*> bag_bus_pairs; // bus_pairs in each bag.
    vector<gravity::node_pairs*> bag_bus_pairs_disjoint; // bus_pairs in each bag.

    map<string, unsigned> indexij;
    map<Bus*, unsigned> indexii;

    for (int c = 0; c < nb_cliques; c++) {
        vector<Bus*> bag_B;
        vector<Gen*> bag_G;
        DebugOff("bag " << c << ": ");
        for (int i = 0; i < grid_augment->_bags[c].size(); i++) {
            Node* n = grid_augment->_bags[c].at(i);
            Bus* b = (Bus*) grid.get_node(n->_name);
            DebugOff(b->_name << " ");
            bag_B.push_back(b);
            if (b->_has_gen) {
                bag_G.insert(bag_G.end(), b->_gen.begin(), b->_gen.end());
            }
        }
        DebugOff(" "<< endl);
        bag_bus.push_back(bag_B);
        bag_gens.push_back(bag_G);
    }


    for (int c = 0; c < nb_cliques; c++) {
        vector<Bus*> bag_B_disjoint;
        vector<Gen*> bag_G_disjoint;

        vector<Line*> bag_A;
        vector<Line*> bag_A_disjoint;

        node_pairs* bag_BP = new node_pairs();
        node_pairs* bag_BP_disjoint = new node_pairs();

        sort(bag_bus[c].begin(), bag_bus[c].end(),node_id_compare);
        for (int i = 0; i < bag_bus[c].size(); i++) {
            Bus* b = bag_bus[c].at(i);
            vector<Node*> N = b->get_neighbours();
            sort(N.begin(), N.end(),node_id_compare);

            if (indexii.find(b)==indexii.end()) {
                bool inclusion = std::includes(bag_bus[c].begin(), bag_bus[c].end(), N.begin(),N.end(),node_id_compare);
                if (inclusion) {
                    indexii.insert(make_pair<>(b,c));
                    bag_B_disjoint.push_back(b);
                    if (b->_has_gen) {
                        bag_G_disjoint.insert(bag_G_disjoint.end(), b->_gen.begin(), b->_gen.end());
                    }
                }
            }
            for (int j = i+1; j < bag_bus[c].size(); j++) {
                std::vector<Arc*> vec_arcs1 = grid.get_arcs(b, bag_bus[c].at(j));
                std::vector<Arc*> vec_arcs2 = grid.get_arcs(bag_bus[c].at(j), b);
                if (vec_arcs1.size() +vec_arcs2.size() > 0) {
                    if (vec_arcs1.size() > 0) {
                        bag_BP->_keys.push_back(new index_pair(index_(b->_name), index_(bag_bus[c].at(j)->_name), true));
                    }
                    else {
                        bag_BP->_keys.push_back(new index_pair(index_(bag_bus[c].at(j)->_name),index_(b->_name), true));
                    }

                    for (auto a: vec_arcs1)
                        bag_A.push_back((Line*)a);
                    for (auto a: vec_arcs2)
                        bag_A.push_back((Line*)a);
                    string key = b->_name + "," +bag_bus[c].at(j)->_name ;
                    string key_inv = bag_bus[c].at(j)->_name + ","+b->_name;
                    if (indexij.find(key) == indexij.end() && indexij.find(key_inv) == indexij.end()) {
                        if (vec_arcs1.size() > 0) {
                            indexij.insert(make_pair<>(key, c));
                            bag_BP_disjoint->_keys.push_back(new index_pair(index_(b->_name), index_(bag_bus[c].at(j)->_name), true));
                        }
                        else {
                            indexij.insert(make_pair<>(key_inv, c));
                            bag_BP_disjoint->_keys.push_back(new index_pair(index_(bag_bus[c].at(j)->_name),index_(b->_name), true));
                        }

                        for (auto a: vec_arcs1) {
                            bag_A_disjoint.push_back((Line*)a);
                        }

                        for (auto a: vec_arcs2) {
                            bag_A_disjoint.push_back((Line*)a);
                        }
                    }
                }
            }
        }

        sort(bag_BP->_keys.begin(), bag_BP->_keys.end(), bus_pair_compare);
        sort(bag_BP_disjoint->_keys.begin(), bag_BP_disjoint->_keys.end(), bus_pair_compare);
        bag_bus_disjoint.push_back(bag_B_disjoint);
        bag_gens_disjoint.push_back(bag_G_disjoint);
        bag_arcs.push_back(bag_A);
        bag_arcs_disjoint.push_back(bag_A_disjoint);
        bag_bus_pairs.push_back(bag_BP);
        bag_bus_pairs_disjoint.push_back(bag_BP_disjoint);
    }

    unsigned  nb_arcs = 0;
    unsigned  nb_buses = 0;
    unsigned  nb_bus_pairs = 0;
    unsigned  nb_gens = 0;

    for (int c =0 ; c < nb_cliques; c++) {
        nb_buses += bag_bus_disjoint[c].size();
        nb_arcs += bag_arcs_disjoint[c].size();
        nb_gens += bag_gens_disjoint[c].size();
        nb_bus_pairs += bag_bus_pairs_disjoint[c]->_keys.size();
    }

    DebugOn("the number of total arcs: " <<  nb_arcs << endl);
    DebugOn("the number of total bus pairs: " <<  nb_bus_pairs << endl);
    DebugOn("the number of total buses: " <<  nb_buses << endl);
    DebugOn("the number of total gens: " <<  nb_gens << endl);
    // intersection_cliques
    for (auto a: cliquetree->arcs){
        std::vector<gravity::index_pair*> v3;
        std::set_intersection(bag_bus_pairs[a->_src->_id]->_keys.begin(), bag_bus_pairs[a->_src->_id]->_keys.end(),
                              bag_bus_pairs[a->_dest->_id]->_keys.begin(), bag_bus_pairs[a->_dest->_id]->_keys.end(),
                              back_inserter(v3), bus_pair_compare);
        if (v3.size() > 0) {
            a->_intersection_clique = v3;
        }
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
        var<Real>  bag_Wii("Wii_" + to_string(c), grid.w_min.in(bag_bus[c]), grid.w_max.in(bag_bus[c]));
        CLT.add_var(bag_Wii^(bag_bus[c].size()));
        bag_Wii.initialize_all(1.001);
        Wii.push_back(bag_Wii);

        if (bag_bus_pairs[c]->_keys.size() > 0){
            var<Real>  bag_R_Wij("R_Wij_" + to_string(c), grid.wr_min.in(bag_bus_pairs[c]->_keys), grid.wr_max.in(bag_bus_pairs[c]->_keys));
            var<Real>  bag_Im_Wij("Im_Wij_" + to_string(c), grid.wi_min.in(bag_bus_pairs[c]->_keys), grid.wi_max.in(bag_bus_pairs[c]->_keys));
            CLT.add_var(bag_R_Wij^(bag_bus_pairs[c]->_keys.size()));
            CLT.add_var(bag_Im_Wij^(bag_bus_pairs[c]->_keys.size()));
            bag_R_Wij.initialize_all(1.0);
            R_Wij.push_back(bag_R_Wij);
            Im_Wij.push_back(bag_Im_Wij);
        }
        else{
            var<Real> empty("empty");
            empty.set_size(0);
            R_Wij.push_back(empty);
            Im_Wij.push_back(empty);
        }
        
        if (bag_gens_disjoint[c].size() > 0){
            var<Real>  bag_Pg("Pg_" + to_string(c), grid.pg_min.in(bag_gens_disjoint[c]), grid.pg_max.in(bag_gens_disjoint[c]));
            var<Real>  bag_Qg ("Qg_" + to_string(c), grid.qg_min.in(bag_gens_disjoint[c]), grid.qg_max.in(bag_gens_disjoint[c]));
            CLT.add_var(bag_Pg^(bag_gens_disjoint[c].size()));
            CLT.add_var(bag_Qg^(bag_gens_disjoint[c].size()));
            Pg.push_back(bag_Pg);
            Qg.push_back(bag_Qg);
        }
        else{
            var<Real> empty("empty");
            empty.set_size(0);
            Pg.push_back(empty);
            Qg.push_back(empty);
        }
    }

    /* Construct the objective function with generations bound constraints */
    func_ obj;
    for (int c = 0; c < nb_cliques; c++) {
        for (auto g:bag_gens_disjoint[c]) {
            if (g->_active) {
                obj += grid.c1(g->_name)*Pg[c](g->_name)+ grid.c2(g->_name)*Pg[c](g->_name)*Pg[c](g->_name)+grid.c0(g->_name);
            }
        }
    }
    CLT.set_objective(min(obj));

    for (int c = 0; c < nb_cliques; c++) {
        if (bag_bus_pairs_disjoint[c]->_keys.size()>0) {
            Constraint SOC("SOC_" + to_string(c));
            SOC =  power(R_Wij[c].in(bag_bus_pairs_disjoint[c]->_keys), 2)
                   + power(Im_Wij[c].in(bag_bus_pairs_disjoint[c]->_keys), 2)
                   - Wii[c].from(bag_bus_pairs_disjoint[c]->_keys)*Wii[c].to(bag_bus_pairs_disjoint[c]->_keys) ;
            CLT.add_constraint(SOC <= 0);
        }

        //* KCL */
        for (auto bus:  bag_bus_disjoint[c]) {
            Constraint KCL_P("KCL_P"+bus->_name);
            Constraint KCL_Q("KCL_Q"+bus->_name);

            for (auto a: bus->get_out()) {
                //assert(std::find(bag_arcs[c].begin(), bag_arcs[c].end(), a) != bag_arcs.end());
                KCL_P += grid.g_ff(a->_name)*Wii[c](a->_src->_name)
                         +grid.g_ft(a->_name)*R_Wij[c](a->_src->_name+","+a->_dest->_name)
                         +grid.b_ft(a->_name)*Im_Wij[c](a->_src->_name+","+a->_dest->_name);

                KCL_Q += -1*grid.b_ff(a->_name)*Wii[c](a->_src->_name)
                         - grid.b_ft(a->_name)*R_Wij[c](a->_src->_name+","+a->_dest->_name)
                         +grid.g_ft(a->_name)*Im_Wij[c](a->_src->_name+","+a->_dest->_name);
            }

            for (auto a: bus->get_in()) {
                //assert(std::find(bag_arcs[c].begin(), bag_arcs[c].end(), a) != bag_arcs.end());
                KCL_P  += grid.g_tt(a->_name)*Wii[c](a->_dest->_name)
                          +grid.g_tf(a->_name)*R_Wij[c](a->_src->_name+","+a->_dest->_name)
                          -grid.b_tf(a->_name)*Im_Wij[c](a->_src->_name+","+a->_dest->_name);

                KCL_Q  -= grid.b_tt(a->_name)*Wii[c](a->_dest->_name)
                          + grid.b_tf(a->_name)*R_Wij[c](a->_src->_name+","+a->_dest->_name)
                          + grid.g_tf(a->_name)*Im_Wij[c](a->_src->_name+","+a->_dest->_name);
            }
            if(bus->_has_gen) {
                KCL_P += bus->pl()- sum(Pg[c].in(bus->_gen));
                KCL_Q += bus->ql()- sum(Qg[c].in(bus->_gen));
            }
            else {
                KCL_P += bus->pl();
                KCL_Q += bus->ql();
            }

            /* Shunts */
            KCL_P += bus->gs()*(Wii[c](bus->_name));
            KCL_Q -= bus->bs()*(Wii[c](bus->_name));

            CLT.add_constraint(KCL_P = 0);
            CLT.add_constraint(KCL_Q = 0);
        }

        /* Phase Angle Bounds constraints */
        if (bag_bus_pairs_disjoint[c]->_keys.size() > 0) {
            Constraint PAD_UB("PAD_UB" + to_string(c));
            PAD_UB = Im_Wij[c].in(bag_bus_pairs_disjoint[c]->_keys);
            PAD_UB -= (grid.tan_th_max).in(bag_bus_pairs_disjoint[c]->_keys)*R_Wij[c].in(bag_bus_pairs_disjoint[c]->_keys);
            CLT.add_constraint(PAD_UB <= 0);

            Constraint PAD_LB("PAD_LB" + to_string(c));
            PAD_LB = Im_Wij[c].in(bag_bus_pairs_disjoint[c]->_keys);
            PAD_LB -= (grid.tan_th_min).in(bag_bus_pairs_disjoint[c]->_keys)*R_Wij[c].in(bag_bus_pairs_disjoint[c]->_keys);
            CLT.add_constraint(PAD_LB >= 0);
        }

        if (bag_arcs_disjoint[c].size() > 0) {
            /* Thermal Limit Constraints */
            Constraint Thermal_Limit_from("Thermal_Limit_from"+to_string(c));
            Thermal_Limit_from += power(grid.g_ff.in(bag_arcs_disjoint[c])*Wii[c].from(bag_arcs_disjoint[c])+
                                        grid.g_ft.in(bag_arcs_disjoint[c])*R_Wij[c].in_pairs(bag_arcs_disjoint[c])
                                        +grid.b_ft.in(bag_arcs_disjoint[c])*Im_Wij[c].in_pairs(bag_arcs_disjoint[c]), 2)
                                  + power(grid.b_ff.in(bag_arcs_disjoint[c])*Wii[c].from(bag_arcs_disjoint[c])
                                          +grid.b_ft.in(bag_arcs_disjoint[c])*R_Wij[c].in_pairs(bag_arcs_disjoint[c])
                                          -grid.g_ft.in(bag_arcs_disjoint[c])*Im_Wij[c].in_pairs(bag_arcs_disjoint[c]), 2);

            Thermal_Limit_from -= power(grid.S_max.in(bag_arcs_disjoint[c]), 2);
            CLT.add_constraint(Thermal_Limit_from <= 0);


            Constraint Thermal_Limit_to("Thermal_Limit_to" + to_string(c));
            Thermal_Limit_to += power(grid.g_tt.in(bag_arcs_disjoint[c])*Wii[c].to(bag_arcs_disjoint[c])
                                      + grid.g_tf.in(bag_arcs_disjoint[c])*R_Wij[c].in_pairs(bag_arcs_disjoint[c])
                                      + grid.b_tf.in(bag_arcs_disjoint[c])*Im_Wij[c].in_pairs(bag_arcs_disjoint[c]), 2)
                                + power(grid.b_tt.in(bag_arcs_disjoint[c])*Wii[c].to(bag_arcs_disjoint[c])
                                        + grid.b_tf.in(bag_arcs_disjoint[c])*R_Wij[c].in_pairs(bag_arcs_disjoint[c])
                                        + grid.g_tf.in(bag_arcs_disjoint[c])*Im_Wij[c].in_pairs(bag_arcs_disjoint[c]), 2);

            Thermal_Limit_to -= power(grid.S_max.in(bag_arcs_disjoint[c]), 2);
            CLT.add_constraint(Thermal_Limit_to <= 0);
        }
    }

// COUPLING CONSTRAINTS
    for (auto a: cliquetree->arcs) {
        Constraint Link_Wii("Link_Wii_" + to_string(a->_id));
        Link_Wii += Wii[a->_src->_id].in(a->_intersection);
        Link_Wii -= Wii[a->_dest->_id].in(a->_intersection);
        CLT.add_constraint(Link_Wii = 0);

        if (a->_intersection_clique.size() > 0) {
            auto v3 = a->_intersection_clique;
            Constraint Link_Im_Wij("Link_Im_Wij_" + to_string(a->_id));
            Link_Im_Wij += Im_Wij[a->_src->_id].in(v3);
            Link_Im_Wij -= Im_Wij[a->_dest->_id].in(v3);
            CLT.add_constraint(Link_Im_Wij = 0);

            Constraint Link_R_Wij("Link_R_Wij_" + to_string(a->_id));
            Link_R_Wij += R_Wij[a->_src->_id].in(v3);
            Link_R_Wij -= R_Wij[a->_dest->_id].in(v3);
            CLT.add_constraint(Link_R_Wij = 0);
//            Constraint Link_Pf("Link_Pf_" + to_string(a->_id));
//            Constraint Link_Qf("Link_Qf_" + to_string(a->_id));
//
//            unsigned s = a->_src->_id;
//            unsigned d = a->_dest->_id;
//            for (auto arc: v3) {
//                string pn = arc->_src->_name +","+ arc->_dest->_name;
//                Link_Pf += (grid.g_ff(arc->_name)*Wii[s](arc->_src->_name) + grid.g_ft(arc->_name)*R_Wij[s](pn)
//                               +grid.b_ft(arc->_name)*Im_Wij[s](pn));
//
//                Link_Pf += (grid.g_tt(arc->_name)*Wii[s](arc->_dest->_name) +grid.g_tf(arc->_name)*R_Wij[s](pn)
//                               + grid.b_tf(arc->_name)*Im_Wij[s](pn));
//
//                Link_Pf -= (grid.g_ff(arc->_name)*Wii[d](arc->_src->_name) + grid.g_ft(arc->_name)*R_Wij[d](pn)
//                               +grid.b_ft(arc->_name)*Im_Wij[d](pn));
//
//                Link_Pf -= (grid.g_tt(arc->_name)*Wii[d](arc->_dest->_name) +grid.g_tf(arc->_name)*R_Wij[d](pn)
//                               + grid.b_tf(arc->_name)*Im_Wij[d](pn));
//
//
//                Link_Qf += (-1*grid.b_ff(arc->_name)*Wii[s](arc->_src->_name) -grid.b_ft(arc->_name)*R_Wij[s](pn)
//                               +grid.g_ft(arc->_name)*Im_Wij[s](pn));
//
//                Link_Qf -= (grid.b_tt(arc->_name)*Wii[s](arc->_dest->_name)+ grid.b_tf(arc->_name)*R_Wij[s](pn)
//                               + grid.g_tf(arc->_name)*Im_Wij[s](pn));
//
//                Link_Qf -= (-1*grid.b_ff(arc->_name)*Wii[d](arc->_src->_name) -grid.b_ft(arc->_name)*R_Wij[d](pn)
//                               +grid.g_ft(arc->_name)*Im_Wij[d](pn));
//
//                Link_Qf += (grid.b_tt(arc->_name)*Wii[d](arc->_dest->_name)+ grid.b_tf(arc->_name)*R_Wij[d](pn)
//                               + grid.g_tf(arc->_name)*Im_Wij[d](pn));
//            }
//
//            CLT.add_constraint(Link_Pf = 0);
//            CLT.add_constraint(Link_Qf = 0);
        }
    }
    //solver SCOPF(CLT, cplex);
    solver SCOPF(CLT, ipopt);
    SCOPF.run();
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
        // fname = "../../data_sets/Power/nesta_case3_lmbd.m";

        //fname = "../../data_sets/Power/nesta_case9241_pegase.m";
        //fname = "../../data_sets/Power/nesta_case2383wp_mp.m";
        //fname = "../../data_sets/Power/nesta_case1354_pegase_api.m";
        //fname = "../../data_sets/Power/nesta_case300_ieee.m";
        fname = "../../data_sets/Power/nesta_case118_ieee.m";
        //fname = "../../data_sets/Power/nesta_case57_ieee.m";

    }
    // ACOPF
    PowerNet grid;
    grid.readgrid(fname);
    cout << "////////////////////////////////////////" << endl;

    //reform(grid);
    //reform1(grid);
    //LP(grid);
    //reform1v(grid);
    //return 0;
    /** Clique tree decomposition **/
    grid.get_tree_decomp_bags();
    auto cliquetree = grid.get_clique_tree_kruskal();
    const unsigned nb_cliques = grid._bags.size();

    vector<vector<Bus*>> bag_bus;
    vector<vector<Bus*>> bag_bus_disjoint;
    vector<vector<Gen*>> bag_gens;
    vector<vector<Gen*>> bag_gens_disjoint;
    vector<vector<Line*>> bag_arcs; //bag_arcs contains the arcs of the power grid while variables associated with W are defined on chordal graph.
    vector<vector<Line*>> bag_arcs_disjoint;
    vector<gravity::node_pairs*> bag_bus_pairs;
    vector<gravity::node_pairs*> bag_bus_pairs_disjoint;
    map<string, unsigned> indexij;
    map<Bus*, unsigned> indexii;

    //for (int c = 0; c < nb_cliques; c++) {
    //    vector<Bus*> bag_B;
    //    vector<Bus*> bag_B_disjoint;
    //    vector<Gen*> bag_G;
    //    vector<Gen*> bag_G_disjoint;
    //    for (int i = 0; i < grid._bags[c].size(); i++) {
    //        Bus* b = (Bus*) grid._bags[c].at(i);
    //        if (b !=nullptr) {
    //            bag_B.push_back(b);
    //            auto pp = indexii.insert(make_pair<>(b, c));
    //            if (pp.second) {
    //                bag_B_disjoint.push_back(b);
    //            }

    //            if (b->_has_gen) {
    //                bag_G.insert(bag_G.end(), b->_gen.begin(), b->_gen.end());
    //            }

    //            if (pp.second && b->_has_gen) {
    //                bag_G_disjoint.insert(bag_G_disjoint.end(), b->_gen.begin(), b->_gen.end());
    //            }
    //        }
    //    }

    //    bag_bus.push_back(bag_B);
    //    bag_gens.push_back(bag_G);
    //    bag_bus_disjoint.push_back(bag_B_disjoint);
    //    bag_gens_disjoint.push_back(bag_G_disjoint);
    //}
    for (int c = 0; c < nb_cliques; c++) {
        vector<Bus*> bag_B;
        vector<Gen*> bag_G;
        DebugOff("bag " << c << ": ");
        for (int i = 0; i < grid._bags[c].size(); i++) {
            Node* n = grid._bags[c].at(i);
            Bus* b = (Bus*) grid.get_node(n->_name);
            DebugOff(b->_name << " ");
            bag_B.push_back(b);
            if (b->_has_gen) {
                bag_G.insert(bag_G.end(), b->_gen.begin(), b->_gen.end());
            }
        }
        DebugOff(" "<< endl);
        bag_bus.push_back(bag_B);
        bag_gens.push_back(bag_G);
    }

    //for (int c = 0; c < nb_cliques; c++) {
    //    vector<Line*> bag_A;
    //    vector<Line*> bag_A_disjoint;

    //    node_pairs* bag_BP = new node_pairs();
    //    node_pairs* bag_BP_disjoint = new node_pairs();

    //    for (int i = 0; i < grid._bags[c].size(); i++) {
    //        Bus* b = bag_bus[c].at(i);
    //        for (int j = i+1; j < grid._bags[c].size(); j++) {
    //            std::vector<Arc*> vec_arcs1 = grid.get_arcs(b, bag_bus[c].at(j));
    //            std::vector<Arc*> vec_arcs2 = grid.get_arcs(bag_bus[c].at(j), b);
    //            if (vec_arcs1.size() +vec_arcs2.size() > 0) {
    //                if (vec_arcs1.size() > 0) {
    //                    bag_BP->_keys.push_back(new index_pair(index_(b->_name), index_(bag_bus[c].at(j)->_name), true));
    //                }
    //                else {
    //                    bag_BP->_keys.push_back(new index_pair(index_(bag_bus[c].at(j)->_name),index_(b->_name), true));
    //                }

    //                for (auto a: vec_arcs1)
    //                    bag_A.push_back((Line*)a);
    //                for (auto a: vec_arcs2)
    //                    bag_A.push_back((Line*)a);
    //                string key = b->_name + "," +bag_bus[c].at(j)->_name ;
    //                string key_inv = bag_bus[c].at(j)->_name + ","+b->_name;
    //                if (indexij.find(key) == indexij.end() && indexij.find(key_inv) == indexij.end()) {
    //                    if (vec_arcs1.size() > 0) {
    //                        indexij.insert(make_pair<>(key, c));
    //                        bag_BP_disjoint->_keys.push_back(new index_pair(index_(b->_name), index_(bag_bus[c].at(j)->_name), true));
    //                    }
    //                    else {
    //                        indexij.insert(make_pair<>(key_inv, c));
    //                        bag_BP_disjoint->_keys.push_back(new index_pair(index_(bag_bus[c].at(j)->_name),index_(b->_name), true));
    //                    }

    //                    for (auto a: vec_arcs1) {
    //                        bag_A_disjoint.push_back((Line*)a);
    //                    }

    //                    for (auto a: vec_arcs2) {
    //                        bag_A_disjoint.push_back((Line*)a);
    //                    }
    //                }
    //            }
    //        }
    //    }


    //    sort(bag_BP->_keys.begin(), bag_BP->_keys.end(), bus_pair_compare);
    //    sort(bag_BP_disjoint->_keys.begin(), bag_BP_disjoint->_keys.end(), bus_pair_compare);
    //    bag_arcs.push_back(bag_A);
    //    bag_arcs_disjoint.push_back(bag_A_disjoint);
    //    bag_bus_pairs.push_back(bag_BP);
    //    bag_bus_pairs_disjoint.push_back(bag_BP_disjoint);
    //}

    for (int c = 0; c < nb_cliques; c++) {
        vector<Bus*> bag_B_disjoint;
        vector<Gen*> bag_G_disjoint;

        vector<Line*> bag_A;
        vector<Line*> bag_A_disjoint;

        node_pairs* bag_BP = new node_pairs();
        node_pairs* bag_BP_disjoint = new node_pairs();

        sort(bag_bus[c].begin(), bag_bus[c].end(),node_id_compare);
        for (int i = 0; i < bag_bus[c].size(); i++) {
            Bus* b = bag_bus[c].at(i);
            vector<Node*> N = b->get_neighbours();
            sort(N.begin(), N.end(),node_id_compare);

            if (indexii.find(b)==indexii.end()) {
                //bool inclusion = std::includes(bag_bus[c].begin(), bag_bus[c].end(), N.begin(),N.end(),node_id_compare);
                //if (inclusion) {
                    indexii.insert(make_pair<>(b,c));
                    bag_B_disjoint.push_back(b);
                    if (b->_has_gen) {
                        bag_G_disjoint.insert(bag_G_disjoint.end(), b->_gen.begin(), b->_gen.end());
                    }
                //}
            }
            for (int j = i+1; j < bag_bus[c].size(); j++) {
                std::vector<Arc*> vec_arcs1 = grid.get_arcs(b, bag_bus[c].at(j));
                std::vector<Arc*> vec_arcs2 = grid.get_arcs(bag_bus[c].at(j), b);
                if (vec_arcs1.size() +vec_arcs2.size() > 0) {
                    if (vec_arcs1.size() > 0) {
                        bag_BP->_keys.push_back(new index_pair(index_(b->_name), index_(bag_bus[c].at(j)->_name), true));
                    }
                    else {
                        bag_BP->_keys.push_back(new index_pair(index_(bag_bus[c].at(j)->_name),index_(b->_name), true));
                    }

                    for (auto a: vec_arcs1)
                        bag_A.push_back((Line*)a);
                    for (auto a: vec_arcs2)
                        bag_A.push_back((Line*)a);
                    string key = b->_name + "," +bag_bus[c].at(j)->_name ;
                    string key_inv = bag_bus[c].at(j)->_name + ","+b->_name;
                    if (indexij.find(key) == indexij.end() && indexij.find(key_inv) == indexij.end()) {
                        if (vec_arcs1.size() > 0) {
                            indexij.insert(make_pair<>(key, c));
                            bag_BP_disjoint->_keys.push_back(new index_pair(index_(b->_name), index_(bag_bus[c].at(j)->_name), true));
                        }
                        else {
                            indexij.insert(make_pair<>(key_inv, c));
                            bag_BP_disjoint->_keys.push_back(new index_pair(index_(bag_bus[c].at(j)->_name),index_(b->_name), true));
                        }

                        for (auto a: vec_arcs1) {
                            bag_A_disjoint.push_back((Line*)a);
                        }

                        for (auto a: vec_arcs2) {
                            bag_A_disjoint.push_back((Line*)a);
                        }
                    }
                }
            }
        }

        sort(bag_BP->_keys.begin(), bag_BP->_keys.end(), bus_pair_compare);
        sort(bag_BP_disjoint->_keys.begin(), bag_BP_disjoint->_keys.end(), bus_pair_compare);
        bag_bus_disjoint.push_back(bag_B_disjoint);
        bag_gens_disjoint.push_back(bag_G_disjoint);
        bag_arcs.push_back(bag_A);
        bag_arcs_disjoint.push_back(bag_A_disjoint);
        bag_bus_pairs.push_back(bag_BP);
        bag_bus_pairs_disjoint.push_back(bag_BP_disjoint);
    }

    unsigned  nb_arcs = 0;
    unsigned  nb_buses = 0;
    unsigned  nb_bus_pairs = 0;
    unsigned  nb_gens = 0;

    for (int c =0 ; c < nb_cliques; c++) {
        nb_buses += bag_bus_disjoint[c].size();
        nb_arcs += bag_arcs_disjoint[c].size();
        nb_gens += bag_gens_disjoint[c].size();
        nb_bus_pairs += bag_bus_pairs_disjoint[c]->_keys.size();
    }

    DebugOn("the number of total arcs: " <<  nb_arcs << endl);
    DebugOn("the number of total bus pairs: " <<  nb_bus_pairs << endl);
    DebugOn("the number of total buses: " <<  nb_buses << endl);
    DebugOn("the number of total gens: " <<  nb_gens << endl);
    // intersection_cliques
    for (auto a: cliquetree->arcs){
        std::vector<gravity::index_pair*> v3;
        std::set_intersection(bag_bus_pairs[a->_src->_id]->_keys.begin(), bag_bus_pairs[a->_src->_id]->_keys.end(),
                              bag_bus_pairs[a->_dest->_id]->_keys.begin(), bag_bus_pairs[a->_dest->_id]->_keys.end(),
                              back_inserter(v3), bus_pair_compare);
        if (v3.size() > 0) {
            a->_intersection_clique = v3;
        }
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
        var<Real>  bag_Wii("Wii_" + to_string(c), grid.w_min.in(bag_bus[c]), grid.w_max.in(bag_bus[c]));
        CLT.add_var(bag_Wii^(bag_bus[c].size()));
        bag_Wii.initialize_all(1.001);
        Wii.push_back(bag_Wii);

        if (bag_bus_pairs[c]->_keys.size() > 0){
            var<Real>  bag_R_Wij("R_Wij_"+ to_string(c), grid.wr_min.in(bag_bus_pairs[c]->_keys), grid.wr_max.in(bag_bus_pairs[c]->_keys));
            var<Real>  bag_Im_Wij("Im_Wij_"+ to_string(c), grid.wi_min.in(bag_bus_pairs[c]->_keys), grid.wi_max.in(bag_bus_pairs[c]->_keys));
            CLT.add_var(bag_R_Wij^(bag_bus_pairs[c]->_keys.size()));
            CLT.add_var(bag_Im_Wij^(bag_bus_pairs[c]->_keys.size()));
            bag_R_Wij.initialize_all(1.0);
            R_Wij.push_back(bag_R_Wij);
            Im_Wij.push_back(bag_Im_Wij);
        }
        else{
            var<Real> empty("empty");
            empty.set_size(0);
            R_Wij.push_back(empty);
            Im_Wij.push_back(empty);
        }
        
        if (bag_gens[c].size() > 0){
            var<Real>  bag_Pg("Pg_" + to_string(c), grid.pg_min.in(bag_gens[c]), grid.pg_max.in(bag_gens[c]));
            var<Real>  bag_Qg ("Qg_" + to_string(c), grid.qg_min.in(bag_gens[c]), grid.qg_max.in(bag_gens[c]));
            CLT.add_var(bag_Pg^(bag_gens_disjoint[c].size()));
            CLT.add_var(bag_Qg^(bag_gens_disjoint[c].size()));
            Pg.push_back(bag_Pg);
            Qg.push_back(bag_Qg);
        }
        else{
            var<Real> empty("empty");
            empty.set_size(0);
            Pg.push_back(empty);
            Qg.push_back(empty);
        }
    }

    /* Construct the objective function with generations bound constraints */
    func_ obj;
    for (int c = 0; c < nb_cliques; c++) {
        for (auto g:bag_gens_disjoint[c]) {
            if (g->_active) {
                obj += grid.c1(g->_name)*Pg[c](g->_name) + grid.c2(g->_name)*Pg[c](g->_name)*Pg[c](g->_name) + grid.c0(g->_name);
            }
        }
    }

    CLT.set_objective(min(obj));

    /** Define constraints */
    for (int c = 0; c < nb_cliques; c++) {
        if (bag_bus_pairs_disjoint[c]->_keys.size()>0) {
            Constraint SOC("SOC_" + to_string(c));
            SOC =  power(R_Wij[c].in(bag_bus_pairs_disjoint[c]->_keys), 2)
                   + power(Im_Wij[c].in(bag_bus_pairs_disjoint[c]->_keys), 2)
                   - Wii[c].from(bag_bus_pairs_disjoint[c]->_keys)*Wii[c].to(bag_bus_pairs_disjoint[c]->_keys) ;
            CLT.add_constraint(SOC <= 0);
        }
    }
   

    for (int c = 0; c < nb_cliques; c++) {
        //* Phase Angle Bounds constraints */
        if (bag_bus_pairs_disjoint[c]->_keys.size() > 0) {
            Constraint PAD_UB("PAD_UB" + to_string(c));
            PAD_UB += Im_Wij[c].in(bag_bus_pairs_disjoint[c]->_keys)- (grid.tan_th_max).in(bag_bus_pairs_disjoint[c]->_keys)*R_Wij[c].in(bag_bus_pairs_disjoint[c]->_keys);
            //PAD_UB -= R_Wij[c].in(bag_bus_pairs_disjoint[c]->_keys);
            //CLT.add_constraint(PAD_UB <= 0);

            Constraint PAD_LB("PAD_L  B" + to_string(c));
            PAD_LB = Im_Wij[c].in(bag_bus_pairs_disjoint[c]->_keys);
            PAD_LB -= (grid.tan_th_min).in(bag_bus_pairs_disjoint[c]->_keys)*R_Wij[c].in(bag_bus_pairs_disjoint[c]->_keys);
            CLT.add_constraint(PAD_LB >= 0);
        }
    }

//        if (bag_arcs_disjoint[c].size() > 0) {
//            /* Thermal Limit Constraints */
//            Constraint Thermal_Limit_from("Thermal_Limit_from" + to_string(c));
//            Thermal_Limit_from += power(grid.g_ff.in(bag_arcs[c])*Wii[c].from(bag_arcs[c])+ grid.g_ft.in(bag_arcs[c])*R_Wij[c].in_pairs(bag_arcs[c])
//                                        + grid.b_ft.in(bag_arcs[c])*Im_Wij[c].in_pairs(bag_arcs[c]), 2)
//                                  + power(grid.g_ft.in(bag_arcs[c])*Im_Wij[c].in_pairs(bag_arcs[c])-grid.b_ff.in(bag_arcs[c])*Wii[c].from(bag_arcs[c])
//                                          - grid.b_ft.in(bag_arcs[c])*R_Wij[c].in_pairs(bag_arcs[c]), 2);
//            Thermal_Limit_from -= power(grid.S_max.in(bag_arcs[c]),2);
//            CLT.add_constraint(Thermal_Limit_from <= 0);
//
//            Constraint Thermal_Limit_to("Thermal_Limit_to" + to_string(c));
//            Thermal_Limit_to += power(grid.g_tt.in(bag_arcs[c])*Wii[c].to(bag_arcs[c]) + grid.g_tf.in(bag_arcs[c])*R_Wij[c].in_pairs(bag_arcs[c])
//                                      + grid.b_tf.in(bag_arcs[c])*Im_Wij[c].in_pairs(bag_arcs[c]), 2)
//                                + power(grid.b_tt.in(bag_arcs[c])*Wii[c].to(bag_arcs[c]) + grid.b_tf.in(bag_arcs[c])*R_Wij[c].in_pairs(bag_arcs[c])
//                                        + grid.g_tf.in(bag_arcs[c])*Im_Wij[c].in_pairs(bag_arcs[c]), 2);
//
//            Thermal_Limit_to -= power(grid.S_max.in(bag_arcs[c]), 2);
//            CLT.add_constraint(Thermal_Limit_to <= 0);
//        }
  //  }

    //// COUPLING CONSTRAINTS
// for (auto b: grid.nodes) {
//        Bus* bus = (Bus*) b;
//        Constraint KCL_P("KCL_P"+bus->_name);
//        Constraint KCL_Q("KCL_Q"+bus->_name);
//        for (int c = 0; c < nb_cliques; c++) {
//            if (std::find(bag_bus_disjoint[c].begin(), bag_bus_disjoint[c].end(), bus) != bag_bus_disjoint[c].end()) {
//                if (bus->_has_gen){
//                    KCL_P += bus->pl()- sum(Pg[c].in(bus->_gen));;
//                    KCL_Q += bus->ql()- sum(Qg[c].in(bus->_gen));
//                }
//                else{
//                    KCL_P += bus->pl();
//                    KCL_Q += bus->ql();
//                }
//                KCL_P +=  bus->gs()*(Wii[c](bus->_name));
//                KCL_Q -=  bus->bs()*(Wii[c](bus->_name));
//            }
//            for (auto a: b->get_out()) {
//                if (std::find(bag_arcs_disjoint[c].begin(), bag_arcs_disjoint[c].end(), a) != bag_arcs_disjoint[c].end()) {
//                    KCL_P  += grid.g_ff(a->_name)*Wii[c](bus->_name)
//                              +grid.g_ft(a->_name)*R_Wij[c](a->_src->_name+","+a->_dest->_name)
//                              +grid.b_ft(a->_name)*Im_Wij[c](a->_src->_name+","+a->_dest->_name);
//
//                    KCL_Q  += -1*grid.b_ff(a->_name)*Wii[c](bus->_name)
//                              - grid.b_ft(a->_name)*R_Wij[c](a->_src->_name+","+a->_dest->_name)
//                              +grid.g_ft(a->_name)*Im_Wij[c](a->_src->_name+","+a->_dest->_name);
//                }
//            }
//
//            for (auto a: b->get_in()) {
//                if (std::find(bag_arcs_disjoint[c].begin(), bag_arcs_disjoint[c].end(), a) != bag_arcs_disjoint[c].end()) {
//                    KCL_P  += grid.g_tt(a->_name)*Wii[c](bus->_name)
//                              +grid.g_tf(a->_name)*R_Wij[c](a->_src->_name+","+a->_dest->_name)-grid.b_tf(a->_name)*Im_Wij[c](a->_src->_name+","+a->_dest->_name);
//
//                    KCL_Q  -= grid.b_tt(a->_name)*Wii[c](bus->_name)
//                              + grid.b_tf(a->_name)*R_Wij[c](a->_src->_name+","+a->_dest->_name)
//                              + grid.g_tf(a->_name)*Im_Wij[c](a->_src->_name+","+a->_dest->_name);
//                }
//            }
//        }
//        CLT.add_constraint(KCL_P = 0);
//        CLT.add_constraint(KCL_Q = 0);
//    }


//    for (auto a: cliquetree->arcs) {
//        Constraint Link_Wii("Link_Wii_" + to_string(a->_id));
//        Link_Wii += Wii[a->_src->_id].in(a->_intersection);
//        Link_Wii -= Wii[a->_dest->_id].in(a->_intersection);
//        CLT.add_constraint(Link_Wii = 0);
//
//        if (a->_intersection_clique.size() > 0) {
//            Constraint Link_Im_Wij("Link_Im_Wij_" + to_string(a->_id));
//            Link_Im_Wij += Im_Wij[a->_src->_id].in(a->_intersection_clique);
//            Link_Im_Wij -= Im_Wij[a->_dest->_id].in(a->_intersection_clique);
//            CLT.add_constraint(Link_Im_Wij = 0);
//
//            if (a->_intersection_clique.size() > 0) {
//                Constraint Link_R_Wij("Link_R_Wij_" + to_string(a->_id));
//                Link_R_Wij += R_Wij[a->_src->_id].in(a->_intersection_clique);
//                Link_R_Wij -= R_Wij[a->_dest->_id].in(a->_intersection_clique);
//                CLT.add_constraint(Link_R_Wij = 0);
//            }
//        }
//    }

    // solver SCOPF(CLT, cplex);
    solver SCOPF(CLT, ipopt);
    SCOPF.run();
}
