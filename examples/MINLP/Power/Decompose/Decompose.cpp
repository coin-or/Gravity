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
//#include "functions.h"
#include <thread>
#include <iomanip>

#ifdef USE_BOOST
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <deque>
#include <iterator>
#endif
struct net_param {
    param<Real> c0, c1, c2; /**< Generation costs */
    param<Real> tan_th_min, tan_th_max;
    param<Real> g_ff, g_ft, g_tt, g_tf, b_ff, b_ft, b_tf, b_tt;
    param<Real> S_max;
};

bool node_id_compare(const Node* n1, const Node* n2) {
    return n1->_id < n2->_id;
}

bool arc_id_compare(const Arc* n1, const Arc* n2) {
    return n1->_id < n2->_id;
}

bool bus_pair_compare(const index_pair* n1, const index_pair* n2) {
    return n1->_name < n2->_name;
}

/** INITIALISE SUBPROBLEM MODEL */
void sub1(net_param np, unsigned c, double rho, var<Real> Pg, var<Real> Qg, var<Real> Wii, var<Real> R_Wij, var<Real> Im_Wij,
          vector<gravity::index_pair*> bag_bus_pairs, vector<Bus*> bag_bus, vector<Bus*>  bag_bus_disjoint,
          vector<Line*> bag_arcs_disjoint, vector<Gen*> bag_gens_disjoint, vector<gravity::index_pair*> bag_bus_pairs_disjoint,
          vector<param<Real>> u, vector<param<Real>> Im_u, vector<param<Real>> R_u,
          param<Real>& Wii_log, param<Real>& R_Wij_log, param<Real>& Im_Wij_log,
          param<Real> Zii_log, param<Real> Im_Zij_log, param<Real> R_Zij_log, double& val)
{

    DebugOff("Solving subproblem associated with maximal clique, "<< c << endl);
    Model Subr("Subr");
    Subr.add_var(Pg);
    Subr.add_var(Qg);
    Subr.add_var(Wii);
    Subr.add_var(R_Wij);
    Subr.add_var(Im_Wij);

    /* Construct the objective function*/
    func_ obj;
    for (auto g:bag_gens_disjoint) {
        if (g->_active) {
            obj += np.c1(g->_name)*Pg(g->_name)+ np.c2(g->_name)*Pg(g->_name)*Pg(g->_name)+np.c0(g->_name);
        }
    }

    // THE QUADRATIC TERM
    for (auto nn: bag_bus) {
        obj += rho*power(Wii(nn->_name) - Zii_log(nn->_name).getvalue() + u[c](nn->_name).getvalue(), 2) - rho*std::pow(u[c](nn->_name).getvalue(), 2);
    }

    for (auto pair: bag_bus_pairs) {
        obj += rho*power(R_Wij(pair->_name) - R_Zij_log(pair->_name).getvalue() + R_u[c](pair->_name).getvalue(), 2) - rho*std::pow(R_u[c](pair->_name).getvalue(), 2);
        obj += rho*power(Im_Wij(pair->_name)- Im_Zij_log(pair->_name).getvalue()+Im_u[c](pair->_name).getvalue(), 2)- rho*std::pow(Im_u[c](pair->_name).getvalue(), 2);
    }

    Subr.set_objective(min(obj));

    if (bag_bus_pairs_disjoint.size() > 0) {
        Constraint SOC("SOC_" + to_string(c));
        SOC =  power(R_Wij.in(bag_bus_pairs_disjoint), 2)
               + power(Im_Wij.in(bag_bus_pairs_disjoint), 2)
               - Wii.from(bag_bus_pairs_disjoint)*Wii.to(bag_bus_pairs_disjoint);
        Subr.add_constraint(SOC <= 0);
    }

    //* KCL */
    for (auto bus:  bag_bus_disjoint) {
        Constraint KCL_P("KCL_P"+bus->_name);
        Constraint KCL_Q("KCL_Q"+bus->_name);

        for (auto a: bus->get_out()) {
            KCL_P += np.g_ff(a->_name)*Wii(a->_src->_name)
                     +np.g_ft(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name)
                     +np.b_ft(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);

            KCL_Q += -1*np.b_ff(a->_name)*Wii(a->_src->_name)
                     - np.b_ft(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name)
                     +np.g_ft(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);
        }

        for (auto a: bus->get_in()) {
            KCL_P  += np.g_tt(a->_name)*Wii(a->_dest->_name)
                      +np.g_tf(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name)
                      -np.b_tf(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);

            KCL_Q  -= np.b_tt(a->_name)*Wii(a->_dest->_name)
                      + np.b_tf(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name)
                      + np.g_tf(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);
        }
        if(bus->_has_gen) {
            KCL_P += bus->pl()- sum(Pg.in(bus->_gen));
            KCL_Q += bus->ql()- sum(Qg.in(bus->_gen));
        }
        else {
            KCL_P += bus->pl();
            KCL_Q += bus->ql();
        }

        /* Shunts */
        KCL_P += bus->gs()*(Wii(bus->_name));
        KCL_Q -= bus->bs()*(Wii(bus->_name));

        Subr.add_constraint(KCL_P = 0);
        Subr.add_constraint(KCL_Q = 0);
    }

    /* Phase Angle Bounds constraints */
    if (bag_bus_pairs_disjoint.size() > 0) {
        Constraint PAD_UB("PAD_UB" + to_string(c));
        PAD_UB = Im_Wij.in(bag_bus_pairs_disjoint);
        PAD_UB -= (np.tan_th_max).in(bag_bus_pairs_disjoint)*R_Wij.in(bag_bus_pairs_disjoint);
        Subr.add_constraint(PAD_UB <= 0);

        Constraint PAD_LB("PAD_LB" + to_string(c));
        PAD_LB = Im_Wij.in(bag_bus_pairs_disjoint);
        PAD_LB -= (np.tan_th_min).in(bag_bus_pairs_disjoint)*R_Wij.in(bag_bus_pairs_disjoint);
        Subr.add_constraint(PAD_LB >= 0);
    }

    if (bag_arcs_disjoint.size() > 0) {
        /* Thermal Limit Constraints */
        Constraint Thermal_Limit_from("Thermal_Limit_from"+to_string(c));
        Thermal_Limit_from += power(np.g_ff.in(bag_arcs_disjoint)*Wii.from(bag_arcs_disjoint)+
                                    np.g_ft.in(bag_arcs_disjoint)*R_Wij.in_pairs(bag_arcs_disjoint)
                                    +np.b_ft.in(bag_arcs_disjoint)*Im_Wij.in_pairs(bag_arcs_disjoint), 2)
                              + power(np.b_ff.in(bag_arcs_disjoint)*Wii.from(bag_arcs_disjoint)
                                      +np.b_ft.in(bag_arcs_disjoint)*R_Wij.in_pairs(bag_arcs_disjoint)
                                      -np.g_ft.in(bag_arcs_disjoint)*Im_Wij.in_pairs(bag_arcs_disjoint), 2);

        Thermal_Limit_from -= power(np.S_max.in(bag_arcs_disjoint), 2);
        Subr.add_constraint(Thermal_Limit_from <= 0);


        Constraint Thermal_Limit_to("Thermal_Limit_to" + to_string(c));
        Thermal_Limit_to += power(np.g_tt.in(bag_arcs_disjoint)*Wii.to(bag_arcs_disjoint)
                                  + np.g_tf.in(bag_arcs_disjoint)*R_Wij.in_pairs(bag_arcs_disjoint)
                                  + np.b_tf.in(bag_arcs_disjoint)*Im_Wij.in_pairs(bag_arcs_disjoint), 2)
                            + power(np.b_tt.in(bag_arcs_disjoint)*Wii.to(bag_arcs_disjoint)
                                    + np.b_tf.in(bag_arcs_disjoint)*R_Wij.in_pairs(bag_arcs_disjoint)
                                    + np.g_tf.in(bag_arcs_disjoint)*Im_Wij.in_pairs(bag_arcs_disjoint), 2);

        Thermal_Limit_to -= power(np.S_max.in(bag_arcs_disjoint), 2);
        Subr.add_constraint(Thermal_Limit_to <= 0);
    }

    /* solve it! */
    solver solve_Subr(Subr, ipopt);
    solve_Subr.run();

    // collect the primal values.
    Wii_log   = (*(var<Real>*)  Subr.get_var("Wii_"+ to_string(c)));
    R_Wij_log  = (*(var<Real>*) Subr.get_var("R_Wij_" + to_string(c)));
    Im_Wij_log = (*(var<Real>*) Subr.get_var("Im_Wij_"+ to_string(c)));

    val = Subr._obj_val;
}


/** INITIALISE SUBPROBLEM MODEL */
void subproblem(net_param np, Net* cliquetree, unsigned c, var<Real> Pg, var<Real> Qg,
                var<Real> Wii, var<Real> R_Wij, var<Real> Im_Wij,
                vector<Bus*>  bag_bus_disjoint, vector<Line*> bag_arcs_disjoint, vector<Gen*> bag_gens_disjoint,
                vector<gravity::index_pair*> bag_bus_pairs_disjoint,
                vector<param<Real>> R_lambda_sep, vector<param<Real>> Im_lambda_sep, vector<param<Real>> lambda_sep,
                param<Real>& Wii_log, param<Real>& R_Wij_log, param<Real>& Im_Wij_log, double& val)
{
    DebugOff("Solving subproblem associated with maximal clique, "<< c << endl);
    Model Subr("Subr");
    Subr.add_var(Pg);
    Subr.add_var(Qg);
    Subr.add_var(Wii);
    Subr.add_var(R_Wij);
    Subr.add_var(Im_Wij);

    /* Construct the objective function*/
    func_ obj;
    for (auto g:bag_gens_disjoint) {
        if (g->_active) {
            obj += np.c1(g->_name)*Pg(g->_name)+ np.c2(g->_name)*Pg(g->_name)*Pg(g->_name)+np.c0(g->_name);
        }
    }

    Node* bag = cliquetree->get_node(to_string(c));
    for (auto arc: bag->get_out()) {
        for (auto nn: arc->_intersection) {
            obj += Wii(nn->_name)*lambda_sep[arc->_id](nn->_name);
        }

        for (auto pair: arc->_intersection_clique) {
            obj += R_Wij(pair->_name)*R_lambda_sep[arc->_id](pair->_name);
            obj += Im_Wij(pair->_name)*Im_lambda_sep[arc->_id](pair->_name);
        }

    }

    for (auto arc: bag->get_in()) {
        for (auto nn: arc->_intersection) {
            obj -= Wii(nn->_name)*lambda_sep[arc->_id](nn->_name);
        }

        for (auto pair: arc->_intersection_clique) {
            obj -= R_Wij(pair->_name)*R_lambda_sep[arc->_id](pair->_name);
            obj -= Im_Wij(pair->_name)*Im_lambda_sep[arc->_id](pair->_name);
        }
    }

    Subr.set_objective(min(obj));

    if (bag_bus_pairs_disjoint.size()>0) {
        Constraint SOC("SOC_" + to_string(c));
        SOC =  power(R_Wij.in(bag_bus_pairs_disjoint), 2)
               + power(Im_Wij.in(bag_bus_pairs_disjoint), 2)
               - Wii.from(bag_bus_pairs_disjoint)*Wii.to(bag_bus_pairs_disjoint) ;
        Subr.add_constraint(SOC <= 0);
    }

    //* KCL */
    for (auto bus:  bag_bus_disjoint) {
        Constraint KCL_P("KCL_P"+bus->_name);
        Constraint KCL_Q("KCL_Q"+bus->_name);

        for (auto a: bus->get_out()) {
            KCL_P += np.g_ff(a->_name)*Wii(a->_src->_name)
                     +np.g_ft(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name)
                     +np.b_ft(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);

            KCL_Q += -1*np.b_ff(a->_name)*Wii(a->_src->_name)
                     - np.b_ft(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name)
                     +np.g_ft(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);
        }

        for (auto a: bus->get_in()) {
            KCL_P  += np.g_tt(a->_name)*Wii(a->_dest->_name)
                      +np.g_tf(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name)
                      -np.b_tf(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);

            KCL_Q  -= np.b_tt(a->_name)*Wii(a->_dest->_name)
                      + np.b_tf(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name)
                      + np.g_tf(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);
        }
        if(bus->_has_gen) {
            KCL_P += bus->pl()- sum(Pg.in(bus->_gen));
            KCL_Q += bus->ql()- sum(Qg.in(bus->_gen));
        }
        else {
            KCL_P += bus->pl();
            KCL_Q += bus->ql();
        }

        /* Shunts */
        KCL_P += bus->gs()*(Wii(bus->_name));
        KCL_Q -= bus->bs()*(Wii(bus->_name));

        Subr.add_constraint(KCL_P = 0);
        Subr.add_constraint(KCL_Q = 0);
    }

    /* Phase Angle Bounds onstraints */
    if (bag_bus_pairs_disjoint.size() > 0) {
        Constraint PAD_UB("PAD_UB" + to_string(c));
        PAD_UB = Im_Wij.in(bag_bus_pairs_disjoint);
        PAD_UB -= (np.tan_th_max).in(bag_bus_pairs_disjoint)*R_Wij.in(bag_bus_pairs_disjoint);
        Subr.add_constraint(PAD_UB <= 0);

        Constraint PAD_LB("PAD_LB" + to_string(c));
        PAD_LB = Im_Wij.in(bag_bus_pairs_disjoint);
        PAD_LB -= (np.tan_th_min).in(bag_bus_pairs_disjoint)*R_Wij.in(bag_bus_pairs_disjoint);
        Subr.add_constraint(PAD_LB >= 0);
    }

    if (bag_arcs_disjoint.size() > 0) {
        /* Thermal Limit Constraints */
        Constraint Thermal_Limit_from("Thermal_Limit_from"+to_string(c));
        Thermal_Limit_from += power(np.g_ff.in(bag_arcs_disjoint)*Wii.from(bag_arcs_disjoint)+
                                    np.g_ft.in(bag_arcs_disjoint)*R_Wij.in_pairs(bag_arcs_disjoint)
                                    +np.b_ft.in(bag_arcs_disjoint)*Im_Wij.in_pairs(bag_arcs_disjoint), 2)
                              + power(np.b_ff.in(bag_arcs_disjoint)*Wii.from(bag_arcs_disjoint)
                                      +np.b_ft.in(bag_arcs_disjoint)*R_Wij.in_pairs(bag_arcs_disjoint)
                                      -np.g_ft.in(bag_arcs_disjoint)*Im_Wij.in_pairs(bag_arcs_disjoint), 2);

        Thermal_Limit_from -= power(np.S_max.in(bag_arcs_disjoint), 2);
        Subr.add_constraint(Thermal_Limit_from <= 0);


        Constraint Thermal_Limit_to("Thermal_Limit_to" + to_string(c));
        Thermal_Limit_to += power(np.g_tt.in(bag_arcs_disjoint)*Wii.to(bag_arcs_disjoint)
                                  + np.g_tf.in(bag_arcs_disjoint)*R_Wij.in_pairs(bag_arcs_disjoint)
                                  + np.b_tf.in(bag_arcs_disjoint)*Im_Wij.in_pairs(bag_arcs_disjoint), 2)
                            + power(np.b_tt.in(bag_arcs_disjoint)*Wii.to(bag_arcs_disjoint)
                                    + np.b_tf.in(bag_arcs_disjoint)*R_Wij.in_pairs(bag_arcs_disjoint)
                                    + np.g_tf.in(bag_arcs_disjoint)*Im_Wij.in_pairs(bag_arcs_disjoint), 2);

        Thermal_Limit_to -= power(np.S_max.in(bag_arcs_disjoint), 2);
        Subr.add_constraint(Thermal_Limit_to <= 0);
    }
    /* solve it! */
    solver solve_Subr(Subr, ipopt);
    solve_Subr.run();

    // collect the primal values.
    Wii_log   = (*(var<Real>*) Subr.get_var("Wii_"+ to_string(c)));
    R_Wij_log  = (*(var<Real>*) Subr.get_var("R_Wij_" + to_string(c)));
    Im_Wij_log = (*(var<Real>*) Subr.get_var("Im_Wij_"+ to_string(c)));

    val = Subr._obj_val;
//    cout << "val: " << val << endl;
}

int inout (PowerNet& grid) {
    cout << "////////////////////////////////////////" << endl;
    Net* grid_augment = grid.clone_undirected();
    DebugOn("number of edges of the augmented graph, " << grid_augment->arcs.size() << endl);
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

// NOTE that grid_augment may NOT be chordal.
    grid_augment->get_tree_decomp_bags();
    auto cliquetree = grid_augment->get_clique_tree_prim(); // Note: it builds on grid_augment.

    const unsigned nb_cliques = grid_augment->_bags.size();

    vector<net_param> bag_param;
    vector<vector<Bus*>> bag_bus;
    vector<vector<Bus*>> bag_bus_disjoint;
    vector<vector<Gen*>> bag_gens;
    vector<vector<Gen*>> bag_gens_disjoint;
    vector<vector<Line*>> bag_arcs; // bag_arcs contains the arcs of the power grid while variables associated with W are defined on chordal graph.
    vector<vector<Line*>> bag_arcs_disjoint;
    vector<vector<gravity::index_pair*>> bag_bus_pairs; // bus_pairs in each bag.
    vector<vector<gravity::index_pair*> > bag_bus_pairs_disjoint; // bus_pairs in each bag.

    map<string, unsigned> indexij;
    map<Bus*, unsigned> indexii;

    for (int c = 0; c < nb_cliques; c++) {
        net_param np;
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
        np.b_ff=grid.b_ff;
        np.b_ft=grid.b_ft;
        np.b_tf=grid.b_tf;
        np.b_tt=grid.b_tt;
        np.g_ff=grid.g_ff;
        np.g_ft=grid.g_ft;
        np.g_tf=grid.g_tf;
        np.g_tt=grid.g_tt;
        np.c0=grid.c0;
        np.c1=grid.c1;
        np.c2=grid.c2;
        np.tan_th_min=grid.tan_th_min;
        np.tan_th_max=grid.tan_th_max;
        np.S_max=grid.S_max;
        bag_param.push_back(np);
    }


    for (int c = 0; c < nb_cliques; c++) {
        vector<Bus*> bag_B_disjoint;
        vector<Gen*> bag_G_disjoint;

        vector<Line*> bag_A;
        vector<Line*> bag_A_disjoint;

        vector<gravity::index_pair*> bag_BP;
        vector<gravity::index_pair*> bag_BP_disjoint;

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
                        bag_BP.push_back(new index_pair(index_(b->_name), index_(bag_bus[c].at(j)->_name), true));
                    }
                    else {
                        bag_BP.push_back(new index_pair(index_(bag_bus[c].at(j)->_name),index_(b->_name), true));
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
                            bag_BP_disjoint.push_back(new index_pair(index_(b->_name), index_(bag_bus[c].at(j)->_name), true));
                        }
                        else {
                            indexij.insert(make_pair<>(key_inv, c));
                            bag_BP_disjoint.push_back(new index_pair(index_(bag_bus[c].at(j)->_name),index_(b->_name), true));
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

        sort(bag_BP.begin(), bag_BP.end(), bus_pair_compare);
        sort(bag_BP_disjoint.begin(), bag_BP_disjoint.end(), bus_pair_compare);
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
        nb_bus_pairs += bag_bus_pairs_disjoint[c].size();
    }

    DebugOn("the number of total arcs: " <<  nb_arcs << endl);
    DebugOn("the number of total bus pairs: " <<  nb_bus_pairs << endl);
    DebugOn("the number of total buses: " <<  nb_buses << endl);
    DebugOn("the number of total gens: " <<  nb_gens << endl);

    // intersection_cliques
    for (auto a: cliquetree->arcs) {
        std::vector<gravity::index_pair*> v3;
        std::set_intersection(bag_bus_pairs[a->_src->_id].begin(), bag_bus_pairs[a->_src->_id].end(),
                              bag_bus_pairs[a->_dest->_id].begin(), bag_bus_pairs[a->_dest->_id].end(),
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
        var<Real>  bag_Wii("Wii_"+ to_string(c), grid.w_min.in(bag_bus[c]), grid.w_max.in(bag_bus[c]));
        CLT.add_var(bag_Wii^(bag_bus[c].size()));
        bag_Wii.initialize_all(1.001);
        Wii.push_back(bag_Wii);

        if (bag_bus_pairs[c].size() > 0) {
            var<Real>  bag_R_Wij("R_Wij_" + to_string(c), grid.wr_min.in(bag_bus_pairs[c]), grid.wr_max.in(bag_bus_pairs[c]));
            var<Real>  bag_Im_Wij("Im_Wij_" + to_string(c), grid.wi_min.in(bag_bus_pairs[c]), grid.wi_max.in(bag_bus_pairs[c]));
            CLT.add_var(bag_R_Wij^(bag_bus_pairs[c].size()));
            CLT.add_var(bag_Im_Wij^(bag_bus_pairs[c].size()));
            bag_R_Wij.initialize_all(1.0);
            R_Wij.push_back(bag_R_Wij);
            Im_Wij.push_back(bag_Im_Wij);
        }
        else {
            var<Real> empty1("R_Wij_"+to_string(c));
            var<Real> empty2("Im_Wij_"+to_string(c));
            empty1.set_size(0);
            empty2.set_size(0);
            R_Wij.push_back(empty1);
            Im_Wij.push_back(empty2);
        }

        if (bag_gens_disjoint[c].size() > 0) {
            var<Real>  bag_Pg("Pg_" + to_string(c), grid.pg_min.in(bag_gens_disjoint[c]), grid.pg_max.in(bag_gens_disjoint[c]));
            var<Real>  bag_Qg ("Qg_" + to_string(c), grid.qg_min.in(bag_gens_disjoint[c]), grid.qg_max.in(bag_gens_disjoint[c]));
            CLT.add_var(bag_Pg^(bag_gens_disjoint[c].size()));
            CLT.add_var(bag_Qg^(bag_gens_disjoint[c].size()));
            Pg.push_back(bag_Pg);
            Qg.push_back(bag_Qg);
        }
        else {
            var<Real> empty("empty");
            empty.set_size(0);
            Pg.push_back(empty);
            Qg.push_back(empty);
        }
    }

///////////////// DEFINE LAGRANGE MULTIPLIERS  ////////////////////////////////
    vector<param<Real>> R_lambda_in;
    vector<param<Real>> Im_lambda_in;
    vector<param<Real>> lambda_in;

    vector<param<Real>> R_lambda_out;
    vector<param<Real>> Im_lambda_out;
    vector<param<Real>> lambda_out;

    vector<param<Real>> R_lambda_sep;
    vector<param<Real>> Im_lambda_sep;
    vector<param<Real>> lambda_sep;

    vector<param<Real>> R_lambda_grad;
    vector<param<Real>> Im_lambda_grad;
    vector<param<Real>> lambda_grad;

    for (auto a: cliquetree->arcs) {
        auto l = a->_id;
        param<Real> lambda_arc_in("R_lambda_arc_in" + to_string(l));
        param<Real> R_lambda_arc_in("R_lambda_arc_in" + to_string(l));
        param<Real> Im_lambda_arc_in("Im_lambda_arc_in" + to_string(l));
        lambda_arc_in^(a->_weight);
        R_lambda_arc_in^(a->_intersection_clique.size());
        Im_lambda_arc_in^(a->_intersection_clique.size());

        param<Real> R_lambda_arc_out("R_lambda_arc_out" + to_string(l));
        param<Real> Im_lambda_arc_out("Im_lambda_arc_out" + to_string(l));
        param<Real> lambda_arc_out("lambda_arc_out" + to_string(l));
        R_lambda_arc_out^(a->_intersection_clique.size());
        Im_lambda_arc_out^(a->_intersection_clique.size());
        lambda_arc_out^(a->_weight);

        param<Real> R_lambda_arc_sep("R_lambda_arc_sep" + to_string(l));
        param<Real> Im_lambda_arc_sep("Im_lambda_arc_sep" + to_string(l));
        param<Real> lambda_arc_sep("lambda_arc_sep" + to_string(l));
        R_lambda_arc_sep^(a->_intersection_clique.size());
        Im_lambda_arc_sep^(a->_intersection_clique.size());
        lambda_arc_sep^(a->_intersection.size());

        param<Real> R_lambda_arc_grad("R_lambda_arc_grad" + to_string(l));
        param<Real> Im_lambda_arc_grad("Im_lambda_arc_grad" + to_string(l));
        param<Real> lambda_arc_grad("lambda_arc_grad" + to_string(l));
        R_lambda_arc_grad^(a->_intersection_clique.size());
        Im_lambda_arc_grad^(a->_intersection_clique.size());
        lambda_arc_grad^(a->_weight);


        R_lambda_arc_in.initialize_all(0);
        Im_lambda_arc_in.initialize_all(0);
        lambda_arc_in.initialize_all(0);

        R_lambda_arc_out.initialize_all(0);
        Im_lambda_arc_out.initialize_all(0);
        lambda_arc_out.initialize_all(0);

        R_lambda_arc_sep.initialize_all(0);
        Im_lambda_arc_sep.initialize_all(0);
        lambda_arc_sep.initialize_all(0);

        R_lambda_in.push_back(R_lambda_arc_in);
        Im_lambda_in.push_back(Im_lambda_arc_in);
        lambda_in.push_back(lambda_arc_in);

        R_lambda_out.push_back(R_lambda_arc_out);
        Im_lambda_out.push_back(Im_lambda_arc_out);
        lambda_out.push_back(lambda_arc_out);


        lambda_sep.push_back(lambda_arc_sep);
        R_lambda_sep.push_back(R_lambda_arc_sep);
        Im_lambda_sep.push_back(Im_lambda_arc_sep);

        R_lambda_grad.push_back(R_lambda_arc_grad);
        Im_lambda_grad.push_back(Im_lambda_arc_grad);
        lambda_grad.push_back(lambda_arc_grad);
    }

/////////////////////////////////// Master Problem ///////////////////////////////////
    Model Master("Master");
    /** param **/
    param<Real> gamma_in("gamma_C_in");
    param<Real> gamma_out("gamma_C_out");
    param<Real> gamma_sep("gamma_C_sep");
    gamma_in^nb_cliques;
    gamma_out^nb_cliques;
    gamma_sep^nb_cliques;

    /** Variables  */
    var<Real> gamma_C("gamma_C");
    Master.add_var(gamma_C^nb_cliques);

    vector<var<Real>> R_lambda_var;
    vector<var<Real>> Im_lambda_var;
    vector<var<Real>> lambda_var;
    for (auto a: cliquetree->arcs) {
        var<Real> lambda("lambda_arc_" + to_string(a->_id));
        var<Real> R_lambda("R_lambda_arc_" + to_string(a->_id));
        var<Real> Im_lambda("Im_lambda_arc_" + to_string(a->_id));
        Master.add_var(lambda^(a->_weight));
        Master.add_var(R_lambda^(a->_intersection_clique.size()));
        Master.add_var(Im_lambda^(a->_intersection_clique.size()));

        lambda_var.push_back(lambda);
        R_lambda_var.push_back(R_lambda);
        Im_lambda_var.push_back(Im_lambda);
    }

    /////////** OBJ*//////////////
    func_ master_obj;
    master_obj += sum(gamma_C);
    Master.set_objective(max(master_obj));
    double bound = 100000;

    Constraint UB;
    UB = sum(gamma_C) - bound;
    Master.add_constraint(UB <= 0);

////////////////  CONVERGENCE INFORMATION /////////////////////////
    unsigned iter_limit;
    cout << "Enter the limit of the number of iterations: ";
    cin >> iter_limit;
    cout << endl;

    double alpha = .5;
    double LBlog[iter_limit];
    double UBlog[iter_limit];

    for(int iter = 0; iter < iter_limit; iter++) {
        LBlog[iter] = 0.0;
        UBlog[iter] = 0.0;
    }

    double LDlog[iter_limit];

    // LOG OF SOLUTIONS
    // Log here means the previous primal and dual solution.
    vector<param<Real>> R_lambda_log;
    vector<param<Real>> Im_lambda_log;
    vector<param<Real>> lambda_log;

    vector<param<Real>> R_Wij_log;
    vector<param<Real>> Im_Wij_log;
    vector<param<Real>> Wii_log;

    for (auto a: cliquetree->arcs) {
        int l = a->_id;
        param<Real> R_lambda_C_log("R_lambda_C_log" + to_string(l));
        param<Real> Im_lambda_C_log("Im_lambda_C_log" + to_string(l));
        param<Real> lambda_C_log("Im_lambda_C_log" + to_string(l));

        lambda_C_log^(a->_weight);
        R_lambda_C_log^(a->_intersection_clique.size());
        Im_lambda_C_log^(a->_intersection_clique.size());

        lambda_log.push_back(lambda_C_log);
        R_lambda_log.push_back(R_lambda_C_log);
        Im_lambda_log.push_back(Im_lambda_C_log);
    }

    for (int c = 0; c < nb_cliques; c++) {
        param<Real> Wii_C_log("Wii_C_log" + to_string(c));
        param<Real> Im_Wij_C_log("Im_Wij_C_log" + to_string(c));
        param<Real> R_Wij_C_log("R_Wij_C_log" + to_string(c));
        Wii_C_log^(bag_bus[c].size());
        R_Wij_C_log^(bag_bus_pairs[c].size());
        Im_Wij_C_log^(bag_bus_pairs[c].size());

        Wii_log.push_back(Wii_C_log);
        R_Wij_log.push_back(R_Wij_C_log);
        Im_Wij_log.push_back(Im_Wij_C_log);
    }

///////////////////////////////// INITIALIZATION ///////////////////////////////////////////
    double wall0 = get_wall_time();
    double cpu0  = get_cpu_time();
    double dual = 0;
    double value_dual[nb_cliques];
    //std::vector<std::thread> threads;
    for (int c = 0; c < nb_cliques; c++) {
        subproblem(bag_param[c], cliquetree, c, Pg[c], Qg[c], Wii[c], R_Wij[c], Im_Wij[c],
                                 bag_bus_disjoint[c], bag_arcs[c], bag_gens_disjoint[c], bag_bus_pairs_disjoint[c],
                                  R_lambda_sep,  Im_lambda_sep, lambda_sep, (Wii_log[c]),
                                  (R_Wij_log[c]), (Im_Wij_log[c]), (value_dual[c]));
        //threads.push_back(std::thread(subproblem, bag_param[c], cliquetree, c, Pg[c], Qg[c], Wii[c], R_Wij[c], Im_Wij[c],
        //                              bag_bus_disjoint[c], bag_arcs[c], bag_gens_disjoint[c], bag_bus_pairs_disjoint[c],
        //                              R_lambda_sep,  Im_lambda_sep, lambda_sep, std::ref(Wii_log[c]),
        //                              std::ref(R_Wij_log[c]), std::ref(Im_Wij_log[c]), std::ref(value_dual[c])));
    }
    //
// join the threads with the main thread
    //for(auto &t: threads) {
    //    t.join();
    //}

    for (int c = 0; c < nb_cliques; c++) {
        dual += value_dual[c] ;
    }

    cout << "Initialization_value,   " << dual <<endl;

    // initialise the in values.
    for (int i= 0; i < nb_cliques; i++) {
        gamma_in(i) = value_dual[i];
        cout << "value_dual[i]: " << value_dual[i]<< endl;
    }

/////////////////// APPEND MORE CONSTRAINTS TO MAIN //////////////////////////////////
    if (iter_limit > 0) {
        for (auto bag: cliquetree->nodes) {
            unsigned c = bag->_id;
            Constraint Concavity("Iter_0_Concavity_" + to_string(c));
            Concavity += gamma_C(c);
            Concavity -= value_dual[c];
            for (auto arc: bag->get_out()) {
                for (auto nn: arc->_intersection) {
                    Concavity -= (lambda_var[arc->_id](nn->_name)-lambda_sep[arc->_id](nn->_name).getvalue())*Wii_log[c](nn->_name).getvalue();
                    DebugOff("Wii_log, " << Wii_log[c](nn->_name).getvalue()<< endl);
                }

                for (auto pair: arc->_intersection_clique) {
                    Concavity -= (R_lambda_var[arc->_id](pair->_name)- R_lambda_sep[arc->_id](pair->_name).getvalue())*R_Wij_log[c](pair->_name).getvalue();
                    Concavity -= (Im_lambda_var[arc->_id](pair->_name)-Im_lambda_sep[arc->_id](pair->_name).getvalue())*Im_Wij_log[c](pair->_name).getvalue();
                    DebugOff("R_Wij_log, " << R_Wij_log[c](pair->_name).getvalue()<< endl);
                    DebugOff("Im_Wij_log, " << Im_Wij_log[c](pair->_name).getvalue()<< endl);
                }
            }

            for (auto arc: bag->get_in()) {
                for (auto nn: arc->_intersection) {
                    Concavity += (lambda_var[arc->_id](nn->_name)-lambda_sep[arc->_id](nn->_name).getvalue())*Wii_log[c](nn->_name).getvalue();
                }
                for (auto pair: arc->_intersection_clique) {
                    Concavity += (R_lambda_var[arc->_id](pair->_name)- R_lambda_sep[arc->_id](pair->_name).getvalue())*R_Wij_log[c](pair->_name).getvalue();
                    Concavity += (Im_lambda_var[arc->_id](pair->_name)-Im_lambda_sep[arc->_id](pair->_name).getvalue())*Im_Wij_log[c](pair->_name).getvalue();
                }
            }
            Master.add_constraint(Concavity <= 0);
        }
    }
    solver solve_Master(Master, ipopt);
    //solver solve_Master(Master, cplex);
    solve_Master.run();

    // initialise the outer point.
    gamma_out = (*(var<Real>*) Master.get_var("gamma_C"));

    for (auto a: cliquetree->arcs) {
        lambda_out[a->_id] = (*(var<Real>*) Master.get_var("lambda_arc_"+ to_string(a->_id)));
        R_lambda_out[a->_id] = (*(var<Real>*) Master.get_var("R_lambda_arc_" + to_string(a->_id)));
        Im_lambda_out[a->_id] = (*(var<Real>*) Master.get_var("Im_lambda_arc_"+ to_string(a->_id)));
    }

    cout << "................  Initialization of Master problem ....................."  <<endl;
    cout << "value: " << Master._obj_val  <<endl;

////////////////////////// BEGIN LAGRANGE ITERATIONS HERE /////////////////////////////////////
    cout << "<<<<<<<<<<< Lagrangian decomposition algorithm >>>>>>>>>"<< endl;
    cout<< setw(15) << left <<"ITERATION" << setw(15) << "LB" << setw(15)  << "UB" << endl;
    for(int itcount = 1; itcount < iter_limit  ; itcount++) {
        //////// CONSTRUCT SEPARATION POINTS
        for (int c = 0; c < nb_cliques; c++) {
            gamma_sep(c) = alpha*gamma_out(c).getvalue() + (1 - alpha)*gamma_in(c).getvalue();
        }

        for (auto a: cliquetree->arcs) {
            int l = a->_id;
            for (auto nn: a->_intersection) {
                lambda_sep[l](nn->_name) = alpha*lambda_out[l](nn->_name).getvalue()
                                           + (1 - alpha)*lambda_in[l](nn->_name).getvalue();
                DebugOff("lambda_sep, " << lambda_sep[a->_id](nn->_name).getvalue() << endl);
            }

            for (auto pair: a->_intersection_clique) {
                R_lambda_sep[l](pair->_name) = alpha*R_lambda_out[l](pair->_name).getvalue()
                                               + (1 - alpha)*R_lambda_in[l](pair->_name).getvalue();
                Im_lambda_sep[l](pair->_name)= alpha*Im_lambda_out[l](pair->_name).getvalue()
                                               + (1 - alpha)*Im_lambda_in[l](pair->_name).getvalue();
            }
        }

        /////////////////SLOVE SUBPROBLEM/////////////////////

        double value_dual[nb_cliques];
        //std::vector<std::thread> threads;
        for (int c = 0; c < nb_cliques; c++) {
             subproblem(bag_param[c], cliquetree, c, Pg[c], Qg[c], Wii[c], R_Wij[c], Im_Wij[c],
                                      bag_bus_disjoint[c], bag_arcs[c], bag_gens_disjoint[c], bag_bus_pairs_disjoint[c],
                                       R_lambda_sep,  Im_lambda_sep, lambda_sep, (Wii_log[c]),
                                       (R_Wij_log[c]), (Im_Wij_log[c]), (value_dual[c]));
            //threads.push_back(std::thread(subproblem, bag_param[c], cliquetree, c, Pg[c], Qg[c], Wii[c], R_Wij[c], Im_Wij[c],
            //                              bag_bus_disjoint[c], bag_arcs[c], bag_gens_disjoint[c], bag_bus_pairs_disjoint[c],
            //                              R_lambda_sep,  Im_lambda_sep, lambda_sep, std::ref(Wii_log[c]),
            //                              std::ref(R_Wij_log[c]), std::ref(Im_Wij_log[c]), std::ref(value_dual[c])));
        }
// join the threads with the main thread
        //for(auto &t: threads) {
        //    t.join();
        //}
        dual =0;
        for (int c = 0; c < nb_cliques; c++) {
            dual += value_dual[c] ;
        }

// UPDATE POINTS of Kelly using in-out algorithm (Ben-Ameur and Neto)
        if (dual- sum(gamma_sep).eval() < 0) {
            for (auto bag: cliquetree->nodes) {
                unsigned c = bag->_id;
                Constraint Concavity("Iter_" + to_string(itcount) + "_Concavity_" + to_string(c));
                Concavity += gamma_C(c);
                Concavity -= value_dual[c];
                for (auto arc: bag->get_out()) {
                    for (auto nn: arc->_intersection) {
                        Concavity -= (lambda_var[arc->_id](nn->_name)-lambda_sep[arc->_id](nn->_name).getvalue())*Wii_log[c](nn->_name).getvalue();
                        DebugOff("Wii_log, " << Wii_log[c](nn->_name).getvalue()<< endl);
                    }
                    for (auto pair: arc->_intersection_clique) {
                        Concavity -= (R_lambda_var[arc->_id](pair->_name)- R_lambda_sep[arc->_id](pair->_name).getvalue())*R_Wij_log[c](pair->_name).getvalue();
                        Concavity -= (Im_lambda_var[arc->_id](pair->_name)-Im_lambda_sep[arc->_id](pair->_name).getvalue())*Im_Wij_log[c](pair->_name).getvalue();
                        DebugOff("R_Wij_log, " << R_Wij_log[c](pair->_name).getvalue()<< endl);
                        DebugOff("Im_Wij_log, " << Im_Wij_log[c](pair->_name).getvalue()<< endl);
                    }
                }

                for (auto arc: bag->get_in()) {
                    for (auto nn: arc->_intersection) {
                        Concavity += (lambda_var[arc->_id](nn->_name)-lambda_sep[arc->_id](nn->_name).getvalue())*Wii_log[c](nn->_name).getvalue();
                        DebugOff("Wii_log, " << Wii_log[c](nn->_name).getvalue()<< endl);
                    }
                    for (auto pair: arc->_intersection_clique) {
                        Concavity += (R_lambda_var[arc->_id](pair->_name)- R_lambda_sep[arc->_id](pair->_name).getvalue())*R_Wij_log[c](pair->_name).getvalue();
                        Concavity += (Im_lambda_var[arc->_id](pair->_name)-Im_lambda_sep[arc->_id](pair->_name).getvalue())*Im_Wij_log[c](pair->_name).getvalue();
                        DebugOff("R_Wij_log, " << R_Wij_log[c](pair->_name).getvalue()<< endl);
                        DebugOff("Im_Wij_log, " << Im_Wij_log[c](pair->_name).getvalue()<< endl);
                    }
                }
                Master.add_constraint(Concavity <= 0);
            }

            if (dual > sum(gamma_in).eval()) {
                for (int c = 0; c < nb_cliques; c++)
                    gamma_in(c) = value_dual[c];

                for (auto a: cliquetree->arcs) {
                    int l = a->_id;
                    for (auto arc: a->_intersection) {
                        lambda_in[l](arc->_name)  = lambda_sep[l](arc->_name).getvalue();
                    }

                    for (auto arc: a->_intersection_clique) {
                        R_lambda_in[l](arc->_name)  = R_lambda_sep[l](arc->_name).getvalue();
                        Im_lambda_in[l](arc->_name)  = Im_lambda_sep[l](arc->_name).getvalue();
                    }
                }
            }
            solver solve_master(Master, ipopt);
            solve_Master.run();
            DebugOff("master problem value: " << Master._obj_val << endl);

            // update the out point.
            gamma_out = (*(var<Real>*) Master.get_var("gamma_C"));
            for (auto a: cliquetree->arcs) {
                lambda_out[a->_id] = (*(var<Real>*) Master.get_var("lambda_arc_"+ to_string(a->_id)));
                R_lambda_out[a->_id] = (*(var<Real>*) Master.get_var("R_lambda_arc_" + to_string(a->_id)));
                Im_lambda_out[a->_id] = (*(var<Real>*) Master.get_var("Im_lambda_arc_"+ to_string(a->_id)));
            }
        }
        else {
            for (int c = 0; c < nb_cliques; c++)
                gamma_in(c) = value_dual[c];
            for (auto a: cliquetree->arcs) {
                int l = a->_id;
                for (auto arc: a->_intersection) {
                    lambda_in[l](arc->_name)  = lambda_sep[l](arc->_name).getvalue();
                }
                for (auto arc: a->_intersection_clique) {
                    R_lambda_in[l](arc->_name)  = R_lambda_sep[l](arc->_name).getvalue();
                    Im_lambda_in[l](arc->_name)  = Im_lambda_sep[l](arc->_name).getvalue();
                }
            }
        }
        LDlog[itcount] = dual;
        LBlog[itcount] = std::max(dual, LBlog[itcount-1]);
        if (itcount > 1)
            UBlog[itcount] = std::min(Master._obj_val, UBlog[itcount - 1]);
        else
            UBlog[itcount] = Master._obj_val;
    }
    double wall1 = get_wall_time();
    double cpu1 = get_cpu_time();
    cout << "CPU time: " << cpu1 -cpu0 << endl;
    cout << "Wall time: " << wall1 -wall0 << endl;
    for(int itcount = 1; itcount < iter_limit  ; itcount++) {
        cout<< setw(15) << left <<itcount << setw(15) << LDlog[itcount]<<setw(15)<< LBlog[itcount] << setw(15) << UBlog[itcount] << endl;
    }
    return 0;
}

int ADMM(PowerNet& grid) {
    cout << "///////////////////// ADMM algorithm //////////////////////////" << endl;
    Net* grid_augment = grid.clone_undirected();
    DebugOn("number of edges of the augmented graph, " << grid_augment->arcs.size() << endl);
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

// NOTE that grid_augment may NOT be chordal.
    grid_augment->get_tree_decomp_bags();
    auto cliquetree = grid_augment->get_clique_tree_prim(); // Note: it builds on grid_augment.

    const unsigned nb_cliques = grid_augment->_bags.size();

    vector<net_param> bag_param;
    vector<vector<Bus*>> bag_bus;
    vector<vector<Bus*>> bag_bus_disjoint;
    vector<vector<Gen*>> bag_gens;
    vector<vector<Gen*>> bag_gens_disjoint;
    vector<vector<Line*>> bag_arcs; // bag_arcs contains the arcs of the power grid while variables associated with W are defined on chordal graph.
    vector<vector<Line*>> bag_arcs_disjoint;
    vector<vector<gravity::index_pair*>> bag_bus_pairs; // bus_pairs in each bag.
    vector<vector<gravity::index_pair*> > bag_bus_pairs_disjoint; // bus_pairs in each bag.

    map<string, unsigned> indexij;
    map<string, vector<unsigned>> indexij_occurance;
    map<Bus*, unsigned> indexii;
    map<Bus*, vector<unsigned>> indexii_occurance;

    for (int c = 0; c < nb_cliques; c++) {
        net_param np;
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
        np.b_ff=grid.b_ff;
        np.b_ft=grid.b_ft;
        np.b_tf=grid.b_tf;
        np.b_tt=grid.b_tt;
        np.g_ff=grid.g_ff;
        np.g_ft=grid.g_ft;
        np.g_tf=grid.g_tf;
        np.g_tt=grid.g_tt;
        np.c0=grid.c0;
        np.c1=grid.c1;
        np.c2=grid.c2;
        np.tan_th_min=grid.tan_th_min;
        np.tan_th_max=grid.tan_th_max;
        np.S_max=grid.S_max;
        bag_param.push_back(np);
    }


    for (int c = 0; c < nb_cliques; c++) {
        vector<Bus*> bag_B_disjoint;
        vector<Gen*> bag_G_disjoint;

        vector<Line*> bag_A;
        vector<Line*> bag_A_disjoint;


        vector<gravity::index_pair*> bag_BP;
        vector<gravity::index_pair*> bag_BP_disjoint;


        sort(bag_bus[c].begin(), bag_bus[c].end(),node_id_compare);
        for (int i = 0; i < bag_bus[c].size(); i++) {
            Bus* b = bag_bus[c].at(i);
            vector<Node*> N = b->get_neighbours();
            sort(N.begin(), N.end(),node_id_compare);

            if (indexii.find(b)==indexii.end()) {
                bool inclusion = std::includes(bag_bus[c].begin(), bag_bus[c].end(), N.begin(),N.end(),node_id_compare);
                if (inclusion) {
                    indexii.insert(make_pair<>(b,c));
                    //indexii_occurance.insert(make_pair<>(b,a));
                    bag_B_disjoint.push_back(b);
                    if (b->_has_gen) {
                        bag_G_disjoint.insert(bag_G_disjoint.end(), b->_gen.begin(), b->_gen.end());
                    }
                }
            }
            vector<unsigned> a;
            a.push_back(c);
            auto pp =indexii_occurance.insert(make_pair<>(b,a));
            if (!pp.second) {
                indexii_occurance[b].push_back(c);
            }
            for (int j = i+1; j < bag_bus[c].size(); j++) {
                std::vector<Arc*> vec_arcs1 = grid.get_arcs(b, bag_bus[c].at(j));
                std::vector<Arc*> vec_arcs2 = grid.get_arcs(bag_bus[c].at(j), b);
                if (vec_arcs1.size() +vec_arcs2.size() > 0) {
                    if (vec_arcs1.size() > 0) {
                        bag_BP.push_back(new index_pair(index_(b->_name), index_(bag_bus[c].at(j)->_name), true));
                    }
                    else {
                        bag_BP.push_back(new index_pair(index_(bag_bus[c].at(j)->_name),index_(b->_name), true));
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
                            vector<unsigned> a;
                            a.push_back(c);
                            indexij_occurance.insert(make_pair<>(key, a));
                            bag_BP_disjoint.push_back(new index_pair(index_(b->_name), index_(bag_bus[c].at(j)->_name), true));
                        }
                        else {
                            indexij.insert(make_pair<>(key_inv, c));
                            vector<unsigned> a;
                            a.push_back(c);
                            indexij_occurance.insert(make_pair<>(key_inv, a));
                            bag_BP_disjoint.push_back(new index_pair(index_(bag_bus[c].at(j)->_name),index_(b->_name), true));
                        }

                        for (auto a: vec_arcs1) {
                            bag_A_disjoint.push_back((Line*)a);
                        }

                        for (auto a: vec_arcs2) {
                            bag_A_disjoint.push_back((Line*)a);
                        }
                    }
                    else {
                        if (vec_arcs1.size() > 0) {
                            indexij_occurance[key].push_back(c);
                        }
                        else {
                            indexij_occurance[key_inv].push_back(c);
                        }
                    }
                }
            }
        }

        sort(bag_BP.begin(), bag_BP.end(), bus_pair_compare);
        sort(bag_BP_disjoint.begin(), bag_BP_disjoint.end(), bus_pair_compare);
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
        nb_bus_pairs += bag_bus_pairs_disjoint[c].size();
    }

    DebugOn("the number of total arcs: " <<  nb_arcs << endl);
    DebugOn("the number of total bus pairs: " <<  nb_bus_pairs << endl);
    DebugOn("the number of total buses: " <<  nb_buses << endl);
    DebugOn("the number of total gens: " <<  nb_gens << endl);

    // intersection_cliques
    for (auto a: cliquetree->arcs) {
        std::vector<gravity::index_pair*> v3;
        std::set_intersection(bag_bus_pairs[a->_src->_id].begin(), bag_bus_pairs[a->_src->_id].end(),
                              bag_bus_pairs[a->_dest->_id].begin(), bag_bus_pairs[a->_dest->_id].end(),
                              back_inserter(v3), bus_pair_compare);
        if (v3.size() > 0) {
            a->_intersection_clique = v3;
        }
    }

    /** build model */
    Model CLT("Clique tree based Model");

    /** Variables */
    // local vars
    vector<var<Real>> R_Wij;
    vector<var<Real>> Im_Wij;
    vector<var<Real>> Wii;
    vector<var<Real>> Pg;
    vector<var<Real>> Qg;
    for (int c = 0; c < nb_cliques; c++) {
        var<Real>  bag_Wii("Wii_"+ to_string(c), grid.w_min.in(bag_bus[c]), grid.w_max.in(bag_bus[c]));
        CLT.add_var(bag_Wii^(bag_bus[c].size()));
        bag_Wii.initialize_all(1.001);
        Wii.push_back(bag_Wii);

        if (bag_bus_pairs[c].size() > 0) {
            var<Real>  bag_R_Wij("R_Wij_" + to_string(c), grid.wr_min.in(bag_bus_pairs[c]), grid.wr_max.in(bag_bus_pairs[c]));
            var<Real>  bag_Im_Wij("Im_Wij_" + to_string(c), grid.wi_min.in(bag_bus_pairs[c]), grid.wi_max.in(bag_bus_pairs[c]));
            CLT.add_var(bag_R_Wij^(bag_bus_pairs[c].size()));
            CLT.add_var(bag_Im_Wij^(bag_bus_pairs[c].size()));
            bag_R_Wij.initialize_all(1.0);
            R_Wij.push_back(bag_R_Wij);
            Im_Wij.push_back(bag_Im_Wij);
        }
        else {
            var<Real> empty1("R_Wij_"+to_string(c));
            var<Real> empty2("Im_Wij_"+to_string(c));
            empty1.set_size(0);
            empty2.set_size(0);
            R_Wij.push_back(empty1);
            Im_Wij.push_back(empty2);
        }

        if (bag_gens_disjoint[c].size() > 0) {
            var<Real>  bag_Pg("Pg_" + to_string(c), grid.pg_min.in(bag_gens_disjoint[c]), grid.pg_max.in(bag_gens_disjoint[c]));
            var<Real>  bag_Qg ("Qg_" + to_string(c), grid.qg_min.in(bag_gens_disjoint[c]), grid.qg_max.in(bag_gens_disjoint[c]));
            CLT.add_var(bag_Pg^(bag_gens_disjoint[c].size()));
            CLT.add_var(bag_Qg^(bag_gens_disjoint[c].size()));
            Pg.push_back(bag_Pg);
            Qg.push_back(bag_Qg);
        }
        else {
            var<Real> empty("empty");
            empty.set_size(0);
            Pg.push_back(empty);
            Qg.push_back(empty);
        }
    }

    // LOG OF SOLUTIONS
    // Log here means the previous primal and dual solution.
    vector<param<Real>> R_Wij_log;
    vector<param<Real>> Im_Wij_log;
    vector<param<Real>> Wii_log;

    param<Real> R_Zij_log("R_Zij_log");
    param<Real> Im_Zij_log("Im_Zij_log");
    param<Real> Zii_log("Zii_log");
    R_Zij_log^grid._bus_pairs._keys.size();
    Im_Zij_log^grid._bus_pairs._keys.size();
    Zii_log^grid.nodes.size();
    for (auto nn: grid.nodes) {
        Zii_log(nn->_name) = 0;
    }
    Im_Zij_log.initialize_all(0);
    R_Zij_log.initialize_all(0);

    for (int c = 0; c < nb_cliques; c++) {
        param<Real> Wii_C_log("Wii_C_log" + to_string(c));
        param<Real> Im_Wij_C_log("Im_Wij_C_log" + to_string(c));
        param<Real> R_Wij_C_log("R_Wij_C_log" + to_string(c));
        Wii_C_log^(bag_bus[c].size());
        R_Wij_C_log^(bag_bus_pairs[c].size());
        Im_Wij_C_log^(bag_bus_pairs[c].size());

        Wii_log.push_back(Wii_C_log);
        R_Wij_log.push_back(R_Wij_C_log);
        Im_Wij_log.push_back(Im_Wij_C_log);
    }

///////////////// DEFINE LAGRANGE MULTIPLIERS  ////////////////////////////////
    unsigned iter_limit;
    cout << "Enter the limit of the number of iterations: ";
    cin >> iter_limit;
    cout << endl;
    double epsilon = 0.0;

    double LRlog[iter_limit];
    double P_res[iter_limit];// primal residual
    double D_res[iter_limit];// dual residual

    vector<param<Real>> R_u;
    vector<param<Real>> Im_u;
    vector<param<Real>> u;

    for (int c = 0; c < nb_cliques; c++) {
        param<Real> u_bag("u_bag" + to_string(c));
        param<Real> Im_u_bag("Im_u_bag" + to_string(c));
        param<Real> R_u_bag("R_u_bag" + to_string(c));

        u_bag^(bag_bus[c].size());
        Im_u_bag^(bag_bus_pairs[c].size());
        R_u_bag^(bag_bus_pairs[c].size());

        u_bag.initialize_all(0);
        Im_u_bag.initialize_all(0);
        R_u_bag.initialize_all(0);

        u.push_back(u_bag);
        Im_u.push_back(Im_u_bag);
        R_u.push_back(R_u_bag);
    }

    //PENALTY PARAMETER
    double rho = 10;
    double mu = 10; // squared version
    double f_incr = 2;
    double f_decr = 2;
    double wall0 = get_wall_time();
    double cpu0  = get_cpu_time();
    //STRUCTURE
    cout << "<<<<<<<<<<< The standard ADMM algorithm starts >>>>>>>>>"<< endl;
    //for(int itcount = 0; itcount < iter_limit  ; itcount++) {
    int itcount =  0;
    P_res[0] = 100;
    D_res[0]=100;
    while ((P_res[itcount] > epsilon || D_res[itcount] < epsilon) && itcount < iter_limit){
        itcount += 1;
        std::vector<std::thread> threads;
        double val[nb_cliques];
        for (int c = 0; c < nb_cliques; c++) {
//        sub1(bag_param[c], c, rho, Pg[c], Qg[c], Wii[c], R_Wij[c], Im_Wij[c], bag_bus_pairs[c],
//                bag_bus[c],bag_bus_disjoint[c], bag_arcs[c], bag_gens_disjoint[c], bag_bus_pairs_disjoint[c],
//                   u, Im_u, R_u,  Wii_log[c], R_Wij_log[c], Im_Wij_log[c], Zii_log,  Im_Zij_log, R_Zij_log, val[c]);
            threads.push_back(std::thread(sub1, bag_param[c],  c, rho, Pg[c], Qg[c],Wii[c], R_Wij[c], Im_Wij[c], bag_bus_pairs[c],
                                          bag_bus[c],bag_bus_disjoint[c], bag_arcs_disjoint[c], bag_gens_disjoint[c], bag_bus_pairs_disjoint[c],
                                          u, Im_u, R_u, std::ref(Wii_log[c]), std::ref(R_Wij_log[c]), std::ref(Im_Wij_log[c]), Zii_log,
                                          Im_Zij_log, R_Zij_log, std::ref(val[c])));
        }

// join the threads with the main thread
        for(auto &t: threads) {
            t.join();
        }

        for (int c = 0; c < nb_cliques; c++) {
            LRlog[itcount] += val[c] ;
        }

//////////////////////////////// update Z//////////////////////////////////////
        D_res[itcount] = 0;
        for (index_pair* pair:grid._bus_pairs._keys) {
            double R=0;
            double Im=0;
            int l = indexij_occurance[pair->_name].size();
            for (int c: indexij_occurance[pair->_name]) {
                R += R_Wij_log[c](pair->_name).getvalue(); // + R_u[c](pair->_name).getvalue();
                Im += Im_Wij_log[c](pair->_name).getvalue(); //+ Im_u[c](pair->_name).getvalue();
            }
            D_res[itcount] += std::pow(R_Zij_log(pair->_name).getvalue() - R/l, 2);
            D_res[itcount] += std::pow(Im_Zij_log(pair->_name).getvalue() - Im/l, 2);
            R_Zij_log(pair->_name) = R/l;
            Im_Zij_log(pair->_name) = Im/l;
        }

        for (auto nn: grid.nodes) {
            double z=0;
            int  l = indexii_occurance[(Bus*)nn].size();
            for (int c: indexii_occurance[(Bus*)nn]) {
                z += Wii_log[c](nn->_name).getvalue(); //+ u[c](nn->_name).getvalue();
            }
            D_res[itcount] += std::pow(Zii_log(nn->_name).getvalue() - z/l, 2);
            Zii_log(nn->_name) = z/l;
        }
        D_res[itcount] = rho*std::sqrt(D_res[itcount]);


///////////////// Convergence information //////////////
        P_res[itcount] = 0;
        for (int c= 0; c < nb_cliques; c++) {
            for (auto nn: bag_bus[c]) {
                P_res[itcount] += std::pow(Wii_log[c](nn->_name).getvalue()-Zii_log(nn->_name).getvalue(), 2);
            }

            for (auto pair:bag_bus_pairs[c]) {
                P_res[itcount] +=std::pow(Im_Wij_log[c](pair->_name).getvalue()-Im_Zij_log(pair->_name).getvalue(), 2);
                P_res[itcount] +=std::pow(R_Wij_log[c](pair->_name).getvalue()-R_Zij_log(pair->_name).getvalue(), 2);
            }
        }
        P_res[itcount] = std::sqrt(P_res[itcount]);

//////////////// update dual multipliers ////////////////////
        for (int c= 0; c < nb_cliques; c++) {
            for (auto nn: bag_bus[c])
                u[c](nn->_name) = u[c](nn->_name).getvalue() + Wii_log[c](nn->_name).getvalue()-Zii_log(nn->_name).getvalue();

            for (auto pair:bag_bus_pairs[c]) {
                Im_u[c](pair->_name) = Im_u[c](pair->_name).getvalue() + Im_Wij_log[c](pair->_name).getvalue()-Im_Zij_log(pair->_name).getvalue();
                R_u[c](pair->_name) = R_u[c](pair->_name).getvalue()+ R_Wij_log[c](pair->_name).getvalue()-R_Zij_log(pair->_name).getvalue();
            }
        }

///////////// update the penalty parameters ///////////////
        if (P_res[itcount] > mu*D_res[itcount]) {
            rho = f_incr*rho;
        }
        else if (D_res[itcount] > mu*P_res[itcount]) {
            rho = rho/f_decr;
        }
        else {
            continue;
        }
    }

    cout<< setw(15) << left <<"ITERATION" << setw(15) << "LR" << setw(15) << "P_res" << setw(15) << "D_res" << endl;
    double wall1 = get_wall_time();
    double cpu1 = get_cpu_time();
    cout << "CPU time: " << cpu1 -cpu0 << endl;
    cout << "Wall time: " << wall1 -wall0 << endl;


    for(int itcount = 0; itcount < iter_limit  ; itcount++) {
        cout<< setw(15) << left <<itcount << setw(15) << LRlog[itcount]<<setw(15) << P_res[itcount]<<setw(15)<< D_res[itcount] << endl;
    }
    return 0;
}

int main (int argc, const char * argv[])
{
    // Decompose
    const char* fname;
    if (argc >= 2) {
        fname = argv[1];
    }
    else {
        //fname = "../../data_sets/Power/nesta_case5_pjm.m";
        //fname = "../../data_sets/Power/nesta_case30_ieee.m";
        //fname = "../../data_sets/Power/nesta_case6_c.m";
        //fname = "../../data_sets/Power/nesta_case5_pjm.m";
        //fname = "../../data_sets/Power/nesta_case3_lmbd.m";
        //fname = "../../data_sets/Power/nesta_case300_ieee.m";
        //fname = "../../data_sets/Power/nesta_case118_ieee.m";
        //fname = "../../data_sets/Power/nesta_case57_ieee.m";
        fname = "../../data_sets/Power/nesta_case14_ieee.m";
        //fname = "../../data_sets/Power/nesta_case57_ieee.m";
    }
    PowerNet grid;
    grid.readgrid(fname);
    cout << "////////////////////////////////////////" << endl;
    inout(grid);
    //ADMM(grid);
    return 0;
}