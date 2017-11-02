// File Name : functions.h
// Time : Sat 28 Oct 14:02:18 2017
// Created By : Guanglei Wang
//
#ifndef functions_h
#define functions_h
#include <stdio.h>
#include <iostream>
#include <string>
#include <stdio.h>
#include <cstring>
#include <fstream>
#include "../PowerNet.h"
#include <gravity/solver.h>

#ifdef USE_BOOST
#include <boost/graph/adjacency_list.hpp>
#include <iterator>
#endif


// CLIQUE BASED REFORMULATION
using namespace std;
using namespace gravity;

void scopf_W(PowerNet& grid, bool include_G)
{
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
    if (include_G) {
        SOCP.add_var(Pg^(nb_gen));
        SOCP.add_var(Qg^(nb_gen));
    }
    else {
        for (auto g:grid.gens) {
            if (g->_active) {
                Constraint Production_P_UB("Production_P_UB" + g->_name);
                Constraint Production_P_LB("Production_P_LB" + g->_name);
                Constraint Production_Q_UB("Production_Q_UB" + g->_name);
                Constraint Production_Q_LB("Production_Q_LB" + g->_name);
                auto b = g->_bus;
                Production_P_UB += b->gs()*Wii(b->_name) + b->pl() - grid.pg_max(g->_name).getvalue();
                Production_P_LB += b->gs()*Wii(b->_name) + b->pl() - grid.pg_min(g->_name).getvalue();
                Production_Q_UB += -b->bs()*Wii(b->_name) + b->ql() - grid.qg_max(g->_name).getvalue();
                Production_Q_LB += -b->bs()*Wii(b->_name) + b->ql() - grid.qg_min(g->_name).getvalue();
                for (auto &a: b->get_out()) {
                    Production_P_UB += grid.g_ff(a->_name)*Wii(b->_name)
                                       + grid.g_ft(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name)
                                       + grid.b_ft(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);

                    Production_P_LB += grid.g_ff(a->_name)*Wii(b->_name)
                                       + grid.g_ft(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name)
                                       + grid.b_ft(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);

                    Production_Q_UB  += -1*grid.b_ff(a->_name)*Wii(b->_name)
                                        -grid.b_ft(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name)
                                        + grid.g_ft(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);

                    Production_Q_LB  += -1*grid.b_ff(a->_name)*Wii(b->_name)
                                        -grid.b_ft(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name)
                                        + grid.g_ft(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);
                }

                for (auto &a: b->get_in()) {
                    Production_P_UB += grid.g_tt(a->_name)*Wii(b->_name)
                                       + grid.g_tf(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name)
                                       - grid.b_tf(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);

                    Production_P_LB += grid.g_tt(a->_name)*Wii(b->_name)
                                       + grid.g_tf(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name)
                                       - grid.b_tf(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);

                    Production_Q_UB -= grid.b_tt(a->_name)*Wii(b->_name)
                                       + grid.b_tf(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name)
                                       + grid.g_tf(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);

                    Production_Q_LB -= grid.b_tt(a->_name)*Wii(b->_name)
                                       + grid.b_tf(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name)
                                       + grid.g_tf(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);
                }
                SOCP.add_constraint(Production_P_UB <= 0);
                SOCP.add_constraint(Production_P_LB >= 0);
                SOCP.add_constraint(Production_Q_UB <= 0);
                SOCP.add_constraint(Production_Q_LB >= 0);
            }
        }
    }

    /** Construct the objective function*/
    func_ obj;
    if (include_G) {
        for (auto g:grid.gens) {
            if (g->_active) {
                obj += grid.c1(g->_name)*Pg(g->_name) + grid.c2(g->_name)*Pg(g->_name)*Pg(g->_name) + grid.c0(g->_name);
            }
        }
    }
    else {
        // NEED TO CONSIDER SECOND ORDER TERMS.
        for (auto g:grid.gens) {
            if (g->_active) {
                auto bus = g->_bus;
                obj += grid.c1(g->_name).getvalue()*bus->pl() + grid.c0(g->_name).getvalue();
                //shunt
                obj  += grid.c1(g->_name)*bus->gs()*Wii(bus->_name);
                for (auto &a: bus->get_out()) {
                    obj  += grid.c1(g->_name)*grid.g_ff(a->_name)*Wii(a->_src->_name)
                            + grid.c1(g->_name)*grid.g_ft(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name)
                            + grid.c1(g->_name)*grid.b_ft(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);
                }

                for (auto &a: bus->get_in()) {
                    obj  += grid.c1(g->_name)*grid.g_tt(a->_name)*Wii(a->_dest->_name)
                            + grid.c1(g->_name)*grid.g_tf(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name)
                            - grid.c1(g->_name)*grid.b_tf(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);
                }
            }
        }
    }
    SOCP.set_objective(min(obj));

    /** Define constraints */
    /* SOCP constraints */
    Constraint SOC("SOC");
    SOC =  power(R_Wij.in(bus_pairs), 2) + power(Im_Wij.in(bus_pairs), 2) - Wii.from(bus_pairs)*Wii.to(bus_pairs) ;
    SOCP.add_constraint(SOC <= 0);

    //KCL
    if (include_G) {
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

            /* Power Conservation */
            //KCL_P += innerproduct(grid.g_ff.in(b->get_out()), Wii.from(b->get_out()))
            //         +innerproduct(grid.g_ft.in(b->get_out()), R_Wij.in_pairs(b->get_out()))
            //         + innerproduct(grid.b_ft.in(b->get_out()), Im_Wij.in_pairs(b->get_out()));

            //KCL_Q += innerproduct(grid.b_ff.in(b->get_out()), Wii.from(b->get_out()))
            //    -innerproduct(grid.b_ft.in(b->get_out()), R_Wij.in_pairs(b->get_out()))
            //    + innerproduct(grid.g_ft.in(b->get_out()), Im_Wij.in_pairs(b->get_out()));

            //KCL_P += innerproduct(grid.g_tt.in(b->get_in()),Wii.to(b->get_in()))
            //    + innerproduct(grid.g_tf.in(b->get_in()), R_Wij.in_pairs(b->get_in()))
            //    - innerproduct(grid.b_tf.in(b->get_in()), Im_Wij.in_pairs(b->get_in()));

            //KCL_Q -= innerproduct(grid.b_tt.in(b->get_in()),Wii.to(b->get_in()))
            //    + innerproduct(grid.b_tf.in(b->get_in()), R_Wij.in_pairs(b->get_in()))
            //    + innerproduct(grid.g_tf.in(b->get_in()), Im_Wij.in_pairs(b->get_in()));

            KCL_P += bus->pl()- sum(Pg.in(bus->_gen));;
            KCL_Q -= bus->ql()- sum(Qg.in(bus->_gen));

            /* Shunts */
            KCL_P +=  bus->gs()*(Wii(bus->_name));
            KCL_Q -=  bus->bs()*(Wii(bus->_name));

            SOCP.add_constraint(KCL_P = 0);
            SOCP.add_constraint(KCL_Q = 0);
        }
    }
    else {
        for (auto b: grid.nodes) {
            Bus* bus = (Bus*) b;
            if (!bus->_has_gen) {
                Constraint KCL_P("KCL_P"+bus->_name);
                Constraint KCL_Q("KCL_Q"+bus->_name);
                for (auto &a: b->get_out()) {
                    KCL_P  += grid.g_ff(a->_name)*Wii(a->_src->_name)
                              +grid.g_ft(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name)
                              +grid.b_ft(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);
                    KCL_Q  += -1*grid.b_ff(a->_name)*Wii(a->_src->_name)
                              - grid.b_ft(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name)
                              +grid.g_ft(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);
                }

                for (auto &a: b->get_in()) {
                    KCL_P  += grid.g_tt(a->_name)*Wii(a->_dest->_name)
                              +grid.g_tf(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name)-grid.b_tf(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);

                    KCL_Q  -= grid.b_tt(a->_name)*Wii(a->_dest->_name)
                              + grid.b_tf(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name)
                              + grid.g_tf(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);
                }

                KCL_P += bus->pl();
                KCL_Q += bus->ql();

                KCL_P +=  bus->gs()*(Wii(bus->_name));
                KCL_Q -=  bus->bs()*(Wii(bus->_name));

                SOCP.add_constraint(KCL_P = 0);
                SOCP.add_constraint(KCL_Q = 0);
            }
        }
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

// clique based reformulation
void OPF_Clique_W(PowerNet& grid)
{
    /** Clique tree decomposition **/
    Net* chordal = grid.get_chordal_extension();
    auto cliquetree = grid.get_clique_tree();
    const unsigned nb_cliques = grid._bags.size();

    vector<vector<Bus*>> bag_bus; 
    vector<vector<Bus*>> bag_bus_disjoint;
    vector<vector<Gen*>> bag_gens;
    vector<vector<Gen*>> bag_gens_disjoint;
    vector<vector<Line*>> bag_arcs; //bag_arcs contains the arcs of the power grid while variables associated with W are defined on chordal graph.
    vector<vector<Line*>> bag_arcs_disjoint;
    vector<gravity::node_pairs*> bag_bus_pairs; // bus_pairs in each bag.
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
        for (int i = 0; i < grid._bags[c].size(); i++) {
            Bus* b = (Bus*) grid.get_node(grid._bags[c].at(i)->_name);
            if (b !=nullptr) {
                DebugOn(b->_name << " ");
                bag_B.push_back(b);
                auto pp = indexii.insert(make_pair<>(b, c));
                if (pp.second) {
                    bag_B_disjoint.push_back(b);
                }

                if (b->_has_gen) {
                    bag_G.insert(bag_G.end(), b->_gen.begin(), b->_gen.end());
                }

                if (pp.second && b->_has_gen) {
                    bag_G_disjoint.insert(bag_G_disjoint.end(), b->_gen.begin(), b->_gen.end());
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
                    }
                }
            }
        }
        
        DebugOn("] " << endl);
        bag_bus.push_back(bag_B);
        bag_gens.push_back(bag_G);
        bag_arcs.push_back(bag_A);
        bag_bus_pairs.push_back(bag_BP);
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
    for (int c = 0; c < nb_cliques; c++) {
        var<Real>  bag_R_Wij("R_Wij_" + to_string(c), grid.wr_min.in(bag_bus_pairs[c]->_keys), grid.wr_max.in(bag_bus_pairs[c]->_keys)); // real part of Wij
        var<Real>  bag_Im_Wij("Im_Wij_" + to_string(c), grid.wi_min.in(bag_bus_pairs[c]->_keys), grid.wi_max.in(bag_bus_pairs[c]->_keys)); // imaginary part of Wij.
        var<Real>  bag_Wii("Wii_" + to_string(c), grid.w_min.in(bag_bus[c]), grid.w_max.in(bag_bus[c]));
        CLT.add_var(bag_Wii^(bag_bus[c].size()));
        CLT.add_var(bag_R_Wij^(bag_bus[c].size()*(bag_bus[c].size()-1)/2));// (Maybe bag_bus[c]*bag_bus[c] -1)/2
        CLT.add_var(bag_Im_Wij^(bag_bus[c].size()*(bag_bus[c].size()-1)/2));

        bag_R_Wij.initialize_all(1.0);
        bag_Wii.initialize_all(1.001);
        R_Wij.push_back(bag_R_Wij);
        Im_Wij.push_back(bag_Im_Wij);
        Wii.push_back(bag_Wii);
    }

    /* Construct the objective function with generations bound constraints */
    func_ obj;
    obj += sum(grid.c0.in(grid.gens));
    for (int c = 0; c < nb_cliques; c++) {
        for (auto g:grid.gens) {
            if (g->_active) {
                auto b = g->_bus;
                if (std::find(bag_bus_disjoint[c].begin(), bag_bus_disjoint[c].end(), b) != bag_bus_disjoint[c].end()) {
                    obj += grid.c1(g->_name).getvalue()*b->pl() + grid.c0(g->_name).getvalue();
                    obj += (grid.c1(g->_name)*b->gs())*Wii[c](b->_name);
                }

                for (auto &a: b->get_out()) {
                    if (std::find(bag_arcs_disjoint[c].begin(), bag_arcs_disjoint[c].end(), a) != bag_arcs_disjoint[c].end()) {
                        obj  += grid.c1(g->_name)*grid.g_ff(a->_name)*Wii[c](a->_src->_name)
                                + grid.c1(g->_name)*grid.g_ft(a->_name)*R_Wij[c](a->_src->_name+","+a->_dest->_name)
                                + grid.c1(g->_name)*grid.b_ft(a->_name)*Im_Wij[c](a->_src->_name+","+a->_dest->_name);
                    }
                }
                for (auto &a: b->get_in()) {
                    if (std::find(bag_arcs_disjoint[c].begin(), bag_arcs_disjoint[c].end(),a) != bag_arcs_disjoint[c].end()) {
                        obj  += grid.c1(g->_name)*grid.g_tt(a->_name)*Wii[c](a->_dest->_name)
                                +grid.c1(g->_name)*grid.g_tf(a->_name)*R_Wij[c](a->_src->_name+","+a->_dest->_name)
                                - grid.c1(g->_name)*grid.b_tf(a->_name)*Im_Wij[c](a->_src->_name+","+a->_dest->_name);
                    }
                }
            }
         }
    }

    CLT.set_objective(min(obj));

    /** Define constraints */
    /* CLT constraints */
    for (int c = 0; c < nb_cliques; c++) {
        Constraint SOC("SOC_" + to_string(c));
        SOC =  power(R_Wij[c].in(bag_bus_pairs[c]->_keys), 2)
               + power(Im_Wij[c].in(bag_bus_pairs[c]->_keys), 2)
               - Wii[c].from(bag_bus_pairs[c]->_keys)*Wii[c].to(bag_bus_pairs[c]->_keys) ;
        CLT.add_constraint(SOC <= 0);
    }


    for (int c = 0; c < nb_cliques; c++) {
        //* Phase Angle Bounds constraints */
        Constraint PAD_UB("PAD_UB" + to_string(c));
        PAD_UB = Im_Wij[c].in(bag_bus_pairs[c]->_keys);
        PAD_UB -= (grid.tan_th_max).in(bag_bus_pairs[c]->_keys)*R_Wij[c].in(bag_bus_pairs[c]->_keys);
        CLT.add_constraint(PAD_UB <= 0);

        Constraint PAD_LB("PAD_LB" + to_string(c));
        PAD_LB = Im_Wij[c].in(bag_bus_pairs[c]->_keys);
        PAD_LB -= (grid.tan_th_min).in(bag_bus_pairs[c]->_keys)*R_Wij[c].in(bag_bus_pairs[c]->_keys);
        CLT.add_constraint(PAD_LB >= 0);

        /* Thermal Limit Constraints */
        Constraint Thermal_Limit_from("Thermal_Limit_from" + to_string(c));
        Thermal_Limit_from += power(grid.g_ff.in(bag_arcs[c])*Wii[c].from(bag_arcs[c])+ grid.g_ft.in(bag_arcs[c])*R_Wij[c].in_pairs(bag_arcs[c])
                                    + grid.b_ft.in(bag_arcs[c])*Im_Wij[c].in_pairs(bag_arcs[c]), 2)
                              + power(grid.g_ft.in(bag_arcs[c])*Im_Wij[c].in_pairs(bag_arcs[c])-grid.b_ff.in(bag_arcs[c])*Wii[c].from(bag_arcs[c])
                                      - grid.b_ft.in(bag_arcs[c])*R_Wij[c].in_pairs(bag_arcs[c]), 2);
        Thermal_Limit_from -= power(grid.S_max.in(bag_arcs[c]),2);
        CLT.add_constraint(Thermal_Limit_from <= 0);

        Constraint Thermal_Limit_to("Thermal_Limit_to" + to_string(c));
        Thermal_Limit_to += power(grid.g_tt.in(bag_arcs[c])*Wii[c].from(bag_arcs[c]) + grid.g_tf.in(bag_arcs[c])*R_Wij[c].in_pairs(bag_arcs[c])
                                  + grid.b_tf.in(bag_arcs[c])*Im_Wij[c].in_pairs(bag_arcs[c]), 2)
                            + power(grid.b_tt.in(bag_arcs[c])*Wii[c].from(bag_arcs[c]) + grid.b_tf.in(bag_arcs[c])*R_Wij[c].in_pairs(bag_arcs[c])
                                    + grid.g_tf.in(bag_arcs[c])*Im_Wij[c].in_pairs(bag_arcs[c]), 2);

        Thermal_Limit_to -= power(grid.S_max.in(bag_arcs[c]), 2);
        CLT.add_constraint(Thermal_Limit_to <= 0);
    }
    
    // COUPLING CONSTRAINTS
    for (auto b: grid.nodes) {
        Bus* bus = (Bus*) b;
        if (!bus->_has_gen) {
            Constraint KCL_P("KCL_P"+bus->_name);
            Constraint KCL_Q("KCL_Q"+bus->_name);
            KCL_P += bus->pl();
            KCL_Q += bus->ql();
            for (int c = 0; c < nb_cliques; c++) {
                if (std::find(bag_bus_disjoint[c].begin(), bag_bus_disjoint[c].end(), bus) != bag_bus_disjoint[c].end()) {
                    KCL_P +=  bus->gs()*(Wii[c](bus->_name));
                    KCL_Q -=  bus->bs()*(Wii[c](bus->_name));
                }
                for (auto &a: b->get_out()) {
                    if (std::find(bag_arcs_disjoint[c].begin(), bag_arcs_disjoint[c].end(), a) != bag_arcs_disjoint[c].end()) {
                        KCL_P  += grid.g_ff(a->_name)*Wii[c](bus->_name)
                                  +grid.g_ft(a->_name)*R_Wij[c](a->_src->_name+","+a->_dest->_name)
                                  +grid.b_ft(a->_name)*Im_Wij[c](a->_src->_name+","+a->_dest->_name);

                        KCL_Q  += -1*grid.b_ff(a->_name)*Wii[c](bus->_name)
                                  - grid.b_ft(a->_name)*R_Wij[c](a->_src->_name+","+a->_dest->_name)
                                  +grid.g_ft(a->_name)*Im_Wij[c](a->_src->_name+","+a->_dest->_name);
                    }
                }

                for (auto &a: b->get_in()) {
                    if (std::find(bag_arcs_disjoint[c].begin(), bag_arcs_disjoint[c].end(), a) != bag_arcs_disjoint[c].end()) {
                        KCL_P  += grid.g_tt(a->_name)*Wii[c](bus->_name)
                                  +grid.g_tf(a->_name)*R_Wij[c](a->_src->_name+","+a->_dest->_name)-grid.b_tf(a->_name)*Im_Wij[c](a->_src->_name+","+a->_dest->_name);

                        KCL_Q  -= grid.b_tt(a->_name)*Wii[c](bus->_name)
                                  + grid.b_tf(a->_name)*R_Wij[c](a->_src->_name+","+a->_dest->_name)
                                  + grid.g_tf(a->_name)*Im_Wij[c](a->_src->_name+","+a->_dest->_name);
                    }
                }
            }
            CLT.add_constraint(KCL_P = 0);
            CLT.add_constraint(KCL_Q = 0);
        }
    }

    for (auto g:grid.gens) {
        if (g->_active) {
            Constraint Production_P_UB("Production_P_UB" + g->_name);
            Constraint Production_P_LB("Production_P_LB" + g->_name);
            Constraint Production_Q_UB("Production_Q_UB" + g->_name);
            Constraint Production_Q_LB("Production_Q_LB" + g->_name);
            auto b = g->_bus;
            Production_P_UB += b->pl() - grid.pg_max(g->_name).getvalue();
            Production_P_LB += b->pl() - grid.pg_min(g->_name).getvalue();
            Production_Q_UB += b->ql() - grid.qg_max(g->_name).getvalue();
            Production_Q_LB += b->ql() - grid.qg_min(g->_name).getvalue();
            for (int c = 0; c < nb_cliques; c++) {
                if (std::find(bag_bus_disjoint[c].begin(), bag_bus_disjoint[c].end(), b) != bag_bus_disjoint[c].end()) {
                    Production_P_UB += b->gs()*Wii[c](b->_name);
                    Production_P_LB += b->gs()*Wii[c](b->_name);
                    Production_Q_UB += -b->bs()*Wii[c](b->_name);
                    Production_Q_LB += -b->bs()*Wii[c](b->_name);
                }
                for (auto &a: b->get_out()) {
                    if (std::find(bag_arcs_disjoint[c].begin(), bag_arcs_disjoint[c].end(), a) != bag_arcs_disjoint[c].end()) {
                        Production_P_UB  += grid.g_ff(a->_name)*Wii[c](a->_src->_name)
                                           + grid.g_ft(a->_name)*R_Wij[c](a->_src->_name+","+a->_dest->_name)
                                           + grid.b_ft(a->_name)*Im_Wij[c](a->_src->_name+","+a->_dest->_name);

                        Production_P_LB  += grid.g_ff(a->_name)*Wii[c](a->_src->_name)
                                           + grid.g_ft(a->_name)*R_Wij[c](a->_src->_name+","+a->_dest->_name)
                                           + grid.b_ft(a->_name)*Im_Wij[c](a->_src->_name+","+a->_dest->_name);

                        Production_Q_UB  +=  -1*grid.b_ff(a->_name)*Wii[c](a->_src->_name)
                                             -grid.b_ft(a->_name)*R_Wij[c](a->_src->_name+","+a->_dest->_name)
                                             + grid.g_ft(a->_name)*Im_Wij[c](a->_src->_name+","+a->_dest->_name);

                        Production_Q_LB  += -1*grid.b_ff(a->_name)*Wii[c](a->_src->_name)
                                            -grid.b_ft(a->_name)*R_Wij[c](a->_src->_name+","+a->_dest->_name)
                                            + grid.g_ft(a->_name)*Im_Wij[c](a->_src->_name+","+a->_dest->_name);
                    }
                }
                for (auto &a: b->get_in()) {
                    if (std::find(bag_arcs_disjoint[c].begin(), bag_arcs_disjoint[c].end(),a) != bag_arcs_disjoint[c].end()) {
                        Production_P_UB = grid.g_tt(a->_name)*Wii[c](a->_dest->_name)
                                         + grid.g_tf(a->_name)*R_Wij[c](a->_src->_name+","+a->_dest->_name)
                                        -1*grid.b_tf(a->_name)*Im_Wij[c](a->_src->_name+","+a->_dest->_name);

                        Production_P_LB = grid.g_tt(a->_name)*Wii[c](a->_dest->_name)
                                         + grid.g_tf(a->_name)*R_Wij[c](a->_src->_name+","+a->_dest->_name)
                                        -1*grid.b_tf(a->_name)*Im_Wij[c](a->_src->_name+","+a->_dest->_name);

                        Production_Q_UB -= grid.b_tt(a->_name)*Wii[c](a->_dest->_name)
                                           + grid.b_tf(a->_name)*R_Wij[c](a->_src->_name+","+a->_dest->_name)
                                           + grid.g_tf(a->_name)*Im_Wij[c](a->_src->_name+","+a->_dest->_name);

                        Production_Q_LB -= grid.b_tt(a->_name)*Wii[c](a->_dest->_name)
                                           + grid.b_tf(a->_name)*R_Wij[c](a->_src->_name+","+a->_dest->_name)
                                           + grid.g_tf(a->_name)*Im_Wij[c](a->_src->_name+","+a->_dest->_name);
                    }
                }
            }
            CLT.add_constraint(Production_P_UB<= 0);
            CLT.add_constraint(Production_P_LB>= 0);
            CLT.add_constraint(Production_Q_UB <= 0);
            CLT.add_constraint(Production_Q_LB >= 0);
        }
    }

    for (auto a: cliquetree->arcs){
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
    
   //
    
    for (int c = 0; c < nb_cliques; c++){
         auto W1 = (*(param<double>*)(CLT.get_var("R_Wij_"+to_string(c))));
        W1.print(true);
    }
    cout << "\n \n" << endl;
    for (int c = 0; c < nb_cliques; c++){
         auto W1 = (*(param<double>*)(CLT.get_var("Im_Wij_"+to_string(c))));
        W1.print(true);
    }

}
#endif
