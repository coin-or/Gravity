/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : functions.h

* Purpose :

* Creation Date : 27-10-2017

* Last Modified : Fri 27 Oct 15:06:01 2017

* Created By : Guanglei Wang
_._._._._._._._._._._._._._._._._._._._._.*/
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
#include <stdio.h>
#include <stdlib.h>

#ifdef USE_BOOST
#include <boost/graph/adjacency_list.hpp>
#include <iterator>
#endif


// CLIQUE BASED REFORMULATION
using namespace std;
using namespace gravity;

void scopf_W(PowerNet* grid, bool include_G)
{
    // Grid Parameters
    auto bus_pairs = grid->get_bus_pairs();
    auto nb_bus_pairs = bus_pairs.size();
    auto nb_gen = grid->get_nb_active_gens();
    auto nb_lines = grid->get_nb_active_arcs();
    auto nb_buses = grid->get_nb_active_nodes();

    /** Clique tree decomposition **/
    Net* chordal = grid->get_chordal_extension();
    grid->get_clique_tree();
    const unsigned nb_cliques = grid->_bags.size();

    vector<vector<Bus*>> bag_bus;
    vector<vector<Gen*>> bag_gens;
    vector<vector<Arc*>> bag_arcs;

    for (int c = 0; c < nb_cliques; c++) {
        vector<Bus*> bag_B;
        vector<Gen*> bag_G;
        vector<Arc*> bag_A;
        for (int i = 0; i < grid->_bags[c].size(); i++) {
            Bus* bus = (Bus*) grid->get_node(grid->_bags[c].at(i)->_name);
            if (bus !=nullptr) {
                bag_B.push_back(bus);
            }
            if (bus->_has_gen) {
                bag_G.insert(bag_G.end(), bus->_gen.begin(), bus->_gen.end());
            }
            for (int j = i+1; j < grid->_bags[c].size(); j++) {
                Arc* a = (Arc*)grid->get_arc(bus, grid->get_node(grid->_bags[c].at(j)->_name));
                if (a != nullptr)
                    bag_A.push_back(a);
            }
        }
        bag_bus.push_back(bag_B);
        bag_gens.push_back(bag_G);
        bag_arcs.push_back(bag_A);
    }
    bag_bus.resize(nb_cliques);
    bag_arcs.resize(nb_cliques);
    bag_gens.resize(nb_cliques);

    Net* cliquetree = grid->get_clique_tree();

    /** build model */
    Model SOCP("SOCP Model");
    /** Variables */
    // lIFTED VARIABLES.
    var<Real>  R_Wij("R_Wij", grid->wr_min.in(bus_pairs), grid->wr_max.in(bus_pairs)); // real part of Wij
    var<Real>  Im_Wij("Im_Wij", grid->wi_min.in(bus_pairs), grid->wi_max.in(bus_pairs)); // imaginary part of Wij.
    var<Real>  Wii("Wii", grid->w_min.in(grid->nodes), grid->w_max.in(grid->nodes));
    SOCP.add_var(Wii^nb_buses);
    SOCP.add_var(R_Wij^nb_bus_pairs);
    SOCP.add_var(Im_Wij^nb_bus_pairs);
    
    R_Wij.initialize_all(1.0);
    Wii.initialize_all(1.001);
    
    // power generation
    var<Real> Pg("Pg", grid->pg_min.in(grid->gens), grid->pg_max.in(grid->gens));
    var<Real> Qg ("Qg", grid->qg_min.in(grid->gens), grid->qg_max.in(grid->gens));
    if (include_G) {
        SOCP.add_var(Pg^(nb_gen));
        SOCP.add_var(Qg^(nb_gen));
    }
    else{
        for (auto g:grid->gens) {
            if (g->_active) {
                Constraint Production_P_UB("Production_P_UB" + g->_name);
                Constraint Production_P_LB("Production_P_LB" + g->_name);
                Constraint Production_Q_UB("Production_Q_UB" + g->_name);
                Constraint Production_Q_LB("Production_Q_LB" + g->_name);
                auto b = g->_bus;
                Production_P_UB += b->gs()*Wii(b->_name) + b->pl() - grid->pg_max(g->_name).getvalue();
                Production_P_LB += b->gs()*Wii(b->_name) + b->pl() - grid->pg_min(g->_name).getvalue();
                Production_Q_UB += b->bs()*Wii(b->_name) + b->ql() - grid->qg_max(g->_name).getvalue();
                Production_Q_LB += -b->bs()*Wii(b->_name) + b->ql() - grid->qg_min(g->_name).getvalue();
                for (auto &a: b->get_out()) {
                    Production_P_UB += grid->g_ff(a->_name)*Wii(a->_src->_name)
                                    + grid->g_ft(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name)
                                    + grid->b_ft(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);

                    Production_P_LB += grid->g_ff(a->_name)*Wii(a->_src->_name)
                                    + grid->g_ft(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name) 
                                    + grid->b_ft(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);

                    Production_Q_UB  += -1*grid->b_ff(a->_name)*Wii(a->_src->_name)
                                      -grid->b_ft(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name) 
                                      + grid->g_ft(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);

                    Production_Q_LB  += -1*grid->b_ff(a->_name)*Wii(a->_src->_name)
                                         -grid->b_ft(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name) 
                                         + grid->g_ft(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);
                }

                for (auto &a: b->get_in()) {
                    Production_P_UB += grid->g_tt(a->_name)*Wii(a->_dest->_name)
                                    + grid->g_tf(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name)
                                    - grid->b_tf(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);

                    Production_P_LB += grid->g_tt(a->_name)*Wii(a->_dest->_name)
                                    + grid->g_tf(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name)
                                    - grid->b_tf(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);

                    Production_Q_UB -= grid->b_tt(a->_name)*Wii(a->_dest->_name) 
                                     + grid->b_tf(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name) 
                                     + grid->g_tf(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);

                    Production_Q_LB -= grid->b_tt(a->_name)*Wii(a->_dest->_name) 
                                     + grid->b_tf(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name) 
                                     + grid->g_tf(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);
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
        for (auto g:grid->gens) {
            if (g->_active) {
                obj += grid->c1(g->_name)*Pg(g->_name) + grid->c2(g->_name)*Pg(g->_name)*Pg(g->_name) + grid->c0(g->_name);
            }
        }
    }
    else {
        // NEED TO CONSIDER SECOND ORDER TERMS. 
        for (auto g:grid->gens) {
            if (g->_active) {
                auto bus = g->_bus;
                obj += grid->c1(g->_name)*bus->pl() + grid->c0(g->_name);
                for (auto &a: bus->get_out()) {
                    obj  += grid->c1(g->_name)*grid->g_ff(a->_name)*Wii(a->_src->_name)
                            + grid->c1(g->_name)*grid->g_ft(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name) 
                            + grid->c1(g->_name)*grid->b_ft(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);
                }

                for (auto &a: bus->get_in()) {
                    obj  += grid->c1(g->_name)*grid->g_tt(a->_name)*Wii(a->_dest->_name)
                            + grid->c1(g->_name)*grid->g_tf(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name)
                            - grid->c1(g->_name)*grid->b_tf(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);
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
        for (auto b: grid->nodes) {
            Bus* bus = (Bus*) b;
            Constraint KCL_P("KCL_P"+bus->_name);
            Constraint KCL_Q("KCL_Q"+bus->_name);

            /* Power Conservation */
            for (auto &a: b->get_out()) {
                KCL_P  += grid->g_ff(a->_name)*Wii(a->_src->_name)
                          +grid->g_ft(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name)+grid->b_ft(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);

                KCL_Q  += -1*grid->b_ff(a->_name)*Wii(a->_src->_name) -grid->b_ft(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name) 
                    + grid->g_ft(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);
            }

            for (auto &a: b->get_in()) {
                KCL_P  += grid->g_tt(a->_name)*Wii(a->_dest->_name)
                          +grid->g_tf(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name)-grid->b_tf(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);

                KCL_Q  -= grid->b_tt(a->_name)*Wii(a->_dest->_name) + grid->b_tf(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name) 
                    + grid->g_tf(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);
            }

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
        for (auto b: grid->nodes) {
            Bus* bus = (Bus*) b;
            Constraint KCL_P("KCL_P"+bus->_name);
            Constraint KCL_Q("KCL_Q"+bus->_name);

            if (!bus->_has_gen) {
                for (auto &a: b->get_out()) {
                    KCL_P  += grid->g_ff(a->_name)*Wii(a->_src->_name)
                              +grid->g_ft(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name)+grid->b_ft(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);
                    KCL_Q  += -1*grid->b_ff(a->_name)*Wii(a->_src->_name)- grid->b_ft(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name)
                              +grid->g_ft(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);
                }

                //KCL_P = sum(grid->g_ff.in(b->get_out()), Wii.from(b->get_out()));

                for (auto &a: b->get_in()) {
                    KCL_P  += grid->g_tt(a->_name)*Wii(a->_dest->_name)
                              +grid->g_tf(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name)-grid->b_tf(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);

                    KCL_Q  -= grid->b_tt(a->_name)*Wii(a->_dest->_name) + grid->b_tf(a->_name)*R_Wij(a->_src->_name+","+a->_dest->_name) 
                             + grid->g_tf(a->_name)*Im_Wij(a->_src->_name+","+a->_dest->_name);
                }

                KCL_P += bus->pl();
                KCL_Q -= bus->ql();

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
    PAD_UB -= (grid->tan_th_max).in(bus_pairs)*R_Wij.in(bus_pairs);
    SOCP.add_constraint(PAD_UB <= 0);

    Constraint PAD_LB("PAD_LB");
    PAD_LB =  Im_Wij.in(bus_pairs);
    PAD_LB -= grid->tan_th_min.in(bus_pairs)*R_Wij.in(bus_pairs);
    SOCP.add_constraint(PAD_LB >= 0);

    /* Thermal Limit Constraints */
    Constraint Thermal_Limit_from("Thermal_Limit_from");
    Thermal_Limit_from += power(grid->g_ff.in(grid->arcs)*Wii.from(grid->arcs)
                                + grid->g_ft.in(grid->arcs)*R_Wij.in_pairs(grid->arcs)+ grid->b_ft.in(grid->arcs)*Im_Wij.in_pairs(grid->arcs), 2)
                          + power(grid->g_ft.in(grid->arcs)*Im_Wij.in_pairs(grid->arcs)-grid->b_ff.in(grid->arcs)*Wii.to(grid->arcs) 
                                 - grid->b_ft.in(grid->arcs)*R_Wij.in_pairs(grid->arcs), 2);
    Thermal_Limit_from -= power(grid->S_max.in(grid->arcs),2);
    SOCP.add_constraint(Thermal_Limit_from <= 0);

    Constraint Thermal_Limit_to("Thermal_Limit_to");
    Thermal_Limit_to += power(grid->g_tt.in(grid->arcs)*Wii.from(grid->arcs) + grid->g_tf.in(grid->arcs)*R_Wij.in_pairs(grid->arcs)
                        + grid->b_tf.in(grid->arcs)*Im_Wij.in_pairs(grid->arcs), 2)
                        + power(grid->b_tt.in(grid->arcs)*Wii.to(grid->arcs) + grid->b_tf.in(grid->arcs)*R_Wij.in_pairs(grid->arcs)
                        + grid->g_tf.in(grid->arcs)*Im_Wij.in_pairs(grid->arcs), 2);

    Thermal_Limit_to -= power(grid->S_max.in(grid->arcs),2);
    SOCP.add_constraint(Thermal_Limit_to <= 0);

    Constraint NL("NL");
    NL = Wii(grid->get_ref_bus())*R_Wij(bus_pairs.front()->_name)*Im_Wij(bus_pairs.front()->_name);
    SOCP.add_constraint(NL <= 0);

    solver SCOPF(SOCP, cplex);
    SCOPF.run();
}

#endif
