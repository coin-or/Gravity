//
//  Decompose.cpp
//
//  Gravity
//
//  Created by Guanglei Wang on 6/9/17.
//
//

#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <stdio.h>
#include <cstring>
#include <fstream>
#include "../PowerNet.h"
#include "SCOPF_W.cpp"
#include <gravity/model.h>
#include <gravity/solver.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;
using namespace gravity;

/** INITIALISE SUBPROBLEM MODEL */
double  subproblem(PowerNet* grid, Net* chordal, unsigned c, Net* cliquetree,
                   vector<Bus*>  bag_bus,  vector<Line*> bag_arcs, vector<Gen*> bag_gens,node_pairs* bag_bus_pairs)
{
    cout << "Solving subproblem associated with maximal clique .........." << c << endl;

    if (bag_arcs.size() == 0) {
        return 0;
    }
    Model Subr("Subr");

    // POWER FLOW
    DebugOn("bag_arcs " << c << " has " << bag_arcs.size() << " lines." << endl);
    
    for (auto &a: bag_arcs) {cout << a->_name << endl;}
    DebugOn("bag_bus  " << c << " has " << bag_bus.size() << " bus." << endl);
    DebugOn("bag_gens " << c << " has " << bag_gens.size() << " gens." << endl);

   
   // FLOW VARIABLES DEFINED ON EDGE SET INVOLVING NODES IN THE CLIQUE.
   var<Real> Pf_from("Pf_from", grid->S_max.in(bag_arcs));
   var<Real> Qf_from("Qf_from", grid->S_max.in(bag_arcs));
   var<Real> Pf_to("Pf_to", grid->S_max.in(bag_arcs));
   var<Real> Qf_to("Qf_to", grid->S_max.in(bag_arcs));
   Subr.add_var(Pf_from^(bag_arcs.size()));
   Subr.add_var(Qf_from^(bag_arcs.size()));
   Subr.add_var(Pf_to^(bag_arcs.size()));
   Subr.add_var(Qf_to^(bag_arcs.size()));
   
   var<Real> Pg("Pg", grid->pg_min.in(bag_gens), grid->pg_max.in(bag_gens));
   var<Real> Qg("Qg", grid->qg_min.in(bag_gens), grid->qg_max.in(bag_gens));
   Subr.add_var(Pg^(bag_gens.size()));
   Subr.add_var(Qg^(bag_gens.size()));

   // LIFTED VARIABLES..
   // Note that there are two parts of W
   // First part: W defined on the clique. (bag_bus_pairs).
   // Second part: W defined across clique sets.
   var<Real>  R_Wij("R_Wij", grid->wr_min.in(bag_bus_pairs->_keys), grid->wr_max.in(bag_bus_pairs->_keys));
   var<Real>  Im_Wij("Im_Wij", grid->wi_min.in(bag_bus_pairs->_keys), grid->wi_max.in(bag_bus_pairs->_keys));
   var<Real>  Wii("Wii", grid->w_min.in(bag_bus), grid->w_max.in(bag_bus));
   Subr.add_var(Wii^(bag_bus.size()));
   Subr.add_var(R_Wij^(bag_bus.size()*(bag_bus.size() - 1)/2));
   Subr.add_var(Im_Wij^(bag_bus.size()*(bag_bus.size() - 1)/2));

   R_Wij.initialize_all(1.0);
   Wii.initialize_all(1.001);

   /* Construct the objective function*/
   func_ obj;
   obj += 0;
   Subr.set_objective(min(obj));

   // POWER GENERATION
   for (auto b: bag_bus) {
       if (!b->_active) {
           continue;
       }
       Constraint KCL_P("KCL_P" + b->_name);
       Constraint KCL_Q("KCL_Q" + b->_name);

       /* Power Conservation */
       KCL_P  = sum(Pf_from.in(b->get_out())) + sum(Pf_to.in(b->get_in())) + b->pl()- sum(Pg.in(b->_gen));
       KCL_Q  = sum(Qf_from.in(b->get_out())) + sum(Qf_to.in(b->get_in())) + b->ql()- sum(Qg.in(b->_gen));

       /* Shunts */
       KCL_P +=  grid->gs(b->_name)*Wii(b->_name);
       KCL_Q -=  grid->bs(b->_name)*Wii(b->_name);

       Subr.add_constraint(KCL_P = 0);
       Subr.add_constraint(KCL_Q = 0);
   }
   /* SOCP constraints */
   Constraint SOC("SOC");
   SOC =  power(R_Wij.in(bag_bus_pairs->_keys), 2) + power(Im_Wij.in(bag_bus_pairs->_keys), 2)
          - Wii.from(bag_bus_pairs->_keys)*Wii.to(bag_bus_pairs->_keys);
   Subr.add_constraint(SOC <= 0);

    
    Wii.print(true);
    R_Wij.print(true);
    Im_Wij.print(true);
    Pf_to.print(true);
    Pf_from.print(true);
    Pg.print(true);
    Qg.print(true);

   //AC Power Flow.
   Constraint Flow_P_From("Flow_P_From");
   Flow_P_From += Pf_from.in(bag_arcs);
   Flow_P_From -= grid->g_ff.in(bag_arcs)*Wii.from(bag_arcs);
   Flow_P_From -= grid->g_ft.in(bag_arcs)*R_Wij.in_pairs(bag_arcs);
   Flow_P_From -= grid->b_ft.in(bag_arcs)*Im_Wij.in_pairs(bag_arcs);
   Subr.add_constraint(Flow_P_From = 0);


   Constraint Flow_P_To("Flow_P_To");
   Flow_P_To += Pf_to.in(bag_arcs);
   Flow_P_To -= grid->g_tt.in(bag_arcs)*Wii.to(bag_arcs);
   Flow_P_To -= grid->g_tf.in(bag_arcs)*R_Wij.in_pairs(bag_arcs);
   Flow_P_To += grid->b_tf.in(bag_arcs)*Im_Wij.in_pairs(bag_arcs);
   Subr.add_constraint(Flow_P_To = 0);
//
//   Constraint Flow_Q_From("Flow_Q_From");
//   Flow_Q_From += Qf_from.in(bag_arcs);
//   Flow_Q_From += grid->b_ff.in(bag_arcs)*Wii.from(bag_arcs);
//   Flow_Q_From += grid->b_ft.in(bag_arcs)*R_Wij.in_pairs(bag_arcs);
//   Flow_Q_From += grid->g_ft.in(bag_arcs)*Im_Wij.in_pairs(bag_arcs);
//   Flow_Q_From = 0;
//   Subr.add_constraint(Flow_Q_From);
//
//   Constraint Flow_Q_To("Flow_Q_To");
//   Flow_Q_To += Qf_to.in(bag_arcs);
//   Flow_Q_To += grid->b_tt.in(bag_arcs)*Wii.to(bag_arcs);
//   Flow_Q_To += grid->b_tf.in(bag_arcs)*R_Wij.in_pairs(bag_arcs);
//   Flow_Q_To -= grid->g_tf.in(bag_arcs)*Im_Wij.in_pairs(bag_arcs);
//   Flow_Q_To = 0;
//   Subr.add_constraint(Flow_Q_To);

   //Constraint PAD_UB("PAD_UB");
   //PAD_UB = Im_Wij.in(bag_bus_pairs->_keys);
   //PAD_UB -= (grid->tan_th_max).in(bag_bus_pairs->_keys)*R_Wij.in(bag_bus_pairs->_keys);
   //Subr.add_constraint(PAD_UB <= 0);

   //Constraint PAD_LB("PAD_LB");
   //PAD_LB =  Im_Wij.in(bag_bus_pairs->_keys);
   //PAD_LB -= grid->tan_th_min.in(bag_bus_pairs->_keys)*R_Wij.in(bag_bus_pairs->_keys);
   //Subr.add_constraint(PAD_LB >= 0);

   ////Thermal Limit Constraints */
   //Constraint Thermal_Limit_from("Thermal_Limit_from");
   //Thermal_Limit_from += power(Pf_from.in(bag_arcs),  2) + power(Qf_from.in(bag_arcs), 2);
   //Thermal_Limit_from -= power(grid->S_max.in(bag_arcs), 2);
   //Subr.add_constraint(Thermal_Limit_from <= 0);

   //Constraint Thermal_Limit_to("Thermal_Limit_to");
   //Thermal_Limit_to += power(Pf_to.in(bag_arcs), 2) + power(Qf_to.in(bag_arcs), 2);
   //Thermal_Limit_to -= power(grid->S_max.in(bag_arcs),2);
   //Subr.add_constraint(Thermal_Limit_to <= 0);

   /* Resolve it! */
   solver solve_Subr(Subr,cplex);
   solve_Subr.run();

   return Subr._obj_val;
}

int main (int argc, const char * argv[])
{
    
        // Decompose
    PowerNet* grid = new PowerNet();
    const char* fname;
    fname = "../../data_sets/Power/nesta_case5_pjm.m";
    //fname = "../../data_sets/Power/nesta_case14_ieee.m";
    grid->readgrid(fname);
    
    scopf_W(grid);
    return 0;
    
    // Grid Parameters
    unsigned nb_gen = grid->get_nb_active_gens();
    unsigned nb_lines = grid->get_nb_active_arcs();
    unsigned nb_buses = grid->get_nb_active_nodes();

    DebugOn("nb gens = " << nb_gen << endl);
    DebugOn("nb lines = " << 2*nb_lines << endl);
    DebugOn("nb buses = " << nb_buses << endl);

    // Schedule
    unsigned T = 1;
    param<Real> rate_ramp("rate_ramp");
    param<Real> rate_switch("rate_switch");
    param<Real> min_up("min_up");
    param<Real> min_down("min_down");
    param<Real> cost_up("cost_up");
    param<Real> cost_down("cost_down");
    
    for (auto g: grid->gens) {
        rate_ramp(g->_name) = max(grid->pg_min(g->_name).getvalue(), 0.25*grid->pg_max(g->_name).getvalue());
        rate_switch(g->_name) = max(grid->pg_min(g->_name).getvalue(), 0.25*grid->pg_max(g->_name).getvalue());
    }
    
    min_up = 1;
    min_down = 1;
    cost_up = 50;
    cost_down = 30;

    //grid->time_expand(T);

    /** Clique tree decomposition **/
    Net* chordal = grid->get_chordal_extension();
    grid->get_clique_tree();
    const unsigned nb_cliques = grid->_bags.size();

    vector<vector<Bus*>> bag_bus; // each clique contains just nodes, not buses! Fixed this by modifying the bag definition.
    vector<vector<Gen*>> bag_gens;
    vector<vector<Line*>> bag_arcs; //bag_arcs contains the arcs of the power grid while variables associated with W are defined on chordal graph.
    vector<vector<Line*>> bag_arcs_S; //arcs across
    vector<gravity::node_pairs*> bag_bus_pairs; // bus_pairs in each bag.

    for (int c = 0; c < nb_cliques; c++) {
        vector<Bus*> bag_B;
        vector<Gen*> bag_G;
        vector<Line*> bag_A;
        vector<Line*> bag_AS;
        node_pairs* bag_BP = new node_pairs();
        for (int i = 0; i < grid->_bags[c].size(); i++) {
            Bus* b = (Bus*) grid->get_node(grid->_bags[c].at(i)->_name);
            if (b !=nullptr) {
                bag_B.push_back(b);
            }
            if (b->_has_gen) {
                bag_G.insert(bag_G.end(), b->_gen.begin(), b->_gen.end());
            }
            for (int j = i+1; j < grid->_bags[c].size(); j++) {
                Line* a = (Line*)grid->get_arc(b, grid->get_node(grid->_bags[c].at(j)->_name));
                if (a != nullptr) {
                    bag_A.push_back(a);
                    bag_BP->_keys.push_back(new index_pair(index_(a->_src->_name), index_(a->_dest->_name), a->_active));
                }
            }
        }
        bag_bus.push_back(bag_B);
        bag_gens.push_back(bag_G);
        bag_arcs.push_back(bag_A);
        bag_bus_pairs.push_back(bag_BP);
    }

    bag_bus.resize(nb_cliques);
    bag_arcs.resize(nb_cliques);
    bag_gens.resize(nb_cliques);

    auto cliquetree = grid->get_clique_tree();


///////////////////////////////// INITIALIZATION ///////////////////////////////////////////
    vector<double> value_dual;
    double val = 0;
    double dual = 0;
//    for (int c = 0; c < nb_cliques; c++) {
    int c = 1;
        val = subproblem(grid, chordal, c, cliquetree,
                         bag_bus[c],  bag_arcs[c],bag_gens[c], bag_bus_pairs[c]);
        value_dual.push_back(val);
        dual +=val;
 //   }
    value_dual.resize(nb_cliques);
    cout << "................  Initialization value:  " << dual <<endl;

    return 0;
}
