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
#include <gravity/model.h>
#include <gravity/solver.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;
using namespace gravity;

/** INITIALISE SUBPROBLEM MODEL */
double  subproblem(PowerNet* grid,Net* chordal, unsigned c, Net* cliquetree,
                   vector<Bus*> bag_bus, vector<Gen*> bag_gens, vector<Arc*> bag_arcs, node_pairs* bag_bus_pairs,
                   vector<param<Real>>& R_lambda_sep, vector<param<Real>>& Im_lambda_sep, vector<param<Real>>& mu_sep,
                   vector<param<Real>>& kappa_sep, vector<param<Real>>& eta_sep,
                   param<Real>& rate_ramp, param<Real>& rate_switch, param<Real>& min_up, param<Real>& min_down,
                   param<Real>& cost_up, param<Real>& cost_down,
                   param<Real>& Wii_log, param<Real>& R_Wij_log, param<Real>& Im_Wij_log,
                   param<Real>& Pg_log, param<Real>& Qg_log, param<Real>& On_off_log )
{
    cout << "Solving subproblem associated with maximal clique .........." << c << endl;

    if (bag_arcs.size() == 0) {
        return 0;
    }

    Model Subr("Subr");
    
    // POWER FLOW
    DebugOn("bag_arcs " << c << " has " << bag_arcs.size() << " lines." << endl);
    DebugOn("bag_bus  " << c << " has " << bag_bus.size() << " bus." << endl);
    DebugOn("bag_gens " << c << " has " << bag_gens.size() << " gens" << endl);

    var<Real> Pf_from("Pf_from", grid->S_max.in(bag_arcs));
    var<Real> Qf_from("Qf_from", grid->S_max.in(bag_arcs));
    var<Real> Pf_to("Pf_to", grid->S_max.in(bag_arcs));
    var<Real> Qf_to("Qf_to", grid->S_max.in(bag_arcs));
    Subr.add_var(Pf_from^(bag_arcs.size()));
    Subr.add_var(Qf_from^(bag_arcs.size()));
    Subr.add_var(Pf_to^(bag_arcs.size()));
    Subr.add_var(Qf_to^(bag_arcs.size()));

    // LIFTED VARIABLES.
    var<Real>  R_Wij("R_Wij", grid->wr_min.in(bag_bus_pairs->_keys), grid->wr_max.in(bag_bus_pairs->_keys));
    var<Real>  Im_Wij("Im_Wij", grid->wi_min.in(bag_bus_pairs->_keys), grid->wi_max.in(bag_bus_pairs->_keys));
    var<Real>  Wii("Wii", grid->w_min.in(bag_bus), grid->w_max.in(bag_bus));
    Subr.add_var(Wii^(bag_bus.size()));
    Subr.add_var(R_Wij^(grid->_bags.size()*(bag_bus.size() - 1)/2));
    Subr.add_var(Im_Wij^(grid->_bags.size()*(bag_bus.size() - 1)/2));

    R_Wij.initialize_all(1.0);
    Wii.initialize_all(1.001);

    // COMMITMENT VARIABLES
    var<Real>  On_off("On_off", 0, 1);
    var<Real>  Start_up("Start_up", 0, 1);
    var<Real>  Shut_down("Shut_down", 0, 1);
    Subr.add_var(On_off^(bag_gens.size()));
    Subr.add_var(Start_up^(bag_gens.size()));
    Subr.add_var(Shut_down^(bag_gens.size()));

    /* Construct the objective function*/
    func_ obj;
    Node* Cr = cliquetree->get_node(to_string(c));
    Arc* arc = nullptr;
    Bus* bus = nullptr;
    for (auto a: Cr->get_out()) {
        Debug("a->_intersection.size " << a->_intersection.size() << endl);
        for (int i = 0; i < a->_intersection.size(); i ++) {
            bus = (Bus*)a->_intersection.at(i);
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
            string name;
            for (int j = i + 1; j < a->_intersection.size(); j ++) {
                arc = chordal->get_arc(bus, a->_intersection.at(j));
                name = arc->_src->_name + "," + arc->_dest->_name;
                obj -= R_lambda_sep[a->_id](name)*R_Wij(name);
                obj -= Im_lambda_sep[a->_id](name)*Im_Wij(name);
            }
        }
    }
    var<Real> Pg("Pg", grid->pg_min.in(bag_gens), grid->pg_max.in(bag_gens));
    var<Real> Qg("Qg", grid->qg_min.in(bag_gens), grid->qg_max.in(bag_gens));
    Subr.add_var(Pg^(bag_gens.size()));
    Subr.add_var(Qg^(bag_gens.size()));

    for (auto g: bag_gens) {
        if (g->_active) {
             obj += grid->c1(g->_name)*Pg(g->_name) + grid->c2(g->_name)*Pg(g->_name)*Pg(g->_name) + grid->c0(g->_name);
              //obj += cost_up.getvalue()*Start_up(g->_name, l)+ cost_down.getvalue()*Shut_down(g->_name, l);
            }
        }
    for (auto a: Cr->get_out()) {
        for (int i = 0; i < a->_intersection.size(); i++) {
            bus = (Bus*) a->_intersection.at(i);
            for (auto g: bus->_gen) {
                if (g->_active) {
                    obj += mu_sep[a->_id](g->_name)*Pg(g->_name);
                    obj += kappa_sep[a->_id](g->_name)*Qg(g->_name);
                    obj += eta_sep[a->_id](g->_name)*On_off(g->_name);
                }
            }
        }
    }

    for (auto a: Cr->get_in()) {
        for (int i = 0; i < a->_intersection.size(); i++) {
            bus = (Bus*) a->_intersection.at(i);
            for (auto g: bus->_gen) {
                if (g->_active) {
                    obj -= mu_sep[a->_id](g->_name)*Pg(g->_name);
                    obj -= kappa_sep[a->_id](g->_name)*Qg(g->_name);
                    obj -= eta_sep[a->_id](g->_name)*On_off(g->_name);
                }
            }
        }
    }
    Subr.set_objective(min(obj));

    // POWER GENERATION
    if (bag_gens.size() > 0) {
        //KCL
       for (auto bus: bag_bus) {
           if (!bus->_active) {
               continue;
           }
           Constraint KCL_P("KCL_P"+bus->_name+ "time_");
           Constraint KCL_Q("KCL_Q"+bus->_name+ "time_");

           /* Power Conservation */
           vector<Arc*> out;
           /* Add comments here, why are we using set_intersection?*/
           /* bus->get_out() includes all out lines of the grid, but here we need the out lines of the bag.*/
           set_intersection(bus->get_out().begin(), bus->get_out().end(), bag_arcs.begin(), bag_arcs.end(), back_inserter(out));
           vector<Arc*> in;
           set_intersection(bus->get_in().begin(), bus->get_in().end(), bag_arcs.begin(), bag_arcs.end(), back_inserter(in));

           KCL_P  = sum(Pf_from.in(out))+ sum(Pf_to.in(in)) + bus->pl() - sum(Pg.in(bus->_gen));
           KCL_Q  = sum(Qf_from.in(out)) + sum(Qf_to.in(in))+ bus->ql() - sum(Qg.in(bus->_gen));

           /* Shunts */
           KCL_P +=  bus->gs()*Wii(bus->_name);
           KCL_Q -=  bus->bs()*Wii(bus->_name);

           Subr.add_constraint(KCL_P = 0);
           Subr.add_constraint(KCL_Q = 0);
       }
        /* Commitment constraints */
        // Inter-temporal constraints
        //for (int t = 1; t < T; t++) {
        //    Constraint MC1("MC1_"+ to_string(t));
        //    Constraint MC2("MC2_"+ to_string(t));
        //    MC1 = On_off.in_at(bag_gens, t) - On_off.in_at(bag_gens, t-1)-  Start_up.in_at(bag_gens, t);
        //    MC2 = On_off.in_at(bag_gens, t-1) - On_off.in_at(bag_gens, t) - Shut_down.in_at(bag_gens, t);
        //    Subr.add_constraint(MC1 <= 0);
        //    Subr.add_constraint(MC2 <= 0);
        //}

        //// Min-up constraints
        //for (int t = 1; t < T; t++) {
        //    Constraint Min_up1("Min_up1_"+ to_string(t));
        //    Min_up1 = On_off.in_at(bag_gens, t) - On_off.in_at(bag_gens, t-1) - Start_up.in_at(bag_gens, t) + Shut_down.in_at(bag_gens, t);
        //    Subr.add_constraint(Min_up1 = 0);
        //}

        //for (int t = min_up.getvalue(); t < T; t++) {
        //    Constraint Min_Up("Min_Up_constraint" + to_string(t));
        //    for (int l = t-min_up.getvalue()+1; l < t +1; l++) {
        //        Min_Up   += Start_up.in_at(bag_gens, l);
        //    }
        //    Min_Up -= On_off.in_at(bag_gens, t);
        //    Subr.add_constraint(Min_Up <= 0);
        //}

        //for (int t = min_down.getvalue(); t < T; t++) {
        //    Constraint Min_Down("Min_Down_constraint" + to_string(t));
        //    for (int l = t-min_down.getvalue()+1; l < t +1; l++) {
        //        Min_Down   += Shut_down.in_at(bag_gens, l);
        //    }
        //    Min_Down -= 1 - On_off.in_at(bag_gens, t);
        //    Subr.add_constraint(Min_Down <= 0);
        //}

        ////Ramp rate
        //Constraint Production_P_LB("Production_P_LB");
        //Constraint Production_P_UB("Production_P_UB");
        //Constraint Production_Q_LB("Production_Q_LB");
        //Constraint Production_Q_UB("Production_Q_UB");

        //Production_P_UB = Pg.in(bag_gens) - grid->pg_max.in(bag_gens)*On_off.in(bag_gens,T);
        //Production_P_LB = Pg.in(bag_gens) - grid->pg_min.in(bag_gens)*On_off.in(bag_gens,T);
        //Subr.add_constraint(Production_P_UB <=0);
        //Subr.add_constraint(Production_P_LB >= 0);

        ////grid->qg_max.print(true);
        ////grid->qg_min.print(true);

        //Production_Q_UB = Qg.in(bag_gens) - grid->qg_max.in(bag_gens)*On_off.in(bag_gens,T);
        //Production_Q_LB = Qg.in(bag_gens) - grid->qg_min.in(bag_gens)*On_off.in(bag_gens,T);
        //Subr.add_constraint(Production_Q_UB <= 0);
        //Subr.add_constraint(Production_Q_LB >= 0);

        //for (int t = 1; t < T; t++) {
        //    Constraint Ramp_up("Ramp_up_constraint" + to_string(t));
        //    Constraint Ramp_down("Ramp_down_constraint" + to_string(t));

        //    Ramp_up = Pg.in_at(bag_gens, t);
        //    Ramp_up -= Pg.in_at(bag_gens, t-1);
        //    Ramp_up -= rate_ramp*On_off.in_at(bag_gens, t-1);
        //    Ramp_up -= rate_switch*(1 - On_off.in_at(bag_gens, t));

        //    Ramp_down = Pg.in_at(bag_gens, t-1);
        //    Ramp_down -= Pg.in_at(bag_gens, t);
        //    Ramp_down -= rate_ramp*On_off.in_at(bag_gens, t);
        //    Ramp_down -= rate_switch*(1 - On_off.in_at(bag_gens, t-1));

        //    Subr.add_constraint(Ramp_up <= 0);
        //    Subr.add_constraint(Ramp_down <= 0);
        //}
    }
    else {
        //KCL
         for (auto bus: bag_bus) {
             if (!bus->_active) {
                 continue;
             }
             Constraint KCL_P("KCL_P"+bus->_name+ "time_");
             Constraint KCL_Q("KCL_Q"+bus->_name+ "time_");

             /* Power Conservation */
             vector<Arc*> out;
             set_intersection(bus->get_out().begin(), bus->get_out().end(), bag_arcs.begin(), bag_arcs.end(), back_inserter(out));
             vector<Arc*> in;
             set_intersection(bus->get_in().begin(), bus->get_in().end(), bag_arcs.begin(), bag_arcs.end(), back_inserter(in));
             KCL_P  = sum(Pf_from.in(out))+ sum(Pf_to.in(in)) + bus->pl();
             KCL_Q  = sum(Qf_from.in(out)) + sum(Qf_to.in(in))+ bus->ql();

             /* Shunts */
             KCL_P +=  bus->gs()*Wii(bus->_name);
             KCL_Q -=  bus->bs()*Wii(bus->_name);

             Subr.add_constraint(KCL_P = 0);
             Subr.add_constraint(KCL_Q = 0);
         }
    }
    /* SOCP constraints */
    Constraint SOC("SOC");
    SOC =  power(R_Wij.in(bag_bus_pairs->_keys), 2) + power(Im_Wij.in(bag_bus_pairs->_keys), 2)
            - Wii.from(bag_bus_pairs->_keys)*Wii.to(bag_bus_pairs->_keys) ;
    Subr.add_constraint(SOC <= 0);

    //AC Power Flow.
    Constraint Flow_P_From("Flow_P_From");
    Flow_P_From += Pf_from.in(bag_arcs);
    Flow_P_From -= grid->g_ff.in(bag_arcs)*Wii.from(bag_arcs);
    Flow_P_From -= grid->g_ft.in(bag_arcs)*R_Wij.in_pairs(bag_arcs);
    Flow_P_From -= grid->b_ft.in(bag_arcs)*Im_Wij.in_pairs(bag_arcs);
    Flow_P_From = 0;
    Subr.add_constraint(Flow_P_From);
//
    Constraint Flow_P_To("Flow_P_To");
    Flow_P_To += Pf_to.in(bag_arcs);
    Flow_P_To -= grid->g_tt.in(bag_arcs)*Wii.to(bag_arcs);
    Flow_P_To -= grid->g_tf.in(bag_arcs)*R_Wij.in_pairs(bag_arcs);
    Flow_P_To += grid->b_tf.in(bag_arcs)*Im_Wij.in_pairs(bag_arcs);
    Flow_P_To = 0;
    Subr.add_constraint(Flow_P_To);
//
    Constraint Flow_Q_From("Flow_Q_From");
    Flow_Q_From += Qf_from.in(bag_arcs);
    Flow_Q_From += grid->b_ff.in(bag_arcs)*Wii.from(bag_arcs);
    Flow_Q_From += grid->b_ft.in(bag_arcs)*R_Wij.in_pairs(bag_arcs);
    Flow_Q_From += grid->g_ft.in(bag_arcs)*Im_Wij.in_pairs(bag_arcs);
    Flow_Q_From = 0;
    Subr.add_constraint(Flow_Q_From);
//
    Constraint Flow_Q_To("Flow_Q_To");
    Flow_Q_To += Qf_to.in(bag_arcs);
    Flow_Q_To += grid->b_tt.in(bag_arcs)*Wii.to(bag_arcs);
    Flow_Q_To += grid->b_tf.in(bag_arcs)*R_Wij.in_pairs(bag_arcs);
    Flow_Q_To -= grid->g_tf.in(bag_arcs)*Im_Wij.in_pairs(bag_arcs);
    Flow_Q_To = 0;
    Subr.add_constraint(Flow_Q_To);

    //NOTE THAT WE SHOULD USE BUS PAIRS!!!
    Constraint PAD_UB("PAD_UB");
    PAD_UB = Im_Wij.in(bag_bus_pairs->_keys);
    PAD_UB -= (grid->tan_th_max).in(bag_bus_pairs->_keys)*R_Wij.in(bag_bus_pairs->_keys);
    Subr.add_constraint(PAD_UB <= 0);

    Constraint PAD_LB("PAD_LB");
    PAD_LB =  Im_Wij.in(bag_bus_pairs->_keys);
    PAD_LB -= grid->tan_th_min.in(bag_bus_pairs->_keys)*R_Wij.in(bag_bus_pairs->_keys);
    Subr.add_constraint(PAD_LB >= 0);

    //Thermal Limit Constraints */
    Constraint Thermal_Limit_from("Thermal_Limit_from");
    Thermal_Limit_from += power(Pf_from.in(bag_arcs),  2) + power(Qf_from.in(bag_arcs), 2);
    Thermal_Limit_from -= power(grid->S_max.in(bag_arcs), 2);
    Subr.add_constraint(Thermal_Limit_from <= 0);

    Constraint Thermal_Limit_to("Thermal_Limit_to");
    Thermal_Limit_to += power(Pf_to.in(bag_arcs), 2) + power(Qf_to.in(bag_arcs), 2);
    Thermal_Limit_to -= power(grid->S_max.in(bag_arcs),2);
    Subr.add_constraint(Thermal_Limit_to <= 0);

    /* Resolve it! */
    solver solve_Subr(Subr,cplex);
    solve_Subr.run();

    Wii_log =  (*(var<Real>*) Subr.get_var("Wii"));
    R_Wij_log =  (*(var<Real>*) Subr.get_var("R_Wij"));
    Im_Wij_log =  (*(var<Real>*) Subr.get_var("Im_Wij"));
    Pg_log = (*(var<Real>*) Subr.get_var("Pg"));
    Qg_log = (*(var<Real>*) Subr.get_var("Qg"));
    On_off_log = (*(var<Real>*) Subr.get_var("On_off"));
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

    vector<vector<Bus*>> bag_bus; //  each clique contains just nodes, not buses! Fixed this by modifying the bag definition.
    vector<vector<Gen*>> bag_gens;
    vector<vector<Arc*>> bag_arcs; //bag_arcs contains the arcs of the power grid while variables associated with W are defined on chordal graph.
    vector<gravity::node_pairs*> bag_bus_pairs; // bus_pairs in each bag.

    for (int c = 0; c < nb_cliques; c++) {
        vector<Bus*> bag_B;
        vector<Gen*> bag_G;
        vector<Arc*> bag_A;
        node_pairs* bag_BP = new node_pairs();
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
                if (a != nullptr){
                    bag_A.push_back(a);
                    bag_BP->_keys.push_back(new index_pair(index_(a->_src->_name), index_(a->_dest->_name), a->_active));
                }
            }
        }
        bag_bus.push_back(bag_B);
        bag_gens.push_back(bag_G);
        bag_arcs.push_back(bag_A);
        bag_bus_pairs.push_back(bag_BP);
        DebugOn("bag " << c << " has " << bag_gens[c].size() << " generators. " << endl;)
        DebugOn("bag " << c << " has " << bag_arcs[c].size() << " line " << endl;)
    }

    bag_bus.resize(nb_cliques);
    bag_arcs.resize(nb_cliques);
    bag_gens.resize(nb_cliques);

    auto cliquetree = grid->get_clique_tree();
    
///////////////// DEFINE LAGRANGE MULTIPLIERS  ////////////////////////////////
    vector<param<Real>> R_lambda_in;
    vector<param<Real>> Im_lambda_in;
    vector<param<Real>> mu_in;
    vector<param<Real>> kappa_in;
    vector<param<Real>> eta_in;

    vector<param<Real>> R_lambda_out;
    vector<param<Real>> Im_lambda_out;
    vector<param<Real>> mu_out;
    vector<param<Real>> kappa_out;
    vector<param<Real>> eta_out;

    vector<param<Real>> R_lambda_sep;
    vector<param<Real>> Im_lambda_sep;
    vector<param<Real>> mu_sep;
    vector<param<Real>> kappa_sep;
    vector<param<Real>> eta_sep;

    vector<param<Real>> R_lambda_grad;
    vector<param<Real>> Im_lambda_grad;
    vector<param<Real>> mu_grad;
    vector<param<Real>> kappa_grad;
    vector<param<Real>> eta_grad;

    for (auto a: cliquetree->arcs) {
        auto l = a->_id;
        param<Real> R_lambda_arc_in("R_lambda_arc_in" + to_string(l));
        param<Real> Im_lambda_arc_in("Im_lambda_arc_in" + to_string(l));
        param<Real> mu_arc_in("R_mu_arc_in" + to_string(l));
        param<Real> kappa_arc_in("kappa_arc_in" + to_string(l));
        param<Real> eta_arc_in("eta_arc_in" + to_string(l));
        R_lambda_arc_in^(a->_weight*(a->_weight-1)/2);
        Im_lambda_arc_in^(a->_weight*(a->_weight-1)/2);
        mu_arc_in^(a->_weight);
        kappa_arc_in^(a->_weight);
        eta_arc_in^(a->_weight);

        param<Real> R_lambda_arc_out("R_lambda_arc_out" + to_string(l));
        param<Real> Im_lambda_arc_out("Im_lambda_arc_out" + to_string(l));
        param<Real> mu_arc_out("mu_arc_out" + to_string(l));
        param<Real> kappa_arc_out("kappa_arc_out" + to_string(l));
        param<Real> eta_arc_out("eta_arc_out" + to_string(l));
        R_lambda_arc_out^(a->_weight*(a->_weight-1)/2);
        Im_lambda_arc_out^(a->_weight*(a->_weight-1)/2);
        mu_arc_out^(a->_weight);
        kappa_arc_out^(a->_weight);
        eta_arc_out^(a->_weight);

        param<Real> R_lambda_arc_sep("R_lambda_arc_sep" + to_string(l));
        param<Real> Im_lambda_arc_sep("Im_lambda_arc_sep" + to_string(l));
        param<Real> mu_arc_sep("mu_arc_sep" + to_string(l));
        param<Real> kappa_arc_sep("kappa_arc_sep" + to_string(l));
        param<Real> eta_arc_sep("eta_arc_sep" + to_string(l));
        R_lambda_arc_sep^(a->_weight*(a->_weight-1)/2);
        Im_lambda_arc_sep^(a->_weight*(a->_weight-1)/2);
        mu_arc_sep^(a->_weight);
        kappa_arc_sep^(a->_weight);
        eta_arc_sep^(a->_weight);

        param<Real> R_lambda_arc_grad("R_lambda_arc_grad" + to_string(l));
        param<Real> Im_lambda_arc_grad("Im_lambda_arc_grad" + to_string(l));
        param<Real> mu_arc_grad("mu_arc_grad" + to_string(l));
        param<Real> kappa_arc_grad("kappa_arc_grad" + to_string(l));
        param<Real> eta_arc_grad("eta_arc_grad" + to_string(l));

        for (int i = 0; i < a->_weight; i++) {
            mu_arc_sep(i) = 0;
        }

        R_lambda_in.push_back(R_lambda_arc_in);
        Im_lambda_in.push_back(Im_lambda_arc_in);
        mu_in.push_back(mu_arc_in);
        kappa_in.push_back(kappa_arc_in);
        eta_in.push_back(eta_arc_in);

        R_lambda_out.push_back(R_lambda_arc_out);
        Im_lambda_out.push_back(Im_lambda_arc_out);
        mu_out.push_back(mu_arc_out);
        kappa_out.push_back(kappa_arc_out);
        eta_out.push_back(eta_arc_out);

        R_lambda_sep.push_back(R_lambda_arc_sep);
        Im_lambda_sep.push_back(Im_lambda_arc_sep);
        mu_sep.push_back(mu_arc_sep);
        kappa_sep.push_back(kappa_arc_sep);
        eta_sep.push_back(eta_arc_sep);

        R_lambda_grad.push_back(R_lambda_arc_grad);
        Im_lambda_grad.push_back(Im_lambda_arc_grad);
        mu_grad.push_back(mu_arc_grad);
        kappa_grad.push_back(kappa_arc_grad);
        eta_grad.push_back(eta_arc_grad);
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
                    mu_sep[l](g->_name) = 0;
                    mu_out[l](g->_name) = 0;
                    mu_in[l] (g->_name) = 0;
                    eta_sep[l](g->_name) = 0;
                    eta_out[l](g->_name) = 0;
                    eta_in[l] (g->_name) = 0;
                    kappa_sep[l](g->_name) = 0;
                    kappa_out[l](g->_name) = 0;
                    kappa_in[l](g->_name) = 0;
                }
            }
        }
    }
    /* Why doing the resize while we have used the ^ operator? */
    // ^ operator is for each parameter; resizing if for each vector; am I missing sth?
    R_lambda_in.resize(cliquetree->arcs.size());
    Im_lambda_in.resize(cliquetree->arcs.size());
    mu_in.resize(cliquetree->arcs.size());
    kappa_in.resize(cliquetree->arcs.size());
    eta_in.resize(cliquetree->arcs.size());
    R_lambda_out.resize(cliquetree->arcs.size());
    Im_lambda_out.resize(cliquetree->arcs.size());
    mu_out.resize(cliquetree->arcs.size());
    kappa_out.resize(cliquetree->arcs.size());
    eta_out.resize(cliquetree->arcs.size());
    R_lambda_sep.resize(cliquetree->arcs.size());
    Im_lambda_sep.resize(cliquetree->arcs.size());
    mu_sep.resize(cliquetree->arcs.size());
    kappa_sep.resize(cliquetree->arcs.size());
    eta_sep.resize(cliquetree->arcs.size());
    R_lambda_grad.resize(cliquetree->arcs.size());
    Im_lambda_grad.resize(cliquetree->arcs.size());
    mu_grad.resize(cliquetree->arcs.size());
    kappa_grad.resize(cliquetree->arcs.size());
    eta_grad.resize(cliquetree->arcs.size());

/////////////////////////////////// INITIALISE MAIN ///////////////////////////////////
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
    vector<var<Real>> R_lambda_var;
    vector<var<Real>> Im_lambda_var;
    vector<var<Real>> mu_var;
    vector<var<Real>> kappa_var;
    vector<var<Real>> eta_var;

    Master.add_var(gamma_C^nb_cliques);
    for (auto a: cliquetree->arcs) {
        var<Real> R_lambda("R_lambda_arc_" + to_string(a->_id));
        var<Real> Im_lambda("Im_lambda_arc_" + to_string(a->_id));
        var<Real> mu("mu_arc_" + to_string(a->_id));
        var<Real> kappa("kappa_arc_" + to_string(a->_id));
        var<Real> eta("eta_arc_" + to_string(a->_id));
        Master.add_var(R_lambda^(a->_weight*(a->_weight -1)/2));
        Master.add_var(Im_lambda^(a->_weight*(a->_weight -1)/2));
        Master.add_var(mu^a->_weight);
        Master.add_var(kappa^a->_weight);
        Master.add_var(eta^a->_weight);

        R_lambda_var.push_back(R_lambda);
        Im_lambda_var.push_back(Im_lambda);
        mu_var.push_back(mu);
        kappa_var.push_back(kappa);
        eta_var.push_back(eta);
    }
    R_lambda_var.resize(cliquetree->arcs.size());
    Im_lambda_var.resize(cliquetree->arcs.size());
    mu_var.resize(cliquetree->arcs.size());
    kappa_var.resize(cliquetree->arcs.size());
    eta_var.resize(cliquetree->arcs.size());

    /////////** OBJ*//////////////
    func_ master_obj = sum(gamma_C);
    Master.set_objective(max(master_obj));
    double bound = 100000000;

    Constraint UB;
    UB = sum(gamma_C) - bound;
    Master.add_constraint(UB <= 0);

////////////////  CONVERGENCE INFORMATION /////////////////////////
    unsigned iter_limit;
    cout << "Enter the limit of the number of iterations: ";
    cin >> iter_limit;
    cout << endl;

    double alpha = 0.5;
    double LBlog[iter_limit];
    double UBlog[iter_limit];

    for(int iter = 0; iter < iter_limit; iter++) {
        LBlog[iter] = 0.0;
    }

    double LDlog[nb_cliques];

    // LOG OF SOLUTIONS
    // Log here means the previous primal and dual solution.
    vector<param<Real>> R_lambda_log;
    vector<param<Real>> Im_lambda_log;
    vector<param<Real>> mu_log;
    vector<param<Real>> kappa_log;
    vector<param<Real>> eta_log;

    vector<param<Real>> R_Wij_log;
    vector<param<Real>> Im_Wij_log;
    vector<param<Real>> Wii_log;
    vector<param<Real>> Pg_log;
    vector<param<Real>> Qg_log;
    vector<param<Real>> On_off_log;

    for (auto a: cliquetree->arcs) {
        int l = a->_id;
        param<Real> R_lambda_C_log("R_lambda_C_log" + to_string(l));
        param<Real> Im_lambda_C_log("Im_lambda_C_log" + to_string(l));
        param<Real> mu_C_log("mu_C_log" + to_string(l));
        param<Real> kappa_C_log("kappa_C_log" + to_string(l));
        param<Real> eta_C_log("eta_C_log" + to_string(l));


        R_lambda_C_log^(a->_weight*(a->_weight-1)/2);
        Im_lambda_C_log^(a->_weight*(a->_weight-1)/2);
        mu_C_log^(a->_weight);
        kappa_C_log^(a->_weight);
        eta_C_log^(a->_weight);

        R_lambda_log.push_back(R_lambda_C_log);
        Im_lambda_log.push_back(Im_lambda_C_log);
        mu_log.push_back(mu_C_log);
        kappa_log.push_back(kappa_C_log);
        eta_log.push_back(eta_C_log);
    }
    R_lambda_log.resize(nb_cliques);
    Im_lambda_log.resize(nb_cliques);
    mu_log.resize(nb_cliques);
    kappa_log.resize(nb_cliques);
    eta_log.resize(nb_cliques);

    for (int c = 0; c < nb_cliques; c++) {
        param<Real> Im_Wij_C_log("Im_Wij_C_log" + to_string(c));
        param<Real> R_Wij_C_log("R_Wij_C_log" + to_string(c));
        param<Real> Wii_C_log("Wii_C_log" + to_string(c));
        param<Real> Pg_C_log("Pg_C_log" + to_string(c));
        param<Real> Qg_C_log("Qg_C_log" + to_string(c));
        param<Real>  On_off_C_log("On_off_C_log" + to_string(c));
        R_Wij_C_log^(T*grid->_bags[c].size()*(bag_bus[c].size() - 1)/2);
        Im_Wij_C_log^(T*grid->_bags[c].size()*(bag_bus[c].size() - 1)/2);
        Wii_C_log^(T*grid->_bags[c].size());
        Pg_C_log^(T*bag_gens[c].size());
        Qg_C_log^(T*bag_gens[c].size());
        On_off_C_log^(T*bag_gens[c].size());

        R_Wij_log.push_back(R_Wij_C_log);
        Im_Wij_log.push_back(Im_Wij_C_log);
        Wii_log.push_back(Wii_C_log);
        Pg_log.push_back(Pg_C_log);
        Qg_log.push_back(Qg_C_log);
        On_off_log.push_back(On_off_C_log);
    }
    
    R_Wij_log.resize(nb_cliques);
    Im_Wij_log.resize(nb_cliques);
    Wii_log.resize(nb_cliques);
    Pg_log.resize(nb_cliques);
    Qg_log.resize(nb_cliques);
    On_off_log.resize(nb_cliques);

///////////////////////////////// INITIALIZATION ///////////////////////////////////////////
    double wall0 = get_wall_time();
    double cpu0  = get_cpu_time();
    vector<double> value_dual;
    double val = 0;
    double dual = 0;
    for (int c = 0; c < nb_cliques; c++) {
        val = subproblem(grid, chordal, c, cliquetree,                          
                         bag_bus[c], bag_gens[c], bag_arcs[c],bag_bus_pairs[c],
                         R_lambda_sep, Im_lambda_sep,  mu_sep, kappa_sep, eta_sep,
                         rate_ramp, rate_switch, min_up, min_down, cost_up, cost_down,
                         Wii_log[c], R_Wij_log[c], Im_Wij_log[c], Pg_log[c], Qg_log[c], On_off_log[c]);
        value_dual.push_back(val);
        gamma_in(c) = val;
        dual +=val;
    }
    value_dual.resize(nb_cliques);
    cout << "................  Initialization value:  " << dual <<endl;

/////////////////// APPEND MORE CONSTRAINTS TO MAIN //////////////////////////////////
    Bus* node = nullptr;
    if (iter_limit > 0) {
        for (int c = 0; c < nb_cliques; c++)
        {
            Constraint Concavity("Iter_0_Concavity_" + to_string(c));
            Concavity += gamma_C(c) - value_dual[c];
            for (auto a: cliquetree->get_node(to_string(c))->get_out()) {
                for (int i = 0; i < a->_intersection.size(); i ++) {
                    bus = (Bus*)a->_intersection.at(i);
                    for (int j = i + 1; j < a->_intersection.size(); j ++) {
                        arc = chordal->get_arc(bus, a->_intersection.at(j)); // we have to use the arc name.
                        Concavity += (R_lambda_sep[a->_id](arc->_name) - R_lambda_var[a->_id](arc->_name))*R_Wij_log[c](arc->_name);
                        Concavity += (Im_lambda_sep[a->_id](arc->_name) - Im_lambda_var[a->_id](arc->_name))*Im_Wij_log[c](arc->_name);
                    }
                    for (auto g: bus->_gen) {
                        if (g->_active) {
                            Concavity +=(mu_sep[a->_id](g->_name) - mu_var[a->_id](g->_name))*Pg_log[c](g->_name);
                            Concavity +=(kappa_sep[a->_id](g->_name) - kappa_var[a->_id](g->_name))*Qg_log[c](g->_name);
                            Concavity +=(eta_sep[a->_id](g->_name) - eta_var[a->_id](g->_name))*On_off_log[c](g->_name);
                        }
                    }
                }
            }
            for (auto a: cliquetree->get_node(to_string(c))->get_in()) {
                for (int i = 0; i < a->_intersection.size(); i ++) {
                    node = (Bus*)a->_intersection.at(i);
                    for (int j = i + 1; j < a->_intersection.size(); j ++) {
                        arc = chordal->get_arc(node, a->_intersection.at(j)); // we have to use the arc name.
                        Concavity += (R_lambda_var[a->_id](arc->_name)- R_lambda_sep[a->_id](arc->_name) )*R_Wij_log[c](arc->_name);
                        Concavity += (Im_lambda_var[a->_id](arc->_name)-Im_lambda_sep[a->_id](arc->_name))*Im_Wij_log[c](arc->_name);
                    }
                }
            }
            Master.add_constraint(Concavity <= 0);
        }
    }

    solver solve_Master(Master, cplex);
    solve_Master.run();

    // initialise the outer point. 
    gamma_out = (*(var<Real>*) Master.get_var("gamma_C"));
    
    for (auto a: cliquetree->arcs) {
        R_lambda_out[a->_id] = (*(var<Real>*) Master.get_var("R_lambda_arc_" + to_string(a->_id)));
        Im_lambda_out[a->_id] = (*(var<Real>*) Master.get_var("Im_lambda_arc_"+ to_string(a->_id)));
        mu_out[a->_id] = (*(var<Real>*) Master.get_var("mu_arc_"+ to_string(a->_id)));
        kappa_out[a->_id] = (*(var<Real>*) Master.get_var("kappa_arc_"+ to_string(a->_id)));
        eta_out[a->_id] = (*(var<Real>*) Master.get_var("eta_arc_"+ to_string(a->_id)));
    }

    cout << "................  Initialization of Master problem ............... "  <<endl;
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
            for (int i = 0; i < a->_intersection.size(); i ++) {
                bus = (Bus*) a->_intersection.at(i);
                for (int j = i + 1; j < a->_intersection.size(); j ++) {
                    arc = chordal->get_arc(bus, a->_intersection.at(j)); // we have to use the arc name.
                    R_lambda_sep[l](arc->_name)  = alpha*R_lambda_in[l](arc->_name).getvalue() + (1 - alpha)*R_lambda_out[l](arc->_name).getvalue();
                    Im_lambda_sep[l](arc->_name) = alpha*Im_lambda_in[l](arc->_name).getvalue() + (1 - alpha)*Im_lambda_out[l](arc->_name).getvalue();
                }

                for (auto g: bus->_gen) {
                    if (g->_active) {
                        mu_sep[l](g->_name) = alpha*mu_out[l](g->_name).getvalue() + (1-alpha)*mu_in[l](g->_name).getvalue();
                        eta_sep[l](g->_name) = alpha*eta_out[l](g->_name).getvalue() + (1 - alpha)*eta_in[l](g->_name).getvalue();
                        kappa_sep[l](g->_name) = alpha*kappa_out[l](g->_name).getvalue() + (1 - alpha)*kappa_in[l](g->_name).getvalue();
                    }
                }
            }
        }

        /////////////////SLOVE SUBPROBLEM/////////////////////
        dual = 0;
        for (int c = 0; c < nb_cliques; c++) {
            val = subproblem(grid, chordal, c, cliquetree, 
                             bag_bus[c], bag_gens[c], bag_arcs[c], bag_bus_pairs[c],
                             R_lambda_sep, Im_lambda_sep,  mu_sep, kappa_sep, eta_sep,
                             rate_ramp, rate_switch, min_up, min_down, cost_up, cost_down,
                             Wii_log[c], R_Wij_log[c], Im_Wij_log[c], Pg_log[c], Qg_log[c], On_off_log[c]);
            value_dual[c] = val;
            dual += val;
        }

// UPDATE POINTS of Kelly using in-out algorithm (Ben-Ameur and Neto)
        if (dual- sum(gamma_sep).eval() < 0) {
            for (int c = 0; c < nb_cliques; c++) {
                Constraint Concavity("Iter_" + to_string(itcount) + "_Concavity_" + to_string(c));
                Concavity += gamma_C(c) - value_dual[c];
                for (auto a: cliquetree->get_node(to_string(c))->get_out()) {
                    for (int i = 0; i < a->_intersection.size(); i ++) {
                        bus = (Bus*)a->_intersection.at(i);
                        for (int j = i + 1; j < a->_intersection.size(); j ++) {
                            arc = chordal->get_arc(bus, a->_intersection.at(j)); // we have to use the arc name.
                            Concavity += (R_lambda_sep[a->_id](arc->_name) - R_lambda_var[a->_id](arc->_name))*R_Wij_log[c](arc->_name);
                            Concavity += (Im_lambda_sep[a->_id](arc->_name) - Im_lambda_var[a->_id](arc->_name))*Im_Wij_log[c](arc->_name);
                        }
                        for (auto g: bus->_gen) {
                            if (g->_active) {
                                Concavity +=(mu_sep[a->_id](g->_name) - mu_var[a->_id](g->_name))*Pg_log[c](g->_name);
                                Concavity +=(kappa_sep[a->_id](g->_name) - kappa_var[a->_id](g->_name))*Qg_log[c](g->_name);
                                Concavity +=(eta_sep[a->_id](g->_name) - eta_var[a->_id](g->_name))*On_off_log[c](g->_name);
                            }
                        }
                    }
                }
            Master.add_constraint(Concavity <= 0);
            }
            if (dual > sum(gamma_in).eval()) {
                for (int c = 0; c < nb_cliques; c++)
                {
                    gamma_in(c) = value_dual[c];
                }
                for (auto a: cliquetree->arcs) {
                    int l = a->_id;
                    for (int i = 0; i < a->_intersection.size(); i ++) {
                        bus = (Bus*) a->_intersection.at(i);
                        for (int j = i + 1; j < a->_intersection.size(); j ++) {
                            arc = chordal->get_arc(bus, a->_intersection.at(j)); // we have to use the arc name.
                            R_lambda_in[l](arc->_name)  = R_lambda_sep[l](arc->_name).getvalue();
                            Im_lambda_in[l](arc->_name)  = Im_lambda_sep[l](arc->_name).getvalue();
                        }

                        for (auto g: bus->_gen) {
                            if (g->_active) {
                                mu_in[l] (g->_name) = mu_sep[l](g->_name).getvalue();
                                eta_in[l] (g->_name) = eta_sep[l](g->_name).getvalue() ;
                                kappa_in[l](g->_name) = kappa_sep[l](g->_name).getvalue();
                            }
                        }
                    }
                }
            }
            solve_Master.run();
            // Update the out-point
            gamma_out = (*(var<Real>*) Master.get_var("gamma_C"));
            for (auto a: cliquetree->arcs) {
                int l = a->_id;
                for (int i = 0; i < a->_intersection.size(); i ++) {
                    bus = (Bus*) a->_intersection.at(i);
                    for (int j = i + 1; j < a->_intersection.size(); j ++) {
                        arc = chordal->get_arc(bus, a->_intersection.at(j)); 
                        R_lambda_out[l](arc->_name)  = R_lambda_sep[l](arc->_name).getvalue();
                        Im_lambda_out[l](arc->_name)  = Im_lambda_sep[l](arc->_name).getvalue();
                    }

                    for (auto g: bus->_gen) {
                        if (g->_active) {
                            mu_out[l] (g->_name) = mu_sep[l](g->_name).getvalue();
                            eta_out[l] (g->_name) = eta_sep[l](g->_name).getvalue() ;
                            kappa_out[l](g->_name) = kappa_sep[l](g->_name).getvalue();
                        }
                    }
                }
            }
        }
        else {
            for (int c = 0; c < nb_cliques; c++)
            {
                gamma_in(c) = value_dual[c];
            }
            for (auto a: cliquetree->arcs) {
                int l = a->_id;
                for (int i = 0; i < a->_intersection.size(); i ++) {
                    bus = (Bus*) a->_intersection.at(i);
                    for (int j = i + 1; j < a->_intersection.size(); j ++) {
                        arc = chordal->get_arc(bus, a->_intersection.at(j)); // we have to use the arc name.
                        R_lambda_in[l](arc->_name)  = R_lambda_sep[l](arc->_name).getvalue();
                        Im_lambda_in[l](arc->_name)  = Im_lambda_sep[l](arc->_name).getvalue();
                    }

                    for (auto g: bus->_gen) {
                        if (g->_active) {
                            mu_in[l] (g->_name) = mu_sep[l](g->_name).getvalue();
                            eta_in[l] (g->_name) = eta_sep[l](g->_name).getvalue() ;
                            kappa_in[l](g->_name) = kappa_sep[l](g->_name).getvalue();
                        }
                    }
                }
            }
        }
        cout<< setw(15) << left <<itcount << setw(15) << LBlog[itcount] << setw(15) << UBlog[itcount] << endl;
    }
    return 0;
}
