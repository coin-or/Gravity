//
//  Decompose.cpp
//  Gravity
//
//  Created by Guanglei Wang on 6/9/17.
//
//

#include <stdio.h>
#include <iostream>
#include <string>
#include <stdio.h>
#include <cstring>
#include <fstream>
#include "../PowerNet.h"
#include <gravity/model.h>
#include <gravity/solver.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef USE_BOOST
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <deque>
#include <iterator>
#endif
using namespace std;
using namespace gravity;

Net* get_cliquetree(Net* grid){
#ifdef USE_BOOST
    /** Note that we also need the edge information of the clique tree **/
    /** boost graph library or implement the expanded version of MCS algorithm by Blair and Peyton */
    typedef boost::adjacency_list <boost::vecS,
            boost::vecS,
            boost::undirectedS,
            boost::no_property,
            boost::property < boost::edge_weight_t, int >
            > Graph;
    typedef boost::graph_traits <Graph>::edge_descriptor Edge;
    typedef boost::graph_traits <Graph>::vertex_descriptor Vertex;

    // BUILD THE INTERSECTION GRAPH OF THE CLIQUES
    typedef std::pair<int, int> E;
    std::vector<E> edges;
    std::vector<int> weights;
    int nb_cliques = grid->_bags.size();
    for (int i = 0; i < nb_cliques; i++) {
        cout << "bag " << i << " has " << grid->_bags[i].size() << " vertices." <<endl;
        sort(grid->_bags[i].begin(), grid->_bags[i].end());
        for (int j = i +1; j < nb_cliques; j++) {
            vector<Node*> v3;
            sort(grid->_bags[j].begin(), grid->_bags[j].end());
            set_intersection(grid->_bags[i].begin(), grid->_bags[i].end(), grid->_bags[j].begin(),
                             grid->_bags[j].end(), back_inserter(v3));
            if (v3.size() > 0) {
                edges.push_back(E(i, j));
                weights.push_back(-v3.size());
            }
        }
    }
    size_t num_edges = edges.size();

#if defined(BOOST_MSVC) && BOOST_MSVC <= 1300
    Graph g(num_nodes);
    boost::property_map<Graph, edge_weight_t>::type weightmap = get(edge_weight, g);
    for (std::size_t j = 0; j < num_edges; ++j) {
        Edge e;
        bool inserted;
        boost::tie(e, inserted) = boost::add_edge(edges[j].first, edges[j].second, g);
        boost::weightmap[e] = weights[j];
    }
#else
    Graph g(edges.begin(), edges.end(), weights.begin(), nb_cliques);
#endif
    boost::property_map < Graph, boost::edge_weight_t >::type weight = get(boost::edge_weight, g);
    std::vector < Edge > spanning_tree;
    boost::kruskal_minimum_spanning_tree(g, std::back_inserter(spanning_tree));

    DebugOn("Print the total " << spanning_tree.size() << " edges in the clique tree:" << endl);

    //////////CLIQUE TREE /////////////////////////////
    Net* cliquetree = new Net();
    Node* node = nullptr;
    Arc*  a = nullptr;
    string name;
    for (int i = 0; i < nb_cliques; i++) {
        node= new Node(to_string(i), i);
        cliquetree->add_node(node);
    }

    for (std::vector < Edge >::iterator ei = spanning_tree.begin();
            ei != spanning_tree.end(); ++ei) {
        int u = source(*ei, g);
        int v = target(*ei, g);
        DebugOn(u << " <--> " << v
                << " with weight of " << -weight[*ei]
                << endl);
        name = (int) cliquetree->arcs.size();
        a = new Arc(name);
        a->_id = cliquetree->arcs.size();

        // intersection
        vector<Node*> v3;
        sort(grid->_bags[u].begin(), grid->_bags[u].end());
        sort(grid->_bags[v].begin(), grid->_bags[v].end());
        set_intersection(grid->_bags[u].begin(), grid->_bags[u].end(),
                         grid->_bags[v].begin(), grid->_bags[v].end(),
                         back_inserter(v3));

        a->_src = cliquetree->get_node(to_string(u));
        a->_dest = cliquetree->get_node(to_string(v));
        a->_weight = -weight[*ei];
        a->_intersection = v3;
        cliquetree->add_arc(a);
        a->connect();
    }
    return cliquetree;
#endif
}

/** INITIALISE SUBPROBLEM MODEL */
// it returns a outer-approximation function object
double  subproblem(PowerNet* grid,Net* chordal, unsigned T, unsigned c, Net* cliquetree,
                   vector<Bus*> bag_bus, vector<Gen*> bag_gens, vector<Arc*> bag_arcs,
                   vector<param<Real>>& R_lambda_sep, vector<param<Real>>& Im_lambda_sep, vector<param<Real>>& mu_sep,
                   vector<param<Real>>& kappa_sep, vector<param<Real>>& eta_sep,
                   param<Real>& rate_ramp, param<Real>& rate_switch, param<Real>& min_up, param<Real>& min_down,
                   param<Real>& cost_up, param<Real>& cost_down,
                   param<Real>& Wii_log, param<Real>& R_Wij_log, param<Real>& Im_Wij_log,
                   param<Real>& Pg_log, param<Real>& Qg_log, param<Real>& On_off_log )
{
//    func_  OA;
    cout << "Solving subproblem associated with maximal clique .........." << c << endl;
    if (bag_arcs.size() == 0) {
        //OA += 0;
        //return OA;
        return 0;
    }
    Model Subr("Subr");
    // POWER FLOW
    DebugOn("bag_arcs " << c << " has " << bag_arcs.size() << " lines." << endl);
    var<Real> Pf_from("Pf_from", grid->S_max.in(bag_arcs, T));
    var<Real> Qf_from("Qf_from", grid->S_max.in(bag_arcs, T));
    var<Real> Pf_to("Pf_to", grid->S_max.in(bag_arcs, T));
    var<Real> Qf_to("Qf_to", grid->S_max.in(bag_arcs, T));
    Subr.add_var(Pf_from^(T*bag_arcs.size()));
    Subr.add_var(Qf_from^(T*bag_arcs.size()));
    Subr.add_var(Pf_to^(T*bag_arcs.size()));
    Subr.add_var(Qf_to^(T*bag_arcs.size()));

    // Lifted variables.
    //var<Real>  R_Wij("R_Wij", grid->wr_min.in(bus_pairs), grid->wr_max.in(bus_pairs)); // real part of Wij
    //var<Real>  Im_Wij("Im_Wij", grid->wi_min.in(bus_pairs), grid->wi_max.in(bus_pairs)); // imaginary part of Wij.
    var<Real>  R_Wij("R_Wij");
    var<Real>  Im_Wij("Im_Wij");
    var<Real>  Wii("Wii", grid->w_min.in(bag_bus, T), grid->w_max.in(bag_bus, T));
    Subr.add_var(Wii^(T*bag_bus.size()));
    Subr.add_var(R_Wij^(T*grid->_bags.size()*(bag_bus.size() - 1)/2));
    Subr.add_var(Im_Wij^(T*grid->_bags.size()*(bag_bus.size() - 1)/2));

    // Commitment variables
    var<Real>  On_off("On_off", 0, 1);
    var<Real>  Start_up("Start_up", 0, 1);
    var<Real>  Shut_down("Shut_down", 0, 1);
    Subr.add_var(On_off^(T*bag_gens.size()));
    Subr.add_var(Start_up^(T*bag_gens.size()));
    Subr.add_var(Shut_down^(T*bag_gens.size()));

    /* Construct the objective function*/
    func_ obj;
    Node* Cr = cliquetree->get_node(to_string(c));
    Arc* arc = nullptr;
    Bus* bus = nullptr;
    for (auto a: Cr->get_out()) {
        Debug("a->_intersection.size " << a->_intersection.size() << endl);
        for (int i = 0; i < a->_intersection.size(); i ++) {
            bus = (Bus*)a->_intersection.at(i);
            for (int j = i + 1; j < a->_intersection.size(); j ++) {
                arc = chordal->get_arc(bus, a->_intersection.at(j)); // we have to use the arc name.
                obj += R_lambda_sep[a->_id](arc->_name)*R_Wij(arc->_name);
                obj += Im_lambda_sep[a->_id](arc->_name)*Im_Wij(arc->_name);
            }
        }
    }

    for (auto a: Cr->get_in()) {
        for (int i = 0; i < a->_intersection.size(); i ++) {
            bus = (Bus*) a->_intersection.at(i);
            for (int j = i + 1; j < a->_intersection.size(); j ++) {
                arc = chordal->get_arc(bus, a->_intersection.at(j));
                obj -= R_lambda_sep[a->_id](arc->_name)*R_Wij(arc->_name);
                obj -= Im_lambda_sep[a->_id](arc->_name)*Im_Wij(arc->_name);
            }
        }
    }
    var<Real> Pg("Pg", grid->pg_min.in(bag_gens, T), grid->pg_max.in(bag_gens, T));
    var<Real> Qg("Qg", grid->qg_min.in(bag_gens, T), grid->qg_max.in(bag_gens, T));
    Subr.add_var(Pg^(T*bag_gens.size()));
    Subr.add_var(Qg^(T*bag_gens.size()));
    // power generation
    if (bag_gens.size() > 0) {
        for (auto g: bag_gens) {
            DebugOn("number of generators: " << bag_gens.size() << endl);
            if (g->_active) {
                DebugOn("generator name: " << g->_name << endl);
                for (int t = 0; t < T; t++) {
                    if (t > 1) {
                        string l = to_string(t);
                        obj += grid->c1(g->_name, l)*Pg(g->_name, l) + grid->c2(g->_name, l)*Pg(g->_name, l)*Pg(g->_name, l) + grid->c0(g->_name, l);
                        //obj += cost_up.getvalue()*Start_up(g->_name, l)+ cost_down.getvalue()*Shut_down(g->_name, l);
                    }
                    else {
                        /* This is weird, Pg should either be indexed by two indices or one, not both.*/
                        // you're right. 
                        obj += grid->c1(g->_name)*Pg(g->_name) + grid->c2(g->_name)*Pg(g->_name)*Pg(g->_name) + grid->c0(g->_name);
                        //obj += cost_up.getvalue()*Start_up(g->_name)+ cost_down.getvalue()*Shut_down(g->_name);
                    }
                }
            }
            obj.print(true);
        }
        for (auto a: Cr->get_out()) {
            for (int i = 0; i < a->_intersection.size(); i ++) {
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
            for (int i = 0; i < a->_intersection.size(); i ++) {
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
        //KCL
        for (int t = 0; t < T; t++) {
            for (auto bus: bag_bus) {
                if (!bus->_active) {
                    continue;
                }
                Constraint KCL_P("KCL_P"+bus->_name+ "time_" + to_string(t));
                Constraint KCL_Q("KCL_Q"+bus->_name+ "time_" + to_string(t));

                /* Power Conservation */
                vector<Arc*> out;
                /* Add comments here, why are we using set_intersection?*/
                /* bus->get_out() includes all out lines of the grid, but here we need the out lines of the bag.*/
                set_intersection(bus->get_out().begin(), bus->get_out().end(), bag_arcs.begin(), bag_arcs.end(), back_inserter(out));
                vector<Arc*> in;
                set_intersection(bus->get_in().begin(), bus->get_in().end(), bag_arcs.begin(), bag_arcs.end(), back_inserter(in));
                /* Pg(*,0) is not defined, in the objective, for t=0 you're only using Pg(*) */
                KCL_P  = sum(Pf_from.in_at(out, t))+ sum(Pf_to.in_at(in, t)) + bus->pl() - sum(Pg.in_at(bus->_gen, t));
                KCL_Q  = sum(Qf_from.in_at(out, t)) + sum(Qf_to.in_at(in, t))+ bus->ql() - sum(Qg.in_at(bus->_gen, t));

                /* Shunts */
                KCL_P +=  bus->gs()*Wii(bus->_name);
                KCL_Q -=  bus->bs()*Wii(bus->_name);

                Subr.add_constraint(KCL_P = 0);
                Subr.add_constraint(KCL_Q = 0);
            }
        }
        /* Commitment constraints */
        // Inter-temporal constraints
        for (int t = 1; t < T; t++) {
            Constraint MC1("MC1_"+ to_string(t));
            Constraint MC2("MC2_"+ to_string(t));
            MC1 = On_off.in_at(bag_gens, t) - On_off.in_at(bag_gens, t-1)-  Start_up.in_at(bag_gens, t);
            MC2 = On_off.in_at(bag_gens, t-1) - On_off.in_at(bag_gens, t) - Shut_down.in_at(bag_gens, t);
            Subr.add_constraint(MC1 <= 0);
            Subr.add_constraint(MC2 <= 0);
        }

        // Min-up constraints
        for (int t = 1; t < T; t++) {
            Constraint Min_up1("Min_up1_"+ to_string(t));
            Min_up1 = On_off.in_at(bag_gens, t) - On_off.in_at(bag_gens, t-1) - Start_up.in_at(bag_gens, t) + Shut_down.in_at(bag_gens, t);
            Subr.add_constraint(Min_up1 = 0);
        }

        for (int t = min_up.getvalue(); t < T; t++) {
            Constraint Min_Up("Min_Up_constraint" + to_string(t));
            for (int l = t-min_up.getvalue()+1; l < t +1; l++) {
                Min_Up   += Start_up.in_at(bag_gens, l);
            }
            Min_Up -= On_off.in_at(bag_gens, t);
            Subr.add_constraint(Min_Up <= 0);
        }

        for (int t = min_down.getvalue(); t < T; t++) {
            Constraint Min_Down("Min_Down_constraint" + to_string(t));
            for (int l = t-min_down.getvalue()+1; l < t +1; l++) {
                Min_Down   += Shut_down.in_at(bag_gens, l);
            }
            Min_Down -= 1 - On_off.in_at(bag_gens, t);
            Subr.add_constraint(Min_Down <= 0);
        }

        //Ramp rate
        Constraint Production_P_LB("Production_P_LB");
        Constraint Production_P_UB("Production_P_UB");
        Constraint Production_Q_LB("Production_Q_LB");
        Constraint Production_Q_UB("Production_Q_UB");

        Production_P_UB = Pg.in(bag_gens, T) - grid->pg_max.in(bag_gens, T)*On_off.in(bag_gens,T);
        Production_P_LB = Pg.in(bag_gens, T) - grid->pg_min.in(bag_gens, T)*On_off.in(bag_gens,T);
        Subr.add_constraint(Production_P_UB <=0);
        Subr.add_constraint(Production_P_LB >= 0);

        //grid->qg_max.print(true);
        //grid->qg_min.print(true);

        Production_Q_UB = Qg.in(bag_gens, T) - grid->qg_max.in(bag_gens, T)*On_off.in(bag_gens,T);
        Production_Q_LB = Qg.in(bag_gens, T) - grid->qg_min.in(bag_gens, T)*On_off.in(bag_gens,T);
        Subr.add_constraint(Production_Q_UB <= 0);
        Subr.add_constraint(Production_Q_LB >= 0);

        for (int t = 1; t < T; t++) {
            Constraint Ramp_up("Ramp_up_constraint" + to_string(t));
            Constraint Ramp_down("Ramp_down_constraint" + to_string(t));

            Ramp_up = Pg.in_at(bag_gens, t);
            Ramp_up -= Pg.in_at(bag_gens, t-1);
            Ramp_up -= rate_ramp*On_off.in_at(bag_gens, t-1);
            Ramp_up -= rate_switch*(1 - On_off.in_at(bag_gens, t));

            Ramp_down = Pg.in_at(bag_gens, t-1);
            Ramp_down -= Pg.in_at(bag_gens, t);
            Ramp_down -= rate_ramp*On_off.in_at(bag_gens, t);
            Ramp_down -= rate_switch*(1 - On_off.in_at(bag_gens, t-1));

            Subr.add_constraint(Ramp_up <= 0);
            Subr.add_constraint(Ramp_down <= 0);
        }
    }
    else {
        Subr.set_objective(min(obj));
        //KCL
        for (int t = 0; t < T; t++)
            for (auto bus: bag_bus) {
                if (!bus->_active) {
                    continue;
                }
                Constraint KCL_P("KCL_P"+bus->_name+ "time_" + to_string(t));
                Constraint KCL_Q("KCL_Q"+bus->_name+ "time_" + to_string(t));

                /* Power Conservation */
                vector<Arc*> out;
                set_intersection(bus->get_out().begin(), bus->get_out().end(), bag_arcs.begin(), bag_arcs.end(), back_inserter(out));
                vector<Arc*> in;
                set_intersection(bus->get_in().begin(), bus->get_in().end(), bag_arcs.begin(), bag_arcs.end(), back_inserter(in));
                KCL_P  = sum(Pf_from.in_at(out, t))+ sum(Pf_to.in_at(in, t)) + bus->pl();
                KCL_Q  = sum(Qf_from.in_at(out, t)) + sum(Qf_to.in_at(in, t))+ bus->ql();

                /* Shunts */
                KCL_P +=  bus->gs()*Wii(bus->_name);
                KCL_Q -=  bus->bs()*Wii(bus->_name);

                Subr.add_constraint(KCL_P = 0);
                Subr.add_constraint(KCL_Q = 0);
            }
    }
    /* SOCP constraints */
    Constraint SOC("SOC");
    SOC =  power(R_Wij.in(bag_arcs, T), 2) + power(Im_Wij.in(bag_arcs, T), 2) - Wii.from(bag_arcs, T)*Wii.to(bag_arcs, T) ;
    Subr.add_constraint(SOC <= 0);

    //AC Power Flow.
    Constraint Flow_P_From("Flow_P_From");
    Flow_P_From += Pf_from.in(bag_arcs, T);
    Flow_P_From -= grid->g_ff.in(bag_arcs, T)*Wii.from(bag_arcs, T);
    Flow_P_From -= grid->g_ft.in(bag_arcs, T)*R_Wij.in(bag_arcs, T);
    Flow_P_From -= grid->b_ft.in(bag_arcs, T)*Im_Wij.in(bag_arcs, T);
    Flow_P_From = 0;
    Subr.add_constraint(Flow_P_From);

    Constraint Flow_P_To("Flow_P_To");
    Flow_P_To += Pf_to.in(bag_arcs, T);
    Flow_P_To -= grid->g_tt.in(bag_arcs, T)*Wii.to(bag_arcs, T);
    Flow_P_To -= grid->g_tf.in(bag_arcs, T)*R_Wij.in(bag_arcs, T);
    Flow_P_To += grid->b_tf.in(bag_arcs, T)*Im_Wij.in(bag_arcs, T);
    Flow_P_To = 0;
    Subr.add_constraint(Flow_P_To);

    Constraint Flow_Q_From("Flow_Q_From");
    Flow_Q_From += Qf_from.in(bag_arcs, T);
    Flow_Q_From += grid->b_ff.in(bag_arcs, T)*Wii.from(bag_arcs, T);
    Flow_Q_From += grid->b_ft.in(bag_arcs, T)*R_Wij.in(bag_arcs, T);
    Flow_Q_From += grid->g_ft.in(bag_arcs, T)*Im_Wij.in(bag_arcs, T);
    Flow_Q_From = 0;
    Subr.add_constraint(Flow_Q_From);

    Constraint Flow_Q_To("Flow_Q_To");
    Flow_Q_To += Qf_to.in(bag_arcs, T);
    Flow_Q_To += grid->b_tt.in(bag_arcs, T)*Wii.to(bag_arcs, T);
    Flow_Q_To += grid->b_tf.in(bag_arcs, T)*R_Wij.in(bag_arcs, T);
    Flow_Q_To -= grid->g_tf.in(bag_arcs, T)*Im_Wij.in(bag_arcs, T);
    Flow_Q_To = 0;
    Subr.add_constraint(Flow_Q_To);

    //Phase Angle Bounds constraints */
    //NOTE THAT WE SHOULD USE BUS PAIRS!!!
    Constraint PAD_UB("PAD_UB");
    PAD_UB = Im_Wij.in(bag_arcs, T);
    PAD_UB -= (grid->tan_th_max).in(bag_arcs, T)*R_Wij.in(bag_arcs, T);
    Subr.add_constraint(PAD_UB <= 0);

    Constraint PAD_LB("PAD_LB");
    PAD_LB =  Im_Wij.in(bag_arcs, T);
    PAD_LB -= grid->tan_th_min.in(bag_arcs, T)*R_Wij.in(bag_arcs, T);
    Subr.add_constraint(PAD_LB >= 0);

    //Thermal Limit Constraints */
    Constraint Thermal_Limit_from("Thermal_Limit_from");
    Thermal_Limit_from += power(Pf_from.in(bag_arcs, T),  2) + power(Qf_from.in(bag_arcs, T), 2);
    Thermal_Limit_from -= power(grid->S_max.in(bag_arcs, T), 2);
    Subr.add_constraint(Thermal_Limit_from <= 0);

    Constraint Thermal_Limit_to("Thermal_Limit_to");
    Thermal_Limit_to += power(Pf_to.in(bag_arcs, T), 2) + power(Qf_to.in(bag_arcs, T), 2);
    Thermal_Limit_to -= power(grid->S_max.in(bag_arcs, T),2);
    Subr.add_constraint(Thermal_Limit_to <= 0);

    /* Resolve it! */
    solver solve_Subr(Subr,ipopt);
    solve_Subr.run();
    // OA = obj.get_outer_app();
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
    //fname = "../../data_sets/Power/nesta_case5_pjm.m";
    //fname = "../../data_sets/Power/nesta_case2383wp_mp.m";
    fname = "../../data_sets/Power/nesta_case14_ieee.m";
    //fname = "../../data_sets/Power/nesta_case1354_pegase.m";

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
        rate_ramp._dim++;
        rate_switch(g->_name) = max(grid->pg_min(g->_name).getvalue(), 0.25*grid->pg_max(g->_name).getvalue());
        rate_switch._dim++;
    }
    min_up = 1;
    min_down = 1;
    cost_up = 50;
    cost_down = 30;

    /* create a function grid->time_expand(T) to do all these below */
    grid->time_expand(T);
    rate_ramp.time_expand(T);
    rate_switch.time_expand(T);

    /** Clique tree decomposition **/
    Net* chordal = grid->get_chordal_extension();
    grid->get_clique_tree();
    const unsigned nb_cliques = grid->_bags.size();

    vector<vector<Bus*>> bag_bus; // Note that each clique contains just nodes, not buses! Fixed this by modifying the bag definition.
    vector<vector<Gen*>> bag_gens;
    vector<vector<Arc*>> bag_arcs; //bag_arcs contains the arcs of the power grid while variables associated with W are defined on chordal graph.

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
        DebugOn("bag " << c << " has " << bag_gens[c].size() << " generators. " << endl;)
        DebugOn("bag " << c << " has " << bag_arcs[c].size() << " line " << endl;)
    }

    bag_bus.resize(nb_cliques);
    bag_arcs.resize(nb_cliques);
    bag_gens.resize(nb_cliques);

    Net* cliquetree = get_cliquetree(grid);
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
        int l = a->_id;
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
        val = subproblem(grid, chordal, T, c, cliquetree,
                         bag_bus[c], bag_gens[c], bag_arcs[c],
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

    solver solve_Master(Master, ipopt);
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
            val = subproblem(grid, chordal, T, c, cliquetree,
                             bag_bus[c], bag_gens[c], bag_arcs[c],
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
                        arc = chordal->get_arc(bus, a->_intersection.at(j)); // we have to use the arc name.
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
