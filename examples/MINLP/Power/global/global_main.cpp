//  Decompose.cpp
//
//  Gravity
//
//  Created by Guanglei Wang on 6/9/17.
//
//
#include "Partition.hpp"
#include "global.hpp"
#include <iomanip>

template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 6)
{
    std::ostringstream out;
    out << std::setprecision(n) << a_value;
    return out.str();
}

bool node_id_compare(const Node* n1, const Node* n2) {
    return n1->_id < n2->_id;
}

bool arc_id_compare(const Arc* n1, const Arc* n2) {
    return n1->_id < n2->_id;
}

bool bus_pair_compare(const index_pair* n1, const index_pair* n2) {
    return n1->_name < n2->_name;
}

struct net_param {
    param<Real> c0, c1, c2; /**< Generation costs */
    param<Real> tan_th_min, tan_th_max;
    param<Real> g_ff, g_ft, g_tt, g_tf, b_ff, b_ft, b_tf, b_tt;
    param<Real> S_max;
};

// Note that time based decomposition needs more constraints relaxed..
// for the ANU project, we neglect min-up and min-down constraints in the original ACUC formulation
// while choose to relax ramp up/down constraints and inter temporal constraints (in total 4 constraints).
int main (int argc, const char * argv[])
{
    // Decompose
    const char* fname;
    int Num_part = 1;
    int T= 1;
    if (argc >= 3) {
        fname = argv[1];
        T = atoi(argv[2]);
        Num_part = atoi(argv[3]);
    }
    else {
         //fname = "../../data_sets/Power/nesta_case5_pjm.m";
        //fname = "../../data_sets/Power/nesta_case30_ieee.m";
        //fname = "../../data_sets/Power/nesta_case6_c.m";
        //fname = "../../data_sets/Power/nesta_case5_pjm.m";
        //fname = "../../data_sets/Power/nesta_case3_lmbd.m";
        //fname = "../../data_sets/Power/nesta_case300_ieee.m";
        //fname = "../../data_sets/Power/nesta_case1354_pegase.m";
        //fname = "../../data_sets/Power/nesta_case14_ieee.m";
        fname = "../../data_sets/Power/nesta_case118_ieee.m";
         //fname = "../../data_sets/Power/nesta_case57_ieee.m";
         Num_part = 1;
         T =2;
    }
    auto grid = new PowerNet();
    grid->readgrid(fname);
    grid->get_tree_decomp_bags();
    grid->update_net();

    //grid->c2.print(true);
    auto nb_bus_pairs = grid->get_nb_active_bus_pairs();
    auto nb_gen = grid->get_nb_active_gens();
    auto nb_lines = grid->get_nb_active_arcs();
    auto nb_buses = grid->get_nb_active_nodes();
    double  total_time_end, total_time;
    double total_time_start = get_cpu_time();
    // Schedule Parameters
    bool include_min_updown = true;
    auto global = new Global(grid, Num_part, T);
    double cst_t = global->getdual_relax_time_(include_min_updown);
    double lr_t = global->LR_bound_time_(include_min_updown);
    cout << "time lr lower bound: " << to_string(lr_t) << endl;
    
    double ub = global->Upper_bound_sequence_(include_min_updown);
    cout << "time upper bound is: " << to_string(ub) << endl;

    // upper bound
    //double cst_s =global->getdual_relax_spatial();
    //double lr_s = global->LR_bound_spatial_();
    //cout << "Spaital lr lower bound: " << to_string(lr_s) << endl;
    total_time_end = get_cpu_time();
    total_time = total_time_end - total_time_start;
    string out = grid->_name + ", " + to_string(nb_buses) + ", " + to_string(nb_lines)
                 +", " + to_string(nb_gen) + ", "+ to_string_with_precision(lr_t, 6)
                    + ", " + to_string_with_precision(ub, 6) +"," + to_string_with_precision(total_time, 6);
   cout << out << endl;
   ofstream outfile("ACUC_MISOCP.txt", ios_base::app);
    if (!outfile)
        cerr << "Oops! Uable to save session data! \n";
    else{
        //      outfile << "Instance,  CPU, Value" << endl;
        outfile << out << endl;
    }

    return 0;
}
