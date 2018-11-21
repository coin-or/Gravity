//
//  Stable_Set.cpp
//  Gravity
//
//  Created by Hijazi, Hassan (Data61, Canberra City) on 6/12/17.
//
//


#include <stdio.h>
#include <iostream>
#include <string>
#include <stdio.h>
#include <cstring>
#include <fstream>
#include <gravity/Net.h>
#include <gravity/model.h>
#include <gravity/solver.h>
#include <DataSet.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;
using namespace gravity;


int main (int argc, char * argv[])
{
    string solver_str ="ipopt", lazy = "no";
    std::cout << "WELCOME, THIS IS AN IMPLEMENTATION OF A SUPPORT VECTOR MACHINE IN GRAVITY\n";

    string fname = string(prj_dir)+"/data_sets/classification/Archive/vowel";
    
    /* Create a OptionParser with options */
    op::OptionParser opt;
    opt.add_option("h", "help",
                   "shows option help"); // no default value means boolean options, which default value is false
    opt.add_option("f", "file", "Input file name", fname);
    opt.add_option("s", "solver", "Solvers: ipopt/cplex/gurobi, default = ipopt", solver_str);
    opt.add_option("lz", "lazy", "Generate constraints in a lazy fashion, default = no", lazy);
    /* Parse the options and verify that all went well. If not, errors and help will be shown */
    bool correct_parsing = opt.parse_options(argc, argv);
    
    if (!correct_parsing) {
        return EXIT_FAILURE;
    }
    
    fname = opt["f"];
    bool has_help = op::str2bool(opt["h"]);
    if (has_help) {
        opt.show_help();
        exit(0);
    }
    solver_str = opt["s"];
    SolverType solv_type = ipopt;
    if (solver_str.compare("Gurobi")==0) {
        solv_type = gurobi;
    }
    else if(solver_str.compare("Cplex")==0) {
        solv_type = cplex;
    }else if(solver_str.compare("Mosek")==0) {
        solv_type = Mosek;
    }
    lazy = opt["lz"];
    bool lazy_bool = false;
    if (lazy.compare("yes")==0) {
        lazy_bool = true;
    }
    
    /* Reading input data */
    DataSet<> training_set;
    training_set.parse(fname);
    training_set.print_stats();
    auto nf = training_set._nb_features;
    auto m = training_set._nb_points;
    auto F = indices(1,nf);
    auto M = indices(1,m);
    auto point = training_set._points[0][0];
    param<> f("f"); /* Feature values */
    /* Model */
    Model SVM;
    /* Variables */
    var<> w("w");
    var<> xi("xi", pos_);
    var<> b("b");
    SVM.add(w.in(F), xi.in(M), b);
    /* Objective function */
    auto Obj= product(0.5,power(w,2)) + 1e+6*sum(xi);
    
    /* Constraints */
    auto InClass = Constraint("InClass");
    InClass = product(f,w) + b + (xi - 1);
    SVM.add(InClass.in(M) >= 0);
    
    auto OutClass = Constraint("OutClass");
    OutClass = product(f,w) + b + (xi - 1);
    SVM.add(OutClass.in(M) <= 0);
    
    /* Start Timers */
    double wall0 = get_wall_time();
    double cpu0  = get_cpu_time();
    double wall1 = get_wall_time();
    double cpu1  = get_cpu_time();
    cout << "Done running the SVM model\n";
    cout << "Wall clock computing time =  " << wall1 - wall0 << "\n";
    cout << "CPU computing time =  " << cpu1 - cpu0 << "\n";
    return 0;
}
