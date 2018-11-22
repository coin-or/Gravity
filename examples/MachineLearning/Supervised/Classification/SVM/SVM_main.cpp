//
//  Stable_Set.cpp
//  Gravity
//
//  Created by Hassan Hijazi on Nov 21 2018
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
    auto total_time_start = get_wall_time();
    string solver_str ="ipopt", lazy_str = "no", mu_str = "1e+6", nb_c_str = "200", output_str = "5";
    std::cout << "WELCOME, THIS IS AN IMPLEMENTATION OF A SUPPORT VECTOR MACHINE IN GRAVITY\n";
    
    string fname = string(prj_dir)+"/data_sets/Classification/Archive/svmguide1";
    
    /* Create a OptionParser with options */
    op::OptionParser opt;
    opt.add_option("h", "help",
                   "shows option help"); // no default value means boolean options, which default value is false
    opt.add_option("f", "file", "Input file name", fname);
    opt.add_option("s", "solver", "Solvers: Ipopt/Cplex/Gurobi, default = Ipopt", solver_str);
    opt.add_option("o", "output", "Output level, default = 5", output_str);
    opt.add_option("lz", "lazy", "Generate constraints in a lazy fashion, default = no", lazy_str);
    opt.add_option("mu", "multiplier", "Value of penalization multiplier, default = 1e+6", mu_str);
    opt.add_option("nb", "nb_init_cstr", "Initial number of constraints to add, default = 200 (will add the constraints corresponding to the first 'nb' points in each class). Setting this parameter to -1 will deactivate constraint generation.", nb_c_str);
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
    lazy_str = opt["lz"];
    bool lazy = false;
    if (lazy_str.compare("yes")==0) {
        lazy = true;
    }
    mu_str = opt["mu"];
    double mu = op::str2double(mu_str);
    
    nb_c_str = opt["nb"];
    int nb_c = op::str2int(nb_c_str);
    
    output_str = opt["o"];
    int output = op::str2int(output_str);
    
    /* Reading input data */
    DataSet<> training_set;
    training_set.parse(fname);
    training_set.print_stats();
    DataSet<> test_set;
    test_set.parse(fname+".t", &training_set);
    test_set.print_stats();
    
    /* Defining parameters ans indices */
    auto nf = training_set._nb_features;
    auto m = training_set._nb_points;
    auto Rnf = indices(1,nf); /* Indices 1 to total number of features */
    auto Rm = indices(1,m); /* Indices 1 to total number of training points */
    auto M1 = indices(1,training_set._class_sizes[0]); /* Indices 1 to number of points in first class */
    auto M2 = indices(1,training_set._class_sizes[1]); /* Indices 1 to number of points in second class */
    auto f = training_set.get_features(); /* Feature values */
    
    /* Model */
    Model SVM;
    /* Variables */
    var<> w("w");
    var<> xi("xi", pos_);
    var<> b("b");
    SVM.add(w.in(Rnf), xi.in(Rm), b.in(R(1)));
    /* Objective function */
    SVM.min(product(0.5,power(w,2)) + mu*sum(xi));
    
    /* Constraints */
    size_t nb_c1 = M1.size();
    if (nb_c >= 0 && nb_c <= M1.size()) { /* Activating constraint generation */
        nb_c1 = nb_c;
    }
    for (auto i = 0; i<nb_c1; i++) {
        auto Class1 = Constraint("Class1_"+to_string(i));
        Class1 = product(f[0][i],w) + b + (xi(i) - 1);
        SVM.add(Class1 >= 0);
    }
    /* Remaining constraints are lazy */
    for (auto i = nb_c1; i<M1.size(); i++) {
        auto Class1 = Constraint("Class1_"+to_string(i));
        Class1 = product(f[0][i],w) + b + (xi(i) - 1);
        SVM.add_lazy(Class1 >= 0);
    }
    
    /* Class 2 constraints */
    size_t nb_c2 = M2.size();
    if (nb_c >= 0 && nb_c < M2.size()) { /* Activating constraint generation */
        nb_c2 = nb_c;
    }
    for (auto i = 0; i<nb_c2; i++) {
        auto Class2 = Constraint("Class2_"+to_string(i));
        Class2 = product(f[1][i],w) + b + (1 - xi(i+M1.size()));
        SVM.add(Class2 <= 0);
    }
    /* Remaining constraints are lazy */
    for (auto i = nb_c2; i<M2.size(); i++) {
        auto Class2 = Constraint("Class2_"+to_string(i));
        Class2 = product(f[1][i],w) + b + (1 - xi(i+M1.size()));
        SVM.add_lazy(Class2 <= 0);
    }
    
    /* Start Timers */
    solver SVM_solver(SVM,solv_type);
    double solver_time_start = get_wall_time();
    SVM_solver.run(output);
    double solver_time_end = get_wall_time();
    double total_time_end = get_wall_time();
    auto solve_time = solver_time_end - solver_time_start;
    auto total_time = total_time_end - total_time_start;
    /* Uncomment line below to print expanded model */
    /* SVM.print(); */
    
    /* Testing */
    unsigned err = 0;
    for (auto c = 0; c<test_set._nb_classes; c++) {
        for (auto i = 0; i<test_set._class_sizes[c]; i++) {
            double lhs = 0;
            auto pt = test_set._points[c][i];
            for (auto j = 0; j<nf; j++) {
                lhs += pt._features[j]*w.eval(j);
            }
            lhs += b.eval();
            if (c == 0 && lhs <= 0) {
                err++;
            }
            else if (c == 1 && lhs > 0) {
                err++;
            }
        }
    }
    cout << "Number of misclassified = " << err << " out of " << test_set._nb_points << endl;
    cout << "Accuracy = " << 100.*(test_set._nb_points-err)/test_set._nb_points << "%" << endl;
    
    /** Terminal output */
    cout << "Done running the SVM model\n";
    cout << "Solve wall clock time =  " << solve_time << "\n";
    cout << "Total wall clock time =  " << total_time << "\n";
    return 0;
}




