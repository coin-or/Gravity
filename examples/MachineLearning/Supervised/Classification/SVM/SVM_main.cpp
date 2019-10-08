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

/* Builds the static SVM model */
unique_ptr<Model<>> build_svm(const DataSet<>& training_set, double mu){
    
    /* Defining parameters ans indices */
    auto nf = training_set._nb_features;
    auto m1 = training_set._class_sizes[0];
    auto m2 = training_set._class_sizes[1];
    auto F = training_set.get_features_matrices(); /* Features matrices, one per-class */
    
    /* Model */
    auto SVM = unique_ptr<Model<>>(new Model<>());
    
    /* Variables */
    var<> w("w");
    var<> xi1("xi1", pos_);
    var<> xi2("xi2", pos_);
    var<> b("b");
    SVM->add(w.in(R(nf)), xi1.in(R(m1)), xi2.in(R(m2)), b.in(R(1)));
    /* Objective function */
    SVM->min(product(0.5,pow(w,2)) + mu*sum(xi1) + mu*sum(xi2));
    
    /* Constraints */
    /* Class 1 constraints */
    auto Class1 = Constraint<>("Class1");
    Class1 = product(F[0],w) + b + (xi1 - 1);
    SVM->add(Class1 >= 0);
    /* Class 2 constraints */
    auto Class2 = Constraint<>("Class2");
    Class2 = product(F[1],w) + b + (1 - xi2);
    SVM->add(Class2 <= 0);
    
    return SVM;
}


/* Builds the dual SVM model */
unique_ptr<Model<>> build_svm_dual(const DataSet<>& training_set, double mu, const string& kernel_type, double gamma, double r, unsigned d){
    
    /* Defining parameters ans indices */
    auto m = training_set._nb_points;
    auto K = training_set.get_kernel_matrix(kernel_type, gamma, r, d);
    auto y = training_set.get_classes();

    /* Model */
    auto SVM = unique_ptr<Model<>>(new Model<>());
    
    /* Variables */
    var<> alpha("𝛂", 0, mu);
    SVM->add(alpha.in(R(m)));
//    SVM->initialize_uniform();
    /* Objective function */
    SVM->min(0.5*alpha.tr()*K*alpha - sum(alpha));
    
    /* Constraints */
    /* Equality constraints */
    auto Equ0 = Constraint<>("Equation");
    Equ0 = y.tr()*alpha;
    SVM->add(Equ0 == 0);
    
    return SVM;
}


unique_ptr<Model<>> build_lazy_svm(const DataSet<>& training_set, int nb_c, double mu){
    
    /* Defining parameters ans indices */
    auto nf = training_set._nb_features;
    auto m = training_set._nb_points;
    auto m1 = training_set._class_sizes[0]; /* Number of points in first class */
    auto m2 = training_set._class_sizes[1]; /* Number of points in second class */
    auto f = training_set.get_features(); /* Feature values */
    
    /* Model */
    unique_ptr<Model<>> SVM = unique_ptr<Model<>>(new Model<>());
    /* Variables */
    var<> w("w");
    var<> xi("xi", pos_);
    var<> b("b");
    SVM->add(w.in(R(nf)), xi.in(R(m)), b.in(R(1)));
    /* Objective function */
    SVM->min(product(0.5,pow(w,2)) + mu*sum(xi));
    
    /* Constraints */
    size_t nb_c1 = m1;
    if (nb_c >= 0 && nb_c <= m1) { /* Activating constraint generation */
        nb_c1 = nb_c;
    }
    for (auto i = 0; i<nb_c1; i++) {
        auto Class1 = Constraint<>("Class1_"+to_string(i));
        Class1 = product(f[0][i],w) + b + (xi(i) - 1);
        SVM->add(Class1 >= 0);
    }
    /* Remaining constraints are lazy */
    for (auto i = nb_c1; i<m1; i++) {
        auto Class1 = Constraint<>("Class1_"+to_string(i));
        Class1 = product(f[0][i],w) + b + (xi(i) - 1);
        SVM->add_lazy(Class1 >= 0);
    }
    
    /* Class 2 constraints */
    size_t nb_c2 = m2;
    if (nb_c >= 0 && nb_c < m2) { /* Activating constraint generation */
        nb_c2 = nb_c;
    }
    for (auto i = 0; i<nb_c2; i++) {
        auto Class2 = Constraint<>("Class2_"+to_string(i));
        Class2 = product(f[1][i],w) + b + (1 - xi(i+m1));
        SVM->add(Class2 <= 0);
    }
    /* Remaining constraints are lazy */
    for (auto i = nb_c2; i<m2; i++) {
        auto Class2 = Constraint<>("Class2_"+to_string(i));
        Class2 = product(f[1][i],w) + b + (1 - xi(i+m1));
        SVM->add_lazy(Class2 <= 0);
    }
    return SVM;
}


int main (int argc, char * argv[])
{
    double tol = 1e-3, r = 0, gamma = 0, mu;
    unsigned d = 3;
    int nbc, output = 5;
    bool dual = false;
    auto total_time_start = get_wall_time();
    string solver_str ="ipopt", dual_str = "no", lazy_str = "no", mu_str = "1e+6", nb_c_str = "200", output_str = "5", kernel = "linear";
    std::cout << "WELCOME, THIS IS AN IMPLEMENTATION OF A SUPPORT VECTOR MACHINE IN GRAVITY\n";
    
    string fname = string(prj_dir)+"/data_sets/Classification/Archive/svmguide1";
    /* Reading input data */
    DataSet<> training_set;
    training_set.parse(fname);
    training_set.print_stats();
    auto nf = training_set._nb_features;
    gamma = 1./nf;
    SolverType solv_type = ipopt;
#ifdef USE_OPT_PARSER
    /* Create a OptionParser with options */
    op::OptionParser opt;
    opt.add_option("h", "help",
                   "shows option help"); // no default value means boolean options, which default value is false
    opt.add_option("f", "file", "Input file name", fname);
    opt.add_option("s", "solver", "Solvers: Ipopt/Cplex/Gurobi, default = Ipopt", solver_str);
    opt.add_option("d", "dual", "Solve dual SVM, yes/no, default = no", dual_str);
    opt.add_option("k", "kernel", "Choose kernel, linear/poly/rbf/sigm, default = linear", kernel);
    opt.add_option("o", "output", "Output level, default = 5", output_str);
    opt.add_option("lz", "lazy", "Generate constraints in a lazy fashion, yes/no, default = no", lazy_str);
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
    if (solver_str.compare("Gurobi")==0) {
        solv_type = gurobi;
    }
    else if(solver_str.compare("Cplex")==0) {
        solv_type = cplex;
    }else if(solver_str.compare("Mosek")==0) {
        solv_type = _mosek;
    }
    lazy_str = opt["lz"];
    bool lazy = false;
    if (lazy_str.compare("yes")==0) {
        lazy = true;
    }
    dual_str = opt["d"];
    if (dual_str.compare("yes")==0) {
        dual = true;
    }
    mu_str = opt["mu"];
    mu = op::str2double(mu_str);
    
    nb_c_str = opt["nb"];
    nb_c = op::str2int(nb_c_str);
    
    output_str = opt["o"];
    output = op::str2int(output_str);
    
    unique_ptr<Model<>> SVM;
    if(dual){
        kernel = opt["k"];
        SVM = build_svm_dual(training_set, mu, kernel,1./training_set._nb_features,0,3);
    }
    else {
        if (lazy && nb_c>0) {
            SVM = build_lazy_svm(training_set, nb_c, mu);
        }
        else {
            SVM = build_svm(training_set, mu);
        }
    }
#else
    auto SVM = build_svm(training_set, mu);
#endif
    SVM->print_symbolic();
//    SVM->print();
    /* Start Timers */
    solver<> SVM_solver(*SVM,solv_type);
    double solver_time_start = get_wall_time();
    SVM_solver.run(output,tol);
    double solver_time_end = get_wall_time();
    double total_time_end = get_wall_time();
    auto solve_time = solver_time_end - solver_time_start;
    auto total_time = total_time_end - total_time_start;
    /* Uncomment line below to print expanded model */
    /* SVM->print(); */
//    SVM->print_solution(false);
    /* Testing */
    unsigned err = 0;
    DataSet<> test_set;
    test_set.parse(fname+".t", &training_set);
    if(dual){
        double b = 0;
        auto alpha = SVM->get_var<double>("𝛂");
        auto y = training_set.get_classes();
        bool point_on_margin = false;
        for (auto i = 0; i<training_set._nb_points; i++) {
            if (alpha.eval(i)>tol && alpha.eval(i)<mu-tol) {//point lies on the margin
                auto pi = training_set._all_points[i];
                b = y.eval(i);
                for (auto k = 0; k<training_set._nb_points; k++) {
                    auto pk = training_set._all_points[k];
                    b -= alpha.eval(k)*y.eval(k)*training_set.K(*pi,*pk, kernel, gamma, r, d);
                }
                point_on_margin = true;
                break;
            }
        }
        if (!point_on_margin) {
            DebugOn("No points found on margin, aborting" << endl);
            return 0;
        }
        DebugOn("b = " << b << endl);        
        for (auto c = 0; c<test_set._nb_classes; c++) {
            for (auto i = 0; i<test_set._class_sizes[c]; i++) {
                double lhs = 0;
                auto pi = &test_set._points[c][i];
                for (auto k = 0; k<training_set._nb_points; k++) {
                    auto pk = training_set._all_points[k];
                    lhs += alpha.eval(k)*y.eval(k)*training_set.K(*pi,*pk, kernel, gamma, r, d);
                }                
                lhs += b;
                if (c == 0 && lhs <= 0) {
                    err++;
                }
                else if (c == 1 && lhs > 0) {
                    err++;
                }
            }
        }
    }
    else {
        auto w = SVM->get_var<double>("w");
        auto b = SVM->get_var<double>("b");
        DebugOn("b = " << b.eval() << endl);
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

    }
    cout << "Number of misclassified = " << err << " out of " << test_set._nb_points << endl;
    cout << "Accuracy = " << 100.*(test_set._nb_points-err)/test_set._nb_points << "%" << endl;
    
    /** Terminal output */
    cout << "Done running the SVM model\n";
    cout << "Solve wall clock time =  " << solve_time << "\n";
    cout << "Total wall clock time =  " << total_time << "\n";
    return 0;
}




