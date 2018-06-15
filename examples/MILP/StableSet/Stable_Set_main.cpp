//
//  Stable_Set.cpp
//  Gravity
//
//  Created by Hijazi, Hassan on 6/12/17.
//
//
#include <stdio.h>
#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <stdlib.h>
#include <gravity/solver.h>

using namespace std;
using namespace gravity;

int main (int argc, const char * argv[])
{
    //  Start Timers
    std::cout << "WELCOME, THIS IS AN IMPLEMENTATION OF THE STABLE SET PROBLEM AND SOME OF ITS RELAXATIONS IN GRAVITY\n";
    
    Net graph;
    const char* fname = "../../data_sets/stable_set/p.3n150.txt";
    graph.read_adjacency_matrix(fname);
    
    Model model;
    unsigned n = graph.nodes.size();
    
    /**  IP model for the stable set problem. **/
    
    /**  Variables **/
    /* Declaring the n-dimensional Real space */
    auto Rn = R(n);
    var<bool> x("x");
    model.add_var(x.in(Rn));
    /**  Objective **/
    model.max(sum(x));
    
    /**  Constraints **/
    Constraint c("Stable_Set");
    c = x.from() + x.to();
    c <= 1;
    c.print();
    model.add_constraint(c.in(graph.arcs));
    
    /**  Solver **/
    SolverType stype = cplex;
    solver s(model,stype);
    s.run();
    return 0;
};
