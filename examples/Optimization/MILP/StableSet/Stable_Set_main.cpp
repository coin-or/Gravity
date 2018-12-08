//
//  Stable_Set.cpp
//  Gravity
//
//  Created by Hijazi, Hassan on 29 Oct 2018.
//
//
#include <iostream>
#include <gravity/solver.h>

using namespace std;
using namespace gravity;

int main (int argc, const char * argv[])
{
    cout << "Welcome, this is an implementation of the Stable Set problem in Gravity" << endl;
    Net graph;
    graph.read_adjacency_matrix(string(prj_dir)+"/data_sets/Stable_set/p.3n150.txt");

    unsigned n = graph.nodes.size();
    Model model;
    /* Declaring the n-dimensional Real space */
    auto Rn = R(n);    
    
    /** Variables **/
    var<double> x("x");
    model.add_var(x.in(Rn));
    /* Or equivalently */
    /* model.add(x.in(graph.nodes)); (graph nodes are indexed in {0,...,n})*/
    
    /** Objective **/
    model.max(sum(x));
    
    /** Constraints **/
    Constraint c("Stable_Set_Const");
    c = x.from() + x.to();
    model.add(c.in(graph.arcs) <= 1);
    
    /** Solver **/
  //solver s(model,cplex);
  solver s(model,Clp);

  s.run();
    
    return 0;
};
