//
//  Sensor_assignment.cpp
//  Gravity
//
//  Created by Hijazi, Hassan on 15 Dec 2021.
//
//
#include <iostream>
#include <gravity/solver.h>

using namespace std;
using namespace gravity;

int main (int argc, const char * argv[])
{
    cout << "Welcome, this is an implementation of the Sensor Assignment problem in Gravity" << endl;
    /* We have n sensors and m space objects, the first n nodes in the graph represent sensors and the remaining ones represent objects */
    int n = 50, m = 10000, degree = 100;
    Net graph;
    bool read_input = false;
    string fname = string(prj_dir)+"/data_sets/sensor/toy.txt";
    if(argc>=2){
        fname = argv[1];
        auto dims =graph.read_pairwise_list(fname);
        n = dims.first;
        m = dims.second;
    }
    else {
        graph.generate_bipartite_random(n,m,degree);
    }
    /* Graph nodes are indexed in {0,...,n+m-1})*/
    indices sensors = range (0,n-1);
    indices objects = range (n,n+m-1);

    assert(n+m==graph.nodes.size());/* Make sure we have the right number of nodes */
    /* Indexing sets */
    indices arcs(graph.arcs);
    param<> w("w");
    w.in(arcs);
    w.initialize_normal(2, 1);
    
    Model<> model("Sensor");
    /* Declaring the n-dimensional Real space */
    
    /** Variables **/
    var<> x("x",0,1);
    model.add(x.in(arcs));
    
    /** Objective **/
    model.max(product(w, x));
    
    /** Constraints **/
    Constraint<> Unique_Obj("Unique_Obj");
    Unique_Obj = x.in_matrix(1, 1);
    model.add(Unique_Obj.in(sensors) == 1);
    
//    Constraint<> Unique_Sensor("Unique_Sensor");
//    Unique_Sensor = x.in_matrix(0, 1);
//    model.add(Unique_Sensor.in(objects) >= 1);
    
//    model.write();
    /** Solver **/
    solver<> s(model,ipopt);
    s.run();
    model.print_solution();
    return 0;
};
