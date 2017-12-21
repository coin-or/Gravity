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
#include <stdio.h>
#include <stdlib.h>

using namespace std;
using namespace gravity;


int main (int argc, const char * argv[])
{
    //  Start Timers
    std::cout << "WELCOME, THIS IS AN IMPLEMENTATION OF THE STABLE SET PROBLEM AND SOME OF ITS RELAXATIONS IN GRAVITY\n";
    
    unsigned i, j;
    
    Net graph;
//    const char* fname = "../../data_sets/stable_set/toy.txt";
    const char* fname = "../../data_sets/stable_set/p.3n150.txt";
    graph.read_adjacency_matrix(fname);
    
    Model model;
    unsigned n = graph.nodes.size();
    
    /**  IP model for the stable set problem. **/
    
    /**  Variables **/
    var<bool> x("x");
//    var<> x("x", 0,1);
    model.add_var(x^n);
    
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
    double wall0 = get_wall_time();
    double cpu0  = get_cpu_time();
    cout << "Running the IP model\n";
    s.run();
    double wall1 = get_wall_time();
    double cpu1  = get_cpu_time();
    cout << "Done running the IP model\n";
    cout << "Wall clock computing time =  " << wall1 - wall0 << "\n";
    cout << "CPU computing time =  " << cpu1 - cpu0 << "\n";

    /* Schriver's SDP relaxation for the stable set problem */
    Model SDP;
    /* Variable declaration */
    var<double> Xii("Xii", 0, 1);
    var<double> Xij("Xij", 0, 1);
    SDP.add_var(Xii^n); /*< Diagonal entries of the matrix */
    SDP.add_var(Xij^(n*(n-1)/2)); /*< Lower left triangular part of the matrix excluding the diagonal*/
    
    /* Constraints declaration */
    ordered_pairs indices(1, n);
    Constraint SOCP("SOCP");
    SOCP =  power(Xij, 2) - Xii.from()*Xii.to() ;
    SDP.add_constraint(SOCP.in(indices._keys) <= 0);
    Constraint diag("diag");
    diag = sum(Xii);
    SDP.add_constraint(diag = 1); // diagonal sum is 1
    
    Constraint zeros("zeros");
    zeros = Xij.in(graph.arcs);
    SDP.add_constraint(zeros = 0); // zero elements
    
    auto Xij_ = Xij.pairs_in(graph._bags, 3);
    auto Xii_ = Xii.in(graph._bags, 3);
    Constraint SDP3("SDP_3D");
    SDP3 = -2*Xij_[0]*Xij_[1]*Xij_[2];
    SDP3 -= Xii_[0]*Xii_[1]*Xii_[2];
    SDP3 += power(Xij_[0],2)*Xii_[2];
    SDP3 += power(Xij_[2],2)*Xii_[1];
    SDP3 += power(Xij_[1],2)*Xii_[0];
//    SDP3.print();
//    SDP.add_constraint(SDP3);
//    set<tuple<int,int,int>> ids;
//        for (i = 0; i < complement_graph._bags.size(); i++){
////        for (i = 0; i < 1; i++){
//            auto bag = complement_graph._bags.at(i);
//            if (bag.size() < 3) {
//                continue;
//            }
////            for (int j = 0; j < bag.size()-2; j++){
//            for (int j = 0; j < 1; j++){
//                i1 = bag[j]->ID;
//                i2 = bag[j+1]->ID;
//                i3 = bag[j+2]->ID;
//                assert(i2>i1 && i3>i2);
//                if(ids.count(make_tuple(i1, i2, i3))==0){
//                    ids.insert(make_tuple(i1, i2, i3));
//                }
//                else {
//                    continue;
//                }
//                Constraint SDP3("SDP3("+to_string(i1)+","+to_string(i2)+","+to_string(i3)+")");
//                SDP3 = -2*Xij(i1,i2)*Xij(i2,i3)*Xij(i1,i3);
//                SDP3 -= Xii(i1)*Xii(i2)*Xii(i3);
//                SDP3 += power(Xij(i1,i2),2)*Xii(i3);
//                SDP3 += power(Xij(i1,i3),2)*Xii(i2);
//                SDP3 += power(Xij(i2,i3),2)*Xii(i1);
//                SDP3.print();
//                SDP.add_constraint(SDP3);
//            }
//        }
    
    /* Objective declaration */
    SDP.max(2*sum(Xij) + sum(Xii));
    

//   solver s1(SDP,ipopt);
    solver s1(SDP,cplex);

    wall0 = get_wall_time();
    cpu0  = get_cpu_time();
    cout << "Running the SDP relaxation\n";
    s1.run();
    wall1 = get_wall_time();
    cpu1  = get_cpu_time();
    cout << "Done running the SDP relaxation\n";
    cout << "\nWall clock computing time =  " << wall1 - wall0 << "\n";
    cout << "CPU computing time =  " << cpu1 - cpu0 << "\n";
    return 0;
    
    /* Outer-Approximation of SOCP and 3d-SDP cuts in Shriver's SDP relaxation for the stable set problem */
    Model OA;
    OA.add_var(Xii);
    OA.add_var(Xij);
    for (auto &cs_p: SDP._cons) {
        if (!cs_p.second->is_linear() && cs_p.second->is_active()){ //Active nonlinear constraint
            DebugOff("Active constraint:" << cs_p.second->to_str() << endl);
            Constraint oa("OA_"+cs_p.second->get_name());
            oa = cs_p.second->get_outer_app();
            //            oa.print();
            OA.add_constraint(oa);
        }
    }
//    OA.add_constraint(diag=1);
    for(auto a: graph.arcs){
        i = (a->_src)->_id;
        j = (a->_dest)->_id;
        Constraint zeros("zeros("+to_string(i)+","+to_string(j)+")");
        zeros = Xij(i,j);
        OA.add_constraint(zeros=0);
    }
    OA.max(2*sum(Xij) + sum(Xii));
    solver s2(OA,cplex);
    wall0 = get_wall_time();
    cpu0  = get_cpu_time();
    cout << "Running the OA-SDP relaxation\n";
    s2.run();
    wall1 = get_wall_time();
    cpu1  = get_cpu_time();
    cout << "Done running the OA-SDP relaxation\n";
    cout << "\nWall clock computing time =  " << wall1 - wall0 << "\n";
    cout << "CPU computing time =  " << cpu1 - cpu0 << "\n";
    return 0;
};
