//
//  model.hpp
//  bilevel_sensor
//
//  Created by Svetlana Riabova on 6/10/22.
//

#pragma once
#include <stdio.h>
#include <random>
#include <gravity/solver.h>

using namespace std;
using namespace gravity;

/*struct Agent {
    int n;
    vector<int> own;
    vector<int> oths;
    indices own_arcs;
    indices oths_arcs;
    vector<double> w;
    
    Agent() { };
    ~Agent() { };
};*/

class myModel {
public:
    Model<> model;
    indices sensors;
    indices objects;
    int M;
    int N;
    int K = 5;
    double e;
    double max_price;
    Net graph;
    indices arcs;
    indices own_arcs;
    indices bought_arcs;
    indices own_sens;
    indices bought_sens;
    indices own_rplc;
    indices oths_rplc;
    indices own_oths_rplc1;
    indices own_oths_rplc2;
    indices jk; //agents X objects
    vector<int> owner;
    
public:
    myModel(){ };
    vector<param<double>> readData(int argc, const char * argv[], int n1, int n2);
    void InitBilevel(param<double> &w0, param<double> &w_own, param<double> &w_bought, double eps);
    void mSolve();
    void saveSolStats();
    void GreedyStart(const param<double> &w0, const param<double> &w_own, const param<double> &w_bought);
    void assignLeader(string &idx, param<double> &wt0, param<double> &wt_own, param<double> &wt_bought);
    void assignOwn(string &idx, param<double> &wt0, param<double> &wt_own, param<double> &wt_bought);
    void assignBought(string &idx, param<double> &wt0, param<double> &wt_own, param<double> &wt_bought);
    double parSum(param<double> w);
    string findMax(const param<double> &w);
    int nthOccurrence(const std::string& str, const std::string& findMe, int nth);
    ~myModel(){ };
};


