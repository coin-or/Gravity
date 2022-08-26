//
//  Sensor_Assign2.hpp
//  sensor_assign
//
//  Created by Svetlana Riabova on 8/18/22.
//

#pragma once
#include <stdio.h>
#include <random>
#include <gravity/solver.h>

using namespace std;
using namespace gravity;

class myModel {
public:
    Model<> model;
    indices sensors;
    indices objects;
    indices agents;
    int M;
    int N;
    int K = 5;
    double e;
    double max_price;
    Net graph;
    indices arcs;
    indices own_arcs;
    indices bought_arcs;
    indices agents_arcs;
    indices own_sens;
    indices bought_sens;
    indices own_rplc;
    indices oths_rplc;
    indices own_oths_rplc1;
    indices own_oths_rplc2;
    indices operations; //agents X sensors
    vector<int> owner;
    
public:
    myModel(){ };
    vector<param<double>> readData(int argc, const char * argv[], int n1, int n2);
    void InitBilevel(param<double> &w0, param<double> &w_own, double eps);
    void mSolve();
    void readGreedySol(string fname);
    //void saveSolStats();
    //void GreedyStart(const param<double> &w0, const param<double> &w_own, const param<double> &w_bought);
    //void assignLeader(string &idx, param<double> &wt0, param<double> &wt_own, param<double> &wt_bought);
    //void assignOwn(string &idx, param<double> &wt0, param<double> &wt_own, param<double> &wt_bought);
    //void assignBought(string &idx, param<double> &wt0, param<double> &wt_own, param<double> &wt_bought);
    //void printAssignment();
    //double parSum(param<double> w);
    //string findMax(const param<double> &w);
    //int nthOccurrence(const std::string& str, const std::string& findMe, int nth);
    ~myModel(){ };
};


