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
    vector<int> owner;/*< vector storing ownership of each sensor */
    vector<double> budget;
    map<string,int> owner_map;/*< map storing ownership of each sensor */
    
public:
    myModel(){ };
    vector<param<double>> readData(int argc, const char * argv[], int n1, int n2);
    vector<param<double>> readHD5(const string& file_name);
    void InitBilevel(param<double> &w0, param<double> &w_own, param<double> &w_bought, double eps);
    void mSolve(bool run_mip=false);
    void saveSolStats();
    void GreedyStart(const param<double> &w0, const param<double> &w_own, const param<double> &w_bought);
    void writeGreedySol();
    void assignLeader(string &idx, param<double> &wt0, param<double> &wt_own, param<double> &wt_bought);
    void assignOwn(string &idx, param<double> &wt0, param<double> &wt_own, param<double> &wt_bought);
    void assignBought(string &idx, param<double> &wt0, param<double> &wt_own, param<double> &wt_bought);
    double parSum(const param<double>& w);
    int nthOccurrence(const std::string& str, const std::string& findMe, int nth);
    ~myModel(){ };
};


