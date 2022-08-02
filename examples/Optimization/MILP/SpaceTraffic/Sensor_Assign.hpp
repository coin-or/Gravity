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
private:
    int N;
    int M;
    int K = 5;
    Net graph;
    indices sensors;
    indices objects;
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
    vector<param<double>> readData(int argc, const char * argv[]);
    void InitBilevel(param<double> w0, param<double> w_own, param<double> w_bought);
    //void GreedyStart(Model<>& model, param<double> w0, param<double> w_own, param<double> w_bought);
    //string findMax(param<double> w);
    int nthOccurrence(const std::string& str, const std::string& findMe, int nth);
    ~myModel(){ };
};


