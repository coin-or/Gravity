//
//  ACOPF.cpp
//  Gravity
//
//  Created by Guanglei Wang on 6/9/17.
//
//
#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <gravity/solver.h>
#include <gravity/csv.h>
#include <stdlib.h>

using namespace std;
using namespace gravity;

size_t nb_cols = 1, nb_conf = 0, nb_spins = 0, tot_nb_samples = 0;
vector<param<int>> configs;
param<int> nb_samples;

bool read_samples(const char* fname){
    int val;
    char* key;
    io::LineReader in0(fname);
    char* line = in0.next_line();
    nb_samples = stoi(strtok(line,","));
    while((key = strtok(NULL,","))!=NULL){
        nb_cols++;
    }
    nb_spins = nb_cols - 1;
    DebugOn("Number of spins = " << nb_spins << endl);
    for (unsigned i = 0; i<nb_spins; i++) {
        configs.push_back(param<int>("configs_spin_"+to_string(i)));
    }
    io::LineReader in(fname);
    unsigned spin = 0;
    while(char* line = in.next_line()){
        nb_samples = stoi(strtok(line,","));
        tot_nb_samples += nb_samples.eval();
        DebugOff("Number of samples = " << nb_samples << endl);
        DebugOff("val = ");
        spin = 0;
        while((key = strtok(NULL,","))!=NULL){
            val = stoi(key);
            DebugOff(val << ", ");
            configs[spin++] = val;
        }
        DebugOff(endl);
        nb_conf++;
    }
    DebugOn("Number of configurations = " << nb_conf << endl);
    DebugOn("Total number of samples = " << tot_nb_samples << endl);
    return true;
}

int main (int argc, const char * argv[])
{
    const char* fname;
    if (argc >= 2) {
        fname = argv[1];
    }
    else {
           fname = "../../data_sets/Ising/samples_bin_sml.csv";
    }
    
    read_samples(fname);
    double regularizor = 0.2;
    double lambda = regularizor*sqrt(log((nb_spins^2)/0.05)/tot_nb_samples);
    DebugOn("Lambda = " << lambda << endl);
    Model Ising("Ising Model");
    /** Variables */
    var<Real> x("x"), z("z", pos_), f("f"), g("g", pos_), obj("obj");
    Ising.add_var(x^nb_spins);
    Ising.add_var(z^nb_spins);
    Ising.add_var(f^nb_conf);
    Ising.add_var(g^nb_conf);
    Ising.add_var(obj);
    Ising.min(obj);
    Constraint Lin("in_lin");
    Lin += f + product(x,configs);
    Ising.add_constraint(Lin=0);
//    solver NLP(Ising,ipopt);
//    NLP.run();
    return 0;
}


