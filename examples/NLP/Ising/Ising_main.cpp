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

int nb_cols = 1, nb_conf = 0, nb_spins = 0, tot_nb_samples = 0;
param<short> configs("configs"), nodal_stat("nodal_stat");
param<int> nb_samples("nb_samples");
param<double> nb_samples_pu("nb_samples_pu");
double regularizor, lambda;


bool read_samples(const char* fname){
    int val;
    char* key;
    io::LineReader in0(fname);
    char* line = in0.next_line();
    strtok(line,",");
    while((key = strtok(NULL,","))!=NULL){
        nb_cols++;
    }
    nb_spins = nb_cols - 1;
    DebugOn("Number of spins = " << nb_spins << endl);
    io::LineReader in(fname);
    unsigned spin = 0;
    while(char* line = in.next_line()){
        nb_samples = stoi(strtok(line,","));
        tot_nb_samples += nb_samples.eval();
        DebugOff("Number of samples = " << nb_samples << endl);
        DebugOff("val = ");
        DebugOff("nodal stat = ");
        spin = 0;
        while((key = strtok(NULL,","))!=NULL){
            val = stoi(key);
            DebugOff(val << ", ");
            configs.set_val(nb_conf, spin, val);
            if (spin==0) {
                nodal_stat.set_val(nb_conf,spin,val);
            }
            else {
                val *= configs.eval(nb_conf,0);
                nodal_stat.set_val(nb_conf,spin,val);
            }
            DebugOff(val << ", ");
            spin++;
        }
        DebugOff(endl);
        nb_conf++;
    }
    DebugOn("Number of configurations = " << nb_conf << endl);
    DebugOn("Total number of samples = " << tot_nb_samples << endl);
    regularizor = 0.2;
    lambda = regularizor*sqrt(log((pow(nb_spins,2))/0.05)/tot_nb_samples);
    DebugOn("Lambda = " << lambda << endl);

    return true;
}

int main (int argc, const char * argv[])
{
    int log_lev = 0;
    bool relax = false;
    const char* fname;
    if (argc >= 2) {
        fname = argv[1];
    }
    else {
           fname = "../../data_sets/Ising/samples_bin_sml.csv";
    }
    
    read_samples(fname);
    for (unsigned i = 0; i<nb_conf; i++) {
        nb_samples_pu = nb_samples.eval(i)/(double)tot_nb_samples;
    }
    Model Ising("Ising Model");
    /** Variables */
    
//    for (unsigned i = 0; i < configs._dim[0]; i++) {
//        for (unsigned j = 0; j < configs._dim[1]; j++) {
//            nodal_stat(i,j) = configs(i,0).eval()*configs(i,j).eval();
//        }
//    }
//    nodal_stat  = [ samples[k,current_row] * (j == current_row ? 1 : samples[k,j]) for k=1:num_conf, j=2:num_row]
    /** Variables */
//    var<Real> x("x"), z("z", pos_), f("f"), g("g", pos_), obj("obj");
    var<Real> x("x"), z("z", pos_), f("f"), obj("obj");
    Ising.add_var(x^nb_spins);
    Ising.add_var(z^nb_spins);
    Ising.add_var(f^nb_conf);
//    Ising.add_var(g^nb_conf);
    Ising.add_var(obj^1);
    Ising.min(obj);
    Constraint Lin("Lin");
    Lin += f + product(nodal_stat,x);
    Ising.add_constraint(Lin=0);
//    for (unsigned i = 0; i<nb_conf; i++) {
//        Constraint Lin("Lin_"+to_string(i));
//        Lin += f(i);
//        for (int j = 0; j<nb_spins; j++) {
//            Lin += x(j) *nodal_stat(i,j);
////            Lin += configs[i][j]*configs[i][0]*x(j);
//        }
////        Lin += configs[i][0]*x(0);
////        Lin += x(0) *configs(i,0).eval();
//        Ising.add_constraint(Lin=0);
////        Constraint Exp("Exp_"+to_string(i));
////        Exp += g(i) - expo(f(i));
////        Ising.add_constraint(Exp>=0);
//    }
    DebugOn("Lin nb instances = " << Lin._nb_instances << endl);
//    Ising.add_constraint(Lin.in(indices._keys)=0);
    
//    Constraint Exp("exp");
//    Exp += g - expo(f);
//    Ising.add_constraint(Exp>=0);
    
    Constraint Absp("Absp");
    Absp += z - x;
    Ising.add_constraint(Absp >= 0);
    Constraint Absn("Absn");
    Absn += z + x;
    Ising.add_constraint(Absn >= 0);
    
    Constraint Obj("Obj");
//    for (unsigned i = 0; i<nb_conf; i++) {
//        Obj -= nb_samples_pu(i)*g(i);
//    }
//    for (unsigned i=1; i<nb_spins; i++) {
//        Obj -= lambda*z(i);
//    }
//    Obj += obj;
//    Obj += obj - product(nb_samples,g);
//    Obj += obj - sum(nb_samples_pu,g) - lambda*sum(z);
//    Obj += obj - sum(nb_samples_pu,expo(f)) - lambda*sum(z);
    
//    Obj += sum(nb_samples_pu,expo(f));
//    Obj += obj - lambda*sum(z);
    Obj += obj - sum(nb_samples_pu,expo(f)) - lambda*sum(z.excl(0));
//    Obj += obj - lambda*sum(z.excl(0));
//    for (unsigned i = 0; i<10; i++) {
//        Obj -= nb_samples_pu.eval(i)*expo(f(i));
//    }
//    Obj += obj - sum(nb_samples_pu,g) - lambda*sum(z.excl(0));
//    Obj += sum(nb_samples_pu,expo(f));
//    Obj += obj - lambda*sum(z) - sum(g);
    Ising.add_constraint(Obj>=0);
    solver NLP(Ising,ipopt);
    NLP.run(log_lev=0,relax=false,"ma57",1e-12);
    return 0;
}

int main_ (int argc, const char * argv[])
{
    int log_lev = 0;
    bool relax = false;
    const char* fname;
    if (argc >= 2) {
        fname = argv[1];
    }
    else {
        fname = "../../data_sets/Ising/samples_bin_sml.csv";
    }
    read_samples(fname);
    for (unsigned i = 0; i<nb_conf; i++) {
        nb_samples_pu = nb_samples.eval(i)/tot_nb_samples;
    }
    Model Ising("Ising Model");
    /** Variables */
    var<Real> x("x"), z("z", pos_), obj("obj");
    Ising.add_var(x^nb_spins);
    Ising.add_var(z^nb_spins);
    Ising.add_var(obj^1);
    Ising.min(obj);
    Constraint Absp("Absp");
    Absp += z - x;
    Ising.add_constraint(Absp >= 0);
    Constraint Absn("Absn");
    Absn += z + x;
    Ising.add_constraint(Absn >= 0);
    Constraint Obj("Obj");
    Obj += obj - sum(nb_samples_pu,expo(-1*product(nodal_stat,x))) - lambda*sum(z);
    Ising.add_constraint(Obj>=0);
    solver NLP(Ising,ipopt);
    NLP.run(log_lev=0,relax=false,"ma57",1e-12);
    return 0;
}



