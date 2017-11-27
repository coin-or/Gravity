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
//        fname = "/users/hh/Downloads/zeros-2x2-2k_2m.csv";
    }
    
    read_samples(fname);
    for (unsigned i = 0; i<nb_conf; i++) {
        nb_samples_pu = nb_samples.eval(i)/(double)tot_nb_samples;
    }
    
    short val;
    for (unsigned main_spin = 0; main_spin<nb_spins; main_spin++) {
    
        for (unsigned conf = 0; conf<nb_conf; conf++) {
            for (unsigned spin = 0; spin<nb_spins; spin++) {
                if (spin==main_spin) {
                    val =configs.eval(conf,spin);
                    nodal_stat.set_val(conf,spin,val);
                }
                else {
                    val = configs.eval(conf,spin)*configs.eval(conf,main_spin);
                    nodal_stat.set_val(conf,spin,val);
                }
            }
        }
        DebugOn(RED << "############ SPIN NUMBER "<< main_spin << " ############\n "<< RESET);
        
        Model Ising("Ising Model");
        /** Variables */
        var<Real> x("x"), z("z", pos_), f("f"), g("g", pos_), obj("obj");
        Ising.add_var(x^nb_spins);
        Ising.add_var(z^nb_spins);
//        Ising.add_var(f^nb_conf);
//        Ising.add_var(g^nb_conf);
        Ising.add_var(obj^1);
        Ising.min(obj);
        
        /** Constraints */
//        Constraint Lin("Lin");
//        Lin += f + product(nodal_stat,x);
//        Ising.add_constraint(Lin=0);
//
//        Constraint Exp("exp");
//        Exp += g - expo(f);
//        Ising.add_constraint(Exp>=0);
        
//        Constraint NLin("NLin");
//        NLin += f - expo(-1*product(nodal_stat,x));
//        Ising.add_constraint(NLin>=0);

        Constraint Absp("Absp");
        Absp += z - x;
        Ising.add_constraint(Absp >= 0);
        Constraint Absn("Absn");
        Absn += z + x;
        Ising.add_constraint(Absn >= 0);
        
        Constraint Obj("Obj");
//        Obj += obj - product(nb_samples_pu,g) - lambda*sum(z.excl(main_spin));
        Obj += obj - product(nb_samples_pu,expo(-1*product(nodal_stat,x))) - lambda*sum(z.excl(main_spin));
        Obj.set_first_derivative(x, (nodal_stat.tr()*(expo(-1*product(nodal_stat,x))).tr())*nb_samples_pu.vec());
        Obj.set_second_derivative(x,x,(nodal_stat.tr()*(expo(-1*product(nodal_stat,x))).tr())*(-1*product(nb_samples_pu,nodal_stat)));
        Ising.add_constraint(Obj>=0);
        
        /** Solver */
        solver NLP(Ising,ipopt);
        NLP.run(log_lev=0,relax=false,"ma57",1e-12,"yes");
    }
    return 0;
}



