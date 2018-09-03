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
#include <thread>
#include <optionParser.hpp>

using namespace std;
using namespace gravity;

vector<double> solver_time;
vector<double> obj_val;
double x_sum;



int nb_cols = 1, nb_conf = 0, nb_spins = 0, tot_nb_samples = 0;
param<short> configs("configs");
param<int> nb_samples("nb_samples");
param<double> nb_samples_pu("nb_samples_pu");
vector<vector<double>> solution;
double regularizor=0.0, lambda;

void write_sol(string input_file, string output_file=""){
    ofstream fs;
    // create a name for the output file
    std::string filename = input_file;
    if (output_file=="") {
        filename = filename.substr(0,filename.size()-4);//remove csv extension
        filename +="_solution.csv";
    }
    else {
        filename = output_file;
    }
    fs.open(filename);
    for (unsigned main_spin = 0; main_spin<nb_spins; main_spin++) {
        for (unsigned spin = 0; spin<nb_spins-1; spin++) {
            fs << solution[main_spin][spin] << ",";
        }
        fs << solution[main_spin][nb_spins-1] << endl;
    }
    fs.close();
}

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
            DebugOff(val << ", ");
            spin++;
        }
        DebugOff(endl);
        nb_conf++;
    }
    DebugOn("Number of configurations = " << nb_conf << endl);
    DebugOn("Total number of samples = " << tot_nb_samples << endl);
    for (unsigned i = 0; i<nb_conf; i++) {
        nb_samples_pu = nb_samples.eval(i)/(double)tot_nb_samples;
    }
    lambda = regularizor*sqrt(log((pow(nb_spins,2))/0.05)/tot_nb_samples);
    DebugOn("Lambda = " << lambda << endl);
    obj_val.resize(nb_spins);
    solver_time.resize(nb_spins);
    return true;
}


void solve_spin(unsigned spin1, unsigned spin2, int log_lev=0, bool relax=false, string mehrotra="yes"){
    for (unsigned main_spin = spin1; main_spin<spin2; main_spin++) {
        param<short> nodal_stat("nodal_stat");
        double val;
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
        var<double> x("x"), z("z", pos_), obj("obj");
        auto Rn = R(nb_spins);
        Ising.add_var(x.in(Rn));
        Ising.add_var(z.in(Rn));
        Ising.add_var(obj);
        Ising.min(obj);
        
        /** Constraints */
        
        if (lambda>0) {
            Constraint Absp("Absp");
            Absp += z - x;
            Ising.add_constraint(Absp >= 0);
            Constraint Absn("Absn");
            Absn += z + x;
            Ising.add_constraint(Absn >= 0);
        }
    
        Constraint Obj("Obj");
        if (lambda > 0) {
            Obj += obj - product(nb_samples_pu,expo(-1*product(nodal_stat,x))) - lambda*sum(z.excl(main_spin));
        }
        else {
            Obj += obj - product(nb_samples_pu,expo(-1*product(nodal_stat,x)));
        }
        Obj.set_first_derivative(x, (nodal_stat.tr()*(expo(-1*product(nodal_stat,x))).tr())*nb_samples_pu.vec());
        Obj.set_second_derivative(x,x,(nodal_stat.tr()*(expo(-1*product(nodal_stat,x))).tr())*(-1*product(nb_samples_pu,nodal_stat)));
        Ising.add_constraint(Obj>=0);
        Obj.print();
        
        /** Solver */
        solver NLP(Ising,ipopt);
        auto solver_time_start = get_wall_time();
        NLP.run(log_lev,relax=false,1e-12,1e-6,"ma27",mehrotra);
//        Ising.print_nl_functions();
        auto solver_time_end = get_wall_time();
        solver_time[main_spin] = solver_time_end - solver_time_start;
        obj_val[main_spin] = Ising._obj_val;
        solution[main_spin].resize(nb_spins);
        for (unsigned spin = 0; spin<nb_spins; spin++) {
            solution[main_spin][spin] = x.eval(spin);
        }
    }
}

void print_sol() {
    cout << "[";
    for (unsigned main_spin = 0; main_spin<nb_spins; main_spin++) {
        for (unsigned spin = 0; spin<nb_spins; spin++) {
            cout << solution[main_spin][spin] << " ";
        }
    }
    cout << "]\n";
}

double sum_x(){
    double res = 0;
    for (unsigned main_spin = 0; main_spin<nb_spins; main_spin++) {
        for (unsigned spin = 0; spin<nb_spins; spin++) {
            res += solution[main_spin][spin];
        }
    }
    return res;
}

int main (int argc, char * argv[])
{
    int output = 5;
    bool relax = false;
    string fname = "../data_sets/Ising/samples_bin_sml.csv", log_level="5";
    string path = argv[0];
    if (path.find("/bin")!=string::npos && path.find("/bin/ising")==string::npos) {//Not running from terminal
        fname = "../" + fname;
    }
    unsigned nr_threads = 1;
    
    // create a OptionParser with options
    op::OptionParser opt;
    opt.add_option("h", "help", "shows option help"); // no default value means boolean options, which default value is false
    opt.add_option("f", "file", "Input file name", fname );
    opt.add_option("t", "threads", "Number of threads to use", "1");
    opt.add_option("r", "regularizer", "Value of regularizer", "0.0");
    opt.add_option("l", "log", "Log level (def. 0)", log_level );
    
    // parse the options and verify that all went well. If not, errors and help will be shown
    bool correct_parsing = opt.parse_options(argc, argv);
    
    if(!correct_parsing){
        return EXIT_FAILURE;
    }
    
    fname = opt["f"];
    nr_threads = op::str2int(opt["t"]);
    regularizor = op::str2double(opt["r"]);
    output = op::str2int(opt["l"]);
    bool has_help = op::str2bool(opt["h"]);
    // show help
    if(has_help) {
        opt.show_help();
        exit(0);
    }
    if (nr_threads<1) {
        nr_threads=1;
    }
    DebugOn("Using " << nr_threads << " threads" << endl);
    DebugOn("Regularizor = " << regularizor << endl);
    string mehrotra="yes";
    if (regularizor==0) {
        mehrotra="no";//IPOPT seems to be failing when mehrotra is on and the reg is zero..
    }
    auto total_time_start = get_wall_time();
    read_samples(fname.c_str());
    solution.resize(nb_spins);
    
    vector<thread> threads;
    /* Split subproblems into nr_threads parts */
    vector<int> limits = bounds(nr_threads, nb_spins);
    /* Launch all threads in parallel */
    for (int i = 0; i < nr_threads; ++i) {
        threads.push_back(thread(solve_spin, limits[i], limits[i+1], output, relax, mehrotra));
    }
    /* Join the threads with the main thread */
    for(auto &t : threads){
        t.join();
    }
//    print_sol();
    write_sol(fname);
    auto total_time_end = get_wall_time();
    auto total_time = total_time_end - total_time_start;
    fname = strtok(strdup(fname.c_str()), "/");
    char* key;
    while((key = strtok(NULL,"/"))!=NULL){
        fname = key;
    }
    string out = "DATA_INV_ISING, " + fname + ", " + to_string(nb_spins) + ", " + to_string(nb_conf) + ", " + to_string(tot_nb_samples) + ", " + to_string(sum_x()) + ", " + to_string(total_time);;
//    for (unsigned spin = 0; spin <nb_spins; spin++) {
//        out += ", " + to_string(obj_val[spin]) + ", " + to_string(solver_time[spin]);
//    }
    DebugOn(out <<endl);
    return 0;
}



