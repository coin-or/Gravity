//
//
//  Gravity
//
//  Created by Hassan Hijazi on August 29 2018.
//
//


#include <iostream>
#include <string>
#include <fstream>
#include <gravity/solver.h>
#include <stdio.h>
#include <stdlib.h>
#include <optionParser.hpp>

#define DEPTH 20
#define HEIGHT 3

using namespace std;
using namespace gravity;


int main (int argc, char * argv[])
{
    #ifdef USE_QPP
    //  Start Timers
    string path = argv[0];
    int output = 0;
    bool relax = false;
    
    double tol = 1e-6;
    string mehrotra = "no", log_level="0";
    
    /** create a OptionParser with options */
    op::OptionParser opt;
    opt.add_option("h", "help", "shows option help"); // no default value means boolean options, which default value is false
    opt.add_option("l", "log", "Log level (def. 0)", log_level );
    
    /** parse the options and verify that all went well. If not, errors and help will be shown */
    bool correct_parsing = opt.parse_options(argc, argv);
    
    if(!correct_parsing){
        return EXIT_FAILURE;
    }
    
    output = op::str2int(opt["l"]);
    output = 5;
    bool has_help = op::str2bool(opt["h"]);
    /** show help */
    if(has_help) {
        opt.show_help();
        exit(0);
    }
    
    /** Parameters **/
    
    /* Qubit number */
    unsigned n = 3;
    unsigned m = pow(2,n);
    /* Circuit depth */
    unsigned d = 1;
    
    /* T and T conjugate transpose gate matrices, one for each row (qubit) and for each column (gate position) */
    param<> T[n][d], Tt[n][d];
    
    /* Hadamard gate matrices */
    param<> H[n][d];
    
    /* Cnot gate matrices */
    vector<param<>> Cnot[d];
    
    /** Initialization **/
    for (auto j = 0; j<d; j++) { // Depth (column) iterator
        for (auto i = 0; i<n; i++) { // Qubit (row) iterator
            T[i][j].set_name("T_"+to_string(i+1)+","+to_string(j+1));
            T[i][j].set_size(2*m, 2*m, 0);/* Representing Complex number as 2x2 real matrices */
            T[i][j].QuantumT(i,n);/* Fill nonzero values corresponding to a T gate */
            Tt[i][j].set_name("Tt_"+to_string(i+1)+","+to_string(j+1));
            Tt[i][j].set_size(2*m, 2*m, 0);
            Tt[i][j].QuantumT(i,n,true);/* Fill nonzero values corresponding to a T conjugate transpose gate */
            H[i][j].set_name("H_"+to_string(i+1)+","+to_string(j+1));
            H[i][j].set_size(2*m, 2*m, 0);
            H[i][j].QuantumH(i,n);/* Fill nonzero values corresponding to an H gate */
            for (auto k = i+1; k<n; k++) {
                Cnot[j].push_back(param<>("Cnot_"+to_string(i+1)+to_string(k+1)+","+to_string(j+1)));
                Cnot[j].back().set_size(2*m, 2*m, 0);
                Cnot[j].back().QuantumCnot(i,k,n);/* Fill nonzero values corresponding to a CNOT gate */
            }
//            exit(-1);
        }
        
    }
    
    
    auto total_time_start = get_wall_time();
    auto total_time_end = get_wall_time();
    auto total_time = total_time_end - total_time_start;
    DebugOn("Total Computing Time = " << total_time << endl);
#else
    cerr << "Error: this version of Gravity "
    "was compiled without QPP support. Please rerun cmake with -DENABLE_QPP=1." << endl;
#endif
    
    return 0;
}
