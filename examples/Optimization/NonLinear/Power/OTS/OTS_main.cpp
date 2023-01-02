//
//  OTS.cpp
//  Gravity
//
//  ACOPF Created by Hassan Hijazi on Dec 7 2017
//  OTS Created by Smitha Gopinath Oct 20 2022
//
#include <stdio.h>
#include <PowerNet.h>
#include <gravity/solver.h>
#ifdef USE_OPT_PARSER
#include <optionParser.hpp>
#endif

using namespace std;
using namespace gravity;



int main (int argc, char * argv[])
{
    string fname = string(prj_dir)+"/data_sets/Power/nesta_case5_pjm.m", mtype = "ACPOL";
    int output = 0;
    double tol = 1e-6;
    string mehrotra = "no", log_level="0";
#ifdef USE_OPT_PARSER
    /** create a OptionParser with options */
    op::OptionParser opt;
    opt.add_option("h", "help", "shows option help"); // no default value means boolean options, which default value is false
    opt.add_option("f", "file", "Input file name (def. ../data_sets/Power/nesta_case5_pjm.m)", fname );
    opt.add_option("l", "log", "Log level (def. 0)", log_level );
    opt.add_option("m", "model", "power flow model: ACPOL/ACRECT (def. ACPOL)", mtype );
    
    /** parse the options and verify that all went well. If not, errors and help will be shown */
    bool correct_parsing = opt.parse_options(argc, argv);
    
    if(!correct_parsing){
        return EXIT_FAILURE;
    }
    
    fname = opt["f"];
    mtype = opt["m"];
    output = op::str2int(opt["l"]);
    output = 5;
    bool has_help = op::str2bool(opt["h"]);
    /** show help */
    if(has_help) {
        opt.show_help();
        exit(0);
    }
#else
    if(argc==2){
        fname=argv[1];
    }
    else{
        fname=string(prj_dir)+"/data_sets/Power/nesta_case5_pjm.m";
    }
#endif
    
   
    
    
    return 0;
}
