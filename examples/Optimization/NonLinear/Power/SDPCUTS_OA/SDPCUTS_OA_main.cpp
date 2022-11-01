//
// Created by kbestuzheva on 12/11/17.
//
#include <stdio.h>
#include <iostream>
#include <string>
#include <stdio.h>
#include <cstring>
#include <fstream>
#include "../PowerNet.h"
#include <gravity/solver.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef USE_OPT_PARSER
#include <optionParser.hpp>
#endif

using namespace std;
using namespace gravity;


/* main */
int main (int argc, char * argv[]) {
    int output = 0;
    bool sdp_cuts = true;
    
    bool current_from = true, llnc=true, current_to=true, loss=true, loss_bounds=true, current;
    
    size_t num_bags = 0;
    string num_bags_s = "100";
    string solver_str = "ipopt";
    string sdp_cuts_s = "yes";
    string current_from_s = "yes";
    string orig_s = "yes";
    string current_to_s="yes";
    string lazy_s = "no";
    bool lazy_bool = false;
    SolverType solv_type =ipopt;
    double tol = 1e-6;
    string mehrotra = "no";
    
    
    string fname = string(prj_dir)+"/data_sets/Power/nesta_case5_pjm.m";
    
    // create a OptionParser with options
    
#ifdef USE_OPT_PARSER
    op::OptionParser opt;
    opt.add_option("h", "help",
                   "shows option help"); // no default value means boolean options, which default value is false
    opt.add_option("f", "file", "Input file name", fname);
    opt.add_option("s", "solver", "Solvers: ipopt/cplex/gurobi, default = ipopt", solver_str);
    opt.add_option("b", "numbags", "Number of bags per iteration", num_bags_s);
    opt.add_option("If", "current_from", "add from current constraints", current_from_s);
    opt.add_option("o", "original", "add original variables and linking constraints", orig_s);
    opt.add_option("It", "current_to", "add to current constraints", current_to_s);
    opt.add_option("lz", "lazy", "Generate 3d SDP cuts in a lazy fashion, default = no", lazy_s);
    // parse the options and verify that all went well. If not, errors and help will be shown
    bool correct_parsing = opt.parse_options(argc, argv);
    
    if (!correct_parsing) {
        return EXIT_FAILURE;
    }
    
    fname = opt["f"];
    bool has_help = op::str2bool(opt["h"]);
    if (has_help) {
        opt.show_help();
        exit(0);
    }
    solver_str = opt["s"];
    if (solver_str.compare("gurobi")==0) {
        solv_type = gurobi;
    }
    else if(solver_str.compare("cplex")==0) {
        solv_type = cplex;
    }else if(solver_str.compare("Mosek")==0) {
        solv_type = _mosek;
    }
    lazy_s = opt["lz"];
    if (lazy_s.compare("no")==0) {
        lazy_bool = false;
    }
    else if(lazy_s.compare("yes")==0) {
        lazy_bool = true;
    }
    
    current_from_s = opt["If"];
    if (current_from_s.compare("no")==0) {
        current_from = false;
    }
    else {
        current_from = true;
    }
    bool add_original = true;
    
    orig_s = opt["o"];
    if (orig_s.compare("no")==0) {
        add_original = false;
    }
    else {
        add_original = true;
    }
    
    
    
    current_to_s = opt["It"];
    if (current_to_s.compare("no")==0) {
        current_to = false;
    }
    else {
        current_to = true;
    }
    num_bags = atoi(opt["b"].c_str());
    
    current_from=true;
    current_to=true;
    loss=true;
#else
    if(argc>=2){
        fname=argv[1];
        DebugOn(fname<<endl);
    }
    else{
        fname=string(prj_dir)+"/data_sets/Power/nesta_case5_pjm.m";
    }
    
    DebugOn(fname<<endl);
    
    
    
#endif
    
    cout << "\nnum bags = " << num_bags << endl;
    
    // double total_time_start = get_wall_time();
    PowerNet grid;
    grid.readgrid(fname);
    //grid.update_ref_bus();
    
    grid.get_tree_decomp_bags();
    //auto bags_3d=grid.decompose_bags_3d();
    
    bool print_bags = false, only_3d_bags = false;
    auto bags_3d=grid.decompose_bags_3d_linear(print_bags, only_3d_bags);

    double upper_bound;
    auto OPF=build_ACOPF(grid, ACRECT);
    solver<> UB_solver(OPF,ipopt);
    UB_solver.run(output = 0, 1e-6,"ma27");
    if(OPF->_status!=0){
        upper_bound=OPF->_obj->_range->second;
    }
    else{
        upper_bound=OPF->get_obj_val();
    }

   
    std::vector<pair<int,std::vector<string>>> _bag_names;
    int count=0;
    for(auto b:grid._bags){
          pair<int,vector<string>> bn;
          bn.first=count++;
//            indices zk_nodes;
//            param<double> midk("midk"+to_string(count));

          for(auto n:b.second){
              bn.second.push_back(n->_name);
         }
        _bag_names.push_back(bn);
    }
    
   
   // int count=0;
//    for(auto b:bags_3d){
//          pair<int,vector<string>> bn;
//          bn.first=count++;
////            indices zk_nodes;
////            param<double> midk("midk"+to_string(count));
//        string bname=b.first;
//        size_t pos=0;
//        while(true){
//            pos= bname.find_first_of(",");
//            if(pos!=string::npos){
//                bn.second.push_back(bname.substr(0, pos));
//                bname=bname.substr(pos+1);
//            }
//            else{
//                bn.second.push_back(bname);
//                break;
//            }
//        }
//
//        _bag_names.push_back(bn);
//    }
    
    
    double solver_time_start;
    solver_time_start = get_wall_time();
    auto SDPa=build_SDPOPF(grid, true, true, true);
    SDPa->print();
    initialize_relaxation(OPF, SDPa, grid, current);
    SDPa->_bag_names=_bag_names;
    SDPa->sdp_dual=false;
    solver<> LBnonlin_solver(SDPa,ipopt);
    LBnonlin_solver.run(output = 5 , 1e-6, "ma57");
    int while_count=0;
    if(solv_type==ipopt){
        count=0;
        while(while_count++<=3000){
        auto res=SDPa->cuts_eigen_bags_primal_complex(1e-6, "Wii", "R_Wij", "Im_Wij");
        if(res.size()>=1){
            for(auto i=0;i<res.size();i++){
                Constraint<> cut("cut"+to_string(count++));
                int j=0;
                for(j=0;j<res[i].size()-1;j+=2){
                    int c=res[i][j];
        
                    for(auto it=SDPa->_vars.begin();it!=SDPa->_vars.end();it++){
                        auto it1=next(it);
                                                if(*it->second->_id<=c && (*it1->second->_id>c || it1==SDPa->_vars.end())){
                                                    auto v= SDPa->get_var<double>(it->second->_name);
                                                    cut+=v(v._indices->_keys->at(c-*it->second->_id))*res[i][j+1];
                                                    break;
                                                }
                                            }
                        
                }
                 cut+=res[i][j];
                SDPa->add(cut<=0);
                //cut.print();
                DebugOff("cut"<<endl);
            }
        }
        else{
            break;
        }
            auto lower_bound = SDPa->get_obj_val();
            
            auto gap = 100*(upper_bound - lower_bound)/upper_bound;
            DebugOn("gap "<<gap<<endl);
            if(gap<=1){
                break;
            }
        //SDPa->reindex();
            LBnonlin_solver.run(output = 0 , 1e-6, "ma57");
        }
    }
//    else
//    SDPOPF.run();

    
    
    double solver_time_end = get_wall_time();
    double solver_time=solver_time_end-solver_time_start;
    double gap=999, lower_bound=999;
    // SDP->print_solution();
    // SDP->print();
    SDPa->print_constraints_stats(tol);
    SDPa->print_nonzero_constraints(tol,true);
    if(SDPa->_status==0)
    {
        lower_bound = SDPa->get_obj_val();
        
        gap = 100*(upper_bound - lower_bound)/upper_bound;
    }
    
    //    auto solve_time = solver_time_end - solver_time_start;
    //
    //    string out = "\nDATA_OPF, " + grid._name + ", " + to_string(nb_buses) + ", " + to_string(nb_lines) +", " + to_string(lower_bound) + ", " + to_string(-numeric_limits<double>::infinity()) + ", " + to_string(solve_time) + ", LocalOptimal, " + to_string(total_time);
    //    DebugOn(out <<endl);
    DebugOn("Final Gap = " << to_string(gap) << "%."<<endl);
    DebugOn("Upper bound = " << to_string(upper_bound) << "."<<endl);
    DebugOn("Lower bound = " << to_string(lower_bound) << "."<<endl);
    //SDP->print();
    string result_name=string(prj_dir)+"/results_SDP/"+grid._name+".txt";
    
    ofstream fout(result_name.c_str());
    fout<<grid._name<<"\t"<<std::fixed<<std::setprecision(5)<<gap<<"\t"<<std::setprecision(5)<<upper_bound<<"\t"<<std::setprecision(5)<<lower_bound<<"\t"<<std::setprecision(5)<<solver_time<<endl;
    fout.close();
    bool not_sdp=true;
   
    
    //    SDP->run_obbt();
    //    SDP->reset_constrs();
    //    solver<> SDPLB1(SDP,solv_type);
    //
    //    auto status = SDPLB1.run(output = 5, tol);
    //    SDP->print_constraints_stats(tol);
    //    bool print_only_relaxed;
    //    SDP->print_nonzero_constraints(tol,print_only_relaxed=true);
    //
    //    //        SDP->print_solution();
    //
    //    //        SDP->print();
    //
    //    if(status==0)
    //    {
    //        total_time_end = get_wall_time();
    //        total_time = total_time_end - total_time_start;
    //        DebugOn("\nResults: " << grid._name << " " << to_string(SDP->get_obj_val()) << " " <<endl);
    //        SDP->print_constraints_stats(tol);
    //
    //        DebugOn("Initial Gap Nonlinear = " << to_string(gap) << "%."<<endl);
    //        lower_bound=SDP->get_obj_val()*upper_bound;
    //        gap = 100*(upper_bound - lower_bound)/upper_bound;
    //        DebugOn("Final Gap = " << to_string(gap) << "%."<<endl);
    //        DebugOn("Upper bound = " << to_string(upper_bound) << "."<<endl);
    //        DebugOn("Lower bound = " << to_string(lower_bound) << "."<<endl);
    //        DebugOn("Time\t"<<total_time<<endl);
    //    }
    //    else {
    //        DebugOn("WARNING: Relaxation did not converge!"<<endl);
    //    }

   // SDP->cuts_eigen_bags_primal_complex(1e-6, "Wii", "R_Wij", "Im_Wij");
    
    
}
