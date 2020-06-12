#include <stdio.h>
#include <iostream>
#include <string>
#include <stdio.h>
#include <cstring>
#include <fstream>
#include <utility>
#include "../PowerNet.h"
#include <gravity/solver.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef USE_OPT_PARSER
#include <optionParser.hpp>
#endif
//#include <math.h>
using namespace std;
using namespace gravity;


/* main */
int main (int argc, char * argv[]) {
#ifdef USE_MPI
    auto err_init = MPI_Init(nullptr,nullptr);
    int worker_id, nb_workers;
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
    auto err_size = MPI_Comm_size(MPI_COMM_WORLD, &nb_workers);
#endif
    int output = 0;
    bool current = false;
    size_t num_bags = 0;
    string num_bags_s = "100";
    string solver_str = "ipopt";
    string sdp_cuts_s = "yes";
    
    string current_s = "yes";
    string time_s = "1000";
    string sdp_kim_s="yes";
    string threads_s="1";
    
    string lazy_s = "no", linearize_s = "no";
    string orig_s = "no";
    bool lazy_bool = false;
    bool add_original=false;
    bool sdp_kim=true;
    SolverType solv_type = ipopt;
    const double tol = 1e-6;
    string mehrotra = "no";
    bool linearize=false;
    
    
    string fname = string(prj_dir)+"/data_sets/Power/nesta_case9_bgm__nco.m";
    
    
#ifdef USE_OPT_PARSER
    
    //create a OptionParser with options
    op::OptionParser opt;
    opt.add_option("h", "help",
                   "shows option help");//No option means Default values which may be seen above using option strings
    opt.add_option("f", "file", "Input file name", fname);
    opt.add_option("t", "time", "Time limit, defaut 60 secs", time_s);
    opt.add_option("threads", "threads", "Number of threads, defaut 1", threads_s);
    opt.add_option("s", "solver", "Solvers: ipopt/cplex/gurobi, default = ipopt", solver_str);
    opt.add_option("b", "numbags", "Number of bags per iteration", num_bags_s);
    opt.add_option("l", "", "add current constraints", current_s); //Adds loss from of true and if true also adds loss_to in a lazy fasion
    opt.add_option("lz", "lazy", "Generate 3d SDP cuts in a lazy fashion, default = no", lazy_s);
    opt.add_option("linear", "linear", "Linearize relaxation using OA cuts, default = no", linearize_s);
    opt.add_option("o", "original", "add original variables and linking constraints", orig_s);
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
    threads_s = opt["threads"];
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
    linearize_s = opt["linear"];
    if (linearize_s.compare("yes")==0) {
        linearize = true;
    }
    
    current_s = opt["l"];
    if (current_s.compare("no")==0) {
        current = false;
    }
    else {
        current = true;
    }
    
    orig_s = opt["o"];
    if (orig_s.compare("no")==0) {
        add_original = false;
    }
    else {
        add_original = true;
    }
    num_bags = atoi(opt["b"].c_str());
    
    auto max_time = op::str2double(opt["t"]);
#else
    if(argc>=4){
        fname=argv[1];
        time_s=argv[2];
        sdp_kim_s=argv[3];
        
    }
    if(argc>4){
        linearize_s=argv[4];
    }
    if (linearize_s.compare("yes")==0) {
        linearize = true;
    }

    current=true;
    
    auto max_time=std::atoi(time_s.c_str());
    if (sdp_kim_s.compare("no")==0) {
        sdp_kim = false;
    }
    else {
        sdp_kim = true;
    }
#endif
    
    current=true;
    auto nb_threads=std::atoi(threads_s.c_str());
    PowerNet grid;
    grid.readgrid(fname);
    grid.get_tree_decomp_bags();
    
    auto nodes = indices(grid.nodes);
    auto arcs = indices(grid.arcs);
    auto v_max = grid.v_max.in(nodes);
    auto bags_3d=grid.decompose_bags_3d();
    auto node_pairs = grid.get_node_pairs();
    auto node_pairs_chord = grid.get_node_pairs_chord(bags_3d);
    
    auto c1 = grid.c1.in(grid.gens);
    auto c2 = grid.c2.in(grid.gens);
    auto c0 = grid.c0.in(grid.gens);
  
    DebugOn("Machine has " << thread::hardware_concurrency() << " threads." << endl);
    
    int nb_total_threads = nb_threads; /** Used when MPI is ON to multipply with the number of workers */
#ifdef USE_MPI
    nb_total_threads *= nb_workers;
#endif
    double lower_bound=numeric_limits<double>::min(), lower_bound_nonlin_init=numeric_limits<double>::min(),total_time=numeric_limits<double>::min();

    auto OPF=build_ACOPF(grid, ACRECT);
    //OPF->print();
    double ub_solver_tol=1e-8, lb_solver_tol=1e-8, range_tol=1e-3, opt_rel_tol=1e-2, opt_abs_tol=1e6;
    int total_iter;
    unsigned max_iter=1e3;
    int oacuts=0, oacuts_init=0;
    SolverType ub_solver_type = ipopt, lb_solver_type = solv_type;
    bool scale_objective=true;
    //linearize=true;
    if(!linearize){
        auto nonlin_obj=true;
        current=true;
        auto SDP= build_SDPOPF(grid, current, nonlin_obj, sdp_kim);
        auto res=OPF->run_obbt(SDP, max_time, max_iter, opt_rel_tol, opt_abs_tol, nb_threads=1, ub_solver_type, lb_solver_type, ub_solver_tol, lb_solver_tol, range_tol, linearize, scale_objective);
        lower_bound = get<6>(res);
        lower_bound_nonlin_init = get<3>(res);
        total_iter=get<1>(res);
        total_time=get<2>(res);
    }
    else{
        current=true;
        auto nonlin_obj=false;
        auto SDP= build_SDPOPF(grid, current, nonlin_obj, sdp_kim);
        auto res=OPF->run_obbt(SDP, max_time, max_iter, opt_rel_tol, opt_abs_tol, nb_threads=1, ub_solver_type, lb_solver_type, ub_solver_tol, lb_solver_tol, range_tol, linearize, scale_objective);
        lower_bound = get<6>(res);
        lower_bound_nonlin_init = get<3>(res);
        total_iter=get<1>(res);
        total_time=get<2>(res);
        oacuts=get<8>(res);
        oacuts_init=get<9>(res);
    }
    string result_name=string(prj_dir)+"/results_obbt/"+grid._name+".txt";
    
    auto upper_bound = OPF->get_obj_val();
    auto gap_init = 100*(upper_bound - lower_bound_nonlin_init)/std::abs(upper_bound);
    auto final_gap = 100*(upper_bound - lower_bound)/std::abs(upper_bound);
#ifdef USE_MPI
    if(worker_id==0){
        ofstream fout(result_name.c_str());
        fout<<grid._name<<"\t"<<std::fixed<<std::setprecision(5)<<gap_init<<"\t"<<std::setprecision(5)<<upper_bound<<"\t"<<std::setprecision(5)<<lower_bound<<"\t"<<std::setprecision(5)<<final_gap<<"\t"<<total_iter<<"\t"<<std::setprecision(5)<<total_time<<"\t"<<oacuts<<"\t"<<oacuts_init<<"\t"<<endl;
        if(lower_bound==numeric_limits<double>::min()){
            fout<<"Lower bound not solved to optimality"<<endl;
        }
        DebugOn("I am worker id "<<worker_id<<" writing to results file "<<endl);
        fout.close();
    }
    MPI_Finalize();
#else
    DebugOn(grid._name<<"\t"<<std::fixed<<std::setprecision(5)<<gap_init<<"\t"<<std::setprecision(5)<<upper_bound<<"\t"<<std::setprecision(5)<<lower_bound<<"\t"<<std::setprecision(5)<<final_gap<<"\t"<<total_iter<<"\t"<<std::setprecision(5)<<total_time<<"\t"<<oacuts<<"\t"<<oacuts_init<<"\t"<<endl);
    
    
    ofstream fout(result_name.c_str());
    fout<<grid._name<<"\t"<<std::fixed<<std::setprecision(5)<<gap_init<<"\t"<<std::setprecision(5)<<upper_bound<<"\t"<<std::setprecision(5)<<lower_bound<<"\t"<<std::setprecision(5)<<final_gap<<"\t"<<total_iter<<"\t"<<std::setprecision(5)<<total_time<<"\t"<<oacuts<<"\t"<<oacuts_init<<"\t"<<endl;
    if(lower_bound==numeric_limits<double>::min()){
        fout<<"Lower bound not solved to optimality"<<endl;
    }
    fout.close();
#endif
    
    return 0;
}

