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

void initialize_relaxation(shared_ptr<Model<double>> OPF, shared_ptr<Model<double>> relax, PowerNet& grid, bool current);
/* Run the OBBT algorithm */
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
    string nb_refine_s="10";
    string nb_root_refine_s="100";
    string init_prim_s="no";
    string viol_obbt_init_s="0.1";
    string viol_root_init_s="0.1";
    
    string lazy_s = "no", linearize_s = "no";
    string orig_s = "no";
    bool lazy_bool = false;
    bool add_original=false;
    bool sdp_kim=true;
    SolverType solv_type = ipopt;
    const double tol = 1e-6;
    string mehrotra = "no";
    bool linearize=true;
    bool initialize_primal=false;
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
    if(argc>1){
        fname=argv[1];
    }
    if(argc>2){
        time_s=argv[2];
    }
    if(argc>3){
        threads_s=argv[3];
    }
    if(argc>4){
        sdp_kim_s=argv[4];
    }
    if(argc>5){
        linearize_s=argv[5];
    }
    if(argc>6){
        solver_str=argv[6];
    }
    if(argc>7){
        nb_refine_s=argv[7];
    }
    if(argc>8){
        nb_root_refine_s=argv[8];
    }
    if(argc>9){
        viol_obbt_init_s=argv[9];
    }
    if(argc>10){
        viol_root_init_s=argv[10];
    }
    if(argc>11){
        init_prim_s=argv[11];
    }
    
    
    if (linearize_s.compare("yes")==0) {
        linearize = true;
    }
    if (init_prim_s.compare("yes")==0) {
        initialize_primal = true;
    }
    current=true;
    
    auto max_time=std::atoi(time_s.c_str());
    if (sdp_kim_s.compare("no")==0) {
        sdp_kim = false;
    }
    else {
        sdp_kim = true;
    }
    if (solver_str.compare("gurobi")==0) {
        solv_type = gurobi;
    }
    else if(solver_str.compare("cplex")==0) {
        solv_type = cplex;
    }else if(solver_str.compare("Mosek")==0) {
        solv_type = _mosek;
    }
#endif
    
    current=true;
    auto nb_threads=std::atoi(threads_s.c_str());
    auto nb_refine=std::atoi(nb_refine_s.c_str());
    auto nb_root_refine=std::atoi(nb_root_refine_s.c_str());
    auto viol_obbt_init=std::stod(viol_obbt_init_s.c_str());
    auto viol_root_init=std::stod(viol_root_init_s.c_str());
    
    PowerNet grid;
    grid.readgrid(fname);
    grid.get_tree_decomp_bags();
    
    
    //    auto bags_3d=grid.decompose_bags_3d();
    
    
    
    DebugOn("Machine has " << thread::hardware_concurrency() << " threads." << endl);
    
    int nb_total_threads = nb_threads; /** Used when MPI is ON to multipply with the number of workers */
#ifdef USE_MPI
    nb_total_threads *= nb_workers;
#endif
    double lower_bound=numeric_limits<double>::min(),upper_bound=numeric_limits<double>::min(), lower_bound_nonlin_init=numeric_limits<double>::min(),total_time=numeric_limits<double>::min();
    auto OPF=build_ACOPF(grid, ACRECT);
    double ub_solver_tol=1e-6, lb_solver_tol=1e-6, range_tol=1e-3, opt_rel_tol=1e-2, opt_abs_tol=1e-2, zero_tol=1e-12;
    int total_iter;
    unsigned max_iter=1e3;
    int oacuts=0, oacuts_init=0, fail=0;
    SolverType ub_solver_type = ipopt, lb_solver_type = solv_type;
    bool scale_objective;
    bool termination=true;
    int status=0;
    bool run_obbt=false;
    if(!linearize){
        auto nonlin_obj=true;
        current=true;
        solver<> UB_solver(OPF,ub_solver_type);
        UB_solver.run(output = 5, ub_solver_tol, 2000,"ma27", 2000);
        if(OPF->_status!=0){
            upper_bound=OPF->_obj->_range->second;
        }
        else{
            upper_bound=OPF->get_obj_val();
        }
        auto time_start=get_wall_time();
        auto SDP= build_SDPOPF(grid, current, nonlin_obj, sdp_kim);
        initialize_relaxation(OPF, SDP, grid, current);
        solver<> LBnonlin_solver(SDP,lb_solver_type);
        LBnonlin_solver.run(output = 5 , 1e-6, 5400, "ma57", 5000);
        auto time_end = get_wall_time();
        SDP->print_constraints_stats(1e-6);
        if(SDP->_status==0)
        {
            lower_bound_nonlin_init = SDP->get_obj_val();
            DebugOn("Initial lower bound = "<<lower_bound_nonlin_init<<endl);
        }
        status=SDP->_status;
        lower_bound=lower_bound_nonlin_init;
        total_time=time_end-time_start;
    }
    else{
        lb_solver_type=ipopt;
        auto nonlin_obj=false;
        current=true;
        std::vector<double> vrbasis;
        std::map<string,double> crbasis;
        lower_bound=0;
        double active_root_tol=1e-6;
        double lb_scale_value=1;
        solver<> UB_solver(OPF,ub_solver_type);
        UB_solver.run(output = 5, ub_solver_tol, 2000,"ma27", 2000);
        if(OPF->_status!=0){
            upper_bound=OPF->_obj->_range->second;
        }
        else{
            upper_bound=OPF->get_obj_val();
        }
        auto time_start=get_wall_time();
        auto SDP= build_SDPOPF(grid, current, nonlin_obj, sdp_kim);
        auto lin_model=SDP->buildOA();
        auto interior_model=lin_model->add_outer_app_solution(*SDP);
        auto close=SDP->root_refine(interior_model, lin_model, lb_solver_type, nb_root_refine, upper_bound, lower_bound, lb_scale_value, lb_solver_tol, active_root_tol, oacuts,  opt_abs_tol, opt_rel_tol, zero_tol, "ma27", 10000, 2000, vrbasis, crbasis, initialize_primal);
        total_time=get_wall_time()-time_start;
    }
    
    string result_name=string(prj_dir)+"/results_obbt/"+grid._name+".txt";
    auto final_gap = 100*(upper_bound - lower_bound)/std::abs(upper_bound+1e-6);
    auto gap_init=final_gap;
    result_name=string(prj_dir)+"/results_obbt/"+grid._name+".txt";
    DebugOn(grid._name<<" & "<<std::fixed<<std::setprecision(5)<<final_gap<<" & "<<std::setprecision(5)<<upper_bound*1e3<<" & "<<std::setprecision(5)<<lower_bound*1e3<<" & "
            <<std::setprecision(5)<<total_time<<endl);
    
    
    ofstream fout(result_name.c_str());
    fout<<grid._name<<" & "<<std::fixed<<std::setprecision(5)<<final_gap<<" & "<<std::setprecision(5)<<upper_bound*1e3<<" & "<<std::setprecision(5)<<lower_bound*1e3<<" & "
    <<std::setprecision(5)<<total_time<<endl;
    if(lower_bound==numeric_limits<double>::min()){
        fout<<"Lower bound not solved to optimality"<<endl;
    }
    fout.close();
    return 0;
}


