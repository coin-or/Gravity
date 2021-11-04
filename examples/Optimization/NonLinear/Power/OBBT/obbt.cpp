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

void initialize_relaxation(shared_ptr<Model<double>> OPF, shared_ptr<Model<double>> relax, PowerNet& grid);
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
    string nb_root_refine_s="10";
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
    bool linearize=false;
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
    //OPF->print();
    double ub_solver_tol=1e-6, lb_solver_tol=1e-8, range_tol=1e-3, opt_rel_tol=1e-2, opt_abs_tol=1e6;
    int total_iter;
    unsigned max_iter=1e3;
    int oacuts=0, oacuts_init=0, fail=0;
    SolverType ub_solver_type = ipopt, lb_solver_type = solv_type;
    bool scale_objective;
    bool termination=true;
    int status=0;
    //linearize=true;
    if(!linearize){
        auto nonlin_obj=true;
        scale_objective=true;
        current=false;
        double lb_scale_value=1.0;
        auto SDP= build_SDPOPF(grid, current, nonlin_obj, sdp_kim);
        //        SDP->print();
        solver<> UB_solver(OPF,ub_solver_type);
        UB_solver.run(output = 5, ub_solver_tol, 2000,"ma27", 2000);
        if(OPF->_status!=0){
            upper_bound=OPF->_obj->_range->second;
        }
        else{
            upper_bound=OPF->get_obj_val();
        }
       
        initialize_relaxation(OPF, SDP, grid);
        
//        OPF->print();
//        SDP->print();
//        SDP->reset();
        //SDP->print_solution();
        //SDP->print_constraints_stats(1e-6);
        solver<> LBnonlin_solver(SDP,lb_solver_type);
        if(scale_objective){
            auto obj = *SDP->_obj/1e5;
            SDP->min(obj);
            SDP->reset();
            lb_scale_value=1;
        }
        SDP->print();
        auto time_start=get_wall_time();
        LBnonlin_solver.run(output = 5 , 1e-6, 1e6, "ma57", 2000);
        auto time_end = get_wall_time();
        if(SDP->_status==0)
        {
            lower_bound_nonlin_init = SDP->get_obj_val()*lb_scale_value;
            DebugOn("Initial lower bound = "<<lower_bound_nonlin_init<<endl);
        }
        status=SDP->_status;
        lower_bound=lower_bound_nonlin_init;
        total_time=time_end-time_start;
    }
    string result_name=string(prj_dir)+"/results_obbt/"+grid._name+".txt";
    auto final_gap = 100*(upper_bound - lower_bound)/std::abs(upper_bound+1e-6);
    auto gap_init=final_gap;
    result_name=string(prj_dir)+"/results_obbt/"+grid._name+".txt";
#ifdef USE_MPI
    if(worker_id==0){
        ofstream fout(result_name.c_str());
        fout<<grid._name<<" & "<<std::fixed<<std::setprecision(5)<<final_gap<<" & "<<std::setprecision(5)<<upper_bound<<" & "<<std::setprecision(5)<<lower_bound<<" & "
        <<std::setprecision(5)<<total_time<<" & "
        <<status<<endl;
        if(lower_bound==numeric_limits<double>::min()){
            fout<<"Lower bound not solved to optimality"<<endl;
        }
        DebugOn("I am worker id "<<worker_id<<" writing to results file "<<endl);
        fout.close();
    }
    MPI_Finalize();
#else
    DebugOn(grid._name<<" & "<<std::fixed<<std::setprecision(5)<<final_gap<<" & "<<std::setprecision(5)<<upper_bound<<" & "<<std::setprecision(5)<<lower_bound<<" & "
            <<std::setprecision(5)<<total_time<<" & "
            <<status<<endl);
    
    
    ofstream fout(result_name.c_str());
    fout<<grid._name<<" & "<<std::fixed<<std::setprecision(5)<<final_gap<<" & "<<std::setprecision(5)<<upper_bound<<" & "<<std::setprecision(5)<<lower_bound<<" & "
    <<std::setprecision(5)<<total_time<<" & "
    <<status<<endl;
    if(lower_bound==numeric_limits<double>::min()){
        fout<<"Lower bound not solved to optimality"<<endl;
    }
    fout.close();
#endif
    return 0;
}
void initialize_relaxation(shared_ptr<Model<double>> OPF, shared_ptr<Model<double>> relax, PowerNet& grid){
    bool current=false;
    OPF->print_constraints_stats(1e-6);
    bool print_bags = false, only_3d_bags = false;
    auto bags_3d=grid.decompose_bags_3d_linear(print_bags, only_3d_bags);
    auto node_pairs = grid.get_node_pairs();
    auto node_pairs_chord = grid.get_node_pairs_chord(bags_3d);
    bool sdp_cuts=true;
    if (grid._tree || !grid.add_3d_nlin || !sdp_cuts) {
        node_pairs_chord = node_pairs;
    }
    auto nodes = indices(grid.nodes);
    auto arcs = indices(grid.arcs);
    auto gens = indices(grid.gens);
    auto b = grid.b.in(arcs);
    auto g = grid.g.in(arcs);
    
    auto vr = OPF->get_var<double>("vr");
    auto vi = OPF->get_var<double>("vi");
    auto pg=  OPF->get_var<double>("Pg");
    auto qg=  OPF->get_var<double>("Qg");
    auto p_from=  OPF->get_var<double>("Pf_from");
    auto q_from=  OPF->get_var<double>("Qf_from");
    auto p_to=  OPF->get_var<double>("Pf_to");
    auto q_to=  OPF->get_var<double>("Qf_to");
   
    auto R_Wij = relax->get_ptr_var<double>("R_Wij");
    auto Im_Wij = relax->get_ptr_var<double>("Im_Wij");
    auto Wii = relax->get_ptr_var<double>("Wii");
    auto Pg=  relax->get_ptr_var<double>("Pg");
    auto Qg=  relax->get_ptr_var<double>("Qg");
    auto Pf_from=  relax->get_ptr_var<double>("Pf_from");
    auto Qf_from=  relax->get_ptr_var<double>("Qf_from");
    auto Pf_to=  relax->get_ptr_var<double>("Pf_to");
    auto Qf_to=  relax->get_ptr_var<double>("Qf_to");
    shared_ptr<var<double>> lij, lji;
    if(current){
        lij=relax->get_ptr_var<double>("lij");
        lji=relax->get_ptr_var<double>("lji");
    }
    
    int count=0;
    for(auto key: *nodes._keys){
        auto value=(pow(vr.eval(key),2)+pow(vi.eval(key),2));
        Wii->_val->at(count)=value;
        count++;
    }
    count=0;
    for(auto key: *node_pairs_chord._keys){
        auto from_key=key.substr(0, key.find_first_of(","));
        auto to_key=key.substr(key.find_first_of(",")+1);
        auto vr_i=vr.eval(from_key);
        auto vi_i=vi.eval(from_key);
        auto vr_j=vr.eval(to_key);
        auto vi_j=vi.eval(to_key);
        R_Wij->_val->at(count)=vr_i*vr_j+vi_i*vi_j;
        Im_Wij->_val->at(count)=vi_i*vr_j-vi_j*vr_i;
        count++;
    }
    count=0;
    for(auto key: *gens._keys){
        Pg->_val->at(count)=pg.eval(key);
        Qg->_val->at(count)=qg.eval(key);
        count++;
    }
    count=0;
    for(auto key: *arcs._keys){
        Pf_from->_val->at(count)=p_from.eval(key);
        Qf_from->_val->at(count)=q_from.eval(key);
        Pf_to->_val->at(count)=p_to.eval(key);
        Qf_to->_val->at(count)=q_to.eval(key);
        count++;
    }
    if(current){
    count=0;
    for(auto key: *arcs._keys){
        auto a_key=key.substr(key.find_first_of(",")+1);
        auto from_key=a_key.substr(0, a_key.find_first_of(",")+1);
        auto to_key=a_key.substr(a_key.find_first_of(",")+1);
        auto vr_i=vr.eval(from_key);
        auto vi_i=vi.eval(from_key);
        auto vr_j=vr.eval(to_key);
        auto vi_j=vi.eval(to_key);
        auto rwij=vr_i*vr_j+vi_i*vi_j;
        auto l=(pow(g.eval(key),2)+pow(b.eval(key),2))*(vi_i*vi_i+vr_i*vr_i+vi_j*vi_j+vr_j*vr_j-2*rwij);
        lij->_val->at(count)=l;
        lji->_val->at(count)=l;
        count++;
    }
    }
}


