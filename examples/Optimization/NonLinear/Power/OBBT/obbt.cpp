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
    string time_s = "100000";
    string threads_s="12";
    
    string lazy_s = "no";
    string orig_s = "no";
    bool lazy_bool = false;
    bool add_original=false;
    SolverType solv_type = ipopt;
    const double tol = 1e-6;
    string mehrotra = "no";
    
    bool nonlin=true;
    
    
    string fname = string(prj_dir)+"/data_sets/Power/nesta_case5_pjm.m";

#ifdef USE_OPT_PARSER
    
    //create a OptionParser with options
    op::OptionParser opt;
    opt.add_option("h", "help",
                   "shows option help");//No option means Default values which may be seen above using option strings
    opt.add_option("f", "file", "Input file name", fname);
    opt.add_option("t", "time", "Time limit, defaut 60 secs", time_s);
    opt.add_option("threads", "threads", "Number of threads, defaut 24", threads_s);
    opt.add_option("s", "solver", "Solvers: ipopt/cplex/gurobi, default = ipopt", solver_str);
    opt.add_option("b", "numbags", "Number of bags per iteration", num_bags_s);
    opt.add_option("l", "", "add current constraints", current_s); //Adds loss from of true and if true also adds loss_to in a lazy fasion
    opt.add_option("lz", "lazy", "Generate 3d SDP cuts in a lazy fashion, default = no", lazy_s);
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


    //    double max_time = 40;

    auto max_time = op::str2double(opt["t"]);
#else
    if(argc==3){
    fname=argv[1];
    time_s=argv[2];
    }
    else{
        fname=string(prj_dir)+"/data_sets/Power/nesta_case5_pjm.m";
        time_s="3600";
    }
    current=true;
    
    auto max_time=std::stod(time_s);
    
    
#endif
    
    auto nb_threads=std::stod(threads_s);
    
    cout << "\nnum bags = " << num_bags << endl;
    
    PowerNet grid;
    grid.readgrid(fname);
    grid.get_tree_decomp_bags();
    
    auto nodes = indices(grid.nodes);
    auto arcs = indices(grid.arcs);
    auto v_max = grid.v_max.in(nodes);
    auto bags_3d=grid.decompose_bags_3d();
    auto bus_pairs = grid.get_bus_pairs();
    auto bus_pairs_chord = grid.get_bus_pairs_chord(bags_3d);
    
    auto c1 = grid.c1.in(grid.gens);
    auto c2 = grid.c2.in(grid.gens);
    auto c0 = grid.c0.in(grid.gens);
    
 

    
    DebugOn("Machine has " << thread::hardware_concurrency() << " threads." << endl);
    

   // int nb_threads = thread::hardware_concurrency();
    nb_threads=12;
    int nb_total_threads = nb_threads; /** Used when MPI is ON to multipply with the number of workers */
#ifdef USE_MPI
    nb_total_threads *= nb_workers;
#endif

    bool lb_solv;
    double gap=999, gapnl=999;
    double lower_bound=-99999, avg=0;
    double solver_time =0;
    int iter=0;
    unsigned max_iter=1000;
    unsigned precision=0;
    
    bool terminate=false, xb_true=false;
    
    auto OPF=build_ACOPF(grid, ACRECT);
    solver<> OPFUB(OPF, solv_type);
    OPFUB.run(output = 0, 1e-6);
//    OPF->print_solution();
    double upper_bound=OPF->get_obj_val();
    
    bool nonlin_obj=false;
   
  //  auto SDP= build_SDPOPF(grid, current, upper_bound, nonlin_obj);
    shared_ptr<Model<double>> SDPO;
    
    double lower_bound_init;
    
//    SDP->print();

  //  solver<> SDPLB(SDP,solv_type);
    //DebugOn("Lower bounding ipopt"<<endl);
   // double solver_time_start=get_wall_time();
   //SDPLB.run(output = 0, tol);
   // double solver_time_end=get_wall_time();
    //double solver_time_lb=solver_time_end-solver_time_start;
  //  SDP->print();
//    SDP->print_solution();
    
  
    
        std::pair<bool,double> ub;
        ub.first=true;
        ub.second=upper_bound;
    nonlin=false;
    if(nonlin){
        //SDPO=SDP->copy();
        
        nonlin_obj=false;
        
        auto SDP= build_SDPOPF(grid, current, upper_bound, nonlin_obj);
//        
//        auto SDPA= build_SDPOPF(grid, current, upper_bound, nonlin_obj);
//
//        auto SDP=SDPA->copy();
        
        SDP->print();
        
        auto res=SDP->run_obbt(max_time, max_iter, ub, precision, OPF, SDP, nonlin);
        if(SDP->_status==0)
        {
            SDP->print();
            
            lower_bound=SDP->get_obj_val()*upper_bound;
            gap=100*(upper_bound - lower_bound)/upper_bound;
            
            terminate=std::get<0>(res);
            iter=std::get<1>(res);
            solver_time=std::get<2>(res);
            lower_bound_init=std::get<3>(res);
            avg=std::get<4>(res);
            xb_true=std::get<5>(res);
            
            
            gapnl = 100*(upper_bound - lower_bound_init)/upper_bound;
            DebugOn("Initial Gap nonlinear = " << to_string(gapnl) << "%."<<endl);
            
            auto SDPOA=SDP->buildOA(15, 10);
              solver<> SDPLB(SDPOA,ipopt);
              SDPLB.run(output = 0, tol);
            
            
            auto gapl = 100*(upper_bound - SDPOA->get_obj_val()*upper_bound)/upper_bound;
            DebugOn("Gap linear at solution of OBBT model = " << to_string(gapl) << "%."<<endl);
            

            
        }
    }
    else{
        
        nonlin_obj=false;
        
        auto SDP= build_SDPOPF(grid, current, upper_bound, nonlin_obj);
//        vector<double> x_sol(SDP->get_nb_vars());
//        SDP->get_solution(x_sol);
         solver<> SDPLB(SDP, ipopt);
//        SDP->print();
        SDPLB.run(output = 0, 1e-6, "ma27");
//        SDP->print_solution();
      //  DebugOn("Objective = " << to_string_with_precision(SDP->get_obj_val(),20) << endl);
        DebugOn("solution of LB"<<endl);
        lb_solv=false;
        if(SDP->_status==0){
            lb_solv=true;
        lower_bound=SDP->get_obj_val()*upper_bound;
        
        gap=100*(upper_bound - lower_bound)/upper_bound;
        DebugOn("Gap "<<gap);
        
         SDPO=SDP->buildOA(10, 10);
       // SDPO->set_solution(x_sol);
     //   SDPO->print();
       //         DebugOn("stats OA-LB"<<endl);
      //  SDPO->print_constraints_stats(1e-10);
      //  SDPO->print();
//          solver<> SDPLin(SDPO, ipopt);
//        SDPLin.run(output = 0, 1e-8);
////        DebugOn("N vars "<<SDPO->_nb_vars<<endl);
////        DebugOn("N cons "<<SDPO->_nb_cons<<endl);
////        DebugOn("SDPO obj value\t"<<SDPO->get_obj_val()<<endl);
//        double gap_lin=100*(upper_bound - SDPO->get_obj_val()*upper_bound)/upper_bound;
//        DebugOn("Gap Linear"<<gap_lin);
        
        
     
        
       
        
        auto res=SDPO->run_obbt(max_time, max_iter, ub, precision, OPF, SDP, nonlin);
        
//        auto SDPO_IIS1=SDPO->build_model_IIS();
//        solver<> IIS_test1(SDPO_IIS1,cplex);
//        IIS_test1.run(output = 0, tol);
     //   SDPO_IIS1->print();
        
     //   SDPO_IIS1->print_solution();
        
        
        
//        solver<> test2(SDPO, cplex);
//        test2.run(output = 5, tol);
//        SDPO->print();
 
        
        if(SDPO->_status==0)
        {
            
            
            
   
            lower_bound=SDPO->get_obj_val()*upper_bound;
            gap=100*(upper_bound - lower_bound)/upper_bound;
            
            terminate=std::get<0>(res);
            iter=std::get<1>(res);
            solver_time=std::get<2>(res);
            lower_bound_init=std::get<3>(res);
            avg=std::get<4>(res);
            xb_true=std::get<5>(res);

            
            gapnl = 100*(upper_bound - lower_bound_init)/upper_bound;
           // DebugOn("Initial Gap= " << to_string(gapnl) << "%."<<endl);
           
        }
        
        
    }
    
    }
    
    

    string result_name=string(prj_dir)+"/results_obbt/"+grid._name+".txt";
#ifdef USE_MPI
    if(worker_id==0){

        
       
        
    	ofstream fout(result_name.c_str());
        fout<<grid._name<<"\t"<<std::fixed<<std::setprecision(5)<<gapnl<<"\t"<<std::setprecision(5)<<upper_bound<<"\t"<<std::setprecision(5)<<lower_bound<<"\t"<<std::setprecision(5)<<gap<<"\t"<<terminate<<"\t"<<iter<<"\t"<<std::setprecision(5)<<solver_time<<"\t"<<std::setprecision(5)<<avg<<"\t"<<xb_true<<"\t"<<lb_solv<<endl;
        DebugOn("I am worker id "<<worker_id<<" writing to results file "<<endl);
        fout.close();
     }
    MPI_Finalize();
#else
    DebugOn(grid._name<<"\t"<<std::fixed<<std::setprecision(5)<<gapnl<<"\t"<<std::setprecision(5)<<upper_bound<<"\t"<<std::setprecision(5)<<lower_bound<<"\t"<<std::setprecision(5)<<gap<<"\t"<<terminate<<"\t"<<iter<<"\t"<<std::setprecision(5)<<solver_time<<"\t"<<std::setprecision(5)<<avg<<"\t"<<xb_true<<endl);

    
	ofstream fout(result_name.c_str());
        fout<<grid._name<<"\t"<<std::fixed<<std::setprecision(5)<<gapnl<<"\t"<<std::setprecision(5)<<upper_bound<<"\t"<<std::setprecision(5)<<lower_bound<<"\t"<<std::setprecision(5)<<gap<<"\t"<<terminate<<"\t"<<iter<<"\t"<<std::setprecision(5)<<solver_time<<"\t"<<std::setprecision(5)<<avg<<"\t"<<xb_true<<endl;
        fout.close();
#endif
    
    return 0;
}

