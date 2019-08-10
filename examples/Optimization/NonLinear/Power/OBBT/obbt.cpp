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
#include <optionParser.hpp>
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
    bool loss_from = false;
    size_t num_bags = 0;
    string num_bags_s = "100";
    string solver_str = "ipopt";
    string sdp_cuts_s = "yes";
    
    string loss_s = "yes";
    string time_s = "60";
    
    string lazy_s = "no";
    string orig_s = "no";
    bool lazy_bool = false;
    bool add_original=false;
    SolverType solv_type = ipopt;
    SolverType solv_type_c = cplex;
    double tol = 1e-6;
    string mehrotra = "no";
    
    
    string fname = string(prj_dir)+"/data_sets/Power/nesta_case5_pjm.m";
    
    // create a OptionParser with options
    op::OptionParser opt;
    opt.add_option("h", "help",
                   "shows option help");//No option means Default values which may be seen above using option strings
    opt.add_option("f", "file", "Input file name", fname);
    opt.add_option("t", "time", "Time limit, defaut 60 secs", time_s);
    opt.add_option("s", "solver", "Solvers: ipopt/cplex/gurobi, default = ipopt", solver_str);
    opt.add_option("b", "numbags", "Number of bags per iteration", num_bags_s);
    opt.add_option("l", "losses", "add loss constraints", loss_s); //Adds loss from of true and if true also adds loss_to in a lazy fasion
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
    
    loss_s = opt["l"];
    if (loss_s.compare("no")==0) {
        loss_from = false;
    }
    else {
        loss_from = true;
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
    
    
    
    cout << "\nnum bags = " << num_bags << endl;
    
    PowerNet grid;
    grid.readgrid(fname);
    grid.get_tree_decomp_bags();
   
    
    
    auto c1 = grid.c1.in(grid.gens);
    auto c2 = grid.c2.in(grid.gens);
    auto c0 = grid.c0.in(grid.gens);
    auto arcs = indices(grid.arcs);
 
    bool sdp_nonlin=true;
    
    DebugOn("Machine has " << thread::hardware_concurrency() << " threads." << endl);
    

    int nb_threads = thread::hardware_concurrency();
    int nb_total_threads = nb_threads; /** Used when MPI is ON to multipply with the number of workers */
#ifdef USE_MPI
    nb_total_threads *= nb_workers;
#endif
    //int nb_threads = thread::hardware_concurrency();
//    int nb_threads =12;
    
    auto OPF=build_ACOPF(grid, ACRECT);
    solver<> OPFUB(OPF, solv_type);
    OPFUB.run(output = 5, tol = 1e-6, "ma57");
    OPF->print_solution();
    double upper_bound=OPF->get_obj_val();
    auto SDPL= build_SDPOPF(grid, loss_from, upper_bound, false, true);
    auto SDP= build_SDPOPF(grid, loss_from, upper_bound, false, sdp_nonlin);
    solver<> SDPLB(SDPL,solv_type);
    DebugOn("Lower bounding ipopt"<<endl);
    SDPLB.run(output = 5, tol = 1e-6);
    SDPL->print();
    double lower_bound=SDPL->get_obj_val();
    const double active_tol=1e-6;
    
    auto gap = 100*(upper_bound - lower_bound)/upper_bound;
    DebugOn("Initial Gap ipopt = " << to_string(gap) << "%."<<endl);
    
//    auto SDP=SDPL;
//    SDP->print();
//    SDP->print_solution();
    
    
    if(!grid._tree && grid.add_3d_nlin && !sdp_nonlin)
    {

    auto SDPcon=SDPL->get_constraint("SDP_3D");

    for(auto i=0;i<SDPcon->get_nb_inst();i++)
        //  for(auto i=0;i<1;i++)
    {
        SDPcon->uneval();
        DebugOn("eval of con "<<SDPcon->eval(i)<<endl);
        SDPcon->uneval();

        if(abs(SDPcon->eval(i))<=active_tol)
        {
            SDPcon->uneval();
            func<> oas=SDPcon->get_outer_app_insti(i);
            oas.eval_all();
            Constraint<> OA_sol("OA_cuts_solution"+to_string(i));
            OA_sol=oas;
            SDP->add(OA_sol>=0);
            oas.uneval();
        }
    }
    }
    
    DebugOn("cplex"<<endl);
    solver<> SDPLB1(SDP,solv_type_c);
    DebugOn("Lower bounding cplex"<<endl);
    SDPLB1.run(output = 7, tol=1e-7);
    SDP->print();
    
    gap = 100*(upper_bound - SDP->get_obj_val())/upper_bound;
    DebugOn("Initial Gap cplex = " << to_string(gap) << "%."<<endl);
    
    
    auto pf_to_min=grid.pf_to_min.in(arcs);
    
    //    auto SDP= build_SDPOPF_QC(grid, loss_from, upper_bound, lower_bound);
    //    solver<> SDPLBI(SDP,solv_type);
    //    SDP->print();
    //    SDPLBI.run(output = 5, tol = 1e-6, "ma57");
    
    vector<shared_ptr<Model<>>> batch_models;
    map<string, bool> fixed_point;
    map<string, double> interval_original, interval_new, ub_original, lb_original;
    string p, pk;
    string vname;
    string mname, mkname, vkname, keyk, dirk;
    string dir_array[2]={"LB", "UB"};
    var<> vark, vk, v;
    int iter=0;
    double boundk1, objk, left, right, mid, temp, tempa;
    bool terminate=false;
    bool infeasible=false;
    
    bool break_flag=false, time_limit = false, lifted_var=false, close=false;
    
    const double upp_low_tol=1e-3, fixed_tol_abs=1e-3, fixed_tol_rel=1e-3, zero_tol=1e-6, range_tol=1e-3, zero_val=1e-6;
    const int max_iter=50,gap_count_int=6;
    
    
    double solver_time_end, solver_time =0, solver_time_start = get_wall_time();
    shared_ptr<map<string,size_t>> p_map;
    //Check if gap is already not zero at root node
    if (upper_bound-lower_bound>=upp_low_tol && (upper_bound-lower_bound)/(upper_bound+zero_tol)>=upp_low_tol)
        
    {
        terminate=false;
        for(auto &it:SDP->_vars_name)
        {
            string vname=it.first;
            v=SDP->get_var<double>(vname);
            auto v_keys=v.get_keys();
            auto v_key_map=v.get_keys_map();
//             if(vname=="Pf_to")
//             {
//                  p_map=pf_to_min.get_keys_map();
//             }
            for(auto &key: *v_keys)
            {
                p=vname+"|"+ key;
                //Do not do OBBT on lifted variables
                if(v._lift){
                    fixed_point[p]=true;
                }
                else{
                    fixed_point[p]=false;
                }
                auto key_pos=v_key_map->at(key);
                
                if(v._off[key_pos]==true)
                {
                    fixed_point[p]=true;
                    DebugOn("Skipping OBBT for "<<vname<<"\t"<<key<<endl);
                }
               
                if(vname=="Pf_to")
                {
                    
            }
//                auto key_pos=p_map->at(key);
//
//                if(pf_to_min._off[key_pos]==true)
//                {
//                   fixed_point[p]=true;
//                    DebugOn("Skipping OBBT for "<<vname<<"\t"<<key<<endl);
//                }
//                }

                interval_original[p]=v.get_ub(key)-v.get_lb(key);
                ub_original[p]=v.get_ub(key);
                lb_original[p]=v.get_lb(key);
                interval_new[p]=v.get_ub(key)-v.get_lb(key);
                
            }
            
        }
        
        solver_time= get_wall_time()-solver_time_start;
         auto v_in_cons=SDP->_v_in_cons;
        while(solver_time<=max_time && !terminate && iter<=max_iter)
        {
            iter++;
            terminate=true;
            for (auto it=SDP->_vars_name.begin(); it!=SDP->_vars_name.end(); it++)
            {
                vname=it->first;
                v = SDP->get_var<double>(vname);
                auto v_keys=v.get_keys();
                for(auto it_key=v.get_keys()->begin(); it_key!=v.get_keys()->end(); it_key++)
                {
                    
                    auto key = *it_key;
                    solver_time_end=get_wall_time();
                    solver_time= solver_time_end-solver_time_start;
                    if(solver_time>=max_time)
                        
                    {
                        break_flag=true;
                        time_limit = true;
                        break;
                    }
                    p=vname+"|"+ key;
                    interval_new[p]=v.get_ub(key)-v.get_lb(key);
                    if(abs(v.get_ub(key)-v.get_lb(key))<=range_tol)
                    {
                        fixed_point[p]=true;
                        
                    }
                    //Either if not fixed point, or if at the last key of the last variable
                    if(fixed_point[p]==false || (next(it)==SDP->_vars_name.end() && next(it_key)==v.get_keys()->end()))
                    {
                        //Loop on directions, upper bound and lower bound
                        for(auto &dir: dir_array)
                        {
                            auto modelk = SDP->copy();
                            mname=vname+"|"+key+"|"+dir;
                            modelk->set_name(mname);
                            
                            vark=modelk->get_var<double>(vname);
                            if(dir=="LB")
                            {
                                modelk->min(vark(key));
                            }
                            else
                            {
                                modelk->max(vark(key));
                                
                            }
                            
                            if(fixed_point[p]==false){
                                batch_models.push_back(modelk);
                            }
                            //When batch models has reached size of nb_threads or when at the last key of last avriable
                            if (batch_models.size()==nb_total_threads || (next(it)==SDP->_vars_name.end() && next(it_key)==v.get_keys()->end() && dir=="UB"))
                            {
                                double batch_time_start = get_wall_time();
#ifdef USE_MPI
                                run_MPI(batch_models,ipopt,1e-6,nb_threads, "ma57",true);
#else
                                run_parallel(batch_models,ipopt,1e-7,nb_threads, "ma57");
                               // run_parallel(batch_models,cplex,1e-7,nb_threads);
#endif
                                double batch_time_end = get_wall_time();
                                auto batch_time = batch_time_end - batch_time_start;
                                DebugOn("Done running batch models, solve time = " << to_string(batch_time) << endl);
                                for (auto model:batch_models)
                                {
                                    mkname=model->get_name();
                                    std::size_t pos = mkname.find("|");
                                    vkname.assign(mkname, 0, pos);
                                    mkname=mkname.substr(pos+1);
                                    pos=mkname.find("|");
                                    keyk.assign(mkname, 0, pos);
                                    dirk=mkname.substr(pos+1);
                                    vk=SDP->get_var<double>(vkname);
                                    pk=vkname+"|"+keyk;
                                    //Update bounds only of the model status is solved to optimal                                }
                                    if(model->_status==0)
                                    {
                                        objk=model->get_obj_val();
                                        if(dirk=="LB")
                                        {
                                            boundk1=vk.get_lb(keyk);
                                            //Uncertainty in objk=obk+-solver_tolerance, here we choose lowest possible value in uncertainty interval
                                            objk=std::max(objk-range_tol, boundk1);
                                        }
                                        else
                                        {
                                            boundk1=vk.get_ub(keyk);
                                            //Uncertainty in objk=obk+-solver_tolerance, here we choose highest possible value in uncertainty interval
                                            objk=std::min(objk+range_tol, boundk1);
                                            
                                        }
                                        if(abs(boundk1-objk) <= fixed_tol_abs || abs((boundk1-objk)/(boundk1+zero_tol))<=fixed_tol_rel)
                                        {//do not close intervals to OBBT before finishing at least one full iteration over all variables
                                            if(iter>1)
                                                fixed_point[pk]=true;
                                            
                                        }
                                        else
                                        {
                                            if(dirk=="LB"){
                                                vk.set_lb(keyk, objk);/* IN MPI this needs to be broadcasted back to the other workers */
                                            }
                                            else{
                                                vk.set_ub(keyk, objk);
                                            }
                                            //If crossover in bounds,just exchange them
                                            if(vk.get_ub(keyk)<vk.get_lb(keyk))
                                            {
                                                fixed_point[pk]=true;
                                                temp=vk.get_ub(keyk);
                                                tempa=vk.get_lb(keyk);
                                                vk.set_ub(keyk, tempa);
                                                vk.set_lb(keyk, temp);
                                                
                                            }
                                            else if(!vk._lift){
                                                fixed_point[pk]=false;
                                                terminate=false;
                                            }
                                            
                                        }
                                        //If interval becomes smaller than range_tol, reset bounds so that interval=range_tol
                                        if(abs(vk.get_ub(keyk)-vk.get_lb(keyk))<range_tol)
                                        {
                                            //If original interval is itself smaller than range_tol, do not have to reset interval
                                            if(interval_original[pk]>=range_tol)
                                            {
                                                DebugOn("Entered reset");
                                                //Mid is the midpoint of interval
                                                mid=(vk.get_ub(keyk)+vk.get_lb(keyk))/2.0;
                                                left=mid-range_tol/2.0;
                                                right=mid+range_tol/2.0;
                                                //If resized interval does not cross original bounds, reset
                                                if(right<=ub_original[pk] && left>=lb_original[pk])
                                                {
                                                    vk.set_ub(keyk, right);
                                                    vk.set_lb(keyk, left);
                                                }
                                                //If resized interval crosses original upperbound, set the new bound to upperbound, and lower bound is expanded to upperbound-range_tolerance
                                                else if(right>ub_original[pk])
                                                {
                                                    
                                                    vk.set_ub(keyk, ub_original[pk]);
                                                    vk.set_lb(keyk, ub_original[pk]-range_tol);
                                                }
                                                //If resized interval crosses original lowerbound, set the new bound to lowerbound, and upper bound is expanded to lowerbound+range_tolerance
                                                else if(left<lb_original[pk])
                                                {
                                                    vk.set_lb(keyk, lb_original[pk]);
                                                    vk.set_ub(keyk, lb_original[pk]+range_tol);
                                                    
                                                }
                                                //In the resized interval both original lower and upper bounds can not be crosses, because original interval is greater
                                                //than range_tol
                                                
                                            }
                                        }
                                    }
                                    else
                                    {
                                        DebugOn("OBBT step has failed in iteration\t"<<iter<<endl);
                                        //                                            model->print();
                                        
                                        //                                        fixed_point[pk]=true;
                                    }
                                }
                                batch_models.clear();
                            }
                        }
                    }
                }
            }
            
            //Check if OBBT has converged, can check every gap_count_int intervals
            if(iter%gap_count_int==0)
            {
                SDP->reset_constrs();
                solver<> SDPLB1(SDP,solv_type);
                SDPLB1.run(output = 5, tol=1e-8);
                if(SDP->_status==0)
                {
                    gap = 100*(upper_bound - (SDP->get_obj_val()))/upper_bound;
                    DebugOn("Gap "<<gap<<endl);
                }
                if (upper_bound-SDP->get_obj_val()<=upp_low_tol || (upper_bound-SDP->get_obj_val())/(upper_bound+zero_tol)<=upp_low_tol)
                {
                    DebugOn("Gap closed at iter "<< iter<<endl);
                    close=true;
                    terminate=true;
                    SDP->print();
                    
                }
            }
            
            if(break_flag==true)
            {
                DebugOn("Maximum Time Exceeded\t"<<max_time<<endl);
                DebugOn("Iterations\t"<<iter<<endl);
                
                break;
            }
            solver_time= get_wall_time()-solver_time_start;
            DebugOn("Solved Fixed Point iteration " << iter << endl);
        }
        vector<double> interval_gap;
        double sum=0, avg, num_var=0.0;
        for(auto &it:SDP->_vars_name)
        {
            string vname=it.first;
            v=SDP->get_var<double>(vname);
            auto v_keys=v.get_keys();
            for(auto &key: *v_keys)
            { num_var++;
                p=vname+"|"+ key;
                interval_gap.push_back((interval_original[p]-interval_new[p])/(interval_original[p]+zero_tol)*100.0);
                sum+=interval_gap.back();
                DebugOn(p<<" " << interval_gap.back()<< " flag = " << fixed_point[p] << endl);
            }
            
        }
        avg=sum/num_var;
        
        DebugOn("Average interval reduction\t"<<avg<<endl);
        
        if(!close)
        {
            
            SDP->reset_constrs();
            solver<> SDPLB1(SDP,solv_type);
            
            SDPLB1.run(output = 5, tol=1e-8);
        }
        
    }
//    avg=sum/num_var;
//
//    DebugOn("Average interval reduction\t"<<avg<<endl);
    
    if(!close)
    {
#ifdef USE_MPI
        if(worker_id==0){
#endif
            SDP->reset_constrs();
            solver<> SDPLB1(SDP,solv_type);
            
            SDPLB1.run(output = 5, tol=1e-6, "ma97");
            SDP->print_constraints_stats(tol);
            bool print_only_relaxed;
            SDP->print_nonzero_constraints(tol,print_only_relaxed=true);
            
            //        SDP->print_solution();
            
            //        SDP->print();
            
            if(SDP->_status==0)
            {
                
                DebugOn("\nResults: " << grid._name << " " << to_string(SDP->get_obj_val()) << " " <<endl);
                DebugOn("Solution Print"<<endl);
                SDP->print();
                //                SDP->print_solution();
                SDP->print_constraints_stats(tol);
                gap = 100*(upper_bound - lower_bound)/upper_bound;
                DebugOn("Initial Gap = " << to_string(gap) << "%."<<endl);
                gap = 100*(upper_bound - (SDP->get_obj_val()))/upper_bound;
                DebugOn("Final Gap = " << to_string(gap) << "%."<<endl);
                DebugOn("Upper bound = " << to_string(upper_bound) << "."<<endl);
                DebugOn("Lower bound = " << to_string((SDP->get_obj_val())) << "."<<endl);
                DebugOn("Time\t"<<solver_time<<endl);
                
            }
            else
            {
                double gap = 100*(upper_bound - lower_bound)/upper_bound;
                DebugOn("Initial Gap = " << to_string(gap) << "%."<<endl);
                DebugOn("Lower bounding problem status = " << SDP->_status <<endl);
                DebugOn("Lower bounding problem not solved to optimality, cannot compute final gap"<<endl);
            }
            if(time_limit){
                DebugOn("Reached Time limit!"<<endl);
            }
            else {
                DebugOn("Terminate\t"<<terminate<<endl);
            }
            
            
            DebugOn("Time\t"<<solver_time<<endl);
            DebugOn("Iterations\t"<<iter<<endl);
#ifdef USE_MPI
        }
#endif
    }
#ifdef USE_MPI
    MPI_Finalize();
#endif
    return 0;
}

