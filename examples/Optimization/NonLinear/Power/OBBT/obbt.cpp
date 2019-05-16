//
// Created by kbestuzheva on 12/11/17.
//
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

using namespace std;
using namespace gravity;


/* main */
int main (int argc, char * argv[]) {
    int output = 0;
    bool loss_from = true;
    size_t num_bags = 0;
    string num_bags_s = "100";
    string solver_str = "ipopt";
    string sdp_cuts_s = "yes";
    string loss_s = "yes";
    string lazy_s = "no";
    string orig_s = "no";
    bool lazy_bool = false;
    bool add_original=false;
    SolverType solv_type = ipopt;
    double tol = 1e-6;
    string mehrotra = "no";
    
    double start_tol=1E-3, interval_tol=1E-3;
    
    //string fname = string(prj_dir)+"/data_sets/Power/nesta_case5_pjm.m";
    
    string fname = string(prj_dir)+"/data_sets/Power/nesta_case5_pjm.m";
    
    // create a OptionParser with options
    op::OptionParser opt;
    opt.add_option("h", "help",
                   "shows option help");//No option means Default values which may be seen above using option strings
    opt.add_option("f", "file", "Input file name", fname);
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
        solv_type = mosek;
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
    
    cout << "\nnum bags = " << num_bags << endl;
    
    PowerNet grid;
    grid.readgrid(fname);
    grid.get_tree_decomp_bags(false,true);
    
    
    
    
    
    DebugOn("Machine has " << thread::hardware_concurrency() << " threads." << endl);
    
    int nb_threads = thread::hardware_concurrency();
    double upper_bound = grid.solve_acopf();
    auto SDP= build_SDPOPF(grid, loss_from, upper_bound);
    solver<> SDPLB(SDP,solv_type);
    SDPLB.run(output = 5, tol = 1e-6);
    double lower_bound=SDP->get_obj_val();
    SDP->print();
    vector<shared_ptr<Model<>>> batch_models;
    map<string, bool> fixed_point;
    map<string, double> interval_original, interval_new;
    string p, pk;
    string vname;
    string mname, mkname, vkname, keyk, dirk;
    string dir_array[2]={"LB", "UB"};
    var<> vark, vk, v;
    int iter=0;
    double boundk1, objk;
    bool terminate=false;
    bool infeasible=false;
    bool break_flag=false;
    const double upp_low_tol=1e-3, fixed_tol=0.001, max_time=100, zero_tol=1e-6, range_tol=1e-6;
    const int dec_round=6;
    double solver_time_end, solver_time, solver_time_start = get_wall_time();
    if ((upper_bound-lower_bound)/(upper_bound+zero_tol)>=upp_low_tol)
    {
        for(auto &it:SDP->_vars_name)
        {
            string vname=it.first;
            v=SDP->get_var<double>(vname);
            auto v_keys=v.get_keys();
            for(auto &key: *v_keys)
            {
                p=vname+"|"+ key;
                fixed_point[p]=false;
                interval_original[p]=v._ub->eval(v.get_keys_map()->at(key))-v._lb->eval(v.get_keys_map()->at(key));
            }
            
        }
        
//        if(SDP->_obj->is_linear())
//        {
//            Constraint<> obj_lower("obj_lower");
//            obj_lower = *(SDP->_obj);
//            SDP->add(obj_lower>=lower_bound);
//        }
        solver_time= get_wall_time()-solver_time_start;
        while(solver_time<=max_time && !terminate)
        {
            iter++;
            terminate=true;
            for(auto &it:SDP->_vars_name)
            {
                vname=it.first;
                v = SDP->get_var<double>(vname);
                auto v_keys=v.get_keys();
                for(auto &key: *(v.get_keys()))
                {
                    solver_time_end=get_wall_time();
                    solver_time= solver_time_end-solver_time_start;
                    if(solver_time>=max_time)
                    {
                        break_flag=true;
                        break;
                    }
                    p=vname+"|"+ key;
                    interval_new[p]=v._ub->eval(v.get_keys_map()->at(key))-v._lb->eval(v.get_keys_map()->at(key));
                    if(abs(v._ub->eval(v.get_keys_map()->at(key))-v._lb->eval(v.get_keys_map()->at(key)))<=range_tol)
                    {
                        fixed_point[p]=true;
                    }
                    if(fixed_point[p]==false)
                    {
                        for(auto &dir: dir_array)
                        {
                            auto modelk = SDP->copy();
                            modelk->reset_constrs();
                            mname=vname+"|"+key+"|"+dir;
                            modelk->set_name(mname);
                            
                            vark=modelk->get_var<double>(vname);
                            if(dir=="LB")
                            {
                                modelk->min(vark(key));
                            }
                            else
                            {
                                modelk->min(vark(key)*(-1));
                                
                            }
                            
                            batch_models.push_back(modelk);
                            if (batch_models.size()==nb_threads || (it==*SDP->_vars_name.end() && key==v.get_keys()->back() && dir=="UB"))
                            {
                                double batch_time_start = get_wall_time();
                                run_parallel(batch_models,ipopt,1e-6,nb_threads,"ma97");
                                double batch_time_end = get_wall_time();
                                auto batch_time = batch_time_end - batch_time_start;
                                DebugOn("Done running batch models, solve time = " << to_string(batch_time) << endl);
                                for (auto model:batch_models)
                                {
                                    //                                    model->print();
                                    mkname=model->get_name();
                                    std::size_t pos = mkname.find("|");
                                    vkname.assign(mkname, 0, pos);
                                    mkname=mkname.substr(pos+1);
                                    pos=mkname.find("|");
                                    keyk.assign(mkname, 0, pos);
                                    dirk=mkname.substr(pos+1);
                                    vk=SDP->get_var<double>(vkname);
                                    pk=vkname+"|"+keyk;
                                    if(model->_status==0)
                                    {
                                        objk=model->get_obj_val();
                                        if(dirk=="LB")
                                            boundk1=vk.get_lb(keyk);
                                        else
                                        {
                                            if(objk!=0)
                                            {
                                            objk*=-1;
                                            }
                                            boundk1=vk.get_ub(keyk);
                                        }
                                        if(abs((boundk1-objk)/(boundk1+zero_tol))<=fixed_tol)
                                        {
                                            fixed_point[pk]=true;
                                        }
                                        else
                                        {
                                            fixed_point[pk]=false;
                                            terminate=false;
                                            if(dirk=="LB")
                                                vk(keyk).set_lb(objk);
                                            else
                                                vk(keyk).set_ub(objk);
                                            if(vk.get_ub(keyk)<vk.get_lb(keyk))
                                            {
                                                fixed_point[pk]=true;
                                                double temp=vk.get_ub(keyk);
                                                vk(keyk).set_ub(vk.get_lb(keyk));
                                                vk(keyk).set_lb(temp);
//                                                int temp=vk.get_ub(keyk)*pow(10, dec_round);
//                                                vk(keyk).set_ub(temp/(1.0*pow(10, dec_round)));
//                                                vk(keyk).set_lb(temp/(1.0*pow(10, dec_round)));
                                            }
                                        }
                                    }
                                    else
                                    {
                                        infeasible=true;
                                        DebugOn("OBBT step has failed in iteration\t"<<iter<<endl);
                                        model->print();
                                        return(-1);
                                    }
                                    
                                    
                                }
                                batch_models.clear();
                            }
                        }
                    }
                }
            }
            if(break_flag==true)
            {
                DebugOn("Maximum Time Exceeded\t"<<max_time<<endl);
                SDP->print();
                break;
            }
            DebugOn("Iteration "<<iter<<" completed"<<endl);
            SDP->print();
        }
    }
    else
    {
        DebugOn("Root node gap is already zero without OBBT\n");
    }
    
    DebugOn("Terminate\t"<<terminate<<endl);
    DebugOn("Time\t"<<solver_time<<endl);
    DebugOn("Iterations\t"<<iter<<endl);
    DebugOn("Variable \t Key \t Interval reduction percentage"<<endl);
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
            DebugOn(p<<interval_gap.back()<<endl);
            
        }
        
    }
    avg=sum/num_var;
    
    DebugOn("Average interval reduction\t"<<avg<<endl);
    
    SDP->reset_constrs();
    SDP->print();
    solver<> SDPLB1(SDP,solv_type);

    SDPLB1.run(output = 5, tol = 1e-6, "ma97");

    SDP->print_solution();

    double gap = 100*(upper_bound - lower_bound)/upper_bound;
    DebugOn("Initial Gap = " << to_string(gap) << "%."<<endl);
    gap = 100*(upper_bound - SDP->get_obj_val())/upper_bound;
    
    DebugOn("Final Gap = " << to_string(gap) << "%."<<endl);
    DebugOn("Upper bound = " << to_string(upper_bound) << "."<<endl);
    DebugOn("Lower bound = " << to_string(SDP->get_obj_val()) << "."<<endl);
    DebugOn("\nResults: " << grid._name << " " << to_string(SDP->get_obj_val()) << " " <<endl);
    return 0;
    
}


