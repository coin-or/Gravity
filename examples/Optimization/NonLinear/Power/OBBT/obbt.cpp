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

double oaa1, oab1;
/* main */
int main (int argc, char * argv[]) {
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

    
    
    
    DebugOn("Machine has " << thread::hardware_concurrency() << " threads." << endl);
    
    //int nb_threads = thread::hardware_concurrency()-1;
    int nb_threads =12;
    
    auto OPF=build_ACOPF(grid, ACRECT);
    solver<> OPFUB(OPF, solv_type);
    OPFUB.run(output = 5, tol = 1e-6, "ma57");
    
    //double upper_bound = grid.solve_acopf();
    //solve_acopf
    //    OPF->print_constraints_stats(tol);
    //    auto rv= OPF->get_var<double>("vr");
    //    DebugOn("W_orig\n");
    //    vector<double> rvv, ivv, w_orig, rrwij, iiwij;
    //
    //    for(auto &k:*rv.get_keys())
    //    {
    //        rvv.push_back(rv.eval(k));
    //        ivv.push_back(OPF->get_var<double>("vi").eval(k));
    //        w_orig.push_back(pow(rvv.back(),2)+pow(ivv.back(),2));
    //        DebugOn(w_orig.back()<<endl);
    //    }
    //    DebugOn("RRwij\tIIwij"<<endl);
    //
    //    for (auto &k:*(grid.bus_pairs._keys))
    //    {
    //        auto k1=k.substr(0, k.find_first_of(","));
    //        auto k2=k.substr(k.find_first_of(",")+1);
    //        auto rva=OPF->get_var<double>("vr").eval(k1);
    //        auto rvb=OPF->get_var<double>("vr").eval(k2);
    //        auto iva=OPF->get_var<double>("vi").eval(k1);
    //        auto ivb=OPF->get_var<double>("vi").eval(k2);
    //
    //        rrwij.push_back(rva*rvb+iva*ivb);
    //        iiwij.push_back(iva*rvb-ivb*rva);
    //        DebugOn(k<<"\t"<<rrwij.back()<<"\t"<<iiwij.back()<<endl);
    //
    //    }
    
    double upper_bound=OPF->get_obj_val();
    auto SDPL= build_SDPOPF(grid, loss_from, upper_bound);
    solver<> SDPLB(SDPL,solv_type);
    SDPLB.run(output = 5, tol = 1e-6, "ma97");
    double lower_bound=SDPL->get_obj_val();
    
    auto SDP= build_SDPOPF_QC(grid, loss_from, upper_bound, lower_bound);
    solver<> SDPLBI(SDP,solv_type);
        SDP->print();
    SDPLBI.run(output = 5, tol = 1e-6, "ma97");
   //auto SDP=SDPL;
    
 //   double lower_bound=SDP->get_obj_val();
    SDP->print_constraints_stats(tol);
    bool print_only_relaxed;
    SDP->print_nonzero_constraints(tol,print_only_relaxed=true);
    double gap = 100*(upper_bound - lower_bound)/upper_bound;
    DebugOn("Initial Gap = " << to_string(gap) << "%."<<endl);
    //    SDP->print_solution();
    vector<shared_ptr<Model<>>> batch_models;
    map<string, bool> fixed_point;
    map<string, double> interval_original, interval_new, ub_original, lb_original;
    string p, pk;
    string vname;
    string mname, mkname, vkname, keyk, dirk;
    string dir_array[2]={"LB", "UB"};
    var<> vark, vk, v;
    int iter=0;
    double boundk1, objk;
    bool terminate=false;
    bool infeasible=false;
    
    bool break_flag=false, time_limit = false, lifted_var=false;
    
    const double upp_low_tol=1e-3, fixed_tol_abs=1e-3, fixed_tol_rel=1e-3, zero_tol=1e-6, range_tol=1e-3, zero_val=1e-6;
    
    double solver_time_end, solver_time =0, solver_time_start = get_wall_time();
    if (upper_bound-lower_bound>=upp_low_tol && (upper_bound-lower_bound)/(upper_bound+zero_tol)>=upp_low_tol)
        
    {
        
        for(auto i = 0; i<1 ;i++){
            terminate=false;
            for(auto &it:SDP->_vars_name)
            {
                string vname=it.first;
                v=SDP->get_var<double>(vname);
                auto v_keys=v.get_keys();
                for(auto &key: *v_keys)
                {
                    p=vname+"|"+ key;
                    fixed_point[p]=false;
                    interval_original[p]=v.get_ub(key)-v.get_lb(key);
                    ub_original[p]=v.get_ub(key);
                    lb_original[p]=v.get_lb(key);
                    //                interval_new[p]=v._ub->eval(v.get_keys_map()->at(key))-v._lb->eval(v.get_keys_map()->at(key));
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
                for (auto it=SDP->_vars_name.begin(); it!=SDP->_vars_name.end(); it++)
                {
                    vname=it->first;
                    v = SDP->get_var<double>(vname);
                    lifted_var=v._lift;
                    if(!lifted_var)
                    {
                        if(vname!="objt")
                        {
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
                            if(fixed_point[p]==false || (next(it)==SDP->_vars_name.end() && next(it_key)==v.get_keys()->end()))
                            {
                                for(auto &dir: dir_array)
                                {
                                    auto modelk = SDP->copy();
//                                    auto Pg=modelk->get_var<double>("Pg");
//                                    Constraint<> obj_UB("obj_UB");
//                                    obj_UB=(product(c1,Pg) + product(c2,pow(Pg,2)) + sum(c0))-upper_bound;
//                                    modelk->add(obj_UB<=0);
//                                    func<double> o=*(SDP->_obj);
//                                    Constraint<> UpperB("UpperB");
//                                    UpperB=o;
//                                    modelk->add(UpperB<=upper_bound);
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
                                    if (batch_models.size()==nb_threads || (next(it)==SDP->_vars_name.end() && next(it_key)==v.get_keys()->end() && dir=="UB"))
                                    {
                                        double batch_time_start = get_wall_time();
                                        run_parallel(batch_models,ipopt,1e-6,nb_threads, "ma57");
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
                                            //                                    if(model->get_name()=="Lift_Im_Vi.from.in(bus_pairs_chordal)_Im_Vi.to.in(bus_pairs_chordal).in(bus_pairs_chordal)|1,2|UB")
                                            //                                        {
                                            //                                            model->print();
                                            //
                                            //
                                            //                                        }
                                            if(model->_status==0)
                                            {
                                                objk=model->get_obj_val();
                                                //                                            if(abs(objk)<=zero_val)
                                                //                                                objk=0.0;
                                                if(dirk=="LB")
                                                {
                                                    boundk1=vk.get_lb(keyk);
                                                    objk=std::max(objk-range_tol, boundk1);
                                                }
                                                else
                                                {
                                                    //                                                if(abs(objk)>=zero_val)
                                                    objk*=-1;
                                                    boundk1=vk.get_ub(keyk);
                                                    objk=std::min(objk+range_tol, boundk1);
                                                    
                                                }
                                                if(abs(boundk1-objk) <= fixed_tol_abs || abs((boundk1-objk)/(boundk1+zero_tol))<=fixed_tol_rel)
                                                {
                                                    if(iter>1)
                                                    fixed_point[pk]=true;
                                                    
                                                }
                                                else
                                                {
                                                    if(dirk=="LB"){
                                                        vk.set_lb(keyk, objk);
                                                    }
                                                    else{
                                                        vk.set_ub(keyk, objk);
                                                    }
                                                    if(vk.get_ub(keyk)<vk.get_lb(keyk))
                                                    {
                                                        fixed_point[pk]=true;
                                                        double temp=vk.get_ub(keyk);
                                                        double tempa=vk.get_lb(keyk);
                                                        vk.set_ub(keyk, tempa);
                                                        vk.set_lb(keyk, temp);
                                                        
                                                    }
                                                    else {
                                                        fixed_point[pk]=false;
                                                        terminate=false;
                                                    }
                                                    
                                                }
                                                if(abs(vk.get_ub(keyk)-vk.get_lb(keyk))<range_tol)
                                                {
//                                                    if(interval_original[pk]>=range_tol && !(abs(vk.get_ub(keyk))<=zero_val && abs(vk.get_lb(keyk))<=zero_val))
                                                                 if(interval_original[pk]>=range_tol)
                                                    {
                                                        DebugOn("Entered reset");
                                                        double mid=(vk.get_ub(keyk)+vk.get_lb(keyk))/2.0;
                                                        
                                                        double left=mid-range_tol/2.0;
                                                        double right=mid+range_tol/2.0;
                                                        DebugOn("UbO"<<ub_original[pk]<<endl);
                                                        DebugOn("LbO"<<lb_original[pk]<<endl);
                                                        if(right<=ub_original[pk] && left>=lb_original[pk])
                                                        {
                                                            vk.set_ub(keyk, right);
                                                            vk.set_lb(keyk, left);
                                                            //  DebugOn("Entered if 1"<<endl);
                                                        }
                                                        else if(right>ub_original[pk])
                                                        {
                                                            
                                                            vk.set_ub(keyk, ub_original[pk]);
                                                            vk.set_lb(keyk, ub_original[pk]-range_tol);
                                                            //  DebugOn("Entered if 2"<<endl);
                                                        }
                                                        else if(left<lb_original[pk])
                                                        {
                                                            vk.set_lb(keyk, lb_original[pk]);
                                                            vk.set_ub(keyk, lb_original[pk]+range_tol);
                                                            //  DebugOn("Entered if 3"<<endl);
                                                            
                                                        }
                                                        
                                                    
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
                        }
                    else
                    {
                        DebugOn("Did not do OBBT for \t" << vname);
                    }
                        
                }
                if(break_flag==true)
                {
                    DebugOn("Maximum Time Exceeded\t"<<max_time<<endl);
                    DebugOn("Iterations\t"<<iter<<endl);
                    //                    SDP->print();
                    break;
                }
                solver_time= get_wall_time()-solver_time_start;
                DebugOn("Solved Fixed Point iteration " << iter << endl);
                //            SDP->reset_constrs();
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
            //SDP->reset_constrs();
            //            auto SDP1=SDP->copy();
            //            //    SDP->print();
            //            solver<> SDPLB1(SDP1,solv_type);
            //
            //            SDPLB1.run(output = 5, tol=1e-6, "ma27");
            //
            //            SDP1->print_solution();
            //            SDP1->print();
            //            if(SDP1->_status==0)
            //            {
            //                double gap = 100*(upper_bound - 1e4*lower_bound)/upper_bound;
            //                DebugOn("Initial Gap = " << to_string(gap) << "%."<<endl);
            //                gap = 100*(upper_bound - 1e4*(SDP1->get_obj_val()))/upper_bound;
            //
            //                DebugOn("Final Gap = " << to_string(gap) << "%."<<endl);
            //                DebugOn("Upper bound = " << to_string(upper_bound) << "."<<endl);
            //                DebugOn("Lower bound = " << to_string(1e4*(SDP1->get_obj_val())) << "."<<endl);
            //                DebugOn("Time\t"<<solver_time<<endl);
            //                DebugOn("\nResults: " << grid._name << " " << to_string(SDP1->get_obj_val()) << " " <<endl);
            //            }
            //            else
            //            {
            //                double gap = 100*(upper_bound - 1e4*lower_bound)/upper_bound;
            //                DebugOn("Initial Gap = " << to_string(gap) << "%."<<endl);
            //                DebugOn("Lower bounding problem status = " << SDP1->_status <<endl);
            //                DebugOn("Lower bounding problem not solved to optimality, cannot compute final gap"<<endl);
            //            }
            
            
            //    SDP->print();
            // SDP->_built=false;
            SDP->reset_constrs();
            solver<> SDPLB1(SDP,solv_type);
            
            SDPLB1.run(output = 5, tol=1e-8);
            
            SDP->print_constraints_stats(tol);
            bool print_only_relaxed;
            SDP->print_nonzero_constraints(tol,print_only_relaxed=true);
            
                       SDP->print_solution();
            
            SDP->print();
            if(SDP->_status==0)
            {
                
                
                
                DebugOn("\nResults: " << grid._name << " " << to_string(SDP->get_obj_val()) << " " <<endl);
                DebugOn("Solution Print"<<endl);
                //                SDP->print_solution();
                SDP->print_constraints_stats(tol);
                double gap = 100*(upper_bound - lower_bound)/upper_bound;
                DebugOn("Initial Gap = " << to_string(gap) << "%."<<endl);
                gap = 100*(upper_bound - (SDP->get_obj_val()))/upper_bound;
                DebugOn("Final Gap = " << to_string(gap) << "%."<<endl);
                DebugOn("Upper bound = " << to_string(upper_bound) << "."<<endl);
                DebugOn("Lower bound = " << to_string((SDP->get_obj_val())) << "."<<endl);
                DebugOn("Time\t"<<solver_time<<endl);
                
                //                DebugOn("Var_map");
                //                for(auto it:SDP->_vars_name)
                //                    DebugOn(it.first<<endl);
                //
                //                DebugOn("SOC"<<endl);
                //                auto c=SDP->get_constraint("SOC");
                //                c->print();
                //                for (auto k=0;k<c->get_nb_instances();k++)
                //                DebugOn(k<<"\t"<<c->eval(k)<<endl);
                //
                //
                //
                //
                //                DebugOn("RealRankType1"<<endl);
                //                c=SDP->get_constraint("Real(RankType1)_lifted");
                //                //c->print();
                //                for (auto k=0;k<c->get_nb_instances();k++)
                //                    DebugOn(c->eval(k)<<endl);
                //
                ////                DebugOn("SDP"<<endl);
                ////                auto c1=SDP->get_constraint("SDP_3D");
                ////                for (auto k=0;k<c1->get_nb_instances();k++)
                ////                    DebugOn(c1->eval(k)<<endl);
                //
                //                 vector<double> rva,rvb,imva,imvb,rwab,imwab, rwij_gap, imwij_gap;
                //                DebugOn("RWij_gap and Im_Wij_gap"<<endl);
                //                auto vw=SDP->get_var<double>("Im_Wij");
                //                for (auto &k:*vw.get_keys())
                //                {
                //                    imwab.push_back(vw.eval(k));
                //                    rwab.push_back(SDP->get_var<double>("R_Wij").eval(k));
                //                    auto k1=k.substr(0, k.find_first_of(","));
                //                    auto k2=k.substr(k.find_first_of(",")+1);
                //                    rva.push_back(SDP->get_var<double>("R_Vi").eval(k1));
                //                    rvb.push_back(SDP->get_var<double>("R_Vi").eval(k2));
                //                    imva.push_back(SDP->get_var<double>("Im_Vi").eval(k1));
                //                    imvb.push_back(SDP->get_var<double>("Im_Vi").eval(k2));
                //                    rwij_gap.push_back(rwab.back()-rva.back()*rvb.back()-imva.back()*imvb.back());
                //                    imwij_gap.push_back(imwab.back()-imva.back()*rvb.back()+imvb.back()*rva.back());
                //                     DebugOn(k<<"\t"<<rwij_gap.back()<<"\t"<<imwij_gap.back()<<endl);
                //
                //                }
                //                DebugOn("R_Wij"<<endl);
                //                auto vr=SDP->get_var<double>("R_Wij");
                //                for (auto &k:*vr.get_keys())
                //                DebugOn(k<<"\t"<<vr.eval(k)<<endl);
                //
                //
                //
                //                vector<double> pf,pt;
                //                DebugOn("Pf_from"<<endl);
                //                vr=SDP->get_var<double>("Pf_from");
                //                for (auto &k:*vr.get_keys())
                //                {
                //                    DebugOn(k<<"\t"<<vr.eval(k)<<endl);
                //                    pf.push_back(vr.eval(k));
                //                }
                //
                //                DebugOn("Pf_to"<<endl);
                //                vr=SDP->get_var<double>("Pf_to");
                //                for (auto &k:*vr.get_keys())
                //                {
                //                    DebugOn(k<<"\t"<<vr.eval(k)<<endl);
                //                      pt.push_back(vr.eval(k));
                //                }
                //                vector<double> qf,qt;
                //                DebugOn("Qf_from"<<endl);
                //                vr=SDP->get_var<double>("Qf_from");
                //                vector<double> b_k;
                //                  vector<double> rvi,rvj,imvi,imvj;
                //                for (auto &k:*vr.get_keys())
                //                {
                //                    DebugOn(k<<"\t"<<vr.eval(k)<<endl);
                //                    qf.push_back(vr.eval(k));
                //                     b_k.push_back(grid.b.eval(k));
                //                    auto k1=k.substr(k.find_first_of(",")+1);
                //                    auto k2=k1.substr(0, k1.find_first_of(","));
                //                    auto k3=k1.substr(k1.find_first_of(",")+1);
                //                    rvi.push_back(SDP->get_var<double>("R_Vi").eval(k2));
                //                    rvj.push_back(SDP->get_var<double>("R_Vi").eval(k3));
                //                    imvi.push_back(SDP->get_var<double>("Im_Vi").eval(k2));
                //                    imvj.push_back(SDP->get_var<double>("Im_Vi").eval(k3));
                //                }
                //
                //                vector<double> Wiii, rviii, imviii;
                //             //   DebugOn("Wii.in(Nodes)"<<endl);
                //                          DebugOn("Wii-rv^2-iv^2"<<endl);
                //                vr=SDP->get_var<double>("Wii");
                //                for (auto &k:*vr.get_keys())
                //                {
                //                        Wiii.push_back(vr.eval(k));
                //                    rviii.push_back(SDP->get_var<double>("R_Vi").eval(k));
                //                    imviii.push_back(SDP->get_var<double>("Im_Vi").eval(k));
                //                    DebugOn(k<<"\t"<<Wiii.back()-pow(rviii.back(),2)-pow(imviii.back(),2)<<endl);
                //                }
                ////                DebugOn("Wii-rv^2-iv^2"<<endl);
                ////                for(auto it=0;it<Wiii.size();it++)
                ////                {
                ////
                ////                }
                //              //  DebugOn("Value of b\t"<<b_k);
                //                DebugOn("Qf_to"<<endl);
                //                vr=SDP->get_var<double>("Qf_to");
                //                for (auto &k:*vr.get_keys())
                //                {
                //                    DebugOn(k<<"\t"<<vr.eval(k)<<endl);
                //                    qt.push_back(vr.eval(k));
                //                }
                //                for (auto it=0;it<pf.size();it++)
                //                {
                //                         DebugOn("RVI\t"<<rvi[it]<<endl);
                //                }
                //                DebugOn("Real_Loss"<<endl);
                //                for (auto it=0;it<pf.size();it++)
                //                {
                //                    DebugOn(pf[it]+pt[it]<<endl);
                //
                //                }
                //                  DebugOn("Imaginary_Loss"<<endl);
                //                for (auto it=0;it<qf.size();it++)
                //                {
                //                    DebugOn(qf[it]+qt[it]<<endl);
                //
                //                }
                //
                //                DebugOn("Lij"<<endl);
                //                 vr=SDP->get_var<double>("lij");
                //                for (auto &k:*vr.get_keys())
                //                {
                //                    DebugOn(k<<"\t"<<vr.eval(k)<<endl);
                //                }
                //                for (auto it=0;it<qf.size();it++)
                //                {
                //                    DebugOn(qf[it]+qt[it]<<endl);
                //
                //                }
                
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
        }
        
        DebugOn("Time\t"<<solver_time<<endl);
        DebugOn("Iterations\t"<<iter<<endl);
    }
    
    //    if(time_limit){
    //        DebugOn("Reached Time limit!"<<endl);
    //    }
    //    else {
    //        DebugOn("Terminate\t"<<terminate<<endl);
    //    }
    //    DebugOn("Time\t"<<solver_time<<endl);
    //    DebugOn("Iterations\t"<<iter<<endl);
    //    DebugOn("Variable \t Key \t Interval reduction percentage"<<endl);
    //    vector<double> interval_gap;
    //    double sum=0, avg, num_var=0.0;
    //    for(auto &it:SDP->_vars_name)
    //    {
    //        string vname=it.first;
    //        v=SDP->get_var<double>(vname);
    //        auto v_keys=v.get_keys();
    //        for(auto &key: *v_keys)
    //        { num_var++;
    //            p=vname+"|"+ key;
    //            interval_gap.push_back((interval_original[p]-interval_new[p])/(interval_original[p]+zero_tol)*100.0);
    //            sum+=interval_gap.back();
    //            DebugOn(p<<" " << interval_gap.back()<< " flag = " << fixed_point[p] << endl);
    //        }
    //
    //    }
    //    avg=sum/num_var;
    //
    //    DebugOn("Average interval reduction\t"<<avg<<endl);
    //
    //    SDP->reset_constrs();
    ////    SDP->print();
    //    solver<> SDPLB1(SDP,solv_type);
    //
    //    SDPLB1.run(output = 5, tol = 1e-6, "ma57");
    //
    //    SDP->print_solution();
    //
    //    double gap = 100*(upper_bound - lower_bound)/upper_bound;
    //    DebugOn("Initial Gap = " << to_string(gap) << "%."<<endl);
    //    gap = 100*(upper_bound - SDP->get_obj_val())/upper_bound;
    //
    //    DebugOn("Final Gap = " << to_string(gap) << "%."<<endl);
    //    DebugOn("Upper bound = " << to_string(upper_bound) << "."<<endl);
    //    DebugOn("Lower bound = " << to_string(SDP->get_obj_val()) << "."<<endl);
    //    DebugOn("Time\t"<<solver_time<<endl);
    //    DebugOn("\nResults: " << grid._name << " " << to_string(SDP->get_obj_val()) << " " <<endl);
    return 0;
    
}
