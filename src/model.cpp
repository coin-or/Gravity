//
//  model.cpp
//  Gravity
//
//  Created by Hijazi, Hassan.
//
//
//
#include <gravity/model.h>
#include <gravity/solver.h>
#include <math.h> //for setting the rounding direction

using namespace std;
namespace gravity {
    
    template <typename type>
    template<typename T,
    typename std::enable_if<is_same<type,double>::value>::type*>
    void Model<type>::run_obbt(double max_time, unsigned max_iter, const pair<bool,double>& upper_bound, unsigned precision) {
                
#ifdef USE_MPI
        auto err_init = MPI_Init(nullptr,nullptr);
        int worker_id, nb_workers;
        auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
        auto err_size = MPI_Comm_size(MPI_COMM_WORLD, &nb_workers);
#endif
        int nb_threads = thread::hardware_concurrency();
        int nb_total_threads = nb_threads; /** Used when MPI is ON to multipply with the number of workers */
#ifdef USE_MPI
        nb_total_threads *= nb_workers;
#endif
        vector<shared_ptr<Model<>>> batch_models;
        int output = 0;
        SolverType solv_type = ipopt;
        double tol = 1e-6;

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
        const int gap_count_int=6;


        double solver_time_end, solver_time =0, solver_time_start = get_wall_time(), gap;
        shared_ptr<map<string,size_t>> p_map;
        //Check if gap is already not zero at root node
        terminate=false;
        for(auto &it:this->_vars_name)
        {
            string vname=it.first;
            v=this->get_var<T>(vname);
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


                interval_original[p]=v.get_ub(key)-v.get_lb(key);
                ub_original[p]=v.get_ub(key);
                lb_original[p]=v.get_lb(key);
                interval_new[p]=v.get_ub(key)-v.get_lb(key);

            }

        }

        solver_time= get_wall_time()-solver_time_start;
        auto v_in_cons=this->_v_in_cons;
        while(solver_time<=max_time && !terminate && iter<=max_iter)
        {
            iter++;
//            terminate=true;
            for (auto it=this->_vars_name.begin(); it!=this->_vars_name.end(); it++)
            {
                vname=it->first;
                v = this->get_var<T>(vname);
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
                    if(fixed_point[p]==false || (next(it)==this->_vars_name.end() && next(it_key)==v.get_keys()->end()))
                    {
                        //Loop on directions, upper bound and lower bound
                        for(auto &dir: dir_array)
                        {
                            auto modelk = this->copy();
                            if(upper_bound.first){/* Add upperbound constraint */
                                Constraint<> obj_UB("obj_UB");
                                obj_UB=*(this->_obj) - upper_bound.second;
                                modelk->add(obj_UB<=0);
                            }
                            mname=vname+"|"+key+"|"+dir;
                            modelk->set_name(mname);
                            vark=modelk->template get_var<T>(vname);
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
                            if (batch_models.size()==nb_total_threads || (next(it)==this->_vars_name.end() && next(it_key)==v.get_keys()->end() && dir=="UB"))
                            {
                                double batch_time_start = get_wall_time();
#ifdef USE_MPI
                                run_MPI(batch_models,ipopt,1e-6,nb_threads, "ma57",true);
#else
                                run_parallel(batch_models,ipopt,1e-6,nb_threads, "ma57");
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
                                    vk=this->get_var<T>(vkname);
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
                                                auto objk_rounded = floor(objk*std::pow(10,precision))/std::pow(10,precision);
//                                                vk.set_lb(keyk, objk);/* IN MPI this needs to be broadcasted back to the other workers */
                                                  vk.set_lb(keyk, objk_rounded);/* IN MPI this needs to be broadcasted back to the other workers */
                                            }
                                            else{
                                                auto objk_rounded = ceil(objk*std::pow(10,precision))/std::pow(10,precision);
                                                vk.set_ub(keyk, objk_rounded);
//                                                vk.set_ub(keyk, objk);
                                            }
                                            //If crossover in bounds,just exchange them
                                            if(vk.get_ub(keyk)<vk.get_lb(keyk))
                                            {
                                                fixed_point[pk]=true;
                                                temp=vk.get_ub(keyk);
                                                tempa=vk.get_lb(keyk);
                                                auto tempa_rounded = ceil(tempa*std::pow(10,precision))/std::pow(10,precision);
                                                auto temp_rounded = floor(temp*std::pow(10,precision))/std::pow(10,precision);
                                                vk.set_ub(keyk, tempa_rounded);
                                                vk.set_lb(keyk, temp_rounded);
//                                                vk.set_ub(keyk, tempa);
//                                                vk.set_lb(keyk, temp);

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
                                                    auto right_rounded = ceil(right*std::pow(10,precision))/std::pow(10,precision);
                                                    auto left_rounded = floor(left*std::pow(10,precision))/std::pow(10,precision);
                                                    vk.set_ub(keyk, right_rounded);
                                                    vk.set_lb(keyk, left_rounded);
//                                                    vk.set_ub(keyk, right);
//                                                    vk.set_lb(keyk, left);
                                                }
                                                //If resized interval crosses original upperbound, set the new bound to upperbound, and lower bound is expanded to upperbound-range_tolerance
                                                else if(right>ub_original[pk])
                                                {
                                                    auto ub_rounded = ceil(ub_original[pk]*std::pow(10,precision))/std::pow(10,precision);
                                                    auto lb_rounded = floor((ub_original[pk]-range_tol)*std::pow(10,precision))/std::pow(10,precision);
                                                    vk.set_ub(keyk, ub_rounded);
                                                    vk.set_lb(keyk, lb_rounded);
//                                                    vk.set_ub(keyk, ub_original[pk]);
//                                                    vk.set_lb(keyk, ub_original[pk]-range_tol);
                                                }
                                                //If resized interval crosses original lowerbound, set the new bound to lowerbound, and upper bound is expanded to lowerbound+range_tolerance
                                                else if(left<lb_original[pk])
                                                {
                                                    auto lb_rounded = floor(lb_original[pk]*std::pow(10,precision))/std::pow(10,precision);
                                                    auto ub_rounded = ceil((lb_original[pk]+range_tol)*std::pow(10,precision))/std::pow(10,precision);
                                                    vk.set_lb(keyk, lb_rounded);
                                                    vk.set_ub(keyk, ub_rounded);
//                                                    vk.set_lb(keyk, lb_original[pk]);
//                                                    vk.set_ub(keyk, lb_original[pk]+range_tol);

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
            //                    if(iter%gap_count_int==0)
            //                    {
            //                        this->reset_constrs();
            //                        solver<> SDPLB1(*this,solv_type);
            //                        SDPLB1.run(output = 5, tol=1e-8);
            //                        if(this->_status==0)
            //                        {
            //                            gap = 100*(upper_bound - (this->get_obj_val()))/upper_bound;
            //                            DebugOn("Gap "<<gap<<endl);
            //                        }
            //                        if (upper_bound-this->get_obj_val()<=upp_low_tol || (upper_bound-this->get_obj_val())/(upper_bound+zero_tol)<=upp_low_tol)
            //                        {
            //                            DebugOn("Gap closed at iter "<< iter<<endl);
            //                            close=true;
            //                            terminate=true;
            //
            //                        }
            //                    }

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
        for(auto &it:this->_vars_name)
        {
            string vname=it.first;
            v=this->get_var<T>(vname);
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

//        if(!close)
//        {
//
//            this->reset_constrs();
//            solver<T> SDPLB1(*this,solv_type);
//
//            SDPLB1.run(output = 5, tol=1e-8);
//        }
//
//        avg=sum/num_var;
//
//        DebugOn("Average interval reduction\t"<<avg<<endl);

        if(!close)
        {
#ifdef USE_MPI
            if(worker_id==0){
#endif
                this->reset_constrs();
                solver<T> SDPLB1(*this,solv_type);

                SDPLB1.run(output = 5, tol=1e-6);
                this->print_constraints_stats(tol);
                bool print_only_relaxed;
                this->print_nonzero_constraints(tol,print_only_relaxed=true);

                //        this->print_solution();

                //        this->print();

                //                if(this->_status==0)
                //                {
                //
                //                    DebugOn("\nResults: " << grid._name << " " << to_string(this->get_obj_val()) << " " <<endl);
                //                    DebugOn("Solution Print"<<endl);
                //                    //                this->print_solution();
                //                    this->print_constraints_stats(tol);
                //                    gap = 100*(upper_bound - lower_bound)/upper_bound;
                //                    DebugOn("Initial Gap = " << to_string(gap) << "%."<<endl);
                //                    gap = 100*(upper_bound - (this->get_obj_val()))/upper_bound;
                //                    DebugOn("Final Gap = " << to_string(gap) << "%."<<endl);
                //                    DebugOn("Upper bound = " << to_string(upper_bound) << "."<<endl);
                //                    DebugOn("Lower bound = " << to_string((this->get_obj_val())) << "."<<endl);
                //                    DebugOn("Time\t"<<solver_time<<endl);
                //
                //                }
                //                else
                //                {
                //                    double gap = 100*(upper_bound - lower_bound)/upper_bound;
                //                    DebugOn("Initial Gap = " << to_string(gap) << "%."<<endl);
                //                    DebugOn("Lower bounding problem status = " << this->_status <<endl);
                //                    DebugOn("Lower bounding problem not solved to optimality, cannot compute final gap"<<endl);
                //                }
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

    }

    template void gravity::Model<double>::run_obbt<double, (void*)0>(double, unsigned int, const pair<bool,double>&, unsigned int);
//    template void Model<double>::run_obbt(double max_time, unsigned max_iter);
//    template func<double> constant<double>::get_real() const;
//    template class Model<double>;
//    template class Model<Cpx>;

}


