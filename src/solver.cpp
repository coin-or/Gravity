//
//  solver.cpp
//  PowerTools++
//
//  Created by Hassan on 30/01/2015.
//  Copyright (c) 2015 NICTA. All rights reserved.
//

#include <gravity/solver.h>

#ifdef USE_BONMIN
#include <coin/BonBonminSetup.hpp>
#include <coin/BonCbc.hpp>
#endif
using namespace gravity;
using namespace std;
namespace {
    void gurobiNotAvailable()
    {
        cerr << "Can't use Gurobi as a solver: this version of Gravity "
        "was compiled without Gurobi support." << endl;
        exit(1);
    }
    void cplexNotAvailable()
    {
        cerr << "Can't use Cplex as a solver: this version of Gravity "
        "was compiled without Cplex support." << endl;
        exit(1);
    }
    
    void bonminNotAvailable()
    {
        cerr << "Can't use Bonmin as a solver: this version of Gravity "
        "was compiled without Bonmin support." << endl;
        exit(1);
    }
    void ipoptNotAvailable()
    {
        cerr << "Can't use Ipopt as a solver: this version of Gravity "
        "was compiled without Ipopt support." << endl;
        exit(1);
    }
    
    void mosekNotAvailable()
    {
        cerr << "Can't use Mosek as a solver: this version of Gravity "
        "was compiled without Mosek support." << endl;
        exit(1);
    }
    
}


solver::solver(Model& model, SolverType stype){
    _stype = stype;
    _model = &model;
    if (_stype==ipopt) {
#ifdef USE_IPOPT
        bool has_int = false;
        //_model->relax();
        //_model->unrelax();
        for (auto &v_p:_model->_vars_name) {
            if (v_p.second->is_integer() || v_p.second->is_binary()) {
                has_int = true;
                auto v = v_p.second;
                auto new_v = new var<double>(v_p.second->_name, 0,1);
                new_v->copy(*v);
                new_v->_is_relaxed = true;
                new_v->_val->resize(new_v->get_dim());
                if (v->get_intype()==integer_) {
                    auto double_var = (var<int>*)v;
                    for (int i = 0; i < double_var->get_dim(); i++) {
                        new_v->_val->at(i) = double_var->_val->at(i);
                    }
                }
                if (v->get_intype()==short_) {
                    auto double_var = (var<short>*)v;
                    for (int i = 0; i < double_var->get_dim(); i++) {
                        new_v->_val->at(i) = double_var->_val->at(i);
                    }
                }
                if (v->get_intype()==binary_) {
                    auto double_var = (var<bool>*)v;
                    _model->_bin_vars[v_p.second->get_vec_id()] = *double_var;
                    for (int i = 0; i < double_var->get_dim(); i++) {
                        new_v->_val->at(i) = double_var->_val->at(i);
                    }
                }
                v_p.second = new_v;
            }
        }
        for (auto &v_p:_model->_vars) {
            if (v_p.second->is_integer() || v_p.second->is_binary()) {
                auto name = v_p.second->_name;
                delete v_p.second;
                v_p.second = _model->get_var(name);
            }
        }
        if(has_int){
            _model->_obj.relax(_model->_vars);
            for (auto &c_p: _model->_cons) {
                c_p.second->relax(_model->_vars);
            }
        }
        

        SmartPtr<IpoptApplication> iapp = IpoptApplicationFactory();
        iapp->RethrowNonIpoptException(true);
        ApplicationReturnStatus status = iapp->Initialize();
        
        if (status != Solve_Succeeded) {
            throw invalid_argument("*** Error during initialization!\n");
        }
        
        if(_model->_objt==maximize){
            _model->_obj *= -1;
        }
        _prog = new IpoptProgram(_model);
#else
        ipoptNotAvailable();
#endif
    }
    else if(_stype==gurobi)
    {
#ifdef USE_GUROBI
    _prog = new GurobiProgram(_model);
#else
        gurobiNotAvailable();
#endif
    }
    else if(_stype==cplex)
    {
#ifdef USE_CPLEX
    _prog = new CplexProgram(_model);
#else
    cplexNotAvailable();
#endif
    }
    else if(_stype == Mosek)
    {
#ifdef USE_MOSEK
    _prog = new MosekProgram(_model);
#else
    mosekNotAvailable();
#endif
    }
    else if(_stype==bonmin) {
#ifdef USE_BONMIN
//        bool has_int = false;
//        //_model->relax();
//        //_model->unrelax();
//        for (auto &v_p:_model->_vars_name) {
//            if (v_p.second->is_integer() || v_p.second->is_binary()) {
//                has_int = true;
//                auto v = v_p.second;
//                auto new_v = new var<double>(v_p.second->_name, 0,1);
//                new_v->copy(*v);
//                new_v->_is_relaxed = true;
//                new_v->_val->resize(new_v->get_dim());
//                if (v->get_intype()==integer_) {
//                    auto double_var = (var<int>*)v;
//                    for (int i = 0; i < double_var->get_dim(); i++) {
//                        new_v->_val->at(i) = double_var->_val->at(i);
//                    }
//                }
//                if (v->get_intype()==short_) {
//                    auto double_var = (var<short>*)v;
//                    for (int i = 0; i < double_var->get_dim(); i++) {
//                        new_v->_val->at(i) = double_var->_val->at(i);
//                    }
//                }
//                if (v->get_intype()==binary_) {
//                    auto double_var = (var<bool>*)v;
//                    _model->_bin_vars[v_p.second->get_vec_id()] = *double_var;
//                    for (int i = 0; i < double_var->get_dim(); i++) {
//                        new_v->_val->at(i) = double_var->_val->at(i);
//                    }
//                }
//                v_p.second = new_v;
//            }
//        }
//        for (auto &v_p:_model->_vars) {
//            if (v_p.second->is_integer() || v_p.second->is_binary()) {
//                auto name = v_p.second->_name;
//                delete v_p.second;
//                v_p.second = _model->get_var(name);
//            }
//        }
//        if(has_int){
//            _model->_obj.relax(_model->_vars);
//            for (auto &c_p: _model->_cons) {
//                c_p.second->relax(_model->_vars);
//            }
//        }
        if(_model->_objt==maximize){
            _model->_obj *= -1;
        }
        _prog = new BonminProgram(_model);
#else
        bonminNotAvailable();
#endif
    }
}

solver::~solver(){
    delete _prog;
}
void run_models(const std::vector<shared_ptr<Model>>& models, int start, int end, SolverType stype, double tol){
    for (auto i = start; i<end; i++) {
        solver(*(models.at(i)),stype).run(5, false, tol, 0.01, "ma27");
        models.at(i)->print_solution(false);
    }
}


void run_parallel(const std::vector<shared_ptr<Model>>& models, SolverType stype, double tol, unsigned nr_threads){
    std::vector<thread> threads;
    std::vector<bool> feasible;
    /* Split models into nr_threads parts */
    std::vector<int> limits = bounds(nr_threads, models.size());
    DebugOff("limits size = " << limits.size() << endl);
    for (int i = 0; i < limits.size(); ++i) {
        DebugOff("limits[" << i << "] = " << limits[i] << endl);
    }
    /* Launch all threads in parallel */
    for (int i = 0; i < nr_threads; ++i) {
        threads.push_back(thread(run_models, ref(models), limits[i], limits[i+1], stype, tol));
    }
    /* Join the threads with the main thread */
    for(auto &t : threads){
        t.join();
    }
}

int solver::run(int print_level, bool relax, double tol, double mipgap, const string& lin_solver, const string& mehrotra, int max_iter){
    int return_status = -1;
    bool violated_constraints = true, optimal = true;
    while(violated_constraints && optimal){
        if (_stype==ipopt) {
    #ifdef USE_IPOPT
//            double mu_init = exp(-1)/exp(10);
            double mu_init = exp(-1)/exp(2);
            SmartPtr<IpoptApplication> iapp = IpoptApplicationFactory();
            iapp->RethrowNonIpoptException(true);
            ApplicationReturnStatus status = iapp->Initialize();
        
            if (status != Solve_Succeeded) {
                std::cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
                return (int) status;
            }

            iapp->Options()->SetStringValue("linear_solver", lin_solver);
            iapp->Options()->SetStringValue("mehrotra_algorithm", mehrotra);
            iapp->Options()->SetNumericValue("tol", tol);
            iapp->Options()->SetIntegerValue("print_level", print_level);
            
            /** Bonmin options */
//            iapp->Options()->SetStringValue("mu_strategy", "adaptive");
//            iapp->Options()->SetStringValue("mu_oracle", "probing");
//            iapp->Options()->SetNumericValue("gamma_phi", 1e-8);
//            iapp->Options()->SetNumericValue("gamma_theta", 1e-4);
//            iapp->Options()->SetNumericValue("bound_push", 1e-12);
//            iapp->Options()->SetNumericValue("bound_frac", 1e-12);
//            iapp->Options()->SetNumericValue("slack_bound_frac", 1e-12);
//            iapp->Options()->SetNumericValue("slack_bound_push", 1e-12);
//            iapp->Options()->SetNumericValue("constr_viol_tol", 1e-8);
//            iapp->Options()->SetNumericValue("dual_inf_tol", 1);
//            iapp->Options()->SetNumericValue("compl_inf_tol", 1e-3);
//            iapp->Options()->SetNumericValue("bound_relax_factor", 1e-9);
//            iapp->Options()->SetNumericValue("bound_relax_factor", 0);
//            iapp->Options()->SetStringValue("derivative_test", "second-order");
            /** Hot start if already solved */
            if (!_model->_first_run) {
//            if (false) {
                DebugOn("Using Hot Start!\n");
                iapp->Options()->SetNumericValue("mu_init", mu_init);
                iapp->Options()->SetStringValue("warm_start_init_point", "yes");
                
//                iapp->Options()->SetNumericValue("warm_start_bound_push", 1e-12);
//                iapp->Options()->SetNumericValue("warm_start_bound_frac", 1e-12);
//                iapp->Options()->SetNumericValue("warm_start_slack_bound_frac", 1e-12);
//                iapp->Options()->SetNumericValue("warm_start_slack_bound_push", 1e-12);
//                iapp->Options()->SetNumericValue("warm_start_mult_bound_push", 1e-12);
            }
            _model->_first_run = false;
            
            
    //mu_init, warm_start_mult_bound_push, warm_start_slack_bound_push, warm_start_bound_push
            
    //                        iapp->Options()->SetStringValue("hessian_approximation", "limited-memory");
    //                        iapp->Options()->SetStringValue("hessian_constant", "yes");
    //                        iapp->Options()->SetStringValue("derivative_test", "only-second-order");
                iapp->Options()->SetNumericValue("ma57_pre_alloc", 10);
//            iapp->Options()->SetIntegerValue("ma57_small_pivot_flag", 1);
            iapp->Options()->SetNumericValue("ma27_liw_init_factor", 20);
            iapp->Options()->SetNumericValue("ma27_la_init_factor", 20);
            iapp->Options()->SetNumericValue("ma27_meminc_factor", 3);
//            iapp->Options()->SetStringValue("ma57_automatic_scaling", "yes");
//                            iapp->Options()->SetStringValue("derivative_test", "second-order");
//                            iapp->Options()->SetNumericValue("derivative_test_perturbation", 1e-6);
    //                        iapp->Options()->SetNumericValue("print_level", 10);

//                            iapp->Options()->SetStringValue("derivative_test", "second-order");
            
                            iapp->Options()->SetIntegerValue("max_iter", max_iter);

    //                        iapp->Options()->SetNumericValue("ma27_liw_init_factor", 50);
    //                        iapp->Options()->SetNumericValue("ma27_pivtol", 1e-6);
    //                        iapp->Options()->SetNumericValue("ma27_la_init_factor", 100);
    //                        iapp->Options()->SetNumericValue("ma27_meminc_factor", 5);
    //                        iapp->Options()->SetStringValue("mu_strategy", "adaptive");
                                iapp->Options()->SetNumericValue("tol", 1e-6);
//                            iapp->Options()->SetNumericValue("dual_inf_tol", 1e-6);
    //            iapp->Options()->SetStringValue("derivative_test", "second-order");
//                            iapp->Options()->SetNumericValue("bound_relax_factor", 1e-14);
    //            iapp.Options()->SetIntegerValue("print_level", 5);

    //                        iapp->Options()->SetStringValue("derivative_test_print_all", "yes");
            _prog->update_model();
            
            SmartPtr<TNLP> tmp = new IpoptProgram(_model);
            status = iapp->OptimizeTNLP(tmp);
                if (status == Solve_Succeeded) {
                    optimal = true;
                    for (auto &v_p:_model->_bin_vars) {
                        auto bin_var = v_p.second;
                        auto double_var = (var<double>*)_model->get_var(v_p.first);
                        for (int i = 0; i < double_var->get_dim(); i++) {
                            if(round(double_var->_val->at(i))==1){
                                bin_var._val->at(i) = true;
                            }
                            else{
                                bin_var._val->at(i) = false;
                            }
                        }
                    }
                    // Retrieve some statistics about the solve
                    
                    //                printf("\n\nSolution of the primal variables:\n");
                    //                _model->print_solution();
    //                return status;
//                    return_status = status;
                }
            else if (status == Solved_To_Acceptable_Level) {
//                return 150;
                optimal = true;
            }
            else {
                optimal = false;
            }
            return_status = optimal ? 100 : -1;
    #else
            ipoptNotAvailable();
    #endif
        }
        else if(_stype==gurobi)
        {
    #ifdef USE_GUROBI
            try{

                auto grb_prog = (GurobiProgram*)(_prog);
                grb_prog->_output = print_level;
    //            prog.grb_prog->reset_model();
                grb_prog->prepare_model();
                optimal = grb_prog->solve(relax,mipgap);
                return_status = optimal ? 100 : -1;
            }catch(GRBException e) {
                cerr << "\nError code = " << e.getErrorCode() << endl;
                cerr << e.getMessage() << endl;
            }
    #else
            gurobiNotAvailable();
    #endif
        }
        else if(_stype==cplex)
        {
    #ifdef USE_CPLEX
            try{
                auto cplex_prog = (CplexProgram*)(_prog);
                cplex_prog->_output = print_level;
                cplex_prog->prepare_model();
                optimal = cplex_prog->solve(relax,mipgap);
                
                return_status = optimal ? 100 : -1;
            }
            catch(IloException e) {
                cerr << e.getMessage() << endl;
            }
    #else
            cplexNotAvailable();
    #endif
        }
        else if(_stype == Mosek)
        {
    #ifdef USE_MOSEK
            try{
                auto mosek_prog = (MosekProgram*)(_prog);
                mosek_prog->_output = print_level;
                mosek_prog->prepare_model();
                bool ok = mosek_prog->solve(relax);
                return_status = ok ? 100 : -1;
            }
            catch(mosek::fusion::FusionException e) {
                cerr << e.toString() << endl;
            }
    #else
            mosekNotAvailable();
    #endif
        }
        else if(_stype==bonmin) {
    #ifdef USE_BONMIN
            
//            using namespace Bonmin;
            BonminSetup bonmin;
            bonmin.initializeOptionsAndJournalist();
            SmartPtr<TMINLP> tmp = new BonminProgram(_model);
            bonmin.initialize(tmp);

            bool ok = false;
    //        bonmin.initialize(prog.bonmin_prog);
            try {
                Bab bb;
                bb(bonmin);
                auto status = bb.mipStatus();
                switch (bb.mipStatus())
                {
                    case Bab::MipStatuses::FeasibleOptimal:
                    case Bab::MipStatuses::Feasible:
                        ok = true;
                        break;
                    default:
                        ok = false;
                }
            }
            catch(TNLPSolver::UnsolvedError *E) {
                //There has been a failure to solve a problem with Bonmin.
                std::cerr << "Bonmin has failed to solve a problem" << std::endl;
                ok = false;
            }
            catch(OsiTMINLPInterface::SimpleError &E) {
                std::cerr << E.className() << "::" << E.methodName() << std::endl << E.message() << std::endl;
                ok = false;
            }
            catch(CoinError &E) {
                std::cerr << E.className() << "::" << E.methodName() << std::endl << E.message() << std::endl;
                ok = false;
            }
            catch(...) {
                std::cerr << "Unknown error: Bonmin failed to solve the problem." << std::endl;
                ok = false;
            }

            if(_model->_objt==maximize){
                _model->_obj_val *= -1;
            }
            
            return_status = ok ? 100 : -1;
    #else
            bonminNotAvailable();
    #endif
        }
        if (optimal) {
            violated_constraints = _model->has_violated_constraints(tol);
            if (violated_constraints) {
                _model->reindex();
            }
        }
        if(_stype!=ipopt) {
            violated_constraints = false;
        }
    }
    _model->_status = return_status;
    return return_status;
}

