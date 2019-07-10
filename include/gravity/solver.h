//
//  solver.h
//  Gravity
//
//  Created by Hassan

#ifndef __Gravity____Solver__
#define __Gravity____Solver__

#include <stdio.h>

#include <gravity/GravityConfig.h>

#include <gravity/model.h>
#ifdef USE_MPI
#include <mpi.h>
#endif
#ifdef USE_IPOPT
#include <gravity/IpoptProgram.h>
//#include "IpoptInterfaceCommon.h"
#include <coin/IpRegOptions.hpp>
#include <coin/IpJournalist.hpp>
#include <coin/IpIpoptApplication.hpp>
#include <coin/IpSolveStatistics.hpp>


using Ipopt::IsValid;
using Ipopt::RegisteredOption;
using Ipopt::EJournalLevel;
using Ipopt::Journal;
using Ipopt::IpoptApplication;
using Ipopt::SmartPtr;
using Ipopt::TNLP;
using Ipopt::ApplicationReturnStatus;
using Ipopt::SolveStatistics;
#endif
#ifdef USE_GUROBI
#include <gravity/GurobiProgram.h>
#endif
#ifdef USE_BONMIN
#include <gravity/BonminProgram.h>
#endif
#ifdef USE_CPLEX
#include <gravity/CplexProgram.h>
#endif
#ifdef USE_CLP
#include <gravity/ClpProgram.h>
#endif
//#ifdef USE_SDPA
//#include "SdpaProgram.h"
//#endif
#ifdef USE_MOSEK
#include "MosekProgram.h"
#endif


void gurobiNotAvailable();

void cplexNotAvailable();

void bonminNotAvailable();

void ipoptNotAvailable();

void mosekNotAvailable();

void ClpNotAvailable();

namespace gravity {

    
    template<typename type = double>
    class solver {
    public:
        shared_ptr<gravity::Model<type>>          _model = nullptr;
        shared_ptr<Program<type>>                 _prog = nullptr;
        SolverType                                _stype;
        type                                      _tol = 1e-6; /*<< Solver tolerance. */
        unsigned                                  _nb_iterations = 0;
        
        /** Constructor */
        //@{
        solver();

        solver(shared_ptr<gravity::Model<type>> model, SolverType stype){
            _stype = stype;
            _model = model;
            _model->_built = true;
            init();
        }
            
        solver(gravity::Model<type>& model, SolverType stype){
            _stype = stype;
            _model = make_shared<gravity::Model<type>>(model);
            _model->_built = true;
            init();
        }
        
        unsigned get_nb_iterations(){
            return _nb_iterations;
        };
        
        void init(){
            if (_stype==ipopt) {
#ifdef USE_IPOPT
                _model->replace_integers();
                SmartPtr<IpoptApplication> iapp = IpoptApplicationFactory();
                iapp->RethrowNonIpoptException(true);
                ApplicationReturnStatus status = iapp->Initialize();
                
                if (status != Solve_Succeeded) {
                    throw invalid_argument("*** Error during initialization!\n");
                }
                
                if(_model->_objt==maximize){
                    _model->_obj->reverse_sign();
                }
                _prog = make_shared<IpoptProgram<type>>(_model);
#else
                ipoptNotAvailable();
#endif
            }
            else if(_stype==gurobi)
            {
#ifdef USE_GUROBI
                _prog = make_shared<GurobiProgram>(_model);
#else
                gurobiNotAvailable();
#endif
            }
            else if(_stype==cplex)
            {
#ifdef USE_CPLEX
                _prog = make_shared<CplexProgram>(_model);
#else
                cplexNotAvailable();
#endif
            }
            else if(_stype == _mosek)
            {
#ifdef USE_MOSEK
                _prog = make_shared<MosekProgram>(_model);
#else
                mosekNotAvailable();
#endif
            }
            else if(_stype==bonmin) {
#ifdef USE_BONMIN
                _model->replace_integers();
                if(_model->_objt==maximize){
                    _model->_obj->reverse_sign();
                }
                _prog = make_shared<BonminProgram>(_model);
#else
                bonminNotAvailable();
#endif
            }
            else if (_stype == clp){
#ifdef USE_CLP
                _model->replace_integers();
                _prog = make_shared<ClpProgram>(_model);
#else
                ClpNotAvailable();
#endif
            }
        }        
        //@}
        void set_model(gravity::Model<type>& m);
        int run(bool relax){
            return run(5, 1e-6, 10000, 1e-6, relax, {false,""}, 1e+6);
        }
        int run(int output, type tol , const string& lin_solver){
            if(lin_solver!="")
                return run(output, tol, 10000, 1e-6, true, {true,lin_solver}, 1e+6);
            else
                return run(output, tol, 10000, 1e-6, true, {false,lin_solver}, 1e+6);
        }
        
        int run(type tol , double time_limit){
            return run(5, tol, 10000, 1e-6, true, {false,""}, time_limit);
        }
        
        //OPF.run(tol,time_limit,"ma97");
        int run(type tol , double time_limit, const string& lin_solver){
            return run(5, tol, 10000, 1e-6, true, {true,lin_solver}, time_limit);
        }
        
        int run(int output=5, type tol=1e-6 , int max_iter=10000){
            return run(output, tol, max_iter, 1e-6, false, {false,""}, 1e+6);
        }
        /* run model */
        int run(int output, type tol , int max_iter, double mipgap, bool relax, pair<bool,string> lin_solver, double time_limit){
            int return_status = -1, dyn_max_iter = 20;
            bool violated_constraints = true, optimal = true;
            unsigned nb_it = 0;
            while(violated_constraints && optimal){
                if (_stype==ipopt) {
#ifdef USE_IPOPT
                    double mu_init = std::exp(5)/std::exp(1);
                    SmartPtr<IpoptApplication> iapp = IpoptApplicationFactory();
                    iapp->RethrowNonIpoptException(true);
                    ApplicationReturnStatus status = iapp->Initialize();
                    
                    if (status != Solve_Succeeded) {
                        std::cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
                        return (int) status;
                    }
                    if(lin_solver.first){
                        iapp->Options()->SetStringValue("linear_solver", lin_solver.second);
                    }
//                    iapp->Options()->SetStringValue("mehrotra_algorithm", mehrotra);
                    iapp->Options()->SetNumericValue("tol", tol);
                    iapp->Options()->SetIntegerValue("print_level", output);
                    
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
//                                iapp->Options()->SetStringValue("derivative_test", "second-order");
                    //            iapp->Options()->SetNumericValue("mu_init", mu_init);
                    //            iapp->Options()->SetNumericValue("obj_scaling_factor", 1e-2);
                    /** Hot start if already solved */
                    if (!_model->_first_run) {
//                        dyn_max_iter *= 2;
//                    if (true) {
//                                    if (false) {
                        mu_init = std::exp(-1)/std::exp(2);
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
                    iapp->Options()->SetIntegerValue("max_iter", max_iter);
                    
                    //mu_init, warm_start_mult_bound_push, warm_start_slack_bound_push, warm_start_bound_push
                    
                    //                        iapp->Options()->SetStringValue("hessian_approximation", "limited-memory");
                    //                        iapp->Options()->SetStringValue("hessian_constant", "yes");
                    //                        iapp->Options()->SetStringValue("derivative_test", "only-second-order");
                    //                iapp->Options()->SetNumericValue("ma57_pre_alloc", 10);
                    //            iapp->Options()->SetIntegerValue("ma57_small_pivot_flag", 1);
                    //            iapp->Options()->SetNumericValue("ma27_liw_init_factor", 20);
                    //            iapp->Options()->SetNumericValue("ma27_la_init_factor", 20);
                    //            iapp->Options()->SetNumericValue("ma27_meminc_factor", 3);
                    //            iapp->Options()->SetStringValue("ma57_automatic_scaling", "yes");
                    //                            iapp->Options()->SetStringValue("derivative_test", "second-order");
                    //                            iapp->Options()->SetNumericValue("derivative_test_perturbation", 1e-6);
                    //                        iapp->Options()->SetNumericValue("print_level", 10);
                    
//                                                iapp->Options()->SetStringValue("derivative_test", "second-order");
                    
                    
                    
                    //                        iapp->Options()->SetNumericValue("ma27_liw_init_factor", 50);
                    //                        iapp->Options()->SetNumericValue("ma27_pivtol", 1e-6);
                    //                        iapp->Options()->SetNumericValue("ma27_la_init_factor", 100);
                    //                        iapp->Options()->SetNumericValue("ma27_meminc_factor", 5);
                    //                        iapp->Options()->SetStringValue("mu_strategy", "adaptive");
                    iapp->Options()->SetNumericValue("tol", tol);
                    //                            iapp->Options()->SetNumericValue("dual_inf_tol", 1e-6);
//                                    iapp->Options()->SetStringValue("derivative_test", "second-order");
                    //                            iapp->Options()->SetNumericValue("bound_relax_factor", 0);
                    //            iapp.Options()->SetIntegerValue("print_level", 5);
                    
                    //                        iapp->Options()->SetStringValue("derivative_test_print_all", "yes");
                    if(!_model->_built){ /* Constraints have been added */
                        _model->reindex();
                    }
                    _model->fill_in_maps();
                    
                    SmartPtr<TNLP> tmp = new IpoptProgram<type>(_model);
                    status = iapp->OptimizeTNLP(tmp);
                    if (IsValid(iapp->Statistics())) {
                        SmartPtr<SolveStatistics> stats = iapp->Statistics();
                        _nb_iterations = stats->IterationCount();                        
                    }
                    DebugOff("Return status = " << status << endl);
                    if (status == Solve_Succeeded) {
                        optimal = true;
                        _model->round_solution();
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
                    return_status = optimal ? 0 : -1;
#else
                    ipoptNotAvailable();
#endif
                }
                else if (_stype == clp){
#ifdef USE_CLP
                    auto clp_prog = (ClpProgram*)(_prog.get());
                    clp_prog->prepare_model();
                    optimal = clp_prog->solve();
                    
#else
                    ClpNotAvailable();
#endif
                }
                else if(_stype==gurobi)
                {
#ifdef USE_GUROBI
                    try{
                        
                        auto grb_prog = (GurobiProgram*)(_prog.get());
                        grb_prog->_output = output;
                        //            prog.grb_prog->reset_model();
                        grb_prog->prepare_model();
                        optimal = grb_prog->solve(relax,mipgap);
                        return_status = optimal ? 0 : -1;
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
                        auto cplex_prog = (CplexProgram*)(_prog.get());
                        cplex_prog->_output = output;
                        cplex_prog->prepare_model();
                        optimal = cplex_prog->solve(relax,mipgap);
                        
                        return_status = optimal ? 0 : -1;
                    }
                    catch(IloException e) {
                        cerr << e.getMessage() << endl;
                    }
#else
                    cplexNotAvailable();
#endif
                }
                else if(_stype == _mosek)
                {
#ifdef USE_MOSEK
                    try{
                        auto mosek_prog = (MosekProgram*)(_prog.get());
                        mosek_prog->_output = output;
                        mosek_prog->prepare_model();
                        bool ok = mosek_prog->solve(relax);
                        return_status = ok ? 0 : -1;
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
                    //            bonmin.options()->SetIntegerValue("max_consecutive_infeasible", 100);
                    //            bonmin.options()->SetIntegerValue("solution_limit", 1);
                    //            bonmin.options()->SetIntegerValue("bb_log_level", 3);
                    //            bonmin.options()->SetNumericValue("allowable_gap", -1e5);
                    //            bonmin.options()->SetNumericValue("allowable_fraction_gap", -1e5);
                    //            bonmin.options()->SetNumericValue("tol", tol);
                    //            bonmin.options()->SetIntegerValue("print_level", 5);
                    bonmin.options()->SetStringValue("tree_search_strategy", "top-node");
                    bonmin.options()->SetStringValue("warm_start", "optimum");
                    
                    //            bonmin.options()->SetStringValue("variable_selection", "random");
                    //            bonmin.options()->SetStringValue("linear_solver", "ma57");
                    //            bonmin.options()->SetIntegerValue("nlp_log_at_root", 5);
                    bonmin.options()->SetNumericValue("cutoff", 0.9);
                    //            bonmin.options()->SetNumericValue("resolve_on_small_infeasibility", INFINITY);
                    bonmin.options()->SetStringValue("heuristic_dive_MIP_fractional", "no");
                    bonmin.options()->SetStringValue("heuristic_dive_MIP_vectorLength", "no");
                    bonmin.options()->SetStringValue("heuristic_dive_fractional", "no");
                    bonmin.options()->SetStringValue("heuristic_dive_vectorLength", "no");
                    bonmin.options()->SetStringValue("heuristic_feasibility_pump", "no");
                    bonmin.options()->SetStringValue("pump_for_minlp", "no");
                    //            bonmin.options()->SetStringValue("algorithm", "B-iFP");
                    bonmin.options()->SetStringValue("node_comparison", "best-guess");
                    //            bonmin.options()->SetStringValue("dynamic_def_cutoff_decr", "yes");
                    bonmin.options()->SetIntegerValue("num_resolve_at_root", 5);
                    
                    _model->fill_in_maps();
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
                    
                    return_status = ok ? 0 : -1;
#else
                    bonminNotAvailable();
#endif
                }
                if (optimal) {
                    violated_constraints = _model->has_violated_constraints(tol);
                    if (violated_constraints) {
                        _model->reindex();
//                        _model->print();
//                        _model->print_symbolic();
                        //                Constraint obj_cstr("obj_cstr_"+to_string(nb_it));
                        //                obj_cstr += _model->_obj - _model->_obj_ub;
                        //                if (_model->_objt==minimize) {
                        //                    _model->add(obj_cstr <= 0);
                        //                }
                        //                else {
                        //                    _model->add(obj_cstr >= 0);
                        //                }
                        //                _model->print();
                    }
                }
                //        if(_stype!=ipopt) {
                //            violated_constraints = false;
                //        }
                nb_it++;
            }
            if (nb_it>1) {
                DebugOn(endl << "####################" << endl);
                DebugOn("Solved in " << nb_it << " constraint generation iterations" << endl);
                DebugOn("####################" << endl);
            }
            _model->_status = return_status;
            if(_model->_status == 0){
                DebugOff("Solved to Optimality/Acceptable "<< endl);
            }
            return return_status;
        }
        
    };
    
    template<typename type>
    int run_models(const std::vector<shared_ptr<Model<type>>>& models, size_t start, size_t end, SolverType stype, type tol, const string& lin_solver=""){
        int return_status = -1;
        for (auto i = start; i<end; i++) {
            return_status = solver<type>((models.at(i)),stype).run(0, tol, lin_solver);
            DebugOn("Return status\t"<<return_status);
            //            models.at(i)->print_solution(24);
        }
        return return_status;
    }
    
    template<typename type>
    void run_parallel(const initializer_list<shared_ptr<gravity::Model<type>>>& models, gravity::SolverType stype = ipopt, type tol = 1e-6, unsigned nr_threads=std::thread::hardware_concurrency(), const string& lin_solver=""){
        run_parallel(vector<shared_ptr<gravity::Model<type>>>(models), stype, tol, nr_threads, lin_solver);
    }
    
    /** Runds models stored in the vector in parallel, using solver of stype and tolerance tol */
    template<typename type>
    
    void run_parallel(const vector<shared_ptr<gravity::Model<type>>>& models, gravity::SolverType stype = ipopt, type tol = 1e-6, unsigned nr_threads=std::thread::hardware_concurrency(), const string& lin_solver=""){
        std::vector<thread> threads;
        std::vector<bool> feasible;
        /* Split models into nr_threads parts */
        std::vector<size_t> limits = bounds(nr_threads, models.size());
        DebugOn("Running on " << nr_threads << " threads." << endl);
        DebugOff("limits size = " << limits.size() << endl);
        for (size_t i = 0; i < limits.size(); ++i) {
            DebugOff("limits[" << i << "] = " << limits[i] << endl);
        }
        /* Launch all threads in parallel */
        auto vec = vector<shared_ptr<gravity::Model<type>>>(models);
        for (size_t i = 0; i < nr_threads; ++i) {
            threads.push_back(thread(run_models<type>, ref(vec), limits[i], limits[i+1], stype, tol, lin_solver));
        }
        /* Join the threads with the main thread */
        for(auto &t : threads){
            t.join();
        }
    }
    
    
#ifdef USE_MPI
    /** Runds models stored in the vector in parallel using MPI
     @models vector of models to run in parallel
     @stype Solver type
     @tol numerical tolerance
     @nb_threads Number of parallel threads
     @lin_solver linear system solver
     */
    template<typename type>
    void run_MPI(const vector<shared_ptr<gravity::Model<type>>>& models, gravity::SolverType stype = ipopt, type tol = 1e-6, unsigned nr_threads=std::thread::hardware_concurrency(), const string& lin_solver=""){
        std::vector<thread> threads;
        std::vector<bool> feasible;
        /* Split models into nr_threads parts */
        std::vector<size_t> limits = bounds(nr_threads, models.size());
        DebugOn("Running on " << nr_threads << " threads." << endl);
        DebugOff("limits size = " << limits.size() << endl);
        for (size_t i = 0; i < limits.size(); ++i) {
            DebugOff("limits[" << i << "] = " << limits[i] << endl);
        }
        /* Launch all threads in parallel */
        auto vec = vector<shared_ptr<gravity::Model<type>>>(models);
        for (size_t i = 0; i < nr_threads; ++i) {
            threads.push_back(thread(run_models<type>, ref(vec), limits[i], limits[i+1], stype, tol, lin_solver));
        }
        /* Join the threads with the main thread */
        for(auto &t : threads){
            t.join();
        }
    }
#endif

}

#endif /* defined(__Gravity____Solver__) */
