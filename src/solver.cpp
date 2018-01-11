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
    if (_stype==ipopt_) {
#ifdef USE_IPOPT
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
    else if(_stype==gurobi_)
    {
#ifdef USE_GUROBI
    _prog = new GurobiProgram(_model);
#else
        gurobiNotAvailable();
#endif
    }
    else if(_stype==cplex_)
    {
#ifdef USE_CPLEX
    _prog = new CplexProgram(_model);
#else
    cplexNotAvailable();
#endif
    }
    else if(_stype == mosek_)
    {
#ifdef USE_MOSEK
    _prog = new MosekProgram(_model);
#else
    mosekNotAvailable();
#endif
    }
    else if(_stype==bonmin_) {
#ifdef USE_BONMIN
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

int solver::run(int output, bool relax, double tol, const string& lin_solver, const string& mehrotra){

    if (_stype==ipopt_) {
#ifdef USE_IPOPT
            SmartPtr<IpoptApplication> iapp = IpoptApplicationFactory();
            iapp->RethrowNonIpoptException(true);
            ApplicationReturnStatus status = iapp->Initialize();
        
            if (status != Solve_Succeeded) {
                std::cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
                return (int) status;
            }

            iapp->Options()->SetStringValue("linear_solver", lin_solver);
            iapp->Options()->SetStringValue("mehrotra_algorithm", mehrotra);

//                        iapp->Options()->SetStringValue("hessian_approximation", "limited-memory");
//                        iapp->Options()->SetStringValue("hessian_constant", "yes");
//                        iapp->Options()->SetStringValue("derivative_test", "only-second-order");
//            iapp->Options()->SetStringValue("derivative_test", "second-order");
//                        iapp->Options()->SetStringValue("derivative_test", "first-order");
//                        iapp->Options()->SetNumericValue("derivative_test_perturbation", 1e-6);
//                        iapp->Options()->SetNumericValue("print_level", 10);

//                        iapp->Options()->SetStringValue("derivative_test", "second-order");
            iapp->Options()->SetNumericValue("tol", tol);
//                        iapp->Options()->SetIntegerValue("max_iter", 50);

//                        iapp->Options()->SetNumericValue("ma27_liw_init_factor", 50);
//                        iapp->Options()->SetNumericValue("ma27_pivtol", 1e-6);
//                        iapp->Options()->SetNumericValue("ma27_la_init_factor", 100);
//                        iapp->Options()->SetNumericValue("ma27_meminc_factor", 5);




//                        iapp->Options()->SetStringValue("mu_strategy", "adaptive");
//                        iapp.Options()->SetNumericValue("tol", 1e-6);
//            iapp->Options()->SetStringValue("derivative_test", "second-order");
//                        iapp->Options()->SetNumericValue("bound_relax_factor", 0);
//            iapp.Options()->SetIntegerValue("print_level", 5);

//                        iapp->Options()->SetStringValue("derivative_test_print_all", "yes");
        _prog->update_model();
        
        SmartPtr<TNLP> tmp = new IpoptProgram(_model);
        status = iapp->OptimizeTNLP(tmp);
//        if(prog.ipopt_prog->_model->_objt==maximize){
//            prog.ipopt_prog->_model->_obj *= -1;
//        }
            if (status == Solve_Succeeded) {
                // Retrieve some statistics about the solve
                
                //                printf("\n\nSolution of the primal variables:\n");
                //                _model->print_solution();
//                return status;
                return 100;
            }
        if (status == Solved_To_Acceptable_Level) {
            return 150;
        }
#else
        ipoptNotAvailable();
#endif
    }
    else if(_stype==gurobi_)
    {
#ifdef USE_GUROBI
        try{

            auto grb_prog = (GurobiProgram*)(_prog);
            grb_prog->_output = output;
//            prog.grb_prog->reset_model();
            grb_prog->prepare_model();
            bool ok = grb_prog->solve(relax);
            return ok ? 100 : -1;
        }catch(GRBException e) {
            cerr << "\nError code = " << e.getErrorCode() << endl;
            cerr << e.getMessage() << endl;
        }
#else
        gurobiNotAvailable();
#endif
    }
    else if(_stype==cplex_)
    {
#ifdef USE_CPLEX
        try{
            auto cplex_prog = (CplexProgram*)(_prog);
            cplex_prog->_output = output;
            cplex_prog->prepare_model();
            bool ok = cplex_prog->solve(relax);
            return ok ? 100 : -1;
        }
        catch(IloException e) {
            cerr << e.getMessage() << endl;
        }
#else
        cplexNotAvailable();
#endif
    }
    else if(_stype == mosek_)
    {
#ifdef USE_MOSEK
        try{
            auto mosek_prog = (MosekProgram*)(_prog);
            mosek_prog->_output = output;
            mosek_prog->prepare_model();
            bool ok = mosek_prog->solve(relax);
            return ok ? 100 : -1;
        }
        catch(mosek::fusion::FusionException e) {
            cerr << e.toString() << endl;
        }
#else
        mosekNotAvailable();
#endif
    }
    else if(_stype==bonmin_) {
#ifdef USE_BONMIN
        
        using namespace Bonmin;
        BonminSetup bonmin;
        bonmin.initializeOptionsAndJournalist();
        SmartPtr<TMINLP> tmp = new BonminProgram(_model);
        bonmin.initialize(tmp);

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
        
        return ok ? 100 : -1;
#else
        bonminNotAvailable();
#endif
    }
    return -1;
}
