//
//  solver.cpp
//  PowerTools++
//
//  Created by Hassan on 30/01/2015.
//  Copyright (c) 2015 NICTA. All rights reserved.
//

#include <Gravity/solver.h>

#ifdef USE_BONMIN
#include <coin/BonBonminSetup.hpp>
#include <coin/BonCbc.hpp>
#endif

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
    
}


solver::solver(Model& model, SolverType stype){
    _stype = stype;
    switch (stype) {
        case cplex:
#ifdef USE_CPLEX
            prog.cplex_prog = new CplexProgram(&model);
#else
            cplexNotAvailable();
#endif
            break;
        case ipopt:
#ifdef USE_IPOPT
            prog.ipopt_prog = new IpoptProgram(&model);
#else
            ipoptNotAvailable();
#endif
            break;
        case gurobi:
#ifdef USE_GUROBI
            try{
                prog.grb_prog = new GurobiProgram(&model);
            }catch(GRBException e) {
                cerr << "\nError code = " << e.getErrorCode() << endl;
                cerr << e.getMessage() << endl;
                exit(1);
            }
#else
            gurobiNotAvailable();
#endif
            break;
        case bonmin:
#ifdef USE_BONMIN
            prog.bonmin_prog = new BonminProgram(model);
#else
            bonminNotAvailable();
#endif
            break;
    }
}

solver::~solver(){
    switch (_stype) {
        case cplex:
#ifdef USE_CPLEX
            delete prog.cplex_prog;
#else
            cplexNotAvailable();
#endif
            break;
        case ipopt:
#ifdef USE_IPOPT
            delete prog.ipopt_prog;
#else
            ipoptNotAvailable();
#endif
            break;
        case gurobi:
#ifdef USE_GUROBI
            delete prog.grb_prog;
#else
            gurobiNotAvailable();
#endif
            break;
        case bonmin:
#ifdef USE_BONMIN
            delete prog.bonmin_prog;
#else
            bonminNotAvailable();
#endif
            break;
    }

}

void solver::set_model(Model& m) {
    
    switch (_stype) {
        case cplex:
#ifdef USE_CPLEX
            prog.cplex_prog->_model = &m;
#else
            cplexNotAvailable();
#endif
        break;
        case ipopt:
#ifdef USE_IPOPT
            prog.ipopt_prog->_model = &m;
#else
            ipoptNotAvailable();
#endif
            break;
        case gurobi:
#ifdef USE_GUROBI
            prog.grb_prog->_model = &m;
#else
            gurobiNotAvailable();
#endif
            break;
        case bonmin:
#ifdef USE_BONMIN
            prog.bonmin_prog->_model = &m;
#else
            bonminNotAvailable();
#endif
            break;
    }
}


int solver::run(int output, bool relax){
    //GurobiProgram* grbprog;
    // Initialize the IpoptApplication and process the options

    if (_stype==ipopt) {
#ifdef USE_IPOPT
        SmartPtr<IpoptApplication> iapp = new IpoptApplication();
            iapp->RethrowNonIpoptException(true);
            ApplicationReturnStatus status = iapp->Initialize();
        
            if (status != Solve_Succeeded) {
                std::cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
                return (int) status;
            }

        if(prog.ipopt_prog->_model->_objt==maximize){
            prog.ipopt_prog->_model->_obj *= -1;
        }
        SmartPtr<TNLP> tmp = new IpoptProgram(prog.ipopt_prog->_model);
//        prog.ipopt_prog;
            //            iapp.Options()->SetStringValue("hessian_constant", "yes");
                        iapp->Options()->SetStringValue("derivative_test", "first-order");
//                        iapp->Options()->SetStringValue("derivative_test", "second-order");
            //            iapp->Options()->SetNumericValue("tol", 1e-6);
//                        iapp.Options()->SetNumericValue("tol", 1e-6);
            //            iapp->Options()->SetStringValue("derivative_test", "second-order");
            //            iapp.Options()->SetNumericValue("bound_relax_factor", 0);
            //            iapp.Options()->SetIntegerValue("print_level", 5);
            
            //            iapp.Options()->SetStringValue("derivative_test_print_all", "yes");
        status = iapp->OptimizeTNLP(tmp);
        if(prog.ipopt_prog->_model->_objt==maximize){
            prog.ipopt_prog->_model->_obj *= -1;
        }
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
    else if(_stype==gurobi)
    {
#ifdef USE_GUROBI
        try{
            //                prog.grbprog = new GurobiProgram();
            prog.grb_prog->_output = output;
//            prog.grb_prog->reset_model();
            prog.grb_prog->prepare_model();
            bool ok = prog.grb_prog->solve(relax);
            return ok ? 100 : -1;
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
            //                prog.grbprog = new GurobiProgram();
            prog.cplex_prog->_output = output;
            //            prog.grb_prog->reset_model();
            prog.cplex_prog->prepare_model();            
            bool ok = prog.cplex_prog->solve(relax);
            return ok ? 100 : -1;
        }
        catch(IloException e) {
            cerr << e.getMessage() << endl;
        }
#else
        cplexNotAvailable();
#endif
    }
    else if(_stype==bonmin) {
#ifdef USE_BONMIN
        if(prog.bonmin_prog->model->_objt==maximize){
            *prog.bonmin_prog->model->_obj *= -1;
        }

        bool ok = true;
        using namespace Bonmin;
        BonminSetup bonmin;
        bonmin.initializeOptionsAndJournalist();
        bonmin.initialize(prog.bonmin_prog);
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

        if(prog.bonmin_prog->model->_objt==maximize){
            *prog.bonmin_prog->model->_obj *= -1;
        }

        return ok ? 100 : -1;
#else
        bonminNotAvailable();
#endif
    }
    return -1;
}
