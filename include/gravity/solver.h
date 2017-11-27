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
#ifdef USE_IPOPT
#include <gravity/IpoptProgram.h>
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
//#ifdef USE_SDPA
//#include "SdpaProgram.h"
//#endif
#ifdef USE_MOSEK
#include "MosekProgram.h"
#endif

namespace gravity {
    class solver {
    public:
        Model*                         _model;
        SolverType                     _stype;
        /** Constructor */
        //@{
        solver();

        solver(Model& model, SolverType stype);
        //@}
        void set_model(Model& m);
        
        /* Destructor */
        ~solver();
        
        int run(int output = 0, bool relax = false, const string& lin_solver = "ma27", double tol = 1e-6, const string& mehrotra = "no");
    };
}
#endif /* defined(__PowerTools____Solver__) */
