//
//  solver.h
//  PowerTools++
//
//  Created by Hassan on 30/01/2015.
//  Copyright (c) 2015 NICTA. All rights reserved.
//

#ifndef __PowerTools____Solver__
#define __PowerTools____Solver__

#include <stdio.h>

#include <Gravity/GravityConfig.h>

#include <Gravity/model.h>
#ifdef USE_IPOPT
#include <Gravity/IpoptProgram.h>
#endif
#ifdef USE_GUROBI
#include <Gravity/GurobiProgram.h>
#endif
#ifdef USE_BONMIN
#include <Gravity/BonminProgram.h>
#endif
#ifdef USE_CPLEX
#include <Gravity/CplexProgram.h>
#endif
#ifdef USE_SDPA
#include "SdpaProgram.h"
#endif

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
    
    int run(int output = 0, bool relax = false);
};
#endif /* defined(__PowerTools____Solver__) */
