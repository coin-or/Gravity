////
////  solver.h
////  Gravity
////
////  Created by Hassan
//
//#ifndef __Gravity____Solver__
//#define __Gravity____Solver__
//
//#include <stdio.h>
//
//#include <gravity/GravityConfig.h>
//
//#include <gravity/model.h>
//#ifdef USE_IPOPT
//#include <gravity/IpoptProgram.h>
//#endif
//#ifdef USE_GUROBI
//#include <gravity/GurobiProgram.h>
//#endif
//#ifdef USE_BONMIN
//#include <gravity/BonminProgram.h>
//#endif
//#ifdef USE_CPLEX
//#include <gravity/CplexProgram.h>
//#endif
//#ifdef USE_CLP
//#include <gravity/ClpProgram.h>
//#endif
////#ifdef USE_SDPA
////#include "SdpaProgram.h"
////#endif
//#ifdef USE_MOSEK
//#include "MosekProgram.h"
//#endif
//
//namespace gravity {
//    class solver {
//    public:
//        gravity::Model*                         _model = nullptr;
//        Program*                       _prog = nullptr;
//        SolverType                     _stype;        
//        double                         _tol = 1e-6; /*<< Solver tolerance. */
//        
//        /** Constructor */
//        //@{
//        solver();
//
//        solver(gravity::Model& model, SolverType stype);
//        //@}
//        void set_model(gravity::Model& m);
//        
//        /* Destructor */
//        ~solver();
//        /* Options: int output = 0, bool relax = false, double tol = 1e-6, double mipgap = 1e-3, const string& lin_solver = "ma57", const string& mehrotra = "no",  int max_iter = 1000 */
//        int run(int output = 0, bool relax = false, double tol = 1e-6, double mipgap = 1e-3, const string& lin_solver = "mumps", const string& mehrotra = "no",  int max_iter = 1000);
//        
//    };
//    
//    void run_parallel(const vector<shared_ptr<gravity::Model>>& models, gravity::SolverType stype, double tol, unsigned nr_threads); /** < Runds models stored in the vector in parallel, using solver of stype and tolerance tol */
//
//}
//
//#endif /* defined(__Gravity____Solver__) */
