//
//  CplexProgram.h
//  Gravity
//
//  Created by Hassan Hijazi on Oct 4 2022.
//
//

#ifndef CplexCallback_h
#define CplexCallback_h

#include <stdio.h>

#ifndef __Gravity____CplexCallback
#define __Gravity____CplexCallback

#ifdef USE_CPLEX
#include <ilcplex/ilocplex.h>
#endif
#include <gravity/model.h>
using namespace gravity;
class Worker
{
public:
    Worker(){}
    ~Worker()
    {
       // Free all memory associated to env
       env.end();
    }
    // This routine separates cuts violated by the current x solution.
    IloBool separate(const IloArray<IloIntVarArray>& x, const IloNumArray2& xSol, IloExpr& cutLhs, IloNum& cutRhs)
    {
       IloBool violatedCutFound = IloFalse;
       

       // A violated cut is available iff the solution status is Unbounded
       if ( cplex.getStatus() == IloAlgorithm::Unbounded ) {

          IloNumVarArray var(env);
          IloNumArray val(env);

          // Get the violated cut as an unbounded ray of the worker LP
          

          // Compute the cut from the violated matrices
          cutLhs.clear();
          cutRhs = 0.0;

          for (IloInt h = 0; h < val.getSize(); ++h) {
             IloInt *index_p = (IloInt*) var[h].getObject();
             IloInt index = *index_p;
//            cutRhs -= val[h];
//            cutLhs += val[h] * x[i][j];
          }
          var.end();
          val.end();
          violatedCutFound = IloTrue;
       }

       return violatedCutFound;
    } // END separate
 private:
    IloEnv env;
    IloCplex cplex;
    IloNumVarArray v;
    IloIntArray vIndex;
    IloNumVarArray u;
    IloIntArray uIndex;
 };

class CplexCallback: public IloCplex::Callback::Function {
private:
    /* Variables to be considered in the callback. */
    IloNumVarArray callback_vars;
    std::vector<Worker*> workers;
    
public:
    
    CplexCallback(IloInt numWorkers=1)
       : workers(numWorkers, 0) {}
    
   
    void invoke(const IloCplex::Callback::Context& context) ILO_OVERRIDE
    {
       int const threadNo = context.getIntInfo(IloCplex::Callback::Context::Info::ThreadId);

       // setup
       if (context.inThreadUp()) {
          delete workers[threadNo];
          workers[threadNo] = new Worker();
          return;
       }

       // teardown
       if (context.inThreadDown()) {
          delete workers[threadNo];
          workers[threadNo] = 0;
          return;
       }
        IloInt numNodes = 0;//_cplex_vars.size()
        IloEnv env = context.getEnv();
        IloNumArray2 xSol(env, numNodes);

       // Get the current x solution
        IloArray<IloIntVarArray> x;
       switch (context.getId()) {
       case IloCplex::Callback::Context::Id::Candidate:
          if ( !context.isCandidatePoint() ) // The model is always bounded
             throw IloCplex::Exception(-1, "Unbounded solution");
          for (IloInt i = 0; i < numNodes; ++i) {
             xSol[i] = IloNumArray(env);
             context.getCandidatePoint(x[i], xSol[i]);
          }
          break;
       case IloCplex::Callback::Context::Id::Relaxation:
          for (IloInt i = 0; i < numNodes; ++i) {
             xSol[i] = IloNumArray(env);
             context.getRelaxationPoint(x[i], xSol[i]);
          }
          break;
       default:
          // Free memory
          for (IloInt i = 0; i < numNodes; ++i) xSol[i].end();
          xSol.end();
          throw IloCplex::Exception(-1, "Unexpected contextID");
       }

       // Get the right worker
       Worker* worker = workers[threadNo];

       // Separate cut
       IloExpr cutLhs(env);
       IloNum cutRhs;
       IloBool sepStat = worker->separate(x, xSol, cutLhs, cutRhs);

       // Free memory
       for (IloInt i = 0; i < numNodes; ++i) xSol[i].end();
       xSol.end();

       if (sepStat) {
          // Add the cut
          IloRange r(env, cutRhs, cutLhs, IloInfinity);

          switch (context.getId()) {
             case IloCplex::Callback::Context::Id::Candidate:
                context.rejectCandidate(r);
                break;
             case IloCplex::Callback::Context::Id::Relaxation:
                context.addUserCut(r,
                                   IloCplex::UseCutPurge,
                                   IloFalse);
                break;
             default:
                r.end();
                throw IloCplex::Exception(-1, "Unexpected contextID");
          }
          r.end();
       }
    }
    
    /// Destructor
    virtual ~CplexCallback();
    
};

#endif /* defined(__Gravity____CplexCallback) */
#endif /* CplexCallback_hpp */
