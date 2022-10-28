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
    Model<>* _model = nullptr;
    Worker(Model<>* m): _model(m){}
    ~Worker()
    {
       // Free all memory associated to env
       env.end();
    }
    // This routine separates cuts violated by the current x solution.
    IloBool separate(const IloNumVarArray& x, const IloNumArray& xSol, IloExpr& cutLhs, IloNum& cutRhs, double tol)
    {
        IloBool violatedCutFound = IloFalse;
        cutLhs.clear();
        cutRhs = 0.0;

       // A violated cut is available iff the solution status is feasible
        if ( cplex.getStatus() == IloAlgorithm::Feasible ) {
            
            IloNumVarArray var(env);
            IloNumArray val(env);
            
            // Get the most violated cut
            auto most_viol = _model->get_most_violated(tol);
            size_t idx = 0, idx_inst = 0, idx1 = 0, idx2 = 0, idx_inst1 = 0, idx_inst2 = 0, nb_inst = 0;
            double cval = 0;
            auto c = most_viol.first;// Most violated symbolic constraint
            auto i = most_viol.second;// Instance of most violated constraint
            if(!c)// No violated cuts found
                return IloFalse;
            IloNumExpr cc(env);
            for (auto& it_qterm: c->get_qterms()) {
                IloNumExpr qterm(env);
                idx1 = it_qterm.second._p->first->get_vec_id();
                idx2 = it_qterm.second._p->second->get_vec_id();
                if (it_qterm.second._p->first->_is_vector) {
                    auto dim = it_qterm.second._p->first->get_dim(i);
                    for (size_t j = 0; j<dim; j++) {
                        qterm += c->eval(it_qterm.second._coef,i,j)*_cplex_vars->at(idx1)[it_qterm.second._p->first->get_id_inst(i,j)]*_cplex_vars->at(idx2)[it_qterm.second._p->second->get_id_inst(i,j)];
                    }
                }
                else {
                    idx_inst1 = it_qterm.second._p->first->get_id_inst(i);
                    idx_inst2 = it_qterm.second._p->second->get_id_inst(i);
                    qterm += c->eval(it_qterm.second._coef, i)*_cplex_vars->at(idx1)[idx_inst1]*_cplex_vars->at(idx2)[idx_inst2];
                }
                if (!it_qterm.second._sign) {
                    qterm *= -1;
                }
                cc += qterm;
                qterm.end();
            }
            
            for (auto& it_lterm: c->get_lterms()) {
                IloNumExpr lterm(env);
                idx = it_lterm.second._p->get_vec_id();
                if (it_lterm.second._p->_is_vector || it_lterm.second._p->is_matrix_indexed() || it_lterm.second._coef->is_matrix()) {
                    auto dim = it_lterm.second._p->get_dim(i);
                    for (int j = 0; j<dim; j++) {
                        lterm += c->eval(it_lterm.second._coef,i,j)*_cplex_vars->at(idx)[it_lterm.second._p->get_id_inst(i,j)];
                    }
                }
                else {
                    idx_inst = it_lterm.second._p->get_id_inst(i);
                    lterm += c->eval(it_lterm.second._coef, i)*_cplex_vars->at(idx)[idx_inst];
                }
                if (!it_lterm.second._sign) {
                    lterm *= -1;
                }
                cc += lterm;
                lterm.end();
            }
            cc += c->eval(c->get_cst(), i);
            
            cutLhs += cc;
            if(c->get_ctype()==leq) {
                cutLhs *= -1;
            }
            var.end();
            val.end();
            return violatedCutFound;
        }
        return IloFalse;
    } // END separate
 private:
    IloEnv env;
    IloCplex cplex;
    shared_ptr<vector<IloNumVarArray>>          _cplex_vars; /*< Mapping variables to Cplex variables */
 };

class CplexCallback: public IloCplex::Callback::Function {
private:
    /* Variables to be considered in the callback. */
    std::vector<Worker*>                        workers;
    shared_ptr<vector<IloNumVarArray>>          _cplex_vars; /*< Mapping variables to Cplex variables */
public:
    Model<>* _model;
    CplexCallback(IloInt numWorkers=1, Model<>* m=nullptr, shared_ptr<vector<IloNumVarArray>> cvars = nullptr)
       : workers(numWorkers, 0), _model(m), _cplex_vars(cvars){}
    
   
    void invoke(const IloCplex::Callback::Context& context) ILO_OVERRIDE
    {
       int const threadNo = context.getIntInfo(IloCplex::Callback::Context::Info::ThreadId);

       // setup
       if (context.inThreadUp()) {
          delete workers[threadNo];
          workers[threadNo] = new Worker(_model);
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
        IloNumArray xSol;

       // Get the current x solution
        IloNumVarArray x(env,_cplex_vars->size());
       switch (context.getId()) {
       case IloCplex::Callback::Context::Id::Candidate:
          if ( !context.isCandidatePoint() ) // The model is always bounded
             throw IloCplex::Exception(-1, "Unbounded solution");
             xSol = IloNumArray(env,_cplex_vars->size());
             context.getCandidatePoint(x, xSol);
          break;
       case IloCplex::Callback::Context::Id::Relaxation:
         xSol = IloNumArray(env,_cplex_vars->size());
         context.getRelaxationPoint(x, xSol);
          break;
       default:
          // Free memory
          xSol.end();
          throw IloCplex::Exception(-1, "Unexpected contextID");
       }

       // Get the right worker
       Worker* worker = workers[threadNo];

       // Separate cut
       IloExpr cutLhs(env);
       IloNum cutRhs;
       IloBool sepStat = worker->separate(x, xSol, cutLhs, cutRhs, 1e-6);

       // Free memory
       xSol.end();

       if (sepStat) {
          // Add the cut
          IloRange r(env, cutRhs, cutLhs, IloInfinity);/* Adds constraint: infinity >= cutLhs >= cutRhs  */

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
