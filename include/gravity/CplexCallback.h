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
    IloEnv* env = nullptr;
    Model<>* _model = nullptr;
    Worker(Model<>* m, shared_ptr<vector<IloNumVarArray>> cvars = nullptr): _model(m), _cplex_vars(cvars){};
    ~Worker(){};
    // This routine separates cuts violated by the current solution.
    IloBool separate(IloExpr& cutLhs, IloNum& cutRhs, double tol)
    {
        IloBool violatedCutFound = IloFalse;
        cutLhs.clear();
        cutRhs = 0.0;
        
        // Get the most violated cut
        auto most_viol = _model->get_most_violated(tol);
        size_t idx = 0, idx_inst = 0, idx1 = 0, idx2 = 0, idx_inst1 = 0, idx_inst2 = 0, nb_inst = 0;
        double cval = 0;
        auto c = most_viol.first;// Most violated symbolic constraint
        auto i = most_viol.second;// Instance of most violated constraint
        if(!c)// No violated cuts found
            return IloFalse;
        IloNumExpr cc(*env);
        for (auto& it_lterm: c->get_lterms()) {
            IloNumExpr lterm(*env);
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
        return IloTrue;
    } // END separate
private:
    shared_ptr<vector<IloNumVarArray>>          _cplex_vars; /**< Cplex variables */
};

class CplexCallback: public IloCplex::Callback::Function {
private:
    /* Variables to be considered in the callback. */
    std::vector<Worker*>                        workers; /**< Cplex Workers */
    shared_ptr<vector<IloNumVarArray>>          _cplex_vars; /**< Cplex variables */
public:
    vector<Model<>*> _models;/**< vector storing nb_threads copies of the original model for parallel cut generation */
    CplexCallback(IloInt numWorkers=1, Model<>* m=nullptr, shared_ptr<vector<IloNumVarArray>> cvars = nullptr)
    : workers(numWorkers, nullptr), _models(numWorkers, nullptr), _cplex_vars(cvars){
        _models[0] = m;
        for(auto i = 1; i<numWorkers; i++){
            _models[i] = new Model<>(*m);
        }
    }
    
    
    void invoke(const IloCplex::Callback::Context& context) ILO_OVERRIDE
    {
        int const threadNo = context.getIntInfo(IloCplex::Callback::Context::Info::ThreadId);
        
        // setup
        if (context.inThreadUp()) {
            delete workers[threadNo];
            workers[threadNo] = new Worker(_models[threadNo],_cplex_vars);
            return;
        }
        
        // teardown
        if (context.inThreadDown()) {
            delete workers[threadNo];
            workers[threadNo] = 0;
            return;
        }
        IloInt numNodes = 0;
        IloEnv env = context.getEnv();
        /** Update the model with current solution
         @note: This should be thread safe since we have multiple copies of the same model.
         */
        switch (context.getId()) {
            case IloCplex::Callback::Context::Id::Candidate:
                if ( !context.isCandidatePoint() ) // The model is unbounded
                    throw IloCplex::Exception(-1, "Unbounded problem");
                for (auto i = 0; i < _cplex_vars->size(); i++) {
                    for (auto j = 0; j < _models[threadNo]->_vars[i]->get_dim(); j++) {
                        _models[threadNo]->_vars[i]->set_double_val(j,context.getCandidatePoint(_cplex_vars->at(i)[j]));
                    }
                }
                break;
            case IloCplex::Callback::Context::Id::Relaxation:
                for (auto i = 0; i < _cplex_vars->size(); i++) {
                    for (auto j = 0; j < _models[threadNo]->_vars[i]->get_dim(); j++) {
                        _models[threadNo]->_vars[i]->set_double_val(j,context.getRelaxationPoint(_cplex_vars->at(i)[j]));
                    }
                }
                break;
            default:
                // Free memory
                throw IloCplex::Exception(-1, "Unexpected contextID");
        }
        
        // Get the right worker
        Worker* worker = workers[threadNo];
        worker->env = &env;
        // Separate cut
        IloExpr cutLhs(env);
        IloNum cutRhs;
        IloBool sepStat = worker->separate(cutLhs, cutRhs, 1e-6);
        
        
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
