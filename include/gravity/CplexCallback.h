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
#ifdef USE_EIGEN
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
using namespace Eigen;
#endif

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
    IloBool separate(IloExprArray& cutLhs, double tol)
    {
        IloBool violatedCutFound = IloFalse;
#ifdef USE_EIGEN
        SelfAdjointEigenSolver<MatrixXd> es;
        es.compute(mat_full);
#endif
        // Get the most violated cuts
        auto violated_cstr = _model->sort_violated_constraints(tol);
        int cstr_id = 0;
        for(cstr_id = 0; cstr_id<std::min(violated_cstr.size(),(size_t)cutLhs.getSize()); cstr_id++){
            auto most_viol = violated_cstr[cstr_id];
            size_t idx = 0, idx_inst = 0, idx1 = 0, idx2 = 0, idx_inst1 = 0, idx_inst2 = 0, nb_inst = 0;
            double cval = 0;
            auto c = get<2>(most_viol);// Most violated symbolic constraint
            auto i = get<3>(most_viol);// Instance of most violated constraint
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
            
            cutLhs[cstr_id] = cc;
            if(c->get_ctype()==leq) {
                cutLhs[cstr_id] *= -1;
            }
            violatedCutFound = true;
            DebugOff("Added " << cstr_id << " lazy cuts\n");
        }
        
        /* SDP CUTS */
        bool hierarc = false, add_bags = true;
        if(add_bags){
            int nb_cuts = 0;
            auto res1=_model->cuts_eigen_bags(1e-6);
            if(res1.size()>=1){
                DebugOff("Adding " << res1.size() << " bag cuts\n");
                for(auto i=0;i<res1.size() && cstr_id<cutLhs.getSize();i++){
                    IloNumExpr cc(*env);
                    size_t j = 0;
                    for(j=0;j<res1[i].size()-1;j++){
                        auto c=res1[i][j];
                        cc -= c.second*_cplex_vars->at(c.first.first)[c.first.second];
                    }
                    cc -= res1[i][j].second;
                    cutLhs[cstr_id++] = cc;
                    nb_cuts++;
                }
                violatedCutFound = true;
            }
            if(res1.size()>=1 && hierarc){
                cutLhs.setSize(cstr_id);
                return IloTrue;
            }
        }
        auto res2=_model->cuts_eigen_full(1e-6);
        if(res2.size()>=1){
            DebugOff("Adding " << res2.size() << " full cuts\n");
            for(auto i=0;i<res2.size() && cstr_id<cutLhs.getSize();i++){
                IloNumExpr cc(*env);
                size_t j = 0;
                for(j=0;j<res2[i].size()-1;j++){
                    auto c=res2[i][j];
                    cc -= c.second*_cplex_vars->at(c.first.first)[c.first.second];
                }
                cc -= res2[i][j].second;
                cutLhs[cstr_id++] = cc;
            }
            violatedCutFound = true;
        }
        DebugOff("Added " << cstr_id << " total cuts\n");
        cutLhs.setSize(cstr_id);
        return violatedCutFound;
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
    int nb_cuts = 10000; /**< @brief number of violated callback constraints to generate each iteration */
    CplexCallback(IloInt numWorkers=1, Model<>* m=nullptr, shared_ptr<vector<IloNumVarArray>> cvars = nullptr)
    : workers(numWorkers, nullptr), _models(numWorkers, nullptr), _cplex_vars(cvars){
        _models[0] = m;
        _models[0]->num_cuts.resize(9,0);
        for(auto i = 1; i<numWorkers; i++){
            _models[i] = new Model<>(*m);
            _models[i]->num_cuts.resize(9,0);
        }
    }
    /**
     @brief Set the number of violated callback constraints to generate each iteration
     */
    void set_nb_cuts(int nb){
        nb_cuts = nb;
    };
    
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
//                _models[threadNo]->print_solution();
//                _models[threadNo]->print_int_solution();
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
        IloExprArray cutLhs(env, nb_cuts);
        IloBool sepStat = worker->separate(cutLhs, 1e-6);
        
        
        if (sepStat) {
            // Add the cut
            IloRangeArray r(env, 0.0, cutLhs, IloInfinity);/* Adds constraint: infinity >= cutLhs >= cutRhs  */
            
            switch (context.getId()) {
                case IloCplex::Callback::Context::Id::Candidate:
                    context.rejectCandidate(r);
                    break;
                case IloCplex::Callback::Context::Id::Relaxation:{
                    for(int cstr_id = 0; cstr_id<r.getSize(); cstr_id++){
                        context.addUserCut(r[cstr_id],
                                           IloCplex::UseCutForce,
                                           IloFalse);
                    }
                    break;}
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
