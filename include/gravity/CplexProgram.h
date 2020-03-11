//
//  CplexProgram.h
//  Gravity
//
//  Created by Hassan Hijazi on 06/06/2017.
//
//

#ifndef CplexProgram_h
#define CplexProgram_h

#include <stdio.h>

#ifndef __Gravity____CplexProgram
#define __Gravity____CplexProgram

#ifdef USE_CPLEX
#include <ilcplex/ilocplex.h>
#endif
#include <gravity/model.h>
//#include <gravity/CplexCallback.h>

using namespace gravity;   

class CplexProgram: public Program<>{
private:
    
    shared_ptr<IloModel> _cplex_model;
    shared_ptr<IloEnv> _cplex_env;
    vector<IloNumVarArray>   _cplex_vars; /** Mapping variables to Cplex variables */
    IloObjective        _cplex_obj;
    
    CPXLONG _cplex_contextmask = 0; /** Context mask for the callback location */
//    CplexCallback _cplex_callback; /** instantiating a callback object */
    
public:
    Model<>* _model;
    int _output;
    CplexProgram();
    CplexProgram(Model<>* m);    
    void update_model();
    void reset_model();
    
    IloModel& get_cplex_model() const{
        return *_cplex_model;
    }
    
    bool solve(int output=0, bool relax=false, double tol=1e-6, double mipgap=0.01);
    void prepare_model();
    void relax_model();
    void warm_start();
    void fill_in_cplex_vars();
    void create_cplex_constraints();
    void set_cplex_objective();
    void create_callback();
    
    void print_constraints();
};


#endif /* defined(__Gravity____CplexProgram) */
#endif /* CplexProgram_hpp */
