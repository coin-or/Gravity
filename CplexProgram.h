//
//  CplexProgram.h
//  Gravity
//
//  Created by Guanglei Wang on 06/06/2017.
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
using namespace gravity;

class CplexProgram: public Program{
private:
    
    IloModel* _cplex_model;
    IloEnv* _cplex_env;
    vector<IloNumVarArray>   _cplex_vars; /** Mapping variables to Cplex variables */
    IloObjective        _cplex_obj;
public:
    Model* _model;
    int _output;
    CplexProgram();
    CplexProgram(Model* m);
    ~CplexProgram();
    void update_model();
    void reset_model();
    
    IloModel& get_cplex_model() const{
        return *_cplex_model;
    }
    
    bool solve(bool relax=false, double mipgap = 0.01);
    void prepare_model();    
    void relax_model();
    
    void fill_in_cplex_vars();
    void create_cplex_constraints();
    void set_cplex_objective();
    
    void print_constraints();
};


#endif /* defined(__Gravity____CplexProgram) */
#endif /* CplexProgram_hpp */
