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

#ifndef __PowerTools____CplexProgram
#define __PowerTools____CplexProgram

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
    vector<IloRangeArray>  _cplex_constraints; /** useful for retrieving dual mulipliers */
   /*  CPLEX does not give us the dual multipliers for quadratic        *
    *   constraints directly. This is because they may not be properly  *
    *   defined at the cone top and deciding whether we are at the cone *
    *   top or not involves (problem specific) tolerance issues. CPLEX  *
    *   instead gives us all the values we need in order to compute the *
    *   dual multipliers if we are not at the cone top.                 */
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
    
    bool solve(bool relax);
    void prepare_model();    
    void relax_model();
    
    void fill_in_cplex_vars();
    void create_cplex_constraints();
    void set_cplex_objective();
    
    void print_constraints();
};


#endif /* defined(__PowerTools____CplexProgram) */
#endif /* CplexProgram_hpp */
