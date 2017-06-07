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
#include <Gravity/model.h>

class CplexProgram {
private:
    
    IloModel* _cplex_model;
    IloEnv* _cplex_env;
    vector<IloNumVar> _cplex_vars; /** Mapping variables to Cplex variables */
public:
    Model* _model;
    int _output;
    CplexProgram();
    CplexProgram(Model* m);
    ~CplexProgram();
    void reset_model();
    
    bool solve(bool relax);
    void prepare_model();
    void update_model();
    void relax_model();
    
    void fill_in_cplex_vars();
    void create_cplex_constraints();
    void set_cplex_objective();
    
    void print_constraints();
};


#endif /* defined(__PowerTools____CplexProgram) */








#endif /* CplexProgram_hpp */
