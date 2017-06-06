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

class CPLEXProgram {
private:
    
    GRBModel* grb_mod;
    GRBEnv* grb_env;
    std::map<string, GRBVar*> _grb_vars; /** Mapping variables to CPLEX variables */
public:
    Model *model;
    int _output;
    CPLEXProgram();
    CPLEXProgram(Model* m);
    ~CPLEXProgram();
    void reset_model();
    
    bool solve(bool relax);
    void prepare_model();
    void update_model();
    void relax_model();
    
    void fill_in_grb_vmap();
    void create_grb_constraints();
    void set_grb_objective();
    
    void print_constraints();
};


#endif /* defined(__PowerTools____CPLEXProgram) */








#endif /* CplexProgram_hpp */
