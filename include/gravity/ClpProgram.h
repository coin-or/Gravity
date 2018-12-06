//
//  ClpProgram.hpp
//  option_parser
//
//  Created by 王光磊 on 2018/12/6.
//

#ifndef ClpProgram_hpp
#define ClpProgram_hpp

#include <stdio.h>

#include <gravity/model.h>
#include "ClpSimplex.hpp"


using namespace gravity; 

class ClpProgram: public Program{
  private:
  ClpModel clp;

  public:
  enum solveType {Simplex, Dual};
    Model* _model;
    int _output;
    ClpProgram();
    ClpProgram(Model* m);
    ~ClpProgram();
    void update_model();
    void reset_model();
  
    bool solve(solveType type,bool relax=false, double mipgap = 0.01);
    void prepare_model();    
    void relax_model();
    
    void fill_in_cplex_vars();
    void create_cplex_constraints();
    void set_cplex_objective();
    
    void print_constraints();
};
#endif /* ClpProgram_hpp */
