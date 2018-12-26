////
////  ClpProgram.hpp
////  
////
////  Created by guanglei Wang on 2018/12/6.
////
//
//#ifndef ClpProgram_hpp
//#define ClpProgram_hpp
//
//#include <stdio.h>
//
//#ifdef USE_CLP
//#include "ClpSimplex.hpp"
//#include "CoinBuild.hpp"
//#endif
//#include <gravity/model.h>
//
//
//using namespace gravity; 
//
//class ClpProgram: public Program{
//  private:
//  ClpModel* _clp = new ClpModel();
//
//  public:
//  enum SolveType {SOLVE_PRIMAL = 0, SOLVE_DUAL = 1};
//    Model* _model;
//    int _output;
//    ClpProgram();
//    ClpProgram(Model* m);
//    ~ClpProgram();
//    void update_model();
//    void reset_model();
//  
//    bool solve(SolveType type = SOLVE_PRIMAL);
//    void prepare_model();    
//    void relax_model();
//    
//    void fill_in_clp_vars();
//    void create_clp_constraints();
//    void set_clp_objective();
//    
//    void print_constraints();
//};
//#endif /* ClpProgram_hpp */
