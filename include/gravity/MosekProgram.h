//
//  MosekProgram.hpp
//  Gravity
//
//  Created by Guanglei Wang on 14/7/17.
//
//   Note mosek fusion API is designed for CONIC Optimization. 

#ifndef MosekProgram_h
#define MosekProgram_h

#include <stdio.h>
#ifdef USE_MOSEK
#include <fusion.h>
#include <monty.h>
#endif
#include <gravity/model.h>

class MosekProgram {
public:
    Model* _model;
    mosek::fusion::Model::t _mosek_model;
    vector<mosek::fusion::Variable::t> _mosek_vars;
    int _output;
    MosekProgram();
    MosekProgram(Model* m);
    ~MosekProgram();
    void reset_model();
    
    
    bool solve(bool relax);
    void prepare_model();
    void update_model();
    void relax_model();
    
    void fill_in_mosek_vars();
    void create_mosek_constraints();
    void set_mosek_objective();
    void print_constraints();
    
   
};

#endif /* MosekProgram_h */
