//
//  model.hpp
//  Gravity
//
//  Created by Guanglei Wang on 19/6/17.
//
//

#ifndef Minkmodel_hpp
#define Minkmodel_hpp

#include <stdio.h>
#include <gravity/Net.h>
#include <gravity/model.h>
#include <gravity/solver.h>
#include <stdlib.h>
#include "armadillo"

typedef enum {MIP,MIP_tree,SDP,SDP_tree, Node_edge} ModelType;

class Minkmodel {
public:
    ModelType   _type;
    gravity::Model       _model;
    gravity::SolverType  _solver;
    double      _K;
    Net*        _graph;
    Net*        _chordal_extension;

    set<tuple<int,int,int>> _ids; 
    shared_ptr<map<std::string,vector<unsigned>>> _cliqueid;
    set<tuple<int,int,int,int>> _ids4; //
    
    arma::mat _eigvec;
    Minkmodel();
    Minkmodel(ModelType type, Net* graph, double K);
    Minkmodel(ModelType type, Net* graph, double K, gravity::SolverType solver);

    ~Minkmodel();
    void reset();
    void build();
    
    /** Variables */
    void add_vars_origin();
    void add_vars_origin_tree();
    void add_vars_lifted();
    void add_vars_lifted_tree();
    
    /** Constraints */
    // different formulations
    //void add_obj();  included in add_vars_*
    //void add_obj_lifted();
    void add_triangle();
    void add_clique();
    void add_clique_tree();
    void add_triangle_lifted();
    void add_triangle_tree();
    void add_triangle_lifted_tree();
    void add_clique_lifted();
    void add_clique_lifted_tree();
    void add_general_clique();
    void add_wheel();
    void add_bicycle();
    void add_3Dcuts();
    void cliquetree_decompose();
    void nchoosek(int, int, int);
    void node_edge_formulation();
    //  post root node relaxation
    bool check_eigenvalues();
    void add_eigcut();
    void construct_fsol();
    
    /** Presolve */
    //    void propagate_bounds();
    /** Solve */
    int solve(int output,bool relax);
    void print();
};
#endif /* model_hpp */
