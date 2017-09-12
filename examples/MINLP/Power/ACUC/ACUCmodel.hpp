//
//  ACUCmodel.hpp
//  Gravity
//
//  Created by Guanglei Wang on 6/9/17.
//
//
#ifndef ACUCmodel_hpp
#define ACUCmodel_hpp
#include <stdio.h>
#include <gravity/Net.h>
#include <gravity/model.h>
#include <gravity/solver.h>
#ifdef __APPLE__
    #include "/usr/local/include/armadillo"
#elif __linux__
    #include "/usr/include/armadillo"
#endif

typedef enum {ACPOL, ACRECT, QC, QC_SDP, OTS, DF, SOCP, SDP, DC, QC_OTS_L, QC_OTS_N, QC_OTS_O, SOCP_OTS, GRB_TEST } PowerModelType;
typedef enum { MinCost, MinLoss, MinDelta, MaxDelta } Obj;

class PowerModel {
public:
    PowerModelType      _type;
    Model*              _model;
    SolverType*         _solver;
    Net*                _net;
    
    PowerModel();
    PowerModel(PowerModelType type, Net* net, SolverType stype);
    ~PowerModel();
    void reset();
    void build();
    //void run_grb_lin_test();
    //void run_grb_quad_test();
    /** Accessors */
    //Function* objective();

    /** Variables */
    void add_AC_gen_vars();
    void add_AC_Pol_vars();
    void add_AC_Rect_vars();
    void add_QC_vars();
    void add_AC_OTS_vars();
    void add_AC_SOCP_vars();
    void add_DC_vars();
    void add_QC_OTS_vars();
    void add_SOCP_OTS_vars();

    /** Constraints */
    void add_Wr_Wi(Arc* a);
    void add_AC_thermal(Arc* a, bool switch_lines);
    void add_AC_Angle_Bounds(Arc* a, bool switch_lines);
    void add_AC_Voltage_Bounds(Node* n);
    void add_AC_Power_Flow(Arc* a, bool polar);
    void add_AC_KCL(Node* n, bool switch_lines);
    void add_Cosine(Arc* a);
    void add_SDP_cuts(int dim);

    /** Operators */
    int get_order(vector<Node*>* bag, int i);
    
    Function sum(vector<Arc*> vec, string var, bool switch_lines);
    Function sum(vector<Gen*> vec, string var);
    
    /** Models */
    void post_AC_Polar();
    void post_AC_Rect();
    void post_QC();
    void post_AC_SOCP();
    void post_AC_OTS();
    void post_DC();
    void post_QC_OTS(bool lin_cos_cuts, bool quad_cos);
    void post_SOCP_OTS();
    /** Presolve */
    void propagate_bounds();

    /** Solve */
    int solve(int output = 1,bool relax = false);
    void min_cost();
    void min_var(var<>& v);
    void max_var(var<>& v);
    void min_cost_load();
    void print();

    void add_SineCutsPos(Arc *a);
    void add_SineCutsNeg(Arc *a);
};
#endif /* defined(__PowerTools____PowerModel__) */
