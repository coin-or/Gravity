#ifndef __Gravity____GurobiProgram
#define __Gravity____GurobiProgram

#ifdef USE_GUROBI
#include <gurobi_c++.h>
#endif
#include <gravity/model.h>
using namespace gravity;


class GurobiProgram {
private:

    GRBModel* grb_mod;
    GRBEnv* grb_env;
    vector<GRBVar> _grb_vars; /** Mapping variables to Gurobi variables */
public:
    Model*  _model;
    int     _output;
    GurobiProgram();
    GurobiProgram(Model* m);
    ~GurobiProgram();
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


#endif /* defined(__Gravity____GurobiProgram) */
