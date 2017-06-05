#ifndef __PowerTools____GurobiProgram
#define __PowerTools____GurobiProgram

#ifdef USE_GUROBI
#include <gurobi_c++.h>
#endif
#include <Gravity/model.h>

class GurobiProgram {
private:

    GRBModel* grb_mod;
    GRBEnv* grb_env;
    std::map<string, GRBVar*> _grb_vars; /** Mapping variables to Gurobi variables */
public:
    Model *model;
    int _output;
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


#endif /* defined(__PowerTools____GurobiProgram) */
