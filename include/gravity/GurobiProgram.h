#ifndef __Gravity____GurobiProgram
#define __Gravity____GurobiProgram

#ifdef USE_GUROBI
#include <gurobi_c++.h>
#endif
#include <gravity/model.h>
using namespace gravity;


class GurobiProgram: public Program<>{
private:

    
    
    vector<GRBVar> _grb_vars; /** Mapping variables to Gurobi variables */
public:
    GRBEnv* grb_env;
    GRBModel* grb_mod;
    Model<>*    _model;
    int         _output;
    GurobiProgram();
    GurobiProgram(Model<>* m);
    GurobiProgram(const shared_ptr<Model<>>& m);
    ~GurobiProgram();
    void reset_model();
    bool solve(bool relax = false, double mipgap = 0.01, double tim_limit = 3600);
    void prepare_model();
    void update_model();
    void update_solution();
    void relax_model();
    void unrelax_model();
    void fill_in_grb_vmap();
    void create_grb_constraints();
    void write_NLCstr(const string &fname);
    void set_grb_objective();
    void print_constraints();
};


#endif /* defined(__Gravity____GurobiProgram) */
