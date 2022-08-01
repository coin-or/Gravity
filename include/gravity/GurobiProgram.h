#ifndef __Gravity____GurobiProgram
#define __Gravity____GurobiProgram

#ifdef USE_GUROBI
#include <gurobi_c++.h>
#endif
#include <gravity/model.h>
using namespace gravity;


class GurobiProgram: public Program<>{
private:

public:
    GRBModel* grb_mod;
    GRBEnv* grb_env;
    vector<GRBVar> _grb_vars; /** Mapping variables to Gurobi variables */

    Model<>*    _model;
    int         _output;
    bool grb_first_run=false;
    GurobiProgram();
    GurobiProgram(Model<>* m);
    GurobiProgram(const shared_ptr<Model<>>& m);
    ~GurobiProgram();
    void reset_model();

    bool solve(int output=0, bool relax = false, double tol=1e-6, double mipgap = 0.01, bool gurobi_crossover=false);
    void prepare_model();
    void update_model();
    void update_solution();
    void relax_model();
    void warm_start();
    void update_grb_constraints(std::map<std::string,std::size_t>);
    void initialize_basis(const std::vector<double>& vbasis, const std::vector<double>& cbasis);
    void get_basis(std::vector<double>& vbasis, std::vector<double>& cbasis);
    void initialize_pstart(const std::vector<double>& vbasis, const std::map<std::string, double>& cbasis);
    void get_pstart(std::vector<double>& vbasis, std::map<std::string, double>& cbasis);

    void fill_in_grb_vmap();
    void create_grb_constraints();
    void set_grb_objective();

    void print_constraints();
};


#endif /* defined(__Gravity____GurobiProgram) */
