#ifndef __Gravity____HiGHSProgram
#define __Gravity____HiGHSProgram

#ifdef USE_HiGHS
#include <Highs.h>
#endif
#include <gravity/model.h>
using namespace gravity;


class HiGHSProgram: public Program<>{
private:

    HighsModel          Highs_mod;/**< Highs model */
    Highs               Highs_inst;/**< Highs instance */
    
public:
    Model<>*            _model;
    int                 _output;
    HiGHSProgram();
    HiGHSProgram(Model<>* m);
    HiGHSProgram(const shared_ptr<Model<>>& m);
    ~HiGHSProgram();
    void reset_model();

    bool solve(bool relax = false, double mipgap = 0.01);
    void prepare_model();
    void update_model();
    void update_solution();
    void relax_model();

    void fill_in_var_map();
    void create_constraints();
    void set_objective();

    void print_constraints();
};


#endif /* defined(__Gravity____HiGHSProgram) */
