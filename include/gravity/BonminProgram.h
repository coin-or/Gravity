#ifndef __Gravity____BonminProgram__
#define __Gravity____BonminProgram__

#ifdef USE_BONMIN
#include <coin/BonTMINLP.hpp>
#include <coin/IpIpoptApplication.hpp>
#include <coin/IpTNLP.hpp>
#endif
#include <gravity/model.h>

using   Ipopt::IpoptApplication;
using namespace  Ipopt;
using namespace Bonmin;
using namespace gravity;

class BonminProgram : public TMINLP, public Program {
public:
    Model* _model = nullptr;
    
    BonminProgram(Model* m);

    /** Pass the type of the variables (INTEGER, BINARY, CONTINUOUS) to the optimizer */
    virtual bool get_variables_types(Index n, VariableType* var_types);

    /** Method to pass some info about the nlp to Ipopt */
    virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                              Index& nnz_h_lag, TNLP::IndexStyleEnum& index_style);

    /** Method to pass the bounds for my problem to Ipopt */
    virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u,
                                 Index m, Number* g_l, Number* g_u);

    /** Method to pass the starting point for the algorithm to Ipopt */
    virtual bool get_starting_point(Index n, bool init_x, Number* x,
                                    bool init_z, Number* z_L, Number* z_U,
                                    Index m, bool init_lambda,

                                    Number* lambda);

    /** Method to return the objective value */
    virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value);

    /** Method to return the gradient of the objective */
    virtual bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f);

    /** Method to return the constraint residuals */
    virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g);

    /** Method to return:
     *   1) The structure of the jacobian (if "values" is NULL)
     *   2) The values of the jacobian (if "values" is not NULL)
     */
    virtual bool eval_jac_g(Index n, const Number* x, bool new_x,
                            Index m, Index nele_jac, Index* iRow, Index *jCol,
                            Number* values);

    /** Method to return:
     *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
     *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
     */
    virtual bool eval_h(Index n, const Number* x, bool new_x,
                        Number obj_factor, Index m, const Number* lambda,
                        bool new_lambda, Index nele_hess, Index* iRow,
                        Index* jCol, Number* values);

    /** Method called by Ipopt */
    virtual void finalize_solution(TMINLP::SolverReturn status,
                                   Index n, const Number* x, Number obj_value);


    virtual bool get_variables_linearity(Index n, TNLP::LinearityType* var_types);
    virtual bool get_constraints_linearity(Index m, TNLP::LinearityType* const_types);

    virtual const SosInfo* sosConstraints() const;
    virtual const BranchingInfo* branchingInfo() const;
};


#endif // __Gravity____BonminProgram__
