//
//  IpoptProgram.h
//

#ifndef __Gravity____IpoptProgram__
#define __Gravity____IpoptProgram__

#include <stdio.h>
#include <assert.h>
#ifdef USE_IPOPT
#include <coin/IpIpoptApplication.hpp>
#include <coin/IpTNLP.hpp>
#endif
#include <gravity/model.h>

using   Ipopt::IpoptApplication;

using namespace Ipopt;
using namespace gravity;

class IpoptProgram : public TNLP, public Program{

    
public:
    Model* _model = nullptr;
    
    IpoptProgram(Model* m);
    
    void update_model();

    /** Method to return some info about the nlp */
    virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                              Index& nnz_h_lag, IndexStyleEnum& index_style);

    /** Method to return the bounds for my problem */
    virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u,
                                 Index m, Number* g_l, Number* g_u);
    
    /** Method to return the starting point for the algorithm */
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
    
    virtual void finalize_solution(
                                   SolverReturn               status    ,
                                   Index                             n         ,
                                   const Number*                     x         ,
                                   const Number*                     z_L       ,
                                   const Number*                     z_U       ,
                                   Index                             m         ,
                                   const Number*                     g         ,
                                   const Number*                     lambda    ,
                                   Number                            obj_value ,
                                   const IpoptData*           ip_data   ,
                                   IpoptCalculatedQuantities* ip_cq
                                   );
    
    virtual bool get_variables_linearity(Index n, LinearityType* var_types);
    virtual bool get_constraints_linearity(Index m, LinearityType* const_types);
    
};



#endif /* defined(__Gravity____IpoptProgram__) */
