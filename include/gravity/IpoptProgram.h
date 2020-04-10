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

template<typename type = double>
class IpoptProgram : public TNLP, public Program<type>{

    
public:
    Model<type>* _model = nullptr;
    
    IpoptProgram(Model<type>* m):_model(m){
    }
    
    
    
    
    void update_model(){
        _model->fill_in_maps();
    }
    
    
    
    
    /** Method to return some info about the nlp */
    bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                                    Index& nnz_h_lag, Ipopt::TNLP::IndexStyleEnum& index_style){
        index_style = Ipopt::TNLP::C_STYLE;
        n = (Index)_model->get_nb_vars();
        m = (Index)_model->get_nb_cons();
        nnz_jac_g = (Index)_model->get_nb_nnz_g();
        DebugOn("Number of nnz in Jacobian = " << nnz_jac_g << endl);
        _model->_jac_vals.resize(nnz_jac_g,0);
        nnz_h_lag = (Index)_model->get_nb_nnz_h();
        DebugOn("Number of nnz in Hessian = " << nnz_h_lag << endl);
        _model->_first_call_jac = true;
        //If quadratic model and we're resolving no need to reset these
        _model->_first_call_hess = true;
        _model->_first_call_grad_obj = true;
        return true;
    }
    
    void finalize_solution(Ipopt::SolverReturn             status    ,
                                         Index                             n         ,
                                         const Number*                     x         ,
                                         const Number*                     z_L       ,
                                         const Number*                     z_U       ,
                                         Index                             m         ,
                                         const Number*                     g         ,
                                         const Number*                     lambda    ,
                                         Number                            obj_value ,
                                         const Ipopt::IpoptData*           ip_data   ,
                                         Ipopt::IpoptCalculatedQuantities* ip_cq
                                         )
    {
        _model->set_x(x);
        for (auto &cp: _model->_cons_name) {
            auto nb_ins = cp.second->get_nb_inst();
            cp.second->_dual.resize(nb_ins);
            auto idx = 0;
            for (size_t inst = 0; inst < nb_ins; inst++) {
                if (!*cp.second->_all_lazy || !cp.second->_lazy[inst]) {
                    cp.second->_dual[inst] = lambda[cp.second->_id + idx++];
                }
            }
        }
        for (auto &vp: _model->_vars) {
            auto nb_inst = vp.second->get_dim();
            vp.second->_u_dual.resize(nb_inst);
            vp.second->_l_dual.resize(nb_inst);
            auto vid = vp.second->get_id();
            for (size_t inst = 0; inst < nb_inst; inst++) {
                vp.second->_u_dual[inst] = z_U[vid + inst];
                vp.second->_l_dual[inst] = z_L[vid + inst];
            }
        }
    }
    
    /** Method to return the bounds for my problem */
    bool get_bounds_info(Index n, Number* x_l, Number* x_u,
                                       Index m, Number* g_l, Number* g_u){
        assert(n==_model->get_nb_vars());
        assert(m==_model->get_nb_cons());
        _model->fill_in_var_bounds(x_l , x_u);
        _model->fill_in_cstr_bounds(g_l , g_u);
        
        return true;
    }
    
    bool get_starting_point(Index n, bool init_x, Number* x,
                             bool init_z, Number* z_L, Number* z_U,
                             Index m, bool init_lambda,
                             Number* lambda){
        return get_starting_point_(n, init_x, x,init_z, z_L, z_U,m, init_lambda,lambda);
    }
    
    /** Method to return the starting point for the algorithm */
    template<typename T=type,typename std::enable_if<is_arithmetic<T>::value>::type* = nullptr>
    bool get_starting_point_(Index n, bool init_x, Number* x,
                                          bool init_z, Number* z_L, Number* z_U,
                                          Index m, bool init_lambda,
                                          Number* lambda){
        assert(n==_model->get_nb_vars());
        assert(m==_model->get_nb_cons());
        
        if (init_x) {
            _model->fill_in_var_init(x);
        }
        if (init_lambda && init_z) {
            _model->fill_in_duals(lambda,z_L,z_U);
        }
        return true;
    }
    
    bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value){
        return eval_f_(n, x, new_x, obj_value);
    }
    
    /** Method to return the objective value */
    template<typename T=type,typename std::enable_if<is_arithmetic<T>::value>::type* = nullptr>
    bool eval_f_(Index n, const Number* x, bool new_x, Number& obj_value){
        
        assert(n==_model->get_nb_vars());
        _model->fill_in_obj(x, obj_value,new_x);
        return true;
    }
    
    bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f){
        return eval_grad_f_(n, x, new_x, grad_f);
    }
    
    /** Method to return the gradient of the objective */
    template<typename T=type,typename std::enable_if<is_arithmetic<T>::value>::type* = nullptr>
    bool eval_grad_f_(Index n, const Number* x, bool new_x, Number* grad_f){
        
        assert(n==_model->get_nb_vars());
        _model->fill_in_grad_obj(x, grad_f, new_x);
        return true;
    }
    
    bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g){
        return eval_g_(n, x, new_x, m, g);
    }
    
    /** Method to return the constraint residuals */
    template<typename T=type,typename std::enable_if<is_arithmetic<T>::value>::type* = nullptr>
    bool eval_g_(Index n, const Number* x, bool new_x, Index m, Number* g){
        
        assert(n==_model->get_nb_vars());
        _model->fill_in_cstr(x, g, new_x);
        return true;
    }
    
    bool eval_jac_g(Index n, const Number* x, bool new_x,
                    Index m, Index nele_jac, Index* iRow, Index *jCol,
                    Number* values){
        return eval_jac_g_(n, x, new_x, m, nele_jac, iRow, jCol, values);
    }
    
    /** Method to return:
     *   1) The structure of the jacobian (if "values" is NULL)
     *   2) The values of the jacobian (if "values" is not NULL)
     */
    template<typename T=type,typename std::enable_if<is_arithmetic<T>::value>::type* = nullptr>
    bool eval_jac_g_(Index n, const Number* x, bool new_x,
                                  Index m, Index nele_jac, Index* iRow, Index *jCol,
                                  Number* values){
        
        assert(n==_model->get_nb_vars());
        assert(m==_model->get_nb_cons());
        assert(nele_jac==_model->get_nb_nnz_g());
        if (values == NULL){
            _model->fill_in_jac_nnz(iRow, jCol);
        } else {
            _model->fill_in_jac(x, values, new_x);
        }
        
        return true;
    }
    
    bool eval_h(Index n, const Number* x, bool new_x,
                Number obj_factor, Index m, const Number* lambda,
                bool new_lambda, Index nele_hess, Index* iRow,
                Index* jCol, Number* values){
        return eval_h_(n, x, new_x, obj_factor, m, lambda, new_lambda, nele_hess, iRow, jCol, values);
    }
    /** Method to return:
     *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
     *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
     */
    template<typename T=type,typename std::enable_if<is_arithmetic<T>::value>::type* = nullptr>
    bool eval_h_(Index n, const Number* x, bool new_x,
                              Number obj_factor, Index m, const Number* lambda,
                              bool new_lambda, Index nele_hess, Index* iRow,
                              Index* jCol, Number* values){
        
        
        assert(n==_model->get_nb_vars());
        assert(m==_model->get_nb_cons());
        assert(nele_hess==_model->get_nb_nnz_h());
        if (values == NULL){
            _model->fill_in_hess_nnz(iRow, jCol);
        } else {
            _model->fill_in_hess(x, obj_factor, lambda, values, new_x);
        }
        
        return true;
    }
    
    bool get_variables_linearity(Index n, LinearityType* var_types){
        assert(n==_model->get_nb_vars());
        return false;
    }
    
    bool get_constraints_linearity(Index m, LinearityType* const_types){
        assert(m==_model->get_nb_cons());
        return false;
    }

};



#endif /* defined(__Gravity____IpoptProgram__) */
