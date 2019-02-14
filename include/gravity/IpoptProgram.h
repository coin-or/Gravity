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
    shared_ptr<Model<type>> _model = nullptr;
    
    IpoptProgram(shared_ptr<Model<type>> m):_model(m){
    }
    
    
    
    
    void update_model(){
//        _model->reset_funcs();
        _model->fill_in_maps();
    }
    
    
    
    
    /** Method to return some info about the nlp */
    bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                                    Index& nnz_h_lag, Ipopt::TNLP::IndexStyleEnum& index_style){
        index_style = Ipopt::TNLP::C_STYLE;
        n = (Index)_model->get_nb_vars();
        m = (Index)_model->get_nb_cons();
        nnz_jac_g = (Index)_model->get_nb_nnz_g();
        _model->_jac_vals.resize(nnz_jac_g,0);
        nnz_h_lag = (Index)_model->get_nb_nnz_h();
        _model->_first_call_jac = true;
        //If quadratic model and we're resolving no need to reset these
        _model->_first_call_hess = true;
        _model->_first_call_gard_obj = true;
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
//        _model->_obj_val = obj_value;
        //    _model->print_solution();
        //    _model->compute_funcs();
        //    _model->check_feasible(x);
        if(_model->_objt==maximize){
            _model->_obj->reverse_sign();
        }
        //    _model->_obj_val = _model->_obj.eval();
        for (auto &cp: _model->_cons) {
            cp.second->_dual.resize(cp.second->_dim[0]);
            auto idx = 0;
            for (size_t inst = 0; inst < cp.second->_dim[0]; inst++) {
                if (!*cp.second->_all_lazy || !cp.second->_lazy[inst]) {
                    cp.second->_dual[inst] = lambda[cp.second->_id + idx++];
                }
            }
        }
//        for (auto &vp: _model->_vars) {
//            auto nb_inst = vp.second->get_dim();
//            vp.second->_u_dual.resize(nb_inst);
//            vp.second->_l_dual.resize(nb_inst);
//            auto vid = vp.second->get_id();
//            for (unsigned inst = 0; inst < nb_inst; inst++) {
//                vp.second->_u_dual[inst] = z_U[vid + vp.second->get_id_inst(inst)];
//                vp.second->_l_dual[inst] = z_L[vid + vp.second->get_id_inst(inst)];
//            }
//        }
//        cout << "\n************** Objective Function Value = " << _model->eval<double>(_model->_obj,0) << " **************" << endl;
    }
    
    /** Method to return the bounds for my problem */
    bool get_bounds_info(Index n, Number* x_l, Number* x_u,
                                       Index m, Number* g_l, Number* g_u){
        //    printf("n = %d;\n", n);
        assert(n==_model->get_nb_vars());
        assert(m==_model->get_nb_cons());
        _model->fill_in_var_bounds(x_l , x_u);
        //    for (int i = 0; i<n; i++) {
        //        if (x_l[i]==x_u[i]) {
        //            printf("%f <= x[%d] <= %f\n",x_l[i], i, x_u[i]);
        //        }
        //    }
        _model->fill_in_cstr_bounds(g_l , g_u);
        //    for (int i = 0; i<m; i++) {
        //        printf("%f <= g[%d] <= %f\n",g_l[i], i, g_u[i]);
        //    }
        
        return true;
    }
    
    /** Method to return the starting point for the algorithm */
    bool get_starting_point(Index n, bool init_x, Number* x,
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
        //    DebugOn("initial point = \n");
        //    DebugOn("x = [ ");
        //    for (int i = 0; i<n; i++) {
        //        DebugOn(to_string(x[i]) << " ");
        //    }
        //    DebugOn("]\n");
        
        return true;
    }
    
    /** Method to return the objective value */
    bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value){
        
        assert(n==_model->get_nb_vars());
        _model->fill_in_obj(x, obj_value,new_x);
        return true;
    }
    
    /** Method to return the gradient of the objective */
    bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f){
        
        assert(n==_model->get_nb_vars());
        _model->fill_in_grad_obj(x, grad_f, new_x);
        return true;
    }
    
    /** Method to return the constraint residuals */
    bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g){
        
        assert(n==_model->get_nb_vars());
        //    if (!new_x)
        _model->fill_in_cstr(x, g, new_x);
        return true;
    }
    
    /** Method to return:
     *   1) The structure of the jacobian (if "values" is NULL)
     *   2) The values of the jacobian (if "values" is not NULL)
     */
    bool eval_jac_g(Index n, const Number* x, bool new_x,
                                  Index m, Index nele_jac, Index* iRow, Index *jCol,
                                  Number* values){
        
        assert(n==_model->get_nb_vars());
        assert(m==_model->get_nb_cons());
        assert(nele_jac==_model->get_nb_nnz_g());
        if (values == NULL){
            _model->fill_in_jac_nnz(iRow, jCol);
        } else {
            //        if (!new_x) {
            _model->fill_in_jac(x, values, new_x);
            //        }
        }
        
        return true;
    }
    
    /** Method to return:
     *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
     *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
     */
    bool eval_h(Index n, const Number* x, bool new_x,
                              Number obj_factor, Index m, const Number* lambda,
                              bool new_lambda, Index nele_hess, Index* iRow,
                              Index* jCol, Number* values){
        
        
        assert(n==_model->get_nb_vars());
        assert(m==_model->get_nb_cons());
        assert(nele_hess==_model->get_nb_nnz_h());
        if (values == NULL){
            _model->fill_in_hess_nnz(iRow, jCol);
        } else {
            //        if(!new_x)
            _model->fill_in_hess(x, obj_factor, lambda, values, new_x);
            //        int nr_threads = 6;
            //        std::vector<std::thread> threads;
            //            //Split constraints into nr_threads parts
            //        std::vector<int> limits = bounds2(nr_threads, (int)n);
            //            //Launch nr_threads threads:
            //        for (int i = 0; i < nr_threads; ++i) {
            //            threads.push_back(std::thread(&_model::fill_in_hess_multithread, _model, x, obj_factor, lambda, values, limits[i], limits[i+1]));
            //        }
            //            //Join the threads with the main thread
            //        for(auto &t : threads){
            //            t.join();
            //        }
        }
        
        return true;
    }
    
    bool get_variables_linearity(Index n, LinearityType* var_types){
        assert(n==_model->get_nb_vars());
        //    _model->fill_in_var_linearity(var_types);
        //    return true;
        return false;
    }
    
    bool get_constraints_linearity(Index m, LinearityType* const_types){
        assert(m==_model->get_nb_cons());
        //    _model->fill_in_cstr_linearity(const_types);
        //    return true;
        return false;
    }

};



#endif /* defined(__Gravity____IpoptProgram__) */
