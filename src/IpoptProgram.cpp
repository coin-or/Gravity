//
//  IpoptProgram.cpp
//  PowerTools++
//
//  Created by Hassan on 30/01/2015.
//  Copyright (c) 2015 NICTA. All rights reserved.
//

#include <Gravity/IpoptProgram.h>

using namespace std;

bool IpoptProgram::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                              Index& nnz_h_lag, Ipopt::TNLP::IndexStyleEnum& index_style){
    index_style = Ipopt::TNLP::C_STYLE;
    
//    _model->fill_in_maps();
    n = (Index)_model->get_nb_vars();
//    printf("n = %d;\n", n);
    m = (Index)_model->get_nb_cons();
    nnz_jac_g = (Index)_model->get_nb_nnz_g();
//    cout << "\n############## CALLING update_hess_link ##############\n";

    nnz_h_lag = (Index)_model->get_nb_nnz_h();
    _model->_nnz_h = nnz_h_lag;
    return true;
}

void IpoptProgram::finalize_solution(    Ipopt::SolverReturn               status    ,
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
    var_* v = NULL;
    var<int>* int_var = NULL;
    var<bool>* bin_var = NULL;
    var<float>* real_var = NULL;
    var<>* long_real_var = NULL;
    for (int i = 0; i<n; i++) {
//        v = _model->getVar(i);
//        switch (v->get_type()) {
//            case real:
//                real_var = (var<float>*)v;
//                real_var->set_val((float)x[i]);
//                break;
//            case longreal:
//                long_real_var = (var<>*)v;
//                long_real_var->set_val(x[i]);
//                break;
//            case integ:
//                int_var = (var<int>*)v;
//                int_var->set_val((int)x[i]);
//                break;
//            case binary:
//                bin_var = (var<bool>*)v;
//                bin_var->set_val((bool)x[i]);
//                break;
//            default:
//                break;
//        } ;
    }
//    if (_model->_store_duals) {
//        auto &cons = _model->get_cons();
//        for (size_t i=0; i<m; i++)
//            cons[i]->_dual = lambda[i];
//        for (int i = 0; i<n; i++) {
//            double dual = z_U[i] - z_L[i]; // upper bound active is +ve dual
//            v = _model->getVar(i);
//            switch (v->get_type()) {
//                case real:
//                    static_cast<var<float> *>(v)->_dual = dual;
//                    break;
//                case longreal:
//                    static_cast<var<> *>(v)->_dual = dual;
//                    break;
//                case integ:
//                    static_cast<var<int> *>(v)->_dual = dual;
//                    break;
//                case binary:
//                    static_cast<var<bool> *>(v)->_dual = dual;
//                    break;
//                default:
//                    break;
//            } ;
//        }
//    }
//    _model->check_feasible(x);
//    _model->_opt = _model->_obj->eval(x);
//    cout << "\n************** Objective Function Value = " << _model->_opt << " **************" << endl;
}

bool IpoptProgram::get_bounds_info(Index n, Number* x_l, Number* x_u,
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

bool IpoptProgram::get_starting_point(Index n, bool init_x, Number* x,
                                    bool init_z, Number* z_L, Number* z_U,
                                    Index m, bool init_lambda,
                                    Number* lambda){
    assert(n==_model->get_nb_vars());
    assert(m==_model->get_nb_cons());
    
    if (init_x) {
        _model->fill_in_var_init(x);
    }
    return true;
}

bool IpoptProgram::eval_f(Index n, const Number* x, bool new_x, Number& obj_value){
    
    assert(n==_model->get_nb_vars());
    _model->fill_in_obj(x, obj_value,new_x);
    return true;
}

bool IpoptProgram::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f){
    
    assert(n==_model->get_nb_vars());
    _model->fill_in_grad_obj(x, grad_f, new_x);
    return true;
}

bool IpoptProgram::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g){
    
    assert(n==_model->get_nb_vars());
//    if (!new_x)
        _model->fill_in_cstr(x, g, new_x);
    return true;
}

bool IpoptProgram::eval_jac_g(Index n, const Number* x, bool new_x,
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

std::vector<int> bounds2(int parts, int mem) {
    std::vector<int>bnd;
    int delta = mem / parts;
    int reminder = mem % parts;
    int N1 = 0, N2 = 0;
    bnd.push_back(N1);
    for (int i = 0; i < parts; ++i) {
        N2 = N1 + delta;
        if (i == parts - 1)
            N2 += reminder;
        bnd.push_back(N2);
        N1 = N2;
    }
    return bnd;
}

bool IpoptProgram::eval_h(Index n, const Number* x, bool new_x,
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

bool IpoptProgram::get_variables_linearity(Index n, LinearityType* var_types){
    assert(n==_model->get_nb_vars());
    _model->fill_in_var_linearity(var_types);
    return true;
//    return false;
}

bool IpoptProgram::get_constraints_linearity(Index m, LinearityType* const_types){
    assert(m==_model->get_nb_cons());
    _model->fill_in_cstr_linearity(const_types);
    return true;
//    return false;
}

