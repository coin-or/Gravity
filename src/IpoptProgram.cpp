//
//  IpoptProgram.cpp
//  PowerTools++
//
//  Created by Hassan on 30/01/2015.
//  Copyright (c) 2015 NICTA. All rights reserved.
//
#define DebugOn(x) cout << x
#include <gravity/IpoptProgram.h>

using namespace std;

bool IpoptProgram::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                              Index& nnz_h_lag, Ipopt::TNLP::IndexStyleEnum& index_style){
    index_style = Ipopt::TNLP::C_STYLE;
    n = (Index)_model->get_nb_vars();
    printf("n = %d;\n", n);
    m = (Index)_model->get_nb_cons();
    printf("m = %d;\n", m);
    nnz_jac_g = (Index)_model->get_nb_nnz_g();
    printf("number of non zeros in Jacobian = %d;\n", nnz_jac_g);
    nnz_h_lag = (Index)_model->get_nb_nnz_h();
    printf("number of non zeros in Hessian = %d;\n", nnz_h_lag);
    return true;
}

void IpoptProgram::finalize_solution(Ipopt::SolverReturn             status    ,
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
    //    _model->check_feasible(x);
    if(_model->_objt==maximize){
                _model->_obj *= -1;
    }
    _model->_obj_val = _model->_obj.eval();
//    for (size_t i=0; i<m; i++) {
//        _model->_cons[i]->_dual = lambda[i];
//    }
    cout << "\n************** Objective Function Value = " << _model->_obj_val << " **************" << endl;
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
    
//    if (init_x) {
        _model->fill_in_var_init(x);
//    DebugOn("initial point = \n");
//    DebugOn("x = [ ");
//    for (int i = 0; i<n; i++) {
//        DebugOn(to_string(x[i]) << " ");
//    }
//    DebugOn("]\n");
//    }
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
