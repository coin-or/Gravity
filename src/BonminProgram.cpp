//#include <Gravity/BonminProgram.h>
//
//using namespace std;
//
//BonminProgram::BonminProgram(Model* m){
//    model = m;
//}
//
//bool BonminProgram::get_variables_types(Index n, VariableType* var_types){
//    assert(n==model->get_nb_vars());
//    model->fill_in_var_types(var_types);
//    return true;
//}
//
//bool BonminProgram::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
//                                 Index& nnz_h_lag, Ipopt::TNLP::IndexStyleEnum& index_style){
//    index_style = Ipopt::TNLP::C_STYLE;
//    n = (Index)model->get_nb_vars();
//    m = (Index)model->get_nb_cons();
//    nnz_jac_g = (Index)model->get_nb_nnz_g();
//    model->update_hess_link();
//    nnz_h_lag = (Index)model->get_nb_nnz_h();
//    model->_nnz_h = nnz_h_lag;
//    return true;
//}
//
//void BonminProgram::finalize_solution(TMINLP::SolverReturn               status    ,
//                                      Index                             n         ,
//                                      const Number*                     x         ,
//                                      Number                            obj_value )
//{
//    if (x == nullptr) return; // Bonmin must have had an issue...
//
//    var_* v = NULL;
//    var<int>* int_var = NULL;
//    var<bool>* bin_var = NULL;
//    var<float>* real_var = NULL;
//    var<>* long_real_var = NULL;
//    for (int i = 0; i<n; i++) {
//        v = model->getVar(i);
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
//    }
////    model->check_feasible(x);
//    model->_opt = model->_obj->eval(x);
//    cout << "\n************** Objective Function Value = " << model->_opt << " **************" << endl;
//}
//
//bool BonminProgram::get_bounds_info(Index n, Number* x_l, Number* x_u,
//                                    Index m, Number* g_l, Number* g_u){
//    assert(n==model->get_nb_vars());
//    assert(m==model->get_nb_cons());
//    model->fill_in_var_bounds(x_l , x_u);
//    for (int i = 0; i<n; i++) {
////        if (x_l[i]==x_u[i]) {
////            printf("%f <= x[%d] <= %f\n",x_l[i], i, x_u[i]);
////        }
//    }
//    model->fill_in_cstr_bounds(g_l , g_u);
////    for (int i = 0; i<m; i++) {
////        printf("%f <= g[%d] <= %f\n",g_l[i], i, g_u[i]);
////    }
//
//    return true;
//}
//
//bool BonminProgram::get_starting_point(Index n, bool init_x, Number* x,
//                                       bool init_z, Number* z_L, Number* z_U,
//                                       Index m, bool init_lambda,
//                                       Number* lambda){
//    assert(n==model->get_nb_vars());
//    assert(m==model->get_nb_cons());
//
//    if (init_x) {
//        model->fill_in_var_init(x);
//    }
//    return true;
//}
//
//bool BonminProgram::eval_f(Index n, const Number* x, bool new_x, Number& obj_value){
//
//    assert(n==model->get_nb_vars());
//    model->fill_in_obj(x, obj_value);
//    return true;
//}
//
//bool BonminProgram::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f){
//
//    assert(n==model->get_nb_vars());
//    model->fill_in_grad_obj(x, grad_f);
//    return true;
//}
//
//bool BonminProgram::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g){
//
//    assert(n==model->get_nb_vars());
////    if (!new_x)
//    model->fill_in_cstr(x, g, new_x);
//    return true;
//}
//
//bool BonminProgram::eval_jac_g(Index n, const Number* x, bool new_x,
//                               Index m, Index nele_jac, Index* iRow, Index *jCol,
//                               Number* values){
//
//    assert(n==model->get_nb_vars());
//    assert(m==model->get_nb_cons());
//    assert(nele_jac==model->get_nb_nnz_g());
//    if (values == NULL){
//        model->fill_in_jac_nnz(iRow, jCol);
//    } else {
////        if (!new_x) {
//        model->fill_in_jac(x, values, new_x);
////        }
//    }
//
//    return true;
//}
//
//bool BonminProgram::eval_h(Index n, const Number* x, bool new_x,
//                           Number obj_factor, Index m, const Number* lambda,
//                           bool new_lambda, Index nele_hess, Index* iRow,
//                           Index* jCol, Number* values){
//
//
//    assert(n==model->get_nb_vars());
//    assert(m==model->get_nb_cons());
//    assert(nele_hess==model->get_nb_nnz_h());
//    if (values == NULL){
//        model->fill_in_hess_nnz(iRow, jCol);
//    } else {
////        if(!new_x)
//        model->fill_in_hess(x, obj_factor, lambda, values, new_x);
////        int nr_threads = 6;
////        std::vector<std::thread> threads;
////            //Split constraints into nr_threads parts
////        std::vector<int> limits = bounds2(nr_threads, (int)n);
////            //Launch nr_threads threads:
////        for (int i = 0; i < nr_threads; ++i) {
////            threads.push_back(std::thread(&Model::fill_in_hess_multithread, model, x, obj_factor, lambda, values, limits[i], limits[i+1]));
////        }
////            //Join the threads with the main thread
////        for(auto &t : threads){
////            t.join();
////        }
//    }
//    return true;
//}
//
//bool BonminProgram::get_variables_linearity(Index n, TNLP::LinearityType* var_types){
//    assert(n==model->get_nb_vars());
//    model->fill_in_var_linearity(var_types);
//    return true;
//}
//
//bool BonminProgram::get_constraints_linearity(Index m, TNLP::LinearityType* const_types){
//    assert(m==model->get_nb_cons());
//    model->fill_in_cstr_linearity(const_types);
//    return true;
//}
