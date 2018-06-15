//
//  MosekProgram.cpp
//  Gravity
//
//  Created by Guanglei Wang on 14/7/17.
//
//
#include <Eigen/Dense>
#include "MosekProgram.h"
using namespace mosek;
using namespace monty;

MosekProgram::MosekProgram() {    
    _mosek_model = new fusion::Model("noname");
    _mosek_model->setLogHandler([=](const std::string & msg){std::cout << msg << std::flush;});

}

MosekProgram::MosekProgram(Model* m) {    
    _mosek_model = new fusion::Model("noname");
    _mosek_model->setLogHandler([=](const std::string & msg){std::cout << msg << std::flush;});
    _model = m;
}

MosekProgram::~MosekProgram() {
    _mosek_model->dispose();
}

void MosekProgram::update_model() {};
    

// remain to do
bool MosekProgram::solve(bool relax) {
    _mosek_model->setSolverParam("log", 10);
    // Set max solution time
//    _mosek_model->setSolverParam("mioMaxTime", 360.0);
    // Set max relative gap (to its default value)
//    _mosek_model->setSolverParam("intpntCoTolRelGap", 1e-10);
//    _mosek_model->setSolverParam("intpntCoTolPfeas", 1);
//    _mosek_model->setSolverParam("intpntQoTolPfeas", 1);
//    _mosek_model->setSolverParam("intpntCoTolMuRed", 1);
    // Set max absolute gap (to its default value)
    _mosek_model->setSolverParam("mioTolRelGap", 1e-4);
//    _mosek_model->setSolverParam("mioTolAbsGap", 0.0);
//    _mosek_model->setSolverParam("numThreads",1);
    if(!_output) _mosek_model->setSolverParam("log", 0);
    if(relax) {
        for (auto &vv: _mosek_vars)
            vv->makeContinuous();
    }
    _mosek_model->solve();

    if (!relax) {
        cout <<_mosek_model->getProblemStatus(fusion::SolutionType::Integer) << endl;
        cout << _mosek_model->getPrimalSolutionStatus() << endl;
    }
    else {
        cout <<_mosek_model->getProblemStatus(fusion::SolutionType::Interior) << endl;
        cout << _mosek_model->getPrimalSolutionStatus() << endl;
    }

    DebugOn("Cost = " << _mosek_model->primalObjValue() << std::endl);
    
    // set the optimal value.
    _model->_obj_val = _mosek_model->primalObjValue();


    // Note that there is only one way to retrieve solutions, i.e.,
    // variable.level(). It only returns double. Thus, we need to cast the
    // solution to required types.
    //std::cout << "Solution = " << std::endl;
    // cout << "dim: " << _mosek_vars.size() << endl;
    //auto sol = _mosek_vars[0]->level();
    //mosek::fusion::Variable::t s= _mosek_vars[0];
    //cout << "size of s: " << s->size() << endl;

    for (auto i = 0; i < _mosek_vars.size(); i++) {
        auto sol = _mosek_vars[i]->level();
        if(i >= _model->_vars.size()) continue;
        if (_model->_vars[i]->get_intype()== binary_ ||_model->_vars[i]->get_intype()== integer_) {
            for (auto j = 0; j < _model->_vars[i]->get_nb_instances(); j++) {
                auto val = (*sol)[j];
                val = round(val);
                poly_set_val(j, val , _model->_vars[i]);
            }
        }
        else if(!_model->_vars[i]->_is_matrix){
            for (auto j = 0; j < _model->_vars[i]->get_nb_instances(); j++) {
                poly_set_val(j, (*sol)[j], _model->_vars[i]);
            }
        }
        else {
            int n = _model->_vars[i]->_dim[0];
            for(auto j1 = 0; j1 < n; j1++) {
                for(auto j2 = j1; j2 < n; j2++) {
                    string key = to_string(j1)+","+to_string(j2);
                    ((param<double>*)_model->_vars[i])->set_val(key, (*sol)[j1*n+j2]);
                }
            }


        }
    }
    return 0;
}

void MosekProgram::fill_in_mosek_vars() {
    param_* v;
    for(auto& v_p: _model->_vars)
    {
        v = v_p.second;
        if(v->get_id()==-1) {
            throw invalid_argument("Variable needs to be added to model first: use add_var(v) function:" + v->get_name());
        }
        switch (v->get_intype()) {
        case float_: {
            if (v->get_type() == var_c) {
                auto real_var = (var<float>*)v;
                //for (int i = 0; i < real_var->get_nb_instances(); i++) {
                //    _mosek_vars.push_back(_mosek_model->variable(real_var->get_name()+"_"+to_string(i),fusion::Domain::inRange(real_var->get_lb(i), real_var->get_ub(i))));
                //}
                auto lb  = new_array_ptr<double,1>(real_var->get_nb_instances());
                auto ub  = new_array_ptr<double,1>(real_var->get_nb_instances());
                for (int i = 0; i < real_var->get_nb_instances(); i++) {
                    (*lb)[i] = real_var->get_lb(i);
                    (*ub)[i] = real_var->get_ub(i);
                }
                if(!real_var->_in_q_cone && !real_var->_psd) _mosek_vars.push_back(_mosek_model->variable(real_var->get_name(), real_var->get_nb_instances(), fusion::Domain::inRange(lb,ub)));
                else if(real_var->_in_q_cone) _mosek_vars.push_back(_mosek_model->variable(real_var->get_name(), fusion::Domain::inQCone(real_var->get_nb_instances())));
                else {
//                    _mosek_vars.push_back(_mosek_model->variable(real_var->get_name(), fusion::Domain::inPSDCone((1 + sqrt(1+8*real_var->get_nb_instances()))/2)));
                    _mosek_vars.push_back(_mosek_model->variable(real_var->get_name(), fusion::Domain::inPSDCone(real_var->_dim[0])));
                }
            }
            else {
                auto sdp_var = (sdpvar<float>*)v;
                _mosek_vars.push_back(_mosek_model->variable(sdp_var->get_name(), fusion::Domain::inPSDCone(sdp_var->_symdim)));
            }
            break;
        }
        case long_: {
            if (v->get_type() == var_c) {
                auto real_var = (var<long double>*)v;
                //for (int i = 0; i < real_var->get_nb_instances(); i++) {
                //    _mosek_vars.push_back(_mosek_model->variable(real_var->get_name()+"_"+to_string(i),fusion::Domain::inRange(real_var->get_lb(i), real_var->get_ub(i))));
                //}
                auto lb  = new_array_ptr<double,1>(real_var->get_nb_instances());
                auto ub  = new_array_ptr<double,1>(real_var->get_nb_instances());
                for (int i = 0; i < real_var->get_nb_instances(); i++) {
                    (*lb)[i] = real_var->get_lb(i);
                    (*ub)[i] = real_var->get_ub(i);
                }
                if(!real_var->_in_q_cone)_mosek_vars.push_back(_mosek_model->variable(real_var->get_name(), real_var->get_nb_instances(), fusion::Domain::inRange(lb,ub)));
                else _mosek_vars.push_back(_mosek_model->variable(real_var->get_name(), fusion::Domain::inQCone(real_var->get_nb_instances())));
            }
            else {
                auto sdp_var = (sdpvar<long double>*)v;
                _mosek_vars.push_back(_mosek_model->variable(sdp_var->get_name(), fusion::Domain::inPSDCone(sdp_var->get_nb_instances())));
            }
            break;
        }
        case double_: {
            if (v->get_type() == var_c) {
                auto real_var = (var<double>*)v;
                //for (int i = 0; i < real_var->get_nb_instances(); i++) {
                //    _mosek_vars.push_back(_mosek_model->variable(real_var->get_name()+"_"+to_string(i),fusion::Domain::inRange(real_var->get_lb(i), real_var->get_ub(i))));
                //}
                auto lb  = new_array_ptr<double,1>(real_var->get_nb_instances());
                auto ub  = new_array_ptr<double,1>(real_var->get_nb_instances());
                for (int i = 0; i < real_var->get_nb_instances(); i++) {
                    (*lb)[i] = real_var->get_lb(i);
                    (*ub)[i] = real_var->get_ub(i);
                }
                if(!real_var->_in_q_cone && !real_var->_psd) _mosek_vars.push_back(_mosek_model->variable(real_var->get_name(), real_var->get_nb_instances(), fusion::Domain::inRange(lb,ub)));
                else if(real_var->_in_q_cone) _mosek_vars.push_back(_mosek_model->variable(real_var->get_name(), fusion::Domain::inQCone(real_var->get_nb_instances())));
                else {
//                    _mosek_vars.push_back(_mosek_model->variable(real_var->get_name(), fusion::Domain::inPSDCone((1 + sqrt(1+8*real_var->get_nb_instances()))/2.)));
                    _mosek_vars.push_back(_mosek_model->variable(real_var->get_name(), fusion::Domain::inPSDCone(real_var->_dim[0])));
                }
            }
            else {
                auto sdp_var = (sdpvar<double>*)v;
                size_t num = sdp_var->_symdim;
                auto c  = _mosek_model->variable(sdp_var->get_name(), fusion::Domain::inPSDCone(num));
                _mosek_vars.push_back(c);
                std::cout << c->toString() << endl;

            }
            break;
        }
        case integer_: {
            if (v->get_type() == var_c) {
                auto real_var = (var<int>*)v;
                auto lb  = new_array_ptr<double,1>(real_var->get_nb_instances());
                auto ub  = new_array_ptr<double,1>(real_var->get_nb_instances());
                for (int i = 0; i < real_var->get_nb_instances(); i++) {
                    (*lb)[i] = real_var->get_lb(i);
                    (*ub)[i] = real_var->get_ub(i);
                }
                _mosek_vars.push_back(_mosek_model->variable(real_var->get_name(), real_var->get_nb_instances(), fusion::Domain::integral(fusion::Domain::inRange(lb,ub))));
            }
            else {
                auto sdp_var = (sdpvar<int>*)v;
                _mosek_vars.push_back(_mosek_model->variable(sdp_var->get_name(), fusion::Domain::inPSDCone(sdp_var->get_nb_instances())));
            }
            break;
        }
        case short_: {
            if (v->get_type() == var_c) {
                auto real_var = (var<short>*)v;
                auto lb  = new_array_ptr<double,1>(real_var->get_nb_instances());
                auto ub  = new_array_ptr<double,1>(real_var->get_nb_instances());
                for (int i = 0; i < real_var->get_nb_instances(); i++) {
                    (*lb)[i] = real_var->get_lb(i);
                    (*ub)[i] = real_var->get_ub(i);
                }
                _mosek_vars.push_back(_mosek_model->variable(real_var->get_name(), real_var->get_nb_instances(), fusion::Domain::integral(fusion::Domain::inRange(lb,ub))));
            }
            else {
                auto sdp_var = (sdpvar<short>*)v;
                _mosek_vars.push_back(_mosek_model->variable(sdp_var->get_name(), fusion::Domain::inPSDCone(sdp_var->get_nb_instances())));
            }
            break;
        }
        case binary_: {
            if (v->get_type() == var_c) {
                auto real_var = (var<bool>*)v;
                //for (unsigned int i = 0; i < real_var->get_nb_instances(); i++) {
                //    _mosek_vars.push_back(_mosek_model->variable(real_var->get_name()+"_"+to_string(i),fusion::Domain::binary()));
                //}
                _mosek_vars.push_back(_mosek_model->variable(real_var->get_name(), real_var->get_nb_instances(), fusion::Domain::binary()));
            }
            else {
                auto sdp_var = (sdpvar<bool>*)v;
                _mosek_vars.push_back(_mosek_model->variable(sdp_var->get_name(), fusion::Domain::inPSDCone(sdp_var->get_nb_instances())));
            }
            break;
        }
        default:
            break;
        }
    }
}

/*void MosekProgram::create_mosek_constraints() {
    //size_t idx = 0;
    //size_t nb_inst = 0;
    size_t idx = 0, idx_inst = 0,idx1 = 0, idx2 = 0, idx_inst1 = 0, idx_inst2 = 0, nb_inst = 0, inst = 0;
    size_t c_idx_inst = 0;
    shared_ptr<Constraint> c;
    for(auto& p: _model->_cons) {
        c = p.second;
//        c->print();
        if (c->is_nonlinear()) {
            cout <<  "We haven't implemented quadratic expressions interface for mosek" << endl;
            throw invalid_argument("Mosek cannot handle nonlinear constraints that are not convex quadratic.\n");
        }
        nb_inst = c->get_nb_instances();
        auto fusion_sums = new ndarray<fusion::Expression::t,1>(shape(nb_inst));
        for (int i = 0; i< nb_inst; i++) {
            // get constant part first.
            monty::rc_ptr<mosek::fusion::Expression > expr= fusion::Expr::constTerm(poly_eval(c->get_cst())); // expr is a pointer to the Expression.
            for (auto& it1: c->get_lterms()) {
                //idx = it1.second._p->get_id();
                idx = it1.second._p->get_vec_id();
                //cout << "get_id: " << it1.second._p->get_id() << endl;
                //cout << "get_vec_id: " << it1.second._p->get_vec_id() << endl;

                CType vartype = it1.second._p->get_type();
                if (vartype == var_c) {
                    if(!it1.second._p->_is_matrix) {

                        if (it1.second._coef->_is_transposed) {
                            auto coefs = new_array_ptr<double, 1>(it1.second._p->get_nb_instances());
                            for (int j = 0; j < it1.second._p->get_nb_instances(); j++) {
                                (*coefs)(j) = poly_eval(it1.second._coef, j);
                            }
                            expr = fusion::Expr::add(expr, fusion::Expr::dot(coefs, _mosek_vars[idx]));
                            //    cout << "expr" << expr->toString() << endl;
                        } else {
                            if (is_indexed(it1.second._p)) {
                                idx_inst = it1.second._p->get_id_inst(i);
                            } else {
                                idx_inst = inst;
                            }
                            if (is_indexed(it1.second._coef)) {
                                c_idx_inst = get_poly_id_inst(it1.second._coef);
                            } else {
                                c_idx_inst = inst;
                            }
                            auto lterm = fusion::Expr::mul(poly_eval(it1.second._coef, c_idx_inst),
                                                           _mosek_vars[idx]->index(idx_inst));

                            if (!it1.second._sign) {
                                lterm = fusion::Expr::mul(-1, lterm);
                            }
                            expr = fusion::Expr::add(expr, lterm);
                        }
                    }else{
                        if (it1.second._coef->_is_transposed) {
                            auto coefs = new_array_ptr<double, 1>(it1.second._p->get_nb_instances());
                            for (int j = 0; j < it1.second._p->get_nb_instances(); j++) {
                                (*coefs)(j) = poly_eval(it1.second._coef, j);
                            }
                            expr = fusion::Expr::add(expr, fusion::Expr::dot(coefs, _mosek_vars[idx]));
                            //    cout << "expr" << expr->toString() << endl;
                        } else {
                            pair<size_t, size_t>  pair;
                            if (is_indexed(it1.second._p)) {
//                                idx_inst = it1.second._p->get_id_inst(i);
                                pair = it1.second._p->get_sdp_inst(i);
                            } else {
                                idx_inst = inst;
                            }
                            if (is_indexed(it1.second._coef)) {
                                c_idx_inst = get_poly_id_inst(it1.second._coef);
                            } else {
                                c_idx_inst = inst;
                            }
                            DebugOff("\nindex: " << pair.first << ", " << pair.second);
                            auto lterm = fusion::Expr::mul(poly_eval(it1.second._coef, c_idx_inst),
                                                           _mosek_vars[idx]->index(pair.first,pair.second));

                            if (!it1.second._sign) {
                                lterm = fusion::Expr::mul(-1, lterm);
                            }
                            expr = fusion::Expr::add(expr, lterm);
                        }
                    }


                }
                else if(vartype == sdpvar_c) {
                    if (it1.second._coef->_is_transposed) {
                        auto coefs = new_array_ptr<double,1>(it1.second._p->get_nb_instances());
                        for (int j = 0; j<it1.second._p->get_nb_instances(); j++) {
                            (*coefs)(j) = poly_eval(it1.second._coef,j);
                        }
                        expr = fusion::Expr::add(expr,fusion::Expr::dot(coefs,_mosek_vars[idx]));
                    }
                    else {
                        pair<size_t, size_t>  pair = make_pair(0, 0);

                        if (is_indexed(it1.second._p)) {
                            pair = it1.second._p->get_sdpid();
                        }

                        if (is_indexed(it1.second._coef)) {
                            c_idx_inst = get_poly_id_inst(it1.second._coef);
                        }
                        else {
                            c_idx_inst = inst;
                        }
                        auto lterm = fusion::Expr::mul(poly_eval(it1.second._coef, c_idx_inst), _mosek_vars[idx]->index(pair.first, pair.second));
                        if (!it1.second._sign) {
                            lterm = fusion::Expr::mul(-1, lterm);
                        }
                        expr = fusion::Expr::add(expr,lterm);
                    }

                }else{ //param
                    if (it1.second._coef->_is_transposed) {
                        auto coefs = new_array_ptr<double,1>(it1.second._p->get_nb_instances());
                        for (int j = 0; j<it1.second._p->get_nb_instances(); j++) {
                            (*coefs)(j) = poly_eval(it1.second._coef,j);
                        }
                        expr = fusion::Expr::add(expr,fusion::Expr::dot(coefs,_mosek_vars[idx]));
                        //    cout << "expr" << expr->toString() << endl;
                    }
                    else { //this results in always taking the first instance? it doesn't use i at all?
                        if (is_indexed(it1.second._p)) {
                            idx_inst = it1.second._p->get_id_inst(i); //fixed here, look for similar errors
                        }
                        else {
                            idx_inst = inst;
                        }
                        if (is_indexed(it1.second._coef)) {
                            c_idx_inst = get_poly_id_inst(it1.second._coef);
                        }
                        else {
                            c_idx_inst = inst;
                        }
                        auto lterm = fusion::Expr::constTerm(poly_eval(it1.second._coef, c_idx_inst) * poly_eval(it1.second._p, i)); // was idx_inst
                        if (!it1.second._sign) {
                            lterm = fusion::Expr::mul(-1, lterm);
                        }
                        expr = fusion::Expr::add(expr,lterm);
                    }
                }
            }

            DebugOff("\nconstr in Mosek: " << expr->toString());
            if(c->get_type()==geq) {
                DebugOff("\n >= " << c->get_rhs());
                _mosek_model->constraint(c->get_name()+to_string(i), expr, fusion::Domain::greaterThan(c->get_rhs()));
            }
            else if(c->get_type()==leq) {
                DebugOff("\n <= " << c->get_rhs());
                _mosek_model->constraint(c->get_name()+to_string(i), expr, fusion::Domain::lessThan(c->get_rhs()));
            }
            else if(c->get_type()==eq) {
                DebugOff("\n" << expr->toString());
                DebugOff(" = " << c->get_rhs());
                _mosek_model->constraint(c->get_name()+to_string(i), expr, fusion::Domain::equalsTo(c->get_rhs()));
            }
            inst++;
        }
    }
}*/

size_t get_num_qterms(map<string, qterm>& qterms) {
    size_t qn = 0;
    for (auto& it1: qterms) {
        if ((it1.second._coef->_is_transposed || it1.second._coef->_is_matrix) && !it1.second._p->first->_is_matrix) {
            size_t dim = it1.second._p->first->get_nb_instances();
            qn += dim;
        }
        else if(it1.second._p->first->_is_transposed){
            size_t dim = it1.second._p->first->get_nb_instances();
            qn += dim;
        }
        else qn++;
    }
    return qn;
}

fusion::Expression::t MosekProgram::create_lin_expr(map<string, lterm>& lt, constant_* cst, size_t inst){
    size_t idx = 0, idx_inst = 0, c_idx_inst = 0;
    auto fusion_lterms = new_array_ptr<fusion::Expression::t,1>(lt.size()+1);
    (*fusion_lterms)[0] = fusion::Expr::constTerm(poly_eval(cst,inst));
    int lterm_idx = 1;

    for (auto& it1: lt) {
        idx = it1.second._p->get_vec_id();

        CType vartype = it1.second._p->get_type();
        if (vartype == var_c) {
//                    cout << "\nvar = " << it1.second._p->get_name();
//                    cout << "\nfull Mosek var: " << _mosek_vars[idx]->toString();
            if(!it1.second._p->_is_matrix) {
                if (it1.second._coef->_is_transposed) {
                    auto coefs = new_array_ptr<double, 1>(it1.second._p->get_nb_instances(inst));
                    auto vars = new_array_ptr<fusion::Variable::t,1>(it1.second._p->get_nb_instances(inst));
                    for (int j = 0; j < it1.second._p->get_nb_instances(inst); j++) {
                        (*coefs)(j) = poly_eval(it1.second._coef,inst,j);
                        idx_inst = it1.second._p->get_id_inst(inst,j);
//                                if(idx_inst >= _mosek_vars[idx]->size()){
//                                    cout << "\nidx_inst >= var size: idx_inst = " << idx_inst << ", var size = " << _mosek_vars[idx]->size();
//                                }
                        (*vars)(j) = _mosek_vars[idx]->index(idx_inst);
//                                cout << "\nj, idx, idx_inst = " << j << ", " << idx << ", " << idx_inst << ", var = " <<  _mosek_vars[idx]->index(idx_inst)->toString();
                    }
                    rc_ptr<fusion::Variable> P;
                    if(it1.second._p->get_nb_instances(inst)!=1) P = fusion::Var::vstack(vars);
                    else P = _mosek_vars[idx]->index(idx_inst);
//                            cout << "\nCoef = ";
//                            for(auto& coef: *coefs)
//                                cout << coef << "; ";
//                            cout << "\nP = " << P->toString();
                    auto lterm = fusion::Expr::dot(coefs, P);
                    if (!it1.second._sign) lterm = fusion::Expr::mul(-1, lterm);
                    (*fusion_lterms)[lterm_idx] = lterm;
                    //expr = fusion::Expr::add(expr, fusion::Expr::dot(coefs, _mosek_vars[idx]));
                } else {
                    if (is_indexed(it1.second._p)) idx_inst = it1.second._p->get_id_inst(inst);
                    else idx_inst = inst;

                    if (is_indexed(it1.second._coef)) c_idx_inst = get_poly_id_inst(it1.second._coef);
                    else c_idx_inst = inst;

                    auto lterm = fusion::Expr::mul(poly_eval(it1.second._coef, inst), //was c_idx_inst instead of i
                                                   _mosek_vars[idx]->index(idx_inst));

                    if (!it1.second._sign) lterm = fusion::Expr::mul(-1, lterm);
                    (*fusion_lterms)[lterm_idx] = lterm;
                    //expr = fusion::Expr::add(expr, lterm);
                }
            }else{
                if (it1.second._coef->_is_transposed) {
                    auto coefs = new_array_ptr<double, 1>(it1.second._p->get_nb_instances());
                    for (int j = 0; j < it1.second._p->get_nb_instances(); j++) {
                        (*coefs)(j) = poly_eval(it1.second._coef, j);
                    }
                    (*fusion_lterms)[lterm_idx] = fusion::Expr::dot(coefs, _mosek_vars[idx]);
                    //expr = fusion::Expr::add(expr, fusion::Expr::dot(coefs, _mosek_vars[idx]));
                } else {
                    pair<size_t, size_t>  pair;
                    if (is_indexed(it1.second._p)) {
//                                idx_inst = it1.second._p->get_id_inst(i);
                        pair = it1.second._p->get_sdp_inst(inst);
                    } else idx_inst = inst;

                    if (is_indexed(it1.second._coef)) c_idx_inst = get_poly_id_inst(it1.second._coef);
                    else c_idx_inst = inst;

                    DebugOff("\nindex: " << pair.first << ", " << pair.second);
                    DebugOff("\nMosek var:" << _mosek_vars[idx]->toString());

                    auto lterm = fusion::Expr::mul(poly_eval(it1.second._coef, c_idx_inst),
                                                   _mosek_vars[idx]->index(pair.first,pair.second));

                    if (!it1.second._sign) lterm = fusion::Expr::mul(-1, lterm);
                    (*fusion_lterms)[lterm_idx] = lterm;
                    //expr = fusion::Expr::add(expr, lterm);
                }
            }

        }
        else{ //param
            if (it1.second._coef->_is_transposed) {
                auto coefs = new_array_ptr<double,1>(it1.second._p->get_nb_instances());
                for (int j = 0; j<it1.second._p->get_nb_instances(); j++) {
                    (*coefs)(j) = poly_eval(it1.second._coef,j);
                }
                (*fusion_lterms)[lterm_idx] = fusion::Expr::dot(coefs,_mosek_vars[idx]);
                //expr = fusion::Expr::add(expr,fusion::Expr::dot(coefs,_mosek_vars[idx]));
            }
            else { //this results in always taking the first instance? it doesn't use i at all?
                if (is_indexed(it1.second._p)) idx_inst = it1.second._p->get_id_inst(inst); //fixed here, look for similar errors
                else idx_inst = inst;

                if (is_indexed(it1.second._coef)) c_idx_inst = get_poly_id_inst(it1.second._coef);
                else c_idx_inst = inst;

                auto lterm = fusion::Expr::constTerm(poly_eval(it1.second._coef, c_idx_inst) * poly_eval(it1.second._p, inst)); // was idx_inst
                if (!it1.second._sign) lterm = fusion::Expr::mul(-1, lterm);

                (*fusion_lterms)[lterm_idx] = lterm;
//                        expr = fusion::Expr::add(expr,lterm);
            }
        }
        lterm_idx++;
    }
    fusion::Expression::t res;
    res = fusion::Expr::add(std::shared_ptr<ndarray<fusion::Expression::t,1>>(fusion_lterms));
    return res;
}

//TODO: check use inst
fusion::Expression::t MosekProgram::form_Fx(map<string, qterm>& qterms, size_t qn, size_t inst){
    size_t aidx = 0, idx1 = 0, idx2 = 0, idx_inst = 0;

    //create a var map
    map<pair<size_t,size_t>,int> qvars;
    int idx_in_A = 0;
    for(auto& it1: qterms) {
        idx1 = it1.second._p->first->get_vec_id();
        idx2 = it1.second._p->second->get_vec_id();
        if ((it1.second._coef->_is_transposed || it1.second._coef->_is_matrix) && !it1.second._p->first->_is_matrix) {
            auto dim = it1.second._p->first->get_nb_instances(inst);
            for (int j = 0; j < dim; j++) {
                idx_inst = it1.second._p->first->get_id_inst(inst, j);
                pair<size_t, size_t> p = make_pair(idx1, idx_inst);
                if (qvars.find(p) == qvars.end()) {
                    qvars.insert(make_pair(p, idx_in_A));
                    idx_in_A++;
                }
                idx_inst = it1.second._p->second->get_id_inst(inst, j);
                p = make_pair(idx2, idx_inst);
                if (qvars.find(p) == qvars.end()) {
                    qvars.insert(make_pair(p, idx_in_A));
                    idx_in_A++;
                }
            }
        } else if (it1.second._p->first->_is_transposed) {
//                auto dim = _p->first->get_nb_instances(i);
//                for (int j = 0; j<dim; j++)
//                    res += poly_eval(_coef,i,j) * poly_eval(_p->first, i,j)* poly_eval(_p->second, i,j);
        } else {
            idx_inst = it1.second._p->first->get_id_inst(inst);
            pair<size_t, size_t> p = make_pair(idx1, idx_inst);
            if (qvars.find(p) == qvars.end()) {
                qvars.insert(make_pair(p, idx_in_A));
                idx_in_A++;
            }
            idx_inst = it1.second._p->second->get_id_inst(inst);
            p = make_pair(idx2, idx_inst);
            if (qvars.find(p) == qvars.end()) {
                qvars.insert(make_pair(p, idx_in_A));
                idx_in_A++;
            }
        }
    }
    auto qvars_arr = new_array_ptr<fusion::Variable::t,1>(qvars.size());
    for(auto& it1: qvars) {
        (*qvars_arr)[it1.second] = _mosek_vars[it1.first.first]->index(it1.first.second);
        DebugOff("\nqvar with index " << it1.second << " is (" << it1.first.first << ", " << it1.first.second << ")");
        DebugOff("\ncorresp mosek vars are " << _mosek_vars[it1.first.first]->index(it1.first.second)->toString());
    }

    Eigen::MatrixXf A = Eigen::MatrixXf::Zero(qvars.size(),qvars.size());
    fusion::Variable::t xi = fusion::Var::vstack(qvars_arr);
//    auto qvarsi = new_array_ptr<fusion::Variable::t,1>(qn);

    size_t j2 = 0;
    for (auto& it1: qterms) {
        double sign;
        if(!it1.second._sign) sign = -2;
        else sign = 2;
        idx1 = it1.second._p->first->get_vec_id();
        idx2 = it1.second._p->second->get_vec_id();
        if ((it1.second._coef->_is_transposed || it1.second._coef->_is_matrix) && !it1.second._p->first->_is_matrix) {
            throw invalid_argument("This type of expression in the constraint is not implemented.");
//            auto dim = it1.second._p->first->get_nb_instances(inst);
//            auto qvarsi = new_array_ptr<fusion::Variable::t,1>(dim);
//            for (int j = 0; j<dim; j++) {
////                    res += poly_eval(_coef,i,j) * poly_eval(_p->first, i,j)* poly_eval(_p->second, i,j);
//                if(it1.second._p->first->get_id_inst(inst,j) != it1.second._p->second->get_id_inst(inst,j))
//                    throw invalid_argument("Bilinear expressions in Mosek objective are not implemented.");
//                A(aidx,aidx) = sign*poly_eval(it1.second._coef,inst,j);
//                aidx++;
//                idx_inst = it1.second._p->first->get_id_inst(inst,j);
////                (*qvarsi)(j) = _mosek_vars[idx1]->index(idx_inst);
//            }
//            xi = fusion::Var::vstack(qvarsi);
        }
        else if(it1.second._p->first->_is_transposed){
//                auto dim = _p->first->get_nb_instances(i);
//                for (int j = 0; j<dim; j++)
//                    res += poly_eval(_coef,i,j) * poly_eval(_p->first, i,j)* poly_eval(_p->second, i,j);
            throw invalid_argument("This type of expression in the constraint is not implemented.");
        }
        else {
//                res = poly_eval(_coef,i) * poly_eval(_p->first, i) * poly_eval(_p->second, i);
            idx_inst = it1.second._p->first->get_id_inst(inst);
            int aidx1 = qvars.find(make_pair(idx1,idx_inst))->second;
            DebugOff("\nvars: " << _mosek_vars[idx1]->index(idx_inst)->toString());
            idx_inst = it1.second._p->second->get_id_inst(inst);
            int aidx2 = qvars.find(make_pair(idx2,idx_inst))->second;
            DebugOff(", " << _mosek_vars[idx2]->index(idx_inst)->toString());
            if(aidx1==aidx2) A(aidx1,aidx2) = sign*poly_eval(it1.second._coef,inst);
            else {
                A(aidx1,aidx2) = 0.5*sign*poly_eval(it1.second._coef,inst);
                A(aidx2,aidx1) = 0.5*sign*poly_eval(it1.second._coef,inst);
            }
//            aidx++;

//            (*qvarsi)(j2) = _mosek_vars[idx1]->index(idx_inst);
//            j2++;
        }
//            if (!it1.second._sign) res *= -1;
    }
    DebugOff("\nA = \n" << A);
    DebugOff("\nx = " << xi->toString());

    Eigen::LLT<Eigen::MatrixXf> lltOfA(A); // compute the Cholesky decomposition of A
    Eigen::MatrixXf L = lltOfA.matrixL();
//    Eigen::MatrixXf L = A.sqrt();
    DebugOff("\nL: " << L);

    std::shared_ptr<ndarray<double, 1>> Farr = new_array_ptr<double,1>(qn*qn);
    for(size_t Fi = 0; Fi < L.rows(); Fi++) {
        for(size_t Fj = 0; Fj < L.cols(); Fj++) {
            (*Farr)[Fi*(L.cols())+Fj] = L(Fi,Fj);
        }
    }
    fusion::Matrix::t F = fusion::Matrix::dense(qvars.size(),qvars.size(),Farr);
    DebugOff("\nF = " << F->toString());
    fusion::Expression::t Fx = fusion::Expr::mul(F, xi);
    DebugOff("\nFx = " << Fx->toString());

    return Fx;
}

void MosekProgram::create_mosek_constraints() {
    size_t idx = 0, idx_inst = 0, idx1 = 0, idx2 = 0, idx_inst1 = 0, idx_inst2 = 0, nb_inst = 0, inst = 0;
    size_t c_idx_inst = 0;
    shared_ptr<Constraint> c;
    for(auto& p: _model->_cons) {
        c = p.second;
        DebugOff("\nconstr " << c->get_name());
        if(c->is_nonlinear() && !c->is_quadratic())
            throw invalid_argument("Mosek cannot handle nonlinear constraints that are not convex quadratic.\n");

        nb_inst = c->get_nb_instances();

        if (c->is_quadratic()) {
            size_t qn = get_num_qterms(c->get_qterms());
            auto fusion_cols = new_array_ptr<fusion::Expression::t,1>(nb_inst);
            for (int i = 0; i< nb_inst; i++) {
                //create matrix F;
                fusion::Expression::t Fx = form_Fx(c->get_qterms(),qn,i);

                //build the linear expression
                fusion::Expression::t lin_expr_i = create_lin_expr(c->get_lterms(), c->get_cst(), i);

                //arrange it all into a row, add the row to an array
                auto Earr = new_array_ptr<fusion::Expression::t,1>(2 + qn);
                (*Earr)[0] = fusion::Expr::constTerm(1);
                (*Earr)[1] = fusion::Expr::mul(lin_expr_i,-1);
                for(size_t Fxi = 0; Fxi < qn; Fxi++) {
                    (*Earr)[Fxi+2] = Fx->index(Fxi);
                }
                fusion::Expression::t qexpr = fusion::Expr::hstack(Earr);
                DebugOff("\nqexpr = " << qexpr->toString());

                (*fusion_cols)[i] = qexpr;
            }

            // use the array of rows to create a matrix
            auto M = fusion::Expr::vstack(fusion_cols);
            DebugOff("\nM = " << M->toString());

            // add a conic constraint
            auto mosek_constr = _mosek_model->constraint(c->get_name(), M, fusion::Domain::inRotatedQCone(nb_inst,2+qn));
            DebugOff("\nConstraint generated: " << mosek_constr->toString());
        } // quadratic expression
        else{
            auto fusion_sums = new_array_ptr<fusion::Expression::t, 1>(nb_inst);
            for (int i = 0; i < nb_inst; i++) {
                fusion::Expression::t lin_expr = create_lin_expr(c->get_lterms(), c->get_cst(), i);
                (*fusion_sums)[i] = lin_expr;
                inst++;
            }
            auto E = fusion::Expr::vstack(fusion_sums);

            DebugOff("\nconstr " << c->get_name() << " in Mosek: " << E->toString());
            if (c->get_type() == geq) {
                DebugOff("\n >= " << c->get_rhs());
                _mosek_model->constraint(c->get_name(), E, fusion::Domain::greaterThan(c->get_rhs()));
            } else if (c->get_type() == leq) {
                DebugOff("\n <= " << c->get_rhs());
                _mosek_model->constraint(c->get_name(), E, fusion::Domain::lessThan(c->get_rhs()));
            } else if (c->get_type() == eq) {
                DebugOff(" = " << c->get_rhs());
                _mosek_model->constraint(c->get_name(), E, fusion::Domain::equalsTo(c->get_rhs()));
            }
        } // linear expression
    } //for(p: _model->_cons)
}

fusion::Expression::t MosekProgram::form_Fx(map<string, qterm>& qterms, size_t qn){
    size_t aidx = 0, idx, idx_inst;
    Eigen::MatrixXf A = Eigen::MatrixXf::Zero(qn,qn);
    fusion::Variable::t xi;
    for (auto& it1: qterms) {
        double sign;
        if(!it1.second._sign) sign = -2;
        else sign = 2;
        idx = it1.second._p->first->get_vec_id();
        if ((it1.second._coef->_is_transposed || it1.second._coef->_is_matrix) && !it1.second._p->first->_is_matrix) {
            auto dim = it1.second._p->first->get_nb_instances();
            auto qvarsi = new_array_ptr<fusion::Variable::t,1>(dim);
            for (int j = 0; j<dim; j++) { //now only supports quadratic (not bilinear) terms
//                    res += poly_eval(_coef,i,j) * poly_eval(_p->first, i,j)* poly_eval(_p->second, i,j);
                if(it1.second._p->first->get_id_inst(j) != it1.second._p->second->get_id_inst(j))
                    throw invalid_argument("Bilinear expressions in Mosek objective are not implemented.");
                A(aidx,aidx) = sign*poly_eval(it1.second._coef,0,j);
                aidx++;
                idx_inst = it1.second._p->first->get_id_inst(j);
                (*qvarsi)(j) = _mosek_vars[idx]->index(idx_inst);
            }
            xi = fusion::Var::vstack(qvarsi);
//                (*qvars)(qtermi) = xi; qtermi++;
        }
        else if(it1.second._p->first->_is_transposed){
//                auto dim = _p->first->get_nb_instances(i);
//                for (int j = 0; j<dim; j++)
//                    res += poly_eval(_coef,i,j) * poly_eval(_p->first, i,j)* poly_eval(_p->second, i,j);
            throw invalid_argument("This type of expression in the objective is not implemented.");
        }
        else {
//                res = poly_eval(_coef,i) * poly_eval(_p->first, i) * poly_eval(_p->second, i);
            throw invalid_argument("This type of expression in the objective is not implemented.");
        }
//            if (!it1.second._sign) res *= -1;
    }
    DebugOff("\nA = " << A);
    DebugOff("\nx = " << xi->toString());

    Eigen::LLT<Eigen::MatrixXf> lltOfA(A); // compute the Cholesky decomposition of A
    Eigen::MatrixXf L = lltOfA.matrixL();
    DebugOff("\nL: " << L);

    std::shared_ptr<ndarray<double, 1>> Farr = new_array_ptr<double,1>(qn*qn);
    for(size_t Fi = 0; Fi < L.rows(); Fi++) {
        for(size_t Fj = 0; Fj < L.cols(); Fj++) {
            (*Farr)[Fi*(L.cols())+Fj] = L(Fi,Fj);
        }
    }
    fusion::Matrix::t F = fusion::Matrix::dense(qn,qn,Farr);
    DebugOff("\nF = " << F->toString());
    fusion::Expression::t Fx = fusion::Expr::mul(F, xi);
    DebugOff("\nFx = " << Fx->toString());

    return Fx;
}

void MosekProgram::set_mosek_objective() {
    size_t idx = 0, idx_inst = 0, idx1 = 0, idx2 = 0, idx_inst1 = 0, idx_inst2 = 0, qn;
    size_t c_idx_inst = 0;
    // initialize with the constant part.
    monty::rc_ptr< ::mosek::fusion::Expression >  expr= fusion::Expr::constTerm(poly_eval(_model->_obj.get_cst())); // expr is a pointer to the Expression.
    if(_model->_obj.get_qterms().empty() == false){
        //cerr << "\nMosek doesn't support quadratic objectives!\n";
        qn = get_num_qterms(_model->_obj.get_qterms());

        _mosek_vars.push_back(_mosek_model->variable("r_obj", 1, fusion::Domain::unbounded()));
        //obj = r + (linear part) (see Mosek modelling cookbook)
        expr = fusion::Expr::add(expr,_mosek_vars[_mosek_vars.size()-1]);
        //auto qvars = new_array_ptr<fusion::Variable::t,1>(_model->_obj.get_qterms().size());
        fusion::Expression::t Fx = form_Fx(_model->_obj.get_qterms(),qn);

        auto Earr = new_array_ptr<fusion::Expression::t,1>(2 + qn);
        (*Earr)[0] = fusion::Expr::constTerm(1);
        (*Earr)[1] = _mosek_vars[_mosek_vars.size()-1]->asExpr();
        for(size_t Fxi = 0; Fxi < qn; Fxi++) {
            (*Earr)[Fxi+2] = Fx->index(Fxi);
        }
        fusion::Expression::t qexpr = fusion::Expr::vstack(Earr);
        DebugOff("\nqexpr = " << qexpr->toString());

        //(1,r,Fx) in Qr
        _mosek_model->constraint(qexpr, fusion::Domain::inRotatedQCone());
    }
    for (auto& it1: _model->_obj.get_lterms()) {
        //idx = it1.second._p->get_id();
        idx = it1.second._p->get_vec_id();
        CType vartype = it1.second._p->get_type();
        if (vartype == var_c) {
            if (it1.second._coef->_is_transposed) {
                auto coefs = new_array_ptr<double,1>((it1.second._p->get_nb_instances()));
                for (int j = 0; j<it1.second._p->get_nb_instances(); j++) {
                    (*coefs)(j)=poly_eval(it1.second._coef, j);
                }
                expr = fusion::Expr::add(expr,fusion::Expr::dot(coefs,_mosek_vars[idx]));
            }
            else {
                // get pos.
                idx_inst = it1.second._p->get_id_inst();
                c_idx_inst = get_poly_id_inst(it1.second._coef);
                //auto t = _mosek_vars[idx]->index
                auto lterm = fusion::Expr::mul(poly_eval(it1.second._coef, c_idx_inst), _mosek_vars[idx]->index(idx_inst));
                if (!it1.second._sign) {
                    lterm = fusion::Expr::mul(-1, lterm);
                }
                expr = fusion::Expr::add(expr,lterm);
            }
        }
        else {
            // sdpvar_c
            if (it1.second._coef->_is_transposed) {
                auto coefs = new_array_ptr<double,1>((it1.second._p->get_nb_instances()));
                for (int j = 0; j<it1.second._p->get_nb_instances(); j++) {
                    (*coefs)(j)=poly_eval(it1.second._coef, j);
                }
                expr = fusion::Expr::add(expr,fusion::Expr::dot(coefs,_mosek_vars[idx]));
            }
            else {
                // retrive the index of the parameter.
                // This is quite different from the general variable definition.
                pair<size_t, size_t>  pair = it1.second._p->get_sdpid();
                c_idx_inst = get_poly_id_inst(it1.second._coef);
                auto lterm = fusion::Expr::mul(poly_eval(it1.second._coef, c_idx_inst), _mosek_vars[idx]->index(pair.first,pair.second));

                if (!it1.second._sign) {
                    lterm = fusion::Expr::mul(-1, lterm);
                }
                expr = fusion::Expr::add(expr,lterm);
            }
        }
    }
    if (_model->_objt == maximize) {
        _mosek_model->objective("obj", mosek::fusion::ObjectiveSense::Maximize, expr);
    }
    else {
        _mosek_model->objective("obj", mosek::fusion::ObjectiveSense::Minimize, expr);
    }
    cout << "\nObj = " << expr->toString() << endl;
}

void MosekProgram::prepare_model() {
    double time_before = get_wall_time();
    fill_in_mosek_vars();
    create_mosek_constraints();
    set_mosek_objective();
    cout << "\nTime spent on building model = " << get_wall_time() - time_before;
    //    print_constraints();
}
