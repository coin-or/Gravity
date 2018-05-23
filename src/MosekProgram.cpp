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

void MosekProgram::create_mosek_constraints() {
    size_t idx = 0, idx_inst = 0,idx1 = 0, idx2 = 0, idx_inst1 = 0, idx_inst2 = 0, nb_inst = 0, inst = 0;
    size_t c_idx_inst = 0;
    shared_ptr<Constraint> c;
    for(auto& p: _model->_cons) {
        c = p.second;
        if (c->is_nonlinear()) {
            cout <<  "We haven't implemented quadratic expressions interface for mosek" << endl;
            throw invalid_argument("Mosek cannot handle nonlinear constraints that are not convex quadratic.\n");
        }

        nb_inst = c->get_nb_instances();
        auto fusion_sums = new_array_ptr<fusion::Expression::t,1>(nb_inst);

        for (int i = 0; i< nb_inst; i++) {
            // get constant part first.
            //monty::rc_ptr<mosek::fusion::Expression > expr= fusion::Expr::constTerm(poly_eval(c->get_cst())); // expr is a pointer to the Expression.
            auto fusion_lterms = new_array_ptr<fusion::Expression::t,1>(c->get_lterms().size()+1);
            (*fusion_lterms)[0] = fusion::Expr::constTerm(poly_eval(c->get_cst(),i));
            int j = 1;

            for (auto& it1: c->get_lterms()) {
                idx = it1.second._p->get_vec_id();

                CType vartype = it1.second._p->get_type();
                if (vartype == var_c) {
//                    cout << "\nvar = " << it1.second._p->get_name();
//                    cout << "\nfull Mosek var: " << _mosek_vars[idx]->toString();
                    if(!it1.second._p->_is_matrix) {
                        if (it1.second._coef->_is_transposed) {
                            auto coefs = new_array_ptr<double, 1>(it1.second._p->get_nb_instances(i));
                            auto vars = new_array_ptr<fusion::Variable::t,1>(it1.second._p->get_nb_instances(i));
                            for (int j = 0; j < it1.second._p->get_nb_instances(i); j++) {
                                (*coefs)(j) = poly_eval(it1.second._coef,i,j);
                                idx_inst = it1.second._p->get_id_inst(i,j);
//                                if(idx_inst >= _mosek_vars[idx]->size()){
//                                    cout << "\nidx_inst >= var size: idx_inst = " << idx_inst << ", var size = " << _mosek_vars[idx]->size();
//                                }
                                (*vars)(j) = _mosek_vars[idx]->index(idx_inst);
//                                cout << "\nj, idx, idx_inst = " << j << ", " << idx << ", " << idx_inst << ", var = " <<  _mosek_vars[idx]->index(idx_inst)->toString();
                            }
                            rc_ptr<fusion::Variable> P;
                            if(it1.second._p->get_nb_instances(i)!=1) P = fusion::Var::vstack(vars);
                            else P = _mosek_vars[idx]->index(idx_inst);
//                            cout << "\nCoef = ";
//                            for(auto& coef: *coefs)
//                                cout << coef << "; ";
//                            cout << "\nP = " << P->toString();
                            auto lterm = fusion::Expr::dot(coefs, P);
                            if (!it1.second._sign) lterm = fusion::Expr::mul(-1, lterm);
                            (*fusion_lterms)[j] = lterm;
                            //expr = fusion::Expr::add(expr, fusion::Expr::dot(coefs, _mosek_vars[idx]));
                        } else {
                            if (is_indexed(it1.second._p)) idx_inst = it1.second._p->get_id_inst(i);
                            else idx_inst = inst;

                            if (is_indexed(it1.second._coef)) c_idx_inst = get_poly_id_inst(it1.second._coef);
                            else c_idx_inst = inst;

                            auto lterm = fusion::Expr::mul(poly_eval(it1.second._coef, i), //was c_idx_inst instead of i
                                                           _mosek_vars[idx]->index(idx_inst));

                            if (!it1.second._sign) lterm = fusion::Expr::mul(-1, lterm);
                            (*fusion_lterms)[j] = lterm;
                            //expr = fusion::Expr::add(expr, lterm);
                        }
                    }else{
                        if (it1.second._coef->_is_transposed) {
                            auto coefs = new_array_ptr<double, 1>(it1.second._p->get_nb_instances());
                            for (int j = 0; j < it1.second._p->get_nb_instances(); j++) {
                                (*coefs)(j) = poly_eval(it1.second._coef, j);
                            }
                            (*fusion_lterms)[j] = fusion::Expr::dot(coefs, _mosek_vars[idx]);
                            //expr = fusion::Expr::add(expr, fusion::Expr::dot(coefs, _mosek_vars[idx]));
                        } else {
                            pair<size_t, size_t>  pair;
                            if (is_indexed(it1.second._p)) {
//                                idx_inst = it1.second._p->get_id_inst(i);
                                pair = it1.second._p->get_sdp_inst(i);
                            } else idx_inst = inst;

                            if (is_indexed(it1.second._coef)) c_idx_inst = get_poly_id_inst(it1.second._coef);
                            else c_idx_inst = inst;

                            DebugOff("\nindex: " << pair.first << ", " << pair.second);

                            auto lterm = fusion::Expr::mul(poly_eval(it1.second._coef, c_idx_inst),
                                                           _mosek_vars[idx]->index(pair.first,pair.second));

                            if (!it1.second._sign) lterm = fusion::Expr::mul(-1, lterm);
                            (*fusion_lterms)[j] = lterm;
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
                        (*fusion_lterms)[j] = fusion::Expr::dot(coefs,_mosek_vars[idx]);
                        //expr = fusion::Expr::add(expr,fusion::Expr::dot(coefs,_mosek_vars[idx]));
                    }
                    else { //this results in always taking the first instance? it doesn't use i at all?
                        if (is_indexed(it1.second._p)) idx_inst = it1.second._p->get_id_inst(i); //fixed here, look for similar errors
                        else idx_inst = inst;

                        if (is_indexed(it1.second._coef)) c_idx_inst = get_poly_id_inst(it1.second._coef);
                        else c_idx_inst = inst;

                        auto lterm = fusion::Expr::constTerm(poly_eval(it1.second._coef, c_idx_inst) * poly_eval(it1.second._p, i)); // was idx_inst
                        if (!it1.second._sign) lterm = fusion::Expr::mul(-1, lterm);

                        (*fusion_lterms)[j] = lterm;
//                        expr = fusion::Expr::add(expr,lterm);
                    }
                }
                j++;
            }
            (*fusion_sums)[i] = fusion::Expr::add(std::shared_ptr<ndarray<fusion::Expression::t,1>>(fusion_lterms));
            inst++;
        }
        auto E = fusion::Expr::vstack(fusion_sums);
        DebugOn("\nconstr in Mosek: " << E->toString());
        if(c->get_type()==geq) {
            DebugOn("\n >= " << c->get_rhs());
            _mosek_model->constraint(c->get_name(), E, fusion::Domain::greaterThan(c->get_rhs()));
        }
        else if(c->get_type()==leq) {
            DebugOn("\n <= " << c->get_rhs());
            _mosek_model->constraint(c->get_name(), E, fusion::Domain::lessThan(c->get_rhs()));
        }
        else if(c->get_type()==eq) {
            DebugOn(" = " << c->get_rhs());
            _mosek_model->constraint(c->get_name(), E, fusion::Domain::equalsTo(c->get_rhs()));
        }
    }
}

void MosekProgram::set_mosek_objective() {
    size_t idx = 0, idx_inst = 0, idx1 = 0, idx2 = 0, idx_inst1 = 0, idx_inst2 = 0, qn = 0;
    size_t c_idx_inst = 0;
    // initialize with the constant part.
    monty::rc_ptr< ::mosek::fusion::Expression >  expr= fusion::Expr::constTerm(poly_eval(_model->_obj.get_cst())); // expr is a pointer to the Expression.
    if(_model->_obj.get_qterms().empty() == false){
        //cerr << "\nMosek doesn't support quadratic objectives!\n";
        for (auto& it1: _model->_obj.get_qterms()) {
            if ((it1.second._coef->_is_transposed || it1.second._coef->_is_matrix) && !it1.second._p->first->_is_matrix) {
                auto dim = it1.second._p->first->get_nb_instances();
                qn += dim;
            }
            else if(it1.second._p->first->_is_transposed){
                auto dim = it1.second._p->first->get_nb_instances();
                qn += dim;
            }
            else qn++;
        }

        _mosek_vars.push_back(_mosek_model->variable("r_obj", 1, fusion::Domain::unbounded()));
        //obj = r + (linear part) (see Mosek modelling cookbook)
        expr = fusion::Expr::add(expr,_mosek_vars[_mosek_vars.size()-1]);
        int aidx = 0, qtermi = 0;
        Eigen::MatrixXf A = Eigen::MatrixXf::Zero(qn,qn);
        rc_ptr<fusion::Variable> x;
        auto qvars = new_array_ptr<fusion::Variable::t,1>(_model->_obj.get_qterms().size());
        rc_ptr<fusion::Variable> xi;
        for (auto& it1: _model->_obj.get_qterms()) {
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
                        cerr << "\nWarning: bilinear expressions in Mosek objective are not implemented.";
                    A(aidx,aidx) = sign*poly_eval(it1.second._coef,0,j);
                    aidx++;
                    idx_inst = it1.second._p->first->get_id_inst(j);
                    (*qvarsi)(j) = _mosek_vars[idx]->index(idx_inst);
                }
                xi = fusion::Var::vstack(qvarsi);
                (*qvars)(qtermi) = xi;
                qtermi++;
            }
            else if(it1.second._p->first->_is_transposed){
//                auto dim = _p->first->get_nb_instances(i);
//                for (int j = 0; j<dim; j++) {
//                    res += poly_eval(_coef,i,j) * poly_eval(_p->first, i,j)* poly_eval(_p->second, i,j);
//                }
                cerr << "\nWarning: this type of expression in the objective is not implemented.";
            }
            else {
//                res = poly_eval(_coef,i) * poly_eval(_p->first, i) * poly_eval(_p->second, i);
                cerr << "\nWarning: this type of expression in the objective is not implemented.";
            }
//            if (!it1.second._sign) {
//                res *= -1;
//            }
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
        auto Fx = fusion::Expr::mul(F, xi);
        DebugOff("\nFx = " << Fx->toString());

        auto Earr = new_array_ptr<fusion::Expression::t,1>(2 + xi->size());
        (*Earr)[0] = fusion::Expr::constTerm(1);
        (*Earr)[1] = _mosek_vars[_mosek_vars.size()-1]->asExpr();
        for(size_t Fxi = 0; Fxi < xi->size(); Fxi++) {
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
