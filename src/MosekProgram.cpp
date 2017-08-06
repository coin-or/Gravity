//
//  MosekProgram.cpp
//  Gravity
//
//  Created by Guanglei Wang on 14/7/17.
//
//
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

// remain to do
bool MosekProgram::solve(bool relax) {
    // Set max solution time
    _mosek_model->setSolverParam("mioMaxTime", 360.0);
    // Set max relative gap (to its default value)
    _mosek_model->setSolverParam("mioTolRelGap", 1e-4);
    // Set max absolute gap (to its default value)
    _mosek_model->setSolverParam("mioTolAbsGap", 0.0);
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

    std::cout << "Cost = " << _mosek_model->primalObjValue() << std::endl;

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
        if (_model->_vars[i]->get_intype()== binary_ ||_model->_vars[i]->get_intype()== integer_) {
            for (auto j = 0; j < _model->_vars[i]->get_dim(); j++) {
                auto val = (*sol)[j];
                val = round(val);
                poly_set_val(j, val , _model->_vars[i]);
            }
        }
        else {
            for (auto j = 0; j < _model->_vars[i]->get_dim(); j++)
                poly_set_val(j,(*sol)[j] , _model->_vars[i]);
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
                //for (int i = 0; i < real_var->get_dim(); i++) {
                //    _mosek_vars.push_back(_mosek_model->variable(real_var->get_name()+"_"+to_string(i),fusion::Domain::inRange(real_var->get_lb(i), real_var->get_ub(i))));
                //}
                auto lb  = new_array_ptr<double,1>(real_var->get_dim());
                auto ub  = new_array_ptr<double,1>(real_var->get_dim());
                for (int i = 0; i < real_var->get_dim(); i++) {
                    (*lb)[i] = real_var->get_lb(i);
                    (*ub)[i] = real_var->get_ub(i);
                }
                _mosek_vars.push_back(_mosek_model->variable(real_var->get_name(), real_var->get_dim(), fusion::Domain::inRange(lb,ub)));
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
                //for (int i = 0; i < real_var->get_dim(); i++) {
                //    _mosek_vars.push_back(_mosek_model->variable(real_var->get_name()+"_"+to_string(i),fusion::Domain::inRange(real_var->get_lb(i), real_var->get_ub(i))));
                //}
                auto lb  = new_array_ptr<double,1>(real_var->get_dim());
                auto ub  = new_array_ptr<double,1>(real_var->get_dim());
                for (int i = 0; i < real_var->get_dim(); i++) {
                    (*lb)[i] = real_var->get_lb(i);
                    (*ub)[i] = real_var->get_ub(i);
                }
                _mosek_vars.push_back(_mosek_model->variable(real_var->get_name(), real_var->get_dim(), fusion::Domain::inRange(lb,ub)));
            }
            else {
                auto sdp_var = (sdpvar<long double>*)v;
                _mosek_vars.push_back(_mosek_model->variable(sdp_var->get_name(), fusion::Domain::inPSDCone(sdp_var->get_dim())));
            }
            break;
        }
        case double_: {
            if (v->get_type() == var_c) {
                auto real_var = (var<double>*)v;
                //for (int i = 0; i < real_var->get_dim(); i++) {
                //    _mosek_vars.push_back(_mosek_model->variable(real_var->get_name()+"_"+to_string(i),fusion::Domain::inRange(real_var->get_lb(i), real_var->get_ub(i))));
                //}
                auto lb  = new_array_ptr<double,1>(real_var->get_dim());
                auto ub  = new_array_ptr<double,1>(real_var->get_dim());
                for (int i = 0; i < real_var->get_dim(); i++) {
                    (*lb)[i] = real_var->get_lb(i);
                    (*ub)[i] = real_var->get_ub(i);
                }
                _mosek_vars.push_back(_mosek_model->variable(real_var->get_name(), real_var->get_dim(), fusion::Domain::inRange(lb,ub)));
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
                auto lb  = new_array_ptr<double,1>(real_var->get_dim());
                auto ub  = new_array_ptr<double,1>(real_var->get_dim());
                for (int i = 0; i < real_var->get_dim(); i++) {
                    (*lb)[i] = real_var->get_lb(i);
                    (*ub)[i] = real_var->get_ub(i);
                }
                _mosek_vars.push_back(_mosek_model->variable(real_var->get_name(), real_var->get_dim(), fusion::Domain::integral(fusion::Domain::inRange(lb,ub))));
            }
            else {
                auto sdp_var = (sdpvar<int>*)v;
                _mosek_vars.push_back(_mosek_model->variable(sdp_var->get_name(), fusion::Domain::inPSDCone(sdp_var->get_dim())));
            }
            break;
        }
        case short_: {
            if (v->get_type() == var_c) {
                auto real_var = (var<short>*)v;
                auto lb  = new_array_ptr<double,1>(real_var->get_dim());
                auto ub  = new_array_ptr<double,1>(real_var->get_dim());
                for (int i = 0; i < real_var->get_dim(); i++) {
                    (*lb)[i] = real_var->get_lb(i);
                    (*ub)[i] = real_var->get_ub(i);
                }
                _mosek_vars.push_back(_mosek_model->variable(real_var->get_name(), real_var->get_dim(), fusion::Domain::integral(fusion::Domain::inRange(lb,ub))));
            }
            else {
                auto sdp_var = (sdpvar<short>*)v;
                _mosek_vars.push_back(_mosek_model->variable(sdp_var->get_name(), fusion::Domain::inPSDCone(sdp_var->get_dim())));
            }
            break;
        }
        case binary_: {
            if (v->get_type() == var_c) {
                auto real_var = (var<bool>*)v;
                //for (unsigned int i = 0; i < real_var->get_dim(); i++) {
                //    _mosek_vars.push_back(_mosek_model->variable(real_var->get_name()+"_"+to_string(i),fusion::Domain::binary()));
                //}
                _mosek_vars.push_back(_mosek_model->variable(real_var->get_name(), real_var->get_dim(), fusion::Domain::binary()));
            }
            else {
                auto sdp_var = (sdpvar<bool>*)v;
                _mosek_vars.push_back(_mosek_model->variable(sdp_var->get_name(), fusion::Domain::inPSDCone(sdp_var->get_dim())));
            }
            break;
        }
        default:
            break;
        }
    }
}

void MosekProgram::create_mosek_constraints() {
    //size_t idx = 0;
    //size_t nb_inst = 0;
    size_t idx = 0, idx_inst = 0,idx1 = 0, idx2 = 0, idx_inst1 = 0, idx_inst2 = 0, nb_inst = 0, inst = 0;
    size_t c_idx_inst = 0;
    Constraint* c;
    for(auto& p: _model->_cons) {
        c = p.second;
        //c->print();
        if (c->is_nonlinear()) {
            cout <<  "We haven't implemented quadratic expressions interface for mosek" << endl;
            throw invalid_argument("Mosek cannot handle nonlinear constraints that are not convex quadratic.\n");
        }
        nb_inst = c->get_nb_instances();
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
                    if (it1.second._coef->_is_transposed) {
                        auto coefs = new_array_ptr<double,1>(it1.second._p->get_dim());
                        for (int j = 0; j<it1.second._p->get_dim(); j++) {
                            (*coefs)(j) = poly_eval(it1.second._coef,j);
                        }
                        expr = fusion::Expr::add(expr,fusion::Expr::dot(coefs,_mosek_vars[idx]));
                        //    cout << "expr" << expr->toString() << endl;
                    }
                    else {
                        if (is_indexed(it1.second._p)) {
                            idx_inst = it1.second._p->get_id_inst();
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
                        auto lterm = fusion::Expr::mul(poly_eval(it1.second._coef, c_idx_inst), _mosek_vars[idx]->index(idx_inst));
                        if (!it1.second._sign) {
                            lterm = fusion::Expr::mul(-1, lterm);
                        }
                        expr = fusion::Expr::add(expr,lterm);
                    }
                }
                else {
                    if (it1.second._coef->_is_transposed) {
                        auto coefs = new_array_ptr<double,1>(it1.second._p->get_dim());
                        for (int j = 0; j<it1.second._p->get_dim(); j++) {
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

                }
            }


            if(c->get_type()==geq) {
                _mosek_model->constraint(c->get_name(), expr, fusion::Domain::greaterThan(c->get_rhs()));
            }
            else if(c->get_type()==leq) {
                _mosek_model->constraint(c->get_name(), expr, fusion::Domain::lessThan(c->get_rhs()));
            }
            else if(c->get_type()==eq) {
                _mosek_model->constraint(c->get_name(), expr, fusion::Domain::equalsTo(c->get_rhs()));
            }
            inst++;
        }
    }
}

void MosekProgram::set_mosek_objective() {
    size_t idx = 0, idx_inst = 0, idx1 = 0, idx2 = 0, idx_inst1 = 0, idx_inst2 = 0;
    size_t c_idx_inst = 0;
    // initialize with the constant part.
    monty::rc_ptr< ::mosek::fusion::Expression >  expr= fusion::Expr::constTerm(poly_eval(_model->_obj.get_cst())); // expr is a pointer to the Expression.
    for (auto& it1: _model->_obj.get_lterms()) {
        //idx = it1.second._p->get_id();
        idx = it1.second._p->get_vec_id();
        CType vartype = it1.second._p->get_type();
        if (vartype == var_c) {
            if (it1.second._coef->_is_transposed) {
                auto coefs = new_array_ptr<double,1>((it1.second._p->get_dim()));
                for (int j = 0; j<it1.second._p->get_dim(); j++) {
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
                auto coefs = new_array_ptr<double,1>((it1.second._p->get_dim()));
                for (int j = 0; j<it1.second._p->get_dim(); j++) {
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
    cout << expr->toString() << endl;
}

void MosekProgram::prepare_model() {
    fill_in_mosek_vars();
    create_mosek_constraints();
    set_mosek_objective();
    //    print_constraints();
}
