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
}

MosekProgram::MosekProgram(Model* m) {
    _mosek_model = new fusion::Model("noname");
    auto _M = monty::finally([&] () {
        _mosek_model->dispose();
    });
    _model = m;
}

MosekProgram::~MosekProgram() {
    delete _mosek_model;
}

// remain to do
bool MosekProgram::solve(bool relax) {
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
            auto real_var = (var<float>*)v;
            for (int i = 0; i < real_var->get_dim(); i++) {
                _mosek_vars.push_back(_mosek_model->variable(v->get_name()+"_"+to_string(i),fusion::Domain::inRange(real_var->get_lb(i), real_var->get_ub(i))));
            }
            break;
        }
        case long_: {
            auto real_var = (var<long double>*)v;
            for (int i = 0; i < real_var->get_dim(); i++) {
                _mosek_vars.push_back(_mosek_model->variable(v->get_name()+"_"+to_string(i),fusion::Domain::inRange(real_var->get_lb(i), real_var->get_ub(i))));
            }
            break;
        }
        case double_: {
            auto real_var = (var<double>*)v;
            for (int i = 0; i < real_var->get_dim(); i++) {
                _mosek_vars.push_back(_mosek_model->variable(v->get_name()+"_"+to_string(i),fusion::Domain::inRange(real_var->get_lb(i), real_var->get_ub(i))));
            }
            break;
        }
        case integer_: {
            auto real_var = (var<int>*)v;
            for (int i = 0; i < real_var->get_dim(); i++) {
                _mosek_vars.push_back(_mosek_model->variable(v->get_name()+"_"+to_string(i),fusion::Domain::integral(fusion::Domain::inRange(real_var->get_lb(i), real_var->get_ub(i)))));
            }
            break;
        }
        case short_: {
            auto real_var = (var<short>*)v;
            for (int i = 0; i < real_var->get_dim(); i++) {
                _mosek_vars.push_back(_mosek_model->variable(v->get_name()+"_"+to_string(i),fusion::Domain::inRange(real_var->get_lb(i), real_var->get_ub(i))));
            }
            break;
        }
        case binary_: {
            auto real_var = (var<bool>*)v;
            for (int i = 0; i < real_var->get_dim(); i++) {
                _mosek_vars.push_back(_mosek_model->variable(v->get_name()+"_"+to_string(i),fusion::Domain::binary()));
            }
            break;
        }
        default:
            break;
        }
    }
}

void MosekProgram::create_mosek_constraints() {
    size_t idx = 0,  idx_inst = 0, inst = 0;
    size_t nb_inst = 0;
    size_t c_idx_inst = 0;
    Constraint* c;
    for(auto& p: _model->_cons) {
        c = p.second;
        if (c->is_nonlinear()) {
            cout <<  "We haven't implemented quadratic expressions interface for mosek" << endl;
            throw invalid_argument("Mosek cannot handle nonlinear constraints that are not convex quadratic.\n");
        }
        nb_inst = c->get_nb_instances();
        inst = 0;
        for (int i = 0; i< nb_inst; i++) {
            fusion::Expression::t expr; // expr is a pointer to the Expression.
            for (auto& it1: c->get_lterms()) {
                idx = it1.second._p->get_id();
                // question: is_transposed and not transposed.
                if (it1.second._coef->_is_transposed) {
                    auto coefs = new_array_ptr<double,1>((it1.second._p->get_dim()));
                    for (int j = 0; j<it1.second._p->get_dim(); j++) {
                        (*coefs)(j) = poly_eval(it1.second._coef,j);
                    }
                    expr = fusion::Expr::add(expr,fusion::Expr::dot(coefs,_mosek_vars[idx]));
                }
                else {
                    if (!it1.second._sign) {
                        expr = fusion::Expr::add(expr,fusion::Expr::mul(poly_eval(it1.second._coef), _mosek_vars[idx]));
                    }
                    else
                        expr = fusion::Expr::add(expr,fusion::Expr::mul(-poly_eval(it1.second._coef), _mosek_vars[idx]));
                }
            }

            if(c->get_type()==geq) {
                _mosek_model->constraint(c->get_name()+"_"+to_string(i), expr, fusion::Domain::greaterThan(c->get_rhs()));
            }
            else if(c->get_type()==leq) {
                _mosek_model->constraint(c->get_name()+"_"+to_string(i), expr, fusion::Domain::lessThan(c->get_rhs()));
            }
            else if(c->get_type()==eq) {
                _mosek_model->constraint(c->get_name()+"_"+to_string(i), expr, fusion::Domain::equalsTo(c->get_rhs()));
            }
            inst++;
        }
    }
}

void MosekProgram::set_mosek_objective() {
    size_t idx = 0;
    fusion::Expression::t expr; // expr is a pointer to the Expression.
    for (auto& it1: _model->_obj.get_lterms()) {
        idx = it1.second._p->get_id();
        if (it1.second._coef->_is_transposed) {
            auto coefs = new_array_ptr<double,1>((it1.second._p->get_dim()));
           // auto coefs = std::shared_ptr<ndarray <double,1>> (new ndarray<double,1> (shape_t<1>(it1.second._p->get_dim())));
            for (int j = 0; j<it1.second._p->get_dim(); j++) {
                (*coefs)(j)=poly_eval(it1.second._coef, j);
            }
            expr = fusion::Expr::add(expr,fusion::Expr::dot(coefs,_mosek_vars[idx]));
        }
        else {
            idx = it1.second._p->get_id();
            if (!it1.second._sign) {
                expr = fusion::Expr::add(expr,fusion::Expr::mul(-poly_eval(it1.second._coef), _mosek_vars[idx]));
            }
            else
                expr = fusion::Expr::add(expr,fusion::Expr::mul(-poly_eval(it1.second._coef), _mosek_vars[idx]));
        }
    }

    if (_model->_objt == maximize) {
        _mosek_model->objective("obj", mosek::fusion::ObjectiveSense::Maximize, expr);
    }
    else {
        _mosek_model->objective("obj", mosek::fusion::ObjectiveSense::Minimize, expr);
    }

}


void MosekProgram::prepare_model() {
    fill_in_mosek_vars();
    create_mosek_constraints();
    set_mosek_objective();
    //    print_constraints();
}
