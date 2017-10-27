//
//  func.cpp
//  Gravity
//
//  Created by Hijazi, Hassan on 24/10/16.
//
//
#include <cmath>
#include <gravity/func.h>

//#define USEDEBUG
#ifdef USEDEBUG
#define Debug(x) cout << x
#else
#define Debug(x)
#endif
#define DebugOn(x) cout << x
#define DebugOff(x)

using namespace std;
namespace gravity{

    bool is_indexed(const constant_* c){
        if (c->is_var() || c->is_param()) {
            return (((param_*)c)->_is_indexed);
        }
        return true;
    }

    size_t get_poly_id(const constant_* c){
        if (c->is_var() || c->is_param()) {
            return ((param_*)c)->get_vec_id();
        }
        return 0;
    }

    size_t get_poly_id_inst(const constant_* c, unsigned inst){
        if (c->is_var() || c->is_param()) {
            return ((param_*)c)->get_id_inst(inst);
        }
        return 0;
    }


    void poly_set_val(unsigned i, double val, param_* p){
        switch (p->get_intype()) {
            case binary_:
                ((param<bool>*)p)->set_val(i, val);
                break;
            case short_:
                ((param<short>*)p)->set_val(i, val);
                break;
            case integer_:
                ((param<int>*)p)->set_val(i, val);
                break;
            case float_:
                ((param<float>*)p)->set_val(i, val);
                break;
            case double_:
                ((param<double>*)p)->set_val(i, val);
                break;
            case long_:
                ((param<long double>*)p)->set_val(i, val);
                break;
            default:
                break;
        }
    }

    double poly_eval(const constant_* c, size_t i){
        if (!c) {
            throw invalid_argument("Cannot evaluate nullptr!");
        }
        switch (c->get_type()) {
            case binary_c: {
                return ((constant<bool>*)(c))->eval();
                break;
            }
            case short_c: {
                return ((constant<short>*)(c))->eval();
                break;
            }
            case integer_c: {
                return ((constant<int>*)(c))->eval();
                break;
            }
            case float_c: {
                return ((constant<float>*)(c))->eval();
                break;
            }
            case double_c: {
                return ((constant<double>*)(c))->eval();
                break;
            }
            case long_c: {
                return (double)((constant<long double>*)(c))->eval();
                break;
            }
            case par_c:{
                auto p_c2 = (param_*)(c);
                switch (p_c2->get_intype()) {
                    case binary_:
                        return ((param<bool>*)p_c2)->eval(i);
                        break;
                    case short_:
                        return ((param<short>*)p_c2)->eval(i);
                        break;
                    case integer_:
                        return ((param<int>*)p_c2)->eval(i);
                        break;
                    case float_:
                        return ((param<float>*)p_c2)->eval(i);
                        break;
                    case double_:
                        return ((param<double>*)p_c2)->eval(i);
                        break;
                    case long_:
                        return (double)((param<long double>*)p_c2)->eval(i);
                        break;
                    default:
                        break;
                }
                break;
            }
            case var_c:{
                auto p_c2 = (param_*)(c);
                switch (p_c2->get_intype()) {
                    case binary_:
                        return ((var<bool>*)p_c2)->eval(i);
                        break;
                    case short_:
                        return ((var<short>*)p_c2)->eval(i);
                        break;
                    case integer_:
                        return ((var<int>*)p_c2)->eval(i);
                        break;
                    case float_:
                        return ((var<float>*)p_c2)->eval(i);
                        break;
                    case double_:
                        return ((var<double>*)p_c2)->eval(i);
                        break;
                    case long_:
                        return (double)((var<long double>*)p_c2)->eval(i);
                        break;
                    default:
                        break;
                }
                break;
            }
            // newly added sdpvar_c (guanglei)    
            case sdpvar_c:{
                auto p_c2 = (param_*)(c);
                switch (p_c2->get_intype()) {
                    case binary_:
                        return ((sdpvar<bool>*)p_c2)->eval(i);
                        break;
                    case short_:
                        return ((sdpvar<short>*)p_c2)->eval(i);
                        break;
                    case integer_:
                        return ((sdpvar<int>*)p_c2)->eval(i);
                        break;
                    case float_:
                        return ((sdpvar<float>*)p_c2)->eval(i);
                        break;
                    case double_:
                        return ((sdpvar<double>*)p_c2)->eval(i);
                        break;
                    case long_:
                        return (double)((sdpvar<long double>*)p_c2)->eval(i);
                        break;
                    default:
                        break;
                }
                break;
            }
            case uexp_c: {
                return ((uexpr*)c)->eval(i);
                break;
            }
            case bexp_c: {
                return ((bexpr*)c)->eval(i);
                break;
            }
            case func_c: {
                return ((func_*)c)->eval(i);
                break;
            }
            default:
                break;
        }
        return 0;
    }
    
    func_ get_poly_derivative(constant_* c, const param_ &v){
        if (!c) {
            throw invalid_argument("Cannot evaluate nullptr!");
        }
        if (c->is_number() || c->is_param()) {
            return func_();
        }
        if (c->is_var()) {
            if ((*(param_*)c)==v) {
                func_() += 1;
            }
            return func_();
        }
        
        if(c->is_uexpr()){
            return ((uexpr*)c)->get_derivative(v);
        }
        if(c->is_bexpr()){
            return ((bexpr*)c)->get_derivative(v);
        }
        if(c->is_function()){
            return ((func_*)c)->get_derivative(v);
        }
        return func_();
    }


    Sign constant_::get_all_sign() const{
        switch (_type) {
            case binary_c: {
                return ((constant<bool>*)this)->get_sign();
                break;
            }
            case short_c: {
                return ((constant<short>*)this)->get_sign();
                break;
            }
            case integer_c: {
                return ((constant<int>*)this)->get_sign();
                break;
            }
            case float_c: {
                return ((constant<float>*)this)->get_sign();
                break;
            }
            case double_c: {
                return ((constant<double>*)this)->get_sign();
                break;
            }
            case long_c: {
                return ((constant<long double>*)this)->get_sign();
                break;
            }
            case par_c:{
                return ((param_*)this)->get_all_sign();
                break;
            }
            case uexp_c: {
                return ((uexpr*)this)->get_all_sign(); // TO UPDATE
                break;
            }
            case bexp_c: {
                return ((bexpr*)this)->get_all_sign(); // TO UPDATE
                break;
            }
            case var_c:{
                return ((param_*)this)->get_all_sign();
                break;
            }
                       //newly added sdpvar (Guanglei)
            case sdpvar_c:{
                return ((param_*)this)->get_all_sign();
                break;
            }
                
            case func_c: {
                return ((func_*)this)->get_all_sign();
                break;
            }
            default:
                break;
        }
        return unknown_;
    }



    Sign constant_::get_sign(int idx) const{
        switch (_type) {
            case binary_c: {
                return ((constant<bool>*)this)->get_sign();
                break;
            }
            case short_c: {
                return ((constant<short>*)this)->get_sign();
                break;
            }
            case integer_c: {
                return ((constant<int>*)this)->get_sign();
                break;
            }
            case float_c: {
                return ((constant<float>*)this)->get_sign();
                break;
            }
            case double_c: {
                return ((constant<double>*)this)->get_sign();
                break;
            }
            case long_c: {
                return ((constant<long double>*)this)->get_sign();
                break;
            }
            case par_c:{
                return ((param_*)this)->get_sign(idx);
                break;
            }
            case uexp_c: {
                return ((uexpr*)this)->get_sign(idx);
                break;
            }
            case bexp_c: {
                return ((bexpr*)this)->get_sign(idx);
                break;
            }
            case var_c:{
                return ((param_*)this)->get_sign(idx);
                break;
            }
                       //newly added sdpvar (guanglei)
            case sdpvar_c:{
                return ((param_*)this)->get_sign(idx);
                break;
            }
                
            case func_c: {
                return ((func_*)this)->get_sign(idx);
                break;
            }
            default:
                break;
        }
        return unknown_;
    }

    Sign param_::get_sign(int idx) const{
        switch (_intype) {
            case binary_:
                if (is_param()) {
                    return ((param<bool>*)this)->get_sign(idx);
                }
                return ((var<bool>*)this)->get_sign(idx);
                break;
            case short_:
                if (is_param()) {
                    return ((param<short>*)this)->get_sign(idx);
                }
                return ((var<short>*)this)->get_sign(idx);
                break;
            case integer_:
                if (is_param()) {
                    return ((param<int>*)this)->get_sign(idx);
                }
                return ((var<int>*)this)->get_sign(idx);
                break;
            case float_:
                if (is_param()) {
                    return ((param<float>*)this)->get_sign(idx);
                }
                return ((var<float>*)this)->get_sign(idx);
                break;
            case double_:
                if (is_param()) {
                    return ((param<double>*)this)->get_sign(idx);
                }
                return ((var<double>*)this)->get_sign(idx);
                break;
            case long_:
                if (is_param()) {
                    return ((param<long double>*)this)->get_sign(idx);
                }
                return ((var<long double>*)this)->get_sign(idx);
                break;
            default:
                break;
        }
    }

    Sign param_::get_all_sign() const{
        switch (_intype) {
            case binary_:
                if (is_param()) {
                    return ((param<bool>*)this)->get_all_sign();
                }
                return ((var<bool>*)this)->get_all_sign();
                break;
            case short_:
                if (is_param()) {
                    return ((param<short>*)this)->get_all_sign();
                }
                return ((var<short>*)this)->get_all_sign();
                break;
            case integer_:
                if (is_param()) {
                    return ((param<int>*)this)->get_all_sign();
                }
                return ((var<int>*)this)->get_all_sign();
                break;
            case float_:
                if (is_param()) {
                    return ((param<float>*)this)->get_all_sign();
                }
                return ((var<float>*)this)->get_all_sign();
                break;
            case double_:
                if (is_param()) {
                    return ((param<double>*)this)->get_all_sign();
                }
                return ((var<double>*)this)->get_all_sign();
                break;
            case long_:
                if (is_param()) {
                    return ((param<long double>*)this)->get_all_sign();
                }
                return ((var<long double>*)this)->get_all_sign();
                break;
            default:
                break;
        }
    }

    constant_* copy(constant_&& c2){/**< Copy c2 into c1 detecting the right class, i.e., constant<>, param<>, uexpr or bexpr. */
        
        switch (c2.get_type()) {
            case binary_c: {
                return new constant<bool>(*(constant<bool>*)move(&c2));
                break;
            }
            case short_c: {
                return new constant<short>(*(constant<short>*)move(&c2));
                break;
            }
            case integer_c: {
                return new constant<int>(*(constant<int>*)move(&c2));
                break;
            }
            case float_c: {
                return new constant<float>(*(constant<float>*)move(&c2));
                break;
            }
            case double_c: {
                return new constant<double>(*(constant<double>*)move(&c2));
                break;
            }
            case long_c: {
                return new constant<long double>(*(constant<long double>*)move(&c2));
                break;
            }
            case par_c:{
                auto p_c2 = (param_*)(&c2);
                switch (p_c2->get_intype()) {
                    case binary_:
                        return new param<bool>(*(param<bool>*)move(p_c2));
                        break;
                    case short_:
                        return new param<short>(*(param<short>*)move(p_c2));
                        break;
                    case integer_:
                        return new param<int>(*(param<int>*)move(p_c2));
                        break;
                    case float_:
                        return new param<float>(*(param<float>*)move(p_c2));
                        break;
                    case double_:
                        return new param<double>(*(param<double>*)move(p_c2));
                        break;
                    case long_:
                        return new param<long double>(*(param<long double>*)move(p_c2));
                        break;
                    default:
                        break;
                }
                break;
            }
            case uexp_c: {
                return new uexpr(*(uexpr*)move(&c2));
                break;
            }
            case bexp_c: {
                return new bexpr(*(bexpr*)move(&c2));
                break;
            }
            case var_c:{
                auto p_c2 = (param_*)(&c2);
                switch (p_c2->get_intype()) {
                    case binary_:
                        return new var<bool>(*(var<bool>*)move(p_c2));
                        break;
                    case short_:
                        return new var<short>(*(var<short>*)move(p_c2));
                        break;
                    case integer_:
                        return new var<int>(*(var<int>*)move(p_c2));
                        break;
                    case float_:
                        return new var<float>(*(var<float>*)move(p_c2));
                        break;
                    case double_:
                        return new var<double>(*(var<double>*)move(p_c2));
                        break;
                    case long_:
                        return new var<long double>(*(var<long double>*)move(p_c2));
                        break;
                    default:
                        break;
                }
                break;
            }
    // new added sdpvar
            case sdpvar_c:{
                auto p_c2 = (param_*)(&c2);
                switch (p_c2->get_intype()) {
                    case binary_:
                        return new sdpvar<bool>(*(sdpvar<bool>*)move(p_c2));
                        break;
                    case short_:
                        return new sdpvar<short>(*(sdpvar<short>*)move(p_c2));
                        break;
                    case integer_:
                        return new sdpvar<int>(*(sdpvar<int>*)move(p_c2));
                        break;
                    case float_:
                        return new sdpvar<float>(*(sdpvar<float>*)move(p_c2));
                        break;
                    case double_:
                        return new sdpvar<double>(*(sdpvar<double>*)move(p_c2));
                        break;
                    case long_:
                        return new sdpvar<long double>(*(sdpvar<long double>*)move(p_c2));
                        break;
                    default:
                        break;
                }
                break;
            }
                
            case func_c: {
                return new func_(move(c2));
                break;
            }
                
            default:
                break;
        }
        return nullptr;
    }

    constant_* copy(const constant_& c2){/**< Copy c2 into c1 detecting the right class, i.e., constant<>, param<>, uexpr or bexpr. */
        
        switch (c2.get_type()) {
            case binary_c: {
                return new constant<bool>(*(constant<bool>*)(&c2));
                break;
            }
            case short_c: {
                return new constant<short>(*(constant<short>*)(&c2));
                break;
            }
            case integer_c: {
                return new constant<int>(*(constant<int>*)(&c2));
                break;
            }
            case float_c: {
                return new constant<float>(*(constant<float>*)(&c2));
                break;
            }
            case double_c: {
                return new constant<double>(*(constant<double>*)(&c2));
                break;
            }
            case long_c: {
                return new constant<long double>(*(constant<long double>*)(&c2));
                break;
            }
            case par_c:{
                auto p_c2 = (param_*)(&c2);
                switch (p_c2->get_intype()) {
                    case binary_:
                        return new param<bool>(*(param<bool>*)p_c2);
                        break;
                    case short_:
                        return new param<short>(*(param<short>*)p_c2);
                        break;
                    case integer_:
                        return new param<int>(*(param<int>*)p_c2);
                        break;
                    case float_:
                        return new param<float>(*(param<float>*)p_c2);
                        break;
                    case double_:
                        return new param<double>(*(param<double>*)p_c2);
                        break;
                    case long_:
                        return new param<long double>(*(param<long double>*)p_c2);
                        break;
                    default:
                        break;
                }
                break;
            }
            case uexp_c: {
                return new uexpr(*(uexpr*)&c2);
                break;
            }
            case bexp_c: {
                return new bexpr(*(bexpr*)&c2);
                break;
            }
            case var_c:{
                auto p_c2 = (param_*)(&c2);
                switch (p_c2->get_intype()) {
                    case binary_:
                        return new var<bool>(*(var<bool>*)p_c2);
                        break;
                    case short_:
                        return new var<short>(*(var<short>*)p_c2);
                        break;
                    case integer_:
                        return new var<int>(*(var<int>*)p_c2);
                        break;
                    case float_:
                        return new var<float>(*(var<float>*)p_c2);
                        break;
                    case double_:
                        return new var<double>(*(var<double>*)p_c2);
                        break;
                    case long_:
                        return new var<long double>(*(var<long double>*)p_c2);
                        break;
                    default:
                        break;
                }
                break;
            }
            case sdpvar_c:{
                auto p_c2 = (param_*)(&c2);
                switch (p_c2->get_intype()) {
                    case binary_:
                        return new sdpvar<bool>(*(sdpvar<bool>*)p_c2);
                        break;
                    case short_:
                        return new sdpvar<short>(*(sdpvar<short>*)p_c2);
                        break;
                    case integer_:
                        return new sdpvar<int>(*(sdpvar<int>*)p_c2);
                        break;
                    case float_:
                        return new sdpvar<float>(*(sdpvar<float>*)p_c2);
                        break;
                    case double_:
                        return new sdpvar<double>(*(sdpvar<double>*)p_c2);
                        break;
                    case long_:
                        return new sdpvar<long double>(*(sdpvar<long double>*)p_c2);
                        break;
                    default:
                        break;
                }
                break;
            }
                
            case func_c: {
                return new func_(c2);
                break;
            }
                
            default:
                break;
        }
        return nullptr;
    }

    double lterm::eval(size_t i) const{
        double res = 0;
        if (_coef->_is_transposed) {
            for (int j = 0; j<_p->_dim; j++) {
                res += poly_eval(_coef,j) * poly_eval(_p, j);
            }
        }
        else {
            res = poly_eval(_coef,i) * poly_eval(_p, i);
        }
        if (!_sign) {
            res *= -1;
        }
        return res;
    }


    double qterm::eval(size_t i) const{
        double res = 0;
        if (_coef->_is_transposed) {
            for (int j = 0; j<_p->first->_dim; j++) {
                res += poly_eval(_coef,j) * poly_eval(_p->first, j)* poly_eval(_p->second, j);
            }
        }
        else if(_p->first->_is_transposed){
            for (int j = 0; j<_p->first->_dim; j++) {
                res += poly_eval(_coef,i) * poly_eval(_p->first, j)* poly_eval(_p->second, j);
            }
        }
        else {
            res = poly_eval(_coef,i) * poly_eval(_p->first, i) * poly_eval(_p->second, i);
        }
        if (!_sign) {
            res *= -1;
        }
        return res;
    }


    double pterm::eval(size_t i) const{
        double res = 0;
        if (_coef->_is_transposed) {
            double pterm = 0;
            for (int j = 0; j<_l->front().first->_dim; j++) {
                pterm = 1;
                for (auto &pair: *_l) {
                    pterm *= pow(poly_eval(pair.first, j), pair.second);
                }
                pterm *= poly_eval(_coef,j);
                res += pterm;
            }
        }// TREAT TRANSPOSED VECTORS IN POLYNOMIAL TERMS HERE
        else {
            res =1;
            for (auto &pair: *_l) {
                res *= pow(poly_eval(pair.first, i), pair.second);
            }
            
            res *= poly_eval(_coef,i);
        }
        if (!_sign) {
            res *= -1;
        }
        return res;
    }


    lterm::lterm(bool sign, constant_* coef, param_* p){
        _coef = coef;
        _p = p;
        _sign = sign;
        if (coef->_is_transposed && p->_is_transposed) {
            throw invalid_argument("Check the transpose operator, there seems to be a dimension issue\n");
        }
        if (coef->is_function()) {
            assert(((func_*)coef)->is_constant());
        }
    };


    lterm& lterm::operator=(const lterm &l){
        if (_coef) {
            delete _coef;
        }
        if (_p) {
            delete _p;
        }
        _coef = copy(*l._coef);
        _p = (param_*)copy(*l._p);
        _sign = l._sign;
        return *this;
    }

    lterm& lterm::operator=(lterm&& l){
        if (_coef) {
            delete _coef;
        }
        if (_p) {
            delete _p;
        }
        _coef = l._coef;
        l._coef = nullptr;
        _p = l._p;
        l._p = nullptr;
        _sign = l._sign;
        return *this;
    }

    qterm& qterm::operator=(const qterm &q){
        if (_coef) {
            delete _coef;
        }
        if (_p) {
            delete _p;
        }
        _coef = copy(*q._coef);
        _p = new pair<param_*, param_*>(make_pair<>((param_*)copy(*q._p->first), (param_*)copy(*q._p->second)));
        _sign = q._sign;
        return *this;
    }


    qterm& qterm::operator=(qterm&& q){
        if (_coef) {
            delete _coef;
        }
        if (_p) {
            delete _p;
        }
        _coef = q._coef;
        q._coef = nullptr;
        _p = q._p;
        q._p = nullptr;
        _sign = q._sign;
        return *this;
    }



    pterm& pterm::operator=(const pterm &p){
        if (_coef) {
            delete _coef;
        }
        if (_l) {
            delete _l;
        }
        _coef = copy(*p._coef);
        _l = new list<pair<param_*, int>>();
        for (auto &pair : *p._l) {
            _l->push_back(make_pair<>((param_*)copy(*pair.first), pair.second));
        }
        _sign = p._sign;
        return *this;
    }


    pterm& pterm::operator=(pterm&& p){
        if (_coef) {
            delete _coef;
        }
        if (_l) {
            delete _l;
        }
        _coef = p._coef;
        p._coef = nullptr;
        _l = p._l;
        p._l = nullptr;
        _sign = p._sign;
        return *this;
    }



    func_::func_(){
        set_type(func_c);
        _params = new map<string, pair<param_*, int>>();
        _vars = new map<string, pair<param_*, int>>();
        _cst = new constant<double>(0);
        _lterms = new map<string, lterm>();
        _qterms = new map<string, qterm>();
        _pterms = new map<string, pterm>();
        _expr = nullptr;
        _DAG = new map<string, expr*>();
        _queue = new queue<expr*>();
        _all_sign = zero_;
        _all_convexity = linear_;
        _all_range = new pair<constant_*, constant_*>(new constant<double>(0), new constant<double>(0));
        _sign = nullptr;
        _convexity = nullptr;
        _range = nullptr;
    };

    func_::func_(constant_&& c){
        if(c.is_function()){
            set_type(func_c);
            auto f = (func_*) &c;
            _ftype = f->_ftype;
            _return_type = f->_return_type;
            _all_convexity = f->_all_convexity;
            _all_sign = f->_all_sign;
            _all_range = f->_all_range;
            f->_all_range = nullptr;
            _lterms = f->_lterms;
            f->_lterms = nullptr;
            _qterms = f->_qterms;
            f->_qterms = nullptr;
            _pterms = f->_pterms;
            f->_pterms = nullptr;
            _expr = f->_expr;
            f->_expr = nullptr;
            _DAG = f->_DAG;
            f->_DAG = nullptr;
            _queue = f->_queue;
            f->_queue = nullptr;
            _vars = f->_vars;
            f->_vars = nullptr;
            _params = f->_params;
            f->_params = nullptr;
            _cst = f->_cst;
            f->_cst = nullptr;
            _range = f->_range;
            f->_range = nullptr;
            _convexity = f->_convexity;
            f->_convexity = nullptr;
            _sign = f->_sign;
            f->_sign = nullptr;
            _is_transposed = f->_is_transposed;
            _is_vector = f->_is_vector;
            _embedded = f->_embedded;
            for (auto &df:f->_dfdx) {
                _dfdx[df.first] = df.second;
            }
            f->_dfdx.clear();
            _nb_instances = f->_nb_instances;
            _nnz_j = f->_nnz_j;
            _nnz_h = f->_nnz_h;
            _hess_link = f->_hess_link;
            
        }
        else {
            set_type(func_c);
            _params = new map<string, pair<param_*, int>>();
            _vars = new map<string, pair<param_*, int>>();
            _cst = nullptr;
            _lterms = new map<string, lterm>();
            _qterms = new map<string, qterm>();
            _pterms = new map<string, pterm>();
            _expr = nullptr;
            _DAG = new map<string, expr*>();
            _queue = new queue<expr*>();
            _all_sign = zero_;
            _all_convexity = linear_;
            _all_range = nullptr;
            _sign = nullptr;
            _convexity = nullptr;
            _range = nullptr;
            
            switch (c.get_type()) {
                case binary_c: {
                    _cst = new constant<bool>(*(constant<bool>*)(&c));
                    _all_sign = ((constant<bool>*)_cst)->get_sign();
                    _all_range = new pair<constant_*, constant_*>(copy(*_cst),copy(*_cst));
                    break;
                }
                case short_c: {
                    _cst = new constant<short>(*(constant<short>*)(&c));
                    _all_sign = ((constant<short>*)_cst)->get_sign();
                    _all_range = new pair<constant_*, constant_*>(copy(*_cst),copy(*_cst));
                    break;
                }
                case integer_c: {
                    _cst = new constant<int>(*(constant<int>*)(&c));
                    _all_sign = ((constant<int>*)_cst)->get_sign();
                    _all_range = new pair<constant_*, constant_*>(copy(*_cst),copy(*_cst));
                    break;
                }
                case float_c: {
                    _cst = new constant<float>(*(constant<float>*)(&c));
                    _all_sign = ((constant<float>*)_cst)->get_sign();
                    _all_range = new pair<constant_*, constant_*>(copy(*_cst),copy(*_cst));
                    break;
                }
                case double_c: {
                    _cst = new constant<double>(*(constant<double>*)(&c));
                    _all_sign = ((constant<double>*)_cst)->get_sign();
                    _all_range = new pair<constant_*, constant_*>(copy(*_cst),copy(*_cst));
                    break;
                }
                case long_c: {
                    _cst = new constant<long double>(*(constant<long double>*)(&c));
                    _all_sign = ((constant<long double>*)_cst)->get_sign();
                    _all_range = new pair<constant_*, constant_*>(copy(*_cst),copy(*_cst));
                    break;
                }
                case par_c:{
                    auto p_c2 = (param_*)copy(move(c));
                    _lterms->insert(make_pair<>(p_c2->get_name(), p_c2));
                    add_param(p_c2);
                    _cst = new constant<double>(0);
                    _all_sign = p_c2->get_all_sign();
                    _all_range = p_c2->get_range();
                    _is_transposed = p_c2->_is_transposed;
                    _is_vector = p_c2->_is_vector;
                    break;
                }
                case var_c:{
                    auto p_c2 = (param_*)copy(move(c));
                    _lterms->insert(make_pair<>(p_c2->get_name(), p_c2));
                    add_var(p_c2);
                    _ftype = lin_;
                    _cst = new constant<double>(0);
                    _all_sign = p_c2->get_all_sign();
                    _all_range = p_c2->get_range();
                    _is_transposed = p_c2->_is_transposed;
                    _is_vector = p_c2->_is_vector;
                    break;
                }
                           // newly added
                case sdpvar_c:{
                    auto p_c2 = (param_*)copy(move(c));
                    _lterms->insert(make_pair<>(p_c2->get_name(), p_c2));
                    add_var(p_c2);
                    _ftype = lin_;
                    _cst = new constant<double>(0);
                    _all_sign = p_c2->get_all_sign();
                    _all_range = p_c2->get_range();
                    _is_transposed = p_c2->_is_transposed;
                    _is_vector = p_c2->_is_vector;
                    break;
                }
                case uexp_c: {
                    _expr = new uexpr(*(uexpr*)move(&c));
                    embed(*_expr);
//                    _DAG->insert(make_pair<>(_expr->get_str(), _expr));
//                    _queue->push(_expr);
                    //sign and convexity
                    break;
                }
                case bexp_c: {
                    _expr = new bexpr(*(bexpr*)move(&c));
                    embed(*_expr);
//                    _DAG->insert(make_pair<>(_expr->get_str(), _expr));
//                    _queue->push(_expr);
                    break;
                }
                default:
                    break;
            }
        }
    }



    func_::func_(const constant_& c){
        set_type(func_c);
        _params = new map<string, pair<param_*, int>>();
        _vars = new map<string, pair<param_*, int>>();
        _cst = nullptr;
        _lterms = new map<string, lterm>();
        _qterms = new map<string, qterm>();
        _pterms = new map<string, pterm>();
        _expr = nullptr;
        _DAG = new map<string, expr*>();
        _queue = new queue<expr*>();
        _all_sign = zero_;
        _all_convexity = linear_;
        _all_range = nullptr;
        _sign = nullptr;
        _convexity = nullptr;
        _range = nullptr;
        _is_transposed = c._is_transposed;
        _is_vector = c._is_vector;

        switch (c.get_type()) {
            case binary_c: {
                _cst = new constant<bool>(*(constant<bool>*)(&c));
                _all_sign = ((constant<bool>*)_cst)->get_sign();
                _all_range = new pair<constant_*, constant_*>(copy(*_cst),copy(*_cst));
                break;
            }
            case short_c: {
                _cst = new constant<short>(*(constant<short>*)(&c));
                _all_sign = ((constant<short>*)_cst)->get_sign();
                _all_range = new pair<constant_*, constant_*>(copy(*_cst),copy(*_cst));
                break;
            }
            case integer_c: {
                _cst = new constant<int>(*(constant<int>*)(&c));
                _all_sign = ((constant<int>*)_cst)->get_sign();
                _all_range = new pair<constant_*, constant_*>(copy(*_cst),copy(*_cst));
                break;
            }
            case float_c: {
                _cst = new constant<float>(*(constant<float>*)(&c));
                _all_sign = ((constant<float>*)_cst)->get_sign();
                _all_range = new pair<constant_*, constant_*>(copy(*_cst),copy(*_cst));
                break;
            }
            case double_c: {
                _cst = new constant<double>(*(constant<double>*)(&c));
                _all_sign = ((constant<double>*)_cst)->get_sign();
                _all_range = new pair<constant_*, constant_*>(copy(*_cst),copy(*_cst));
                break;
            }
            case long_c: {
                _cst = new constant<long double>(*(constant<long double>*)(&c));
                _all_sign = ((constant<long double>*)_cst)->get_sign();
                _all_range = new pair<constant_*, constant_*>(copy(*_cst),copy(*_cst));
                break;
            }
            case par_c:{
                auto p_c2 = (param_*)copy(c);
                _lterms->insert(make_pair<>(p_c2->get_name(), p_c2));
                add_param(p_c2);
                _cst = new constant<double>(0);
                _all_sign = p_c2->get_all_sign();
                _all_range = p_c2->get_range();
                break;
            }
            case var_c:{
                auto p_c2 = (param_*)copy(c);
                _lterms->insert(make_pair<>(p_c2->get_name(), p_c2));
                add_var(p_c2);
                _ftype = lin_;
                _cst = new constant<double>(0);
                _all_sign = p_c2->get_all_sign();
                _all_range = p_c2->get_range();
                break;
            }
            case sdpvar_c:{
                auto p_c2 = (param_*)copy(c);
                _lterms->insert(make_pair<>(p_c2->get_name(), p_c2));
                add_var(p_c2);
                _ftype = lin_;
                _cst = new constant<double>(0);
                _all_sign = p_c2->get_all_sign();
                _all_range = p_c2->get_range();
                break;
            }
            case uexp_c: {
                _cst = new constant<double>(0);
                auto ue = (uexpr*)(&c);
                switch (ue->_otype) {
                    case sin_:
                        _all_range = new pair<constant_*, constant_*>(new constant<double>(-1),new constant<double>(1)); // TO UPDATE
                        _all_sign = unknown_;
                        break;
                    case cos_:
                        _all_range = new pair<constant_*, constant_*>(new constant<double>(-1),new constant<double>(1)); // TO UPDATE
                        _all_sign = unknown_;
                        break;
                    case sqrt_:
                        _all_range = new pair<constant_*, constant_*>(new constant<double>(0),new constant<double>(numeric_limits<double>::max())); // TO UPDATE
                        _all_sign = non_neg_;
                        break;
                    case exp_:
                        _all_range = new pair<constant_*, constant_*>(new constant<double>(numeric_limits<double>::lowest()),new constant<double>(numeric_limits<double>::max())); // TO UPDATE
                        _all_sign = pos_;
                        break;
                    case log_:
                        _all_range = new pair<constant_*, constant_*>(new constant<double>(numeric_limits<double>::lowest()),new constant<double>(numeric_limits<double>::max())); // TO UPDATE
                        _all_sign = unknown_;
                        break;
                    default:
                        break;
                }
                _expr = new uexpr(*ue);
                embed(*_expr);
                if (!_vars->empty()) {
                    _ftype = nlin_;
                }
//                _DAG->insert(make_pair<>(_expr->get_str(), _expr));
//                _queue->push(_expr);
                //sign and convexity
                break;
            }
            case bexp_c: {
                auto be = (bexpr*)&c;
                _cst = new constant<double>(0);
                _expr = new bexpr(*(bexpr*)&c);
                _all_range = new pair<constant_*, constant_*>(new constant<double>(numeric_limits<double>::lowest()),new constant<double>(numeric_limits<double>::max())); // TO UPDATE
                _all_sign = be->get_all_sign();
                embed(*_expr);
                if (!_vars->empty()) {
                    _ftype = nlin_;
                }
//                _DAG->insert(make_pair<>(_expr->get_str(), _expr));
//                _queue->push(_expr);
                break;
            }
            case func_c: {
                *this = *(func_*)&c;
            }
            default:
                break;
        }
    }

    func_::func_(func_&& f){
        set_type(func_c);
        _ftype = f._ftype;
        _return_type = f._return_type;
        _all_convexity = f._all_convexity;
        _all_sign = f._all_sign;
        _all_range = f._all_range;
        f._all_range = nullptr;
        _lterms = f._lterms;
        f._lterms = nullptr;
        _qterms = f._qterms;
        f._qterms = nullptr;
        _pterms = f._pterms;
        f._pterms = nullptr;
        _expr = f._expr;
        f._expr = nullptr;
        _DAG = f._DAG;
        f._DAG = nullptr;
        _queue = f._queue;
        f._queue = nullptr;
        _vars = f._vars;
        f._vars = nullptr;
        _params = f._params;
        f._params = nullptr;
        _cst = f._cst;
        f._cst = nullptr;
        _range = f._range;
        f._range = nullptr;
        _convexity = f._convexity;
        f._convexity = nullptr;
        _sign = f._sign;
        f._sign = nullptr;
        _is_transposed = f._is_transposed;
        _is_vector = f._is_vector;
        _embedded = f._embedded;
        for (auto &df:f._dfdx) {
            _dfdx[df.first] = df.second;
        }
        f._dfdx.clear();
        _nnz_j = f._nnz_j;
        _nnz_h = f._nnz_h;
        _hess_link = f._hess_link;
        _nb_instances = f._nb_instances;
    }

    func_::func_(const func_& f){
        set_type(func_c);
        _params = new map<string, pair<param_*, int>>();
        _vars = new map<string, pair<param_*, int>>();
        _lterms = new map<string, lterm>();
        _qterms = new map<string, qterm>();
        _pterms = new map<string, pterm>();
        _expr = nullptr;
        _DAG = new map<string, expr*>();
        _queue  = new queue<expr*>();
        _ftype = f._ftype;
        _return_type = f._return_type;
        _all_convexity = f._all_convexity;
        _all_sign = f._all_sign;
        _is_transposed = f._is_transposed;
        _is_vector = f._is_vector;
        _embedded = f._embedded;
        _cst = copy(*f._cst);
        _all_range = new pair<constant_*, constant_*>(copy(*f._all_range->first), copy(*f._all_range->second));
        _sign = nullptr;
        _convexity = nullptr;
        _range = nullptr;
        
        if(f._sign){
            _sign = new vector<Sign>();
            *_sign = *f._sign;
        }
        if (f._convexity) {
            _convexity = new vector<Convexity>();
            *_convexity = *f._convexity;
        }
        if (f._range) {
            _range= new vector<pair<constant_*, constant_*>>();
            *_range = *f._range;// Make sure this creates new pointers inside each pair, otherwise use below.
    //        for (auto &elem: *f._range) {
    //            _range->push_back(make_pair<>(copy(elem.first), copy(elem.second)));
    //        }
        }
        
        for (auto &pair:*f._lterms) {
            insert(pair.second);
        }
        for (auto &pair:*f._qterms) {
            insert(pair.second);
        }
        for (auto &pair:*f._pterms) {
            insert(pair.second);
        }
        if (f._expr) {
            _expr = (expr*)copy(*f._expr);
            embed(*_expr);
//            _DAG->insert(make_pair<>(_expr->get_str(), _expr));
//            _queue->push(_expr);
        }
        for (auto &df:f._dfdx) {
            _dfdx[df.first] = new func_(*df.second);
        }
        _nnz_j = f._nnz_j;
        _nnz_h = f._nnz_h;
        _hess_link = f._hess_link;
        _nb_instances = f._nb_instances;

    }

    bool func_::operator==(const func_& f) const{
        if (_ftype!=f._ftype || _all_sign!=f._all_sign || _vars->size()!=f._vars->size() || _pterms->size() != f._pterms->size()|| _qterms->size() != f._qterms->size() || _lterms->size() != f._lterms->size()) {
            return false;
        }
        if (this->to_str()!=f.to_str()) {
            return false;
        }
        return true;
    }

    func_& func_::operator=(const func_& f){
        if (_all_range) {
            delete _all_range->first;
            delete _all_range->second;
        }
        delete _all_range;
        if (_vars) {
            if (!_embedded) {
                for (auto &elem: *_vars) {
                    delete elem.second.first;
                }
            }
            delete _vars;
        }
        if (_params) {
            if (!_embedded) {
                for (auto &elem: *_params) {
                    delete elem.second.first;
                }
            }
            delete _params;
        }
        if (_range) {
            for (auto &elem: *_range) {
                delete elem.first;
                delete elem.second;
            }
        }
        for (auto &f_p: _dfdx) {
            delete f_p.second;
        }
        _dfdx.clear();
        delete _range;
        delete _convexity;
        delete _sign;
        delete _lterms;
        delete _qterms;
        delete _pterms;
        delete _expr;
        delete _DAG;
        delete _queue;
        delete _cst;
        
        set_type(func_c);
        _params = new map<string, pair<param_*, int>>();
        _vars = new map<string, pair<param_*, int>>();
        _lterms = new map<string, lterm>();
        _qterms = new map<string, qterm>();
        _pterms = new map<string, pterm>();
        _DAG = new map<string, expr*>();
        _queue = new queue<expr*>();
        _ftype = f._ftype;
        _return_type = f._return_type;
        _all_convexity = f._all_convexity;
        _all_sign = f._all_sign;
        _is_transposed = f._is_transposed;
        _is_vector = f._is_vector;
        _embedded = f._embedded;
        _return_type = f._return_type;
        _cst = copy(*f._cst);
        _all_range = new pair<constant_*, constant_*>(copy(*f._all_range->first), copy(*f._all_range->second));
        _sign = nullptr;
        _convexity = nullptr;
        _range = nullptr;
        
        if(f._sign){
            _sign = new vector<Sign>();
            *_sign = *f._sign;
        }
        if (f._convexity) {
            _convexity = new vector<Convexity>();
            *_convexity = *f._convexity;
        }
        if (f._range) {
            _range= new vector<pair<constant_*, constant_*>>();
            *_range = *f._range;// Make sure this creates new pointers inside each pair, otherwise use below.
            //        for (auto &elem: *f._range) {
            //            _range->push_back(make_pair<>(copy(elem.first), copy(elem.second)));
            //        }
        }

        for (auto &pair:*f._lterms) {
            insert(pair.second);
        }
        for (auto &pair:*f._qterms) {
            insert(pair.second);
        }
        for (auto &pair:*f._pterms) {
            insert(pair.second);
        }
        if (f._expr) {
            _expr = (expr*)copy(*f._expr);
            embed(*_expr);
//            _DAG->insert(make_pair<>(_expr->get_str(), _expr));
//            _queue->push(_expr);
        }
        for (auto &df:_dfdx) {
            _dfdx[df.first] = new func_(*df.second);
        }
        _nnz_j = f._nnz_j;
        _nnz_h = f._nnz_h;
        _hess_link = f._hess_link;
        _nb_instances = f._nb_instances;
        return *this;
    }


    func_& func_::operator=(func_&& f){
        if (_all_range) {
            delete _all_range->first;
            delete _all_range->second;
        }
        delete _all_range;
        if (_vars) {
            if (!_embedded) {
                for (auto &elem: *_vars) {
                    delete elem.second.first;
                }
            }
            delete _vars;
        }
        if (_params) {
            if (!_embedded) {
                for (auto &elem: *_params) {
                    delete elem.second.first;
                }
            }
            delete _params;
        }
        if (_range) {
            for (auto &elem: *_range) {
                delete elem.first;
                delete elem.second;
            }
        }
        for (auto &f_p: _dfdx) {
            delete f_p.second;
        }
        _dfdx.clear();
        delete _range;
        delete _convexity;
        delete _sign;
        delete _lterms;
        delete _qterms;
        delete _pterms;
        delete _expr;
        delete _DAG;
        delete _queue;
        delete _cst;
        set_type(func_c);
        _ftype = f._ftype;
        _return_type = f._return_type;
        _all_convexity = f._all_convexity;
        _all_sign = f._all_sign;
        _all_range = f._all_range;
        f._all_range = nullptr;
        _lterms = f._lterms;
        f._lterms = nullptr;
        _qterms = f._qterms;
        f._qterms = nullptr;
        _pterms = f._pterms;
        f._pterms = nullptr;
        _expr = f._expr;
        f._expr = nullptr;
        _DAG = f._DAG;
        f._DAG = nullptr;
        _queue = f._queue;
        f._queue = nullptr;
        _vars = f._vars;
        f._vars = nullptr;
        _params = f._params;
        f._params = nullptr;
        _cst = f._cst;
        f._cst = nullptr;
        _range = f._range;
        f._range = nullptr;
        _convexity = f._convexity;
        f._convexity = nullptr;
        _sign = f._sign;
        f._sign = nullptr;
        _is_transposed = f._is_transposed;
        _is_vector = f._is_vector;
        _embedded = f._embedded;
        for (auto &df:f._dfdx) {
            _dfdx[df.first] = df.second;
        }
        f._dfdx.clear();
        _nnz_j = f._nnz_j;
        _nnz_h = f._nnz_h;
        _hess_link = f._hess_link;
        _nb_instances = f._nb_instances;
        return *this;
    }



    func_::~func_(){
        if (_all_range) {
            delete _all_range->first;
            delete _all_range->second;
        }
        delete _all_range;
        if (_vars) {
            if (!_embedded) {
                for (auto &elem: *_vars) {
                    delete elem.second.first;
                }
            }
            delete _vars;
        }
        if (_params) {
            if (!_embedded) {
                for (auto &elem: *_params) {
                    delete elem.second.first;
                }
            }
            delete _params;
        }
        if (_range) {
            for (auto &elem: *_range) {
                delete elem.first;
                delete elem.second;
            }
        }
        for (auto &f_p: _dfdx) {
            delete f_p.second;
        }    
        delete _range;
        delete _convexity;
        delete _sign;    
        delete _lterms;
        delete _qterms;
        delete _pterms;
        delete _expr;
        delete _DAG;
        delete _queue;
        delete _cst;
    };

    bool all_zeros(const string& s){
        auto it = s.begin();
        while (it != s.end()) {
            if ((*it)!='0' && (*it)!='.') {
                return false;
            }
            it++;
        }
        return true;
    }

    bool constant_::is_zero() const{ /**< Returns true if constant equals 0 */
//        auto a = poly_to_str(this);
        if (is_number() && all_zeros(poly_to_str(this))){
            return true;
        }
        if (is_param() || is_var()) {
            auto p_c = (param_*)this;
            return p_c->get_all_sign()==zero_;
        }
        if (is_uexpr() || is_bexpr()) {
            auto e_p = (expr*)this;
            return e_p->get_all_sign()==zero_;
        }
        if (is_function()) {
            return ((func_*)(this))->is_zero();
        }
        return false;
    }

    bool constant_::is_unit() const{ /**< Returns true if constant equals 1 */
        if(is_number() && poly_to_str(this)=="1"){
            return true;
        }
        if (is_param()) {
            auto p_c = (param_*)this;
            switch (p_c->get_intype()) {
                case binary_:
                    return ((param<bool>*)p_c)->is_unit();
                    break;
                case short_:
                    return ((param<short>*)p_c)->is_unit();
                    break;
                case integer_:
                    return ((param<int>*)p_c)->is_unit();
                    break;
                case float_:
                    return ((param<float>*)p_c)->is_unit();
                    break;
                case double_:
                    return ((param<double>*)p_c)->is_unit();
                    break;
                case long_:
                    return ((param<long double>*)p_c)->is_unit();
                    break;
                default:
                    break;
            }
        }
        return false;
    }


    bool constant_::is_positive() const{
        if (get_all_sign()==pos_) {
            return true;
        }
        return false;
    }


    bool constant_::is_non_positive() const{
        if (get_all_sign()==non_pos_) {
            return true;
        }
        return false;
    }

    bool constant_::is_non_negative() const{
        if (get_all_sign()==non_neg_) {
            return true;
        }
        return false;
    }

    bool constant_::is_negative() const{
        if (get_all_sign()==neg_) {
            return true;
        }
        return false;
    }

     void func_::update_sign(){
         
     }


    func_& func_::operator+=(const constant_& c){
        if (c.is_zero()) {
            return *this;
        }
        if (c.is_number() || (!is_constant() && c.is_param())) {
            _cst = add(_cst, c);
            update_sign(*_cst);
            return *this;
        }
        if (c.is_param() || c.is_var()) {
            this->insert(true, constant<double>(1), *(param_*)&c);
        }
        if (c.is_function()) {
            func_* f = (func_*)&c;
            if (!is_constant() && f->is_constant()) {//check _expr
                _cst = add(_cst, c);
                return *this;
            }
            if (!f->_cst->is_zero()) {
                _cst = add(_cst, *f->_cst);
            }
            for (auto &pair:*f->_lterms) {
                this->insert(pair.second);
            }
            for (auto &pair:*f->_qterms) {
                this->insert(pair.second);
            }
            for (auto &pair:*f->_pterms) {
                this->insert(pair.second);
            }
            
            if (_expr && f->_expr) {
                _expr = new bexpr(plus_, _expr, copy(*f->_expr));
                embed(*_expr);
//                _DAG->insert(make_pair<>(_expr->get_str(), _expr));
//                _queue->push(_expr);
            }
            if (!_expr && f->_expr) {
                _expr = (expr*)copy(*f->_expr);
                embed(*_expr);
//                _DAG->insert(make_pair<>(_expr->get_str(), _expr));
//                _queue->push(_expr);
                if (!_vars->empty()) {
                    _ftype = nlin_;
                }
            }
            update_sign(*f);
            update_convexity();
            return *this;
        }
        return *this;
    }

    func_& func_::operator-=(const constant_& c){
        if (c.is_zero()) {
            return *this;
        }
        if (c.is_number() || (!is_constant() && c.is_param())) {
            _cst = substract(_cst, c);
            update_sign(*_cst);
            return *this;
        }
        if (c.is_param() || c.is_var()) {
            this->insert(false, constant<double>(1), *(param_*)&c);
        }
        if (c.is_function()) {
            func_* f = (func_*)&c;
            if (!is_constant() && f->is_constant()) {
                _cst = substract(_cst, c);
                return *this;
            }
            _cst = substract(_cst, *f->_cst);
            for (auto &pair:*f->_lterms) {
                this->insert(!pair.second._sign, *pair.second._coef, *pair.second._p);
            }
            for (auto &pair:*f->_qterms) {
                this->insert(!pair.second._sign, *pair.second._coef, *pair.second._p->first, *pair.second._p->second);
            }
            for (auto &pair:*f->_pterms) {
                this->insert(!pair.second._sign, *pair.second._coef, *pair.second._l);
            }
            if (_expr && f->_expr) {
                _expr = new bexpr(minus_, _expr, copy(*f->_expr));
                embed(*_expr);
//                _DAG->insert(make_pair<>(_expr->get_str(), _expr));
//                _queue->push(_expr);
            }
            if (!_expr && f->_expr) {
                _expr = (expr*)copy(-1*(*f->_expr));
                embed(*_expr);
//                _DAG->insert(make_pair<>(_expr->get_str(), _expr));
//                _queue->push(_expr);
                if (!_vars->empty()) {
                    _ftype = nlin_;
                }
            }
            update_sign(*f);
            update_convexity();
            return *this;
        }
        return *this;
    }


    func_& func_::operator*=(const constant_& c){
        if (is_zero()) {
            return *this;
        }
    //    if (c.is_unit()) {
    //        return *this;
    //    }
        if (c.is_zero()) {
            reset();
            return *this;
        }
        /* Case where c is a number */
        if (c.is_number()){
            if (!_cst->is_zero()) {
                _cst = multiply(_cst, c);
            }
            for (auto &pair:*_lterms) {
                pair.second._coef = multiply(pair.second._coef, c);
            }
            for (auto &pair:*_qterms) {
                pair.second._coef = multiply(pair.second._coef, c);
            }
            for (auto &pair:*_pterms) {
                pair.second._coef = multiply(pair.second._coef, c);
            }
            if (c.is_negative()) {
                reverse_sign();
            }
            if (_expr) {
                auto be = new bexpr(product_, _expr, copy(c));
                _expr = be;
                embed(*_expr);
//                _DAG->insert(make_pair<>(_expr->get_str(), _expr));
//                _queue->push(_expr);
            }
            return *this;
        }
        if (is_nonlinear() || (c.is_function() && ((func_*)&c)->is_nonlinear())) {
            auto be = bexpr(product_, copy(*this), copy(c));
            *this = func_(be);
            return *this;
        }
        /* Case where the current function is not constant and the other operand is */
         if(!is_constant() && (c.is_param() || (c.is_function() && ((func_*)&c)->is_constant()))) {
             if (!_cst->is_zero()) {
                 _cst = multiply(_cst, c);
             }
            for (auto &pair:*_lterms) {
                pair.second._coef = multiply(pair.second._coef, c);
                if (c._is_transposed) {
                    pair.second._p->_is_vector = true;
                }
            }
            for (auto &pair:*_qterms) {
                pair.second._coef = multiply(pair.second._coef, c);
                if (c._is_transposed) {
                    pair.second._p->first->_is_vector = true;
                    pair.second._p->second->_is_vector = true;
                }

            }
            for (auto &pair:*_pterms) {
                pair.second._coef = multiply(pair.second._coef, c);
                if (c._is_transposed) {
                    for (auto& p: (*pair.second._l)) {
                        p.first->_is_vector = true;
                    }
                }

            }
            if (_expr) {
                auto be = new bexpr(product_, _expr, copy(c));
                _expr = be;
                embed(*_expr);
//                _DAG->insert(make_pair<>(_expr->get_str(), _expr));
//                _queue->push(_expr);
            }
            if (c.is_negative()) {
                reverse_sign();
            }
            if (c.get_all_sign()==unknown_) {
                _all_sign = unknown_;
                if (!_qterms->empty()) {
                    _all_convexity = undet_;
                }
            }
            return *this;
        }
        /* Case where the current function is constant and the other operand is not (we go to previous case) */
        if (is_constant() && (c.is_var() || (c.is_function() && !((func_*)&c)->is_constant()))) {
            func_ f(c);
            f *= *this;
            *this = move(f);
            return *this;
        }
        if (c.is_param() || c.is_var()) {
            func_ f(c);
            *this *= f;
            return *this;
        }
        /* Case where the multiplication invlolves multiplying variables/parameters together, i.e., they are both parametric or both include variables  */
        if (c.is_function()) {
            func_* f = (func_*)&c;
            constant_* coef;
            vector<bool>* is_sum = nullptr;
            func_ res;
            for (auto& t1: *_pterms) {
                if (t1.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial function), we cannot factor the coefficients. Just create a binary expression and return it.
                    bexpr e;
                    e += *this;
                    e *= c;
                    *this = e;
                    return *this;
                }
                for (auto& t2: *f->_pterms) {
                    is_sum = nullptr;
                    if (t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial function), see comment above.
                        bexpr e;
                        e += *this;
                        e *= c;
                        *this = e;
                        return *this;
                    }
                    auto newl(*t1.second._l);
                    for (auto& it: *t2.second._l) {
                        newl.push_back(make_pair<>(it.first, it.second));
                    }
                    coef = copy(*t1.second._coef);
                    coef = multiply(coef, *t2.second._coef);
                    res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
                    delete coef;
                }
                for (auto& t2: *f->_qterms) {
                    if (t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial function), see comment above.
                        bexpr e;
                        e += *this;
                        e *= c;
                        *this = e;
                        return *this;
                    }
                    auto newl(*t1.second._l);
                    newl.push_back(make_pair<>((t2.second._p->first), 1));
                    newl.push_back(make_pair<>((t2.second._p->second), 1));
                    coef = copy(*t1.second._coef);
                    coef = multiply(coef, *t2.second._coef);
                    res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
                    delete coef;
                }
                for (auto& t2: *f->_lterms) {
                    if (t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial function)
                        bexpr e;
                        e += *this;
                        e *= c;
                        *this = e;
                        return *this;
                    }
                    auto newl(*t1.second._l);
                    newl.push_back(make_pair<>((t2.second._p), 1));
                    coef = copy(*t1.second._coef);
                    coef = multiply(coef, *t2.second._coef);
                    res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
                    delete coef;
                }
                if (!f->_cst->is_zero()) {
                    auto newl(*t1.second._l);
                    coef = copy(*f->_cst);
                    coef = multiply(coef, *t1.second._coef);
                    res.insert(t1.second._sign, *coef, newl);
                    delete coef;
                }
            }

            for (auto& t1: *_qterms) {
                if (t1.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(Quadratic term)
                    bexpr e;
                    e += *this;
                    e *= c;
                    *this = e;
                    return *this;
                }
                for (auto& t2: *f->_pterms) {
                    is_sum = nullptr;
                    if (t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial term)
                        bexpr e;
                        e += *this;
                        e *= c;
                        *this = e;
                        return *this;
                    }
                    auto newl(*t2.second._l);
                    newl.push_front(make_pair<>(t1.second._p->first, 1));
                    newl.push_front(make_pair<>(t1.second._p->second, 1));
                    coef = copy(*t1.second._coef);
                    coef = multiply(coef, *t2.second._coef);
                    res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
                    delete coef;
                }
                for (auto& t2: *f->_qterms) {
                    if (t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial term)
                        bexpr e;
                        e += *this;
                        e *= c;
                        *this = e;
                        return *this;
                    }
                    coef = copy(*t1.second._coef);
                    coef = multiply(coef, *t2.second._coef);
                    list<pair<param_*, int>> newl;
                    newl.push_back(make_pair<>(t1.second._p->first, 1));
                    newl.push_back(make_pair<>(t1.second._p->second, 1));
                    newl.push_back(make_pair<>(t2.second._p->first, 1));
                    newl.push_back(make_pair<>(t2.second._p->second, 1));
                    res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
                    delete coef;
                }
                for (auto& t2: *f->_lterms) {
                    is_sum = nullptr;
                    if (t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial term)
                        bexpr e;
                        e += *this;
                        e *= c;
                        *this = e;
                        return *this;
                    }
                    coef = copy(*t1.second._coef);
                    coef = multiply(coef, *t2.second._coef);
                    list<pair<param_*, int>> newl;
                    newl.push_back(make_pair<>(t1.second._p->first, 1));
                    newl.push_back(make_pair<>(t1.second._p->second, 1));
                    newl.push_back(make_pair<>(t2.second._p, 1));
                    res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
                    delete coef;
                }
                if (!f->_cst->is_zero()) {
                    coef = copy(*t1.second._coef);
                    coef = multiply(coef, *f->_cst);
                    res.insert(t1.second._sign, *coef, *t1.second._p->first, *t1.second._p->second);
                    delete coef;
                }
                
            }
            for (auto& t1: *_lterms) {
                if (t1.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(Quadratic term)
                    bexpr e;
                    e += *this;
                    e *= c;
                    *this = e;
                    return *this;
                }
                for (auto& t2: *f->_pterms) {
                    is_sum = nullptr;
                    if (t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial term)
                        bexpr e;
                        e += *this;
                        e *= c;
                        *this = e;
                        return *this;
                    }
                    auto newl(*t2.second._l);
                    newl.push_front(make_pair<>((t1.second._p), 1));
                    coef = copy(*t1.second._coef);
                    coef = multiply(coef, *t2.second._coef);
                    res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
                    delete coef;
                }
                for (auto& t2: *f->_qterms) {
                    if (t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial term)
                        bexpr e;
                        e += *this;
                        e *= c;
                        *this = e;
                        return *this;
                    }
                    coef = copy(*t1.second._coef);
                    coef = multiply(coef, *t2.second._coef);
                    list<pair<param_*, int>> newl;
                    newl.push_back(make_pair<>(t1.second._p, 1));
                    newl.push_back(make_pair<>(t2.second._p->first, 1));
                    newl.push_back(make_pair<>(t2.second._p->second, 1));
                    res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
                    delete coef;
                }
                for (auto& t2: *f->_lterms) {
                    if (t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial term)
                        bexpr e;
                        e += *this;
                        e *= c;
                        *this = e;
                        return *this;
                    }
                    coef = copy(*t1.second._coef);
                    coef = multiply(coef, *t2.second._coef);
                    res.insert(!(t1.second._sign^t2.second._sign), *coef, *t1.second._p, *t2.second._p);
                    delete coef;
                }
                if (!f->_cst->is_zero()) {
                    coef = copy(*t1.second._coef);
                    coef = multiply(coef, *f->_cst);
                    res.insert(t1.second._sign, *coef, *t1.second._p);
                    delete coef;
                }
            }
            if (!_cst->is_zero()) {
                for (auto& t2: *f->_pterms) {
                    if (t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial term)
                        bexpr e;
                        e += *this;
                        e *= c;
                        *this = e;
                        return *this;
                    }
                    coef = copy(*_cst);
                    coef = multiply(coef, *t2.second._coef);
                    res.insert(t2.second._sign, *coef, *t2.second._l);
                    delete coef;
                }
                for (auto& t2: *f->_qterms) {
                    if (t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial term)
                        bexpr e;
                        e += *this;
                        e *= c;
                        *this = e;
                        return *this;
                    }
                    coef = copy(*_cst);
                    coef = multiply(coef, *t2.second._coef);
                    res.insert(t2.second._sign, *coef, *t2.second._p->first, *t2.second._p->second);
                    delete coef;
                }
                for (auto& t2: *f->_lterms) {
                    if (t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial term)
                        bexpr e;
                        e += *this;
                        e *= c;
                        *this = e;
                        return *this;
                    }
                    coef = copy(*_cst);
                    coef = multiply(coef, *t2.second._coef);
                    res.insert(t2.second._sign, *coef, *t2.second._p);
                    delete coef;
                }
                if (!f->_cst->is_zero()) {
                    coef = copy(*_cst);
                    coef = multiply(coef, *f->_cst);
                    delete res._cst;
                    res._cst = coef;
                }
            }
            
            *this = move(res);
        }
        return *this;
    }

    func_& func_::operator/=(const constant_& c){
        if (is_zero()) {
            return *this;
        }
        if (c.is_unit()) {
            return *this;
        }
        if (c.is_zero()) {
            throw invalid_argument("dividing by zero!\n");
        }
        /* Case where c is a number */
        if (c.is_number()){
            if (!_cst->is_zero()) {
                _cst = divide(_cst, c);
            }
            for (auto &pair:*_lterms) {
                pair.second._coef = divide(pair.second._coef, c);
            }
            for (auto &pair:*_qterms) {
                pair.second._coef = divide(pair.second._coef, c);
            }
            for (auto &pair:*_pterms) {
                pair.second._coef = divide(pair.second._coef, c);
            }
            if (_expr) {
                auto be = new bexpr(div_, _expr, copy(c));
                _expr = be;
                embed(*_expr);
//                _DAG->insert(make_pair<>(_expr->get_str(), _expr));
//                _queue->push(_expr);
            }
            if (c.is_negative()) {
                reverse_convexity();
                reverse_sign();
            }
        }
        /* Case where the current function is not constant and the other operand is */
        if(!is_constant() && (c.is_param() || (c.is_function() && ((func_*)&c)->is_constant()))) {
            if (!_cst->is_zero()) {
                _cst = divide(_cst, c);
            }
            for (auto &pair:*_lterms) {
                pair.second._coef = divide(pair.second._coef, c);
            }
            for (auto &pair:*_qterms) {
                pair.second._coef = divide(pair.second._coef, c);
            }
            for (auto &pair:*_pterms) {
                pair.second._coef = divide(pair.second._coef, c);
            }
            if (c.is_negative()) {
                reverse_convexity();
                reverse_sign();
            }
            if (c.get_all_sign()==unknown_) {
                _all_sign = unknown_;
                if (!_qterms->empty()) {
                    _all_convexity = undet_;
                }
            }
            if (_expr) {
                auto be = new bexpr(div_, _expr, copy(c));
                _expr = be;
                embed(*_expr);
//                _DAG->insert(make_pair<>(_expr->get_str(), _expr));
//                _queue->push(_expr);
            }
            return *this;
        }
    //    /* Case where the current function is constant and the other operand is not (we go to previous case) */
        if (is_constant() && (c.is_var() || (c.is_function() && !((func_*)&c)->is_constant()))) {
            func_ f(c);
            f /= *this;
            *this = move(f);
            return *this;
        }
        if (c.is_param() || c.is_var()) {
            func_ f(c);
            *this /= f;
            return *this;
        }
        /* Case where the multiplication invlolves multiplying variables/parameters together, i.e., they are both parametric or both include variables  */
        if (c.is_function()) {
            auto be = bexpr(div_, copy(*this), copy(c));
            *this = func_(be);
            return *this;
        }
        return *this;
    }

    //func_& func_::operator+=(const constant_& c){
    //    };
    //
    //func_& func_::operator-=(const constant_& c){
    //    return *this;
    //};
    //
    //func_& func_::operator*=(const constant_& c){
    //    return *this;
    //};

    //func_::func_(const func_& f){
    //    set_type(func_c);
    //    _ftype = f._ftype;
    //    _return_type = f._return_type;
    //    _convex = f._convex;
    //    _params = new map<string, pair<param_*, int>>();
    //    _vars = new map<string, pair<param_*, int>>();
    //    _lterms = new map<string, lterm>();
    //    _qterms = new map<string, qterm>();;
    //    _pterms = new map<string, pterm>();;
    //    param_* p_new;
    //    constant_* c_new;
    //    string str;
    //    for (auto &pair:*f._lterms) {
    //        p_new = (param_*)copy(pair.second._p);
    //        c_new = copy(pair.second._coef);
    //        str = p_new->get_name();
    //        _lterms->insert(make_pair<>(str, lterm(pair.second._sign, c_new, p_new)));
    //        if (p_new->is_var()) {
    //            _vars->insert(make_pair<>(str, make_pair<>(p_new, 1)));
    //        }
    //        else {
    //            _params->insert(make_pair<>(str, make_pair<>(p_new, 1)));
    //        }
    //    }
    //    for (auto &pair:*f._qterms) {
    //        insert(pair.second);
    //    }
    //    _cst = copy(f._cst);
    //}

    // TODO revisit embed and make it robust
    void func_::embed(expr& e){
        switch (e.get_type()) {
            case uexp_c:{
                auto ue = (uexpr*)&e;
                if (ue->_son->is_function()) {
                    auto f = (func_*)ue->_son;
                    embed(*f);
                }
                else if(ue->_son->is_expr()){
                    embed(*(expr*)ue->_son);
                }
                else if (ue->_son->is_var()){
                    if (_vars->count(((param_*)ue->_son)->get_name())==0) {
                        add_var((param_*)copy(*ue->_son));
                    }
                }
                break;
            }
            case bexp_c:{
                auto be = (bexpr*)&e;
                if (be->_lson->is_function()) {
                    auto f = (func_*)be->_lson;
                    embed(*f);
                }
                else if(be->_lson->is_expr()){
                    embed(*(expr*)be->_lson);
                }
                else if (be->_lson->is_var()){
                    if (_vars->count(((param_*)be->_lson)->get_name())==0) {
                        add_var((param_*)copy(*be->_lson));
                    }
                }
                if (be->_rson->is_function()) {
                    auto f = (func_*)be->_rson;
                    embed(*f);
                }
                else if(be->_rson->is_expr()){
                    embed(*(expr*)be->_rson);
                }
                else if (be->_rson->is_var()){
                    if (_vars->count(((param_*)be->_rson)->get_name())==0) {
                        add_var((param_*)copy(*be->_rson));
                    }
                }
                break;
            }
            default:
                break;
        }
    }

    void func_::embed(func_& f){
        f._embedded = true;
        param_* p = nullptr;
        param_* p1 = nullptr;
        param_* p2 = nullptr;
        for (auto &pair:*f._lterms) {
            p = pair.second._p;
            if (p->is_var()) {
                auto it = _vars->find(p->get_name());
                if (it==_vars->end()) {
                    add_var(p);
                }
                else{
                    p = it->second.first;
                    pair.second._p = p;
                    it->second.second++;
                }
            }
            else {
                auto it = _params->find(p->get_name());
                if (it==_params->end()) {
                    add_param(p);
                }
                else{
                    p = it->second.first;
                    pair.second._p = p;
                    it->second.second++;
                }
            }
        }
        for (auto &pair:*f._qterms) {
            p1 = pair.second._p->first;
            p2 = pair.second._p->second;
            if (p1->is_var()) {
                auto it1 = _vars->find(p1->get_name());
                if (it1==_vars->end()) {
                    add_var(p1);
                }
                else{
                    p1 = it1->second.first;
                    pair.second._p->first = p1;
                    it1->second.second++;
                }
                auto it2 = _vars->find(p2->get_name());
                if (it2==_vars->end()) {
                    add_var(p2);
                }
                else{
                    p2 = it2->second.first;
                    pair.second._p->second = p2;
                    it2->second.second++;
                }
            }
            else {
                auto it1 = _params->find(p1->get_name());
                if (it1==_params->end()) {
                    add_param(p1);
                }
                else{
                    p1 = it1->second.first;
                    pair.second._p->first = p1;
                    it1->second.second++;
                }
                auto it2 = _params->find(p2->get_name());
                if (it2==_params->end()) {
                    add_param(p2);
                }
                else{
                    p2 = it2->second.first;
                    pair.second._p->second = p2;
                    it2->second.second++;
                }
            }
        }
        for (auto &pair:*f._pterms) {
            auto list = pair.second._l;
            for (auto &ppi: *list) {
                p = ppi.first;
                if (p->is_var()) {
                    auto it = _vars->find(p->get_name());
                    if (it==_vars->end()) {
                        add_var(p);
                    }
                    else{
                        p = it->second.first;
                        ppi.first = p;
                        it->second.second++;
                    }
                }
                else {
                    auto it = _params->find(p->get_name());
                    if (it==_params->end()) {
                        add_param(p);
                    }
                    else{
                        p = it->second.first;
                        ppi.first = p;
                        it->second.second++;
                    }
                }
            }
        }
        if (f._expr) {
            embed(*f._expr);
        }
        auto old_vars = *f._vars;
        for (auto &vp: old_vars) {
            auto vv = (*_vars)[vp.first].first;
            if (vv != vp.second.first) {
                f._vars->erase(vp.first);
                f._vars->insert(make_pair<>(vp.first, make_pair<>(vv, 1)));
            }
        }
        auto old_params = *f._params;
        for (auto &pp: old_params) {
            auto p = (*_params)[pp.first].first;
            if (p != pp.second.first) {            
                f._params->erase(pp.first);
                f._params->insert(make_pair<>(pp.first, make_pair<>(p, 1)));
            }
        }
    }

    void func_::reset(){
        if (_all_range) {
            delete _all_range->first;
            delete _all_range->second;
        }
        delete _all_range;
        _all_range = new pair<constant_*, constant_*>(new constant<double>(0), new constant<double>(0));
        if (_vars) {
            if (!_embedded) {
                for (auto &elem: *_vars) {
                    delete elem.second.first;
                }
            }
            _vars->clear();
        }
        if (_params) {
            if (!_embedded) {
                for (auto &elem: *_params) {
                    delete elem.second.first;
                }
            }
            _params->clear();
        }
        if (_range) {
            for (auto &elem: *_range) {
                delete elem.first;
                delete elem.second;
            }
        }
        for (auto &df:_dfdx) {
            delete df.second;
        }
        _dfdx.clear();
        _hess_link.clear();
        delete _range;
        _range = nullptr;
        delete _convexity;
        _convexity = nullptr;
        delete _sign;
        _sign = nullptr;
        delete _expr;
        _expr = nullptr;
        delete _DAG;
        _DAG = nullptr;
        delete _queue;
        _queue = nullptr;
        set_type(func_c);
        _ftype = const_;
        _return_type = integer_;
        _all_convexity = linear_;
        _all_sign = zero_;
        _is_transposed = false;
        _is_vector = false;
        _lterms->clear();
        _qterms->clear();
        _pterms->clear();
        delete _cst;
        _cst = new constant<double>(0);
        _nb_instances = 1;
        _nnz_h = 0;
        _nnz_j = 0;
    };

    //void func_::reverse_sign(){ /*<< Reverse the sign of all terms in the function */
    //    for (auto &pair: *_lterms) {
    //        pair.second.reverse_sign();
    //    }
    //    for (auto &pair: *_qterms) {
    //        pair.second.reverse_sign();
    //    }
    //    ::reverse_sign(_cst);
    //    reverse_convexity();
    //    if (_sign==neg_) {
    //        _sign=pos_;
    //    }
    //    else if (_sign==pos_) {
    //        _sign=neg_;
    //    }
    //    else if(_sign==non_neg_) {
    //        _sign=non_pos_;
    //    }
    //    else if(_sign==non_pos_) {
    //        _sign=non_neg_;
    //    }
    //}

    void func_::reverse_convexity(){
        if (_all_convexity==convex_) {
            _all_convexity=concave_;
        }
        else if (_all_convexity==concave_) {
            _all_convexity=convex_;
        }
    }

    void func_::reverse_sign(){
        reverse_convexity();
        if (_all_sign==pos_) {
            _all_sign = neg_;
        }
        else if (_all_sign==neg_) {
            _all_sign = pos_;
        }
        else if (_all_sign== non_neg_){
            _all_sign = non_pos_;
        }
        else if (_all_sign== non_pos_){
            _all_sign = non_neg_;
        }
    }

    void func_::update_sign(const constant_& c){
        Sign sign = c.get_all_sign();
        if (sign==unknown_ || ((_all_sign==non_neg_ || _all_sign==pos_) && sign!=non_neg_ && sign!=pos_)) {
            _all_sign = unknown_;
        }
        else if((_all_sign==non_pos_ || _all_sign==neg_) && sign!=non_pos_ && sign!=neg_){
            _all_sign = unknown_;
        }
        else if(_all_sign==zero_ || _all_sign==pos_ || _all_sign==neg_){// take weaker sign
            _all_sign = sign;
        }
    }

    void func_::update_sign(const lterm& l){
        Sign sign = get_all_sign(l);
        if (sign==unknown_ || ((_all_sign==non_neg_ || _all_sign==pos_) && sign!=non_neg_ && sign!=pos_)) {
            _all_sign = unknown_;
        }
        else if((_all_sign==non_pos_ || _all_sign==neg_) && sign!=non_pos_ && sign!=neg_){
            _all_sign = unknown_;
        }
        else if(_all_sign==zero_ || _all_sign==pos_ || _all_sign==neg_){// take weaker sign
            _all_sign = sign;
        }
    }

    void func_::update_sign(const qterm& q){
        Sign sign = get_all_sign(q);
        if (sign==unknown_ || ((_all_sign==non_neg_ || _all_sign==pos_) && sign!=non_neg_ && sign!=pos_)) {
            _all_sign = unknown_;
        }
        else if((_all_sign==non_pos_ || _all_sign==neg_) && sign!=non_pos_ && sign!=neg_){
            _all_sign = unknown_;
        }
        else if(_all_sign==zero_ || _all_sign==pos_ || _all_sign==neg_){// take weaker sign
            _all_sign = sign;
        }
    }

    void func_::update_sign(const pterm& p){
        Sign sign = get_all_sign(p);
        if (sign==unknown_ || ((_all_sign==non_neg_ || _all_sign==pos_) && sign!=non_neg_ && sign!=pos_)) {
            _all_sign = unknown_;
        }
        else if((_all_sign==non_pos_ || _all_sign==neg_) && sign!=non_pos_ && sign!=neg_){
            _all_sign = unknown_;
        }
        else if(_all_sign==zero_ || _all_sign==pos_ || _all_sign==neg_){// take weaker sign
            _all_sign = sign;
        }
    }

    void func_::update_convexity(const qterm& q){
        Convexity conv = get_convexity(q);
        if (_all_convexity==undet_ || conv ==undet_ || (_all_convexity==convex_ && conv==concave_) || (_all_convexity==concave_ && conv==convex_)) {
            _all_convexity = undet_;
        }
        else {
            _all_convexity = conv;
        }
    }

    bool func_::insert(bool sign, const constant_& coef, const param_& p){/**< Adds coef*p to the linear function. Returns true if added new term, false if only updated coef of p */
        param_* p_new;
        auto pname = p.get_name();
        
        auto pair_it = _lterms->find(pname);
        if (_ftype == const_ && p.is_var()) {
            _ftype = lin_;
        }
        
        if (pair_it == _lterms->end()) {
            auto c_new = copy(coef);
            if (p.is_var()) {
                p_new = get_var(pname);
                if (!p_new) {
                    p_new = (param_*)copy(p);
                    add_var(p_new);
                }
                else {
                    incr_occ_var(pname);
                }
            }
            else {
                p_new = get_param(pname);
                if (!p_new) {
                    p_new = (param_*)copy(p);
                    add_param(p_new);
                }
                else {
                    incr_occ_param(pname);
                }
            }
            lterm l(sign, c_new, p_new);
            update_sign(l);
            _lterms->insert(make_pair<>(pname, move(l)));
            return true;
        }
        else {
            if (pair_it->second._sign == sign) {
                pair_it->second._coef = add(pair_it->second._coef, coef);
            }
            else{
                pair_it->second._coef = substract(pair_it->second._coef, coef);
            }
            if (pair_it->second._coef->is_zero()) {            
                if (p.is_var()) {
                    decr_occ_var(pname);
                }
                else{
                    decr_occ_param(pname);
                }
                _lterms->erase(pair_it);
                //update_sign();
            }
            else {
                update_sign(pair_it->second);
            }
            return false;
        }
    };

    void func_::insert(const lterm& term){
        insert(term._sign, *term._coef, *term._p);
    }

    bool func_::insert(bool sign, const constant_& coef, const param_& p1, const param_& p2){/**< Adds coef*p1*p2 to the function. Returns true if added new term, false if only updated coef of p1*p2 */
        auto ps1 = p1.get_name();
        auto ps2 = p2.get_name();
        auto qname = ps1+","+ps2;
        auto pair_it = _qterms->find(qname);
        param_* p_new1;
        param_* p_new2;
        
        if (_ftype <= lin_ && p1.is_var()) {
            _ftype = quad_;
        }
        
        if (pair_it == _qterms->end()) {
            if (p1.is_var()) {
                p_new1 = (param_*)get_var(ps1);
                if (!p_new1) {
                    p_new1 = (param_*)copy(p1);
                    add_var(p_new1);
                }
                else {
                    incr_occ_var(ps1);
                }
            }
            else {
                p_new1 = (param_*)get_param(ps1);
                if (!p_new1) {
                    p_new1 = (param_*)copy(p1);
                    add_param(p_new1);
                }
                else {
                    incr_occ_param(ps1);
                }
                
            }
            if (p2.is_var()) {
                p_new2 = get_var(ps2);
                if (!p_new2) {
                    p_new2 = (param_*)copy(p2);
                    add_var(p_new2);
                }
                else {
                    incr_occ_var(ps2);
                }
            }
            else {
                p_new2 = get_param(ps2);
                if (!p_new2) {
                    p_new2 = (param_*)copy(p2);
                    add_param(p_new2);
                }
                else {
                    incr_occ_param(ps2);
                }
            }
            auto c_new = copy(coef);
            qterm q(sign, c_new, p_new1, p_new2);
            update_sign(q);
            update_convexity(q);
            _qterms->insert(make_pair<>(qname, move(q)));
            return true;
        }
        else {
            if (pair_it->second._sign == sign) {
                pair_it->second._coef = add(pair_it->second._coef, coef);
            }
            else{
                pair_it->second._coef = substract(pair_it->second._coef, coef);
            }
            if (pair_it->second._coef->is_zero()) {
                if (p1.is_var()) {
                    decr_occ_var(ps1);
                }
                else {
                    decr_occ_param(ps1);
                }
                if (p2.is_var()) {
                    decr_occ_var(ps2);
                }
                else {
                    decr_occ_param(ps2);
                }
                _qterms->erase(pair_it);
                update_sign();
                update_convexity();
            }
            else {
                update_sign(pair_it->second);
                update_convexity(pair_it->second);
            }
            return false;
        }
    };

    void func_::insert(const qterm& term){
        insert(term._sign, *term._coef, *term._p->first, *term._p->second);
    }


    bool func_::insert(bool sign, const constant_& coef, const list<pair<param_*, int>>& l){/**< Adds polynomial term to the function. Returns true if added new term, false if only updated corresponding coef */
        _all_convexity = undet_;
        string name;
        string s;
        bool newv = true;
//        int i = 0;
        for (auto &pair:l) {
            name += pair.first->get_name();
            name += "^"+to_string(pair.second);
            name += ",";
        }
        auto pair_it = _pterms->find(name);
        param_* p = l.begin()->first;
        param_* pnew;
        if (_ftype <= quad_ && p->is_var()) {
            _ftype = pol_;
        }    
        if (pair_it == _pterms->end()) {
            auto newl = new list<pair<param_*, int>>();
//            i = 1;
            for (auto &pair:l) {
                p = pair.first;
                s = p->get_name();
                if (p->is_var()) {
                    pnew = (param_*)get_var(s);
                    if (!pnew) {
                        pnew = (param_*)copy(*p);
                        add_var(pnew,pair.second);
                    }
                    else {
                        incr_occ_var(s);
                    }
                }
                else {
                    pnew = (param_*)get_param(s);
                    if (!pnew) {
                        pnew = (param_*)copy(*p);
                        add_param(pnew);
                    }
                    else {
                        incr_occ_param(s);
                    }
                }
                newv = true;
                for (auto& p_it:*newl) {
                    if (p_it.first->get_name()==s) {
                        p_it.second++;
                        newv = false;
                        break;
                    }
                }
                if (newv) {
                    newl->push_back(make_pair<>(pnew, pair.second));
                }
            }
            auto c_new = copy(coef);
            pterm p(sign, c_new, newl);
            update_sign(p);
            _pterms->insert(make_pair<>(name, move(p)));
            return true;
        }
        else {
            if (pair_it->second._sign == sign) {
                pair_it->second._coef = add(pair_it->second._coef, coef);
            }
            else{
                pair_it->second._coef = substract(pair_it->second._coef, coef);
            }
            
            if (pair_it->second._coef->is_zero()) {
                for (auto& it:*pair_it->second._l) {
                    p = it.first;
                    s = p->get_name();
                    if (p->is_var()) {
                        decr_occ_var(s,it.second);
                    }
                    else {
                        decr_occ_param(s,it.second);
                    }
                }
                _pterms->erase(pair_it);
                update_sign();
                update_convexity();
            }
            else {
                update_sign(pair_it->second);            
            }
            return false;
        }

    }

    void func_::insert(const pterm& term){
        insert(term._sign, *term._coef, *term._l);
    }

    void func_::insert(expr& e){
    //    insert(term._sign, *term._coef, *term._l);
        auto name = e.get_str();
        auto pair_it = _DAG->find(name);
        
        _ftype = nlin_;
        
    //    if (pair_it == _DAG->end()) {
    //        update_sign(e);
    //        update_convexity(e);
    //        _DAG->insert(make_pair<>(name, e));
    //    }
    //    else {
    //        if (pair_it->second._sign == e.sign) {
    ////            pair_it->second._coef = add(pair_it->second._coef, coef);
    ////            if (!pair_it->second._sign) { // both negative
    ////                pair_it->second._sign = true;
    ////            }
    //        }
    //        else{
    ////            pair_it->second._coef = substract(pair_it->second._coef, coef);
    //        }
    ////        if (pair_it->second._coef->is_zero()) {
    ////            if (p1.is_var()) {
    ////                decr_occ_var(s1);
    ////            }
    ////            else {
    ////                decr_occ_param(s1);
    ////            }
    ////            if (p2.is_var()) {
    ////                decr_occ_var(s2);
    ////            }
    ////            else {
    ////                decr_occ_param(s2);
    ////            }
    ////            _qterms->erase(pair_it);
    ////            update_sign();
    ////            update_convexity();
    ////        }
    ////        else {
    ////            update_sign(pair_it->second);
    ////            update_convexity(pair_it->second);
    ////        }
    ////        return false;
    //    }
    ////    _DAG->insert(<#const value_type &__v#>)
    }


    string expr::get_str(){
        if (_to_str=="" || _to_str=="null") {
            if (is_uexpr()) {
                return _to_str = ((uexpr*)this)->to_str();
            }
            else {
                return _to_str = ((bexpr*)this)->to_str();
            }
        }
        else {
            return _to_str;
        }
    }

    uexpr::uexpr(const uexpr& exp){
        _otype = exp._otype;
        _son = copy(*exp._son);
        _to_str = exp._to_str;
        _type = uexp_c;
    }

    uexpr::uexpr(uexpr&& exp){
        _otype = exp._otype;
        _son = move(exp._son);
        exp._son = nullptr;
        _to_str = exp._to_str;
        _type = uexp_c;
    }

    uexpr& uexpr::operator=(const uexpr& exp){
        delete _son;
        _son = copy(*exp._son);
        _otype = exp._otype;
        _to_str = exp._to_str;
        return *this;
    }

    uexpr& uexpr::operator=(uexpr&& exp){
        delete _son;
        _son = move(exp._son);
        exp._son = nullptr;
        _otype = exp._otype;
        _to_str = exp._to_str;
        return *this;
    }


    bool uexpr::contains(const constant_* c) const{
        if (!_son) {
            return false;
        }
        if (equals(_son, c)) {
            return true;
        }
        
        if (_son->get_type()==uexp_c) {
            return ((uexpr*)_son)->contains(c);
        }
        if (_son->get_type()==bexp_c) {
            return ((bexpr*)_son)->contains(c);
        }
        return false;
    };


    bool uexpr::operator==(const uexpr &c)const{
        return (_otype == c._otype && equals(_son,c._son));
        
    }



    func_ cos(const constant_& c){
        uexpr e;
        e._otype = cos_;
        if(!c.is_function()){
            e._son = new func_(c);
        }
        else {
            e._son = copy(c);
        }
        func_ res;
        res._expr = new uexpr(move(e));
        res.embed(*res._expr);
//        res._DAG->insert(make_pair<>(res._expr->get_str(), res._expr));
//        res._queue->push(res._expr);
        if (!res._vars->empty()) {
            res._ftype = nlin_;
        }
        return res;
    };

    func_ cos(constant_&& c){
        uexpr e;
        e._otype = cos_;
        if(!c.is_function()){
            e._son = new func_(move(c));
        }
        else {
            e._son = copy(move(c));
        }
        func_ res;
        res._expr = new uexpr(move(e));
        res.embed(*res._expr);
//        res._DAG->insert(make_pair<>(res._expr->get_str(), res._expr));
//        res._queue->push(res._expr);
        if (!res._vars->empty()) {
            res._ftype = nlin_;
        }
        return res;
    };


    func_ sin(const constant_& c){
        uexpr e;
        e._otype = sin_;
        if(!c.is_function()){
            e._son = new func_(c);
        }
        else {
            e._son = copy(c);
        }
        func_ res;
        res._expr = new uexpr(move(e));
        res.embed(*res._expr);
//        res._DAG->insert(make_pair<>(res._expr->get_str(), res._expr));
//        res._queue->push(res._expr);
        if (!res._vars->empty()) {
            res._ftype = nlin_;
        }
        return res;
    };

    func_ sin(constant_&& c){
        uexpr e;
        e._otype = sin_;
        if(!c.is_function()){
            e._son = new func_(move(c));
        }
        else {
            e._son = copy(move(c));
        }
        func_ res;
        res._expr = new uexpr(move(e));
        res.embed(*res._expr);
//        res._DAG->insert(make_pair<>(res._expr->get_str(), res._expr));
//        res._queue->push(res._expr);
        if (!res._vars->empty()) {
            res._ftype = nlin_;
        }
        return res;
    };


    func_ sqrt(const constant_& c){
        uexpr e;
        e._otype = sqrt_;
        if(!c.is_function()){
            e._son = new func_(c);
        }
        else {
            e._son = copy(c);
        }
        func_ res;
        res._expr = new uexpr(move(e));
        res.embed(*res._expr);
//        res._DAG->insert(make_pair<>(res._expr->get_str(), res._expr));
//        res._queue->push(res._expr);
        if (!res._vars->empty()) {
            res._ftype = nlin_;
        }
        return res;
    };

    func_ sqrt(constant_&& c){
        uexpr e;
        e._otype = sqrt_;
        if(!c.is_function()){
            e._son = new func_(move(c));
        }
        else {
            e._son = copy(move(c));
        }
        func_ res;
        res._expr = new uexpr(move(e));
        res.embed(*res._expr);
//        res._DAG->insert(make_pair<>(res._expr->get_str(), res._expr));
//        res._queue->push(res._expr);
        if (!res._vars->empty()) {
            res._ftype = nlin_;
        }
        return res;
    };


    func_ expo(const constant_& c){
        uexpr e;
        e._otype = exp_;
        if(!c.is_function()){
            e._son = new func_(c);
        }
        else {
            e._son = copy(c);
        }
        func_ res;
        res._expr = new uexpr(move(e));
        res.embed(*res._expr);
//        res._DAG->insert(make_pair<>(res._expr->get_str(), res._expr));
//        res._queue->push(res._expr);
        if (!res._vars->empty()) {
            res._ftype = nlin_;
        }
        return res;
    };


    func_ expo(constant_&& c){
        uexpr e;
        e._otype = exp_;
        if(!c.is_function()){
            e._son = new func_(move(c));
        }
        else {
            e._son = copy(move(c));
        }
        func_ res;
        res._expr = new uexpr(move(e));
        res.embed(*res._expr);
//        res._DAG->insert(make_pair<>(res._expr->get_str(), res._expr));
//        res._queue->push(res._expr);
        if (!res._vars->empty()) {
            res._ftype = nlin_;
        }
        return res;
    };



    func_ log(const constant_& c){
        uexpr e;
        e._otype = log_;
        if(!c.is_function()){
            e._son = new func_(c);
        }
        else {
            e._son = copy(c);
        }
        func_ res;
        res._expr = new uexpr(move(e));
        res.embed(*res._expr);
//        res._DAG->insert(make_pair<>(res._expr->get_str(), res._expr));
//        res._queue->push(res._expr);
        if (!res._vars->empty()) {
            res._ftype = nlin_;
        }
        return res;
    };

    func_ log(constant_&& c){
        uexpr e;
        e._otype = log_;
        if(!c.is_function()){
            e._son = new func_(move(c));
        }
        else {
            e._son = copy(move(c));
        }
        func_ res;
        res._expr = new uexpr(move(e));
        res.embed(*res._expr);
//        res._DAG->insert(make_pair<>(res._expr->get_str(), res._expr));
//        res._queue->push(res._expr);
        if (!res._vars->empty()) {
            res._ftype = nlin_;
        }

        return res;
    };

    double uexpr::eval(size_t i) const{
        if (!_son) {
            throw invalid_argument("Cannot evaluate empty expression!");
        }
        switch (_otype) {
            case cos_:
                return std::cos(poly_eval(_son,i));
                break;
            case sin_:
                return std::sin(poly_eval(_son,i));
                break;
            case sqrt_:
                return std::sqrt(poly_eval(_son,i));
                break;
            case log_:
                return std::log(poly_eval(_son,i));
                break;
            case exp_:
                return exp(poly_eval(_son,i));
                break;
            default:
                throw invalid_argument("Unsupported unary operator");
                break;
        }
        
    }


    /* BINARY EXPRESSIONS */

    bool bexpr::operator==(const bexpr &c)const{
        return (_otype == c._otype && equals(_lson,c._lson) && equals(_rson,c._rson));
    }


    bool bexpr::contains(const constant_* c) const{
        if(_lson){
            if (equals(_lson, c)) {
                return true;
            }
            if (_lson->get_type()==uexp_c && ((uexpr*)_lson)->contains(c)) {
                return true;
            }
            if (_lson->get_type()==bexp_c && ((bexpr*)_lson)->contains(c)) {
                return true;
            }
        }
        if(_rson) {
            if (equals(_rson, c)) {
                return true;
            }
            if (_rson->get_type()==uexp_c && ((uexpr*)_rson)->contains(c)) {
                return true;
            }
            if (_rson->get_type()==bexp_c && ((bexpr*)_rson)->contains(c)) {
                return true;
            }
        }
        return false;
    };

    bexpr::bexpr(){
        _otype = id_;
        _lson = nullptr;
        _rson = nullptr;
        _type = bexp_c;
        _to_str = "null";
    }

    bexpr::bexpr(OperatorType otype, constant_* lson, constant_* rson){
        _otype = otype;
        _lson = lson;
        _rson = rson;
        _type = bexp_c;
        _to_str = poly_to_str(_lson);
        switch (otype) {
            case plus_:
                _to_str += " + ";
                break;
            case minus_:
                _to_str += " - ";
                break;
            case product_:
                _to_str += " * ";
                break;
            case div_:
                _to_str += " / ";
                break;
            default:
                break;
        }
        _to_str += poly_to_str(_rson);
    };

    bexpr::bexpr(const bexpr& exp){ /**< Copy constructor from binary expression tree */
        _otype = exp._otype;
        _lson = copy(*exp._lson);
        _rson =  copy(*exp._rson);
        _type = bexp_c;
        _to_str = exp._to_str;
    };

    bexpr::bexpr(bexpr&& exp){ /**< Move constructor from binary expression tree */
        _otype = exp._otype;
        _lson = move(exp._lson);
        _rson = move(exp._rson);
        _type = bexp_c;
        _to_str = exp._to_str;
    };

    bexpr& bexpr::operator=(const bexpr& exp){
        delete _lson;
        delete _rson;
        _lson = copy(*exp._lson);
        _rson =  copy(*exp._rson);
        _otype = exp._otype;
        _to_str = exp._to_str;
        return *this;
    }

    bexpr& bexpr::operator=(bexpr&& exp){
        delete _lson;
        delete _rson;
        _lson = move(exp._lson);
        _rson = move(exp._rson);
        _otype = exp._otype;
        _to_str = exp._to_str;
        return *this;
    }



    template<typename other_type> bexpr& bexpr::operator+=(const other_type& v){
        bexpr res(plus_, copy(*this), copy(v));
        delete _lson;
        delete _rson;
        return *this=res;
    }


    template<typename other_type> bexpr& bexpr::operator-=(const other_type& v){
        bexpr res(minus_, copy(*this), copy(v));
        delete _lson;
        delete _rson;
        return *this=res;
    }


    template<typename other_type> bexpr& bexpr::operator*=(const other_type& v){
        bexpr res(product_, copy(*this), copy(v));
        delete _lson;
        delete _rson;
        return *this=res;
    }

    template<typename other_type> bexpr& bexpr::operator/=(const other_type& v){
        bexpr res(div_, copy(*this), copy(v));
        delete _lson;
        delete _rson;
        return *this=res;
    }



    /* Polymorphic functions */

    //void reverse_sign(constant_* c){ /**< Reverses the sign of the constant. */
    //    switch (c->get_type()) {
    //        case binary_c: {
    //            ((constant<bool>*)c)->set_val(!((constant<bool>*)c)->eval());
    //            break;
    //        }
    //        case short_c: {
    //            ((constant<short>*)c)->set_val(-1*((constant<short>*)c)->eval());
    //            break;
    //        }
    //        case integer_c: {
    //            ((constant<int>*)c)->set_val(-1*((constant<int>*)c)->eval());
    //            break;
    //        }
    //        case float_c: {
    //            ((constant<float>*)c)->set_val(-1*((constant<float>*)c)->eval());
    //            break;
    //        }
    //        case double_c: {
    //            ((constant<double>*)c)->set_val(-1*((constant<double>*)c)->eval());
    //            break;
    //        }
    //        case long_c: {
    //            ((constant<long double>*)c)->set_val(-1*((constant<long double>*)c)->eval());
    //            break;
    //        }
    //        default:
    //            throw invalid_argument("Cannot reverse sign of non-numeric constant");
    //            break;
    //    }
    //    
    //}

    void poly_print(const constant_* c){/**< printing c, detecting the right class, i.e., constant<>, param<>, uexpr or bexpr. */
        cout << poly_to_str(c);
    }


    string poly_to_str(const constant_* c){/**< printing c, detecting the right class, i.e., constant<>, param<>, uexpr or bexpr. */
        
        if (!c) {
            return "null";
        }    
        switch (c->get_type()) {
            case binary_c: {
                return ((constant<bool>*)(c))->to_str();
                break;
            }
            case short_c: {
                return ((constant<short>*)(c))->to_str();
                break;
            }
            case integer_c: {
                return ((constant<int>*)(c))->to_str();
                break;
            }
            case float_c: {
                return ((constant<float>*)(c))->to_str();
                break;
            }
            case double_c: {
                return ((constant<double>*)(c))->to_str();
                break;
            }
            case long_c: {
                return ((constant<long double>*)(c))->to_str();
                break;
            }
            case par_c:{
                auto p_c = (param_*)(c);
                switch (p_c->get_intype()) {
                    case binary_:
                        return ((param<bool>*)p_c)->get_name();
                        break;
                    case short_:
                        return ((param<short>*)p_c)->get_name();
                        break;
                    case integer_:
                        return ((param<int>*)p_c)->get_name();
                        break;
                    case float_:
                        return ((param<float>*)p_c)->get_name();
                        break;
                    case double_:
                        return ((param<double>*)p_c)->get_name();
                        break;
                    case long_:
                        return ((param<long double>*)p_c)->get_name();
                        break;
                    default:
                        break;
                }
                break;
            }
            case uexp_c: {
                return ((uexpr*)c)->get_str();
                break;
            }
            case bexp_c: {
                return ((bexpr*)c)->get_str();
                break;
            }
            case var_c: {
                auto p_c = (param_*)(c);
                switch (p_c->get_intype()) {
                    case binary_:
                        return ((var<bool>*)p_c)->get_name();
                        break;
                    case short_:
                        return ((var<short>*)p_c)->get_name();
                        break;
                    case integer_:
                        return ((var<int>*)p_c)->get_name();
                        break;
                    case float_:
                        return ((var<float>*)p_c)->get_name();
                        break;
                    case double_:
                        return ((var<double>*)p_c)->get_name();
                        break;
                    case long_:
                        return ((var<long double>*)p_c)->get_name();
                        break;
                    default:
                        break;
                }
                break;
            }

            case sdpvar_c: {
                auto p_c = (param_*)(c);
                switch (p_c->get_intype()) {
                    case binary_:
                        return ((sdpvar<bool>*)p_c)->get_name();
                        break;
                    case short_:
                        return ((sdpvar<short>*)p_c)->get_name();
                        break;
                    case integer_:
                        return ((sdpvar<int>*)p_c)->get_name();
                        break;
                    case float_:
                        return ((sdpvar<float>*)p_c)->get_name();
                        break;
                    case double_:
                        return ((sdpvar<double>*)p_c)->get_name();
                        break;
                    case long_:
                        return ((sdpvar<long double>*)p_c)->get_name();
                        break;
                    default:
                        break;
                }
                break;
            }
            case func_c: {
                return ((func_*)c)->to_str();
                break;
            }
            default:
                break;
        }
        return "null";
    }

    bool equals(const constant_* c1, const constant_* c2){/**< Checks if c2 equals c1 detecting the right class, i.e., constant<>, param<>, uexpr or bexpr. */
        if ((!c1 && !c2) || (c1==c2)) {
            return true;
        }
        if ((c1 && !c2) || (!c1 && c2)) {
            return false;
        }
        
        switch (c2->get_type()) {
            case binary_c: {
                return (c1->is_binary() && *(constant<bool>*)c1 == *(constant<bool>*)c2);
                break;
            }
            case short_c: {
                return (c1->is_short() && *(constant<short>*)c1 == *(constant<short>*)c2);
                break;
            }
            case integer_c: {
                return (c1->is_integer() && *(constant<int>*)c1 == *(constant<int>*)c2);
                break;
            }
            case float_c: {
                return (c1->is_float() && *(constant<float>*)c1 == *(constant<float>*)c2);
                break;
            }
            case double_c: {
                return (c1->is_double() && *(constant<double>*)c1 == *(constant<double>*)c2);
                break;
            }
            case long_c: {
                return (c1->is_long() && *(constant<long double>*)c1 == *(constant<long double>*)c2);
                break;
            }
            case par_c:{
                return (c1->is_param() && *(param_ *)c1 == *(param_ *)c2);
                break;
            }
            case uexp_c: {
                return (c1->is_uexpr() && *(uexpr *)c1 == *(uexpr *)c2);
                break;
            }
            case bexp_c: {
                return (c1->is_bexpr() && *(bexpr *)c1 == *(bexpr *)c2);
                break;
            }
            case var_c:{
                return (c1->is_var() && *(param_ *)c1 == *(param_ *)c2);
                break;
            }
            case sdpvar_c:{
                return (c1->is_sdpvar() && *(param_ *)c1 == *(param_ *)c2);
                break;
            }
            case func_c:{
                if (c1->is_function()){
                    
                    auto f1 = (func_ *)c1;
                    auto f2 = (func_ *)c2;
                    return f1==f2;
                }
                break;
            }
            default:
                break;
        }
        return false;
    }



//    template<typename type> type  eval(const constant_* c1){
//        return eval<type>(0,c1);
//    };

    double  bexpr::eval(ind i) const{
        if (!_lson || !_rson) {
            throw invalid_argument("Cannot evaluate empty expression!");
        }
        switch (_otype) {
            case plus_:
                return poly_eval(_lson,i) + poly_eval(_rson,i);
                break;
            case minus_:
                return poly_eval(_lson,i) - poly_eval(_rson,i);
                break;
            case product_:
                return poly_eval(_lson,i) * poly_eval(_rson,i);
                break;
            case div_:
                return poly_eval(_lson,i)/poly_eval(_rson,i);
                break;
            case power_:
                return powl(poly_eval(_lson,i),poly_eval(_rson,i));
                break;
            default:
                throw invalid_argument("Unsupported binary operator");
                break;
        }
        
    }

    func_ operator+(const constant_& c1, const constant_& c2){
        return func_(c1) += c2;
    }

//    func_ operator+(func_&& f, const constant_& c){
//        return f += c;
//    }
//
//    func_ operator+(const constant_& c, func_&& f){
//        return f += c;
//    }


    func_ operator-(const constant_& c1, const constant_& c2){
        return func_(c1) -= c2;
    }

//    func_ operator-(func_&& f, const constant_& c){
//        return f -= c;
//    }
//
//    func_ operator-(const constant_& c, func_&& f){
//        return (f *= -1) += c;
//    }
//

    func_ operator*(const constant_& c1, const constant_& c2){// Rewrite this to change res after the multiplication is done, make sure both vars are now vecs.
        if(c1.is_number()) {
            if (c1._is_transposed) {
                auto new_c2 = copy(c2);
                new_c2->_is_vector = true;
                auto res = func_(c1) *= move(*new_c2);
                res._is_vector = false;
                res._is_transposed = false;
                delete new_c2;
                return res;
            }
            else {
                return func_(c2) *= c1;
            }
        }
        else {
            if (c1._is_transposed) {
                auto new_c2 = copy(c2);
                new_c2->_is_vector = true;
                auto res = func_(c1) *= move(*new_c2);
                delete new_c2;
                return res;
            }
            else {
                return func_(c1) *= c2;
            }
        }
    }

//    func_ operator*(func_&& f, const constant_& c){
//        return f *= c;
//    }
//
//    func_ operator*(const constant_& c, func_&& f){
//        return f *= c;
//    }

    func_ operator/(const constant_& c1, const constant_& c2){
        return func_(c1) /= c2;
    }

//    func_ operator/(func_&& f, const constant_& c){
//        return f /= c;
//    }

    //
    //func_ operator+(const func_& f1, const func_& f2){
    //    return func_(f1) += f2;
    //}
    //
    //func_ operator+(func_&& f1, const func_& f2){
    //    return f1 += f2;
    //}
    //
    //func_ operator+(const func_& f1, func_&& f2){
    //    return f2 += f1;
    //}
    //
    //
    //func_ operator-(const func_& f1, const func_& f2){
    //    return func_(f1) -= f2;
    //}
    //
    //func_ operator-(func_&& f1, const func_& f2){
    //    return f1 -= f2;
    //}
    //
    //func_ operator-(const func_& f1, func_&& f2){
    //    return (f2 *= -1) += f1;
    //}
    //
    constant_* add(constant_* c1, const func_& f){
        switch (c1->get_type()) {
            case binary_c: {
                auto res = new func_(f);
                *res += ((constant<bool>*)c1)->eval();
                delete c1;
                return c1 = res;
                break;
            }
            case short_c: {
                auto res = new func_(f);
                *res += ((constant<short>*)c1)->eval();
                delete c1;
                return c1 = res;
                break;
            }
            case integer_c: {
                auto res = new func_(f);
                *res += ((constant<int>*)c1)->eval();
                delete c1;
                return c1 = res;
                break;
            }
            case float_c: {
                auto res = new func_(f);
                *res += ((constant<float>*)c1)->eval();
                delete c1;
                return c1 = res;
                break;
            }
            case double_c: {
                auto res = new func_(f);
                *res += ((constant<double>*)c1)->eval();
                delete c1;
                return c1 = res;
                break;
            }
            case long_c: {
                auto res = new func_(f);
                *res += ((constant<long double>*)c1)->eval();
                delete c1;
                return c1 = res;
                break;
            }
            case par_c:{
                auto res = new func_(f);
                if (f.is_constant()) {
                    res->insert(true, constant<double>(1), *(param_*)c1);
                }
                else{
                    auto cst = res->get_cst();
                    cst = add(cst, *c1);
                }
                delete c1;
                return c1 = res;
                break;
            }
            case var_c:{
                auto res = new func_(f);
                delete c1;
                res->insert(true, constant<double>(1), *(param_*)c1);
                return c1 = res;
                break;
            }
            case sdpvar_c:{
                auto res = new func_(f);
                delete c1;
                res->insert(true, constant<double>(1), *(param_*)c1);
                return c1 = res;
                break;
            }
            case uexp_c: {
                auto res = new bexpr(plus_, copy(*c1), copy(f));
                delete c1;
                c1 = (constant_*)res;
                return c1;
                break;
            }
            case bexp_c: {
                auto res = new bexpr(plus_, copy(*c1), copy(f));
                delete c1;
                c1 = (constant_*)res;
                return c1;
                break;
            }
            case func_c: {
                auto res = new func_(f);
                *res += *c1;
                delete c1;
                return c1 = res;
                break;        }
                
            default:
                break;
        }
        return c1;
    }

    constant_* add(constant_* c1, const constant_& c2){ /**< adds c2 to c1, updates its type and returns the result **/
        switch (c2.get_type()) {
            case binary_c: {
                return add(c1, *(constant<bool>*)&c2);
                break;
            }
            case short_c: {
                return add(c1, *(constant<short>*)&c2);
                break;
            }
            case integer_c: {
                return add(c1, *(constant<int>*)&c2);
                break;
            }
            case float_c: {
                return add(c1, *(constant<float>*)&c2);
                break;
            }
            case double_c: {
                return add(c1, *(constant<double>*)&c2);
                break;
            }
            case long_c: {
                return add(c1, *(constant<long double>*)&c2);
                break;
            }
            case par_c:{
                auto pc2 = (param_*)(&c2);
                switch (pc2->get_intype()) {
                    case binary_:
                        return add(c1, *(param<bool>*)pc2);
                        break;
                    case short_:
                        return add(c1, *(param<short>*)pc2);
                        break;
                    case integer_:
                        return add(c1, *(param<int>*)pc2);
                        break;
                    case float_:
                        return add(c1, *(param<float>*)pc2);
                        break;
                    case double_:
                        return add(c1, *(param<double>*)pc2);
                        break;
                    case long_:
                        return add(c1, *(param<long double>*)pc2);
                        break;
                    default:
                        break;
                }
                break;
            }
            case var_c:{
                auto pc2 = (param_*)(&c2);
                switch (pc2->get_intype()) {
                    case binary_:
                        return add(c1, *(var<bool>*)pc2);
                        break;
                    case short_:
                        return add(c1, *(var<short>*)pc2);
                        break;
                    case integer_:
                        return add(c1, *(var<int>*)pc2);
                        break;
                    case float_:
                        return add(c1, *(var<float>*)pc2);
                        break;
                    case double_:
                        return add(c1, *(var<double>*)pc2);
                        break;
                    case long_:
                        return add(c1, *(var<long double>*)pc2);
                        break;
                    default:
                        break;
                }
                break;
            }
            case sdpvar_c:{
                auto pc2 = (param_*)(&c2);
                switch (pc2->get_intype()) {
                    case binary_:
                        return add(c1, *(sdpvar<bool>*)pc2);
                        break;
                    case short_:
                        return add(c1, *(sdpvar<short>*)pc2);
                        break;
                    case integer_:
                        return add(c1, *(sdpvar<int>*)pc2);
                        break;
                    case float_:
                        return add(c1, *(sdpvar<float>*)pc2);
                        break;
                    case double_:
                        return add(c1, *(sdpvar<double>*)pc2);
                        break;
                    case long_:
                        return add(c1, *(sdpvar<long double>*)pc2);
                        break;
                    default:
                        break;
                }
                break;
            }
            case uexp_c: {
                auto res = new bexpr(plus_, copy(*c1), copy(c2));
                delete c1;
                c1 = (constant_*)res;
                return c1;
                break;
            }
            case bexp_c: {
                auto res = new bexpr(plus_, copy(*c1), copy(c2));
                delete c1;
                c1 = (constant_*)res;
                return c1;
                break;
            }
            case func_c: {
                return add(c1, *(func_*)&c2);
                break;
            }
            default:
                break;
        }
        return nullptr;
    }

    constant_* substract(constant_* c1, const constant_& c2){ /**< adds c2 to c1, updates its type and returns the result **/
        switch (c2.get_type()) {
            case binary_c: {
                return substract(c1, *(constant<bool>*)&c2);
                break;
            }
            case short_c: {
                return substract(c1, *(constant<short>*)&c2);
                break;
            }
            case integer_c: {
                return substract(c1, *(constant<int>*)&c2);
                break;
            }
            case float_c: {
                return substract(c1, *(constant<float>*)&c2);
                break;
            }
            case double_c: {
                return substract(c1, *(constant<double>*)&c2);
                break;
            }
            case long_c: {
                return substract(c1, *(constant<long double>*)&c2);
                break;
            }
            case par_c:{
                auto pc2 = (param_*)(&c2);
                switch (pc2->get_intype()) {
                    case binary_:
                        return substract(c1, *(param<bool>*)pc2);
                        break;
                    case short_:
                        return substract(c1, *(param<short>*)pc2);
                        break;
                    case integer_:
                        return substract(c1, *(param<int>*)pc2);
                        break;
                    case float_:
                        return substract(c1, *(param<float>*)pc2);
                        break;
                    case double_:
                        return substract(c1, *(param<double>*)pc2);
                        break;
                    case long_:
                        return substract(c1, *(param<long double>*)pc2);
                        break;
                    default:
                        break;
                }
                break;
            }
            case var_c:{
                auto pc2 = (param_*)(&c2);
                switch (pc2->get_intype()) {
                    case binary_:
                        return substract(c1, *(var<bool>*)pc2);
                        break;
                    case short_:
                        return substract(c1, *(var<short>*)pc2);
                        break;
                    case integer_:
                        return substract(c1, *(var<int>*)pc2);
                        break;
                    case float_:
                        return substract(c1, *(var<float>*)pc2);
                        break;
                    case double_:
                        return substract(c1, *(var<double>*)pc2);
                        break;
                    case long_:
                        return substract(c1, *(var<long double>*)pc2);
                        break;
                    default:
                        break;
                }
                break;
            }
            case sdpvar_c:{
                auto pc2 = (param_*)(&c2);
                switch (pc2->get_intype()) {
                    case binary_:
                        return substract(c1, *(sdpvar<bool>*)pc2);
                        break;
                    case short_:
                        return substract(c1, *(sdpvar<short>*)pc2);
                        break;
                    case integer_:
                        return substract(c1, *(sdpvar<int>*)pc2);
                        break;
                    case float_:
                        return substract(c1, *(sdpvar<float>*)pc2);
                        break;
                    case double_:
                        return substract(c1, *(sdpvar<double>*)pc2);
                        break;
                    case long_:
                        return substract(c1, *(sdpvar<long double>*)pc2);
                        break;
                    default:
                        break;
                }
                break;
            }
            case uexp_c: {
                auto res = new bexpr(minus_, copy(*c1), copy(c2));
                delete c1;
                c1 = (constant_*)res;
                return c1;
                break;
            }
            case bexp_c: {
                auto res = new bexpr(minus_, copy(*c1), copy(c2));
                delete c1;
                c1 = (constant_*)res;
                return c1;
                break;
            }
            case func_c: {
                auto f = new func_(*c1);
                delete c1;
                *f -= c2;
                return c1 = (constant_*)f;
                break;
                break;
            }
        }
        return nullptr;
    }

    constant_* multiply(constant_* c1, const constant_& c2){ /**< adds c2 to c1, updates its type and returns the result **/
        switch (c2.get_type()) {
            case binary_c: {
                return multiply(c1, *(constant<bool>*)&c2);
                break;
            }
            case short_c: {
                return multiply(c1, *(constant<short>*)&c2);
                break;
            }
            case integer_c: {
                return multiply(c1, *(constant<int>*)&c2);
                break;
            }
            case float_c: {
                return multiply(c1, *(constant<float>*)&c2);
                break;
            }
            case double_c: {
                return multiply(c1, *(constant<double>*)&c2);
                break;
            }
            case long_c: {
                return multiply(c1, *(constant<long double>*)&c2);
                break;
            }
            case uexp_c: {
                auto res = new bexpr(product_, copy(*c1), copy(c2));
                delete c1;
                c1 = (constant_*)res;
                return c1;
                break;
            }
            case bexp_c: {
                auto res = new bexpr(product_, copy(*c1), copy(c2));
                delete c1;
                c1 = (constant_*)res;
                return c1;
                break;
            }
            case par_c:{
                auto pc2 = (param_*)(&c2);
                switch (pc2->get_intype()) {
                    case binary_:
                        return multiply(c1, *(param<bool>*)pc2);
                        break;
                    case short_:
                        return multiply(c1, *(param<short>*)pc2);
                        break;
                    case integer_:
                        return multiply(c1, *(param<int>*)pc2);
                        break;
                    case float_:
                        return multiply(c1, *(param<float>*)pc2);
                        break;
                    case double_:
                        return multiply(c1, *(param<double>*)pc2);
                        break;
                    case long_:
                        return multiply(c1, *(param<long double>*)pc2);
                        break;
                    default:
                        break;
                }
                break;
            }
            case var_c:{
                auto pc2 = (param_*)(&c2);
                switch (pc2->get_intype()) {
                    case binary_:
                        return multiply(c1, *(var<bool>*)pc2);
                        break;
                    case short_:
                        return multiply(c1, *(var<short>*)pc2);
                        break;
                    case integer_:
                        return multiply(c1, *(var<int>*)pc2);
                        break;
                    case float_:
                        return multiply(c1, *(var<float>*)pc2);
                        break;
                    case double_:
                        return multiply(c1, *(var<double>*)pc2);
                        break;
                    case long_:
                        return multiply(c1, *(var<long double>*)pc2);
                        break;
                    default:
                        break;
                }
                break;
            }
            case sdpvar_c:{
                auto pc2 = (param_*)(&c2);
                switch (pc2->get_intype()) {
                    case binary_:
                        return multiply(c1, *(sdpvar<bool>*)pc2);
                        break;
                    case short_:
                        return multiply(c1, *(sdpvar<short>*)pc2);
                        break;
                    case integer_:
                        return multiply(c1, *(sdpvar<int>*)pc2);
                        break;
                    case float_:
                        return multiply(c1, *(sdpvar<float>*)pc2);
                        break;
                    case double_:
                        return multiply(c1, *(sdpvar<double>*)pc2);
                        break;
                    case long_:
                        return multiply(c1, *(sdpvar<long double>*)pc2);
                        break;
                    default:
                        break;
                }
                break;
            }
            case func_c: {
                auto f = new func_(c2);
                *f *= *c1;
                delete c1;
                return c1 = (constant_*)f;
                break;
            }
            default:
                break;
        }
        return nullptr;
    }
        
        
        constant_* divide(constant_* c1, const constant_& c2){ /**< adds c2 to c1, updates its type and returns the result **/
            switch (c2.get_type()) {
                case binary_c: {
                    return divide(c1, *(constant<bool>*)&c2);
                    break;
                }
                case short_c: {
                    return divide(c1, *(constant<short>*)&c2);
                    break;
                }
                case integer_c: {
                    return divide(c1, *(constant<int>*)&c2);
                    break;
                }
                case float_c: {
                    return divide(c1, *(constant<float>*)&c2);
                    break;
                }
                case double_c: {
                    return divide(c1, *(constant<double>*)&c2);
                    break;
                }
                case long_c: {
                    return divide(c1, *(constant<long double>*)&c2);
                    break;
                }
                case uexp_c: {
                    auto res = new bexpr(div_, copy(*c1), copy(c2));
                    delete c1;
                    c1 = (constant_*)res;
                    return c1;
                    break;
                }
                case bexp_c: {
                    auto res = new bexpr(div_, copy(*c1), copy(c2));
                    delete c1;
                    c1 = (constant_*)res;
                    return c1;
                    break;
                }
                default:{
                    auto f = new func_(*c1);
                    *f /= func_(c2);
                    delete c1;
                    return c1 = (constant_*)f;
                    break;
                }
            }
            return nullptr;
        }
        
    string pterm::to_str(int ind) const{
        string str;
        constant_* c_new = _coef;
        if (c_new->is_number()){
            string v = poly_to_str(c_new);
            if (_sign) {
                if (v=="-1") {
                    str += " - ";
                }
                else if (ind>0) {
                    str += " + ";
                    if(v!="1") {
                        str += v;
                    }
                }
                else if(v!="1") {
                    str += v;
                }
            }
            if(!_sign) {
                if (v == "-1" && ind>0) {
                    str += " + ";
                }
                else if (v.front()=='-'){
                    if (ind > 0) {
                        str += " + ";
                    }
                    str += v.substr(1);
                }
                else if (v=="1"){
                    str += " - ";
                }
                else if(v!="-1"){
                    str += " - " + v;
                }
            }
        }
        else{
            if (!_sign) {
                str += " - ";
            }
            if(ind > 0 && _sign) {
                str += " + ";
            }
            str += "(";
            str += poly_to_str(c_new);
            str += ")";
        }
        for (auto& p: *_l) {
            str += poly_to_str(p.first);
            if (p.second != 1) {
                str += "^" + to_string(p.second);
            }
        }
        return str;
    }

    void pterm::print(int ind) const{
        cout << this->to_str(ind);
    }

    string qterm::to_str(int ind) const {
        string str;
        constant_* c_new = _coef;
        param_* p_new1 = (param_*)_p->first;
        param_* p_new2 = (param_*)_p->second;
        if (c_new->is_number()){
            string v = poly_to_str(c_new);
            if (_sign) {
                if (v=="-1") {
                    str += " - ";
                }
                else if (ind>0) {
                    str += " + ";
                    if(v!="1") {
                        str += v;
                    }
                }
                else if(v!="1") {
                    str += v;
                }
            }
            if(!_sign) {
                if (v == "-1" && ind>0) {
                    str += " + ";
                }
                else if (v.front()=='-'){
                    if (ind > 0) {
                        str += " + ";
                    }
                    str += v.substr(1);
                }
                else if (v=="1"){
                    str += " - ";
                }
                else if(v!="-1"){
                    str += " - " + v;
                }
            }
        }
        else{
            if (!_sign) {
                str += " - ";
            }
            if(ind > 0 && _sign) {
                str += " + ";
            }
            str += "(";
            str += poly_to_str(c_new);
            str += ")";
        }
        str += poly_to_str(p_new1);
        if (p_new1==p_new2) {
            str += "^2";
        }
        else {
            str += "*";
            str += poly_to_str(p_new2);
        }
        return str;
    }


    void qterm::print(int ind) const {
        cout << this->to_str(ind);
    }


    string lterm::to_str(int ind) const{
        string str;
        constant_* c_new = _coef;
        param_* p_new = (param_*)_p;
        if (c_new->is_number()){
            string v = poly_to_str(c_new);
            if (_sign) {
                if (v=="-1") {
                    str += " - ";
                }
                else if (ind>0) {
                    str += " + ";
                    if(v!="1") {
                        str += v;
                    }
                }
                else if(v!="1") {
                    str += v;
                }
            }
            if(!_sign) {
                if (v == "-1" && ind>0) {
                    str += " + ";
                }
                else if (v.front()=='-'){
                    if (ind > 0) {
                        str += " + ";
                    }
                    str += v.substr(1);
                }
                else if (v=="1"){
                    str += " - ";
                }
                else if(v!="-1"){
                    str += " - " + v;
                }
            }
        }
        else{
            if (!_sign) {
                str += " - ";
            }
            if(ind > 0 && _sign) {
                str += " + ";
            }
            str += "(";
            str += poly_to_str(c_new);
            str += ")";
        }
        str += poly_to_str(p_new);
        return str;
    }

    void lterm::print(int ind) const{
        cout << this->to_str(ind);
    }

    func_* func_::get_stored_derivative(const unique_id& vid) const{
        auto it = _dfdx.find(vid);
        if (it!=_dfdx.end()) {
            return it->second;
        }
        else {
            throw invalid_argument("No derivatives stored!\n");
        }
    }

     func_* func_::compute_derivative(const param_ &v){
         auto vid = v._unique_id;
        if(_dfdx.count(vid)!=0){
            return _dfdx[vid];
        }
         
        auto df = new func_(get_derivative(v));
    //    embed(*df);
//         if(_dfdx.count(v._unique_id)!=0){
//             throw invalid_argument("storing same df twice! " + poly_to_str(&v) + this->to_str());
//         }
        _dfdx[vid] = df;
        Debug( "First derivative with respect to " << v.get_name() << " = ");
    //    df->print();
        return df;
    }

    void func_::compute_derivatives(){ /**< Computes and stores the derivative of f with respect to all variables. */
        size_t vid = 0, vjd = 0;
        param_* vi;
        param_* vj;
        for (auto &vp: *_vars) {
            vi = vp.second.first;
            vid = vi->get_id();
            auto df = compute_derivative(*vi);
            for (auto &vp2: *_vars) {
                vj = vp2.second.first;
                vjd = vj->get_id();
                if (df->has_var(*vj)) { //only store lower left part of hessian matrix since it is symmetric.
                    auto d2f = df->compute_derivative(*vj);
                    Debug( "Second derivative with respect to " << vp2.first << " and " << vp.first << " = " << d2f.to_str());
    //                d2f->print();
                }
            }
            
        }
    };

    bool func_::has_var(const param_& v) const{
        return _vars->count(v.get_name())>0;
    }

    func_ func_::get_derivative(const param_ &v) const{
        func_ res;
        for (auto &lt: *_lterms) {
            if (*lt.second._p == v) {
                if(lt.second._coef->_is_transposed){
                    lt.second._coef->_dim = v._dim;
                    constant<double> ones(1);
                    if(lt.second._sign){
                        res += ones.tr()*(*lt.second._coef);
                    }
                    else {
                        res -= ones.tr()*(*lt.second._coef);
                    }
                }
                else {
                    if(lt.second._sign){
                        res += (*lt.second._coef);
                    }
                    else {
                        res -= (*lt.second._coef);
                    }
                }
            }
        }
        for (auto &lt: *_qterms) {
            if (*lt.second._p->first == v) {
                if(lt.second._sign) {
                    res += *lt.second._coef*(*lt.second._p->second);
                }
                else {
                    res -= *lt.second._coef*(*lt.second._p->second);
                }
            }
            if (*lt.second._p->second == v) {
                if(lt.second._sign) {
                    res += *lt.second._coef*(*lt.second._p->first);
                }
                else {
                    res -= *lt.second._coef*(*lt.second._p->first);
                }
            }
        }
        for (auto &lt: *_pterms) {
            for (auto &p: *lt.second._l) {
                if (*p.first == v) {
                    func_ pterm = constant<>(1);
                    if (!lt.second._sign) {
                        pterm = constant<>(-1);
                    }
                    auto expo = p.second;
                    if (expo > 1) {
                        pterm *= expo;
                        pterm *= *lt.second._coef;
                        pterm *= (*p.first);
                        for (int i = 1; i<expo-1; i++) {
                            pterm *= *p.first;
                        }
                    }
                    for (auto &p2: *lt.second._l) {
                        if (p2!=p) {
                            func_ pterm2(*p2.first);
                            for (int i = 1; i<p2.second; i++) {
                                pterm2 *= *p2.first;
                            }
                            pterm *= pterm2;
                        }
                    }
                    res += *lt.second._coef*pterm;
                }
            }
        }
        if (!_expr) {
            return res;
        }
        else { // f is a composition of functions
            return res += get_poly_derivative(_expr, v);
        }
        return res;
    }

    double func_::eval(size_t i) const{
        double res = 0;
        for (auto &pair:*_pterms) {
            res += pair.second.eval(i);
        }
        for (auto &pair:*_qterms) {
            res += pair.second.eval(i);
        }
        for (auto &pair:*_lterms) {
            res += pair.second.eval(i);
        }
        res += poly_eval(_cst,i);
        if (_expr) {
            res += poly_eval(_expr, i);
        }
        return res;
    }

    string func_::to_str(bool display_input) const{
        string str;
        int ind = 0;
        string sign = " + ";
//        for (int inst = 0; inst < _nb_instances ; inst++) {
            ind = 0;
            if (display_input) {
                if (!_embedded && !is_constant() && _all_convexity==convex_) {
                    str += "Convex function: ";
                }
                if (!_embedded && !is_constant() && _all_convexity==concave_) {
                    str += "Concave function: ";
                }
                if (!_embedded && !is_constant()) {
                    str += "f(";
                    for (auto pair_it = _vars->begin(); pair_it != _vars->end();) {
        //                if (!pair_it->second.first->_is_vector) {
//                        str += pair_it->second.first->get_name();
                            str += pair_it->second.first->get_name()+"[";
                        str += to_string(pair_it->second.first->get_id_inst())+"]";
                            if (next(pair_it) != _vars->end()) {
                                str += ",";
                            }
        //                }
                        pair_it++;
                    }
                    str += ") = ";
                }
            }
            for (auto &pair:*_pterms) {
                str += pair.second.to_str(ind++);
            }
            if (!_pterms->empty() && (!_qterms->empty() || !_lterms->empty())) {
                str += " + ";
            }
            ind = 0;
            for (auto &pair:*_qterms) {
                str += pair.second.to_str(ind++);
            }
            if (!_qterms->empty() && !_lterms->empty()) {
                str += " + ";
            }
            ind = 0;
            for (auto &pair:*_lterms) {
                str += pair.second.to_str(ind++);
            }
            if (_cst->is_number()) {
                auto val = poly_to_str(_cst);
                if (val.front()=='-') {
                    str += " - " + val.substr(1);
                }
                else if (val != "0"){
                    if (!_pterms->empty() || !_qterms->empty() || !_lterms->empty()) {
                        str += " + ";
                    }
                    str += val;
                }
            }
            else {
                if (!_pterms->empty() || !_qterms->empty() || !_lterms->empty()) {
                    str += " + ";
                }
                str += "(";
                str += poly_to_str(_cst);
                str += ")";
            }
            if (!_queue->empty() && (!_pterms->empty() || !_qterms->empty() || !_lterms->empty() || !_cst->is_zero())) {
                str += " + ";
            }
            if (!_queue->empty()) {
                str += _expr->get_str();
            }
            if (_is_vector) {
                str = "[" + str +"]";
            }
            if (_is_transposed) {
                str += "^T";
            }
//            str += "\n";
//        }
        return str;
    }


    void func_::print(bool endline, bool display_input) const{
        cout << this->to_str(display_input);
        if (endline)
            cout << endl;
    }


    /* UNARY EXPRESSIONS */

    uexpr::uexpr(){
        _otype = id_;
        _son = nullptr;
        _type = uexp_c;
        _to_str = "null";
    }

    string uexpr::to_str() const{
        string str;
        switch (_otype) {
            case log_:
                str += "log(";
                str += poly_to_str(_son);
                str += ")";
                break;
                
            case exp_:
                str += "exp(";
                str += poly_to_str(_son);
                str += ")";
                break;
                
            case cos_:
                str += "cos(";
                str += poly_to_str(_son);
                str += ")";
                break;
                
            case sin_:
                str += "sin(";
                str += poly_to_str(_son);
                str += ")";
                break;
                
            case sqrt_:
                str += "sqrt(";
                str += poly_to_str(_son);
                str += ")";
                break;
            default:
                break;
        }
        return str;
    }

    void uexpr::print(bool endline) const{
        cout << _to_str;
        if(endline) {
            cout << endl;
        }

    }
    
    vector<param_*> uexpr::get_nl_vars() const{
        vector<param_*> res;
        if (_son->is_function()) {
            auto vars = ((func_*)_son)->get_vars();
            for (auto &pairs: vars){
                res.push_back(pairs.second.first);
            }
        }
        else if(_son->is_uexpr()) {
            return ((uexpr*)_son)->get_nl_vars();
        }
        else if (_son->is_bexpr()){
            return ((bexpr*)_son)->get_nl_vars();
        }
        return res;
    }
    
    func_ uexpr::get_derivative(const param_ &v) const{
       // Unary operators
       // f(g(x))' = f'(g(x))*g'(x).
       switch (_otype) {
           case cos_:
               return -1*get_poly_derivative(_son,v)* sin(*_son);
               break;
           case sin_:
               return get_poly_derivative(_son,v)*cos(*_son);
               break;
           case sqrt_:
               return get_poly_derivative(_son,v)/(2*sqrt(*_son));
               break;
           case exp_:
               return get_poly_derivative(_son,v)*expo(*_son);
               break;
           case log_:
               return get_poly_derivative(_son,v)/(*_son);
               break;
           default:
               std::cerr << "ok unsupported unary operation";
               exit(-1);
               break;
       }
        return func_();
    }
    
    vector<param_*> bexpr::get_nl_vars() const{
        vector<param_*> res;
        if (_lson->is_function()) {
            auto vars = ((func_*)_lson)->get_vars();
            for (auto &pairs: vars){
                res.push_back(pairs.second.first);
            }
        }
        else if(_lson->is_uexpr()) {
            res = ((uexpr*)_lson)->get_nl_vars();
        }
        else if (_lson->is_bexpr()){
            res = ((bexpr*)_lson)->get_nl_vars();
        }
        if (_rson->is_function()) {
            auto vars = ((func_*)_rson)->get_vars();
            for (auto &pairs: vars){
                res.push_back(pairs.second.first);
            }
        }
        else if(_rson->is_uexpr()) {
            auto vars = ((uexpr*)_rson)->get_nl_vars();
            res.insert(res.end(), vars.begin(), vars.end() );
        }
        else if (_rson->is_bexpr()){
            auto vars = ((bexpr*)_rson)->get_nl_vars();
            res.insert(res.end(), vars.begin(), vars.end() );
        }
        return res;
    }
    
    func_ bexpr::get_derivative(const param_ &v) const{
        switch (_otype) {
            case plus_:
                return get_poly_derivative(_lson,v) + get_poly_derivative(_rson,v);
                break;
            case minus_:
                return get_poly_derivative(_lson,v) - get_poly_derivative(_rson,v);
                break;
            case product_:
                return get_poly_derivative(_lson,v)*(*_rson) + (*_lson)*get_poly_derivative(_rson,v);
                // f'g + fg'
                break;
            case div_:
                return (get_poly_derivative(_lson,v)*(*_rson) - (*_lson)*get_poly_derivative(_rson,v))/((*_rson)*(*_rson));
                // (f'g - fg')/g^2
                break;
            case power_:
                if (!_rson->is_number()) {
                    throw invalid_argument("Function in exponent not supported yet.\n");
                }
//                auto exponent = poly_eval(_rson);
//                return (exponent*get_poly_derivative(_lson,v)*power(*_lson, exponent-1));// nf'f^n-1
                break;
            default:
                throw invalid_argument("unsupported operation");
                break;
        }
        return func_();
    }

    string bexpr::to_str() const{
        string str;
        if((_otype==product_ || _otype==div_) && (_lson->get_type()==uexp_c || _lson->get_type()==bexp_c)) {
            str += "(";
            str+= poly_to_str(_lson);
            str += ")";
        }
        else
            str+= poly_to_str(_lson);

        if (_otype==plus_) {
            str+= " + ";
        }
        if (_otype==minus_) {
            str+= " - ";
        }
        if (_otype==product_) {
            str+= " * ";
        }
        if (_otype==div_) {
            str+= "/";
        }

        if (_otype==power_) {
            str+= "^";
        }

        if (_otype==plus_ || (_rson->get_type()!=uexp_c && _rson->get_type()!=bexp_c)) {
            str+= poly_to_str(_rson);
        }
        else {
            str+= "(";
            str+= poly_to_str(_rson);
            str+= ")";
        }
        return str;
    }

    void bexpr::print(bool endline) const {
        cout << _to_str;
    //    if((_otype==product_ || _otype==div_) && (_lson->get_type()==uexp_c || _lson->get_type()==bexp_c)) {
    //        cout << "(";
    //        poly_print(_lson);
    //        cout << ")";
    //    }
    //    else
    //        poly_print(_lson);
    //    
    //    if (_otype==plus_) {
    //        cout << " + ";
    //    }
    //    if (_otype==minus_) {
    //        cout << " - ";
    //    }
    //    if (_otype==product_) {
    //        cout << " * ";
    //    }
    //    if (_otype==div_) {
    //        cout << "/";
    //    }
    //    
    //    if (_otype==power_) {
    //        cout << "^";
    //    }
    //    
    //    if (_otype==plus_ || (_rson->get_type()!=uexp_c && _rson->get_type()!=bexp_c)) {
    //        poly_print(_rson);
    //    }
    //    else {
    //        cout << "(";
    //        poly_print(_rson);
    //        cout << ")";
    //    }
        if(endline)
            cout << endl;
    }

    size_t func_::get_nb_vars() const{
        size_t n = 0;
        for (auto &p: *_vars) {
            if (p.second.first->_is_vector) {
                n += p.second.first->get_dim();
            }
            else {
                n += 1;
            }
        }
        return n;
    }

    size_t func_::get_nb_instances() const{
        return _nb_instances;
    }


    constant_* func_::get_cst() {
        return _cst;
    }

    param_* func_::get_var(string name){
        if (_vars->empty()) {
            return nullptr;
        }
        auto pair_it = _vars->find(name);
        if (pair_it==_vars->end()) {
            return nullptr;
        }
        else {
            return get<1>(*pair_it).first;
        }
    }

    param_* func_::get_param(string name){
        auto pair_it = _params->find(name);
        if (pair_it==_params->end()) {
            return nullptr;
        }
        else {
            return get<1>(*pair_it).first;
        }
    }


    void func_::add_var(param_* v, int nb){/**< Inserts the variable in this function input list. nb represents the number of occurences v has. WARNING: Assumes that v has not been added previousely!*/
        if (_vars->count(v->get_name())!=0) {
            v->set_id(_vars->size());
        }
        _vars->insert(make_pair<>(v->get_name(), make_pair<>(v, nb)));
        if (!v->_is_vector) {
            _nb_instances = max(_nb_instances, v->get_nb_instances());
        }
    }


    void func_::add_param(param_* p){/**< Inserts the parameter in this function input list. WARNING: Assumes that p has not been added previousely!*/
        if (_params->count(p->get_name())!=0) {
            p->set_id(_params->size());
        }
        _params->insert(make_pair<>(p->get_name(), make_pair<>(p, 1)));
        if (!p->_is_vector) {
            _nb_instances = max(_nb_instances, p->get_nb_instances());
        }
    }



    void func_::delete_var(const string& vid){
        auto vit = _vars->find(vid);
        if (vit==_vars->end()) {
            return;
        }
        _vars->erase(vit);
    }

    void func_::delete_param(const string& vid){
        auto vit = _params->find(vid);
        if (vit==_params->end()) {
            return;
        }
        _params->erase(vit);
    }




    int func_::nb_occ_var(string name) const{/**< Returns the number of occurences the variable has in this function. */
        auto pair_it = _vars->find(name);
        if (pair_it==_vars->end()) {
            return 0;
        }
        else {
            return get<1>(*pair_it).second;
        }
    }

    int func_::nb_occ_param(string name) const{/**< Returns the number of occurences the parameter has in this function. */
        auto pair_it = _params->find(name);
        if (pair_it==_params->end()) {
            return 0;
        }
        else {
            return get<1>(*pair_it).second;
        }
    }

    void func_::incr_occ_var(string str){/**< Increases the number of occurences the variable has in this function. */
        auto pair_it = _vars->find(str);
        if (pair_it==_vars->end()) {
            throw invalid_argument("Non-existing variable in function!\n");
        }
        else {
            get<1>(*pair_it).second++;
        }
    }

    void func_::incr_occ_param(string str){/**< Increases the number of occurences the parameter has in this function. */
        auto pair_it = _params->find(str);
        if (pair_it==_params->end()) {
            throw invalid_argument("Non-existing variable in function!\n");
        }
        else {
            get<1>(*pair_it).second++;
        }
    }

    void func_::decr_occ_var(string str, int nb){/**< Decreases the number of occurences the variable has in this function by nb. */
        auto pair_it = _vars->find(str);
        if (pair_it==_vars->end()) {
            return;
        }
        else {
            get<1>(*pair_it).second-=nb;
            if (get<1>(*pair_it).second==0) {
                delete get<1>(*pair_it).first;
                _vars->erase(pair_it);
            }
        }
    }

    void func_::decr_occ_param(string str, int nb){/**< Decreases the number of occurences the parameter has in this function by nb. */
        auto pair_it = _params->find(str);
        if (pair_it==_params->end()) {
            return;
        }
        else {
            get<1>(*pair_it).second -= nb;
            if (get<1>(*pair_it).second==0) {
                delete get<1>(*pair_it).first;                
                _params->erase(pair_it);
            }
        }
    }



    bool func_::is_convex(int idx) const{
        return (_convexity->at(idx)==convex_ || _convexity->at(idx)==linear_);
    };

    bool func_::is_concave(int idx) const{
        return (_convexity->at(idx)==concave_ || _convexity->at(idx)==linear_);
    };


    bool func_::is_number() const{
        return (_vars->empty() && _params->empty());
    }

    bool func_::is_constant() const{
        return (_ftype==const_);
    }

    bool func_::is_linear() const{
        return (_ftype==lin_);
    };

    bool func_::is_quadratic() const{
        return (_ftype==quad_);
    };

    bool func_::is_polynomial() const{
        return (_ftype==pol_);
    };

    bool func_::is_nonlinear() const{
        return (_ftype==nlin_);
    };

    bool func_::is_zero(){/*<< A function is zero if it is constant and equals zero or if it is a sum of zero valued parameters */
        if (_ftype==const_ && _cst->is_zero()){
            for (auto& it:*_params) {
                if (!it.second.first->is_zero()) {
                    return false;
                }
            }
            return true;
        }
        return false;
    }

    bool func_::is_transposed() const {
        return _is_transposed;
    }

    FType func_::get_ftype() const { return _ftype;}


    qterm* func_::get_square(param_* p){ /**< Returns the quadratic term containing a square of p or nullptr if none exists. **/
        for (auto pair_it = _qterms->begin(); pair_it != _qterms->end(); pair_it++) {
            if (pair_it->second._p->first==p && pair_it->second._p->second==p) {
                return &pair_it->second;
            }
        }
        return nullptr;
    }

    func_ func_::get_outer_app(){ /**< Returns an outer-approximation of the function using the current value of the variables **/
        func_ res; // res = gradf(x*)*(x-x*) + f(x*)
        param_* v;
        for(auto &it: *_vars){
            v = it.second.first;
            res += (get_stored_derivative(v->_unique_id)->eval())*((*v) - poly_eval(v));
        }
        res += eval();
        return res;
    }

    Sign func_::get_all_sign() const{
        return _all_sign;
    }

    Sign func_::get_sign(int idx) const{
        return _sign->at(idx);
    }

    Sign func_::get_all_sign(const lterm& l) {
        if (l._coef->is_zero()) {
            return zero_;
        }
        if (l._coef->get_all_sign()==unknown_ || l._p->get_all_sign()==unknown_) {
            return unknown_;
        }
        auto s = l._coef->get_all_sign() * l._p->get_all_sign();
        if(s == 1 || s == 2) {
            if (l._sign) {
                return non_neg_;
            }
            else {
                return non_pos_;
            }
        }
        if(s == 4) {
            if (l._sign) {
                return pos_;
            }
            else {
                return neg_;
            }
        }
        if(s == -1 || s == -2) {
            if (l._sign) {
                return non_pos_;
            }
            else{
                return non_neg_;
            }
        }
        if(s == -4) {
            if (l._sign) {
                return neg_;
            }
            else {
                return pos_;
            }
        }
        return unknown_;
    }

    Sign func_::get_all_sign(const qterm& l) {
        if (l._coef->is_zero()) {
            return zero_;
        }
        if (l._coef->get_all_sign()==unknown_ || l._p->first->get_all_sign()==unknown_ || l._p->second->get_all_sign()==unknown_) {
            return unknown_;
        }
        auto s = l._coef->get_all_sign() * l._p->first->get_all_sign() * l._p->second->get_all_sign();
        if(s == 1 || s == 2 || s == 4) {
            if (l._sign) {
                return non_neg_;
            }
            else {
                return non_pos_;
            }
        }
        if(s == 8) {
            if (l._sign) {
                return pos_;
            }
            else {
                return neg_;
            }
        }
        if(s == -1 || s == -2 || s == -4) {
            if (l._sign) {
                return non_pos_;
            }
            else{
                return non_neg_;
            }
        }
        if(s == -8) {
            if (l._sign) {
                return neg_;
            }
            else {
                return pos_;
            }
        }
        return unknown_;
    }

    Sign func_::get_all_sign(const pterm& l) {
        if (l._coef->is_zero()) {
            return zero_;
        }
        //        if (l._coef->get_all_sign()==unknown_ || l._p->first->get_all_sign()==unknown_ || l._p->second->get_all_sign()==unknown_) {
        //            return unknown_;
        //        }
        //        if (l._sign) {
        //            if(l._coef->get_all_sign() * l._p->first->get_all_sign() * l._p->second->get_all_sign() == 2) {
        //                return non_neg_;
        //            }
        //            if(l._coef->get_all_sign() * l._p->first->get_all_sign() * l._p->second->get_all_sign() == 4) {
        //                return pos_;
        //            }
        //            if(l._coef->get_all_sign() * l._p->first->get_all_sign() * l._p->second->get_all_sign() == -2) {
        //                return non_pos_;
        //            }
        //            if(l._coef->get_all_sign() * l._p->first->get_all_sign() * l._p->second->get_all_sign() == -4) {
        //                return neg_;
        //            }
        //        }
        //        else {
        //            if(l._coef->get_all_sign() * l._p->first->get_all_sign() * l._p->second->get_all_sign() == 2) {
        //                return non_pos_;
        //            }
        //            if(l._coef->get_all_sign() * l._p->first->get_all_sign() * l._p->second->get_all_sign() == 4) {
        //                return neg_;
        //            }
        //            if(l._coef->get_all_sign() * l._p->first->get_all_sign() * l._p->second->get_all_sign() == -2) {
        //                return non_neg_;
        //            }
        //            if(l._coef->get_all_sign() * l._p->first->get_all_sign() * l._p->second->get_all_sign() == -4) {
        //                return pos_;
        //            }
        //        }
        return unknown_;
    }



    Convexity func_::get_convexity(const qterm& q) {
        if(q._p->first == q._p->second){
            if (q._sign && (q._coef->is_positive() || q._coef->is_non_negative())) {
                return convex_;
            }
            if (q._sign && (q._coef->is_negative() || q._coef->is_non_positive())) {
                return concave_;
            }
            if (!q._sign && (q._coef->is_negative() || q._coef->is_non_positive())) {
                return convex_;
            }
            if (!q._sign && (q._coef->is_negative() || q._coef->is_non_positive())) {
                return concave_;
            }
        }
        // At this stage, we know that q._p->first !=q._p->second
        // Checking if the product can be factorized
        auto sqr1 = get_square(q._p->first);
        auto sqr2 = get_square(q._p->second);
        if (sqr1 && sqr2){
            auto c1 = sqr1->_coef;
            auto c2 = sqr2->_coef;
            if (!(sqr1->_sign^c1->is_positive())==!(sqr2->_sign^c2->is_positive())) {// && c0->is_at_least_half(c1) && c0->is_at_least_half(c2)
                if (!(sqr1->_sign^c1->is_positive())) {
                    return convex_;
                }
                return concave_;
            }
        }
        return undet_;
    }

    template<typename type>
    func_ power(const param<type>& v, unsigned p){
        func_ res(v);
        for (int i = 1; i < p; i++) {
            res *= v;
        }
        return res;
    }
    template func_ power<double>(const param<double>& v, unsigned p);
    template func_ power<float>(const param<float>& v, unsigned p);
    template func_ power<int>(const param<int>& v, unsigned p);
    template func_ power<long double>(const param<long double>& v, unsigned p);
    template func_ power<short>(const param<short>& v, unsigned p);
    template func_ power<bool>(const param<bool>& v, unsigned p);
    
//    template<typename type>
//    func_ sum(const var<type>& v){
//        func_ res;
//        if (v.get_dim()==0) {
//            return res;
//        }
//        return constant<double>(1).tr()*v;
//    }
    
    template<typename type>
    func_ sum(const param<type>& p){
        func_ res;
        if (p.get_dim()==0 || p.is_zero()) {
            return res;
        }
        return constant<double>(1).tr()*p;
    }
    
    template<typename type>
    func_ sum(const param<type>& p1, const func_& f){
        func_ res;
        if (p1.get_dim()==0 || p1.is_zero() || f.constant_::is_zero()) {
            return res;
        }
        return p1.tr()*f;
    }
    
    template func_ sum<double>(const param<double>& v);
    template func_ sum<float>(const param<float>& v);
    template func_ sum<long double>(const param<long double>& v);
    template func_ sum<int>(const param<int>& v);
    template func_ sum<short>(const param<short>& v);
    template func_ sum<bool>(const param<bool>& v);

    template func_ sum<double>(const param<double>& v1, const func_& f);
    template func_ sum<float>(const param<float>& v1, const func_& f);
    template func_ sum<long double>(const param<long double>& v1, const func_& f);
    template func_ sum<int>(const param<int>& v1, const func_& f);
    template func_ sum<short>(const param<short>& v1, const func_& f);
    template func_ sum<bool>(const param<bool>& v1, const func_& f);
    
    void func_::update_convexity(){
        if (!_pterms->empty()) {
            _all_convexity = undet_;
            return;
        }
        if (_qterms->empty()) {
            _all_convexity = linear_;
            return;
        }
        _all_convexity = get_convexity(_qterms->begin()->second);
        for (auto pair_it = next(_qterms->begin()); pair_it != _qterms->end(); pair_it++) {
            Convexity conv = get_convexity(pair_it->second);
            if (_all_convexity==undet_ || conv ==undet_ || (_all_convexity==convex_ && conv==concave_) || (_all_convexity==concave_ && conv==convex_)) {
                _all_convexity = undet_;
                return;
            }
            else {
                _all_convexity = conv;
            }
        }
    }

}
