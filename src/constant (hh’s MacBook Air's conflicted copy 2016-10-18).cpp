//
// Created by Hassan on 19/11/2015.
//

#include <math.h>
#include <Gravity/constant.h>
#include <Gravity/param.h>
#include <sstream>

/* Polymorphic functions */

bool equals(const constant_* c1, const constant_* c2){/**< Copy c2 into c1 detecting the right class, i.e., constant<>, param<>, uexpr or bexpr. */
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
        case parameter:{
            return (c1->is_param() && *(param_ *)c1 == *(param_ *)c2);
            break;
        }
        case unary_exp: {
            return (c1->is_uexpr() && *(uexpr *)c1 == *(uexpr *)c2);
            break;
        }
        case binary_exp: {
            return (c1->is_bexpr() && *(bexpr *)c1 == *(bexpr *)c2);
            break;
        }
        default:
            break;
    }
    
}


void copy(constant_* c1, constant_* c2){/**< Copy c2 into c1 detecting the right class, i.e., constant<>, param<>, uexpr or bexpr. */
    if (!c1 && !c2) {
        return;
    }
    if (c1) {
        delete c1;
        c1 = nullptr;
    }
    if (!c2) {
        return;
    }
    switch (c2->get_type()) {
        case binary_c: {
            c1 = new constant<bool>(((constant<bool>*)(c2))->get_val());
            break;
        }
        case short_c: {
            c1 = new constant<short>(((constant<short>*)(c2))->get_val());
            break;
        }
        case integer_c: {
            c1 = new constant<int>(((constant<int>*)(c2))->get_val());
            break;
        }
        case float_c: {
            c1 = new constant<float>(((constant<float>*)(c2))->get_val());
            break;
        }
        case double_c: {
            c1 = new constant<double>(((constant<double>*)(c2))->get_val());
            break;
        }
        case long_c: {
            c1 = new constant<long double>(((constant<long double>*)(c2))->get_val());
            break;
        }
        case parameter:{
            auto p_c2 = (param_*)(c2);
            switch (p_c2->get_intype()) {
                case binary_p:
                    c1 = new param<bool>(*(param<bool>*)p_c2);
                    break;
                case short_p:
                    c1 = new param<short>(*(param<short>*)p_c2);
                    break;
                case integer_p:
                    c1 = new param<int>(*(param<int>*)p_c2);
                    break;
                case float_p:
                    c1 = new param<float>(*(param<float>*)p_c2);
                    break;
                case double_p:
                    c1 = new param<double>(*(param<double>*)p_c2);
                    break;
                case long_p:
                    c1 = new param<long double>(*(param<long double>*)p_c2);
                    break;
                default:
                    break;
            }
            break;
        }
        case unary_exp: {
            c1 = new uexpr(*(uexpr*)c2);
            break;
        }
        case binary_exp: {
            c1 = new bexpr(*(bexpr*)c2);
            break;
        }
        default:
            break;
    }
    
}

/* UNARY EXPRESSIONS */

uexpr::uexpr(const uexpr& exp){
    _otype = exp._otype;
    copy(_son, exp._son);
}

uexpr::uexpr(uexpr&& exp){
    _otype = exp._otype;
    _son = move(exp._son);
    exp._son = nullptr;
}


bool uexpr::contains(const constant_* c) const{
    if (!_son) {
        return false;
    }
    if (equals(_son, c)) {
        return true;
    }
    
    if (_son->get_type()==unary_exp) {
        return ((uexpr*)_son)->contains(c);
    }
    if (_son->get_type()==binary_exp) {
        return ((bexpr*)_son)->contains(c);
    }
    return false;
};


bool uexpr::operator==(const uexpr &c)const{
    return (_otype == c._otype && equals(_son,c._son));
    
}


uexpr cos(const constant_& c){
    uexpr res;
    res._otype = cos_;
    copy(res._son,(constant_*)&c);
    return res;
};


uexpr sin(const constant_& c){
    uexpr res;
    res._otype = sin_;
    copy(res._son,(constant_*)&c);
    return res;
};


uexpr sqrt(const constant_& c){
    uexpr res;
    res._otype = sqrt_;
    copy(res._son,(constant_*)&c);
    return res;
};

uexpr expo(const constant_& c){
    uexpr res;
    res._otype = exp_;
    copy(res._son,(constant_*)&c);
    return res;
};

uexpr log(const constant_& c){
    uexpr res;
    res._otype = log_;
    copy(res._son,(constant_*)&c);
    return res;
};

/* BINARY EXPRESSIONS */

bool bexpr::operator==(const bexpr &c)const{
    return (_otype == c._otype && equals(_lson,c._lson) && equals(_rson,c._rson));
}


bool bexpr::contains(const constant_* c) const{
    if(_lson){
        if (equals(_lson, c)) {
            return true;
        }
        if (_lson->get_type()==unary_exp && ((uexpr*)_lson)->contains(c)) {
            return true;
        }
        if (_lson->get_type()==binary_exp && ((bexpr*)_lson)->contains(c)) {
            return true;
        }
    }
    if(_rson) {
        if (equals(_rson, c)) {
            return true;
        }
        if (_rson->get_type()==unary_exp && ((uexpr*)_rson)->contains(c)) {
            return true;
        }
        if (_rson->get_type()==binary_exp && ((bexpr*)_rson)->contains(c)) {
            return true;
        }
    }
    return false;
};


bexpr::bexpr(const bexpr& exp){ /**< Constructor from unary expression tree */
    _otype = exp._otype;
    copy(_lson, exp._lson);
    copy(_rson, exp._rson);
    _type = binary_exp;
};

bexpr::bexpr(bexpr&& exp){ /**< Move constructor from unary expression tree */
    _otype = exp._otype;
    delete _lson;
    _lson = move(exp._lson);
    exp._lson = nullptr;
    delete _rson;
    _rson = move(exp._rson);
    _type = binary_exp;
};


//constant::constant(constant&& c){ /**< Move constructor */
//    _content = move(c._content);
//    c._content = nullptr;
//}
//
//constant::constant(const constant& cc){ /**< Copy constructor */
//    switch (cc.get_type()) {
//        case binary_c: {
//            auto c = static_cast<content<bool>*>(cc._content);
//            _content = new content<bool>(c->get_val());
//            break;
//        }
//        case short_c: {
//            auto c = static_cast<content<short>*>(cc._content);
//            _content = new content<short>(c->get_val());
//            break;
//        }
//        case integer_c: {
//            auto c = static_cast<content<int>*>(cc._content);
//            _content = new content<short>(c->get_val());
//            break;
//        }
//        case float_c: {
//            auto c = static_cast<content<float>*>(cc._content);
//            _content = new content<short>(c->get_val());
//            break;
//        }
//        case double_c: {
//            auto c = static_cast<content<double>*>(cc._content);
//            _content = new content<short>(c->get_val());
//            break;
//        }
//        case long_c: {
//            auto c = static_cast<content<long double>*>(cc._content);
//            _content = new content<short>(c->get_val());
//            break;
//        }
//        case parameter: {
//            auto c = static_cast<content<param_*>*>(cc._content);
//            auto arg = c->get_arg();
//            switch (arg->get_type()) {
//                case binary_p: {
//                    auto c = static_cast<param<bool>*>(arg);
//                    _content = new content<param_*>(new param<bool>(c->get_val()));
//                    break;
//                }
//                case short_p: {
//                    auto c = static_cast<param<short>*>(arg);
//                    _content = new content<param_*>(new param<short>(c->get_val()));
//                    break;
//                }
//                case integer_p: {
//                    auto c = static_cast<param<int>*>(arg);
//                    _content = new content<param_*>(new param<int>(c->get_val()));
//                    break;
//                }
//                case float_p: {
//                    auto c = static_cast<param<float>*>(arg);
//                    _content = new content<param_*>(new param<float>(c->get_val()));
//                    break;
//                }
//                case double_p: {
//                    auto c = static_cast<param<double>*>(arg);
//                    _content = new content<param_*>(new param<double>(c->get_val()));
//                    break;
//                }
//                case long_p: {
//                    auto c = static_cast<param<long double>*>(arg);
//                    _content = new content<param_*>(new param<long double>(c->get_val()));
//                    break;
//                }
//                default:
//                    throw bad_typeid();
//                    break;
//            }
//            break;
//        }
//        case unary_exp: {
//            auto c = static_cast<content<uexp*>*>(cc._content);
//            auto exp = new uexp();
//            exp->_otype = c->get_arg()->_otype;
//            exp->_son = new constant(c->get_arg()->_son);
//            _content = new content<uexp*>(exp);
//            break;
//        }
//        case binary_exp: {
//            auto c = static_cast<content<bexp*>*>(cc._content);
//            auto exp = new bexp();
//            exp->_otype = c->get_arg()->_otype;
//            exp->_lson = new constant(c->get_arg()->_lson);
//            exp->_rson = new constant(c->get_arg()->_rson);
//            _content = new content<uexp*>(exp);
//            break;
//        }
//        default:
//            throw bad_typeid();
//            break;
//    }
//}
//
//
//constant::constant(bool c){
//    _content = new content<bool>(c);
//}
//
//
//constant::constant(short c){
//    _content = new content<short>(c);
//}
//
//constant::constant(int c){
//    _content = new content<int>(c);
//}
//
//constant::constant(float c){
//    _content = new content<float>(c);
//}
//
//constant::constant(double c){
//    _content = new content<double>(c);
//}
//
//constant::constant(long double c){
//    _content = new content<long double>(c);
//}
//
//constant::constant(const char* s){
//    _content = new content<param_*>(new param<int>(s));
//}
//
//
//constant&  constant::operator=(const constant& c){
//    if (_content) {
//        delete _content;
//    }
//    _content = c._content;
//    return *this;
//}
//
//constant& constant::operator=(bool v) {
//    if (is_binary()) {
//        auto in_c = static_cast<content<bool>*>(_content);
//        in_c->set_val(v);
//        return *this;
//    }
//    if (is_param()) {
//        auto in_c = static_cast<content<param_*>*>(_content);
//        if (in_c->get_arg()->get_type()==binary_p) {
//            auto in_p = static_cast<param<bool>*>(in_c->get_arg());
//            in_p->set_val(v);
//            return *this;
//        }
//    }
//    throw bad_function_call();
//}
//
//
//constant& constant::operator=(int v) {
//    if (is_integer()) {
//        auto in_c = static_cast<content<int>*>(_content);
//        in_c->set_val(v);
//        return *this;
//    }
//    if (is_float()) {
//        auto in_c = static_cast<content<float>*>(_content);
//        in_c->set_val(v);
//        return *this;
//    }
//    if (is_double()) {
//        auto in_c = static_cast<content<double>*>(_content);
//        in_c->set_val(v);
//        return *this;
//    }
//    if (is_long()) {
//        auto in_c = static_cast<content<long double>*>(_content);
//        in_c->set_val(v);
//        return *this;
//    }
//    if (is_param()) {
//        auto in_c = static_cast<content<param_*>*>(_content);
//        if (in_c->get_arg()->get_type()==integer_p) {
//            auto in_p = static_cast<param<int>*>(in_c->get_arg());
//            in_p->set_val(v);
//            return *this;
//        }
//        if (in_c->get_arg()->get_type()==float_p) {
//            auto in_p = static_cast<param<float>*>(in_c->get_arg());
//            in_p->set_val(v);
//            return *this;
//        }
//        if (in_c->get_arg()->get_type()==double_p) {
//            auto in_p = static_cast<param<double>*>(in_c->get_arg());
//            in_p->set_val(v);
//            return *this;
//        }
//        if (in_c->get_arg()->get_type()==long_p) {
//            auto in_p = static_cast<param<long double>*>(in_c->get_arg());
//            in_p->set_val(v);
//            return *this;
//        }
//    }
//    throw bad_function_call();
//}
//
//constant& constant::operator=(float v) {
//    if (is_float()) {
//        auto in_c = static_cast<content<float>*>(_content);
//        in_c->set_val(v);
//        return *this;
//    }
//    if (is_double()) {
//        auto in_c = static_cast<content<double>*>(_content);
//        in_c->set_val(v);
//        return *this;
//    }
//    if (is_long()) {
//        auto in_c = static_cast<content<long double>*>(_content);
//        in_c->set_val(v);
//        return *this;
//    }
//    if (is_param()) {
//        auto in_c = static_cast<content<param_*>*>(_content);
//        if (in_c->get_arg()->get_type()==float_p) {
//            auto in_p = static_cast<param<float>*>(in_c->get_arg());
//            in_p->set_val(v);
//            return *this;
//        }
//        if (in_c->get_arg()->get_type()==double_p) {
//            auto in_p = static_cast<param<double>*>(in_c->get_arg());
//            in_p->set_val(v);
//            return *this;
//        }
//        if (in_c->get_arg()->get_type()==long_p) {
//            auto in_p = static_cast<param<long double>*>(in_c->get_arg());
//            in_p->set_val(v);
//            return *this;
//        }
//    }
//    throw bad_function_call();
//}
//
//constant& constant::operator=(double v) {
//    if (is_double()) {
//        auto in_c = static_cast<content<double>*>(_content);
//        in_c->set_val(v);
//        return *this;
//    }
//    if (is_long()) {
//        auto in_c = static_cast<content<long double>*>(_content);
//        in_c->set_val(v);
//        return *this;
//    }
//    if (is_param()) {
//        auto in_c = static_cast<content<param_*>*>(_content);
//        if (in_c->get_arg()->get_type()==double_p) {
//            auto in_p = static_cast<param<double>*>(in_c->get_arg());
//            in_p->set_val(v);
//            return *this;
//        }
//        if (in_c->get_arg()->get_type()==long_p) {
//            auto in_p = static_cast<param<long double>*>(in_c->get_arg());
//            in_p->set_val(v);
//            return *this;
//        }
//    }
//    throw bad_function_call();
//}
//
//bool constant::contains(const constant* c) const{
//    if (get_type()==unary_exp) {
//        auto cc = static_cast<content<uexp*>*>(_content);
//        auto son = cc ->get_arg()->_son;
//        return (son==c || son->contains(c));
//    }
//    if (get_type()==binary_exp) {
//        auto cc = static_cast<content<bexp*>*>(_content);
//        auto lson = cc ->get_arg()->_lson;
//        auto rson = cc ->get_arg()->_rson;
//        if (lson == c || lson->contains(c)) {
//            return true;
//        }
//        if (rson == c || rson->contains(c)) {
//            return true;
//        }
//    }
//    return false;
//}
//
//bool constant::is_lson(const constant& c) const{
//    if (get_type()==binary_exp) {
//        auto cc = static_cast<content<bexp*>*>(_content);
//        auto lson = cc ->get_arg()->_lson;
//        if (lson == c) {
//            return true;
//        }
//    }
//    return false;
//}
//
//constant* constant::get_lson() const{
//    if (get_type()==binary_exp) {
//        auto cc = static_cast<content<bexp*>*>(_content);
//        return cc ->get_arg()->_lson;
//    }
//    return nullptr;
//}
//
//constant* constant::get_rson() const{
//    if (get_type()==binary_exp) {
//        auto cc = static_cast<content<bexp*>*>(_content);
//        return cc ->get_arg()->_rson;
//    }
//    return nullptr;
//}
//
//
//OperatorType constant::get_otype() const{
//    if (get_type()==unary_exp) {
//        auto cc = static_cast<content<uexp*>*>(_content);
//        return (cc ->get_arg()->_otype);
//    }
//    if (get_type()==binary_exp) {
//        auto cc = static_cast<content<bexp*>*>(_content);
//        return (cc ->get_arg()->_otype);
//    }
//    throw bad_typeid();
//}
//
//
//void constant::reset(){
//    delete _content;
//    _content = nullptr;
//}
//
//bool constant::operator==(const constant &cc) const {
//    if (get_type()!=cc.get_type()) {
//        return false;
//    }
//    if (_content==cc._content){
//        return true;
//    }
//    switch (get_type()) {
//        case binary_c: {
//            return get_val<bool>()==cc->get_val<bool>();
//            break;
//        }
//        case short_c: {
//            return get_val<short>()==cc->get_val<short>();
//            break;
//        }
//        case integer_c: {
//            return get_val<int>()==cc->get_val<int>();
//            break;
//        }
//        case float_c: {
//            return get_val<float>()==cc->get_val<float>();
//            break;
//        }
//        case double_c: {
//            return get_val<double>()==cc->get_val<double>();
//            break;
//        }
//        case long_c: {
//            return get_val<long double>()==cc->get_val<long double>();
//            break;
//        }
//        case parameter: {
//            auto c = static_cast<content<param_*>*>(cc._content);
//            auto this_c = static_cast<content<param_*>*>(_content);
//            return (c->get_arg()->get_name() == this_c->get_arg()->get_name() && c->get_arg()->get_type() == this_c->get_arg()->get_type());
//            break;
//        }
//        case unary_exp: {
//            auto c = static_cast<content<uexp*>*>(cc._content);
//            auto this_c = static_cast<content<uexp*>*>(_content);
//            if (c->get_arg()->_otype!=this_c->get_arg()->_otype) {
//                return false;
//            }
//            if (c->get_arg()->_son != this_c->get_arg()->_son && *c->get_arg()->_son != *this_c->get_arg()->_son) {
//                return false;
//            }
//            return true;
//            break;
//        }
//        case binary_exp: {
//            auto c = static_cast<content<bexp*>*>(cc._content);
//            auto this_c = static_cast<content<uexp*>*>(_content);
//            if (c->get_arg()->_otype!=this_c->get_arg()->_otype) {
//                return false;
//            }
//            if (c->get_arg()->_lson!=this_c->get_arg()->_lson && *c->get_arg()->_lson != *this_c->get_arg()->_lson) {
//                return false;
//            }
//            if (c->get_arg()->_rson!=this_c->get_arg()->_rson && *c->get_arg()->_rson != *this_c->get_arg()->_rson) {
//                return false;
//            }
//            return true;
//            break;
//        }
//        default:
//            throw bad_typeid();
//            break;
//    }
//}
//
//
//
////
////double constant::eval(Function* f, int i) const{
////    if (is_number()) {
////        return _cst;
////    }
////    if (_otype==id_) {
////        return _coef*get_concretization(f, i)->eval(f) + _cst;
////    }
////    
////    switch (_otype) {
////        case plus_:
////            return _coef * (_lson->eval(f,i) + _rson->eval(f,i)) + _cst;
////            break;
////        case minus_:
////            return _coef * (_lson->eval(f,i) - _rson->eval(f,i)) + _cst;
////            break;
////        case product_:
////            return _coef * (_lson->eval(f,i) * _rson->eval(f,i)) + _cst;
////            break;
////        case div_:
////            return _coef * (_lson->eval(f,i) / _rson->eval(f,i)) + _cst;
////            break;
////        case power_:
////            return _coef * pow(_lson->eval(f,i),_rson->eval(f,i)) + _cst;
////            break;
////            
////        default:
////            break;
////    }
////    cerr << "ERROR in Quadratic functions evaluation\n";
////    exit(-1);
////}
////
////double constant::eval(Function* f) const{
////    return eval(f, 0);
////}
////
////
////
////
////bool constant::operator!=(const constant &c) const {
////    return !(*this==c);
////}
////
////bool constant::operator==(double c) const{
////    return (_otype==number_ && _cst==c);
////}
////
////bool constant::factor_of(const constant& c) const{
////    if (_otype==id_) {
////        if (c._otype!=id_ || _name.compare(c._name)!=0) {
////            return false;
////        }
////    }
////    if(_otype != c._otype || _cst != c._cst){
////        return false;
////    }
////    if((!_lson && c._lson) || (!_rson && c._rson)){
////        return false;
////    }
////    if (_lson) {
////        if(!c._lson || *_lson!=*c._lson)
////            return false;
////    }
////    if (_rson) {
////        if(!c._rson || *_rson!=*c._rson)
////            return false;
////    }
////    return true;
////}
////
////bool constant::inverse_factor_of(const constant& c) const{
////    if(_name.compare(c._name)!=0 || _otype != div_ || _otype != c._otype)
////        return false;
////    if((!_lson && c._lson) || (!_rson && c._rson))
////        return false;
////    if (_lson) {
////        if(!c._rson || _lson!=c._rson)
////            return false;
////    }
////    if (_rson) {
////        if(!c._lson || _rson!=c._lson)
////            return false;
////    }
////    return true;
////}
////
////
////void constant::set_coef(double c){
////    if (c==0) {
////        if (_lson &&_lson->is_number()) {
////            delete _lson;
////        }
////        _lson = nullptr;
////        if (_rson && _rson->is_number()) {
////            delete _rson;
////        }
////        _rson = nullptr;
////        _otype = number_;
////        _coef = 0;
////    }
////    else{
////        _coef = c;
////    }
////    
////}
////
////void constant::set_lson(constant *c){
////    if (_lson && _lson->is_number()) {
////        delete _lson;
////    }
////    _lson = c;
////}
////
////void constant::set_rson(constant *c){
////    if (_rson && _rson->is_number()) {
////        delete _rson;
////    }
////    _rson = c;
////}
////
////
////
////constant& constant::operator+=(double v){
////    _cst += v;
////    return *this;
////}
////
////
////constant& constant::operator-=(double v){
////    _cst -= v;
////    return *this;
////}
////
////
////constant& constant::operator*=(double v){
////    if (v==0) {
////        reset();
////        return *this;
////    }
////    _coef *= v;
////    _cst *= v;
////    return *this;
////}
////
////
////constant& constant::operator/=(double v){
////    if(v==0){
////        cerr << "ERROR: Constant" << _name << " dividing by zero!"<< endl;
////        exit(-1);
////    }
////    _coef /= v;
////    _cst /= v;
////    return *this;
////}
////
////
////constant operator+(const constant& c1, const constant &c2) {
////    if (c1.is_number()) {
////        return constant(c2) += c1.get_cst();
////    }
////    if (c2.is_number()) {
////        return constant(c1) += c2.get_cst();
////    }
////    constant res;
////    if ((c1 - c1._cst).factor_of(c2 - c2._cst)) {
////        double coef = c1._coef + c2._coef;
////        if (coef != 0) {
////            res = c1;
////        }
////        res.set_coef(coef);
////        res._cst = c1._cst + c2._cst;
////        return res;
////    }
////    res._coef=1;
////    res._lson = (constant *) &c1;
////    res._otype = plus_;
////    res._rson = (constant *) &c2;
////    return res;
////}
////
////
////constant operator-(const constant& c1, const constant &c2) {
////    if (c1.is_number()) {
////        return constant(-1*c2) += c1.get_cst();
////    }
////    if (c2.is_number()) {
////        return constant(c1) -= c2.get_cst();
////    }
////    constant res;
////    if ((c1 - c1._cst).factor_of(c2 - c2._cst)) {
////        double coef = c1._coef - c2._coef;
////        if (coef != 0) {
////            res = c1;
////        }
////        res.set_coef(coef);
////        res._cst = c1._cst - c2._cst;
////        return res;
////    }
////    res._coef=1;
////    res._lson = (constant *) &c1;
////    res._otype = minus_;
////    res._rson = (constant *) &c2;
////    return res;
////}
////
////
////constant operator*(const constant& c1, const constant &c2) {
////    if (c1.is_number()) {
////        return constant(c2) *= c1.get_cst();
////    }
////    if (c2.is_number()) {
////        return constant(c1) *= c2.get_cst();
////    }
////    if (c1._cst == 0 && c2._cst==0 && c1.inverse_factor_of(c2)) {
////        return constant(c1._coef*c2._coef);
////    }
////    if (c2._cst==0 && c2._cst==0 && c1.power_factor(c2)) {
////        constant res;
////        res._otype = power_;
////        res._coef = c1._coef * c2._coef;
////        if (c1._otype==power_ && c2._otype==power_) {
////            if (c1._lson->_coef > 0 && c2._lson->_coef > 0) {
////                res._lson = c1._lson;
////            }
////            else if (c1._lson->_coef < 0 && c2._lson->_coef < 0){
////                res._lson = c1._lson;                
////            }
////            else if (c1._lson->_coef > 0 && c2._lson->_coef < 0){
////                res._lson = c1._lson;
////                if ((int)c2._rson->_cst%2) { // both odd powers with negative coefficients
////                    res._coef *= -1;
////                }
////                
////            }
////            else if (c1._lson->_coef < 0 && c2._lson->_coef > 0){
////                res._lson = c2._lson;
////                if ((int)c1._rson->_cst%2) { // both odd powers with negative coefficients
////                    res._coef *= -1;
////                }
////                
////            }
////            res._rson = new constant (c1._rson->get_cst() + c2._rson->get_cst());
////        }
////        else if (c1._otype==power_ && c2._otype==id_) {
////            if (c1._lson->_coef >0 && c2._coef > 0) {
////                res._lson = c1._lson;
////            }
////            else if (c1._lson->_coef < 0 && c2._coef < 0){
////                res._lson = c1._lson;
////                res._coef *= -1;
////            }
////            else if (c1._lson->_coef > 0 && c2._coef < 0){
////                res._lson = c1._lson;
////            }
////            else if (c1._lson->_coef < 0 && c2._coef > 0){
////                res._lson = c1._lson;
////                if ((int)(c1._rson->_cst+1)%2) {
////                    res._coef *= -1;
////                }
////            }
////            res._rson = new constant (c1._rson->get_cst() + 1);
////        }
////        else if (c1._otype==id_ && c2._otype==power_) {
////            if (c2._lson->_coef >0 && c1._coef > 0) {
////                res._lson = c2._lson;
////            }
////            else if (c2._lson->_coef < 0 && c1._coef < 0){
////                res._lson = c2._lson;
////                res._coef *= -1;
////            }
////            else if (c2._lson->_coef > 0 && c1._coef < 0){
////                res._lson = c2._lson;
////            }
////            else if (c2._lson->_coef < 0 && c1._coef > 0){
////                res._lson = c1._lson;
////                if ((int)(c1._rson->_cst+1)%2) {
////                    res._coef *= -1;
////                }
////            }
////            res._rson = new constant (c2._rson->get_cst() + 1);
////        }
////        else if (c1._otype==id_ && c2._otype==id_) {
////            if (c2._coef >0 && c1._coef > 0) {
////                res._lson = (constant*)&c1;
////            }
////            else if (c2._coef < 0 && c1._coef < 0){
////                res._lson = (constant*)&c1;
////                res._coef *= -1;
////            }
////            else if (c2._coef > 0 && c1._coef < 0){
////                res._lson = (constant*)&c2;
////            }
////            else if (c2._lson->_coef < 0 && c1._coef > 0){
////                res._lson = (constant*)&c1;
////            }
////            res._rson = new constant (2);
////        }
////        return res;
////    }
////    constant res;
////    res._coef=1;
////    res._lson = (constant *) &c1;
////    res._otype = product_;
////    res._rson = (constant *) &c2;
////    return res;
////}
////
////bool constant::power_factor(const constant& c) const{
////    if (_otype == id_ && c._otype == id_) {
////        return (this->factor_of(c));
////    }
////    if (c._otype == power_ && _otype== power_) {
////        return ((*c._lson == *_lson) || (*c._lson==-1*(*_lson)));
////    }
////    if (c._otype == power_ && _otype == id_) {
////        return (c._lson->factor_of(*this) && (c._lson->_coef==1 || c._lson->_coef==-1));
////    }
////    if (_otype == power_ && c._otype == id_) {
////        return (_lson->factor_of(c) && (_lson->_coef==1 || _lson->_coef==-1));
////    }
////    return false;
////}
////
////constant operator/(const constant& c1, const constant &c2) {
////    if (c2.is_number()) {
////        return constant(c1) /= c2.get_cst();
////    }
////    if (c1.is_number() && c1._cst == 0) {
////        return constant(c1);
////    }
////    if (c1._cst == 0 && c2._cst == 0 && c1.factor_of(c2)) {
////        return constant(c1._coef/c2._coef);
////    }
////    
////    if (c1.is_number()){
////        return c1.get_cst()/c2;
////    }
////    constant res;
////    res._coef=1;
////    res._lson = (constant *) &c1;
////    res._otype = div_;
////    res._rson = (constant *) &c2;
////    return res;
////}
////
////
////constant operator^(const constant& c, int p) {
////    constant res;
////    res._coef=1;
////    res._lson = (constant*)&c;
////    res._otype = power_;
////    res._rson = new constant(p);
////    return res;
////}
////
////
////constant operator^(const constant& c1, const constant& c2) {
////    constant res;
////    res._coef=1;
////    res._lson = (constant *) &c1;
////    res._otype = power_;
////    res._rson = (constant *) &c2;
////    return res;
////}
////
////constant cos(constant &c) {
////    constant res;
////    res._coef=1;
////    res._lson = (constant *) &c;
////    res._otype = cos_;
////    res._rson = nullptr;
////    return res;
////}
////
////constant sin(constant &c) {
////    constant res;
////    res._coef=1;
////    res._lson = (constant *) &c;
////    res._otype = sin_;
////    res._rson = nullptr;
////    return res;
////}
////
////constant sqrt(constant &c) {
////    constant res;
////    res._coef=1;
////    res._lson = (constant *) &c;
////    res._otype = sqrt_;
////    res._rson = nullptr;
////    return res;
////}
////
////constant expo(constant &c) {
////    constant res;
////    res._coef=1;
////    res._lson = (constant *) &c;
////    res._otype = exp_;
////    res._rson = nullptr;
////    return res;
////}
////
////constant log(constant &c) {
////    constant res;
////    res._coef=1;
////    res._lson = (constant *) &c;
////    res._otype = log_;
////    res._rson = nullptr;
////    return res;
////    
////}
////
////
////constant operator*(double cst, const constant &c) {
////    constant res(c);
////    return move(res *= cst);
////}
////
////constant operator*(const constant &c, double cst) {
////    constant res(c);
////    return move(res *= cst);
////}
////
////
////
////constant operator+(double cst, const constant &c) {
////    constant res(c);
////    return move(res += cst);
////}
////
////constant operator+(const constant &c, double cst) {
////    constant res(c);
////    return move(res += cst);
////}
////
////
////
////constant operator-(double cst, const constant &c) {
////    constant res(c);
////    return move((res *= -1) += cst);
////}
////
////constant operator-(const constant &c, double cst) {
////    constant res(c);
////    return move(res -= cst);
////}
////
////
////constant operator/(double cst, const constant &c) {
////    constant res;
////    res._coef=1;
////    res._lson = new constant(cst);
////    res._otype = div_;
////    res._rson = (constant*)&c;
////    return res;
////}
////
////constant operator/(const constant &c, double cst) {
////    constant res(c);
////    return move(res /= cst);
////}
////
////
////
////constant operator*(double cst, constant &&c) {
////    return move(c *= cst);
////}
////
////constant operator*(constant &&c, double cst) {
////    return move(c *= cst);
////}
////
////constant operator+(double cst, constant &&c) {
////    return move(c += cst);
////}
////
////constant operator+(constant &&c, double cst) {
////    return move(c += cst);
////}
////
////
////constant operator-(double cst, constant &&c) {
////    return move((c *= -1) += cst);
////}
////
////constant operator-(constant &&c, double cst) {
////    return move(c -= cst);
////}
////
////
////constant operator/(constant &&c, double cst) {
////    return move(c /= cst);
////}
////
////
////
////string constant::to_string() const {
////    
////    string s = _name;
////    if(_otype==id_){
////        s  +=  "[";
////        if(_conc) {
////            for (unsigned int i = 0; i < _conc->size(); i++) {
////                ostringstream strs;
////                strs << _conc->at(i);
////                s  +=  strs.str();
////                if (i < _conc->size()-1)
////                    s  +=  " ,";
////            }
////        }
////        s  +=  "]";
////        return s;
////    }
////    
////    
////    if(!_rson){
////        if(_otype==cos_) {
////            s  +=  "cos";
////        }
////        if(_otype==sin_) {
////            s  +=  "sin";
////        }
////        if (_otype==exp_) {
////            s  +=  "exp";
////        }
////        if (_otype==log_) {
////            s  +=  "log";
////        }
////        if (_otype==sqrt_) {
////            s  +=  "sqrt";
////        }
////        s  +=  "(";
////        s += _lson->to_string();
////        s  +=  ")";
////        s  +=  "\n";
////        return s;
////        
////    }
////    
////    if(_otype!=plus_ && _otype!=minus_ && _lson->_otype!=id_) {
////        s  +=  "(";
////        s += _lson->to_string();
////        s  +=  ")";
////    }
////    else
////        s += _lson->to_string();
////    
////    if (_otype==plus_) {
////        s  +=  " + ";
////    }
////    if (_otype==minus_) {
////        s  +=  " - ";
////    }
////    if (_otype==product_) {
////        if(!_lson || _lson->_otype!=id_)
////            s  +=  " * ";
////    }
////    if (_otype==div_) {
////        s  +=  "/";
////    }
////    
////    if (_otype==power_) {
////        s  +=  "^";
////    }
////    
////    if (_otype==plus_ || _rson->_otype==id_ || _rson->_otype==number_) {
////        s += _rson->to_string();
////    }
////    else {
////        s  +=  "(";
////        s += _rson->to_string();
////        s  +=  ")";
////    }
////    
////    return s;
////    
////}
////
////void constant::print_tree() const {
////    if (_otype==number_) {
////        cout << _cst;
////        return;
////    }
////    
////    if(_otype==id_){
////        if (_coef != 1 && _coef!=-1) {
////            cout << _coef;
////        }
////        if (_coef==-1){
////            cout << "-";
////        }
////        cout << _name;
////        if (_cst !=0) {
////            if (_cst > 0 && _coef!=1) {
////                cout << "+";
////            }
////            cout << _cst;
////        }
////        return;
////    }
////    
////    if (_cst !=0) {
////        cout << "(";
////    }
////    if (_coef != 1 && _coef!=-1) {
////        cout << _coef << "(";
////    }
////    if (_coef==-1){
////        cout << "-(";
////    }
////    
////    if(!_rson){
////        if(_otype==cos_) {
////            cout << "cos";
////        }
////        if(_otype==sin_) {
////            cout << "sin";
////        }
////        if (_otype==exp_) {
////            cout << "exp";
////        }
////        if (_otype==log_) {
////            cout << "log";
////        }
////        if (_otype==sqrt_) {
////            cout << "sqrt";
////        }
////        cout << "(";
////        _lson->print_tree();
////        cout << ")";
////        if (_coef!=1) {
////            cout << ")";
////        }
////        if (_cst !=0) {
////            if (_cst > 0) {
////                cout << "+";
////            }
////            cout << _cst << ")";
////        }
////        return;
////        
////    }
////    
////    if((_otype == power_) || (_otype!=plus_ && _otype!=minus_ && _lson->_otype!=id_ && _lson->_otype!=number_)) {
////        cout << "(";
////        _lson->print_tree();
////        cout << ")";
////    }
////    else
////        _lson->print_tree();
////    
////    if (_otype==plus_) {
////        cout << " + ";
////    }
////    if (_otype==minus_) {
////        cout << " - ";
////    }
////    if (_otype==product_) {
////        if(_lson->_otype!=id_ || _rson->_coef < 0 || (_rson->_coef==0 && _rson->_cst < 0))
////            cout << " * ";
////    }
////    if (_otype==div_) {
////        cout << "/";
////    }
////    
////    if (_otype==power_) {
////        cout << "^";
////    }
////    
////    if (_otype==plus_ || _rson->_otype==id_ || _rson->_otype==number_) {
////        _rson->print_tree();
////    }
////    else {
////        cout << "(";
////        _rson->print_tree();
////        cout << ")";
////    }
////    if (_coef!=1) {
////        cout << ")";
////    }
////    if (_cst !=0) {
////        if (_cst > 0) {
////            cout << "+";
////        }
////        cout << _cst << ")";
////    }
////    
////}
////
////
////void constant::print() const {
////    
////    
////    if(_otype==id_){
////        if (_coef != 1 && _coef!=-1) {
////            cout << _coef;
////        }
////        if (_coef==-1){
////            cout << "-";
////        }
////        cout << _name;
////        if (_cst >0) {
////            cout << "+" << _cst;
////        }
////        if (_cst <0) {
////            cout << _cst;
////        }
////    }
////    else {
////        cout << _name << " = ";
////        print_tree();
////    }
////    if (_conc) {
////        cout << " = ";
////        cout << "[";
////        for (unsigned int i= 0; i< _conc->size()-1; i++) {
////            _conc->at(i)->print_();
////            cout << ",";
////        }
////        _conc->at(_conc->size()-1)->print_();
////        cout << "];";
////    }
////    cout << endl;
////}
////
////
////void constant::print_() const {
////    
////    if (is_number()) {
////        cout << _cst;
////        return;
////    }
////    if(_otype==id_){
////        if (_coef != 1 && _coef!=-1) {
////            cout << _coef;
////        }
////        if (_coef==-1){
////            cout << "-";
////        }
////        cout << _name;
////        if (_cst >0) {
////            cout << "+" << _cst;
////        }
////        if (_cst <0) {
////            cout << _cst;
////        }
////    }
////    else {
////        print_tree();
////    }
////}
//
////
////
////
////unsigned long constant::assign(double v) {
////    if (!is_letter()) {
////        cerr << "ERROR: Cannot assign value to non-atomic constants!\n";
////        exit(-1);
////    }
////    if (!_conc) {
////        _conc = shared_ptr<vector<constant*>>(new vector<constant*>(1));
////        _conc->at(0) = new constant(v);
////        return 0;
////    }
////    else {
////        _conc->push_back(new constant(v));
////        return _conc->size()-1;
////    }
////}
////
////
////
////constant* constant::get_concretization(Function *f, int i) const{
////    if (!_conc) {
////        throw bad_alloc();
////    }
////    for (auto & p: *_conc){
////        if (p.first==f) {
////            return p.second.at(i);
////        }
////    }
////    throw invalid_argument("no such concretization");
////}
////
////void constant::concretize(Function *f, constant* c){
////    vector<constant*> vec;
////    if (!_conc) {
////        _conc = shared_ptr<forward_list<pair<Function*, vector<constant*>>>>(new forward_list<pair<Function*, vector<constant*>>>());
////        vec.push_back(c);
////        _conc->push_front(make_pair(f,vec));
////    }
////    else {
////        for (auto& p: *_conc){
////            if (p.first==f) {
////                p.second.push_back(c);
////                return;
////            }
////        }
////        vec.push_back(c);
////        _conc->push_front(make_pair(f,vec));
////    }
////}