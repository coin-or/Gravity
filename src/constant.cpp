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



constant_* copy(constant_* c2){/**< Copy c2 into c1 detecting the right class, i.e., constant<>, param<>, uexpr or bexpr. */
    
    if (!c2) {
        return nullptr;
    }
    switch (c2->get_type()) {
        case binary_c: {
        return new constant<bool>(((constant<bool>*)(c2))->eval());
            break;
        }
        case short_c: {
            return new constant<short>(((constant<short>*)(c2))->eval());
            break;
        }
        case integer_c: {
            return new constant<int>(((constant<int>*)(c2))->eval());
            break;
        }
        case float_c: {
            return new constant<float>(((constant<float>*)(c2))->eval());
            break;
        }
        case double_c: {
            return new constant<double>(((constant<double>*)(c2))->eval());
            break;
        }
        case long_c: {
            return new constant<long double>(((constant<long double>*)(c2))->eval());
            break;
        }
        case parameter:{
            auto p_c2 = (param_*)(c2);
            switch (p_c2->get_intype()) {
                case binary_p:
                    return new param<bool>(*(param<bool>*)p_c2);
                    break;
                case short_p:
                    return new param<short>(*(param<short>*)p_c2);
                    break;
                case integer_p:
                    return new param<int>(*(param<int>*)p_c2);
                    break;
                case float_p:
                    return new param<float>(*(param<float>*)p_c2);
                    break;
                case double_p:
                    return new param<double>(*(param<double>*)p_c2);
                    break;
                case long_p:
                    return new param<long double>(*(param<long double>*)p_c2);
                    break;
                default:
                    break;
            }
            break;
        }
        case unary_exp: {
            return new uexpr(*(uexpr*)c2);
            break;
        }
        case binary_exp: {
            return new bexpr(*(bexpr*)c2);
            break;
        }
        default:
            break;
    }
    
}

double eval(int i, const constant_* c){
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
        case parameter:{
            auto p_c2 = (param_*)(c);
            switch (p_c2->get_intype()) {
                case binary_p:
                    return ((param<bool>*)p_c2)->eval(i);
                    break;
                case short_p:
                    return ((param<short>*)p_c2)->eval(i);
                    break;
                case integer_p:
                    return ((param<int>*)p_c2)->eval(i);
                    break;
                case float_p:
                    return ((param<float>*)p_c2)->eval(i);
                    break;
                case double_p:
                    return ((param<double>*)p_c2)->eval(i);
                    break;
                case long_p:
                    return (double)((param<long double>*)p_c2)->eval(i);
                    break;
                default:
                    break;
            }
            break;
        }
        case unary_exp: {
            return ((uexpr*)c)->eval(i);
            break;
        }
        case binary_exp: {
            return ((bexpr*)c)->eval(i);
            break;
        }
        default:
            break;
    }
    return 0;
}


void poly_print(const constant_* c2){/**< Copy c2 into c1 detecting the right class, i.e., constant<>, param<>, uexpr or bexpr. */
    
    if (!c2) {
        return;
    }
    switch (c2->get_type()) {
        case binary_c: {
            ((constant<bool>*)(c2))->print();
            break;
        }
        case short_c: {
            ((constant<short>*)(c2))->print();
            break;
        }
        case integer_c: {
            ((constant<int>*)(c2))->print();
            break;
        }
        case float_c: {
            ((constant<float>*)(c2))->print();
            break;
        }
        case double_c: {
            ((constant<double>*)(c2))->print();
            break;
        }
        case long_c: {
            ((constant<long double>*)(c2))->print();
            break;
        }
        case parameter:{
            auto p_c2 = (param_*)(c2);
            switch (p_c2->get_intype()) {
                case binary_p:
                    ((param<bool>*)p_c2)->print();
                    break;
                case short_p:
                    ((param<short>*)p_c2)->print();
                    break;
                case integer_p:
                    ((param<int>*)p_c2)->print();
                    break;
                case float_p:
                    ((param<float>*)p_c2)->print();
                    break;
                case double_p:
                    ((param<double>*)p_c2)->print();
                    break;
                case long_p:
                    ((param<long double>*)p_c2)->print();
                    break;
                default:
                    break;
            }
            break;
        }
        case unary_exp: {
            ((uexpr*)c2)->print(false);
            break;
        }
        case binary_exp: {
            ((bexpr*)c2)->print(false);
            break;
        }
        default:
            break;
    }
    
}


/* UNARY EXPRESSIONS */

uexpr::uexpr(){
    _otype = id_;
    _son = nullptr;
    _type = unary_exp;
}

uexpr::uexpr(const uexpr& exp){
    _otype = exp._otype;
    _son = copy(exp._son);
    _type = unary_exp;
}

uexpr::uexpr(uexpr&& exp){
    _otype = exp._otype;
    _son = move(exp._son);
    exp._son = nullptr;
    _type = unary_exp;
}

uexpr& uexpr::operator=(const uexpr& exp){
    delete _son;
    _son = copy(exp._son);
    _otype = exp._otype;
    return *this;
}

uexpr& uexpr::operator=(uexpr&& exp){
    delete _son;
    _son = move(exp._son);
    exp._son = nullptr;
    _otype = exp._otype;
    return *this;
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
    res._son = copy((constant_*)&c);
    return res;
};


uexpr sin(const constant_& c){
    uexpr res;
    res._otype = sin_;
    res._son = copy((constant_*)&c);
    return res;
};


uexpr sqrt(const constant_& c){
    uexpr res;
    res._otype = sqrt_;
    res._son = copy((constant_*)&c);
    return res;
};

uexpr expo(const constant_& c){
    uexpr res;
    res._otype = exp_;
    res._son = copy((constant_*)&c);
    return res;
};

uexpr log(const constant_& c){
    uexpr res;
    res._otype = log_;
    res._son = copy((constant_*)&c);
    return res;
};

double uexpr::eval(unsigned int i) const{
    if (!_son) {
        throw invalid_argument("Cannot evaluate empty expression!");
    }
    switch (_otype) {
        case cos_:
            return cos(::eval(i,_son));
            break;
        case sin_:
            return sin(::eval(i,_son));
            break;
        case sqrt_:
            return sqrt(::eval(i,_son));
            break;
        case log_:
            return log(::eval(i,_son));
            break;
        case exp_:
            return exp(::eval(i,_son));
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

bexpr::bexpr(){
    _otype = id_;
    _lson = nullptr;
    _rson = nullptr;
    _type = binary_exp;
}

bexpr::bexpr(const bexpr& exp){ /**< Copy constructor from binary expression tree */
    _otype = exp._otype;
    _lson = copy(exp._lson);
    _rson =  copy(exp._rson);
    _type = binary_exp;
};

bexpr::bexpr(bexpr&& exp){ /**< Move constructor from binary expression tree */
    _otype = exp._otype;
    _lson = move(exp._lson);
    _rson = move(exp._rson);
    _type = binary_exp;
};

bexpr& bexpr::operator=(const bexpr& exp){
    delete _lson;
    delete _rson;
    _lson = copy(exp._lson);
    _rson =  copy(exp._rson);
    _otype = exp._otype;
    return *this;
}

bexpr& bexpr::operator=(bexpr&& exp){
    delete _lson;
    delete _rson;
    _lson = move(exp._lson);
    _rson = move(exp._rson);
    _otype = exp._otype;
    return *this;
}



template<typename other_type> bexpr& bexpr::operator+=(const other_type& v){
    bexpr res;
    constant_* old = _lson;
    _lson = copy(this);
    delete old;
    if (is_arithmetic<other_type>::value) {
        _rson = new constant<other_type>(v);
    }
    else {
        delete _rson;
        _rson = copy((constant_*)&v);
    }
    _otype = plus_;
    return *this;
}


template<typename other_type> bexpr& bexpr::operator-=(const other_type& v){
    constant_* old = _lson;
    _lson = copy(this);
    delete old;
    if (is_arithmetic<other_type>::value) {
        _rson = new constant<other_type>(v);
    }
    else {
        delete _rson;
        _rson = copy((constant_*)&v);
    }
    _otype = minus_;
    return *this;
}


template<typename other_type> bexpr& bexpr::operator*=(const other_type& v){
    constant_* old = _lson;
    _lson = copy(this);
    delete old;
    if (is_arithmetic<other_type>::value) {
        _rson = new constant<other_type>(v);
    }
    else {
        delete _rson;
        _rson = copy((constant_*)&v);
    }
    
    _otype = product_;
    return *this;
}

template<typename other_type> bexpr& bexpr::operator/=(const other_type& v){
    constant_* old = _lson;
    _lson = copy(this);
    delete old;
    if (is_arithmetic<other_type>::value) {
        _rson = new constant<other_type>(v);
    }
    else {
        delete _rson;
        _rson = copy((constant_*)&v);
    }    
    _otype = div_;
    return *this;
}





//bexpr operator+(const constant_& c1, const constant_& c2){
//    bexpr res;
//    res._otype = plus_;
//    copy(res._lson,(constant_*)&c1);
//    copy(res._rson, (constant_*)&c2);
//    return res;
//}
//
//
//bexpr operator-(const constant_& c1, const constant_& c2){
//    bexpr res;
//    res._otype = minus_;
//    copy(res._lson,(constant_*)&c1);
//    copy(res._rson, (constant_*)&c2);
//    return res;
//}
//
//bexpr operator*(const constant_& c1, const constant_& c2){
//    bexpr res;
//    res._otype = product_;
//    copy(res._lson,(constant_*)&c1);
//    copy(res._rson, (constant_*)&c2);
//    return res;
//}
//
//bexpr operator/(const constant_& c1, const constant_& c2){
//    bexpr res;
//    res._otype = div_;
//    copy(res._lson,(constant_*)&c1);
//    copy(res._rson, (constant_*)&c2);
//    return res;
//}







double bexpr::eval(unsigned int i) const{
    if (!_lson || !_rson) {
        throw invalid_argument("Cannot evaluate empty expression!");
    }
    switch (_otype) {
        case plus_:
            return ::eval(i,_lson) + ::eval(i,_rson);
            break;
        case minus_:
            return ::eval(i,_lson) + ::eval(i,_rson);
            break;
        case product_:
            return ::eval(i,_lson) * ::eval(i,_rson);
            break;
        case div_:
            return ::eval(i,_lson)/::eval(i,_rson);
            break;
        case power_:
            return powl(::eval(i,_lson),::eval(i,_rson));
            break;
        default:
            throw invalid_argument("Unsupported binary operator");
            break;
    }
    
}

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


//    if(_otype==cos_) {
//        cout << "cos";
//        }
//        if(_otype==sin_) {
//            cout << "sin";
//        }
//        if (_otype==exp_) {
//            cout << "exp";
//        }
//        if (_otype==log_) {
//            cout << "log";
//        }
//        if (_otype==sqrt_) {
//            cout << "sqrt";
//        }
//        cout << "(";
//        _lson->print_tree();
//        cout << ")";
//        if (_coef!=1) {
//            cout << ")";
//        }
//        if (_cst !=0) {
//            if (_cst > 0) {
//                cout << "+";
//            }
//            cout << _cst << ")";
//        }
//        return;
//        
//        }

void uexpr::print(bool endline) const{
    switch (_otype) {
        case log_:
            cout << "log(";
            poly_print(_son);
            cout << ")";
            break;
            
        case exp_:
            cout << "exp(";
            poly_print(_son);
            cout << ")";
            break;
            
        case cos_:
            cout << "cos(";
            poly_print(_son);
            cout << ")";
            break;
            
        case sin_:
            cout << "sin(";
            poly_print(_son);
            cout << ")";
            break;
            
        case sqrt_:
            cout << "sqrt(";
            poly_print(_son);
            cout << ")";
            break;
        default:
            break;
    }
    if(endline)
        cout << endl;
}

void bexpr::print(bool endline) const {
    if((_otype==product_ || _otype==div_) && (_lson->get_type()==unary_exp || _lson->get_type()==binary_exp)) {
        cout << "(";
        poly_print(_lson);
        cout << ")";
    }
    else
        poly_print(_lson);

    if (_otype==plus_) {
        cout << " + ";
    }
    if (_otype==minus_) {
        cout << " - ";
    }
    if (_otype==product_) {
        cout << " * ";
    }
    if (_otype==div_) {
        cout << "/";
    }

    if (_otype==power_) {
        cout << "^";
    }
    
    if (_otype==plus_ || (_rson->get_type()!=unary_exp && _rson->get_type()!=binary_exp)) {
        poly_print(_rson);
    }
    else {
        cout << "(";
        poly_print(_rson);
        cout << ")";
    }
    if(endline)
        cout << endl;
}

//    
//    if(_otype==id_){
//        if (_coef != 1 && _coef!=-1) {
//            cout << _coef;
//        }
//        if (_coef==-1){
//            cout << "-";
//        }
//        cout << _name;
//        if (_cst >0) {
//            cout << "+" << _cst;
//        }
//        if (_cst <0) {
//            cout << _cst;
//        }
//    }
//    else {
//        cout << _name << " = ";
//        print_tree();
//    }
//    if (_conc) {
//        cout << " = ";
//        cout << "[";
//        for (unsigned int i= 0; i< _conc->size()-1; i++) {
//            _conc->at(i)->print_();
//            cout << ",";
//        }
//        _conc->at(_conc->size()-1)->print_();
//        cout << "];";
//    }
//    cout << endl;
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
