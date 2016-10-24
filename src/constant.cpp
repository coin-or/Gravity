//
// Created by Hassan on 19/11/2015.
//

#include <math.h>
#include <Gravity/var.h>
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
        case par_:{
            return (c1->is_param() && *(param_ *)c1 == *(param_ *)c2);
            break;
        }
        case uexp_: {
            return (c1->is_uexpr() && *(uexpr *)c1 == *(uexpr *)c2);
            break;
        }
        case bexp_: {
            return (c1->is_bexpr() && *(bexpr *)c1 == *(bexpr *)c2);
            break;
        }
        default:
            break;
    }
    return false;
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
        case par_:{
            auto p_c2 = (param_*)(c2);
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
        case uexp_: {
            return new uexpr(*(uexpr*)c2);
            break;
        }
        case bexp_: {
            return new bexpr(*(bexpr*)c2);
            break;
        }
        case var_:{
            auto p_c2 = (param_*)(c2);
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
    
        default:
            break;
    }
    return nullptr;
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
        case par_:{
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
        case uexp_: {
            return ((uexpr*)c)->eval(i);
            break;
        }
        case bexp_: {
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
        case par_:{
            auto p_c2 = (param_*)(c2);
            switch (p_c2->get_intype()) {
                case binary_:
                    ((param<bool>*)p_c2)->print();
                    break;
                case short_:
                    ((param<short>*)p_c2)->print();
                    break;
                case integer_:
                    ((param<int>*)p_c2)->print();
                    break;
                case float_:
                    ((param<float>*)p_c2)->print();
                    break;
                case double_:
                    ((param<double>*)p_c2)->print();
                    break;
                case long_:
                    ((param<long double>*)p_c2)->print();
                    break;
                default:
                    break;
            }
            break;
        }
        case uexp_: {
            ((uexpr*)c2)->print(false);
            break;
        }
        case bexp_: {
            ((bexpr*)c2)->print(false);
            break;
        }
        case var_: {
            auto p_c2 = (param_*)(c2);
            switch (p_c2->get_intype()) {
                case binary_:
                    ((var<bool>*)p_c2)->print();
                    break;
                case short_:
                    ((var<short>*)p_c2)->print();
                    break;
                case integer_:
                    ((var<int>*)p_c2)->print();
                    break;
                case float_:
                    ((var<float>*)p_c2)->print();
                    break;
                case double_:
                    ((var<double>*)p_c2)->print();
                    break;
                case long_:
                    ((var<long double>*)p_c2)->print();
                    break;
                default:
                    break;
            }
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
    _type = uexp_;
}

uexpr::uexpr(const uexpr& exp){
    _otype = exp._otype;
    _son = copy(exp._son);
    _type = uexp_;
}

uexpr::uexpr(uexpr&& exp){
    _otype = exp._otype;
    _son = move(exp._son);
    exp._son = nullptr;
    _type = uexp_;
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
    
    if (_son->get_type()==uexp_) {
        return ((uexpr*)_son)->contains(c);
    }
    if (_son->get_type()==bexp_) {
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
        if (_lson->get_type()==uexp_ && ((uexpr*)_lson)->contains(c)) {
            return true;
        }
        if (_lson->get_type()==bexp_ && ((bexpr*)_lson)->contains(c)) {
            return true;
        }
    }
    if(_rson) {
        if (equals(_rson, c)) {
            return true;
        }
        if (_rson->get_type()==uexp_ && ((uexpr*)_rson)->contains(c)) {
            return true;
        }
        if (_rson->get_type()==bexp_ && ((bexpr*)_rson)->contains(c)) {
            return true;
        }
    }
    return false;
};

bexpr::bexpr(){
    _otype = id_;
    _lson = nullptr;
    _rson = nullptr;
    _type = bexp_;
}

bexpr::bexpr(const bexpr& exp){ /**< Copy constructor from binary expression tree */
    _otype = exp._otype;
    _lson = copy(exp._lson);
    _rson =  copy(exp._rson);
    _type = bexp_;
};

bexpr::bexpr(bexpr&& exp){ /**< Move constructor from binary expression tree */
    _otype = exp._otype;
    _lson = move(exp._lson);
    _rson = move(exp._rson);
    _type = bexp_;
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
    if((_otype==product_ || _otype==div_) && (_lson->get_type()==uexp_ || _lson->get_type()==bexp_)) {
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
    
    if (_otype==plus_ || (_rson->get_type()!=uexp_ && _rson->get_type()!=bexp_)) {
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
