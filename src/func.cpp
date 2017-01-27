//
//  func.cpp
//  Gravity
//
//  Created by Hijazi, Hassan on 24/10/16.
//
//
#include <math.h>
#include <Gravity/func.h>

constant_* copy(const constant_* c2){/**< Copy c2 into c1 detecting the right class, i.e., constant<>, param<>, uexpr or bexpr. */
    
    if (!c2) {
        return nullptr;
    }
    switch (c2->get_type()) {
        case binary_c: {
            return new constant<bool>(*(constant<bool>*)(c2));
            break;
        }
        case short_c: {
            return new constant<short>(*(constant<short>*)(c2));
            break;
        }
        case integer_c: {
            return new constant<int>(*(constant<int>*)(c2));
            break;
        }
        case float_c: {
            return new constant<float>(*(constant<float>*)(c2));
            break;
        }
        case double_c: {
            return new constant<double>(*(constant<double>*)(c2));
            break;
        }
        case long_c: {
            return new constant<long double>(*(constant<long double>*)(c2));
            break;
        }
        case par_c:{
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
        case uexp_c: {
            return new uexpr(*(uexpr*)c2);
            break;
        }
        case bexp_c: {
            return new bexpr(*(bexpr*)c2);
            break;
        }
        case var_c:{
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
            
        case func_c: {
            return new func_(*c2);
            break;
        }
            
        default:
            break;
    }
    return nullptr;
}



func_::func_(){
    set_type(func_c);
    _params = new map<string, pair<param_*, int>>();
    _vars = new map<string, pair<param_*, int>>();
    _lterms = new map<string, lterm>();
    _qterms = new map<string, qterm>();;
    _pterms = new map<string, pterm>();;
    _cst = new constant<int>(0);
};

func_::func_(const constant_& c){
    set_type(func_c);
    _params = new map<string, pair<param_*, int>>();
    _vars = new map<string, pair<param_*, int>>();
    _lterms = new map<string, lterm>();
    _qterms = new map<string, qterm>();
    _pterms = new map<string, pterm>();
    switch (c.get_type()) {
        case binary_c: {
            _cst = new constant<bool>(*(constant<bool>*)(&c));
            break;
        }
        case short_c: {
            _cst = new constant<short>(*(constant<short>*)(&c));
            break;
        }
        case integer_c: {
            _cst = new constant<int>(*(constant<int>*)(&c));
            break;
        }
        case float_c: {
            _cst = new constant<float>(*(constant<float>*)(&c));
            break;
        }
        case double_c: {
            _cst = new constant<double>(*(constant<double>*)(&c));
            break;
        }
        case long_c: {
            _cst = new constant<long double>(*(constant<long double>*)(&c));
            break;
        }
        case par_c:{
            auto p_c2 = (param_*)copy(&c);
            _lterms->insert(make_pair<>(p_c2->get_name(), p_c2));
            add_param(p_c2);
            _cst = new constant<int>(0);
            break;
        }
        case uexp_c: {
//            return new uexpr(*(uexpr*)(&c2));
            break;
        }
        case bexp_c: {
//            return new bexpr(*(bexpr*)&c2);
            break;
        }
        case var_c:{
            auto p_c2 = (param_*)copy(&c);
            _lterms->insert(make_pair<>(p_c2->get_name(), p_c2));
            add_var(p_c2);
            _ftype = lin_;
            _cst = new constant<int>(0);
            break;
        }
        case func_c: {
            auto f = (func_*)&c;
            _ftype = f->_ftype;
            _return_type = f->_return_type;
            _convex = f->_convex;
            for (auto &pair:*f->_lterms) {
                insert(pair.second);
            }
            for (auto &pair:*f->_qterms) {
                insert(pair.second);
            }
            _cst = copy(f->_cst);
            break;
        }
        default:
            break;
    }
}

func_::func_(func_&& f){
    set_type(func_c);
    _ftype = f._ftype;
    _return_type = f._return_type;
    _convex = f._convex;
    _lterms = f._lterms;
    f._lterms = nullptr;
    _qterms = f._qterms;
    f._qterms = nullptr;
    _pterms = f._pterms;
    f._pterms = nullptr;
    _vars = f._vars;
    f._vars = nullptr;
    _params = f._params;
    f._vars = nullptr;
    _cst = f._cst;
    f._cst = nullptr;
}

func_::func_(const func_& f){
    set_type(func_c);
    _params = new map<string, pair<param_*, int>>();
    _vars = new map<string, pair<param_*, int>>();
    _lterms = new map<string, lterm>();
    _qterms = new map<string, qterm>();
    _pterms = new map<string, pterm>();
    _ftype = f._ftype;
    _return_type = f._return_type;
    _convex = f._convex;
    for (auto &pair:*f._lterms) {
        insert(pair.second);
    }
    for (auto &pair:*f._qterms) {
        insert(pair.second);
    }
    _cst = copy(f._cst);
}


func_::~func_(){
    if (_vars) {
        for (auto &elem: *_vars) {
            delete elem.second.first;
        }
        delete _vars;
    }
    if (_params) {
        for (auto &elem: *_params) {
            delete elem.second.first;
        }
        delete _params;
    }
    delete _lterms;
    delete _qterms;
    delete _pterms;
    delete _cst;
};

bool constant_::is_zero() const{ /**< Returns true if constant equals 0 */
    return (is_number() && eval(this)==0);
}

func_& func_::operator+=(double v){
    add(_cst, constant<double>(v));
    return *this;
}

func_& func_::operator-=(double v){
    substract(_cst, constant<double>(v));
    return *this;
}

func_& func_::operator*=(double v){
    multiply(_cst, constant<double>(v));
    return *this;
}

func_& func_::operator/=(double v){
    divide(_cst, constant<double>(v));
    return *this;
}


func_& func_::operator+=(const func_& f){
    if (!is_constant() && f.is_constant()) {
        _cst = add(_cst, f);
        return *this;
    }
    _cst = add(_cst, *f._cst);
    for (auto &pair:*f._lterms) {
        this->insert(pair.second._sign, *pair.second._coef, *pair.second._p);
    }
    for (auto &pair:*f._qterms) {
        this->insert(pair.second._sign, *pair.second._coef, *pair.second._p->first, *pair.second._p->second);
    }
//    update_convexity();
    return *this;
}

func_& func_::operator-=(const func_& f){
    if (!is_constant() && f.is_constant()) {
        _cst = substract(_cst, f);
        return *this;
    }
    _cst = substract(_cst, *f._cst);
    for (auto &pair:*f._lterms) {
        this->insert(!pair.second._sign, *pair.second._coef, *pair.second._p);
    }
    for (auto &pair:*f._qterms) {
        this->insert(!pair.second._sign, *pair.second._coef, *pair.second._p->first, *pair.second._p->second);
    }
    //    update_convexity();
    return *this;
}


func_& func_::operator*=(const func_& f){
    if (f.is_constant()) {
        _cst = multiply(_cst, f);
        return *this;
    }
    return *this;
}

func_& func_::operator/=(const func_& f){
    
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





void func_::reset(){
    set_type(func_c);
    _ftype = const_;
    _return_type = integer_;
    _convex = linear_;
    _lterms->clear();
    _qterms->clear();
    _pterms->clear();
    if (!_in_model) {
        for (auto &elem: *_vars) {
            delete elem.second.first;
        }
        for (auto &elem: *_params) {
            delete elem.second.first;
        }
    }
    _vars->clear();
    _params->clear();
    delete _cst;
    _cst = new constant<int>(0);
};

bool func_::insert(bool sign, const constant_& coef, const param_& p){/**< Adds coef*p to the linear function. Returns true if added new term, false if only updated coef of p */
    param_* p_new;
    auto name = p.get_name();
    auto pair_it = _lterms->find(name);
    if (_ftype == const_ && p.is_var()) {
        _ftype = lin_;
    }
    
    if (pair_it == _lterms->end()) {
        auto c_new = copy(&coef);
        if (p.is_var()) {
            p_new = get_var(name);
            if (!p_new) {
                p_new = (param_*)copy(&p);
                add_var(p_new);
            }
            else {
                incr_occ_var(name);
            }
        }
        else {
            p_new = get_param(name);
            if (!p_new) {
                p_new = (param_*)copy(&p);
                add_param(p_new);
            }
            else {
                incr_occ_param(name);
            }
        }
        _lterms->insert(make_pair<>(name, lterm(sign, c_new, p_new)));
        return true;
    }
    else {
        if (pair_it->second._sign == sign) {
            pair_it->second._coef = add(pair_it->second._coef, coef);
            if (!pair_it->second._sign) { // both negative
                pair_it->second._sign = true;
            }
        }
        else{
            pair_it->second._coef = substract(pair_it->second._coef, coef);
        }
        if (pair_it->second._coef->is_zero()) {
            _lterms->erase(pair_it);
            if (p.is_var()) {
                decr_occ_var(name);
            }
            else{
                decr_occ_param(name);
            }
            
        }
        return false;
    }
};

void func_::insert(const lterm& term){
    insert(term._sign, *term._coef, *term._p);
}

bool func_::insert(bool sign, const constant_& coef, const param_& p1, const param_& p2){/**< Adds coef*p1*p2 to the function. Returns true if added new term, false if only updated coef of p1*p2 */
    auto s1 = p1.get_name();
    auto s2 = p2.get_name();
    auto name = s1+","+s2;
    auto pair_it = _qterms->find(name);
    param_* p_new1;
    param_* p_new2;
    
    if (_ftype == const_ && p1.is_var()) {
        _ftype = lin_;
    }
    if (_ftype == lin_ && p1.is_var() && p2.is_var()) {
        _ftype = quad_;
    }
    
    if (pair_it == _qterms->end()) {
        if (p1.is_var()) {
            p_new1 = (param_*)get_var(s1);
            if (!p_new1) {
                p_new1 = (param_*)copy(&p1);
                add_var(p_new1);
            }
            else {
                incr_occ_var(s1);
            }
        }
        else {
            p_new1 = (param_*)get_param(s1);
            if (!p_new1) {
                p_new1 = (param_*)copy(&p1);
                add_param(p_new1);
            }
            else {
                incr_occ_param(s1);
            }
            
        }
        if (p2.is_var()) {
            p_new2 = get_var(s2);
            if (!p_new2) {
                p_new2 = (param_*)copy(&p2);
                add_var(p_new2);
            }
            else {
                incr_occ_var(s2);
            }
        }
        else {
            p_new2 = get_param(s2);
            if (!p_new2) {
                p_new2 = (param_*)copy(&p2);
                add_param(p_new2);
            }
            else {
                incr_occ_param(s2);
            }
        }
        auto c_new = copy(&coef);
        _qterms->insert(make_pair<>(name, qterm(sign, c_new, p_new1, p_new2)));
        return true;
    }
    else {
        if (pair_it->second._sign == sign) {
            pair_it->second._coef = add(pair_it->second._coef, coef);
            if (!pair_it->second._sign) { // both negative
                pair_it->second._sign = true;
            }
        }
        else{
            pair_it->second._coef = substract(pair_it->second._coef, coef);
        }
        
        if (pair_it->second._coef->is_zero()) {
            _qterms->erase(pair_it);
            if (p1.is_var()) {
                decr_occ_var(s1);
            }
            else {
                decr_occ_param(s1);
            }
            if (p2.is_var()) {
                decr_occ_var(s2);
            }
            else {
                decr_occ_param(s2);
            }
//            update_convexity();
        }
        return false;
    }
};

void func_::insert(const qterm& term){
    insert(term._sign, *term._coef, *term._p->first, *term._p->second);
}



uexpr::uexpr(const uexpr& exp){
    _otype = exp._otype;
    _son = copy(exp._son);
    _type = uexp_c;
}

uexpr::uexpr(uexpr&& exp){
    _otype = exp._otype;
    _son = move(exp._son);
    exp._son = nullptr;
    _type = uexp_c;
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

double uexpr::eval(ind i) const{
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
}

bexpr::bexpr(const bexpr& exp){ /**< Copy constructor from binary expression tree */
    _otype = exp._otype;
    _lson = copy(exp._lson);
    _rson =  copy(exp._rson);
    _type = bexp_c;
};

bexpr::bexpr(bexpr&& exp){ /**< Move constructor from binary expression tree */
    _otype = exp._otype;
    _lson = move(exp._lson);
    _rson = move(exp._rson);
    _type = bexp_c;
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
    delete _rson;
    _rson = copy((constant_*)&v);
    _otype = plus_;
    return *this;
}


template<typename other_type> bexpr& bexpr::operator-=(const other_type& v){
    constant_* old = _lson;
    _lson = copy(this);
    delete old;
    delete _rson;
    _rson = copy((constant_*)&v);
    _otype = minus_;
    return *this;
}


template<typename other_type> bexpr& bexpr::operator*=(const other_type& v){
    constant_* old = _lson;
    _lson = copy(this);
    delete old;
    delete _rson;
    _rson = copy((constant_*)&v);
    _otype = product_;
    return *this;
}

template<typename other_type> bexpr& bexpr::operator/=(const other_type& v){
    constant_* old = _lson;
    _lson = copy(this);
    delete old;
    delete _rson;
    _rson = copy((constant_*)&v);
    _otype = div_;
    return *this;
}




void poly_print(const constant_* c){/**< Copy c into c1 detecting the right class, i.e., constant<>, param<>, uexpr or bexpr. */
    
    if (!c) {
        return;
    }
    switch (c->get_type()) {
        case binary_c: {
            ((constant<bool>*)(c))->print();
            break;
        }
        case short_c: {
            ((constant<short>*)(c))->print();
            break;
        }
        case integer_c: {
            ((constant<int>*)(c))->print();
            break;
        }
        case float_c: {
            ((constant<float>*)(c))->print();
            break;
        }
        case double_c: {
            ((constant<double>*)(c))->print();
            break;
        }
        case long_c: {
            ((constant<long double>*)(c))->print();
            break;
        }
        case par_c:{
            auto p_c = (param_*)(c);
            switch (p_c->get_intype()) {
                case binary_:
                    ((param<bool>*)p_c)->print();
                    break;
                case short_:
                    ((param<short>*)p_c)->print();
                    break;
                case integer_:
                    ((param<int>*)p_c)->print();
                    break;
                case float_:
                    ((param<float>*)p_c)->print();
                    break;
                case double_:
                    ((param<double>*)p_c)->print();
                    break;
                case long_:
                    ((param<long double>*)p_c)->print();
                    break;
                default:
                    break;
            }
            break;
        }
        case uexp_c: {
            ((uexpr*)c)->print(false);
            break;
        }
        case bexp_c: {
            ((bexpr*)c)->print(false);
            break;
        }
        case var_c: {
            auto p_c = (param_*)(c);
            switch (p_c->get_intype()) {
                case binary_:
                    ((var<bool>*)p_c)->print();
                    break;
                case short_:
                    ((var<short>*)p_c)->print();
                    break;
                case integer_:
                    ((var<int>*)p_c)->print();
                    break;
                case float_:
                    ((var<float>*)p_c)->print();
                    break;
                case double_:
                    ((var<double>*)p_c)->print();
                    break;
                case long_:
                    ((var<long double>*)p_c)->print();
                    break;
                default:
                    break;
            }
            break;
        }
        case func_c: {
            auto f = (func_*)c;
            f->print(false);
            break;
        }
        default:
            break;
    }
    
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

double eval(ind i, const constant_* c){
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
        case uexp_c: {
            return ((uexpr*)c)->eval(i);
            break;
        }
        case bexp_c: {
            return ((bexpr*)c)->eval(i);
            break;
        }
        default:
            break;
    }
    return 0;
}

double eval(const constant_* c1){
    return eval(0,c1);
};

double bexpr::eval(ind i) const{
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


func_ operator+(const func_& f1, const func_& f2){
    return func_(f1) += f2;
}

func_ operator+(func_&& f1, const func_& f2){
    return f1 += f2;
}

func_ operator+(const func_& f1, func_&& f2){
    return f2 += f1;
}


func_ operator-(const func_& f1, const func_& f2){
    return func_(f1) -= f2;
}

func_ operator-(func_&& f1, const func_& f2){
    return f1 -= f2;
}

func_ operator-(const func_& f1, func_&& f2){
    return (f2 *= -1) += f1;
}

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
                res->insert(true, constant<int>(1), *(param_*)c1);
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
            res->insert(true, constant<int>(1), *(param_*)c1);
            return c1 = res;
            break;
        }
        case uexp_c: {
//            auto res = new bexpr(*(uexpr*)c1 + c2);
//            delete c1;
//            c1 = (constant_*)res;
//            return c1;
//            break;
        }
        case bexp_c: {
//            auto res = new bexpr(*(bexpr*)c1 + c2);
//            delete c1;
//            c1 = (constant_*)res;
//            return c1;
//            break;
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
        case uexp_c: {
            auto res = new bexpr(*(uexpr*)c1 + c2);
            delete c1;
            c1 = (constant_*)res;
            return c1;
            break;
        }
        case bexp_c: {
            auto res = new bexpr(*(bexpr*)c1 + c2);
            delete c1;
            c1 = (constant_*)res;
            return c1;
            break;
        }
        case func_c: {
            return add(c1, (func_)c2);
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
        case uexp_c: {
            auto res = new bexpr(*(uexpr*)c1 - c2);
            delete c1;
            c1 = (constant_*)res;
            return c1;
            break;
        }
        case bexp_c: {
            auto res = new bexpr(*(bexpr*)c1 - c2);
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
            auto res = new bexpr(*(uexpr*)c1 * c2);
            delete c1;
            c1 = (constant_*)res;
            return c1;
            break;
        }
        case bexp_c: {
            auto res = new bexpr(*(bexpr*)c1 * c2);
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
        case func_c: {
            auto f = new func_(*c1);
            delete c1;
            *f *= c2;
            return c1 = (constant_*)f;
            break;
        }
        default:
            break;
    }
    return nullptr;
}



void const func_::print(bool endline){
    param_* p_new1;
    param_* p_new2;
    constant_* c_new;
    int ind = 0;
    string sign = " + ";
    if (_convex==convex_) {
        cout << "Convex function: ";
    }
    if (_convex==concave_) {
        cout << "Concave function: ";
    }
    cout << "f(";
    for (auto pair_it = _vars->begin(); pair_it != _vars->end();) {
        poly_print(pair_it->second.first);
        if (++pair_it != _vars->end()) {
            cout << ",";
        }
    }
    cout << ") = ";
    for (auto &pair:*_qterms) {
        c_new = pair.second._coef;
        p_new1 = (param_*)pair.second._p->first;
        p_new2 = (param_*)pair.second._p->second;
        if (c_new->is_number()){
            auto v = eval(c_new);
            if (pair.second._sign) {
                if (v==-1) {
                    cout << " - ";
                }
                else if (ind>0) {
                    cout << " + ";
                    if(v!=1) {
                        cout << v;
                    }
                }
                else if(v!=1) {
                    cout << v;
                }
            }
            if(!pair.second._sign) {
                if (v == -1 && ind>0) {
                    cout << " + ";
                }
                else if (v < 0 && ind>0){
                    cout << " + " << -1*v;
                }
                else if (v==1){
                    cout << " - ";
                }
                else if(v!=-1){
                    cout << " - " << v;
                }
            }
        }
        else{
            if(ind > 0) {
                if (!pair.second._sign) {
                    cout << " - ";
                }
                else {
                    cout << " + ";
                }
            }
            cout << "(";
            poly_print(c_new);
            cout << ")";
        }
        poly_print(p_new1);
        cout << ".";
        poly_print(p_new2);
        ind++;
    }
    if (!_qterms->empty() && !_lterms->empty()) {
        cout << " + ";
    }
    param_* p_new;
    for (auto &pair:*_lterms) {
        c_new = pair.second._coef;
        p_new = (param_*)pair.second._p;
        if (c_new->is_number()){
            auto v = eval(c_new);
            if (pair.second._sign) {
                if (v==-1) {
                    cout << " - ";
                }
                else if (ind>0) {
                    cout << " + ";
                    if(v!=1) {
                        cout << v;
                    }
                }
                else if(v!=1) {
                    cout << v;
                }
            }
            if(!pair.second._sign) {
                if (v == -1 && ind>0) {
                    cout << " + ";
                }
                else if (v < 0 && ind>0){
                    cout << " + " << -1*v;
                }
                else if (v==1){
                    cout << " - ";
                }
                else if(v!=-1){
                    cout << " - " << v;
                }
            }
        }
        else{
            if(ind > 0) {
                if (!pair.second._sign) {
                    cout << " - ";
                }
                else {
                    cout << " + ";
                }
            }
            cout << "(";
            poly_print(c_new);
            cout << ")";
        }
        poly_print(p_new);
        ind++;
    }
    if (_cst->is_number()) {
        auto val = eval(_cst);
        if (val < 0) {
            cout <<" - " << -val;
        }
        else if (val > 0){
            cout <<" + " << val;
        }
    }
    else {
        cout <<" + (";
        poly_print(_cst);
        cout <<")";
    }
    if (endline)
        cout << endl;
}