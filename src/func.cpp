//
//  func.cpp
//  Gravity
//
//  Created by Hijazi, Hassan on 24/10/16.
//
//
#include <math.h>
#include <Gravity/func.h>



//pair<ind,param_*> param_::operator[](ind i){
//    return make_pair(i, (param_*)copy((constant_*)this));
//}

lin::lin(constant_* coef, param_* p):lin(){
    _lterms->insert(make_pair<>(p->get_name(), lterm(coef, p)));
    if (p->is_var()) {
        _vars->insert((var_*)p);
    }
};

lin::lin(param_* p):lin(){
    _lterms->insert(make_pair<>(p->get_name(), lterm(p)));
    if (p->is_var()) {
        _vars->insert((var_*)p);
    }
};

lin::~lin(){
    delete _lterms;
    delete _vars;
    delete _cst;
};


bool lin::insert(bool sign, constant_* coef, param_* p){/**< Adds coef*p to the linear function. Returns true if added new term, false if only updated coef of p */
    auto name = p->get_name();
    auto pair_it = _lterms->find(name);
    if (pair_it == _lterms->end()) {
        auto p_new = (param_*)copy(p);
        auto c_new = copy(coef);
        _lterms->insert(make_pair<>(name, lterm(sign, c_new, p_new)));
        if (p->is_var()) {
            _vars->insert((var_*)p_new);
        }
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
            _vars->erase((var_*)pair_it->second._p);
            _lterms->erase(pair_it);
        }
        return false;
    }
};


void lin::reverse_sign(){ /*<< Reverse the sign of all linear terms and constant in the function */
    for (auto &pair: *_lterms) {
        pair.second.reverse_sign();
    }
    ::reverse_sign(_cst);
}

bool lterm::operator==(const lterm& l) const{
    return (equals(_coef, l._coef) && equals(_p, l._p));
}

bool lterm::operator!=(const lterm& l) const{
    return (!equals(_coef, l._coef) || !equals(_p, l._p));
}

bool lin::operator==(const lin& l) const{
    if(!equals(_cst, l._cst)){
        return false;
    }
    auto it = l._lterms->begin();
    for (auto &pair: *_lterms) {
        it = l._lterms->find(pair.first);
        if (it==l._lterms->end()) {
            return false;
        }
        else if (it->second != pair.second){
            return false;
        }
    }
    return true;
}

bool lin::operator!=(const lin& l) const{
    return !(*this==l);
}


lin& lin::operator=(const lin& l){
    func_::set_ftype(lin_);
    _convex = l._convex;
    delete _lterms;
    _lterms = new map<string, lterm>();
    delete _vars;
    _vars = new set<var_*>;
    param_* p_new;
    constant_* c_new;
    for (auto &pair:*l._lterms) {
        p_new = (param_*)copy(pair.second._p);
        c_new = copy(pair.second._coef);
        _lterms->insert(make_pair<>(p_new->get_name(), lterm(pair.second._sign, c_new, p_new)));
        if (p_new->is_var()) {
            _vars->insert((var_*) p_new);
        }
    }
    delete _cst;
    _cst = copy(l._cst);
    return *this;
};


lin& lin::operator=(lin&& l){
    delete _lterms;
    _lterms = l._lterms;
    l._lterms = nullptr;
    delete _vars;
    _vars = l._vars;
    l._vars = nullptr;
    delete _cst;
    _cst = l._cst;
    l._cst = nullptr;
    return *this;
};



bool constant_::is_unit() const{ /**< Returns true if constant equals 1 */
    return (is_number() && eval(this)==1);
}

bool constant_::is_zero() const{ /**< Returns true if constant equals 0 */
    return (is_number() && eval(this)==0);
}


bool constant_::is_neg_unit() const{ /**< Returns true if constant equals -1 */
    return (is_number() && eval(this)==-1);
}


//template <typename T> lin operator*(const lin& l, const constant<T>& c){
//    return lin(l) *= c;
//}
//
//template lin operator*(const lin& l, const constant<short>& c);
//template lin operator*(const lin& l, const constant<int>& c);
//template lin operator*(const lin& l, const constant<float>& c);
//template lin operator*(const lin& l, const constant<double>& c);
//template lin operator*(const lin& l, const constant<long double>& c);


//missing long double

void const lin::print(bool endline){
    param_* p_new;
    constant_* c_new;
    int ind = 0;
    string sign = " + ";
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


lterm& lterm::operator=(const lterm& l){
    delete _coef;
    _coef = copy(l._coef);
    delete _p;
    _p = (param_*)copy(l._p);
    _sign = l._sign;
    return *this;
};

lterm& lterm::operator=(lterm&& l){
    delete _coef;
    _coef = l._coef;
    l._coef = nullptr;
    delete _p;
    _p = l._p;
    l._p = nullptr;
    _sign = l._sign;
    return *this;
};

/* Polymorphic functions */

void reverse_sign(constant_* c){ /**< Reverses the sign of the constant. */
    switch (c->get_type()) {
        case binary_c: {
            ((constant<bool>*)c)->set_val(!((constant<bool>*)c)->eval());
            break;
        }
        case short_c: {
            ((constant<short>*)c)->set_val(-1*((constant<short>*)c)->eval());
            break;
        }
        case integer_c: {
            ((constant<int>*)c)->set_val(-1*((constant<int>*)c)->eval());
            break;
        }
        case float_c: {
            ((constant<float>*)c)->set_val(-1*((constant<float>*)c)->eval());
            break;
        }
        case double_c: {
            ((constant<double>*)c)->set_val(-1*((constant<double>*)c)->eval());
            break;
        }
        case long_c: {
            ((constant<long double>*)c)->set_val(-1*((constant<long double>*)c)->eval());
            break;
        }
        default:
            throw invalid_argument("Cannot reverse sign of non-numeric constant");
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
                if (f1->_return_type != f2->_return_type || f1->get_ftype() != f2->get_ftype() || f1->_convex != f2->_convex) {
                    return false;
                }
                switch (f1->get_ftype()) {
                    case lin_:{
                        auto l1 = (lin*)c1;
                        auto l2 = (lin*)c2;
                        return *l1 == *l2;
                        break;
                    }
                    case quad_: {
                        
                        break;
                    }
                    case pol_: {
                        
                        break;
                    }
                    case nlin_: {
                        
                        break;
                    }
                    default:
                        break;
                }
            }
            break;
        }
        default:
            break;
    }
    return false;
}



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
            auto f = (func_*)c2;
            switch (f->get_ftype()) {
                case lin_:
                    return new lin(*(lin*)f);
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
            switch (f->get_ftype()) {
                case lin_:
                    ((lin*)f)->print(false);
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



constant_* add(constant_* c1, const lin& l1){
    switch (c1->get_type()) {
        case binary_c: {
            auto val = ((constant<bool>*)c1);
            delete c1;
            auto l = new lin(l1);
            *l += *val;
            c1 = (constant_*)(l);
            return c1;
            break;
        }
        case short_c: {
            auto val = ((constant<short>*)c1);
            delete c1;
            auto l = new lin(l1);
            *l += *val;
            c1 = (constant_*)(l);
            return c1;
            break;
        }
        case integer_c: {
            auto val = ((constant<int>*)c1);
            delete c1;
            auto l = new lin(l1);
            *l += *val;
            c1 = (constant_*)(l);
            return c1;
            break;
        }
        case float_c: {
            auto val = ((constant<float>*)c1);
            delete c1;
            auto l = new lin(l1);
            *l += *val;
            c1 = (constant_*)(l);
            return c1;
            break;
        }
        case double_c: {
            auto val = ((constant<double>*)c1);
            delete c1;
            auto l = new lin(l1);
            *l += *val;
            c1 = (constant_*)(l);
            return c1;
            break;
        }
        case long_c: {
            auto val = ((constant<long double>*)c1);
            delete c1;
            auto l = new lin(l1);
            *l += *val;
            c1 = (constant_*)(l);
            return c1;
            break;
        }
        case par_c:{
            auto pc1 = (param_*)(c1);
            auto l = new lin(pc1);
            *l += l1;
            if (l->is_constant()) {
                c1 = copy(l->get_cst());
                delete l;
            }
            else {
                c1 =(constant_*)(l);
            }
            return c1;
            break;
        }
        case uexp_c: {
            auto res = new bexpr(*(uexpr*)c1 + l1);
            delete c1;
            c1 = (constant_*)res;
            return c1;
            break;
        }
        case bexp_c: {
            auto res = new bexpr(*(bexpr*)c1 + l1);
            delete c1;
            c1 = (constant_*)res;
            return c1;
            break;
        }
        case func_c: {
            switch (((func_*)c1)->get_ftype()) {
                case lin_: {
                    auto res = new lin(*(lin*)c1 + l1);
                    delete c1;
                    if (res->is_constant()) {
                        c1 = copy(res->get_cst());
                        delete res;
                        return c1;
                    }
                    else {
                        return c1 = res;
                    }
                    break;
                }
                default:
                    break;
            }
        }
        default:
            break;
    }
    return nullptr;
}


constant_* multiply(constant_* c1, const lin& l1){
    switch (c1->get_type()) {
        case binary_c: {
            auto val = ((constant<bool>*)c1);
            delete c1;
            auto l = new lin(l1);
            *l *= *val;
            c1 = (constant_*)(l);
            return c1;
            break;
        }
        case short_c: {
            auto val = ((constant<short>*)c1);
            delete c1;
            auto l = new lin(l1);
            *l *= *val;
            c1 = (constant_*)(l);
            return c1;
            break;
        }
        case integer_c: {
            auto val = ((constant<int>*)c1);
            delete c1;
            auto l = new lin(l1);
            *l *= *val;
            c1 = (constant_*)(l);
            return c1;
            break;
        }
        case float_c: {
            auto val = ((constant<float>*)c1);
            delete c1;
            auto l = new lin(l1);
            *l *= *val;
            c1 = (constant_*)(l);
            return c1;
            break;
        }
        case double_c: {
            auto val = ((constant<double>*)c1);
            delete c1;
            auto l = new lin(l1);
            *l *= *val;
            c1 = (constant_*)(l);
            return c1;
            break;
        }
        case long_c: {
            auto val = ((constant<long double>*)c1);
            delete c1;
            auto l = new lin(l1);
            *l *= *val;
            c1 = (constant_*)(l);
            return c1;
            break;
        }
        case par_c:{
            auto pc1 = (param_*)(c1);
            auto l = new lin(l1);
            switch (pc1->get_intype()) {
                case binary_:
                    *l *= *((param<bool>*)pc1);
                    delete c1;
                    return c1 = l;
                    break;
                case short_:
                    *l *= *((param<short>*)pc1);
                    delete c1;
                    return c1 = l;
                    break;
                case integer_:
                    *l *= *((param<int>*)pc1);
                    delete c1;
                    return c1 = l;
                    
                    break;
                case float_:
                    *l *= *((param<float>*)pc1);
                    delete c1;
                    return c1 = l;
                    break;
                case double_:
                    *l *= *((param<double>*)pc1);
                    delete c1;
                    return c1 = l;
                    break;
                case long_:
                    *l *= *((param<long double>*)pc1);
                    delete c1;
                    return c1 = l;
                    break;
                default:
                    break;
            }
            break;
        }
        case uexp_c: {
            auto res = new bexpr(*(uexpr*)c1 * l1);
            delete c1;
            c1 = (constant_*)res;
            return c1;
            break;
        }
        case bexp_c: {
            auto res = new bexpr(*(bexpr*)c1 * l1);
            delete c1;
            c1 = (constant_*)res;
            return c1;
            break;
        }
        case func_c: {
            switch (((func_*)c1)->get_ftype()) {
                case lin_: {
//                    auto res = new lin(*(lin*)c1 + l1);
//                    delete c1;
//                    return c1 = res;
                    cerr << "Unsupported yet!\n";
                    break;
                }
                default:
                    break;
            }
        }
        default:
            break;
    }
    return nullptr;
}


constant_* substract(constant_* c1, const lin& l1){
    switch (c1->get_type()) {
        case binary_c: {
            auto val = ((constant<bool>*)c1);
            delete c1;
            auto l = new lin(l1);
            *l -= *val;
            c1 = (constant_*)(l);
            return c1;
            break;
        }
        case short_c: {
            auto val = ((constant<short>*)c1);
            delete c1;
            auto l = new lin(l1);
            *l -= *val;
            c1 = (constant_*)(l);
            return c1;
            break;
        }
        case integer_c: {
            auto val = ((constant<int>*)c1);
            delete c1;
            auto l = new lin(l1);
            *l -= *val;
            c1 = (constant_*)(l);
            return c1;
            break;
        }
        case float_c: {
            auto val = ((constant<float>*)c1);
            delete c1;
            auto l = new lin(l1);
            *l -= *val;
            c1 = (constant_*)(l);
            return c1;
            break;
        }
        case double_c: {
            auto val = ((constant<double>*)c1);
            delete c1;
            auto l = new lin(l1);
            *l -= *val;
            c1 = (constant_*)(l);
            return c1;
            break;
        }
        case long_c: {
            auto val = ((constant<long double>*)c1);
            delete c1;
            auto l = new lin(l1);
            *l -= *val;
            c1 = (constant_*)(l);
            return c1;
            break;
        }
        case par_c:{
            auto pc1 = (param_*)(c1);
            auto l = new lin(pc1);
            *l -= l1;
            c1 =(constant_*)(l);
            return c1;
            break;
        }
        case uexp_c: {
            auto res = new bexpr(*(uexpr*)c1 - l1);
            delete c1;
            c1 = (constant_*)res;
            return c1;
            break;
        }
        case bexp_c: {
            auto res = new bexpr(*(bexpr*)c1 - l1);
            delete c1;
            c1 = (constant_*)res;
            return c1;
            break;
        }
        case func_c: {
            switch (((func_*)c1)->get_ftype()) {
                case lin_: {
                    auto res = new lin(*(lin*)c1 - l1);
                    delete c1;
                    return c1 = res;
                    break;
                }
                default:
                    break;
            }
        }
        default:
            break;
            
    }
    return nullptr;
}


constant_* add(constant_* c1, constant_* c2){ /**< adds c2 to c1, updates its type and returns the result **/
    switch (c2->get_type()) {
        case binary_c: {
            return add(c1, *(constant<bool>*)c2);
            break;
        }
        case short_c: {
            return add(c1, *(constant<short>*)c2);
            break;
        }
        case integer_c: {
            return add(c1, *(constant<int>*)c2);
            break;
        }
        case float_c: {
            return add(c1, *(constant<float>*)c2);
            break;
        }
        case double_c: {
            return add(c1, *(constant<double>*)c2);
            break;
        }
        case long_c: {
            return add(c1, *(constant<long double>*)c2);
            break;
        }
        case par_c:{
            auto pc2 = (param_*)(c2);
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
            switch (((func_*)c2)->get_ftype()) {
                case lin_:
                    return add(c1, (*(lin*)c2));
                    break;
                    
                default:
                    break;
            }
        }
        default:
            break;
    }
    return nullptr;
}

constant_* substract(constant_* c1, constant_* c2){ /**< adds c2 to c1, updates its type and returns the result **/
    switch (c2->get_type()) {
        case binary_c: {
            return substract(c1, *(constant<bool>*)c2);
            break;
        }
        case short_c: {
            return substract(c1, *(constant<short>*)c2);
            break;
        }
        case integer_c: {
            return substract(c1, *(constant<int>*)c2);
            break;
        }
        case float_c: {
            return substract(c1, *(constant<float>*)c2);
            break;
        }
        case double_c: {
            return substract(c1, *(constant<double>*)c2);
            break;
        }
        case long_c: {
            return substract(c1, *(constant<long double>*)c2);
            break;
        }
        case par_c:{
            auto pc2 = (param_*)(c2);
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
            switch (((func_*)c2)->get_ftype()) {
                case lin_:
                    return substract(c1, (*(lin*)c2));
                    break;
                    
                default:
                    break;
            }
        }
        default:
            break;
    }
    return nullptr;
}

constant_* multiply(constant_* c1, constant_* c2){ /**< adds c2 to c1, updates its type and returns the result **/
    switch (c2->get_type()) {
        case binary_c: {
            return multiply(c1, *(constant<bool>*)c2);
            break;
        }
        case short_c: {
            return multiply(c1, *(constant<short>*)c2);
            break;
        }
        case integer_c: {
            return multiply(c1, *(constant<int>*)c2);
            break;
        }
        case float_c: {
            return multiply(c1, *(constant<float>*)c2);
            break;
        }
        case double_c: {
            return multiply(c1, *(constant<double>*)c2);
            break;
        }
        case long_c: {
            return multiply(c1, *(constant<long double>*)c2);
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
        case func_c: {
            switch (((func_*)c2)->get_ftype()) {
                case lin_:
                    return multiply(c1, (*(lin*)c2));
                    break;
                    
                default:
                    break;
            }
        }
        default:
            break;
    }
    return nullptr;
}

lin operator+(const lin& l1, const lin& l2){
    return lin(l1) += l2;
}

lin operator+(lin&& l1, const lin& l2){
    return l1 += l2;
}


lin operator-(const lin& l1, const lin& l2){
    return lin(l1) -= l2;
}

lin operator-(lin&& l1, const lin& l2){
    return l1 -= l2;
}

template <typename T> lterm& lterm::operator*=(const param<T>& c){
    if (_coef->is_number()) {
        _coef = multiply(_coef, c);
    }
    else if (_coef->is_param()) {
//        _coef = (constant_*)new lin(new constant<T>(c), (param_*)_coef);
        cerr << "Unsupported yet;\n";
    }
    else if (_coef->is_function()) {
        auto f = (func_*)_coef;
        switch (f->get_ftype()) {
            case lin_: {
//                auto l = (lin*) f;
//                *l *= c;
                cerr << "Unsupported yet;\n";
                break;
            }
            default:
                break;
        }
    }
    return *this;
}

template <typename T> lterm& lterm::operator+=(const constant<T>& c){
    if (c.is_zero()) {
        return *this;
    }
    if (_coef->is_number()) {
        _coef = add(_coef, c);
    }
    else if (_coef->is_param()) {
        auto l = new lin((param_*)_coef);
        *l += c;
        _coef = (constant_*)l;
    }
    else if (_coef->is_function()) {
        auto f = (func_*)_coef;
        switch (f->get_ftype()) {
            case lin_: {
                auto l = (lin*) f;
                *l += c;
                break;
            }
            default:
                break;
        }
    }
    return *this;
}


template <typename T> lterm& lterm::operator*=(const constant<T>& c){
    if (c.is_zero()) {
        delete _coef;
        _coef = new constant<int>(0);
        _sign = true;
        cerr << "\nWARNING multiplying lterm by zero!\n";
        return *this;
    }
    if (c.is_unit()) {
        return *this;
    }
    
    if (_coef->is_number()) {
        _coef = multiply(_coef, c);
    }
    else if (_coef->is_param()) {
        _coef = (constant_*)new lin(new constant<T>(c), (param_*)_coef);                
    }
    else if (_coef->is_function()) {
        auto f = (func_*)_coef;
        switch (f->get_ftype()) {
            case lin_: {
                auto l = (lin*) f;
                *l *= c;
                break;
            }
            default:
                break;
        }
    }
    return *this;
}


template <typename T> lterm& lterm::operator/=(const constant<T>& c){
    if (c.is_zero()) {
        throw invalid_argument("\nERROT dividing lterm by zero!\n");
        return *this;
    }
    if (c.is_unit()) {
        return *this;
    }
    if (_coef->is_number()) {
        _coef = divide(_coef, c);
    }
    else if (_coef->is_param()) {
        auto l = new lin((param_*)_coef);
        *l /= c;
        _coef = (constant_*)l;
    }
    else if (_coef->is_function()) {
        auto f = (func_*)_coef;
        switch (f->get_ftype()) {
            case lin_: {
                auto l = (lin*) f;
                *l /= c;
                break;
            }
            default:
                break;
        }
    }
    return *this;
}






lin& lin::operator+=(const lin& l){
    _cst = add(_cst, l._cst);
    for (auto &pair:*l._lterms) {
        this->insert(pair.second._sign, pair.second._coef, pair.second._p);
    }
    return *this;
}

lin& lin::operator-=(const lin& l){
    return *this += -1*l;
}

//template<typename T> lin& lin::operator+=(const constant<T>& c){
//    if (c.is_zero()) {
//        return *this;
//    }
//    _cst = add(_cst, c);
//    return *this;
//}

//template lin& lin::operator+=(const constant<bool>& c);
//template lin& lin::operator+=(const constant<short>& c);
//template lin& lin::operator+=(const constant<int>& c);
//template lin& lin::operator+=(const constant<float>& c);
//template lin& lin::operator+=(const constant<double>& c);
//template lin& lin::operator+=(const constant<long double>& c);


//template<typename T> lin& lin::operator-=(const constant<T>& c){
//    if (c.is_zero()) {
//        return *this;
//    }
//    _cst = substract(_cst, c);
//    return *this;
//}
//
//template lin& lin::operator-=(const constant<bool>& c);
//template lin& lin::operator-=(const constant<short>& c);
//template lin& lin::operator-=(const constant<int>& c);
//template lin& lin::operator-=(const constant<float>& c);
//template lin& lin::operator-=(const constant<double>& c);
//template lin& lin::operator-=(const constant<long double>& c);

//template<typename T> lin& lin::operator+=(const param<T>& p){
//    
//    auto pair_it = _lterms->find(p.get_name());
//    if (pair_it == _lterms->end()) {
//        auto p_new = (param_*)copy((constant_*)&p);
//        _lterms->insert(make_pair<>(p_new->get_name(), lterm(p_new)));
//        if (p.is_var()) {
//            _vars->insert((var_*)p_new);
//        }
//    }
//    else {
//        auto lterm = pair_it->second;
//        if (lterm._coef->is_neg_unit()) {
//            _lterms->erase(pair_it);
//        }
//        else{
//            lterm += constant<bool>(1);
//        }
//    }
//    return *this;
//}
//
//template lin& lin::operator+=(const param<bool>& c);
//template lin& lin::operator+=(const param<short>& c);
//template lin& lin::operator+=(const param<int>& c);
//template lin& lin::operator+=(const param<float>& c);
//template lin& lin::operator+=(const param<double>& c);
//template lin& lin::operator+=(const param<long double>& c);

//template<typename T> lin& lin::operator-=(const param<T>& p){
//    
//    auto pair_it = _lterms->find(p.get_name());
//    if (pair_it == _lterms->end()) {
//        auto p_new = (param_*)copy((constant_*)&p);
//        _lterms->insert(make_pair<>(p_new->get_name(), lterm(new constant<short>(-1), p_new)));
//        if (p.is_var()) {
//            _vars->insert((var_*)p_new);
//        }
//    }
//    else {
//        auto lterm = pair_it->second;
//        if (lterm._coef->is_unit()) {
//            _lterms->erase(pair_it);
//        }
//        else{
//            lterm += constant<short>(-1);
//        }
//    }
//    return *this;
//}
//
//template lin& lin::operator-=(const param<bool>& c);
//template lin& lin::operator-=(const param<short>& c);
//template lin& lin::operator-=(const param<int>& c);
//template lin& lin::operator-=(const param<float>& c);
//template lin& lin::operator-=(const param<double>& c);
//template lin& lin::operator-=(const param<long double>& c);

//template<typename T> lin& lin::operator*=(const param<T>& p){
//    if (!_cst->is_zero()) {
//        auto pair_it = _lterms->find(p.get_name());
//        if (pair_it == _lterms->end()) {
//            auto p_new = (param_*)copy((constant_*)&p);
//            _lterms->insert(make_pair<>(p_new->get_name(), lterm(_cst, p_new)));
//        }
//        else {
//            auto lterm = pair_it->second;
//            switch (_cst->get_type()) {
//                case binary_c: {
//                    lterm += (*(constant<bool>*)_cst);
//                    break;
//                }
//                case short_c: {
//                    lterm += (*(constant<short>*)_cst);
//                    break;
//                }
//                case integer_c: {
//                    lterm += (*(constant<int>*)_cst);
//                    break;
//                }
//                case float_c: {
//                    lterm += (*(constant<float>*)_cst);
//                    break;
//                }
//                case double_c: {
//                    lterm += (*(constant<double>*)_cst);
//                    break;
//                }
//                case long_c: {
//                    lterm += (*(constant<long double>*)_cst);
//                    break;
//                }
//                default:
//                    break;
//            }
//        }
//        _cst = new constant<bool>(false);
//    }
//    for (auto &pair:*_lterms) {
//        pair.second *= p;
//    }
//    return *this;
//}
//
//template lin& lin::operator*=(const param<bool>& c);
//template lin& lin::operator*=(const param<short>& c);
//template lin& lin::operator*=(const param<int>& c);
//template lin& lin::operator*=(const param<float>& c);
//template lin& lin::operator*=(const param<double>& c);
//template lin& lin::operator*=(const param<long double>& c);





//template<typename T> lin& lin::operator*=(const constant<T>& c){
//    if (c.is_unit()) {
//        return *this;
//    }
//    if (c.is_zero()) {
//        reset();
//        return *this;
//    }
//    _cst = multiply(_cst, c);
//    
//    for (auto &pair:*_lterms) {
//        pair.second *= c;
//    }
//    return *this;
//}
//
//template lin& lin::operator*=(const constant<short>& c);
//template lin& lin::operator*=(const constant<int>& c);
//template lin& lin::operator*=(const constant<float>& c);
//template lin& lin::operator*=(const constant<double>& c);
//template lin& lin::operator*=(const constant<long double>& c);


//template<typename T> lin& lin::operator/=(const constant<T>& c){
//    if (c.is_zero()) {
//        throw invalid_argument("\nERROT dividing linear function by zero!\n");
//        return *this;
//    }
//    if (c.is_unit()) {
//        return *this;
//    }
//    _cst = divide(_cst, c);
//    
//    for (auto &pair:*_lterms) {
//        pair.second /= c;
//    }
//    return *this;
//}

lin operator+(const lin& l, bool c){
    return lin(l) += constant<bool>(c);
}

lin operator+(const lin& l, int c){
    return lin(l) += constant<int>(c);
}

lin operator+(const lin& l, short c){
    return lin(l) += constant<short>(c);
}

lin operator+(const lin& l, float c){
    return lin(l) += constant<float>(c);
}

lin operator+(const lin& l, double c){
    return lin(l) += constant<double>(c);
}

lin operator+(const lin& l, long double c){
    return lin(l) += constant<long double>(c);
}

lin operator-(const lin& l, bool c){
    return lin(l) -= constant<bool>(c);
}

lin operator-(const lin& l, int c){
    return lin(l) -= constant<int>(c);
}

lin operator-(const lin& l, short c){
    return lin(l) -= constant<short>(c);
}

lin operator-(const lin& l, float c){
    return lin(l) -= constant<float>(c);
}

lin operator-(const lin& l, double c){
    return lin(l) -= constant<double>(c);
}

lin operator-(const lin& l, long double c){
    return lin(l) -= constant<long double>(c);
}

lin operator*(const lin& l, bool c){
    return lin(l) *= constant<bool>(c);
}

lin operator*(const lin& l, int c){
    return lin(l) *= constant<int>(c);
}

lin operator*(const lin& l, short c){
    return lin(l) *= constant<short>(c);
}

lin operator*(const lin& l, float c){
    return lin(l) *= constant<float>(c);
}

lin operator*(const lin& l, double c){
    return lin(l) *= constant<double>(c);
}

lin operator*(const lin& l, long double c){
    return lin(l) *= constant<long double>(c);
}

lin operator/(const lin& l, bool c){
    return lin(l) /= constant<bool>(c);
}

lin operator/(const lin& l, int c){
    return lin(l) /= constant<int>(c);
}

lin operator/(const lin& l, short c){
    return lin(l) /= constant<short>(c);
}

lin operator/(const lin& l, float c){
    return lin(l) /= constant<float>(c);
}

lin operator/(const lin& l, double c){
    return lin(l) /= constant<double>(c);
}

lin operator/(const lin& l, long double c){
    return lin(l) /= constant<long double>(c);
}

lin operator+(lin&& l, bool c){
    return l += constant<bool>(c);
}

lin operator+(lin&& l, int c){
    return l += constant<int>(c);
}

lin operator+(lin&& l, short c){
    return l += constant<short>(c);
}

lin operator+(lin&& l, float c){
    return l += constant<float>(c);
}

lin operator+(lin&& l, double c){
    return l += constant<double>(c);
}

lin operator+(lin&& l, long double c){
    return l += constant<long double>(c);
}

lin operator-(lin&& l, bool c){
    return l -= constant<bool>(c);
}

lin operator-(lin&& l, int c){
    return l -= constant<int>(c);
}

lin operator-(lin&& l, short c){
    return l -= constant<short>(c);
}

lin operator-(lin&& l, float c){
    return l -= constant<float>(c);
}

lin operator-(lin&& l, double c){
    return l -= constant<double>(c);
}

lin operator-(lin&& l, long double c){
    return l -= constant<long double>(c);
}

lin operator*(lin&& l, bool c){
    return l *= constant<bool>(c);
}

lin operator*(lin&& l, int c){
    return l *= constant<int>(c);
}

lin operator*(lin&& l, short c){
    return l *= constant<short>(c);
}

lin operator*(lin&& l, float c){
    return l *= constant<float>(c);
}

lin operator*(lin&& l, double c){
    return l *= constant<double>(c);
}

lin operator*(lin&& l, long double c){
    return l *= constant<long double>(c);
}

lin operator/(lin&& l, bool c){
    return l /= constant<bool>(c);
}

lin operator/(lin&& l, int c){
    return l /= constant<int>(c);
}

lin operator/(lin&& l, short c){
    return l /= constant<short>(c);
}

lin operator/(lin&& l, float c){
    return l /= constant<float>(c);
}

lin operator/(lin&& l, double c){
    return l /= constant<double>(c);
}

lin operator/(lin&& l, long double c){
    return l /= constant<long double>(c);
}


lin operator+(bool c, const lin& l){
    return lin(l) += constant<bool>(c);
}

lin operator+(int c, const lin& l){
    return lin(l) += constant<int>(c);
}

lin operator+(short c, const lin& l){
    return lin(l) += constant<short>(c);
}

lin operator+(float c, const lin& l){
    return lin(l) += constant<float>(c);
}

lin operator+(double c, const lin& l){
    return lin(l) += constant<double>(c);
}

lin operator+(long double c, const lin& l){
    return lin(l) += constant<long double>(c);
}

lin operator-(bool c, const lin& l){
    return lin(l) -= constant<bool>(c);
}

lin operator-(int c, const lin& l){
    return lin(l) -= constant<int>(c);
}

lin operator-(short c, const lin& l){
    return lin(l) -= constant<short>(c);
}

lin operator-(float c, const lin& l){
    return lin(l) -= constant<float>(c);
}

lin operator-(double c, const lin& l){
    return lin(l) -= constant<double>(c);
}

lin operator-(long double c, const lin& l){
    return lin(l) -= constant<long double>(c);
}

lin operator*(bool c, const lin& l){
    return lin(l) *= constant<bool>(c);
}

lin operator*(int c, const lin& l){
    return lin(l) *= constant<int>(c);
}

lin operator*(short c, const lin& l){
    return lin(l) *= constant<short>(c);
}

lin operator*(float c, const lin& l){
    return lin(l) *= constant<float>(c);
}

lin operator*(double c, const lin& l){
    return lin(l) *= constant<double>(c);
}

lin operator*(long double c, const lin& l){
    return lin(l) *= constant<long double>(c);
}

lin operator/(bool c, const lin& l){
    return lin(l) /= constant<bool>(c);
}

lin operator/(int c, const lin& l){
    return lin(l) /= constant<int>(c);
}

lin operator/(short c, const lin& l){
    return lin(l) /= constant<short>(c);
}

lin operator/(float c, const lin& l){
    return lin(l) /= constant<float>(c);
}

lin operator/(double c, const lin& l){
    return lin(l) /= constant<double>(c);
}

lin operator/(long double c, const lin& l){
    return lin(l) /= constant<long double>(c);
}

lin operator+(bool c, lin&& l){
    return l += constant<bool>(c);
}

lin operator+(int c, lin&& l){
    return l += constant<int>(c);
}

lin operator+(short c, lin&& l){
    return l += constant<short>(c);
}

lin operator+(float c, lin&& l){
    return l += constant<float>(c);
}

lin operator+(double c, lin&& l){
    return l += constant<double>(c);
}

lin operator+(long double c, lin&& l){
    return l += constant<long double>(c);
}

lin operator-(bool c, lin&& l){
    return l -= constant<bool>(c);
}

lin operator-(int c, lin&& l){
    return l -= constant<int>(c);
}

lin operator-(short c, lin&& l){
    return l -= constant<short>(c);
}

lin operator-(float c, lin&& l){
    return l -= constant<float>(c);
}

lin operator-(double c, lin&& l){
    return l -= constant<double>(c);
}

lin operator-(long double c, lin&& l){
    return l -= constant<long double>(c);
}

lin operator*(bool c, lin&& l){
    return l *= constant<bool>(c);
}

lin operator*(int c, lin&& l){
    return l *= constant<int>(c);
}

lin operator*(short c, lin&& l){
    return l *= constant<short>(c);
}

lin operator*(float c, lin&& l){
    return l *= constant<float>(c);
}

lin operator*(double c, lin&& l){
    return l *= constant<double>(c);
}

lin operator*(long double c, lin&& l){
    return l *= constant<long double>(c);
}

lin operator/(bool c, lin&& l){
    return l /= constant<bool>(c);
}

lin operator/(int c, lin&& l){
    return l /= constant<int>(c);
}

lin operator/(short c, lin&& l){
    return l /= constant<short>(c);
}

lin operator/(float c, lin&& l){
    return l /= constant<float>(c);
}

lin operator/(double c, lin&& l){
    return l /= constant<double>(c);
}

lin operator/(long double c, lin&& l){
    return l /= constant<long double>(c);
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



//template bexpr operator+(expr const&, constant<bool> const&);
//template bexpr operator/(expr const&, constant<short> const&);
//template bexpr operator/(expr const&, constant<int> const&);
//template bexpr operator/(expr const&, constant<float> const&);
//template bexpr operator/(expr const&, constant<double> const&);
//template bexpr operator/(expr const&, constant<long double> const&);
//template bexpr operator/(constant<bool> const&, expr const&);
//template bexpr operator/(constant<short> const&, expr const&);
//template bexpr operator/(constant<int> const&, expr const&);
//template bexpr operator/(constant<float> const&, expr const&);
//template bexpr operator/(constant<double> const&, expr const&);
//template bexpr operator/(constant<long double> const&, expr const&);
//
//
//
//template bexpr operator/(expr const&, constant<bool> const&);
//template bexpr operator/(expr const&, constant<short> const&);
//template bexpr operator/(expr const&, constant<int> const&);
//template bexpr operator/(expr const&, constant<float> const&);
//template bexpr operator/(expr const&, constant<double> const&);
//template bexpr operator/(expr const&, constant<long double> const&);
//template bexpr operator/(constant<bool> const&, expr const&);
//template bexpr operator/(constant<short> const&, expr const&);
//template bexpr operator/(constant<int> const&, expr const&);
//template bexpr operator/(constant<float> const&, expr const&);
//template bexpr operator/(constant<double> const&, expr const&);
//template bexpr operator/(constant<long double> const&, expr const&);
