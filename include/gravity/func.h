//
//  func.h
//  Gravity
//
//  Created by Hijazi, Hassan on 21/10/16.
//
//

#ifndef func_h
#define func_h

#include <Gravity/var.h>
#include <stdio.h>
#include <map>
#include <list>
#include <set>

using namespace std;

void reverse_sign(constant_* c); /**< Reverses the sign of the constant. */
constant_* copy(const constant_* c2); /**< Copy c2 into a new constant_* detecting the right class, i.e., constant<>, param<>, uexpr or bexpr. */
bool equals(const constant_* c1, const constant_* c2);
double eval(ind i, const constant_* c1);
double eval(const constant_* c1);


/** Backbone class for function */
class func_ : public constant_{
private:
    FType                   _ftype = const_;
public:
    NType                   _return_type;
    bool                    _convex = true;
    set<var_*>*             _vars;/**< Set of variables appearing in this function **/
    
    virtual ~func_(){};
    pair<ind,func_*> operator[](ind i);
    bool is_convex() const{
        return (_convex);
    };
    
    bool is_constant() const{
        return (_ftype==const_ || (_ftype==lin_ && _vars->empty()));
    }
    
    bool is_linear() const{
        return (_ftype==lin_);
    };
    
    bool is_quadratic() const{
        return (_ftype==quad_);
    };
    
    bool is_polynome() const{
        return (_ftype==pol_);
    };
    
    bool is_nonlinear() const{
        return (_ftype==nlin_);
    };
    
    FType get_ftype() const { return _ftype;}
    
     void set_ftype(FType t){ _ftype = t;}
    
};


/** A function can return a bool, a short, an int, a float or a double or any user specified number type. The return type is set by default to int. */
template<typename type = int>
class func: public func_{
protected:
    vector<type>                _vals;
    
    
public:
    
    func(){
        update_type();
    };
    
    ~func(){
        
    };
    
   
    
    void update_type() {
        set_ftype(const_);
        if(typeid(type)==typeid(bool)){
            _return_type = binary_;
            return;
        }
        if(typeid(type)==typeid(short)) {
            _return_type = short_;
            return;
        }
        if(typeid(type)==typeid(int)) {
            _return_type = integer_;
            return;
        }
        if(typeid(type)==typeid(float)) {
            _return_type = float_;
            return;
        }
        if(typeid(type)==typeid(double)) {
            _return_type = double_;
            return;
        }
        if(typeid(type)==typeid(long double)) {
            _return_type = long_;
            return;
        }
        throw invalid_argument("Unsupported type");
    }

    
    /** Querries */
    
    

};


//typedef pair<map<ind, ind>*, var_*> ind_var; /**< pair <map, var> to store indices of a given variable */

/** A class to represent a linear term, e.g. 2x. */
class lterm{
    
public:
    constant_*              _coef;
    param_*                 _p;
    bool                    _sign = true; /**< True if +, flase if - */
    
    lterm(){
        _coef = nullptr;
        _p = nullptr;
    }
    
    lterm(const lterm& t){
        _coef = copy(t._coef);
        _p = (param_*)copy(t._p);
        _sign = t._sign;
    };
    
    lterm(lterm&& t){
        _coef = t._coef;
        t._coef = nullptr;
        _p = t._p;
        t._p = nullptr;
        _sign = t._sign;
    };
    
    
    lterm(param_* p){
        _coef = new constant<int>(1);
        _p = p;
    };
    
    
    lterm(constant_* coef, param_* p){
        _coef = coef;
        _p = p;        
    };
    
    lterm(bool sign, param_* p){
        _coef = new constant<int>(1);
        _p = p;
        _sign = sign;
    };
    
    lterm(bool sign, constant_* coef, param_* p){
        _coef = coef;
        _p = p;
        _sign = sign;
    };
    
    void reverse_sign() {
        _sign = ! _sign;
    }
    
    ~lterm(){
        delete _coef;
        delete _p;
    };
    
    bool operator==(const lterm& l) const;
    bool operator!=(const lterm& l) const;
    
    lterm& operator=(const lterm& l);
    lterm& operator=(lterm&& l);
    
    template<typename T> lterm& operator+=(const constant<T>& c);
    
    template <typename T> lterm& operator-=(const constant<T>& c){
        return *this += constant<T>(c.eval()*-1);
    }
    
    template<typename T> lterm& operator*=(const constant<T>& c);
    template<typename T> lterm& operator/=(const constant<T>& c);
    
    template<typename T> lterm& operator+=(const param<T>& p);
    template<typename T> lterm& operator-=(const param<T>& p);
    template<typename T> lterm& operator*=(const param<T>& p);
    template<typename T> lterm& operator/=(const param<T>& p);
    
};



/** A linear function */
class lin: public func_{
protected:

    constant_*                 _cst;/**< Constant part of the linear function */
    map<string, lterm>*        _lterms;

public:
    lin(){
        func_::set_ftype(lin_);
        set_type(func_c);
        _convex = true;
        _lterms = new map<string, lterm>();
        _vars = new set<var_*>;
        _cst = new constant<int>(0);
    };

    lin(constant_* coef, param_* p);
    
    bool insert(bool sign, constant_* coef, param_* p);/**< Adds coef*p to the linear function. Returns true if added new term, false if only updated coef of p */
    
    lin(param_* p);
    
    lin(const lin& l){
        func_::set_ftype(lin_);
        set_type(func_c);
        _convex = l._convex;
        _lterms = new map<string, lterm>();
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
        _cst = copy(l._cst);
    }
    
    lin(lin&& l){
        func_::set_ftype(lin_);
        set_type(func_c);
        _convex = l._convex;
        _lterms = l._lterms;
        l._lterms = nullptr;
        _vars = l._vars;
        l._vars = nullptr;
        _cst = l._cst;
        l._cst = nullptr;
    }
    
    void reset(){
        func_::set_ftype(lin_);
        set_type(func_c);
        _convex = true;
        delete _lterms;
        _lterms = new map<string, lterm>();
        delete _vars;
        _vars = new set<var_*>;
        delete _cst;
        _cst = new constant<int>(0);
        
    };
    
    ~lin();
    
    void reverse_sign(); /*<< Reverse the sign of all linear terms + constant in the function */
    
    lin& operator=(const lin& l);
    
    lin& operator=(lin&& l);
    
    bool operator==(const lin& l) const;
    bool operator!=(const lin& l) const;
    bool is_constant() const{
        return (_lterms->empty());
    }

    constant_* get_cst() {
        return _cst;
    }
    
    lin& operator+=(const lin& l);
    lin& operator-=(const lin& l);
    
    
    lin& operator*=(bool n){
        return *this *= constant<bool>(n);
    };
    lin& operator*=(short n){
        return *this *= constant<short>(n);
    };
    lin& operator*=(int n){
        return *this *= constant<int>(n);
    };
    lin& operator*=(float n){
        return *this *= constant<float>(n);
    };
    lin& operator*=(double n){
        return *this *= constant<double>(n);
    };
    lin& operator*=(long double n){
        return *this *= constant<long double>(n);
    };
    
    template<typename T> lin& operator+=(const constant<T>& c){
        if (c.is_zero()) {
            return *this;
        }
        _cst = add(_cst, c);
        return *this;
    }
    template<typename T> lin& operator-=(const constant<T>& c){
        if (c.is_zero()) {
            return *this;
        }
        _cst = substract(_cst, c);
        return *this;
    }
    
    template<typename T> lin& operator*=(const constant<T>& c){
        if (c.is_unit()) {
            return *this;
        }
        if (c.is_zero()) {
            reset();
            return *this;
        }
        _cst = multiply(_cst, c);
        
        for (auto &pair:*_lterms) {
            pair.second *= c;
        }
        return *this;
    }

    template<typename T> lin& operator/=(const constant<T>& c){
        if (c.is_zero()) {
            throw invalid_argument("\nERROT dividing linear function by zero!\n");
            return *this;
        }
        if (c.is_unit()) {
            return *this;
        }
        _cst = divide(_cst, c);
        
        for (auto &pair:*_lterms) {
            pair.second /= c;
        }
        return *this;
    }

    template<typename T> lin& operator+=(const param<T>& p){
        
        auto pair_it = _lterms->find(p.get_name());
        if (pair_it == _lterms->end()) {
            auto p_new = (param_*)copy((constant_*)&p);
            _lterms->insert(make_pair<>(p_new->get_name(), lterm(p_new)));
            if (p.is_var()) {
                _vars->insert((var_*)p_new);
            }
        }
        else {
            if (pair_it->second._coef->is_neg_unit()) {
                _lterms->erase(pair_it);
            }
            else{
                pair_it->second += constant<int>(1);
            }
        }
        return *this;
    }
    
    template<typename T> lin& operator-=(const param<T>& p){
        
        auto pair_it = _lterms->find(p.get_name());
        if (pair_it == _lterms->end()) {
            auto p_new = (param_*)copy((constant_*)&p);
            _lterms->insert(make_pair<>(p_new->get_name(), lterm(new constant<int>(-1), p_new)));
            if (p.is_var()) {
                _vars->insert((var_*)p_new);
            }
        }
        else {
            if (pair_it->second._coef->is_unit()) {
                _lterms->erase(pair_it);
            }
            else{
                pair_it->second -= constant<int>(1);
            }
        }
        return *this;
    }

    template<typename T> lin& operator*=(const param<T>& p){
        if (!_cst->is_zero()) {
            auto pair_it = _lterms->find(p.get_name());
            if (pair_it == _lterms->end()) {
                auto p_new = (param_*)copy((constant_*)&p);
                _lterms->insert(make_pair<>(p_new->get_name(), lterm(_cst, p_new)));
            }
            else {
                switch (_cst->get_type()) {
                    case binary_c: {
                        pair_it->second += (*(constant<bool>*)_cst);
                        break;
                    }
                    case short_c: {
                        pair_it->second += (*(constant<short>*)_cst);
                        break;
                    }
                    case integer_c: {
                        pair_it->second += (*(constant<int>*)_cst);
                        break;
                    }
                    case float_c: {
                        pair_it->second += (*(constant<float>*)_cst);
                        break;
                    }
                    case double_c: {
                        pair_it->second += (*(constant<double>*)_cst);
                        break;
                    }
                    case long_c: {
                        pair_it->second += (*(constant<long double>*)_cst);
                        break;
                    }
                    default:
                        break;
                }
            }
            _cst = new constant<int>(0);
        }
        for (auto &pair:*_lterms) {
            pair.second *= p;
        }
        return *this;
    }
//    template<typename T> lin& operator+=(const constant<T>& c);
//    template<typename T> lin& operator-=(const constant<T>& c);
//    template<typename T> lin& operator+=(const param<T>& c);
//    template<typename T> lin& operator-=(const param<T>& c);
//    template<typename T> lin& operator*=(const constant<T>& c);
//    template<typename T> lin& operator/=(const constant<T>& c);
//    template<typename T> lin& operator*=(const param<T>& c);
//    template<typename T> lin& operator/=(const param<T>& c);
    
    void const print(bool endline=true);

};



/** A quadratic function */
class quad: public lin{
public:
    quad();
    ~quad();
    
};

//typedef pair<ind_var*, ind> pvar; /**< pair <ind_var, exponent> to represent x^p */
//typedef pair<constant_*, map<ind, pvar*>*> pterm; /**< pair <front coeffcient, list of pvars sorted by var ids> to represent a polynomial term, e.g., 2x^2y^3z */
//


/** A polynomial function */
class polynome: public quad{
    
public:
//    shared_ptr<vector<type>>    _lb; /**< Lower Bound */
//    shared_ptr<vector<type>>    _ub; /**< Upper Bound */
//    
//    
//    /* Constructors */
//    
//    //@{
//    /** Unbounded variable constructor */
//    var();
//    var(const char* name);
//    var(const var<type>& v);
//    //@}
//    
//    
//    //@{
//    /** Bounded variable constructor */
//    var(const char* name, type lb, type ub);
//    //@}
//    
//    
//    /* Querries */
//    
//    type    get_lb() const;
//    type    get_ub() const;
//
//    
//    bool is_bounded_above(int i = 0) const{
//        return (_ub!=nullptr &&  i < _ub->size() && _ub->at(i)!=numeric_limits<type>::max());
//    };
//
//    bool is_bounded_below(int i = 0) const{
//        return (_lb!=nullptr &&  i < _lb->size() && _lb->at(i)!=numeric_limits<type>::min());
//    };
//    
//    bool is_constant(int i=0) const{
//        return (is_bounded_below() && is_bounded_above() && _lb->at(i)==_ub->at(i));
//    };
//    
//    
//    /* Modifiers */
//    void    reserve(int size);
//    
//    void    add_bounds(type lb, type ub);
//    void    add_lb_only(type v); /**< Adds a new lower bound and infinity in the corresponding upper bound*/
//    void    add_ub_only(type v); /**< Adds a new upper bound and -infinity in the corresponding lower bound*/
//    
//    void    set_lb(int i, type v);
//    void    set_ub(int i, type v);
//
//    
//    /* Operators */
//    bool operator==(const var& v) const;
//    bool operator!=(const var& v) const;
//    
//    
//    /* Output */
//    void print(bool bouds=false) const;
    

    
};

constant_* add(constant_* c1, constant_* c2);
constant_* substract(constant_* c1, constant_* c2);

template<typename T> constant_* add(constant_* c1, const constant<T>& c2){ /**< adds c2 to c1, updates its type and returns the result **/
    switch (c1->get_type()) {
        case binary_c: {
            if (c2.is_binary() ) {
                *(constant<bool>*)c1 += c2.eval();
            }
            else {
                auto val = ((constant<bool>*)c1)->eval();
                delete c1;
                c1 = new constant<T>(c2 + val);
            }
            return c1;
            break;
        }
        case short_c: {
            if (c2.get_type() <= short_c) {
                *((constant<short>*)c1) += c2.eval();
            }
            else {
                auto val = ((constant<short>*)c1)->eval();
                delete c1;
                c1 = new constant<T>(c2 + val);
            }
            break;
        }
        case integer_c: {
            if (c2.get_type() <= integer_c) {
                *((constant<int>*)c1) += c2.eval();
            }
            else {
                auto val = ((constant<int>*)c1)->eval();
                delete c1;
                c1 = new constant<T>(c2 + val);
            }
            break;
        }
        case float_c: {
            if (c2.get_type() <= float_c) {
                *((constant<float>*)c1) += c2.eval();
            }
            else {
                auto val = ((constant<float>*)c1)->eval();
                delete c1;
                c1 = new constant<T>(c2 + val);
            }
            break;
        }
        case double_c: {
            if (c2.get_type() <= double_c) {
                *((constant<double>*)c1) += c2.eval();
            }
            else {
                auto val = ((constant<double>*)c1)->eval();
                delete c1;
                c1 = new constant<T>(c2 + val);
            }
            break;
        }
        case long_c: {
            if (c2.get_type() <= long_c) {
                *((constant<long double>*)c1) += c2.eval();
            }
            else {
                auto val = ((constant<long double>*)c1)->eval();
                delete c1;
                c1 = new constant<T>(c2 + val);
            }
            break;
        }
        case par_c:{
            auto pc1 = (param_*)(c1);
            auto l = new lin(pc1);
            *l += c2;
            c1 =(constant_*)(l);
            return c1;
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
            switch (((func_*)c1)->get_ftype()) {
                case lin_: {
                    auto res = new lin(*(lin*)c1 + c2);
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
    return c1;
}


template<typename T> constant_* add(constant_* c1, const param<T>& c2){ /**< adds c2 to c1, updates its type and returns the result **/
    switch (c1->get_type()) {
        case binary_c: {
            auto val = ((constant<bool>*)c1);
            delete c1;
            auto l = new lin(new param<T>(c2));
            *l += *val;
            c1 = (constant_*)(l);
            return c1;
            break;
        }
        case short_c: {
            auto val = ((constant<short>*)c1);
            delete c1;
            auto l = new lin(new param<T>(c2));
            *l += *val;
            c1 = (constant_*)(l);
            return c1;
            break;
        }
        case integer_c: {
            auto val = ((constant<int>*)c1);
            delete c1;
            auto l = new lin(new param<T>(c2));
            *l += *val;
            c1 = (constant_*)(l);
            return c1;
            break;
        }
        case float_c: {
            auto val = ((constant<float>*)c1);
            delete c1;
            auto l = new lin(new param<T>(c2));
            *l += *val;
            c1 = (constant_*)(l);
            return c1;
            break;
        }
        case double_c: {
            auto val = ((constant<double>*)c1);
            delete c1;
            auto l = new lin(new param<T>(c2));
            *l += *val;
            c1 = (constant_*)(l);
            return c1;
            break;
        }
        case long_c: {
            auto val = ((constant<long double>*)c1);
            delete c1;
            auto l = new lin(new param<T>(c2));
            *l += *val;
            c1 = (constant_*)(l);
            return c1;
            break;
        }
        case par_c:{
            auto pc1 = (param_*)(c1);
            auto l = new lin(pc1);
            *l += c2;
            c1 =(constant_*)(l);
            return c1;
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
            switch (((func_*)c1)->get_ftype()) {
                case lin_: {
                    auto res = new lin(*(lin*)c1 + c2);
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
    return c1;
}


template<typename T> constant_* substract(constant_* c1, const param<T>& c2){ /**< adds c2 to c1, updates its type and returns the result **/
    switch (c1->get_type()) {
        case binary_c: {
            auto val = ((constant<bool>*)c1);
            delete c1;
            auto l = new lin(new param<T>(c2));
            *l -= *val;
            c1 = (constant_*)(l);
            return c1;
            break;
        }
        case short_c: {
            auto val = ((constant<short>*)c1);
            delete c1;
            auto l = new lin(new param<T>(c2));
            *l -= *val;
            c1 = (constant_*)(l);
            return c1;
            break;
        }
        case integer_c: {
            auto val = ((constant<int>*)c1);
            delete c1;
            auto l = new lin(new param<T>(c2));
            *l -= *val;
            c1 = (constant_*)(l);
            return c1;
            break;
        }
        case float_c: {
            auto val = ((constant<float>*)c1);
            delete c1;
            auto l = new lin(new param<T>(c2));
            *l -= *val;
            c1 = (constant_*)(l);
            return c1;
            break;
        }
        case double_c: {
            auto val = ((constant<double>*)c1);
            delete c1;
            auto l = new lin(new param<T>(c2));
            *l -= *val;
            c1 = (constant_*)(l);
            return c1;
            break;
        }
        case long_c: {
            auto val = ((constant<long double>*)c1);
            delete c1;
            auto l = new lin(new param<T>(c2));
            *l -= *val;
            c1 = (constant_*)(l);
            return c1;
            break;
        }
        case par_c:{
            auto pc1 = (param_*)(c1);
            auto l = new lin(pc1);
            *l -= c2;
            c1 =(constant_*)(l);
            return c1;
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
            switch (((func_*)c1)->get_ftype()) {
                case lin_: {
                    auto res = new lin(*(lin*)c1 - c2);
                    delete c1;
                    if (res->is_constant()) {
                        c1 = copy(res->get_cst());
                        delete res;
                        return c1;
                    }
                    else {
                        return c1 = res;
                    }
                    break;                }
                default:
                    break;
            }
        }
        default:
            break;
    }
    return c1;
}



template<typename T> constant_* multiply(constant_* c1, const param<T>& c2){ /**< adds c2 to c1, updates its type and returns the result **/
    switch (c1->get_type()) {
        case binary_c: {
            auto val = ((constant<bool>*)c1);
            delete c1;
            auto l = new lin(new param<T>(c2));
            *l *= *val;
            c1 = (constant_*)(l);
            return c1;
            break;
        }
        case short_c: {
            auto val = ((constant<short>*)c1);
            delete c1;
            auto l = new lin(new param<T>(c2));
            *l *= *val;
            c1 = (constant_*)(l);
            return c1;
            break;
        }
        case integer_c: {
            auto val = ((constant<int>*)c1);
            delete c1;
            auto l = new lin(new param<T>(c2));
            *l *= *val;
            c1 = (constant_*)(l);
            return c1;
            break;
        }
        case float_c: {
            auto val = ((constant<float>*)c1);
            delete c1;
            auto l = new lin(new param<T>(c2));
            *l *= *val;
            c1 = (constant_*)(l);
            return c1;
            break;
        }
        case double_c: {
            auto val = ((constant<double>*)c1);
            delete c1;
            auto l = new lin(new param<T>(c2));
            *l *= *val;
            c1 = (constant_*)(l);
            return c1;
            break;
        }
        case long_c: {
            auto val = ((constant<long double>*)c1);
            delete c1;
            auto l = new lin(new param<T>(c2));
            *l *= *val;
            c1 = (constant_*)(l);
            return c1;
            break;
        }
        case par_c:{
            //            auto pc1 = (param_*)(c1);
            //            auto l = new lin(pc1);
            //            *l *= c2;
            //            c1 =(constant_*)(l);
            //            return c1;
            cerr << "Unsupported yet;\n";
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
            switch (((func_*)c1)->get_ftype()) {
                case lin_: {
                    cerr << "Unsupported yet;\n";
                    //                    auto res = new lin(*(lin*)c1 + c2);
                    //                    delete c1;
                    //                    return c1 = res;
                    break;
                }
                default:
                    break;
            }
        }
            
        default:
            break;
    }
    return c1;
}




template<typename T> constant_* substract(constant_* c1, const constant<T>& c2){ /**< adds c2 to c1, updates its type and returns the result **/
    switch (c1->get_type()) {
        case binary_c: {
            if (c2.is_binary() ) {
                *(constant<bool>*)c1 -= c2.eval();
            }
            else {
                auto val = ((constant<bool>*)c1)->eval();
                delete c1;
                c1 = new constant<T>(c2 - val);
            }
            return c1;
            break;
        }
        case short_c: {
            if (c2.get_type() <= short_c) {
                *((constant<short>*)c1) -= c2.eval();
            }
            else {
                auto val = ((constant<short>*)c1)->eval();
                delete c1;
                c1 = new constant<T>(c2 - val);
            }
            break;
        }
        case integer_c: {
            if (c2.get_type() <= integer_c) {
                *((constant<int>*)c1) -= c2.eval();
            }
            else {
                auto val = ((constant<int>*)c1)->eval();
                delete c1;
                c1 = new constant<T>(c2 - val);
            }
            break;
        }
        case float_c: {
            if (c2.get_type() <= float_c) {
                *((constant<float>*)c1) -= c2.eval();
            }
            else {
                auto val = ((constant<float>*)c1)->eval();
                delete c1;
                c1 = new constant<T>(c2 - val);
            }
            break;
        }
        case double_c: {
            if (c2.get_type() <= double_c) {
                *((constant<double>*)c1) -= c2.eval();
            }
            else {
                auto val = ((constant<double>*)c1)->eval();
                delete c1;
                c1 = new constant<T>(c2 - val);
            }
            break;
        }
        case long_c: {
            if (c2.get_type() <= long_c) {
                *((constant<long double>*)c1) -= c2.eval();
            }
            else {
                auto val = ((constant<long double>*)c1)->eval();
                delete c1;
                c1 = new constant<T>(c2 - val);
            }
            break;
        }
        case par_c:{
            auto pc1 = (param_*)(c1);
            auto l = new lin(pc1);
            *l -= c2;
            c1 =(constant_*)(l);
            return c1;
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
            switch (((func_*)c1)->get_ftype()) {
                case lin_: {
                    auto res = new lin(*(lin*)c1 - c2);
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
    return c1;
}

template<typename T> constant_* multiply(constant_* c1, const constant<T>& c2){ /**< adds c2 to c1, updates its type and returns the result **/
    if (eval(c1)==0) {
        return c1;
    }
    if (c2.eval()==1) {
        return c1;
    }
    if (c2.eval()==0) {
        delete c1;
        c1 = new constant<int>(0);
        return c1;
    }
    switch (c1->get_type()) {
        case binary_c: {
            if (c2.is_binary() ) {
                *(constant<bool>*)c1 *= c2.eval();
            }
            else {
                auto val = ((constant<bool>*)c1)->eval();
                delete c1;
                c1 = new constant<T>(c2 * val);
            }
            return c1;
            break;
        }
        case short_c: {
            if (c2.get_type() <= short_c) {
                *((constant<short>*)c1) *= c2.eval();
            }
            else {
                auto val = ((constant<short>*)c1)->eval();
                delete c1;
                c1 = new constant<T>(c2 * val);
            }
            break;
        }
        case integer_c: {
            if (c2.get_type() <= integer_c) {
                *((constant<int>*)c1) *= c2.eval();
            }
            else {
                auto val = ((constant<int>*)c1)->eval();
                delete c1;
                c1 = new constant<T>(c2 * val);
            }
            break;
        }
        case float_c: {
            if (c2.get_type() <= float_c) {
                *((constant<float>*)c1) *= c2.eval();
            }
            else {
                auto val = ((constant<float>*)c1)->eval();
                delete c1;
                c1 = new constant<T>(c2 * val);
            }
            break;
        }
        case double_c: {
            if (c2.get_type() <= double_c) {
                *((constant<double>*)c1) *= c2.eval();
            }
            else {
                auto val = ((constant<double>*)c1)->eval();
                delete c1;
                c1 = new constant<T>(c2 * val);
            }
            break;
        }
        case long_c: {
            if (c2.get_type() <= long_c) {
                *((constant<long double>*)c1) *= c2.eval();
            }
            else {
                auto val = ((constant<long double>*)c1)->eval();
                delete c1;
                c1 = new constant<T>(c2 * val);
            }
            break;
        }
        case par_c:{
            auto pc1 = (param_*)(c1);
            auto l = new lin(pc1);
            *l *= c2;
            c1 =(constant_*)(l);
            return c1;
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
            switch (((func_*)c1)->get_ftype()) {
                case lin_: {
                    auto res = new lin(*(lin*)c1 * c2);
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
    return c1;
}

template<typename T> constant_* divide(constant_* c1, const constant<T>& c2){ /**< adds c2 to c1, updates its type and returns the result **/
    if (c2.eval()==0) {
        throw invalid_argument("dividing by zero!\n");
    }
    switch (c1->get_type()) {
        case binary_c: {
            if (c2.is_binary() ) {
                *(constant<bool>*)c1 /= c2.eval();
            }
            else {
                auto val = ((constant<bool>*)c1)->eval();
                delete c1;
                c1 = new constant<T>(c2 / val);
            }
            return c1;
            break;
        }
        case short_c: {
            if (c2.get_type() <= short_c) {
                *((constant<short>*)c1) /= c2.eval();
            }
            else {
                auto val = ((constant<short>*)c1)->eval();
                delete c1;
                c1 = new constant<T>(c2 / val);
            }
            break;
        }
        case integer_c: {
            if (c2.get_type() <= integer_c) {
                *((constant<int>*)c1) /= c2.eval();
            }
            else {
                auto val = ((constant<int>*)c1)->eval();
                delete c1;
                c1 = new constant<T>(c2 / val);
            }
            break;
        }
        case float_c: {
            if (c2.get_type() <= float_c) {
                *((constant<float>*)c1) /= c2.eval();
            }
            else {
                auto val = ((constant<float>*)c1)->eval();
                delete c1;
                c1 = new constant<T>(c2 / val);
            }
            break;
        }
        case double_c: {
            if (c2.get_type() <= double_c) {
                *((constant<double>*)c1) /= c2.eval();
            }
            else {
                auto val = ((constant<double>*)c1)->eval();
                delete c1;
                c1 = new constant<T>(c2 / val);
            }
            break;
        }
        case long_c: {
            if (c2.get_type() <= long_c) {
                *((constant<long double>*)c1) /= c2.eval();
            }
            else {
                auto val = ((constant<long double>*)c1)->eval();
                delete c1;
                c1 = new constant<T>(c2 / val);
            }
            break;
        }
        case par_c:{
            auto pc1 = (param_*)(c1);
            auto l = new lin(pc1);
            *l /= c2;
            c1 =(constant_*)(l);
            return c1;
            break;
        }
        case uexp_c: {
            auto res = new bexpr(*(uexpr*)c1 / c2);
            delete c1;
            c1 = (constant_*)res;
            return c1;
            break;
        }
        case bexp_c: {
            auto res = new bexpr(*(bexpr*)c1 / c2);
            delete c1;
            c1 = (constant_*)res;
            return c1;
            break;
        }
        case func_c: {
            switch (((func_*)c1)->get_ftype()) {
                case lin_: {
                    cerr << "Unsupported yet;\n";
                    //                    auto res = new lin(*(lin*)c1 * 1/c2);
                    //                    delete c1;
                    //                    return c1 = res;
                    break;
                }
                default:
                    break;
            }
        }
        default:
            break;
    }
    return c1;
}


lin operator+(const lin& l1, const lin& l2);
lin operator+(lin&& l1, const lin& l2);

lin operator-(const lin& l1, const lin& l2);
lin operator-(lin&& l1, const lin& l2);



template<typename T1, typename T2> lin operator+(const param<T1>& v1, const param<T2>& v2){
    if (v1.is_var()) {
        auto vv1 = (var<T1>*)&v1;
        return lin(new var<T1>(*vv1))+=v2;
    }
    return lin(new param<T1>(v1))+=v2;
}

template<typename T> lin operator+(const constant<T>& c, const lin& l){
    return lin(l) += c;
}

template<typename T> lin operator+(const constant<T>& c, lin&& l){
    return l += c;
}

template<typename T> lin operator+(const lin& l, const constant<T>& c){
    return lin(l) += c;
}

template<typename T> lin operator+(lin&& l, const constant<T>& c){
    return l += c;
}

template<typename T> lin operator+(const param<T>& v, const lin& l){
    return lin(l)+=v;
}

template<typename T> lin operator+(const param<T>& v, lin&& l){
    return l+=v;
}

template<typename T> lin operator+(const lin& l, const param<T>& v){
    return lin(l)+=v;
}


template<typename T> lin operator+(lin&& l, const param<T>& v){
    return l+=v;
}


template<typename T1, typename T2> lin operator+(const param<T1>& v, const constant<T2>& c){
    return lin(c) += v;
}

template<typename T1, typename T2> lin operator+(const constant<T2>& c, const var<T1>& v){
    return lin(c) += v;
}


template<typename T1, typename T2> lin operator-(const param<T1>& v1, const param<T2>& v2){
    if (v1.is_var()) {
        auto vv1 = (var<T1>*)&v1;
        return lin(new var<T1>(*vv1))-=v2;
    }
    return lin(new param<T1>(v1))-=v2;
}

template<typename T> lin operator-(const constant<T>& c, const lin& l){
    return  (lin(l) *= -1) += c;
}

template<typename T> lin operator-(const constant<T>& c, lin&& l){
    return  (l *= -1) += c;
}

template<typename T> lin operator-(const lin& l, const constant<T>& c){
    return lin(l) -= c;
}

template<typename T> lin operator-(lin&& l, const constant<T>& c){
    return l -= c;
}

template<typename T> lin operator-(const param<T>& v, const lin& l){
    return  (lin(l) *= -1) += v;
}

template<typename T> lin operator-(const param<T>& v, lin&& l){
    return  (l *= -1) += v;
}

template<typename T> lin operator-(const lin& l, const param<T>& v){
    return lin(l)-=v;
}


template<typename T> lin operator-(lin&& l, const param<T>& v){
    return l-=v;
}


template<typename T1, typename T2> lin operator-(const param<T1>& v, const constant<T2>& c){
    return lin(v) -= c;
}

template<typename T1, typename T2> lin operator-(const constant<T2>& c, const param<T1>& v){
    return lin(c) -= v;
}


template<typename T1, typename T2> lin operator*(const param<T1>& v1, const param<T2>& v2){
    if (v1.is_var()) {
        auto vv1 = (var<T1>*)&v1;
        return lin(new param<T2>(v2), new var<T1>(*vv1));
    }
    else{
        auto vv2 = (var<T2>*)&v2;
        return lin(new param<T1>(v1), new var<T2>(*vv2));
    }
    throw invalid_argument("multpilying two variables does not give a linear function");
}

template<typename T1, typename T2> quad operator*(const var<T1>& v1, const var<T2>& v2){
    return quad();
}


template<typename T1, typename T2> lin operator*(const param<T1>& v, const constant<T2>& c){
    if (v.is_var()) {
        auto vv = (var<T1>*)&v;
        return lin(new constant<T2>(c), new var<T1>(*vv));
    }
    else{
        return lin(new constant<T2>(c), new param<T1>(v));
    }
}

template<typename T1, typename T2> lin operator*(const constant<T2>& c, const param<T1>& v){
    if (v.is_var()) {
        auto vv = (var<T1>*)&v;
        return lin(new constant<T2>(c), new var<T1>(*vv));
    }
    else{
        return lin(new constant<T2>(c), new param<T1>(v));
    }

}

template<typename T> lin operator*(const int c, const param<T>& v){
    if (v.is_var()) {
        auto vv = (var<T>*)&v;
        return lin(new constant<int>(c), new var<T>(*vv));
    }
    return lin(new constant<int>(c), new param<T>(v));
}

template<typename T> lin operator*(const param<T>& v, const int c){
    if (v.is_var()) {
        auto vv = (var<T>*)&v;
        return lin(new constant<int>(c), new var<T>(*vv));
    }
    return lin(new constant<int>(c), new param<T>(v));
}


template<typename T> lin operator*(const double c, const param<T>& v){
    if (v.is_var()) {
        auto vv = (var<T>*)&v;
        return lin(new constant<double>(c), new var<T>(*vv));
    }
    return lin(new constant<double>(c), new param<T>(v));
}


template<typename T> lin operator*(const param<T>& v, const double c){
    if (v.is_var()) {
        auto vv = (var<T>*)&v;
        return lin(new constant<double>(c), new var<T>(*vv));
    }
    return lin(new constant<double>(c), new param<T>(v));
}


template<typename T> lin operator*(const float c, const param<T>& v){
    if (v.is_var()) {
        auto vv = (var<T>*)&v;
        return lin(new constant<float>(c), new var<T>(*vv));
    }
    return lin(new constant<float>(c), new param<T>(v));
}

template<typename T> lin operator*(const param<T>& v, const float c){
    if (v.is_var()) {
        auto vv = (var<T>*)&v;
        return lin(new constant<float>(c), new var<T>(*vv));
    }
    return lin(new constant<float>(c), new param<T>(v));
}


template <typename T> lin operator*(const lin& l, const constant<T>& c){
    return lin(l) *= c;
}

//template <typename T> lin operator*(const lin& l, const constant<T>& c);


template<typename T1, typename T2> lin operator*(const lin& l, const constant<T2>& c){
    return lin(l) *= c;
}

template<typename T1, typename T2> lin operator*(lin&& l, const constant<T2>& c){
    return l *= c;
}

template<typename T1, typename T2> lin operator*(const constant<T2>& c, const lin& l){
    return lin(l) *= c;
}

template<typename T1, typename T2> lin operator*(const constant<T2>& c, lin&& l){
    return l *= c;
}


template<typename T1, typename T2> lin operator/(const var<T1>& v, const constant<T2>& c);




template<typename T> quad operator+(const constant<T>& c, const quad& q);
template<typename T> quad operator+(const quad& q, const constant<T>& c);
template<typename T> quad operator+(const var<T>& v, const quad& q);
template<typename T> quad operator+(const quad& q, const var<T>& v);
template<typename T> quad operator+(const constant<T>& c, const quad& q);


template<typename T> quad operator-(const constant<T>& c, const quad& q);
template<typename T> quad operator-(const quad& q, const constant<T>& c);
template<typename T> quad operator-(const var<T>& v, const quad& q);
template<typename T> quad operator-(const quad& q, const var<T>& v);
template<typename T> quad operator-(const constant<T>& c, const quad& q);


template<typename T> quad operator*(const quad& q, const constant<T>& c);
template<typename T> quad operator*(const constant<T>& c, const quad& q);
template<typename T> quad operator*(const lin& l, const var<T>& v);
template<typename T> quad operator*(const var<T>& v, const lin& l);
quad operator*(const lin& l1, const lin& l2);

template<typename T> quad operator/(const quad& q, const constant<T>& c);


template<typename T> polynome operator+(const constant<T>& c, const polynome& p);
template<typename T> polynome operator+(const polynome& p, const constant<T>& c);
template<typename T> polynome operator+(const var<T>& v, const polynome& p);
template<typename T> polynome operator+(const polynome& p, const var<T>& v);
template<typename T> polynome operator+(const constant<T>& c, const polynome& p);


template<typename T> polynome operator-(const constant<T>& c, const polynome& p);
template<typename T> polynome operator-(const polynome& p, const constant<T>& c);
template<typename T> polynome operator-(const var<T>& v, const polynome& p);
template<typename T> polynome operator-(const polynome& p, const var<T>& v);
template<typename T> polynome operator-(const constant<T>& c, const polynome& p);


template<typename T> polynome operator*(const polynome& p, const constant<T>& c);
template<typename T> polynome operator*(const constant<T>& c, const polynome& p);
template<typename T> polynome operator*(const quad& q, const var<T>& v);
template<typename T> polynome operator*(const var<T>& v, const quad& q);
polynome operator*(const quad& q1, const quad& q2);


template<typename T> polynome operator/(const polynome& p, const constant<T>& c);


lin operator+(const lin& l, bool c);
lin operator+(const lin& l, int c);
lin operator+(const lin& l, short c);
lin operator+(const lin& l, float c);
lin operator+(const lin& l, double c);
lin operator+(const lin& l, long double c);
lin operator-(const lin& l, bool c);
lin operator-(const lin& l, int c);
lin operator-(const lin& l, short c);
lin operator-(const lin& l, float c);
lin operator-(const lin& l, double c);
lin operator-(const lin& l, long double c);
lin operator*(const lin& l, bool c);
lin operator*(const lin& l, int c);
lin operator*(const lin& l, short c);
lin operator*(const lin& l, float c);
lin operator*(const lin& l, double c);
lin operator*(const lin& l, long double c);
lin operator/(const lin& l, bool c);
lin operator/(const lin& l, int c);
lin operator/(const lin& l, short c);
lin operator/(const lin& l, float c);
lin operator/(const lin& l, double c);
lin operator/(const lin& l, long double c);
lin operator+(lin&& l, bool c);
lin operator+(lin&& l, int c);
lin operator+(lin&& l, short c);
lin operator+(lin&& l, float c);
lin operator+(lin&& l, double c);
lin operator+(lin&& l, long double c);
lin operator-(lin&& l, bool c);
lin operator-(lin&& l, int c);
lin operator-(lin&& l, short c);
lin operator-(lin&& l, float c);
lin operator-(lin&& l, double c);
lin operator-(lin&& l, long double c);
lin operator*(lin&& l, bool c);
lin operator*(lin&& l, int c);
lin operator*(lin&& l, short c);
lin operator*(lin&& l, float c);
lin operator*(lin&& l, double c);
lin operator*(lin&& l, long double c);
lin operator/(lin&& l, bool c);
lin operator/(lin&& l, int c);
lin operator/(lin&& l, short c);
lin operator/(lin&& l, float c);
lin operator/(lin&& l, double c);
lin operator/(lin&& l, long double c);
lin operator+(bool c, const lin& l);
lin operator+(int c, const lin& l);
lin operator+(short c, const lin& l);
lin operator+(float c, const lin& l);
lin operator+(double c, const lin& l);
lin operator+(long double c, const lin& l);
lin operator-(bool c, const lin& l);
lin operator-(int c, const lin& l);
lin operator-(short c, const lin& l);
lin operator-(float c, const lin& l);
lin operator-(double c, const lin& l);
lin operator-(long double c, const lin& l);
lin operator*(bool c, const lin& l);
lin operator*(int c, const lin& l);
lin operator*(short c, const lin& l);
lin operator*(float c, const lin& l);
lin operator*(double c, const lin& l);
lin operator*(long double c, const lin& l);
lin operator/(bool c, const lin& l);
lin operator/(int c, const lin& l);
lin operator/(short c, const lin& l);
lin operator/(float c, const lin& l);
lin operator/(double c, const lin& l);
lin operator/(long double c, const lin& l);
lin operator+(bool c, lin&& l);
lin operator+(int c, lin&& l);
lin operator+(short c, lin&& l);
lin operator+(float c, lin&& l);
lin operator+(double c, lin&& l);
lin operator+(long double c, lin&& l);
lin operator-(bool c, lin&& l);
lin operator-(int c, lin&& l);
lin operator-(short c, lin&& l);
lin operator-(float c, lin&& l);
lin operator-(double c, lin&& l);
lin operator-(long double c, lin&& l);
lin operator*(bool c, lin&& l);
lin operator*(int c, lin&& l);
lin operator*(short c, lin&& l);
lin operator*(float c, lin&& l);
lin operator*(double c, lin&& l);
lin operator*(long double c, lin&& l);
lin operator/(bool c, lin&& l);
lin operator/(int c, lin&& l);
lin operator/(short c, lin&& l);
lin operator/(float c, lin&& l);
lin operator/(double c, lin&& l);
lin operator/(long double c, lin&& l);




template<typename other_type> bexpr operator+(const other_type& c1, const expr& c2){
    bexpr res;
    res._otype = plus_;
    res._lson = copy((constant_*)&c1);
    res._rson =  copy((constant_*)&c2);
    return res;
}

template<typename other_type> bexpr operator+(const expr& c1, const other_type& c2){
    bexpr res;
    res._otype = plus_;
    res._lson = copy((constant_*)&c1);
    res._rson =  copy((constant_*)&c2);
    return res;
}


template<typename other_type> bexpr operator-(const other_type& c1, const expr& c2){
    bexpr res;
    res._otype = minus_;
    res._lson = copy((constant_*)&c1);
    res._rson =  copy((constant_*)&c2);
    return res;
}

template<typename other_type> bexpr operator-(const expr& c1, const other_type& c2){
    bexpr res;
    res._otype = minus_;
    res._lson = copy((constant_*)&c1);
    res._rson =  copy((constant_*)&c2);
    return res;
}


template<typename other_type> bexpr operator*(const other_type& c1, const expr& c2){
    bexpr res;
    res._otype = product_;
    res._lson = copy((constant_*)&c1);
    res._rson =  copy((constant_*)&c2);
    return res;
}

template<typename other_type> bexpr operator*(const expr& c1, const other_type& c2){
    bexpr res;
    res._otype = product_;
    res._lson = copy((constant_*)&c1);
    res._rson =  copy((constant_*)&c2);
    return res;
}



template<typename other_type> bexpr operator/(const other_type& c1, const expr& c2){
    bexpr res;
    res._otype = div_;
    res._lson = copy((constant_*)&c1);
    res._rson =  copy((constant_*)&c2);
    return res;
}
template<typename other_type> bexpr operator/(const expr& c1, const other_type& c2){
    bexpr res;
    res._otype = div_;
    res._lson = copy((constant_*)&c1);
    res._rson =  copy((constant_*)&c2);
    return res;
}
#endif /* func_h */

