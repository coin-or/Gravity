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
    FType                                  _ftype = const_;
public:
    NType                                  _return_type;
    bool                                   _convex = true;
    map<string, pair<param_*, int>>*       _vars;/**< map <variable name, <variable pointer, number of times it appears in function>>**/
    
    param_* get_var(string name){
        auto pair_it = _vars->find(name);
        if (pair_it==_vars->end()) {
            return nullptr;
        }
        else {
            return get<1>(*pair_it).first;
        }
    }
    
    void add_var(param_* p){/**< Inserts the variable in this function input list. WARNING: Assumes that p has not been added previousely!*/
        assert(_vars->count(p->get_name())==0);
        _vars->insert(make_pair<>(p->get_name(), make_pair<>(p, 1)));
    }
    
    int nb_occ(string name) const{/**< Returns the number of occurences the variable has in this function. */
        auto pair_it = _vars->find(name);
        if (pair_it==_vars->end()) {
            return 0;
        }
        else {
            return get<1>(*pair_it).second;
        }
    }
    
    void incr_occ(string str){/**< Increases the number of occurences the variable has in this function. */
        auto pair_it = _vars->find(str);
        if (pair_it==_vars->end()) {
            throw invalid_argument("Non-existing variable in function!\n");
        }
        else {
            get<1>(*pair_it).second++;
        }
    }
    
    void decr_occ(string str){/**< Decreases the number of occurences the variable has in this function. */
        auto pair_it = _vars->find(str);
        if (pair_it==_vars->end()) {
            return;
        }
        else {
            get<1>(*pair_it).second--;
            if (get<1>(*pair_it).second==0) {
                delete get<1>(*pair_it).first;
                _vars->erase(pair_it);
            }
        }
    }
    
    virtual ~func_(){};
    pair<ind,func_*> operator[](ind i);
    bool is_convex() const{
        return (_convex);
    };
    
    bool is_constant() const{
        return (_ftype==const_);
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
//        delete _p;
    };
    
    bool operator==(const lterm& l) const;
    bool operator!=(const lterm& l) const;
    
    lterm& operator=(const lterm& l);
    lterm& operator=(lterm&& l);
    
    
    template <typename T> lterm& operator+=(const constant<T>& c){
        if (c.is_zero()) {
            return *this;
        }
        _coef = add(_coef, c);
        return *this;
    }
    
    template <typename T> lterm& operator-=(const constant<T>& c){
        return *this += constant<T>(c.eval()*-1);
    }

    
    template <typename T> lterm& operator*=(const constant<T>& c){
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
        _coef = multiply(_coef, c);
        return *this;
    }
    
    template <typename T> lterm& operator/=(const constant<T>& c){
        if (c.is_zero()) {
            throw invalid_argument("\nERROT dividing lterm by zero!\n");
            return *this;
        }
        if (c.is_unit()) {
            return *this;
        }
        _coef = divide(_coef, c);
        return *this;
    }
    
    
    template <typename T> lterm& operator+=(const param<T>& p){
        _coef = add(_coef, p);
        return *this;
    }
    
    template <typename T> lterm& operator-=(const param<T>& p){
        _coef = substract(_coef, p);
        return *this;
    }
    
    template <typename T> lterm& operator*=(const param<T>& p){
        _coef = multiply(_coef, p);
        return *this;
    }

    template<typename T> lterm& operator/=(const param<T>& p){
        _coef = divide(_coef, p);
        return *this;
    }
    
};


/** A class to represent a quadratic term, e.g. 2xy or 3x^2. */
class qterm{
    
public:
    constant_*                  _coef;
    pair<param_*,param_*>*      _p;
    bool                        _sign = true; /**< True if +, flase if - */
    
    qterm(){
        _coef = nullptr;
        _p = nullptr;
    }
    
    qterm(const qterm& t){
        _coef = copy(t._coef);
        _p = new pair<param_*, param_*>(make_pair((param_*)copy(t._p->first), (param_*)copy(t._p->second)));
        _sign = t._sign;
    };
    
    qterm(qterm&& t){
        _coef = t._coef;
        t._coef = nullptr;
        _p = t._p;
        t._p = nullptr;
        _sign = t._sign;
    };
    
    
    qterm(param_* p1, param_* p2){
        _coef = new constant<int>(1);
        _p = new pair<param_*, param_*>(make_pair(p1,p2));
    };
    
    
    qterm(constant_* coef, param_* p1, param_* p2){
        _coef = coef;
        _p = new pair<param_*, param_*>(make_pair(p1,p2));
    };
    
    qterm(bool sign, param_* p1, param_* p2){
        _coef = new constant<int>(1);
        _p = new pair<param_*, param_*>(make_pair(p1,p2));
        _sign = sign;
    };
    
    qterm(bool sign, constant_* coef, param_* p1, param_* p2){
        _coef = coef;
        _p = new pair<param_*, param_*>(make_pair(p1,p2));
        _sign = sign;
    };
    
    void reverse_sign() {
        _sign = ! _sign;
    }
    
    ~qterm(){
        delete _coef;
        if (_p) {
            delete _p;
        }
    };
    
    bool operator==(const qterm& l) const;
    bool operator!=(const qterm& l) const;
    
    qterm& operator=(const qterm& l);
    qterm& operator=(qterm&& l);
    
    template <typename T> qterm& operator+=(const constant<T>& c){
        if (c.is_zero()) {
            return *this;
        }
        _coef = add(_coef, c);
        return *this;
    }
    
    template <typename T> qterm& operator-=(const constant<T>& c){
        return *this += constant<T>(c.eval()*-1);
    }
    
    
    template <typename T> qterm& operator*=(const constant<T>& c){
        if (c.is_zero()) {
            delete _coef;
            _coef = new constant<int>(0);
            _sign = true;
            cerr << "\nWARNING multiplying qterm by zero!\n";
            return *this;
        }
        if (c.is_unit()) {
            return *this;
        }
        _coef = multiply(_coef, c);
        return *this;
    }
    
    template <typename T> qterm& operator/=(const constant<T>& c){
        if (c.is_zero()) {
            throw invalid_argument("\nERROT dividing qterm by zero!\n");
            return *this;
        }
        if (c.is_unit()) {
            return *this;
        }
        _coef = divide(_coef, c);
        return *this;
    }
    
    
    template <typename T> qterm& operator+=(const param<T>& p){
        _coef = add(_coef, p);
        return *this;
    }
    
    template <typename T> qterm& operator-=(const param<T>& p){
        _coef = substract(_coef, p);
        return *this;
    }
    
    template<typename T> qterm& operator/=(const param<T>& p){
        _coef = divide(_coef, p);
        return *this;
    }

    
};



/** A linear function */
class lin: public func_{
protected:


public:
    constant_*                 _cst;/**< Constant part of the linear function */
    map<string, lterm>*        _lterms;
    
    lin(){
        func_::set_ftype(lin_);
        set_type(func_c);
        _convex = true;
        _lterms = new map<string, lterm>();
        _vars = new map<string, pair<param_*, int>>();
        _cst = new constant<int>(0);
    };

    lin(constant_* coef, param_* p);
    
    bool insert(bool sign, const constant_& coef, const param_& p);/**< Adds coef*p to the linear function. Returns true if added new term, false if only updated coef of p */
    
    lin(param_* p);
    
    lin(const lin& l){
        func_::set_ftype(lin_);
        set_type(func_c);
        _convex = l._convex;
        _lterms = new map<string, lterm>();
        _vars = new map<string, pair<param_*, int>>();
        param_* p_new;
        constant_* c_new;
        string str;
        for (auto &pair:*l._lterms) {
            p_new = (param_*)copy(pair.second._p);
            c_new = copy(pair.second._coef);
            str = p_new->get_name();
            _lterms->insert(make_pair<>(str, lterm(pair.second._sign, c_new, p_new)));
//            if (p_new->is_var()) {
                _vars->insert(make_pair<>(str, make_pair<>(p_new, 1)));
//            }
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
        _lterms->clear();
        for (auto &elem: *_vars) {
            delete elem.second.first;
        }
        _vars->clear();
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
        insert(true, constant<int>(1), p);
        return *this;
    }
    
    template<typename T> lin& operator-=(const param<T>& p){
       insert(false, constant<int>(1), p);
        return *this;
    }
    
    void const print(bool endline=true);

};



/** A quadratic function */
class quad: public lin{
protected:
    map<string, qterm>*        _qterms;

public:
    quad():lin(){
        func_::set_ftype(quad_);
        _qterms = new map<string, qterm>();
    };
    
    quad(constant_* coef, param_* p1, param_* p2);
    
    bool insert(bool sign, const constant_& coef, const param_& p);/**< Adds coef*p to the linear function. Returns true if added new term, false if only updated coef of p */
    
    bool insert(bool sign, const constant_& coef, const param_& p1, const param_& p2);/**< Adds coef*p1*p2 to the linear function. Returns true if added new term, false if only updated coef of p */
    
    quad(param_* p1, param_* p2);
    
//    quad(const lin& l):lin(l){
//        _qterms = new map<string, qterm>();
//    };
    
    quad(const quad& q):lin(q){
        func_::set_ftype(quad_);
        _qterms = new map<string, qterm>();
        param_* p_new1;
        param_* p_new2;
        string s1;
        string s2;
        constant_* c_new;
        for (auto &pair:*q._qterms) {
            auto p1 = pair.second._p->first;
            auto p2 = pair.second._p->second;
            s1 = p1->get_name();
            p_new1 = get_var(s1);
            if (!p_new1) {
                p_new1 = (param_*)copy(p1);
                add_var(p_new1);
            }
            else {
                incr_occ(s1);
            }
            s2 = p2->get_name();
            p_new2 = get_var(s2);
            if (!p_new2) {
                p_new2 = (param_*)copy(p2);
                add_var(p_new2);
            }
            else {
                incr_occ(s2);
            }
            c_new = copy(pair.second._coef);
            _qterms->insert(make_pair<>(s1+","+s2, qterm(pair.second._sign, c_new, p_new1, p_new2)));
        }
    }
    
    quad(quad&& q){
        func_::set_ftype(quad_);
        set_type(func_c);
        _convex = q._convex;
        _lterms = q._lterms;
        q._lterms = nullptr;
        _qterms = q._qterms;
        q._qterms = nullptr;
        _vars = q._vars;
        q._vars = nullptr;
        _cst = q._cst;
        q._cst = nullptr;
    }
    
    void reset(){
        func_::set_ftype(quad_);
        set_type(func_c);
        _convex = true;
        _lterms->clear();
        _qterms->clear();
        for (auto &elem: *_vars) {
            delete elem.second.first;
        }
        _vars->clear();
        delete _cst;
        _cst = new constant<int>(0);
    };
    
    ~quad();
    
    void reverse_sign(); /*<< Reverse the sign of all quadratic and linear terms + constant in the function */
    
    quad& operator=(const quad& q);
    
    quad& operator=(quad&& q);
    
    bool operator==(const quad& q) const;
    bool operator!=(const quad& q) const;
    
    bool is_constant() const{
        return (_lterms->empty() && _qterms->empty());
    }
    
    bool is_linear() const{
        return (!_lterms->empty() && _qterms->empty());
    }
    
    bool is_quadratic() const{
        return (!_qterms->empty());
    }
    
    constant_* get_cst() {
        return _cst;
    }
    
    quad& operator+=(const lin& l);
    
    quad& operator+=(const quad& q);
    quad& operator-=(const quad& q);
//    quad& operator-=(quad&& q);
    
    
    quad& operator*=(bool n){
        return *this *= constant<bool>(n);
    };
    quad& operator*=(short n){
        return *this *= constant<short>(n);
    };
    quad& operator*=(int n){
        return *this *= constant<int>(n);
    };
    quad& operator*=(float n){
        return *this *= constant<float>(n);
    };
    quad& operator*=(double n){
        return *this *= constant<double>(n);
    };
    quad& operator*=(long double n){
        return *this *= constant<long double>(n);
    };
    
    template<typename T> quad& operator+=(const constant<T>& c){
        if (c.is_zero()) {
            return *this;
        }
        _cst = add(_cst, c);
        return *this;
    }
    template<typename T> quad& operator-=(const constant<T>& c){
        if (c.is_zero()) {
            return *this;
        }
        _cst = substract(_cst, c);
        return *this;
    }
    
    template<typename T> quad& operator*=(const constant<T>& c){
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
        for (auto &pair:*_qterms) {
            pair.second *= c;
        }
        return *this;
    }
    
    template<typename T> quad& operator/=(const constant<T>& c){
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
        for (auto &pair:*_qterms) {
            pair.second /= c;
        }
        return *this;
    }
    
    template<typename T> quad& operator+=(const param<T>& p){
        insert(true, constant<int>(1), p);
        return *this;
    }
    
    template<typename T> quad& operator-=(const param<T>& p){
        string name = p.get_name();
        insert(false, constant<int>(1), p);
        return *this;
    }
    
    
    
    void const print(bool endline=true);
    
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

constant_* add(constant_* c1, const constant_& c2);

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
            if (c2.is_var()) {
                auto l = new lin(new var<T>(*(var<T>*)&c2));
                *l += *val;
                c1 = (constant_*)(l);
            }
            else{
                auto l = new lin(new param<T>(c2));
                *l += *val;
                c1 = (constant_*)(l);
            }
            return c1;
            break;
        }
        case short_c: {
            auto val = ((constant<short>*)c1);
            delete c1;
            if (c2.is_var()) {
                auto l = new lin(new var<T>(*(var<T>*)&c2));
                *l += *val;
                c1 = (constant_*)(l);
            }
            else{
                auto l = new lin(new param<T>(c2));
                *l += *val;
                c1 = (constant_*)(l);
            }
            return c1;
            break;
        }
        case integer_c: {
            auto val = ((constant<int>*)c1);
            delete c1;
            if (c2.is_var()) {
                auto l = new lin(new var<T>(*(var<T>*)&c2));
                *l += *val;
                c1 = (constant_*)(l);
            }
            else{
                auto l = new lin(new param<T>(c2));
                *l += *val;
                c1 = (constant_*)(l);
            }
            return c1;
            break;
        }
        case float_c: {
            auto val = ((constant<float>*)c1);
            delete c1;
            if (c2.is_var()) {
                auto l = new lin(new var<T>(*(var<T>*)&c2));
                *l += *val;
                c1 = (constant_*)(l);
            }
            else{
                auto l = new lin(new param<T>(c2));
                *l += *val;
                c1 = (constant_*)(l);
            }
            return c1;
            break;
        }
        case double_c: {
            auto val = ((constant<double>*)c1);
            delete c1;
            if (c2.is_var()) {
                auto l = new lin(new var<T>(*(var<T>*)&c2));
                *l += *val;
                c1 = (constant_*)(l);
            }
            else{
                auto l = new lin(new param<T>(c2));
                *l += *val;
                c1 = (constant_*)(l);
            }
            return c1;
            break;
        }
        case long_c: {
            auto val = ((constant<long double>*)c1);
            delete c1;
            if (c2.is_var()) {
                auto l = new lin(new var<T>(*(var<T>*)&c2));
                *l += *val;
                c1 = (constant_*)(l);
            }
            else{
                auto l = new lin(new param<T>(c2));
                *l += *val;
                c1 = (constant_*)(l);
            }
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

constant_* substract(constant_* c1, const constant_& c2);


template<typename T> constant_* substract(constant_* c1, const param<T>& c2){ /**< adds c2 to c1, updates its type and returns the result **/
    switch (c1->get_type()) {
        case binary_c: {
            auto val = ((constant<bool>*)c1);
            delete c1;
            if (c2.is_var()) {
                auto l = new lin(new var<T>(*(var<T>*)&c2));
                *l -= *val;
                c1 = (constant_*)(l);
            }
            else{
                auto l = new lin(new param<T>(c2));
                *l -= *val;
                c1 = (constant_*)(l);
            }
            return c1;
            break;
        }
        case short_c: {
            auto val = ((constant<short>*)c1);
            delete c1;
            if (c2.is_var()) {
                auto l = new lin(new var<T>(*(var<T>*)&c2));
                *l -= *val;
                c1 = (constant_*)(l);
            }
            else{
                auto l = new lin(new param<T>(c2));
                *l -= *val;
                c1 = (constant_*)(l);
            }
            return c1;
            break;
        }
        case integer_c: {
            auto val = ((constant<int>*)c1);
            delete c1;
            if (c2.is_var()) {
                auto l = new lin(new var<T>(*(var<T>*)&c2));
                *l -= *val;
                c1 = (constant_*)(l);
            }
            else{
                auto l = new lin(new param<T>(c2));
                *l -= *val;
                c1 = (constant_*)(l);
            }
            return c1;
            break;
        }
        case float_c: {
            auto val = ((constant<float>*)c1);
            delete c1;
            if (c2.is_var()) {
                auto l = new lin(new var<T>(*(var<T>*)&c2));
                *l -= *val;
                c1 = (constant_*)(l);
            }
            else{
                auto l = new lin(new param<T>(c2));
                *l -= *val;
                c1 = (constant_*)(l);
            }
            return c1;
            break;
        }
        case double_c: {
            auto val = ((constant<double>*)c1);
            delete c1;
            if (c2.is_var()) {
                auto l = new lin(new var<T>(*(var<T>*)&c2));
                *l -= *val;
                c1 = (constant_*)(l);
            }
            else{
                auto l = new lin(new param<T>(c2));
                *l -= *val;
                c1 = (constant_*)(l);
            }
            return c1;
            break;
        }
        case long_c: {
            auto val = ((constant<long double>*)c1);
            delete c1;
            if (c2.is_var()) {
                auto l = new lin(new var<T>(*(var<T>*)&c2));
                *l -= *val;
                c1 = (constant_*)(l);
            }
            else{
                auto l = new lin(new param<T>(c2));
                *l -= *val;
                c1 = (constant_*)(l);
            }
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
            if (c2.is_var()) {
                auto l = new lin(new var<T>(*(var<T>*)&c2));
                *l *= *val;
                c1 = (constant_*)(l);
            }
            else{
                auto l = new lin(new param<T>(c2));
                *l *= *val;
                c1 = (constant_*)(l);
            }
            return c1;
            break;
        }
        case short_c: {
            auto val = ((constant<short>*)c1);
            delete c1;
            if (c2.is_var()) {
                auto l = new lin(new var<T>(*(var<T>*)&c2));
                *l *= *val;
                c1 = (constant_*)(l);
            }
            else{
                auto l = new lin(new param<T>(c2));
                *l *= *val;
                c1 = (constant_*)(l);
            }
            return c1;
            break;
        }
        case integer_c: {
            auto val = ((constant<int>*)c1);
            delete c1;
            if (c2.is_var()) {
                auto l = new lin(new var<T>(*(var<T>*)&c2));
                *l *= *val;
                c1 = (constant_*)(l);
            }
            else{
                auto l = new lin(new param<T>(c2));
                *l *= *val;
                c1 = (constant_*)(l);
            }
            return c1;
            break;
        }
        case float_c: {
            auto val = ((constant<float>*)c1);
            delete c1;
            if (c2.is_var()) {
                auto l = new lin(new var<T>(*(var<T>*)&c2));
                *l *= *val;
                c1 = (constant_*)(l);
            }
            else{
                auto l = new lin(new param<T>(c2));
                *l *= *val;
                c1 = (constant_*)(l);
            }
            return c1;
            break;
        }
        case double_c: {
            auto val = ((constant<double>*)c1);
            delete c1;
            if (c2.is_var()) {
                auto l = new lin(new var<T>(*(var<T>*)&c2));
                *l *= *val;
                c1 = (constant_*)(l);
            }
            else{
                auto l = new lin(new param<T>(c2));
                *l *= *val;
                c1 = (constant_*)(l);
            }
            return c1;
            break;
        }
        case long_c: {
            auto val = ((constant<long double>*)c1);
            delete c1;
            if (c2.is_var()) {
                auto l = new lin(new var<T>(*(var<T>*)&c2));
                *l *= *val;
                c1 = (constant_*)(l);
            }
            else{
                auto l = new lin(new param<T>(c2));
                *l *= *val;
                c1 = (constant_*)(l);
            }
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
    if (c1->is_number() && eval(c1)==0) {
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

template<typename T1, typename T2> lin operator+(const constant<T2>& c, const param<T1>& v){
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



template<typename T1, typename T2> quad operator*(const param<T1>& v1, const param<T2>& v2){
    return ((lin() += v1) * v2);
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


template<typename T1, typename T2> lin operator/(const param<T1>& v, const constant<T2>& c);




template<typename T> quad operator+(const constant<T>& c, const quad& q){
    return quad(q) += c;
};

template<typename T> quad operator+(const constant<T>& c, quad&& q){
    return q += c;
};

template<typename T> quad operator+(const quad& q, const constant<T>& c){
    return quad(q) += c;
};

template<typename T> quad operator+(quad&& q, const constant<T>& c){
    return q += c;
};


template<typename T> quad operator+(const param<T>& v, const quad& q){
    return quad(q) += v;
};

template<typename T> quad operator+(const param<T>& v, quad&& q){
    return q += v;
};

template<typename T> quad operator+(const quad& q, const param<T>& v){
    return quad(q) += v;
};

template<typename T> quad operator+(quad&& q, const param<T>& v){
    return q +=v;
};



template<typename T> quad operator-(const constant<T>& c, const quad& q){
    return quad(c) -= q;
};

template<typename T> quad operator-(const constant<T>& c, quad&& q){
    return (q *= -1) += c;
};

template<typename T> quad operator-(const quad& q, const constant<T>& c){
    return quad(q) -= c;
};

template<typename T> quad operator-(quad&& q, const constant<T>& c){
    return q -= c;
};


template<typename T> quad operator-(const param<T>& v, const quad& q){
    return quad(v) -= q;
};

template<typename T> quad operator-(const param<T>& v, quad&& q){
    return (q *= -1) += v;
};

template<typename T> quad operator-(const quad& q, const param<T>& v){
    return quad(q) -= v;
};

template<typename T> quad operator-(quad&& q, const param<T>& v){
    return q -= v;
};



template<typename T> quad operator*(const quad& q, const constant<T>& c);
template<typename T> quad operator*(const constant<T>& c, const quad& q);
template<typename T> quad operator*(const lin& l, const param<T>& v){
    quad res;
    if (!l._cst->is_zero()) {
        res.insert(true, *l._cst, v);
    }
    param_* p_new1;
    param_* p_new2 = (param_*) &v;
    constant_* c_new;
    for (auto &pair:*l._lterms) {
        c_new = pair.second._coef;
        p_new1 = pair.second._p;
        res.insert(pair.second._sign, *c_new, *p_new1, *p_new2);
    }
    return res;
}

template<typename T> quad operator*(const param<T>& v, const lin& l){
    return l*v;
};

quad operator*(const lin& l1, const lin& l2);

template<typename T> quad operator/(const quad& q, const constant<T>& c);


template<typename T> polynome operator+(const constant<T>& c, const polynome& p);
template<typename T> polynome operator+(const polynome& p, const constant<T>& c);
template<typename T> polynome operator+(const param<T>& v, const polynome& p);
template<typename T> polynome operator+(const polynome& p, const param<T>& v);
template<typename T> polynome operator+(const constant<T>& c, const polynome& p);


template<typename T> polynome operator-(const constant<T>& c, const polynome& p);
template<typename T> polynome operator-(const polynome& p, const constant<T>& c);
template<typename T> polynome operator-(const param<T>& v, const polynome& p);
template<typename T> polynome operator-(const polynome& p, const param<T>& v);
template<typename T> polynome operator-(const constant<T>& c, const polynome& p);


template<typename T> polynome operator*(const polynome& p, const constant<T>& c);
template<typename T> polynome operator*(const constant<T>& c, const polynome& p);
template<typename T> polynome operator*(const quad& q, const param<T>& v);
template<typename T> polynome operator*(const param<T>& v, const quad& q);
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

quad operator+(quad&& q, int c);




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

