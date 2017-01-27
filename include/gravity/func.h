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


/** A class to represent a polynomial term, e.g. 2xyz or 3x^2y. */
class pterm{
    
public:
    constant_*                      _coef;
    vector<pair<param_*, int>>*     _p; /**< A polynomial term is represented as a vector of pairs <param_*,int> where the first element points to the parameter and the second indicates the exponent */
    bool                            _sign = true; /**< True if +, flase if - */
    
    pterm(){
        _coef = nullptr;
        _p = nullptr;
    }
    
    pterm(const pterm& t){
        _coef = copy(t._coef);
        _p = new vector<pair<param_*, int>>();
        for (auto &pair : *t._p) {
            _p->push_back(make_pair<>((param_*)copy(pair.first), pair.second));
        }
        _sign = t._sign;
    };
    
    pterm(pterm&& t){
        _coef = t._coef;
        t._coef = nullptr;
        _p = t._p;
        t._p = nullptr;
        _sign = t._sign;
    };
    
    void reverse_sign() {
        _sign = ! _sign;
    }
    
    ~pterm(){
        delete _coef;
        if (_p) {
            delete _p;
        }
    };
    
    bool operator==(const pterm& l) const;
    bool operator!=(const pterm& l) const;
    
    pterm& operator=(const pterm& l);
    pterm& operator=(pterm&& l);
    
    template <typename T> pterm& operator+=(const constant<T>& c){
        if (c.is_zero()) {
            return *this;
        }
        _coef = add(_coef, c);
        return *this;
    }
    
    template <typename T> pterm& operator-=(const constant<T>& c){
        return *this += constant<T>(c.eval()*-1);
    }
    
    
    template <typename T> pterm& operator*=(const constant<T>& c){
        if (c.is_zero()) {
            delete _coef;
            _coef = new constant<int>(0);
            _sign = true;
            cerr << "\nWARNING multiplying pterm by zero!\n";
            return *this;
        }
        if (c.is_unit()) {
            return *this;
        }
        _coef = multiply(_coef, c);
        return *this;
    }
    
    template <typename T> pterm& operator/=(const constant<T>& c){
        if (c.is_zero()) {
            throw invalid_argument("\nERROT dividing pterm by zero!\n");
            return *this;
        }
        if (c.is_unit()) {
            return *this;
        }
        _coef = divide(_coef, c);
        return *this;
    }
    
    
    template <typename T> pterm& operator+=(const param<T>& p){
        _coef = add(_coef, p);
        return *this;
    }
    
    template <typename T> pterm& operator-=(const param<T>& p){
        _coef = substract(_coef, p);
        return *this;
    }
    
    template<typename T> pterm& operator/=(const param<T>& p){
        _coef = divide(_coef, p);
        return *this;
    }
    
    
};

/** Backbone class for function */
class func_ : public constant_{
private:
    FType                                  _ftype = const_; /**< Function type, e.g., constant, linear, quadratic... >>**/
    NType                                  _return_type = integer_; /**< Return type, e.g., bool, integer, complex... >>**/
    Convexity                              _convex = linear_; /**< Convexity type, i.e., linear, convex, concave or indeterminate. >>**/
    bool                                   _in_model = false; /**< If the function is in a mathematical model, the latter is responsible for memory management. >>**/
    map<string, pair<param_*, int>>*       _params;/**< map <parameter name, <paramter pointer, number of times it appears in function>>**/
    map<string, pair<param_*, int>>*       _vars;/**< map <variable name, <variable pointer, number of times it appears in function>>**/
    constant_*                             _cst;/**< Constant part of the linear function */
    map<string, lterm>*                    _lterms; /**< Map of linear terms */
    map<string, qterm>*                    _qterms; /**< Map of quadratic terms */
    map<string, pterm>*                    _pterms; /**< Map of polynomial terms */

public:
    
    func_();
    
    func_(const constant_& c);
    
    func_(const func_& f); /**< Copy constructor */
    
    func_(func_&& f); /**< Move constructor */

    ~func_();

    bool insert(bool sign, const constant_& coef, const param_& p);/**< Adds coef*p to the function. Returns true if added new term, false if only updated coef of p */
    bool insert(bool sign, const constant_& coef, const param_& p1, const param_& p2);/**< Adds coef*p1*p2 to the function. Returns true if added new term, false if only updated coef of p1*p2 */
        
    
    void insert(const lterm& term);
    
    void insert(const qterm& term);
    
    constant_* get_cst() {
        return _cst;
    }
    
    param_* get_var(string name){
        auto pair_it = _vars->find(name);
        if (pair_it==_vars->end()) {
            return nullptr;
        }
        else {
            return get<1>(*pair_it).first;
        }
    }
    
    param_* get_param(string name){
        auto pair_it = _params->find(name);
        if (pair_it==_params->end()) {
            return nullptr;
        }
        else {
            return get<1>(*pair_it).first;
        }
    }
    
    
    void add_var(param_* v){/**< Inserts the variable in this function input list. WARNING: Assumes that v has not been added previousely!*/
        assert(_vars->count(v->get_name())==0);
        _vars->insert(make_pair<>(v->get_name(), make_pair<>(v, 1)));
    }
    
    
    void add_param(param_* p){/**< Inserts the parameter in this function input list. WARNING: Assumes that p has not been added previousely!*/
        assert(_params->count(p->get_name())==0);
        _params->insert(make_pair<>(p->get_name(), make_pair<>(p, 1)));
    }
    
    
    int nb_occ_var(string name) const{/**< Returns the number of occurences the variable has in this function. */
        auto pair_it = _vars->find(name);
        if (pair_it==_vars->end()) {
            return 0;
        }
        else {
            return get<1>(*pair_it).second;
        }
    }
    
    int nb_occ_param(string name) const{/**< Returns the number of occurences the parameter has in this function. */
        auto pair_it = _params->find(name);
        if (pair_it==_params->end()) {
            return 0;
        }
        else {
            return get<1>(*pair_it).second;
        }
    }
    
    void incr_occ_var(string str){/**< Increases the number of occurences the variable has in this function. */
        auto pair_it = _vars->find(str);
        if (pair_it==_vars->end()) {
            throw invalid_argument("Non-existing variable in function!\n");
        }
        else {
            get<1>(*pair_it).second++;
        }
    }
    
    void incr_occ_param(string str){/**< Increases the number of occurences the parameter has in this function. */
        auto pair_it = _params->find(str);
        if (pair_it==_params->end()) {
            throw invalid_argument("Non-existing variable in function!\n");
        }
        else {
            get<1>(*pair_it).second++;
        }
    }
    
    void decr_occ_var(string str){/**< Decreases the number of occurences the variable has in this function. */
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
    
    void decr_occ_param(string str){/**< Decreases the number of occurences the parameter has in this function. */
        auto pair_it = _params->find(str);
        if (pair_it==_params->end()) {
            return;
        }
        else {
            get<1>(*pair_it).second--;
            if (get<1>(*pair_it).second==0) {
                delete get<1>(*pair_it).first;
                _params->erase(pair_it);
            }
        }
    }
    
    
    pair<ind,func_*> operator[](ind i);
    
    bool is_convex() const{
        return (_convex==convex_ || _convex==linear_);
    };
    
    bool is_concave() const{
        return (_convex==concave_ || _convex==linear_);
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
    
    void reset();
    
    void reverse_sign(); /*<< Reverse the sign of all terms in the function */
    
    func_& operator=(const func_& f);
    
    func_& operator=(func_&& f);
    
    bool operator==(const func_& f) const;
    
    bool operator!=(const func_& f) const;
    
//    template<typename T> func_& operator+=(const T& v){
//        if (is_arithmetic<T>::value || typeid(T)==typeid(constant<bool>) || typeid(T)==typeid(constant<int>) || typeid(T)==typeid(constant<short>)) {
//            _cst = add(_cst, v);
//            return *this;
//        }
//        return *this;
//    };
//    template<typename T> func_& operator-=(const T& v){
//        return *this;
//    };
//    template<typename T> func_& operator*=(const T& v){
//        return *this;
//    };
//    template<typename T> func_& operator/=(const T& v){
//        return *this;
//    };
//    
//    func_& operator+=(const constant_& c){
//        return *this;
//    };
    
//    template<typename T> func_& operator+=(const constant<T>& c){
//        switch (_cst->get_type()) {
//            case binary_c: {
//                *(constant<bool>*)_cst += c;
//                return *this;
//                break;
//            }
//            case short_c: {
//                *(constant<short>*)_cst += c;
//                return *this;
//                break;
//            }
//            case integer_c: {
//                *(constant<int>*)_cst += c;
//                return *this;
//                break;
//            }
//            case float_c: {
//                *(constant<float>*)_cst += c;
//                return *this;
//                break;
//            }
//            case double_c: {
//                *(constant<double>*)_cst += c;
//                return *this;
//                break;
//            }
//            case long_c: {
//                *(constant<long double>*)_cst += c;
//                return *this;
//                break;
//            }
//            case par_c:{
//                auto f = new func_(*(param_*)_cst + c);
//                delete _cst;
//                _cst = (constant_*)f;
//                break;
//            }
//            case uexp_c: {
////                auto res = new bexpr(*(uexpr*)c1 + c2);
////                delete c1;
////                c1 = (constant_*)res;
////                return c1;
////                break;
//            }
//            case bexp_c: {
////                auto res = new bexpr(*(bexpr*)c1 + c2);
////                delete c1;
////                c1 = (constant_*)res;
////                return c1;
////                break;
//            }
//            case func_c: {
//                auto f = new func_(*(func_*)_cst + c);
//                delete _cst;
//                _cst = (constant_*)f;
//                break;
////                switch (((func_*)&c2)->get_ftype()) {
////                    case lin_:
////                        return add(c1, (*(lin*)&c2));
////                        break;
////                        
////                    default:
////                        break;
////                }
//            }
//            default:
//                break;
//        }
//        return *this;
//    }
//    template<typename T> func_& operator-=(const constant<T>& c);
//    template<typename T> func_& operator*=(const constant<T>& c);
//    template<typename T> func_& operator/=(const constant<T>& c);
//
    
    func_& operator+=(const func_& f);
    func_& operator-=(const func_& f);
    func_& operator*=(const func_& f);
    func_& operator/=(const func_& f);
    
    func_& operator+=(double v);
    func_& operator-=(double v);
    func_& operator*=(double v);
    func_& operator/=(double v);
    //    template<typename T> func_& operator+=(const param<T>& c){
//        switch (_cst->get_type()) {
//            case binary_c: {
//                auto f = new func_(*(constant<bool>*)_cst + c);
//                delete _cst;
//                _cst = (constant_*)f;
//                return *this;
//                break;
//            }
//            case short_c: {
//                auto f = new func_(*(constant<short>*)_cst + c);
//                delete _cst;
//                _cst = (constant_*)f;
//                return *this;
//                break;
//            }
//            case integer_c: {
//                auto f = new func_(*(constant<int>*)_cst + c);
//                delete _cst;
//                _cst = (constant_*)f;
//                return *this;
//                break;
//            }
//            case float_c: {
//                auto f = new func_(*(constant<float>*)_cst + c);
//                delete _cst;
//                _cst = (constant_*)f;
//
//                return *this;
//                break;
//            }
//            case double_c: {
//                auto f = new func_(*(constant<double>*)_cst + c);
//                delete _cst;
//                _cst = (constant_*)f;
//
//                return *this;
//                break;
//            }
//            case long_c: {
//                auto f = new func_(*(constant<long double>*)_cst + c);
//                delete _cst;
//                _cst = (constant_*)f;
//                return *this;
//                break;
//            }
//            case par_c:{
//                auto f = new func_(*(param_*)_cst + c);
//                delete _cst;
//                _cst = (constant_*)f;
//                break;
//            }
//            case uexp_c: {
//                //                auto res = new bexpr(*(uexpr*)c1 + c2);
//                //                delete c1;
//                //                c1 = (constant_*)res;
//                //                return c1;
//                //                break;
//            }
//            case bexp_c: {
//                //                auto res = new bexpr(*(bexpr*)c1 + c2);
//                //                delete c1;
//                //                c1 = (constant_*)res;
//                //                return c1;
//                //                break;
//            }
//            case func_c: {
//                auto f = new func_(*(func_*)_cst + c);
//                delete _cst;
//                _cst = (constant_*)f;
//                break;
//                //                switch (((func_*)&c2)->get_ftype()) {
//                //                    case lin_:
//                //                        return add(c1, (*(lin*)&c2));
//                //                        break;
//                //                        
//                //                    default:
//                //                        break;
//                //                }
//            }
//            default:
//                break;
//        }
//        return *this;
//    };
//    
//    template<typename T> func_& operator-=(const param<T>& c);
//    template<typename T> func_& operator*=(const param<T>& c);
//    template<typename T> func_& operator/=(const param<T>& c);
//    
//    template<typename T> func_& operator+=(const func_& c);
//    template<typename T> func_& operator-=(const func_& c);
//    template<typename T> func_& operator*=(const func_& c);
//    template<typename T> func_& operator/=(const func_& c);
    
    
    qterm* get_square(param_* p){ /**< Returns the quadratic term containing a square of p or nullptr if none exists. **/
        for (auto pair_it = _qterms->begin(); pair_it != _qterms->end(); pair_it++) {
            if (pair_it->second._p->first==p && pair_it->second._p->second==p) {
                return &pair_it->second;
            }
        }
        return nullptr;
    }
    
//    void update_convexity(){// WRONG!!
//        if (_qterms->empty() && _pterms->empty()) {
//            _convex = linear_;
//            return;
//        }
//        _convex = get_convexity(_qterms->begin()->second);
//        for (auto pair_it = next(_qterms->begin()); pair_it != _qterms->end(); pair_it++) {
//            if (_convex != get_convexity(pair_it->second)) {
//                _convex = indet_;
//                return;
//            }
//        }
//    }
    
    void const print(bool endline=true);
    
};


func_ operator+(const func_& f1, const func_& f2);
func_ operator+(func_&& f1, const func_& f2);

func_ operator-(const func_& f1, const func_& f2);
func_ operator-(func_&& f1, const func_& f2);


//template<typename T> func_ operator+(func_&& f, const T& v){
//    return f += v;
//};
//
//template<typename T1, typename T2> func_ operator-(const T1& v1, const T2& v2){
//    return func_();
//};
//template<typename T1, typename T2> func_ operator*(const T1& v1, const T2& v2){
//    return func_();
//};
//template<typename T1, typename T2> func_ operator/(const T1& v1, const T2& v2){
//    return func_();
//};


/** A function can return a bool, a short, an int, a float or a double or any user specified number type. The return type is set by default to int. */
//template<typename type = int>
//class func: public func_{
//protected:
//    vector<type>                _vals;
//    
//    
//public:
//    
//    func(){
////        update_type();
//    };
//    
//
//   
//    
////    void update_type() {
////        set_ftype(const_);
////        if(typeid(type)==typeid(bool)){
////            _return_type = binary_;
////            return;
////        }
////        if(typeid(type)==typeid(short)) {
////            _return_type = short_;
////            return;
////        }
////        if(typeid(type)==typeid(int)) {
////            _return_type = integer_;
////            return;
////        }
////        if(typeid(type)==typeid(float)) {
////            _return_type = float_;
////            return;
////        }
////        if(typeid(type)==typeid(double)) {
////            _return_type = double_;
////            return;
////        }
////        if(typeid(type)==typeid(long double)) {
////            _return_type = long_;
////            return;
////        }
////        throw invalid_argument("Unsupported type");
////    }
//
//    
//    /** Querries */
//    
//    
//
//};


//typedef pair<map<ind, ind>*, var_*> ind_var; /**< pair <map, var> to store indices of a given variable */






    

    
//Convexity get_convexity(const qterm& q) {
//    if(q._p->first == q._p->second){
//        if (q._sign) {
//            return convex_;
//        }
//        return concave_;
//    }
//    //q._p->first !=q._p->second
//    auto sqr1 = get_square(q._p->first);
//    auto sqr2 = get_square(q._p->second);
//    if (sqr1 && sqr2 && sqr1->_sign==sqr2->_sign) {
//        if (sqr1->_sign) {
//            return convex_;
//        }
//        return concave_;
//    }
//    return indet_;
//}


    
    
constant_* add(constant_* c1, const constant_& c2); /**< adds c2 to c1, updates its type and returns the result **/
//
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
            auto f = new func_(*c1);
            *f += c2;
            c1 =(constant_*)(f);
            return c1;
            break;
        }
        case var_c:{
            auto f = new func_(*c1);
            *f += c2;
            c1 =(constant_*)(f);
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
            auto res = new func_((*(func_*)c1) + c2);
            delete c1;
            return c1 = (constant_*)res;
            break;
        }
            
        default:
            break;
    }
    return c1;
}

//constant_* add(constant_* c1, const func_& f);

template<typename T> constant_* add(constant_* c1, const param<T>& c2){ /**< adds c2 to c1, updates its type and returns the result **/
    switch (c1->get_type()) {
        case binary_c: {
            auto val = ((constant<bool>*)c1)->eval();
            delete c1;
            auto f = new func_(c2);
            *f += val;
            c1 = (constant_*)(f);
            return c1;
            break;
        }
        case short_c: {
            auto val = ((constant<short>*)c1)->eval();
            delete c1;
            auto f = new func_(c2);
            *f += val;
            c1 = (constant_*)(f);
            return c1;
            break;
        }
        case integer_c: {
            auto val = ((constant<int>*)c1)->eval();
            delete c1;
            auto f = new func_(c2);
            *f += val;
            c1 = (constant_*)(f);
            return c1;
            break;
        }
        case float_c: {
            auto val = ((constant<float>*)c1)->eval();
            delete c1;
            auto f = new func_(c2);
            *f += val;
            c1 = (constant_*)(f);
            return c1;
            break;
        }
        case double_c: {
            auto val = ((constant<double>*)c1)->eval();
            delete c1;
            auto f = new func_(c2);
            *f += val;
            c1 = (constant_*)(f);
            return c1;
            break;
        }
        case long_c: {
            auto val = ((constant<long double>*)c1)->eval();
            delete c1;
            auto f = new func_(c2);
            *f += val;
            c1 = (constant_*)(f);
            return c1;
            break;
        }
        case par_c:{
            auto res = new func_(*c1);
            delete c1;
            res->insert(true, constant<int>(1), c2);
            return c1 = res;
            break;
        }
        case var_c:{
            auto res = new func_(*c1);
            delete c1;
            if (c2.is_var()) {
                res->insert(true, constant<int>(1), c2);
            }
            else {
                auto cst = res->get_cst();
                cst = add(cst, c2);
            }
            return c1 = res;
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
            auto res = new func_(*c1);
            delete c1;
            *res += c2;
            return c1 = (constant_*)res;
            break;
        }
        default:
            break;
    }
    return c1;
}

constant_* substract(constant_* c1, const constant_& c2);


template<typename T> constant_* substract(constant_* c1, const param<T>& c2){ /**< Substracts c2 from c1, updates its type and returns the result **/
    switch (c1->get_type()) {
        case binary_c: {
            auto val = ((constant<bool>*)c1)->eval();
            delete c1;
            auto f = new func_(c2);
            (*f *= -1) += val;
            c1 = (constant_*)(f);
            return c1;
            break;
        }
        case short_c: {
            auto val = ((constant<short>*)c1)->eval();
            delete c1;
            auto f = new func_(c2);
            (*f *= -1) += val;
            c1 = (constant_*)(f);
            return c1;
            break;
        }
        case integer_c: {
            auto val = ((constant<int>*)c1)->eval();
            delete c1;
            auto f = new func_(c2);
            (*f *= -1) += val;
            c1 = (constant_*)(f);
            return c1;
            break;
        }
        case float_c: {
            auto val = ((constant<float>*)c1)->eval();
            delete c1;
            auto f = new func_(c2);
            (*f *= -1) += val;
            c1 = (constant_*)(f);
            return c1;
            break;
        }
        case double_c: {
            auto val = ((constant<double>*)c1)->eval();
            delete c1;
            auto f = new func_(c2);
            (*f *= -1) += val;
            c1 = (constant_*)(f);
            return c1;
            break;
        }
        case long_c: {
            auto val = ((constant<long double>*)c1)->eval();
            delete c1;
            auto f = new func_(c2);
            (*f *= -1) += val;
            c1 = (constant_*)(f);
            return c1;
            break;
        }
        case par_c:{
            auto res = new func_(*c1);
            delete c1;
            res->insert(false, constant<int>(1), c2);
            return c1 = res;
            break;
        }
        case var_c:{
            auto res = new func_(*c1);
            delete c1;
            if (c2.is_var()) {
                res->insert(false, constant<int>(1), c2);
            }
            else {
                auto cst = res->get_cst();
                cst = add(cst, c2);
            }
            return c1 = res;
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
            auto res = new func_(*c1);
            delete c1;
            *res -= c2;
            return c1 = res;
        }
        default:
            break;
    }
    return c1;
}



template<typename T> constant_* multiply(constant_* c1, const param<T>& c2){ /**< adds c2 to c1, updates its type and returns the result **/
    switch (c1->get_type()) {
        case binary_c: {
            auto val = ((constant<bool>*)c1)->eval();
            delete c1;
            auto f = new func_(c2);
            *f *= val;
            c1 = (constant_*)(f);
            return c1;
            break;
        }
        case short_c: {
            auto val = ((constant<short>*)c1)->eval();
            delete c1;
            auto f = new func_(c2);
            *f *= val;
            c1 = (constant_*)(f);
            return c1;
            break;
        }
        case integer_c: {
            auto val = ((constant<int>*)c1)->eval();
            delete c1;
            auto f = new func_(c2);
            *f *= val;
            c1 = (constant_*)(f);
            return c1;
            break;
        }
        case float_c: {
            auto val = ((constant<float>*)c1)->eval();
            delete c1;
            auto f = new func_(c2);
            *f *= val;
            c1 = (constant_*)(f);
            return c1;
            break;
        }
        case double_c: {
            auto val = ((constant<double>*)c1)->eval();
            delete c1;
            auto f = new func_(c2);
            *f *= val;
            c1 = (constant_*)(f);
            return c1;
            break;
        }
        case long_c: {
            auto val = ((constant<long double>*)c1)->eval();
            delete c1;
            auto f = new func_(c2);
            *f *= val;
            c1 = (constant_*)(f);
            return c1;
            break;
        }
        case par_c:{
            auto f = new func_(*c1);
            delete c1;
            *f *= c2;
            c1 =(constant_*)(f);
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
            auto res = new func_(*c1);
            delete c1;
            *res *= c2;
            return c1 = res;
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
            auto f = new func_(*c1);
            *f -= c2;
            c1 =(constant_*)(f);
            return c1;
            break;
        }
        case var_c:{
            auto f = new func_(*c1);
            *f += c2;
            c1 =(constant_*)(f);
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
            auto res = new func_(*(func_*)c1 + c2);
            delete c1;
            return c1 = (constant_*)res;
        }        default:
            break;
    }
    return c1;
}

constant_* multiply(constant_* c1, const constant_& c2);

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
            auto l = new func_(*pc1);
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
//            auto res = new func_(*(func_*)c1 * c2);
//            delete c1;
//            return c1 = res;
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
            auto l = new func_(*pc1);
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
                    //                    auto res = new func_(*(func_*)c1 * 1/c2);
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

