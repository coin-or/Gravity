//
// Created by Hassan on 19/11/2015.
//

#ifndef GRAVITY_CONSTANT_H
#define GRAVITY_CONSTANT_H
#include <iostream>
#include <vector>
#include <forward_list>
#include <assert.h>
#include <string>
#include <memory>
#include <Gravity/Type.h>


class param_;
class expr;
class uexpr;
class bexpr;
using namespace std;



/** Backbone class for constant */
class constant_{
protected:
    ConstType _type;
public:

    virtual ~constant_(){};
    ConstType get_type() const { return _type;}
    void set_type(ConstType type){ _type = type;}
    
    /** Querries */
    bool is_binary() const{
        return (_type==binary_c);
    };
    
    bool is_short() const{
        return (_type==short_c);
    }
    
    bool is_integer() const{
        return (_type==integer_c);
    };
    
    bool is_float() const{
        return (_type==float_c);
    };
    
    bool is_double() const{
        return (_type==double_c);
    };
    
    bool is_long() const{
        return (_type==long_c);
    };
    
    bool is_param() const{
        return (_type==parameter);
    };

    bool is_uexpr() const{
        return (_type==unary_exp);
    };

    bool is_bexpr() const{
        return (_type==binary_exp);
    };
    
    

    
//        friend constant sin(constant& c);
//        friend constant sqrt(constant& c);
//        friend constant expo(constant& c);
//        friend constant log(constant& c);
//    
//        friend constant cos(constant&& c);
//        friend constant sin(constant&& c);
//        friend constant sqrt(constant&& c);
//        friend constant expo(constant&& c);
//        friend constant log(constant&& c);

};

/** Polymorphic class constant, can store an arithmetic number (int. float, double..).*/
template<typename type = short>
class constant: public constant_{
protected:
    type        _val;
public:
    
    /** Constructors */
    constant(){
        if(typeid(type)==typeid(bool)){
            set_type(binary_c);
//            _val = false;
            return;
        }
        if(typeid(type)==typeid(short)) {
            set_type(short_c);
//            _val = 0;
            return;
        }
        if(typeid(type)==typeid(int)) {
            set_type(integer_c);
//            _val = 0;
            return;
        }
        if(typeid(type)==typeid(float)) {
            set_type(float_c);
//            _val = 0;
            return;
        }
        if(typeid(type)==typeid(double)) {
            set_type(double_c);
//            _val = 0;
            return;
        }
        if(typeid(type)==typeid(long double)) {
            set_type(long_c);
//            _val = 0;
            return;
        }
        throw invalid_argument("Unknown constant type.");
    }
    
    constant(const constant& c){ /**< Copy constructor */
        _type = c._type;
        _val = c._val;
    };

    constant(type val):constant(){
        _val = val;
    };
    


    ~constant(){};
    
    type eval() const { return _val;}
    
    void set_val(type val) {
        _val = val;
    }
    
    /** Operators */
    bool operator==(const constant& c) const {
        return (_type==c._type && _val==c._val);
    }
    
    bool operator==(const type& v) const{
        return _val==v;
    }

    constant& operator=(const type& val){
        _val = val;
        return *this;
    }
    
    constant& operator+=(const type& v){
        _val += v;
        return *this;
    }
    
    constant& operator-=(const type& v){
        _val -= v;
        return *this;
    }
    
    constant& operator*=(const type& v){
        _val *= v;
        return *this;
    }
    
    constant& operator/=(const type& v){
        _val /= v;
        return *this;
    }
    
    friend constant operator+(const constant& c1, const constant& c2){
        return constant(c1._val + c2._val);
    }
    
    friend constant operator-(const constant& c1, const constant& c2){
        return constant(c1._val - c2._val);
    }
    
    friend constant operator/(const constant& c1, const constant& c2){
        return constant(c1._val / c2._val);
    }
    
    friend constant operator*(const constant& c1, const constant& c2){
        return constant(c1._val * c2._val);
    }
    
    friend constant operator^(const constant& c1, const constant& c2){
        return constant(pow(c1._val,c2._val));
    }


    friend constant operator+(const constant& c, type cst){
        return constant(c._val + cst);
    }
    
    friend constant operator-(const constant& c, type cst){
        return constant(c._val - cst);
    }
    
    friend constant operator*(const constant& c, type cst){
        return constant(c._val * cst);
    }

    
    friend constant operator/(const constant& c, type cst){
        return constant(c._val / cst);
    }

    friend constant operator+(type cst, const constant& c){
        return constant(c._val + cst);
    }
    
    friend constant operator-(type cst, const constant& c){
        return constant(cst - c._val);
    }
    
    friend constant operator*(type cst, const constant& c){
        return constant(c._val * cst);
    }
    
    
    friend constant operator/(type cst, const constant& c){
        return constant(cst / c._val);
    }

    friend constant cos(const constant& c){
        return constant(cos(c._val));
    }
    
    friend constant sin(const constant& c){
        return constant(sin(c._val));
    }
    
    friend constant sqrt(const constant& c){
        return constant(sqrt(c._val));
    }
    
    friend constant expo(const constant& c){
        return constant(exp(c._val));
    }
    
    friend constant log(const constant& c){
        return constant(log(c._val));
    }

    
    /** Output */
    void print() const{
        cout << _val;
    }

    
};

/** Polymorphic class uexpr (unary expression), stores a unary expression tree. */

class expr{

public:
virtual ~expr(){};

};

class uexpr: public constant_, public expr{

public:
    OperatorType    _otype;
    constant_*      _son;
    
    uexpr();
    uexpr(const uexpr& exp);
    uexpr(uexpr&& exp);
    
    uexpr& operator=(const uexpr& e);
    uexpr& operator=(uexpr&& e);
    
    
    ~uexpr(){
            delete _son;        
    };
    
    bool contains(const constant_* c) const;
    
    
    void reset(){
        delete _son;
        _son = nullptr;

    };
    
    
    OperatorType get_otype() const{
        return _otype;
    };
    
    /** Operators */
    
    
    bool operator==(const uexpr &c)const;
    
    bool operator!=(const uexpr& c) const{
        return !(*this==c);
    };

    double eval(unsigned int i) const;
    
    double eval() const{
        return eval(0);
    }
    
    void print(bool endline = true) const;

    //
//    expr& operator-=(const type& v){
//        *_arg -= v;
//        return *this;
//    }
//    
//    expr& operator*=(const type& v){
//        *_arg *= v;
//        return *this;
//    }
//    
//    expr& operator/=(const type& v){
//        *_arg /= v;
//        return *this;
//    }

};


class bexpr: public constant_, public expr{
private:
    
public:
    OperatorType    _otype;
    constant_*      _lson;
    constant_*      _rson;
    
    bexpr();
    
    bexpr(const bexpr& exp);
    bexpr(bexpr&& exp);
    
    bexpr& operator=(const bexpr& e);
    
    bexpr& operator=(bexpr&& e);
    
    ~bexpr(){
        delete _lson;
        delete _rson;
    }
    
    void reset(){
        _otype = id_;
        delete _lson;
        _lson = nullptr;
        delete _rson;
        _rson = nullptr;
    };

    bool is_lson(const constant_* c) const{
        return (_lson == c);
    };

    bool is_rson(const constant_* c) const{
        return (_rson == c);
    };

    constant_* get_lson() const{
        return _lson;
    };

    constant_* get_rson() const{
        return _rson;
    };
        
    void set_lson(constant_* c){
        delete _lson;
        _lson = c;
    };
        
    void set_rson(constant_* c){
        delete _rson;
        _rson = c;
    };
        
    OperatorType get_otype() const {
        return _otype;
    };

    bool contains(const constant_* c) const;
    
    bool operator==(const bexpr &c)const;
    
    bool operator!=(const bexpr& c) const{
        return !(*this==c);
    };
    
    template<typename other_type> bexpr& operator+=(const other_type& v);
    template<typename other_type> bexpr& operator-=(const other_type& v);
    template<typename other_type> bexpr& operator*=(const other_type& v);
    template<typename other_type> bexpr& operator/=(const other_type& v);
    
    
    void print(bool endline = true) const;
    
    void print_tree() const;

    double eval(unsigned int i) const;
    
};





    //            switch (_arg->get_type()) {
    //                case binary_p: {
    //                    auto c = static_cast<param<bool>*>(_arg);
    //                    *c += v;
    //                    break;
    //                }
    //                case short_p: {
    //                    auto c = static_cast<param<short>*>(_arg);
    //                    *c += v;
    //                    break;
    //                }
    //                case integer_p: {
    //                    auto c = static_cast<param<int>*>(_arg);
    //                    *c += v;
    //                    break;
    //                }
    //                case float_p: {
    //                    auto c = static_cast<param<float>*>(_arg);
    //                    *c += v;
    //                    break;
    //                }
    //                case double_p: {
    //                    auto c = static_cast<param<double>*>(_arg);
    //                    *c += v;
    //                    break;
    //                }
    //                case long_p: {
    //                    auto c = static_cast<param<long double>*>(_arg);
    //                    *c += v;
    //                    break;
    //                }
    //                default:
    //                    throw bad_typeid();
    //                    break;
    //            }

    
//    friend constant operator+(const constant& c1, const constant& c2){
//        return constant(c1._val + c2._val);
//    }
//    
//    friend constant operator-(const constant& c1, const constant& c2){
//        return constant(c1._val - c2._val);
//    }
//    
//    friend constant operator/(const constant& c1, const constant& c2){
//        return constant(c1._val / c2._val);
//    }
//    
//    friend constant operator*(const constant& c1, const constant& c2){
//        return constant(c1._val * c2._val);
//    }
//    
//    friend constant operator^(const constant& c1, const constant& c2){
//        return constant(pow(c1._val,c2._val));
//    }
//    
//    
//    friend constant operator+(const constant& c, type cst){
//        return constant(c._val + cst);
//    }
//    
//    friend constant operator-(const constant& c, type cst){
//        return constant(c._val - cst);
//    }
//    
//    friend constant operator*(const constant& c, type cst){
//        return constant(c._val * cst);
//    }
//    
//    
//    friend constant operator/(const constant& c, type cst){
//        return constant(c._val / cst);
//    }
//    
//    friend constant operator+(type cst, const constant& c){
//        return constant(c._val + cst);
//    }
//    
//    friend constant operator-(type cst, const constant& c){
//        return constant(cst - c._val);
//    }
//    
//    friend constant operator*(type cst, const constant& c){
//        return constant(c._val * cst);
//    }
//    
//    
//    friend constant operator/(type cst, const constant& c){
//        return constant(cst / c._val);
//    }
//    
//    friend constant cos(const constant& c){
//        return constant(cos(c._val));
//    }
//    
//    friend constant sin(const constant& c){
//        return constant(sin(c._val));
//    }
//    
//    friend constant sqrt(const constant& c){
//        return constant(sqrt(c._val));
//    }
//    
//    friend constant expo(const constant& c){
//        return constant(exp(c._val));
//    }
//    
//    friend constant log(const constant& c){
//        return constant(log(c._val));
//    }
//    
//    
//    /** Output */
//    void print() const{
//        cout << _val;
//    }
//
//
//};
//
//
//
//
//public:
//    /** Constructors */
//
//
//
//
//
//
//    bool power_factor(const constant& c) const;
//    bool factor_of(const constant& c) const;
//    bool inverse_factor_of(const constant& c) const;
//
//    /** Accessors */
//
//    
//
//
//
//
//
//
//    /** Modifiers */
//    void shallow_copy(const constant& c);// Ignore _lson and _rson in Copy
//
//
//    /* Operators */
//
//
//
//
//
//
//    friend constant operator+(constant&& c, double cst);
//    friend constant operator-(constant&& c, double cst);
//    friend constant operator*(constant&& c, double cst);
//    friend constant operator/(constant&& c, double cst);
//
//
//    friend constant operator^(const constant& c, int p);
//
//    friend constant operator+(double cst, const constant& c);
//    friend constant operator-(double cst, const constant& c);
//    friend constant operator*(double cst, const constant& c);
//    friend constant operator/(double cst, const constant& c);
//
//    friend constant operator+(double cst, constant&& c);
//    friend constant operator-(double cst, constant&& c);
//    friend constant operator*(double cst, constant&& c);
//
//    friend constant cos(constant& c);
//    friend constant sin(constant& c);
//    friend constant sqrt(constant& c);
//    friend constant expo(constant& c);
//    friend constant log(constant& c);
//
//    friend constant cos(constant&& c);
//    friend constant sin(constant&& c);
//    friend constant sqrt(constant&& c);
//    friend constant expo(constant&& c);
//    friend constant log(constant&& c);
//
//
//
//    /** Output */

//    void print_() const;
//    void print_tree() const;
//    string to_string() const;
//
//
//};


//class meta_constant: public constant{
//protected:
//    shared_ptr<forward_list<pair<Function*, vector<constant*>>>>        _conc;/**< Pointer to a list of vectors pointing to concretizations of this constant in every function it appears.*/
//
//public:
////        unsigned long get_nbc() const { /**< A meta_constant can have multiple concretisations, return the number of concretisations */
////            return _conc->size();
////        }
//
//        unsigned long assign(double v); /**< assign new concretisation, return the number of concretisations */
//        constant* get_concretization(Function *f, int i) const;
//        void concretize(Function *f, constant* c);
//        double eval(Function* f, int i) const;
//        double eval(Function* f) const;
//
//};

//uexpr(const uexp& exp){ /**< Constructor from binary expression tree */
//    _arg = new bexp();
//    _arg->_otype = exp._otype;
//    switch (_arg->_lson->_type) {
//        case binary_c: {
//            _arg->_lson = new constant<bool>(((constant<bool>*)(_arg->_lson))->get_val());
//            break;
//        }
//        case short_c: {
//            _arg->_lson = new constant<short>(((constant<short>*)(_arg->_lson))->get_val());
//            break;
//        }
//        case integer_c: {
//            _arg->_lson = new constant<int>(((constant<int>*)(_arg->_lson))->get_val());
//            break;
//        }
//        case float_c: {
//            _arg->_lson = new constant<float>(((constant<float>*)(_arg->_lson))->get_val());
//            break;
//        }
//        case double_c: {
//            _arg->_lson = new constant<double>(((constant<double>*)(_arg->_lson))->get_val());
//            break;
//        }
//        case long_c: {
//            _arg->_lson = new constant<long double>(((constant<long double>*)(_arg->_lson))->get_val());
//            break;
//        }
//        default:
//            break;
//    }
//    switch (_arg->_rson->_type) {
//        case binary_c: {
//            _arg->_rson = new constant<bool>(((constant<bool>*)(_arg->_rson))->get_val());
//            break;
//        }
//        case short_c: {
//            _arg->_rson = new constant<short>(((constant<short>*)(_arg->_rson))->get_val());
//            break;
//        }
//        case integer_c: {
//            _arg->_rson = new constant<int>(((constant<int>*)(_arg->_rson))->get_val());
//            break;
//        }
//        case float_c: {
//            _arg->_rson = new constant<float>(((constant<float>*)(_arg->_rson))->get_val());
//            break;
//        }
//        case double_c: {
//            _arg->_rson = new constant<double>(((constant<double>*)(_arg->_rson))->get_val());
//            break;
//        }
//        case long_c: {
//            _arg->_rson = new constant<long double>(((constant<long double>*)(_arg->_rson))->get_val());
//            break;
//        }
//        default:
//            break;
//    }
//    
//    _type = binary_exp;
//};

constant_* copy(constant_* c2); /**< Copy c2 into a new constant_* detecting the right class, i.e., constant<>, param<>, uexpr or bexpr. */
bool equals(const constant_* c1, const constant_* c2);
double eval(int i, const constant_* c1);

uexpr cos(const constant_& c);

uexpr sin(const constant_& c);


uexpr sqrt(const constant_& c);

uexpr expo(const constant_& c);

uexpr log(const constant_& c);

//bexpr operator+(const constant_& c1, const constant_& c2);
//bexpr operator-(const constant_& c1, const constant_& c2);
//bexpr operator*(const constant_& c1, const constant_& c2);
//bexpr operator/(const constant_& c1, const constant_& c2);


template<typename T1, typename T2> bexpr operator+(const T1& c1, const T2& c2){
    bexpr res;
    res._otype = plus_;
    if (is_arithmetic<T1>::value) {
        res._lson = new constant<T1>(c1);
    }
    else{
        res._lson = copy((constant_*)&c1);
    }
    if (is_arithmetic<T2>::value) {
        res._rson = new constant<T2>(c2);
    }
    else{
        res._rson =  copy((constant_*)&c2);
    }
    return res;
}

template<typename T1, typename T2> bexpr operator-(const T1& c1, const T2& c2){
    bexpr res;
    res._otype = minus_;
    if (is_arithmetic<T1>::value) {
        res._lson = new constant<T1>(c1);
    }
    else{
        res._lson = copy((constant_*)&c1);
    }
    if (is_arithmetic<T2>::value) {
        res._rson = new constant<T2>(c2);
    }
    else{
        res._rson =  copy((constant_*)&c2);
    }
    return res;
}

template<typename T1, typename T2> bexpr operator*(const T1& c1, const T2& c2){
    bexpr res;
    res._otype = product_;
    if (is_arithmetic<T1>::value) {
        res._lson = new constant<T1>(c1);
    }
    else{
        res._lson = copy((constant_*)&c1);
    }
    if (is_arithmetic<T2>::value) {
        res._rson = new constant<T2>(c2);
    }
    else{
        res._rson =  copy((constant_*)&c2);
    }
    return res;
}

template<typename T1, typename T2> bexpr operator/(const T1& c1, const T2& c2){
    bexpr res;
    res._otype = div_;
    if (is_arithmetic<T1>::value) {
        res._lson = new constant<T1>(c1);
    }
    else{
        res._lson = copy((constant_*)&c1);
    }
    if (is_arithmetic<T2>::value) {
        res._rson = new constant<T2>(c2);
    }
    else{
        res._rson =  copy((constant_*)&c2);
    }
    return res;
}

template<typename T1, typename T2> bexpr operator^(const T1& c1, const T2& c2){
    bexpr res;
    res._otype = power_;
    if (is_arithmetic<T1>::value) {
        res._lson = new constant<T1>(c1);
    }
    else{
        res._lson = copy((constant_*)&c1);
    }
    if (is_arithmetic<T2>::value) {
        res._rson = new constant<T2>(c2);
    }
    else{
        res._rson =  copy((constant_*)&c2);
    }
    return res;
}


void poly_print(const constant_* c);

#endif //GRAVITY_CONSTANT_H
