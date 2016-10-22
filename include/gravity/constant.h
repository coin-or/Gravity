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



using namespace std;



/** Backbone class for constant */
class constant_{
protected:
    CType _type;
public:

    virtual ~constant_(){};
    CType get_type() const { return _type;}
    void set_type(CType type){ _type = type;}
    
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
        return (_type==par_);
    };

    bool is_uexpr() const{
        return (_type==uexp_);
    };

    bool is_bexpr() const{
        return (_type==bexp_);
    };
    
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
            return;
        }
        if(typeid(type)==typeid(short)) {
            set_type(short_c);
            return;
        }
        if(typeid(type)==typeid(int)) {
            set_type(integer_c);
            return;
        }
        if(typeid(type)==typeid(float)) {
            set_type(float_c);
            return;
        }
        if(typeid(type)==typeid(double)) {
            set_type(double_c);
            return;
        }
        if(typeid(type)==typeid(long double)) {
            set_type(long_c);
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


constant_* copy(constant_* c2); /**< Copy c2 into a new constant_* detecting the right class, i.e., constant<>, param<>, uexpr or bexpr. */
bool equals(const constant_* c1, const constant_* c2);
double eval(int i, const constant_* c1);

uexpr cos(const constant_& c);

uexpr sin(const constant_& c);


uexpr sqrt(const constant_& c);

uexpr expo(const constant_& c);

uexpr log(const constant_& c);


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
