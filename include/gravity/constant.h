//
// Created by Hassan on 19/11/2015.
//

#ifndef GRAVITY_CONSTANT_H
#define GRAVITY_CONSTANT_H
#include <iostream>
#include <iomanip>
#include <vector>
#include <forward_list>
#include <assert.h>
#include <string>
#include <map>
#include <complex>
#include <memory>
#include <typeinfo>
#include <limits>
#include <gravity/types.h>
#include <gravity/utils.h>


using namespace std;

namespace gravity {

    /**
     Transform a scalar to a string with user-specified precision.
     @param[in] a_value number to be transformed.
     @param[in] n number of decimals in transformation.
     @return a string with the specified precision.
     */
    template<class T, class = typename enable_if<is_arithmetic<T>::value>::type>
    string to_string_with_precision(const T a_value, const int n)
    {
        std::ostringstream out;
        if(std::numeric_limits<T>::is_specialized && a_value==numeric_limits<T>::lowest()){
            return "−∞";
        }
        if(std::numeric_limits<T>::is_specialized && a_value==numeric_limits<T>::max()){
            return "+∞";
        }
        out << std::setprecision(n) << a_value;
        return out.str();
    }
    /**
     Transform a complex number to a string with user-specified precision.
     @param[in] a_value complex number to be transformed.
     @param[in] n number of decimals in transformation.
     @return a string with the specified precision.
     */
    string to_string_with_precision(const Cpx& a_value, const int n);
    
    /** Backbone class for constant */
    class constant_{
    protected:
        CType                           _type; /**< Constant type: { binary_c, short_c, integer_c, float_c, double_c, long_c, complex_c, par_c, uexp_c, bexp_c, var_c, func_c}*/
        void set_type(CType type){ _type = type;}
        
    public:
        bool                            _is_transposed = false; /**< True if the constant is transposed */
        bool                            _is_vector = false; /**< True if the constant is a vector or matrix */
        size_t                          _dim[2] = {1,1}; /*< dimension of current object */
        
        virtual ~constant_(){};
        CType get_type() const { return _type;}
        
        
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
        
        bool is_complex() const{
            return (_type==complex_c);
        };

        
        virtual bool is_number() const{
            return (_type!=par_c && _type!=uexp_c && _type!=bexp_c && _type!=var_c && _type!=func_c);
        }
        
        bool is_param() const{
            return (_type==par_c);
        };

        /**
         @return true if current object is a unitary expression.
         */
        bool is_uexpr() const{
            return (_type==uexp_c);
        };

        /**
         @return true if current object is a binary expression.
         */
        bool is_bexpr() const{
            return (_type==bexp_c);
        };
        
        /**
         @return true if current object is a unitary or binary expression.
         */
        bool is_expr() const{
            return (_type==uexp_c || _type==bexp_c);
        };

        
        bool is_var() const{
            return (_type==var_c);
        };
        
        bool is_matrix() const{
            return (_dim[0]>1 && _dim[1]>1);
        }
        
        bool is_function() const{
            return (_type==func_c);
        };
        
        
        virtual size_t get_dim(size_t i) const {
            if (i>1) {
                throw invalid_argument("In function: size_t constant_::get_dim(size_t i) const, i is out of range!\n");
            }
            return _dim[i];
        }
        
        size_t get_dim() const {
            size_t dim = _dim[0];
            if(_dim[1]>0){
                dim*= _dim[1];
            }
            return dim;
        }
        
        
        void vec(){
            _is_vector = true;
        }
        
        virtual void transpose(){
            _is_transposed = !_is_transposed;
            _is_vector = true;
            auto temp = _dim[0];
            _dim[0] = _dim[1];
            _dim[1] = temp;
        }
        
        /**
         Update the dimensions of current object after it is multiplied with c2.
         @param[in] c2 object multiplying this.
         @return true if dimensions were updated, false otherwise.
         */
        bool update_dot_dim(const constant_& c2){
            return update_dot_dim(*this, c2);
        }
        
        /**
         Update the dimensions of current object to correspond to c1.c2.
         @param[in] c1 first element in product.
         @param[in] c2 second element in product.
         @return true if dimensions were updated, false otherwise.
         */
        bool update_dot_dim(const constant_& c1, const constant_& c2){
            if(c1._is_vector || c2._is_vector){/* If both c1 and c2 are scalars, no need to update dimensions */
                if(c1.is_matrix() && (!c2._is_vector || (c2._is_vector && !c2._is_transposed))){/* If multiplying a scalar/column vector with a matrix */
                    _dim[0] = c2._dim[0];
                    _dim[1] = c2._dim[1];
                    _is_vector = true;
                    return true;
                }
                if(c2.is_matrix() && (!c1._is_vector || (c1._is_vector && !c1._is_transposed))){/* If multiplying a scalar/column vector with a matrix */
                    _dim[0] = c1._dim[0];
                    _dim[1] = c1._dim[1];
                    _is_vector = true;
                    return true;
                }
                /* Otherwise, we have a dot product */
                _dim[0] = c1._dim[0];
                _dim[1] = c2._dim[1];
                if(get_dim()>1){
                    _is_vector = true;
                }
                _is_transposed = false;/* The result of a product with a vector/matrix is not transposed */
                return true;
            }
            return false;
        }
        
        
        bool is_zero() const; /**< Returns true if constant equals 0 */
        bool is_unit() const; /**< Returns true if constant equals 1 */
        bool is_neg_unit() const; /**< Returns true if constant equals -1 */
        bool is_positive() const; /**< Returns true if constant is positive */
        bool is_negative() const; /**< Returns true if constant is negative */
        bool is_non_positive() const; /**< Returns true if constant is non positive */
        bool is_non_negative() const; /**< Returns true if constant is non negative */
    };

    template<typename type>
    class param;

    /** Polymorphic class constant, can store an arithmetic or a complex number.*/
    template<typename type = double>
    class constant: public constant_{
    protected:
        type        _val;/**< value of current constant */
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
            if(typeid(type)==typeid(complex<double>)) {
                set_type(complex_c);
                return;
            }

            throw invalid_argument("Unknown constant type.");
        }
        
        constant(const constant& c){ /**< Copy constructor */
            _type = c._type;
            _val = c._val;
            _is_transposed = c._is_transposed;
            _is_vector = c._is_vector;
            _dim[0] = c._dim[0];
            _dim[1] = c._dim[1];
        };

        constant(type val):constant(){
            _val = val;
        };
        


        ~constant(){};
        
        constant tr(){
            auto newc(*this);
            newc.transpose();
            return newc;
        };
                
        
        type eval() const { return _val;}
        
        void set_val(type val) {
            _val = val;
        }
        
        template<class T=type, class = typename enable_if<is_same<T, Cpx>::value>::type> Sign get_sign() const{
            if (_val == Cpx(0,0)) {
                return zero_;
            }
            if ((_val.real() < 0 && _val.imag() < 0)) {
                return neg_;
            }
            if ((_val.real() > 0 && _val.imag() > 0)) {
                return pos_;
            }
            if (_val.real() >= 0 && _val.imag() >= 0) {
                return non_neg_;
            }
            if (_val.real() <= 0 && _val.imag() <= 0) {
                return non_pos_;
            }
            return unknown_;
        }
        
        template<typename T=type,
        typename std::enable_if<is_arithmetic<T>::value>::type* = nullptr> Sign get_sign() const{
            if (_val==0) {
                return zero_;
            }
            if (_val > 0) {
                return pos_;
            }
            if (_val < 0) {
                return neg_;
            }
            return unknown_;
        }
        
        
        
        
        /** Operators */
        bool is_negative() const {
            return get_sign()==neg_;
        }
        
        bool is_positive() const {
            return get_sign()==pos_;
        }
        
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
        void print(int prec = 10) const{
            cout << to_str(prec);
        }
        
        void println(int prec = 10) const{
            cout << to_str(prec) << endl;
        }
        
        string to_str(int prec = 10) const{
            return to_string_with_precision(_val,prec);
        }
        
    };
    /**
     Returns the conjugate of cst.
     @param[in] cst complex number.
     @return the conjugate of cst.
     */
    constant<Cpx> conj(const constant<Cpx>& cst);
    constant<double> real(const constant<Cpx>& cst);
    constant<double> imag(const constant<Cpx>& cst);
    /**
     Returns the square magnitude of cst.
     @param[in] cst complex number.
     @return the square magnitude of cst.
     */
    constant<double> sqrmag(const constant<Cpx>& cst);
    constant<double> angle(const constant<Cpx>& cst);

}

#endif //GRAVITY_CONSTANT_H
