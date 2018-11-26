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

    template<class T, class = typename enable_if<is_arithmetic<T>::value>::type>
    std::string to_string_with_precision(const T a_value, const int n)
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
    
    std::string to_string_with_precision(const Cpx& a_value, const int n);
    
    /** Backbone class for constant */
    class constant_{
    protected:
        CType                           _type;
        
        
    public:
        bool                            _is_transposed = false; /**< True if the constant is transposed */
        bool                            _is_vector = false; /**< True if the constant is a vector */
        bool                            _is_conjugate = false; /**< True if the constant is a complex number and is conjugated */
        bool                            _is_sqrmag = false; /**< True if the constant is the magnitude squared of a complex number */
        bool                            _is_angle = false; /**< True if the constant is the angle of a complex number */
        bool                            _is_real = false; /**< True if the constant is the real part of a complex number */
        bool                            _is_imag = false; /**< True if the constant is the imaginary part of a complex number */
        size_t                          _dim[2] = {1,1}; /*< dimension of current object */
        
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
        
        bool is_complex() const{
            return (_type==complex_c);
        };

        
        virtual bool is_number() const{
            return (_type!=par_c && _type!=uexp_c && _type!=bexp_c && _type!=var_c && _type!=func_c);
        }
        bool is_param() const{
            return (_type==par_c);
        };

        bool is_uexpr() const{
            return (_type==uexp_c);
        };

        bool is_bexpr() const{
            return (_type==bexp_c);
        };
        
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
            if (i>0) {
                return 1;
                throw invalid_argument("In Function: size_t get_dim(size_t i) const, i is out of range!\n");
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
        
//        virtual size_t get_nb_instances() const {
//            if(_is_vector && _is_transposed && !is_matrix()){
//                return 1;
//            }
//            return _dim[0];
//        }
        
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
        
        bool update_dot(const constant_& c2){
            if(_is_transposed || c2._is_vector){
                _is_transposed = false;
                if(!is_matrix() && c2.is_matrix()){
                    _dim[0] = c2._dim[0];
                    _dim[1] = c2._dim[1];
                }
                else if(!is_matrix() || c2.is_matrix() || !c2._is_transposed){
                    _dim[1] = c2._dim[1];
                }
                if(is_matrix()){
                    _is_vector = true;
                }
                return true;
            }
            return false;
        }
        
        bool update_dot_dim(const constant_& c1, const constant_& c2){
            _is_transposed = false;
            if(c1._is_vector || c2._is_vector){
                _dim[0] = c1._dim[0];
                _dim[1] = c2._dim[1];
                /* Instructions above work for matrix and vector dot products, below we check if it's a component-wise vector,matrix product */
                if(!c1.is_matrix() && c2.is_matrix()){
                    _dim[0] = c2._dim[0];
                }
                if(c1.is_matrix() && !c2.is_matrix() && c2._is_transposed){
                    _dim[1] = c1._dim[1];
                }
                if(is_matrix()){
                    _is_vector = true;
                }
                return true;
            }
            return false;
        }
        /* Update the dimensions based on the product (this*c2) */
//        bool update_dot_dim(const constant_& c2) {
//            if (_is_vector || c2._is_vector) {
//                if(!is_matrix() || c2.is_matrix()){
//                    _dim[1] = c2._dim[1];
//                }
//                if (is_matrix()) {
//                    _is_transposed = false;
//                }
//                else {//this is a vector
//                    if (_is_transposed) {
//                        if (!c2._is_transposed && !c2.is_matrix()) {//if c2 is not transposed and is a vector, we have a scalar product
//                            _is_vector = false;
//                            _is_transposed = false;
//                        }
//                        else {//c2 is either transposed or a matrix at this stage
//                            _is_transposed = true;
//                            _is_vector = true;
//                        }
//                    }
//                    else {//this is a column vector
//                        if (!c2.is_matrix()) {//if c2 is not a matrix, the result is component wise vector product
//                            _is_vector = true;
//                            _is_transposed = false;
//                        }
//                        else {//c2 a matrix, the result is a vector matrix columnwise product
//                            _is_transposed = c2._is_transposed;
//                        }
//                    }
//
//                }
//                if(!is_matrix() && c2.is_matrix()){
//                    _dim[0] = c2._dim[0];
//                }
//                return true;
//            }
//            return false;
//        }
        
        
        Sign get_all_sign() const;
        Sign get_sign(size_t idx=0) const;
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

    /** Polymorphic class constant, can store an arithmetic number (int. float, double..).*/
    template<typename type = double>
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
            _is_angle = c._is_angle;
            _is_sqrmag = c._is_sqrmag;
            _is_conjugate = c._is_conjugate;
            _is_real = c._is_real;
            _is_imag = c._is_imag;
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
        
        Sign get_sign() const{
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
            return _val < 0;
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
    
    constant<Cpx> conj(const constant<Cpx>& cst);
    constant<double> real(const constant<Cpx>& cst);
    constant<double> imag(const constant<Cpx>& cst);

}

#endif //GRAVITY_CONSTANT_H
