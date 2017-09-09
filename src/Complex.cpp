//
//  Complex.cpp
//  PowerTools++
//
//  Created by Hassan on 30/01/2015.
//  Copyright (c) 2015 NICTA. All rights reserved.
//

#include <gravity/Complex.h>
#include <iostream>
#include <string>
#include <stdio.h>
#include <cstring>
using namespace std;

/** Constructor */
//@{
//Complex::Complex<>():_idx(-1), _polar(true), _lifted(false){};
//
//Complex::Complex(Complex& c): _name(c._name), _real(c._real), _imag(c._imag), _magni(c._magni), _angle(c._angle), _idx(c._idx), _polar(c._polar), _lifted(c._lifted){}
//
//Complex::Complex(string name, var<>* real, var<>* imag):_name(name), _real(*real), _imag(*imag), _polar(false), _lifted(false){
//    _idx = -1;
//    _magni += (*real^2) + (*imag^2);
//};
//
//Complex::Complex(string name, var<>* real, var<>* imag, var<>* smag):Complex(name, real, imag){
//    _smagni = Function(*smag);
//};
//
//
//Complex::Complex(string name, var<>* angle):_name(name), _angle(*angle), _polar(true), _lifted(false){
//    _idx = -1;
//};
//
//
//Complex::Complex(string name, var<>* real, var<>* imag, var<>* angle, var<>* magni):_name(name), _real(*real), _imag(*imag), _angle(*angle), _magni(*magni), _polar(true), _lifted(false){
//};

//Complex::Complex(string name, Function* real, Function* imag, Function* angle):Complex(name, 0, real, imag, angle){
//};
//
//
//Complex::Complex(string name, int idx, Function* real, Function* imag, Function* angle):_name(new string(name)), _idx(idx), _real(real), _imag(imag), _angle(angle){
////    Function* r = new Function(_real);
////    Function* i = new Function(_imag);
////    _magni = &(((*r)^2) + ((*i)^2));
//};
//
//Complex::Complex(string name, int idx, var<>* real, var<>* imag, var<>* angle, var<>* magni):_name(new string(name)), _idx(idx), _real(new Function(real)), _imag(new Function(imag)), _angle(new Function(angle)), _magni(new Function(magni)){
//    
//}

/* Accessors */
//Function* Complex::get_real(){
//    return &_real;
//};
//Function* Complex::get_imag(){
//    return &_imag;
//};
//Function* Complex::get_angle(){
//    return &_angle;
//};
//
//Function Complex::square_magnitude(){
//    if (_lifted) {
//        return (_smagni);
//    }
//    if (_polar) {
//        return (_magni^2);
//    }
//    else{
////        Function res = ((_real*_real) + (_imag*_imag));
////        res += ((_real^2) + (_imag^2));
////        res += (_real*_real);
////        res += (_imag^2);
////        res.print(true);
//        return ((_real^2) + (_imag^2));
////        return res;
//    }
//};
//
//
////@}
//
///* Destructor */
//Complex::~Complex(){    
//};


 
/* Operators */

//Complex& Complex::operator+=(Complex& c){
//    _real += c._real;
//    _imag += c._imag;
//    return *this;
//}
//
//Complex& Complex::operator-=(Complex& c){
//    _real -= c._real;
//    _imag -= c._imag;
//    return *this;
//}

//    Function rr(_real);
//    rr *= c._real;
//    Function ii(_imag);
//    ii *= c._imag;
//    _real += rr;
//    _real -= ii;
//    Function ri(_real);
//    ri *= c._imag;
//    Function ir(_imag);
//    ir *= c._real;
//    _real += ri;
//    _real += ir;
//Complex& Complex::operator*=(Complex& c){
//    _real = (_real*c._real) - (_imag*c._imag);
//    _imag = (_imag*c._real) + (_real*c._imag);
//    return *this;
//}
//
//Complex& Complex::operator/=(Complex& c){
//    _real = (_real*c._real) + (_imag*c._imag);
//    _real /= ((c._real^2) + (c._imag^2));
//    _imag = (_imag*c._real) + (_real*c._imag);
//    _imag /= ((c._real^2) + (c._imag^2));
//    return *this;
//}
//
//
//Complex operator+(Complex c1, Complex& c2){
//    return c1+=c2;
//};
//
//Complex operator-(Complex c1, Complex& c2){
//    return c1-=c2;
//};
//
///** Conjugate operator */
//Complex Complex::operator--(){
//    Complex res(*this);
//    res._imag *= -1;
//    return res;
//};
//
//void Complex::print() const{
//    cout << _name << " = ";
//    _real.print(false);
//    cout << " + i(";
//    _imag.print(false);
//    cout << ");\n";
//}
