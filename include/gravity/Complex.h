//
//  Complex.h
//  PowerTools++
//
//  Created by Hassan on 30/01/2015.
//  Copyright (c) 2015 NICTA. All rights reserved.
//

#ifndef Complex_h
#define Complex_h

#include <stdio.h>
#include <gravity/constant.h>
#include <gravity/var.h>


/** Polymorphic class template, can store an arithmetic number (int. float, double..).*/
template<typename type = float> class Complex {
    
protected:
    var<type>       _real;
    var<type>       _imag;
    /** Angle */
    var<type>       _angle;
    /** Magnitude  = sqrt(_real^2 + _imag^2)*/
    var<type>       _magni;
    /** Magnitude sqaured*/
    var<type>       _smagni;
public:
    
    int             _idx;
    std::string    _name;
    bool           _polar; /**< True if we are using a polar representation of this Complex number, i.e. _magnitude(cos(angle) + i.sin(angle)). False indicates we are using a rectangular representation, i.e. _real +i._imag **/
    bool           _lifted; /**< True if we are working in the lifted space W = VV^* **/
    
    /** Constructor */
    //@{
    Complex();
    Complex(std::string name, var<>* angle);
    Complex(std::string name, var<>* real, var<>* imag);
    Complex(std::string name, var<>* real, var<>* imag, var<>* smag);
    Complex(std::string name, var<>* real, var<>* imag, var<>* angle, var<>* magni);
    Complex(std::string name, int idx, var<>* real, var<>* imag, var<>* angle, var<>* magni);
    Complex(Complex& c);
    
//    Complex(std::string name, Function* real, Function* imag, Function* angle);
//    Complex(std::string name, int idx, Function* real, Function* imag, Function* angle);
    /* Accessors */
    constant<type>* get_real();
    constant<type>* get_imag();
    constant<type>* get_angle();
    constant<type>* square_magnitude();
    //@}
    /* Modifiers */
    void lift(){_lifted = true;};
    /* Destructor */
    ~Complex();
    
    /* Operators */
    Complex& operator+=(Complex& c2);
    Complex& operator-=(Complex& c2);
    Complex& operator*=(Complex& c2);
    Complex& operator/=(Complex& c2);
    
    friend Complex operator+(Complex c1, const Complex& c2);
    friend Complex operator-(Complex c1, const Complex& c2);
    friend Complex operator*(Complex c1, const Complex& c2);
    friend Complex operator/(Complex c1, const Complex& c2);
    /** Conjugate operator */
    Complex operator--();
//    bool operator==(Complex& f);
    void print() const;
};

#endif /* defined(__PowerTools____Complex__) */
