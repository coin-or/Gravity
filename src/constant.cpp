//
// Created by Hassan on 19/11/2015.
//

#include <gravity/constant.h>
#include <gravity/func.h>

using namespace std;
using namespace gravity;


/**
 Transform a complex number to a string with user-specified precision.
 @param[in] a_value complex number to be transformed.
 @param[in] n number of decimals in transformation.
 @return a string with the specified precision.
 */
string gravity::to_string_with_precision(const Cpx& a_value, const int n){
    std::ostringstream out;
    out << std::setprecision(n) << a_value;
    return out.str();
}


constant<Cpx> gravity::conj(const constant<Cpx>& cst){
    constant<Cpx> newc(cst);
    newc.set_val(conj(newc.eval()));
    return newc;
}

constant<double> gravity::real(const constant<Cpx>& cst){
    return real(cst.eval());
}

constant<double> gravity::sqrmag(const constant<Cpx>& cst){
    return std::pow(std::abs(cst.eval()),2);
}

constant<double> gravity::angle(const constant<Cpx>& cst){
    return arg(cst.eval());
}

constant<double> gravity::imag(const constant<Cpx>& cst){
    return imag(cst.eval());
}



