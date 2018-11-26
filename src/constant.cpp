//
// Created by Hassan on 19/11/2015.
//

#include <gravity/constant.h>

using namespace std;
using namespace gravity;
constant<Cpx> conj(const constant<Cpx>& cst){
    constant<Cpx> newc(cst);
    newc.set_val(conj(newc.eval()));
    return newc;
}

constant<double> real(const constant<Cpx>& cst){
    return real(cst.eval());
}

constant<double> imag(const constant<Cpx>& cst){
    return imag(cst.eval());
}

string gravity::to_string_with_precision(const Cpx& val, const int n) {
    std::ostringstream out;
    out << std::setprecision(n) << val;
    return out.str();
}
