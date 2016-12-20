//
// Created by Hassan on 19/11/2015.
//

#include <math.h>
#include <Gravity/constant.h>
//#include <Gravity/func.h>
#include <sstream>



/* UNARY EXPRESSIONS */

uexpr::uexpr(){
    _otype = id_;
    _son = nullptr;
    _type = uexp_c;
}


void uexpr::print(bool endline) const{
    switch (_otype) {
        case log_:
            cout << "log(";
            poly_print(_son);
            cout << ")";
            break;
            
        case exp_:
            cout << "exp(";
            poly_print(_son);
            cout << ")";
            break;
            
        case cos_:
            cout << "cos(";
            poly_print(_son);
            cout << ")";
            break;
            
        case sin_:
            cout << "sin(";
            poly_print(_son);
            cout << ")";
            break;
            
        case sqrt_:
            cout << "sqrt(";
            poly_print(_son);
            cout << ")";
            break;
        default:
            break;
    }
    if(endline)
        cout << endl;
}

void bexpr::print(bool endline) const {
    if((_otype==product_ || _otype==div_) && (_lson->get_type()==uexp_c || _lson->get_type()==bexp_c)) {
        cout << "(";
        poly_print(_lson);
        cout << ")";
    }
    else
        poly_print(_lson);

    if (_otype==plus_) {
        cout << " + ";
    }
    if (_otype==minus_) {
        cout << " - ";
    }
    if (_otype==product_) {
        cout << " * ";
    }
    if (_otype==div_) {
        cout << "/";
    }

    if (_otype==power_) {
        cout << "^";
    }
    
    if (_otype==plus_ || (_rson->get_type()!=uexp_c && _rson->get_type()!=bexp_c)) {
        poly_print(_rson);
    }
    else {
        cout << "(";
        poly_print(_rson);
        cout << ")";
    }
    if(endline)
        cout << endl;
}
