//
//  Line.cpp
//  Cycle_Basis_PF
//
//  Created by Sumiran on 17/06/2014.
//  Copyright (c) 2014 NICTA. All rights reserved.
//

#include "Line.h"
#include <iostream>

using namespace std;

Line::Line(){
    ch = 0;
    tr = 0;
    as = 0;
    limit = 0;
    r = 0;
    x = 0;
    g = 0;
    b = 0;
    cc = 0;
    dd = 0;
    tbound.min = 0;
    tbound.max = 0;
    status = 0;
    smax = 0;
}

Line::Line(const string& name):Line(){
    _name = name;    
}


Line::~Line(){
}

//void Line::init_complex(){
//    _Si_ = Complex("S", &pi, &qi);
//    _Sj_ = Complex("S", &pj, &qj);
//    _Si_._name.append(_name);
//    _Sj_._name.append(_name);
//}

void Line::print(){
    cout << _name << " : ";
//    cout << "(" << _src->_name<< "," <<_dest->_name << "): ";
    cout << "Expansion cost = " << cost;
    cout << "; limit = " << smax;
    cout << "; r = " << r;
    cout << "; x = " << x;
    cout << "; charge = " << ch;
    cout << "; tap ration = " << tr;
    cout << "; phase shift = " << as << std::endl;
}
