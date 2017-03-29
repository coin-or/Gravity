//
//  param.cpp
//  
//
//  Created by Hassan on 20/05/2016.
//
//

#include <Gravity/param.h>


string param_::get_name(bool indices) const {
    string name = _name;
    int nb = _indices->size() - 1;
    if (indices && nb >= 0) {
        name += "(";
        auto iter = _indices->begin();
        for (auto i = 0; i < nb; i++) {
            name += (*iter++).first;
            name += ",";
        }
        name += (*iter).first;
        name += ")";
    }
    return name;
};


