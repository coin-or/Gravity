//
//  param.cpp
//  
//
//  Created by Hassan on 20/05/2016.
//
//

#include <Gravity/param.h>


string param_::get_name() const {
    string name = _name;
    int nb = _indices->size() - 1;
    if (nb >= 0) {        
        name += "(";
        for (auto i = 0; i < nb; i++) {
            name += to_string((*_indices)[i]);
            name += ",";
        }
        name += to_string((*_indices)[nb]);
        name += ")";
    }
    return name;
};


