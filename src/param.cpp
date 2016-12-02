//
//  param.cpp
//  
//
//  Created by Hassan on 20/05/2016.
//
//

#include <Gravity/param.h>

pair<ind,param_*> param_::operator[](ind i){
    return make_pair(i, this);
}


