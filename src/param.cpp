//
//  param.cpp
//  
//
//  Created by Hassan on 20/05/2016.
//
//
#include <gravity/param.h>
namespace gravity{
    string param_::get_name(bool indices) const {
        string name = _name;
        if (_is_vector) {
            name = "[" + name + "]";
        }
        if (_is_transposed) {
            name += "^T";
        }
        return name;
    };

    pair<constant_*, constant_*>* param_::get_range() const{
        switch (get_intype()) {
            case binary_:
                return new pair<constant_*,constant_*>(new constant<bool>(((param<bool>*)this)->_range.first), new constant<bool>(((param<bool>*)this)->_range.second));
                break;
            case short_:
                return new pair<constant_*,constant_*>(new constant<short>(((param<short>*)this)->_range.first), new constant<short>(((param<short>*)this)->_range.second));
                break;
            case integer_:
                return new pair<constant_*,constant_*>(new constant<int>(((param<int>*)this)->_range.first), new constant<int>(((param<int>*)this)->_range.second));
                break;
            case float_:
                return new pair<constant_*,constant_*>(new constant<float>(((param<float>*)this)->_range.first), new constant<float>(((param<float>*)this)->_range.second));
                break;
            case double_:
                return new pair<constant_*,constant_*>(new constant<double>(((param<double>*)this)->_range.first), new constant<double>(((param<double>*)this)->_range.second));
                break;
            case long_:
                return new pair<constant_*,constant_*>(new constant<long double>(((param<long double>*)this)->_range.first), new constant<long double>(((param<long double>*)this)->_range.second));
                break;
            default:
                break;
        }
        
    }
}
