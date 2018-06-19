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
//        if (_is_transposed) {
//            name += "^T";
//        }
        return name;
    };

    pair<Real, Real>* param_::get_range() const{
        switch (get_intype()) {
            case binary_:
                return new pair<Real,Real>(((param<bool>*)this)->_range->first, ((param<bool>*)this)->_range->second);
                break;
            case short_:
                return new pair<Real,Real>(((param<short>*)this)->_range->first, ((param<short>*)this)->_range->second);
                break;
            case integer_:
                return new pair<Real,Real>(((param<int>*)this)->_range->first, ((param<int>*)this)->_range->second);
                break;
            case float_:
                return new pair<Real,Real>(((param<float>*)this)->_range->first, ((param<float>*)this)->_range->second);
            case double_:
                return new pair<Real,Real>(((param<double>*)this)->_range->first, ((param<double>*)this)->_range->second);
                break;
            case long_:
                return new pair<Real,Real>(((param<long double>*)this)->_range->first, ((param<long double>*)this)->_range->second);
                break;
            default:
                break;
        }
        
    }
    std::vector<index_> indices(unsigned p1 ,unsigned p2){
        std::vector<index_> _keys;
        for (int i = p1; i <= p2; i++){
            _keys.push_back(index_(to_string(i)));
        }
        return _keys;
    }
}


