//
//  var.cpp
//  Gravity
//
//  Created by Hijazi, Hassan on 21/10/16.
//  note: Sdpvar needs to be tested (Guanglei). 
//
#include <gravity/var.h>

template<typename type> var<type>::var():param<type>(){
    param<type>::set_type(var_c);
};

template<typename type> var<type>::var(const string& name):param<type>(name){
    param<type>::set_type(var_c);
    _lb = make_shared<vector<type>>();
    _ub = make_shared<vector<type>>();
};

template<typename type> var<type>::var(const var<type>& v):param<type>(v){
    param<type>::set_type(var_c);
    _lb = v._lb;
    _ub = v._ub;
};

template<typename type> var<type>::var(var<type>&& v):param<type>(v){    
    param<type>::set_type(var_c);
    _lb = v._lb;
    _ub = v._ub;
};

template<typename type> var<type>::var(const string& name, type lb, type ub):var(name){
    _lb->push_back(lb);
    _ub->push_back(ub);
    if (lb < param<type>::_range.first) {
        param<type>::_range.first = lb;
    }
    if (ub > param<type>::_range.second) {
        param<type>::_range.second = ub;
    }
};

/* Modifiers */
template<typename type> void   var<type>::set_size(size_t s, type val){
    param<type>::set_size(s,val);
    if (_lb->empty()) {
        _lb->resize(s, numeric_limits<type>::lowest());
    }
    else {
        _lb->resize(s, _lb->at(0));
    }
    if (_ub->empty()) {
        _ub->resize(s, numeric_limits<type>::max());
    }
    else {
        _ub->resize(s, _ub->at(0));
    }
};


template<typename type> void  var<type>::add_bounds(type lb, type ub){
    if (ub<lb) {
        throw invalid_argument("add_bounds(lb, ub): Invalid bounds!");
    }
    _lb->push_back(lb);
    _ub->push_back(ub);
    if (lb < param<type>::_range.first) {
        param<type>::_range.first = lb;
    }
    if (ub > param<type>::_range.second) {
        param<type>::_range.second = ub;
    }
}

template<typename type> void   var<type>::add_lb_only(type val){
    _lb->push_back(val);
    _ub->push_back(numeric_limits<type>::max());
    if (val < param<type>::_range.first) {
        param<type>::_range.first = val;
    }
    param<type>::_range.second = numeric_limits<type>::max();
}

template<typename type> void   var<type>::add_ub_only(type val){
    _lb->push_back(numeric_limits<type>::lowest());
    _ub->push_back(val);
    param<type>::_range.first = numeric_limits<type>::lowest();
    if (val > param<type>::_range.second) {
        param<type>::_range.second = val;
    }
}

template<typename type> void   var<type>::set_lb(int i, type val){
    if (_lb->size() <= i) {
        throw out_of_range("set_lb(int i, type val)");
    }
    _lb->at(i) = val;
    if (val < param<type>::_range.first) {
        param<type>::_range.first = val;
    }
}

template<typename type> void   var<type>::set_ub(int i, type val){
    if (_ub->size() <= i) {
        throw out_of_range("set_lb(int i, type val)");
    }
    _ub->at(i) = val;
    if (val > param<type>::_range.second) {
        param<type>::_range.second = val;
    }
}

/* Operators */
template<typename type> bool var<type>::operator==(const var& v) const{
    return (param<type>::operator==(v));
//    return (param<type>::operator==(v) && *_lb==*v._lb && *_ub==*v._ub);
};

template<typename type> bool var<type>::operator!=(const var& v) const{
    return !(*this==v);
}

/* Output */
template<typename type> void var<type>::print(bool bounds) const{
    param<type>::print(false);
    if (!bounds) {
        return;
    }
    size_t idx = 0;
    cout << " in ";

    if(_lb==nullptr && _ub!=nullptr){
        if (param_::_is_indexed) {
            idx = param_::get_id_inst();
            cout << "(" << idx << ") = ";
            cout << " [ -inf , " << (*_ub)[idx] << " ]";
        }
        else {
            for(int i = 0 ; i < param_::get_dim(); i++){
                cout << "(" << idx << ") = ";
                cout << " [ -inf , " << (*_ub)[i] << " ]\n";
            }
        }
        cout << ";\n";
    }
    if(_lb!=nullptr && _ub==nullptr){
        if (param_::_is_indexed) {
            idx = param_::get_id_inst();
            cout << " [ " << (*_lb)[idx] << ", inf ]";
        }
        else {
            for(int i = 0 ; i < param_::get_dim(); i++){
                cout << "(" << i << ") = ";
                cout << " [ " << (*_lb)[i] << ", inf ]\n";
            }
        }
        cout << ";\n";
    }
    if(_lb!=nullptr && _ub!=nullptr){
        if (param_::_is_indexed) {
            idx = param_::get_id_inst();
            cout << " [ " << (*_lb)[idx] << ", " << (*_ub)[idx] << "]";
        }
        else {
            for(int i = 0 ; i < param_::get_dim(); i++){
                cout << "(" << i << ") = ";
                cout << " [ " << (*_lb)[i] << ", " << (*_ub)[i] << "]\n";
            }
        }
        cout << ";\n";
    }    
};

template<typename type> sdpvar<type>::sdpvar():var<type>(){
    var<type>::set_type(sdpvar_c);
};

template<typename type> sdpvar<type>::sdpvar(const string& name):var<type>(name){
    var<type>::set_type(sdpvar_c);
};

template<typename type> sdpvar<type>::sdpvar(const sdpvar<type>& v):var<type>(v){
    var<type>::set_type(sdpvar_c);
    _symdim = v._symdim;
};

template<typename type> sdpvar<type>::sdpvar(sdpvar<type>&& v):var<type>(v){    
    var<type>::set_type(sdpvar_c);
    _symdim = v._symdim;
};

// sdpvar 

template<typename type> sdpvar<type>::sdpvar(const string& name, type lb, type ub):sdpvar<type>(name){
    var<type>(name, lb, ub);
};

//template<typename type> void   sdpvar<type>::set_size(size_t s, type val){
//    param<type>::set_size(s,val);
//};

// Operators
template<typename type> bool sdpvar<type>::operator==(const sdpvar& v) const{
    return (param<type>::operator==(v));
};

template<typename type> bool sdpvar<type>::operator!=(const sdpvar& v) const{
    return !(*this==v);
}

/* output */
template<typename type> void sdpvar<type>::print() const{
    param<type>::print(false);
};

template class var<bool>;
template class var<short>;
template class var<int>;
template class var<float>;
template class var<double>;
template class var<long double>;

template class sdpvar<bool>;
template class sdpvar<short>;
template class sdpvar<int>;
template class sdpvar<float>;
template class sdpvar<double>;
template class sdpvar<long double>;
