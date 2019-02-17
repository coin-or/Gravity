//
//  var.cpp
//  Gravity
//
//  Created by Hijazi, Hassan on 21/10/16.
//
#include <gravity/var.h>
#include <gravity/func.h>
#include <gravity/Net.h>

namespace gravity {

    
    
    
//template<typename type> var<type>::var(const string& name, Sign s):param<type>(name) {
//    constant_::set_type(var_c);
//    _lb = make_shared<func<type>>();
//    _ub = make_shared<func<type>>();
//    if (s==non_neg_ || s==pos_) {
//        if(numeric_limits<type>::is_specialized){
//            add_lb_only(0);
//        }
//    }
//    else if (s==non_pos_ || s==neg_) {
//        if(numeric_limits<type>::is_specialized){
//            add_ub_only(0);
//        }
//    }
//};

template<typename type> var<type>::var(const var<type>& v){
    *this = v;
};

template<typename type> var<type>::var(var<type>&& v){
    *this = move(v);
};

    


template<typename type> var<type>& var<type>::operator=(const var<type>& v) {
    this->param<type>::operator=(v);
    constant_::set_type(var_c);
    _lb = v._lb;
    _ub = v._ub;
    return *this;
};

template<typename type> var<type>& var<type>::operator=(var<type>&& v) {
    this->param<type>::operator=(move(v));
    constant_::set_type(var_c);
    _lb = move(v._lb);
    _ub = move(v._ub);
    return *this;
};

/* Modifiers */
    
    /* Create a vector of variables indexed as pair of nodes from bags of size bag_size (WARNING assumes bags are unique in bags)*/
    template<typename type>
    vector<var<type>> var<type>::pairs_in_bags(const vector<vector<Node*>>& bags, size_t bag_size){
    vector<var> res;
    vector<indices> ids_vec;
    string key;
    res.resize(bag_size);
    ids_vec.resize(bag_size);
//    set<vector<Node*>> unique_bags;
    for (auto i = 0; i<bag_size; i++) {
        ids_vec[i] = *param_::_indices;
        ids_vec[i]._name = "pairs_"+to_string(i);
        if(ids_vec[i]._ids){
            ids_vec[i]._ids = nullptr;
        }
    }
    for (auto &bag: bags) {
        /* Make sure it's a bag with size=bag_size */
//        if (bag.size() == bag_size && unique_bags.insert(bag).second) {
        if (bag.size() == bag_size) {
            for (size_t i = 0; i< bag_size-1; i++) {
                key = bag[i]->_name + "," + bag[i+1]->_name;
                ids_vec[i].add(key);
            }
            /* Loop back pair */
            key = bag[0]->_name + "," + bag[bag_size-1]->_name;
            ids_vec[bag_size-1].add(key);
        }
    }
    for (auto i = 0; i<bag_size; i++) {
        res[i] = this->in(ids_vec[i]);
    }
    return res;
}
    
/* Create a vector of variables indexed based on nodes from bags of size bag_size */
template<typename type>
vector<var<type>> var<type>::in_bags(const vector<vector<Node*>>& bags, size_t bag_size){
    vector<var> res;
    vector<indices> ids_vec;
    res.resize(bag_size);
    ids_vec.resize(bag_size);
    set<vector<Node*>> unique_bags;
    for (auto i = 0; i<bag_size; i++) {
        ids_vec[i] = *param_::_indices;
        ids_vec[i]._name = "nodes_"+to_string(i);
    }
    for (auto &bag: bags) {
        /* Make sure it's a new bag with size=bag_size */
        if (bag.size() == bag_size && unique_bags.insert(bag).second) {
            for (size_t i = 0; i< bag_size; i++) {
                ids_vec[i].add(bag[i]->_name);
            }
        }
    }
    for (auto i = 0; i<bag_size; i++) {
        res[i] = this->in(ids_vec[i]);
    }
    return res;
}
    
template<typename type> void   var<type>::set_size(vector<size_t> dims) {
    param<type>::set_size(dims);
};

template<typename type> void   var<type>::set_size(size_t s) {
    param<type>::set_size(s);
};

template<typename type>
type    var<type>::get_lb(size_t i) const {
    if (_lb->is_number()) {
        return _lb->_val->at(0);
    }
    return _lb->eval(i);
};

template<typename type>
type    var<type>::get_ub(size_t i) const {
    if (_ub->is_number()) {
        return _ub->_val->at(0);
    }
    return _ub->eval(i);
};
    
    template <typename type>
    template<typename T,
    typename std::enable_if<is_arithmetic<T>::value>::type*>
    bool var<type>::is_bounded_above(size_t i) const {
        return (_ub->eval(i)!=numeric_limits<type>::max());
    };
    
    template <typename type>
    template<typename T,
    typename std::enable_if<is_arithmetic<T>::value>::type*>
    bool var<type>::is_bounded_below(size_t i) const {
        return (_lb->eval(i)!=numeric_limits<type>::lowest());
    };
    
    template <typename type>
    template<typename T,
    typename std::enable_if<is_arithmetic<T>::value>::type*>
    bool var<type>::is_constant(size_t i) const {
        return (is_bounded_below() && is_bounded_above() && _lb->eval(i)==_ub->eval(i));
    };
//
    

    template <typename type>
    template<typename T,
    typename std::enable_if<is_arithmetic<T>::value>::type*>
    Sign var<type>::get_sign(size_t idx) const{
        if (_lb->eval(idx) == 0 && _ub->eval(idx) == 0) {
            return zero_;
        }
        if (_ub->eval(idx) < 0) {
            return neg_;
        }
        if (_lb->eval(idx) > 0) {
            return pos_;
        }
        if (_ub->eval(idx) == 0) {
            return non_pos_;
        }
        if (_lb->eval(idx) == 0) {
            return non_neg_;
        }
        return unknown_;
    }
    
//    template Sign var<bool>::get_sign(size_t t) const;
//    template Sign var<short>::get_sign(size_t t) const;
//    template Sign var<int>::get_sign(size_t t) const;
//    template Sign var<float>::get_sign(size_t t) const;
//    template Sign var<double>::get_sign(size_t t) const;
//    template Sign var<long double>::get_sign(size_t t) const;
    
//template<typename type> void  var<type>::add_bounds(type lb, type ub) {
//    if (ub<lb) {
//        throw invalid_argument("add_bounds(lb, ub): Invalid bounds!");
//    }
//    _lb->set_val(lb);
//    _ub->push_back(ub);
//    if (lb < param<type>::_range->first) {
//        param<type>::_range->first = lb;
//    }
//    if (ub > param<type>::_range->second) {
//        param<type>::_range->second = ub;
//    }
//}

template<typename type> void   var<type>::add_lb_only(type val) {
    *_lb = constant<type>(val);
    param<type>::_range->first = val;
    if(numeric_limits<type>::is_specialized){
        *_ub = constant<type>(numeric_limits<type>::max());
        param<type>::_range->second = numeric_limits<type>::max();
    }
    
    
}

template<typename type> void   var<type>::add_ub_only(type val) {
    if(numeric_limits<type>::is_specialized){
        *_lb = constant<type>(numeric_limits<type>::lowest());
        param<type>::_range->first = numeric_limits<type>::lowest();
    }
    *_ub = constant<type>(val);
    param<type>::_range->second = val;
}

//template<typename type> void   var<type>::set_lb(int i, type val) {
//    *_lb = constant<type>(val);
//    param<type>::_range->first = val;
//}
//
//template<typename type> void   var<type>::set_ub(int i, type val) {
//    *_ub = constant<type>(val);
//    param<type>::_range->second = val;
//}

    
template<typename type> void   var<type>::set_lb(type val) {
    param<type>::_range->first = val;
}

template<typename type> void   var<type>::set_ub(type val) {
    param<type>::_range->second = val;
}

template<typename type> void   var<type>::in_q_cone() {
        _in_q_cone = true;
    };

//    template<typename type> void  var<type>::set_lb(string name, type val){
//        auto it = param_::_indices->find(name);
//        if (it !=param_::_indices->end()) {
//            _ub->at(it->second) = val;
//            if (val > param<type>::_range->first) {
//                param<type>::_range->second = val;
//            }
//        }
//        else{
//            throw out_of_range("set_lb(int i, type val)");
//        }
//    }
//
//    template<typename type> void var<type>::set_ub(string name, type val){
//        auto it = param_::_indices->find(name);
//        if (it !=param_::_indices->end()) {
//            _ub->at(it->second) = val;
//            if (val > param<type>::_range->first) {
//                param<type>::_range->second = val;
//            }
//        }
//        else{
//            throw out_of_range("set_lb(int i, type val)");
//        }
//    }

/* Operators */
template<typename type> bool var<type>::operator==(const var& v) const {
    return (param<type>::operator==(v));
    //    return (param<type>::operator==(v) && *_lb==*v._lb && *_ub==*v._ub);
};

template<typename type> bool var<type>::operator!=(const var& v) const {
    return !(*this==v);
}

/* Output */
template<typename type> string var<type>::to_str_bounds(bool bounds, int prec) const {
    string str = param<type>::to_str_vals(false, prec);
    if (!bounds) {
        return str;
    }
    if(_lb->func_is_number() && _ub->func_is_number()){
        str += " ∈ [" + _lb->to_str(0,3) +"," + _ub->to_str(0,3) +"]^" + to_string(this->get_dim()) + "\n";
        return str;
    }
    str += " : ";
    auto space_size = str.size();
    if(this->_indices) {
        for (size_t i = 0; i < this->_dim[0]; i++) {
            if(i>0){
                str.insert(str.end(), space_size, ' ');
            }
            auto idx = this->get_id_inst(i);
            str += "(" + this->_indices->_keys->at(idx) + ") ∈ ";
            str += " [" + _lb->to_str(i,prec) + "," + _ub->to_str(i,prec) + "]\n";
            str += " \n";
        }
    }
    else {
        for (size_t idx = 0; idx < this->_dim[0]; idx++) {
            str += "["+to_string(idx) + "] = ";
            str += " [" + _lb->to_str(idx,prec) + "," + _ub->to_str(idx,prec) + "]\n";
            str += " \n";
        }
    }
    str+= ";\n";
    return str;
}

template<typename type> void var<type>::print() {
    cout << to_str_bounds(true,10);
}

    
template<typename type> void var<type>::print_bounds(bool bounds, int prec) {
    cout << to_str_bounds(bounds,prec);
}
    
var<Cpx> conj(const var<Cpx>& v){
    var<Cpx> newv(v);
    if(newv._is_conjugate){
        newv._name = newv._name.substr(newv._name.find("("),newv._name.find(")"));
    }
    else {
        newv._name = "conj("+newv._name+")";
    }
    newv._is_conjugate = !newv._is_conjugate;
    return newv;
}

var<Cpx> ang(const var<Cpx>& v){
    var<Cpx> newv(v);
    newv._is_angle = true;
    newv._name = "ang("+newv._name+")";
    return newv;
}

var<Cpx> sqrmag(const var<Cpx>& v){
    var<Cpx> newv(v);
    newv._is_sqrmag = true;
    newv._name = "|"+newv._name+"|²";
    return newv;
}

var<Cpx> real(const var<Cpx>& v){
    var<Cpx> newv(v);
    newv._is_real = true;
    newv._name = "real("+newv._name+")";
    return newv;
}
    
var<Cpx> imag(const var<Cpx>& v){
    var<Cpx> newv(v);
    newv._is_imag = true;
    newv._name = "imag("+newv._name+")";
    return newv;
}
    
template class var<bool>;
template class var<short>;
template class var<int>;
template class var<float>;
template class var<double>;
template class var<long double>;
template class var<complex<double>>;


    
}
