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

template<typename type>
var<type> var<type>::deep_copy() const{
    var<type> res;
    res.param<type>::operator=(this->param<type>::deep_copy());
    res.constant_::set_type(var_c);
    res._num_partns = make_shared<int>(*_num_partns);
    res._lb = make_shared<func<type>>();
    res._ub = make_shared<func<type>>();
    res._lb->deep_copy(*_lb);
    res._ub->deep_copy(*_ub);
    res._lift=_lift;
    return res;
}

template<typename type> var<type>& var<type>::operator=(const var<type>& v) {
    this->param<type>::operator=(v);
    constant_::set_type(var_c);
    _lb = v._lb;
    _ub = v._ub;
    _lift=v._lift;
    _lift_lb = v._lift_lb;
    _lift_ub = v._lift_ub;
    _in_SOC_partn=v._in_SOC_partn;
    _num_partns = v._num_partns;
    _cur_partn = v._cur_partn;
    _original_vars = v._original_vars;
    return *this;
};

template<typename type> var<type>& var<type>::operator=(var<type>&& v) {
    this->param<type>::operator=(move(v));
    constant_::set_type(var_c);
    _lb = move(v._lb);
    _ub = move(v._ub);
    _lift=v._lift;
    _lift_lb = v._lift_lb;
    _lift_ub = v._lift_ub;
    _in_SOC_partn=v._in_SOC_partn;
    _num_partns = v._num_partns;
    _cur_partn = v._cur_partn;
    _original_vars = move(v._original_vars);
    return *this;
};

/* Modifiers */
    
    /* Create a vector of variables indexed as pair of nodes from bags of size bag_size (WARNING assumes bags are unique in bags)*/
    template<typename type>
    vector<var<type>> var<type>::pairs_in_bags(const vector<pair<string,vector<Node*>>>& bags, size_t bag_size){
    vector<var<type>> res;
    vector<indices> ids_vec;
    string key, rev_key;
    res.resize(bag_size);
    ids_vec.resize(bag_size);
//    set<vector<Node*>> unique_bags;
    for (auto i = 0; i<bag_size; i++) {
        ids_vec[i] = *param_::_indices;
        ids_vec[i].set_name("pairs_"+to_string(i));
        if(ids_vec[i]._ids){
            ids_vec[i]._ids = nullptr;
        }
    }
    for (auto &bag: bags) {
        /* Make sure it's a bag with size=bag_size */
//        if (bag.size() == bag_size && unique_bags.insert(bag).second) {
        if (bag.second.size() == bag_size) {
            for (size_t i = 0; i< bag_size-1; i++) {
                key = bag.second[i]->_name + "," + bag.second[i+1]->_name;
                rev_key = bag.second[i+1]->_name + "," + bag.second[i]->_name;
                if(key.find("x")!=std::string::npos){
                    auto key1 = bag.second[i]->_name.substr(bag.second[i]->_name.find("[")+1);
                    key1 = key1.substr(0,key1.size()-1);
                    auto key2 = bag.second[i+1]->_name.substr(bag.second[i+1]->_name.find("[")+1);
                    key2 = key2.substr(0,key2.size()-1);
                    key = key1 + "," + key2;
                    rev_key = key2 + "," + key1;
                }
                if(ids_vec[i]._keys_map->count(key)!=0) {
                    ids_vec[i].add_ref(key);
                }
                else if(ids_vec[i]._keys_map->count(rev_key)!=0){
                    ids_vec[i].add_ref(rev_key);
                }
                else{
                    throw invalid_argument("key not found");
                }
            }
            /* Loop back pair */
            key = bag.second[0]->_name + "," + bag.second[bag_size-1]->_name;
            rev_key =  bag.second[bag_size-1]->_name + "," + bag.second[0]->_name;
            if(key.find("x")!=std::string::npos){
                auto key1 = bag.second[0]->_name.substr(bag.second[0]->_name.find("[")+1);
                key1 = key1.substr(0,key1.size()-1);
                auto key2 = bag.second[bag_size-1]->_name.substr(bag.second[bag_size-1]->_name.find("[")+1);
                key2 = key2.substr(0,key2.size()-1);
                key = key1 + "," + key2;
                rev_key = key2 + "," + key1;
            }
            if(ids_vec[bag_size-1]._keys_map->count(key)!=0){
                ids_vec[bag_size-1].add_ref(key);
            }
            else if(ids_vec[bag_size-1]._keys_map->count(rev_key)!=0){
                ids_vec[bag_size-1].add_ref(rev_key);
            }
            else{
                throw invalid_argument("key not found");
            }
        }
    }
    for (auto i = 0; i<bag_size; i++) {
        res[i] = this->in(ids_vec[i]);
        if(res[i].is_complex()){
            if(res[i]._real){
                res[i]._real->_name = "R_"+res[i]._name;
                res[i]._real->_indices = res[i]._indices;
            }
            if(res[i]._imag){
                res[i]._imag->_name = "Im_"+res[i]._name;;
                res[i]._imag->_indices = res[i]._indices;
            }
            if(res[i]._ang){
                res[i]._ang->_name += "pairs_"+to_string(i);
                res[i]._ang->_indices = res[i]._indices;
            }
            if(res[i]._mag){
                res[i]._mag->_name += "pairs_"+to_string(i);
                res[i]._mag->_indices = res[i]._indices;
            }
        }

    }
    return res;
}
    
/* Create a vector of variables indexed based on nodes from bags of size bag_size */
template<typename type>
vector<var<type>> var<type>::in_bags(const vector<pair<string,vector<Node*>>>& bags, size_t bag_size){
    vector<var> res;
    vector<indices> ids_vec;
    res.resize(bag_size);
    ids_vec.resize(bag_size);
    map<string,vector<Node*>> unique_bags;
    for (auto i = 0; i<bag_size; i++) {
        ids_vec[i] = *param_::_indices;
        ids_vec[i].set_name("nodes_"+to_string(i));
    }
    for (auto &bag: bags) {
        /* Make sure it's a new bag with size=bag_size */
        if (bag.second.size() == bag_size && unique_bags.insert(bag).second) {
            for (size_t i = 0; i< bag_size; i++) {
                auto key = bag.second[i]->_name;
                if(key.find("x")!=std::string::npos){
                    key = key.substr(key.find("[")+1);
                    key = key.substr(0,key.size()-1);
                }
                ids_vec[i].add_ref(key);
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
    if (_lb->func_is_number()) {
        return _lb->_val->at(0);
    }
    return _lb->eval(i);
};

template<typename type>
type    var<type>::get_ub(size_t i) const {
    if (_ub->func_is_number()) {
        return _ub->_val->at(0);
    }
    return _ub->eval(i);
};
    
    template<typename type>
    type    var<type>::get_lb(const string& key) const {
        auto i = this->_indices->_keys_map->at(key);
        if (_lb->func_is_number()) {
            return _lb->_val->at(0);
        }
        return _lb->eval(i);
    };
    
    template<typename type>
    type    var<type>::get_ub(const string& key) const {
        auto i = this->_indices->_keys_map->at(key);
        if (_ub->func_is_number()) {
            return _ub->_val->at(0);
        }
        return _ub->eval(i);
    };
    
    template<typename type>
    param<type>    var<type>::get_lb() const {
        if(!_lift)
            return *static_pointer_cast<param<type>>(_lb->_params->begin()->second.first);
        param<type> lb(this->_name+"-lb");
        _lb->eval_all();
        if(!_lb->func_is_number()){
            lb.index_in(*this->_indices);
        }
        lb._val = _lb->_val;
        lb._range = _lb->_range;
        lb._dim[0] = _lb->_dim[0];
        lb._dim[1] = _lb->_dim[1];
        return lb;
    };


/* If this is a lifted variable lifted(xy)= xy, return the lowerbound on x*/
    template<typename type>
    shared_ptr<param<type>>    var<type>::get_bilinear_lb1() const{
        assert(_lift);
        auto lson = static_pointer_cast<func<type>>(_lb->_expr->get_lson());
        auto rson = static_pointer_cast<func<type>>(lson->_expr->get_rson());
        /* prod_b2 = (lb1*ub2);*/
        for (auto p_it: *rson->_params) {
            if (p_it.first.find("lb") != string::npos) {
                return static_pointer_cast<param<type>>(p_it.second.first);
            }
        }
        return nullptr;
    }

/* If this is a lifted variable lifted(x^2)= x^2, return the lowerbound on x*/
//ub = gravity::max(gravity::max(prod_b1,prod_b2).in(unique_ids),prod_b3);
template<typename type>
shared_ptr<param<type>>    var<type>::get_square_lb() const{
    assert(_lift);
    auto lson = static_pointer_cast<func<type>>(_ub->_expr->get_lson());
    auto lson2 = static_pointer_cast<func<type>>(lson->_expr->get_lson());
    return static_pointer_cast<param<type>>(lson2->_params->begin()->second.first);
}

template<typename type>
shared_ptr<param<type>>    var<type>::get_square_ub() const{
    assert(_lift);
    auto rson = static_pointer_cast<func<type>>(_ub->_expr->get_rson());
    return static_pointer_cast<param<type>>(rson->_params->begin()->second.first);
}

/* If this is a lifted variable lifted(xy)= xy, return the lowerbound on y*/
    template<typename type>
    shared_ptr<param<type>> var<type>::get_bilinear_lb2() const{
        assert(_lift);
        auto rson = static_pointer_cast<func<type>>(_lb->_expr->get_rson());
        auto lson = static_pointer_cast<func<type>>(rson->_expr->get_lson());
        /* prod_b3 = (ub1*lb2);*/
        for (auto p_it: *lson->_params) {
            if (p_it.first.find("lb") != string::npos) {
                return static_pointer_cast<param<type>>(p_it.second.first);
            }
        }
        return nullptr;
    }

/* If this is a lifted variable lifted(xy)= xy, return the upperbound on x*/
    template<typename type>
    shared_ptr<param<type>>    var<type>::get_bilinear_ub1() const{
        assert(_lift);
        auto rson = static_pointer_cast<func<type>>(_lb->_expr->get_rson());
        auto lson = static_pointer_cast<func<type>>(rson->_expr->get_lson());
        /* prod_b3 = (ub1*lb2);*/
        for (auto p_it: *lson->_params) {
            if (p_it.first.find("ub") != string::npos) {
                return static_pointer_cast<param<type>>(p_it.second.first);
            }
        }
        return nullptr;
    }

/* If this is a lifted variable lifted(xy)= xy, return the upperbound on y*/
    template<typename type>
    shared_ptr<param<type>>    var<type>::get_bilinear_ub2() const{
        assert(_lift);
        auto lson = static_pointer_cast<func<type>>(_lb->_expr->get_lson());
        auto rson = static_pointer_cast<func<type>>(lson->_expr->get_rson());
        /* prod_b2 = (lb1*ub2);*/
        for (auto p_it: *rson->_params) {
            if (p_it.first.find("ub") != string::npos) {
                return static_pointer_cast<param<type>>(p_it.second.first);
            }
        }
        return nullptr;
    }
    
    template<typename type>
    param<type>    var<type>::get_ub() const {
        if(!_lift)
            return *static_pointer_cast<param<type>>(_ub->_params->begin()->second.first);
        param<type> ub(this->_name+"-ub");
        _ub->eval_all();
        if(!_ub->func_is_number()){
            ub.index_in(*this->_indices);
        }
        ub._val = _ub->_val;
        ub._range = _ub->_range;
        ub._dim[0] = _ub->_dim[0];
        ub._dim[1] = _ub->_dim[1];
        return ub;
    };
    
//    template <typename type>
//    template<typename T,
//    typename std::enable_if<is_arithmetic<T>::value>::type*>
//    bool var<type>::is_bounded_above(size_t i) const {
//        return (_ub->eval(i)!=numeric_limits<type>::max());
//    };
    
//    template <typename type>
//    template<typename T,
//    typename std::enable_if<is_arithmetic<T>::value>::type*>
//    bool var<type>::is_bounded_below(size_t i) const {
//        return (_lb->eval(i)!=numeric_limits<type>::lowest());
//    };
    
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
    if(this->is_indexed()){
        _lb->set_val(this->get_id_inst(),val);
        _lb->update_range(val);
        this->update_range(val);
        _lb->_evaluated = true;
    }
    else {
        _lb->set_val(val);
        param<type>::_range->first = val;
    }
    _lb->_evaluated = true;
}
    
template<typename type> void  var<type>::set_lb(const string& key, type val){
    auto key_it = this->_indices->_keys_map->find(key);
    /* First check if variable has this key */
    if(key_it== this->_indices->_keys_map->end()){
        throw invalid_argument("in set_lb(string, val), unknown key" + key);
    }
    /* Also check if lower-bound has this key */
    key_it = _lb->_indices->_keys_map->find(key);
    if(key_it== _lb->_indices->_keys_map->end()){
        throw invalid_argument("in set_ub(string, val), unknown key " + key);
    }
    _lb->eval_all();
    _lb->_val->at(key_it->second) = val;
    _lb->update_range(val);
    this->update_range(val);
}

template<typename type> void  var<type>::set_ub(const string& key, type val){
    auto key_it = this->_indices->_keys_map->find(key);
    /* First check if variable has this key */
    if(key_it== this->_indices->_keys_map->end()){
        throw invalid_argument("in set_ub(string, val), unknown key " + key);
    }
    /* Also check if upper-bound has this key */
    key_it = _ub->_indices->_keys_map->find(key);
    if(key_it== _ub->_indices->_keys_map->end()){
        throw invalid_argument("in set_ub(string, val), unknown key " + key);
    }
    _ub->eval_all();
    _ub->_val->at(key_it->second) = val;
    _ub->update_range(val);
    this->update_range(val);
}

    
template<typename type> void   var<type>::set_ub(type val) {
    if(this->is_indexed()){
        _ub->set_val(this->get_id_inst(),val);
        _ub->update_range(val);
        this->update_range(val);
        _ub->_evaluated = true;
    }
    else {
        _ub->set_val(val);
        param<type>::_range->second = val;
    }
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
template<typename type> string var<type>::to_str_bounds(bool bounds, int prec) {
    string str = param<type>::to_str_vals(false, prec);
    if (!bounds) {
        return str;
    }
    if(_lb->func_is_number() && _ub->func_is_number()){
        if(this->_is_relaxed)
            str += " ∈ {" + _lb->to_str(0,3) +"," + _ub->to_str(0,3) +"}^" + to_string(this->get_dim()) + "\n";
        else
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
            if(this->_is_relaxed)
                str += " {" + _lb->to_str(i,prec) + "," + _ub->to_str(i,prec) + "}\n";
            else
                str += " [" + _lb->to_str(i,prec) + "," + _ub->to_str(i,prec) + "]\n";
            str += " \n";
        }
    }
    else {
        for (size_t idx = 0; idx < this->_dim[0]; idx++) {
            str += "["+to_string(idx) + "] = ";
            if(this->_is_relaxed)
                str += " {" + _lb->to_str(idx,prec) + "," + _ub->to_str(idx,prec) + "}\n";
            else
                str += " [" + _lb->to_str(idx,prec) + "," + _ub->to_str(idx,prec) + "]\n";
            str += " \n";
        }
    }
    str+= ";\n";
    return str;
}

template<typename type> void var<type>::print() {
    print(16);
}

    
    
template<typename type> void var<type>::print(int prec) {
    cout << to_str_bounds(true,prec);
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
