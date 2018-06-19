//
//  var.cpp
//  Gravity
//
//  Created by Hijazi, Hassan on 21/10/16.
//  note: Sdpvar needs to be tested (Guanglei).
//
#include <gravity/var.h>
#include <gravity/func.h>
#include <gravity/Net.h>
#define DebugOn(x) cout << x
#define DebugOff(x)

namespace gravity {

template<typename type> var<type>::var():param<type>() {
    param<type>::set_type(var_c);
};

template<typename type> var<type>::var(const string& name):param<type>(name) {
    param<type>::set_type(var_c);
    if (this->is_binary()) {
        _lb = make_shared<func_>(constant<>(0));
        _ub = make_shared<func_>(constant<>(1));
    }
    else {
        _lb = make_shared<func_>(constant<type>(numeric_limits<type>::lowest()));
        _ub = make_shared<func_>(constant<type>(numeric_limits<type>::max()));
    }
};
    
template<typename type> var<type>::var(const string& name, Sign s):param<type>(name) {
    param<type>::set_type(var_c);
    _lb = make_shared<func_>();
    _ub = make_shared<func_>();
    if (s==non_neg_ || s==pos_) {
        add_lb_only(0);
    }
    else if (s==non_pos_ || s==neg_) {
        add_ub_only(0);
    }
};

template<typename type> var<type>::var(const var<type>& v):param<type>(v) {
    param<type>::set_type(var_c);
    _lb = v._lb;
    _ub = v._ub;
    _in_q_cone = v._in_q_cone;
    _psd = v._psd;
    param<type>::_is_relaxed = v._is_relaxed;
};

template<typename type> var<type>::var(var<type>&& v):param<type>(v) {
    param<type>::set_type(var_c);
    _lb = move(v._lb);
    _ub = move(v._ub);
    param<type>::_is_relaxed = v._is_relaxed;
};

template<typename type> var<type>::var(const string& name, type lb, type ub):var(name) {
    _lb = make_shared<func_>(constant<type>(lb));
    _ub = make_shared<func_>(constant<type>(ub));
    if (lb < param<type>::_range->first) {
        param<type>::_range->first = lb;
    }
    if (ub > param<type>::_range->second) {
        param<type>::_range->second = ub;
    }
};
    
template<typename type> var<type>::var(const string& name, func_&& lb, func_&& ub):var(name) {
    _lb = make_shared<func_>(move(lb));
    _ub = make_shared<func_>(move(ub));
}
    
template<typename type> var<type>::var(const string& name, const func_& lb, const func_& ub):var(name) {
    _lb = make_shared<func_>(lb);
    _ub = make_shared<func_>(ub);
}
    
//template<typename type> var<type>::var(const string& name, const param<type>& lb, const param<type>& ub):var(name) {
//    _lb = make_shared<func_>(lb);
//    _ub = make_shared<func_>(ub);
//    auto nb_inst = _lb->_nb_instances;
//    _lb->_val->resize(nb_inst);
//    _ub->_val->resize(nb_inst);
////    this->_val->resize(nb_inst);
//    this->_dim[0] = nb_inst;
//    assert(nb_inst == _ub->_nb_instances);
//    param<type>::_range->first = lb._range->first;
//    param<type>::_range->second = ub._range->second;
////    unsigned i = 0;    
////    for (auto &p: *lb.get_indices()) {
////    this->_indices = lb.get_indices();
////    this->_rev_indices = lb.get_rev_indices();
////    for(unsigned i = 0;i<nb_inst; i++){
////        _lb->_val->at(i) = lb._val->at(i);
////        _ub->_val->at(i) = ub._val->at(i);
////    }
////    _lb->_evaluated = true;
////    _ub->_evaluated = true;
//}
    
    
    
//template<typename type> var<type>::var(const string& name, func_&& lb, func_&& ub):var(name) {
//    _lb = make_shared<func_>(move(lb));
//    _ub = make_shared<func_>(move(ub));
//    param<type>::_range->first = lb.get_all_range()->first;
//    param<type>::_range->second = ub.get_all_range()->second;
//}

template<typename type> var<type>::var(const string& name, const param<type>& sb):var(name) {
    _lb = make_shared<func_>(-1*sb);
    _ub = make_shared<func_>(sb);
    auto nb_inst = _lb->_nb_instances;
    this->_dim[0] = nb_inst;
    _lb->_val->resize(nb_inst);
    _ub->_val->resize(nb_inst);
    param<type>::_range->first = min(-1*sb._range->first, -1*sb._range->second);
    param<type>::_range->second = sb._range->second;
//    unsigned i = 0;
//    for(unsigned i = 0;i<nb_inst; i++){
//        //        auto index = param_::_indices->size();
//        //        param_::_indices->insert(make_pair<>(p.first, index));
//        //        param_::_rev_indices->resize(max(param_::_rev_indices->size(),index+1));
//        //        param_::_rev_indices->at(index) = p.first;
//        _lb->_val->at(i) = -1*sb._val->at(i);
//        _ub->_val->at(i) = sb._val->at(i);
//    }
//    _lb->_evaluated = true;
//    _ub->_evaluated = true;
};
    
//    template<typename type> var<type>::var(const string& name, func_&& sb):var(name) {
//        _lb = make_shared<func_>(sb*-1);
//        _ub = make_shared<func_>(move(sb));
//        param<type>::_range->first = min(-1*sb.get_all_range()->first,-1*sb.get_all_range()->second);
//        param<type>::_range->second = sb.get_all_range()->second;
//        param<type>::_ids = sb._ids;
//    };



template<typename type> var<type>& var<type>::operator=(const var<type>& v) {
    this->param<type>::operator=(v);
    param<type>::set_type(var_c);
    _lb = v._lb;
    _ub = v._ub;
    param<type>::_is_relaxed = v._is_relaxed;
    return *this;
};

template<typename type> var<type>& var<type>::operator=(var<type>&& v) {
    this->param<type>::operator=(move(v));
    param<type>::set_type(var_c);
    _lb = v._lb;
    _ub = v._ub;
    param<type>::_is_relaxed = v._is_relaxed;
    return *this;
};

/* Modifiers */
    
template<typename type> void   var<type>::set_size(vector<size_t> dims) {
    param<type>::set_size(dims);
};

template<typename type> void   var<type>::set_size(size_t s, type val) {
    param<type>::set_size(s,val);
};

template<typename type>
type    var<type>::get_lb(size_t i) const {
    if (_lb->is_number()) {
        return _lb->_val->at(0);
    }
    return _lb->eval(i);
//    unsigned index = 0;
//    if (_lb->_ids->at(0).empty()) {
//        index = i;
//    }
//    else {
//        index = _lb->_ids->at(0).at(i);
//    }
//    return poly_eval(_lb->_params->begin()->second.first, i);
};

template<typename type>
type    var<type>::get_ub(size_t i) const {
    if (_ub->is_number()) {
        return _ub->_val->at(0);
    }
//    unsigned index = 0;
//    if (_ub->_ids->at(0).empty()) {
//        index = i;
//    }
//    else {
//        index = _ub->_ids->at(0).at(i);
//    }
//    return _ub->get_val(index);
    return _ub->eval(i);
};
    
    template<typename type>
    bool var<type>::is_bounded_above(int i) const {
        return (_ub!=nullptr && _ub->eval(i)!=numeric_limits<type>::max());
    };
    
    template<typename type>
    bool var<type>::is_bounded_below(int i) const {
        return (_lb!=nullptr && _lb->eval(i)!=numeric_limits<type>::lowest());
    };
    
    template<typename type>
    bool var<type>::is_constant(int i) const {
        return (is_bounded_below() && is_bounded_above() && _lb->eval(i)==_ub->eval(i));
    };
    
    template<typename type>
    Sign var<type>::get_sign(int idx) const {
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
    *_ub = constant<type>(numeric_limits<type>::max());
    param<type>::_range->first = val;
    param<type>::_range->second = numeric_limits<type>::max();
}

template<typename type> void   var<type>::add_ub_only(type val) {
    *_lb = constant<type>(numeric_limits<type>::lowest());
    *_ub = constant<type>(val);
    param<type>::_range->first = numeric_limits<type>::lowest();
    param<type>::_range->second = val;
}

template<typename type> void   var<type>::set_lb(int i, type val) {
    *_lb = constant<type>(val);
    param<type>::_range->first = val;
}

template<typename type> void   var<type>::set_ub(int i, type val) {
    *_ub = constant<type>(val);
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
template<typename type> void var<type>::print(bool bounds) const {
    param<type>::print(false);
    if (!bounds) {
        return;
    }
    size_t idx = 0;
    cout << " in: \n";

    if(_lb == nullptr && _ub != nullptr) {
        if (param_::_is_indexed) {
            idx = param_::get_id_inst();
            cout << "(" << idx << ") = ";
            cout << " [ -inf , " << (*_ub->_val)[idx] << " ]";
        }
        else {
            for(int i = 0 ; i < param_::get_dim(); i++) {
                cout << "(" << idx << ") = ";
                cout << " [ -inf , " << (*_ub->_val)[i] << " ]\n";
            }
        }
        cout << ";\n";
    }
    if(_lb != nullptr && _ub == nullptr) {
        if (param_::_is_indexed) {
            idx = param_::get_id_inst();
            cout << " [ " << (*_lb->_val)[idx] << ", inf ]";
        }
        else {
            for(int i = 0 ; i < param_::get_dim(); i++) {
                cout << "(" << i << ") = ";
                cout << " [ " << (*_lb->_val)[i] << ", inf ]\n";
            }
        }
        cout << ";\n";
    }
    if(_lb != nullptr && _ub != nullptr) {
//        if (param_::_is_indexed) {
//            if (_lb->_val->size() >1){
//                idx = param_::get_id_inst();
//                cout << " [ " << (*_lb->_val)[idx] << ", " << (*_ub->_val)[idx] << "]";
//            }
//            else
//                cout << " [ " << (*_lb->_val)[0] << ", " << (*_ub->_val)[0] << "]";
//        }
//        else {
            for(int i = 0 ; i < param_::get_dim(); i++) {
                cout << "(" << param_::_rev_indices->at(i) << ") = ";
                cout << " [ " << _lb->eval(i) << ", " << _ub->eval(i) << "]\n";
            }
//        }
        cout << ";\n";
    }
};

template<typename type> sdpvar<type>::sdpvar():var<type>() {
    var<type>::set_type(sdpvar_c);
};

template<typename type> sdpvar<type>::sdpvar(const string& name):var<type>(name) {
    var<type>::set_type(sdpvar_c);
};

template<typename type> sdpvar<type>::sdpvar(const sdpvar<type>& v):var<type>(v) {
    var<type>::set_type(sdpvar_c);
    _symdim = v._symdim;
};

template<typename type> sdpvar<type>::sdpvar(sdpvar<type>&& v):var<type>(v) {
    var<type>::set_type(sdpvar_c);
    _symdim = v._symdim;
};

// sdpvar
template<typename type> sdpvar<type>::sdpvar(const string& name, type lb, type ub):sdpvar<type>(name) {
    var<type>(name, lb, ub);
};


// Operators
template<typename type> bool sdpvar<type>::operator==(const sdpvar& v) const {
    return (param<type>::operator==(v));
};

template<typename type> bool sdpvar<type>::operator!=(const sdpvar& v) const {
    return !(*this==v);
}

/* output */
template<typename type> void sdpvar<type>::print(bool vals) const {
    param<type>::print(vals);
};

template<typename type>vector<var<type>> var<type>::pairs_in_directed(Net& net, const std::vector<std::vector<Node*>>& bags, unsigned size){
    vector<var> res;
    string key;
    res.resize(size);
    for (int i = 0; i<size; i++) {
        res[i]._id = this->_id;
        res[i]._vec_id = this->_vec_id;
        res[i]._is_vector = this->_is_vector;
        res[i]._intype = this->_intype;
        res[i]._range = this->_range;
        res[i]._val = this->_val;
        res[i]._lb = this->_lb;
        res[i]._ub = this->_ub;
        res[i]._name = this->_name+to_string(i);
        res[i]._unique_id = make_tuple<>(res[i]._id,in_,typeid(type).hash_code(), 0, i);
        res[i]._is_indexed = true;
        res[i]._rev_indices = this->_rev_indices;
        res[i]._indices = this->_indices;
    }
    set<vector<unsigned>> ids;
    for (auto &bag: bags) {
        if (bag.size() != size) {
            continue;
        }
        vector<unsigned> ids_bag;
        for (int i = 0; i<size; i++) {
            ids_bag.push_back(bag[i]->_id);
        }
        if(ids.count(ids_bag)==0) {
            ids.insert(ids_bag);
        }
        else {
            continue;
        }
        for (int i = 0; i< size-1; i++) {
            if(net.get_directed_arc(bag[i]->_name, bag[i+1]->_name)!=nullptr) {
                key = bag[i]->_name + "," + bag[i+1]->_name;
            }
            else {
                key = bag[i+1]->_name + "," + bag[i]->_name;
            }
            auto index = param_::_indices->size();
            auto pp = param_::_indices->insert(make_pair<>(key,index));
            if(pp.second) { //new index inserted
                this->_val->resize(max(this->_val->size(),index+1));
                this->_dim[0] = max(this->_dim[0],this->_val->size());
                param_::_rev_indices->resize(max(param_::_rev_indices->size(),index+1));
                param_::_rev_indices->at(index) = key;
                res[i]._ids->at(0).push_back(index);
            }
            else {
                res[i]._ids->at(0).push_back(pp.first->second);
            }
        }
        /* Loop back pair */
        if(net.get_directed_arc(bag[0]->_name, bag[size-1]->_name)!=nullptr) {
            key = bag[0]->_name + "," + bag[size-1]->_name;
        }
        else{
            key = bag[size-1]->_name + "," + bag[0]->_name;
        }
        auto index = param_::_indices->size();
        auto pp = param_::_indices->insert(make_pair<>(key,index));
        if(pp.second) { //new index inserted
            this->_val->resize(max(this->_val->size(),index+1));
            this->_dim[0] = max(this->_dim[0],this->_val->size());
            param_::_rev_indices->resize(max(param_::_rev_indices->size(),index+1));
            param_::_rev_indices->at(index) = key;
            res[size-1]._ids->at(0).push_back(index);
        }
        else {
            res[size-1]._ids->at(0).push_back(pp.first->second);
        }
    }
    for (int i = 0; i<size; i++) {
        res[i]._dim[0]=res[i]._ids->at(0).size();
    }
    return res;
}


template<typename type>vector<var<type>> var<type>::pairs_in(const std::vector<std::vector<Node*>>& bags, unsigned size) {
    vector<var> res;
    string key;
    res.resize(size,(this->_name));
    for (int i = 0; i<size; i++) {
        res[i]._id = this->_id;
        res[i]._vec_id = this->_vec_id;
        res[i]._intype = this->_intype;
        res[i]._range = this->_range;
        res[i]._val = this->_val;
        res[i]._lb = this->_lb;
        res[i]._ub = this->_ub;
        res[i]._name += "_in_bags_"+to_string(i);
        res[i]._unique_id = make_tuple<>(res[i]._id,in_,typeid(type).hash_code(), 0, i);
        res[i]._is_indexed = true;
    }
    set<vector<unsigned>> ids;
    for (auto &bag: bags) {
        if (bag.size() != size) {
            continue;
        }
        vector<unsigned> ids_bag;
        for (int i = 0; i<size; i++) {
            ids_bag.push_back(bag[i]->_id);
        }
        if(ids.count(ids_bag)==0) {
            ids.insert(ids_bag);
        }
        else {
            continue;
        }
        for (int i = 0; i< size-1; i++) {
            key = bag[i]->_name + "," + bag[i+1]->_name;
//            assert(bag[i+2]->ID > bag[i]->ID);
            auto index = param_::_indices->size();
            auto pp = param_::_indices->insert(make_pair<>(key,index));
            if(pp.second) { //new index inserted
                param_::_rev_indices->resize(max(param_::_rev_indices->size(),index+1));
                param_::_rev_indices->at(index) = key;
                if(res[i]._indices->insert(make_pair<>(key,index)).second) {
                    res[i]._dim[0]++;
                }
                res[i]._ids->at(0).push_back(index);
            }
            else {
                if(res[i]._indices->insert(make_pair<>(key,pp.first->second)).second) {
                    res[i]._dim[0]++;
                }
                res[i]._ids->at(0).push_back(pp.first->second);
            }
        }
        /* Loop back pair */
        key = bag[0]->_name + "," + bag[size-1]->_name;
//        assert(bag[size-1]->ID > bag[0]->ID);
        auto index = param_::_indices->size();
        auto pp = param_::_indices->insert(make_pair<>(key,index));
        if(pp.second) { //new index inserted
            param_::_rev_indices->resize(max(param_::_rev_indices->size(),index+1));
            param_::_rev_indices->at(index) = key;
            if(res[size-1]._indices->insert(make_pair<>(key,index)).second) {
                res[size-1]._dim[0]++;
            }
            res[size-1]._ids->at(0).push_back(index);
        }
        else {
            if(res[size-1]._indices->insert(make_pair<>(key,pp.first->second)).second) {
                res[size-1]._dim[0]++;
            }
            res[size-1]._ids->at(0).push_back(pp.first->second);
        }
    }
    return res;
}

template<typename type>vector<var<type>> var<type>::in_bags(const std::vector<std::vector<Node*>>& bags, unsigned size) {
    vector<var> res;
    string key;
    res.resize(size);
    for (int i = 0; i<size; i++) {
        res[i]._id = this->_id;
        res[i]._vec_id = this->_vec_id;
        res[i]._intype = this->_intype;
        res[i]._is_vector = this->_is_vector;
        res[i]._range = this->_range;
        res[i]._val = this->_val;
        res[i]._lb = this->_lb;
        res[i]._ub = this->_ub;
        res[i]._name = this->_name+to_string(i);//+"_in_bags_"+to_string(i);
        res[i]._unique_id = make_tuple<>(res[i]._id,in_,typeid(type).hash_code(), 0, i);
        res[i]._is_indexed = true;
        res[i]._rev_indices = this->_rev_indices;
        res[i]._indices = this->_indices;
    }
    set<vector<unsigned>> ids;
    for (auto &bag: bags) {
        if (bag.size() != size) {
            continue;
        }
        vector<unsigned> ids_bag;
        for (int i = 0; i<size; i++) {
            ids_bag.push_back(bag[i]->_id);
        }
        if(ids.count(ids_bag)==0) {
            ids.insert(ids_bag);
        }
        else {
            continue;
        }

        for (int i = 0; i<size; i++) {
            key = bag[i]->_name;
            auto index = param_::_indices->size();
            auto pp = param_::_indices->insert(make_pair<>(key,index));
            if(pp.second) { //new index inserted
                this->_val->resize(max(this->_val->size(),index+1));
                this->_dim[0] = max(this->_dim[0],this->_val->size());
                param_::_rev_indices->resize(max(param_::_rev_indices->size(),index+1));
                param_::_rev_indices->at(index) = key;
                res[i]._ids->at(0).push_back(index);
            }
            else {
                res[i]._ids->at(0).push_back(pp.first->second);
            }
        }
    }
    for (int i = 0; i<size; i++) {
        res[i]._dim[0]=res[i]._ids->at(0).size();
    }
    return res;
}


    
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
}
