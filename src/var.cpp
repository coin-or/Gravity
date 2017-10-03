//
//  var.cpp
//  Gravity
//
//  Created by Hijazi, Hassan on 21/10/16.
//  note: Sdpvar needs to be tested (Guanglei). 
//
#include <gravity/var.h>
#define DebugOn(x) cout << x
#define DebugOff(x)

namespace gravity{

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
    
    template<typename type> var<type>::var(const string& name, param<type> lb, param<type> ub):var(name){
//        for (int i= 0; i<lb.get_dim(); i++) {
//            _lb->push_back(lb.eval(i));
//            _ub->push_back(ub.eval(i));
//            if (_lb->back() < param<type>::_range.first) {
//                param<type>::_range.first = _lb->back() ;
//            }
//            if (_ub->back() > param<type>::_range.second) {
//                param<type>::_range.second = _ub->back();
//            }
//        }
        _lb = lb.get_vals();
        _ub = ub.get_vals();
        param<type>::_ids = unique_ptr<vector<unsigned>>(new vector<unsigned>(lb.get_ids()));
        //TODO: update range first and second
    };
    
    template<typename type> var<type>::var(const string& name, param<type> sb):var(name){
        _ub = sb.get_vals();
        for (int i= 0; i<_ub->size(); i++) {
            _lb->push_back(-1.*_ub->at(i));
            DebugOff(_lb->at(i) << endl);
            if (_lb->back() < param<type>::_range.first) {
                param<type>::_range.first = _lb->back() ;
            }
        }
        DebugOff("Done" << endl);
        param<type>::_ids = unique_ptr<vector<unsigned>>(new vector<unsigned>(sb.get_ids()));
        //TODO: update range first and second
    };

    template<typename type> var<type>& var<type>::operator=(const var<type>& v){
        this->param<type>::operator=(v);
        param<type>::set_type(var_c);
        _lb = v._lb;
        _ub = v._ub;
        return *this;
    };
    
    template<typename type> var<type>& var<type>::operator=(var<type>&& v){
        this->param<type>::operator=(move(v));
        param<type>::set_type(var_c);
        _lb = v._lb;
        _ub = v._ub;
        return *this;
    };

    /* Modifiers */
    template<typename type> void   var<type>::set_size(size_t s, type val){
        param<type>::set_size(s,val);
        if (_lb->empty()) {
            _lb->resize(s, numeric_limits<type>::lowest());
        }
        else if(_lb->size() < s){
            _lb->resize(s, _lb->at(0));
        }
        if (_ub->empty()) {
            _ub->resize(s, numeric_limits<type>::max());
        }
        else if(_ub->size() < s){
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

        if(_lb == nullptr && _ub != nullptr){
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
        if(_lb != nullptr && _ub == nullptr){
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
        if(_lb != nullptr && _ub != nullptr){
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

    template<typename type>var<type> var<type>::from(const vector<Arc*>& arcs){
        var<type> res(this->_name);
        res.param<type>::operator=(param<type>::from(arcs));
        res.param<type>::set_type(var_c);
        res._lb = this->_lb;
        res._ub = this->_ub;        
        return res;
    }


    template<typename type>var<type> var<type>::to(const vector<Arc*>& arcs){
        var<type> res(this->_name);
        res.param<type>::operator=(param<type>::to(arcs));
        res.param<type>::set_type(var_c);
        res._lb = this->_lb;
        res._ub = this->_ub;
        return res;
    }

    template<typename type>var<type> var<type>::in(const vector<Arc*>& arcs){
        var<type> res(this->_name);
        res.param<type>::operator=(param<type>::in(arcs));
        res.param<type>::set_type(var_c);
        res._lb = this->_lb;
        res._ub = this->_ub;
        return res;
    }

    
    template<typename type>var<type> var<type>::in(const ordered_pairs& pairs){
        var<type> res(this->_name);
        res.param<type>::operator=(param<type>::in(pairs));
        res.param<type>::set_type(var_c);
        res._lb = this->_lb;
        res._ub = this->_ub;
        return res;
    }
    
    template<typename type>var<type> var<type>::from(const ordered_pairs& pairs){
        var<type> res(this->_name);
        res.param<type>::operator=(param<type>::from(pairs));
        res.param<type>::set_type(var_c);
        res._lb = this->_lb;
        res._ub = this->_ub;
        return res;
    }
    
    template<typename type>var<type> var<type>::to(const ordered_pairs& pairs){
        var<type> res(this->_name);
        res.param<type>::operator=(param<type>::to(pairs));
        res.param<type>::set_type(var_c);
        res._lb = this->_lb;
        res._ub = this->_ub;
        return res;
    }
    
    template<typename type>vector<var<type>> var<type>::pairs_in(const std::vector<std::vector<Node*>>& bags, unsigned size){
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
            res[i]._unique_id = make_tuple<>(res[i]._id,in_bags_, 0, i);
            res[i]._is_indexed = true;
        }
        set<vector<unsigned>> ids;
        for (auto &bag: bags){
            if (bag.size() < size) {
                continue;
            }
            vector<unsigned> ids_bag;
            for (int i = 0; i<size; i++) {
                ids_bag.push_back(bag[i]->ID);
            }
            if(ids.count(ids_bag)==0){
                ids.insert(ids_bag);
            }
            else {
                continue;
            }
            for (int i = 0; i<size-1; i++) {
                key = bag[i]->_name + "," + bag[i+1]->_name;
                assert(bag[i+2]->ID > bag[i]->ID);
                auto pp = param_::_indices->insert(make_pair<>(key,param_::_indices->size()));
                if(pp.second){//new index inserted
                    if(res[i]._indices->insert(make_pair<>(key,param_::_indices->size()-1)).second){
                        res[i]._dim++;
                    }
                    res[i]._ids->push_back(param_::_indices->size()-1);
                }
                else {
                    if(res[i]._indices->insert(make_pair<>(key,pp.first->second)).second){
                        res[i]._dim++;
                    }
                    res[i]._ids->push_back(pp.first->second);
                }
            }
            /* Loop back pair */
            key = bag[0]->_name + "," + bag[size-1]->_name;
            assert(bag[size-1]->ID > bag[0]->ID);
            auto pp = param_::_indices->insert(make_pair<>(key,param_::_indices->size()));
            if(pp.second){//new index inserted
                if(res[size-1]._indices->insert(make_pair<>(key,param_::_indices->size()-1)).second){
                    res[size-1]._dim++;
                }
                res[size-1]._ids->push_back(param_::_indices->size()-1);
            }
            else {
                if(res[size-1]._indices->insert(make_pair<>(key,pp.first->second)).second){
                    res[size-1]._dim++;
                }
                res[size-1]._ids->push_back(pp.first->second);
            }
        }
        return res;
    }
    
    template<typename type>vector<var<type>> var<type>::in(const std::vector<std::vector<Node*>>& bags, unsigned size){
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
            res[i]._unique_id = make_tuple<>(res[i]._id,in_bags_, 0, i);
            res[i]._is_indexed = true;
        }
        set<vector<unsigned>> ids;
        for (auto &bag: bags){
            if (bag.size() < size) {
                continue;
            }
            vector<unsigned> ids_bag;
            for (int i = 0; i<size; i++) {
                ids_bag.push_back(bag[i]->ID);
            }
            if(ids.count(ids_bag)==0){
                ids.insert(ids_bag);
            }
            else {
                continue;
            }

            for (int i = 0; i<size; i++) {
                key = bag[i]->_name;
                auto pp = param_::_indices->insert(make_pair<>(key,param_::_indices->size()));
                if(pp.second){//new index inserted
                    if(res[i]._indices->insert(make_pair<>(key,param_::_indices->size()-1)).second){
                        res[i]._dim++;
                    }
                    res[i]._ids->push_back(param_::_indices->size()-1);
                }
                else {
                    if(res[i]._indices->insert(make_pair<>(key,pp.first->second)).second){
                        res[i]._dim++;
                    }
                    res[i]._ids->push_back(pp.first->second);
                }
            }
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
