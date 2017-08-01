//
//  var.h
//  Gravity
//
//  Created by Hijazi, Hassan on 21/10/16.
//
//

#ifndef var_h
#define var_h

#include <gravity/param.h>
#include <stdio.h>
#include <string>
#include <set>
#include <list>
#include <limits>

//class func;
using namespace std;

/** Backbone class for parameter */
class var_ {
public:
    virtual ~var_(){};
};

/** A variable can be a bool, a short, an int, a float or a double*/
template<typename type = double>

// define variable as a parameter with bounds
class var: public param<type>, public var_{
    
public:
    shared_ptr<vector<type>>    _lb; /**< Lower Bound */
    shared_ptr<vector<type>>    _ub; /**< Upper Bound */
    /* Constructors */
    //@{
    /** Unbounded variable constructor */
    var();
    ~var(){};
    var(const string& name);
    var(const var<type>& v);
    var(var<type>&& v);
    //@}
    
    //@{
    /** Bounded variable constructor */
    var(const string& name, type lb, type ub);
    //@}
        
    template<typename... Args>
    var operator()(size_t t1, Args&&... args){
        var res(this->_name);
        res._id = this->_id;
        res._vec_id = this->_vec_id;
        res._intype = this->_intype;
        res._range = this->_range;
        res._val = this->_val;
        res._lb = this->_lb;
        res._ub = this->_ub;
        list<size_t> indices;
        indices = {forward<size_t>(args)...};
        //cout << "indices: "  << endl;
        //for(list<size_t>::iterator iter = indices.begin(); iter != indices.end(); iter++){
        //    cout<<*iter<<endl;
        //}

        indices.push_front(t1);
        string key;
        auto it = indices.begin();
        for (size_t i= 0; i<indices.size(); i++) {
            key += to_string(*it);
            if (i<indices.size()-1) {
                key += ",";
            }
            it++;
        }
        auto it2 = param_::_indices->find(key);
        if (it2 == param_::_indices->end()) {
            res._indices->insert(make_pair<>(key,param_::_indices->size()));
            param_::_indices->insert(make_pair<>(key,param_::_indices->size()));
            //            throw invalid_argument("Index " + key + "not found in" + param<type>::to_str()+".\n");
        }
        else {
            size_t idx = param_::_indices->at(key);
            res._indices->insert(make_pair<>(key,idx));
            res._dim = 1;
        }
        res._name += "["+key+"]";
        res._is_indexed = true;
        return res;
    }
    
    /* Querries */
    
    type    get_lb(size_t i = 0) const{
        if (_lb->size() <= i) {
            throw out_of_range("get_lb(size_t i, index: " + to_string(i) + ")\n");
        }
        return _lb->at(i);
    };
    
    type    get_ub(size_t i = 0) const{
        if (_ub->size() <= i) {
            throw out_of_range("get_ub(size_t i), index: " + to_string(i)+ ")\n");
        }
        return _ub->at(i);
    };

    
    bool is_bounded_above(int i = 0) const{
        return (_ub!=nullptr &&  i < _ub->size() && _ub->at(i)!=numeric_limits<type>::max());
    };

    bool is_bounded_below(int i = 0) const{
        return (_lb!=nullptr &&  i < _lb->size() && _lb->at(i)!=numeric_limits<type>::min());
    };
    
    bool is_constant(int i=0) const{
        return (is_bounded_below() && is_bounded_above() && _lb->at(i)==_ub->at(i));
    };
    
    Sign get_sign(int idx = 0) const{
        if (_lb->at(idx) == 0 && _ub->at(idx) == 0) {
            return zero_;
        }
        if (_ub->at(idx) < 0) {
            return neg_;
        }
        if (_lb->at(idx) > 0) {
            return pos_;
        }
        if (_ub->at(idx) == 0){
            return non_pos_;
        }
        if (_lb->at(idx) == 0) {
            return non_neg_;
        }
        return unknown_;
    }
    
    /* Modifiers */
    void    set_size(size_t s, type val = 0);
    
    void    add_bounds(type lb, type ub);
    void    add_lb_only(type v); /**< Adds a new lower bound and infinity in the corresponding upper bound*/
    void    add_ub_only(type v); /**< Adds a new upper bound and -infinity in the corresponding lower bound*/
    
    void    set_lb(int i, type v);
    void    set_ub(int i, type v);
    
    /* Operators */
    var& operator=(type v){
        param<type>::_val->push_back(v);
        param<type>::update_range(v);
        param<type>::_dim++;
        return *this;
    }
    
    bool operator==(const var& v) const;
    bool operator!=(const var& v) const;
    var& operator^(size_t d){
        set_size(d);
        return *this;
    }
    
    var tr() const{
        auto v = var(*this);
        v._is_transposed = true;
        return v;
    }
    
    /* Output */
    void print(bool bounds=false) const;
    
};

template<typename type>
var<type> all(const var<type>& p){
    auto pp = var<type>(p);
    pp._is_vector = true;
    return pp;
}


// In contrast to a general variable, indices of an SDP variable should be
// recorded using _sdpindices, moreover the size of an SDP variable is d.  
template<typename type = double>
//class sdpvar: public param<type>, public var_{
class sdpvar: public var<type>{

public:
    size_t _symdim=0;
    //@{
    /** Unbounded sdp-variable constructor */
    sdpvar();
    ~sdpvar(){};

    sdpvar(const string& name);
    sdpvar(const sdpvar<type>& v);
    sdpvar(sdpvar<type>&& v);
    //@}

    /** bounded sdp-variable constructor */
    sdpvar(const string& name, type lb, type ub);
    
    template<typename... Args>
    sdpvar operator()(size_t t1, Args&&... args){
        sdpvar res(this->_name);
        res._id = this->_id;
        res._vec_id = this->_vec_id;
        res._intype = this->_intype;
        res._val = this->_val;
        res._range = this->_range;
        res._lb = this->_lb;
        res._ub = this->_ub;
        list<size_t> indices;
        indices = {forward<size_t>(args)...};
        indices.push_front(t1);
        string key;
        auto it = indices.begin();

        for (size_t i= 0; i<indices.size(); i++) {
            key += to_string(*it);
            if (i<indices.size()-1) {
                key += ",";
            }
            it++;
        }

        auto it2 = param_::_sdpindices->find(key);
        if (it2 == param_::_sdpindices->end()) {
            //not defined before.  
            auto temp = make_pair<>(t1, (*(++indices.begin())));
            res._sdpindices->insert(make_pair<>(key,temp));
            param_::_sdpindices->insert(make_pair<>(key,temp));
        }
        else {
            auto temp = param_::_sdpindices->at(key);
            res._sdpindices->insert(make_pair<>(key,temp));
            res._dim = 1;
            res._symdim = 1;
        }
        res._name += "["+key+"]";
        res._is_indexed = true;
        return res;
    }

    /* Modifiers */
    //void    set_size(size_t s, type val = 0);

    /* Operators */
    sdpvar& operator=(type v){
        param<type>::_val->push_back(v);
        param<type>::_dim++;
        return *this;
    }

    bool operator==(const sdpvar& v) const;
    bool operator>=(const sdpvar& v) const; 
    bool operator!=(const sdpvar& v) const;
    sdpvar& operator^(size_t d){
        // the upper/lower triangular part. 
        param<type>::set_size(d*(d+1)/2);
        _symdim = d;
        return *this;
    }

    sdpvar tr() const{
        auto v = sdpvar(*this);
        v._is_transposed = true;
        return v;
    }
    /* Output */
    void print() const;
};
#endif /* sdpvar_h */
