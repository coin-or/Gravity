//
//  var.h
//  Gravity
//
//  Created by Hijazi, Hassan on 21/10/16.
//
//

#ifndef var_h
#define var_h

#include <Gravity/param.h>
#include <stdio.h>
#include <string>
#include <set>
#include <list>

//class func;

using namespace std;

/** Backbone class for parameter */
class var_ {
public:
    virtual ~var_(){};
    
};

/** A variable can be a bool, a short, an int, a float or a double*/
template<typename type = int>

// define variable as a parameter with bounds
class var: public param<type>, public var_{
    
public:
    ind                         _id;
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
        auto res(*this);
//        res._indices->clear();
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
        auto it2 = param_::_indices->find(key);
        if (it2 == param_::_indices->end()) {
            res._indices->insert(make_pair<>(key,param_::_indices->size()));
            param_::_indices->insert(make_pair<>(key,param_::_indices->size()));
        }
        else {
            size_t idx = param_::_indices->at(key);
            res._indices->insert(make_pair<>(key,idx));
        }
        return res;
    }
    
//    var operator()(int idx...){
//        auto res(*this);
//        res._indices.insert(idx);
//        va_list arg_list;
//        va_start(arg_list, idx);
//        int i = 0;
//        while (*idx != '\0') {
//            if (*idx == 'd') {
//                i = va_arg(arg_list, int);
//                res._indices.insert(i);
//            }
//            else {
//                throw invalid_argument("indices can only be integers");
//            }
//            ++idx;
//        }
//        va_end(arg_list);
//    }
    
    
    /* Querries */
    
    type    get_lb(size_t i) const{
        if (_lb->size() <= i) {
            throw out_of_range("get_lb(size_t i, index: " + to_string(i) + ")\n");
        }
        return _lb->at(i);
    };
    
    type    get_ub(size_t i) const{
        if (_ub->size() <= i) {
            throw out_of_range("get_ub(size_t i), index: " + to_string(i)+ ")\n");
        }
        return _lb->at(i);
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
    bool operator==(const var& v) const;
    bool operator!=(const var& v) const;
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

#endif /* var_h */
