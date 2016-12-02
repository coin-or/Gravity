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
#include <set>


//class func;

using namespace std;

/** Backbone class for parameter */
class var_ {
public:
    virtual ~var_(){};
    
};




/** A variable can be a bool, a short, an int, a float or a double*/
template<typename type = int>
class var: public var_, public param<type>{
    
public:
    ind                         _id;
    shared_ptr<vector<type>>    _lb; /**< Lower Bound */
    shared_ptr<vector<type>>    _ub; /**< Upper Bound */
    
    
    /* Constructors */
    
    //@{
    /** Unbounded variable constructor */
    var();
    var(const char* name);
    var(const var<type>& v);
    //@}
    
    
    //@{
    /** Bounded variable constructor */
    var(const char* name, type lb, type ub);
    //@}
    
    
    /* Querries */
    
    type    get_lb() const;
    type    get_ub() const;

    
    bool is_bounded_above(int i = 0) const{
        return (_ub!=nullptr &&  i < _ub->size() && _ub->at(i)!=numeric_limits<type>::max());
    };

    bool is_bounded_below(int i = 0) const{
        return (_lb!=nullptr &&  i < _lb->size() && _lb->at(i)!=numeric_limits<type>::min());
    };
    
    bool is_constant(int i=0) const{
        return (is_bounded_below() && is_bounded_above() && _lb->at(i)==_ub->at(i));
    };
    
    
    /* Modifiers */
    void    set_size(int s);
    
    void    add_bounds(type lb, type ub);
    void    add_lb_only(type v); /**< Adds a new lower bound and infinity in the corresponding upper bound*/
    void    add_ub_only(type v); /**< Adds a new upper bound and -infinity in the corresponding lower bound*/
    
    void    set_lb(int i, type v);
    void    set_ub(int i, type v);

    
    /* Operators */
    bool operator==(const var& v) const;
    bool operator!=(const var& v) const;
    
    
    /* Output */
    void print(bool bouds=false) const;
    

    
};

#endif /* var_h */
