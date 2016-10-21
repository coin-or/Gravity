//
//  var.hpp
//  Gravity
//
//  Created by Hijazi, Hassan (Data61, Canberra City) on 21/10/16.
//
//

#ifndef var_hpp
#define var_hpp

#include <Gravity/param.h>
#include <stdio.h>

/** A variable can be a bool, a short, an int, a float or a double*/
template<typename type = int>
class var: public param<type>{
    
public:
    shared_ptr<vector<type>>                _lb; /**< Lower Bound */
    shared_ptr<vector<type>>                _ub; /**< Upper Bound */
    
    
    /* Constructors */
    
    //@{
    /** Unbounded variable constructor */
    var();
    var(string name);
    var(const var<type>& v);
    //@}
    
    
    //@{
    /** Bounded variable constructor */
    var(string name, type lb, type ub);
    var(string name, type lb, type ub, type lb_off, type ub_off);
    //@}
    
    /* Querries */
//    bool is_bounded_above() const{
//        return _bounded_up;
//    };
//    
//    bool is_bounded_below() const{
//        return _bounded_down;
//    };
    
//    bool is_constant() const{
//        return (_type==fixed_);
//    };
    
    
//    bool is_int() const{
//        return (_type==integer_);
//    };
//    
//    bool is_short() const{
//        return (_type==short_);
//    };
//    
//    
//    bool is_binary() const{
//        return (_type==binary_);
//    };
//    
//    bool is_float() const{
//        return (_type==float_);
//    };
//    
//    bool is_double() const{
//        return (_type==double_);
//    };
//    
//    bool is_long() const{
//        return (_type==long_);
//    };

};

#endif /* var_hpp */
