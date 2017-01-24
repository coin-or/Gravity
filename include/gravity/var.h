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
    var(const char* name);
    var(const var<type>& v);
    //@}
    
    
    //@{
    /** Bounded variable constructor */
    var(const char* name, type lb, type ub);
    //@}
    
    
//    template<class... Args>
//    var operator()(Args&&... args){
//        auto res(*this);
//        //        va_list arg_list;
//        va_start(args, sizeof...(args));
//        for (int i = 0; i < sizeof...(args); ++i) {
//            res._indices.insert(va_arg(args, int));
//        }
//        va_end(args);
//    }
    
    template<class... Args>
    var operator()(Args&&... args){
        auto res(*this);
        *res._indices = {forward<Args>(args)...};
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
    void print(bool bounds=false) const;
    

    
};

#endif /* var_h */
