//
//  constraint.hpp
//  Gravity
//
//  Created by Hijazi, Hassan (Data61, Canberra City) on 6/5/17.
//
//

#ifndef constraint_hpp
#define constraint_hpp

#include <stdio.h>
#include <Gravity/func.h>

class Constraint :public func_{
    
protected:
    string                      _name = "no_name";
    
public:
    
    ConstraintType              _ctype; /**< Constraint type: leq, geq or eq */
    double _rhs;
    /** Constructor */
    //@{
    Constraint();
    Constraint(const Constraint& c);
    Constraint(std::string name);
    Constraint(std::string name, ConstraintType ctype);
    //@}
    
    
    
    /* Destructor */
    ~Constraint();
    
    
    /* Boolean Requests */
    
    
    /* Operators */
    Constraint& operator<=(double rhs);
    Constraint& operator>=(double rhs);
    Constraint& operator=(double rhs);
    Constraint& operator=(const func_& f);
    //
    //    Constraint& operator<=(int rhs);
    //    Constraint& operator>=(int rhs);
    //    Constraint& operator=(int rhs);
    
    /* Accessors */
    string get_name() const;
    int get_type() const;
    double get_rhs() const;
    
    
    
    /* Modifiers */
    
    /* Output */
    void print() const;
    
    
};

#endif /* constraint_hpp */
