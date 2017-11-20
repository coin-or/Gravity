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
#include <gravity/func.h>

namespace gravity {
    class Constraint :public func_{
        
    protected:
        string                      _name = "no_name";
        
    public:        
        
        unsigned                    _id = -1;
        ConstraintType              _ctype = leq; /**< Constraint type: leq, geq or eq */
        double                      _rhs = 0;
        double                      _dual = 0; /**< Lagrange multipliers at a KKT point */
        
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
        Constraint& operator <=(double rhs);
        Constraint& operator >=(double rhs);
        Constraint& operator =(double rhs);
        Constraint& operator =(const func_& f);
        
        /* Accessors */
        string get_name() const;
        int get_type() const;
        double get_rhs() const;
        bool is_active() const;
        size_t get_id_inst(size_t ind) const;
        
        
        /* Modifiers */
        
        template<typename Tobj> Constraint& in(const vector<Tobj*>& vec){
            this->func_::in(vec);
            return *this;
        };
        
        /* Output */
        void print() const;
        
        
    };
}
#endif /* constraint_hpp */
