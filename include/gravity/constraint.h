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
        vector<double>              _dual ; /**< Lagrange multipliers at a KKT point */
        
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
        Constraint& operator=(const Constraint& c);
        Constraint& operator=(Constraint&& c);
        
        Constraint& operator <=(double rhs);
        Constraint& operator >=(double rhs);        
        Constraint& operator ==(double rhs);
        Constraint& operator <=(const func_& rhs);
        Constraint& operator >=(const func_& rhs);
        Constraint& operator ==(const func_& rhs);
        Constraint& operator =(const func_& rhs);
        
        /* Accessors */
        string get_name() const;
        int get_type() const;
        double get_rhs() const;
        bool is_active(unsigned inst = 0) const;
        bool is_convex() const;
        bool is_concave() const;
        
        size_t get_id_inst(size_t ind) const;
        
        
        /* Modifiers */
        
        Constraint& in(const node_pairs& np){
            this->func_::in(np);
            return *this;
        };
        
        template<typename Tobj> Constraint& in(const vector<Tobj*>& vec){
            this->func_::in(vec);
            return *this;
        };
        
        template<typename Tobj> Constraint& in(const vector<Tobj>& vec){
            this->func_::in(vec);
            return *this;
        };
        
        /* Output */
        void print_expanded();
        void print(unsigned);
        void print();
        
        
    };
}
#endif /* constraint_hpp */
