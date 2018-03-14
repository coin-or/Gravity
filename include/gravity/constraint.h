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
        shared_ptr<vector<string>>  _rev_indices = nullptr; /*<< A vector storing all the indices this parameter has */
        
    public:        
        unsigned                    _jac_cstr_idx; /* Firt index of the corresponding non-zero values in the Jacobian */
        unsigned                    _id = -1;
        ConstraintType              _ctype = leq; /**< Constraint type: leq, geq or eq */
        double                      _rhs = 0;
        vector<double>              _dual ; /**< Lagrange multipliers at a KKT point */
        bool                        _all_active = true;
        vector<bool>                _active;
        shared_ptr<bool>            _all_lazy;
        vector<bool>                _lazy;
        bool                        _all_satisfied = true;
        vector<bool>                _violated;

        
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
        size_t get_nb_instances() const;
        string get_name() const;
        int get_type() const;
        double get_rhs() const;
        bool is_active(unsigned inst = 0) const;
        bool is_convex() const;
        bool is_concave() const;
        bool is_ineq() const;
        
        size_t get_id_inst(size_t ind) const;
        
        
        /* Modifiers */
        
        void make_lazy() {
            *_all_lazy = true;
            _lazy.resize(_nb_instances,true);
        }
        
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

        template<typename Tobj> Constraint& in_at(const vector<Tobj>& vec, const unsigned t){
//            std::vector<string> keys;
//            string key;
//            for(auto it = vec.begin(); it!= vec.end(); it++) {
//                key = (*it)->_name;
//                key += "," + to_string(t);
//                keys.push_back(key);
//            }
            this->func_::in_at(vec, t);
            return *this;
        };

        template<typename Tobj> Constraint& in(const vector<Tobj*>& vec, const unsigned T){
            this->func_::in(vec, T);
            return *this;
        };
        
        template<typename Tobj> Constraint& in_at(const vector<Tobj*>& vec, const unsigned T){
            this->func_::in_at(vec, T);
            return *this;
        };
        
        /* Output */
        void print_expanded();
        void print(unsigned);
        void print();
        
        
    };
}
#endif /* constraint_hpp */
