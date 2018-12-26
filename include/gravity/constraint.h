////
////  constraint.hpp
////  Gravity
////
////  Created by Hijazi, Hassan (Data61, Canberra City) on 6/5/17.
////
////
//
//#ifndef constraint_hpp
//#define constraint_hpp
//
//#include <stdio.h>
//#include <gravity/func.h>
//
//namespace gravity {
//    class Constraint_ :public func_{
//        
//    public:
//        size_t                    _jac_cstr_idx; /* Firt index of the corresponding non-zero values in the Jacobian */
//        size_t                      _id = 0;
//        ConstraintType              _ctype = leq; /**< Constraint type: leq, geq or eq */
//        vector<double>              _dual ; /**< Lagrange multipliers at a KKT point */
//        bool                        _all_active = true;
//        vector<bool>                _active;
//        shared_ptr<bool>            _all_lazy;
//        vector<bool>                _lazy;
//        bool                        _all_satisfied = true;
//        vector<bool>                _violated;
//
//        
//        /** Constructors */
//        //@{
//        Constraint_();
//        Constraint_(const Constraint_& c);
//        Constraint_(std::string name);
//        Constraint_(std::string name, ConstraintType ctype);
//        //@}
//        
//        
//        /* Destructor */
//        ~Constraint_();
//        
//        
//        /* Boolean Requests */
//        
//        /* Operators */
//        Constraint_& operator=(const Constraint_& c);
//        Constraint_& operator=(Constraint_&& c);
//        
//        Constraint_& operator <=(const func_& rhs);
//        Constraint_& operator >=(const func_& rhs);
//        Constraint_& operator ==(const func_& rhs);
//        Constraint_& operator =(const func_& rhs);
//        
//        /* Accessors */
//        size_t get_nb_instances() const;
//        string get_name() const;
//        int get_type() const;
//        bool is_convex() const;
//        bool is_concave() const;
//        bool is_ineq() const;
//        
//        size_t get_id_inst(size_t ind) const;
//        
//        
//        /* Modifiers */
//        
//        void make_lazy() {
//            *_all_lazy = true;
//            _lazy.resize(get_dim(),true);
//        }
//        
////
////
////        Constraint_& in(const node_pairs& np){
////            this->func_::in(np);
////            return *this;
////        };
////
//        Constraint_& in(const vector<Node*>& vec) {
//            this->func_::in(vec);
//            return *this;
//        }
//        
//        Constraint_& in(const indices& ids){
//            if(ids.empty()){
//                _dim[0] = 0;
//                return *this;
//            }
//            this->func_::in(ids);
//            return *this;
//        };
////
//        
//        /* Output */
//        void print();
//        void print(size_t);
//        void print_symbolic();
//        
//        
//    };
//    
//    template<typename type = double>
//    class Constraint: public Constraint_, public func<type>{
//        void print(){
//            auto nb_inst = _dim[0];
//            allocate_mem();
//            for (size_t inst = 0; inst<nb_inst; inst++) {
//                if (*_all_lazy && _lazy[inst]) {
//                    continue;
//                }
//                print(inst);
//            }
//        }
//        bool is_active(size_t inst = 0, double tol = 1e-6) const{
//            return fabs(var<type>::_val->at(inst)) < tol;
//        }
//    };
//}
//#endif /* constraint_hpp */
