//
//  model.hpp
//  Gravity
//
//  Created by Hijazi, Hassan (Data61, Canberra City) on 6/5/17.
//
//

#ifndef model_hpp
#define model_hpp

#include <stdio.h>
#include <gravity/constraint.h>
#include <gravity/utils.h>
#include <map>
#include <unordered_set>
#include <math.h>
#include <vector>
#include <deque>
#include <thread>
#ifdef USE_IPOPT
#include <coin/IpIpoptApplication.hpp>
#include <coin/IpTNLP.hpp>
#endif
#ifdef USE_GUROBI
#include <gurobi_c++.h>
#endif
#ifdef USE_BONMIN
#include <coin/BonTMINLP.hpp>
#endif

using namespace std;

namespace gravity {
    
    class Model {
        
    protected:
        string                          _name;
        set<pair<size_t,size_t>>        _hess; /* Pairs of variables linked in the hessian, storing Ipopt indices here. */
        deque<shared_ptr<func_>>        _nl_funcs;
        map<string,shared_ptr<func_>>   _nl_funcs_map;
        void add_var(param_* v);        //Add variables without reallocating memory
        void add_param(param_* v);      //Add variables without reallocating memory
        
    public:
        bool                            _has_lazy = false; /*<< Has lazy constraints. */
        bool                            _built = false; /* Indicates if this model has been already built */
        bool                            _first_run = true; /* Indicates if this model has been already */
                                                            
        bool                            _first_call_gard_obj = true; /* Indicates if this is the first call to fill_in_grad_obj */
        bool                            _first_call_jac = true; /* Indicates if this is the first call to fill_in_jac */
        bool                            _first_call_hess = true; /* Indicates if this is the first call to fill_in_hess */
        bool                            _is_convex = true; /* Indicates if this model is convex */
        MType                           _type = lin_m; /* Model type, e.g., linar, quadratic, polynomial, NLP.. */
        size_t                          _nb_vars = 0;
        size_t                          _nb_params = 0;
        size_t                          _nb_cons = 0;
        size_t                          _nnz_g = 0; /* Number of non zeros in the Jacobian */
        size_t                          _nnz_h = 0; /* Number of non zeros in the Hessian */
        size_t                          _nnz_g_obj = 0; /* Number of non zeros in the Objective gradient */
        vector<double>                  _jac_vals; /* Jacobian values stored in sparse format */        
        vector<double>                  _obj_grad_vals; /* Objective gradient values stored in sparse format */
        vector<double>                  _hess_vals; /* Hessian values stored in sparse format */
        map<unsigned, param_*>          _params; /**< Sorted map pointing to all parameters contained in this model */
        map<unsigned, param_*>          _vars; /**< Sorted map pointing to all variables contained in this model. Note that a variable is a parameter with a bounds attribute. */
        map<string, param_*>            _params_name; /**< Sorted map pointing to all parameters contained in this model */
        map<string, param_*>            _vars_name; /**< Sorted map pointing to all variables contained in this model. Note that a variable is a parameter with a bounds attribute. */
        vector<shared_ptr<Constraint>>              _cons_vec;
        map<unsigned, shared_ptr<Constraint>>       _cons; /**< Sorted map (increasing index) pointing to all constraints contained in this model */
        map<string, shared_ptr<Constraint>>         _cons_name; /**< Sorted map (increasing index) pointing to all constraints contained in this model */
        map<unique_id, set<Constraint*>>  _v_in_cons; /**< Set of constraints where each variable appears */

        /**< Set of variables linked to one another in the hessian, indexed by pairs
         * of variable ids, a pair contains the metavar id and the instanciated var
         * id. The last set contains pairs of indices of constraints where both
         * variables appear, the pair stores the metaconstraint id and the
         * instanciated constraint id 
         * */

        map<pair<string, string>,set<pair<func_*,func_*>>>            _hess_link; /* for each pair of variables appearing in the hessian, storing the set of constraints they appear together */

        func_                           _obj; /** Objective function */
        ObjectiveType                   _objt; /** Minimize or maximize */
        double                          _obj_val = 0;/** Objective function value */

        /** Constructor */
        //@{
        Model();
        Model(const string& name){
            _name = name;
        };
        //@}
        
        /* Destructor */
        ~Model();
        
        /* Accessors */
        
        size_t get_nb_vars() const;
        
        size_t get_nb_cons() const;
        
        size_t get_nb_ineq() const;
        
        size_t get_nb_nnz_g();
        
        size_t get_nb_nnz_h();
        
        bool is_convex() const{
            return _is_convex;
        }
        
        bool is_nonconvex() const{
            return !_is_convex;
        }
        
        param_* get_var(const string& vname) const;
        
        shared_ptr<Constraint> get_constraint(const string& name) const;
        
        void add_indices(const string& constr_name, const node_pairs& indices);
        
        bool has_var(const param_& v) const{
            return (_vars.count(v.get_vec_id())!=0);
        };

        
        
        /* Modifiers */
        
        void set_x(const double* x); // Assign values to all variables based on array x.
        void compute_funcs();
        
        template <typename type>
        void add_var(var<type>& v){//Add variables by copy
            
            if (v._is_indexed) {
                auto nb_ind = v.get_nb_instances();
                for (unsigned i = 0; i<nb_ind; i++) {
                    v._lb->eval(i);
                    v._ub->eval(i);
                }
            }
            if (_vars_name.count(v._name)==0) {
                v.set_id(_nb_vars);
                v.set_vec_id(_vars.size());
                param_* newv;
                if (v._dim[0]==0) {
                    newv = (param_*)copy(v.in(R(1)));
                }
                else {
                    newv = (param_*)copy(v);
                }
                _vars_name[v._name] = newv;
                _vars[v.get_vec_id()] = newv;        
                _nb_vars += v.get_dim();        
            }
        };
        
        
        template <typename type>
        void add_var(var<type>&& v){//Add variables by copy
            if (v._is_indexed) {
                auto nb_ind = v.get_nb_instances();
                for (unsigned i = 0; i<nb_ind; i++) {
                    v._lb->eval(i);
                    v._ub->eval(i);
                }
//                shared_ptr<vector<type>> lb = make_shared<vector<type>>(nb_ind);
//                shared_ptr<vector<type>> ub = make_shared<vector<type>>(nb_ind);
//                unsigned i = 0;
////                for (auto &p: *v.get_rev_indices()) {//TODO fix this
////                    lb->at(i) = v._lb->eval(v._lb->get_indices()->at(p));
////                    ub->at(i) = v._ub->eval(v._ub->get_indices()->at(p));
////                    i++;
////                }
//                v._lb->_val = move(lb);
//                v._ub->_val = move(ub);
//                v._lb->_evaluated = true;
//                v._ub->_evaluated = true;
//                v._lb->_nb_instances = nb_ind;
//                v._lb->_dim[0] = nb_ind;
//                v._ub->_nb_instances = nb_ind;
//                v._ub->_dim[0] = nb_ind;                
            }
            if (_vars_name.count(v._name)==0) {
                v.set_id(_nb_vars);
                v.set_vec_id(_vars.size());
                param_* newv;
                if (v._dim[0]==0) {
                    newv = (param_*)copy(v.in(R(1)));
                }
                else {
                    newv = (param_*)copy(v);
                }
                _vars_name[v._name] = newv;
                _vars[v.get_vec_id()] = newv;
                _nb_vars += v.get_dim();
            }
        };
        
        void del_var(const param_& v);
        
        
        void add_param(param_& v); //Add variables by copying variable
        void del_param(const param_& v);
        
        
        void add_lazy(const Constraint& c);
        void add(const Constraint& c);
        shared_ptr<Constraint> add_constraint(const Constraint& c);
        
        shared_ptr<func_> embed(shared_ptr<func_> f);/**<  Transfer all variables and parameters to the model, useful for a centralized memory management. */
        void embed(expr& e);/**<  Transfer all variables and parameters to the model, useful for a centralized memory management. */
        void del_constraint(const Constraint& c);
        void min(const func_& f);
        void max(const func_& f);
        void set_objective(const func_& f, ObjectiveType t = minimize);
        void set_objective(pair<func_*, ObjectiveType> p);
        void set_objective_type(ObjectiveType);
        void init_indices();// Initialize the indices of all variables involved in the model
        void reindex(); /*<< Reindexes the constraints after violated ones have been detected and aded to the formulation */
        bool has_violated_constraints(double tol); /*<< Returns true if some constraints are violated by the current solution with tolerance tol */
        void check_feasible(const double* x);
        void reset_funcs();
        void fill_in_maps();/*< Fill the _hess and _v_in_ maps to link variables with their constraints and compute the Jacobian & Hessian matrices */
        void fill_in_var_bounds(double* x_l ,double* x_u);
        void fill_in_var_init(double* x);
        void fill_in_duals(double* lambda,double* z_l, double* z_u);
        void fill_in_cstr_bounds(double* g_l ,double* g_u);
        void fill_in_obj(const double* x , double& res, bool new_x);
        void fill_in_grad_obj(const double* x , double* res, bool new_x);
        void fill_in_cstr(const double* x , double* res, bool new_x);
        void fill_in_jac(const double* x , double* res, bool new_x);
        void fill_in_jac_nnz(int* iRow , int* jCol);
        void fill_in_hess(const double* x , double obj_factor, const double* lambda, double* res, bool new_x);
        void eval_funcs_parallel(const double* x , int start, int end);
        void fill_in_hess_nnz(int* iRow , int* jCol);
    #ifdef USE_IPOPT
        void fill_in_var_linearity(Ipopt::TNLP::LinearityType* param_types);
        void fill_in_cstr_linearity(Ipopt::TNLP::LinearityType* const_types);
    #endif
        
    #ifdef USE_BONMIN
        void fill_in_var_types(Bonmin::TMINLP::VariableType* var_types);
    #endif        
        void add_on_off(const Constraint& c, var<bool>& on);
        void add_on_off(var<>& v, var<bool>& on);
        
        void add_on_off_McCormick(string name, var<>& v, var<>& v1, var<>& v2, var<bool>& on);
        void add_McCormick(string name, var<>& v, var<>& v1, var<>& v2);

        
        
        /* Operators */
        
        
        
        /* Output */
        void print_nl_functions() const;
        void print() const;
        void print_expanded() const;
        void print_solution() const;
        void print_constraints() const;
        
    };
    
//    void compute_constrs(vector<Constraint*>& v, double* res, unsigned i, unsigned j);

    pair<func_*, ObjectiveType> max(const func_& f);
    pair<func_*, ObjectiveType> min(const func_& f);

    class Program{
    public:
        virtual void update_model(){};
        virtual ~Program(){};
    };
}
#endif /* model_hpp */
