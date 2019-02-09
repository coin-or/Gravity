//
//  model.hpp
//  Gravity
//
//  Created by Hijazi, Hassan
//
//

#ifndef model_hpp
#define model_hpp

#include <stdio.h>
#include <gravity/constraint.h>
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
    
    template<typename type = double>
    class Model {
        
    protected:
        string                                  _name; /**< Model name. */
        set<pair<size_t,size_t>>                _hess; /**< Pairs of variables linked in the hessian, storing Ipopt indices here. */
        deque<shared_ptr<func<type>>>           _nl_funcs; /**< Queue of all the nonlinear functions appearing in the model. */
        map<string,shared_ptr<func<type>>>      _nl_funcs_map;/**< Map of all the nonlinear functions appearing in the model. */
        
        void add_var(shared_ptr<param_> v);        /**< Add variables without reallocating memory */
        void add_param(shared_ptr<param_> v);      /**< Add variables without reallocating memory */
        
    public:
        bool                            _has_lazy = false; /**< Has lazy constraints. */
        bool                            _built = false; /**< Indicates if this model has been already built. */
        bool                            _first_run = true; /**< Indicates if a solver was ran on this model. */
        
        bool                            _first_call_gard_obj = true; /**< Indicates if this is the first call to fill_in_grad_obj */
        bool                            _first_call_jac = true; /**< Indicates if this is the first call to fill_in_jac */
        bool                            _first_call_hess = true; /**< Indicates if this is the first call to fill_in_hess */
        Convexity                       _convexity = linear_; /**< Indicates the convexity type of the current model */
        MType                           _type = lin_m; /**< Model type, e.g., linar, quadratic, polynomial, NLP.. */
        size_t                          _nb_vars = 0; /**< Number of variables. */
        size_t                          _nb_params = 0; /**< Number of parameters. */
        size_t                          _nb_cons = 0; /**< Number of constraints. */
        size_t                          _nnz_g = 0; /**< Number of non zeros in the Jacobian */
        size_t                          _nnz_h = 0; /**< Number of non zeros in the Hessian */
        size_t                          _nnz_g_obj = 0; /**< Number of non zeros in the Objective gradient */
        
        /* Ipopt data structures */
        vector<double>                  _jac_vals; /* Jacobian values stored in sparse format */
        vector<double>                  _obj_grad_vals; /* Objective gradient values stored in sparse format */
        vector<double>                  _hess_vals; /* Hessian values stored in sparse format */
        /* */
        map<size_t, shared_ptr<param_>>                     _params; /**< Sorted map pointing to all parameters contained in this model. */
        map<size_t, shared_ptr<param_>>                     _vars; /**< Sorted map pointing to all variables contained in this model. */
        map<size_t, shared_ptr<var<bool>>>                  _bin_vars; /**< Sorted map pointing to all binary variables contained in this model. */
        map<string, shared_ptr<param_>>                     _params_name; /**< Sorted map (by name) pointing to all parameters contained in this model. */
        map<string, shared_ptr<param_>>                     _vars_name; /**< Sorted map (by name) pointing to all variables contained in this model. */
        vector<shared_ptr<Constraint<type>>>                _cons_vec; /**< vector pointing to all constraints contained in this model. */
        map<size_t, shared_ptr<Constraint<type>>>           _cons; /**< Sorted map (increasing index) pointing to all constraints contained in this model. */
        map<string, shared_ptr<Constraint<type>>>           _cons_name; /**< Sorted map (by name) pointing to all constraints contained in this model. */
        map<unique_id, set<shared_ptr<Constraint<type>>>>   _v_in_cons; /**< Set of constraints where each variable appears. */
        shared_ptr<func_>                                   _obj = nullptr; /**< Pointer to objective function */
        ObjectiveType                                       _objt = minimize; /**< Minimize or maximize */
        int                                                 _status = -1;/**< status when last solved */
        map<pair<string, string>,set<pair<shared_ptr<func<type>>,shared_ptr<func<type>>>>>            _hess_link; /* for each pair of variables appearing in the hessian, storing the set of constraints they appear together in */

        /** Constructor */
        //@{
        Model();
        Model(const string& name){
            _name = name;
        };
        //@}
        
        template<typename T>
        void update_convexity(const func<T>& f);
        
        template<typename T>
        void update_convexity(const Constraint<T>& f);
        
        /* Accessors */
        
        string get_name() const {return _name;}
        
        size_t get_nb_vars() const;
        
        size_t get_nb_cons() const;
        
        size_t get_nb_ineq() const;
        
        size_t get_nb_nnz_g();
        
        size_t get_nb_nnz_h();
        
        bool is_convex() const{
            return _convexity==convex_;
        }
        
        bool is_concave() const{
            return _convexity==concave_;
        }
        
        shared_ptr<param_> get_var_ptr(const string& vname) const;
        
        shared_ptr<param_> get_var_ptr(size_t idx) const;
        
        
        shared_ptr<Constraint<type>> get_constraint(const string& name) const;
        
        
        bool has_var(const param_& v) const{
            return (_vars.count(v.get_vec_id())!=0);
        };

        
        
        /* Modifiers */
        
        void set_x(const double* x); // Assign values to all variables based on array x.
        void compute_funcs();
        
        template <typename T>
        void add_var(var<T>& v){//Add variables by copy
            if (_vars_name.count(v._name)==0) {
                v.set_id(_nb_vars);
                v.set_vec_id(_vars.size());
                shared_ptr<param_> newv;
                if (v._val->empty()) {
                    Warning("WARNING adding unindexed variable to model, treating it as a one dimensional Real.\n");
                    newv = (v.in(R(1))).copy();
                }
                else {
                    newv = v.copy();
                }
                _vars_name[v._name] = newv;
                _vars[v.get_vec_id()] = newv;
                _nb_vars += v.get_dim();
            }
        };
        
        template <typename T>
        void add(var<T>& v){//Add variables by copy
            add_var(v);
        }
        
        template <typename T, typename... Args>
        void add(var<T>&& v, Args&&... args){
            list<var<T>> vars;
            vars = {forward<var<T>>(args)...};
            vars.push_front(move(v));
            for (auto &v:vars) {
                add_var(move(v));
            }
        }
        
        
        void del_var(const param_& v);
        
        
        void add_integrality_cuts();
        void add_lazy(const Constraint<type>& c);
        void add(const Constraint<type>& c);
        void project();/**<  Use the equations where at least one variable appear linearly to express it as a function of other variables in the problem */
        void replace(shared_ptr<param_> v, func<type>& f);/**<  Replace v with function f everywhere it appears */
        shared_ptr<Constraint<type>> add_constraint(const Constraint<type>& c);
        
        
        shared_ptr<Model> build_McCormick();/**< Build the sequential McCormick relaxation for polynomial programs */
        
        shared_ptr<func<type>> embed(shared_ptr<func<type>> f);/**<  Transfer all variables and parameters to the model. */
        void embed(expr& e);/**<  Transfer all variables and parameters to the model. */
        void del_constraint(const Constraint<type>& c);
        template<typename T>
        void min(const func<T>& f);
        template<typename T>
        void max(const func<T>& f);
        void set_objective(const func<type>& f, ObjectiveType t = minimize);
        void set_objective(pair<shared_ptr<func_>, ObjectiveType> p);
        void set_objective_type(ObjectiveType);
        void init_indices();// Initialize the indices of all variables involved in the model
        void reindex(); /*<< Reindexes the constraints after violated ones have been detected and aded to the formulation */
        bool has_violated_constraints(type tol); /*<< Returns true if some constraints are violated by the current solution with tolerance tol */
        bool is_feasible(type tol);
        void reset_funcs();
        void fill_in_maps();/*< Fill the _hess and _v_in_ maps to link variables with their constraints and compute the Jacobian & Hessian matrices */
        /** IPOPT functions */
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
        void add_on_off(const Constraint<type>& c, var<bool>& on);
        void add_on_off(var<>& v, var<bool>& on);
        
        void add_on_off_McCormick(string name, var<>& v, var<>& v1, var<>& v2, var<bool>& on);
        
        void add_McCormick(string name, shared_ptr<param_> vlift, shared_ptr<param_> v1, shared_ptr<param_> v2);
        
        void replace_integers();/*< Replace internal type of integer variables so that continuous relaxations can be computed */

        
        /* Evaluation Operators */
        
        type eval(shared_ptr<func_> f, size_t i){
            return func<type>::eval(f,i);
        }
        
        type eval(shared_ptr<func_> f, size_t i, size_t j){
            return func<type>::eval(f,i,j);
        }
        
        /* Output */
        void print_nl_functions() const;
        void print_symbolic();
        void print();
        void print_solution() const;
        void round_solution();
        void add_round_solution_cuts();
        void add_round_solution_obj(bool balance_obj = true);
        void add_integrality();
        void print_constraints() const;
        
    };
    
//    void compute_constrs(vector<Constraint*>& v, double* res, unsigned i, unsigned j);

    template<typename T>
    pair<shared_ptr<func_>, ObjectiveType> max(const func<T>& f);
    template<typename T>
    pair<shared_ptr<func_>, ObjectiveType> min(const func<T>& f);

    class Program{
    public:
        virtual void update_model(){};
        virtual ~Program(){};
        string _status;
    };
    
    
}



#endif /* model_hpp */
