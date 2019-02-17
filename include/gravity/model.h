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
    
    /**
     Parallel computation of the constraints stored in v[i] to v[j]
     @param[in] v vector of constraints
     @param[in] res vector storing the computation results
     @param[in] i starting index in v
     @param[in] j ending index in v
     */
    template<typename type>
    void compute_constrs(vector<shared_ptr<Constraint<type>>>& v, double* res, size_t i, size_t j){
        DebugOff("Calling compute_constrts with i =  " << i << "and j = "<< j << endl);
        for (size_t idx = i; idx < j; idx++) {
            auto c = v[idx];
            size_t nb_ins = c->_dim[0];
            size_t ind = 0;
            for (size_t inst = 0; inst< nb_ins; inst++){
                if (!*c->_all_lazy || !c->_lazy[inst]) {
                    res[c->_id+ind++] = c->eval(inst);
                    DebugOff("Accessing res at position " << c->_id+inst << endl);
                    //                _cons_vals[index++] = res[c->_id+inst];
                    DebugOff("g[" << to_string(c->_id+inst) << "] = " << to_string(res[c->_id+inst]) << endl);
                }
            }
        }
    }
    
    /**
     Parallel computation of the jacobian of the constraints stored in v[i] to v[j]
     @param[in] v vector of constraints
     @param[in] res vector storing the computation results
     @param[in] i starting index in v
     @param[in] j ending index in v
     @param[in] first_call true if first_call to this function, false otherwise
     @param[in] jac_vals Gravity's internal vector for storing the results
     */
    template<typename type>
    void compute_jac(vector<shared_ptr<Constraint<type>>>& vec, double* res, size_t i, size_t j, bool first_call, vector<double>& jac_vals){
        size_t cid = 0, id_inst = 0;
        string vid;
        shared_ptr<Constraint<type>> c = NULL;
        param_* v = NULL;
        shared_ptr<func<type>> dfdx;
        auto idx = vec[i]->_jac_cstr_idx;
        for (size_t ind = i; ind < j; ind++) {
            c = vec[ind];
            auto nb_ins = c->_dim[0];
            id_inst = 0;
            if (c->is_linear() && !first_call) {
                //        if (false) {
                DebugOff("Linear constraint, using stored jacobian!\n");
                for (size_t i = 0; i<nb_ins; i++) {
                    if (!*c->_all_lazy || !c->_lazy[i]) {
                        for (size_t j = 0; j<c->get_nb_vars(i); j++) {
                            res[idx] = jac_vals[idx];
                            idx++;
                        }
                    }
                }
            }
            else {
                for (auto &v_p: c->get_vars()){
                    v = v_p.second.first.get();
                    vid = v->_name;
                    dfdx = c->get_stored_derivative(vid);
                    id_inst = 0;
                    for (size_t inst = 0; inst< nb_ins; inst++){
                        if (!*c->_all_lazy || !c->_lazy[inst]) {
                            cid = c->_id+id_inst++;
                            if (v->_is_vector) {
                                auto dim = v->get_dim(inst);
                                for (size_t j = 0; j<dim; j++) {
                                    res[idx] = dfdx->eval(inst,j);//TODO: res[idx] += .. (account for vectors with repeated identical entries)
                                    jac_vals[idx] = res[idx];
                                    DebugOff("jac_val["<< idx <<"] = " << jac_vals[idx] << endl);
                                    idx++;
                                }
                            }
                            else {
                                res[idx] = dfdx->eval(inst);
                                jac_vals[idx] = res[idx];
                                idx++;
                            }
                        }
                    }
                }
            }
        }
    }
    
    template<typename type = double>
    class Model {
        
    protected:
        string                                  _name; /**< Model name. */
        set<pair<size_t,size_t>>                _hess; /**< Pairs of variables linked in the hessian, storing Ipopt indices here. */
        deque<shared_ptr<func<type>>>           _nl_funcs; /**< Queue of all the nonlinear functions appearing in the model. */
        map<string,shared_ptr<func<type>>>      _nl_funcs_map;/**< Map of all the nonlinear functions appearing in the model. */
        
        
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
        map<size_t, shared_ptr<param_>>                     _int_vars; /**< Sorted map pointing to all binary variables contained in this model. */
        map<string, shared_ptr<param_>>                     _params_name; /**< Sorted map (by name) pointing to all parameters contained in this model. */
        map<string, shared_ptr<param_>>                     _vars_name; /**< Sorted map (by name) pointing to all variables contained in this model. */
        vector<shared_ptr<Constraint<type>>>                _cons_vec; /**< vector pointing to all constraints contained in this model. */
        map<size_t, shared_ptr<Constraint<type>>>           _cons; /**< Sorted map (increasing index) pointing to all constraints contained in this model. */
        map<string, shared_ptr<Constraint<type>>>           _cons_name; /**< Sorted map (by name) pointing to all constraints contained in this model. */
        map<unique_id, set<shared_ptr<Constraint<type>>>>   _v_in_cons; /**< Set of constraints where each variable appears. */
        shared_ptr<func<type>>                              _obj = nullptr; /**< Pointer to objective function */
        ObjectiveType                                       _objt = minimize; /**< Minimize or maximize */
        int                                                 _status = -1;/**< status when last solved */
        map<pair<string, string>,set<pair<shared_ptr<func<type>>,shared_ptr<func<type>>>>>            _hess_link; /* for each pair of variables appearing in the hessian, storing the set of constraints they appear together in */

        /** Constructor */
        //@{
        Model(){
            _obj = make_shared<func<type>>();
        };
        Model(const string& name){
            _name = name;
            _obj = make_shared<func<type>>();
        };
        //@}

        
        /* Accessors */
        
        string get_name() const {return _name;}
        
        bool is_linear() const{
            return _type==lin_m;
        }
        
        bool is_convex() const{
            return _convexity==convex_;
        }
        
        bool is_concave() const{
            return _convexity==concave_;
        }
        
        bool has_var(const param_& v) const{
            return (_vars.count(v.get_vec_id())!=0);
        };

        
        
        /* Modifiers */
        
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
        void add_var(var<T>&& v){//Add variables by copy
            if (_vars_name.count(v._name)==0) {
                v.set_id(_nb_vars);
                v.set_vec_id(_vars.size());
                shared_ptr<param_> newv;
                if (v._val->empty()) {
                    Warning("WARNING adding unindexed variable to model, treating it as a one dimensional Real.\n");
                    newv = make_shared<var<type>>(move((v.in(R(1)))));
                }
                else {
                    newv = make_shared<var<type>>(move(v));
                }
                _vars_name[v._name] = newv;
                _vars[v.get_vec_id()] = newv;
                _nb_vars += newv->get_dim();
            }
        };
        
        template <typename T>
        void add(var<T>& v){//Add variables by copy
            add_var(v);
        }
        
        template <typename T, typename... Args>
        void add(var<T>&& v, Args&&... args){
            list<var<type>> vars;
            vars = {forward<var<type>>(args)...};
            vars.push_front(move(v));
            for (auto &v:vars) {
                add_var(move(v));
            }
        }
        
        
        /* Output */
        
        
        /* Accessors */
        
        size_t get_nb_vars() const{
            return _nb_vars;
        };
        
        
        size_t get_nb_cons() const{
            size_t n = 0;
            for (auto &cp:_cons) {
                n += cp.second->get_nb_instances();
            }
            return n;
        };
        
        
        size_t get_nb_ineq() const{
            size_t n = 0;
            for (auto &cp:_cons) {
                if (cp.second->is_ineq()) {
                    n += cp.second->get_nb_instances();
                }
            }
            return n;
        };
        
        size_t get_nb_eq() const{
            size_t n = 0;
            for (auto &cp:_cons) {
                if (cp.second->is_eq()) {
                    n += cp.second->get_nb_instances();
                }
            }
            return n;
        };
        
        
        size_t get_nb_nnz_g(){
            _nnz_g = 0;
            for (auto &cp:_cons) {
                auto c = cp.second;
                auto nb_inst = c->_dim[0];
                for (size_t inst = 0; inst<nb_inst; inst++) {
                    if (!*c->_all_lazy || !c->_lazy[inst]) {
                        _nnz_g += c->get_nb_vars(inst);
                    }
                }
            }
            return _nnz_g;
        };
        
        
        /* Return the number of nonzeros in the lower left part of the hessian */
        
        size_t get_nb_nnz_h(){
            size_t idx = 0;
            bool idx_inc = false;
            shared_ptr<Constraint<type>> c = nullptr;
            for (auto &pairs: _hess_link) {
                for (auto &f_pair:pairs.second) {
                    idx_inc = false;
                    auto f = f_pair.first;
                    if (f->_is_constraint) {
                        c = static_pointer_cast<Constraint<type>>(f);
                    }
                    auto d2f = f_pair.second;
                    size_t nb_inst = f->_dim[0];
                    for (size_t inst = 0; inst<nb_inst; inst++) {
                        if (!(f->_is_constraint && *c->_all_lazy && c->_lazy[inst])) {
                            if (d2f->is_matrix()) {
                                for (size_t i = 0; i < d2f->_dim[0]; i++) {
                                    for (size_t j = i; j < d2f->_dim[1]; j++) {
                                        idx++;
                                    }
                                }
                            }
                            else if(d2f->_is_vector){
                                for (size_t j = 0; j < d2f->_dim[0]; j++) {
                                    idx++;
                                }
                            }
                            else {
                                idx++;
                            }
                        }
                    }
                }
            }
            _nnz_h = idx;
            return _nnz_h;
        };
        
        
        
        shared_ptr<Constraint<type>> get_constraint(const string& cname) const{
            return _cons_name.at(cname);
        }
        
        
        shared_ptr<param_> get_var_ptr(const string& vname) const{
            auto it = _vars_name.find(vname);
            if (it==_vars_name.end()) {
                return nullptr;
            }
            return it->second;
        }
        
        shared_ptr<param_> get_var_ptr(size_t idx) const{
            return _vars.at(idx);
        }
        
        
        /* Modifiers */
        
        /** Reindexes the constraints after violated ones have been detected and added to the formulation */
        void reindex(){
            size_t cid = 0, new_cid = 0, nb_inst = 0;
            shared_ptr<Constraint<type>> c = nullptr;
            map<size_t, shared_ptr<Constraint<type>>>  new_cons;
            for(auto& c_p: _cons_name)
            {
                c = c_p.second;
                nb_inst = c->get_nb_instances();
                if (nb_inst==0) {
                    continue;
                }
                if (_type==lin_m && c->is_quadratic()) {
                    _type = quad_m;
                }
                else if ((_type==lin_m || _type==quad_m) && c->is_polynomial()) {
                    _type = pol_m;
                }
                else if (c->is_nonlinear()) {
                    _type = nlin_m;
                }
                update_convexity(*c);
                cid = c->_id;
                if (cid!=new_cid) {
                    c->_id = new_cid;
                }
                new_cons[c->_id] = c;
                new_cid = c->_id+nb_inst;
            }
            _cons = new_cons;
            _nb_cons = get_nb_cons();
        }
        
        
        void init_indices(){/**< Initialize the indices of all variables involved in the model */
            shared_ptr<param_> v= nullptr;
            size_t idx = 0, vec_idx = 0;
            for(auto& v_p: _vars)
            {
                v = v_p.second;
                v->set_vec_id(vec_idx++);
                v->set_id(idx);
                for (size_t i = 0; i < v->get_dim(); i++) {
                    idx++;
                }
            }
        }
        
        
        void del_var(const param_& v){
            auto it = _vars.find(v.get_id());
            if (it!=_vars.end()) {
                _nb_vars -= v.get_dim();
                delete it->second;
                _vars.erase(it);
            }
            init_indices();
        };
        
        
        
        void add(const Constraint<type>& c){
            if (c.get_dim()==0) {
                return;
            }
            add_constraint(c);
        }
        
        
        void add_lazy(Constraint<type>& c){
            if (c.get_dim()==0) {
                return;
            }
            c.make_lazy();
            add_constraint(c);
            _has_lazy = true;
        }
        
        template<typename T>
        void replace(shared_ptr<param_> v, func<T>& f){/**<  Replace v with function f everywhere it appears */
            for (auto &c_p: _cons_name) {
                auto c = c_p.second;
                if (!c->has_var(*v)) {
                    continue;
                }
                c->replace(v, f);
            }
            _vars_name.erase(v->_name);
            auto vid = *v->_vec_id;
            delete _vars.at(vid);
            _vars.erase(vid);
            init_indices();
        }
        
        
        void project() {/**<  Use the equations where at least one variable appear linearly to express it as a function of other variables in the problem */
            for (auto& c_pair:_cons_name) {
                if (!c_pair.second->is_ineq()) {
                    auto &lterms = c_pair.second->get_lterms();
                    if (!lterms.empty()) {
                        auto first = lterms.begin();
                        auto v = first->second._p;
                        if (v->_is_vector) {
                            continue;
                        }
                        auto f = *c_pair.second;
                        if (first->second._sign) {
                            //                    f -= *v;
                            //                    f *= -1;
                        }
                        else {
                            //                    f += *v;
                        }
                        DebugOff(f.to_str());
                        _cons.erase(c_pair.second->_id);
                        _cons_name.erase(c_pair.first);
                        replace(v,f);
                        //                project();
                        return;
                    }
                }
            }
        }
        
        
        shared_ptr<Constraint<type>> add_constraint(const Constraint<type>& c){
            if (c.get_dim()==0) {
                return nullptr;
            }
            if (_cons_name.count(c.get_name())==0) {
                auto newc = make_shared<Constraint<type>>(c);
                newc->check_soc();
                newc->check_rotated_soc();
                if (newc->is_constant()) {
                    switch (newc->_ctype) {
                        case leq:
                            if (newc->is_positive()) {
                                throw invalid_argument("Adding violated constant constraint!\n");
                            }
                            break;
                        case geq:
                            if (newc->is_negative()) {
                                throw invalid_argument("Adding violated constant constraint!\n");
                            }
                            break;
                        case eq:
                            if (!newc->is_zero()) {
                                throw invalid_argument("Adding violated constant constraint!\n");
                            }
                            break;
                        default:
                            break;
                    }
                    Warning("WARNING: Adding redundant constant constraint, Gravity will be ignoring it.\n");
                    return newc;
                }
                newc->update_str();
                embed(newc);
                update_convexity(*newc);
                newc->_violated.resize(newc->_dim[0],true);
                _cons_name[c.get_name()] = newc;
                if(*newc->_all_lazy){
                    newc->_lazy.resize(newc->_dim[0],true);
                    newc->allocate_mem();
                    return newc;
                }
                size_t nb_inst = c.get_nb_instances();
                if (nb_inst>0) {
                    newc->_id = _nb_cons;
                    _cons[newc->_id] = newc;
                    _cons_vec.push_back(newc);
                    _built = false;
                    _nb_cons += nb_inst;
                    for (size_t inst = 0; inst<nb_inst; inst++) {
                        _nnz_g += c.get_nb_vars(inst);
                    }
                    if (_type==lin_m && c.is_quadratic()) {
                        _type = quad_m;
                    }
                    else if ((_type==lin_m || _type==quad_m) && c.is_polynomial()) {
                        _type = pol_m;
                    }
                    else if(c.is_nonlinear()){
                        _type = nlin_m;
                    }
                    newc->allocate_mem();
                }
                return newc;
            }
            else {
                throw invalid_argument("rename constraint as this name has been used by another one: " + c.get_name());
            }
        };
        
        
        
        void del_constraint(const Constraint<type>& c){
            _cons_name.erase(c.get_name());
            _cons.erase(c._id);
            _built = false;
            _cons_vec.clear();
            _cons_vec.resize(_cons.size());
            size_t i = 0;
            for(auto& c_p: _cons)
            {
                _cons_vec[i++] = c_p.second;
            }
            /* TODO: make sure other infos in model are updated */
        };
        
        /**
         Update the convexity of the current model after adding f (either objective or constraint)
         @param[in] f function/constraint to be added to the model
         */
        template<typename T1>
        void update_convexity(const func<T1>& f){
            if(_convexity==linear_){
                if ((_objt==minimize && f.is_convex()) || (_objt==maximize && f.is_concave())) {
                    _convexity = convex_;
                }
                else if ((_objt==minimize && f.is_concave()) || (_objt==maximize && f.is_convex())) {
                    _convexity = concave_;
                }
            }
            else if(_convexity==convex_){
                if (!((_objt==minimize && f.is_convex()) || (_objt==maximize && f.is_concave()))) {
                    _convexity = undet_;
                }
            }
            else if(_convexity==concave_){
                if (!((_objt==minimize && f.is_concave()) || (_objt==maximize && f.is_convex()))) {
                    _convexity = undet_;
                }
            }
            else {
                _convexity = undet_;
            }
        }
        
        
        template<typename T>
        void set_objective(const func<T>& f, ObjectiveType t) {
            *_obj = f;
            _objt = t;
            update_convexity(f);
            embed(_obj);
        }
        
        template<typename T1>
        void min(const func<T1>& f){
            set_objective(f, minimize);
        }
        
        template<typename T1>
        void max(const func<T1>& f){
            set_objective(f, maximize);
        }
        
        
        void set_objective_type(ObjectiveType t) {
            _objt = t;
        }
        
        
        bool has_violated_constraints(type tol){/*<< Returns true if some constraints are violated by the current solution with tolerance tol */
            //    if (!_has_lazy) {
            //        return false;
            //    }
            //    int cid = 0;
            size_t nb_inst = 0, nb_viol = 0, nb_viol_all = 0;
            size_t nb_active = 0, nb_active_all = 0;
            double diff = 0;
            shared_ptr<Constraint<type>> c = nullptr;
            bool violated = false;
            for(auto& c_p: _cons_name)
            {
                c = c_p.second;
                //        cid = c->_id;
                nb_inst = c->_dim[0];
                nb_viol = 0;
                nb_active = 0;
                c->_all_satisfied = true;
                c->_violated.resize(nb_inst);
                c->_active.resize(nb_inst);
                switch (c->get_ctype()) {
                    case eq:
                        for (size_t inst=0; inst<nb_inst; inst++) {
                            diff = abs(c->eval(inst));
                            if(diff > tol) {
                                DebugOff("Violated equation: ");
                                //                        c->print(inst);
                                DebugOff(", violation = "<< diff << endl);
                                nb_viol++;
                                //                        violated = true;
                                if (*c->_all_lazy) {
                                    c->_all_satisfied = false;
                                    c->_violated[inst] = true;
                                    violated = true;
                                    c->_lazy[inst] = false;
                                }
                                else {
                                    //                            throw runtime_error("Non-lazy constraint is violated, solution declared optimal by solver!\n" + c->to_str(inst));
                                }
                                //                        c->_violated[inst] = true;
                            }
                            else {
                                //                        c->_violated[inst] = false;
                            }
                            //                    nb_active++;
                        }
                        break;
                    case leq:
                        for (size_t inst=0; inst<nb_inst; inst++) {
                            c->_violated[inst] = false;
                            diff = c->eval(inst);
                            if(diff > tol) {
                                DebugOn("Violated inequality: ");
                                c->print(inst);
                                DebugOn(", violation = "<< diff << endl);
                                nb_viol++;
                                //                        violated = true;
                                if (*c->_all_lazy) {
                                    *c->_all_lazy = false;
                                    c->_all_satisfied = false;
                                    c->_violated[inst] = true;
                                    violated = true;
                                    c->_lazy[inst] = false;
                                }
                                else {
                                    //                            throw runtime_error("Non-lazy constraint is violated, solution declared optimal by solver!\n" + c->to_str(inst));
                                }
                            }
                            else if (abs(diff)>tol) {
                                c->_active[inst] = false;
                                //                        if (*c->_all_lazy) {
                                //                            c->_lazy[inst] = true;
                                //                        }
                            }
                            else {
                                nb_active++;
                            }
                        }
                        break;
                    case geq:
                        for (size_t inst=0; inst<nb_inst; inst++) {
                            c->_violated[inst] = false;
                            diff = c->eval(inst);
                            if(diff < -tol) {
                                DebugOff("Violated inequality: ");
                                //                        c->print(inst);
                                DebugOff(", violation = "<< diff << endl);
                                nb_viol++;
                                //                        violated = true;
                                if (*c->_all_lazy) {
                                    *c->_all_lazy = false;
                                    c->_all_satisfied = false;
                                    c->_violated[inst] = true;
                                    violated = true;
                                    c->_lazy[inst] = false;
                                }
                                else {
                                    //                            throw runtime_error("Non-lazy constraint is violated, solution declared optimal by solver!\n" + c->to_str(inst));
                                }
                            }
                            else if (abs(diff)> tol) {
                                c->_active[inst] = false;
                                //                        if (*c->_all_lazy) {
                                //                            c->_lazy[inst] = true;
                                //                        }
                            }
                            else {
                                nb_active++;
                            }
                        }
                        break;
                        
                    default:
                        break;
                }
                //        *c->_all_lazy = false;
                nb_viol_all += nb_viol;
                nb_active_all += nb_active;
                if (nb_viol>0) {
                    DebugOn("Percentage of violated constraints for " << c->get_name() << " = " << to_string(100.*nb_viol/nb_inst) << "%\n");
                }
                if (c->get_ctype()!=eq) {
                    DebugOff("Percentage of active constraints for " << c->get_name() << " = " << to_string(100.*nb_active/nb_inst) << "%\n");
                }
            }
            DebugOn("Total percentage of violated constraints = " << to_string(100.*nb_viol_all/_nb_cons) << "%\n");
            auto nb_ineq = get_nb_ineq();
            DebugOn("Total percentage of active constraints = " << to_string(100.*nb_active_all/nb_ineq) << "%\n");
            return violated;
        }
        
        
        
        /**
         Returns true if the current solution satisfies bounds and constraints upt to tolerance tol
         @param[in] tol tolerance for constraint satisfaction
         */
        bool is_feasible(type tol){            
            shared_ptr<param_> v;
            double viol = 0;
            bool feasible = true;
            for(auto& v_p: _vars)
            {
                v = v_p.second;
                for (size_t i = 0; i < v->get_dim(); i++) {
                    viol = v->get_lb_violation(i);
                    if(viol  > tol){
                        clog << "Violated lower bound for var: " << v->_name << ", index: "<< i << ", violation = " << viol << endl;
                        feasible = false;
                        
                    }
                    viol = v->get_ub_violation(i);
                    if(viol > tol){
                        clog << "Violated upper bound for var: " << v->_name << ", index: "<< i << ", violation = " << viol << endl;
                        feasible = false;
                    }
                }
            }
            feasible = feasible && !has_violated_constraints(tol);
            return feasible;
        }
        
        
#ifdef USE_BONMIN
        void fill_in_var_types(Bonmin::TMINLP::VariableType* var_types){
            size_t vid;
            param_* v;
            for(auto& v_p: _vars)
            {
                v = v_p.second;
                vid = v->get_id();
                auto bonmin_type = Bonmin::TMINLP::CONTINUOUS;
                auto type = v->get_intype();
                if(type==short_ || type==integer_){
                    bonmin_type = Bonmin::TMINLP::INTEGER;
                }
                else if(type==binary_){
                    bonmin_type = Bonmin::TMINLP::BINARY;
                }
                for (size_t i = 0; i < v->get_dim(); i++) {
                    var_types[vid+i] = bonmin_type;
                }
            }
                
        }
#endif
        
        /**
         Fill the lower and upper bound values in x_l and x_u
         @param[out] x_l lower bound values to fill
         @param[out] x_u upper bound values to fill
         */
        void fill_in_var_bounds(double* x_l ,double* x_u) {
            for(auto &v_p: _vars)
            {
                v_p.second->set_double_lb(x_l);
                v_p.second->set_double_ub(x_u);
            }
        }
        
        /**
         Initialize the model variables using values from x
         @param[in] x values to initialize to
         */
        void set_x(const double* x){
            for(auto &v_p: _vars)
            {
                v_p.second->get_double_val(x);
            }
        }
        
        /**
         Evaluate all nonlinear functions at current point
         */
        void compute_funcs() {
            auto it = _nl_funcs.begin();
            while (it!=_nl_funcs.end()) {
                auto f = (*it++);
                DebugOff(f->to_str() << endl);
                if (f->is_constant() && f->_evaluated) {
                    continue;
                }
                if (!f->is_constant() && (_type!=nlin_m)) {//no need to precompute these for polyomial programs
                    continue;
                }
                if (!f->is_matrix()) {
                    DebugOff(f->to_str()<<endl);
                    for (size_t inst = 0; inst < f->get_dim(); inst++) {
                        f->eval(inst);
                    }
                }
                else {
                    DebugOff(f->to_str()<<endl);
                    f->eval_matrix();
                }
                if (f->is_constant()) {
                    f->_evaluated = true;
                }
            }
        }
        
        void fill_in_obj(const double* x , double& res, bool new_x){
            if (new_x) {
                set_x(x);
                compute_funcs();
            }
            res = _obj->eval();
            DebugOff("Objective = " << to_string(res) << endl);
        }
        
        
        void fill_in_grad_obj(const double* x , double* res, bool new_x){
            param_* v;
            shared_ptr<func<type>> df;
            size_t vid, vid_inst, index = 0;
            unique_id v_unique;
            if (new_x) {
                set_x(x);
                compute_funcs();
            }
            for (size_t i = 0; i<_nb_vars; i++) {
                res[i] = 0;
            }
            if (_first_call_gard_obj) {
                _obj_grad_vals.resize(_obj->get_nb_vars());
                _first_call_gard_obj = false;
            }
            else if (_obj->is_linear()) {
                //    else if (false) { /* No need to recompute jacobian for linear objectives */
                for(auto& vi_p: _obj->get_vars())
                {
                    v = vi_p.second.first.get();
                    vid = v->get_id();
                    if (v->_is_vector) {
                        for (size_t i = 0; i < v->get_dim(); i++) {
                            vid_inst = vid + v->get_id_inst(i);
                            res[vid_inst] = _obj_grad_vals[index++];
                        }
                    }
                    else {
                        vid_inst = vid + v->get_id_inst();
                        res[vid_inst] = _obj_grad_vals[index++];
                    }
                }
                return;
            }
            for(auto& vi_p: _obj->get_vars()) {
                v = vi_p.second.first.get();
                vid = v->get_id();
                df = _obj->get_stored_derivative(v->_name);
                if (v->is_matrix()) {
                    //Unsupported yet.
                    throw invalid_argument("Matrices in the objective unsupported");
                }
                else if (v->_is_vector) {
                    for (size_t i = 0; i < v->get_dim(); i++) {
                        vid_inst = vid + v->get_id_inst(i);
                        res[vid_inst] = df->eval(i);
                        _obj_grad_vals[index++] =res[vid_inst];
                    }
                }
                else {
                    vid_inst = vid + v->get_id_inst();
                    res[vid_inst] = df->eval();
                    _obj_grad_vals[index++] =res[vid_inst];
                }
            }
        }
        
        //
        //template<typename T1>
        //void add_McCormick(std::string name, const var<T1>& vlift, const var<T1>& v1, const var<T1>& v2) {
        //    Constraint<type> MC1(name+"_McCormick1");
        //    auto lb1 = v1.get_lb(v1.get_id_inst());
        //    auto lb2 = v2.get_lb(v2.get_id_inst());
        //    auto ub1 = v1.get_ub(v1.get_id_inst());
        //    auto ub2 = v2.get_ub(v2.get_id_inst());
        //    bool template_cstr = v1->_dim[0]>1;
        //    MC1 += vlift;
        //    if(template_cstr){//Template constraint
        //        MC1 -= (*v1._lb)*v2 + (*v2._lb)*v1 - (*v1._lb)*(*v2._lb);
        //    }
        //    else {
        //        MC1 -= lb1*v2 + lb2*v1 - lb1*lb2;
        //    }
        //    MC1 >= 0;
        //    add(MC1);
        //    //    MC1.print();
        //    Constraint<type> MC2(name+"_McCormick2");
        //    MC2 += vlift;
        //    if(template_cstr){//Template constraint
        //        MC2 -= (*v1._ub)*v2 + (*v2._ub)*v1 - (*v1._ub)*(*v2._ub);
        //    }
        //    else {
        //        MC2 -= ub1*v2 + ub2*v1 - ub1*ub2;
        //    }
        //    MC2 >= 0;
        //    add(MC2);
        ////    //    MC2.print();
        //    Constraint<type> MC3(name+"_McCormick3");
        //    MC3 += vlift;
        //    if(template_cstr){//Template constraint
        //        MC3 -= (*v1._lb)*v2 + (*v2._ub)*v1 - (*v1._lb)*(*v2._ub);
        //    }
        //    else {
        //        MC3 -= lb1*v2 + ub2*v1 - lb1*ub2;
        //    }
        //    MC3 <= 0;
        //    add(MC3);
        ////    //    MC3.print();
        //    Constraint<type> MC4(name+"_McCormick4");
        //    MC4 += vlift;
        //    if(template_cstr){//Template constraint
        //        MC4 -= (*v1._ub)*v2 + (*v2._lb)*v1 - (*v1._ub)*(*v2._lb);
        //    }
        //    else{
        //        MC4 -= ub1*v2 + lb2*v1 - ub1*lb2;
        //    }
        //    MC4 <= 0;
        //    add(MC4);
        //    //    MC4.print();
        //}
        
        
        /** Build the sequential McCormick relaxation for polynomial programs **/
        
        shared_ptr<Model<type>> build_McCormick(){
            replace_integers();
            if (_type==nlin_m) {
                cerr << "Can only build a McCormick relaxation for polynomial programs, returning null" << endl;
                return nullptr;
            }
            if (_type==lin_m) {
                cerr << "No need to build a McCormick relaxation for a linear program, returning null" << endl;
                return nullptr;
            }
            shared_ptr<Model> Mc = make_shared<Model>("McCormick Relaxation");
            //    shared_ptr<Constraint<type>> cstr;
            //    param_* v;
            //    for (auto &var_p:_vars) {
            //        v = var_p.second;
            //
            //        switch (v->get_intype()) {// TODO check for other types
            //            case double_:
            //                Mc->add(*(var<double>*)v);
            //                break;
            //            case integer_:
            //                Mc->add(*(var<int>*)v);
            //                break;
            //            case binary_:
            //                Mc->add(*(var<bool>*)v);
            //                break;
            //            default:
            //                break;
            //        }
            //    }
            //    if ((_obj->is_convex() && _objt==minimize) || (_obj->is_concave() && _objt==maximize)) {
            //        Mc->_obj = _obj;
            //    }
            //    else{
            //        func<double> new_obj = *_obj->get_cst();
            //        Mc->_objt = _objt;
            //        for(auto &lt : _obj->get_lterms()){
            //            new_obj->insert(lt.second);
            //        }
            //        for(auto &qt : _obj->get_qterms()){
            //
            //            if ((_objt==maximize && _obj->get_convexity(qt.second)==concave_) || (_objt==minimize && _obj->get_convexity(qt.second)==convex_)) {
            //                new_obj->insert(qt.second);
            //            }
            //            else {//Lift products
            //                string new_name;
            //                if(qt.second._p->first==qt.second._p->second){
            //                    new_name  = qt.second._p->first->_name+"²_lifted";
            //                }
            //                else {
            //                    new_name = qt.second._p->first->_name+"_"+qt.second._p->second->_name+"_lifted";
            //                }
            //                var<> v(new_name);
            //                auto new_v = Mc->get_var_ptr(new_name);
            //                if(new_v == nullptr){
            //                    Mc->add(v);
            //                    new_v = (param_*)&v;
            //                }
            //                Mc->add_McCormick(new_v->_name, new_v, qt.second._p->first, qt.second._p->second);
            //                if(qt.second._sign) {
            //                    new_obj += *new_v;
            //                }
            //                else {
            //                    new_obj -= *new_v;
            //                }
            //            }
            //        }
            //    }
            //    for (auto &cst_p:_cons) {
            //        cstr = cst_p.second;
            //        if (cstr->is_convex()) {
            //            Mc->add(*cstr);
            //        }
            //        else if(cstr->func_::is_convex()){
            //            Constraint new_cstr(*cstr);
            //            new_cstr._ctype = leq;
            //            Mc->add(new_cstr);
            //        }
            //        else if(cstr->func_::is_concave()){
            //            Constraint new_cstr(*cstr);
            //            new_cstr._ctype = geq;
            //            Mc->add(new_cstr);
            //        }
            //        else{
            //            Constraint new_cstr(cstr->get_name()+"_lifted");
            //            new_cstr += *cstr->get_cst();
            //            new_cstr._ctype = cstr->_ctype;
            //            new_cstr._rhs = cstr->get_rhs();
            //            for(auto &lt : cstr->get_lterms()){
            //                new_cstr.insert(lt.second);
            //            }
            //            for(auto &qt : cstr->get_qterms()){
            //
            //                if ((cstr->_ctype==geq && cstr->get_convexity(qt.second)==concave_) || (cstr->_ctype==leq && cstr->get_convexity(qt.second)==convex_)) {
            //                    new_cstr.insert(qt.second);
            //                }
            //                else {//Lift products
            //                    string new_name;
            //                    if(qt.second._p->first==qt.second._p->second){
            //                        new_name  = qt.second._p->first->_name+"²_lifted";
            //                    }
            //                    else {
            //                        new_name = qt.second._p->first->_name+"_"+qt.second._p->second->_name+"_lifted";
            //                    }
            //                    var<> v(new_name);
            //                    auto new_v = Mc->get_var_ptr(new_name);
            //                    if(new_v == nullptr){
            //                        Mc->add(v);
            //                        new_v = (param_*)&v;
            //                    }
            //                    Mc->add_McCormick(new_v->_name, new_v, qt.second._p->first, qt.second._p->second);
            //                    if(qt.second._sign) {
            //                        new_cstr += *new_v;
            //                    }
            //                    else {
            //                        new_cstr -= *new_v;
            //                    }
            //                }
            //            }
            //            //            for(auto &qt : cstr->get_qterms()){
            //            //TODO polynomial part
            //            Mc->add(new_cstr);
            //        }
            //    }
            return Mc;
        }
        
        
        
        
        
        void fill_in_cstr(const double* x , double* res, bool new_x){
            if (new_x) {
                set_x(x);
                compute_funcs();
            }
            unsigned nr_threads = std::thread::hardware_concurrency()/2;
            if (nr_threads==0) {
                nr_threads = 1;
            }
            vector<thread> threads;
            /* Split cons into nr_threads parts */
            vector<size_t> limits = bounds(nr_threads, _cons_vec.size());
            /* Launch all threads in parallel */
            for (unsigned i = 0; i < nr_threads; ++i) {
                threads.push_back(thread(compute_constrs<type>, ref(_cons_vec), res, limits[i], limits[i+1]));
            }
            /* Join the threads with the main thread */
            for(auto &t : threads){
                t.join();
            }
        }
        
        
        
        void fill_in_jac_nnz(int* iRow , int* jCol){
            size_t idx=0, id = 0;
            size_t cid = 0;
            size_t vid = 0;
            Constraint<type>* c = NULL;
            param_* v = NULL;
            /* return the structure of the jacobian */
            for(auto& c_p :_cons)
            {
                c = c_p.second.get();
                c->_jac_cstr_idx = idx;
                auto nb_ins = c->_dim[0];
                for (auto &v_p: c->get_vars()){
                    v = v_p.second.first.get();
                    vid = v->get_id();
                    id = 0;
                    for (size_t inst = 0; inst< nb_ins; inst++){
                        if (!*c->_all_lazy || !c->_lazy[inst]) {
                            cid = c->_id+id++;
                            if (v->_is_vector) {
                                auto dim = v->get_dim(inst);
                                for (size_t j = 0; j<dim; j++) {
                                    iRow[idx] = cid;
                                    jCol[idx] = vid + v->get_id_inst(inst, j);
                                    idx++;
                                }
                            }
                            else {
                                iRow[idx] = cid;
                                jCol[idx] = vid + v->get_id_inst(inst);
                                idx++;
                            }
                        }
                    }
                }
            }
            if (idx!=_nnz_g) {
                throw invalid_argument("idx!=_nnz_g");
            }
        }
        
        
        
        /* Fill the nonzero values in the jacobian */
        
        void fill_in_jac(const double* x , double* res, bool new_x){
            //    if (!_first_call_jac && (!new_x || _type==lin_m)) { /* No need to recompute jacobian for linear models */
            if (new_x) {
                set_x(x);
                compute_funcs();
            }
            
            if (!_first_call_jac && (_type==lin_m)) { /* No need to recompute jacobian for linear models */
                for (size_t i = 0; i< _nnz_g; i++) {
                    res[i] = _jac_vals[i];
                }
                return;
            }
            size_t idx=0;
            size_t cid = 0;
            string vid;
            Constraint<type>* c = NULL;
            param_* v = NULL;
            shared_ptr<func<type>> dfdx;
            //    vector<Constraint*> cons;
            if (_type!=nlin_m) {//Polynomial, Quadratic or Linear
                _cons_vec.clear();
                for(auto& c_p: _cons)
                {
                    c_p.second->_new = false;
                    _cons_vec.push_back(c_p.second);
                }
                
                //        compute_jac(_cons_vec, res, 0, _cons_vec.size(), _first_call_jac, _jac_vals);
                unsigned nr_threads = std::thread::hardware_concurrency()/2;
                if (nr_threads==0) {
                    nr_threads = 1;
                }
                vector<thread> threads;
                /* Split cons into nr_threads parts */
                vector<size_t> limits = bounds(nr_threads, _cons_vec.size());
                
                /* Launch all threads in parallel */
                for (unsigned i = 0; i < nr_threads; ++i) {
                    threads.push_back(thread(compute_jac<type>, ref(_cons_vec), res, limits[i], limits[i+1], _first_call_jac, ref(_jac_vals)));
                }
                /* Join the threads with the main thread */
                for(auto &t : threads){
                    t.join();
                }
                _first_call_jac = false;
                return;
            }
            for(auto& c_p :_cons)
            {
                c = c_p.second.get();
                auto nb_ins = c->_dim[0];
                size_t id = 0;
                if (c->is_linear() && !_first_call_jac) {
                    //        if (false) {
                    DebugOff("Linear constraint, using stored jacobian!\n");
                    for (size_t i = 0; i<nb_ins; i++) {
                        if (!*c->_all_lazy || !c->_lazy[i]) {
                            for (size_t j = 0; j<c->get_nb_vars(i); j++) {
                                res[idx] = _jac_vals[idx];
                                idx++;
                            }
                        }
                    }
                }
                else {
                    //            if (_type==nlin_m) {
                    for (auto &v_p: c->get_vars()){
                        v = v_p.second.first.get();
                        vid = v->_name;
                        dfdx = c->get_stored_derivative(vid);
                        if (dfdx->is_number()) {
                            for (size_t inst = 0; inst< nb_ins; inst++){
                                if (!*c->_all_lazy || !c->_lazy[inst]) {
                                    cid = c->_id+id++;
                                    if (v->_is_vector) {
                                        auto dim = v->get_dim(inst);
                                        for (size_t j = 0; j<dim; j++) {
                                            res[idx] = c->eval(dfdx,0);
                                            _jac_vals[idx] = res[idx];
                                            idx++;
                                        }
                                    }
                                    else {
                                        res[idx] = dfdx->eval();
                                        _jac_vals[idx] = res[idx];
                                        idx++;
                                    }
                                }
                            }
                        }
                        else {
                            for (size_t inst = 0; inst< nb_ins; inst++){
                                if (!*c->_all_lazy || !c->_lazy[inst]) {
                                    cid = c->_id+id++;
                                    if (v->_is_vector) {
                                        auto dim = v->get_dim(inst);
                                        if (dfdx->is_matrix()) {
                                            for (size_t j = 0; j<dim; j++) {
                                                res[idx] = dfdx->eval(j,inst);
                                                _jac_vals[idx] = res[idx];
                                                idx++;
                                            }
                                        }
                                        else {
                                            for (size_t j = 0; j<dim; j++) {
                                                res[idx] = dfdx->eval(j);//TODO check double indexed funcs
                                                _jac_vals[idx] = res[idx];
                                                DebugOff("jac_val["<< idx <<"] = " << _jac_vals[idx] << endl);
                                                idx++;
                                            }
                                        }
                                    }
                                    else {
                                        res[idx] = dfdx->eval(inst);
                                        _jac_vals[idx] = res[idx];
                                        idx++;
                                    }
                                }
                            }
                        }
                    }
                    //            }
                    //            else {
                    //                for (auto &v_p: c->get_vars()){
                    //                    v = v_p.second.first.get();
                    //                    vid = v->_name;
                    //                    dfdx = c->get_stored_derivative(vid);
                    //                    for (size_t inst = 0; inst< nb_ins; inst++){
                    //                        cid = c->_id+inst;
                    //                        if (v->_is_vector) {
                    //                            auto dim = v->get_dim(inst);
                    //                            for (size_t j = 0; j<dim; j++) {
                    //                                res[idx] = dfdx->eval(inst,j);
                    //                                _jac_vals[idx] = res[idx];
                    //                                idx++;
                    //                            }
                    //                        }
                    //                        else {
                    //                            res[idx] = dfdx->eval(inst);
                    //                            _jac_vals[idx] = res[idx];
                    //                            idx++;
                    //                        }
                    //                    }
                    //                }
                    //            }
                }
            }
            _first_call_jac = false;
        }
        
        
        
        
#ifdef USE_IPOPT
        
        void fill_in_var_linearity(Ipopt::TNLP::LinearityType* param_types){
            size_t vid = 0;
            bool linear = true;
            for(auto& vi: _vars)
            {
                vid = vi.second->get_id();
                for (size_t i = 0; i < vi.second->get_dim(); i++) {
                    //            linear = true;
                    //            for(auto &c: _v_in_cons[vid + i])
                    //            {
                    //                if (!c->is_linear()) {
                    //                    linear=false;
                    //                }
                    //            }
                    if (linear) param_types[vid + i]=Ipopt::TNLP::LINEAR;
                    else param_types[vid + i] = Ipopt::TNLP::NON_LINEAR;
                }
            }
        }
        
        
        void fill_in_cstr_linearity(Ipopt::TNLP::LinearityType* const_types){
            Constraint<type>* c = nullptr;
            bool lin = false;
            size_t cid = 0, id = 0;
            for(auto& c_p :_cons)
            {
                c = c_p.second.get();
                id = 0;
                if (c->is_linear() || c->is_constant()) {
                    lin = true;
                }
                else {
                    lin = false;
                }
                auto nb_ins = c->_dim[0];
                for (size_t i = 0; i< nb_ins; i++){
                    if (!*c->_all_lazy || !c->_lazy[i]) {
                        cid = c->_id+id++;
                        if (lin) {
                            const_types[cid]=Ipopt::TNLP::LINEAR;
                        }
                        else {
                            const_types[cid] = Ipopt::TNLP::NON_LINEAR;
                        }
                    }
                }
            }
        }
#endif
        
        
        void fill_in_hess_nnz(int* iRow , int* jCol){
            size_t idx = 0, idx_all=0, idx_pair=0, vid, vjd;
            string vi_name, vj_name;
            shared_ptr<param_> vi;
            shared_ptr<param_> vj;
            shared_ptr<Constraint<type>> c;
            for (auto &pairs: _hess_link) {
                vi_name = pairs.first.first;
                vj_name = pairs.first.second;
                vi = (pairs.second.begin())->first->get_var(vi_name);
                vj = (pairs.second.begin())->first->get_var(vj_name);
                if (vi_name.compare(vj_name) > 0) {//ONLY STORE LOWER TRIANGULAR PART OF HESSIAN
                    throw invalid_argument("SHOULD BE SORTED CORRECTLY IN FILL_MAPS");
                }
                vid = vi->get_id();
                vjd = vj->get_id();
                idx_pair = idx;
                //        auto max_f_idx = 0;
                for (auto &f_pair:pairs.second) {
                    //            auto f_idx = 0;
                    //            idx = idx_pair;
                    //        auto f_pair = *pairs.second.begin();
                    auto f = f_pair.first;
                    if (f->_is_constraint) {
                        c = static_pointer_cast<Constraint<type>>(f);
                    }
                    auto d2f = f_pair.second;
                    size_t nb_inst = f->_dim[0];
                    for (size_t inst = 0; inst<nb_inst; inst++) {
                        if (!(f->_is_constraint && *c->_all_lazy && c->_lazy[inst])) {
                            
                            if (d2f->is_matrix()) {
                                for (size_t i = 0; i < d2f->_dim[0]; i++) {
                                    for (size_t j = i; j < d2f->_dim[1]; j++) {
                                        idx_all++;
                                        iRow[idx] = vid + vi->get_id_inst(i);
                                        jCol[idx] = vjd + vj->get_id_inst(j);
                                        idx++;
                                        //                                f_idx++;
                                    }
                                }
                            }
                            else if(d2f->_is_vector){
                                //                    if (d2f->_dim[0] != d2f->_dim[0]) {
                                //                        throw invalid_argument("error");
                                //                    }
                                for (size_t j = 0; j < d2f->_dim[0]; j++) {
                                    idx_all++;
                                    iRow[idx] = vid + vi->get_id_inst(j);
                                    jCol[idx] = vjd + vj->get_id_inst(j);
                                    idx++;
                                    //                            f_idx++;
                                }
                            }
                            else {
                                idx_all++;
                                iRow[idx] = vid + vi->get_id_inst(inst);
                                jCol[idx] = vjd + vj->get_id_inst(inst);
                                idx++;
                                //                        f_idx++;
                            }
                        }
                    }
                    //            if(max_f_idx < f_idx){
                    //                max_f_idx = f_idx;
                    //            }
                }
                //        idx = idx_pair+max_f_idx;
            }
            if (idx!=_nnz_h) {
                throw invalid_argument("idx!=_nnz_h");
            }
            _hess_vals.resize(idx_all);
        }
        
        
        void fill_in_hess(const double* x , double obj_factor, const double* lambda, double* res, bool new_x){
            size_t idx = 0, idx_in = 0, c_inst = 0, idx_pair=0;
            shared_ptr<Constraint<type>> c;
            bool idx_inc = false;
            double hess = 0;
            for (size_t i = 0; i<_nnz_h; i++) {
                res[i] = 0;
            }
            if (new_x) {
                set_x(x);
                compute_funcs();
            }
            if (_first_call_hess) {
                for (auto &pairs: _hess_link) {
                    idx_pair = idx;
                    //            auto max_f_idx = 0;
                    for (auto &f_pair:pairs.second) {
                        //                auto f_idx = 0;
                        //                idx = idx_pair;
                        idx_inc = false;
                        auto f = f_pair.first;
                        if (f->_is_constraint) {
                            c = static_pointer_cast<Constraint<type>>(f);
                        }
                        auto d2f = f_pair.second;
                        size_t nb_inst = f->_dim[0];
                        size_t id_inst = 0;
                        for (size_t inst = 0; inst<nb_inst; inst++) {
                            if (f->_is_constraint) {
                                if (!*c->_all_lazy || !c->_lazy[inst]) {
                                    if (c->is_nonlinear()) {
                                        idx_inc = false;
                                        c_inst = c->get_id_inst(id_inst++);
                                        //                            if(f->is_nonlinear()){
                                        if (d2f->is_matrix()) {
                                            for (size_t i = 0; i < d2f->_dim[0]; i++) {
                                                for (size_t j = i; j < d2f->_dim[1]; j++) {
                                                    hess = c->eval(d2f,i,j);
                                                    _hess_vals[idx_in++] = hess;
                                                    res[idx++] += lambda[c->_id + c_inst] * hess;
                                                    //                                        f_idx++;
                                                    idx_inc = true;
                                                }
                                            }
                                        }
                                        else if(d2f->_is_vector){
                                            for (size_t j = 0; j < d2f->_dim[0]; j++) {
                                                hess = d2f->get_val(j);
                                                _hess_vals[idx_in++] = hess;
                                                res[idx] += lambda[c->_id + c_inst] * hess;
                                                idx++;
                                                //                                    f_idx++;
                                                idx_inc = true;
                                                //                    }
                                            }
                                        }
                                        else {
                                            if (d2f->is_number()) {
                                                hess = d2f->_val->at(0);
                                            }
                                            else {
                                                hess = d2f->_val->at(inst);
                                            }
                                            _hess_vals[idx_in++] = hess;
                                            res[idx] += lambda[c->_id + c_inst] * hess;
                                        }
                                    }
                                    else {
                                        idx_inc = false;
                                        c_inst = c->get_id_inst(id_inst++);
                                        //                            if(f->is_nonlinear()){
                                        if (d2f->is_matrix()) {
                                            for (size_t i = 0; i < d2f->_dim[0]; i++) {
                                                for (size_t j = i; j < d2f->_dim[1]; j++) {
                                                    hess = d2f->eval(i,j);
                                                    _hess_vals[idx_in++] = hess;
                                                    res[idx++] += lambda[c->_id + c_inst] * hess;
                                                    //                                        f_idx++;
                                                    idx_inc = true;
                                                }
                                            }
                                        }
                                        else if(d2f->_is_vector){
                                            for (size_t j = 0; j < d2f->_dim[0]; j++) {
                                                hess = d2f->eval(j);
                                                _hess_vals[idx_in++] = hess;
                                                res[idx] += lambda[c->_id + c_inst] * hess;
                                                idx++;
                                                //                                    f_idx++;
                                                idx_inc = true;
                                                //                    }
                                            }
                                        }
                                        else {
                                            if (d2f->is_number()) {
                                                hess = d2f->_val->at(0);
                                            }
                                            else {
                                                hess = d2f->eval(inst);
                                            }
                                            _hess_vals[idx_in++] = hess;
                                            res[idx] += lambda[c->_id + c_inst] * hess;
                                        }
                                    }
                                }
                                else {
                                    idx_inc = true;
                                }
                            }
                            else {
                                if (d2f->is_matrix()) {
                                    for (size_t i = 0; i < d2f->_dim[0]; i++) {
                                        for (size_t j = i; j < d2f->_dim[1]; j++) {
                                            hess = d2f->eval(i,j);
                                            _hess_vals[idx_in++] = hess;
                                            res[idx++] += obj_factor * hess;
                                            //                                        f_idx++;
                                            idx_inc = true;
                                        }
                                    }
                                }
                                else if(d2f->_is_vector){
                                    for (size_t j = 0; j < d2f->_dim[0]; j++) {
                                        hess = d2f->eval(j);
                                        _hess_vals[idx_in++] = hess;
                                        res[idx] += obj_factor * hess;
                                        idx++;
                                        //                                f_idx++;
                                        idx_inc = true;
                                    }
                                }
                                else {
                                    hess = d2f->eval(0);
                                    _hess_vals[idx_in++] = hess;
                                    res[idx] += obj_factor * hess;
                                }
                            }
                            if (!idx_inc) {
                                idx++;
                                //                        f_idx++;
                            }
                        }
                        //                if(max_f_idx < f_idx){
                        //                    max_f_idx = f_idx;
                        //                }
                    }
                    //            idx = idx_pair+max_f_idx;
                }
                _first_call_hess = false;
                return;
            }
            if ((_type==lin_m || _type==quad_m)) { /* No need to recompute Hessian for quadratic models, used stored values */
                size_t id_inst = 0;
                for (auto &pairs: _hess_link) {
                    idx_pair = idx;
                    //            auto max_f_idx = 0;
                    for (auto &f_pair:pairs.second) {
                        //                auto f_idx = 0;
                        //                idx = idx_pair;
                        idx_inc = false;
                        auto f = f_pair.first;
                        if (f->_is_constraint) {
                            c = static_pointer_cast<Constraint<type>>(f);
                        }
                        auto d2f = f_pair.second;
                        size_t nb_inst = f->_dim[0];
                        id_inst = 0;
                        for (size_t inst = 0; inst<nb_inst; inst++) {
                            if (f->_is_constraint) {
                                if (!*c->_all_lazy || !c->_lazy[inst]) {
                                    idx_inc = false;
                                    c_inst = c->get_id_inst(id_inst++);
                                    if (d2f->is_matrix()) {
                                        for (size_t i = 0; i < d2f->_dim[0]; i++) {
                                            for (size_t j = i; j < d2f->_dim[1]; j++) {
                                                res[idx] += lambda[c->_id + c_inst] * _hess_vals[idx_in++];
                                                idx++;
                                                //                                        f_idx++;
                                                idx_inc = true;
                                            }
                                        }
                                    }
                                    else if(d2f->_is_vector){
                                        for (size_t j = 0; j < d2f->_dim[0]; j++) {
                                            res[idx] += lambda[c->_id + c_inst] * _hess_vals[idx_in++];
                                            idx++;
                                            //                                    f_idx++;
                                            idx_inc = true;
                                        }
                                    }
                                    else {
                                        res[idx] += lambda[c->_id + c_inst] * _hess_vals[idx_in++];
                                    }
                                }
                                else {
                                    idx_inc = true;
                                }
                            }
                            else {
                                if(d2f->_is_vector){
                                    for (size_t j = 0; j < d2f->_dim[0]; j++) {
                                        res[idx] += obj_factor * _hess_vals[idx_in++];
                                        idx++;
                                        //                                f_idx++;
                                        idx_inc = true;
                                    }
                                }
                                else {
                                    res[idx] += obj_factor * _hess_vals[idx_in++];
                                }
                            }
                            if (!idx_inc) {
                                idx++;
                                //                        f_idx++;
                            }
                        }
                        //                if(max_f_idx < f_idx){
                        //                    max_f_idx = f_idx;
                        //                }
                    }
                    //            idx = idx_pair+max_f_idx;
                }
                return;
            }
            for (auto &pairs: _hess_link) {
                idx_pair = idx;
                //        auto max_f_idx = 0;
                for (auto &f_pair:pairs.second) {
                    //            auto f_idx = 0;
                    //            idx = idx_pair;
                    idx_inc = false;
                    auto f = f_pair.first;
                    if (f->_is_constraint) {
                        c = static_pointer_cast<Constraint<type>>(f);
                    }
                    auto d2f = f_pair.second;
                    size_t nb_inst = f->_dim[0];
                    size_t id_inst = 0;
                    for (size_t inst = 0; inst<nb_inst; inst++) {
                        if (f->_is_constraint) {
                            if (!*c->_all_lazy || !c->_lazy[inst]) {
                                idx_inc = false;
                                c_inst = c->get_id_inst(id_inst++);
                                if (c->is_quadratic()) {
                                    if (d2f->is_matrix()) {
                                        for (size_t i = 0; i < d2f->_dim[0]; i++) {
                                            for (size_t j = i; j < d2f->_dim[1]; j++) {
                                                res[idx++] += lambda[c->_id + c_inst] * _hess_vals[idx_in++];
                                                //                                        f_idx++;
                                                idx_inc = true;
                                            }
                                        }
                                    }
                                    else if(d2f->_is_vector){
                                        for (size_t j = 0; j < d2f->_dim[0]; j++) {
                                            res[idx] += lambda[c->_id + c_inst] * _hess_vals[idx_in++];
                                            idx++;
                                            //                                    f_idx++;
                                            idx_inc = true;
                                        }
                                    }
                                    else {
                                        res[idx] += lambda[c->_id + c_inst] * _hess_vals[idx_in++];
                                    }
                                }
                                else if (d2f->is_matrix()) {
                                    if (c->is_nonlinear()) {
                                        for (size_t i = 0; i < d2f->_dim[0]; i++) {
                                            for (size_t j = i; j < d2f->_dim[1]; j++) {
                                                hess = d2f->eval(i,j);
                                                _hess_vals[idx_in++] = hess;
                                                res[idx++] += lambda[c->_id + c_inst] * hess;
                                                //                                    f_idx++;
                                                idx_inc = true;
                                            }
                                        }
                                    }
                                    else {
                                        for (size_t i = 0; i < d2f->_dim[0]; i++) {
                                            for (size_t j = i; j < d2f->_dim[1]; j++) {
                                                hess = d2f->eval(i,j);
                                                _hess_vals[idx_in++] = hess;
                                                res[idx++] += lambda[c->_id + c_inst] * hess;
                                                //                                    f_idx++;
                                                idx_inc = true;
                                            }
                                        }
                                    }
                                }
                                else if(d2f->_is_vector){
                                    if (c->is_nonlinear()) {
                                        for (size_t j = 0; j < d2f->_dim[0]; j++) {
                                            hess = d2f->get_val(j);
                                            _hess_vals[idx_in++] = hess;
                                            res[idx] += lambda[c->_id + c_inst] * hess;
                                            idx++;
                                            //                                f_idx++;
                                            idx_inc = true;
                                        }
                                    }
                                    else {
                                        for (size_t j = 0; j < d2f->_dim[0]; j++) {
                                            hess = d2f->eval(j);
                                            _hess_vals[idx_in++] = hess;
                                            res[idx] += lambda[c->_id + c_inst] * hess;
                                            idx++;
                                            //                                f_idx++;
                                            idx_inc = true;
                                        }
                                        
                                    }
                                }
                                else {
                                    if (d2f->is_number()) {
                                        hess = d2f->_val->at(0);
                                    }
                                    else {
                                        if (c->is_nonlinear()) {
                                            hess = d2f->_val->at(inst);
                                        }
                                        else {
                                            hess = d2f->eval(inst);
                                        }
                                    }
                                    _hess_vals[idx_in++] = hess;
                                    res[idx] += lambda[c->_id + c_inst] * hess;
                                }
                            }
                            else {
                                idx_inc = true;
                            }
                        }
                        else {
                            if (_obj->is_quadratic()) {
                                if(d2f->_is_vector){
                                    for (size_t j = 0; j < d2f->_dim[0]; j++) {
                                        res[idx] += obj_factor * _hess_vals[idx_in++];
                                        idx++;
                                        //                                f_idx++;
                                        idx_inc = true;
                                    }
                                }
                                else {
                                    res[idx] += obj_factor * _hess_vals[idx_in++];
                                }
                            }
                            else if(d2f->_is_vector){
                                for (size_t j = 0; j < d2f->_dim[0]; j++) {
                                    //                    for (size_t j = i; j < (pairs.second.begin())->second->_dim[1]; j++) {
                                    hess = d2f->eval(j);
                                    _hess_vals[idx_in++] = hess;
                                    res[idx] += obj_factor * hess;
                                    idx++;
                                    //                            f_idx++;
                                    idx_inc = true;
                                    //                    }
                                }
                            }
                            else {
                                hess = d2f->eval(0);
                                _hess_vals[idx_in++] = hess;
                                res[idx] += obj_factor * hess;
                            }
                        }
                        if (!idx_inc) {
                            idx++;
                            //                    f_idx++;
                        }
                    }
                    //            if(max_f_idx < f_idx){
                    //                max_f_idx = f_idx;
                    //            }
                }
                //        idx = idx_pair+max_f_idx;
            }
        }
        
        
        
        
        
        
//        void reset_funcs() {
//            for (auto& f:_nl_funcs) {
//                f->reset_val();
//            }
//            for (auto& c:_cons) {
//                c.second->reset_val();
//            }
//            //    _obj->reset_val();
//        }
        
        void fill_in_maps() {/*< Fill the _hess and _v_in_ maps to link variables with their constraints and compute the Jacobian & Hessian matrices */
            string vi_name, vj_name;
            param_* vi;
            param_* vj;
            
            _built = true;
            
            
            if (_obj->_new) {
                _obj->compute_derivatives();
                if (!_obj->is_linear()) {
                    for (auto &vi_p: _obj->get_vars()) {
                        vi = vi_p.second.first.get();
                        vi_name = vi_p.first;
                        //            vid = vi->get_id();
                        auto df = _obj->get_stored_derivative(vi->_name);
                        for (auto &vj_p: df->get_vars()) {
                            vj = vj_p.second.first.get();
                            vj_name = vj_p.first;
                            //                vjd = vj->get_id();
                            if (vi_name.compare(vj_name) < 0) {//ONLY STORE LOWER TRIANGULAR PART OF HESSIAN
                                _hess_link[make_pair<>(vi_name,vj_name)].insert(make_pair<>(_obj,_obj->get_stored_derivative(vi->_name)->get_stored_derivative(vj->_name)));
                            }
                            else {
                                _hess_link[make_pair<>(vj_name,vi_name)].insert(make_pair<>(_obj,_obj->get_stored_derivative(vj->_name)->get_stored_derivative(vi->_name)));
                            }
                            //                for (int inst = 0; inst<vi->get_dim(); inst++) {
                            //                    vid = vi->get_id();
                            //                    vid_inst = vid + vi->get_id_inst(inst);
                            //                    vjd = vj->get_id();
                            //                    vjd_inst = vjd + vj->get_id_inst(inst);
                            //                    _hess.insert(make_pair<>(vid_inst, vjd_inst));
                            //                }
                        }
                    }
                }
            }
            shared_ptr<Constraint<type>> c = nullptr;
            for(auto& c_p :_cons)
            {
                c = c_p.second;
                //        c->print();
                if (c->_new) {
                    c->compute_derivatives();
                    //            if (_type==nlin_m) {
                    for (auto &df_p:*c->get_dfdx()) {
                        auto df = static_pointer_cast<func<type>>(df_p.second);
                        DebugOff(df->to_str() << endl);
                        if (df->get_expr() || _type==nlin_m) {
                            df_p.second = embed(df);
                        }
                        else {
                            embed(df);
                        }
                        for (auto &df2_p:*df_p.second->get_dfdx()) {
                            //                            if (df2_p.second->get_expr()) {
                            df2_p.second = embed(df2_p.second);
                            //                            }
                        }
                    }
                    //            }
                    if (!c->is_linear()) {
                        for (auto &vi_p: c->get_vars()) {
                            vi = vi_p.second.first.get();
                            vi_name = vi_p.first;
                            auto df = c->get_stored_derivative(vi->_name);
                            for (auto &vj_p: df->get_vars()) {
                                vj = vj_p.second.first.get();
                                vj_name = vj_p.first;
                                if (vi_name.compare(vj_name) <= 0) {//ONLY STORE LOWER TRIANGULAR PART OF HESSIAN
                                    _hess_link[make_pair<>(vi_name,vj_name)].insert(make_pair<>(c, c->get_stored_derivative(vi->_name)->get_stored_derivative(vj->_name)));
                                }
                                else {
                                    _hess_link[make_pair<>(vj_name,vi_name)].insert(make_pair<>(c, c->get_stored_derivative(vj->_name)->get_stored_derivative(vi->_name)));
                                }
                            }
                        }
                    }
                    DebugOff(c->to_str() << endl);
                }
            }
            //    print_nl_functions();
        }
        
        
        void fill_in_duals(double* lambda, double* z_L, double* z_U){
            for (auto &cp: _cons) {
                size_t idx = 0;
                //        for (size_t inst = 0; inst < cp.second->_dim[0]; inst++) {
                //            if (!*cp.second->_all_lazy || !cp.second->_lazy[inst]) {
                //                lambda[cp.second->_id + idx++] = 0;
                //            }
                //
                //        }
                idx = 0;
                for (size_t inst = 0; inst < cp.second->_dual.size(); inst++) {
                    if (!*cp.second->_all_lazy || !cp.second->_lazy[inst]) {
                        auto index = cp.second->_id + idx++;
                        lambda[index] = cp.second->_dual[inst];
                    }
                    //            else
                    //            {
                    //                lambda[index] = 100;
                    //            }
                }
            }
            for (auto &vp: _vars) {
                auto nb_inst = vp.second->get_dim();
                auto vid = vp.second->get_id();
                for (size_t inst = 0; inst < nb_inst; inst++) {
                    auto id_inst = vp.second->get_id_inst(inst);
                    //            z_L[vid + id_inst] = vp.second->_l_dual[inst];
                    //            z_U[vid + id_inst] = vp.second->_u_dual[inst];
                    //            z_L[vp.second->get_id() + vp.second->get_id_inst(inst)] = 0;
                    //            z_U[vp.second->get_id() + vp.second->get_id_inst(inst)] = 0;
                }
            }
            
        }
        
        template<typename T=type, typename=enable_if<is_arithmetic<T>::value>>
        void fill_in_var_init(double* x) {
            for(auto& v_p: _vars){
                v_p.second->set_double_val(x);
            }
        }
        
        
        void fill_in_cstr_bounds(double* g_l ,double* g_u) {
            size_t cid = 0;
            Constraint<type>* c = NULL;
            for(auto& c_p :_cons)
            {
                c = c_p.second.get();
                switch (c->get_ctype()) {
                    case eq:{
                        auto nb_ins = c->_dim[0];
                        size_t idx= 0;
                        for (size_t i = 0; i< nb_ins; i++){
                            if (!*c->_all_lazy || !c->_lazy[i]) {
                                cid = c->_id+idx++;
                                g_l[cid] = 0;
                                g_u[cid] = 0;
                            }
                        }
                        break;
                    }
                    case leq:{
                        auto nb_ins = c->_dim[0];
                        size_t idx= 0;
                        for (size_t i = 0; i< nb_ins; i++){
                            if (!*c->_all_lazy || !c->_lazy[i]) {
                                cid = c->_id+idx++;
                                g_l[cid] = numeric_limits<double>::lowest();
                                g_u[cid] = 0;
                            }
                        }
                        break;
                    }
                    case geq:{
                        auto nb_ins = c->_dim[0];
                        size_t idx= 0;
                        for (size_t i = 0; i< nb_ins; i++){
                            if (!*c->_all_lazy || !c->_lazy[i]) {
                                cid = c->_id+idx++;
                                g_l[cid] = 0;
                                g_u[cid] = numeric_limits<double>::max();
                            }
                        }
                        break;
                    }
                    default:
                        throw invalid_argument("Undefined constraint type!\n");
                        exit(-1);
                        break;
                }
            }
        }
        
        
        void embed(expr<type>& e){/**<  Transfer all variables and parameters to the model. */
            switch (e.get_type()) {
                case uexp_c:{
                    auto ue = (uexpr<type>*)&e;
                    if (ue->_son->is_function()) {
                        auto f = static_pointer_cast<func<type>>(ue->_son);
                        bool found_cpy = false;
                        auto name = f->to_str();
                        if (name.back()=='T') {
                            name = name.substr(0,name.size()-2);
                            if (_nl_funcs_map.count(name)>0) {
                                auto cpy = _nl_funcs_map.at(name);
                                f->_val = cpy->_val;
                                f->_evaluated = true;
                                found_cpy = true;
                            }
                        }
                        else {
                            auto name1 = "[["+name+"]]\u1D40";
                            if (_nl_funcs_map.count(name1)>0) {
                                auto cpy = _nl_funcs_map.at(name1);
                                f->_val = cpy->_val;
                                f->_evaluated = true;
                                found_cpy = true;
                            }
                            auto name2 = name+"\u1D40";
                            if (_nl_funcs_map.count(name2)>0) {
                                auto cpy = _nl_funcs_map.at(name2);
                                f->_val = cpy->_val;
                                f->_evaluated = true;
                                found_cpy = true;
                            }
                        }
                        if (!found_cpy) {
                            auto f_p = _nl_funcs_map.insert(make_pair<>(f->to_str(), f));
                            if (f_p.second) {
                                embed(f);
                                _nl_funcs.push_back(f);
                                DebugOff(f->to_str() << endl);
                                //                f->_val = make_shared<vector<double>>();
                                //                f->_val->resize(f->_dim[0]);
                            }
                            else {
                                ue->_son = f_p.first->second;
                            }
                        }
                    }
                    break;
                }
                case bexp_c:{
                    auto be = (bexpr<type>*)&e;
                    if (be->_lson->is_function()) {
                        auto f = static_pointer_cast<func<type>>(be->_lson);
                        bool found_cpy = false;
                        auto name = f->to_str();
                        if (name.back()=='T') {
                            name = name.substr(0,name.size()-2);
                            if (_nl_funcs_map.count(name)>0) {
                                auto cpy = _nl_funcs_map.at(name);
                                f->_val = cpy->_val;
                                f->_evaluated = true;
                                found_cpy = true;
                            }
                        }
                        else {
                            auto name1 = "[["+name+"]]\u1D40";
                            if (_nl_funcs_map.count(name1)>0) {
                                auto cpy = _nl_funcs_map.at(name1);
                                f->_val = cpy->_val;
                                f->_evaluated = true;
                                found_cpy = true;
                            }
                            auto name2 = name+"\u1D40";
                            if (_nl_funcs_map.count(name2)>0) {
                                auto cpy = _nl_funcs_map.at(name2);
                                f->_val = cpy->_val;
                                f->_evaluated = true;
                                found_cpy = true;
                            }
                        }
                        if (!found_cpy) {
                            auto f_p = _nl_funcs_map.insert(make_pair<>(f->to_str(), f));
                            if (f_p.second) {
                                embed(f);
                                DebugOff(f->to_str() << endl);
                                _nl_funcs.push_back(f);
                                //                f->_val = make_shared<vector<double>>();
                                //                f->_val->resize(f->_dim[0]);
                            }
                            else {
                                be->_lson = f_p.first->second;
                            }
                        }
                    }
                    if (be->_rson->is_function()) {
                        auto f = static_pointer_cast<func<type>>(be->_rson);
                        auto found_cpy = false;
                        auto name = f->to_str();
                        if (name.back()=='T') {
                            name = name.substr(0,name.size()-2);
                            if (_nl_funcs_map.count(name)>0) {
                                auto cpy = _nl_funcs_map.at(name);
                                f->_val = cpy->_val;
                                f->_evaluated = true;
                                found_cpy = true;
                            }
                        }
                        else {
                            auto name1 = "[["+name+"]]\u1D40";
                            if (_nl_funcs_map.count(name1)>0) {
                                auto cpy = _nl_funcs_map.at(name1);
                                f->_val = cpy->_val;
                                f->_evaluated = true;
                                found_cpy = true;
                            }
                            auto name2 = name+"\u1D40";
                            if (_nl_funcs_map.count(name2)>0) {
                                auto cpy = _nl_funcs_map.at(name2);
                                f->_val = cpy->_val;
                                f->_evaluated = true;
                                found_cpy = true;
                            }
                        }
                        if (!found_cpy) {
                            auto f_p = _nl_funcs_map.insert(make_pair<>(f->to_str(), f));
                            if (f_p.second) {
                                embed(f);
                                DebugOff(f->to_str() << endl);
                                _nl_funcs.push_back(f);
                                //                f->_val = make_shared<vector<double>>();
                                //                f->_val->resize(f->_dim[0]);
                            }
                            else {
                                be->_rson = f_p.first->second;
                            }
                        }
                    }
                    break;
                }
                default:
                    break;
            }
        }
        
        
        shared_ptr<func<type>> embed(shared_ptr<func<type>> f){/**<  Transfer all variables and parameters to the model. */
            DebugOff(f->to_str() << endl);
            for (auto &p_t: f->get_lterms()) {
                if (p_t.second._coef->is_function()) {
                    auto cf = static_pointer_cast<func<type>>(p_t.second._coef);
                    auto exp = cf->get_expr();
                    if (exp) {
                        embed(*exp);
                    }
                    if(cf->is_matrix()){
                        auto newf = embed(make_shared<func<type>>(*cf));
                        cf->_val = newf->_val;
                    }
                }
            }
            for (auto &p_t: f->get_qterms()) {
                if (p_t.second._coef->is_function()) {
                    auto cf = static_pointer_cast<func<type>>(p_t.second._coef);
                    auto exp = cf->get_expr();
                    if (exp) {
                        embed(*exp);
                    }
                    if(cf->is_matrix()){
                        auto newf = embed(make_shared<func<type>>(*cf));
                        cf->_val = newf->_val;
                    }
                }
            }
            for (auto &p_t: f->get_pterms()) {
                if (p_t.second._coef->is_function()) {
                    auto cf = static_pointer_cast<func<type>>(p_t.second._coef);
                    auto exp = cf->get_expr();
                    if (exp) {
                        embed(*exp);
                    }
                }
            }
            if (f->get_cst()->is_function()) {
                auto c = static_pointer_cast<func<type>>(f->get_cst());
                auto exp = c->get_expr();
                if (exp) {
                    embed(*exp);
                }
                if(c->is_matrix()){
                    auto newf = embed(make_shared<func<type>>(*c));
                    c->_val = newf->_val;
                }
            }
            if (f->get_expr()) {
                embed(*f->get_expr());
            }
            bool found_cpy = false;
            auto name = f->to_str();
            if (name.back()=='T') {
                name = name.substr(0,name.size()-2);
                if (_nl_funcs_map.count(name)>0) {
                    auto cpy = _nl_funcs_map.at(name);
                    f->_val = cpy->_val;
                    f->_evaluated = true;
                    found_cpy = true;
                }
            }
            else {
                auto name1 = "[["+name+"]]\u1D40";
                if (_nl_funcs_map.count(name1)>0) {
                    auto cpy = _nl_funcs_map.at(name1);
                    f->_val = cpy->_val;
                    f->_evaluated = true;
                    found_cpy = true;
                }
                auto name2 = name+"\u1D40";
                if (_nl_funcs_map.count(name2)>0) {
                    auto cpy = _nl_funcs_map.at(name2);
                    f->_val = cpy->_val;
                    f->_evaluated = true;
                    found_cpy = true;
                }
            }
            if (!found_cpy) {
                auto f_p = _nl_funcs_map.insert(make_pair<>(f->to_str(), f));
                if (f_p.second) {
                    _nl_funcs.push_back(f_p.first->second);
                    DebugOff(f->to_str() << endl);
                    f->allocate_mem();
                    return f;
                    //        f_p.first->second->_val = make_shared<vector<double>>();
                    //        f_p.first->second->_val->resize(f_p.first->second->_dim[0]);
                }
                //        if (f->_new) {
                //            f_p.first->second = f;
                //            return f;
                //        }
                if (f->_dim[0] > f_p.first->second->_dim[0]) {
                    *f_p.first->second = *f;
                }
                else if (f->_dfdx->size()>0) {
                    *f_p.first->second = *f;
                }
                //        f_p.first->second->allocate_mem();
                return f_p.first->second;
            }
            return f;
        }
        
        
        
        void print_nl_functions() const{
            cout << "Number of atomic functions = " << _nl_funcs.size();
            cout << endl;
            for (auto& f: _nl_funcs){
                f->print_symbolic();
                f->print();
                cout << endl;
            }
            cout << endl;
        }
        
        
        void print_solution(int prec=5) const{
            for (auto &v_pair:_vars) {
                auto v = v_pair.second;
                v->print(true,prec);
            }
        }
        
        shared_ptr<param_> get_int_var(size_t idx){
            return _int_vars.at(idx);
        }
        
        void round_solution(){
            for (auto &v_pair:_vars) {
                if(v_pair.second->_is_relaxed){
                    v_pair.second->round_vals();
                    auto int_var = get_int_var(v_pair.first);
                    int_var->copy_vals(*v_pair.second);
                }
            }
        }
        
        
        void print_symbolic(){
            cout << "-------------------------" << endl;
            print_properties();
            _obj->print_symbolic();
            for(auto& p: _cons){
                p.second->print_symbolic();
            }
            for(auto& v: _vars){
                v.second->print_symbolic();
            }
            cout << "-------------------------" << endl;
        }
        
        size_t print_properties() const{
            string str = "\n";
            if(is_linear()){
                str += "Linear ";
            }
            else if(is_convex()){
                str += "Convex ";
            }
            else if(is_concave()){
                str += "Concave ";
            }
            else {
                str += "Nonconvex ";
            }
            str += "Model: " + _name+"\n";;
            auto size_header = str.size()-1;
            str.insert(0,size_header,'-');
            str.append(size_header,'-');
            cout << str << endl;
            cout << "Number of variables = " << get_nb_vars() << endl;
            cout << "Number of constraints = " << get_nb_cons() << " (" << get_nb_ineq() << " inequalities, " << get_nb_eq() << " equations)" << endl;
            //    compute_funcs();
            cout << "Objective: ";
            if(_objt==minimize){
                cout << "Min ";
            }
            else{
                cout << "Max ";
            }
            return size_header;
        }
        
        void print(){
            auto size_header = print_properties();
            _obj->print();
            cout << "s.t." << endl;
            for(auto& p: _cons){
                p.second->print();
            }
            for(auto& v: _vars){
                v.second->print();
            }
            string tail;
            tail.insert(0,size_header,'-');
            cout << tail << endl;
        }
        
        
        
        
        void print_constraints() const{
            for(auto& p: _cons){
                p.second->print_symbolic();
            }
        }
        
        
        
        
        void replace_integers(){/*< Replace internal type of integer variables so that continuous relaxations can be computed */
            bool has_int = false;
            for (auto &v_p:this->_vars_name) {
                if (v_p.second->is_integer() || v_p.second->is_binary()) {
                    has_int = true;
                    auto v = v_p.second;
                    auto new_v = make_shared<var<double>>(v_p.second->_name);
                    new_v->shallow_copy(*v);
                    new_v->_is_relaxed = true;                    
                    new_v->copy_vals(*v);
                    new_v->copy_bounds(*v);
                    v_p.second = new_v;
                }
            }
            for (auto &v_p:this->_vars) {
                if (v_p.second->is_integer() || v_p.second->is_binary()) {
                    auto name = v_p.second->_name;
                    v_p.second = this->get_var_ptr(name);
                }
            }
            if(has_int){
                this->_obj->relax(this->_vars);
                for (auto &c_p: this->_cons) {
                    c_p.second->relax(this->_vars);
                }
            }
        }
        
        
        void add_on_off(const Constraint<type>& c, var<bool>& on){
            if (c.get_ftype() != lin_) {
                cerr << "Nonlinear constraint.\n";
                exit(-1);
            }
            Constraint<type> res(c.get_name() + "_on/off");
            //    double b;
            //    for(auto it: orig_q->_coefs) {
            //        v = getparam_<double>(it.first);
            //        if (!v->is_bounded_below() || !v->is_bounded_above()) {
            //            cerr << "Variable " << v->_name << " in constraint " << c._name << " does not have finite bounds.\n";
            //            exit(1);
            //        }
            //        if (c.get_type() == leq || c.get_type() == eq) {
            //            if (it.second < 0) res -= it.second*v->get_lb_off()*(1-on);
            //            else res -= it.second*v->get_ub_off()*(1-on);
            //        }
            //        else{ // geq
            //            if (it.second < 0) res -= it.second*v->get_ub_off()*(1-on);
            //            else res -= it.second*v->get_lb_off()*(1-on);
            //        }
            //    }
            //    if (c.get_type() == eq) {
            //        Constraint res2(c.get_name() + "_on/off2");
            //        for(auto it: orig_q->_coefs) {
            //            v = getparam_<double>(it.first);
            //            if (it.second < 0) res2 -= it.second*v->get_ub_off()*(1-on);
            //            else res2 -= it.second*v->get_lb_off()*(1-on);
            //        }
            //        res2 += *orig_q;
            //        res2 -= b*on;
            //        res2 >= 0;
            //        addConstraint(res2);
            //    }
            //    res += *orig_q;
            //    res -= orig_q->get_const();
            //    res -= b*on;
            //    if (c.get_type() == eq or c.get_type() == leq) res <= 0;
            //    else res >= 0;
            add_constraint(res);
        }
        
        
        void add_on_off(var<>& v, var<bool>& on){
            //    if(v.get_ub() != v.get_ub_off()) {
            //        Constraint UB(v._name + "_UB_on/off");
            //        UB += v - v.get_ub() * on - (1 - on) * v.get_ub_off();
            //        UB <= 0;
            //        addConstraint(UB);
            //    }
            //    if(v.get_lb() != v.get_lb_off()) {
            //        Constraint LB(v._name + "_LB_on/off");
            //        LB += v - v.get_lb() * on - (1 - on) * v.get_lb_off();
            //        LB >= 0;
            //        addConstraint(LB);
            //    }
        }
        
        
        
        void add_on_off_McCormick(std::string name, var<>& v, var<>& v1, var<>& v2, var<bool>& on) {
            //    Constraint MC1(name+"_McCormick1");
            //    MC1 += v;
            //    MC1 -= v1.get_lb()*v2 + v2.get_lb()*v1 - v1.get_lb()*v2.get_lb();
            //    MC1 >= 0;
            //    add_on_off(MC1, on);
            //    Constraint MC2(name+"_McCormick2");
            //    MC2 += v;
            //    MC2 -= v1.get_ub()*v2 + v2.get_ub()*v1 - v1.get_ub()*v2.get_ub();
            //    MC2 >= 0;
            //    add_on_off(MC2, on);
            //    Constraint MC3(name+"_McCormick3");
            //    MC3 += v;
            //    MC3 -= v1.get_lb()*v2 + v2.get_ub()*v1 - v1.get_lb()*v2.get_ub();
            //    MC3 <= 0;
            //    add_on_off(MC3, on);
            //    Constraint MC4(name+"_McCormick4");
            //    MC4 += v;
            //    MC4 -= v1.get_ub()*v2 + v2.get_lb()*v1 - v1.get_ub()*v2.get_lb();
            //    MC4 <= 0;
            //    add_on_off(MC4, on);
        }
        
        
        

        
    };
    
//    void compute_constrs(vector<Constraint*>& v, double* res, unsigned i, unsigned j);

    template<typename T>
    pair<shared_ptr<func_>, ObjectiveType> max(const func<T>& f){
        auto fcpy = f.copy();
        f->_val->resize(1);
        return make_pair<>(f,maximize);
    };
    
    template<typename T>
    pair<shared_ptr<func_>, ObjectiveType> min(const func<T>& f){
        auto fcpy = f.copy();
        f->_val->resize(1);
        return make_pair<>(f,minimize);
    };

    template<typename type = double>
    class Program{
    public:
        virtual void update_model(){};
        string _status;
    };
    
    
}



#endif /* model_hpp */
