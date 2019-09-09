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
#include <bitset>
#include <gravity/constraint.h>
#include <map>
#include <unordered_set>
#include <math.h>
#include <vector>
#include <deque>
#include <thread>
#include <iostream>
#include <functional>
#ifdef USE_MPI
#include <mpi.h>
#endif
#ifdef USE_IPOPT
#define HAVE_STDDEF_H
#include <coin/IpIpoptApplication.hpp>
#include <coin/IpTNLP.hpp>
#undef HAVE_STDDEF_H
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
    template<typename type,typename std::enable_if<is_arithmetic<type>::value>::type* = nullptr>
    void compute_constrs(vector<shared_ptr<Constraint<type>>>& v, double* res, size_t i, size_t j){
        DebugOff("Calling compute_constrts with i =  " << i << "and j = "<< j << endl);
        for (size_t idx = i; idx < j; idx++) {
            auto c = v[idx];
            //            c->print_symbolic();
            //            if(c->_name == "Real(Linking_V_mag)_lifted"){
            //                cout << "ok";
            //            }
            c->_new = false;
            c->_evaluated = false;
            size_t nb_ins = c->get_nb_inst();
            size_t ind = 0;
            for (size_t inst = 0; inst< nb_ins; inst++){
                if (!*c->_all_lazy || !c->_lazy[inst]) {
                    res[c->_id+ind++] = c->eval(inst);
                    DebugOff("Accessing res at position " << c->_id+inst << endl);
                    //                _cons_vals[index++] = res[c->_id+inst];
                    DebugOff("g[" << to_string(c->_id+inst) << "] = " << to_string(res[c->_id+inst]) << endl);
                    //                    if(c->_id+inst==15){
                    //                        cout << "ok";
                    //                    }
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
    template<typename type,typename std::enable_if<is_arithmetic<type>::value>::type* = nullptr>
    void compute_jac(vector<shared_ptr<Constraint<type>>>& vec, double* res, size_t i, size_t j, bool first_call, vector<double>& jac_vals){
        size_t cid = 0, id_inst = 0;
        string vid;
        shared_ptr<Constraint<type>> c = NULL;
        param_* v = NULL;
        shared_ptr<func<type>> dfdx;
        auto idx = vec[i]->_jac_cstr_idx;
        for (size_t ind = i; ind < j; ind++) {
            c = vec[ind];
            auto nb_ins = c->get_nb_inst();
            id_inst = 0;
            if (c->is_linear() && !first_call) {
                //                        if (false) {
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
                    if(!dfdx->is_constant()){
                        dfdx->_evaluated=false;
                    }
                    id_inst = 0;
                    for (size_t inst = 0; inst< nb_ins; inst++){
                        if (!*c->_all_lazy || !c->_lazy[inst]) {
                            cid = c->_id+id_inst++;
                            
                            if (v->_is_vector || v->is_matrix_indexed()) {
                                auto dim = v->get_dim(inst);
                                for (size_t j = 0; j<dim; j++) {
                                    res[idx] += dfdx->eval(inst,j);
                                    jac_vals[idx] = res[idx];
                                    DebugOff("jac_val["<< idx <<"] = " << jac_vals[idx] << endl);
                                    idx++;
                                }
                            }
                            else {
                                res[idx] += dfdx->eval(inst);
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
        
        bool                            _has_callback = false; /**< Has callback option. */
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
        map<pair<size_t,size_t>, size_t>                    _nnz_pairs;
        vector<size_t>                                      _idx_it;
        map<size_t, shared_ptr<param_>>                     _params; /**< Sorted map pointing to all parameters contained in this model. */
        map<size_t, shared_ptr<param_>>                     _vars; /**< Sorted map pointing to all variables contained in this model. */
        map<size_t, shared_ptr<param_>>                     _int_vars; /**< Sorted map pointing to all binary variables contained in this model. */
        map<string, shared_ptr<param_>>                     _params_name; /**< Sorted map (by name) pointing to all parameters contained in this model. */
        map<string, shared_ptr<param_>>                     _vars_name; /**< Sorted map (by name) pointing to all variables contained in this model. */
        vector<shared_ptr<Constraint<type>>>                _cons_vec; /**< vector pointing to all constraints contained in this model. */
        map<size_t, shared_ptr<Constraint<type>>>           _cons; /**< Sorted map (increasing index) pointing to all constraints contained in this model. */
        map<string, shared_ptr<Constraint<type>>>           _cons_name; /**< Sorted map (by name) pointing to all constraints contained in this model. */
        map<string, set<shared_ptr<Constraint<type>>>>   _v_in_cons; /**< Set of constraints where each variable appears. */
        shared_ptr<func<type>>                              _obj = nullptr; /**< Pointer to objective function */
        ObjectiveType                                       _objt = minimize; /**< Minimize or maximize */
        int                                                 _status = -1;/**< status when last solved */
        map<pair<string, string>,map<int,pair<shared_ptr<func<type>>,shared_ptr<func<type>>>>>            _hess_link; /* for each pair of variables appearing in the hessian, storing the set of constraints they appear together in */
        
        void merge_vars(const shared_ptr<expr<type>>& e){/**<  Transfer all variables and parameters to the model. */
            switch (e->get_type()) {
                case uexp_c:{
                    auto ue = (uexpr<type>*)e.get();
                    if (ue->_son->is_function()) {
                        auto f = static_pointer_cast<func<type>>(ue->_son);
                        merge_vars(f);
                    }
                    break;
                }
                case bexp_c:{
                    auto be = (bexpr<type>*)e.get();
                    if (be->_lson->is_function()) {
                        auto f = static_pointer_cast<func<type>>(be->_lson);
                        merge_vars(f);
                    }
                    if (be->_rson->is_function()) {
                        auto f = static_pointer_cast<func<type>>(be->_rson);
                        merge_vars(f);
                    }
                    break;
                }
                default:
                    break;
            }
        }
        
        /**
         Subfunction of embed(func_&& f). Merge variables and parameters with f. If a variable x in f exists in the current funtion, x will now point to the same variable appearing in current function.
         @param[in] f function to merge variables and parameters with.
         */
        void merge_vars(const shared_ptr<func<type>>& f){
            for (auto &pair:*f->_lterms) {
                auto p = pair.second._p;
                if (p->is_var()) {
                    auto pid = *p->_vec_id;
                    p->share_vals(_vars.at(pid));
                }
            }
            for (auto &pair:*f->_qterms) {
                auto coef = pair.second._coef;
                auto p1 = pair.second._p->first;
                auto p2 = pair.second._p->second;
                if (p1->is_var()) {
                    auto pid1 = *p1->_vec_id;
                    p1->share_vals(_vars.at(pid1));
                }
                if (p2->is_var()) {
                    auto pid2 = *p2->_vec_id;
                    p2->share_vals(_vars.at(pid2));
                }
            }
            for (auto &pair:*f->_pterms) {
                auto list = pair.second._l;
                for (auto &ppi: *list) {
                    auto p = ppi.first;
                    if (p->is_var()) {
                        auto pid = *p->_vec_id;
                        ppi.first->share_vals(_vars.at(pid));
                    }
                }
            }
            if (f->_expr) {
                merge_vars(f->_expr);
            }
        }
        
        shared_ptr<Model<type>> copy() const{
            shared_ptr<Model<type>> cpy = make_shared<Model<type>>();
            cpy->_name = _name;
            for(auto &vp: _vars){
                switch (vp.second->get_intype()) {
                    case binary_: {
                        auto vv = *static_pointer_cast<var<bool>>(vp.second);
                        cpy->add(vv.deep_copy());
                        break;
                    }
                    case short_: {
                        auto vv = *static_pointer_cast<var<short>>(vp.second);
                        cpy->add(vv.deep_copy());
                        break;
                    }
                    case integer_: {
                        auto vv = *static_pointer_cast<var<int>>(vp.second);
                        cpy->add(vv.deep_copy());
                        break;
                    }
                    case float_: {
                        auto vv = *static_pointer_cast<var<float>>(vp.second);
                        cpy->add(vv.deep_copy());
                        break;
                    }
                    case double_: {
                        auto vv = *static_pointer_cast<var<double>>(vp.second);
                        cpy->add(vv.deep_copy());
                        break;
                    }
                    case long_: {
                        auto vv = *static_pointer_cast<var<long double>>(vp.second);
                        cpy->add(vv.deep_copy());
                        break;
                    }
                    case complex_: {
                        auto vv = *static_pointer_cast<var<Cpx>>(vp.second);
                        cpy->add(vv.deep_copy());
                        break;
                    }
                }
            }
            for(auto &cp: _cons_name){
                cpy->add(*cp.second);
                cpy->merge_vars(cpy->_cons_vec.back());
                cpy->_cons_vec.back()->uneval();
            }
            cpy->set_objective(*_obj, _objt);
            return cpy;
        }
        
        //        Model& operator=(const Model& m){
        //            _name = m._name;
        //            _hess = m._hess;
        //            _nl_funcs =m._nl_funcs;
        //            _nl_funcs_map = m._nl_funcs_map;
        //            _has_lazy = m._has_lazy;
        //            _built = m._built;
        //            _first_run = m._first_run;
        //            _first_call_gard_obj = m._first_call_gard_obj;
        //            _first_call_jac = m._first_call_jac;
        //            _first_call_hess = m._first_call_hess;
        //            _convexity = m._convexity;
        //            _type = m._type;
        //            _nb_vars = m._nb_vars;
        //            _nb_params = m._nb_params;
        //            _nb_cons = m._nb_cons;
        //            _nnz_g = m._nnz_g;
        //            _nnz_h = m._nnz_h;
        //            _nnz_g_obj = m._nnz_g_obj;
        //            _jac_vals = m._jac_vals;
        //            _obj_grad_vals = m._obj_grad_vals;
        //            _hess_vals = m._hess_vals;
        //            _params = m._params;
        //            _vars = m._vars;
        //            _int_vars = m._int_vars;
        //            _params_name = m._params_name;
        //            _vars_name = m._vars_name;
        //            _cons_vec = m._cons_vec;
        //            _cons = m._cons;
        //            _cons_name = m._cons_name;
        //            _v_in_cons = m._v_in_cons;
        //            _obj = m._obj->copy();
        //            _objt = m._objt;
        //            _status = m._status;
        //            _hess_link = m._hess_link;
        //            return *this;
        //        }
        
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
        
        void set_name(const string& name){
            _name = name;
        }
        
        void add_var(const shared_ptr<param_>& v){
            switch (v->get_intype()) {
                case binary_:
                    add(*static_pointer_cast<var<bool>>(v));
                    break;
                case short_:
                    add(*static_pointer_cast<var<short>>(v));
                    break;
                case integer_:
                    add(*static_pointer_cast<var<int>>(v));
                    break;
                case float_:
                    add(*static_pointer_cast<var<float>>(v));
                    break;
                case double_:
                    add(*static_pointer_cast<var<double>>(v));
                    break;
                case long_:
                    add(*static_pointer_cast<var<long double>>(v));
                    break;
                case complex_:
                    add(*static_pointer_cast<var<Cpx>>(v));
                    break;
                default:
                    break;
            }
        }
        
        template <typename T>
        void add_var(var<T>& v){//Add variables by copy
            auto name = v._name.substr(0,v._name.find_first_of("."));
            //            auto name = v._name;
            v._name = name;
            v._lb->eval_all();
            v._ub->eval_all();
            if (_vars_name.count(v._name)==0) {
                v.set_id(_nb_vars);
                v.set_vec_id(_vars.size());
                shared_ptr<param_> newv;
                if (v._val->empty()) {
                    Warning("WARNING adding unindexed variable to model, treating it as a one dimensional Real.\n");
                    newv = (v.in(R(1))).pcopy();
                }
                else {
                    newv = v.pcopy();
                }
                _vars_name[name] = newv;
                _vars[v.get_vec_id()] = newv;
                _nb_vars += v.get_dim();
            }
            else {
                throw invalid_argument("adding variable with same name, please rename: " + v._name);
            }
        };
        
        
        template <typename T>
        void add_var(var<T>&& v){//Add variables by copy
            if(v.get_dim()==0)
                return;
            auto name = v._name.substr(0,v._name.find_first_of("."));
            //            auto name = v._name;
            v._name = name;
            v._lb->eval_all();
            v._ub->eval_all();
            if (_vars_name.count(v._name)==0) {
                v.set_id(_nb_vars);
                v.set_vec_id(_vars.size());
                shared_ptr<param_> newv;
                if (v._val->empty()) {
                    Warning("WARNING adding unindexed variable to model, treating it as a one dimensional Real.\n");
                    newv = make_shared<var<T>>(move((v.in(R(1)))));
                }
                else {
                    newv = make_shared<var<T>>(move(v));
                }
                _vars_name[name] = newv;
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
            list<var<T>> vars;
            vars = {forward<var<T>>(args)...};
            vars.push_front(move(v));
            for (auto &v:vars) {
                add_var(move(v));
            }
        }
        
        
        /* Output */
        
        
        /* Accessors */
        
        void update_nb_vars(){
            size_t n = 0;
            for (auto &vp:_vars) {
                n += vp.second->get_dim();
            }
            _nb_vars = n;
        };
        
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
                auto nb_inst = c->get_nb_inst();
                for (size_t inst = 0; inst<nb_inst; inst++) {
                    if (!*c->_all_lazy || !c->_lazy[inst]) {
                        _nnz_g += c->get_nb_vars(inst);
                    }
                }
            }
            return _nnz_g;
        };
        
        template <typename T>
        var<T> get_var(const string& vname) const{
            auto it = _vars_name.find(vname);
            if (it==_vars_name.end()) {
                throw invalid_argument("In function: Model::get_var(const string& vname) const, unable to find variable with given name");
            }
            auto v = dynamic_pointer_cast<var<T>>(it->second);
            if(v){
                return *v;
            }
            throw invalid_argument("In function: Model::get_var<T>(const string& vname) const, cannot cast variable, make sure to use the right numerical type T");
        }
        
        /* Return the number of nonzeros in the lower left part of the hessian */
        
        size_t get_nb_nnz_h(){
            size_t idx = 0, vid, vjd;
            string vi_name, vj_name;
            shared_ptr<param_> vi;
            shared_ptr<param_> vj;
            shared_ptr<Constraint<type>> c;
            for (auto &pairs: _hess_link) {
                vi_name = pairs.first.first;
                vj_name = pairs.first.second;
                vi = (pairs.second.begin())->second.first->get_var(vi_name);
                vj = (pairs.second.begin())->second.first->get_var(vj_name);
                if (vi_name.compare(vj_name) > 0) {//ONLY STORE LOWER TRIANGULAR PART OF HESSIAN
                    throw invalid_argument("SHOULD BE SORTED CORRECTLY IN FILL_MAPS");
                }
                vid = vi->get_id();
                vjd = vj->get_id();
                
                
                
                for (auto &f_pair:pairs.second) {
                    auto f = f_pair.second.first;
                    if (f->_is_constraint) {
                        c = static_pointer_cast<Constraint<type>>(f);
                    }
                    auto d2f = f_pair.second.second;
                    size_t nb_inst = f->get_nb_inst();
                    for (size_t inst = 0; inst<nb_inst; inst++) {
                        if (!(f->_is_constraint && *c->_all_lazy && c->_lazy[inst])) {
                            if(d2f->is_matrix_indexed()){
                                auto dim = d2f->get_dim(inst);
                                for (size_t j = 0; j<dim; j++) {
                                    if(_nnz_pairs.insert({{vid + vi->get_id_inst(inst,j),vjd+vj->get_id_inst(inst,j)}, idx}).second){
                                        idx++;
                                    };
                                    _idx_it.push_back(_nnz_pairs.at({vid + vi->get_id_inst(inst,j),vjd+vj->get_id_inst(inst,j)}));
                                }
                            }
                            else if (d2f->is_matrix()) {
                                for (size_t i = 0; i < d2f->get_nb_inst(); i++) {
                                    for (size_t j = i; j < d2f->_dim[1]; j++) {
                                        if(_nnz_pairs.insert({{vid + vi->get_id_inst(i),vjd+vj->get_id_inst(j)}, idx}).second){
                                            idx++;
                                        };
                                        _idx_it.push_back(_nnz_pairs.at({vid + vi->get_id_inst(i),vjd+vj->get_id_inst(j)}));
                                    }
                                }
                            }
                            else if(d2f->_is_vector){
                                for (size_t j = 0; j < d2f->get_nb_inst(); j++) {
                                    if(_nnz_pairs.insert({{vid + vi->get_id_inst(j),vjd+vj->get_id_inst(j)}, idx}).second){
                                        idx++;
                                    };
                                    _idx_it.push_back(_nnz_pairs.at({vid + vi->get_id_inst(j),vjd+vj->get_id_inst(j)}));
                                }
                            }
                            else {
                                if(_nnz_pairs.insert({{vid + vi->get_id_inst(inst),vjd+vj->get_id_inst(inst)}, idx}).second){
                                    idx++;
                                };
                                _idx_it.push_back(_nnz_pairs.at({vid + vi->get_id_inst(inst),vjd+vj->get_id_inst(inst)}));
                            }
                        }
                    }
                }
            }
            _nnz_h = _nnz_pairs.size();
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
            _idx_it.clear();
            _cons_vec.clear();
            _nnz_pairs.clear();
            _jac_vals.clear();
            _hess_vals.clear();
            _first_call_jac = true;
            _first_call_hess = true;
            _first_call_gard_obj = true;
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
                _cons_vec.push_back(c);
                new_cid = c->_id+nb_inst;
            }
            _cons = new_cons;
            _nb_cons = get_nb_cons();
        }
        
        
        void reindex_vars(){/**< Re-index all variables involved in the model and update total number of variables */
            shared_ptr<param_> v= nullptr;
            size_t idx = 0, vec_idx = 0;
            for(auto& v_p: _vars)
            {
                v = v_p.second;
                v->set_vec_id(vec_idx++);
                v->set_id(idx);
                idx += v->get_dim();
            }
            _nb_vars = idx;
        }
        
        
        void del_var(const param_& v){
            auto it = _vars.find(v.get_id());
            if (it!=_vars.end()) {
                _nb_vars -= v.get_dim();
                _vars.erase(it);
            }
            reindex_vars();
        };
        
        //generic add constraint
        //INPUT: a constraint, lifting option boolean, disjunctive union methods from ("on/off", "lambda_II", "lambda_III")
        //OUTPUT: addition of this constraint in the model, if lift is selected, the convex relaxation and the partitioning of the individual variables will be also included to the mathematical formulation
        template<typename T=type,typename std::enable_if<is_arithmetic<T>::value>::type* = nullptr>
        void add(const Constraint<Cpx>& c, bool convexify = false, string method_type = "on/off"){
            if (c.get_dim()==0) {
                return;
            }
            auto real_imag = get_real_imag(c);
            Constraint<type> c_real;
            c_real += real_imag.first;
            c_real._name = "Real(" + c._name + ")";
            c_real._ctype = c._ctype;
            c_real._indices = c._indices;
            c_real._dim[0] = c._dim[0];
            Constraint<type> c_imag;
            c_imag += real_imag.second;
            c_imag._name = "Imag(" + c._name + ")";
            c_imag._ctype = c._ctype;
            c_imag._indices = c._indices;
            c_imag._dim[0] = c._dim[0];
            if(convexify){
                c_real.check_soc();c_real.check_rotated_soc();
                auto lifted_real = lift(c_real, method_type);
                c_imag.check_soc();c_imag.check_rotated_soc();
                auto lifted_imag = lift(c_imag, method_type);
                lifted_real._ctype = c._ctype;
                lifted_real._indices = c._indices;
                lifted_real._dim[0] = c._dim[0];
                lifted_imag._ctype = c._ctype;
                lifted_imag._indices = c._indices;
                lifted_imag._dim[0] = c._dim[0];
                add_constraint(lifted_real);
                add_constraint(lifted_imag);
                if(c_real.func<type>::is_convex() && c_real._ctype==eq){
                    DebugOn("Convex left hand side of equation detected, splitting constraint into <= and ==" << endl);
                    Constraint<type> c_real_cvx(c_real._name+"_convex");
                    c_real_cvx = c_real;
                    c_real_cvx._relaxed = true;
                    add_constraint(c_real_cvx <= 0);
                }
                if(c_real.func<type>::is_concave() && c_real._ctype==eq){
                    DebugOn("Concave left hand side of equation detected, splitting constraint into >= and ==" << endl);
                    Constraint<type> c_real_ccve(c_real._name+"_concave");
                    c_real_ccve = c_real;
                    c_real_ccve._relaxed = true;
                    add_constraint(c_real_ccve >= 0);
                }
                if(c_imag.func<type>::is_convex() && c_real._ctype==eq){
                    DebugOn("Convex left hand side of equation detected, splitting constraint into <= and ==" << endl);
                    Constraint<type> c_imag_cvx(c_imag._name+"_convex");
                    c_imag_cvx = c_imag;
                    c_imag_cvx._relaxed = true;
                    add_constraint(c_imag_cvx <= 0);
                }
                if(c_imag.func<type>::is_concave() && c_real._ctype==eq){
                    DebugOn("Concave left hand side of equation detected, splitting constraint into >= and ==" << endl);
                    Constraint<type> c_imag_ccve(c_imag._name+"_concave");
                    c_imag_ccve = c_imag;
                    c_imag_ccve._relaxed = true;
                    add_constraint(c_imag_ccve >= 0);
                }
            }
            else {
                add_constraint(c_real);
                add_constraint(c_imag);
            }
        }
        
        /** Reformulate rotated SOCs into standard SOCs and add linking constraints to the model.
         @param[in] c: constraint to reformulate
         @return the linearized constraint
         @note This function will add constraints linking the lifted variables to the original ones, if a variable's partition is greater than 1, it will also add the disjunctive constraints corresponding to the partitionning of the variables.
         **/
        template<typename T=type,typename std::enable_if<is_arithmetic<type>::value>::type* = nullptr>
        Constraint<type> get_standard_SOC(Constraint<type>& c){
            if(!c.is_rotated_soc()){
                return c;
            }
            Constraint<type> newc_standard(c._name+"_standard");
            if (!c.get_cst()->is_zero()) { //here we extract the constant part of the constraint
                if (c.get_cst()->is_number()) {
                    auto f_cst = static_pointer_cast<constant<type>>(c.get_cst());
                    newc_standard.add_cst(*f_cst);
                }
                else if (c.get_cst()->is_param()) {
                    auto f_cst = static_pointer_cast<param<type>>(c.get_cst());
                    newc_standard.add_cst(*f_cst);
                }
                else {
                    auto f_cst = static_pointer_cast<func<type>>(c.get_cst());
                    newc_standard.add_cst(*f_cst);
                }
                if (newc_standard._cst->is_function()) {
                    newc_standard.embed(*static_pointer_cast<func<type>>(newc_standard._cst));
                }
            }
            for (auto &pair:*c._lterms) { //here we extract the linear terms appearing in the constraint
                auto term = pair.second;
                newc_standard.insert(term);
            }
            //go over the quadratic terms (including the bilinears)
            //here we will introduce two auxiliary variables for the bilinear term
            for (auto &pair:*c._qterms) {
                if (pair.second._p->first!=pair.second._p->second) { //means it is bilinear term
                    qterm lt1;
                    qterm lt2;
                    lt1._sign = !(pair.second._sign);
                    lt2._sign = pair.second._sign;
                    if (pair.second._coef->is_function()) {
                        auto coef = *static_pointer_cast<func<type>>(pair.second._coef);
                        lt1._coef = func<type>(coef).copy();
                        lt2._coef = func<type>(coef).copy();
                    }
                    else if(pair.second._coef->is_param()) {
                        auto coef = *static_pointer_cast<param<type>>(pair.second._coef);
                        lt1._coef = param<type>(coef).copy();
                        lt2._coef = param<type>(coef).copy();
                    }
                    else if(pair.second._coef->is_number()) {
                        auto coef = *static_pointer_cast<constant<type>>(pair.second._coef);
                        lt1._coef = constant<type>(coef).copy();
                        lt2._coef = constant<type>(coef).copy();
                    }
                    
                    auto v1 = *static_pointer_cast<var<type>>(pair.second._p->first); //assign the pointers to the variables
                    auto v2 = *static_pointer_cast<var<type>>(pair.second._p->second);
                    
                    if((v1._name > v2._name)){    //get the variables in the alphabetical order
                        v2 = *static_pointer_cast<var<type>>(pair.second._p->first);
                        v1 = *static_pointer_cast<var<type>>(pair.second._p->second);
                    }
                    auto ids = combine(*v1._indices,*v2._indices); //get the combined index set
                    
                    // I am not sure that we need the following part or not!! ***************************************************************************************************************************************
                    auto unique_ids = ids.get_unique_keys(); /* In case of an indexed variable, keep the unique keys only */
                    auto v1_ids = *v1._indices;
                    auto v2_ids = *v2._indices;
                    if(unique_ids.size()!=ids.size()){/* If some keys are repeated, remove them from the refs of o1 and o2 */
                        auto keep_refs = ids.get_unique_refs();
                        v1_ids.filter_refs(keep_refs);
                        v2_ids.filter_refs(keep_refs);
                    }
                    // ***************************************************************************************************************************************
                    
                    param<type> lb1("lb1"), ub1("ub1");
                    lb1.in(unique_ids);ub1.in(unique_ids);
                    param<type> lb2("lb2"), ub2("ub2");
                    lb2.in(unique_ids);ub2.in(unique_ids);
                    
                    //get the bounds for the auxiliary variables y1 and y2
                    for (int i=0; i<unique_ids.size(); i++) {
                        size_t id1;
                        size_t id2;
                        if(v1_ids._ids == nullptr){
                            id1 = i;
                        }
                        else id1 = v1_ids._ids->at(0).at(i);
                        if(v2_ids._ids == nullptr){
                            id2 = i;
                        }
                        else id2 = v2_ids._ids->at(0).at(i);
                        auto key1 = v1_ids._keys->at(id1);
                        auto key2 = v2_ids._keys->at(id2);
                        
                        auto sum_b1 = v1.get_lb(key1) - v2.get_ub(key2);
                        auto sum_b2 = v1.get_ub(key1) - v2.get_lb(key2);
                        auto sum_b3 = v1.get_lb(key1) + v2.get_lb(key2);
                        auto sum_b4 = v1.get_ub(key1) + v2.get_ub(key2);
                        
                        lb1.set_val(key1+","+key2, sum_b1/2);
                        ub1.set_val(key1+","+key2, sum_b2/2);
                        
                        lb2.set_val(key1+","+key2, sum_b3/2);
                        ub2.set_val(key1+","+key2, sum_b4/2);
                    }
                    //define the names of auxiliary variables to check their prior existence in the model
                    string aux1_name = "aux1("+v1.get_name(true,true)+v2.get_name(true,true)+")";
                    string aux2_name = "aux2("+v1.get_name(true,true)+v2.get_name(true,true)+")";
                    
                    auto it1 = _vars_name.find(aux1_name);
                    auto it2 = _vars_name.find(aux2_name);
                    
                    if(it1==_vars_name.end()){ //means the variables have not been included in the model, so create new ones, and add them to the model
                        //define variables
                        var<type> y1(aux1_name, lb1, ub1);
                        add(y1.in(unique_ids));
                        y1._num_partns = v1._num_partns + v2._num_partns;
                        lt1._p = make_shared<gravity::pair< shared_ptr<param_>,shared_ptr<param_> >>(make_pair(make_shared<var<type>>(y1.in(ids)), make_shared<var<type>>(y1.in(ids))));
                        var<type> y2(aux2_name, lb2, ub2);
                        add(y2.in(unique_ids));
                        lt2._p = make_shared<gravity::pair< shared_ptr<param_>,shared_ptr<param_> >>(make_pair(make_shared<var<type>>(y2.in(ids)), make_shared<var<type>>(y2.in(ids))));
                        //add constraints
                        Constraint<type> link1(pair.first+"_link1");
                        link1 = y1.in(unique_ids) - (v1.in(v1_ids) - v2.in(v2_ids))/2;
                        add(link1.in(unique_ids) == 0);
                        
                        Constraint<type> link2(pair.first+"_link2");
                        link2 = y2.in(unique_ids) - (v1.in(v1_ids) + v2.in(v2_ids))/2;
                        add(link2.in(unique_ids) == 0);
                    }
                    else{ //if the variables are already included in the model before, access the correct pointer and add the new indices/linkage
                        //get variables
                        auto y1 = static_pointer_cast<var<type>>(it1->second);
                        auto added1 = y1->add_bounds(lb1,ub1);
                        y1->_num_partns = v1._num_partns + v2._num_partns;
                        lt1._p = make_shared<gravity::pair< shared_ptr<param_>,shared_ptr<param_> >>(make_pair(make_shared<var<type>>(y1->in(ids)), make_shared<var<type>>(y1->in(ids))));
                        if(!added1.empty()){
                            assert(v1._indices->size()==v2._indices->size());
                            if(added1.size()!=v1._indices->size()){/* If some keys are repeated, remove them from the refs of o1 and o2 */
                                auto keep_refs = ids.diff_refs(added1);
                                v1_ids.filter_refs(keep_refs);
                                v2_ids.filter_refs(keep_refs);
                            }
                            reindex_vars();
                        }
                        auto y2 = static_pointer_cast<var<type>>(it2->second);
                        auto added2 = y2->add_bounds(lb2,ub2);
                        lt2._p = make_shared<gravity::pair< shared_ptr<param_>,shared_ptr<param_> >>(make_pair(make_shared<var<type>>(y2->in(ids)), make_shared<var<type>>(y2->in(ids))));
                        if(!added2.empty()){
                            assert(v1._indices->size()==v2._indices->size());
                            if(added2.size()!=v1._indices->size()){/* If some keys are repeated, remove them from the refs of o1 and o2 */
                                auto keep_refs = ids.diff_refs(added2);
                                v1_ids.filter_refs(keep_refs);
                                v2_ids.filter_refs(keep_refs);
                            }
                            reindex_vars();
                        }
                        
                        //create constraints for linking the auxiliary variables to the bilinear term
                        Constraint<type> link1(pair.first+"_link1");
                        link1 = y1->in(added1) - (v1.in(v1_ids) - v2.in(v2_ids))/2;
                        add(link1.in(unique_ids) == 0);
                        
                        Constraint<type> link2(pair.first+"_link2");
                        link2 = y2->in(added2) - (v1.in(v1_ids) + v2.in(v2_ids))/2;
                        add(link2.in(unique_ids) == 0);
                        
                    }
                    //insert the standardized terms into the new constraint
                    newc_standard.insert(lt1);
                    newc_standard.insert(lt2);
                }
                else { /* simply insert the squared term (not bilinear) */
                    newc_standard.insert(pair.second);
                }
            }
            //set the other properties of the standardized constraint properly
            newc_standard._range = c._range;
            newc_standard._all_convexity = c._all_convexity;
            newc_standard._all_sign = c._all_sign;
            newc_standard._ftype = c._ftype;
            newc_standard._ctype = c._ctype;
            newc_standard._indices = c._indices;
            newc_standard._dim[0] = c._dim[0];
            newc_standard._dim[1] = c._dim[1];
            return newc_standard;
        }
        
        
        /** Lift and linearize the nonlinear constraint c, return the linearized form and add linking constraints to the model.
         @param[in] c: constraint to linearize
         @param[in] partition_model: formulation used for partitionning the nonconvex parts of the constraint
         @return the linearized constraint
         @note This function will add constraints linking the lifted variables to the original ones, if a variable's partition is greater than 1, it will also add the disjunctive constraints corresponding to the partitionning of the variables.
         **/
        template<typename T=type,typename std::enable_if<is_arithmetic<type>::value>::type* = nullptr>
        Constraint<type> lift(Constraint<type>& c, string model_type){
            if(c.is_constant() || c.is_linear()){
                return c;
            }
            if(c.is_nonlinear() || c.is_polynomial()){
                throw invalid_argument("lift can only be called on quadratic constraints for now");
            }
            /* Lambda models are taken from Padberg's paper as they are described in type II and type III */
            if((model_type != "on/off") && (model_type != "lambda_II") && (model_type != "lambda_III")){
                throw invalid_argument("model_type can only be one of the following: 'on/off', 'lambda_II', 'lambda_III' ");
            }
            Constraint<type> lifted(c._name+"_lifted");
            if (!c.get_cst()->is_zero()) {
                if (c.get_cst()->is_number()) {
                    auto f_cst = static_pointer_cast<constant<type>>(c.get_cst());
                    lifted.add_cst(*f_cst);
                }
                else if (c.get_cst()->is_param()) {
                    auto f_cst = static_pointer_cast<param<type>>(c.get_cst());
                    lifted.add_cst(*f_cst);
                }
                else {
                    auto f_cst = static_pointer_cast<func<type>>(c.get_cst());
                    lifted.add_cst(*f_cst);
                }
                if (lifted._cst->is_function()) {
                    lifted.embed(*static_pointer_cast<func<type>>(lifted._cst));
                }
            }
            for (auto &pair:*c._lterms) {
                auto term = pair.second;
                if (term._coef->is_function()) {
                    auto coef = *static_pointer_cast<func<type>>(term._coef);
                    term._coef = func<type>(coef).copy();
                }
                else if(term._coef->is_param()) {
                    auto coef = *static_pointer_cast<param<type>>(term._coef);
                    term._coef = param<type>(coef).copy();
                }
                else if(term._coef->is_number()) {
                    auto coef = *static_pointer_cast<constant<type>>(term._coef);
                    term._coef = constant<type>(coef).copy();//TODO if T2==type no need to cast
                }
                lifted.insert(term);
            }
            bool lift_sign; /* create lift_sign for correct lower/upper bounding of the variables */
            for (auto &pair:*c._qterms) {
                auto term = pair.second;
                lterm lt;
                lt._sign = term._sign;
                if (term._coef->is_function()) {
                    auto coef = *static_pointer_cast<func<type>>(term._coef);
                    lt._coef = func<type>(coef).copy();
                }
                else if(term._coef->is_param()) {
                    auto coef = *static_pointer_cast<param<type>>(term._coef);
                    lt._coef = param<type>(coef).copy();
                }
                else if(term._coef->is_number()) {
                    auto coef = *static_pointer_cast<constant<type>>(term._coef);
                    lt._coef = constant<type>(coef).copy();
                    lift_sign = (term._sign ^ coef.is_negative()); //TODO: update prod_sign in other cases of coef type. Don't know how to do!
                }
                
                if (c.func<type>::is_concave()) //reverse the sign if the constraint is concave
                {
                    DebugOn("Changing the sign of the lifted variable." << endl);
                     lift_sign = !lift_sign;
                }
                else{
                    DebugOn("Keeping the sign of the lifted variable same." << endl);
                }
                
                //arrange the variables so that if they have the same base name, use them ordered in name
                auto o1 = *static_pointer_cast<var<type>>(term._p->first);
                auto o2 = *static_pointer_cast<var<type>>(term._p->second);
                if((o1 != o2) && (o1.get_name(true,true) == o2.get_name(true,true)) && (o1._name > o2._name) ){
                    o2 = *static_pointer_cast<var<type>>(term._p->first);
                    o1 = *static_pointer_cast<var<type>>(term._p->second);
                    DebugOff("O1 name "<< o1._name << endl);
                    DebugOff("O2 name "<< o2._name << endl);
                }
                
                string name;
                indices ids;
                if(o1==o2){
                    name = "Lift("+o1.get_name(true,true)+"^2)";
                    ids = *o1._indices;
                }
                else {
                    name = "Lift("+o1.get_name(true,true)+o2.get_name(true,true)+")";
                    ids = combine(*o1._indices,*o2._indices);
                }
                auto unique_ids = ids.get_unique_keys(); /* In case of an indexed variable, keep the unique keys only */
                auto o1_ids = *o1._indices;
                auto o2_ids = *o2._indices;
                if(unique_ids.size()!=ids.size()){/* If some keys are repeated, remove them from the refs of o1 and o2 */
                    auto keep_refs = ids.get_unique_refs();
                    o1_ids.filter_refs(keep_refs);
                    o2_ids.filter_refs(keep_refs);
                }
                
                // collect the number of partitions of each variable
                int num_partns1 = o1._num_partns;
                int num_partns2 = o2._num_partns;
                
                param<type> lb("lb"), ub("ub");
                lb.in(unique_ids);ub.in(unique_ids);
                
                //calculate the tightest valid bounds
                if(o1==o2) //if variables are same, calculate the bounds more efficiently
                {
                    for (int i=0; i<unique_ids.size(); i++) {
                        //calculate all the possibilities and assign the worst case
                        size_t id1;
                        if(o1_ids._ids == nullptr){
                            id1 = i;
                        }
                        else id1 = o1_ids._ids->at(0).at(i);
                        auto key1 = o1_ids._keys->at(id1);
                        
                        auto prod_b1 = o1.get_lb(key1)*o1.get_lb(key1);
                        auto prod_b2 = o1.get_lb(key1)*o1.get_ub(key1);
                        auto prod_b3 = o1.get_ub(key1)*o1.get_ub(key1);
                        
                        lb.set_val(key1, std::max(std::min(std::min(prod_b1,prod_b2), prod_b3), (type)0 ));
                        ub.set_val(key1, std::max(std::max(prod_b1,prod_b2),prod_b3));
                    }
                }
                else //if variables are different, need to check all four combinations
                {
                    for (int i=0; i<unique_ids.size(); i++) {
                        //calculate all the possibilities and assign the worst case
                        size_t id1;
                        size_t id2;
                        if(o1_ids._ids == nullptr){
                            id1 = i;
                        }
                        else id1 = o1_ids._ids->at(0).at(i);
                        if(o2_ids._ids == nullptr){
                            id2 = i;
                        }
                        else id2 = o2_ids._ids->at(0).at(i);
                        auto key1 = o1_ids._keys->at(id1);
                        auto key2 = o2_ids._keys->at(id2);
                        
                        auto prod_b1 = o1.get_lb(key1)*o2.get_lb(key2);
                        auto prod_b2 = o1.get_lb(key1)*o2.get_ub(key2);
                        auto prod_b3 = o1.get_ub(key1)*o2.get_lb(key2);
                        auto prod_b4 = o1.get_ub(key1)*o2.get_ub(key2);
                        
                        lb.set_val(key1+","+key2, std::min(std::min(prod_b1,prod_b2),std::min(prod_b3,prod_b4)));
                        ub.set_val(key1+","+key2, std::max(std::max(prod_b1,prod_b2),std::max(prod_b3,prod_b4)));
                    }
                }
                
                auto it = _vars_name.find(name);
                
                auto name1 = o1.get_name(true,true);
                auto name2 = o2.get_name(true,true);
                
                if(it==_vars_name.end()){
                    
                    // If some keys are repeated in individual indices, remove them from the refs of o1 and o2
                    auto o1_ids_uq = o1_ids;
                    auto o2_ids_uq = o2_ids;
                    auto keep_refs1 = o1_ids_uq.get_unique_refs();
                    auto keep_refs2 = o2_ids_uq.get_unique_refs();
                    o1_ids_uq.filter_refs(keep_refs1);
                    o2_ids_uq.filter_refs(keep_refs2);
                    reindex_vars();
                    
                    //create the lifted variable with proper lower and upper bounds
                    var<type> vlift(name, lb, ub);
                    vlift._lift = true;
                    add(vlift.in(unique_ids));
                    lt._p = make_shared<var<type>>(vlift.in(ids));
                    
                    //check the sign of the lift and the correspoinding bounding functions
                    if(c.check_soc() && c.is_eq()){
                    if(lift_sign){
                        vlift._lift_ub = true;
                        vlift._lift_lb = false;
                    }
                    else{
                        vlift._lift_ub = false;
                        vlift._lift_lb = true;
                    }
                    }
                    else{
                        vlift._lift_ub = true;
                        vlift._lift_lb = true;
                    }
                        
                        
                    if((num_partns1 > 1) || (num_partns2 > 1)) {
                        if (o1 == o2) //if the variables are same add 1d partition
                        {

                            DebugOn("<<<<<<<<<< THIS IS NOT SEEN BOTH -> SINGLE <<<<<<<<<<<" << endl);
                            //create the binary variables for the partitions
                            var<int> on(name1+"_binary",0,1);
                            
                            //create the proper indices and add the binary variables to the model
                            indices partns("partns");
                            for (int i = 0; i < num_partns1 ; ++i)
                            {
                                partns.add(name1+ "{" +to_string(i+1) + "}");
                            }
                            auto inst_partition = indices(unique_ids,partns);
                            add(on.in(inst_partition));
                            
                            //collect the number of entries in each of the index set
                            auto nb_entries_v1 = o1_ids.get_nb_entries();
                            auto nb_entries = unique_ids.get_nb_entries();
                            auto total_entries = inst_partition.get_nb_entries();
                            
                            //add the binary assignment constraint
                            Constraint<> onSum(pair.first + "_binarySum");
                            onSum += sum(on.in_matrix(nb_entries,total_entries-nb_entries));
                            add(onSum.in(unique_ids) == 1);
                            
                            //if the model type is selected as on/off, call on_off formulation for activation of individual constraints
                            if(model_type == "on/off"){
                                add_on_off_McCormick_refined(pair.first, vlift.in(unique_ids), o1.in(o1_ids), o2.in(o2_ids), on);
                            }
                            
                            else{ //means it is one of the lambda formulations
                                
                                //difference is this has one more partition index
                                indices partns_lambda("partns_lambda");
                                for (int i = 0; i < num_partns1+1 ; ++i)
                                {
                                    partns_lambda.add(name1+ "{" +to_string(i+1) + "}");
                                }
                                auto inst_partition_lambda = indices(unique_ids,partns_lambda);
                                
                                // Convex combination variables
                                var<> lambda(name1+"_lambda",pos_);
                                add(lambda.in(inst_partition_lambda));
                                
                                /** Parameters */
                                // Bounds on variable v1 & v2
                                param<> bounds(name1+"_bounds");
                                bounds.in(inst_partition_lambda);
                                
                                // Function values on the extreme points
                                param<> EP(name1+name2+"_grid_values");
                                EP.in(inst_partition_lambda);
                                
                                size_t nb_ins = vlift.in(unique_ids).get_nb_inst();
                                auto o1_global_lb = o1.get_lb();
                                auto increment = (o1.get_ub() - o1_global_lb)/num_partns1;
                                
                                // fill bounds and function values
                                for (int i=0 ; i<num_partns1+1; ++i) {
                                    auto bound_partn = o1_global_lb + increment*i;
                                    bound_partn.eval_all();
                                    for (size_t inst = 0; inst< nb_ins; inst++){
                                        auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                        auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                        string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"}";
                                        bounds.set_val(cur_idx,bound_partn.eval(inst));
                                        EP.set_val(cur_idx,(bound_partn.eval(inst)*bound_partn.eval(inst)));
                                    }
                                }
                                
                                // Lambda coefficient matrix when linking with partition variables
                                param<> lambda_coef(name1+"_lambda_linking_coefficients");
                                
                                // Partition coefficient matrix when linking with lambda variables
                                param<> on_coef(name1+"_partition_linking_coefficients");
                                
                                // create constraint indices
                                indices const_idx("const_idx");
                                
                                if(model_type == "lambda_II"){
                                    
                                    //fill constraint indices
                                    for (int i=0; i<num_partns1+1; ++i){
                                        const_idx.add(to_string(i+1));
                                    }
                                    
                                    // Lambda coefficient matrix when linking with partition variables
                                    lambda_coef.in(indices(inst_partition_lambda, const_idx));
                                    
                                    // Partition coefficient matrix when linking with lambda variables
                                    on_coef.in(indices(inst_partition, const_idx));
                                    
                                    // fill lambda_coef
                                    for (size_t inst = 0; inst< nb_ins; inst++){
                                        auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                        auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                        for (int i=0 ; i<num_partns1+1; ++i) {
                                            string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"}," + to_string(i+1);
                                            lambda_coef.set_val(cur_idx,1);
                                        }
                                    }
                                    
                                    // fill on_coef
                                    for (size_t inst = 0; inst< nb_ins; inst++){
                                        auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                        auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                        string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"}," + to_string(1);
                                        on_coef.set_val(cur_idx,1);
                                        for (int i=1 ; i<num_partns1; ++i) {
                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(i)+"}," + to_string(i+1);
                                            on_coef.set_val(cur_idx,1);
                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"}," + to_string(i+1);
                                            on_coef.set_val(cur_idx,1);
                                        }
                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(num_partns1)+"}," + to_string(num_partns1+1);
                                        on_coef.set_val(cur_idx,1);
                                    }
                                }
                                
                                else /*means model_type == "lambda_III" */{
                                    
                                    //fill constraint indices
                                    for (int i=0; i<(num_partns1-2)*2+2; ++i){
                                        const_idx.add(to_string(i+1));
                                    }
                                    
                                    // Lambda coefficient matrix when linking with partition variables
                                    lambda_coef.in(indices(inst_partition_lambda, const_idx));
                                    
                                    // Partition coefficient matrix when linking with lambda variables
                                    on_coef.in(indices(inst_partition, const_idx));
                                    
                                    // fill lambda_coef
                                    for (size_t inst = 0; inst< nb_ins; inst++){
                                        auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                        auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                        string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"}," + to_string(1);
                                        lambda_coef.set_val(cur_idx,1);
                                        for (int i=1 ; i<(num_partns1-2)*2+1; i=i+2) {
                                            for (int j=(i-1)/2 + 2; j<num_partns1+1; ++j) {
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"}," + to_string(i+1) ;
                                                lambda_coef.set_val(cur_idx,1);
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"}," + to_string(i+2);
                                                lambda_coef.set_val(cur_idx,-1);
                                            }
                                        }
                                        cur_idx =  cur_var_idx+","+name1+"{"+to_string(num_partns1+1)+"}," + to_string((num_partns1-2)*2+2);
                                        lambda_coef.set_val(cur_idx,1);
                                    }
                                    
                                    // fill on_coef
                                    for (size_t inst = 0; inst< nb_ins; inst++){
                                        auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                        auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                        string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"}," + to_string(1);
                                        on_coef.set_val(cur_idx,1);
                                        
                                        for (int i=1; i<num_partns1; ++i) {
                                            cur_idx =  cur_var_idx+","+name1+"{"+to_string(i+1)+"}," + to_string(2);
                                            on_coef.set_val(cur_idx, 1);
                                        }
                                        
                                        for (int i=2 ; i<(num_partns1-2)*2+2; i=i+2) {
                                            for (int j=i/2+1; j<num_partns1; ++j) {
                                                cur_idx =  cur_var_idx +","+name1+"{"+to_string(j+1)+"}," +  to_string(i+1);
                                                on_coef.set_val(cur_idx,-1);
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"}," + to_string(i+2) ;
                                                on_coef.set_val(cur_idx,1);
                                            }
                                        }
                                    }
                                    
                                }
                                
                                
                                /** Constraints */
                                if (vlift._lift_ub){
                                // Representation of the quadratic term with secant
                                Constraint<> quad_ub(pair.first+"_quad_ub");
                                quad_ub = EP.in_matrix(nb_entries,total_entries-nb_entries)*lambda.in_matrix(nb_entries,total_entries-nb_entries) - vlift.in(unique_ids);
                                add(quad_ub.in(unique_ids) >= 0); /*using it as the upper bound to be valid*/
                                }
                                if (vlift._lift_lb){
                                Constraint<> quad_lb(pair.first+"_quad_lb");
                                quad_lb = pow(o1.from_ith(0,unique_ids),2) - vlift.in(unique_ids);
                                quad_lb._relaxed = true;
                                add(quad_lb.in(unique_ids) <= 0); /*using it as the lower bound to be valid*/
                                }
                                
                                // Representation of o1 with convex combination
                                Constraint<> o1_rep(pair.first+"_o1_rep");
                                o1_rep = bounds.from_ith(0,inst_partition_lambda).in_matrix(nb_entries, 1) * lambda.in_matrix(nb_entries,total_entries-nb_entries) - o1.in(o1_ids);
                                add(o1_rep.in(unique_ids) == 0);
                                
                                // Linking partition variables with lambda for both lambda formulations
                                if(model_type == "lambda_II"){
                                    Constraint<> on_link_lambda(pair.first+"_on_link_lambda_II");
                                    on_link_lambda = lambda.in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef.in_matrix(nb_entries,total_entries-nb_entries) - on.in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,on_coef.get_matrix_ids(nb_entries,total_entries-nb_entries)) * on_coef.in_matrix(nb_entries,total_entries-nb_entries);
                                    add(on_link_lambda.in(indices(unique_ids,const_idx)) <= 0);

                                }
                                else{ //means model_type == "lambda_III"
                                    Constraint<> on_link_lambda(pair.first+"_on_link_lambda_III");
                                    on_link_lambda = lambda.in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef.in_matrix(nb_entries,total_entries-nb_entries) -
                                    on.in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,on_coef.get_matrix_ids(nb_entries,total_entries-nb_entries)) * on_coef.in_matrix(nb_entries,total_entries-nb_entries);
                                    add(on_link_lambda.in(indices(unique_ids,const_idx)) <= 0);

                                }
                                
                                
                                // sum over lambda
                                Constraint<> lambdaSum(pair.first+"_lambdaSum");
                                lambdaSum = sum(lambda.in_matrix(nb_entries,total_entries-nb_entries));
                                add(lambdaSum.in(unique_ids) == 1);
                            }
                            
                        }
                        else{ //else add 2d partition
                            
                            auto binvar_ptr1 = _vars_name.find(name1+"_binary");
                            auto binvar_ptr2 = _vars_name.find(name2+"_binary");
                            
                            if(binvar_ptr1 !=_vars_name.end()){ //means v1 has been partitioned before
                                
                                //if the variables are same core name (means they are same symbolic variable with different indexing)
                                if(name1 == name2){
                                    DebugOn("<<<<<<<<<< THIS IS NOT SEEN BOTH -> SEEN FIRST -> SAME VARS <<<<<<<<<<<" << endl);
                                    
                                    //create the proper index set for partitions
                                    indices partns1("partns1");
                                    for (int i = 0; i < num_partns1 ; ++i)
                                    {
                                        partns1.add(name1+ "{" +to_string(i+1) + "}");
                                    }
                                    
                                    indices partns2("partns2");
                                    for (int i = 0; i < num_partns2 ; ++i)
                                    {
                                        partns2.add(name2+ "{" +to_string(i+1) + "}");
                                    }
                                    
                                    //cast the variable pointer for further use
                                    auto binvar1 = static_pointer_cast<var<int>>(binvar_ptr1->second);
                                    
                                    //define the lower and upper bounds
                                    param<int> lb1("lb1"), ub1("ub1");
                                    lb1.in(union_ids(o1_ids_uq, o2_ids_uq),partns1);
                                    ub1.in(union_ids(o1_ids_uq, o2_ids_uq),partns1);
                                    lb1.set_val(0), ub1.set_val(1);
                                    auto added1 = binvar1->add_bounds(lb1,ub1);
                                    reindex_vars();
                                    
                                    //collect the number of entries in each index set
                                    auto nb_entries_v1 = o1_ids.get_nb_entries();
                                    auto nb_entries_v2 = o2_ids.get_nb_entries();
                                    auto nb_entries = unique_ids.get_nb_entries();
                                    
                                    //if there are new indices for the previously defined variable add the corresponding constraint for partitions
                                    if(!added1.empty()){
                                        Constraint<> onSum1(o1._name+"_binarySum");
                                        onSum1 = sum(binvar1->in(added1).in_matrix(nb_entries_v1,1));
                                        auto vset1 = added1.from_ith(0,nb_entries_v1);
                                        vset1.filter_refs(vset1.get_unique_refs());
                                        add(onSum1.in(vset1) == 1);
                                    }
                                    
                                    //if the on/off formulation is chosen for activating constraint
                                    if(model_type == "on/off"){
                                        var<int> on(name1+name2+"_binary",0,1);
                                        
                                        indices partns("partns");
                                        partns = indices(partns1,partns2);
                                        auto inst_partition = indices(unique_ids,partns);
                                        add(on.in(inst_partition));
                                        auto total_entries = inst_partition.get_nb_entries();
                                        
                                        Constraint<> onLink1(pair.first+"_binaryLink1");
                                        onLink1 = binvar1->from_ith(0,inst_partition.ignore_ith(nb_entries_v1, nb_entries_v2)) - on;
                                        add(onLink1.in(inst_partition) >= 0);
                                        
                                        Constraint<> onLink2(pair.first+"_binaryLink2");
                                        onLink2 = binvar1->in_ignore_ith(nb_entries_v2,1,inst_partition.ignore_ith(0,nb_entries_v1)) - on;
                                        add(onLink2.in(inst_partition) >= 0);
                                        
                                        Constraint<> onLink3(pair.first+"_binaryLink3");
                                        onLink3 = binvar1->from_ith(0,inst_partition.ignore_ith(nb_entries_v1, nb_entries_v2)) + binvar1->in_ignore_ith(nb_entries_v2,1,inst_partition.ignore_ith(0,nb_entries_v1)) - 1 - on;
                                        add(onLink3.in(inst_partition) <= 0);
                                        
                                        Constraint<> onSumComb(pair.first+"_binarySum");
                                        onSumComb = sum(on.in_matrix(nb_entries,total_entries-nb_entries));
                                        add(onSumComb.in(unique_ids) == 1);
                                        
                                        add_on_off_McCormick_refined(pair.first, vlift.in(unique_ids), o1.in(o1_ids), o2.in(o2_ids), on);
                                    }
                                    
                                    else{ //means it is one of the lambda formulations
                                        
                                        //difference is this has one more partition index
                                        indices partns1_lambda("partns1_lambda");
                                        for (int i = 0; i < num_partns1+1; ++i)
                                        {
                                            partns1_lambda.add(name1+ "{" +to_string(i+1) + "}");
                                        }
                                        
                                        indices partns2_lambda("partns2_lambda");
                                        for (int i = 0; i < num_partns2+1; ++i)
                                        {
                                            partns2_lambda.add(name2+ "{" +to_string(i+1) + "}");
                                        }
                                        
                                        indices partns_lambda("partns_lambda");
                                        partns_lambda = indices(partns1_lambda,partns2_lambda);
                                        auto inst_partition_lambda = indices(unique_ids,partns_lambda);
                                        auto inst_partition_bounds1 = indices(unique_ids,partns1_lambda);
                                        auto inst_partition_bounds2 = indices(unique_ids,partns2_lambda);
                                        
                                        // Convex combination variables
                                        var<> lambda(name1+name2+"_lambda",pos_);
                                        add(lambda.in(inst_partition_lambda));
                                        
                                        /** Parameters */
                                        // Bounds on variable v1 & v2
                                        param<> bounds1(name1+"_bounds1");
                                        bounds1.in(inst_partition_bounds1);
                                        
                                        param<> bounds2(name2+"_bounds2");
                                        bounds2.in(inst_partition_bounds2);
                                        
                                        // Function values on the extreme points
                                        param<> EP(name1+name2+"_grid_values");
                                        EP.in(inst_partition_lambda);
                                        auto total_entries = inst_partition_lambda.get_nb_entries();
                                        
                                        size_t nb_ins = vlift.in(unique_ids).get_nb_inst();
                                        auto o1_global_lb = o1.get_lb();
                                        auto increment1 = (o1.get_ub() - o1_global_lb)/num_partns1;
                                        
                                        auto o2_global_lb = o2.get_lb();
                                        auto increment2 = (o2.get_ub() - o2_global_lb)/num_partns2;
                                        
                                        // fill bounds and function values
                                        for (int i=0 ; i<num_partns1+1; ++i) {
                                            auto bound_partn1 = o1_global_lb + increment1*i;
                                            bound_partn1.eval_all();
                                            auto bound_partn2 = o2_global_lb + increment2*i;
                                            bound_partn2.eval_all();
                                            for (size_t inst = 0; inst< nb_ins; inst++){
                                                auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                                auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                                string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"}";
                                                bounds1.set_val(cur_idx,bound_partn1.eval(inst));
                                                bounds2.set_val(cur_idx,bound_partn2.eval(inst));
                                                for(int j=0; j<num_partns2+1; ++j){
                                                    auto bound_partn2_temp = o2_global_lb + increment2*j;
                                                    bound_partn2_temp.eval_all();
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"}";
                                                    EP.set_val(cur_idx,(bound_partn1.eval(inst)*bound_partn2_temp.eval(inst)));
                                                }
                                            }
                                        }
                                        
                                        // Lambda coefficient matrix when linking with partition variables
                                        param<> lambda_coef1(name1+"_lambda_linking_coefficients1");
                                        param<> lambda_coef2(name2+"_lambda_linking_coefficients2");
                                        
                                        // Partition coefficient matrix when linking with lambda variables
                                        param<> on_coef1(name1+"_partition_linking_coefficients1");
                                        param<> on_coef2(name2+"_partition_linking_coefficients2");
                                        
                                        // create constraint indices
                                        indices const_idx1("const_idx1");
                                        indices const_idx2("const_idx2");
                                        
                                        if(model_type == "lambda_II"){
                                            
                                            //fill constraint indices
                                            for (int i=0; i<num_partns1+1; ++i){
                                                const_idx1.add(to_string(i+1));
                                            }
                                            
                                            //fill constraint indices
                                            for (int i=0; i<num_partns2+1; ++i){
                                                const_idx2.add(to_string(i+1));
                                            }
                                            
                                            // Lambda coefficient matrix when linking with partition variables
                                            lambda_coef1.in(indices(inst_partition_lambda, const_idx1));
                                            lambda_coef2.in(indices(inst_partition_lambda, const_idx2));
                                            
                                            // Lambda coefficient matrix when linking with partition variables
                                            on_coef1.in(indices(unique_ids, partns1, const_idx1));
                                            on_coef2.in(indices(unique_ids, partns2, const_idx2));
                                            
                                            // fill lambda_coef1 and lambda_coef2
                                            for (size_t inst = 0; inst< nb_ins; inst++){
                                                auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                                auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                                for (int i=0 ; i<num_partns1+1; ++i) {
                                                    for (int j=0 ; j<num_partns2+1; ++j) {
                                                        string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                        lambda_coef1.set_val(cur_idx,1);
                                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(j+1);
                                                        lambda_coef2.set_val(cur_idx,1);
                                                    }
                                                }
                                            }
                                            
                                            // fill on_coef1 and on_coef2
                                            for (size_t inst = 0; inst< nb_ins; inst++){
                                                auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                                auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                                string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+to_string(1);
                                                on_coef1.set_val(cur_idx,1);
                                                cur_idx = cur_var_idx+","+name2+"{"+to_string(1)+"},"+to_string(1);
                                                on_coef2.set_val(cur_idx,1);
                                                for (int i=1 ; i<num_partns1; ++i) {
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(i)+"},"+to_string(i+1);
                                                    on_coef1.set_val(cur_idx,1);
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+to_string(i+1);
                                                    on_coef1.set_val(cur_idx,1);
                                                }
                                                for (int i=1 ; i<num_partns2; ++i) {
                                                    cur_idx = cur_var_idx+","+name2+"{"+to_string(i)+"},"+to_string(i+1);
                                                    on_coef2.set_val(cur_idx,1);
                                                    cur_idx = cur_var_idx+","+name2+"{"+to_string(i+1)+"},"+to_string(i+1);
                                                    on_coef2.set_val(cur_idx,1);
                                                }
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(num_partns1)+"},"+to_string(num_partns1+1);
                                                on_coef1.set_val(cur_idx,1);
                                                cur_idx = cur_var_idx+","+name2+"{"+to_string(num_partns2)+"},"+to_string(num_partns2+1);
                                                on_coef2.set_val(cur_idx,1);
                                            }
                                        }
                                        
                                        
                                        else /*means model_type == "lambda_III" */{
                                            
                                            //fill constraint indices
                                            for (int i=0; i<(num_partns1-2)*2+2; ++i){
                                                const_idx1.add(to_string(i+1));
                                            }
                                            
                                            //fill constraint indices
                                            for (int i=0; i<(num_partns2-2)*2+2; ++i){
                                                const_idx2.add(to_string(i+1));
                                            }
                                            
                                            // Lambda coefficient matrix when linking with partition variables
                                            lambda_coef1.in(indices(inst_partition_lambda, const_idx1));
                                            lambda_coef2.in(indices(inst_partition_lambda, const_idx2));
                                            
                                            // Partition coefficient matrix when linking with lambda variables
                                            on_coef1.in(indices(unique_ids, partns1, const_idx1));
                                            on_coef2.in(indices(unique_ids, partns2, const_idx2));
                                            
                                            // fill lambda_coef1 and lambda_coef2
                                            for (size_t inst = 0; inst< nb_ins; inst++){
                                                auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                                auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                                for (int j=0; j<num_partns2+1; ++j) {
                                                    string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(1);
                                                    lambda_coef1.set_val(cur_idx,1);
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(num_partns1+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string((num_partns1-2)*2+2);
                                                    lambda_coef1.set_val(cur_idx,1);
                                                }
                                                
                                                for (int i=1 ; i<(num_partns1-2)*2+1; i=i+2) {
                                                    for (int j=(i-1)/2 + 2; j<num_partns1+1; ++j) {
                                                        for(int k=0; k<num_partns2+1; ++k){
                                                            string cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+name2+"{"+to_string(k+1)+"},"+to_string(i+1);
                                                            lambda_coef1.set_val(cur_idx,1);
                                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+name2+"{"+to_string(k+1)+"},"+to_string(i+2);
                                                            lambda_coef1.set_val(cur_idx,-1);
                                                        }
                                                    }
                                                }
                                                
                                                for (int i=0; i<num_partns1+1; ++i) {
                                                    string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(1)+"},"+to_string(1);
                                                    lambda_coef2.set_val(cur_idx,1);
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(num_partns2+1)+"},"+to_string((num_partns2-2)*2+2);
                                                    lambda_coef2.set_val(cur_idx,1);
                                                }
                                                
                                                for (int i=1 ; i<(num_partns2-2)*2+1; i=i+2) {
                                                    for (int j=(i-1)/2 + 2; j<num_partns2+1; ++j) {
                                                        for(int k=0; k<num_partns1+1; ++k){
                                                            string cur_idx = cur_var_idx+","+name1+"{"+to_string(k+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                            lambda_coef2.set_val(cur_idx,1);
                                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(k+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                            lambda_coef2.set_val(cur_idx,-1);
                                                        }
                                                    }
                                                }
                                            }
                                            
                                            
                                            
                                            // fill on_coef1 and on_coef2
                                            for (size_t inst = 0; inst< nb_ins; inst++){
                                                auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                                auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                                string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+to_string(1);
                                                on_coef1.set_val(cur_idx,1);
                                                
                                                for (int i=1; i<num_partns1; ++i) {
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+to_string(2);
                                                    on_coef1.set_val(cur_idx, 1);
                                                }
                                                
                                                for (int i=2 ; i<(num_partns1-2)*2+2; i=i+2) {
                                                    for (int j=i/2+1; j<num_partns1; ++j) {
                                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                        on_coef1.set_val(cur_idx,-1);
                                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                        on_coef1.set_val(cur_idx,1);
                                                    }
                                                }
                                                
                                                cur_idx = cur_var_idx+","+name2+"{"+to_string(1)+"},"+to_string(1);
                                                on_coef2.set_val(cur_idx,1);
                                                for (int i=1; i<num_partns2; ++i) {
                                                    cur_idx = cur_var_idx+","+name2+"{"+to_string(i+1)+"},"+to_string(2);
                                                    on_coef2.set_val(cur_idx, 1);
                                                }
                                                
                                                for (int i=2 ; i<(num_partns2-2)*2+2; i=i+2) {
                                                    for (int j=i/2+1; j<num_partns2; ++j) {
                                                        cur_idx = cur_var_idx+","+name2+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                        on_coef2.set_val(cur_idx,-1);
                                                        cur_idx = cur_var_idx+","+name2+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                        on_coef2.set_val(cur_idx,1);
                                                    }
                                                }
                                            }
                                        }
                                        
                                        
                                        /** Constraints */
                                        // Representation of the bilinear term with convex combination
                                        Constraint<> bln_rep(pair.first+"_bln_rep");
                                        bln_rep = EP.in_matrix(nb_entries,total_entries-nb_entries)*lambda.in_matrix(nb_entries,total_entries-nb_entries) - vlift.in(unique_ids);
                                        add(bln_rep.in(unique_ids) == 0);
                                        
                                        // Representation of o1 with convex combination
                                        Constraint<> o1_rep(pair.first+"_o1_rep");

                                        o1_rep = bounds1.from_ith(0,inst_partition_lambda).in_matrix(nb_entries, 1) * lambda.in_matrix(nb_entries,total_entries-nb_entries) - o1.in(o1_ids);

                                        add(o1_rep.in(unique_ids) == 0);
                                        
                                        // Representation of o2 with convex combination
                                        Constraint<> o2_rep(pair.first+"_o2_rep");
              
                                        o2_rep = bounds2.in_ignore_ith(nb_entries, 1, inst_partition_lambda).in_matrix(nb_entries,1) * lambda.in_matrix(nb_entries,total_entries-nb_entries) - o2.in(o2_ids);
                                        add(o2_rep.in(unique_ids) == 0);
                                        
                                        // Linking partition variables1 with lambda
                                        if(model_type == "lambda_II"){
                                            Constraint<> on_link_lambda1(pair.first+"_on_link_lambda1_II");
                   
                                            on_link_lambda1 = lambda.in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef1.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef1.in_matrix(nb_entries,total_entries-nb_entries) - binvar1->in_ignore_ith(nb_entries_v1,nb_entries_v2,on_coef1.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef1.in_matrix(nb_entries,1);
                                            add(on_link_lambda1.in(indices(unique_ids,const_idx1)) <= 0);
                                            
                                            Constraint<> on_link_lambda2(pair.first+"_on_link_lambda2_II");
                                            on_link_lambda2 = lambda.in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef2.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef2.in_matrix(nb_entries,total_entries-nb_entries) - binvar1->in_ignore_ith(0,nb_entries_v1,on_coef2.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef2.in_matrix(nb_entries,1);
                                            add(on_link_lambda2.in(indices(unique_ids,const_idx2)) <= 0);

                                        }
                                        else{
                                            Constraint<> on_link_lambda1(pair.first+"_on_link_lambda1_III");

                                            on_link_lambda1 = lambda.in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef1.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef1.in_matrix(nb_entries,total_entries-nb_entries) - binvar1->in_ignore_ith(nb_entries_v1,nb_entries_v2,on_coef1.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef1.in_matrix(nb_entries,1);
                                            add(on_link_lambda1.in(indices(unique_ids,const_idx1)) <= 0);
                                            
                                            Constraint<> on_link_lambda2(pair.first+"_on_link_lambda2_III");
                                            on_link_lambda2 = lambda.in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef2.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef2.in_matrix(nb_entries,total_entries-nb_entries) - binvar1->in_ignore_ith(0,nb_entries_v1,on_coef2.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef2.in_matrix(nb_entries,1);
                                            add(on_link_lambda2.in(indices(unique_ids,const_idx2)) <= 0);
                                        }
                                        // sum over lambda
                                        Constraint<> lambdaSum(pair.first+"_lambdaSum");
                                        lambdaSum = sum(lambda.in_matrix(nb_entries,total_entries-nb_entries));
                                        add(lambdaSum.in(unique_ids) == 1);
                                    }
                                }
                                
                                else{
                                    DebugOn("<<<<<<<<<< THIS IS NOT SEEN BOTH -> SEEN FIRST -> DIFF VARS <<<<<<<<<<<" << endl);
                                    
                                    auto binvar1 = static_pointer_cast<var<int>>(binvar_ptr1->second);
                                    indices partns1("partns1");
                                    for (int i = 0; i < num_partns1 ; ++i)
                                    {
                                        partns1.add(name1+ "{" +to_string(i+1) + "}");
                                    }
                                    
                                    param<int> lb1("lb1"), ub1("ub1");
                                    lb1.in(o1_ids_uq,partns1);
                                    ub1.in(o1_ids_uq,partns1);
                                    lb1.set_val(0), ub1.set_val(1);
                                    auto added1 = binvar1->add_bounds(lb1,ub1);
                                    reindex_vars();
                                    
                                    var<int> on2(name2+"_binary",0,1);
                                    indices partns2("partns2");
                                    for (int i = 0; i < num_partns2 ; ++i)
                                    {
                                        partns2.add(name2+ "{" + to_string(i+1) + "}");
                                    }
                                    
                                    
                                    add(on2.in(o2_ids_uq,partns2));
                                    
                                    auto nb_entries_v1 = o1_ids.get_nb_entries();
                                    auto nb_entries_v2 = o2_ids.get_nb_entries();
                                    auto nb_entries = unique_ids.get_nb_entries();
                                    
                                    if(!added1.empty()){
                                        Constraint<> onSum1(o1._name+"_binarySum");
                                        onSum1 = sum(binvar1->in(added1).in_matrix(nb_entries_v1,1));
                                        auto vset1 = added1.from_ith(0,nb_entries_v1);
                                        vset1.filter_refs(vset1.get_unique_refs());
                                        add(onSum1.in(vset1) == 1);
                                    }
                                    
                                    Constraint<> onSum2(o2._name+"_binarySum");
                                    onSum2 = sum(on2.in_matrix(nb_entries_v2,1));
                                    add(onSum2.in(o2_ids_uq) == 1);
                                    
                                    if(model_type == "on/off"){ //if on/off is chosen
                                        
                                        var<int> on(name1+name2+"_binary",0,1);
                                        
                                        indices partns("partns");
                                        partns = indices(partns1,partns2);
                                        auto inst_partition = indices(unique_ids,partns);
                                        add(on.in(inst_partition));
                                        auto total_entries = inst_partition.get_nb_entries();
                                        
                                        Constraint<> onLink1(pair.first+"_binaryLink1");
                                        onLink1 = binvar1->from_ith(0,inst_partition.ignore_ith(nb_entries_v1, nb_entries_v2)) - on;
                                        add(onLink1.in(inst_partition) >= 0);
                                        
                                        Constraint<> onLink2(pair.first+"_binaryLink2");
                                        onLink2 = on2.in_ignore_ith(nb_entries_v2,1,inst_partition.ignore_ith(0,nb_entries_v1)) - on;
                                        add(onLink2.in(inst_partition) >= 0);
                                        
                                        Constraint<> onLink3(pair.first+"_binaryLink3");
                                        onLink3 = binvar1->from_ith(0,inst_partition.ignore_ith(nb_entries_v1, nb_entries_v2)) + on2.in_ignore_ith(nb_entries_v2,1,inst_partition.ignore_ith(0,nb_entries_v1)) - 1 - on;
                                        add(onLink3.in(inst_partition) <= 0);
                                        
                                        Constraint<> onSumComb(pair.first+"_binarySum");
                                        onSumComb = sum(on.in_matrix(nb_entries,total_entries-nb_entries));
                                        add(onSumComb.in(unique_ids) == 1);
                                        
                                        add_on_off_McCormick_refined(pair.first, vlift.in(unique_ids), o1.in(o1_ids), o2.in(o2_ids), on);
                                    }
                                    
                                    else{ //means it is one of the lambda formulations
                                        
                                        //difference is this has one more partition index
                                        indices partns1_lambda("partns1_lambda");
                                        for (int i = 0; i < num_partns1+1; ++i)
                                        {
                                            partns1_lambda.add(name1+ "{" +to_string(i+1) + "}");
                                        }
                                        
                                        indices partns2_lambda("partns2_lambda");
                                        for (int i = 0; i < num_partns2+1; ++i)
                                        {
                                            partns2_lambda.add(name2+ "{" +to_string(i+1) + "}");
                                        }
                                        
                                        indices partns_lambda("partns_lambda");
                                        partns_lambda = indices(partns1_lambda,partns2_lambda);
                                        auto inst_partition_lambda = indices(unique_ids,partns_lambda);
                                        auto inst_partition_bounds1 = indices(unique_ids,partns1_lambda);
                                        auto inst_partition_bounds2 = indices(unique_ids,partns2_lambda);
                                        
                                        // Convex combination variables
                                        var<> lambda(name1+name2+"_lambda",pos_);
                                        add(lambda.in(inst_partition_lambda));
                                        
                                        /** Parameters */
                                        // Bounds on variable v1 & v2
                                        param<> bounds1(name1+"_bounds1");
                                        bounds1.in(inst_partition_bounds1);
                                        
                                        param<> bounds2(name2+"_bounds2");
                                        bounds2.in(inst_partition_bounds2);
                                        
                                        // Function values on the extreme points
                                        param<> EP(name1+name2+"_grid_values");
                                        EP.in(inst_partition_lambda);
                                        auto total_entries = inst_partition_lambda.get_nb_entries();
                                        
                                        size_t nb_ins = vlift.in(unique_ids).get_nb_inst();
                                        auto o1_global_lb = o1.get_lb();
                                        auto increment1 = (o1.get_ub() - o1_global_lb)/num_partns1;
                                        
                                        auto o2_global_lb = o2.get_lb();
                                        auto increment2 = (o2.get_ub() - o2_global_lb)/num_partns2;
                                        
                                        // fill bounds1 and function values
                                        for (int i=0 ; i<num_partns1+1; ++i) {
                                            auto bound_partn1 = o1_global_lb + increment1*i;
                                            bound_partn1.eval_all();
                                            for (size_t inst = 0; inst< nb_ins; inst++){
                                                auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                                auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                                string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"}";
                                                bounds1.set_val(cur_idx,bound_partn1.eval(inst));
                                                for(int j=0; j<num_partns2+1; ++j){
                                                    auto bound_partn2 = o2_global_lb + increment2*j;
                                                    bound_partn2.eval_all();
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"}";
                                                    EP.set_val(cur_idx,(bound_partn1.eval(inst)*bound_partn2.eval(inst)));
                                                }
                                            }
                                        }
                                        // fill bounds2
                                        for (int i=0 ; i<num_partns2+1; ++i) {
                                            auto bound_partn2 = o2_global_lb + increment2*i;
                                            bound_partn2.eval_all();
                                            for (size_t inst = 0; inst< nb_ins; inst++){
                                                auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                                auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                                string cur_idx = cur_var_idx+","+name2+"{"+to_string(i+1)+"}";
                                                bounds2.set_val(cur_idx,bound_partn2.eval(inst));
                                            }
                                        }
                                        
                                        // Lambda coefficient matrix when linking with partition variables
                                        param<> lambda_coef1(name1+"_lambda_linking_coefficients1");
                                        param<> lambda_coef2(name2+"_lambda_linking_coefficients2");
                                        
                                        // Partition coefficient matrix when linking with lambda variables
                                        param<> on_coef1(name1+"_partition_linking_coefficients1");
                                        param<> on_coef2(name2+"_partition_linking_coefficients2");
                                        
                                        // create constraint indices
                                        indices const_idx1("const_idx1");
                                        indices const_idx2("const_idx2");
                                        
                                        if(model_type == "lambda_II"){
                                            
                                            //fill constraint indices
                                            for (int i=0; i<num_partns1+1; ++i){
                                                const_idx1.add(to_string(i+1));
                                            }
                                            
                                            //fill constraint indices
                                            for (int i=0; i<num_partns2+1; ++i){
                                                const_idx2.add(to_string(i+1));
                                            }
                                            
                                            // Lambda coefficient matrix when linking with partition variables
                                            if(num_partns1 > 1) lambda_coef1.in(indices(inst_partition_lambda, const_idx1));
                                            if(num_partns2 > 1) lambda_coef2.in(indices(inst_partition_lambda, const_idx2));
                                            
                                            // Lambda coefficient matrix when linking with partition variables
                                            if(num_partns1 > 1) on_coef1.in(indices(unique_ids, partns1, const_idx1));
                                            if(num_partns2 > 1) on_coef2.in(indices(unique_ids, partns2, const_idx2));
                                            
                                            // fill lambda_coef1 and lambda_coef2
                                            for (size_t inst = 0; inst< nb_ins; inst++){
                                                auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                                auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                                for (int i=0 ; i<num_partns1+1; ++i) {
                                                    for (int j=0 ; j<num_partns2+1; ++j) {
                                                        string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                        if(num_partns1 > 1) lambda_coef1.set_val(cur_idx,1);
                                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(j+1);
                                                        if(num_partns2 > 1) lambda_coef2.set_val(cur_idx,1);
                                                    }
                                                }
                                            }
                                            
                                            // fill on_coef1 and on_coef2
                                            for (size_t inst = 0; inst< nb_ins; inst++){
                                                auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                                auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                                string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+to_string(1);
                                                if(num_partns1 > 1) on_coef1.set_val(cur_idx,1);
                                                cur_idx = cur_var_idx+","+name2+"{"+to_string(1)+"},"+to_string(1);
                                                if(num_partns2 > 1) on_coef2.set_val(cur_idx,1);
                                                if(num_partns1 > 1){
                                                    for (int i=1 ; i<num_partns1; ++i) {
                                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(i)+"},"+to_string(i+1);
                                                        on_coef1.set_val(cur_idx,1);
                                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+to_string(i+1);
                                                        on_coef1.set_val(cur_idx,1);
                                                    }
                                                }
                                                if(num_partns2 > 1){
                                                    for (int i=1 ; i<num_partns2; ++i) {
                                                        cur_idx = cur_var_idx+","+name2+"{"+to_string(i)+"},"+to_string(i+1);
                                                        on_coef2.set_val(cur_idx,1);
                                                        cur_idx = cur_var_idx+","+name2+"{"+to_string(i+1)+"},"+to_string(i+1);
                                                        on_coef2.set_val(cur_idx,1);
                                                    }
                                                }
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(num_partns1)+"},"+to_string(num_partns1+1);
                                                if(num_partns1 > 1) on_coef1.set_val(cur_idx,1);
                                                cur_idx = cur_var_idx+","+name2+"{"+to_string(num_partns2)+"},"+to_string(num_partns2+1);
                                                if(num_partns2 > 1) on_coef2.set_val(cur_idx,1);
                                            }
                                        }
                                        
                                        
                                        else /*means model_type == "lambda_III" */{
                                            
                                            //fill constraint indices
                                            for (int i=0; i<(num_partns1-2)*2+2; ++i){
                                                const_idx1.add(to_string(i+1));
                                            }
                                            
                                            //fill constraint indices
                                            for (int i=0; i<(num_partns2-2)*2+2; ++i){
                                                const_idx2.add(to_string(i+1));
                                            }
                                            
                                            // Lambda coefficient matrix when linking with partition variables
                                            if(num_partns1 > 1) lambda_coef1.in(indices(inst_partition_lambda, const_idx1));
                                            if(num_partns2 > 1) lambda_coef2.in(indices(inst_partition_lambda, const_idx2));
                                            
                                            // Partition coefficient matrix when linking with lambda variables
                                            if(num_partns1 > 1) on_coef1.in(indices(unique_ids, partns1, const_idx1));
                                            if(num_partns2 > 1) on_coef2.in(indices(unique_ids, partns2, const_idx2));
                                            
                                            // fill lambda_coef1 and lambda_coef2
                                            for (size_t inst = 0; inst< nb_ins; inst++){
                                                auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                                auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                                if(num_partns1 > 1) {
                                                    for (int j=0; j<num_partns2+1; ++j) {
                                                        string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(1);
                                                        lambda_coef1.set_val(cur_idx,1);
                                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(num_partns1+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string((num_partns1-2)*2+2);
                                                        lambda_coef1.set_val(cur_idx,1);
                                                    }
                                                    
                                                    for (int i=1 ; i<(num_partns1-2)*2+1; i=i+2) {
                                                        for (int j=(i-1)/2 + 2; j<num_partns1+1; ++j) {
                                                            for(int k=0; k<num_partns2+1; ++k){
                                                                string cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+name2+"{"+to_string(k+1)+"},"+to_string(i+1);
                                                                lambda_coef1.set_val(cur_idx,1);
                                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+name2+"{"+to_string(k+1)+"},"+to_string(i+2);
                                                                lambda_coef1.set_val(cur_idx,-1);
                                                            }
                                                        }
                                                    }
                                                }
                                                
                                                if(num_partns2 > 1){
                                                    for (int i=0; i<num_partns1+1; ++i) {
                                                        string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(1)+"},"+to_string(1);
                                                        lambda_coef2.set_val(cur_idx,1);
                                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(num_partns2+1)+"},"+to_string((num_partns2-2)*2+2);
                                                        lambda_coef2.set_val(cur_idx,1);
                                                    }
                                                    
                                                    for (int i=1 ; i<(num_partns2-2)*2+1; i=i+2) {
                                                        for (int j=(i-1)/2 + 2; j<num_partns2+1; ++j) {
                                                            for(int k=0; k<num_partns1+1; ++k){
                                                                string cur_idx = cur_var_idx+","+name1+"{"+to_string(k+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                                lambda_coef2.set_val(cur_idx,1);
                                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(k+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                                lambda_coef2.set_val(cur_idx,-1);
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                            
                                            
                                            
                                            // fill on_coef1 and on_coef2
                                            for (size_t inst = 0; inst< nb_ins; inst++){
                                                auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                                auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                                string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+to_string(1);
                                                if(num_partns1 > 1) {
                                                    on_coef1.set_val(cur_idx,1);
                                                    
                                                    
                                                    for (int i=1; i<num_partns1; ++i) {
                                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+to_string(2);
                                                        on_coef1.set_val(cur_idx, 1);
                                                    }
                                                    
                                                    for (int i=2 ; i<(num_partns1-2)*2+2; i=i+2) {
                                                        for (int j=i/2+1; j<num_partns1; ++j) {
                                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                            on_coef1.set_val(cur_idx,-1);
                                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                            on_coef1.set_val(cur_idx,1);
                                                        }
                                                    }
                                                }
                                                
                                                if(num_partns2 > 1) {
                                                    cur_idx = cur_var_idx+","+name2+"{"+to_string(1)+"},"+to_string(1);
                                                    on_coef2.set_val(cur_idx,1);
                                                    for (int i=1; i<num_partns2; ++i) {
                                                        cur_idx = cur_var_idx+","+name2+"{"+to_string(i+1)+"},"+to_string(2);
                                                        on_coef2.set_val(cur_idx, 1);
                                                    }
                                                    
                                                    for (int i=2 ; i<(num_partns2-2)*2+2; i=i+2) {
                                                        for (int j=i/2+1; j<num_partns2; ++j) {
                                                            cur_idx = cur_var_idx+","+name2+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                            on_coef2.set_val(cur_idx,-1);
                                                            cur_idx = cur_var_idx+","+name2+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                            on_coef2.set_val(cur_idx,1);
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                        
                                        
                                        /** Constraints */
                                        // Representation of the bilinear term with convex combination
                                        Constraint<> bln_rep(pair.first+"_bln_rep");
                                        bln_rep = EP.in_matrix(nb_entries,total_entries-nb_entries)*lambda.in_matrix(nb_entries,total_entries-nb_entries) - vlift.in(unique_ids);
                                        add(bln_rep.in(unique_ids) == 0);
                                        
                                        // Representation of o1 with convex combination
                                        Constraint<> o1_rep(pair.first+"_o1_rep");
                                        o1_rep = bounds1.from_ith(0,inst_partition_lambda).in_matrix(nb_entries, 1) * lambda.in_matrix(nb_entries,total_entries-nb_entries) - o1.in(o1_ids);
                                        add(o1_rep.in(unique_ids) == 0);
                                        
                                        // Representation of o2 with convex combination
                                        Constraint<> o2_rep(pair.first+"_o2_rep");

                                        o2_rep = bounds2.in_ignore_ith(nb_entries, 1, inst_partition_lambda).in_matrix(nb_entries,1) * lambda.in_matrix(nb_entries,total_entries-nb_entries) - o2.in(o2_ids);

                                        add(o2_rep.in(unique_ids) == 0);
                                        
                                        // Linking partition variables1 with lambda
                                        if(model_type == "lambda_II"){
                                            if(num_partns1 > 1) {
                                                Constraint<> on_link_lambda1(pair.first+"_on_link_lambda1_II");
 
                                                on_link_lambda1 = lambda.in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef1.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef1.in_matrix(nb_entries,total_entries-nb_entries) - binvar1->in_ignore_ith(nb_entries_v1,nb_entries_v2,on_coef1.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef1.in_matrix(nb_entries,1);
                                                add(on_link_lambda1.in(indices(unique_ids,const_idx1)) <= 0);
                                            }
                                            if(num_partns2 > 1) {
                                                Constraint<> on_link_lambda2(pair.first+"_on_link_lambda2_II");
                                              on_link_lambda2 = lambda.in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef2.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef2.in_matrix(nb_entries,total_entries-nb_entries) - on2.in_ignore_ith(0,nb_entries_v1,on_coef2.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef2.in_matrix(nb_entries,1);
                                                add(on_link_lambda2.in(indices(unique_ids,const_idx2)) <= 0);

                                            }
                                        }
                                        else{
                                            if(num_partns1 > 1) {
                                                Constraint<> on_link_lambda1(pair.first+"_on_link_lambda1_III");

                                                on_link_lambda1 = lambda.in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef1.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef1.in_matrix(nb_entries,total_entries-nb_entries) - binvar1->in_ignore_ith(nb_entries_v1,nb_entries_v2,on_coef1.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef1.in_matrix(nb_entries,1);
                                                add(on_link_lambda1.in(indices(unique_ids,const_idx1)) <= 0);
                                                

                                            }
                                            if(num_partns2 > 1) {
                                                Constraint<> on_link_lambda2(pair.first+"_on_link_lambda2_III");

                                                on_link_lambda2 = lambda.in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef2.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef2.in_matrix(nb_entries,total_entries-nb_entries) - on2.in_ignore_ith(0,nb_entries_v1,on_coef2.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef2.in_matrix(nb_entries,1);
                                                add(on_link_lambda2.in(indices(unique_ids,const_idx2)) <= 0);
                                            }
                                        }
                                        // sum over lambda
                                        Constraint<> lambdaSum(pair.first+"_lambdaSum");
                                        lambdaSum = sum(lambda.in_matrix(nb_entries,total_entries-nb_entries));
                                        add(lambdaSum.in(unique_ids) == 1);
                                    }
                                }
                            }
                            
                            else if(binvar_ptr2 !=_vars_name.end()){ //means v2 has been partitioned before)
                                DebugOn("<<<<<<<<<< THIS IS NOT SEEN BOTH -> SEEN SECOND -> DIFF VARS <<<<<<<<<<<" << endl);
                                
                                var<int> on1(name1+"_binary",0,1);
                                indices partns1("partns1");
                                for (int i = 0; i < num_partns1 ; ++i)
                                {
                                    partns1.add(name1+ "{" +to_string(i+1) + "}");
                                }
                                add(on1.in(o1_ids_uq,partns1));
                                
                                indices partns2("partns2");
                                for (int i = 0; i < num_partns2 ; ++i)
                                {
                                    partns2.add(name2+ "{" + to_string(i+1) + "}");
                                }
                                
                                auto binvar2 = static_pointer_cast<var<int>>(binvar_ptr2->second);
                                param<int> lb2("lb2"), ub2("ub2");
                                lb2.in(o2_ids_uq,partns2);
                                ub2.in(o2_ids_uq,partns2);
                                lb2.set_val(0), ub2.set_val(1);
                                auto added2 = binvar2->add_bounds(lb2,ub2);
                                reindex_vars();
                                
                                auto nb_entries_v1 = o1_ids.get_nb_entries();
                                auto nb_entries_v2 = o2_ids.get_nb_entries();
                                auto nb_entries = unique_ids.get_nb_entries();
                                
                                if(!added2.empty()){
                                    Constraint<> onSum2(o2._name+"_binarySum");
                                    onSum2 = sum(binvar2->in(added2).in_matrix(nb_entries_v2,1));
                                    auto vset2 = added2.from_ith(0,nb_entries_v2);
                                    vset2.filter_refs(vset2.get_unique_refs());
                                    add(onSum2.in(vset2) == 1);
                                }
                                
                                Constraint<> onSum1(o1._name+"_binarySum");
                                onSum1 = sum(on1.in_matrix(nb_entries_v1,1));
                                add(onSum1.in(o1_ids_uq) == 1);
                                
                                if(model_type == "on/off"){//if on/off is chosen
                                    var<int> on(name1+name2+"_binary",0,1);
                                    
                                    indices partns("partns");
                                    partns = indices(partns1,partns2);
                                    auto inst_partition = indices(unique_ids,partns);
                                    add(on.in(inst_partition));
                                    auto total_entries = inst_partition.get_nb_entries();
                                    
                                    Constraint<> onLink1(pair.first+"_binaryLink1");
                                    onLink1 = on1.from_ith(0,inst_partition.ignore_ith(nb_entries_v1, nb_entries_v2)) - on;
                                    add(onLink1.in(inst_partition) >= 0);
                                    
                                    Constraint<> onLink2(pair.first+"_binaryLink2");
                                    onLink2 = binvar2->in_ignore_ith(nb_entries_v2,1,inst_partition.ignore_ith(0,nb_entries_v1)) - on;
                                    add(onLink2.in(inst_partition) >= 0);
                                    
                                    Constraint<> onLink3(pair.first+"_binaryLink3");
                                    onLink3 = on1.from_ith(0,inst_partition.ignore_ith(nb_entries_v1, nb_entries_v2)) + binvar2->in_ignore_ith(nb_entries_v2,1,inst_partition.ignore_ith(0,nb_entries_v1)) - 1 - on;
                                    add(onLink3.in(inst_partition) <= 0);
                                    
                                    Constraint<> onSumComb(pair.first+"_binarySum");
                                    onSumComb = sum(on.in_matrix(nb_entries,total_entries-nb_entries));
                                    add(onSumComb.in(unique_ids) == 1);
                                    
                                    add_on_off_McCormick_refined(pair.first, vlift.in(unique_ids), o1.in(o1_ids), o2.in(o2_ids), on);
                                }
                                
                                
                                else{ //means it is one of the lambda formulations
                                    
                                    //difference is this has one more partition index
                                    indices partns1_lambda("partns1_lambda");
                                    for (int i = 0; i < num_partns1+1; ++i)
                                    {
                                        partns1_lambda.add(name1+ "{" +to_string(i+1) + "}");
                                    }
                                    
                                    indices partns2_lambda("partns2_lambda");
                                    for (int i = 0; i < num_partns2+1; ++i)
                                    {
                                        partns2_lambda.add(name2+ "{" +to_string(i+1) + "}");
                                    }
                                    
                                    indices partns_lambda("partns_lambda");
                                    partns_lambda = indices(partns1_lambda,partns2_lambda);
                                    auto inst_partition_lambda = indices(unique_ids,partns_lambda);
                                    auto inst_partition_bounds1 = indices(unique_ids,partns1_lambda);
                                    auto inst_partition_bounds2 = indices(unique_ids,partns2_lambda);
                                    
                                    // Convex combination variables
                                    var<> lambda(name1+name2+"_lambda",pos_);
                                    add(lambda.in(inst_partition_lambda));
                                    
                                    /** Parameters */
                                    // Bounds on variable v1 & v2
                                    param<> bounds1(name1+"_bounds1");
                                    bounds1.in(inst_partition_bounds1);
                                    
                                    param<> bounds2(name2+"_bounds2");
                                    bounds2.in(inst_partition_bounds2);
                                    
                                    // Function values on the extreme points
                                    param<> EP(name1+name2+"_grid_values");
                                    EP.in(inst_partition_lambda);
                                    auto total_entries = inst_partition_lambda.get_nb_entries();
                                    
                                    size_t nb_ins = vlift.in(unique_ids).get_nb_inst();
                                    auto o1_global_lb = o1.get_lb();
                                    auto increment1 = (o1.get_ub() - o1_global_lb)/num_partns1;
                                    
                                    auto o2_global_lb = o2.get_lb();
                                    auto increment2 = (o2.get_ub() - o2_global_lb)/num_partns2;
                                    
                                    // fill bounds1 and function values
                                    for (int i=0 ; i<num_partns1+1; ++i) {
                                        auto bound_partn1 = o1_global_lb + increment1*i;
                                        bound_partn1.eval_all();
                                        for (size_t inst = 0; inst< nb_ins; inst++){
                                            auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                            auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                            string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"}";
                                            bounds1.set_val(cur_idx,bound_partn1.eval(inst));
                                            for(int j=0; j<num_partns2+1; ++j){
                                                auto bound_partn2 = o2_global_lb + increment2*j;
                                                bound_partn2.eval_all();
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"}";
                                                EP.set_val(cur_idx,(bound_partn1.eval(inst)*bound_partn2.eval(inst)));
                                            }
                                        }
                                    }
                                    // fill bounds2
                                    for (int i=0 ; i<num_partns2+1; ++i) {
                                        auto bound_partn2 = o2_global_lb + increment2*i;
                                        bound_partn2.eval_all();
                                        for (size_t inst = 0; inst< nb_ins; inst++){
                                            auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                            auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                            string cur_idx = cur_var_idx+","+name2+"{"+to_string(i+1)+"}";
                                            bounds2.set_val(cur_idx,bound_partn2.eval(inst));
                                        }
                                    }
                                    
                                    // Lambda coefficient matrix when linking with partition variables
                                    param<> lambda_coef1(name1+"_lambda_linking_coefficients1");
                                    param<> lambda_coef2(name2+"_lambda_linking_coefficients2");
                                    
                                    // Partition coefficient matrix when linking with lambda variables
                                    param<> on_coef1(name1+"_partition_linking_coefficients1");
                                    param<> on_coef2(name2+"_partition_linking_coefficients2");
                                    
                                    // create constraint indices
                                    indices const_idx1("const_idx1");
                                    indices const_idx2("const_idx2");
                                    
                                    if(model_type == "lambda_II"){
                                        
                                        //fill constraint indices
                                        for (int i=0; i<num_partns1+1; ++i){
                                            const_idx1.add(to_string(i+1));
                                        }
                                        
                                        //fill constraint indices
                                        for (int i=0; i<num_partns2+1; ++i){
                                            const_idx2.add(to_string(i+1));
                                        }
                                        
                                        // Lambda coefficient matrix when linking with partition variables
                                        if(num_partns1 > 1) lambda_coef1.in(indices(inst_partition_lambda, const_idx1));
                                        if(num_partns2 > 1) lambda_coef2.in(indices(inst_partition_lambda, const_idx2));
                                        
                                        // Lambda coefficient matrix when linking with partition variables
                                        if(num_partns1 > 1)  on_coef1.in(indices(unique_ids, partns1, const_idx1));
                                        if(num_partns2 > 1)  on_coef2.in(indices(unique_ids, partns2, const_idx2));
                                        
                                        // fill lambda_coef1 and lambda_coef2
                                        for (size_t inst = 0; inst< nb_ins; inst++){
                                            auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                            auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                            for (int i=0 ; i<num_partns1+1; ++i) {
                                                for (int j=0 ; j<num_partns2+1; ++j) {
                                                    string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                    if(num_partns1 > 1)  lambda_coef1.set_val(cur_idx,1);
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(j+1);
                                                    if(num_partns2 > 1)  lambda_coef2.set_val(cur_idx,1);
                                                }
                                            }
                                        }
                                        
                                        // fill on_coef1 and on_coef2
                                        for (size_t inst = 0; inst< nb_ins; inst++){
                                            auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                            auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                            string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+to_string(1);
                                            if(num_partns1 > 1)  on_coef1.set_val(cur_idx,1);
                                            cur_idx = cur_var_idx+","+name2+"{"+to_string(1)+"},"+to_string(1);
                                            if(num_partns2 > 1)  on_coef2.set_val(cur_idx,1);
                                            
                                            if(num_partns1 > 1) {
                                                for (int i=1 ; i<num_partns1; ++i) {
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(i)+"},"+to_string(i+1);
                                                    on_coef1.set_val(cur_idx,1);
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+to_string(i+1);
                                                    on_coef1.set_val(cur_idx,1);
                                                }
                                            }
                                            if(num_partns2 > 1) {
                                                for (int i=1 ; i<num_partns2; ++i) {
                                                    cur_idx = cur_var_idx+","+name2+"{"+to_string(i)+"},"+to_string(i+1);
                                                    on_coef2.set_val(cur_idx,1);
                                                    cur_idx = cur_var_idx+","+name2+"{"+to_string(i+1)+"},"+to_string(i+1);
                                                    on_coef2.set_val(cur_idx,1);
                                                }
                                            }
                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(num_partns1)+"},"+to_string(num_partns1+1);
                                            if(num_partns1 > 1)  on_coef1.set_val(cur_idx,1);
                                            cur_idx = cur_var_idx+","+name2+"{"+to_string(num_partns2)+"},"+to_string(num_partns2+1);
                                            if(num_partns1 > 2)  on_coef2.set_val(cur_idx,1);
                                        }
                                    }
                                    
                                    
                                    else /*means model_type == "lambda_III" */{
                                        
                                        //fill constraint indices
                                        for (int i=0; i<(num_partns1-2)*2+2; ++i){
                                            const_idx1.add(to_string(i+1));
                                        }
                                        
                                        //fill constraint indices
                                        for (int i=0; i<(num_partns2-2)*2+2; ++i){
                                            const_idx2.add(to_string(i+1));
                                        }
                                        
                                        // Lambda coefficient matrix when linking with partition variables
                                        if(num_partns1 > 1)  lambda_coef1.in(indices(inst_partition_lambda, const_idx1));
                                        if(num_partns2 > 1)  lambda_coef2.in(indices(inst_partition_lambda, const_idx2));
                                        
                                        // Partition coefficient matrix when linking with lambda variables
                                        if(num_partns1 > 1)  on_coef1.in(indices(unique_ids, partns1, const_idx1));
                                        if(num_partns2 > 1)  on_coef2.in(indices(unique_ids, partns2, const_idx2));
                                        
                                        // fill lambda_coef1 and lambda_coef2
                                        for (size_t inst = 0; inst< nb_ins; inst++){
                                            auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                            auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                            if(num_partns1 > 1){
                                                for (int j=0; j<num_partns2+1; ++j) {
                                                    string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(1);
                                                    lambda_coef1.set_val(cur_idx,1);
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(num_partns1+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string((num_partns1-2)*2+2);
                                                    lambda_coef1.set_val(cur_idx,1);
                                                }
                                                
                                                for (int i=1 ; i<(num_partns1-2)*2+1; i=i+2) {
                                                    for (int j=(i-1)/2 + 2; j<num_partns1+1; ++j) {
                                                        for(int k=0; k<num_partns2+1; ++k){
                                                            string cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+name2+"{"+to_string(k+1)+"},"+to_string(i+1);
                                                            lambda_coef1.set_val(cur_idx,1);
                                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+name2+"{"+to_string(k+1)+"},"+to_string(i+2);
                                                            lambda_coef1.set_val(cur_idx,-1);
                                                        }
                                                    }
                                                }
                                            }
                                            if(num_partns2 > 1) {
                                                for (int i=0; i<num_partns1+1; ++i) {
                                                    string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(1)+"},"+to_string(1);
                                                    lambda_coef2.set_val(cur_idx,1);
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(num_partns2+1)+"},"+to_string((num_partns2-2)*2+2);
                                                    lambda_coef2.set_val(cur_idx,1);
                                                }
                                                
                                                for (int i=1 ; i<(num_partns2-2)*2+1; i=i+2) {
                                                    for (int j=(i-1)/2 + 2; j<num_partns2+1; ++j) {
                                                        for(int k=0; k<num_partns1+1; ++k){
                                                            string cur_idx = cur_var_idx+","+name1+"{"+to_string(k+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                            lambda_coef2.set_val(cur_idx,1);
                                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(k+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                            lambda_coef2.set_val(cur_idx,-1);
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                        
                                        
                                        
                                        // fill on_coef1 and on_coef2
                                        for (size_t inst = 0; inst< nb_ins; inst++){
                                            auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                            auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                            string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+to_string(1);
                                            if(num_partns1 > 1) {
                                                on_coef1.set_val(cur_idx,1);
                                                
                                                for (int i=1; i<num_partns1; ++i) {
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+to_string(2);
                                                    on_coef1.set_val(cur_idx, 1);
                                                }
                                                
                                                for (int i=2 ; i<(num_partns1-2)*2+2; i=i+2) {
                                                    for (int j=i/2+1; j<num_partns1; ++j) {
                                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                        on_coef1.set_val(cur_idx,-1);
                                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                        on_coef1.set_val(cur_idx,1);
                                                    }
                                                }
                                            }
                                            
                                            if(num_partns2 > 1) {
                                                cur_idx = cur_var_idx+","+name2+"{"+to_string(1)+"},"+to_string(1);
                                                on_coef2.set_val(cur_idx,1);
                                                for (int i=1; i<num_partns2; ++i) {
                                                    cur_idx = cur_var_idx+","+name2+"{"+to_string(i+1)+"},"+to_string(2);
                                                    on_coef2.set_val(cur_idx, 1);
                                                }
                                                
                                                for (int i=2 ; i<(num_partns2-2)*2+2; i=i+2) {
                                                    for (int j=i/2+1; j<num_partns2; ++j) {
                                                        cur_idx = cur_var_idx+","+name2+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                        on_coef2.set_val(cur_idx,-1);
                                                        cur_idx = cur_var_idx+","+name2+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                        on_coef2.set_val(cur_idx,1);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    
                                    
                                    /** Constraints */
                                    // Representation of the bilinear term with convex combination
                                    Constraint<> bln_rep(pair.first+"_bln_rep");
                                    bln_rep = EP.in_matrix(nb_entries,total_entries-nb_entries)*lambda.in_matrix(nb_entries,total_entries-nb_entries) - vlift.in(unique_ids);
                                    add(bln_rep.in(unique_ids) == 0);
                                    
                                    // Representation of o1 with convex combination
                                    Constraint<> o1_rep(pair.first+"_o1_rep");

                                    o1_rep = bounds1.from_ith(0,inst_partition_lambda).in_matrix(nb_entries, 1) * lambda.in_matrix(nb_entries,total_entries-nb_entries) - o1.in(o1_ids);
                                    add(o1_rep.in(unique_ids) == 0);
                                    
                                    // Representation of o2 with convex combination
                                    Constraint<> o2_rep(pair.first+"_o2_rep");
 
                                    o2_rep = bounds2.in_ignore_ith(nb_entries, 1, inst_partition_lambda).in_matrix(nb_entries,1) * lambda.in_matrix(nb_entries,total_entries-nb_entries) - o2.in(o2_ids);
                                    add(o2_rep.in(unique_ids) == 0);
                                    
                                    // Linking partition variables1 with lambda
                                    if(model_type == "lambda_II"){
                                        if(num_partns1 > 1) {
                                            Constraint<> on_link_lambda1(pair.first+"_on_link_lambda1_II");
                                            on_link_lambda1 = lambda.in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef1.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef1.in_matrix(nb_entries,total_entries-nb_entries) - on1.in_ignore_ith(nb_entries_v1,nb_entries_v2,on_coef1.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef1.in_matrix(nb_entries,1);
                                            add(on_link_lambda1.in(indices(unique_ids,const_idx1)) <= 0);
                                        }
                                        if(num_partns2 > 1) {
                                            Constraint<> on_link_lambda2(pair.first+"_on_link_lambda2_II");
                                            on_link_lambda2 = lambda.in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef2.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef2.in_matrix(nb_entries,total_entries-nb_entries) - binvar2->in_ignore_ith(0,nb_entries_v1,on_coef2.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef2.in_matrix(nb_entries,1);
                                            add(on_link_lambda2.in(indices(unique_ids,const_idx2)) <= 0);

                                        }
                                    }
                                    else{
                                        if(num_partns1 > 1) {
                                            Constraint<> on_link_lambda1(pair.first+"_on_link_lambda1_III");
                                            on_link_lambda1 = lambda.in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef1.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef1.in_matrix(nb_entries,total_entries-nb_entries) - on1.in_ignore_ith(nb_entries_v1,nb_entries_v2,on_coef1.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef1.in_matrix(nb_entries,1);
                                            add(on_link_lambda1.in(indices(unique_ids,const_idx1)) <= 0);
                                        }
                                        if(num_partns2 > 1) {
                                            Constraint<> on_link_lambda2(pair.first+"_on_link_lambda2_III");
                                            on_link_lambda2 = lambda.in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef2.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef2.in_matrix(nb_entries,total_entries-nb_entries) - binvar2->in_ignore_ith(0,nb_entries_v1,on_coef2.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef2.in_matrix(nb_entries,1);
                                            add(on_link_lambda2.in(indices(unique_ids,const_idx2)) <= 0);
                                        }
                                    }
                                    // sum over lambda
                                    Constraint<> lambdaSum(pair.first+"_lambdaSum");
                                    lambdaSum = sum(lambda.in_matrix(nb_entries,total_entries-nb_entries));
                                    add(lambdaSum.in(unique_ids) == 1);
                                }
                                
                            }
                            else{ //means both variables v1 and v2 haven't been partitioned
                                if(name1==name2){
                                    DebugOn("<<<<<<<<<< THIS IS NOT SEEN BOTH -> DOUBLE -> SAME VARS <<<<<<<<<<<" << endl);
                                    
                                    var<int> on1(name1+"_binary",0,1);
                                    indices partns1("partns1");
                                    for (int i = 0; i < num_partns1 ; ++i)
                                    {
                                        partns1.add(name1+ "{" +to_string(i+1) + "}");
                                    }
                                    add(on1.in(union_ids(o1_ids_uq, o2_ids_uq),partns1));
                                    indices partns2("partns2");
                                    for (int i = 0; i < num_partns2 ; ++i)
                                    {
                                        partns2.add(name2+ "{" +to_string(i+1) + "}");
                                    }
                                    
                                    auto nb_entries_v1 = o1_ids.get_nb_entries();
                                    auto nb_entries_v2 = o2_ids.get_nb_entries();
                                    auto nb_entries = unique_ids.get_nb_entries();
                                    
                                    Constraint<> onSum1(o1._name+"_binarySum");
                                    onSum1 = sum(on1.in_matrix(nb_entries_v1,1));
                                    add(onSum1.in(union_ids(o1_ids_uq,o2_ids_uq)) == 1);
                                    
                                    if(model_type == "on/off")
                                    {
                                        var<int> on(name1+name2+"_binary",0,1);
                                        
                                        indices partns("partns");
                                        partns = indices(partns1,partns2);
                                        auto inst_partition = indices(unique_ids,partns);
                                        add(on.in(inst_partition));
                                        auto total_entries = inst_partition.get_nb_entries();
                                        
                                        Constraint<> onLink1(pair.first+"_binaryLink1");
                                        onLink1 = on1.from_ith(0,inst_partition.ignore_ith(nb_entries_v1, nb_entries_v2)) - on;
                                        add(onLink1.in(inst_partition) >= 0);
                                        
                                        Constraint<> onLink2(pair.first+"_binaryLink2");
                                        onLink2 = on1.in_ignore_ith(nb_entries_v1,1,inst_partition.ignore_ith(0,nb_entries_v1)) - on;
                                        add(onLink2.in(inst_partition) >= 0);
                                        
                                        Constraint<> onLink3(pair.first+"_binaryLink3");
                                        onLink3 = on1.from_ith(0,inst_partition.ignore_ith(nb_entries_v1, nb_entries_v2)) + on1.in_ignore_ith(nb_entries_v1,1,inst_partition.ignore_ith(0,nb_entries_v1)) - 1 - on;
                                        add(onLink3.in(inst_partition) <= 0);
                                        
                                        Constraint<> onSumComb(pair.first+"_binarySum");
                                        onSumComb = sum(on.in_matrix(nb_entries,total_entries-nb_entries));
                                        add(onSumComb.in(unique_ids) == 1);
                                        
                                        add_on_off_McCormick_refined(pair.first, vlift.in(unique_ids), o1.in(o1_ids), o2.in(o2_ids), on);
                                    }
                                    
                                    else{ //means it is one of the lambda formulations
                                        
                                        //difference is this has one more partition index
                                        indices partns1_lambda("partns1_lambda");
                                        for (int i = 0; i < num_partns1+1; ++i)
                                        {
                                            partns1_lambda.add(name1+ "{" +to_string(i+1) + "}");
                                        }
                                        
                                        indices partns2_lambda("partns2_lambda");
                                        for (int i = 0; i < num_partns2+1; ++i)
                                        {
                                            partns2_lambda.add(name2+ "{" +to_string(i+1) + "}");
                                        }
                                        
                                        indices partns_lambda("partns_lambda");
                                        partns_lambda = indices(partns1_lambda,partns2_lambda);
                                        auto inst_partition_lambda = indices(unique_ids,partns_lambda);
                                        auto inst_partition_bounds1 = indices(unique_ids,partns1_lambda);
                                        auto inst_partition_bounds2 = indices(unique_ids,partns2_lambda);
                                        
                                        // Convex combination variables
                                        var<> lambda(name1+name2+"_lambda",pos_);
                                        add(lambda.in(inst_partition_lambda));
                                        
                                        /** Parameters */
                                        // Bounds on variable v1 & v2
                                        param<> bounds1(name1+"_bounds1");
                                        bounds1.in(inst_partition_bounds1);
                                        
                                        param<> bounds2(name2+"_bounds2");
                                        bounds2.in(inst_partition_bounds2);
                                        
                                        // Function values on the extreme points
                                        param<> EP(name1+name2+"_grid_values");
                                        EP.in(inst_partition_lambda);
                                        auto total_entries = inst_partition_lambda.get_nb_entries();
                                        
                                        size_t nb_ins = vlift.in(unique_ids).get_nb_inst();
                                        
                                        auto o1_global_lb = o1.get_lb();
                                        auto increment1 = (o1.get_ub() - o1_global_lb)/num_partns1;
                                        
                                        auto o2_global_lb = o2.get_lb();
                                        auto increment2 = (o2.get_ub() - o2_global_lb)/num_partns2;
                                        
                                        // fill bounds and function values
                                        for (int i=0 ; i<num_partns1+1; ++i) {
                                            auto bound_partn1 = o1_global_lb + increment1*i;
                                            bound_partn1.eval_all();
                                            auto bound_partn2 = o2_global_lb + increment2*i;
                                            bound_partn2.eval_all();
                                            for (size_t inst = 0; inst< nb_ins; inst++){
                                                auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                                auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                                string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"}";
                                                bounds1.set_val(cur_idx,bound_partn1.eval(inst));
                                                bounds2.set_val(cur_idx,bound_partn2.eval(inst));
                                                for(int j=0; j<num_partns2+1; ++j){
                                                    auto bound_partn2_temp = o2_global_lb + increment2*j;
                                                    bound_partn2_temp.eval_all();
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"}";
                                                    EP.set_val(cur_idx,(bound_partn1.eval(inst)*bound_partn2_temp.eval(inst)));
                                                }
                                            }
                                        }
                                        
                                        // Lambda coefficient matrix when linking with partition variables
                                        param<> lambda_coef1(name1+"_lambda_linking_coefficients1");
                                        param<> lambda_coef2(name2+"_lambda_linking_coefficients2");
                                        
                                        // Partition coefficient matrix when linking with lambda variables
                                        param<> on_coef1(name1+"_partition_linking_coefficients1");
                                        param<> on_coef2(name2+"_partition_linking_coefficients2");
                                        
                                        // create constraint indices
                                        indices const_idx1("const_idx1");
                                        indices const_idx2("const_idx2");
                                        
                                        if(model_type == "lambda_II"){
                                            
                                            //fill constraint indices
                                            for (int i=0; i<num_partns1+1; ++i){
                                                const_idx1.add(to_string(i+1));
                                            }
                                            
                                            //fill constraint indices
                                            for (int i=0; i<num_partns2+1; ++i){
                                                const_idx2.add(to_string(i+1));
                                            }
                                            
                                            // Lambda coefficient matrix when linking with partition variables
                                            lambda_coef1.in(indices(inst_partition_lambda, const_idx1));
                                            lambda_coef2.in(indices(inst_partition_lambda, const_idx2));
                                            
                                            // Lambda coefficient matrix when linking with partition variables
                                            on_coef1.in(indices(unique_ids, partns1, const_idx1));
                                            on_coef2.in(indices(unique_ids, partns2, const_idx2));
                                            
                                            // fill lambda_coef1 and lambda_coef2
                                            for (size_t inst = 0; inst< nb_ins; inst++){
                                                auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                                auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                                for (int i=0 ; i<num_partns1+1; ++i) {
                                                    for (int j=0 ; j<num_partns2+1; ++j) {
                                                        string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                        lambda_coef1.set_val(cur_idx,1);
                                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(j+1);
                                                        lambda_coef2.set_val(cur_idx,1);
                                                    }
                                                }
                                            }
                                            
                                            // fill on_coef1 and on_coef2
                                            for (size_t inst = 0; inst< nb_ins; inst++){
                                                auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                                auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                                string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+to_string(1);
                                                on_coef1.set_val(cur_idx,1);
                                                cur_idx = cur_var_idx+","+name2+"{"+to_string(1)+"},"+to_string(1);
                                                on_coef2.set_val(cur_idx,1);
                                                for (int i=1 ; i<num_partns1; ++i) {
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(i)+"},"+to_string(i+1);
                                                    on_coef1.set_val(cur_idx,1);
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+to_string(i+1);
                                                    on_coef1.set_val(cur_idx,1);
                                                }
                                                for (int i=1 ; i<num_partns2; ++i) {
                                                    cur_idx = cur_var_idx+","+name2+"{"+to_string(i)+"},"+to_string(i+1);
                                                    on_coef2.set_val(cur_idx,1);
                                                    cur_idx = cur_var_idx+","+name2+"{"+to_string(i+1)+"},"+to_string(i+1);
                                                    on_coef2.set_val(cur_idx,1);
                                                }
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(num_partns1)+"},"+to_string(num_partns1+1);
                                                on_coef1.set_val(cur_idx,1);
                                                cur_idx = cur_var_idx+","+name2+"{"+to_string(num_partns2)+"},"+to_string(num_partns2+1);
                                                on_coef2.set_val(cur_idx,1);
                                            }
                                        }
                                        
                                        
                                        else /*means model_type == "lambda_III" */{
                                            
                                            //fill constraint indices
                                            for (int i=0; i<(num_partns1-2)*2+2; ++i){
                                                const_idx1.add(to_string(i+1));
                                            }
                                            
                                            //fill constraint indices
                                            for (int i=0; i<(num_partns2-2)*2+2; ++i){
                                                const_idx2.add(to_string(i+1));
                                            }
                                            
                                            // Lambda coefficient matrix when linking with partition variables
                                            lambda_coef1.in(indices(inst_partition_lambda, const_idx1));
                                            lambda_coef2.in(indices(inst_partition_lambda, const_idx2));
                                            
                                            // Partition coefficient matrix when linking with lambda variables
                                            on_coef1.in(indices(unique_ids, partns1, const_idx1));
                                            on_coef2.in(indices(unique_ids, partns2, const_idx2));
                                            
                                            // fill lambda_coef1 and lambda_coef2
                                            for (size_t inst = 0; inst< nb_ins; inst++){
                                                auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                                auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                                for (int j=0; j<num_partns2+1; ++j) {
                                                    string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(1);
                                                    lambda_coef1.set_val(cur_idx,1);
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(num_partns1+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string((num_partns1-2)*2+2);
                                                    lambda_coef1.set_val(cur_idx,1);
                                                }
                                                
                                                for (int i=1 ; i<(num_partns1-2)*2+1; i=i+2) {
                                                    for (int j=(i-1)/2 + 2; j<num_partns1+1; ++j) {
                                                        for(int k=0; k<num_partns2+1; ++k){
                                                            string cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+name2+"{"+to_string(k+1)+"},"+to_string(i+1);
                                                            lambda_coef1.set_val(cur_idx,1);
                                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+name2+"{"+to_string(k+1)+"},"+to_string(i+2);
                                                            lambda_coef1.set_val(cur_idx,-1);
                                                        }
                                                    }
                                                }
                                                
                                                for (int i=0; i<num_partns1+1; ++i) {
                                                    string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(1)+"},"+to_string(1);
                                                    lambda_coef2.set_val(cur_idx,1);
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(num_partns2+1)+"},"+to_string((num_partns2-2)*2+2);
                                                    lambda_coef2.set_val(cur_idx,1);
                                                }
                                                
                                                for (int i=1 ; i<(num_partns2-2)*2+1; i=i+2) {
                                                    for (int j=(i-1)/2 + 2; j<num_partns2+1; ++j) {
                                                        for(int k=0; k<num_partns1+1; ++k){
                                                            string cur_idx = cur_var_idx+","+name1+"{"+to_string(k+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                            lambda_coef2.set_val(cur_idx,1);
                                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(k+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                            lambda_coef2.set_val(cur_idx,-1);
                                                        }
                                                    }
                                                }
                                            }
                                            
                                            
                                            
                                            // fill on_coef1 and on_coef2
                                            for (size_t inst = 0; inst< nb_ins; inst++){
                                                auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                                auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                                string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+to_string(1);
                                                on_coef1.set_val(cur_idx,1);
                                                
                                                for (int i=1; i<num_partns1; ++i) {
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+to_string(2);
                                                    on_coef1.set_val(cur_idx, 1);
                                                }
                                                
                                                for (int i=2 ; i<(num_partns1-2)*2+2; i=i+2) {
                                                    for (int j=i/2+1; j<num_partns1; ++j) {
                                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                        on_coef1.set_val(cur_idx,-1);
                                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                        on_coef1.set_val(cur_idx,1);
                                                    }
                                                }
                                                
                                                cur_idx = cur_var_idx+","+name2+"{"+to_string(1)+"},"+to_string(1);
                                                on_coef2.set_val(cur_idx,1);
                                                for (int i=1; i<num_partns2; ++i) {
                                                    cur_idx = cur_var_idx+","+name2+"{"+to_string(i+1)+"},"+to_string(2);
                                                    on_coef2.set_val(cur_idx, 1);
                                                }
                                                
                                                for (int i=2 ; i<(num_partns2-2)*2+2; i=i+2) {
                                                    for (int j=i/2+1; j<num_partns2; ++j) {
                                                        cur_idx = cur_var_idx+","+name2+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                        on_coef2.set_val(cur_idx,-1);
                                                        cur_idx = cur_var_idx+","+name2+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                        on_coef2.set_val(cur_idx,1);
                                                    }
                                                }
                                            }
                                        }
                                        
                                        
                                        /** Constraints */
                                        // Representation of the bilinear term with convex combination
                                        Constraint<> bln_rep(pair.first+"_bln_rep");
                                        bln_rep = EP.in_matrix(nb_entries,total_entries-nb_entries)*lambda.in_matrix(nb_entries,total_entries-nb_entries) - vlift.in(unique_ids);
                                        add(bln_rep.in(unique_ids) == 0);
                                        
                                        // Representation of o1 with convex combination
                                        Constraint<> o1_rep(pair.first+"_o1_rep");
                                        o1_rep = bounds1.from_ith(0,inst_partition_lambda).in_matrix(nb_entries, 1) * lambda.in_matrix(nb_entries,total_entries-nb_entries) - o1.in(o1_ids);
                                        add(o1_rep.in(unique_ids) == 0);
                                        
                                        // Representation of o2 with convex combination
                                        Constraint<> o2_rep(pair.first+"_o2_rep");
                                        o2_rep = bounds2.in_ignore_ith(nb_entries, 1, inst_partition_lambda).in_matrix(nb_entries,1) * lambda.in_matrix(nb_entries,total_entries-nb_entries) - o2.in(o2_ids);
                                        add(o2_rep.in(unique_ids) == 0);
                                        
                                        // Linking partition variables1 with lambda
                                        if(model_type == "lambda_II"){
                                            Constraint<> on_link_lambda1(pair.first+"_on_link_lambda1_II");
                                            on_link_lambda1 = lambda.in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef1.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef1.in_matrix(nb_entries,total_entries-nb_entries) - on1.in_ignore_ith(nb_entries_v1,nb_entries_v2,on_coef1.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef1.in_matrix(nb_entries,1);
                                            add(on_link_lambda1.in(indices(unique_ids,const_idx1)) <= 0);
                                            
                                            Constraint<> on_link_lambda2(pair.first+"_on_link_lambda2_II");
                                            on_link_lambda2 = lambda.in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef2.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef2.in_matrix(nb_entries,total_entries-nb_entries) - on1.in_ignore_ith(0,nb_entries_v1,on_coef2.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef2.in_matrix(nb_entries,1);
                                            add(on_link_lambda2.in(indices(unique_ids,const_idx2)) <= 0);
                                        }
                                        else{
                                            Constraint<> on_link_lambda1(pair.first+"_on_link_lambda1_III");
                                            on_link_lambda1 = lambda.in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef1.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef1.in_matrix(nb_entries,total_entries-nb_entries) - on1.in_ignore_ith(nb_entries_v1,nb_entries_v2,on_coef1.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef1.in_matrix(nb_entries,1);
                                            add(on_link_lambda1.in(indices(unique_ids,const_idx1)) <= 0);
                                            
                                            Constraint<> on_link_lambda2(pair.first+"_on_link_lambda2_III");
                                            on_link_lambda2 = lambda.in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef2.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef2.in_matrix(nb_entries,total_entries-nb_entries) - on1.in_ignore_ith(0,nb_entries_v1,on_coef2.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef2.in_matrix(nb_entries,1);
                                            add(on_link_lambda2.in(indices(unique_ids,const_idx2)) <= 0);
                                        }
                                        // sum over lambda
                                        Constraint<> lambdaSum(pair.first+"_lambdaSum");
                                        lambdaSum = sum(lambda.in_matrix(nb_entries,total_entries-nb_entries));
                                        add(lambdaSum.in(unique_ids) == 1);
                                    }
                                }
                                else{
                                    DebugOn("<<<<<<<<<< THIS IS NOT SEEN BOTH -> DOUBLE -> DIFF VARS <<<<<<<<<<<" << endl);
                                    
                                    var<int> on1(name1+"_binary",0,1);
                                    indices partns1("partns1");
                                    for (int i = 0; i < num_partns1 ; ++i)
                                    {
                                        partns1.add(name1+ "{" + to_string(i+1) + "}");
                                    }
                                    add(on1.in(o1_ids_uq,partns1));
                                    
                                    var<int> on2(name2+"_binary",0,1);
                                    indices partns2("partns2");
                                    for (int i = 0; i < num_partns2 ; ++i)
                                    {
                                        partns2.add(name2+ "{" + to_string(i+1) + "}");
                                    }
                                    add(on2.in(o2_ids_uq,partns2));
                                    
                                    auto nb_entries_v1 = o1_ids.get_nb_entries();
                                    auto nb_entries_v2 = o2_ids.get_nb_entries();
                                    auto nb_entries = unique_ids.get_nb_entries();
                                    
                                    Constraint<> onSum1(o1._name+"_binarySum");
                                    onSum1 = sum(on1.in_matrix(nb_entries_v1,1));
                                    add(onSum1.in(o1_ids_uq) == 1);
                                    
                                    Constraint<> onSum2(o2._name+"_binarySum");
                                    onSum2 = sum(on2.in_matrix(nb_entries_v2,1));
                                    add(onSum2.in(o2_ids_uq) == 1);
                                    
                                    if(model_type == "on/off"){//if on/off is chosen
                                        var<int> on(name1+name2+"_binary",0,1);
                                        
                                        indices partns("partns");
                                        partns = indices(partns1,partns2);
                                        auto inst_partition = indices(unique_ids,partns);
                                        add(on.in(inst_partition));
                                        auto total_entries = inst_partition.get_nb_entries();
                                        
                                        Constraint<> onLink1(pair.first+"_binaryLink1");
                                        onLink1 = on1.from_ith(0,inst_partition.ignore_ith(nb_entries_v1, nb_entries_v2)) - on;
                                        add(onLink1.in(inst_partition) >= 0);
                                        
                                        Constraint<> onLink2(pair.first+"_binaryLink2");
                                        onLink2 = on2.in_ignore_ith(nb_entries_v2,1,inst_partition.ignore_ith(0,nb_entries_v1)) - on;
                                        add(onLink2.in(inst_partition) >= 0);
                                        
                                        Constraint<> onLink3(pair.first+"_binaryLink3");
                                        onLink3 = on1.from_ith(0,inst_partition.ignore_ith(nb_entries_v1, nb_entries_v2)) + on2.in_ignore_ith(nb_entries_v2,1,inst_partition.ignore_ith(0,nb_entries_v1)) - 1 - on;
                                        add(onLink3.in(inst_partition) <= 0);
                                        
                                        Constraint<> onSumComb(pair.first+"_binarySum");
                                        onSumComb = sum(on.in_matrix(nb_entries,total_entries-nb_entries));
                                        add(onSumComb.in(unique_ids) == 1);
                                        
                                        add_on_off_McCormick_refined(pair.first, vlift.in(unique_ids), o1.in(o1_ids), o2.in(o2_ids), on);
                                    }
                                    
                                    else{//means it is one of the lambda formulations
                                        
                                        //difference is this has one more partition index
                                        indices partns1_lambda("partns1_lambda");
                                        for (int i = 0; i < num_partns1+1; ++i)
                                        {
                                            partns1_lambda.add(name1+ "{" +to_string(i+1) + "}");
                                        }
                                        
                                        indices partns2_lambda("partns2_lambda");
                                        for (int i = 0; i < num_partns2+1; ++i)
                                        {
                                            partns2_lambda.add(name2+ "{" +to_string(i+1) + "}");
                                        }
                                        
                                        indices partns_lambda("partns_lambda");
                                        partns_lambda = indices(partns1_lambda,partns2_lambda);
                                        auto inst_partition_lambda = indices(unique_ids,partns_lambda);
                                        auto inst_partition_bounds1 = indices(unique_ids,partns1_lambda);
                                        auto inst_partition_bounds2 = indices(unique_ids,partns2_lambda);
                                        
                                        // Convex combination variables
                                        var<> lambda(name1+name2+"_lambda",pos_);
                                        add(lambda.in(inst_partition_lambda));
                                        
                                        /** Parameters */
                                        // Bounds on variable v1 & v2
                                        param<> bounds1(name1+"_bounds1");
                                        bounds1.in(inst_partition_bounds1);
                                        
                                        param<> bounds2(name2+"_bounds2");
                                        bounds2.in(inst_partition_bounds2);
                                        
                                        // Function values on the extreme points
                                        param<> EP(name1+name2+"_grid_values");
                                        EP.in(inst_partition_lambda);
                                        auto total_entries = inst_partition_lambda.get_nb_entries();
                                        
                                        size_t nb_ins = vlift.in(unique_ids).get_nb_inst();
                                        auto o1_global_lb = o1.get_lb();
                                        auto increment1 = (o1.get_ub() - o1_global_lb)/num_partns1;
                                        
                                        auto o2_global_lb = o2.get_lb();
                                        auto increment2 = (o2.get_ub() - o2_global_lb)/num_partns2;
                                        
                                        // fill bounds1 and function values
                                        for (int i=0 ; i<num_partns1+1; ++i) {
                                            auto bound_partn1 = o1_global_lb + increment1*i;
                                            bound_partn1.eval_all();
                                            for (size_t inst = 0; inst< nb_ins; inst++){
                                                auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                                auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                                string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"}";
                                                bounds1.set_val(cur_idx,bound_partn1.eval(inst));
                                                for(int j=0; j<num_partns2+1; ++j){
                                                    auto bound_partn2 = o2_global_lb + increment2*j;
                                                    bound_partn2.eval_all();
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"}";
                                                    EP.set_val(cur_idx,(bound_partn1.eval(inst)*bound_partn2.eval(inst)));
                                                }
                                            }
                                        }
                                        // fill bounds2
                                        for (int i=0 ; i<num_partns2+1; ++i) {
                                            auto bound_partn2 = o2_global_lb + increment2*i;
                                            bound_partn2.eval_all();
                                            for (size_t inst = 0; inst< nb_ins; inst++){
                                                auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                                auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                                string cur_idx = cur_var_idx+","+name2+"{"+to_string(i+1)+"}";
                                                bounds2.set_val(cur_idx,bound_partn2.eval(inst));
                                            }
                                        }
                                        
                                        // Lambda coefficient matrix when linking with partition variables
                                        param<> lambda_coef1(name1+"_lambda_linking_coefficients1");
                                        param<> lambda_coef2(name2+"_lambda_linking_coefficients2");
                                        
                                        // Partition coefficient matrix when linking with lambda variables
                                        param<> on_coef1(name1+"_partition_linking_coefficients1");
                                        param<> on_coef2(name2+"_partition_linking_coefficients2");
                                        
                                        // create constraint indices
                                        indices const_idx1("const_idx1");
                                        indices const_idx2("const_idx2");
                                        
                                        if(model_type == "lambda_II"){
                                            
                                            //fill constraint indices
                                            for (int i=0; i<num_partns1+1; ++i){
                                                const_idx1.add(to_string(i+1));
                                            }
                                            
                                            //fill constraint indices
                                            for (int i=0; i<num_partns2+1; ++i){
                                                const_idx2.add(to_string(i+1));
                                            }
                                            
                                            // Lambda coefficient matrix when linking with partition variables
                                            if(num_partns1 > 1) lambda_coef1.in(indices(inst_partition_lambda, const_idx1));
                                            if(num_partns2 > 1) lambda_coef2.in(indices(inst_partition_lambda, const_idx2));
                                            
                                            // Lambda coefficient matrix when linking with partition variables
                                            if(num_partns1 > 1) on_coef1.in(indices(unique_ids, partns1, const_idx1));
                                            if(num_partns2 > 1) on_coef2.in(indices(unique_ids, partns2, const_idx2));
                                            
                                            // fill lambda_coef1 and lambda_coef2
                                            for (size_t inst = 0; inst< nb_ins; inst++){
                                                auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                                auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                                for (int i=0 ; i<num_partns1+1; ++i) {
                                                    for (int j=0 ; j<num_partns2+1; ++j) {
                                                        string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                        if(num_partns1 > 1) lambda_coef1.set_val(cur_idx,1);
                                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(j+1);
                                                        if(num_partns2 > 1) lambda_coef2.set_val(cur_idx,1);
                                                    }
                                                }
                                            }
                                            
                                            // fill on_coef1 and on_coef2
                                            for (size_t inst = 0; inst< nb_ins; inst++){
                                                auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                                auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                                string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+to_string(1);
                                                if(num_partns1 > 1) on_coef1.set_val(cur_idx,1);
                                                cur_idx = cur_var_idx+","+name2+"{"+to_string(1)+"},"+to_string(1);
                                                if(num_partns2 > 1) on_coef2.set_val(cur_idx,1);
                                                if(num_partns1 > 1) {
                                                    for (int i=1 ; i<num_partns1; ++i) {
                                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(i)+"},"+to_string(i+1);
                                                        on_coef1.set_val(cur_idx,1);
                                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+to_string(i+1);
                                                        on_coef1.set_val(cur_idx,1);
                                                    }
                                                }
                                                if(num_partns2 > 1) {
                                                    for (int i=1 ; i<num_partns2; ++i) {
                                                        cur_idx = cur_var_idx+","+name2+"{"+to_string(i)+"},"+to_string(i+1);
                                                        on_coef2.set_val(cur_idx,1);
                                                        cur_idx = cur_var_idx+","+name2+"{"+to_string(i+1)+"},"+to_string(i+1);
                                                        on_coef2.set_val(cur_idx,1);
                                                    }
                                                }
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(num_partns1)+"},"+to_string(num_partns1+1);
                                                if(num_partns1 > 1) on_coef1.set_val(cur_idx,1);
                                                cur_idx = cur_var_idx+","+name2+"{"+to_string(num_partns2)+"},"+to_string(num_partns2+1);
                                                if(num_partns2 > 1) on_coef2.set_val(cur_idx,1);
                                            }
                                        }
                                        
                                        
                                        else /*means model_type == "lambda_III" */{
                                            
                                            //fill constraint indices
                                            for (int i=0; i<(num_partns1-2)*2+2; ++i){
                                                const_idx1.add(to_string(i+1));
                                            }
                                            
                                            //fill constraint indices
                                            for (int i=0; i<(num_partns2-2)*2+2; ++i){
                                                const_idx2.add(to_string(i+1));
                                            }
                                            
                                            // Lambda coefficient matrix when linking with partition variables
                                            if(num_partns1 > 1) lambda_coef1.in(indices(inst_partition_lambda, const_idx1));
                                            if(num_partns2 > 1) lambda_coef2.in(indices(inst_partition_lambda, const_idx2));
                                            
                                            // Partition coefficient matrix when linking with lambda variables
                                            if(num_partns1 > 1) on_coef1.in(indices(unique_ids, partns1, const_idx1));
                                            if(num_partns2 > 1) on_coef2.in(indices(unique_ids, partns2, const_idx2));
                                            
                                            // fill lambda_coef1 and lambda_coef2
                                            for (size_t inst = 0; inst< nb_ins; inst++){
                                                auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                                auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                                if(num_partns1 > 1) {
                                                    for (int j=0; j<num_partns2+1; ++j) {
                                                        string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(1);
                                                        lambda_coef1.set_val(cur_idx,1);
                                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(num_partns1+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string((num_partns1-2)*2+2);
                                                        lambda_coef1.set_val(cur_idx,1);
                                                    }
                                                    
                                                    for (int i=1 ; i<(num_partns1-2)*2+1; i=i+2) {
                                                        for (int j=(i-1)/2 + 2; j<num_partns1+1; ++j) {
                                                            for(int k=0; k<num_partns2+1; ++k){
                                                                string cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+name2+"{"+to_string(k+1)+"},"+to_string(i+1);
                                                                lambda_coef1.set_val(cur_idx,1);
                                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+name2+"{"+to_string(k+1)+"},"+to_string(i+2);
                                                                lambda_coef1.set_val(cur_idx,-1);
                                                            }
                                                        }
                                                    }
                                                }
                                                if(num_partns2 > 1) {
                                                    for (int i=0; i<num_partns1+1; ++i) {
                                                        string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(1)+"},"+to_string(1);
                                                        lambda_coef2.set_val(cur_idx,1);
                                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(num_partns2+1)+"},"+to_string((num_partns2-2)*2+2);
                                                        lambda_coef2.set_val(cur_idx,1);
                                                    }
                                                    for (int i=1 ; i<(num_partns2-2)*2+1; i=i+2) {
                                                        for (int j=(i-1)/2 + 2; j<num_partns2+1; ++j) {
                                                            for(int k=0; k<num_partns1+1; ++k){
                                                                string cur_idx = cur_var_idx+","+name1+"{"+to_string(k+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                                lambda_coef2.set_val(cur_idx,1);
                                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(k+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                                lambda_coef2.set_val(cur_idx,-1);
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                            
                                            
                                            // fill on_coef1 and on_coef2
                                            for (size_t inst = 0; inst< nb_ins; inst++){
                                                auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                                auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                                string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+to_string(1);
                                                if(num_partns1 > 1) {
                                                    on_coef1.set_val(cur_idx,1);
                                                    
                                                    for (int i=1; i<num_partns1; ++i) {
                                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+to_string(2);
                                                        on_coef1.set_val(cur_idx, 1);
                                                    }
                                                    
                                                    for (int i=2 ; i<(num_partns1-2)*2+2; i=i+2) {
                                                        for (int j=i/2+1; j<num_partns1; ++j) {
                                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                            on_coef1.set_val(cur_idx,-1);
                                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                            on_coef1.set_val(cur_idx,1);
                                                        }
                                                    }
                                                }
                                                if(num_partns2 > 1) {
                                                    cur_idx = cur_var_idx+","+name2+"{"+to_string(1)+"},"+to_string(1);
                                                    on_coef2.set_val(cur_idx,1);
                                                    for (int i=1; i<num_partns2; ++i) {
                                                        cur_idx = cur_var_idx+","+name2+"{"+to_string(i+1)+"},"+to_string(2);
                                                        on_coef2.set_val(cur_idx, 1);
                                                    }
                                                    
                                                    for (int i=2 ; i<(num_partns2-2)*2+2; i=i+2) {
                                                        for (int j=i/2+1; j<num_partns2; ++j) {
                                                            cur_idx = cur_var_idx+","+name2+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                            on_coef2.set_val(cur_idx,-1);
                                                            cur_idx = cur_var_idx+","+name2+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                            on_coef2.set_val(cur_idx,1);
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                        
                                        
                                        /** Constraints */
                                        // Representation of the bilinear term with convex combination
                                        Constraint<> bln_rep(pair.first+"_bln_rep");
                                        bln_rep = EP.in_matrix(nb_entries,total_entries-nb_entries)*lambda.in_matrix(nb_entries,total_entries-nb_entries) - vlift.in(unique_ids);
                                        add(bln_rep.in(unique_ids) == 0);
                                        
                                        // Representation of o1 with convex combination
                                        Constraint<> o1_rep(pair.first+"_o1_rep");
                                        o1_rep = bounds1.from_ith(0,inst_partition_lambda).in_matrix(nb_entries, 1) * lambda.in_matrix(nb_entries,total_entries-nb_entries) - o1.in(o1_ids);
                                        add(o1_rep.in(unique_ids) == 0);
                                        
                                        // Representation of o2 with convex combination
                                        Constraint<> o2_rep(pair.first+"_o2_rep");
                                        o2_rep = bounds2.in_ignore_ith(nb_entries, 1, inst_partition_lambda).in_matrix(nb_entries,1) * lambda.in_matrix(nb_entries,total_entries-nb_entries) - o2.in(o2_ids);
                                        add(o2_rep.in(unique_ids) == 0);
                                        
                                        // Linking partition variables1 with lambda
                                        if(model_type == "lambda_II"){
                                            if(num_partns1 > 1) {
                                                Constraint<> on_link_lambda1(pair.first+"_on_link_lambda1_II");
                                                on_link_lambda1 = lambda.in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef1.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef1.in_matrix(nb_entries,total_entries-nb_entries) - on1.in_ignore_ith(nb_entries_v1,nb_entries_v2,on_coef1.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef1.in_matrix(nb_entries,1);
                                                add(on_link_lambda1.in(indices(unique_ids,const_idx1)) <= 0);
                                            }
                                            if(num_partns2 > 1) {
                                                Constraint<> on_link_lambda2(pair.first+"_on_link_lambda2_II");
                                                on_link_lambda2 = lambda.in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef2.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef2.in_matrix(nb_entries,total_entries-nb_entries) - on2.in_ignore_ith(0,nb_entries_v1,on_coef2.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef2.in_matrix(nb_entries,1);
                                                add(on_link_lambda2.in(indices(unique_ids,const_idx2)) <= 0);
                                            }
                                        }
                                        else{
                                            if(num_partns1 > 1) {
                                                Constraint<> on_link_lambda1(pair.first+"_on_link_lambda1_III");
                                                on_link_lambda1 = lambda.in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef1.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef1.in_matrix(nb_entries,total_entries-nb_entries) - on1.in_ignore_ith(nb_entries_v1,nb_entries_v2,on_coef1.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef1.in_matrix(nb_entries,1);
                                                add(on_link_lambda1.in(indices(unique_ids,const_idx1)) <= 0);
                                            }
                                            if(num_partns2 > 1) {
                                                Constraint<> on_link_lambda2(pair.first+"_on_link_lambda2_III");
                                                on_link_lambda2 = lambda.in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef2.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef2.in_matrix(nb_entries,total_entries-nb_entries) - on2.in_ignore_ith(0,nb_entries_v1,on_coef2.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef2.in_matrix(nb_entries,1);
                                                add(on_link_lambda2.in(indices(unique_ids,const_idx2)) <= 0);
                                            }
                                        }
                                        // sum over lambda
                                        Constraint<> lambdaSum(pair.first+"_lambdaSum");
                                        lambdaSum = sum(lambda.in_matrix(nb_entries,total_entries-nb_entries));
                                        add(lambdaSum.in(unique_ids) == 1);
                                    }
                                }
                            }
                        }
                    }
                    else {
                        add_McCormick(pair.first, vlift.in(unique_ids), o1.in(o1_ids), o2.in(o2_ids));
                    }
                }
                else {
                    auto vlift = static_pointer_cast<var<type>>(it->second);
                    auto added = vlift->add_bounds(lb,ub);
                    lt._p = make_shared<var<type>>(vlift->in(ids));
                    if(!added.empty()){
                        assert(o1._indices->size()==o2._indices->size());
                        if(added.size()!=o1._indices->size()){/* If some keys are repeated, remove them from the refs of o1 and o2 */
                            auto keep_refs = ids.diff_refs(added);
                            o1_ids.filter_refs(keep_refs);
                            o2_ids.filter_refs(keep_refs);
                        }
                        reindex_vars();
                        // If some keys are repeated in individual indices, remove them from the refs of o1 and o2
                        auto o1_ids_uq = o1_ids;
                        auto o2_ids_uq = o2_ids;
                        auto keep_refs1 = o1_ids_uq.get_unique_refs();
                        auto keep_refs2 = o2_ids_uq.get_unique_refs();
                        o1_ids_uq.filter_refs(keep_refs1);
                        o2_ids_uq.filter_refs(keep_refs2);
                        reindex_vars();
                        
                        //check the sign of the lift and the correspoinding boudning functions
                        if(c.check_soc() && c.is_eq()){
                        if(lift_sign){
                            vlift->_lift_ub = true;
                            vlift->_lift_lb = false;
                        }
                        else{
                            vlift->_lift_ub = false;
                            vlift->_lift_lb = true;
                        }
                        }
                        else{
                            vlift->_lift_ub = true;
                            vlift->_lift_lb = true;
                        }
                        if((num_partns1 > 1) || (num_partns2 > 1)) {
                            if (o1 == o2) //if the variables are same add 1d partition
                            {
                                DebugOn("<<<<<<<<<< THIS IS SEEN BOTH -> SINGLE <<<<<<<<<<<" << endl);
                                auto binvar_ptr1 = _vars_name.find(name1+"_binary");
                                auto binvar1 = static_pointer_cast<var<int>>(binvar_ptr1->second);
                                
                                indices partns("partns");
                                for (int i = 0; i < num_partns1 ; ++i)
                                {
                                    partns.add(name1+  "{" + to_string(i+1) + "}");
                                }
                                auto inst_partition = indices(added,partns);
                                
                                param<int> lb1("lb1"), ub1("ub1");
                                lb1.in(added,partns);
                                ub1.in(added,partns);
                                lb1.set_val(0), ub1.set_val(1);
                                auto added1 = binvar1->add_bounds(lb1,ub1);
                                reindex_vars();
                                
                                auto nb_entries_v1 = o1_ids.get_nb_entries();
                                auto nb_entries = added.get_nb_entries();
                                auto total_entries = inst_partition.get_nb_entries();
                                
                                Constraint<> onSumComb(pair.first+"_binarySum");
                                onSumComb = sum((binvar1->in(added1)).in_matrix(nb_entries,total_entries-nb_entries));
                                add(onSumComb.in(added) == 1);
                                
                                if(model_type == "on/off"){//if on/off is chosen
                                    add_on_off_McCormick_refined(pair.first, vlift->in(added), o1.in(o1_ids), o2.in(o2_ids), binvar1->in(added1));
                                }
                                
                                else{ //means it is one of the lambda formulations
                                    
                                    //difference is this has one more partition index
                                    indices partns_lambda("partns_lambda");
                                    for (int i = 0; i < num_partns1+1 ; ++i)
                                    {
                                        partns_lambda.add(name1+ "{" +to_string(i+1) + "}");
                                    }
                                    auto inst_partition_lambda = indices(added,partns_lambda);
                                    
                                    // Convex combination variables
                                    auto lambda_ptr = _vars_name.find(name1+"_lambda");
                                    auto lambda = static_pointer_cast<var<double>>(lambda_ptr->second);
                                    param<double> lb_lambda("lb_lambda"), ub_lambda("ub_lambda");
                                    lb_lambda.in(added,partns_lambda);
                                    ub_lambda.in(added,partns_lambda);
                                    lb_lambda.set_val(0), ub_lambda.set_val(1);
                                    auto added_lambda = lambda->add_bounds(lb_lambda,ub_lambda);
                                    reindex_vars();
                                    
                                    /** Parameters */
                                    // Bounds on variable v1 & v2
                                    param<> bounds(name1+"_bounds");
                                    bounds.in(inst_partition_lambda);
                                    
                                    // Function values on the extreme points
                                    param<> EP(name1+name2+"_grid_values");
                                    EP.in(inst_partition_lambda);
                                    
                                    size_t nb_ins = vlift->in(added).get_nb_inst();
                                    auto o1_global_lb = o1.get_lb();
                                    auto increment = (o1.get_ub() - o1_global_lb)/num_partns1;
                                    
                                    // fill bounds and function values
                                    for (int i=0 ; i<num_partns1+1; ++i) {
                                        auto bound_partn = o1_global_lb + increment*i;
                                        bound_partn.eval_all();
                                        for (size_t inst = 0; inst< nb_ins; inst++){
                                            auto cur_var_id = vlift->get_id_inst(inst);
                                            auto cur_var_idx = added._keys->at(cur_var_id);
                                            string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"}";
                                            bounds.set_val(cur_idx,bound_partn.eval(inst));
                                            EP.set_val(cur_idx,(bound_partn.eval(inst)*bound_partn.eval(inst)));
                                        }
                                    }
                                    
                                    // Lambda coefficient matrix when linking with partition variables
                                    param<> lambda_coef(name1+"_lambda_linking_coefficients");
                                    
                                    // Partition coefficient matrix when linking with lambda variables
                                    param<> on_coef(name1+"_partition_linking_coefficients");
                                    
                                    // create constraint indices
                                    indices const_idx("const_idx");
                                    
                                    if(model_type == "lambda_II"){
                                        
                                        //fill constraint indices
                                        for (int i=0; i<num_partns1+1; ++i){
                                            const_idx.add(to_string(i+1));
                                        }
                                        
                                        // Lambda coefficient matrix when linking with partition variables
                                        lambda_coef.in(indices(inst_partition_lambda, const_idx));
                                        
                                        // Partition coefficient matrix when linking with lambda variables
                                        on_coef.in(indices(inst_partition, const_idx));
                                        
                                        // fill lambda_coef
                                        for (size_t inst = 0; inst< nb_ins; inst++){
                                            auto cur_var_id = vlift->get_id_inst(inst);
                                            auto cur_var_idx = added._keys->at(cur_var_id);
                                            for (int i=0 ; i<num_partns1+1; ++i) {
                                                string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+to_string(i+1);
                                                lambda_coef.set_val(cur_idx,1);
                                            }
                                        }
                                        
                                        // fill on_coef
                                        for (size_t inst = 0; inst< nb_ins; inst++){
                                            auto cur_var_id = vlift->get_id_inst(inst);
                                            auto cur_var_idx = added._keys->at(cur_var_id);
                                            string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+to_string(1);
                                            on_coef.set_val(cur_idx,1);
                                            for (int i=1 ; i<num_partns1; ++i) {
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(i)+"},"+to_string(i+1);
                                                on_coef.set_val(cur_idx,1);
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+to_string(i+1);
                                                on_coef.set_val(cur_idx,1);
                                            }
                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(num_partns1)+"},"+to_string(num_partns1+1);
                                            on_coef.set_val(cur_idx,1);
                                        }
                                    }
                                    
                                    else /*means model_type == "lambda_III" */{
                                        
                                        //fill constraint indices
                                        for (int i=0; i<(num_partns1-2)*2+2; ++i){
                                            const_idx.add(to_string(i+1));
                                        }
                                        
                                        // Lambda coefficient matrix when linking with partition variables
                                        lambda_coef.in(indices(inst_partition_lambda, const_idx));
                                        
                                        // Partition coefficient matrix when linking with lambda variables
                                        on_coef.in(indices(inst_partition, const_idx));
                                        
                                        // fill lambda_coef
                                        for (size_t inst = 0; inst< nb_ins; inst++){
                                            auto cur_var_id = vlift->get_id_inst(inst);
                                            auto cur_var_idx = added._keys->at(cur_var_id);
                                            string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+to_string(1);
                                            lambda_coef.set_val(cur_idx,1);
                                            for (int i=1 ; i<(num_partns1-2)*2+1; i=i+2) {
                                                for (int j=(i-1)/2 + 2; j<num_partns1+1; ++j) {
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                    lambda_coef.set_val(cur_idx,1);
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                    lambda_coef.set_val(cur_idx,-1);
                                                }
                                            }
                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(num_partns1+1)+"},"+to_string((num_partns1-2)*2+2);
                                            lambda_coef.set_val(cur_idx,1);
                                        }
                                        
                                        // fill on_coef
                                        for (size_t inst = 0; inst< nb_ins; inst++){
                                            auto cur_var_id = vlift->get_id_inst(inst);
                                            auto cur_var_idx = added._keys->at(cur_var_id);
                                            string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+to_string(1);
                                            on_coef.set_val(cur_idx,1);
                                            
                                            for (int i=1; i<num_partns1; ++i) {
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+to_string(2);
                                                on_coef.set_val(cur_idx, 1);
                                            }
                                            
                                            for (int i=2 ; i<(num_partns1-2)*2+2; i=i+2) {
                                                for (int j=i/2+1; j<num_partns1; ++j) {
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                    on_coef.set_val(cur_idx,-1);
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                    on_coef.set_val(cur_idx,1);
                                                }
                                            }
                                        }
                                        
                                    }
                                    
                                    
                                    /** Constraints */
                                    
                                    if (vlift->_lift_ub){
                                    // Representation of the quadratic term with secant
                                    Constraint<> quad_ub(pair.first+"_quad_ub");
                                    quad_ub = EP.in_matrix(nb_entries,total_entries-nb_entries)*lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries) - vlift->in(added);
                                    add(quad_ub.in(added) >= 0); /*using it as the upper bound to be valid*/
                                    }
                                    
                                    if (vlift->_lift_lb){
                                    Constraint<> quad_lb(pair.first+"_quad_lb");
                                    quad_lb = o1.from_ith(0,added)*o2.from_ith(nb_entries_v1,added) - vlift->in(added);
                                    quad_lb._relaxed = true;
                                    add(quad_lb.in(added) <= 0); /*using it as the lower bound to be valid*/
                                    }
                                    
                                    // Representation of o1 with convex combination
                                    Constraint<> o1_rep(pair.first+"_o1_rep");
                                    o1_rep = bounds.from_ith(0,inst_partition_lambda).in_matrix(nb_entries, 1) * lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries) - o1.in(o1_ids);
                                    add(o1_rep.in(added) == 0);
                                    
                                    // Linking partition variables with lambda
                                    if(model_type == "lambda_II"){
                                        Constraint<> on_link_lambda(pair.first+"_on_link_lambda_II");
                                        on_link_lambda = lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef.in_matrix(nb_entries,total_entries-nb_entries) - binvar1->in(added1).in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,on_coef.get_matrix_ids(nb_entries,total_entries-nb_entries)) * on_coef.in_matrix(nb_entries,total_entries-nb_entries);
                                        add(on_link_lambda.in(indices(added,const_idx)) <= 0);
                                    }
                                    else{
                                        Constraint<> on_link_lambda(pair.first+"_on_link_lambda_III");
                                        on_link_lambda = lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef.in_matrix(nb_entries,total_entries-nb_entries) - binvar1->in(added1).in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,on_coef.get_matrix_ids(nb_entries,total_entries-nb_entries)) * on_coef.in_matrix(nb_entries,total_entries-nb_entries);
                                        add(on_link_lambda.in(indices(added,const_idx)) <= 0);
                                    }
                                    
                                    
                                    // sum over lambda
                                    Constraint<> lambdaSum(pair.first+"_lambdaSum");
                                    lambdaSum = sum(lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries));
                                    add(lambdaSum.in(added) == 1);
                                }
                            }
                            else{ //else add 2d partition
                                
                                if(name1 == name2){
                                    DebugOn("<<<<<<<<<< THIS IS SEEN BOTH -> DOUBLE -> SAME VARS <<<<<<<<<<<" << endl);
                                    
                                    indices partns1("partns1");
                                    for (int i = 0; i < num_partns1 ; ++i)
                                    {
                                        partns1.add(name1+ "{" + to_string(i+1) + "}");
                                    }
                                    
                                    auto binvar_ptr1 = _vars_name.find(name1+"_binary");
                                    auto binvar1 = static_pointer_cast<var<int>>(binvar_ptr1->second);
                                    param<int> lb1("lb1"), ub1("ub1");
                                    lb1.in(union_ids(o1_ids_uq,o2_ids_uq),partns1);
                                    ub1.in(union_ids(o1_ids_uq,o2_ids_uq),partns1);
                                    lb1.set_val(0), ub1.set_val(1);
                                    auto added1 = binvar1->add_bounds(lb1,ub1);
                                    reindex_vars();
                                    
                                    indices partns2("partns2");
                                    for (int i = 0; i < num_partns2 ; ++i)
                                    {
                                        partns2.add(name2+ "{" + to_string(i+1) + "}");
                                    }
                                    
                                    auto nb_entries_v1 = o1_ids.get_nb_entries();
                                    auto nb_entries_v2 = o2_ids.get_nb_entries();
                                    auto nb_entries = added.get_nb_entries();
                                    
                                    if(!added1.empty()){
                                        Constraint<> onSum1(o1._name+"_binarySum");
                                        onSum1 = sum(binvar1->in(added1).in_matrix(nb_entries_v1,1));
                                        auto vset1 = added1.from_ith(0,nb_entries_v1);
                                        vset1.filter_refs(vset1.get_unique_refs());
                                        add(onSum1.in(vset1) == 1);
                                    }
                                    
                                    if(model_type == "on/off"){//if on/off is chosen
                                        
                                        indices partns("partns");
                                        partns = indices(partns1,partns2);
                                        auto inst_partition = indices(added,partns);
                                        auto total_entries = inst_partition.get_nb_entries();
                                        
                                        auto binvar_ptr3 = _vars_name.find(name1+name2+"_binary");
                                        auto binvar3 = static_pointer_cast<var<int>>(binvar_ptr3->second);
                                        param<int> lb3("lb3"), ub3("ub3");
                                        lb3.in(added,partns);
                                        ub3.in(added,partns);
                                        lb3.set_val(0), ub3.set_val(1);
                                        auto added3 = binvar3->add_bounds(lb3,ub3);
                                        reindex_vars();
                                        
                                        Constraint<> onLink1(pair.first+"_binaryLink1");
                                        onLink1 = binvar1->from_ith(0,inst_partition.ignore_ith(nb_entries_v1, nb_entries_v1)) - binvar3->in(inst_partition);
                                        add(onLink1.in(inst_partition) >= 0);
                                        
                                        Constraint<> onLink2(pair.first+"_binaryLink2");
                                        onLink2 = binvar1->in_ignore_ith(nb_entries_v1,1,inst_partition.ignore_ith(0,nb_entries_v1)) - binvar3->in(inst_partition);
                                        add(onLink2.in(inst_partition) >= 0);
                                        
                                        Constraint<> onLink3(pair.first+"_binaryLink3");
                                        onLink3 = binvar1->from_ith(0,inst_partition.ignore_ith(nb_entries_v1, nb_entries_v1)) + binvar1->in_ignore_ith(nb_entries_v1,1,inst_partition.ignore_ith(0,nb_entries_v1)) - 1 - binvar3->in(inst_partition);
                                        add(onLink3.in(inst_partition) <= 0);
                                        
                                        Constraint<> onSumComb(pair.first+"_binarySum");
                                        onSumComb = sum((binvar3->in(added3)).in_matrix(nb_entries,total_entries-nb_entries));
                                        add(onSumComb.in(added) == 1);
                                        
                                        add_on_off_McCormick_refined(pair.first, vlift->in(added), o1.in(o1_ids), o2.in(o2_ids), binvar3->in(added3));
                                    }
                                    
                                    else{//means it is one of the lambda formulations
                                        
                                        //difference is this has one more partition index
                                        indices partns1_lambda("partns1_lambda");
                                        for (int i = 0; i < num_partns1+1; ++i)
                                        {
                                            partns1_lambda.add(name1+ "{" +to_string(i+1) + "}");
                                        }
                                        
                                        indices partns2_lambda("partns2_lambda");
                                        for (int i = 0; i < num_partns2+1; ++i)
                                        {
                                            partns2_lambda.add(name2+ "{" +to_string(i+1) + "}");
                                        }
                                        
                                        indices partns_lambda("partns_lambda");
                                        partns_lambda = indices(partns1_lambda,partns2_lambda);
                                        auto inst_partition_lambda = indices(added,partns_lambda);
                                        auto inst_partition_bounds1 = indices(added,partns1_lambda);
                                        auto inst_partition_bounds2 = indices(added,partns2_lambda);
                                        
                                        // Convex combination variables
                                        auto lambda_ptr = _vars_name.find(name1+name2+"_lambda");
                                        auto lambda = static_pointer_cast<var<double>>(lambda_ptr->second);
                                        param<double> lb_lambda("lb_lambda"), ub_lambda("ub_lambda");
                                        lb_lambda.in(added,partns_lambda);
                                        ub_lambda.in(added,partns_lambda);
                                        lb_lambda.set_val(0), ub_lambda.set_val(1);
                                        auto added_lambda = lambda->add_bounds(lb_lambda,ub_lambda);
                                        reindex_vars();
                                        
                                        /** Parameters */
                                        // Bounds on variable v1 & v2
                                        param<> bounds1(name1+"_bounds1");
                                        bounds1.in(inst_partition_bounds1);
                                        
                                        param<> bounds2(name2+"_bounds2");
                                        bounds2.in(inst_partition_bounds2);
                                        
                                        // Function values on the extreme points
                                        param<> EP(name1+name2+"_grid_values");
                                        EP.in(inst_partition_lambda);
                                        auto total_entries = inst_partition_lambda.get_nb_entries();
                                        
                                        size_t nb_ins = vlift->in(added).get_nb_inst();
                                        auto o1_global_lb = o1.get_lb();
                                        auto increment1 = (o1.get_ub() - o1_global_lb)/num_partns1;
                                        
                                        auto o2_global_lb = o2.get_lb();
                                        auto increment2 = (o2.get_ub() - o2_global_lb)/num_partns2;
                                        
                                        // fill bounds and function values
                                        for (int i=0 ; i<num_partns1+1; ++i) {
                                            auto bound_partn1 = o1_global_lb + increment1*i;
                                            bound_partn1.eval_all();
                                            auto bound_partn2 = o2_global_lb + increment2*i;
                                            bound_partn2.eval_all();
                                            for (size_t inst = 0; inst< nb_ins; inst++){
                                                auto cur_var_id = vlift->get_id_inst(inst);
                                                auto cur_var_idx = added._keys->at(cur_var_id);
                                                string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"}";
                                                bounds1.set_val(cur_idx,bound_partn1.eval(inst));
                                                bounds2.set_val(cur_idx,bound_partn2.eval(inst));
                                                for(int j=0; j<num_partns2+1; ++j){
                                                    auto bound_partn2_temp = o2_global_lb + increment2*j;
                                                    bound_partn2_temp.eval_all();
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"}";
                                                    EP.set_val(cur_idx,(bound_partn1.eval(inst)*bound_partn2_temp.eval(inst)));
                                                }
                                            }
                                        }
                                        
                                        // Lambda coefficient matrix when linking with partition variables
                                        param<> lambda_coef1(name1+"_lambda_linking_coefficients1");
                                        param<> lambda_coef2(name2+"_lambda_linking_coefficients2");
                                        
                                        // Partition coefficient matrix when linking with lambda variables
                                        param<> on_coef1(name1+"_partition_linking_coefficients1");
                                        param<> on_coef2(name2+"_partition_linking_coefficients2");
                                        
                                        // create constraint indices
                                        indices const_idx1("const_idx1");
                                        indices const_idx2("const_idx2");
                                        
                                        if(model_type == "lambda_II"){
                                            
                                            //fill constraint indices
                                            for (int i=0; i<num_partns1+1; ++i){
                                                const_idx1.add(to_string(i+1));
                                            }
                                            
                                            //fill constraint indices
                                            for (int i=0; i<num_partns2+1; ++i){
                                                const_idx2.add(to_string(i+1));
                                            }
                                            
                                            // Lambda coefficient matrix when linking with partition variables
                                            lambda_coef1.in(indices(inst_partition_lambda, const_idx1));
                                            lambda_coef2.in(indices(inst_partition_lambda, const_idx2));
                                            
                                            // Lambda coefficient matrix when linking with partition variables
                                            on_coef1.in(indices(added, partns1, const_idx1));
                                            on_coef2.in(indices(added, partns2, const_idx2));
                                            
                                            // fill lambda_coef1 and lambda_coef2
                                            for (size_t inst = 0; inst< nb_ins; inst++){
                                                auto cur_var_id = vlift->get_id_inst(inst);
                                                auto cur_var_idx = added._keys->at(cur_var_id);
                                                for (int i=0 ; i<num_partns1+1; ++i) {
                                                    for (int j=0 ; j<num_partns2+1; ++j) {
                                                        string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                        lambda_coef1.set_val(cur_idx,1);
                                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(j+1);
                                                        lambda_coef2.set_val(cur_idx,1);
                                                    }
                                                }
                                            }
                                            
                                            // fill on_coef1 and on_coef2
                                            for (size_t inst = 0; inst< nb_ins; inst++){
                                                auto cur_var_id = vlift->get_id_inst(inst);
                                                auto cur_var_idx = added._keys->at(cur_var_id);
                                                string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+to_string(1);
                                                on_coef1.set_val(cur_idx,1);
                                                cur_idx = cur_var_idx+","+name2+"{"+to_string(1)+"},"+to_string(1);
                                                on_coef2.set_val(cur_idx,1);
                                                for (int i=1 ; i<num_partns1; ++i) {
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(i)+"},"+to_string(i+1);
                                                    on_coef1.set_val(cur_idx,1);
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+to_string(i+1);
                                                    on_coef1.set_val(cur_idx,1);
                                                }
                                                for (int i=1 ; i<num_partns2; ++i) {
                                                    cur_idx = cur_var_idx+","+name2+"{"+to_string(i)+"},"+to_string(i+1);
                                                    on_coef2.set_val(cur_idx,1);
                                                    cur_idx = cur_var_idx+","+name2+"{"+to_string(i+1)+"},"+to_string(i+1);
                                                    on_coef2.set_val(cur_idx,1);
                                                }
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(num_partns1)+"},"+to_string(num_partns1+1);
                                                on_coef1.set_val(cur_idx,1);
                                                cur_idx = cur_var_idx+","+name2+"{"+to_string(num_partns2)+"},"+to_string(num_partns2+1);
                                                on_coef2.set_val(cur_idx,1);
                                            }
                                        }
                                        
                                        
                                        else /*means model_type == "lambda_III" */{
                                            
                                            //fill constraint indices
                                            for (int i=0; i<(num_partns1-2)*2+2; ++i){
                                                const_idx1.add(to_string(i+1));
                                            }
                                            
                                            //fill constraint indices
                                            for (int i=0; i<(num_partns2-2)*2+2; ++i){
                                                const_idx2.add(to_string(i+1));
                                            }
                                            
                                            // Lambda coefficient matrix when linking with partition variables
                                            lambda_coef1.in(indices(inst_partition_lambda, const_idx1));
                                            lambda_coef2.in(indices(inst_partition_lambda, const_idx2));
                                            
                                            // Partition coefficient matrix when linking with lambda variables
                                            on_coef1.in(indices(added, partns1, const_idx1));
                                            on_coef2.in(indices(added, partns2, const_idx2));
                                            
                                            // fill lambda_coef1 and lambda_coef2
                                            for (size_t inst = 0; inst< nb_ins; inst++){
                                                auto cur_var_id = vlift->get_id_inst(inst);
                                                auto cur_var_idx = added._keys->at(cur_var_id);
                                                for (int j=0; j<num_partns2+1; ++j) {
                                                    string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(1);
                                                    lambda_coef1.set_val(cur_idx,1);
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(num_partns1+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string((num_partns1-2)*2+2);
                                                    lambda_coef1.set_val(cur_idx,1);
                                                }
                                                
                                                for (int i=1 ; i<(num_partns1-2)*2+1; i=i+2) {
                                                    for (int j=(i-1)/2 + 2; j<num_partns1+1; ++j) {
                                                        for(int k=0; k<num_partns2+1; ++k){
                                                            string cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+name2+"{"+to_string(k+1)+"},"+to_string(i+1);
                                                            lambda_coef1.set_val(cur_idx,1);
                                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+name2+"{"+to_string(k+1)+"},"+to_string(i+2);
                                                            lambda_coef1.set_val(cur_idx,-1);
                                                        }
                                                    }
                                                }
                                                
                                                for (int i=0; i<num_partns1+1; ++i) {
                                                    string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(1)+"},"+to_string(1);
                                                    lambda_coef2.set_val(cur_idx,1);
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(num_partns2+1)+"},"+to_string((num_partns2-2)*2+2);
                                                    lambda_coef2.set_val(cur_idx,1);
                                                }
                                                
                                                for (int i=1 ; i<(num_partns2-2)*2+1; i=i+2) {
                                                    for (int j=(i-1)/2 + 2; j<num_partns2+1; ++j) {
                                                        for(int k=0; k<num_partns1+1; ++k){
                                                            string cur_idx = cur_var_idx+","+name1+"{"+to_string(k+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                            lambda_coef2.set_val(cur_idx,1);
                                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(k+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                            lambda_coef2.set_val(cur_idx,-1);
                                                        }
                                                    }
                                                }
                                            }
                                            
                                            
                                            
                                            // fill on_coef1 and on_coef2
                                            for (size_t inst = 0; inst< nb_ins; inst++){
                                                auto cur_var_id = vlift->get_id_inst(inst);
                                                auto cur_var_idx = added._keys->at(cur_var_id);
                                                string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+to_string(1);
                                                on_coef1.set_val(cur_idx,1);
                                                
                                                for (int i=1; i<num_partns1; ++i) {
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+to_string(2);
                                                    on_coef1.set_val(cur_idx, 1);
                                                }
                                                
                                                for (int i=2 ; i<(num_partns1-2)*2+2; i=i+2) {
                                                    for (int j=i/2+1; j<num_partns1; ++j) {
                                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                        on_coef1.set_val(cur_idx,-1);
                                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                        on_coef1.set_val(cur_idx,1);
                                                    }
                                                }
                                                
                                                cur_idx = cur_var_idx+","+name2+"{"+to_string(1)+"},"+to_string(1);
                                                on_coef2.set_val(cur_idx,1);
                                                for (int i=1; i<num_partns2; ++i) {
                                                    cur_idx = cur_var_idx+","+name2+"{"+to_string(i+1)+"},"+to_string(2);
                                                    on_coef2.set_val(cur_idx, 1);
                                                }
                                                
                                                for (int i=2 ; i<(num_partns2-2)*2+2; i=i+2) {
                                                    for (int j=i/2+1; j<num_partns2; ++j) {
                                                        cur_idx = cur_var_idx+","+name2+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                        on_coef2.set_val(cur_idx,-1);
                                                        cur_idx = cur_var_idx+","+name2+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                        on_coef2.set_val(cur_idx,1);
                                                    }
                                                }
                                            }
                                        }
                                        
                                        
                                        /** Constraints */
                                        // Representation of the bilinear term with convex combination
                                        Constraint<> bln_rep(pair.first+"_bln_rep");
                                        bln_rep = EP.in_matrix(nb_entries,total_entries-nb_entries)*lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries) - vlift->in(added);
                                        add(bln_rep.in(added) == 0);
                                        
                                        // Representation of o1 with convex combination
                                        Constraint<> o1_rep(pair.first+"_o1_rep");
                                        o1_rep = bounds1.from_ith(0,inst_partition_lambda).in_matrix(nb_entries, 1) * lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries) - o1.in(o1_ids);
                                        add(o1_rep.in(added) == 0);
                                        
                                        // Representation of o2 with convex combination
                                        Constraint<> o2_rep(pair.first+"_o2_rep");
                                        o2_rep = bounds2.in_ignore_ith(nb_entries, 1, inst_partition_lambda).in_matrix(nb_entries,1) * lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries) - o2.in(o2_ids);
                                        add(o2_rep.in(added) == 0);
                                        
                                        // Linking partition variables1 with lambda
                                        if(model_type == "lambda_II"){
                                            Constraint<> on_link_lambda1(pair.first+"_on_link_lambda1_II");
                                            on_link_lambda1 = lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef1.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef1.in_matrix(nb_entries,total_entries-nb_entries) - binvar1->in_ignore_ith(nb_entries_v1,nb_entries_v2,on_coef1.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef1.in_matrix(nb_entries,1);
                                            add(on_link_lambda1.in(indices(added,const_idx1)) <= 0);
                                            
                                            Constraint<> on_link_lambda2(pair.first+"_on_link_lambda2_II");
                                            on_link_lambda2 = lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef2.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef2.in_matrix(nb_entries,total_entries-nb_entries) - binvar1->in_ignore_ith(0,nb_entries_v1,on_coef2.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef2.in_matrix(nb_entries,1);
                                            add(on_link_lambda2.in(indices(added,const_idx2)) <= 0);
                                        }
                                        else{
                                            Constraint<> on_link_lambda1(pair.first+"_on_link_lambda1_III");
                                            on_link_lambda1 = lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef1.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef1.in_matrix(nb_entries,total_entries-nb_entries) - binvar1->in_ignore_ith(nb_entries_v1,nb_entries_v2,on_coef1.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef1.in_matrix(nb_entries,1);
                                            add(on_link_lambda1.in(indices(added,const_idx1)) <= 0);
                                            
                                            Constraint<> on_link_lambda2(pair.first+"_on_link_lambda2_III");
                                            on_link_lambda2 = lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef2.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef2.in_matrix(nb_entries,total_entries-nb_entries) - binvar1->in_ignore_ith(0,nb_entries_v1,on_coef2.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef2.in_matrix(nb_entries,1);
                                            add(on_link_lambda2.in(indices(added,const_idx2)) <= 0);
                                        }
                                        // sum over lambda
                                        Constraint<> lambdaSum(pair.first+"_lambdaSum");
                                        lambdaSum = sum(lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries));
                                        add(lambdaSum.in(added) == 1);
                                    }
                                    
                                }
                                else{
                                    DebugOn("<<<<<<<<<< THIS IS SEEN BOTH -> DOUBLE -> DIFF VARS <<<<<<<<<<<" << endl);
                                    
                                    auto binvar_ptr1 = _vars_name.find(name1+"_binary");
                                    auto binvar1 = static_pointer_cast<var<int>>(binvar_ptr1->second);
                                    indices partns1("partns1");
                                    for (int i = 0; i < num_partns1 ; ++i)
                                    {
                                        partns1.add(name1+ "{" +to_string(i+1) + "}");
                                    }
                                    param<int> lb1("lb1"), ub1("ub1");
                                    lb1.in(o1_ids_uq,partns1);
                                    ub1.in(o1_ids_uq,partns1);
                                    lb1.set_val(0), ub1.set_val(1);
                                    auto added1 = binvar1->add_bounds(lb1,ub1);
                                    reindex_vars();
                                    
                                    auto binvar_ptr2 = _vars_name.find(name2+"_binary");
                                    auto binvar2 = static_pointer_cast<var<int>>(binvar_ptr2->second);
                                    indices partns2("partns2");
                                    for (int i = 0; i < num_partns2 ; ++i)
                                    {
                                        partns2.add(name2+ "{" + to_string(i+1) + "}");
                                    }
                                    param<int> lb2("lb2"), ub2("ub2");
                                    lb2.in(o2_ids_uq,partns2);
                                    ub2.in(o2_ids_uq,partns2);
                                    lb2.set_val(0), ub2.set_val(1);
                                    auto added2 = binvar2->add_bounds(lb2,ub2);
                                    reindex_vars();
                                    
                                    auto nb_entries_v1 = o1_ids.get_nb_entries();
                                    auto nb_entries_v2 = o2_ids.get_nb_entries();
                                    auto nb_entries = added.get_nb_entries();
                                    
                                    if(!added1.empty()){
                                        Constraint<> onSum1(o1._name+"_binarySum");
                                        onSum1 = sum(binvar1->in(added1).in_matrix(nb_entries_v1,1));
                                        auto vset1 = added1.from_ith(0,nb_entries_v1);
                                        vset1.filter_refs(vset1.get_unique_refs());
                                        add(onSum1.in(vset1) == 1);
                                    }
                                    
                                    if(!added2.empty()){
                                        Constraint<> onSum2(o2._name+"_binarySum");
                                        onSum2 = sum(binvar2->in(added2).in_matrix(nb_entries_v2,1));
                                        auto vset2 = added2.from_ith(0,nb_entries_v2);
                                        vset2.filter_refs(vset2.get_unique_refs());
                                        add(onSum2.in(vset2) == 1);
                                    }
                                    
                                    if(model_type == "on/off"){//if on/off is chosen
                                        auto binvar_ptr3 = _vars_name.find(name1+name2+"_binary");
                                        auto binvar3 = static_pointer_cast<var<int>>(binvar_ptr3->second);
                                        
                                        indices partns("partns");
                                        partns = indices(partns1,partns2);
                                        auto inst_partition = indices(added,partns);
                                        auto total_entries = inst_partition.get_nb_entries();
                                        
                                        param<int> lb3("lb3"), ub3("ub3");
                                        lb3.in(added,partns);
                                        ub3.in(added,partns);
                                        lb3.set_val(0), ub3.set_val(1);
                                        auto added3 = binvar3->add_bounds(lb3,ub3);
                                        reindex_vars();
                                        
                                        Constraint<> onLink1(pair.first+"_binaryLink1");
                                        onLink1 = binvar1->from_ith(0,inst_partition.ignore_ith(nb_entries_v1, nb_entries_v2)) - binvar3->in(inst_partition);
                                        add(onLink1.in(inst_partition) >= 0);
                                        
                                        Constraint<> onLink2(pair.first+"_binaryLink2");
                                        onLink2 = binvar2->in_ignore_ith(nb_entries_v2,1,inst_partition.ignore_ith(0,nb_entries_v1)) - binvar3->in(inst_partition);
                                        add(onLink2.in(inst_partition) >= 0);
                                        
                                        Constraint<> onLink3(pair.first+"_binaryLink3");
                                        onLink3 = binvar1->from_ith(0,inst_partition.ignore_ith(nb_entries_v1, nb_entries_v2)) + binvar2->in_ignore_ith(nb_entries_v2,1,inst_partition.ignore_ith(0,nb_entries_v1)) - 1 - binvar3->in(inst_partition);
                                        add(onLink3.in(inst_partition) <= 0);
                                        
                                        Constraint<> onSumComb(pair.first+"_binarySum");
                                        onSumComb = sum((binvar3->in(added3)).in_matrix(nb_entries,total_entries-nb_entries));
                                        add(onSumComb.in(added) == 1);
                                        
                                        add_on_off_McCormick_refined(pair.first, vlift->in(added), o1.in(o1_ids), o2.in(o2_ids), binvar3->in(added3));
                                    }
                                    
                                    else{//means it is one of the lambda formulations
                                        
                                        //difference is this has one more partition index
                                        indices partns1_lambda("partns1_lambda");
                                        for (int i = 0; i < num_partns1+1; ++i)
                                        {
                                            partns1_lambda.add(name1+ "{" +to_string(i+1) + "}");
                                        }
                                        
                                        indices partns2_lambda("partns2_lambda");
                                        for (int i = 0; i < num_partns2+1; ++i)
                                        {
                                            partns2_lambda.add(name2+ "{" +to_string(i+1) + "}");
                                        }
                                        
                                        indices partns_lambda("partns_lambda");
                                        partns_lambda = indices(partns1_lambda,partns2_lambda);
                                        auto inst_partition_lambda = indices(added,partns_lambda);
                                        auto inst_partition_bounds1 = indices(added,partns1_lambda);
                                        auto inst_partition_bounds2 = indices(added,partns2_lambda);
                                        
                                        // Convex combination variables
                                        auto lambda_ptr = _vars_name.find(name1+name2+"_lambda");
                                        auto lambda = static_pointer_cast<var<double>>(lambda_ptr->second);
                                        param<double> lb_lambda("lb_lambda"), ub_lambda("ub_lambda");
                                        lb_lambda.in(added,partns_lambda);
                                        ub_lambda.in(added,partns_lambda);
                                        lb_lambda.set_val(0), ub_lambda.set_val(1);
                                        auto added_lambda = lambda->add_bounds(lb_lambda,ub_lambda);
                                        reindex_vars();
                                        
                                        /** Parameters */
                                        // Bounds on variable v1 & v2
                                        param<> bounds1(name1+"_bounds1");
                                        bounds1.in(inst_partition_bounds1);
                                        
                                        param<> bounds2(name2+"_bounds2");
                                        bounds2.in(inst_partition_bounds2);
                                        
                                        // Function values on the extreme points
                                        param<> EP(name1+name2+"_grid_values");
                                        EP.in(inst_partition_lambda);
                                        auto total_entries = inst_partition_lambda.get_nb_entries();
                                        
                                        size_t nb_ins = vlift->in(added).get_nb_inst();
                                        auto o1_global_lb = o1.get_lb();
                                        auto increment1 = (o1.get_ub() - o1_global_lb)/num_partns1;
                                        
                                        auto o2_global_lb = o2.get_lb();
                                        auto increment2 = (o2.get_ub() - o2_global_lb)/num_partns2;
                                        
                                        // fill bounds1 and function values
                                        for (int i=0 ; i<num_partns1+1; ++i) {
                                            auto bound_partn1 = o1_global_lb + increment1*i;
                                            bound_partn1.eval_all();
                                            for (size_t inst = 0; inst< nb_ins; inst++){
                                                auto cur_var_id = vlift->get_id_inst(inst);
                                                auto cur_var_idx = added._keys->at(cur_var_id);
                                                string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"}";
                                                bounds1.set_val(cur_idx,bound_partn1.eval(inst));
                                                for(int j=0; j<num_partns2+1; ++j){
                                                    auto bound_partn2 = o2_global_lb + increment2*j;
                                                    bound_partn2.eval_all();
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"}";
                                                    EP.set_val(cur_idx,(bound_partn1.eval(inst)*bound_partn2.eval(inst)));
                                                }
                                            }
                                        }
                                        // fill bounds2
                                        for (int i=0 ; i<num_partns2+1; ++i) {
                                            auto bound_partn2 = o2_global_lb + increment2*i;
                                            bound_partn2.eval_all();
                                            for (size_t inst = 0; inst< nb_ins; inst++){
                                                auto cur_var_id = vlift->get_id_inst(inst);
                                                auto cur_var_idx = added._keys->at(cur_var_id);
                                                string cur_idx = cur_var_idx+","+name2+"{"+to_string(i+1)+"}";
                                                bounds2.set_val(cur_idx,bound_partn2.eval(inst));
                                            }
                                        }
                                        
                                        // Lambda coefficient matrix when linking with partition variables
                                        param<> lambda_coef1(name1+"_lambda_linking_coefficients1");
                                        param<> lambda_coef2(name2+"_lambda_linking_coefficients2");
                                        
                                        // Partition coefficient matrix when linking with lambda variables
                                        param<> on_coef1(name1+"_partition_linking_coefficients1");
                                        param<> on_coef2(name2+"_partition_linking_coefficients2");
                                        
                                        // create constraint indices
                                        indices const_idx1("const_idx1");
                                        indices const_idx2("const_idx2");
                                        
                                        if(model_type == "lambda_II"){
                                            
                                            //fill constraint indices
                                            for (int i=0; i<num_partns1+1; ++i){
                                                const_idx1.add(to_string(i+1));
                                            }
                                            
                                            //fill constraint indices
                                            for (int i=0; i<num_partns2+1; ++i){
                                                const_idx2.add(to_string(i+1));
                                            }
                                            
                                            // Lambda coefficient matrix when linking with partition variables
                                            if(num_partns1 > 1) lambda_coef1.in(indices(inst_partition_lambda, const_idx1));
                                            if(num_partns2 > 1) lambda_coef2.in(indices(inst_partition_lambda, const_idx2));
                                            
                                            // Lambda coefficient matrix when linking with partition variables
                                            if(num_partns1 > 1) on_coef1.in(indices(added, partns1, const_idx1));
                                            if(num_partns2 > 1) on_coef2.in(indices(added, partns2, const_idx2));
                                            
                                            // fill lambda_coef1 and lambda_coef2
                                            for (size_t inst = 0; inst< nb_ins; inst++){
                                                auto cur_var_id = vlift->get_id_inst(inst);
                                                auto cur_var_idx = added._keys->at(cur_var_id);
                                                for (int i=0 ; i<num_partns1+1; ++i) {
                                                    for (int j=0 ; j<num_partns2+1; ++j) {
                                                        string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                        if(num_partns1 > 1) lambda_coef1.set_val(cur_idx,1);
                                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(j+1);
                                                        if(num_partns2 > 1) lambda_coef2.set_val(cur_idx,1);
                                                    }
                                                }
                                            }
                                            
                                            // fill on_coef1 and on_coef2
                                            for (size_t inst = 0; inst< nb_ins; inst++){
                                                auto cur_var_id = vlift->get_id_inst(inst);
                                                auto cur_var_idx = added._keys->at(cur_var_id);
                                                string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+to_string(1);
                                                if(num_partns1 > 1) on_coef1.set_val(cur_idx,1);
                                                cur_idx = cur_var_idx+","+name2+"{"+to_string(1)+"},"+to_string(1);
                                                if(num_partns2 > 1) on_coef2.set_val(cur_idx,1);
                                                if(num_partns1 > 1) {
                                                    for (int i=1 ; i<num_partns1; ++i) {
                                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(i)+"},"+to_string(i+1);
                                                        on_coef1.set_val(cur_idx,1);
                                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+to_string(i+1);
                                                        on_coef1.set_val(cur_idx,1);
                                                    }
                                                }
                                                if(num_partns2 > 1) {
                                                    for (int i=1 ; i<num_partns2; ++i) {
                                                        cur_idx = cur_var_idx+","+name2+"{"+to_string(i)+"},"+to_string(i+1);
                                                        on_coef2.set_val(cur_idx,1);
                                                        cur_idx = cur_var_idx+","+name2+"{"+to_string(i+1)+"},"+to_string(i+1);
                                                        on_coef2.set_val(cur_idx,1);
                                                    }
                                                }
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(num_partns1)+"},"+to_string(num_partns1+1);
                                                if(num_partns1 > 1) on_coef1.set_val(cur_idx,1);
                                                cur_idx = cur_var_idx+","+name2+"{"+to_string(num_partns2)+"},"+to_string(num_partns2+1);
                                                if(num_partns2 > 1) on_coef2.set_val(cur_idx,1);
                                            }
                                        }
                                        
                                        
                                        else /*means model_type == "lambda_III" */{
                                            
                                            //fill constraint indices
                                            for (int i=0; i<(num_partns1-2)*2+2; ++i){
                                                const_idx1.add(to_string(i+1));
                                            }
                                            
                                            //fill constraint indices
                                            for (int i=0; i<(num_partns2-2)*2+2; ++i){
                                                const_idx2.add(to_string(i+1));
                                            }
                                            
                                            // Lambda coefficient matrix when linking with partition variables
                                            if(num_partns1 > 1) lambda_coef1.in(indices(inst_partition_lambda, const_idx1));
                                            if(num_partns2 > 1) lambda_coef2.in(indices(inst_partition_lambda, const_idx2));
                                            
                                            // Partition coefficient matrix when linking with lambda variables
                                            if(num_partns1 > 1) on_coef1.in(indices(added, partns1, const_idx1));
                                            if(num_partns2 > 1) on_coef2.in(indices(added, partns2, const_idx2));
                                            
                                            // fill lambda_coef1 and lambda_coef2
                                            for (size_t inst = 0; inst< nb_ins; inst++){
                                                auto cur_var_id = vlift->get_id_inst(inst);
                                                auto cur_var_idx = added._keys->at(cur_var_id);
                                                if(num_partns1 > 1) {
                                                    for (int j=0; j<num_partns2+1; ++j) {
                                                        string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(1);
                                                        lambda_coef1.set_val(cur_idx,1);
                                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(num_partns1+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string((num_partns1-2)*2+2);
                                                        lambda_coef1.set_val(cur_idx,1);
                                                    }
                                                    
                                                    for (int i=1 ; i<(num_partns1-2)*2+1; i=i+2) {
                                                        for (int j=(i-1)/2 + 2; j<num_partns1+1; ++j) {
                                                            for(int k=0; k<num_partns2+1; ++k){
                                                                string cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+name2+"{"+to_string(k+1)+"},"+to_string(i+1);
                                                                lambda_coef1.set_val(cur_idx,1);
                                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+name2+"{"+to_string(k+1)+"},"+to_string(i+2);
                                                                lambda_coef1.set_val(cur_idx,-1);
                                                            }
                                                        }
                                                    }
                                                }
                                                if(num_partns2 > 1) {
                                                    for (int i=0; i<num_partns1+1; ++i) {
                                                        string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(1)+"},"+to_string(1);
                                                        lambda_coef2.set_val(cur_idx,1);
                                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(num_partns2+1)+"},"+to_string((num_partns2-2)*2+2);
                                                        lambda_coef2.set_val(cur_idx,1);
                                                    }
                                                    
                                                    for (int i=1 ; i<(num_partns2-2)*2+1; i=i+2) {
                                                        for (int j=(i-1)/2 + 2; j<num_partns2+1; ++j) {
                                                            for(int k=0; k<num_partns1+1; ++k){
                                                                string cur_idx = cur_var_idx+","+name1+"{"+to_string(k+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                                lambda_coef2.set_val(cur_idx,1);
                                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(k+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                                lambda_coef2.set_val(cur_idx,-1);
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                            
                                            
                                            
                                            // fill on_coef1 and on_coef2
                                            for (size_t inst = 0; inst< nb_ins; inst++){
                                                auto cur_var_id = vlift->get_id_inst(inst);
                                                auto cur_var_idx = added._keys->at(cur_var_id);
                                                string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+to_string(1);
                                                if(num_partns1 > 1) {
                                                    on_coef1.set_val(cur_idx,1);
                                                    
                                                    for (int i=1; i<num_partns1; ++i) {
                                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+to_string(2);
                                                        on_coef1.set_val(cur_idx, 1);
                                                    }
                                                    
                                                    for (int i=2 ; i<(num_partns1-2)*2+2; i=i+2) {
                                                        for (int j=i/2+1; j<num_partns1; ++j) {
                                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                            on_coef1.set_val(cur_idx,-1);
                                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                            on_coef1.set_val(cur_idx,1);
                                                        }
                                                    }
                                                }
                                                if(num_partns2 > 1) {
                                                    cur_idx = cur_var_idx+","+name2+"{"+to_string(1)+"},"+to_string(1);
                                                    on_coef2.set_val(cur_idx,1);
                                                    for (int i=1; i<num_partns2; ++i) {
                                                        cur_idx = cur_var_idx+","+name2+"{"+to_string(i+1)+"},"+to_string(2);
                                                        on_coef2.set_val(cur_idx, 1);
                                                    }
                                                    
                                                    for (int i=2 ; i<(num_partns2-2)*2+2; i=i+2) {
                                                        for (int j=i/2+1; j<num_partns2; ++j) {
                                                            cur_idx = cur_var_idx+","+name2+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                            on_coef2.set_val(cur_idx,-1);
                                                            cur_idx = cur_var_idx+","+name2+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                            on_coef2.set_val(cur_idx,1);
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                        
                                        
                                        /** Constraints */
                                        // Representation of the bilinear term with convex combination
                                        Constraint<> bln_rep(pair.first+"_bln_rep");
                                        bln_rep = EP.in_matrix(nb_entries,total_entries-nb_entries)*lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries) - vlift->in(added);
                                        add(bln_rep.in(added) == 0);
                                        
                                        // Representation of o1 with convex combination
                                        Constraint<> o1_rep(pair.first+"_o1_rep");
                                        o1_rep = bounds1.from_ith(0,inst_partition_lambda).in_matrix(nb_entries, 1) * lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries) - o1.in(o1_ids);
                                        add(o1_rep.in(added) == 0);
                                        
                                        // Representation of o2 with convex combination
                                        Constraint<> o2_rep(pair.first+"_o2_rep");
                                        o2_rep = bounds2.in_ignore_ith(nb_entries, 1, inst_partition_lambda).in_matrix(nb_entries,1) * lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries) - o2.in(o2_ids);
                                        add(o2_rep.in(added) == 0);
                                        
                                        // Linking partition variables1 with lambda
                                        if(model_type == "lambda_II"){
                                            if(num_partns1 > 1) {
                                                Constraint<> on_link_lambda1(pair.first+"_on_link_lambda1_II");
                                                on_link_lambda1 = lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef1.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef1.in_matrix(nb_entries,total_entries-nb_entries) - binvar1->in_ignore_ith(nb_entries_v1,nb_entries_v2,on_coef1.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef1.in_matrix(nb_entries,1);
                                                add(on_link_lambda1.in(indices(added,const_idx1)) <= 0);
                                            }
                                            if(num_partns2 > 1) {
                                                Constraint<> on_link_lambda2(pair.first+"_on_link_lambda2_II");
                                                on_link_lambda2 = lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef2.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef2.in_matrix(nb_entries,total_entries-nb_entries) - binvar2->in_ignore_ith(0,nb_entries_v1,on_coef2.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef2.in_matrix(nb_entries,1);
                                                add(on_link_lambda2.in(indices(added,const_idx2)) <= 0);
                                            }
                                        }
                                        else{
                                            if(num_partns1 > 1) {
                                                Constraint<> on_link_lambda1(pair.first+"_on_link_lambda1_III");
                                                on_link_lambda1 = lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef1.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef1.in_matrix(nb_entries,total_entries-nb_entries) - binvar1->in_ignore_ith(nb_entries_v1,nb_entries_v2,on_coef1.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef1.in_matrix(nb_entries,1);
                                                add(on_link_lambda1.in(indices(added,const_idx1)) <= 0);
                                            }
                                            if(num_partns2 > 1) {
                                                Constraint<> on_link_lambda2(pair.first+"_on_link_lambda2_III");
                                                on_link_lambda2 = lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef2.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef2.in_matrix(nb_entries,total_entries-nb_entries) - binvar2->in_ignore_ith(0,nb_entries_v1,on_coef2.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef2.in_matrix(nb_entries,1);
                                                add(on_link_lambda2.in(indices(added,const_idx2)) <= 0);
                                            }
                                        }
                                        // sum over lambda
                                        Constraint<> lambdaSum(pair.first+"_lambdaSum");
                                        lambdaSum = sum(lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries));
                                        add(lambdaSum.in(added) == 1);
                                    }
                                }
                            }
                        }
                        else {
                            add_McCormick(pair.first, vlift->in(added), o1.in(o1_ids), o2.in(o2_ids));
                        }
                    }
                }
                
                lifted.insert(lt);
            }
            for (auto &pair:*c._pterms) {
                auto term = pair.second;
                lterm lt;
                lt._sign = term._sign;
                if (term._coef->is_function()) {
                    auto coef = *static_pointer_cast<func<type>>(term._coef);
                    lt._coef = func<type>(coef).copy();
                }
                else if(term._coef->is_param()) {
                    auto coef = *static_pointer_cast<param<type>>(term._coef);
                    lt._coef = param<type>(coef).copy();
                }
                else if(term._coef->is_number()) {
                    auto coef = *static_pointer_cast<constant<type>>(term._coef);
                    lt._coef = constant<type>(coef).copy();
                }
                func<type> prod = 1;
                string prod_name = "Lift(";
                auto list = pair.second._l;
                for (auto &ppi: *list) {
                    auto p = ppi.first;
                    auto orig_var = *static_pointer_cast<var<type>>(p);
                    if(ppi.second>1){
                        prod_name += orig_var.get_name(true,true)+"("+orig_var._indices->get_name()+")^"+to_string(ppi.second);
                        //TODO Lift univarite power function
                    }
                    else{
                        prod_name += orig_var.get_name(true,true)+"("+orig_var._indices->get_name()+")";
                    }
                    prod *= pow(orig_var,ppi.second);
                }
                prod_name += ")";
                
                auto ids = *c._indices;
                param<type> lb("lb"), ub("ub");
                lb.in(ids);ub.in(ids);
                lb.set_val(prod._range->first);
                ub.set_val(prod._range->second);
                var<type> vlift(prod_name, lb, ub);
                auto it = _vars_name.find(prod_name);
                if(it==_vars_name.end()){
                    vlift._lift=true;
                    add(vlift.in(ids));
                    lt._p = make_shared<var<type>>(vlift);
                }
                else {
                    vlift = *static_pointer_cast<var<type>>(it->second);
                    lt._p = make_shared<var<type>>(vlift);
                }
                lifted.insert(lt);
            }
            lifted._range = c._range;
            lifted._all_convexity = linear_;
            lifted._all_sign = c._all_sign;
            lifted._ftype = lin_;
            lifted._ctype = c._ctype;
            lifted._indices = c._indices;
            lifted._dim[0] = c._dim[0];
            lifted._dim[1] = c._dim[1];
            return lifted;
        }
        
        template<typename T=type,typename std::enable_if<is_arithmetic<T>::value>::type* = nullptr>
        void add_real(const Constraint<Cpx>& c){
            if (c.get_dim()==0) {
                return;
            }
            auto real_imag = get_real_imag(c);
            Constraint<type> c_real;
            c_real += real_imag.first;
            c_real._name = "Real(" + c._name + ")";
            c_real._ctype = c._ctype;
            c_real._indices = c._indices;
            c_real._dim[0] = c._dim[0];
            add_constraint(c_real);
        }
        
        template<typename T=type,typename std::enable_if<is_arithmetic<T>::value>::type* = nullptr>
        void add_imag(const Constraint<Cpx>& c){
            if (c.get_dim()==0) {
                return;
            }
            auto real_imag = get_real_imag(c);
            Constraint<type> c_imag;
            c_imag += real_imag.second;
            c_imag._name = "Imag(" + c._name + ")";
            c_imag._ctype = c._ctype;
            c_imag._indices = c._indices;
            c_imag._dim[0] = c._dim[0];
            add_constraint(c_imag);
        }
        
        template<typename T=type,typename std::enable_if<is_arithmetic<T>::value>::type* = nullptr>
        void add(Constraint<Cpx>&& c, bool convexify = false, string model_type = "on/off"){
            if (c.get_dim()==0) {
                return;
            }
            auto real_imag = get_real_imag(c);
            Constraint<type> c_real = real_imag.first;
            c_real._name = "Real(" + c._name + ")";
            c_real._ctype = c._ctype;
            c_real._indices = c._indices;
            Constraint<type> c_imag = real_imag.first;
            c_imag._name = "Imag(" + c._name + ")";
            c_imag._ctype = c._ctype;
            c_imag._indices = c._indices;
            if(convexify){
                add(lift(c_real, model_type));
                add(lift(c_imag, model_type));
            }
            else {
                add_constraint(c_real);
                add_constraint(c_imag);
            }
        }
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        void add(Constraint<type>& c, bool convexify = false, string method_type = "on/off"){
            if (c.get_dim()==0) {
                return;
            }
            add_constraint(c,convexify, method_type);
        }
        
        
        void add_lazy(Constraint<type>& c, bool convexify = false){
            if (c.get_dim()==0) {
                return;
            }
            c.make_lazy();
            add_constraint(c, convexify);
            _has_lazy = true;
        }
        
        void add_callback(Constraint<type>& c){/**<  Adds a callback to the model with given constraints (maybe parameters and variables are needed!). */
            if (c.get_dim()==0) {
                return;
            }
            //            c.make_lazy();   //Should we do a similar thing like this?
            add_constraint(c);
            _has_callback = true;
        }
        
        //        template<typename T>
        //        void replace(const shared_ptr<param_>& v, func<T>& f){/**<  Replace v with function f everywhere it appears */
        //            for (auto &c_p: _cons_name) {
        //                auto c = c_p.second;
        //                if (!c->has_var(*v)) {
        //                    continue;
        //                }
        //                c->replace(v, f);
        //            }
        //            _vars_name.erase(v->_name);
        //            auto vid = *v->_vec_id;
        //            delete _vars.at(vid);
        //            _vars.erase(vid);
        //            reindex_vars();
        //        }
        
        
        //        void project() {/**<  Use the equations where at least one variable appear linearly to express it as a function of other variables in the problem */
        //            for (auto& c_pair:_cons_name) {
        //                if (!c_pair.second->is_ineq()) {
        //                    auto &lterms = c_pair.second->get_lterms();
        //                    if (!lterms.empty()) {
        //                        auto first = lterms.begin();
        //                        auto v = first->second._p;
        //                        if (v->_is_vector) {
        //                            continue;
        //                        }
        //                        auto f = *c_pair.second;
        //                        if (first->second._sign) {
        //                            //                    f -= *v;
        //                            //                    f *= -1;
        //                        }
        //                        else {
        //                            //                    f += *v;
        //                        }
        //                        DebugOff(f.to_str());
        //                        _cons.erase(c_pair.second->_id);
        //                        _cons_name.erase(c_pair.first);
        //                        replace(v,f);
        //                        //                project();
        //                        return;
        //                    }
        //                }
        //            }
        //        }
        
        
        /** Add constraint to model
         @param[in] c: constraint to add
         @param[in] lift_flag: if true, add a linearized version of this constraint and the constraints linking the lifted variables to the original ones.
         @return a pointer to the added constraint
         @note If lift = true, this function will add constraints linking the lifted variables to the original ones, if a variable's partition is greater than 1, it will also add the disjunctive constraints corresponding to the partitionning of the variables. Note also that if this is an equation f(x) = 0 s.t. f(x)<=0 or f(x)>=0 is convex, will add the convex inequality to the model.
         **/
        shared_ptr<Constraint<type>> add_constraint(Constraint<type>& c, bool lift_flag = false, string method_type = "on/off"){
            if (c.get_dim()==0) {
                return nullptr;
            }
            if (_cons_name.count(c.get_name())==0) {
                auto newc = make_shared<Constraint<type>>(c);
                for (auto &vp: *newc->_vars) {
                    _v_in_cons[vp.second.first->_name].insert(newc);
                }
                newc->_val = c._val;
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
                            if (newc->is_positive() || newc->is_negative()) {
                                throw invalid_argument("Adding violated constant equation!\n");
                            }
                            break;
                        default:
                            break;
                    }
                    Warning("WARNING: Adding redundant constant constraint, Gravity will be ignoring it.\n");
                    newc->print();
                    return newc;
                }
                newc->update_str();
                if(lift_flag){
                    if(newc->func<type>::is_convex() && newc->_ctype==eq){
                        DebugOn("Convex left hand side of equation detected, splitting constraint into <= and ==" << endl);
                        Constraint<type> c_cvx(*newc);
                        c_cvx._name = newc->_name+"_convex";
                        c_cvx._relaxed = true;
                        add_constraint(c_cvx <= 0);
                    }
                    if(newc->func<type>::is_concave() && newc->_ctype==eq){
                        DebugOn("Concave left hand side of equation detected, splitting constraint into >= and ==" << endl);
                        Constraint<type> c_ccve(*newc);
                        c_ccve._name = newc->_name+"_concave";
                        c_ccve._relaxed = true;
                        add_constraint(c_ccve >= 0);
                    }
                    
                    /* Call the lift function and add corresponding constraints */
                    auto lifted = lift(*newc, method_type);
                    return add_constraint(lifted);
                }
                    
                embed(newc, false);
                update_convexity(*newc);
                newc->_violated.resize(newc->get_nb_inst(),true);
                _cons_name[c.get_name()] = newc;
                if(*newc->_all_lazy && newc->_lazy.size()==0){
                    newc->_lazy.resize(newc->get_nb_inst(),true);
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
                    c.allocate_mem();
                }
                return newc;
            }
            else {
                DebugOn("constraint with same name: " << c.get_name() << endl);
                c.update_str();
                if(!_cons_name.at(c.get_name())->equal(c)) {
                    throw invalid_argument("rename constraint as this name has been used by another one: " + c.get_name());
                }
                else{
                    DebugOn("Both constraints are identical, ignoring this one." << endl);
                    return nullptr;
                }
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
                else {
                    _convexity = undet_;
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
            _obj->_indices = nullptr;
            update_convexity(f);
            if (_type==lin_m && f.is_quadratic()) {
                _type = quad_m;
            }
            else if ((_type==lin_m || _type==quad_m) && f.is_polynomial()) {
                _type = pol_m;
            }
            else if(f.is_nonlinear()){
                _type = nlin_m;
            }
            embed(_obj);
        }
        
        template<typename T1>
        void min(const param<T1>& p){
            set_objective(func<T1>(p), minimize);
        }
        
        template<typename T1>
        void min(T1 c){
            set_objective(func<T1>(c), minimize);
        }
        
        template<typename T1>
        void min(const var<T1>& v){
            set_objective(func<T1>(v), minimize);
        }
        
        template<typename T1>
        void min(const func<T1>& f){
            set_objective(f, minimize);
        }
        
        
        template<typename T1>
        void max(const param<T1>& p){
            set_objective(func<T1>(p), maximize);
        }
        
        template<typename T1>
        void max(T1 c){
            set_objective(func<T1>(c), maximize);
        }
        
        template<typename T1>
        void max(const var<T1>& v){
            set_objective(func<T1>(v), maximize);
        }
        
        template<typename T1>
        void max(const func<T1>& f){
            set_objective(f, maximize);
        }
        
        
        void set_objective_type(ObjectiveType t) {
            _objt = t;
        }
        
        /** Prints the constraints that have a non-zero value at the current solution where zero <= tol
         @param[tol] zero value tolerance
         @param[only_relaxed] if set to true, only print non-zero values of inequalities that were generated by relaxing non-convex equations
         */
        template<typename T=type,typename std::enable_if<is_arithmetic<T>::value>::type* = nullptr>
        void print_nonzero_constraints(double tol, bool only_relaxed = false) const{
            size_t nb_inst = 0;
            shared_ptr<Constraint<type>> c = nullptr;
            for(auto& c_p: _cons_name)
            {
                c = c_p.second;
                if(only_relaxed && !c->_relaxed){
                    continue;
                }
                nb_inst = c->get_nb_inst();
                switch (c->get_ctype()) {
                    case eq:
                        for (size_t inst=0; inst<nb_inst; inst++) {
                            auto diff = std::abs(c->eval(inst));
                            if(diff>tol){
                                DebugOn(c->_name << " Non-zero equation: " << to_string(inst) << ", value = "<< diff << endl);
                            }
                        }
                        break;
                    case leq:
                        for (size_t inst=0; inst<nb_inst; inst++) {
                            auto diff = c->eval(inst);
                            if(diff < -tol) {
                                DebugOn(c->_name << " Non-zero <= inequality: " << to_string(inst) << ", value = "<< diff << endl);
                            }
                        }
                        break;
                    case geq:
                        for (size_t inst=0; inst<nb_inst; inst++) {
                            auto diff = c->eval(inst);
                            if(diff > tol) {
                                DebugOn(c->_name << " Non-zero >= inequality: " << to_string(inst) << ", value = "<< diff << endl);
                            }
                        }
                        break;
                        
                    default:
                        break;
                }
            }
        }
        
        template<typename T=type,typename std::enable_if<is_arithmetic<T>::value>::type* = nullptr>
        vector<tuple<double, int, int>> sorted_nonzero_constraints(double tol, bool only_relaxed = false, bool print_name = false) const{
            // the tuple has the following form <value, constraint_id, instance_id>
            
            vector<tuple<double, int, int>> v;
            size_t nb_inst = 0;
            shared_ptr<Constraint<type>> c = nullptr;
            for(auto& c_p: _cons_name)
            {
                c = c_p.second;
                if(only_relaxed && !c->_relaxed){
                    continue;
                }
                nb_inst = c->get_nb_inst();
                switch (c->get_ctype()) {
                    case eq:
                        for (size_t inst=0; inst<nb_inst; inst++) {
                            auto diff = std::abs(c->eval(inst));
                            if(diff>tol){
                                v.push_back(make_tuple(std::abs(diff), c->_id, inst));
                                if(print_name) DebugOn(" Non-zero >= inequality: " << c->_name << " instance: " << to_string(inst) << ", value = "<< std::abs(diff) << endl);
                            }
                        }
                        break;
                    case leq:
                        for (size_t inst=0; inst<nb_inst; inst++) {
                            auto diff = c->eval(inst);
                            if(diff < -tol) {
                                v.push_back(make_tuple(std::abs(diff), c->_id, inst));
                                if(print_name) DebugOn(" Non-zero >= inequality: " << c->_name << " instance: " << to_string(inst) << ", value = "<< std::abs(diff) << endl);
                            }
                        }
                        break;
                    case geq:
                        for (size_t inst=0; inst<nb_inst; inst++) {
                            auto diff = c->eval(inst);
                            if(diff > tol) {
                                v.push_back(make_tuple(std::abs(diff), c->_id, inst));
                                if(print_name) DebugOn(" Non-zero >= inequality: " << c->_name << " instance: " << to_string(inst) << ", value = "<< std::abs(diff) << endl);
                            }
                        }
                        break;
                        
                    default:
                        break;
                }
            }
            sort(v.begin(), v.end(), std::greater<tuple<double,int,int>>());
            return v;
        }
        
        indices sorted_nonzero_constraint_indices(double tol, bool print_name, string constraint_name) const{
            // returns the indices of the constraint given a constraint_name
            
            vector<tuple<double, string>> v; //violation amount & the index of a given constraint
            size_t nb_inst = 0;
            shared_ptr<Constraint<type>> c = nullptr;
            for(auto& c_p: _cons_name)
            {
                c = c_p.second;
                if(c->_name != constraint_name){
                    continue;
                }
                
                nb_inst = c->get_nb_inst();
                switch (c->get_ctype()) {
                    case eq:
                        for (size_t inst=0; inst<nb_inst; inst++) {
                            auto diff = std::abs(c->eval(inst));
                            if(diff>tol){
                                v.push_back(make_tuple(std::abs(diff), c->_indices->_keys->at(inst)));
                                if(print_name) DebugOn(" Non-zero >= inequality: " << c->_name << " instance: " << to_string(inst) << ", value = "<< std::abs(diff) << endl);
                            }
                        }
                        break;
                    case leq:
                        for (size_t inst=0; inst<nb_inst; inst++) {
                            auto diff = c->eval(inst);
                            if(diff < -tol) {
                                v.push_back(make_tuple(std::abs(diff), c->_indices->_keys->at(inst)));
                                if(print_name) DebugOn(" Non-zero >= inequality: " << c->_name << " instance: " << to_string(inst) << ", value = "<< std::abs(diff) << endl);
                            }
                        }
                        break;
                    case geq:
                        for (size_t inst=0; inst<nb_inst; inst++) {
                            auto diff = c->eval(inst);
                            if(diff > tol) {
                                v.push_back(make_tuple(std::abs(diff), c->_indices->_keys->at(inst)));
                                if(print_name) DebugOn(" Non-zero >= inequality: " << c->_name << " instance: " << to_string(inst) << ", value = "<< std::abs(diff) << endl);
                            }
                        }
                        break;
                        
                    default:
                        break;
                }
            }
            sort(v.begin(), v.end(), std::greater<tuple<double,string>>());
            
            // HERE IS THE PART TO COLLECT THE SORTED INDEX SET
            indices nonzero_idx("nonzero_idx"); //the indices of the nonzero_constraint instances
            for (int i = 0; i < v.size(); i++)
                nonzero_idx.add(get<1>(v[i]));
            
            return nonzero_idx;
        }
        
        
        template<typename T=type,typename std::enable_if<is_arithmetic<T>::value>::type* = nullptr>
        bool print_constraints_stats(type tol){/*<< Prints stats on constraints active status with the current solution and tolerance tol */
            size_t nb_inst = 0, nb_viol = 0, nb_viol_all = 0;
            size_t nb_active = 0, nb_active_all = 0;
            double diff = 0;
            shared_ptr<Constraint<type>> c = nullptr;
            bool violated = false;
            for(auto& c_p: _cons_name)
            {
                c = c_p.second;
                //        cid = c->_id;
                nb_inst = c->get_nb_inst();
                nb_viol = 0;
                nb_active = 0;
                c->_all_satisfied = true;
                c->_violated.resize(nb_inst);
                c->_active.resize(nb_inst);
                switch (c->get_ctype()) {
                    case eq:
                        for (size_t inst=0; inst<nb_inst; inst++) {
                            diff = std::abs(c->eval(inst));
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
                                //                                c->print(inst);
                                DebugOn(", violation = "<< diff << endl);
                                nb_viol++;
                                //                        violated = true;
                                if (*c->_all_lazy) {
                                    //                                    *c->_all_lazy = false;
                                    c->_all_satisfied = false;
                                    c->_violated[inst] = true;
                                    violated = true;
                                    c->_lazy[inst] = false;
                                }
                                else {
                                    //                            throw runtime_error("Non-lazy constraint is violated, solution declared optimal by solver!\n" + c->to_str(inst));
                                }
                            }
                            else if (std::abs(diff)>tol) {
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
                                DebugOn("Violated inequality: ");
                                //                        c->print(inst);
                                DebugOn(", violation = "<< diff << endl);
                                nb_viol++;
                                //                        violated = true;
                                if (*c->_all_lazy) {
                                    //                                    *c->_all_lazy = false;
                                    c->_all_satisfied = false;
                                    c->_violated[inst] = true;
                                    violated = true;
                                    c->_lazy[inst] = false;
                                }
                                else {
                                    //                            throw runtime_error("Non-lazy constraint is violated, solution declared optimal by solver!\n" + c->to_str(inst));
                                }
                            }
                            else if (std::abs(diff)> tol) {
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
                if (nb_viol>0 && c->get_ctype()!=eq) {
                    DebugOff("Percentage of violated constraints for " << c->get_name() << " = (" << nb_viol << "/" << nb_inst << ") " << to_string_with_precision(100.*nb_viol/nb_inst,3) << "%\n");
                }
                if (c->get_ctype()!=eq) {
                    DebugOn("Percentage of active constraints for " << c->get_name() << " = (" << nb_active << "/" << nb_inst << ") " << to_string_with_precision(100.*nb_active/nb_inst,3) << "%\n");
                }
            }
            auto nb_ineq = get_nb_ineq();
            DebugOn("Total percentage of violated constraints = (" << nb_viol_all << "/" << nb_ineq << ") " << to_string_with_precision(100.*nb_viol_all/nb_ineq,3) << "%\n");
            DebugOn("Total percentage of active constraints = (" << nb_active_all << "/" << nb_ineq << ") "  << to_string_with_precision(100.*nb_active_all/nb_ineq,3) << "%\n");
            return violated;
        }
        
        template<typename T=type,typename std::enable_if<is_arithmetic<T>::value>::type* = nullptr>
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
                nb_inst = c->get_nb_inst();
                nb_viol = 0;
                nb_active = 0;
                c->_all_satisfied = true;
                c->_violated.resize(nb_inst);
                c->_active.resize(nb_inst);
                switch (c->get_ctype()) {
                    case eq:
                        for (size_t inst=0; inst<nb_inst; inst++) {
                            diff = std::abs(c->eval(inst));
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
                                DebugOff("Violated inequality: ");
                                //                                c->print(inst);
                                DebugOff(", violation = "<< diff << endl);
                                nb_viol++;
                                //                        violated = true;
                                if (*c->_all_lazy) {
                                    //                                    *c->_all_lazy = false;
                                    c->_all_satisfied = false;
                                    c->_violated[inst] = true;
                                    violated = true;
                                    c->_lazy[inst] = false;
                                }
                                else {
                                    //                            throw runtime_error("Non-lazy constraint is violated, solution declared optimal by solver!\n" + c->to_str(inst));
                                }
                            }
                            else if (std::abs(diff)>tol) {
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
                                    //                                    *c->_all_lazy = false;
                                    c->_all_satisfied = false;
                                    c->_violated[inst] = true;
                                    violated = true;
                                    c->_lazy[inst] = false;
                                }
                                else {
                                    //                            throw runtime_error("Non-lazy constraint is violated, solution declared optimal by solver!\n" + c->to_str(inst));
                                }
                            }
                            else if (std::abs(diff)> tol) {
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
                if (nb_viol>0 && c->get_ctype()!=eq) {
                    DebugOff("Percentage of violated constraints for " << c->get_name() << " = (" << nb_viol << "/" << nb_inst << ") " << to_string_with_precision(100.*nb_viol/nb_inst,3) << "%\n");
                }
                if (c->get_ctype()!=eq) {
                    DebugOff("Percentage of active constraints for " << c->get_name() << " = (" << nb_active << "/" << nb_inst << ") " << to_string_with_precision(100.*nb_active/nb_inst,3) << "%\n");
                }
            }
            auto nb_ineq = get_nb_ineq();
            DebugOff("Total percentage of violated constraints = (" << nb_viol_all << "/" << nb_ineq << ") " << to_string_with_precision(100.*nb_viol_all/nb_ineq,3) << "%\n");
            DebugOff("Total percentage of active constraints = (" << nb_active_all << "/" << nb_ineq << ") "  << to_string_with_precision(100.*nb_active_all/nb_ineq,3) << "%\n");
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
                v_p.second->get_double_lb(x_l);
                v_p.second->get_double_ub(x_u);
            }
        }
        
        
        /**
         Initialize the model variables using values from x
         @param[in] x values to initialize to
         */
        void set_x(const double* x){
            for(auto &v_p: _vars)
            {
                v_p.second->set_double_val(x);
            }
        }
        
        /**
         Initialize the model variables using values from x
         @param[in] x values to initialize to
         */
        void set_solution(const vector<double>& x){
            for(auto &v_p: _vars)
            {
                v_p.second->set_solution(x);
            }
        }
        
        /**
         Evaluate all nonlinear functions at current point
         */
        void compute_funcs() {
            if(_type!=nlin_m){
                return;
            }
            //            size_t tot_evals = 0;
            auto it = _nl_funcs.begin();
            while (it!=_nl_funcs.end()) {
                auto f = (*it++);
                DebugOff(f->to_str() << endl);
                if (f->is_constant() && f->_evaluated) {
                    continue;
                }
                if (f->is_matrix_indexed()) {
                    
                    
                    f->_evaluated = false;
                    DebugOff(f->to_str()<<endl);
                    
                    //                    cout << " | nb_instances = " << f->get_dim() << endl;
                    for (size_t i = 0; i < f->_indices->_ids->size(); i++) {
                        
                        for (size_t j = 0; j < f->_indices->_ids->at(i).size(); j++) {
                            f->eval(i,j);
                        }
                        //                        tot_evals++;
                    }
                }
                else if (!f->is_matrix()) {
                    f->_evaluated = false;
                    DebugOff(f->to_str()<<endl);
                    //                    cout << " | nb_instances = " << f->get_dim() << endl;
                    for (size_t inst = 0; inst < f->get_nb_inst(); inst++) {
                        f->eval(inst);
                        //                        tot_evals++;
                    }
                }
                else {
                    DebugOff(f->to_str()<<endl);
                    f->_evaluated = false;
                    f->eval_matrix();
                }
                //                if (f->is_constant()) {
                f->_evaluated = true;
                
                //                }
            }
            //            cout << "tot_evals = " << tot_evals << endl;
        }
        
        template<typename T=type,typename std::enable_if<is_arithmetic<T>::value>::type* = nullptr>
        void fill_in_obj(const double* x , double& res, bool new_x){
            if (new_x) {
                set_x(x);
                compute_funcs();
            }
            if(!_obj->is_constant()){
                _obj->_evaluated = false;
            }
            res = _obj->eval();
            _obj->_new = false;
            DebugOff("Objective = " << to_string(res) << endl);
        }
        
        template<typename T=type,typename std::enable_if<is_arithmetic<T>::value>::type* = nullptr>
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
                if(!df->is_constant()){
                    df->_evaluated = false;
                }
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
        
        template<typename T1>
        void partition(std::string name, std::string model_type, const var<T1>& vlift, const var<T1>& v1, const var<T1>& v2, const vector<T1>& p1, const vector<T1>& p2) {
            // *************SHOULD I USE SOMETHING LIKE THIS??? add(MC4.in(*vlift._indices)); ***************//
            if(v1._name!=v2._name)
            {
                if ((model_type == "Model_II") || (model_type == "Model_III")) {
                    
                    // Obtain the number of partitions
                    int num_partition1 = p1.size()-1;  // number of partition on v1
                    int num_partition2 = p2.size()-1;  // number of partition on v2
                    
                    /** Parameters */
                    // Bounds on variable v1
                    param<double> bounds1(name+"_bounds on variable v1");
                    bounds1.in(R(num_partition1+1));
                    for (int i=0; i<num_partition1+1; ++i) {
                        bounds1.set_val(i,p1[i]);
                    }
                    // Bounds on variable v2
                    param<double> bounds2(name+"_bounds on variable v2");
                    bounds2.in(R(num_partition2+1));
                    for (int i=0; i<num_partition2+1; ++i) {
                        bounds2.set_val(i,p2[i]);
                    }
                    // Function values on the extreme points
                    param<double> EP(name+"_function values in the grid");
                    EP.in(R((num_partition1+1)*(num_partition2+1)));
                    for (int j=0; j<num_partition2+1; ++j) {
                        for (int i=0; i<num_partition1+1; ++i) {
                            EP.set_val(i+j*(num_partition1+1),p1[i]*p2[j]);
                        }
                    }
                    
                    // Lambda coefficient matrix when linking with partition variables
                    param<type> lambdaCOEF(name+"_lambda linking coefficients");
                    
                    // Partition coefficient matrix when linking with lambda variables
                    param<type> zCOEF(name+"_partition linking coefficients");
                    
                    // Partition assignment coefficients
                    param<type> zASGNCOEF(name+"_partition assignment coefficients");
                    
                    if (model_type == "Model_II"){
                        // check if both of the variables involve a partition
                        if ((num_partition1>1) && (num_partition2>1)) {
                            
                            // Lambda coefficient matrix when linking with partition variables
                            lambdaCOEF.set_size(num_partition1+num_partition2+2,((num_partition1+1)*(num_partition2+1)));
                            
                            // Coefficients related to partition on x2
                            for (int i=0; i<num_partition2+1; ++i) {
                                for (int j=0; j<num_partition1+1; ++j) {
                                    lambdaCOEF.set_val(i,j+i*(num_partition1+1),1);
                                }
                            }
                            
                            // Collect the corresponding indices of lambda related to partition on v1
                            vector<int> myIndexVector(num_partition2+1); //initial index vector for lambdas bounding partition variables for v1
                            for (int i = 0; i<myIndexVector.size(); ++i) {
                                myIndexVector[i] = (i*(num_partition1+1));
                            }
                            
                            // Coefficients related to partition on v1
                            for (int i=0; i<num_partition1+1; ++i) {
                                for (int j=0; j<num_partition2+1; ++j) {
                                    lambdaCOEF.set_val(i+num_partition2+1,myIndexVector[j]+i,1);
                                }
                            }
                            
                            // Partition coefficient matrix when linking with lambda variables
                            zCOEF.set_size(num_partition1+num_partition2+2, num_partition1+num_partition2);
                            
                            // Coefficients related to partition on v2
                            zCOEF.set_val(0, num_partition1, 1);
                            for (int i=1; i<num_partition2; ++i) {
                                zCOEF.set_val(i,num_partition1+i-1,1);
                                zCOEF.set_val(i,num_partition1+i, 1);
                            }
                            zCOEF.set_val(num_partition2, num_partition1+num_partition2-1, 1);
                            
                            // Coefficients related to partition on v1
                            zCOEF.set_val(num_partition2+1, 0, 1);
                            for (int i=2; i<num_partition1+1; ++i) {
                                zCOEF.set_val(i+num_partition2,i-2,1);
                                zCOEF.set_val(i+num_partition2,i-1,1);
                            }
                            zCOEF.set_val(num_partition1+num_partition2+1, num_partition1-1, 1);
                            
                            // Partition assignment coefficients
                            zASGNCOEF.set_size(2,num_partition1+num_partition2);
                            
                            // Coefficients related to partition on v1
                            for (int i=0; i<num_partition1; ++i) {
                                zASGNCOEF.set_val(0,i,1);
                            }
                            
                            // Coefficients related to partition on v2
                            for (int i=0; i<num_partition2; ++i) {
                                zASGNCOEF.set_val(1,i+num_partition1,1);
                            }
                        }
                        // check if the first variable v1 involves a partition
                        else if (num_partition1>1){
                            // Lambda coefficient matrix when linking with partition variables
                            lambdaCOEF.set_size(num_partition1+1,((num_partition1+1)*(num_partition2+1)));
                            
                            // Collect the corresponding indices of lambda related to partition on v1
                            vector<int> myIndexVector(num_partition2+1); //initial index vector for lambdas bounding partition variables for x1
                            for (int i = 0; i<myIndexVector.size(); ++i) {
                                myIndexVector[i] = (i*(num_partition1+1));
                            }
                            
                            // Coefficients related to partition on v1
                            for (int i=0; i<num_partition1+1; ++i) {
                                for (int j=0; j<num_partition2+1; ++j) {
                                    lambdaCOEF.set_val(i,myIndexVector[j]+i,1);
                                }
                            }
                            
                            // Partition coefficient matrix when linking with lambda variables
                            zCOEF.set_size(num_partition1+1, num_partition1);
                            
                            // Coefficients related to partition on v1
                            zCOEF.set_val(0, 0, 1);
                            for (int i=1; i<num_partition1; ++i) {
                                zCOEF.set_val(i,i-1,1);
                                zCOEF.set_val(i,i,1);
                            }
                            zCOEF.set_val(num_partition1, num_partition1-1, 1);
                            
                            // Partition assignment coefficients
                            zASGNCOEF.set_size(num_partition1);
                            
                            // Coefficients related to partition on v1
                            for (int i=0; i<num_partition1; ++i) {
                                zASGNCOEF.set_val(i,1);
                            }
                        }
                        // check if the first variable v2 involves a partition
                        else if (num_partition2>1){
                            // Lambda coefficient matrix when linking with partition variables
                            lambdaCOEF.set_size(num_partition2+1,((num_partition1+1)*(num_partition2+1)));
                            
                            // Coefficients related to partition on v2
                            for (int i=0; i<num_partition2+1; ++i) {
                                for (int j=0; j<num_partition1+1; ++j) {
                                    lambdaCOEF.set_val(i,j+i*(num_partition1+1),1);
                                }
                            }
                            
                            // Partition coefficient matrix when linking with lambda variables
                            zCOEF.set_size(num_partition2+1, num_partition2);
                            
                            // Coefficients related to partition on v2
                            zCOEF.set_val(0, 0, 1);
                            for (int i=1; i<num_partition2; ++i) {
                                zCOEF.set_val(i,i-1,1);
                                zCOEF.set_val(i,i, 1);
                            }
                            zCOEF.set_val(num_partition2, num_partition2-1, 1);
                            
                            // Partition assignment coefficients
                            zASGNCOEF.set_size(num_partition2);
                            
                            // Coefficients related to partition on v2
                            for (int i=0; i<num_partition2; ++i) {
                                zASGNCOEF.set_val(i,1);
                            }
                        }
                    }
                    else if (model_type == "Model_III"){
                        // check if both of the variables involve a partition
                        if ((num_partition1>1) && (num_partition2>1)) {
                            
                            // Lambda coefficient matrix when linking with partition variables
                            lambdaCOEF.set_size((num_partition1-2)*2+2+(num_partition2-2)*2+2,((num_partition1+1)*(num_partition2+1)));
                            
                            // Coefficients related to partition on v2
                            for (int i=0; i<num_partition1+1; i++) {
                                lambdaCOEF.set_val(0,i,1);
                            }
                            
                            for (int i=1; i<(num_partition2-2)*2+1; i=i+2) {
                                for (int j=(i-1)/2 + 2; j<num_partition2+1; ++j) {
                                    for(int k=0; k<num_partition1+1; ++k){
                                        lambdaCOEF.set_val(i,k+j*(num_partition1+1),1);
                                        lambdaCOEF.set_val(i+1,k+j*(num_partition1+1),-1);
                                    }
                                }
                            }
                            
                            for (int i=0; i<num_partition1+1; i++) {
                                lambdaCOEF.set_val((num_partition2-2)*2+1,i+(num_partition2)*(num_partition1+1),1);
                            }
                            
                            
                            // Collect the corresponding indices of lambda related to partition on v1
                            vector<int> myIndexVector(num_partition2+1); //initial index vector for lambdas bounding partition variables for x1
                            for (int i = 0; i<myIndexVector.size(); ++i) {
                                myIndexVector[i] = (i*(num_partition1+1));
                            }
                            
                            // Coefficients related to partition on v1
                            for (int i=0; i<num_partition2+1; i++) {
                                lambdaCOEF.set_val((num_partition2-2)*2+2,myIndexVector[i],1);
                            }
                            
                            for (int i=1; i<(num_partition1-2)*2+1; i=i+2) {
                                for (int j=(i-1)/2 + 2; j<num_partition1+1; ++j) {
                                    for(int k=0; k<num_partition2+1; ++k){
                                        lambdaCOEF.set_val(i+(num_partition2-2)*2+2,myIndexVector[k]+j,1);
                                        lambdaCOEF.set_val(i+(num_partition2-2)*2+3,myIndexVector[k]+j,-1);
                                    }
                                }
                            }
                            
                            for (int i=0; i<num_partition2+1; i++) {
                                lambdaCOEF.set_val((num_partition1-2)*2+2+(num_partition2-2)*2+1,myIndexVector[i]+num_partition1,1);
                            }
                            
                            // Partition coefficient matrix when linking with lambda variables
                            zCOEF.set_size((num_partition1-2)*2+2+(num_partition2-2)*2+2, num_partition1+num_partition2);
                            
                            // Coefficients related to partition on v2
                            zCOEF.set_val(0, num_partition1, 1);
                            for (int i=1; i<num_partition2; ++i) {
                                zCOEF.set_val(1, num_partition1+i, 1);
                            }
                            
                            for (int i=2; i<(num_partition2-2)*2+2; i=i+2) {
                                for (int j=i/2+1; j<num_partition2; ++j) {
                                    zCOEF.set_val(i, num_partition1+j, -1);
                                    zCOEF.set_val(i+1, num_partition1+j, 1);
                                }
                            }
                            
                            // Coefficients related to partition on v1
                            zCOEF.set_val((num_partition2-2)*2+2, 0, 1);
                            for (int i=1; i<num_partition1; ++i) {
                                zCOEF.set_val((num_partition2-2)*2+3, i, 1);
                            }
                            
                            for (int i=2; i<(num_partition1-2)*2+2; i=i+2) {
                                for (int j=i/2+1; j<num_partition1; ++j) {
                                    zCOEF.set_val(i+(num_partition2-2)*2+2, j, -1);
                                    zCOEF.set_val(i+(num_partition2-2)*2+3, j, 1);
                                }
                            }
                            
                            // Partition assignment coefficients
                            zASGNCOEF.set_size(2,num_partition1+num_partition2);
                            
                            // Coefficients related to partition on v1
                            for (int i=0; i<num_partition1; ++i) {
                                zASGNCOEF.set_val(0,i,1);
                            }
                            
                            // Coefficients related to partition on v2
                            for (int i=0; i<num_partition2; ++i) {
                                zASGNCOEF.set_val(1,i+num_partition1,1);
                            }
                        }
                        // check if the first variable v1 involves a partition
                        else if (num_partition1>1){
                            // Lambda coefficient matrix when linking with partition variables
                            lambdaCOEF.set_size((num_partition1-2)*2+2,((num_partition1+1)*(num_partition2+1)));
                            
                            // Collect the corresponding indices of lambda related to partition on v1
                            vector<int> myIndexVector(num_partition2+1); //initial index vector for lambdas bounding partition variables for x1
                            for (int i = 0; i<myIndexVector.size(); ++i) {
                                myIndexVector[i] = (i*(num_partition1+1));
                            }
                            
                            // Coefficients related to partition on v1
                            for (int i=0; i<num_partition2+1; i++) {
                                lambdaCOEF.set_val(0,myIndexVector[i],1);
                            }
                            
                            for (int i=1; i<(num_partition1-2)*2+1; i=i+2) {
                                for (int j=(i-1)/2 + 2; j<num_partition1+1; ++j) {
                                    for(int k=0; k<num_partition2+1; ++k){
                                        lambdaCOEF.set_val(i,myIndexVector[k]+j,1);
                                        lambdaCOEF.set_val(i+1,myIndexVector[k]+j,-1);
                                    }
                                }
                            }
                            
                            for (int i=0; i<num_partition2+1; i++) {
                                lambdaCOEF.set_val((num_partition1-2)*2+1,myIndexVector[i]+num_partition1,1);
                            }
                            
                            // Partition coefficient matrix when linking with lambda variables
                            zCOEF.set_size((num_partition1-2)*2+2, num_partition1);
                            
                            // Coefficients related to partition on v1
                            zCOEF.set_val(0, 0, 1);
                            for (int i=1; i<num_partition1; ++i) {
                                zCOEF.set_val(1, i, 1);
                            }
                            
                            for (int i=2; i<(num_partition1-2)*2+2; i=i+2) {
                                for (int j=i/2+1; j<num_partition1; ++j) {
                                    zCOEF.set_val(i, j, -1);
                                    zCOEF.set_val(i+1, j, 1);
                                }
                            }
                            
                            // Partition assignment coefficients
                            zASGNCOEF.set_size(num_partition1);
                            
                            // Coefficients related to partition on v1
                            for (int i=0; i<num_partition1; ++i) {
                                zASGNCOEF.set_val(i,1);
                            }
                        }
                        // check if the first variable v2 involves a partition
                        else if (num_partition2>1){
                            // Lambda coefficient matrix when linking with partition variables
                            lambdaCOEF.set_size((num_partition2-2)*2+2,((num_partition1+1)*(num_partition2+1)));
                            
                            // Coefficients related to partition on v2
                            for (int i=0; i<num_partition1+1; i++) {
                                lambdaCOEF.set_val(0,i,1);
                            }
                            
                            for (int i=1; i<(num_partition2-2)*2+1; i=i+2) {
                                for (int j=(i-1)/2 + 2; j<num_partition2+1; ++j) {
                                    for(int k=0; k<num_partition1+1; ++k){
                                        lambdaCOEF.set_val(i,k+j*(num_partition1+1),1);
                                        lambdaCOEF.set_val(i+1,k+j*(num_partition1+1),-1);
                                    }
                                }
                            }
                            
                            for (int i=0; i<num_partition1+1; i++) {
                                lambdaCOEF.set_val((num_partition2-2)*2+1,i+(num_partition2)*(num_partition1+1),1);
                            }
                            
                            
                            // Partition coefficient matrix when linking with lambda variables
                            zCOEF.set_size((num_partition2-2)*2+2, num_partition2);
                            
                            // Coefficients related to partition on v2
                            zCOEF.set_val(0, 0, 1);
                            for (int i=1; i<num_partition2; ++i) {
                                zCOEF.set_val(1, i, 1);
                            }
                            
                            for (int i=2; i<(num_partition2-2)*2+2; i=i+2) {
                                for (int j=i/2+1; j<num_partition2; ++j) {
                                    zCOEF.set_val(i, j, -1);
                                    zCOEF.set_val(i+1, j, 1);
                                }
                            }
                            
                            // Partition assignment coefficients
                            zASGNCOEF.set_size(num_partition2);
                            
                            // Coefficients related to partition on v2
                            for (int i=0; i<num_partition2; ++i) {
                                zASGNCOEF.set_val(i,1);
                            }
                        }
                    }
                    
                    
                    /** Variables */
                    // Convex combination variables
                    var<type> lambda(name+"_lambda",pos_);
                    lambda.in(R((num_partition1+1)*(num_partition2+1)));
                    add(lambda);
                    
                    /** Constraints */
                    // Representation of the bilinear term with convex combination
                    Constraint<type> BLNREP(name+"_Representation of bilinear term");
                    BLNREP = EP.tr()*lambda - vlift;
                    add(BLNREP == 0);
                    
                    // Representation of v1 with convex combination
                    Constraint<type> v1REP(name+"_Representation of v1");
                    for(int i=0;i<num_partition2+1;i++)
                    {
                        v1REP += bounds1.tr()*lambda.in(range(i*(num_partition1+1),i*(num_partition1+1)+num_partition1));
                    }
                    v1REP -= v1;
                    add(v1REP == 0);
                    
                    // Representation of v2 with convex combination
                    Constraint<type> v2REP(name+"_Representation of v2");
                    for(int i=0;i<num_partition2+1;i++)
                    {
                        v2REP += bounds2.eval(i)*sum(lambda.in(range(i*(num_partition1+1),i*(num_partition1+1)+num_partition1)));
                    }
                    v2REP -= v2;
                    add(v2REP == 0);
                    
                    // check if any of the variables involve a partition
                    if ((num_partition1>1) || (num_partition2>1)){
                        // Partition variables
                        var<int> z(name+"_z",0,1);
                        z.in(R((num_partition1)*(num_partition1>1)+(num_partition2)*(num_partition2>1)));
                        add(z);
                        
                        // Partition Assignment
                        Constraint<type> zASGN(name+"_Partition assignment");
                        zASGN = product(zASGNCOEF,z);
                        add(zASGN == 1);
                        
                        // Linking partition variables with lambda
                        Constraint<type> zLINKlambda(name+"_Linking the partition variables to lambda");
                        zLINKlambda = product(lambdaCOEF,lambda) - product(zCOEF,z);
                        add(zLINKlambda<= 0);
                    }
                    
                    // Summation over lambda is 1
                    Constraint<type> lambdaSUM(name+"_Summation over lambda");
                    lambdaSUM = sum(lambda);
                    add(lambdaSUM == 1);
                    
                    
                }
                else{
                    throw invalid_argument("Model type is invalid. The options are 'Model_II' and 'Model_III'!\n");
                }
            }
            else{ //this case means the two variables are same and the approximation is now on quadratic term
                if ((model_type == "Model_II") || (model_type == "Model_III")) {
                    
                    if (p1 == p2) //the partition vectors must be same since the variables are same
                    {
                        // Obtain the number of partitions
                        int num_partition1 = p1.size()-1;  // number of partition on v1 & v2 (they are same)
                        
                        /** Parameters */
                        // Bounds on variable v1 & v2
                        param<double> bounds1(name+"_bounds on variables v1&v2");
                        bounds1.in(R(num_partition1+1));
                        for (int i=0; i<num_partition1+1; ++i) {
                            bounds1.set_val(i,p1[i]);
                        }
                        
                        // Function values on the extreme points
                        param<double> EP(name+"_function values in the grid");
                        EP.in(R((num_partition1+1)));
                        for (int i=0; i<num_partition1+1; ++i) {
                            EP.set_val(i,p1[i]*p1[i]);
                        }
                        
                        // Lambda coefficient matrix when linking with partition variables
                        param<type> lambdaCOEF(name+"_lambda linking coefficients");
                        
                        // Partition coefficient matrix when linking with lambda variables
                        param<type> zCOEF(name+"_partition linking coefficients");
                        
                        // Partition assignment coefficients
                        param<type> zASGNCOEF(name+"_partition assignment coefficients");
                        
                        if (model_type == "Model_II"){
                            // check if the the variables involve a partition
                            if (num_partition1>1){
                                // Lambda coefficient matrix when linking with partition variables
                                lambdaCOEF.set_size(num_partition1+1,num_partition1+1);
                                
                                // Coefficients related to partition on variables v1 & v2
                                for (int i=0; i<num_partition1+1; ++i) {
                                    lambdaCOEF.set_val(i,i,1);
                                }
                                
                                // Partition coefficient matrix when linking with lambda variables
                                zCOEF.set_size(num_partition1+1, num_partition1);
                                
                                // Coefficients related to partition on variables v1 & v2
                                zCOEF.set_val(0, 0, 1);
                                for (int i=1; i<num_partition1; ++i) {
                                    zCOEF.set_val(i,i-1,1);
                                    zCOEF.set_val(i,i,1);
                                }
                                zCOEF.set_val(num_partition1, num_partition1-1, 1);
                                
                                // Partition assignment coefficients
                                zASGNCOEF.set_size(num_partition1);
                                
                                // Coefficients related to partition on variables v1 & v2
                                for (int i=0; i<num_partition1; ++i) {
                                    zASGNCOEF.set_val(i,1);
                                }
                            }
                        }
                        else if (model_type == "Model_III"){
                            // check if both of the variables involve a partition
                            if (num_partition1>1){
                                // Lambda coefficient matrix when linking with partition variables
                                lambdaCOEF.set_size((num_partition1-2)*2+2,num_partition1+1);
                                
                                // Coefficients related to partition on variables v1 & v2
                                lambdaCOEF.set_val(0,0,1);
                                
                                for (int i=1; i<(num_partition1-2)*2+1; i=i+2) {
                                    for (int j=(i-1)/2 + 2; j<num_partition1+1; ++j) {
                                        lambdaCOEF.set_val(i,j,1);
                                        lambdaCOEF.set_val(i+1,j,-1);
                                    }
                                }
                                lambdaCOEF.set_val((num_partition1-2)*2+1,num_partition1,1);
                                
                                // Partition coefficient matrix when linking with lambda variables
                                zCOEF.set_size((num_partition1-2)*2+2, num_partition1);
                                
                                // Coefficients related to partition on v1
                                zCOEF.set_val(0, 0, 1);
                                for (int i=1; i<num_partition1; ++i) {
                                    zCOEF.set_val(1, i, 1);
                                }
                                
                                for (int i=2; i<(num_partition1-2)*2+2; i=i+2) {
                                    for (int j=i/2+1; j<num_partition1; ++j) {
                                        zCOEF.set_val(i, j, -1);
                                        zCOEF.set_val(i+1, j, 1);
                                    }
                                }
                                
                                // Partition assignment coefficients
                                zASGNCOEF.set_size(num_partition1);
                                
                                // Coefficients related to partition on v1
                                for (int i=0; i<num_partition1; ++i) {
                                    zASGNCOEF.set_val(i,1);
                                }
                            }
                        }
                        
                        /** Variables */
                        // Convex combination variables
                        var<type> lambda(name+"_lambda",pos_);
                        lambda.in(R((num_partition1+1)));
                        add(lambda);
                        
                        /** Constraints */
                        // Representation of the bilinear term with secant
                        Constraint<type> BLNREP_UB(name+"_Representation of quadratic term (UB)");
                        BLNREP_UB = EP.tr()*lambda - vlift;
                        add(BLNREP_UB >= 0); /*using it as the upper bound to be valid*/
                        
                        Constraint<type> BLNREP_LB(name+"_Representation of quadratic term (LB)");
                        BLNREP_LB = pow(v1,2) - vlift;
                        add(BLNREP_LB <= 0); /*using it as the lower bound to be valid*/
                        
                        // Representation of v1=v2 with convex combination
                        Constraint<type> v1REP(name+"_Representation of v1");
                        v1REP == bounds1.tr()*lambda - v1;
                        add(v1REP == 0);
                        
                        // check if any of the variables involve a partition
                        if (num_partition1>1) {
                            // Partition variables
                            var<int> z(name+"_z",0,1);
                            z.in(R((num_partition1)));
                            add(z);
                            
                            // Partition Assignment
                            Constraint<type> zASGN(name+"_Partition assignment");
                            zASGN = product(zASGNCOEF,z);
                            add(zASGN == 1);
                            
                            // Linking partition variables with lambda
                            Constraint<type> zLINKlambda(name+"_Linking the partition variables to lambda");
                            zLINKlambda = product(lambdaCOEF,lambda) - product(zCOEF,z);
                            add(zLINKlambda<= 0);
                        }
                        
                        // Summation over lambda is 1
                        Constraint<type> lambdaSUM(name+"_Summation over lambda");
                        lambdaSUM = sum(lambda);
                        add(lambdaSUM == 1);
                    }
                    else
                    {
                        throw invalid_argument("Partition bounds must be same since the two varibles are same.\n");
                    }
                }
                else{
                    throw invalid_argument("Model type is invalid. The options are 'Model_II' and 'Model_III'!\n");
                }
            }
        }
        
        
        template<typename T1>
        void add_McCormick(std::string name, const var<T1>& vlift, const var<T1>& v1, const var<T1>& v2) {
            Constraint<type> MC1(name+"_McCormick1");
            param<T1> lb1_ = v1.get_lb(), lb2_ = v2.get_lb(), ub1_ = v1.get_ub(), ub2_= v2.get_ub();
            if(!lb1_.func_is_number()){
                lb1_ = v1.get_lb().in(*v1._indices);
            }
            if(!lb2_.func_is_number()){
                lb2_ = v2.get_lb().in(*v2._indices);
            }
            if(!ub1_.func_is_number()){
                ub1_ = v1.get_ub().in(*v1._indices);
            }
            if(!ub2_.func_is_number()){
                ub2_ = v2.get_ub().in(*v2._indices);
            }
            bool var_equal=false;
            if(v1._name==v2._name)
                var_equal=true;
            if(!var_equal)
            {
                Constraint<type> MC1(name+"_McCormick1");
                MC1 += vlift;
                MC1 -= lb1_*v2 + lb2_*v1 - lb1_*lb2_;
                MC1 >= 0;
                //                MC1._relaxed = true; /* MC1 is a relaxation of a non-convex constraint */
                add(MC1.in(*vlift._indices));
                //    MC1.print();
                Constraint<type> MC2(name+"_McCormick2");
                MC2 += vlift;
                MC2 -= ub1_*v2 + ub2_*v1 - ub1_*ub2_;
                MC2 >= 0;
                //                MC2._relaxed = true; /* MC2 is a relaxation of a non-convex constraint */
                add(MC2.in(*vlift._indices));
                
                //    //    MC2.print();
                Constraint<type> MC3(name+"_McCormick3");
                MC3 += vlift;
                MC3 -= lb1_*v2 + ub2_*v1 - lb1_*ub2_;
                MC3 <= 0;
                //                MC3._relaxed = true; /* MC3 is a relaxation of a non-convex constraint */
                add(MC3.in(*vlift._indices));
                //    //    MC3.print();
                Constraint<type> MC4(name+"_McCormick4");
                MC4 += vlift;
                MC4 -= ub1_*v2 + lb2_*v1 - ub1_*lb2_;
                MC4 <= 0;
                //                MC4._relaxed = true; /* MC4 is a relaxation of a non-convex constraint */
                add(MC4.in(*vlift._indices));
            }
            else {
                if (vlift._lift_ub){
                Constraint<type> MC4(name+"_Secant");
                MC4 += vlift;
                MC4 -= (ub1_+lb1_)*v1 - ub1_*lb1_;
                MC4 <= 0;
                add(MC4.in(*vlift._indices));
                }
                
                if (vlift._lift_lb){
                Constraint<type> MC5(name+"_McCormick_squared");
                if(var_equal)
                {
                    MC5 += vlift;
                    
                    MC5 -= v1*v1;
                    
                    MC5 >= 0;
                    MC5._relaxed = true; /* MC5 is a relaxation of a non-convex constraint */
                    add(MC5.in(*vlift._indices));
                }
                }
            }
            //    MC4.print();
        }
        
        //        template<typename T1>
        //        void add_McCormick(std::string name, const var<T1>& vlift, const var<T1>& v1, const var<T1>& v2) {
        //            Constraint<type> MC1(name+"_McCormick1");
        //            auto lb1 = v1.get_lb(v1.get_id_inst());
        //            auto lb2 = v2.get_lb(v2.get_id_inst());
        //            auto ub1 = v1.get_ub(v1.get_id_inst());
        //            auto ub2 = v2.get_ub(v2.get_id_inst());
        //            auto lb1_ = v1.get_lb().in(*v1._indices);
        //            auto lb2_ = v2.get_lb().in(*v2._indices);
        //            auto ub1_ = v1.get_ub().in(*v1._indices);
        //            auto ub2_ = v2.get_ub().in(*v2._indices);
        //            bool template_cstr = v1._dim[0]>1;
        //            MC1 += vlift;
        //            if(template_cstr){//Template constraint
        //                MC1 -= lb1_*v2 + lb2_*v1 - lb1_*lb2_;
        //            }
        //            else {
        //                MC1 -= lb1*v2 + lb2*v1 - lb1*lb2;
        //            }
        //            MC1 >= 0;
        //            add(MC1.in(*vlift._indices));
        //            //    MC1.print();
        //            Constraint<type> MC2(name+"_McCormick2");
        //            MC2 += vlift;
        //            if(template_cstr){//Template constraint
        //                MC2 -= ub1_*v2 + ub2_*v1 - ub1_*ub2_;
        //            }
        //            else {
        //                MC2 -= ub1*v2 + ub2*v1 - ub1*ub2;
        //            }
        //            MC2 >= 0;
        //            add(MC2.in(*vlift._indices));
        //            //    //    MC2.print();
        //            Constraint<type> MC3(name+"_McCormick3");
        //            MC3 += vlift;
        //            if(template_cstr){//Template constraint
        //                MC3 -= lb1_*v2 + ub2_*v1 - lb1_*ub2_;
        //            }
        //            else {
        //                MC3 -= lb1*v2 + ub2*v1 - lb1*ub2;
        //            }
        //            MC3 <= 0;
        //            add(MC3.in(*vlift._indices));
        //            //    //    MC3.print();
        //            Constraint<type> MC4(name+"_McCormick4");
        //            MC4 += vlift;
        //            if(template_cstr){//Template constraint
        //                MC4 -= ub1_*v2 + lb2_*v1 - ub1_*lb2_;
        //            }
        //            else{
        //                MC4 -= ub1*v2 + lb2*v1 - ub1*lb2;
        //            }
        //            MC4 <= 0;
        //            add(MC4.in(*vlift._indices));
        //            //    MC4.print();
        //        }
        
        
        /** Build the sequential McCormick relaxation for polynomial programs **/
        
        shared_ptr<Model<type>> build_McCormick(){
            //            replace_integers();
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
        
        type get_obj_val() const{
            _obj->allocate_mem();
            return _obj->eval();
        }
        
        void print_obj_val(int prec = 5) const{
            _obj->allocate_mem();
            cout << "Objective = " << to_string_with_precision(_obj->eval(),prec) << endl;
        }
        
        template<typename T=type,typename std::enable_if<is_arithmetic<T>::value>::type* = nullptr>
        void fill_in_cstr(const double* x , double* res, bool new_x){
            if (new_x) {
                set_x(x);
                compute_funcs();
            }
            compute_constrs<type>(_cons_vec, res, 0, _cons_vec.size());
            return;
            unsigned nr_threads = std::thread::hardware_concurrency();
            if(nr_threads>_cons_vec.size()){
                nr_threads=_cons_vec.size();
            }
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
                auto nb_ins = c->get_nb_inst();
                for (auto &v_p: c->get_vars()){
                    v = v_p.second.first.get();
                    vid = v->get_id();
                    id = 0;
                    for (size_t inst = 0; inst< nb_ins; inst++){
                        if (!*c->_all_lazy || !c->_lazy[inst]) {
                            cid = c->_id+id++;
                            if (v->_is_vector || v->is_matrix_indexed()) {
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
        template<typename T=type,typename std::enable_if<is_arithmetic<T>::value>::type* = nullptr>
        void fill_in_jac(const double* x , double* res, bool new_x){
            //    if (!_first_call_jac && (!new_x || _type==lin_m)) { /* No need to recompute jacobian for linear models */
            if (new_x) {
                set_x(x);
                compute_funcs();
            }
            for (size_t i = 0; i<_nnz_g; i++) {
                res[i] = 0;
            }
            if (!_first_call_jac && (_type==lin_m)) { /* No need to recompute jacobian for linear models */
                //            if (false) { /* No need to recompute jacobian for linear models */
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
                
                //                compute_jac(_cons_vec, res, 0, _cons_vec.size(), _first_call_jac, _jac_vals);
                //                return;
                unsigned nr_threads = std::thread::hardware_concurrency();
                if(nr_threads>_cons_vec.size()){
                    nr_threads=_cons_vec.size();
                }
                if (nr_threads==0) {
                    nr_threads = 1;
                }
#ifdef USE_MPI
                nr_threads = 1;
#endif
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
                auto nb_ins = c->get_nb_inst();
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
                        if (dfdx->func_is_number()) {
                            
                            
                            
                            
                            for (size_t inst = 0; inst< nb_ins; inst++){
                                if (!*c->_all_lazy || !c->_lazy[inst]) {
                                    cid = c->_id+id++;
                                    
                                    if (v->_is_vector || v->is_matrix_indexed()) {
                                        auto dim = v->get_dim(inst);
                                        for (size_t j = 0; j<dim; j++) {
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            res[idx] = dfdx->_val->at(0);
                                            _jac_vals[idx] = res[idx];
                                            idx++;
                                        }
                                    }
                                    else {
                                        res[idx] = dfdx->_val->at(0);
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
                                    if (v->_is_vector || v->is_matrix_indexed()) {
                                        auto dim = v->get_dim(inst);
                                        if (dfdx->is_matrix()) {
                                            for (size_t j = 0; j<dim; j++) {
                                                res[idx] += dfdx->eval(j,inst);
                                                _jac_vals[idx] = res[idx];
                                                idx++;
                                            }
                                        }
                                        else {
                                            for (size_t j = 0; j<dim; j++) {
                                                res[idx] += dfdx->eval(inst,j);
                                                //                                                res[idx] += dfdx->_val->at(j);//TODO check double indexed funcs
                                                _jac_vals[idx] = res[idx];
                                                DebugOff("jac_val["<< idx <<"] = " << _jac_vals[idx] << endl);
                                                idx++;
                                            }
                                        }
                                    }
                                    else {
                                        res[idx] = dfdx->_val->at(inst);
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
                auto nb_ins = c->get_nb_inst();
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
                vi = (pairs.second.begin())->second.first->get_var(vi_name);
                vj = (pairs.second.begin())->second.first->get_var(vj_name);
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
                    auto f = f_pair.second.first;
                    if (f->_is_constraint) {
                        c = static_pointer_cast<Constraint<type>>(f);
                    }
                    auto d2f = f_pair.second.second;
                    size_t nb_inst = f->get_nb_inst();
                    for (size_t inst = 0; inst<nb_inst; inst++) {
                        if (!(f->_is_constraint && *c->_all_lazy && c->_lazy[inst])) {
                            if(d2f->is_matrix_indexed()){
                                auto dim = d2f->get_dim(inst);
                                for (size_t j = 0; j<dim; j++) {
                                    idx_all++;
                                    iRow[_idx_it[idx]] = vid + vi->get_id_inst(inst,j);
                                    jCol[_idx_it[idx]] = vjd + vj->get_id_inst(inst,j);
                                    idx++;
                                }
                            }
                            else  if (d2f->is_matrix()) {
                                for (size_t i = 0; i < d2f->get_nb_inst(); i++) {
                                    for (size_t j = i; j < d2f->_dim[1]; j++) {
                                        idx_all++;
                                        iRow[_idx_it[idx]] = vid + vi->get_id_inst(i);
                                        jCol[_idx_it[idx]] = vjd + vj->get_id_inst(j);
                                        idx++;
                                        //                                f_idx++;
                                    }
                                }
                            }
                            else if(d2f->_is_vector){
                                //                    if (d2f->get_nb_inst() != d2f->get_nb_inst()) {
                                //                        throw invalid_argument("error");
                                //                    }
                                for (size_t j = 0; j < d2f->get_nb_inst(); j++) {
                                    idx_all++;
                                    iRow[_idx_it[idx]] = vid + vi->get_id_inst(j);
                                    jCol[_idx_it[idx]] = vjd + vj->get_id_inst(j);
                                    idx++;
                                    //                            f_idx++;
                                }
                            }
                            else {
                                idx_all++;
                                iRow[_idx_it[idx]] = vid + vi->get_id_inst(inst);
                                jCol[_idx_it[idx]] = vjd + vj->get_id_inst(inst);
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
            //            if (idx!=_nnz_h) {
            //                throw invalid_argument("idx!=_nnz_h");
            //            }
            _hess_vals.resize(idx_all);
        }
        
        template<typename T=type,typename std::enable_if<is_arithmetic<T>::value>::type* = nullptr>
        void fill_in_hess(const double* x , double obj_factor, const double* lambda, double* res, bool new_x){
            size_t idx = 0, idx_in = 0, c_inst = 0, idx_pair=0;
            shared_ptr<Constraint<type>> c;
            double hess = 0;
            if (new_x) {
                set_x(x);
                compute_funcs();
            }
            for (size_t i = 0; i<_nnz_h; i++) {
                res[i] = 0;
            }
            if (_first_call_hess) {
                for (auto &pairs: _hess_link) {
                    idx_pair = idx;
                    for (auto &f_pair:pairs.second) {
                        auto f = f_pair.second.first;
                        if (f->_is_constraint) {
                            c = static_pointer_cast<Constraint<type>>(f);
                        }
                        auto d2f = f_pair.second.second;
                        size_t nb_inst = f->get_nb_inst();
                        size_t id_inst = 0;
                        for (size_t inst = 0; inst<nb_inst; inst++) {
                            if (f->_is_constraint) {
                                if (!*c->_all_lazy || !c->_lazy[inst]) {
                                    if (c->is_nonlinear()) {
                                        c_inst = c->get_id_inst(id_inst++);
                                        if(d2f->is_matrix_indexed()){
                                            auto dim = d2f->get_dim(inst);
                                            for (size_t j = 0; j<dim; j++) {
                                                hess = d2f->eval(inst,j);
                                                _hess_vals[idx_in++] = hess;
                                                res[_idx_it[idx++]] += lambda[c->_id + c_inst] * hess;
                                            }
                                        }
                                        else if (d2f->is_matrix()) {
                                            for (size_t i = 0; i < d2f->get_nb_inst(); i++) {
                                                for (size_t j = i; j < d2f->_dim[1]; j++) {
                                                    hess = d2f->eval(i,j);
                                                    _hess_vals[idx_in++] = hess;
                                                    res[_idx_it[idx++]] += lambda[c->_id + c_inst] * hess;
                                                }
                                            }
                                        }
                                        else if(d2f->_is_vector){
                                            for (size_t j = 0; j < d2f->get_nb_inst(); j++) {
                                                if (d2f->func_is_number()) {
                                                    hess = d2f->_val->at(0);
                                                }
                                                else {
                                                    hess = d2f->_val->at(j);
                                                }
                                                _hess_vals[idx_in++] = hess;
                                                res[_idx_it[idx++]] += lambda[c->_id + c_inst] * hess;
                                            }
                                        }
                                        else {
                                            if (d2f->func_is_number()) {
                                                hess = d2f->_val->at(0);
                                            }
                                            else {
                                                hess = d2f->_val->at(inst);
                                            }
                                            _hess_vals[idx_in++] = hess;
                                            res[_idx_it[idx++]] += lambda[c->_id + c_inst] * hess;
                                        }
                                    }
                                    else {
                                        if(!d2f->is_constant()){
                                            d2f->_evaluated=false;
                                        }
                                        c_inst = c->get_id_inst(id_inst++);
                                        if(d2f->is_matrix_indexed()){
                                            auto dim = d2f->get_dim(inst);
                                            for (size_t j = 0; j<dim; j++) {
                                                hess = d2f->eval(inst,j);
                                                _hess_vals[idx_in++] = hess;
                                                res[_idx_it[idx++]] += lambda[c->_id + c_inst] * hess;
                                            }
                                        }
                                        else if (d2f->is_matrix()) {
                                            for (size_t i = 0; i < d2f->get_nb_inst(); i++) {
                                                for (size_t j = i; j < d2f->_dim[1]; j++) {
                                                    hess = d2f->eval(i,j);
                                                    _hess_vals[idx_in++] = hess;
                                                    res[_idx_it[idx++]] += lambda[c->_id + c_inst] * hess;
                                                }
                                            }
                                        }
                                        else if(d2f->_is_vector){
                                            for (size_t j = 0; j < d2f->get_nb_inst(); j++) {
                                                if (d2f->func_is_number()) {
                                                    hess = d2f->eval(0);
                                                }
                                                else {
                                                    hess = d2f->eval(j);
                                                }
                                                _hess_vals[idx_in++] = hess;
                                                res[_idx_it[idx++]] += lambda[c->_id + c_inst] * hess;
                                            }
                                        }
                                        else {
                                            if (d2f->func_is_number()) {
                                                hess = d2f->eval(0);
                                            }
                                            else {
                                                hess = d2f->eval(inst);
                                            }
                                            _hess_vals[idx_in++] = hess;
                                            res[_idx_it[idx++]] += lambda[c->_id + c_inst] * hess;
                                        }
                                    }
                                }
                            }
                            else {
                                if (d2f->is_matrix()) {
                                    for (size_t i = 0; i < d2f->get_nb_inst(); i++) {
                                        for (size_t j = i; j < d2f->_dim[1]; j++) {
                                            hess = d2f->eval(i,j);
                                            _hess_vals[idx_in++] = hess;
                                            res[_idx_it[idx++]] += obj_factor * hess;
                                        }
                                    }
                                }
                                else if(d2f->_is_vector){
                                    for (size_t j = 0; j < d2f->get_nb_inst(); j++) {
                                        if (d2f->func_is_number()) {
                                            hess = d2f->eval(0);
                                        }
                                        else {
                                            hess = d2f->eval(j);
                                        }
                                        _hess_vals[idx_in++] = hess;
                                        res[_idx_it[idx++]] += obj_factor * hess;
                                    }
                                }
                                else {
                                    hess = d2f->eval(0);
                                    _hess_vals[idx_in++] = hess;
                                    res[_idx_it[idx++]] += obj_factor * hess;
                                }
                            }
                        }
                    }
                }
                _first_call_hess = false;
                return;
            }
            if ((_type==lin_m || _type==quad_m)) { /* No need to recompute Hessian for quadratic models, used stored values */
                size_t id_inst = 0;
                for (auto &pairs: _hess_link) {
                    idx_pair = idx;
                    for (auto &f_pair:pairs.second) {
                        auto f = f_pair.second.first;
                        if (f->_is_constraint) {
                            c = static_pointer_cast<Constraint<type>>(f);
                        }
                        auto d2f = f_pair.second.second;
                        size_t nb_inst = f->get_nb_inst();
                        id_inst = 0;
                        for (size_t inst = 0; inst<nb_inst; inst++) {
                            if (f->_is_constraint) {
                                if (!*c->_all_lazy || !c->_lazy[inst]) {
                                    c_inst = c->get_id_inst(id_inst++);
                                    if(d2f->is_matrix_indexed()){
                                        auto dim = d2f->get_dim(inst);
                                        for (size_t j = 0; j<dim; j++) {
                                            res[_idx_it[idx++]] += lambda[c->_id + c_inst] * _hess_vals[idx_in++];
                                        }
                                    }
                                    else if (d2f->is_matrix()) {
                                        for (size_t i = 0; i < d2f->get_nb_inst(); i++) {
                                            for (size_t j = i; j < d2f->_dim[1]; j++) {
                                                res[_idx_it[idx++]] += lambda[c->_id + c_inst] * _hess_vals[idx_in++];
                                            }
                                        }
                                    }
                                    else if(d2f->_is_vector){
                                        for (size_t j = 0; j < d2f->get_nb_inst(); j++) {
                                            res[_idx_it[idx++]] += lambda[c->_id + c_inst] * _hess_vals[idx_in++];
                                        }
                                    }
                                    else {
                                        res[_idx_it[idx++]] += lambda[c->_id + c_inst] * _hess_vals[idx_in++];
                                    }
                                }
                            }
                            else {
                                if(d2f->_is_vector){
                                    for (size_t j = 0; j < d2f->get_nb_inst(); j++) {
                                        res[_idx_it[idx++]] += obj_factor * _hess_vals[idx_in++];
                                    }
                                }
                                else {
                                    res[_idx_it[idx++]] += obj_factor * _hess_vals[idx_in++];
                                }
                            }
                        }
                    }
                }
                return;
            }
            for (auto &pairs: _hess_link) {
                idx_pair = idx;
                for (auto &f_pair:pairs.second) {
                    auto f = f_pair.second.first;
                    if (f->_is_constraint) {
                        c = static_pointer_cast<Constraint<type>>(f);
                    }
                    auto d2f = f_pair.second.second;
                    if(!d2f->is_constant()){
                        d2f->_evaluated=false;
                    }
                    size_t nb_inst = f->get_nb_inst();
                    size_t id_inst = 0;
                    for (size_t inst = 0; inst<nb_inst; inst++) {
                        if (f->_is_constraint) {
                            if (!*c->_all_lazy || !c->_lazy[inst]) {
                                c_inst = c->get_id_inst(id_inst++);
                                if (c->is_quadratic()) {
                                    if(d2f->is_matrix_indexed()){
                                        auto dim = d2f->get_dim(inst);
                                        for (size_t j = 0; j<dim; j++) {
                                            res[_idx_it[idx++]] += lambda[c->_id + c_inst] * _hess_vals[idx_in++];
                                        }
                                    }
                                    else if (d2f->is_matrix()) {
                                        for (size_t i = 0; i < d2f->get_nb_inst(); i++) {
                                            for (size_t j = i; j < d2f->_dim[1]; j++) {
                                                res[_idx_it[idx++]] += lambda[c->_id + c_inst] * _hess_vals[idx_in++];
                                            }
                                        }
                                    }
                                    else if(d2f->_is_vector){
                                        for (size_t j = 0; j < d2f->get_nb_inst(); j++) {
                                            res[_idx_it[idx++]] += lambda[c->_id + c_inst] * _hess_vals[idx_in++];
                                        }
                                    }
                                    else {
                                        res[_idx_it[idx++]] += lambda[c->_id + c_inst] * _hess_vals[idx_in++];
                                    }
                                }
                                else if(d2f->is_matrix_indexed()){
                                    auto dim = d2f->get_dim(inst);
                                    for (size_t j = 0; j<dim; j++) {
                                        hess = d2f->eval(inst,j);
                                        _hess_vals[idx_in++] = hess;
                                        res[_idx_it[idx++]] += lambda[c->_id + c_inst] * hess;
                                    }
                                }
                                else if (d2f->is_matrix()) {
                                    for (size_t i = 0; i < d2f->get_nb_inst(); i++) {
                                        for (size_t j = i; j < d2f->_dim[1]; j++) {
                                            hess = d2f->eval(i,j);
                                            _hess_vals[idx_in++] = hess;
                                            res[_idx_it[idx++]] += lambda[c->_id + c_inst] * hess;
                                        }
                                    }
                                }
                                else if(d2f->_is_vector){
                                    for (size_t j = 0; j < d2f->get_nb_inst(); j++) {
                                        if (d2f->func_is_number()) {
                                            hess = d2f->_val->at(0);
                                        }
                                        else {
                                            hess = d2f->_val->at(j);
                                        }
                                        _hess_vals[idx_in++] = hess;
                                        res[_idx_it[idx++]] += lambda[c->_id + c_inst] * hess;
                                    }
                                }
                                else {
                                    if (d2f->func_is_number()) {
                                        hess = d2f->eval(0);
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
                                    res[_idx_it[idx++]] += lambda[c->_id + c_inst] * hess;
                                }
                            }
                        }
                        else {
                            if (_obj->is_quadratic()) {
                                if(d2f->_is_vector){
                                    for (size_t j = 0; j < d2f->get_nb_inst(); j++) {
                                        res[_idx_it[idx++]] += obj_factor * _hess_vals[idx_in++];
                                    }
                                }
                                else {
                                    res[_idx_it[idx++]] += obj_factor * _hess_vals[idx_in++];
                                }
                            }
                            else if(d2f->_is_vector){
                                for (size_t j = 0; j < d2f->get_nb_inst(); j++) {
                                    hess = d2f->eval(j);
                                    _hess_vals[idx_in++] = hess;
                                    res[_idx_it[idx++]] += obj_factor * hess;
                                }
                            }
                            else {
                                hess = d2f->eval(0);
                                _hess_vals[idx_in++] = hess;
                                res[_idx_it[idx++]] += obj_factor * hess;
                            }
                        }
                    }
                    
                    
                    
                }
                
            }
        }
        
        
        void initialize_zero(){
            for (auto &vp: _vars) {
                vp.second->initialize_zero();
            }
        }
        
        void initialize_uniform(){
            for (auto &vp: _vars) {
                vp.second->initialize_uniform();
            }
        }
        
        
        void fill_in_maps() {/*< Fill the _hess and _v_in_ maps to link variables with their constraints and compute the Jacobian & Hessian matrices */
            string vi_name, vj_name;
            param_* vi;
            param_* vj;
            _built = true;
            if (_obj->_new) {
                _obj->allocate_mem();
                _obj->compute_derivatives();
                for (auto &df_p:*_obj->get_dfdx()) {
                    auto df = static_pointer_cast<func<type>>(df_p.second);
                    DebugOff(df->to_str() << endl);
                    for (auto &df2_p:*df_p.second->get_dfdx()) {
                        if (df2_p.second->get_expr() || _type==nlin_m) {
                            df2_p.second = embed(df2_p.second);
                        }
                    }
                    if (df->get_expr() || _type==nlin_m) {
                        df_p.second = embed(df);
                    }
                    else {
                        embed(df);
                    }
                }
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
                                _hess_link[make_pair<>(vi_name,vj_name)][-1] = (make_pair<>(_obj,_obj->get_stored_derivative(vi->_name)->get_stored_derivative(vj->_name)));
                            }
                            else {
                                _hess_link[make_pair<>(vj_name,vi_name)][-1] = (make_pair<>(_obj,_obj->get_stored_derivative(vj->_name)->get_stored_derivative(vi->_name)));
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
                        for (auto &df2_p:*df_p.second->get_dfdx()) {
                            if (df2_p.second->get_expr() || _type==nlin_m) {
                                df2_p.second = embed(df2_p.second);
                            }
                        }
                        if (df->get_expr() || _type==nlin_m) {
                            df_p.second = embed(df);
                        }
                        else {
                            embed(df);
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
                                    _hess_link[make_pair<>(vi_name,vj_name)][c->_id] = (make_pair<>(c, c->get_stored_derivative(vi->_name)->get_stored_derivative(vj->_name)));
                                }
                                else {
                                    _hess_link[make_pair<>(vj_name,vi_name)][c->_id] = (make_pair<>(c, c->get_stored_derivative(vj->_name)->get_stored_derivative(vi->_name)));
                                }
                            }
                        }
                    }
                    DebugOff(c->to_str() << endl);
                }
            }
            //                print_nl_functions();
        }
        
        
        void fill_in_duals(double* lambda, double* z_L, double* z_U){
            for (auto &cp: _cons) {
                if(cp.second->_new){
                    continue;
                }
                size_t idx = 0;
                //        for (size_t inst = 0; inst < cp.second->get_nb_inst(); inst++) {
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
                if(vp.second->_new){
                    continue;
                }
                auto nb_inst = vp.second->get_dim();
                auto vid = vp.second->get_id();
                for (size_t inst = 0; inst < nb_inst; inst++) {
                    auto id_inst = vp.second->get_id_inst(inst);
                    z_L[vid + id_inst] = vp.second->_l_dual[inst];
                    z_U[vid + id_inst] = vp.second->_u_dual[inst];
                }
            }
            
        }
        /**
         Initialize x with model variables values
         @param[out] x values to initialize
         */
        template<typename T=type,typename std::enable_if<is_arithmetic<T>::value>::type* = nullptr>
        void fill_in_var_init(double* x) {
            for(auto& v_p: _vars){
                v_p.second->get_double_val(x);
            }
        }
        
        /**
         Initialize x with model variables values
         @param[out] x values to initialize
         */
        template<typename T=type,typename std::enable_if<is_arithmetic<T>::value>::type* = nullptr>
        void get_solution(vector<double>& x) const{
            for(auto& v_p: _vars){
                v_p.second->get_solution(x);
            }
        }
        
        void reset() {
            for(auto& c_p :_cons)
            {
                c_p.second->uneval();
            }
            for(auto &v_p: _vars)
            {
                v_p.second->reset_bounds();
            }
        }
        
        void reset_constrs() {
            for(auto& c_p :_cons)
            {
                c_p.second->uneval();
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
                        auto nb_ins = c->get_nb_inst();
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
                        auto nb_ins = c->get_nb_inst();
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
                        auto nb_ins = c->get_nb_inst();
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
                                //                f->_val->resize(f->get_nb_inst());
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
                                //                f->_val->resize(f->get_nb_inst());
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
                                //                f->_val->resize(f->get_nb_inst());
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
        
        
        shared_ptr<func<type>> embed(const shared_ptr<func<type>>& f, bool insert_in_map = true){/**<  Transfer all variables and parameters to the model. */
            f->allocate_mem();
            merge_vars(f);
            //            return f;
            DebugOff(f->to_str() << endl);
            for (auto &p_t: f->get_lterms()) {
                if (p_t.second._coef->is_function() && !p_t.second._coef->func_is_number()) {
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
                if (p_t.second._coef->is_function() && !p_t.second._coef->func_is_number()) {
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
                if (p_t.second._coef->is_function() && !p_t.second._coef->func_is_number()) {
                    auto cf = static_pointer_cast<func<type>>(p_t.second._coef);
                    auto exp = cf->get_expr();
                    if (exp) {
                        embed(*exp);
                    }
                }
            }
            if (f->get_cst()->is_function() && !f->get_cst()->func_is_number()) {
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
            if (!found_cpy && insert_in_map) {
                auto f_p = _nl_funcs_map.insert(make_pair<>(f->to_str(), f));
                if (f_p.second) {
                    _nl_funcs.push_back(f_p.first->second);
                    DebugOff(f->to_str() << endl);
                    f->allocate_mem();
                    return f;
                    //        f_p.first->second->_val = make_shared<vector<double>>();
                    //        f_p.first->second->_val->resize(f_p.first->second->get_nb_inst());
                }
                //        if (f->_new) {
                //            f_p.first->second = f;
                //            return f;
                //        }
                //                if (f->get_nb_inst() > f_p.first->second->get_nb_inst()) {
                //                    *f_p.first->second = *f;
                //                }
                
                //                else if (f->_dfdx->size()>0) {
                //                    *f_p.first->second = *f;
                
                //                }
                //        f_p.first->second->allocate_mem();
                return f_p.first->second;
            }
            return f;
        }
        
        
        
        void print_nl_functions() const{
            cout << "Number of atomic functions = " << _nl_funcs.size();
            cout << endl;
            for (auto& f: _nl_funcs){
                cout << f->_to_str;
                //                f->print();
                cout << endl;
            }
            cout << endl;
        }
        
        
        void print_solution(int prec=5) const{
            print_obj_val(prec);
            for (auto &v_pair:_vars) {
                auto v = v_pair.second;
                v->print_vals(prec);
            }
        }
        
        shared_ptr<param_> get_int_var(size_t idx){
            return _int_vars.at(idx);
        }
        
        void round_solution(){
            for (auto &v_pair:_vars) {
                if(v_pair.second->_is_relaxed){
                    v_pair.second->round_vals();
                    auto int_var = get_int_var(v_pair.second->get_id());
                    int_var->copy_vals(v_pair.second);
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
        
        void print(int prec = 10){
            auto size_header = print_properties();
            _obj->print(prec);
            cout << "s.t." << endl;
            for(auto& p: _cons){
                p.second->print(prec);
            }
            for(auto& v: _vars){
                v.second->print(prec);
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
                    _int_vars[v->get_id()] = v;
                    auto new_v = make_shared<var<double>>(v_p.second->_name);
                    new_v->shallow_copy(*v);
                    new_v->_is_relaxed = true;
                    new_v->copy_vals(v);
                    new_v->copy_bounds(v);
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
        
        
        
        /*** adding on/off option to a constraint within an interval ***/
        template<typename T=type,typename std::enable_if<is_arithmetic<type>::value>::type* = nullptr>
        void add_on_off_univariate(const Constraint<type>& c, var<bool>& on, int num_partns, int cur_partn){
            if (c.get_ftype() != lin_) {
                cerr << "Nonlinear constraint.\n";
                exit(-1);
            }
            
            param<type> M1 ("M1");
            param<type> M2 ("M2");
            
            size_t nb_ins = c.get_nb_inst();
            
            for (size_t inst = 0; inst<nb_ins; inst++)
            {
                for (auto &pair:*c._lterms) {
                    auto term = pair.second;
                    type coef_val = 0;
                    if (term._coef->is_function()) {
                        auto coef = static_pointer_cast<func<type>>(term._coef);
                        coef_val = coef->eval(inst);//this will give you the value of this instance
                    }
                    else if(term._coef->is_param()) {
                        auto coef = static_pointer_cast<param<type>>(term._coef);
                        coef_val = coef->eval(inst);//this will give you the value of this instance
                    }
                    else { /*means (term._coef->is_number())*/
                        auto coef = static_pointer_cast<constant<type>>(term._coef);
                        coef_val = coef->eval();
                    }
                    
                    auto LB = (term._p->get_double_lb(inst));
                    auto UB = (term._p->get_double_ub(inst));
                    
                    if (c.get_ctype() == eq) {
                        double LB_partn = (LB*(num_partns - cur_partn + 1) + UB*(cur_partn - 1))/num_partns;
                        double UB_partn = (LB*(num_partns - cur_partn) + UB*(cur_partn))/num_partns;
                        
                        if (coef_val < 0){
                            M1.add_val(coef_val * LB_partn);
                            M2.add_val(coef_val * UB_partn);
                        }
                        else {
                            M1.add_val(coef_val * UB_partn);
                            M2.add_val(coef_val * LB_partn);
                        }
                        
                    }
                    
                    else if (c.get_ctype() == leq) {
                        if (coef_val < 0){
                            type LB_partn = (LB*(num_partns - cur_partn + 1) + UB*(cur_partn - 1))/num_partns;
                            M1.add_val(coef_val * LB_partn);
                            
                        }
                        else {
                            type UB_partn = (LB*(num_partns - cur_partn) + UB*(cur_partn))/num_partns;
                            M1.add_val(coef_val * UB_partn);
                        }
                    }
                    
                    else {
                        if (coef_val < 0){
                            auto UB_partn = (LB*(num_partns - cur_partn) + UB*(cur_partn))/num_partns;
                            M2.add_val(coef_val * UB_partn);
                        }
                        else {
                            auto LB_partn = (LB*(num_partns - cur_partn + 1) + UB*(cur_partn - 1))/num_partns;
                            M2.add_val(coef_val * LB_partn);
                            
                        }
                    }
                }
            }
            
            
            if (c.get_ctype() == eq){
                Constraint<type> res1(c.get_name() + "_on/off");
                res1 = c - M1*(1-on);
                add_constraint(res1<=0);
                
                Constraint<type> res2(c.get_name() + "_on/off2");
                res2 = c - M2*(1-on);
                add_constraint(res2>=0);
            }
            
            if (c.get_ctype() == leq){
                Constraint<type> res1(c.get_name() + "_on/off");
                res1 = c - M1*(1-on);
                add_constraint(res1<=0);
            }
            
            if (c.get_ctype() == geq){
                Constraint<type> res2(c.get_name() + "_on/off2");
                res2 = c - M2*(1-on);
                add_constraint(res2>=0);
            }
        }
        
        
        //this function calculates on-off coefficients for representing the convex hull of the disjunctive union of constraints
        //INPUT: a constraint to be activated based on some binary variables, and a constraint type for making sure the coefficients are gathered correctly. Follows the formulation described by Hijazi et. al.
        void get_on_off_coefficients(Constraint<type>& c, const ConstraintType c_type){
            if (c.get_ftype() != lin_) {
                cerr << "Nonlinear constraint.\n";
                exit(-1);
            }
            //TODO : only consider the non-lazy instances to add on-off constraint
            
            /*allocate the coefficient vectors and the sum values to update them*/
            type M1sum_off, M2sum_off;
            type M1sum_on, M2sum_on;
            
            /*allocate the bound values*/
            type LB_off,UB_off;
            type LB_on,UB_on;
            
            size_t nb_ins = c.get_nb_inst();
            
            //go over all the instances and fill the coefficient values
            for (size_t inst = 0; inst<nb_ins; inst++)
            {
                M1sum_off = 0;
                M2sum_off = 0;
                M1sum_on = 0;
                M2sum_on = 0;
                
                c.eval_all();
                
                //cast a correct pointer for the constant part in the constraint
                if (!c.get_cst()->is_zero()) {
                    if (c.get_cst()->is_number()) {
                        auto f_cst = static_pointer_cast<constant<type>>(c.get_cst());
                        M1sum_on -= f_cst->eval();
                        M2sum_on -= f_cst->eval();
                    }
                    else if (c.get_cst()->is_param()) {
                        auto f_cst = static_pointer_cast<param<type>>(c.get_cst());
                        M1sum_on -= f_cst->eval(inst);
                        M2sum_on -= f_cst->eval(inst);
                    }
                    else {
                        auto f_cst = static_pointer_cast<func<type>>(c.get_cst());
                        M1sum_on -= f_cst->eval(inst);
                        M2sum_on -= f_cst->eval(inst);
                    }
                }
                //collect the instance index as a string
                auto partition_info = c._indices->_keys->at(inst);
                
                //go over the linear terms and cast pointers
                for (auto &pair:*c._lterms) {
                    auto term = pair.second;
                    
                    auto in_S = *term._p->_in; //collect that the lterm is in S or not
                    type coef_val = 0;
                    if (term._coef->is_function()) {
                        auto coef = static_pointer_cast<func<type>>(term._coef);
                        coef_val = coef->eval(inst);//this will give you the value of this instance
                    }
                    else if(term._coef->is_param()) {
                        auto coef = static_pointer_cast<param<type>>(term._coef);
                        coef_val = coef->eval(inst);//this will give you the value of this instance
                    }
                    else { /*means (term._coef->is_number())*/
                        auto coef = static_pointer_cast<constant<type>>(term._coef);
                        coef_val = coef->eval();
                    }
                    
                    auto inst_id = term._p->get_id_inst(inst);
                    auto num_partns = term._p->get_num_partns();
                    
                    /* update the coef_val as coef_val * sign */
                    if (!term._sign) coef_val = -coef_val;
                    
                    //set the global bounds
                    LB_off = (term._p->get_double_lb(inst_id));
                    UB_off = (term._p->get_double_ub(inst_id));
                    
                    auto lifted = term._p->get_lift();
                    if (lifted){ //if lifted to LB_on values should be the global bounds since the number of partitions is 1
                        LB_on = (term._p->get_double_lb(inst_id));
                        UB_on = (term._p->get_double_ub(inst_id));
                    }
                    else {
                        //collect the cur_partn number from the instance index (this is not the info _cur_partn stored in the variable, it is stored in the indices of the constraint)
                        auto name1 = term._p->get_name(true,true);
                        auto loc1 = partition_info.find(name1) + name1.length() +1 ;
                        auto loc2 = partition_info.find_first_of('}', loc1);
                        int cur_partn = stoi(partition_info.substr(loc1,loc2-loc1));
                        
                        if (cur_partn > num_partns) throw invalid_argument("Current partition is out of range (larger than the number of partitions)");
                        
                        LB_on = (LB_off*(num_partns - cur_partn + 1) + UB_off*(cur_partn - 1))/num_partns;
                        UB_on = (LB_off*(num_partns - cur_partn) + UB_off*(cur_partn))/num_partns;
                    }
                    
                    // update the coefficients based on the constraint type leq or geq
                    if (c_type == leq) {
                        if (coef_val < 0){
                            if(in_S){
                                M1sum_on -= coef_val * UB_on;
                            }
                            else {
                                M1sum_off += coef_val * LB_off;
                            }
                        }
                        else {
                            if(in_S){
                                M1sum_on -= coef_val * LB_on;
                            }
                            else {
                                M1sum_off += coef_val * UB_off;
                            }
                        }
                    }
                    // if the constraint is of type geq update the coefficients accordingly
                    else if (c_type == geq){
                        if (coef_val < 0){
                            if(in_S){
                                M2sum_on -= coef_val * LB_on;
                            }
                            else {
                                M2sum_off += coef_val * UB_off;
                            }
                        }
                        else {
                            if(in_S){
                                M2sum_on -= coef_val * UB_on;
                            }
                            else {
                                M2sum_off += coef_val * LB_off;
                            }
                        }
                    }
                    else { //throw an error if the constraint is of some other type (i.e. eq)
                        throw invalid_argument("Only leq and geq types are allowed. If you want to get coefficients for eq, use leq and geq consecutively.");
                    }
                }
                //set the coefficient values in the attribute belonging to constraint object
                if (c_type == leq){
                    c._offCoef.set_val(inst,M1sum_off);
                    c._onCoef.set_val(inst,M1sum_on);
                }
                else {
                    c._offCoef.set_val(inst,M2sum_off);
                    c._onCoef.set_val(inst,M2sum_on);
                }
            }
            
        }
        
        //this function calculates on-off coefficients for representing the convex hull of the disjunctive union of constraints
        //INPUT: a constraint to be activated based on some binary variables, assuming the constraint is of standard type "leq (<=) ". Follows the formulation described by Hijazi et. al.
        void get_on_off_coefficients_standard(Constraint<type>& c){
            if (c.get_ftype() != lin_) {
                cerr << "Nonlinear constraint.\n";
                exit(-1);
            }
            
            //TODO : only consider the non-lazy instances to add on-off constraint
            
            /*allocate the coefficient vectors and the sum values to update them*/
            type M1sum_off, M1sum_on;
            
            /*allocate the bound values*/
            type LB_off,UB_off;
            type LB_on,UB_on;
            
            size_t nb_ins = c.get_nb_inst();
            
            for (size_t inst = 0; inst<nb_ins; inst++)
            {
                string prev_name = "";
                M1sum_off = 0;
                M1sum_on = 0;
                
                c.eval_all();
                //collect the constant part in the constraint by casting the correct pointer
                if (!c.get_cst()->is_zero()) {
                    if (c.get_cst()->is_number()) {
                        auto f_cst = static_pointer_cast<constant<type>>(c.get_cst());
                        M1sum_on -= f_cst->eval();
                    }
                    else if (c.get_cst()->is_param()) {
                        auto f_cst = static_pointer_cast<param<type>>(c.get_cst());
                        M1sum_on -= f_cst->eval(inst);
                    }
                    else {
                        auto f_cst = static_pointer_cast<func<type>>(c.get_cst());
                        M1sum_on -= f_cst->eval(inst);
                    }
                }
                //collect the instance index as a string
                auto partition_info = c._indices->_keys->at(inst);
                
                //go over the lterms and calculate the proper coefficients
                for (auto &pair:*c._lterms) {
                    auto in_S = *pair.second._p->_in; //collect that the lterm is in S or not
                    type coef_val = 0;
                    if (pair.second._coef->is_function()) {
                        auto coef = static_pointer_cast<func<type>>(pair.second._coef);
                        coef_val = coef->eval(inst);//this will give you the value of this instance
                    }
                    else if(pair.second._coef->is_param()) {
                        auto coef = static_pointer_cast<param<type>>(pair.second._coef);
                        coef_val = coef->eval(inst);//this will give you the value of this instance
                    }
                    else { /*means (term._coef->is_number())*/
                        auto coef = static_pointer_cast<constant<type>>(pair.second._coef);
                        coef_val = coef->eval();
                    }
                    
                    auto inst_id = pair.second._p->get_id_inst(inst);
                    auto num_partns = pair.second._p->get_num_partns();
                    auto in_SOC_partn = pair.second._p->get_in_SOC_partn();
                    
                    /* update the coef_val as coef_val * sign */
                    if (!pair.second._sign) coef_val = -coef_val;
                    
                    //set the global bounds
                    LB_off = (pair.second._p->get_double_lb(inst_id));
                    UB_off = (pair.second._p->get_double_ub(inst_id));
                    
                    auto lifted = pair.second._p->get_lift();
                    if (lifted || (num_partns == 1) || in_SOC_partn){ //if lifted to LB_on values should be the global bounds since the number of partitions is 1, similarly if number of partitions is 1
                        LB_on = (pair.second._p->get_double_lb(inst_id));
                        UB_on = (pair.second._p->get_double_ub(inst_id));
                    }
                    else {
                        //collect the cur_partn number from the instance index (this is not the info _cur_partn stored in the variable, it is stored in the indices of the constraint)
                        auto name1 = pair.second._p->get_name(true,true);
                        int cur_partn;
                        //this uses the variable names to access the correct number of partitions. The convention is if there are two variables with same base name, the names will be ordered alphabetically, and accessed similarly as well
                        if(prev_name == name1){
                            auto loc1 = partition_info.rfind(name1) + name1.length() +1 ;
                            auto loc2 = partition_info.find_first_of('}', loc1);
                            cur_partn = stoi(partition_info.substr(loc1,loc2-loc1));
                            
                        }
                        else{
                            auto loc1 = partition_info.find(name1) + name1.length() +1 ;
                            auto loc2 = partition_info.find_first_of('}', loc1);
                            cur_partn = stoi(partition_info.substr(loc1,loc2-loc1));
                        }
                        
                        //set the prev_name for variables having the same name
                        prev_name = name1;
                        
                        if (cur_partn > num_partns) throw invalid_argument("Current partition is out of range (larger than the number of partitions)");
                        
                        LB_on = (LB_off*(num_partns - cur_partn + 1) + UB_off*(cur_partn - 1))/num_partns;
                        UB_on = (LB_off*(num_partns - cur_partn) + UB_off*(cur_partn))/num_partns;
                    }
                    if (coef_val < 0){
                        if(in_S){
                            M1sum_on -= coef_val * UB_on;
                        }
                        else {
                            M1sum_off += coef_val * LB_off;
                        }
                    }
                    else {
                        if(in_S){
                            M1sum_on -= coef_val * LB_on;
                        }
                        else {
                            M1sum_off += coef_val * UB_off;
                        }
                    }
                }
                //set the coefficient values in the attribute belonging to constraint object
                c._offCoef.set_val(inst,M1sum_off);
                c._onCoef.set_val(inst,M1sum_on);
            }
            
        }
        
        //This function adds on-off version of a given linear constraint and the binary variables to activate. The boolean option handles all the facet definining inequalities of the convex hull (if true), else it only adds the Big_M version of the constraint
        //INPUT: linear constraint to be activated, corresponding binary variables to form a disjunctive union, big_M version of the constraint or the whole convex hull
        //OUTPUT: disjunctive union of the constraints provided by "c" and linked by "on"
        void add_on_off_multivariate_refined(Constraint<type>& c, const var<int>& on, bool big_M = false){
            if (c.get_ftype() != lin_) {
                cerr << "Nonlinear constraint.\n";
                exit(-1);
            }
            
            // TODO: consider only the not lazy ones in on-off
            // TODO: maybe do this somehow in the constructor
            c._onCoef.in(range(1,c.get_nb_inst()));
            c._offCoef.in(range(1,c.get_nb_inst()));
            
            //use bitset vector to represent S efficiently
            //currently this function can handle up to 64 linear variables in the constraint, and throws and error otherwise
            auto n_terms = c._lterms->size();
            if (n_terms > 64){
                throw invalid_argument("Currently we can not handle more than 64 linear terms in an on/off constraint. Please do not use partitioning or decrease the number of linear terms.");
            }
            std::bitset<64> S;
            //decide on the subset selection limit
            int num_subset;
            
            //if there are less than 3 variables appearing in the linear constraint, the function will only add the big-M version, since that with the bounds on the variables defines the convex hull
            if ((n_terms <= 2) || (big_M)) num_subset = 1;
            else num_subset = std::pow(2,n_terms) -1;
            
            for (int i = 0 ; i< num_subset ; ++i) { //should be num_subset
                S = i;
                int j = 0;
                
                shared_ptr<pair<type,type>> term_range;
                func<type> LHS; //to handle the left hand side of the constaint
                for (auto &lt:*c._lterms) { //set the _in values and create LHS
                    *(lt.second._p->_in) = S[j];
                    if (!S[j]){ //only if not in S
                        //cast the correct pointer to the coef part of the linear term and add that to the LHS depending on their _in values
                        auto coef = lt.second._coef->copy();
                        if (coef->is_function()) {
                            auto f_cst = *((func<type>*)(coef.get()));
                            auto var_range = make_shared<pair<type,type>>(c.get_range(lt.second._p));
                            term_range = get_product_range(f_cst._range,var_range);
                            LHS.insert(lt.second._sign, f_cst, *lt.second._p);
                        }
                        else if(coef->is_param()) {
                            auto p_cst = *((param<type>*)(coef.get()));
                            auto var_range = make_shared<pair<type,type>>(c.get_range(lt.second._p));
                            term_range = get_product_range(p_cst._range,var_range);
                            LHS.insert(lt.second._sign, p_cst, *lt.second._p);
                        }
                        else if(coef->is_number()) {
                            auto p_cst = *((constant<type>*)(coef.get()));
                            auto var_range = make_shared<pair<type,type>>(c.get_range(lt.second._p));
                            term_range = get_product_range(make_shared<pair<type,type>>(p_cst.eval(),p_cst.eval()),var_range);
                            LHS.insert(lt.second._sign, p_cst, *lt.second._p);
                        }
                        //update the ranges of the function appropriately
                        if(lt.second._sign){
                            LHS._range = get_plus_range(LHS._range, term_range);
                        }
                        else {
                            LHS._range = get_minus_range(LHS._range, term_range);
                        }
                    }
                    j++;
                }
                
                // all the cases are standardized into the leq form, then for each case, we obtain the _onCoef and _offCoef for the constraint and add it to the model with the proper LHS value for each subset S
                if (c.get_ctype() == eq) {
                    get_on_off_coefficients_standard(c);
                    auto offCoef1 = c._offCoef.deep_copy();
                    auto onCoef1 = c._onCoef.deep_copy();
                    Constraint<type> res1(c.get_name() + "_" + to_string(i) + "_on/off");
                    res1 = LHS - offCoef1*(1-on) - onCoef1*on;
                    add_constraint(res1.in(*c._indices)<=0);
                    Constraint<type> n_c(c);
                    n_c *= -1;
                    get_on_off_coefficients_standard(n_c);
                    auto offCoef2 = n_c._offCoef.deep_copy();
                    auto onCoef2 = n_c._onCoef.deep_copy();
                    Constraint<type> res2(c.get_name() +  "_" + to_string(i) + "_on/off2");
                    res2 = -1 * LHS - offCoef2*(1-on) - onCoef2*on;
                    add_constraint(res2.in(*c._indices)<=0);
                    
                }
                
                else if (c.get_ctype() == leq) {
                    get_on_off_coefficients_standard(c);
                    auto offCoef = c._offCoef.deep_copy();
                    auto onCoef = c._onCoef.deep_copy();
                    Constraint<type> res1(c.get_name() +  "_" + to_string(i) + "_on/off");
                    res1 = LHS - offCoef*(1-on) - onCoef*on;
                    add_constraint(res1.in(*c._indices)<=0);
                }
                
                else { //if c.get_ctype() == geq
                    Constraint<type> n_c(c);
                    n_c *= -1;
                    n_c._ctype = leq;
                    get_on_off_coefficients_standard(n_c);
                    auto offCoef = n_c._offCoef.deep_copy();
                    auto onCoef = n_c._onCoef.deep_copy();
                    Constraint<type> res2(c.get_name() +  "_" + to_string(i) + "_on/off2");
                    res2 = -1 * LHS - offCoef*(1-on) - onCoef*on;
                    add_constraint(res2.in(*c._indices)<=0);
                }
            }
        }
        
        //this function adds the on/off version of the McCormick relaxations for bilinear terms, and on/off version of the secants for the quadratic terms
        //INPUT: base name for the constraints to be added, lifted variable, correspoinding variables v1 and v2, and the partitioning indicator variable on (binary)
        //OUTPUT: the on/off version of the McCormick relaxation (for every partition and the disjunctive union of that) for bilinears and the secant relaxation (for every partition and the disjunctive union of that) for quadratic variables
        template<typename T1>
        void add_on_off_McCormick_refined(std::string name, var<T1>&& vlift, var<T1>&& v1, var<T1>&& v2, const var<int>& on) {
            
            if (!v1.is_bounded_below() || !v2.is_bounded_below() || !vlift.is_bounded_below() || !v1.is_bounded_above() || !v2.is_bounded_above() || !vlift.is_bounded_above()){
                throw invalid_argument("Variables have to be bounded. Please set bounds for all!");}
            
            if (!(vlift._lift)){
                throw invalid_argument("You forgot to set _lift to true for the lifted variable.");
            }
            
            //collect the number of partitions for each variable
            int num_partns1 = v1._num_partns;
            int num_partns2 = v2._num_partns;
            
            //collect the base name of each variable
            auto name1 = v1.get_name(true,true);
            auto name2 = v2.get_name(true,true);
            
            //if the variables are different apply McCormick relaxations for each partition pair
            if(v1._name!=v2._name)
            {
                //a dimension check to prevent any errors
                if (on.get_dim() != v1.get_dim() * num_partns1 * num_partns2){
                    throw invalid_argument("Number of on variables are not conforming with the given number of partitions");}
                
                //create the proper index sets for partitions
                indices partns1("partns1");
                for (int i = 0; i < num_partns1 ; ++i)
                {
                    partns1.add(name1+ "{" +to_string(i+1) + "}");
                }
                
                indices partns2("partns2");
                for (int i = 0; i < num_partns2 ; ++i)
                {
                    partns2.add(name2+ "{"+ to_string(i+1) + "}");
                }
                
                indices partns("partns");
                //                partns = indices(range(1,num_partns1),range(1,num_partns2));
                partns = indices(partns1,partns2);
                
                //combine the index sets for instances
                auto var_indices = combine(*v1._indices,*v2._indices);
                auto inst_partition = indices(var_indices,partns);
                
                // Create the parameters for the McCormick constraints (there are 4 differenc MC constraints)
                param<type> V1par_MC1("V1par_MC1");
                V1par_MC1.in(inst_partition);
                param<type> V2par_MC1("V2par_MC1");
                V2par_MC1.in(inst_partition);
                param<type> Cpar_MC1("Cpar_MC1");
                Cpar_MC1.in(inst_partition);
                
                param<type> V1par_MC2("V1par_MC2");
                V1par_MC2.in(inst_partition);
                param<type> V2par_MC2("V2par_MC2");
                V2par_MC2.in(inst_partition);
                param<type> Cpar_MC2("Cpar_MC2");
                Cpar_MC2.in(inst_partition);
                
                param<type> V1par_MC3("V1par_MC3");
                V1par_MC3.in(inst_partition);
                param<type> V2par_MC3("V2par_MC3");
                V2par_MC3.in(inst_partition);
                param<type> Cpar_MC3("Cpar_MC3");
                Cpar_MC3.in(inst_partition);
                
                param<type> V1par_MC4("V1par_MC4");
                V1par_MC4.in(inst_partition);
                param<type> V2par_MC4("V2par_MC4");
                V2par_MC4.in(inst_partition);
                param<type> Cpar_MC4("Cpar_MC4");
                Cpar_MC4.in(inst_partition);
                
                //create the parameters for the bounds on the variables (for on and off cases)
                param<type> v1_on_LB("v1_on_LB");
                v1_on_LB.in(inst_partition);
                param<type> v1_off_LB("v1_off_LB");
                v1_off_LB.in(inst_partition);
                
                param<type> v1_on_UB("v1_on_UB");
                v1_on_UB.in(inst_partition);
                param<type> v1_off_UB("v1_off_UB");
                v1_off_UB.in(inst_partition);
                
                param<type> v2_on_LB("v2_on_LB");
                v2_on_LB.in(inst_partition);
                param<type> v2_off_LB("v2_off_LB");
                v2_off_LB.in(inst_partition);
                
                param<type> v2_on_UB("v2_on_UB");
                v2_on_UB.in(inst_partition);
                param<type> v2_off_UB("v2_off_UB");
                v2_off_UB.in(inst_partition);
                
                //fill the parameters for the variable on-off bounds and the corresponding function values
                size_t nb_ins = v1.get_nb_inst();
                
                auto v1_global_lb = v1.get_lb();
                auto v1_global_ub = v1.get_ub();
                auto v2_global_lb = v2.get_lb();
                auto v2_global_ub = v2.get_ub();
                auto increment1 = (v1_global_ub - v1_global_lb)/num_partns1;
                auto increment2 = (v2_global_ub - v2_global_lb)/num_partns2;
                
                //for each partition pair and instance of the original constraint, set the values of all of the parameters
                for (int i=0 ; i<num_partns1; ++i) {
                    auto LB_partn1 = v1_global_lb + increment1*i;
                    auto UB_partn1 = LB_partn1 + increment1;
                    LB_partn1.eval_all();
                    UB_partn1.eval_all();
                    for (int j=0 ; j<num_partns2; ++j) {
                        auto LB_partn2 = v2_global_lb + increment2*j;
                        auto UB_partn2 = LB_partn2 + increment2;
                        LB_partn2.eval_all();
                        UB_partn2.eval_all();
                        for (size_t inst = 0; inst< nb_ins; inst++){
                            //collect the proper indexing combining all the information available
                            auto cur_var_idx = var_indices._keys->at(inst);
                            string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"}";
                            
                            v1_off_LB.set_val(cur_idx,v1_global_lb.eval(inst));
                            v1_off_UB.set_val(cur_idx,v1_global_ub.eval(inst));
                            v1_on_LB.set_val(cur_idx,LB_partn1.eval(inst));
                            v1_on_UB.set_val(cur_idx,UB_partn1.eval(inst));
                            
                            v2_off_LB.set_val(cur_idx,v2_global_lb.eval(inst));
                            v2_off_UB.set_val(cur_idx,v2_global_ub.eval(inst));
                            v2_on_LB.set_val(cur_idx,LB_partn2.eval(inst));
                            v2_on_UB.set_val(cur_idx,UB_partn2.eval(inst));
                            
                            V2par_MC1.set_val(cur_idx,LB_partn1.eval(inst));
                            V1par_MC1.set_val(cur_idx,LB_partn2.eval(inst));
                            Cpar_MC1.set_val(cur_idx,LB_partn1.eval(inst)*LB_partn2.eval(inst));
                            
                            V2par_MC2.set_val(cur_idx,UB_partn1.eval(inst));
                            V1par_MC2.set_val(cur_idx,UB_partn2.eval(inst));
                            Cpar_MC2.set_val(cur_idx,UB_partn1.eval(inst)*UB_partn2.eval(inst));
                            
                            Cpar_MC3.set_val(cur_idx,LB_partn1.eval(inst)*UB_partn2.eval(inst));
                            
                            Cpar_MC4.set_val(cur_idx,UB_partn1.eval(inst)*LB_partn2.eval(inst));
                            
                        }
                    }
                    
                }
                //add the McCormick constraints (all four of them) by using the previously calculated parameters, original variables and the indicator variables for the partitions
                auto nb_entries_v1 = v1._indices->get_nb_entries();
                Constraint<type> MC1(name+"_McCormick1");
                MC1 = vlift.from_ith(0,inst_partition) - V1par_MC1*v1.from_ith(0,inst_partition) - V2par_MC1*v2.from_ith(nb_entries_v1,inst_partition) + Cpar_MC1;
                MC1.in(inst_partition) >= 0;
                add_on_off_multivariate_refined(MC1, on);
                
                Constraint<type> MC2(name+"_McCormick2");
                MC2 = vlift.from_ith(0,inst_partition) - V1par_MC2*v1.from_ith(0,inst_partition) - V2par_MC2*v2.from_ith(nb_entries_v1,inst_partition) + Cpar_MC2;
                MC2.in(inst_partition) >= 0;
                add_on_off_multivariate_refined(MC2, on);
                
                Constraint<type> MC3(name+"_McCormick3");
                MC3 = vlift.from_ith(0,inst_partition) - V1par_MC2*v1.from_ith(0,inst_partition) - V2par_MC1*v2.from_ith(nb_entries_v1,inst_partition) + Cpar_MC3;
                MC3.in(inst_partition) <= 0;
                add_on_off_multivariate_refined(MC3, on);
                
                Constraint<type> MC4(name+"_McCormick4");
                MC4 = vlift.from_ith(0,inst_partition) - V1par_MC1*v1.from_ith(0,inst_partition) - V2par_MC2*v2.from_ith(nb_entries_v1,inst_partition) + Cpar_MC4;
                MC4.in(inst_partition) <= 0;
                add_on_off_multivariate_refined(MC4, on);
                
                //add the on-off bound constraints on original the variables
                Constraint<type> v1_on_off_LB(name+"_v1_on_off_LB");
                v1_on_off_LB = v1.from_ith(0,inst_partition) - on*v1_on_LB - (1-on)*v1_off_LB;
                v1_on_off_LB.in(inst_partition) >= 0;
                add(v1_on_off_LB);
                
                Constraint<type> v1_on_off_UB(name+"_v1_on_off_UB");
                v1_on_off_UB = v1.from_ith(0,inst_partition) - on*v1_on_UB - (1-on)*v1_off_UB;
                v1_on_off_UB.in(inst_partition) <= 0;
                add(v1_on_off_UB);
                
                Constraint<type> v2_on_off_LB(name+"_v2_on_off_LB");
                v2_on_off_LB = v2.from_ith(nb_entries_v1,inst_partition) - on*v2_on_LB - (1-on)*v2_off_LB;
                v2_on_off_LB.in(inst_partition) >= 0;
                add(v2_on_off_LB);
                
                Constraint<type> v2_on_off_UB(name+"_v2_on_off_UB");
                v2_on_off_UB = v2.from_ith(nb_entries_v1,inst_partition) - on*v2_on_UB - (1-on)*v2_off_UB;
                v2_on_off_UB.in(inst_partition) <= 0;
                add(v2_on_off_UB);
            }
            
            
            else{ //if the variables are same lift the quadratic via on/off secant
                if (num_partns1 != num_partns2) throw invalid_argument("Partition numbers must be same since the two varibles are same.\n");
                if (on.get_dim() != v1.get_dim() * num_partns1){
                    throw invalid_argument("Number of on variables are not conforming with the given number of partitions");}
                
                //create the index set for the partitions and combine them with the instance indices
                indices partns("partns");
                for (int i = 0; i < num_partns1 ; ++i)
                {
                    partns.add(name1+"{"+to_string(i+1) + "}");
                }
                auto var_indices = *v1._indices;
                auto inst_partition = indices(var_indices,partns);
                
                //Create the parameters for on/off constraints (including secant parameters and the on/off bounds on the variables)
                param<type> Vpar("Vpar");
                Vpar.in(inst_partition);
                param<type> Cpar("Cpar");
                Cpar.in(inst_partition);
                
                param<type> v1_on_LB("v1_on_LB");
                v1_on_LB.in(inst_partition);
                param<type> v1_off_LB("v1_off_LB");
                v1_off_LB.in(inst_partition);
                
                param<type> v1_on_UB("v1_on_UB");
                v1_on_UB.in(inst_partition);
                param<type> v1_off_UB("v1_off_UB");
                v1_off_UB.in(inst_partition);
                
                //fill the parameters for the variable on-off bounds and the corresponding function values
                size_t nb_ins = v1.get_nb_inst();
                
                //collect the bounds on the variables for efficient on/off bound calculation (considering uniform partitioning in the domain)
                auto v1_global_lb = v1.get_lb();
                auto v1_global_ub = v1.get_ub();
                auto increment = (v1_global_ub - v1_global_lb)/num_partns1;
                
                //for every partition and instance, fill the parameter values accordingly
                for (int i=0 ; i<num_partns1; ++i) {
                    auto LB_partn = v1_global_lb + increment*i;
                    auto UB_partn = LB_partn + increment;
                    LB_partn.eval_all();
                    UB_partn.eval_all();
                    for (size_t inst = 0; inst< nb_ins; inst++){
                        auto cur_var_id = v1.get_id_inst(inst);
                        auto cur_var_idx = var_indices._keys->at(cur_var_id);
                        string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"}";
                        v1_off_LB.set_val(cur_idx,v1_global_lb.eval(inst));
                        v1_off_UB.set_val(cur_idx,v1_global_ub.eval(inst));
                        v1_on_LB.set_val(cur_idx,LB_partn.eval(inst));
                        v1_on_UB.set_val(cur_idx,UB_partn.eval(inst));
                        Vpar.set_val(cur_idx,LB_partn.eval(inst)+UB_partn.eval(inst));
                        Cpar.set_val(cur_idx,LB_partn.eval(inst)*UB_partn.eval(inst));
                        
                    }
                }
                //add the secant constraint if the quadratic term appears on the LHS
                if (vlift._lift_ub){
                Constraint<type> MC_secant(name+"_secant");
                MC_secant = vlift.from_ith(0,inst_partition) - Vpar*v1.from_ith(0,inst_partition) + Cpar;
                MC_secant.in(inst_partition) <= 0;
                add_on_off_multivariate_refined(MC_secant, on);
                }
                
                //add the global lb constraint if the quadratic term appears on the RHS
                if (vlift._lift_lb){
                Constraint<type> MC_squared(name+"_McCormick_squared");
                MC_squared += vlift;
                MC_squared -= v1*v1;
                MC_squared >= 0;
                MC_squared._relaxed = true; /* MC_squared is a relaxation of a non-convex constraint */
                add(MC_squared.in(*vlift._indices));
                }
                
                //add the on-off bound constraints on the variables
                Constraint<type> v1_on_off_LB(name+"_v1_on_off_LB");
                v1_on_off_LB = v1.from_ith(0,inst_partition) - on*v1_on_LB - (1-on)*v1_off_LB;
                v1_on_off_LB.in(inst_partition) >= 0;
                add(v1_on_off_LB);
                
                Constraint<type> v1_on_off_UB(name+"_v1_on_off_UB");
                v1_on_off_UB = v1.from_ith(0,inst_partition) - on*v1_on_UB - (1-on)*v1_off_UB;
                v1_on_off_UB.in(inst_partition) <= 0;
                add(v1_on_off_UB);
                
            }
        }
        
        /** Outer approximation of model. Throws exception if model has nonlinear equality constraints
         @param[in] nb_discr: number of OA cuts per nonlinear constraint
         @return Model with OA cuts. OA cuts are added to the model (for all func instances) in an uniform grid (nb_discr)
         **/
        shared_ptr<Model<>> buildOA(int nb_discr)
        {
            auto OA=make_shared<Model<>>(_name+"-OA Model");
            for (auto &it: _vars)
            {
                auto v = it.second;
                if(!OA->has_var(*v)){
                    OA->add_var(v);
                }
            }
            auto obj=*_obj;
            if(_objt==minimize){
                OA->min(obj);
            }
            else {
                OA->max(obj);
            }
            string cname;
            for (auto &con: _cons_vec)
            {
                if(!con->is_linear()) {
                    if(con->_ctype==eq)
                    {
                        DebugOn("Exception: Equality constraint is not currently supported"<<endl);
                        DebugOn("Throw exception" <<endl);
                        
                    }
                    else
                    {
                        OA->add_outer_app_uniform(nb_discr, *con);
                    }
                }
                else
                {
                    OA->add(*con);
                }
            }
            return OA;
        }
        
        
        /** Outer approximation of model. Throws exception if model has nonlinear equality constraints
         @param[in] nb_discr:
         @param[in] nb_perturb:
         @param[in] xinterior:
         @return Model with OA cuts. OA cuts are added to the model (for all func instances) in an uniform grid (nb_discr) and at the solution and by perturbing the solution
         **/
        shared_ptr<Model<>> buildOA(int nb_discr, int nb_perturb, bool interior, vector<double> xinterior)
        {
            auto OA= buildOA(nb_discr);
            vector<double> xsolution(_nb_vars);
            set_solution(xsolution);
            OA->add_outer_app_active(*this, nb_perturb, interior, xinterior);
            set_solution(xsolution);
            return OA;
        }
        
        /** Outer approximation of active (nonlinear) constraints of the model
         @param[in] nonlin: original nonlinear model at whose solution (at the active point) OA cuts are added:
         @param[in] nb_perturb:
         @return void. OA cuts are added to the model that calls the function (for all func instances) at the solution and by perturbing the solution
         Assumption: nonlinear constraint to be perturbed does not have any vector variables
         **/
        void add_outer_app_active(Model<> nonlin, int nb_perturb, bool interior, vector<double> xinterior)
        {
            const double active_tol=1e-6, perturb_dist=1e-3;
            vector<double> xsolution(_nb_vars);
            vector<double> xactive, xcurrent;
            double fk;
            bool outer;
            int counter=0;
            size_t posv;
            for (auto &con: nonlin._cons_vec)
            {
                if(!con->is_linear() && con->_name!="obj_UB") {
                    
                    con->uneval();
                    for(auto i=0;i<con->get_nb_inst();i++){
                        if(std::abs(con->eval(i))<=active_tol || (con->is_convex() && !con->is_rotated_soc() && !con->check_soc())){
                            Constraint<> OA_sol("OA_cuts_solution "+con->_name+to_string(i));
                            OA_sol=con->get_outer_app_insti(i);
                            if(con->_ctype==leq) {
                                add(OA_sol<=0);
                            }
                            else {
                                add(OA_sol>=0);
                            }
                        }
                    }
                }
                else if(con->_name=="obj_UB"){
                    add(*con);
                }
            }
            if(interior)
            {
                get_solution(xsolution);
                for (auto &con: nonlin._cons_vec)
                {
                    if(!con->is_linear()) {
                        for(auto i=0;i<con->get_nb_inst();i++){
                            con->uneval();
                            if(std::abs(con->eval(i))<=active_tol && (!con->is_convex() || con->is_rotated_soc() || con->check_soc())){
                                set_solution(xinterior);
                                xinterior=con->get_x(i);
                                set_solution(xsolution);
                                xcurrent=con->get_x(i);
                                for(auto j=1;j<=nb_perturb;j++)
                                {
                                    counter=0;
                                    for(auto &it: *con->_vars)
                                    {
                                        auto v = it.second.first;
                                        if(v->_is_vector)
                                        {
                                            DebugOn("Exception: Vector variables are not currently supported"<<endl);
                                            DebugOn("Throw exception" <<endl);
                                            break;
                                        }
                                        else
                                        {
                                            outer=false;
                                            for(auto k=-1;k<=1;k+=2)
                                            {
                                                posv=v->get_id_inst(i);
                                                v->set_double_val(posv, xcurrent[counter]*(1+k*j*perturb_dist));
                                                con->uneval();
                                                fk=con->eval(i);
                                                if((fk>active_tol && con->_ctype==leq) || (fk<(active_tol*(-1)) && con->_ctype==geq)){
                                                    outer=true;
                                                    break;
                                                }
                                            }
                                            if(outer)
                                            {
                                                auto res_search=con->linesearchbinary(xinterior, i, con->_ctype);
                                                if(res_search){
                                                    Constraint<> OA_active("OA_active "+con->_name+to_string(i)+to_string(j)+v->_name);
                                                    OA_active=con->get_outer_app_insti(i);
                                                    if(con->_ctype==leq) {
                                                        add(OA_active<=0);
                                                    }
                                                    else {
                                                        add(OA_active>=0);
                                                    }
                                                    
                                                }
                                                
                                            }
                                        }
                                        con->set_x(i, xcurrent);
                                        counter++;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        //this function partitions a given SOC constraint to given number of uniform regions and construct hyperplanes in order to satisfy the SOC constraint at equality with an inner approximation as a convex relaxation (which is originally a non-convex constraint)
        //INPUT: an SOC constraint to satisfy at equality, number of desired partitions, another number of partitions if the original SOC constraint involves more than 3 variables (we need to seperate that into different SOC constraints), use_lambda option for using the lambda formulation to activate the hyperplanes where the default is the on/off formulation
        template<typename T=type>
        void SOC_partition(Constraint<type>& c, int num_SOC_partitions1 = 10, int num_SOC_partitions2 = 10, bool use_lambda = false) { //currently the function asssumes there are only qterms in the original SOC constraint
            // TODO: incorporate linear terms in the constraint, or convert those into the standard format as well!
            
            auto is_rotated_SOC = c.check_rotated_soc(); //check the constraint is either an SOC or a rotated SOC
            auto is_SOC = c.check_soc();
            
            if (!is_rotated_SOC && !is_SOC) // if the constraint is not an SOC, throw an error
            {
                throw invalid_argument("SOC partition can only be applied to second-order cones! \n");
            }
            if (!c._lterms->empty()) //currently, there is no support for linear terms in the SOC constraint
            {
                throw invalid_argument("Current SOC partition version can only handle quadratic & bilinear terms! \n");
            }
            if (c.get_ctype() == eq){ //the original constraint should be either leq or geq since we need to also have the convex relaxation of the original constraint as the upper bound on the cone itself
                throw invalid_argument("Please provide the constraint in the convex format (either <= or >=), and make sure you add the convex constraint to the model as well!");
            }
            //convert the constraint into standard format if the constraint is given as a >= constraint
            Constraint<type> n_c(c);
            if (c.get_ctype() == geq){
                n_c *= -1;
                n_c._ctype = leq;
            }
            
            //create the variables for bound calculations on the auxiliary variable for the case when there are more than 3 terms in the given SOC constraint
            var<> v1;
            var<> v2;
            
            if (is_SOC) //if the original constraint is an SOC
            {
                unsigned num_qterms = n_c._qterms->size();
                if (num_qterms >= 4){ //we need to split the constraint into two
                    if (num_qterms > 4){ //current version do not support more than 4 quadratic terms in the constraint
                        //TODO: create a generic scheme to further divide the constraint into smaller constraints and add them recursively
                        throw invalid_argument("Current SOC partition version can only up to four quadratic terms! \n");
                    }
                    //create the two constraints for the divided version of the original constraint
                    Constraint<type> SOC_1(n_c._name + "_SOC_1"); //to split the constraint (first half)
                    Constraint<type> SOC_2(n_c._name + "_SOC_2"); //to split the constraint (second half)
                    
                    //in order to decide which constraint +t^2 should be added (the auxiliary variable to connect the separated SOC constraints)
                    bool first = false;
                    auto aux_idx = *n_c._indices;
                    
                    //flag for assignment when assigning the variable pointers
                    bool first_occupied = false;
                    
                    //go over the qterms
                    unsigned counter = 0; // counter is for dividing the constraint into two with 2 quadratic terms each
                    for (auto &qt_pair: *n_c._qterms) {
                        if (!qt_pair.second._p->first->is_double()) {
                            throw invalid_argument("Current hyperplanes only support double type variables!");
                        }
                        if (!qt_pair.second._coef->is_number()) { /*means coef is not a number*/
                            throw invalid_argument("Current hyperplanes only support constant coefficients for the variables");
                        }
                        auto sign = qt_pair.second._sign;
                        auto coef = static_pointer_cast<constant<type>>(qt_pair.second._coef);
                        if (counter < 2) {
                            SOC_1.insert(qt_pair.second);
                            if (sign ^ (coef->is_positive())) { //update the first to true if the variable appearing has a negative sign
                                first = true;
                            }
                        }
                        else SOC_2.insert(qt_pair.second);
                        ++counter;
                    }
                    
                    //TODO: remove the ranges of extra terms in SOC_1 and SOC_2
                    //set the ranges of the functions properly
                    SOC_1._range = n_c._range;
                    SOC_2._range = n_c._range;
                    
                    if (first){ //if the t^2 is added to the first constraint, colect the other variables' pointers to obtain proper bounds for t
                        for (auto &qt_pair: *SOC_2._qterms) {
                            if (!first_occupied){
                                v1 = *static_pointer_cast<var<double>>(qt_pair.second._p->first);
                                first_occupied = true;
                            }
                            else v2 = *static_pointer_cast<var<double>>(qt_pair.second._p->first);
                        }
                    }
                    else{
                        for (auto &qt_pair: *SOC_1._qterms) {
                            if (!first_occupied){
                                v1 = *static_pointer_cast<var<double>>(qt_pair.second._p->first);
                                first_occupied = true;
                            }
                            else v2 = *static_pointer_cast<var<double>>(qt_pair.second._p->first);
                        }
                    }

                    //update the upper bound of t properly based on the bounds on the other variables appering the reduced SOC constraint
                    param<type> ub("ub");
                    ub.in(aux_idx);
                    
                    //go over all the instances and calculate the upper bound for t
                    for (int i=0; i<aux_idx.size(); i++) {
                        //calculate all the possibilities and assign the worst case
                        size_t id1;
                        size_t id2;
                        
                        if(v1._indices->_ids == nullptr){
                            id1 = i;
                        }
                        else id1 = v1._indices->_ids->at(0).at(i);
                        auto key1 = v1._indices->_keys->at(id1);
                        
                        if(v2._indices->_ids == nullptr){
                            id2 = i;
                        }
                        else id2 = v2._indices->_ids->at(0).at(i);
                        auto key2 = v2._indices->_keys->at(id2);
                        
                        //the largest possible values for the squared variables
                        auto prod_b1 = std::max(v1.get_lb(key1)*v1.get_lb(key1), v1.get_ub(key1)*v1.get_ub(key1));
                        auto prod_b2 = std::max(v2.get_lb(key2)*v2.get_lb(key2), v2.get_ub(key2)*v2.get_ub(key2));
                        
                        size_t id3;
                        if(aux_idx._ids == nullptr){
                            id3 = i;
                        }
                        else id3 = aux_idx._ids->at(0).at(i);
                        auto key3 = aux_idx._keys->at(id3);
                        
                        ub.set_val(key3, std::sqrt(prod_b1+prod_b2)); //since it is a cone constraint t can be at most this
                    }
                    
                    var<type> t("t_" + n_c._name, 0, ub); //create the auxilary variable
                    add(t.in(aux_idx));
                    //TODO: consider the case where there can be multiple negative terms!
                    if (first) { //add constraints accordingly
                        SOC_1 += pow(t.in(aux_idx),2);
                        add(SOC_1.in(aux_idx) <= 0);
                        
                        SOC_2 -= pow(t.in(aux_idx),2);
                        add(SOC_2.in(aux_idx) <= 0);
                    }
                    else{
                        SOC_1 -= pow(t.in(aux_idx),2);
                        add(SOC_1.in(aux_idx) <= 0);
                        
                        SOC_2 += pow(t.in(aux_idx),2);
                        add(SOC_2.in(aux_idx) <= 0);
                    }
                    
                    //call the hyperplane function to generate the disjunctive union of hyperplanes, if use_lambda is specified use the lambda formulation, otherwise use the on_off formulation
                    if(use_lambda){
                        add_lambda_SOC_hyperplanes_symmetric(SOC_1, num_SOC_partitions1);
                        add_lambda_SOC_hyperplanes_symmetric(SOC_2, num_SOC_partitions2);
                    }
                    else{
                    add_on_off_SOC_hyperplanes(SOC_1, num_SOC_partitions1);
                    add_on_off_SOC_hyperplanes(SOC_2, num_SOC_partitions2);
                    }
                }
                else { //if the number of qterms is less than 4, we can directly generate the hyperplanes
                    if(use_lambda){
                        add_lambda_SOC_hyperplanes_symmetric(n_c, num_SOC_partitions1); //this uses the first int as the number of partitions not the total one
                    }
                    else{
                        add_on_off_SOC_hyperplanes(n_c, num_SOC_partitions1); //this uses the first int as the number of partitions not the total one
                    }
                }
            }
            else { //means the original constraint is a rotated SOC (then we should take care of the bilinear terms)
                
                unsigned num_qterms = 0;
                unsigned num_blnterms = 0;
                for (auto &qt_pair: *n_c._qterms) { //go over the terms and count the standardized number of quadratic terms (which means a bilinear term will have 2 terms since it can be standardized with two quadratic terms)
                    if (qt_pair.second._p->first!=qt_pair.second._p->second) ++num_blnterms;
                    else ++num_qterms;
                }
                if (num_blnterms > 1){ //current procedure does not allow two bilinear terms
                    throw invalid_argument("Current SOC partition version can only allow one bilinear term! \n");
                }
                if (num_qterms + 2*num_blnterms >= 4){ //we need to split the constraint into two
                    
                    if (num_qterms + 2*num_blnterms > 4){ //current procedure does not allow more than 4 standard quadratic terms
                        throw invalid_argument("Current SOC partition version can only up to four quadratic terms! \n");
                    }
                    //create the individual SOC constraint for the separated halves
                    Constraint<type> SOC_1(n_c._name + "_SOC_1"); //to split the constraint (first half)
                    Constraint<type> SOC_2(n_c._name + "_SOC_2"); //to split the constraint (second half)
                    
                    auto aux_idx = *n_c._indices; /*** use this or n_c._indices->_keys->at(inst)? ***/
                    
                    //flag for assignment when assigning the variable pointers
                    bool first_occupied = false;
                    
                    //go over the qterms
                    for (auto &qt_pair: *n_c._qterms) {
                        if (!qt_pair.second._p->first->is_double()) {
                            throw invalid_argument("Current hyperplanes only support double type variables!");
                        }
                        if (qt_pair.second._p->first!=qt_pair.second._p->second){ //if the term is bilinear insert that into the first constraint by default
                            SOC_1.insert(qt_pair.second);
                        }
                        else { //if the term is quadratic, automatically insert that into the second part, and get the variables for calculating the bound of the auxiliary variable t
                            if (!first_occupied) {
                                v1 = *static_pointer_cast<var<double>>(qt_pair.second._p->first);
                                first_occupied = true;
                            }
                            else {
                                v2 = *static_pointer_cast<var<double>>(qt_pair.second._p->first);
                            }
                            SOC_2.insert(qt_pair.second);
                        }
                    }
                    //update the range of the functions appropriately
                    //TODO: remove the ranges of extra terms in SOC_1 and SOC_2
                    SOC_1._range = n_c._range;
                    SOC_2._range = n_c._range;
                    
                    
                    //create the parameter for the upper bound of t
                    param<type> ub("ub");
                    ub.in(aux_idx);
                    
                    //go over all the instances and calculate the best possible bound for the aux variable
                    for (int i=0; i<aux_idx.size(); i++) {
                        //calculate all the possibilities and assign the worst case
                        size_t id1;
                        size_t id2;
                        
                        if(v1._indices->_ids == nullptr){
                            id1 = i;
                        }
                        else id1 = v1._indices->_ids->at(0).at(i);
                        auto key1 = v1._indices->_keys->at(id1);
                        
                        if(v2._indices->_ids == nullptr){
                            id2 = i;
                        }
                        else id2 = v2._indices->_ids->at(0).at(i);
                        auto key2 = v2._indices->_keys->at(id2);

                        //these terms are the largest possible values that the squared variables can take
                        auto prod_b1 = std::max(v1.get_lb(key1)*v1.get_lb(key1), v1.get_ub(key1)*v1.get_ub(key1));
                        auto prod_b2 = std::max(v2.get_lb(key2)*v2.get_lb(key2), v2.get_ub(key2)*v2.get_ub(key2));
                        
                        size_t id3;
                        if(aux_idx._ids == nullptr){
                            id3 = i;
                        }
                        else id3 = aux_idx._ids->at(0).at(i);
                        auto key3 = aux_idx._keys->at(id3);
                        
                        ub.set_val(key3, std::sqrt(prod_b1+prod_b2)); //since the constraint is SOC, this is the maximum value t can take
                    }
                    
                    var<type> t("t_" + n_c._name, 0, ub); //create the auxilary variable
                    add(t.in(aux_idx));
                    //TODO: consider the case where there can be multiple negative terms!
                    //add constraints accordingly with including the auxiliary term
                    SOC_1 += pow(t.in(aux_idx),2);
                    add(SOC_1.in(aux_idx) <= 0);
                    
                    SOC_2 -= pow(t.in(aux_idx),2);
                    add(SOC_2.in(aux_idx) <= 0);
                    
                    //call the hyperplane function to generate the disjunctive union of hyperplanes (either with using lambda formulation or on/off formulation)
                    if(use_lambda){
                        add_lambda_SOC_hyperplanes_symmetric(SOC_1,num_SOC_partitions1);
                        add_lambda_SOC_hyperplanes_symmetric(SOC_2,num_SOC_partitions2);
                    }
                    else{
                    add_on_off_SOC_hyperplanes(SOC_1,num_SOC_partitions1);
                    add_on_off_SOC_hyperplanes(SOC_2,num_SOC_partitions2);
                    }
                }
                
                
                else{ //if the total number of standard terms are less than 4, we can directly generate the hyperplanes by standardizing the bilinear term
                    if(use_lambda){
                        add_lambda_SOC_hyperplanes_symmetric(n_c,num_SOC_partitions1);
                    }
                    else{
                        add_on_off_SOC_hyperplanes(n_c,num_SOC_partitions1);
                    }
                }
            }
        }

        /** Returns an interior point of a model
         @param[in] nonlin: model for which interior point with respect to nonlinear constraints is computed
         Assuming model has no nonlinear equality constraints
         **/
        Model<> build_model_interior()
        {
            Model<> Interior(_name+"Interior");
            vector<double> xinterior(_nb_vars);
            
            
            for (auto &it: _vars)
            {
                auto v = it.second;
                if(!Interior.has_var(*v)){
                    Interior.add_var(v);
                }
            }
            var<> eta_int("eta_int", -1, 0);
            
            Interior.add(eta_int.in(range(0,0)));
            auto obj=eta_int;
            
            Interior.min(obj);
            
            
            
            for (auto &con: _cons_vec)
            {
                if(!con->is_linear()) {
                    Constraint<> Inter_con(*con);
                    if(con->_ctype==leq)
                    {
                        Inter_con -= eta_int;
                        Interior.add(Inter_con<=0);
                    }
                    else
                    {
                        Inter_con += eta_int;
                        Interior.add(Inter_con>=0);
                    }
                }
                else
                {
                    Interior.add(*con);
                }
            }
            
            return Interior;
        }
        
        
        
        /** Discretizes Constraint con and adds OA cuts to the model that calls it. Discretization of squared constraint only currently implemented
         @param[in] nb_discr:
         @param[in] con:
         @return void. OA cuts are added to the model that calls the function (for all func instances)
         **/
        void add_outer_app_uniform(int nb_discr, Constraint<> con)
        {
            
            func<> res;
            double lb,ub;
            size_t posv;
            if(con.is_quadratic() && con._lterms->size()==1 && con._qterms->size()==1 && con._qterms->begin()->second._p->first==con._qterms->begin()->second._p->second) //This if is specific to constraints of the form ay- x^2 or x^2-ay
            {
                
                auto x=con._qterms->begin()->second._p->first;
                
                for(auto d=0;d<nb_discr;d++)
                {
                    for(auto i=0;i<con.get_nb_inst();i++)
                    {
                        posv=x->get_id_inst(i);
                        lb=x->get_double_lb(posv);
                        ub=x->get_double_ub(posv);
                        x->set_double_val(posv, lb+d*(ub-lb)/nb_discr);
                    }
                    Constraint<> OA_uniform("OA_cuts_uniform "+con._name+to_string(d));
                    OA_uniform=con.get_outer_app_squared();
                    if(con._ctype==leq) {
                        add(OA_uniform<=0);
                    }
                    else {
                        add(OA_uniform>=0);
                    }
                    
                }
            }/*TODO Else (discretization for general constraint)*/
        }

        // INPUT: an SOC type constraint, and total number of binary variables
        // OUTPUT: disjunctive union of hyperplanes as an inner approximation to the SOC, where the disjunctive union is made by lambda formulation
        // IMPORTANT NOTE: we also utilize the symmetric nature of the formulation to make the formulation more efficient
        void add_lambda_SOC_hyperplanes_symmetric(Constraint<type>& c, int num_SOC_partitions){ //currently this is not fully correct
            //TODO: scale the coefficients properly
            
            DebugOn("SOC_hyperplane function!" << endl);
            c.print();
            
            auto is_rotated_SOC = c.check_rotated_soc(); //collect the information about the cone
            auto is_SOC = c.check_soc();
            
            //create hyperplane indices
            indices hyper_idx("hyper_idx");
            for (int i=0; i<num_SOC_partitions; ++i) {
                hyper_idx.add(to_string(i+1));
            }
            //get the combined index set
            auto inst_hyper = indices(*c._indices,hyper_idx);
            
            //create lambda indices
            indices lambda_idx("lambda_idx");
            for (int i=0; i<num_SOC_partitions+2; ++i) {
                lambda_idx.add(to_string(i+1));
            }
            //get the combined index set
            auto inst_lambda = indices(*c._indices,lambda_idx);
            
            //create on_link_lambda indices
            indices on_link_lambda_idx("on_link_lambda_idx");
            for (int i=0; i<num_SOC_partitions+1; ++i) {
                on_link_lambda_idx.add(to_string(i+1));
            }
            
            Constraint<type> v1_rep_pos(c._name + "_v1_rep_pos"); //create the v1_rep constraint (positive side)
            Constraint<type> v1_rep_neg(c._name + "_v1_rep_neg"); //create the v1_rep constraint (negative side)
            Constraint<type> v2_rep_pos(c._name + "_v2_rep_pos"); //create the v2_rep constraint (positive side)
            Constraint<type> v2_rep_neg(c._name + "_v2_rep_neg"); //create the v2_rep constraint (negative side)
            Constraint<type> v3_rep(c._name + "_v3_rep"); //create the v3_rep constraint
            var<int> on(c._name + "_binary",0,1); //create the partition variable
            Constraint<type> onSum(c._name + "_binarySum"); //create the partition assignment constraint
            var<> lambda(c._name + "_lambda",0,1); //create the lambda variable
            Constraint<type> lambdaSum(c._name + "_lambdaSum"); //create the lambda constraint
            Constraint<type> on_link_lambda(c._name + "_on_link_lambda"); // lambda linking with the partition variables
            var<int> v1_sign(c._name + "_v1_sign",0,1); //create the sign indicator for v1 (1 if negative 0 positive)
            var<int> v2_sign(c._name + "_v2_sign",0,1); //create the sign indicator for v2 (1 if negative 0 positive)
            
            if (is_SOC){ //this will follow the standard creation of the hyperplane
                
                //create the variables
                var<> lhs_first_var;
                var<> lhs_second_var;
                var<> rhs_var;
                
                //create the scaling factors for the variables
                double lhs_first_scale;
                double lhs_second_scale;
                double rhs_scale;
                
                //flag for assignment
                bool first_occupied = false;
                
                //go over the quadratic terms and collect the variables and their multipliers (scales)
                for (auto &qt_pair: *c._qterms) {
                    if (!qt_pair.second._p->first->is_double()) {
                        throw invalid_argument("Current hyperplanes only support double type variables!");
                    }
                    if (!qt_pair.second._coef->is_number()) { /*means coef is not a number*/
                        throw invalid_argument("Current hyperplanes only support constant coefficients for the variables");
                    }
                    auto sign = qt_pair.second._sign;
                    auto coef = static_pointer_cast<constant<type>>(qt_pair.second._coef);
                    if (sign ^ (coef->is_negative())) {
                        if (!first_occupied){
                            lhs_first_var = *static_pointer_cast<var<double>>(qt_pair.second._p->first);
                            if(coef->is_positive()) lhs_first_scale = std::sqrt(coef->eval());
                            else lhs_first_scale = std::sqrt((-1)*coef->eval());
                            first_occupied = true;
                        }
                        else{
                            lhs_second_var = *static_pointer_cast<var<double>>(qt_pair.second._p->first);
                            if(coef->is_positive()) lhs_second_scale = std::sqrt(coef->eval());
                            else lhs_second_scale = std::sqrt((-1)*coef->eval());
                            
                        }
                    }
                    else{
                        rhs_var = *static_pointer_cast<var<double>>(qt_pair.second._p->first);
                        if(coef->is_positive()) rhs_scale = std::sqrt(coef->eval());
                        else rhs_scale = std::sqrt((-1)*coef->eval());
                    }
                }
                
                //combine all the indices
                //How to combine safely when the index sets are same? ***********************************************************************************************************************
                auto inst_combined  = combine(*lhs_first_var._indices,*lhs_second_var._indices,*rhs_var._indices);
                auto inst_combined_partn = indices(inst_combined, hyper_idx);
                auto inst_combined_lambda = indices(inst_combined, lambda_idx);
                
                //collect the number of entries
                auto nb_entries_v1 = lhs_first_var._indices->get_nb_entries();
                auto nb_entries_v2 = lhs_second_var._indices->get_nb_entries();
                auto nb_entries_v3 = rhs_var._indices->get_nb_entries();
                
                //create the multipliers for the hyperplane
                param<double> v1_coef("v1_coef");
                param<double> v2_coef("v2_coef");
                param<double> v3_coef("v3_coef");
                
                v1_coef.in(inst_combined_lambda);
                v2_coef.in(inst_combined_lambda);
                v3_coef.in(inst_combined_lambda);
                
                // Lambda coefficient matrix when linking with partition variables
                param<> lambda_coef(c._name+"_lambda_linking_coefficients");
                // Partition coefficient matrix when linking with lambda variables
                param<> on_coef(c._name+"_partition_linking_coefficients");
                
                // we can use constraint_idx as the hyper_idx since there will be exactly |num_SOC_partitions| many linking constraints
                lambda_coef.in(indices(inst_combined_lambda, on_link_lambda_idx));
                on_coef.in(indices(inst_combined_partn, on_link_lambda_idx));
                
                //on variable constraint and definition
                add(on.in(inst_combined_partn));
                
                //create the summation constraint for partition
                onSum = sum(on.in_matrix(nb_entries_v1+nb_entries_v2+nb_entries_v3,1));
                add(onSum.in(*c._indices) == 1);
                
                //lambda variable constraint and definition
                add(lambda.in(inst_combined_lambda));
                
                //create the summation constraint for partition
                lambdaSum = sum(lambda.in_matrix(nb_entries_v1+nb_entries_v2+nb_entries_v3,1));
                add(lambdaSum.in(*c._indices) == 1);
                
                add(v1_sign.in(inst_combined));
                add(v2_sign.in(inst_combined));
                
                size_t nb_ins = lhs_first_var.get_nb_inst(); //get the number of instances
                
                // fill lambda_coef
                for (size_t inst = 0; inst< nb_ins; inst++){
                    // lhs_first
                    auto cur_var_id_lhs_first = lhs_first_var.get_id_inst(inst);
                    auto cur_var_idx_lhs_first = lhs_first_var._indices->_keys->at(cur_var_id_lhs_first);
                    
                    // lhs_second
                    auto cur_var_id_lhs_second = lhs_second_var.get_id_inst(inst);
                    auto cur_var_idx_lhs_second = lhs_second_var._indices->_keys->at(cur_var_id_lhs_second);
                    
                    // rhs
                    auto cur_var_id_rhs = rhs_var.get_id_inst(inst);
                    auto cur_var_idx_rhs = rhs_var._indices->_keys->at(cur_var_id_rhs);
                    
                    //for constraints
                    for (int i=0 ; i<num_SOC_partitions+1; ++i) {
                        string cur_idx = cur_var_idx_lhs_first+","+cur_var_idx_lhs_second+","+cur_var_idx_rhs+","+to_string(i+2)+","+to_string(i+1);
                        lambda_coef.set_val(cur_idx,1);
                    }
                }
                
                // fill on_coef
                for (size_t inst = 0; inst< nb_ins; inst++){
                    // lhs_first
                    auto cur_var_id_lhs_first = lhs_first_var.get_id_inst(inst);
                    auto cur_var_idx_lhs_first = lhs_first_var._indices->_keys->at(cur_var_id_lhs_first);
                    
                    // lhs_second
                    auto cur_var_id_lhs_second = lhs_second_var.get_id_inst(inst);
                    auto cur_var_idx_lhs_second = lhs_second_var._indices->_keys->at(cur_var_id_lhs_second);
                    
                    // rhs
                    auto cur_var_id_rhs = rhs_var.get_id_inst(inst);
                    auto cur_var_idx_rhs = rhs_var._indices->_keys->at(cur_var_id_rhs);
                    
                    //first and last constraint
                    string cur_idx = cur_var_idx_lhs_first+","+cur_var_idx_lhs_second+","+cur_var_idx_rhs+","+to_string(1)+","+to_string(1);
                    on_coef.set_val(cur_idx,1);
                    cur_idx = cur_var_idx_lhs_first+","+cur_var_idx_lhs_second+","+cur_var_idx_rhs+","+to_string(num_SOC_partitions)+","+to_string(num_SOC_partitions+1);
                    on_coef.set_val(cur_idx,1);
                    
                    //for constraints
                    for (int i=0 ; i<num_SOC_partitions-1; ++i) {
                        cur_idx = cur_var_idx_lhs_first+","+cur_var_idx_lhs_second+","+cur_var_idx_rhs+","+to_string(i+1)+","+to_string(i+2);
                        on_coef.set_val(cur_idx,1);
                        cur_idx = cur_var_idx_lhs_first+","+cur_var_idx_lhs_second+","+cur_var_idx_rhs+","+to_string(i+2)+","+to_string(i+2);
                        on_coef.set_val(cur_idx,1);
                    }
                }
                
                // add the link constraint between partition variables and lambda variables
                on_link_lambda = lambda.in_matrix(nb_entries_v1+nb_entries_v2+nb_entries_v3,1).from_ith(0,lambda_coef.get_matrix_ids(nb_entries_v1+nb_entries_v2+nb_entries_v3,1))*lambda_coef.in_matrix(nb_entries_v1+nb_entries_v2+nb_entries_v3,1) - on.in_matrix(nb_entries_v1+nb_entries_v2+nb_entries_v3,1).from_ith(0,on_coef.get_matrix_ids(nb_entries_v1+nb_entries_v2+nb_entries_v3,1)) * on_coef.in_matrix(nb_entries_v1+nb_entries_v2+nb_entries_v3,1);
                add(on_link_lambda.in(indices(*c._indices,on_link_lambda_idx)) <= 0);
                

                // fill lambda representation coefficients
                for (size_t inst = 0; inst< nb_ins; inst++){
                    
                    size_t id1;
                    if(lhs_first_var._indices->_ids == nullptr){
                        id1 = inst;
                    }
                    else id1 = lhs_first_var._indices->_ids->at(0).at(inst);
                    auto key1 = lhs_first_var._indices->_keys->at(id1);
                    
                    size_t id2;
                    if(lhs_second_var._indices->_ids == nullptr){
                        id2 = inst;
                    }
                    else id2 = lhs_second_var._indices->_ids->at(0).at(inst);
                    auto key2 = lhs_second_var._indices->_keys->at(id2);
                    
                    size_t id3;
                    if(rhs_var._indices->_ids == nullptr){
                        id3 = inst;
                    }
                    else id3 = rhs_var._indices->_ids->at(0).at(inst);
                    auto key3 = rhs_var._indices->_keys->at(id3);
                    
                    double radius1 = std::max(std::abs(lhs_first_var.get_lb(key1)), std::abs(lhs_first_var.get_ub(key1)))*lhs_first_scale;
                    double radius2 = std::max(std::abs(lhs_second_var.get_lb(key2)), std::abs(lhs_first_var.get_ub(key2)))*lhs_second_scale;
                    double radius3 = rhs_var.get_ub(key3)*rhs_scale;
                    
                    double radius = std::min( std::sqrt(std::pow(radius1,2)+std::pow(radius2,2)), radius3);
                    
                    // lhs_first
                    auto cur_var_id_lhs_first = lhs_first_var.get_id_inst(inst);
                    auto cur_var_idx_lhs_first = lhs_first_var._indices->_keys->at(cur_var_id_lhs_first);
                    
                    // lhs_second
                    auto cur_var_id_lhs_second = lhs_second_var.get_id_inst(inst);
                    auto cur_var_idx_lhs_second = lhs_second_var._indices->_keys->at(cur_var_id_lhs_second);
                    
                    // rhs
                    auto cur_var_id_rhs = rhs_var.get_id_inst(inst);
                    auto cur_var_idx_rhs = rhs_var._indices->_keys->at(cur_var_id_rhs);
                    
                    string cur_idx = cur_var_idx_lhs_first+","+cur_var_idx_lhs_second+","+cur_var_idx_rhs+","+to_string(1);
                    v1_coef.set_val(cur_idx,0);
                    v2_coef.set_val(cur_idx,0);
                    v3_coef.set_val(cur_idx,0);
                    
                    for (int i=0 ; i<num_SOC_partitions+1; ++i) {
                        //calculate the hyperplane coefficients
                        auto v1_val = radius*std::cos(M_PI*i/(2*num_SOC_partitions));
                        auto v2_val = radius*std::sin(M_PI*i/(2*num_SOC_partitions));
                        auto v3_val = radius;
                        cur_idx = cur_var_idx_lhs_first+","+cur_var_idx_lhs_second+","+cur_var_idx_rhs+","+to_string(i+2);
                        v1_coef.set_val(cur_idx,v1_val/lhs_first_scale);
                        v2_coef.set_val(cur_idx,v2_val/lhs_second_scale);
                        v3_coef.set_val(cur_idx,v3_val/rhs_scale);
                    }
                }
                
                //add the constraints involving the individual variables (tying lambda to the variables as lower and upper bounds)
                v1_rep_pos = lhs_first_var + 100 * v1_sign - v1_coef.in_matrix(nb_entries_v1+nb_entries_v2+nb_entries_v3,1) * lambda.in_matrix(nb_entries_v1+nb_entries_v2+nb_entries_v3,1);
                v1_rep_pos.in(*c._indices) >= 0;
                add(v1_rep_pos);
                
                v1_rep_neg = lhs_first_var - 100 * (1 - v1_sign) + v1_coef.in_matrix(nb_entries_v1+nb_entries_v2+nb_entries_v3,1) * lambda.in_matrix(nb_entries_v1+nb_entries_v2+nb_entries_v3,1);
                v1_rep_neg.in(*c._indices) <= 0;
                add(v1_rep_neg);
                
                v2_rep_pos = lhs_second_var + 100 * v2_sign - v2_coef.in_matrix(nb_entries_v1+nb_entries_v2+nb_entries_v3,1) * lambda.in_matrix(nb_entries_v1+nb_entries_v2+nb_entries_v3,1);
                v2_rep_pos.in(*c._indices) >= 0;
                add(v2_rep_pos);
                
                v2_rep_neg = lhs_second_var - 100 * (1 - v2_sign) + v2_coef.in_matrix(nb_entries_v1+nb_entries_v2+nb_entries_v3,1) * lambda.in_matrix(nb_entries_v1+nb_entries_v2+nb_entries_v3,1);
                v2_rep_neg.in(*c._indices) <= 0;
                add(v2_rep_neg);
                
                //this is the variable on the rhs
                v3_rep = rhs_var - v3_coef.in_matrix(nb_entries_v1+nb_entries_v2+nb_entries_v3,1) * lambda.in_matrix(nb_entries_v1+nb_entries_v2+nb_entries_v3,1);
                v3_rep.in(*c._indices) <= 0;
                add(v3_rep);
            }
            
            else if (is_rotated_SOC){ //this will follow bilinear scheme
                
                //create the variables
                var<> quad_var;
                var<> bln_first_var;
                var<> bln_second_var;
                
                //create the scaling factors for the variables
                double quad_scale;
                double bln_scale;
                
                for (auto &qt_pair: *c._qterms) {
                    if (!qt_pair.second._p->first->is_double()) {
                        throw invalid_argument("Current hyperplanes only support double type variables!");
                    }
                    if (!qt_pair.second._coef->is_number()) { /*means coef is not a number*/
                        throw invalid_argument("Current hyperplanes only support constant coefficients for the variables");
                    }
                    if (qt_pair.second._p->first!=qt_pair.second._p->second) {
                        if (!qt_pair.second._p->second->is_double()) {
                            throw invalid_argument("Current hyperplanes only support double type variables!");
                        }
                        //we do not check the sign of the bilinear term, since current assumption is that the sign is negative
                        bln_first_var = *static_pointer_cast<var<double>>(qt_pair.second._p->first);
                        bln_second_var = *static_pointer_cast<var<double>>(qt_pair.second._p->second);
                        auto coef = static_pointer_cast<constant<type>>(qt_pair.second._coef);
                        if(coef->is_positive()) bln_scale = std::sqrt(coef->eval())/2; //dividing to two since the standard form has 1/4 as the multiplier
                        else bln_scale = std::sqrt((-1)*coef->eval())/(2);
                    }
                    else{
                        quad_var = *static_pointer_cast<var<double>>(qt_pair.second._p->first);
                        auto coef = static_pointer_cast<constant<type>>(qt_pair.second._coef);
                        if(coef->is_positive()) quad_scale = std::sqrt(coef->eval());
                        else quad_scale = std::sqrt((-1)*coef->eval());
                    }
                }
                
                //combine all the indices
                auto inst_combined  = combine(*quad_var._indices,*bln_first_var._indices,*bln_second_var._indices);
                auto inst_combined_partn = indices(inst_combined,hyper_idx);
                auto inst_combined_lambda = indices(inst_combined, lambda_idx);
                
                //collect the number of entries
                auto nb_entries_v1 = quad_var._indices->get_nb_entries();
                auto nb_entries_v2 = bln_first_var._indices->get_nb_entries();
                auto nb_entries_v3 = bln_second_var._indices->get_nb_entries();
                
                //create the multipliers for the hyperplane
                param<double> quad_coef("quad_coef");
                param<double> first_minus_second_coef("first_minus_second_coef");
                param<double> first_plus_second_coef("first_plus_second_coef");
                
                quad_coef.in(inst_combined_lambda);
                first_minus_second_coef.in(inst_combined_lambda);
                first_plus_second_coef.in(inst_combined_lambda);
                
                // Lambda coefficient matrix when linking with partition variables
                param<> lambda_coef(c._name+"_lambda_linking_coefficients");
                // Partition coefficient matrix when linking with lambda variables
                param<> on_coef(c._name+"_partition_linking_coefficients");
                
                // we can use constraint_idx as the hyper_idx since there will be exactly |num_SOC_partitions| many linking constraints
                lambda_coef.in(indices(inst_combined_lambda, on_link_lambda_idx));
                on_coef.in(indices(inst_combined_partn, on_link_lambda_idx));
                
                //on variable constraint and definition
                add(on.in(inst_combined_partn));
                
                //create the summation constraint for partition
                onSum = sum(on.in_matrix(nb_entries_v1+nb_entries_v2+nb_entries_v3,1));
                add(onSum.in(*c._indices) == 1);
                
                //lambda variable constraint and definition
                add(lambda.in(inst_combined_lambda));
                
                //create the summation constraint for partition
                lambdaSum = sum(lambda.in_matrix(nb_entries_v1+nb_entries_v2+nb_entries_v3,1));
                add(lambdaSum.in(*c._indices) == 1);
                
                add(v1_sign.in(inst_combined));
                add(v2_sign.in(inst_combined));
                
                size_t nb_ins = quad_var.get_nb_inst(); //get the number of instances
                
                // fill lambda_coef
                for (size_t inst = 0; inst< nb_ins; inst++){
                    // quad
                    auto cur_var_id_quad = quad_var.get_id_inst(inst);
                    auto cur_var_idx_quad = quad_var._indices->_keys->at(cur_var_id_quad);
                    
                    // bln_first
                    auto cur_var_id_bln_first = bln_first_var.get_id_inst(inst);
                    auto cur_var_idx_bln_first = bln_first_var._indices->_keys->at(cur_var_id_bln_first);
                    
                    // bln_second
                    auto cur_var_id_bln_second = bln_second_var.get_id_inst(inst);
                    auto cur_var_idx_bln_second = bln_second_var._indices->_keys->at(cur_var_id_bln_second);
                    
                    //for constraints
                    for (int i=0 ; i<num_SOC_partitions+1; ++i) {
                        string cur_idx = cur_var_idx_quad+","+cur_var_idx_bln_first+","+cur_var_idx_bln_second+","+to_string(i+2)+","+to_string(i+1);
                        lambda_coef.set_val(cur_idx,1);
                    }
                }
                
                // fill on_coef
                for (size_t inst = 0; inst< nb_ins; inst++){
                    // quad
                    auto cur_var_id_quad = quad_var.get_id_inst(inst);
                    auto cur_var_idx_quad = quad_var._indices->_keys->at(cur_var_id_quad);
                    
                    // bln_first
                    auto cur_var_id_bln_first = bln_first_var.get_id_inst(inst);
                    auto cur_var_idx_bln_first = bln_first_var._indices->_keys->at(cur_var_id_bln_first);
                    
                    // bln_second
                    auto cur_var_id_bln_second = bln_second_var.get_id_inst(inst);
                    auto cur_var_idx_bln_second = bln_second_var._indices->_keys->at(cur_var_id_bln_second);
                    
                    //first and last constraints
                    string cur_idx = cur_var_idx_quad+","+cur_var_idx_bln_first+","+cur_var_idx_bln_second+","+to_string(1)+","+to_string(1);
                    on_coef.set_val(cur_idx,1);
                    cur_idx = cur_var_idx_quad+","+cur_var_idx_bln_first+","+cur_var_idx_bln_second+","+to_string(num_SOC_partitions)+","+to_string(num_SOC_partitions+1);
                    on_coef.set_val(cur_idx,1);
                    
                    //for constraints
                    for (int i=0 ; i<num_SOC_partitions-1; ++i) {
                        cur_idx = cur_var_idx_quad+","+cur_var_idx_bln_first+","+cur_var_idx_bln_second+","+to_string(i+1)+","+to_string(i+2);
                        on_coef.set_val(cur_idx,1);
                        cur_idx = cur_var_idx_quad+","+cur_var_idx_bln_first+","+cur_var_idx_bln_second+","+to_string(i+2)+","+to_string(i+2);
                        on_coef.set_val(cur_idx,1);
                    }
                }
                
                // add the link constraint between partition variables and lambda variables
                on_link_lambda = lambda.in_matrix(nb_entries_v1+nb_entries_v2+nb_entries_v3,1).from_ith(0,lambda_coef.get_matrix_ids(nb_entries_v1+nb_entries_v2+nb_entries_v3,1))*lambda_coef.in_matrix(nb_entries_v1+nb_entries_v2+nb_entries_v3,1) - on.in_matrix(nb_entries_v1+nb_entries_v2+nb_entries_v3,1).from_ith(0,on_coef.get_matrix_ids(nb_entries_v1+nb_entries_v2+nb_entries_v3,1)) * on_coef.in_matrix(nb_entries_v1+nb_entries_v2+nb_entries_v3,1);
                add(on_link_lambda.in(indices(*c._indices,on_link_lambda_idx)) <= 0);
                
                // fill lambda representation coefficients
                for (size_t inst = 0; inst< nb_ins; inst++){
                    
                    size_t id1;
                    if(quad_var._indices->_ids == nullptr){
                        id1 = inst;
                    }
                    else id1 = quad_var._indices->_ids->at(0).at(inst);
                    auto key1 = quad_var._indices->_keys->at(id1);
                    
                    size_t id2;
                    if(bln_first_var._indices->_ids == nullptr){
                        id2 = inst;
                    }
                    else id2 = bln_first_var._indices->_ids->at(0).at(inst);
                    auto key2 = bln_first_var._indices->_keys->at(id2);
                    
                    size_t id3;
                    if(bln_second_var._indices->_ids == nullptr){
                        id3 = inst;
                    }
                    else id3 = bln_second_var._indices->_ids->at(0).at(inst);
                    auto key3 = bln_second_var._indices->_keys->at(id3);
                    
                    double radius1 = std::max(std::abs(quad_var.get_lb(key1)), std::abs(quad_var.get_ub(key1)))*quad_scale;
                    double radius2 = std::max(std::abs(bln_first_var.get_lb(key2) - bln_second_var.get_ub(key3)), std::abs(bln_first_var.get_ub(key2) - bln_second_var.get_lb(key3)))*bln_scale;
                    double radius3 = (bln_first_var.get_ub(key2) + bln_second_var.get_ub(key3))*bln_scale;
                    
                    double radius = std::min( std::sqrt(std::pow(radius1,2)+std::pow(radius2,2)), radius3);
                    
                    // quad
                    auto cur_var_id_quad = quad_var.get_id_inst(inst);
                    auto cur_var_idx_quad = quad_var._indices->_keys->at(cur_var_id_quad);
                    
                    // bln_first
                    auto cur_var_id_bln_first = bln_first_var.get_id_inst(inst);
                    auto cur_var_idx_bln_first = bln_first_var._indices->_keys->at(cur_var_id_bln_first);
                    
                    // bln_second
                    auto cur_var_id_bln_second = bln_second_var.get_id_inst(inst);
                    auto cur_var_idx_bln_second = bln_second_var._indices->_keys->at(cur_var_id_bln_second);
                    
                    string cur_idx = cur_var_idx_quad+","+cur_var_idx_bln_first+","+cur_var_idx_bln_second+","+to_string(1);
                    quad_coef.set_val(cur_idx,0);
                    first_minus_second_coef.set_val(cur_idx,0);
                    first_plus_second_coef.set_val(cur_idx,0);
                    
                    for (int i=0 ; i<num_SOC_partitions+1; ++i) {
                        //calculate the hyperplane coefficients
                        auto quad_val = radius*std::cos(M_PI*i/(2*num_SOC_partitions));
                        auto first_minus_second_val = radius*std::sin(M_PI*i/(2*num_SOC_partitions));
                        auto first_plus_second_val = radius;
                        cur_idx = cur_var_idx_quad+","+cur_var_idx_bln_first+","+cur_var_idx_bln_second+","+to_string(i+2);
                        quad_coef.set_val(cur_idx,quad_val/quad_scale);
                        first_minus_second_coef.set_val(cur_idx,first_minus_second_val/bln_scale);
                        first_plus_second_coef.set_val(cur_idx,first_plus_second_val/bln_scale);
                    }
                }
                
                //add the constraints involving the individual variables (tying lambda to the variables as lower and upper bounds)
                v1_rep_pos = quad_var + 100 * v1_sign - quad_coef.in_matrix(nb_entries_v1+nb_entries_v2+nb_entries_v3,1) * lambda.in_matrix(nb_entries_v1+nb_entries_v2+nb_entries_v3,1);
                v1_rep_pos.in(*c._indices) >= 0;
                add(v1_rep_pos);
                
                v1_rep_neg = quad_var - 100 * (1 - v1_sign) + quad_coef.in_matrix(nb_entries_v1+nb_entries_v2+nb_entries_v3,1) * lambda.in_matrix(nb_entries_v1+nb_entries_v2+nb_entries_v3,1);
                v1_rep_neg.in(*c._indices) <= 0;
                add(v1_rep_neg);
                
                v2_rep_pos = (bln_first_var - bln_second_var) + 100 * v2_sign - first_minus_second_coef.in_matrix(nb_entries_v1+nb_entries_v2+nb_entries_v3,1) * lambda.in_matrix(nb_entries_v1+nb_entries_v2+nb_entries_v3,1);
                v2_rep_pos.in(*c._indices) >= 0;
                add(v2_rep_pos);
                
                v2_rep_neg = (bln_first_var - bln_second_var)  - 100 * (1 - v2_sign) + first_minus_second_coef.in_matrix(nb_entries_v1+nb_entries_v2+nb_entries_v3,1) * lambda.in_matrix(nb_entries_v1+nb_entries_v2+nb_entries_v3,1);
                v2_rep_neg.in(*c._indices) <= 0;
                add(v2_rep_neg);
                
                //this is the second part of the bilinear, which ends up in the rhs of the SOC constraint as the bounding term in the standard format
                v3_rep = (bln_first_var + bln_second_var) - first_plus_second_coef.in_matrix(nb_entries_v1+nb_entries_v2+nb_entries_v3,1) * lambda.in_matrix(nb_entries_v1+nb_entries_v2+nb_entries_v3,1);
                v3_rep.in(*c._indices) <= 0;
                add(v3_rep);
                
            }
            
        }
        
        // INPUT: an SOC type constraint, and total number of partitions(binary variables)
        // OUTPUT: disjunctive union of hyperplanes as an inner approximation to the SOC, where the disjunctive union is made by on/off-formulation by Hijazi et. al.
        template<typename T=type> //function for creating hyperplanes to have an inner approximation of an SOC
        void add_on_off_SOC_hyperplanes(Constraint<type>& c, int num_SOC_partitions) { //currently the function asssumes there are only qterms!
            
            DebugOn("SOC_hyperplane function!" << endl);
            c.print();
            
            auto is_rotated_SOC = c.check_rotated_soc(); //collect the information about the cone
            auto is_SOC = c.check_soc();
            
            //create hyperplane indices
            indices hyper_idx("hyper_idx");
            for (int i=0; i<num_SOC_partitions; ++i) {
                hyper_idx.add(to_string(i+1));
            }
            //get the combined index set
            auto inst_hyper = indices(*c._indices,hyper_idx);

            Constraint<type> SOC_hyperplanes(c._name + "_hyperplane"); //create the hyperplane constraint
            var<int> on(c._name + "_binary",0,1); //create the partition variable
            Constraint<type> onSum(c._name + "_binarySum"); //create the partition assignment constraint
            
            if (is_SOC){ //this will follow the standard creation of the hyperplane
                
                //create the variables
                var<> lhs_first_var;
                var<> lhs_second_var;
                var<> rhs_var;
                
                //create the scaling factors for the variables
                double lhs_first_scale;
                double lhs_second_scale;
                double rhs_scale;
                
                //flag for assignment
                bool first_occupied = false;
                
                for (auto &qt_pair: *c._qterms) {
                    if (!qt_pair.second._p->first->is_double()) {
                        throw invalid_argument("Current hyperplanes only support double type variables!");
                    }
                    if (!qt_pair.second._coef->is_number()) { /*means coef is not a number*/
                        throw invalid_argument("Current hyperplanes only support constant coefficients for the variables");
                    }
                    auto sign = qt_pair.second._sign;
                    auto coef = static_pointer_cast<constant<type>>(qt_pair.second._coef);
                    if (sign ^ (coef->is_negative())) {
                        if (!first_occupied){
                            lhs_first_var = *static_pointer_cast<var<double>>(qt_pair.second._p->first);
                            if(coef->is_positive()) lhs_first_scale = std::sqrt(coef->eval());
                            else lhs_first_scale = std::sqrt((-1)*coef->eval());
                            first_occupied = true;
                        }
                        else{
                            lhs_second_var = *static_pointer_cast<var<double>>(qt_pair.second._p->first);
                            if(coef->is_positive()) lhs_second_scale = std::sqrt(coef->eval());
                            else lhs_second_scale = std::sqrt((-1)*coef->eval());
                            
                        }
                    }
                    else{
                        rhs_var = *static_pointer_cast<var<double>>(qt_pair.second._p->first);
                        if(coef->is_positive()) rhs_scale = std::sqrt(coef->eval());
                        else rhs_scale = std::sqrt((-1)*coef->eval());
                    }
                }
               
                
                //combine all the indices
                auto inst_combined  = combine(*lhs_first_var._indices,*lhs_second_var._indices,*rhs_var._indices);
                auto inst_combined_partn = indices(inst_combined,hyper_idx);
                
                //collect the number of entries
                auto nb_entries_v1 = lhs_first_var._indices->get_nb_entries();
                auto nb_entries_v2 = lhs_second_var._indices->get_nb_entries();
                auto nb_entries_v3 = rhs_var._indices->get_nb_entries();
                
                //create the multipliers for the hyperplane
                param<double> lhs_first_coef("lhs_first_coef");
                param<double> lhs_second_coef("lhs_second_coef");
                param<double> rhs_coef("rhs_coef");
                
                lhs_first_coef.in(inst_combined_partn);
                lhs_second_coef.in(inst_combined_partn);
                rhs_coef.in(inst_combined_partn);
                
                add(on.in(inst_combined_partn));

                //create the summation constraint for partition
                onSum = sum(on.in_matrix(nb_entries_v1+nb_entries_v2+nb_entries_v3,1));
                add(onSum.in(*c._indices) == 1);
                
                size_t nb_ins = lhs_first_var.get_nb_inst(); //get the number of instances
                
                // fill hyperplane coefficients
                for (int i=0 ; i<num_SOC_partitions; ++i) {
                    //calculate the hyperplane coefficients
                    auto rhs_val = std::cos(2*M_PI*i/num_SOC_partitions)*std::sin(2*M_PI*(i+1)/num_SOC_partitions) - std::cos(2*M_PI*(i+1)/num_SOC_partitions)*std::sin(2*M_PI*i/num_SOC_partitions);
                    auto lhs_first_val = std::sin(2*M_PI*i/num_SOC_partitions)-std::sin(2*M_PI*(i+1)/num_SOC_partitions);
                    auto lhs_second_val = std::cos(2*M_PI*(i+1)/num_SOC_partitions)-std::cos(2*M_PI*i/num_SOC_partitions);
                    if (rhs_val < 0)
                    {
                        rhs_val = -rhs_val;
                        lhs_first_val = -lhs_first_val;
                        lhs_second_val = -lhs_second_val;
                    }
                    
                    for (size_t inst = 0; inst< nb_ins; inst++){
                        // lhs_first
                        auto cur_var_id_lhs_first = lhs_first_var.get_id_inst(inst);
                        auto cur_var_idx_lhs_first = lhs_first_var._indices->_keys->at(cur_var_id_lhs_first);
                        
                        // lhs_second
                        auto cur_var_id_lhs_second = lhs_second_var.get_id_inst(inst);
                        auto cur_var_idx_lhs_second = lhs_second_var._indices->_keys->at(cur_var_id_lhs_second);
                        
                        // rhs
                        auto cur_var_id_rhs = rhs_var.get_id_inst(inst);
                        auto cur_var_idx_rhs = rhs_var._indices->_keys->at(cur_var_id_rhs);
                        
                        string cur_idx = cur_var_idx_lhs_first+","+cur_var_idx_lhs_second+","+cur_var_idx_rhs+","+to_string(i+1);
                        lhs_first_coef.set_val(cur_idx,lhs_first_val*lhs_first_scale);
                        lhs_second_coef.set_val(cur_idx,lhs_second_val*lhs_second_scale);
                        rhs_coef.set_val(cur_idx,rhs_val*rhs_scale);
                    }
                }
                
                // set the _in_SOC_partn to true
                lhs_first_var._in_SOC_partn = true;
                lhs_second_var._in_SOC_partn = true;
                rhs_var._in_SOC_partn = true;
              
                SOC_hyperplanes = lhs_first_var.from_ith(0, inst_combined_partn) * lhs_first_coef + lhs_second_var.from_ith(nb_entries_v1, inst_combined_partn) * lhs_second_coef + rhs_var.from_ith(nb_entries_v1 + nb_entries_v2, inst_combined_partn) * rhs_coef;
                SOC_hyperplanes.in(inst_hyper) <= 0;
                //the third argument for this function decides full facets or only Big-M type of version
                add_on_off_multivariate_refined(SOC_hyperplanes, on, false);
                
                // set the _in_SOC_partn to false back for removing the confusion in get_on_off_coefficients
                lhs_first_var._in_SOC_partn = false;
                lhs_second_var._in_SOC_partn = false;
                rhs_var._in_SOC_partn = false;
            }
            
            else if (is_rotated_SOC){ //this will follow bilinear scheme
                
                //create the variables
                var<> quad_var;
                var<> bln_first_var;
                var<> bln_second_var;
                
                //create the scaling factors for the variables
                double quad_scale;
                double bln_scale;
                
                for (auto &qt_pair: *c._qterms) {
                    if (!qt_pair.second._p->first->is_double()) {
                        throw invalid_argument("Current hyperplanes only support double type variables!");
                    }
                    if (!qt_pair.second._coef->is_number()) { /*means coef is not a number*/
                        throw invalid_argument("Current hyperplanes only support constant coefficients for the variables");
                    }
                    if (qt_pair.second._p->first!=qt_pair.second._p->second) {
                        if (!qt_pair.second._p->second->is_double()) {
                            throw invalid_argument("Current hyperplanes only support double type variables!");
                        }
                            //we do not check the sign of the bilinear term, since current assumption is that the sign is negative
                            bln_first_var = *static_pointer_cast<var<double>>(qt_pair.second._p->first);
                            bln_second_var = *static_pointer_cast<var<double>>(qt_pair.second._p->second);
                            auto coef = static_pointer_cast<constant<type>>(qt_pair.second._coef);
                            if(coef->is_positive()) bln_scale = std::sqrt(coef->eval())/2; //dividing to two since the standard form has 1/4 as the multiplier
                            else bln_scale = std::sqrt((-1)*coef->eval())/(2);
                    }
                    else{
                        quad_var = *static_pointer_cast<var<double>>(qt_pair.second._p->first);
                        auto coef = static_pointer_cast<constant<type>>(qt_pair.second._coef);
                        if(coef->is_positive()) quad_scale = std::sqrt(coef->eval());
                        else quad_scale = std::sqrt((-1)*coef->eval());
                    }
                }
                
                
                //combine all the indices
                auto inst_combined  = combine(*quad_var._indices,*bln_first_var._indices,*bln_second_var._indices);
                auto inst_combined_partn = indices(inst_combined,hyper_idx);
                
                //collect the number of entries
                auto nb_entries_v1 = quad_var._indices->get_nb_entries();
                auto nb_entries_v2 = bln_first_var._indices->get_nb_entries();
                auto nb_entries_v3 = bln_second_var._indices->get_nb_entries();
                
                //create the multipliers for the hyperplane
                param<double> quad_coef("quad_coef");
                param<double> first_minus_second_coef("first_minus_second_coef");
                param<double> first_plus_second_coef("first_plus_second_coef");
                
                quad_coef.in(inst_combined_partn);
                first_minus_second_coef.in(inst_combined_partn);
                first_plus_second_coef.in(inst_combined_partn);
                
                
                add(on.in(inst_combined_partn));
                //create the summation constraint for partition
                onSum = sum(on.in_matrix(nb_entries_v1+nb_entries_v2+nb_entries_v3,1));
                add(onSum.in(*c._indices) == 1);
                
                size_t nb_ins = quad_var.get_nb_inst(); //get the number of instances
                
                // fill hyperplane coefficients
                for (int i=0 ; i<num_SOC_partitions; ++i) {
                    //calculate the hyperplane coefficients
                    auto first_plus_second_val = std::cos(2*M_PI*i/num_SOC_partitions)*std::sin(2*M_PI*(i+1)/num_SOC_partitions) - std::cos(2*M_PI*(i+1)/num_SOC_partitions)*std::sin(2*M_PI*i/num_SOC_partitions);
                    auto quad_val = std::sin(2*M_PI*i/num_SOC_partitions)-std::sin(2*M_PI*(i+1)/num_SOC_partitions);
                    auto first_minus_second_val = std::cos(2*M_PI*(i+1)/num_SOC_partitions)-std::cos(2*M_PI*i/num_SOC_partitions);
                    if (first_plus_second_val < 0)
                    {
                        first_plus_second_val = -first_plus_second_val;
                        quad_val = -quad_val;
                        first_minus_second_val = -first_minus_second_val;
                    }
                    
                    for (size_t inst = 0; inst< nb_ins; inst++){
                        // quad
                        auto cur_var_id_quad = quad_var.get_id_inst(inst);
                        auto cur_var_idx_quad = quad_var._indices->_keys->at(cur_var_id_quad);
                        
                        // bln_first
                        auto cur_var_id_bln_first = bln_first_var.get_id_inst(inst);
                        auto cur_var_idx_bln_first = bln_first_var._indices->_keys->at(cur_var_id_bln_first);
                        
                        // bln_second
                        auto cur_var_id_bln_second = bln_second_var.get_id_inst(inst);
                        auto cur_var_idx_bln_second = bln_second_var._indices->_keys->at(cur_var_id_bln_second);
                        
                        string cur_idx = cur_var_idx_quad+","+cur_var_idx_bln_first+","+cur_var_idx_bln_second+","+to_string(i+1);
                        quad_coef.set_val(cur_idx,quad_val*quad_scale);
                        first_minus_second_coef.set_val(cur_idx,first_minus_second_val*bln_scale);
                        first_plus_second_coef.set_val(cur_idx,first_plus_second_val*bln_scale);
                    }
                }
                
                // set the _in_SOC_partn to true
                quad_var._in_SOC_partn = true;
                bln_first_var._in_SOC_partn = true;
                bln_second_var._in_SOC_partn = true;
                
                SOC_hyperplanes = quad_var.from_ith(0, inst_combined_partn) * quad_coef + (bln_first_var.from_ith(nb_entries_v1, inst_combined_partn) - bln_second_var.from_ith(nb_entries_v1+nb_entries_v2, inst_combined_partn) ) * first_minus_second_coef + (bln_first_var.from_ith(nb_entries_v1, inst_combined_partn) + bln_second_var.from_ith(nb_entries_v1+nb_entries_v2, inst_combined_partn) ) * first_plus_second_coef;
                SOC_hyperplanes.in(inst_hyper) <= 0;
                //the third argument for this function decides full facets or only Big-M type of version
                add_on_off_multivariate_refined(SOC_hyperplanes, on, false);
                
                // set the _in_SOC_partn to false back for removing the confusion in get_on_off_coefficients
                quad_var._in_SOC_partn = false;
                bln_first_var._in_SOC_partn = false;
                bln_second_var._in_SOC_partn = false;
            }
        }
        
        //optimality based bound tightening procedure for a given mathematical formulation
        //INPUT: a given mathematical model, tolerances, maximum number of iterations, max amount of CPU time, and an upper bound for the current formulation to further tighten the bounds
        template<typename T=type,
        typename std::enable_if<is_same<type,double>::value>::type* = nullptr>
        void run_obbt(double max_time = 1000, unsigned max_iter=1e4, const pair<bool,double>& upper_bound = make_pair<bool,double>(false,0), unsigned precision=6);
        
        
//        void add_on_off(var<>& v, var<bool>& on){
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
//        }
        
        
        
        //        void add_on_off_McCormick(std::string name, var<>& v, var<>& v1, var<>& v2, var<bool>& on) {
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
        //        }
        
        
        //functions for the evaluation of the objective function
        template<class T=type>
        inline type eval(const shared_ptr<constant_>& c, size_t i=0) {
            return _obj->eval(c,i);
        }
        
        template<class T=type>
        inline type eval(const shared_ptr<constant_>& c, size_t i, size_t j){
            return _obj->eval(c,i,j);
        }
        
        
    };
    
    //    void compute_constrs(vector<Constraint*>& v, double* res, unsigned i, unsigned j);
    
    //    template<typename T>
    //    pair<shared_ptr<func_>, ObjectiveType> max(const func<T>& f){
    //        auto fcpy = f.copy();
    //        fcpy->allocate_mem();
    //        return make_pair<>(fcpy,maximize);
    //    };
    //
    //    template<typename T>
    //    pair<shared_ptr<func_>, ObjectiveType> min(const func<T>& f){
    //        auto fcpy = f.copy();
    //        f->_val->resize(1);
    //        return make_pair<>(fcpy,minimize);
    //    };
    //
    //    template<typename T>
    //    pair<shared_ptr<func_>, ObjectiveType> min(func<T>&& f){
    //        f->_val->resize(1);
    //        auto fcpy = move(f.copy());
    //        return make_pair<>(fcpy,minimize);
    //    };
    
    template<typename type = double>
    class Program{
    public:
        //        virtual void update_model(){};
        string _status;
    };
    
    
}




#endif /* model_hpp */
