//
//  func.h
//  Gravity
//
//  Created by Hijazi, Hassan on 21/10/16.
//
//

#ifndef func_h
#define func_h


#include <gravity/expr.h>
#include <gravity/poly.h>
#include <gravity/var.h>
#include <gravity/Auxiliary.h>
#include <stdio.h>
#include <map>
#include <math.h>
#include <iterator>
#include <queue>
#include <list>
#include <limits>
#include <set>
//
using namespace std;

namespace gravity {
    //
    //
    //
    /** Backbone class for function */
    class func_ : public constant_{
    private:
        
        
    public:
        FType                                                             _ftype = const_; /**< Function type, e.g., constant, linear, quadratic... >>**/
        NType                                                             _return_type = double_; /**< Return type, e.g., bool, integer, complex... >>**/
        
        shared_ptr<map<string, pair<shared_ptr<param_>, unsigned>>>       _params = nullptr;/**< Set of parameters in current function, stored as a map <parameter name, <paramter pointer, number of times it appears in function>>**/
        shared_ptr<map<string, pair<shared_ptr<param_>, unsigned>>>       _vars = nullptr;/**< Set of variables in current function, stored as a map <variable name, <variable pointer, number of times it appears in function>>**/
        
        shared_ptr<constant_>                                             _cst = nullptr;/**< Constant part of the function */
        shared_ptr<map<string, lterm>>                                    _lterms = nullptr; /**< Set of linear terms, stored as a map <string describing term, term>. */
        shared_ptr<map<string, qterm>>                                    _qterms = nullptr; /**< Set of quadratic terms, stored as a map <string describing term, term>.  */
        shared_ptr<map<string, pterm>>                                    _pterms = nullptr; /**< Set of polynomial terms, stored as a map <string describing term, term>.  */
        //        map<string, expr*>*                    _DAG = nullptr; /**< Map of experssions stored in the expression tree (a Directed Acyclic Graph) */
        //        deque<shared_ptr<expr>>*               _queue = nullptr; /**< A queue storing the expression tree from the leaves to the root (the root is stored at the end of the queue)*/
        Convexity                                                         _all_convexity = linear_; /**< If all instances of this function have the same convexity type, it stores it here, i.e. linear, convex, concave, otherwise it stores unknown. >>**/
        Sign                                                              _all_sign = zero_; /**< If all instances of this function have the same sign, it stores it here, otherwise it stores unknown. >>**/
        
        shared_ptr<vector<Convexity>>                                     _convexity = nullptr; /**< Vector of convexity types, i.e., linear, convex, concave or unknown. This is a vector since a function can have multiple instances (different constants coefficients, and bounds, but same structure) >>**/
        shared_ptr<vector<Sign>>                                          _sign = nullptr; /**< vector storing the sign of return value if known. >>**/
        shared_ptr<map<size_t, set<size_t>>>                              _hess_link = nullptr; /**< Set of variables linked to one another in the hessian, stored by variable ids  */
        
        bool                                                              _new = true; /**< Will become false once this function is added to a program. Can be useful for iterative model solving. */
        bool                                                              _is_constraint = false;
        bool                                                              _is_hessian = false;
        bool                                                              _embedded = false; /**< If the function is embedded in a mathematical model or in another function, this is used for memory management. >>**/
        bool                                                              _evaluated = true;/**< If the function has already been evaluated, useful for constant funcs */
        string                                                            _to_str = "";/**< A string representation of the expression */
        
        size_t                                                            _nb_vars = 0; /**< Number of variables */
        
        size_t                                                            _nnz_j = 0; /**< Number of nonzeros in the Jacobian **/
        size_t                                                            _nnz_h = 0; /**< Number of nonzeros in the Hessian **/
        
        string                                         _name = "noname";
        shared_ptr<indices>                            _indices = nullptr; /*< If indexed, point to the indexing set */
        /** Accessors */
        FType get_ftype() const;
        inline NType get_return_type() const{ return _return_type;};
        
        map<size_t, set<size_t>>& get_hess_link() { return *_hess_link;};
        map<string, pair<shared_ptr<param_>, unsigned>>& get_vars() { return *_vars;};
        map<string, pair<shared_ptr<param_>, unsigned>>& get_params() { return *_params;};
        
        /** true/false statements */
        bool is_indexed() const{
            return (_indices && _indices->_ids);
        }
        bool has_var(const param_& v) const;
        bool has_var(const string& name) const;
        bool is_binary() const {
            return (_return_type==binary_);
        };
        
        bool is_integer() const {
            return (_return_type==integer_);
        };
        
        bool is_float() const {
            return (_return_type==float_);
        };
        
        bool is_double() const {
            return (_return_type==double_);
        };
        
        bool is_long() const {
            return (_return_type==long_);
        };
        
        bool is_complex() const {
            return (_return_type==complex_);
        };
        virtual bool is_zero() const{ return false;};
        virtual bool is_convex() const{
            return (_all_convexity==convex_ || _all_convexity==linear_);
        }
        
        virtual bool is_concave() const{
            return (_all_convexity==concave_ || _all_convexity==linear_);
        }
        
        bool check_soc();
        bool check_rotated_soc();
        
        bool is_convex(size_t idx) const;
        bool is_concave(size_t idx) const;
        bool is_linear() const;
        bool is_quadratic() const;
        bool is_polynomial() const;
        bool is_nonlinear() const;
        bool is_transposed() const;
        bool func_is_number() const{
            return (_vars->empty() && _params->empty());
        }
        
        bool is_unitary() const{
            return (_vars->size()==1);
        }
        
        bool has_square() const {
            for (auto &qt: *_qterms) {
                if (qt.second._p->first==qt.second._p->second && !qt.second._p->first->_is_transposed && !qt.second._coef_p1_tr) {
                    return true;
                }
            }
            return false;
        }
        
        
        /* Virtual functions */
        
        virtual shared_ptr<func_> fcopy() const{return nullptr;};
        
        /* Modifiers */
        
        
        
        
        
        /**
         Copy and embed derivatives of f.
         @param[in] f function to copy derivatives from.
         */
        void copy_derivatives(const func_& f);
        
        
        unsigned nb_occ_var(string name) const;/**< Returns the number of occurences the variable has in this function. */
        
        unsigned nb_occ_param(string name) const;/**< Returns the number of occurences the parameter has in this function. */
        
        void incr_occ_var(string str);/**< Increases the number of occurences the variable has in this function. */
        
        void incr_occ_param(string str);/**< Increases the number of occurences the parameter has in this function. */
        
        void decr_occ_var(string str, int nb=1);/**< Decreases the number of occurences the variable has in this function by nb. */
        
        void decr_occ_param(string str, int nb=1);/**< Decreases the number of occurences the parameter has in this function by nb. */
        
        
        map<string, lterm>& get_lterms() const{
            return *_lterms;
        }
        
        map<string, qterm>& get_qterms() const{
            return *_qterms;
        }
        
        map<string, pterm>& get_pterms() const{
            return *_pterms;
        }
        
        
        
        
        
        
        /**
         Returns a vector of monomials of degree d using the variables in the current function
         @param[in] d degree of monomials
         @return a vector of monomials of degree d using the variables in the current function
         */
        vector<pterm> get_monomials(unsigned d);
        
        
        qterm* get_square(shared_ptr<param_> p); /**< Returns the quadratic term containing a square of p or nullptr if none exists. **/
        
        
        
        
        
        /**
         Returns the number of variables per-instance.
         @param[in] instance number.
         @return number of variables per-instance.
         */
        size_t get_nb_vars(unsigned inst) const{
            size_t n = 0;
            for (auto &vp:*_vars) {
                if(vp.second.first->_is_vector || vp.second.first->is_matrix_indexed()){
                    n += vp.second.first->get_dim(inst);
                }
                else {
                    n += 1;
                }
            }
            return n;
        };
        
        size_t get_nb_vars() const{
            size_t n = 0;
            for (auto &vp:*_vars) {
                if(vp.second.first->_is_vector){
                    n += vp.second.first->get_dim();
                }
                else {
                    n += 1;
                }
            }
            return n;
        };
        
        size_t get_id_inst(size_t inst = 0) const {
            if (is_indexed()) {
                if(_indices->_ids->at(0).size() <= inst){
                    throw invalid_argument("func_::get_id_inst(size_t inst) inst is out of range");
                }
                return _indices->_ids->at(0).at(inst);
            }
            auto dim = get_dim();
            if(inst > dim-1){
                throw invalid_argument("func_::get_id_inst(size_t inst) inst is out of range");
            }
            return inst;
        };
        
        size_t get_id_inst(size_t inst1, size_t inst2) const {
            if (is_matrix_indexed()) {
                if (_indices->_ids->size()<=inst1) {
                    throw invalid_argument("get_id_inst(size_t inst1, size_t inst2) inst1 out of range\n");
                }
                if (_indices->_ids->at(inst1).size()<=inst2) {
                    throw invalid_argument("get_id_inst(size_t inst1, size_t inst2) inst2 out of range\n");
                }
                return _indices->_ids->at(inst1).at(inst2);
            }
            throw invalid_argument("Calling get_id_inst(size_t inst1, size_t inst2) on a non-indexed param\n");
        };
        
        
        
        /**
         Returns a pointer to the constant part of the function.
         @return a pointer to the constant part of the function.
         */
        shared_ptr<constant_> get_cst() const;
        
        /**
         Returns a pointer to the variable matching the name provided.
         @param[in] name variable name.
         @return a pointer to the variable matching the name provided.
         */
        shared_ptr<param_> get_var(const string& name) const;
        
        /**
         Returns a pointer to the variable matching the index provided.
         @param[in] idx variable index.
         @return a pointer to the variable matching the index provided.
         */
        shared_ptr<param_> get_var(size_t idx) const;
        
        /**
         Returns a pointer to the parameter matching the name provided.
         @param[in] name variable name.
         @return a pointer to the parameter matching the name provided.
         */
        shared_ptr<param_> get_param(string name) const;
        
        void add_var(shared_ptr<param_> v, int nb = 1);/**< Inserts the variable in this function input list. nb represents the number of occurences v has. WARNING: Assumes that v has not been added previousely!*/
        
        void add_param(shared_ptr<param_> v, int nb = 1);/**< Inserts the parameter in this function input list. nb represents the number of occurences v has. WARNING: Assumes that v has not been added previousely!*/
        
        
        /**
         Reverse the convexity property of the current function
         */
        void reverse_convexity();
        
        
        
        
        void update_sign_add(const constant_& c);
        void update_sign_multiply(const constant_& c);
        
        void update_quad_convexity();
        
        void update_convexity_add(Convexity c){
            if(_all_convexity==linear_){
                _all_convexity = c;
            }
            else if(c!=linear_){
                if(_all_convexity!=c){
                    _all_convexity = undet_;
                }
            }
        }
        virtual void print(){};
        virtual void print(int prec){};
        virtual void print(size_t index, int prec = 10) {};
        virtual void print(size_t i, size_t j, int prec = 10) {};
        
        
        
        void print_symbolic(bool endline = true, bool display_input = true);
    };
    
    
    
    /*
     From Table 1.1 in http://www2.math.uni-wuppertal.de/~xsc/preprints/prep_01_4.pdf
     Replacing R* with +INF or -INF
     */
    template<class T,typename enable_if<is_arithmetic<T>::value>::type* = nullptr>
    T extended_plus(T x, T y){
        
        if(x==numeric_limits<T>::max() && y==numeric_limits<T>::lowest()){
            throw invalid_argument("In function extended_plus cannot add +inf to -inf");
        }
        if(x==numeric_limits<T>::lowest() && y==numeric_limits<T>::max()){
            throw invalid_argument("In function extended_plus cannot add -inf to +inf");
        }
        //+INF
        if(x==numeric_limits<T>::max() || y==numeric_limits<T>::max()){
            return numeric_limits<T>::max();
        }
        if(x==numeric_limits<T>::lowest() || y==numeric_limits<T>::lowest()){
            return numeric_limits<T>::lowest();
        }
        auto res = (x+y);
        if (res>numeric_limits<T>::max()){
            return numeric_limits<T>::max();
        }
        if (res<numeric_limits<T>::lowest()){
            return numeric_limits<T>::lowest();
        }
        return res;
    }
    
    /*
     From Table 1.2 in http://www2.math.uni-wuppertal.de/~xsc/preprints/prep_01_4.pdf
     Replacing R* with +INF or -INF
     */
    template<class T,typename enable_if<is_arithmetic<T>::value>::type* = nullptr>
    T extended_minus(T x, T y){
        
        if(x==numeric_limits<T>::max() && y==numeric_limits<T>::max()){
            return numeric_limits<T>::max();
            throw invalid_argument("In function extended_minus cannot substract +inf to +inf");
        }
        if(x==numeric_limits<T>::lowest() && y==numeric_limits<T>::lowest()){
            return numeric_limits<T>::lowest();
            throw invalid_argument("In function extended_minus cannot substract -inf to -inf");
        }
        //+INF
        if(x==numeric_limits<T>::max() || y==numeric_limits<T>::lowest()){
            return numeric_limits<T>::max();
        }
        if(x==numeric_limits<T>::lowest() || y==numeric_limits<T>::max()){
            return numeric_limits<T>::lowest();
        }
        auto res = (x-y);
        if (res>numeric_limits<T>::max()){
            return numeric_limits<T>::max();
        }
        if (res<numeric_limits<T>::lowest()){
            return numeric_limits<T>::lowest();
        }
        return res;
    }
    
    /*
     From Table 1.3 in http://www2.math.uni-wuppertal.de/~xsc/preprints/prep_01_4.pdf
     Replacing R* with +INF or -INF
     */
    template<class T,typename enable_if<is_arithmetic<T>::value>::type* = nullptr>
    T extended_mult(T x, T y){
        //+INF
        if(x==numeric_limits<T>::lowest() && y==numeric_limits<T>::lowest()){
            return numeric_limits<T>::max();
        }
        if(x==numeric_limits<T>::lowest() && y<0){
            return numeric_limits<T>::max();
        }
        if(y==numeric_limits<T>::lowest() && x<0){
            return numeric_limits<T>::max();
        }
        if(x==numeric_limits<T>::lowest() && y==0){
            return numeric_limits<T>::lowest();
        }
        if(y==numeric_limits<T>::lowest() && x==0){
            return numeric_limits<T>::lowest();
        }
        if(x==numeric_limits<T>::max() && y==numeric_limits<T>::max()){
            return numeric_limits<T>::max();
        }
        if(x==numeric_limits<T>::max() && y>=0){
            return numeric_limits<T>::max();
        }
        if(y==numeric_limits<T>::max() && x>=0){
            return numeric_limits<T>::max();
        }
        //-INF
        if(x==numeric_limits<T>::lowest() && y==numeric_limits<T>::max()){
            return numeric_limits<T>::lowest();
        }
        if(x==numeric_limits<T>::lowest() && y>=0){
            return numeric_limits<T>::lowest();
        }
        if(y==numeric_limits<T>::lowest() && x>=0){
            return numeric_limits<T>::lowest();
        }
        if(x==numeric_limits<T>::max() && y==numeric_limits<T>::lowest()){
            return numeric_limits<T>::lowest();
        }
        if(x==numeric_limits<T>::max() && y<0){
            return numeric_limits<T>::lowest();
        }
        if(y==numeric_limits<T>::max() && x<0){
            return numeric_limits<T>::lowest();
        }
        if(x==numeric_limits<T>::max() && y==0){
            return numeric_limits<T>::max();
        }
        if(y==numeric_limits<T>::max() && x==0){
            return numeric_limits<T>::max();
        }
        if(y==0 && x==0){
            return 0;
        }
        auto res = (x*y);
        if (res>numeric_limits<T>::max()){
            return numeric_limits<T>::max();
        }
        if (res<numeric_limits<T>::lowest()){
            return numeric_limits<T>::lowest();
        }
        return res;
    }
    
    template<class T,typename enable_if<is_same<T,Cpx>::value>::type* = nullptr>
    T extended_minus(T x, T y){
        T res;
        res.real(extended_minus(x.real(),y.real()));
        res.imag(extended_minus(x.imag(),y.imag()));
        return res;
    }
    
    template<class T,typename enable_if<is_same<T,Cpx>::value>::type* = nullptr>
    T extended_plus(T x, T y){
        T res;
        res.real(extended_plus(x.real(),y.real()));
        res.imag(extended_plus(x.imag(),y.imag()));
        return res;
    }
    
    template<class T,typename enable_if<is_same<T,Cpx>::value>::type* = nullptr>
    T extended_mult(T x, T y){
        T res;
        res.real(extended_minus(extended_mult(x.real(),y.real()),extended_mult(x.imag(),y.imag())));
        res.imag(extended_plus(extended_mult(x.real(),y.imag()),extended_mult(x.imag(),y.real())));
        return res;
    }
    
    template<class T1, class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) <= sizeof(T1)>::type* = nullptr>
    shared_ptr<pair<T1,T1>> get_product_range(shared_ptr<pair<T1,T1>> x, shared_ptr<pair<T2,T2>> y){
        shared_ptr<pair<T1,T1>> res = make_shared<pair<T1,T1>>();
        T1 x1 = x->first;
        T1 x2 = x->second;
        T1 y1 = y->first;
        T1 y2 = y->second;
        T1 min1 = gravity::min(extended_mult(x1,y1), extended_mult(x1,y2));
        T1 min2 = gravity::min(extended_mult(x2,y1), extended_mult(x2,y2));
        T1 max1 = gravity::max(extended_mult(x1,y1), extended_mult(x1,y2));
        T1 max2 = gravity::max(extended_mult(x2,y1), extended_mult(x2,y2));
        res->first = gravity::min(min1,min2);
        res->second = gravity::max(max1,max2);
        return res;
    }
    
    template<class T1, class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T1) < sizeof(T2)>::type* = nullptr>
    shared_ptr<pair<T2,T2>> get_product_range(shared_ptr<pair<T1,T1>> x, shared_ptr<pair<T2,T2>> y){
        shared_ptr<pair<T2,T2>> res = make_shared<pair<T2,T2>>();
        T2 x1 = x->first;
        T2 x2 = x->second;
        T2 y1 = y->first;
        T2 y2 = y->second;
        T2 min1 = gravity::min(extended_mult(x1,y1), extended_mult(x1,y2));
        T2 min2 = gravity::min(extended_mult(x2,y1), extended_mult(x2,y2));
        T2 max1 = gravity::max(extended_mult(x1,y1), extended_mult(x1,y2));
        T2 max2 = gravity::max(extended_mult(x2,y1), extended_mult(x2,y2));
        res->first = gravity::min(min1,min2);
        res->second = gravity::max(max1,max2);
        return res;
    }
    
    template<class T1, class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) <= sizeof(T1)>::type* = nullptr>
    shared_ptr<pair<T1,T1>> get_div_range(shared_ptr<pair<T1,T1>> range1, shared_ptr<pair<T2,T2>> range2){
        if(range2->first==numeric_limits<T2>::lowest() || range2->second==numeric_limits<T2>::max()
           || range1->first==numeric_limits<T1>::lowest()|| range1->second==numeric_limits<T1>::max()){
            shared_ptr<pair<T1,T1>> res = make_shared<pair<T1,T1>>();
            res->first =numeric_limits<T1>::lowest();
            res->second =numeric_limits<T1>::max();
            return res;
        }
        auto inv_range2 = make_shared<pair<T2,T2>>(*range2);
        inv_range2->first = 1./inv_range2->first;
        inv_range2->second = 1./inv_range2->second;
        return get_product_range<T1, T1, nullptr>(range1, inv_range2);
    }
    
    template<class T1, class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T1) < sizeof(T2)>::type* = nullptr>
    shared_ptr<pair<T2,T2>> get_div_range(shared_ptr<pair<T1,T1>> range1, shared_ptr<pair<T2,T2>> range2){
        shared_ptr<pair<T2,T2>> res = make_shared<pair<T2,T2>>();
        if(range2->first==numeric_limits<T2>::lowest() || range2->second==numeric_limits<T2>::max()
           || range1->first==numeric_limits<T1>::lowest()|| range1->second==numeric_limits<T1>::max()){
            res->first =numeric_limits<T2>::lowest();
            res->second =numeric_limits<T2>::max();
            return res;
        }
        auto inv_range2 = make_shared<pair<T2,T2>>(*range2);
        inv_range2->first = 1./inv_range2->first;
        inv_range2->second = 1./inv_range2->second;
        return get_product_range<T2, T2, nullptr>(range1, inv_range2);
    }
    
    template<class T1, class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) <= sizeof(T1)>::type* = nullptr>
    shared_ptr<pair<T1,T1>> get_plus_range(shared_ptr<pair<T1,T1>> x, shared_ptr<pair<T2,T2>> y){
        shared_ptr<pair<T1,T1>> res = make_shared<pair<T1,T1>>();
        T1 x1 = x->first;
        T1 x2 = x->second;
        T1 y1 = y->first;
        T1 y2 = y->second;
        res->first = extended_plus(x1,y1);
        res->second = extended_plus(x2,y2);
        return res;
    }
    
    template<class T1, class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T1) < sizeof(T2)>::type* = nullptr>
    shared_ptr<pair<T2,T2>> get_plus_range(shared_ptr<pair<T1,T1>> x, shared_ptr<pair<T2,T2>> y){
        shared_ptr<pair<T2,T2>> res = make_shared<pair<T2,T2>>();
        T2 x1 = x->first;
        T2 x2 = x->second;
        T2 y1 = y->first;
        T2 y2 = y->second;
        res->first = extended_plus(x1,y1);
        res->second = extended_plus(x2,y2);
        return res;
    }
    
    
    template<class T1, class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) <= sizeof(T1)>::type* = nullptr>
    shared_ptr<pair<T1,T1>> get_minus_range(shared_ptr<pair<T1,T1>> x, shared_ptr<pair<T2,T2>> y){
        shared_ptr<pair<T1,T1>> res = make_shared<pair<T1,T1>>();
        T1 x1 = x->first;
        T1 x2 = x->second;
        T1 y1 = y->first;
        T1 y2 = y->second;
        res->first = extended_minus(x1,y2);
        res->second = extended_minus(x2,y1);
        return res;
    }
    
    template<class T1, class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T1) < sizeof(T2)>::type* = nullptr>
    shared_ptr<pair<T2,T2>> get_minus_range(shared_ptr<pair<T1,T1>> x, shared_ptr<pair<T2,T2>> y){
        shared_ptr<pair<T2,T2>> res = make_shared<pair<T2,T2>>();
        T2 x1 = x->first;
        T2 x2 = x->second;
        T2 y1 = y->first;
        T2 y2 = y->second;
        res->first = extended_minus(x1,y2);
        res->second = extended_minus(x2,y1);
        return res;
    }
    
    
    template<typename type>
    class func: public func_{
    private:
        /** Computes and stores the derivative of f with respect to variable v. Returns a pointer to the stored function. */
        shared_ptr<func> compute_derivative(const param_& v){
            auto vid = v._name;
            if(_dfdx->count(vid)!=0){
                return _dfdx->at(vid);
            }
            auto df = make_shared<func>(get_derivative(v));
            if(_is_vector){
                df->_is_vector=true;//TODO check this.
            }
            df->_evaluated = false;
            df->allocate_mem();
            (*_dfdx)[vid] = df;
            DebugOff( "First derivative with respect to " << v.get_name(false,false) << " = " << df->to_str() << endl);
            return df;
        }
        
    public:
        shared_ptr<expr<type>>                                  _expr = nullptr; /**< Nonlinear part of the function, this points to the root node in _DAG */
        
        shared_ptr<map<string,shared_ptr<func>>>                _dfdx = nullptr;/**< A map storing the derivatives indexed by variables' names */
        shared_ptr<vector<type>>                                _val = nullptr; /**< vector of values **/
        shared_ptr<pair<type,type>>                             _range = nullptr; /**< (Min,Max) values in vals **/
        shared_ptr<vector<pair<type,type>>>                     _all_range = nullptr; /**< Vector of (Min,Max) values for each instance of this func **/
        
        void update_range(){
            _range = make_shared<pair<type,type>>(make_pair<>(zero<type>().eval(), zero<type>().eval()));
        }
        
        
        func(){
            update_type();
            update_range();
            _cst = make_shared<constant<type>>();
            _lterms = make_shared<map<string, lterm>>();
            _qterms = make_shared<map<string, qterm>>();
            _pterms = make_shared<map<string, pterm>>();
            _vars = make_shared<map<string, pair<shared_ptr<param_>, unsigned>>>();
            _params = make_shared<map<string, pair<shared_ptr<param_>, unsigned>>>();
            _dfdx = make_shared<map<string,shared_ptr<func<type>>>>();
            _val = make_shared<vector<type>>();
        };
        
        
        void update_type() {
            _type = func_c;
            if(typeid(type)==typeid(bool)) {
                _return_type = binary_;
                return;
            }
            if(typeid(type)==typeid(short)) {
                _return_type = short_;
                return;
            }
            if(typeid(type)==typeid(int)) {
                _return_type = integer_;
                return;
            }
            if(typeid(type)==typeid(float)) {
                _return_type = float_;
                return;
            }
            if(typeid(type)==typeid(double)) {
                _return_type = double_;
                return;
            }
            if(typeid(type)==typeid(long double)) {
                _return_type = long_;
                return;
            }
            if(typeid(type)==typeid(Cpx)) {
                _return_type = complex_;
                return;
            }
            throw invalid_argument("Unsupported numerical function type");
        }
        
        
        bool is_constant() const{
            return (_vars->empty());
        }
        
        
        
        /**
         Index the function and its variables/parameters using the indices in ids
         @param[in] ids indices
         @return current function
         */
        func& in(const indices& ids){//TODO assert each var and param that is not transposed has same dim as ids
            //            index_in(ids);
            //            _nb_vars = 0;
            //            string key;
            //            auto iter = _vars->begin();
            //            while (iter!=_vars->end()) {
            //                auto pair = (*iter++);
            //                auto v = pair.second.first;
            //                if(!v->is_indexed()){
            //                    v->index_in(ids);
            //                }
            //                if (!v->_is_vector) {// i.e., it is not transposed
            //                    _nb_vars++;
            //                }
            //                else {
            //                    _nb_vars += v->get_dim();
            //                }
            //            }
            //            iter = _params->begin();
            //            while (iter!=_params->end()) {
            //                auto pair = (*iter++);
            //                auto p = pair.second.first;
            //                if(!p->is_indexed()){
            //                    p->index_in(ids);
            //                }
            //            }
            _indices = make_shared<indices>(ids);
            _dim[0] = ids.size();
            //            if(_expr){
            //                _expr->in(ids);
            //            }
            return *this;
        }
        
        
        /** Operators */
        bool operator==(const func& p) const {
            if (_type!=p._type || _return_type!=p._return_type || _dim[0]!=p._dim[0] || _dim[1]!=p._dim[1] || _to_str!=p._to_str) return false;
            if(_indices==p._indices) return true; /* accounts for both being nullptr */
            if((_indices && !p._indices) || (p._indices && !_indices) || (*_indices != *p._indices)) return false;
            return true;
        }
        
        func vec() const{
            auto f(*this);
            if(func_is_param()){
                auto vi = f._vars->begin()->second.first;
                if(!vi->_is_vector){
                    vi->_is_vector = true;
                    vi->_name = "["+vi->_name+"]";
                }
            }
            f._is_vector = true;
            return f;
        }
        
        shared_ptr<expr<type>> get_expr() const{
            return _expr;
        }
        
        void transpose(){
            _is_transposed = !_is_transposed;
            _is_vector = true;
            auto temp = _dim[0];
            _dim[0] = _dim[1];
            _dim[1] = temp;
            for (auto &vp: *_vars) {
                if(vp.second.first->_is_transposed){
                    vp.second.first->transpose();
                }
            }
            for (auto &pp: *_params) {
                if(pp.second.first->_is_transposed){
                    pp.second.first->transpose();
                }
            }
        }
        
        func get_derivative(shared_ptr<constant_> ex, const param_& v) const{
            auto name = v.get_name(false,false);
            if(ex->is_var()){
                auto vv = static_pointer_cast<param_>(ex);
                if(vv->get_name(false,false)==name){
                    return unit<type>();
                }
            }
            else if(ex->is_function()){
                auto f = static_pointer_cast<func>(ex);
                return f->get_derivative(v);
            }
            else if(ex->is_uexpr()){
                func son;
                auto uexp = static_pointer_cast<uexpr<type>>(ex);
                if (uexp->_son->is_function()) {
                    auto f = static_pointer_cast<func>(uexp->_son);
                    son = (*f);
                }
                else if(uexp->_son->is_var()) {
                    auto vv = static_pointer_cast<param_>(uexp->_son);
                    if(vv->get_name(false,false)==name){
                        son.insert(v);
                    }
                    else{
                        return func();
                    }
                    
                }
                else {
                    return func();
                }
                // f(g(x))' = f'(g(x))*g'(x).
                switch (uexp->_otype) {
                    case cos_:
                        return uexp->_coef*-1.*get_derivative(uexp->_son,v)*sin(son);
                        break;
                    case sin_:
                        return uexp->_coef*get_derivative(uexp->_son,v)*cos(son);
                        break;
                    case sqrt_:
                        return uexp->_coef*get_derivative(uexp->_son,v)/(2.*sqrt(son));
                        break;
                    case exp_:{
                        return uexp->_coef*get_derivative(uexp->_son,v)*exp(son);
                        break;
                    }
                    case log_:
                        return uexp->_coef*get_derivative(uexp->_son,v)/son;
                        break;
                    case relu_:
                        return uexp->_coef*unit_step(son)*get_derivative(uexp->_son,v);
                        break;
                    case unit_step_:
                        return zero<type>();
                        break;
                    default:
                        throw invalid_argument("Unsupported unary operation");
                        break;
                }
            }
            else if(ex->is_bexpr()){
                func lson, rson;
                auto bexp = static_pointer_cast<bexpr<type>>(ex);
                if (bexp->_lson->is_function()) {
                    auto f = static_pointer_cast<func>(bexp->_lson);
                    lson = *f;
                }
                else if(bexp->_lson->is_var()) {
                    auto vv = static_pointer_cast<param_>(bexp->_lson);
                    if(vv->get_name(false,false)==name){
                        lson.insert(v);
                    }
                }
                else {
                    return func();
                }
                if (bexp->_rson->is_function()) {
                    auto f = static_pointer_cast<func>(bexp->_rson);
                    rson = *f;
                }
                else if(bexp->_rson->is_var()) {
                    auto vv = static_pointer_cast<param_>(bexp->_rson);
                    if(vv->get_name(false,false)==name){
                        rson.insert(v);
                    }
                }
                else {
                    return func();
                }
                switch (bexp->_otype) {
                    case plus_:
                        return bexp->_coef*(lson.get_derivative(v) + rson.get_derivative(v));
                        break;
                    case minus_:
                        return bexp->_coef*(lson.get_derivative(v) - rson.get_derivative(v));
                        break;
                    case product_:{
                        return bexp->_coef*(lson.get_derivative(v)*(rson) + (lson)*rson.get_derivative(v));
                        // f'g + fg'
                        break;
                    }
                    case div_:
                        return bexp->_coef*((lson.get_derivative(v)*(rson) - (lson)*rson.get_derivative(v))/(rson*rson));
                        // (f'g - fg')/g^2
                        break;
                    case power_:
                        if (!rson.is_number()) {
                            throw invalid_argument("Function in exponent not supported yet.\n");
                        }
                        //                auto exponent = poly_eval(_rson);
                        //                return (exponent*get_poly_derivative(_lson,v)*power(*_lson, exponent-1));// nf'f^n-1
                        break;
                    default:
                        throw invalid_argument("unsupported operation");
                        break;
                }
            }
            return func();
        }
        
        shared_ptr<func> get_stored_derivative(const string& vid) const{ /**< Returns the stored derivative with respect to variable v. */
            auto it = _dfdx->find(vid);
            if (it!=_dfdx->end()) {
                return it->second;
            }
            else {
                throw invalid_argument("No derivatives stored!\n");
            }
        }
        
        //        func_ get_derivative(const param_& v) const; /**< Computes and returns the derivative with respect to variable v. */
        
        //        func_ get_dfdx(const param_& v); /**< Computes all derivatives and returns a copy of the derivative with respect to variable v. */
        
        
        template<typename T=type,
        typename std::enable_if<is_arithmetic<T>::value>::type* = nullptr>
        void init_range() {
            _range = make_shared<pair<type,type>>(make_pair<>(numeric_limits<type>::max(), numeric_limits<type>::lowest()));
        }
        
        
        template<class T=type, class = typename enable_if<is_same<T, Cpx>::value>::type>
        void init_range() {
            _range = make_shared<pair<type,type>>(make_pair<>(Cpx(numeric_limits<double>::max(), numeric_limits<double>::max()), Cpx(numeric_limits<double>::lowest(), numeric_limits<double>::lowest())));
        }
        
        /**
         Recompute range based on stored values.
         */
        void reset_range(){
            init_range();
            if(is_matrix_indexed()){
                for(auto i = 0; i<_indices->_ids->size();i++){
                    for(auto j = 0; j<_indices->_ids->at(i).size();j++){
                        auto idx = _indices->_ids->at(i).at(j);
                        auto v = _val->at(idx);
                        if(_range->first > v){
                            _range->first = v;
                        }
                        if(_range->second  < v){
                            _range->second = v;
                        }
                    }
                }
            }
            else if(is_indexed()){
                for(auto i = 0; i<_indices->_ids->at(0).size();i++){
                    auto idx = _indices->_ids->at(0).at(i);
                    auto v = _val->at(idx);
                    if(_range->first > v){
                        _range->first = v;
                    }
                    if(_range->second  < v){
                        _range->second = v;
                    }
                }
            }
            else {
                for (auto v:*_val) {
                    if(_range->first > v){
                        _range->first = v;
                    }
                    if(_range->second  < v){
                        _range->second = v;
                    }
                }
            }
        }
        
        func<type> get_outer_app(){ /**< Returns an outer-approximation of the function using the current value of the variables **/
            func<type> res; // res = gradf(x*)*(x-x*) + f(x*)
            param<type> f_xstar("f_xstar");
            f_xstar = *this;
            for(auto &it: *_vars){
                auto v = it.second.first;
                param<type> xstar("xstar_"+v->_name);
                xstar.in(*v->_indices);
                xstar.copy_vals(v);
                param<type> df_xstar("df_xstar"+v->_name);
                auto df = *compute_derivative(*v);
                //df.uneval();
                df.eval_all();
                df_xstar = df;
                auto ids = v->_indices->get_unique_keys();
                df_xstar._indices = make_shared<indices>(ids);
                res.insert(true, df_xstar, *v);
                res -= df_xstar*xstar;
            }
            res += f_xstar;
            res._indices = this->_indices;
            merge_vars(res);
            return res;
        }
        
        func<type> get_outer_app_insti(size_t nb_inst){ /**< Returns an outer-approximation of the function using the current value of the variables **/
            func<type> res; // res = gradf(x*)*(x-x*) + f(x*)
            double f_xstar, xv, dfv;
            vector<double> xcurrent, dfvector;
            uneval();
            f_xstar=eval(nb_inst);
            DebugOn("F_xstar in func.h\t"<<f_xstar<<endl);
            for(auto &it: *_vars){
                auto v = it.second.first;
                size_t posv=v->get_id_inst(nb_inst);
                v->get_double_val(posv, xv);
                xcurrent.push_back(xv);
                auto df = *compute_derivative(*v);
                df.uneval();
                dfv=df.eval(nb_inst);
                dfvector.push_back(dfv);
                indices ids("ids");
                auto key=v->_indices->_keys;
                ids.add((*key)[posv]);
                
                
                
                param<type> df_xstar("df_xstar"+v->_name);
                df_xstar.in(ids);
                df_xstar.set_val(dfv);
                
                switch (v->get_intype()) {
                    case binary_:
                        res += df_xstar*(*static_pointer_cast<param<bool>>(v)).in(ids);
                        break;
                    case short_:
                        res += df_xstar*(*static_pointer_cast<param<short>>(v)).in(ids);
                        break;
                    case integer_:
                        res += df_xstar*(*static_pointer_cast<param<int>>(v)).in(ids);
                        break;
                    case float_:
                        res += df_xstar*(*static_pointer_cast<param<float>>(v)).in(ids);
                        break;
                    case double_:
                        res += df_xstar*(*static_pointer_cast<param<double>>(v)).in(ids);
                        break;
//                    case long_:
//                        res += df_xstar*(*static_pointer_cast<param<long double>>(v)).in(ids);
//                        break;
//                    case complex_:
//                        res += df_xstar*(*static_pointer_cast<param<Cpx>>(v)).in(ids);
//                        break;
                    default:
                        break;
                }
//
//                auto v1=v->pcopy();
//                v1->_indices=make_shared<gravity::indices>(ids);
//                res.insert(true, df_xstar, *v1);
                
                
                //Alterntaively tried the below as well, both forms give correct functional form of OA cut in the absence of merge_vars
                
                //                auto indcopy=v->_indices;
                //                v->_indices=make_shared<gravity::indices>(ids);
                //                res.insert(true, df_xstar, *v);
                //                merge_vars(res);
                //                v->_indices=indcopy;
                res -= dfv*xv;
            }
            
            
            res += f_xstar;
            
            res.print();
            
            
            indices res_ind("res_ind");
            res_ind.add(to_string(nb_inst));
            res._indices=make_shared<gravity::indices>(res_ind);
            
            res.eval_all();
            res.uneval();
//            res.merge_vars(*this);
            DebugOn("Eval of OA_cut in get_outer_app_insti\t"<<res.eval(0)<<endl);
            
           
            
            
            //res.print();
            //            DebugOn("Xcurrent from get_outer_app_insti"<<endl);
            //            for(auto i=0;i<xcurrent.size();i++)
            //                DebugOn(xcurrent[i]<<"\t");
            //            DebugOn(endl);
            //            DebugOn("DF at Xcurrent from get_outer_app_insti"<<endl);
            //            for(auto i=0;i<xcurrent.size();i++)
            //                DebugOn(dfvector[i]<<"\t");
            //            DebugOn(endl);
            return res;
        }
        
        
        
        double l2norm(vector<double> x)
        {
            double res=0;
            for(auto i=0;i<x.size();i++)
            {
                res+=x[i]*x[i];
            }
            res=sqrt(res);
            return(res);
        }
        
        //x_start is an interior point and x_end is an outer point.
        //Interior and outer point clasification depends on constraint type (\geq 0 or \leq 0) as input by ctype
        pair<vector<double>,bool> linesearchbinary(vector<double> x_start, vector<double> x_end, size_t nb_inst, ConstraintType ctype)
        {
            pair<vector<double>,bool> res;
            const double int_tol=1e-6, zero_tol=1e-6;
            const int max_iter=1000;
            vector<double> x_f, x_t, xcurrent, interval, mid;
            double  f_a,f_b,f_f, f_t, f_mid, interval_norm, xv;
            bool solution_found=false;
            int iter=0;
            for(auto i=0;i<x_start.size();i++)
            {
                interval.push_back(x_end[i]-x_start[i]);
                mid.push_back((x_end[i]+x_start[i])*0.5);
            }
            int counter=0;
            for(auto &it: *_vars)
            {
                auto v = it.second.first;
                size_t posv=v->get_id_inst(nb_inst);
                v->get_double_val(posv, xv);
                xcurrent.push_back(xv);
                v->set_double_val(posv, x_start[counter++]);
            }
            uneval();
            f_a=eval(nb_inst);
            
            counter=0;
            for(auto &it: *_vars)
            {
                auto v = it.second.first;
                size_t posv=v->get_id_inst(nb_inst);
                v->set_double_val(posv, x_end[counter++]);
            }
            uneval();
            f_b=eval(nb_inst);
            if(ctype==leq)
            {
                f_f=f_a;
                f_t=f_b;
                x_f=x_start;
                x_t=x_end;
            }
            else
            {
                f_f=f_b;
                f_t=f_a;
                x_f=x_end;
                x_t=x_start;
            }
            interval_norm=l2norm(interval);
            
            if(f_f<=0 && f_t>=0 )
            {
                while(interval_norm>int_tol && iter<=max_iter)
                {
                    for(auto i=0;i<x_start.size();i++)
                    {
                        mid[i]=(x_f[i]+x_t[i])*0.5;
                    }
                    counter=0;
                    for(auto &it: *_vars)
                    {
                        auto v = it.second.first;
                        size_t posv=v->get_id_inst(nb_inst);
                        v->set_double_val(posv, mid[counter++]);
                    }
                    uneval();
                    f_mid=eval(nb_inst);
                    // DebugOn("F_mid "<< f_mid<<endl);
                    //DebugOn("xf\t xt\t xmid"<<endl);
                    for(auto i=0;i<x_start.size();i++)
                    {
                        //  DebugOn(x_f[i]<<"\t"<<x_t[i]<<"\t"<<mid[i]<<endl);
                        
                    }
                    
                    if(f_mid>=zero_tol && f_mid<=f_t)
                    {
                        x_t=mid;
                    }
                    else if(f_mid<=zero_tol*(-1) && f_mid>=f_f)
                    {
                        x_f=mid;
                    }
                    else
                    {
                        //DebugOn("Reached answer"<<endl);
                        solution_found=true;
                        break;
                    }
                    for(auto i=0;i<x_start.size();i++)
                    {
                        interval[i]=x_t[i]-x_f[i];
                    }
                    
                    interval_norm=l2norm(interval);
                    iter++;
                }
                //
                //            DebugOn("F_mid "<<f_mid<<endl);
                //            DebugOn("Interval_Norm "<<interval_norm<<endl);
                //            DebugOn("Iter "<<iter<<endl);
            }
            
            res.first=mid;
            res.second=solution_found;
            if(res.second)
            {
                // DebugOn("Solution to line search found"<<endl);
                for(auto i=0;i<res.first.size();i++)
                    //DebugOn(res.first[i]<<endl);
                    counter=0;
                for(auto &it: *_vars)
                {
                    auto v = it.second.first;
                    size_t posv=v->get_id_inst(nb_inst);
                    v->set_double_val(posv, mid[counter++]);
                }
                uneval();
                //  DebugOn("Function value at pos "<<nb_inst<<" at solution of line search "<<eval(nb_inst));
                
            }
            counter=0;
            for(auto &it: *_vars)
            {
                auto v = it.second.first;
                size_t posv=v->get_id_inst(nb_inst);
                v->set_double_val(posv, xcurrent[counter++]);
            }
            return res;
        }
        
        
        
        
        /** Finds a vector of outer points perturbing along each direction */
        //Algorithm finds an outer point for each index of each variable if available
        //First, if available,the outer point is at least at a distance perturb_distance greater than original value of variable
        //Else, if available, the algorithm returns any outer point produced by perturbing variable
        //Else, the algorithm does not return anything
        //Interior and outer point clasification depends on constraint type (\geq 0 or \leq 0) as input by ctype
        vector<vector<double> > get_outer_point(size_t nb_inst, ConstraintType ctype)
        {
            vector<vector<double> > res(_nb_vars);
            vector<double> xcurrent, ub_v, lb_v;
            const int max_iter=1000;
            const double step_tol=1e-6, step_init=1e-3, perturb_dist=1e-3, zero_tol=1e-6;
            double step, f_start, xv=0,xv_p=0,f,ub,lb, fnew, dfdv;
            int count=0, iter, sign, iter_dir;
            bool perturb=true, dir;
            f_start=eval(nb_inst);
            for(auto &it: *_vars)
            {
                auto v = it.second.first;
                size_t posv=v->get_id_inst(nb_inst);
                v->get_double_val(posv, xv);
                xcurrent.push_back(xv);
            }
            int res_count=0;
            
            
            //No backtracking
            
            
            //Once feasible direction is found algorithm does not reverse direction. So shall work from any current point only for convex function and will work to identify one outer point, not necessarily at greater than perturb_dist from an active point for any nonconvex function
            
            //Perturb so that distance between new point and current point is greater than perturb dist
            for(auto &it: *_vars)
            {
                perturb=true;
                dir=false;
                iter=0;
                auto v = it.second.first;
                size_t posv=v->get_id_inst(nb_inst);
                xv=xcurrent[count];
                step=step_tol;
                iter_dir=0;
                ub=v->get_double_ub(posv);
                lb=v->get_double_lb(posv);
                auto df = *compute_derivative(*v);
                df.uneval();
                dfdv=df.eval(nb_inst);
                // if interval zero do not perturb, perturb=false, else if x at upper bound (within perturb_dist) do not step out but set sign=-1,else if x at lower bound set sign=1, else go to white loop
                if((ub-lb)<=perturb_dist)
                {
                    dir=false;
                }
                else if((xv-lb)<=perturb_dist)
                {
                    sign=1;
                    dir=true;
                    
                    if(ctype==leq && dfdv<0)
                    {
                        dir=false;
                    }
                    else if(ctype==geq && dfdv>0)
                    {
                        dir=false;
                    }
                    
                }
                else if((ub-xv)<=perturb_dist)
                {
                    sign=-1;
                    dir=true;
                    if(ctype==leq && dfdv<0)
                    {
                        dir=false;
                    }
                    else if(ctype==geq && dfdv>0)
                    {
                        dir=false;
                    }
                    
                }
                else
                {
                    if(ctype==leq)
                    {
                        if(dfdv>0)
                        {
                            dir=true;
                            sign=1;
                            
                        }
                        else if(dfdv<0)
                        {
                            dir=true;
                            sign=-1;
                        }
                    }
                    else if(ctype==geq)
                    {
                        if(dfdv<0)
                        {
                            dir=true;
                            sign=1;
                            
                        }
                        else if(dfdv>0)
                        {
                            dir=true;
                            sign=-1;
                        }
                        
                    }
                }
                if(dir)
                {
                    step=dfdv;
                    perturb=false;
                    f=f_start;
                    while(!perturb && iter<=max_iter)
                    {
                        if(sign==1)
                        {
                            xv=std::min(xv*(1+step), ub);
                            v->set_double_val(posv, xv);
                        }
                        else
                        {
                            xv=std::max(xv*(1-step), lb);
                            v->set_double_val(posv, xv);
                        }
                        uneval();
                        fnew=eval(nb_inst);
                        if(ctype==leq)
                        {
                            if(fnew>zero_tol && abs(xv-xcurrent[count])>=perturb_dist)
                            {
                                perturb=true;
                                xv_p=xv;
                                break;
                            }
                            else if(fnew>f)
                            {
                                f=fnew;
                                xv_p=xv;
                            }
                            else if(fnew<=f)
                            {
                                perturb=false;
                                break;
                            }
                        }
                        if(ctype==geq)
                        {
                            if(fnew<(zero_tol*(-1)) && abs(xv-xcurrent[count])>=perturb_dist)
                            {
                                perturb=true;
                                xv_p=xv;
                                break;
                            }
                            else if(fnew<f)
                            {
                                f=fnew;
                                xv_p=xv;
                            }
                            else if(fnew>=f)
                            {
                                perturb=false;
                                break;
                            }
                        }
                        iter++;
                    }
                    if(perturb==true || (f>zero_tol && ctype==leq) || (f<(zero_tol*(-1)) && ctype==geq) )
                    {
                        for(auto i=0;i<count;i++)
                            res[res_count].push_back(xcurrent[i]);
                        res[res_count].push_back(xv_p);
                        for(auto i=count+1;i<_nb_vars;i++)
                            res[res_count].push_back(xcurrent[i]);
                        res_count++;
                        
                        for(auto &it: *_vars)
                        {
                            auto v = it.second.first;
                            size_t posv=v->get_id_inst(nb_inst);
                            v->get_double_val(posv, xv);
                            //DebugOn("Xvalues of Outer point\t"<<xv<<endl);
                        }
                        uneval();
                        //DebugOn("fvalue at pos "<<nb_inst<<" at the outer point\t"<<eval(nb_inst)<<endl);
                        
                    }
                }
                v->set_double_val(posv, xcurrent[count]);
                f=f_start;
                count++;
            }
            return(res);
        }
        
     
        pair<bool, vector<double>> newton_raphson(vector<double> x, string vname, size_t nb_inst, ConstraintType ctype)
        {
            pair<bool, vector<double>> res;
            vector<double> xcurrent, xk, xsolution;
            double xvk, xvk1, fk, dfdvk, xv,ub,lb;
            const int max_iter=100;
            const double active_tol=1e-6,zero_tol=1e-8;

            int counter=0,iter=0;
            bool solution=false;
            for(auto &it: *_vars)
            {
                auto v = it.second.first;
                size_t posv=v->get_id_inst(nb_inst);
                v->get_double_val(posv, xv);
               // DebugOn("Name\t"<<name<<"Posv\t"<<posv<<"XV\t"<<xv<<endl);
                xcurrent.push_back(xv);
                v->set_double_val(posv, x[counter++]);
                
            }
            xk=x;
            
            
            fk=eval(nb_inst);
            auto vk=_vars->at(vname).first;
            size_t posvk=vk->get_id_inst(nb_inst);
            vk->get_double_val(posvk, xvk);
            auto df = *compute_derivative(*vk);
            df.uneval();
            dfdvk=df.eval(nb_inst);
            ub=vk->get_double_ub(posvk);
            lb=vk->get_double_lb(posvk);

            
            
             while(iter<=max_iter && !solution)
                       {
            vk->set_double_val(posvk, xvk);
                           uneval();
                           fk=eval(nb_inst);
                           if(abs(fk)<=active_tol)
                           {
                               solution=true;
                           break;
                       }
                                               auto df = *compute_derivative(*vk);
                                               df.uneval();
                                               dfdvk=df.eval(nb_inst);
                           if(abs(dfdvk)>=zero_tol)
                           xvk1=xvk-fk/dfdvk;
                           else
                               break;
                        
                           if(abs(xvk1-xvk)<=zero_tol)
                               break;
                           if(xvk1>=ub)
                               xvk=ub;
                           else if (xvk1<=lb)
                               xvk=lb;
                           else
                               xvk=xvk1;
                       
                           
                           iter++;
            
        }
            
   

            for(auto &it: *_vars)
            {
  
                
                auto v = it.second.first;
                size_t posv=v->get_id_inst(nb_inst);
                
                v->get_double_val(posv, xv);
                
                // DebugOn("Name\t"<<name<<"Posv\t"<<posv<<"XV\t"<<xv<<endl);
                xsolution.push_back(xv);
                
                //Reset xcurrent
                
            }
            res.first=solution;
            res.second=xsolution;
            counter=0;
            for(auto &it: *_vars)
            {
                auto v = it.second.first;
                size_t posv=v->get_id_inst(nb_inst);
                v->set_double_val(posv, xcurrent[counter++]);
                
            }
            return res;
        }

        vector<vector<double> > get_active_point(size_t nb_inst, ConstraintType ctype)
        {
            vector<vector<double> > res(_nb_vars);
            vector<double> xcurrent;
            int res_count=0;
            double xv;
            
            for(auto &it: *_vars)
            {
                auto v = it.second.first;
                size_t posv=v->get_id_inst(nb_inst);
                v->get_double_val(posv, xv);
                // DebugOn("Name\t"<<name<<"Posv\t"<<posv<<"XV\t"<<xv<<endl);
                xcurrent.push_back(xv);
                
            }
            
            for(auto &it: *_vars)
            {
                auto vname=it.first;
                
                auto res_nr=newton_raphson(xcurrent, vname, nb_inst, ctype);
                
                if(res_nr.first==true)
                {
                    
                    for(auto i=0;i<xcurrent.size();i++)
                        res[res_count].push_back(res_nr.second[i]);
                res_count++;
                }
            
        }
            return res;
        }
        
        /** Computes and stores the derivative of f with respect to all variables. */
        void compute_derivatives(){
            size_t vid = 0, vjd = 0;
            param_* vi;
            param_* vj;
            DebugOff( "Computing derivatives for " << to_str() << endl);
            for (auto &vp: *_vars) {
                vi = vp.second.first.get();
                vid = vi->get_id();
                auto vi_name = vp.first;
                auto df = compute_derivative(*vi);
                //            if (is_nonlinear()) {
                DebugOff( "First derivative with respect to " << vp.first << " = " << df->to_str() << endl);
                //                df->print();
                //            }
                for (auto &vp2: *df->_vars) {
                    vj = vp2.second.first.get();
                    vjd = vj->get_id();
                    auto vj_name = vp2.first;
                    if (vi_name.compare(vj_name) <= 0) { //only store lower left part of hessian matrix since it is symmetric.
                        auto d2f = df->compute_derivative(*vj);
                        DebugOff( "Second derivative with respect to " << vp.first << " and " << vp2.first << " = " << d2f->to_str() << endl);
                        //                        d2f->print();
                    }
                }
                
            }
        }
        
        shared_ptr<map<string,shared_ptr<func>>> get_dfdx() const{
            return _dfdx;
        };
        
        void set_first_derivative(const param_& v, func&& f){
            DebugOff(f.to_str()<<endl);
            (*_dfdx)[v._name] = make_shared<func>(move(f));
        }
        
        void set_second_derivative(const param_& v1, const param_& v2, func&& f){
            DebugOff(f.to_str()<<endl);
            (*_dfdx)[v1._name]->_dfdx->insert(make_pair<>(v2._name, make_shared<func>(move(f))));
        }
        
        pair<type,type> get_range(shared_ptr<list<pair<shared_ptr<param_>, int>>> pterm) const{
            pair<type,type> res = {unit<type>().eval(),unit<type>().eval()};
            auto it = pterm->begin();
            while(it!=pterm->end()){
                auto it2 = next(it);
                auto range1 = get_range(*it);
                if(it2!=pterm->end()){
                    auto range2 = get_range(*it2);
                    auto prod = get_product_range(make_shared<pair<type,type>>(range1), make_shared<pair<type,type>>(range2));
                    res = *get_product_range(make_shared<pair<type,type>>(res),prod);
                    advance(it,2);
                }
                else {
                    res = *get_product_range(make_shared<pair<type,type>>(res),make_shared<pair<type,type>>(range1));
                    it++;
                }
            }
            return res;
        }
        
        pair<type,type> get_range(const pair<shared_ptr<param_>, int>& p) const{
            auto range = get_range(p.first);
            pair<type,type> res;
            res.first=pow(range.first,p.second);
            res.second=pow(range.second,p.second);
            if(p.second%2==0){
                res.first=zero<type>().eval();
                if(p.first->is_positive()|| p.first->is_negative()){
                    res.first=pow(range.first,p.second);
                }
            }
            return res;
        }
        
        
        template<typename T=type,
        typename std::enable_if<is_arithmetic<T>::value>::type* = nullptr>
        pair<type,type> get_range(shared_ptr<param_> p) const{
            switch (p->get_intype()) {
                case binary_:
                    return *static_pointer_cast<param<bool>>(p)->_range;
                    break;
                case short_:
                    return *static_pointer_cast<param<short>>(p)->_range;
                    break;
                case integer_:
                    return *static_pointer_cast<param<int>>(p)->_range;
                    break;
                case float_:
                    return *static_pointer_cast<param<float>>(p)->_range;
                    break;
                case double_:
                    return *((param<double>*)(p.get()))->_range;
                    break;
                case long_:
                    return *static_pointer_cast<param<long double>>(p)->_range;
                    break;
                default:
                    break;
            }
            return pair<type,type>();
        };
        
        template<class T=type, typename enable_if<is_same<T, Cpx>::value>::type* = nullptr>
        pair<type,type> get_range(shared_ptr<param_> p) const{
            switch (p->get_intype()) {
                case binary_:
                    return *static_pointer_cast<param<bool>>(p)->_range;
                    break;
                case short_:
                    return *static_pointer_cast<param<short>>(p)->_range;
                    break;
                case integer_:
                    return *static_pointer_cast<param<int>>(p)->_range;
                    break;
                case float_:
                    return *static_pointer_cast<param<float>>(p)->_range;
                    break;
                case double_:
                    return *((param<double>*)(p.get()))->_range;
                    break;
                case long_:
                    return *static_pointer_cast<param<long double>>(p)->_range;
                    break;
                case complex_:
                    return *static_pointer_cast<param<Cpx>>(p)->_range;
                    break;
                default:
                    break;
            }
            return pair<type,type>();
        };
        
        
        func get_derivative(const param_& v) const{
            func res;
            shared_ptr<pair<type,type>> term_range;
            if(!has_var(v)){
                return res;
            }
            if(v._is_vector){
                res._is_vector=true;
            }
            auto name = v.get_name(false,false);
            for (auto &lt: *_qterms) {
                if (lt.second._p->first->get_name(false,false) == name) {
                    res.set_max_dim(v);
                    auto coef = lt.second._coef->copy();
                    if (coef->_is_transposed) {
                        coef->transpose();//TODO is this needed?
                        coef->_is_vector = false;
                    }
                    if (coef->is_function()) {
                        auto f_cst = *((func<type>*)(coef.get()));
                        auto var_range = make_shared<pair<type,type>>(get_range(lt.second._p->second));
                        term_range = get_product_range(f_cst._range,var_range);
                        res.insert(lt.second._sign, f_cst, *lt.second._p->second);
                    }
                    else if(coef->is_param()) {
                        auto p_cst = *((param<type>*)(coef.get()));
                        auto var_range = make_shared<pair<type,type>>(get_range(lt.second._p->second));
                        term_range = get_product_range(p_cst._range,var_range);
                        res.insert(lt.second._sign, p_cst, *lt.second._p->second);
                        if(p_cst.is_indexed()){
                            res._indices = p_cst._indices;
                        }
                    }
                    else if(coef->is_number()) {
                        auto p_cst = *((constant<type>*)(coef.get()));
                        auto var_range = make_shared<pair<type,type>>(get_range(lt.second._p->second));
                        term_range = get_product_range(make_shared<pair<type,type>>(p_cst.eval(),p_cst.eval()),var_range);
                        res.insert(lt.second._sign, p_cst, *lt.second._p->second);
                    }
                    if (lt.second._coef->_is_vector || lt.second._coef->is_matrix()) {
                        res._is_vector = true;
                    }
                    if(lt.second._sign){
                        res._range = get_plus_range(res._range, term_range);
                    }
                    else {
                        res._range = get_minus_range(res._range, term_range);
                    }
                }
                if (lt.second._p->second->get_name(false,false) == name) {
                    res.set_max_dim(v);
                    auto coef = lt.second._coef->copy();
                    if (coef->_is_transposed) {
                        coef->transpose();//TODO is this needed?
                    }
                    if (coef->is_function()) {
                        auto f_cst = *((func<type>*)(coef.get()));
                        auto var_range = make_shared<pair<type,type>>(get_range(lt.second._p->second));
                        term_range = get_product_range(f_cst._range,var_range);
                        res.insert(lt.second._sign, f_cst, *lt.second._p->first);
                    }
                    else if(coef->is_param()) {
                        auto p_cst = *((param<type>*)(coef.get()));
                        auto var_range = make_shared<pair<type,type>>(get_range(lt.second._p->second));
                        term_range = get_product_range(p_cst._range,var_range);
                        res.insert(lt.second._sign, p_cst, *lt.second._p->first);
                        if(p_cst.is_indexed()){
                            res._indices = p_cst._indices;
                        }
                    }
                    else if(coef->is_number()) {
                        auto p_cst = *((constant<type>*)(coef.get()));
                        auto var_range = make_shared<pair<type,type>>(get_range(lt.second._p->second));
                        term_range = get_product_range(make_shared<pair<type,type>>(p_cst.eval(),p_cst.eval()),var_range);
                        res.insert(lt.second._sign, p_cst, *lt.second._p->first);
                    }
                    if (lt.second._coef->_is_vector || lt.second._coef->is_matrix()) {
                        res._is_vector = true;
                    }
                    if(lt.second._sign){
                        res._range = get_plus_range(res._range, term_range);
                    }
                    else {
                        res._range = get_minus_range(res._range, term_range);
                    }
                }
            }
            int expo = 0;
            for (auto &lt: *_pterms) {
                auto coef = lt.second._coef->copy();
                if (coef->_is_transposed) {
                    coef->transpose();
                    coef->_is_vector = false;
                }
                bool has_v = false;
                shared_ptr<list<pair<shared_ptr<param_>, int>>> newl = make_shared<list<pair<shared_ptr<param_>,int>>>();
                auto it = lt.second._l->begin();
                while(it!=lt.second._l->end()){
                    auto vv = it->first;
                    if (vv->get_name(false,false) == name) {
                        has_v = true;
                        expo = it->second;
                        if(expo!=1){
                            newl->push_back(make_pair<>(vv, expo-1));
                            res.set_max_dim(*vv);
                        }
                    }
                    else {
                        newl->push_back(*it);
                    }
                    it++;
                }
                if(has_v){//TODO fix range propagation for square terms
                    res.set_max_dim(v);
                    auto pterm_range = get_range(newl);
                    pterm_range.first *= expo;
                    pterm_range.second *= expo;
                    if (coef->is_function()) {
                        auto f_cst = *static_pointer_cast<func<type>>(lt.second._coef);
                        if(newl->size()==1 && newl->front().second==2){
                            res.insert(lt.second._sign, expo*f_cst, *newl->front().first,*newl->front().first);
                        }
                        else if(newl->size()==2 && newl->front().second==1 && newl->back().second==1){
                            res.insert(lt.second._sign, expo*f_cst, *newl->front().first,*newl->back().first);
                        }
                        else{
                            res.insert(lt.second._sign, expo*f_cst, *newl);
                        }
                        pterm_range = *get_product_range(make_shared<pair<type,type>>(pterm_range), f_cst._range);
                    }
                    else if(coef->is_param()) {
                        auto p_cst = *static_pointer_cast<param<type>>(lt.second._coef);
                        if(newl->size()==1 && newl->front().second==2){
                            res.insert(lt.second._sign, expo*p_cst, *newl->front().first,*newl->front().first);
                        }
                        else if(newl->size()==2 && newl->front().second==1 && newl->back().second==1){
                            res.insert(lt.second._sign, expo*p_cst, *newl->front().first,*newl->back().first);
                        }
                        else{
                            res.insert(lt.second._sign, expo*p_cst, *newl);
                        }
                        pterm_range = *get_product_range(make_shared<pair<type,type>>(pterm_range), p_cst._range);
                    }
                    else if(coef->is_number()) {
                        auto p_cst = *static_pointer_cast<constant<type>>(lt.second._coef);
                        if(newl->size()==1 && newl->front().second==2){
                            res.insert(lt.second._sign, expo*p_cst, *newl->front().first,*newl->front().first);
                        }
                        else if(newl->size()==2 && newl->front().second==1 && newl->back().second==1){
                            res.insert(lt.second._sign, expo*p_cst, *newl->front().first,*newl->back().first);
                        }
                        else{
                            res.insert(lt.second._sign, expo*p_cst, *newl);
                        }
                        pterm_range.first = gravity::min(pterm_range.first*p_cst.eval(), pterm_range.second*p_cst.eval());
                        pterm_range.second = gravity::max(pterm_range.first*p_cst.eval(), pterm_range.second*p_cst.eval());
                    }
                    if(lt.second._sign){
                        res._range = get_plus_range(res._range, make_shared<pair<type,type>>(pterm_range));
                    }
                    else {
                        res._range = get_minus_range(res._range, make_shared<pair<type,type>>(pterm_range));
                    }
                }
            }
            if (_expr) {
                res += get_derivative(_expr,v);
            }
            for (auto &lt: *_lterms) {
                if (lt.second._p->get_name(false,false) == name) {
                    auto coef = lt.second._coef->copy();
                    if (coef->_is_transposed) {
                        coef->transpose();
                    }
                    if (coef->is_function()) {
                        if(lt.second._sign){
                            res += *((func<type>*)(coef.get()));
                        }
                        else {
                            res -= *((func<type>*)(coef.get()));
                        }
                        
                    }
                    else if(coef->is_param()) {
                        auto p_cst = *((param<type>*)(coef.get()));
                        if(lt.second._sign){
                            res += p_cst;
                        }
                        else {
                            res -= p_cst;
                        }
                        if(p_cst.is_indexed()){
                            res._indices = p_cst._indices;
                        }
                    }
                    else if(coef->is_number()) {
                        if(lt.second._sign){
                            res += *((constant<type>*)(coef.get()));
                        }
                        else {
                            res -= *((constant<type>*)(coef.get()));
                        }
                    }
                    break;
                }
            }
            
            res.update_double_index();
            return res;
        }
        
        
        bool insert(bool sign, const constant_& coef, const param_& p){/**< Adds coef*p to the linear function. Returns true if added new term, false if only updated coef of p */
            shared_ptr<param_> p_new = p.pcopy();
            bool transpose = false;
            _evaluated=false;
            string pname;
            if(coef._is_transposed && !p_new->_is_vector){
                p_new->_is_vector=true;
                p_new->_name = "["+p_new->_name+"]";
            }
            if(coef.get_dim()>1 && p._is_transposed){// Situation where p^T*coef is transformed into coef^T*p
                if(coef._is_transposed){
                    throw invalid_argument("In  bool insert(bool sign, const constant_& coef, const param_& p), both coef and p are transposed.");
                }
                p_new->transpose();
                transpose = true;
            }
            pname = p_new->get_name(false,false);
            auto pair_it = _lterms->find(pname);
            if (pair_it != _lterms->end() && pair_it->second._p->get_type() != p.get_type()) {
                throw invalid_argument("param and var with same name: " + pname);
            }
            _evaluated = false;
            if (_ftype == const_ && p.is_var()) {
                _ftype = lin_;
            }
            
            if (pair_it == _lterms->end()) {
                auto c_new = coef.copy();
                if(transpose){
                    c_new->transpose();
                }
                if (c_new->is_function()) {
                    embed(*static_pointer_cast<func>(c_new));
                }
                if (p.is_var()) {
                    auto p_exist = get_var(pname);
                    if (!p_exist) {
                        add_var(p_new);
                    }
                    else {
                        incr_occ_var(pname);
                    }
                }
                else {
                    auto p_exist = get_var(pname);
                    if (!p_exist) {
                        add_param(p_new);
                    }
                    else {
                        incr_occ_param(pname);
                    }
                }
                _lterms->insert(make_pair<>(pname, lterm(sign, c_new, p_new)));
                return true;
            }
            else {
                if (pair_it->second._sign == sign) {
                    if (coef.is_function()) {
                        auto coef2 = *(func<type>*)(&coef);
                        pair_it->second._coef = add(pair_it->second._coef,coef2);
                    }
                    else if(coef.is_param()) {
                        auto coef2 = *(param<type>*)(&coef);
                        pair_it->second._coef = add(pair_it->second._coef,coef2);
                    }
                    else if(coef.is_number()) {
                        auto coef2 = *(constant<type>*)(&coef);
                        pair_it->second._coef = add(pair_it->second._coef,coef2);
                    }
                }
                else{
                    if (coef.is_function()) {
                        auto coef2 = *(func<type>*)(&coef);
                        pair_it->second._coef = subtract(pair_it->second._coef,coef2);
                    }
                    else if(coef.is_param()) {
                        auto coef2 = *(param<type>*)(&coef);
                        pair_it->second._coef = subtract(pair_it->second._coef,coef2);
                    }
                    else if(coef.is_number()) {
                        auto coef2 = *(constant<type>*)(&coef);
                        pair_it->second._coef = subtract(pair_it->second._coef,coef2);
                    }
                }
                if (pair_it->second._coef->is_zero()) {
                    if (p.is_var()) {
                        decr_occ_var(pname);
                    }
                    else{
                        decr_occ_param(pname);
                    }
                    _lterms->erase(pair_it);
                    if(is_constant()){
                        _ftype = const_;
                        _val->resize(1);
                    }
                }
                return false;
            }
        };
        
        bool is_rotated_soc(){
            if (_qterms->empty() || !_pterms->empty() || _expr) {
                return false;
            }
            unsigned nb_bilinear = 0, nb_quad = 0;
            Sign bilinear_sign = unknown_, quadratic_sign = unknown_, var1_sign = unknown_, var2_sign = unknown_;
            for (auto &qt_pair: *_qterms) {
                if (qt_pair.second._p->first!=qt_pair.second._p->second) {
                    bilinear_sign = qt_pair.second.get_all_sign();
                    var1_sign = qt_pair.second._p->first->get_all_sign();
                    var2_sign = qt_pair.second._p->second->get_all_sign();
                    if (bilinear_sign==unknown_ || var1_sign==neg_ || var2_sign==neg_) {
                        return false;
                    }
                    nb_bilinear++;
                    if (nb_bilinear > 1) {
                        return false;
                    }
                }
                else{
                    nb_quad++;
                    auto sign = qt_pair.second.get_all_sign();
                    if (quadratic_sign!=unknown_ && quadratic_sign!=sign) {
                        return false;
                    }
                    if (quadratic_sign!=unknown_ && quadratic_sign==bilinear_sign) {
                        return false;
                    }
                    else {
                        quadratic_sign = sign;
                    }
                }
            }
            if(nb_quad==0){
                return false;
            }
            if (bilinear_sign==pos_) {
                _all_convexity = concave_;
                return true;
            }
            else if(bilinear_sign==neg_) {
                _all_convexity = convex_;
                return true;
            }
            return false;
        };
        
        bool insert(const constant_& coef, const param_& p){
            return insert(true, coef, p);
        }
        
        bool insert(const param_& p){
            return insert(true, unit<type>(), p);
        }
        
        bool insert(const lterm& term){return insert(term._sign, *term._coef, *term._p);};
        
        bool insert(const param_& p1, const param_& p2, bool coef_p1_tr=false){
            return insert(true, unit<type>(), p1, p2, coef_p1_tr);
        };
        
        bool insert(const param_& p2, int exp){
            return insert(true, unit<type>(), p2, exp);
        };
        
        /**
         Reverse the sign of all terms in the function, also reverses convexity.
         */
        void reverse_sign(){
            _cst->reverse_sign();
            for (auto &pair: *_lterms) {
                pair.second.reverse_sign();
            }
            for (auto &pair: *_qterms) {
                pair.second.reverse_sign();
            }
            for (auto &pair: *_pterms) {
                pair.second.reverse_sign();
            }
            if(_expr){
                _expr->reverse_sign();
            }
            if(_evaluated){
                for (auto i = 0; i<_val->size(); i++) {
                    _val->at(i) = -1.*eval(i);
                }
            }
            reverse_convexity();
            reverse_all_sign();
            reverse_range();
        }
        
        void reverse_range(){
            auto temp = _range->first;
            _range->first = -1.*_range->second;
            _range->second = -1.*temp;
        }
        
        
        template<class T=type, typename enable_if<is_same<T, Cpx>::value>::type* = nullptr> void update_all_sign(){
            if (_range->first == Cpx(0,0) && _range->second == Cpx(0,0)) {
                _all_sign = zero_;
            }
            else if ((_range->second.real() < 0 && _range->second.imag() < 0)) {
                _all_sign = neg_;
            }
            else if ((_range->second.real() > 0 && _range->second.imag() > 0)) {
                _all_sign = pos_;
            }
            else if (_range->second == Cpx(0,0)) {
                _all_sign = non_pos_;
            }
            else if (_range->first == Cpx(0,0)) {
                _all_sign = non_neg_;
            }
            else {
                _all_sign = unknown_;
            }
        }
        
        template<typename T=type,
        typename std::enable_if<is_arithmetic<T>::value>::type* = nullptr> void update_all_sign() {
            if (_range->first == 0 && _range->second == 0) {
                _all_sign = zero_;
            }
            else if (_range->second < 0  && _range->first < 0) {
                _all_sign = neg_;
            }
            else if (_range->first > 0 && _range->second > 0) {
                _all_sign = pos_;
            }
            else if (_range->second == 0   && _range->first < 0) {
                _all_sign = non_pos_;
            }
            else if (_range->first == 0  && _range->second > 0) {
                _all_sign = non_neg_;
            }
            else {
                _all_sign = unknown_;
            }
        }
        
        void reverse_all_sign(){
            if(_all_sign==neg_){
                _all_sign=pos_;
            }
            else if(_all_sign==pos_){
                _all_sign=neg_;
            }
            else if(_all_sign==non_pos_){
                _all_sign=non_neg_;
            }
            else if(_all_sign==non_neg_){
                _all_sign=non_pos_;
            }
        }
        //        shared_ptr<constant_> multiply(const func_& f){
        //
        //            //        switch (f.get_return_type()) {
        //            //            case binary_: {
        //            //                auto newf = static_cast<func<bool>>(f);
        //            //                newf *= *this;
        //            //                return make_shared<func<bool>>(vv*c2);
        //            //                break;
        //            //            }
        //            //            case short_: {
        //            //                auto vv = static_pointer_cast<var<short>>(c1);
        //            //                return make_shared<func<short>>(vv*c2);
        //            //                break;
        //            //            }
        //            //            case integer_: {
        //            //                auto vv = static_pointer_cast<var<int>>(c1);
        //            //                return make_shared<func<int>>(vv*c2);
        //            //                break;
        //            //            }
        //            //            case float_: {
        //            //                auto vv = static_pointer_cast<var<float>>(c1);
        //            //                return make_shared<func<float>>(vv*c2);
        //            //                break;
        //            //            }
        //            //            case double_: {
        //            //                auto vv = static_pointer_cast<var<double>>(c1);
        //            //                return make_shared<func<double>>(vv*c2);
        //            //                break;
        //            //            }
        //            //            case long_: {
        //            //                auto vv = static_pointer_cast<var<long double>>(c1);
        //            //                return make_shared<func<long double>>(vv*c2);
        //            //                break;
        //            //            }
        //            //            case complex_: {
        //            //                auto vv = static_pointer_cast<var<Cpx>>(c1);
        //            //                return make_shared<func<Cpx>>(vv*c2);
        //            //                break;
        //            //            }
        //            //            default:
        //            //                break;
        //        }
        
        //        Sign get_all_sign() const{ /**< If all instances of the current parameter/variable have the same sign, it returns it, otherwise, it returns unknown. **/
        //            return get_all_sign();
        //        };
        //        Sign get_sign(size_t idx = 0) const{ /**< returns the sign of one instance of the current parameter/variable. **/
        //            return get_sign(idx);
        //        }
        
        //        template<typename... Args>
        //        func in(const indices& vec1, Args&&... args) {
        //            func<type> res(*this);
        ////            res.operator=(in(vec1, forward<Args>(args)...));
        //            return res;
        //        }
        
        template<typename... Args>
        void index_in(const indices& ids1, Args&&... args) {
            _indices = make_shared<indices>(ids1, args...);
        }
        
        void print() {
            allocate_mem();
            print(10);
        }
        
        void print(size_t index, int prec = 10) {
            cout << to_str(index,prec);
        }
        
        void print(size_t i, size_t j, int prec = 10) {
            cout << to_str(i,j,prec);
        }
        
        string to_str() {
            string str;
            for (auto &pair:*_pterms) {
                str += pair.second.to_str();
            }
            for (auto &pair:*_qterms) {
                str += pair.second.to_str();
            }
            for (auto &pair:*_lterms) {
                str += pair.second.to_str();
            }
            if(!_cst->is_zero()){
                auto cst_str = _cst->to_str();
                if (cst_str.front()=='-'){
                    str += " - " + cst_str.substr(1);
                }
                else {
                    str += " + " + cst_str;
                }
            }
            if (_expr) {
                str += " + ";
                str += _expr->to_str();
            }
            if (str.size() > 2 && str.at(1)=='+') {
                str = str.substr(3);
            }
            if (_is_vector) {
                str = "[" + str +"]";
            }
            if (_is_transposed) {
                str += "\u1D40";
            }
            if(str.size()==0){
                str = "0";
            }
            _to_str = str;
            return str;
        }
        
        void update_str(){
            _to_str = to_str();
        }
        
        void update_double_index() {
            if(_expr){
                _expr->update_double_index();
            }
            for (auto &v_p:*_vars) {
                auto v = v_p.second.first;
                if(v->is_indexed() && !v->_is_transposed){
                    _indices = v->_indices;
                    _dim[0] = v->_dim[0];
                }
                if(v->is_matrix_indexed()){
                    _indices = v->_indices;
                    return;
                }
            }
            for (auto &p_p:*_params) {
                auto p = p_p.second.first;
                if(p->is_indexed()&& !p->_is_transposed){
                    _indices = p->_indices;
                    _dim[0] = p->_dim[0];
                }
                if(p->is_matrix_indexed()){
                    _indices = p->_indices;
                    return;
                }
            }
        }
        
        bool is_matrix_indexed() const{
            return (_indices && _indices->_ids && _indices->_ids->size()>1);
        }
        
        string to_str(size_t index, int prec) {
            if (is_constant() && !this->is_matrix_indexed()) {
                return to_string_with_precision(eval(index),prec);
            }
            string str;
            for (auto &pair:*_pterms) {
                str += pair.second.to_str(index, prec);
            }
            for (auto &pair:*_qterms) {
                str += pair.second.to_str(index, prec);
            }
            for (auto &pair:*_lterms) {
                str += pair.second.to_str(index, prec);
            }
            if (_expr) {
                str += _expr->to_str(index, prec);
            }
            if(!_cst->is_zero()){
                auto val = _cst->to_str(index, prec);
                if (val.front()=='-') {
                    str += " - " + val.substr(1);
                }
                else {
                    str += " + ";
                    str += val;
                }
            }
            if (str.size() > 2 && str.at(1)=='+') {
                str = str.substr(3);
            }
            //            if (_is_vector && (is_number() || _vars->size()>1 || _params->size()>1)) {
            //                str = "[" + str +"]";
            //            }
            //            if (_is_transposed && (is_number() || _vars->size()>1 || _params->size()>1)) {
            //                str += "\u1D40";
            //            }
            return str;
        }
        
        size_t get_max_cell_size(){
            auto max_size = 0;
            for (size_t i = 0; i<_dim[0]; i++) {
                for (size_t j = 0; j<_dim[1]; j++) {
                    eval(i,j);
                    auto cell = to_str(i,j,5);
                    if(max_size < cell.size()){
                        max_size = cell.size();
                    }
                }
            }
            return max_size;
        }
        
        size_t get_dim() const{
            return constant_::get_dim();
        }
        
        size_t get_dim(size_t i) const{
            if(is_matrix_indexed())
                return _indices->_ids->at(i).size();
            if (is_indexed()) {
                return _indices->_ids->at(0).size();
            }
            return constant_::get_dim(i);
        }
        
        size_t get_nb_inst() const{
            if(is_matrix_indexed())
                return _indices->_ids->size();
            if(is_indexed() && !_is_transposed){
                return _indices->_ids->at(0).size();
            }
            return this->_dim[0];
        }
        
        void print(int prec){
            string str;
            if (is_constant()) {
                str += " (Constant";
            }
            else if (is_linear()) {
                str += " (Linear";
            }
            else if (is_convex()) {
                str += " (Convex";
            }
            else if (is_concave()){
                str += " (Concave";
            }
            else {
                str += " (Unknown";
            }
            if (is_complex()) {
                str += " Complex) : ";
            }
            else {
                str += ") : ";
            }
            if (!_embedded && !is_constant()) {
                str += "f(";
                for (auto pair_it = _vars->begin(); pair_it != _vars->end();) {
                    str += pair_it->second.first->get_name(false,true);
                    if (next(pair_it) != _vars->end()) {
                        str += ",";
                    }
                    pair_it++;
                }
                str += ") = ";
            }
            auto space_size = str.size();
            auto nb_inst = this->get_nb_inst();
            allocate_mem();
            if (is_matrix()) {
                auto max_cell_size = get_max_cell_size();
                for (size_t i = 0; i<_dim[0]; i++) {
                    if (i>0) {
                        str.insert(str.end(), space_size, ' ');
                    }
                    str += "|";
                    for (size_t j = 0; j<_dim[1]; j++) {
                        auto cell = to_str(i,j,prec);
                        auto cell_size = cell.size();
                        cell.insert(0, floor((max_cell_size - cell_size)/2.), ' ');
                        cell.append(ceil((max_cell_size - cell_size)/2.), ' ');
                        str += cell;
                        if(j!=_dim[1]-1){
                            str += " ";
                        }
                    }
                    str += "|\n";
                }
            }
            else {
                for (size_t inst = 0; inst<nb_inst; inst++) {
                    eval(inst);
                    if (inst>0) {
                        str.insert(str.end(), space_size, ' ');
                    }
                    str += to_str(inst,prec);
                    str += "\n";
                }
            }
            cout << str;
        }
        
        
        string to_str(size_t index1, size_t index2, int prec) {
            if (is_constant()) {
                return to_string_with_precision(eval(index1,index2),prec);
            }
            string str;
            for (auto &pair:*_pterms) {
                str += pair.second.to_str(index1,index2, prec);
            }
            for (auto &pair:*_qterms) {
                str += pair.second.to_str(index1,index2, prec);
            }
            for (auto &pair:*_lterms) {
                str += pair.second.to_str(index1,index2, prec);
            }
            if(!_cst->is_zero()){
                auto val = _cst->to_str(index1,index2, prec);
                if (val.front()=='-') {
                    str += " - " + val.substr(1);
                }
                else {
                    str += " + ";
                    str += val;
                }
            }
            if (_expr) {
                str += " + ";
                str += _expr->to_str(index1,index2, prec);
            }
            if (str.size() > 2 && str.at(1)=='+') {
                str = str.substr(3);
            }
            //            if (_is_vector && (func_is_number() || _vars->size()>1 || _params->size()>1)) {
            //                str = "[" + str +"]";
            //            }
            //            if (_is_transposed && (func_is_number() || _vars->size()>1 || _params->size()>1)) {
            //                str += "\u1D40";
            //            }
            return str;
        }
        
        
        void propagate_dim(size_t d){
            if (is_matrix()) {
                return;
            }
            if(_is_transposed){
                _dim[1] = d;
            }
            else {
                _dim[0] = d;
            }
            for (auto &pair:*_lterms) {
                auto coef = pair.second._coef;
                coef->propagate_dim(d);
            }
            for (auto &pair:*_qterms) {
                auto coef = pair.second._coef;
                coef->propagate_dim(d);
            }
            for (auto &pair:*_pterms) {
                auto coef = pair.second._coef;
                coef->propagate_dim(d);
            }
            _cst->propagate_dim(d);
            if (_expr) {
                _expr->propagate_dim(d);
            }
        }
        
        void allocate_mem(){
            //            _evaluated = false;
            if(is_matrix_indexed()){
                for(auto i = 0; i<_indices->_ids->size();i++){
                    for(auto j = 0; j<_indices->_ids->at(i).size();j++){
                        _dim[0] = max(_dim[0],_indices->_ids->at(i).at(j)+1);
                    }
                }
            }
            _val->resize(get_dim());
            for (auto &pair:*_lterms) {
                auto coef = pair.second._coef;
                coef->allocate_mem();
            }
            for (auto &pair:*_qterms) {
                auto coef = pair.second._coef;
                coef->allocate_mem();
            }
            for (auto &pair:*_pterms) {
                auto coef = pair.second._coef;
                coef->allocate_mem();
            }
            _cst->allocate_mem();
            if (_expr) {
                _expr->allocate_mem();
            }
        }
        
        /** Returns true if the current function is just a wrapper to a parameter */
        bool func_is_param() const{
            if(_vars->size()==1 && _params->size()==0 && _vars->begin()->second.first->get_intype()==_return_type)
                return true;
            if(_vars->size()==0 && _params->size()==1 && _params->begin()->second.first->get_intype()==_return_type)
                return true;
            return false;
        }
        
        param<type> func_get_param() const{
            auto p = *(param<type>*)(_lterms->begin()->second._p.get());
            if(_is_transposed && ! p._is_transposed){
                p.transpose();
            }
            return p;
        }
        
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        shared_ptr<constant_> multiply(shared_ptr<constant_> coef, const constant<T2>& c){
            if (coef->is_function()) {
                auto f_cst = *((func<type>*)(coef.get()));
                //                if(c.is_unit() && f_cst.func_is_param()){
                //                    return f_cst.func_get_param().copy();
                //                }
                f_cst *= func<type>(c);
                return f_cst.copy();
            }
            else if(coef->is_param()) {
                auto p_cst = *((param<type>*)(coef.get()));
                if(c.is_unit()){
                    return p_cst.pcopy();
                }
                auto new_cst = p_cst * c;
                return new_cst.copy();
            }
            else if(coef->is_number()) {
                auto p_cst = *((constant<type>*)(coef.get()));
                auto new_cst = p_cst * c;
                return new_cst.copy();
            }
            return nullptr;
        }
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        shared_ptr<constant_> multiply(shared_ptr<constant_> coef, const param<T2>& p){
            if (coef->is_function()) {
                auto f_cst = *((func<type>*)(coef.get()));
                f_cst *= func<type>(p);
                return f_cst.copy();
            }
            else if(coef->is_param()) {
                auto p_cst = *((param<type>*)(coef.get()));
                auto new_cst = p_cst * p;
                return new_cst.copy();
            }
            else if(coef->is_number()) {
                auto p_cst = *((constant<type>*)(coef.get()));
                auto new_cst = p_cst * p;
                return new_cst.copy();
            }
            return nullptr;
        }
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        shared_ptr<constant_> multiply(shared_ptr<constant_> coef, const func<T2>& f){
            if (coef->is_function()) {
                auto f_cst = *((func<type>*)(coef.get()));
                f_cst *= func<type>(f);
                embed(f_cst);
                return f_cst.copy();
            }
            else if(coef->is_param()) {
                auto p_cst = *((param<type>*)(coef.get()));
                auto new_cst = p_cst * f;
                return new_cst.copy();
            }
            else if(coef->is_number()) {
                auto p_cst = *((constant<type>*)(coef.get()));
                if(f.func_is_number()){
                    return make_shared<constant<type>>(p_cst *= eval(f.copy(),0));
                }
                return make_shared<func<type>>(func(p_cst) *= f);
            }
            return nullptr;
        }
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        shared_ptr<constant_> add(shared_ptr<constant_> coef, const constant<T2>& c){
            if (coef->is_function()) {
                auto f_cst = *((func<type>*)(coef.get()));
                f_cst += func<type>(c);
                return f_cst.copy();
            }
            else if(coef->is_param()) {
                auto p_cst = *((param<type>*)(coef.get()));
                auto new_cst = p_cst + c;
                return new_cst.copy();
            }
            else if(coef->is_number()) {
                auto p_cst = *((constant<type>*)(coef.get()));
                auto new_cst = p_cst + c;
                return new_cst.copy();
            }
            return nullptr;
        }
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        shared_ptr<constant_> add(shared_ptr<constant_> coef, const param<T2>& p){
            if (coef->is_function()) {
                auto f_cst = *((func<type>*)(coef.get()));
                f_cst += func<type>(p);
                return f_cst.copy();
            }
            else if(coef->is_param()) {
                auto p_cst = *((param<type>*)(coef.get()));
                auto new_cst = p_cst + p;
                return new_cst.copy();
            }
            else if(coef->is_number()) {
                auto p_cst = *((constant<type>*)(coef.get()));
                auto new_cst = p_cst + p;
                return new_cst.copy();
            }
            return nullptr;
        }
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        shared_ptr<constant_> add(shared_ptr<constant_> coef, const func<T2>& f){
            if (coef->is_function()) {
                auto f_cst = *((func<type>*)(coef.get()));
                f_cst += func<type>(f);
                embed(f_cst);
                return f_cst.copy();
            }
            else if(coef->is_param()) {
                auto p_cst = *((param<type>*)(coef.get()));
                auto new_cst = p_cst + f;
                return new_cst.copy();
            }
            else if(coef->is_number()) {
                auto p_cst = (*((constant<type>*)(coef.get())));
                if(f.func_is_number()){
                    return make_shared<constant<type>>(p_cst += eval(f.copy(),0));
                }
                return make_shared<func<type>>(func(p_cst) += f);
            }
            return nullptr;
        }
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        shared_ptr<constant_> subtract(shared_ptr<constant_> coef, const constant<T2>& c){
            if (coef->is_function()) {
                auto f_cst = *((func<type>*)(coef.get()));
                f_cst -= func<type>(c);
                return f_cst.copy();
            }
            else if(coef->is_param()) {
                auto p_cst = *((param<type>*)(coef.get()));
                if(c.is_unit()){
                    return p_cst.pcopy();
                }
                auto new_cst = p_cst - c;
                return new_cst.copy();
            }
            else if(coef->is_number()) {
                auto p_cst = *((constant<type>*)(coef.get()));
                auto new_cst = p_cst - c;
                return new_cst.copy();
            }
            return nullptr;
        }
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        shared_ptr<constant_> subtract(shared_ptr<constant_> coef, const param<T2>& p){
            if (coef->is_function()) {
                auto f_cst = *((func<type>*)(coef.get()));
                f_cst -= func<type>(p);
                return f_cst.copy();
            }
            else if(coef->is_param()) {
                auto p_cst = *((param<type>*)(coef.get()));
                auto new_cst = p_cst - p;
                return new_cst.copy();
            }
            else if(coef->is_number()) {
                auto p_cst = *((constant<type>*)(coef.get()));
                auto new_cst = p_cst - p;
                return new_cst.copy();
            }
            return nullptr;
        }
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        shared_ptr<constant_> subtract(shared_ptr<constant_> coef, const func<T2>& f){
            if (coef->is_function()) {
                auto f_cst = *((func<type>*)(coef.get()));
                f_cst -= func<type>(f);
                embed(f_cst);
                return f_cst.copy();
            }
            else if(coef->is_param()) {
                auto p_cst = *((param<type>*)(coef.get()));
                auto new_cst = p_cst - f;
                return new_cst.copy();
            }
            else if(coef->is_number()) {
                auto p_cst = (*((constant<type>*)(coef.get())));
                if(f.func_is_number()){
                    return make_shared<constant<type>>(p_cst -= eval(f.copy(),0));
                }
                return make_shared<func<type>>(func(p_cst) -= f);
            }
            return nullptr;
        }
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        void add_cst(const constant<T2>& f){
            if (_cst->is_function()) {
                auto f_cst = *static_pointer_cast<func<type>>(_cst);
                if(f_cst.func_is_number()){
                    _cst = make_shared<constant<type>>(eval(f_cst.copy(),0) + eval(f.copy(),0));
                }
                else {
                    f_cst += func<type>(f);
                    _cst = make_shared<func<type>>(move(f_cst));
                }
            }
            else if(_cst->is_param()) {
                auto p_cst = *static_pointer_cast<param<type>>(_cst);
                auto new_cst = f + p_cst;
                _cst = make_shared<func<type>>(move(new_cst));
            }
            else if(_cst->is_number()) {
                auto p_cst = *static_pointer_cast<constant<type>>(_cst);
                auto new_cst = f + p_cst;
                _cst = make_shared<constant<type>>(new_cst);
            }
        }
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        void add_cst(const param<T2>& f){
            if (_cst->is_function()) {
                auto f_cst = *static_pointer_cast<func<type>>(_cst);
                f_cst += func<type>(f);
                embed(f_cst);
                _cst = make_shared<func<type>>(move(f_cst));
            }
            else if(_cst->is_param()) {
                auto p_cst = *static_pointer_cast<param<type>>(_cst);
                auto f_cst = f + p_cst;
                embed(f_cst);
                _cst = make_shared<func<type>>(move(f_cst));
            }
            else if(_cst->is_number()) {
                auto p_cst = *static_pointer_cast<constant<type>>(_cst);
                auto f_cst = f + p_cst;
                embed(f_cst);
                _cst = make_shared<func<type>>(move(f_cst));
            }
        }
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        void add_cst(const func<T2>& f){
            if (_cst->is_function()) {
                auto f_cst = *static_pointer_cast<func<type>>(_cst);
                if(f_cst.func_is_number() && f.func_is_number()){
                    _cst = make_shared<constant<type>>(eval(f_cst.copy(),0) + eval(f.copy(),0));
                }
                else {
                    f_cst += f;
                    embed(f_cst);
                    _cst = make_shared<func<type>>(move(f_cst));
                }
            }
            else if(_cst->is_param()) {
                auto p_cst = *static_pointer_cast<param<type>>(_cst);
                auto f_cst = f + func<type>(p_cst);
                embed(f_cst);
                _cst = make_shared<func<type>>(move(f_cst));
            }
            else if(_cst->is_number()) {
                auto p_cst = *static_pointer_cast<constant<type>>(_cst);
                if(f.func_is_number()){
                    _cst = make_shared<constant<type>>(p_cst += eval(f.copy(),0));
                }
                else {
                    auto f_cst = f + func<type>(p_cst);
                    embed(f_cst);
                    _cst = make_shared<func<type>>(move(f_cst));
                }
            }
        }
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        func(const uexpr<T2>& ue):func(){
            _expr = make_shared<uexpr<type>>(ue);
            embed(_expr);
            if (!is_constant()) {
                _ftype = nlin_;
            }
            _dim[0] = ue._dim[0];
            _dim[1] = ue._dim[1];
            _evaluated = false;
            _range->first = ue._range->first;
            _range->second = ue._range->second;
            _all_convexity = ue._all_convexity;
            _all_sign = ue._all_sign;
            
        };
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        func(const bexpr<T2>& be):func(){
            _expr = make_shared<bexpr<type>>(be);
            embed(_expr);
            if (!is_constant()) {
                _ftype = nlin_;
            }
            _dim[0] = be._dim[0];
            _dim[1] = be._dim[1];
            _evaluated = false;
            _range->first = be._range->first;
            _range->second = be._range->second;
            _all_convexity = be._all_convexity;
            _all_sign = be._all_sign;
        };
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        func(T2 c):func(){
            *this = constant<T2>(c);
        }
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        func(const constant<T2>& c):func(){
            *this = c;
        }
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        func(const param<T2>& c):func(){
            *this = c;
        }
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        func(const var<T2>& c):func(){
            *this = c;
        }
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) < sizeof(type)>::type* = nullptr>
        func(const func<T2>& f): func(){
            *this = f;
        }
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        func& operator=(const constant<T2>& c){
            reset();
            static_pointer_cast<constant<type>>(_cst)->set_val(c.eval());
            _all_sign = _cst->get_sign();
            _val->resize(1);
            _val->at(0) = c.eval();
            set_range(_val->at(0));
            _all_sign = c.get_all_sign();
            _is_vector = c._is_vector;
            _is_transposed = c._is_transposed;
            _dim[0] = c._dim[0];
            _dim[1] = c._dim[1];
            _evaluated = true;
            return *this;
        }
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        func& operator=(const param<T2>& c){
            reset();
            insert(true,unit<type>(),c);
            _dim[0] = c.get_nb_inst();
            _dim[1] = c._dim[1];
            _is_transposed = c._is_transposed;
            _is_vector = c._is_vector;
            _val->clear();
            _range->first = c._range->first;
            _range->second = c._range->second;
            _all_sign = c.get_all_sign();
            _evaluated = false;
            if(c._indices){
                _indices = make_shared<indices>(*c._indices);
            }
            //            if(c.is_matrix_indexed()){
            //                _indices = c._indices;
            //            }
            return *this;
        }
        
        
        
        func(func&& f){
            *this = move(f);
        }
        
        func(const func& f){
            *this = f;
        }
        
        
        
        shared_ptr<func_> fcopy() const{return make_shared<func>(*this);};
        
        shared_ptr<constant_> copy()const{return make_shared<func>(*this);};
        
        void deep_copy(const func& f){
            constant_::_type = f._type;
            _ftype = f._ftype;
            _return_type = f._return_type;
            _to_str = f._to_str;
            _all_convexity = f._all_convexity;
            _all_sign = f._all_sign;
            _cst = f._cst->copy();
            _val = make_shared<vector<type>>();
            _range = make_shared<pair<type,type>>();
            _vars = make_shared<map<string, pair<shared_ptr<param_>, unsigned>>>();
            _params = make_shared<map<string, pair<shared_ptr<param_>, unsigned>>>();
            _lterms = make_shared<map<string, lterm>>();
            _qterms = make_shared<map<string, qterm>>();
            _pterms = make_shared<map<string, pterm>>();
            for (auto &pair:*f._lterms) {
                this->insert(pair.second);
            }
            for (auto &pair:*f._qterms) {
                this->insert(pair.second);
            }
            for (auto &pair:*f._pterms) {
                this->insert(pair.second);
            }
            if(f._expr){
                if (f._expr->is_uexpr()) {
                    _expr = make_shared<uexpr<type>>(*static_pointer_cast<uexpr<type>>(f._expr));
                }
                else {
                    _expr = make_shared<bexpr<type>>(*static_pointer_cast<bexpr<type>>(f._expr));
                }
                embed(_expr);
            }
            else {
                _expr = nullptr;
            }
            if(f._indices){
                _indices = make_shared<indices>(*f._indices);
            }
            else {
                _indices = nullptr;
            }
            _range->first = f._range->first;
            _range->second = f._range->second;
            //            _val->clear();
            _val->resize(f._val->size());
            for(auto i = 0; i< f._val->size(); i++){
                _val->at(i) = f._val->at(i);
            }
            *_convexity = *f._convexity;
            _sign = f._sign;
            constant_::_is_transposed = f._is_transposed;
            constant_::_is_vector = f._is_vector;
            if(f._is_constraint)
                _is_constraint = true;
            _is_hessian = f._is_hessian;
            constant_::_dim[0] = f._dim[0];
            constant_::_dim[1] = f._dim[1];
            _embedded = f._embedded;
            _dfdx = make_shared<map<string,shared_ptr<func<type>>>>();
            //            copy_derivatives(f);
            if(f._hess_link){
                _hess_link = make_shared<map<size_t, set<size_t>>>(*f._hess_link);
            }
            else {
                _hess_link = nullptr;
            }
            _nnz_j = f._nnz_j;
            _nnz_h = f._nnz_h;
            _hess_link = f._hess_link;
            _nb_vars = f._nb_vars;
            _evaluated = f._evaluated;
        }
        
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) < sizeof(type)>::type* = nullptr>
        void deep_copy(const func<T2>& f){
            _to_str = f._to_str;
            _all_convexity = f._all_convexity;
            _all_sign = f._all_sign;
            if (f._cst->is_function()) {
                auto coef = *static_pointer_cast<func<T2>>(f._cst);
                _cst = func(coef).copy();
            }
            else if(f._cst->is_param()) {
                auto coef = *static_pointer_cast<param<T2>>(f._cst);
                _cst = param<type>(coef).copy();
            }
            else if(f._cst->is_number()) {
                auto coef = *static_pointer_cast<constant<T2>>(f._cst);
                _cst = constant<type>(coef.eval()).copy();
            }
            _val = make_shared<vector<type>>();
            _range = make_shared<pair<type,type>>();
            _vars = make_shared<map<string, pair<shared_ptr<param_>, unsigned>>>();
            _params = make_shared<map<string, pair<shared_ptr<param_>, unsigned>>>();
            _lterms = make_shared<map<string, lterm>>();
            _qterms = make_shared<map<string, qterm>>();
            _pterms = make_shared<map<string, pterm>>();
            for (auto &pair:*f._lterms) {
                auto term = pair.second;
                if (term._coef->is_function()) {
                    auto coef = *static_pointer_cast<func<T2>>(term._coef);
                    term._coef = func(coef).copy();
                }
                else if(term._coef->is_param()) {
                    auto coef = *static_pointer_cast<param<T2>>(term._coef);
                    term._coef = param<type>(coef).copy();
                }
                else if(term._coef->is_number()) {
                    auto coef = *static_pointer_cast<constant<T2>>(term._coef);
                    term._coef = constant<type>(coef).copy();
                }
                this->insert(term);
            }
            for (auto &pair:*f._qterms) {
                auto term = pair.second;
                if (term._coef->is_function()) {
                    auto coef = *static_pointer_cast<func<T2>>(term._coef);
                    term._coef = func(coef).copy();
                }
                else if(term._coef->is_param()) {
                    auto coef = *static_pointer_cast<param<T2>>(term._coef);
                    term._coef = param<type>(coef).copy();
                }
                else if(term._coef->is_number()) {
                    auto coef = *static_pointer_cast<constant<T2>>(term._coef);
                    term._coef = constant<type>(coef.eval()).copy();
                }
                this->insert(term);
            }
            for (auto &pair:*f._pterms) {
                auto term = pair.second;
                if (term._coef->is_function()) {
                    auto coef = *static_pointer_cast<func<T2>>(term._coef);
                    term._coef = func(coef).copy();
                }
                else if(term._coef->is_param()) {
                    auto coef = *static_pointer_cast<param<T2>>(term._coef);
                    term._coef = param<type>(coef).copy();
                }
                else if(term._coef->is_number()) {
                    auto coef = *static_pointer_cast<constant<T2>>(term._coef);
                    term._coef = constant<type>(coef.eval()).copy();
                }
                this->insert(term);
            }
            if(f._expr){
                if (f._expr->is_uexpr()) {
                    auto uexp = make_shared<uexpr<type>>(*static_pointer_cast<uexpr<T2>>(f._expr));
                    if (uexp->_son->is_function()) {
                        auto f = static_pointer_cast<func<T2>>(uexp->_son);
                        uexp->_son = make_shared<func>(*f);
                    }
                    _expr = uexp;
                }
                else {
                    auto bexp = make_shared<bexpr<type>>(*static_pointer_cast<bexpr<T2>>(f._expr));
                    if (bexp->_lson->is_function()) {
                        auto f = static_pointer_cast<func<T2>>(bexp->_lson);
                        bexp->_lson = make_shared<func>(*f);
                    }
                    if (bexp->_rson->is_function()) {
                        auto f = static_pointer_cast<func<T2>>(bexp->_rson);
                        bexp->_rson = make_shared<func>(*f);
                    }
                    _expr = bexp;
                }
            }
            else {
                _expr = nullptr;
            }
            if(f._indices){
                _indices = make_shared<indices>(*f._indices);
            }
            else {
                _indices = nullptr;
            }
            _range->first = f._range->first;
            _range->second = f._range->second;
            _val->clear();
            _val->resize(f._val->size());
            for(auto i = 0; i< f._val->size(); i++){
                _val->at(i) = f._val->at(i);
            }
            *_convexity = *f._convexity;
            _sign = f._sign;
            constant_::_is_transposed = f._is_transposed;
            constant_::_is_vector = f._is_vector;
            if(f._is_constraint)
                _is_constraint = true;
            _is_hessian = f._is_hessian;
            constant_::_dim[0] = f._dim[0];
            constant_::_dim[1] = f._dim[1];
            _embedded = f._embedded;
            _dfdx = make_shared<map<string,shared_ptr<func>>>();
            //            copy_derivatives(f);
            if(f._hess_link){
                _hess_link = make_shared<map<size_t, set<size_t>>>(*f._hess_link);
            }
            else {
                _hess_link = nullptr;
            }
            _nnz_j = f._nnz_j;
            _nnz_h = f._nnz_h;
            _hess_link = f._hess_link;
            _nb_vars = f._nb_vars;
            _evaluated = f._evaluated;
        }
        
        func& operator=(const func& f){
            deep_copy(f);
            return  *this;
        }
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        func& operator=(const func<T2>& f){
            deep_copy(f);
            return  *this;
        }
        
        func& operator=(func&& f){
            constant_::_type = f._type;
            _ftype = f._ftype;
            _to_str = f._to_str;
            _return_type = f._return_type;
            _all_convexity = f._all_convexity;
            _all_sign = f._all_sign;
            _lterms = move(f._lterms);
            _qterms = move(f._qterms);
            _pterms = move(f._pterms);
            _expr = move(f._expr);
            _vars = move(f._vars);
            _params = move(f._params);
            _cst = move(f._cst);
            _indices = move(f._indices);
            _range = move(f._range);
            _val = move(f._val);
            _convexity = move(f._convexity);
            _sign = f._sign;
            f._sign = nullptr;
            constant_::_is_transposed = f._is_transposed;
            constant_::_is_vector = f._is_vector;
            if(f._is_constraint)
                _is_constraint = true;
            _is_hessian = f._is_hessian;
            constant_::_dim[0] = f._dim[0];
            constant_::_dim[1] = f._dim[1];
            _embedded = f._embedded;
            _dfdx = move(f._dfdx);
            _nnz_j = f._nnz_j;
            _nnz_h = f._nnz_h;
            _hess_link = f._hess_link;
            _nb_vars = f._nb_vars;
            _evaluated = f._evaluated;
            return *this;
        }
        
        /* Modifiers */
        void    set_size(vector<size_t> dims){
            if (dims.size()==1) {
                set_size(dims[0]);
            }
            else if (dims.size()==2){
                set_size(dims[0],dims[1]);
            }
            else {
                throw invalid_argument("In Function set_size(vector<size_t> dims), dims.size() should be less or equal 2. \n");
            }
        }
        
        void   set_size(size_t s1, size_t s2) {
            _dim[0] = s1;
            _dim[1] = s2;
            auto dim = _dim[0]*_dim[1];
            _val->resize(dim);
            if (is_matrix()) {
                _is_vector = true;
            }
        };
        
        void   set_size(size_t s) {
            _val->resize(s);
            _dim[0] = s;
        };
        
        
        void add_val(type val) {
            if(is_matrix()){
                throw invalid_argument("Cannot call func::add_val(type val) on matrix");
            }
            _val->push_back(val);
            update_range(val);
            _dim[0] = _val->size();
        }
        
        void set_range(type val) {
            _range->first = val;
            _range->second = val;
        }
        
        void update_range(type val) {
            if (val <= _range->first) {
                _range->first = val;
            }
            if (val >= _range->second) {
                _range->second = val;
            }
        }
        
        void add_val(size_t i, type val) {
            if(is_matrix()){
                throw invalid_argument("Cannot call func::add_val(type val) on matrix");
            }
            _dim[0] = gravity::max(_dim[0],i+1);
            _val->resize(gravity::max(_val->size(),i+1));
            _val->at(i) = val;
            update_range(val);
        }
        
        void set_val(size_t i, size_t j, type val) {
            if(!is_matrix()){
                throw invalid_argument("Function set_val(size_t i, size_t j, type val) should be called on a matrix");
            }
            if(_dim[0] <= i || _dim[1] <= j){
                throw invalid_argument("In Function set_val(size_t i, size_t j, type val), i or j are out of bounds");
            }
            if (_is_transposed) {
                _val->at(_dim[0]*j+i) = val;
            }
            _val->at(_dim[1]*i+j) = val;
            update_range(val);
        }
        
        
        
        size_t set_val(const string& key, type val) {
            auto it = _indices->_keys_map->find(key);
            if (it == _indices->_keys_map->end()){
                throw invalid_argument("in Function size_t set_val(const string& key, type val), unknown key" + key);
            }
            _val->at(it->second) = val;
            update_range(val);
            return it->second;
        }
        
        size_t add_val(const string& key, type val) {
            if(!_indices){
                _indices = make_shared<indices>();
            }
            auto index = _indices->size();
            auto pp = _indices->_keys_map->insert(make_pair<>(key,index));
            if (pp.second) {//new index inserted
                _val->resize(gravity::max(_val->size(),index+1));
                _dim[0] = gravity::max(_dim[0],_val->size());
                _indices->_keys->resize(_val->size());
                _indices->_keys->at(index) = key;
                _val->at(index) = val;
                update_range(val);
                return index;
            }
            else {
                Warning("WARNING: calling add_val(const string& key, T val) with an existing key, overriding existing value" << endl);
                _val->at(pp.first->second) = val;
                update_range(val);
                return pp.first->second;
            }
        }
        
        void add_val(size_t i, size_t j, type val) {
            _is_vector = true;
            _dim[0] = gravity::max(_dim[0],i+1);
            _dim[1] = gravity::max(_dim[1],j+1);
            auto index = _dim[1]*i+j;
            _val->resize(gravity::max(_val->size(),index+1));
            _val->at(index) = val;
            update_range(val);
        }
        
        void set_val(size_t i, type val) {
            if(is_matrix()){
                throw invalid_argument("set_val(size_t i, type val) should be called with double index here\n");
            }
            if (is_indexed()) {
                if (_indices->_ids->size()>1) {
                    throw invalid_argument("set_val(size_t i, type val) should be called with double index here\n");
                }
                if (_val->size()<=_indices->_ids->at(0).at(i)){
                    throw invalid_argument("Param set_val(size_t i, type val) out of range");
                }
                _val->at(_indices->_ids->at(0).at(i)) = val;
            }
            if (_val->size()<=i){
                throw invalid_argument("Param set_val(size_t i, type val) out of range");
            }
            _val->at(i) = val;
            update_range(val);
        }
        
        void set_val(type val) {
            if(is_indexed()){
                for(auto &idx: _indices->_ids->at(0)){
                    _val->at(idx) = val;
                }
            }
            else {
                for (auto i = 0; i<_val->size() ;i++) {
                    _val->at(i) = val;
                }
            }
            update_range(val);
        }
        
        template<typename T=type,
        typename enable_if<is_arithmetic<T>::value>::type* = nullptr> Sign get_sign(size_t idx) const{
            if (_val->at(idx)==0) {
                return zero_;
            }
            if (_val->at(idx)< 0) {
                return neg_;
            }
            if (_val->at(idx)> 0) {
                return pos_;
            }
            return unknown_;
        }
        
        template<class T=type, typename enable_if<is_same<T, Cpx>::value>::type* = nullptr> Sign get_sign(size_t idx) const{
            if (_val->at(idx) == Cpx(0,0)) {
                return zero_;
            }
            if ((_val->at(idx).real() < 0 && _val->at(idx).imag() < 0)) {
                return neg_;
            }
            if ((_val->at(idx).real() > 0 && _val->at(idx).imag() > 0)) {
                return pos_;
            }
            if ((_val->at(idx).real() <= 0 && _val->at(idx).imag() <= 0)) {
                return non_pos_;
            }
            if ((_val->at(idx).real() >= 0 && _val->at(idx).imag() >= 0)) {
                return non_neg_;
            }
            return unknown_;
        }
        
        
        //        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        //        shared_ptr<constant_> add(shared_ptr<param<T2>> c1){
        //
        //        }
        
        
        
        Sign get_all_sign() const {
            return _all_sign;
        }
        
        
        
        type get_val() const {
            if (is_indexed()) {
                return _val->at(_indices->_ids->at(0).back());
            }
            return _val->back();
        }
        
        type get_val(size_t i) const{
            if(func_is_number()){
                return _val->at(0);
            }
            auto idx = get_id_inst(i);
            if (is_indexed()) {
                if (_indices->_ids->size()>1) {
                    throw invalid_argument("eval() should be called with double index here\n");
                }
                if (_val->size()<=idx){
                    throw invalid_argument("Param eval out of range");
                }
                return _val->at(idx);
            }
            if (_val->size()<=idx){
                throw invalid_argument("Param eval out of range");
            }
            return _val->at(idx);
        }
        
        type get_val(size_t i, size_t j) const{
            auto idx = get_id_inst(i,j);
            if (is_indexed()) {
                //                if (_indices->_ids->size()>1) {
                //                    throw invalid_argument("eval() should be called with double index here\n");
                //                }
                if (_val->size()<=idx){
                    throw invalid_argument("Param eval out of range");
                }
                return _val->at(idx);
            }
            if (_val->size()<=idx){
                throw invalid_argument("Param eval out of range");
            }
            return _val->at(idx);
        }
        
        void eval_matrix() {
            for (size_t i = 0; i < _dim[0]; i++) {
                for (size_t j = 0; j < _dim[1]; j++) {
                    set_val(i,j,eval(i,j));
                }
            }
        }
        
        
        
        void uneval() {
            _evaluated = false;
            _cst->uneval();
            for (auto &pair:*_lterms) {
                auto coef = pair.second._coef;
                if(coef->is_function()){
                    static_pointer_cast<func>(coef)->uneval();
                }
            }
            for (auto &pair:*_qterms) {
                auto coef = pair.second._coef;
                if (coef->is_function()){
                    static_pointer_cast<func>(coef)->uneval();
                }
            }
            for (auto &pair:*_pterms) {
                auto coef = pair.second._coef;
                if(coef->is_function()){
                    static_pointer_cast<func>(coef)->uneval();
                }
            }
            if(_expr){
                _expr->uneval();
            }
        }
        
        void eval_all(){
            //            if(_val->size()==0){
            allocate_mem();
            //            }
            auto nb_inst = get_nb_inst();
            for (size_t inst = 0; inst<nb_inst; inst++) {
                eval(inst);
            }
            _evaluated = true;
        }
        
        type eval(size_t i=0) {
            if(is_zero()){
                if (func_is_number()){
                    assert(_val->size()>0);
                    _val->at(0) = this->_range->first;
                    return _val->at(0);
                }
                assert(_val->size()>i);
                _val->at(i) = this->_range->first;
                return _range->first;
            }
            //            if (is_constant() && _evaluated) {
            if (_evaluated) {
                if (func_is_number()){
                    assert(_val->size()>0);
                    return _val->at(0);
                }
                assert(_val->size()>i);
                return _val->at(i);
            }
            type res = 0;
            if(!_cst->is_zero())
                res += eval_cst(i);
            if(!_lterms->empty()){
                for (auto &pair:*_lterms) {
                    if (pair.second._coef->_is_transposed || pair.second._coef->is_matrix() || pair.second._p->is_matrix_indexed()) {
                        auto dim = pair.second._p->get_dim(i);
                        if (pair.second._sign) {
                            for (size_t j = 0; j<dim; j++) {
                                res += eval_coef(pair.second._coef,i,j) * eval(pair.second._p,i,j);
                            }
                        }
                        else {
                            for (size_t j = 0; j<dim; j++) {
                                res -= eval_coef(pair.second._coef,i,j) * eval(pair.second._p,i,j);
                            }
                        }
                    }
                    else {
                        if (pair.second._sign) {
                            res += eval_coef(pair.second._coef,i) * eval(pair.second._p,i);
                        }
                        else {
                            res -= eval_coef(pair.second._coef,i) * eval(pair.second._p,i);
                        }
                    }
                }
            }
            //                res += eval_lterms(i);
            if(!_qterms->empty()){
                for (auto &pair:*_qterms) {
                    type qval = 0;
                    if(pair.second._p->second->is_matrix_indexed()){
                        auto dim = pair.second._p->first->get_dim(i);
                        if (pair.second._sign) {
                            for (size_t j = 0; j<dim; j++) {
                                res += eval_coef(pair.second._coef,i,j) * eval(pair.second._p->first,i,j) * eval(pair.second._p->second,i,j);
                            }
                        }
                        else {
                            for (size_t j = 0; j<dim; j++) {
                                res -= eval_coef(pair.second._coef,i,j) * eval(pair.second._p->first,i,j) * eval(pair.second._p->second,i,j);
                            }
                        }
                    }
                    else if (pair.second._coef_p1_tr) { // qterm = (coef*p1)^T*p2
                        //                        assert(pair.second._p->first->_dim[1]==1 && pair.second._coef->_dim[0]==pair.second._p->second->_dim[0]);
                        for (auto i = 0; i<pair.second._p->first->_dim[0]; i++) {
                            for (auto j = 0; j<pair.second._p->first->_dim[0]; j++) {
                                qval += eval_coef(pair.second._coef,i,j) * eval(pair.second._p->first,j) * eval(pair.second._p->second,i);
                            }
                        }
                    }
                    else if (pair.second._p->first->is_matrix() && !pair.second._p->second->is_matrix() && !pair.second._p->second->_is_transposed) {//matrix * vect
                        for (size_t j = 0; j<pair.second._p->second->_dim[0]; j++) {
                            qval += eval(pair.second._p->first,i,j) * eval(pair.second._p->second,j);
                        }
                        qval *= eval_coef(pair.second._coef,i);
                    }
                    else if (!pair.second._p->first->is_matrix() && pair.second._p->first->_is_transposed && pair.second._p->second->is_matrix() ) {//transposed vect * matrix
                        for (size_t j = 0; j<pair.second._p->first->_dim[0]; j++) {
                            qval += eval(pair.second._p->first,j) * eval(pair.second._p->second,j,i);
                        }
                        qval *= eval_coef(pair.second._coef,i);
                    }
                    else if (!pair.second._p->first->is_matrix() && pair.second._p->first->_is_transposed && !pair.second._p->second->is_matrix() && i==0) {//transposed vect * vec, a dot product of two vectors
                        for (size_t j = 0; j<pair.second._p->first->_dim[1]; j++) {
                            qval += eval(pair.second._p->first,j) * eval(pair.second._p->second,j);
                        }
                        qval *= eval_coef(pair.second._coef,i);
                    }
                    else if (!pair.second._coef->is_matrix() && pair.second._coef->_is_transposed && !pair.second._p->first->is_matrix()) {//transposed vect * vec, a dot product of two vectors
                        for (size_t j = 0; j<pair.second._p->first->_dim[0]; j++) {
                            qval += eval_coef(pair.second._coef,j) * eval(pair.second._p->first,j) * eval(pair.second._p->second,j);
                        }
                    }
                    else {
                        qval += eval_coef(pair.second._coef,i) * eval(pair.second._p->first,i) * eval(pair.second._p->second,i);
                    }
                    if (!pair.second._sign) {
                        qval *= -1.;
                    }
                    res += qval;
                }
            }
            //                res += eval_qterms(i);
            if(!_pterms->empty()){
                for (auto &pair:*_pterms) {
                    if (pair.second._coef->_is_transposed) {//transposed vect * vec, a dot product of two vectors
                        for (size_t j = 0; j<pair.second._coef->_dim[1]; j++) {
                            res += eval_pterm(pair.second,j);
                        }
                    }
                    else {
                        type pval = unit<type>().eval();
                        for (auto &vpair: *pair.second._l) {
                            pval *= std::pow(eval(vpair.first, i), vpair.second);
                        }
                        pval *= eval_coef(pair.second._coef,i);
                        if (!pair.second._sign) {
                            pval *= -1.;
                        }
                        res += pval;
                    }
                }
                //                res += eval_pterms(i);
            }
            if(_expr)
                res += eval_expr(_expr,i);
            if (func_is_number()) {
                assert(_val->size()>0);
                _val->at(0) = res;
                _evaluated = true;
            }
            else {
                //                if (is_constant() && i==_val->size()-1) {
                //                if (i==_val->size()-1) {
                //                    _evaluated = true;
                //                }
                assert(_val->size()>i);
                _val->at(i) = res;
            }
            return res;
        }
        
        inline type eval_cst(size_t i) {
            return eval_coef(_cst, i);
        }
        
        inline type eval_cst(size_t i, size_t j) {
            return eval_coef(_cst, i, j);
        }
        
        template<class T=type, typename enable_if<is_same<T, Cpx>::value>::type* = nullptr>
        inline type get_val(const shared_ptr<constant_>& c, size_t i=0) {
            switch (c->get_type()) {
                case binary_c:
                    return static_pointer_cast<constant<bool>>(c)->eval();
                    break;
                case short_c:
                    return static_pointer_cast<constant<short>>(c)->eval();
                    break;
                case integer_c:
                    return static_pointer_cast<constant<int>>(c)->eval();
                    break;
                case float_c:
                    return static_pointer_cast<constant<float>>(c)->eval();
                    break;
                case double_c:
                    return ((constant<double>*)(c.get()))->eval();
                    break;
                case long_c:
                    return static_pointer_cast<constant<long double>>(c)->eval();
                    break;
                case complex_c:
                    return static_pointer_cast<constant<Cpx>>(c)->eval();
                    break;
                case func_c:{
                    auto f = ((func_*)(c.get()));
                    switch (f->get_return_type()) {
                        case complex_:
                            return ((func<Cpx>*)(f))->_val->at(i);
                            break;
                        default:
                            break;
                    }
                    break;
                }
                case uexp_c:{
                    return eval_uexpr((uexpr<type>*)(c.get()),i);
                    break;
                }
                case bexp_c:{
                    return eval_bexpr((bexpr<type>*)(c.get()),i);
                    break;
                }
                default:{
                    auto p = ((param_*)(c.get()));
                    switch (p->get_intype()) {
                        case binary_:
                            return ((param<bool>*)(p))->eval(i);
                            break;
                        case short_:
                            return ((param<short>*)(p))->eval(i);
                            break;
                        case integer_:
                            return ((param<int>*)(p))->eval(i);
                            break;
                        case float_:
                            return ((param<float>*)(p))->eval(i);
                            break;
                        case double_:
                            return ((param<double>*)(p))->eval(i);
                            break;
                        case long_:
                            return ((param<long double>*)(p))->eval(i);
                            break;
                        case complex_:
                            return ((param<Cpx>*)(p))->eval(i);
                            break;
                        default:
                            break;
                    }
                    break;
                }
            }
            throw invalid_argument("Unsupported type");
        }
        
        template<class T=type, typename enable_if<is_arithmetic<T>::value>::type* = nullptr>
        inline type get_val(const shared_ptr<constant_>& c, size_t i=0) {
            switch (c->get_type()) {
                case binary_c:
                    return static_pointer_cast<constant<bool>>(c)->eval();
                    break;
                case short_c:
                    return static_pointer_cast<constant<short>>(c)->eval();
                    break;
                case integer_c:
                    return static_pointer_cast<constant<int>>(c)->eval();
                    break;
                case float_c:
                    return static_pointer_cast<constant<float>>(c)->eval();
                    break;
                case double_c:
                    return ((constant<double>*)(c.get()))->eval();
                    break;
                case long_c:
                    return static_pointer_cast<constant<long double>>(c)->eval();
                    break;
                case func_c:{
                    auto f = ((func_*)(c.get()));
                    switch (f->get_return_type()) {
                        case binary_:
                            if(f->func_is_number()){
                                return ((func<bool>*)(f))->_val->at(0);
                            }
                            return ((func<bool>*)(f))->_val->at(i);
                            break;
                        case short_:
                            if(f->func_is_number()){
                                return ((func<short>*)(f))->_val->at(0);
                            }
                            return ((func<short>*)(f))->_val->at(i);
                            break;
                        case integer_:
                            if(f->func_is_number()){
                                return ((func<int>*)(f))->_val->at(0);
                            }
                            return ((func<int>*)(f))->_val->at(i);
                            break;
                        case float_:
                            if(f->func_is_number()){
                                return ((func<float>*)(f))->_val->at(0);
                            }
                            return ((func<float>*)(f))->_val->at(i);
                            break;
                        case double_:
                            if(f->func_is_number()){
                                return ((func<double>*)(f))->_val->at(0);
                            }
                            return ((func<double>*)(f))->_val->at(i);
                            break;
                        case long_:
                            if(f->func_is_number()){
                                return ((func<long double>*)(f))->_val->at(0);
                            }
                            return ((func<long double>*)(f))->_val->at(i);
                            break;
                        default:
                            break;
                    }
                    break;
                }
                case uexp_c:{
                    return eval_uexpr(((uexpr<type>*)c.get()),i);
                    break;
                }
                case bexp_c:{
                    return eval_bexpr(((bexpr<type>*)c.get()),i);
                    break;
                }
                default:{
                    auto p = ((param_*)(c.get()));
                    switch (p->get_intype()) {
                        case binary_:
                            return ((param<bool>*)(p))->eval(i);
                            break;
                        case short_:
                            return ((param<short>*)(p))->eval(i);
                            break;
                        case integer_:
                            return ((param<int>*)(p))->eval(i);
                            break;
                        case float_:
                            return ((param<float>*)(p))->eval(i);
                            break;
                        case double_:
                            return ((param<double>*)(p))->eval(i);
                            break;
                        case long_:
                            return ((param<long double>*)(p))->eval(i);
                            break;
                        default:
                            break;
                    }
                    break;
                }
            }
            throw invalid_argument("Unsupported type");
        }
        
        template<class T=type, typename enable_if<is_same<T, Cpx>::value>::type* = nullptr>
        inline type get_val(const shared_ptr<constant_>& c, size_t i, size_t j) {
            switch (c->get_type()) {
                case binary_c:
                    return static_pointer_cast<constant<bool>>(c)->eval();
                    break;
                case short_c:
                    return static_pointer_cast<constant<short>>(c)->eval();
                    break;
                case integer_c:
                    return static_pointer_cast<constant<int>>(c)->eval();
                    break;
                case float_c:
                    return static_pointer_cast<constant<float>>(c)->eval();
                    break;
                case double_c:
                    return ((constant<double>*)(c.get()))->eval();
                    break;
                case long_c:
                    return static_pointer_cast<constant<long double>>(c)->eval();
                    break;
                case complex_c:
                    return static_pointer_cast<constant<Cpx>>(c)->eval();
                    break;
                case func_c:{
                    auto f = ((func_*)(c.get()));
                    switch (f->get_return_type()) {
                        case complex_:
                            return ((func<Cpx>*)(f))->get_val(i,j);
                            break;
                        default:
                            break;
                    }
                    break;
                }
                case uexp_c:{
                    return eval_uexpr((uexpr<type>*)(c.get()),i,j);
                    break;
                }
                case bexp_c:{
                    return eval_bexpr((bexpr<type>*)(c.get()),i,j);
                    break;
                }
                default:{
                    auto p = ((param_*)(c.get()));
                    switch (p->get_intype()) {
                        case binary_:
                            return ((param<bool>*)(p))->eval(i,j);
                            break;
                        case short_:
                            return ((param<short>*)(p))->eval(i,j);
                            break;
                        case integer_:
                            return ((param<int>*)(p))->eval(i,j);
                            break;
                        case float_:
                            return ((param<float>*)(p))->eval(i,j);
                            break;
                        case double_:
                            return ((param<double>*)(p))->eval(i,j);
                            break;
                        case long_:
                            return ((param<long double>*)(p))->eval(i,j);
                            break;
                        case complex_:
                            return ((param<Cpx>*)(p))->eval(i,j);
                            break;
                        default:
                            break;
                    }
                    break;
                }
            }
            throw invalid_argument("Unsupported type");
        }
        
        template<class T=type, typename enable_if<is_arithmetic<T>::value>::type* = nullptr>
        inline type get_val(const shared_ptr<constant_>& c, size_t i, size_t j) {
            switch (c->get_type()) {
                case binary_c:
                    return static_pointer_cast<constant<bool>>(c)->eval();
                    break;
                case short_c:
                    return static_pointer_cast<constant<short>>(c)->eval();
                    break;
                case integer_c:
                    return static_pointer_cast<constant<int>>(c)->eval();
                    break;
                case float_c:
                    return static_pointer_cast<constant<float>>(c)->eval();
                    break;
                case double_c:
                    return ((constant<double>*)(c.get()))->eval();
                    break;
                case long_c:
                    return static_pointer_cast<constant<long double>>(c)->eval();
                    break;
                case func_c:{
                    auto f = ((func_*)(c.get()));
                    switch (f->get_return_type()) {
                        case binary_:
                            return ((func<bool>*)(f))->get_val(i,j);
                            break;
                        case short_:
                            return ((func<short>*)(f))->get_val(i,j);
                            break;
                        case integer_:
                            return ((func<int>*)(f))->get_val(i,j);
                            break;
                        case float_:
                            return ((func<float>*)(f))->get_val(i,j);
                            break;
                        case double_:
                            return ((func<double>*)(f))->get_val(i,j);
                            break;
                        case long_:
                            return ((func<long double>*)(f))->get_val(i,j);
                            break;
                        default:
                            break;
                    }
                    break;
                }
                case uexp_c:{
                    return eval_uexpr(((uexpr<type>*)c.get()),i,j);
                    break;
                }
                case bexp_c:{
                    return eval_bexpr(((bexpr<type>*)c.get()),i,j);
                    break;
                }
                default:{
                    auto p = ((param_*)(c.get()));
                    switch (p->get_intype()) {
                        case binary_:
                            return ((param<bool>*)(p))->eval(i,j);
                            break;
                        case short_:
                            return ((param<short>*)(p))->eval(i,j);
                            break;
                        case integer_:
                            return ((param<int>*)(p))->eval(i,j);
                            break;
                        case float_:
                            return ((param<float>*)(p))->eval(i,j);
                            break;
                        case double_:
                            return ((param<double>*)(p))->eval(i,j);
                            break;
                        case long_:
                            return ((param<long double>*)(p))->eval(i,j);
                            break;
                        default:
                            break;
                    }
                    break;
                }
            }
            throw invalid_argument("Unsupported type");
        }
        
        template<class T=type, typename enable_if<is_arithmetic<T>::value>::type* = nullptr>
        inline type eval(const shared_ptr<constant_>& c, size_t i=0) {
            switch (c->get_type()) {
                case binary_c:
                    return static_pointer_cast<constant<bool>>(c)->eval();
                    break;
                case short_c:
                    return static_pointer_cast<constant<short>>(c)->eval();
                    break;
                case integer_c:
                    return static_pointer_cast<constant<int>>(c)->eval();
                    break;
                case float_c:
                    return static_pointer_cast<constant<float>>(c)->eval();
                    break;
                case double_c:
                    return ((constant<double>*)(c.get()))->eval();
                    break;
                case long_c:
                    return static_pointer_cast<constant<long double>>(c)->eval();
                    break;
                case func_c:{
                    auto f = ((func_*)(c.get()));
                    switch (f->get_return_type()) {
                        case binary_:
                            return ((func<bool>*)(f))->eval(i);
                            break;
                        case short_:
                            return ((func<short>*)(f))->eval(i);
                            break;
                        case integer_:
                            return ((func<int>*)(f))->eval(i);
                            break;
                        case float_:
                            return ((func<float>*)(f))->eval(i);
                            break;
                        case double_:
                            return ((func<double>*)(f))->eval(i);
                            break;
                        case long_:
                            return ((func<long double>*)(f))->eval(i);
                            break;
                        default:
                            break;
                    }
                    break;
                }
                case uexp_c:{
                    return eval_uexpr(((uexpr<type>*)c.get()),i);
                    break;
                }
                case bexp_c:{
                    return eval_bexpr(((bexpr<type>*)c.get()),i);
                    break;
                }
                default:{
                    auto p = ((param_*)(c.get()));
                    switch (p->get_intype()) {
                        case binary_:
                            return ((param<bool>*)(p))->eval(i);
                            break;
                        case short_:
                            return ((param<short>*)(p))->eval(i);
                            break;
                        case integer_:
                            return ((param<int>*)(p))->eval(i);
                            break;
                        case float_:
                            return ((param<float>*)(p))->eval(i);
                            break;
                        case double_:
                            return ((param<double>*)(p))->eval(i);
                            break;
                        case long_:
                            return ((param<long double>*)(p))->eval(i);
                            break;
                        default:
                            break;
                    }
                    break;
                }
            }
            throw invalid_argument("Unsupported type");
        }
        
        
        template<class T=type, typename enable_if<is_arithmetic<T>::value>::type* = nullptr>
        inline type eval(const shared_ptr<constant_>& c, size_t i, size_t j) {
            switch (c->get_type()) {
                case binary_c:
                    return static_pointer_cast<constant<bool>>(c)->eval();
                    break;
                case short_c:
                    return static_pointer_cast<constant<short>>(c)->eval();
                    break;
                case integer_c:
                    return static_pointer_cast<constant<int>>(c)->eval();
                    break;
                case float_c:
                    return static_pointer_cast<constant<float>>(c)->eval();
                    break;
                case double_c:
                    return ((constant<double>*)(c.get()))->eval();
                    break;
                case long_c:
                    return static_pointer_cast<constant<long double>>(c)->eval();
                    break;
                case func_c:{
                    auto f = ((func_*)(c.get()));
                    switch (f->get_return_type()) {
                        case binary_:
                            return ((func<bool>*)(f))->eval(i,j);
                            break;
                        case short_:
                            return ((func<short>*)(f))->eval(i,j);
                            break;
                        case integer_:
                            return ((func<int>*)(f))->eval(i,j);
                            break;
                        case float_:
                            return ((func<float>*)(f))->eval(i,j);
                            break;
                        case double_:
                            return ((func<double>*)(f))->eval(i,j);
                            break;
                        case long_:
                            return ((func<long double>*)(f))->eval(i,j);
                            break;
                        default:
                            break;
                    }
                    break;
                }
                case uexp_c:{
                    return eval_uexpr((uexpr<type>*)(c.get()),i,j);
                    break;
                }
                case bexp_c:{
                    return eval_bexpr((bexpr<type>*)(c.get()),i,j);
                    break;
                }
                default:{
                    auto p = ((param_*)(c.get()));
                    switch (p->get_intype()) {
                        case binary_:
                            return ((param<bool>*)(p))->eval(i,j);
                            break;
                        case short_:
                            return ((param<short>*)(p))->eval(i,j);
                            break;
                        case integer_:
                            return ((param<int>*)(p))->eval(i,j);
                            break;
                        case float_:
                            return ((param<float>*)(p))->eval(i,j);
                            break;
                        case double_:
                            return ((param<double>*)(p))->eval(i,j);
                            break;
                        case long_:
                            return ((param<long double>*)(p))->eval(i,j);
                            break;
                        default:
                            break;
                    }
                    break;
                }
            }
            throw invalid_argument("Unsupported type");
        }
        
        template<class T=type, typename enable_if<is_same<T, Cpx>::value>::type* = nullptr>
        inline type eval(const shared_ptr<constant_>& c, size_t i=0) {
            switch (c->get_type()) {
                case binary_c:
                    return static_pointer_cast<constant<bool>>(c)->eval();
                    break;
                case short_c:
                    return static_pointer_cast<constant<short>>(c)->eval();
                    break;
                case integer_c:
                    return static_pointer_cast<constant<int>>(c)->eval();
                    break;
                case float_c:
                    return static_pointer_cast<constant<float>>(c)->eval();
                    break;
                case double_c:
                    return ((constant<double>*)(c.get()))->eval();
                    break;
                case long_c:
                    return static_pointer_cast<constant<long double>>(c)->eval();
                    break;
                case complex_c:
                    return static_pointer_cast<constant<Cpx>>(c)->eval();
                    break;
                case func_c:{
                    auto f = ((func_*)(c.get()));
                    switch (f->get_return_type()) {
                        case binary_:
                            return ((func<bool>*)(f))->eval(i);
                            break;
                        case short_:
                            return ((func<short>*)(f))->eval(i);
                            break;
                        case integer_:
                            return ((func<int>*)(f))->eval(i);
                            break;
                        case float_:
                            return ((func<float>*)(f))->eval(i);
                            break;
                        case double_:
                            return ((func<double>*)(f))->eval(i);
                            break;
                        case long_:
                            return ((func<long double>*)(f))->eval(i);
                            break;
                        case complex_:
                            return ((func<Cpx>*)(f))->eval(i);
                            break;
                        default:
                            break;
                    }
                    break;
                }
                case uexp_c:{
                    return eval_uexpr((uexpr<type>*)(c.get()),i);
                    break;
                }
                case bexp_c:{
                    return eval_bexpr((bexpr<type>*)(c.get()),i);
                    break;
                }
                default:{
                    auto p = ((param_*)(c.get()));
                    switch (p->get_intype()) {
                        case binary_:
                            return ((param<bool>*)(p))->eval(i);
                            break;
                        case short_:
                            return ((param<short>*)(p))->eval(i);
                            break;
                        case integer_:
                            return ((param<int>*)(p))->eval(i);
                            break;
                        case float_:
                            return ((param<float>*)(p))->eval(i);
                            break;
                        case double_:
                            return ((param<double>*)(p))->eval(i);
                            break;
                        case long_:
                            return ((param<long double>*)(p))->eval(i);
                            break;
                        case complex_:
                            return ((param<Cpx>*)(p))->eval(i);
                            break;
                        default:
                            break;
                    }
                    break;
                }
            }
            throw invalid_argument("Unsupported type");
        }
        
        template<class T=type, typename enable_if<is_same<T, Cpx>::value>::type* = nullptr>
        inline type eval(const shared_ptr<constant_>& c, size_t i, size_t j) {
            switch (c->get_type()) {
                case binary_c:
                    return static_pointer_cast<constant<bool>>(c)->eval();
                    break;
                case short_c:
                    return static_pointer_cast<constant<short>>(c)->eval();
                    break;
                case integer_c:
                    return static_pointer_cast<constant<int>>(c)->eval();
                    break;
                case float_c:
                    return static_pointer_cast<constant<float>>(c)->eval();
                    break;
                case double_c:
                    return ((constant<double>*)(c.get()))->eval();
                    break;
                case long_c:
                    return static_pointer_cast<constant<long double>>(c)->eval();
                    break;
                case complex_c:
                    return static_pointer_cast<constant<Cpx>>(c)->eval();
                    break;
                case func_c:{
                    auto f = ((func_*)(c.get()));
                    switch (f->get_return_type()) {
                        case binary_:
                            return ((func<bool>*)(f))->eval(i,j);
                            break;
                        case short_:
                            return ((func<short>*)(f))->eval(i,j);
                            break;
                        case integer_:
                            return ((func<int>*)(f))->eval(i,j);
                            break;
                        case float_:
                            return ((func<float>*)(f))->eval(i,j);
                            break;
                        case double_:
                            return ((func<double>*)(f))->eval(i,j);
                            break;
                        case long_:
                            return ((func<long double>*)(f))->eval(i,j);
                            break;
                        case complex_:
                            return ((func<Cpx>*)(f))->eval(i,j);
                            break;
                        default:
                            break;
                    }
                    break;
                }
                case uexp_c:{
                    return eval(static_pointer_cast<uexpr<type>>(c),i,j);
                    break;
                }
                case bexp_c:{
                    return eval(static_pointer_cast<bexpr<type>>(c),i,j);
                    break;
                }
                default:{
                    auto p = ((param_*)(c.get()));
                    switch (p->get_intype()) {
                        case binary_:
                            return ((param<bool>*)(p))->eval(i,j);
                            break;
                        case short_:
                            return ((param<short>*)(p))->eval(i,j);
                            break;
                        case integer_:
                            return ((param<int>*)(p))->eval(i,j);
                            break;
                        case float_:
                            return ((param<float>*)(p))->eval(i,j);
                            break;
                        case double_:
                            return ((param<double>*)(p))->eval(i,j);
                            break;
                        case long_:
                            return ((param<long double>*)(p))->eval(i,j);
                            break;
                        case complex_:
                            return ((param<Cpx>*)(p))->eval(i,j);
                            break;
                        default:
                            break;
                    }
                    break;
                }
            }
            throw invalid_argument("Unsupported type");
        }
        
        inline type eval_coef(const shared_ptr<constant_>& coef, size_t i) {
            auto coef_type = coef->get_type();
            if (coef_type==func_c) {
                auto f_cst = ((func<type>*)(coef.get()));
                return f_cst->eval(i);
            }
            else if(coef_type==par_c || coef_type==var_c) {
                auto p_cst = ((param<type>*)(coef.get()));
                return p_cst->eval(i);
            }
            else {
                auto p_cst = ((constant<type>*)(coef.get()));
                return p_cst->eval();
            }
            throw invalid_argument("in function eval_coef(shared_ptr<constant_> coef, size_t i), coef should be a constant");
        }
        
        inline type eval_coef(const shared_ptr<constant_>& coef, size_t i, size_t j) {
            auto coef_type = coef->get_type();
            if (coef_type==func_c) {
                auto f_cst = ((func<type>*)(coef.get()));
                return f_cst->eval(i,j);
            }
            else if(coef_type==par_c || coef_type==var_c) {
                auto p_cst = ((param<type>*)(coef.get()));
                return p_cst->eval(i,j);
            }
            else {
                auto p_cst = ((constant<type>*)(coef.get()));
                return p_cst->eval();
            }
            throw invalid_argument("in function eval_coef(shared_ptr<constant_> coef, size_t i), coef should be a constant");
        }
        
        type eval_lterm(const lterm& lt, size_t i){
            type res = zero<type>().eval();
            if ((lt._coef->_is_transposed || lt._coef->is_matrix() || (lt._p->is_indexed() && lt._p->_indices->_ids->size()>1)) && !lt._p->is_matrix()) {
                auto dim = lt._p->get_dim(i);
                if (lt._sign) {
                    for (size_t j = 0; j<dim; j++) {
                        res += eval_coef(lt._coef,i,j) * eval(lt._p,i,j);
                    }
                }
                else {
                    for (size_t j = 0; j<dim; j++) {
                        res -= eval_coef(lt._coef,i,j) * eval(lt._p,i,j);
                    }
                }
            }
            else {
                if (lt._sign) {
                    res += eval_coef(lt._coef,i) * eval(lt._p, i);
                }
                else {
                    res -= eval_coef(lt._coef,i) * eval(lt._p, i);
                }
            }
            return res;
        }
        
        type eval_lterm(const lterm& lt, size_t i, size_t j){
            type res = zero<type>().eval();
            if (lt._coef->is_matrix() && lt._p->is_matrix()) {
                //matrix product
                if(lt._sign){
                    for (size_t col = 0; col<lt._coef->_dim[1]; col++) {
                        res += eval_coef(lt._coef, i,col) * eval(lt._p,col,j);
                    }
                }
                else {
                    for (size_t col = 0; col<lt._coef->_dim[1]; col++) {
                        res -= eval_coef(lt._coef, i,col) * eval(lt._p,col,j);
                    }
                }
                return res;
            }
            if (lt._coef->is_matrix() && !lt._p->is_matrix() && lt._p->_is_transposed) {//matrix * transposed vect
                if(lt._sign){
                    return eval_coef(lt._coef, i,j) * eval(lt._p,j);
                }
                else {
                    return res -= eval_coef(lt._coef, i,j) * eval(lt._p,j);
                }
            }
            
            if (!lt._coef->is_matrix() && !lt._coef->_is_transposed && lt._p->is_matrix()) {//vect * matrix
                if(lt._sign) {
                    return eval_coef(lt._coef, i) * eval(lt._p,i,j);
                }
                else {
                    return res -= eval_coef(lt._coef, i) * eval(lt._p,i,j);
                }
            }
            if (lt._coef->is_matrix() && lt._p->_is_vector) {//matrix*vect
                if(lt._sign) {
                    return eval_coef(lt._coef, i,j) * eval(lt._p,i);
                }
                else {
                    return res -= eval_coef(lt._coef, i,j) * eval(lt._p,i);
                }
            }
            if(lt._sign) {
                return eval_coef(lt._coef, i,j) * eval(lt._p,i,j);
            }
            else {
                return res -= eval_coef(lt._coef, i,j) * eval(lt._p,i,j);
            }
        }
        
        
        type eval_qterm(const qterm& qt, size_t i){
            type res = zero<type>().eval();
            if (qt._coef_p1_tr) { // qterm = (coef*p1)^T*p2
                assert(qt._p->first->_dim[1]==1 && qt._coef->_dim[0]==qt._p->second->_dim[0]);
                for (auto i = 0; i<qt._p->first->_dim[0]; i++) {
                    for (auto j = 0; j<qt._p->first->_dim[0]; j++) {
                        res += eval_coef(qt._coef,i,j) * eval(qt._p->first,i) * eval(qt._p->second,j);
                    }
                }
                if (!qt._sign) {
                    res *= -1;
                }
                return res;
                
            }
            if (qt._p->first->is_matrix() && !qt._p->second->is_matrix() && !qt._p->second->_is_transposed) {//matrix * vect
                for (size_t j = 0; j<qt._p->second->_dim[0]; j++) {
                    res += eval(qt._p->first,i,j) * eval(qt._p->second,j);
                }
                res *= eval_coef(qt._coef,i);
            }
            else if (!qt._p->first->is_matrix() && qt._p->first->_is_transposed && qt._p->second->is_matrix() ) {//transposed vect * matrix
                for (size_t j = 0; j<qt._p->first->_dim[0]; j++) {
                    res += eval(qt._p->first,j) * eval(qt._p->second,j,i);
                }
                res *= eval_coef(qt._coef,i);
            }
            else if (!qt._p->first->is_matrix() && qt._p->first->_is_transposed && !qt._p->second->is_matrix() && i==0) {//transposed vect * vec, a dot product of two vectors
                for (size_t j = 0; j<qt._p->first->_dim[1]; j++) {
                    res += eval(qt._p->first,j) * eval(qt._p->second,j);
                }
                res *= eval_coef(qt._coef,i);
            }
            else if (!qt._coef->is_matrix() && qt._coef->_is_transposed && !qt._p->first->is_matrix()) {//transposed vect * vec, a dot product of two vectors
                for (size_t j = 0; j<qt._p->first->_dim[0]; j++) {
                    res += eval_coef(qt._coef,j) * eval(qt._p->first,j) * eval(qt._p->second,j);
                }
            }
            else {
                res = eval_coef(qt._coef,i) * eval(qt._p->first,i) * eval(qt._p->second,i);
            }
            if (!qt._sign) {
                res *= -1;
            }
            return res;
        }
        
        type eval_qterm(const qterm& qt, size_t i, size_t j){
            type res = zero<type>().eval();
            if (qt._p->first->is_matrix() && qt._p->second->is_matrix()) {
                //matrix product
                for (size_t col = 0; col<qt._p->first->_dim[1]; col++) {
                    res += eval(qt._p->first,i,col) * eval(qt._p->second,col,j);
                }
                res *= eval_coef(qt._coef,i);
            }
            else if (qt._p->first->is_matrix() && !qt._p->second->is_matrix() && qt._p->second->_is_transposed) {//matrix * transposed vect
                res = eval_coef(qt._coef,i) * eval(qt._p->first,i,j) * eval(qt._p->second,j);
            }
            else if (!qt._p->first->is_matrix() && !qt._p->first->_is_transposed && qt._p->second->is_matrix() ) {//vect * matrix
                res = eval_coef(qt._coef,i) * eval(qt._p->first,i) * eval(qt._p->second,i,j);
            }
            else {
                throw invalid_argument("eval(i,j) on non-matrix function");
            }
            if (!qt._sign) {
                res *= -1;
            }
            return res;
        }
        
        type eval_pterm(const pterm& pt, size_t i){
            type res = zero<type>().eval();
            if (pt._coef->_is_transposed) {
                throw invalid_argument("Unspported operation\n");
            }// TREAT TRANSPOSED VECTORS IN POLYNOMIAL TERMS HERE
            else {
                res += 1;
                for (auto &pair: *pt._l) {
                    res *= pow(eval(pair.first, i), pair.second);
                }
                res *= eval_coef(pt._coef,i);
            }
            if (!pt._sign) {
                res *= -1;
            }
            return res;
        }
        
        type eval_pterm(const pterm& pt, size_t i, size_t j){
            type res = zero<type>().eval();
            res += 1;
            for (auto &pair: *pt._l) {
                res *= pow(eval(pair.first,i,j), pair.second);
            }
            
            res *= eval_coef(pt._coef,i,j);
            if (!pt._sign) {
                res *= -1;
            }
            return res;
        }
        
        type eval_lterms(size_t i) {
            type res = zero<type>().eval();
            auto it = _lterms->begin();
            while(it!=_lterms->end()){
                res += eval_lterm(it->second,i);
                it++;
            }
            return res;
        }
        
        type eval_qterms(size_t i) {
            type res = zero<type>().eval();
            auto it = _qterms->begin();
            while(it!=_qterms->end()){
                res += eval_qterm(it->second,i);
                it++;
            }
            return res;
        }
        
        type eval_pterms(size_t i) {
            type res = zero<type>().eval();
            auto it = _pterms->begin();
            while(it!=_pterms->end()){
                res += eval_pterm(it->second,i);
                it++;
            }
            return res;
        }
        
        bool is_evaluated() const{
            return _evaluated;
        }
        
        void evaluate(bool v){
            _evaluated = v;
        }
        
        inline type eval_expr(const shared_ptr<expr<type>>& exp, size_t i) {
            if (exp->is_uexpr()) {
                return eval_uexpr(((uexpr<type>*)exp.get()),i);
            }
            return eval_bexpr(((bexpr<type>*)exp.get()),i);
        }
        
        inline type eval_expr(const shared_ptr<expr<type>>& exp, size_t i, size_t j) {
            if (exp->is_uexpr()) {
                return eval_uexpr(((uexpr<type>*)exp.get()),i, j);
            }
            return eval_bexpr(((bexpr<type>*)exp.get()),i,j);
        }
        
        template<class T=type, typename enable_if<is_arithmetic<T>::value>::type* = nullptr>
        inline T eval_uexpr(uexpr<type>* exp, size_t i) {
            if (exp->_son->is_constant() && !exp->_son->is_evaluated()) {
                exp->_son->eval_all();
            }
            T res = get_val(exp->_son,i);
            switch (exp->_otype) {
                case cos_:
                    return exp->_coef*std::cos(res);
                    break;
                case sin_:
                    return exp->_coef*std::sin(res);
                    break;
                case tan_:
                    return exp->_coef*std::tan(res);
                    break;
                case atan_:
                    return exp->_coef*std::atan(res);
                    break;
                case acos_:
                    return exp->_coef*std::acos(res);
                    break;
                case asin_:
                    return exp->_coef*std::asin(res);
                    break;
                case sqrt_:
                    return exp->_coef*std::sqrt(res);
                    break;
                case log_:
                    return exp->_coef*std::log(res);
                    break;
                case exp_:
                    return exp->_coef*std::exp(res);
                    break;
                case relu_:{
                    if(res < 0)
                        res = 0;
                    return exp->_coef*res;
                }
                case unit_step_:{
                    if(res <= 0)
                        return 0;
                    return exp->_coef*1;
                }
                    break;
                default:
                    throw invalid_argument("Unsupported unary operator");
                    break;
            }
        }
        template<class T=type, typename enable_if<is_same<T, Cpx>::value>::type* = nullptr>
        Cpx eval_uexpr(uexpr<T>* exp, size_t i) {
            if (exp->_son->is_constant() && !exp->_son->is_evaluated()) {
                exp->_son->eval_all();
            }
            Cpx res = get_val(exp->_son,i);
            switch (exp->_otype) {
                case cos_:
                    return exp->_coef*std::cos(res);
                    break;
                case sin_:
                    return exp->_coef*std::sin(res);
                    break;
                case tan_:
                    return exp->_coef*std::tan(res);
                    break;
                case atan_:
                    return exp->_coef*std::atan(res);
                    break;
                case acos_:
                    return exp->_coef*std::acos(res);
                    break;
                case asin_:
                    return exp->_coef*std::asin(res);
                    break;
                case sqrt_:
                    return exp->_coef*std::sqrt(res);
                    break;
                case log_:
                    return exp->_coef*std::log(res);
                    break;
                case exp_:
                    return exp->_coef*std::exp(res);
                    break;
                case relu_:{
                    if(res.real() < 0)
                        res.real(0);
                    if(res.imag() < 0)
                        res.imag(0);
                    return exp->_coef*res;
                }
                case unit_step_:{
                    if(res.real() <= 0 || res.imag() <= 0)
                        return Cpx(0,0);
                    return exp->_coef*Cpx(1,0);
                }
                    break;
                default:
                    throw invalid_argument("Unsupported unary operator");
                    break;
            }
        }
        
        template<class T=type, typename enable_if<is_same<T, Cpx>::value>::type* = nullptr>
        inline T  eval_bexpr(bexpr<type>* exp, size_t i){
            if (exp->_lson->is_constant() && !exp->_lson->is_evaluated()) {
                exp->_lson->eval_all();
            }
            if (exp->_rson->is_constant() && !exp->_rson->is_evaluated()) {
                exp->_rson->eval_all();
            }
            if(exp->_otype==product_ && (exp->_lson->is_matrix_indexed() || exp->_rson->is_matrix_indexed()))
            {
                auto dim = exp->_lson->get_dim(i);
                if(exp->_rson->is_matrix_indexed()){
                    dim = exp->_rson->get_dim(i);
                }
                if(dim==0){
                    return 0;
                }
                T res = 0.;
                for (auto idx = 0; idx <dim; idx++) {
                    T lval = get_val(exp->_lson,i,idx);
                    T rval = get_val(exp->_rson,i,idx);
                    res += exp->_coef*(lval*rval);
                }
                return res;
            }
            T lval = get_val(exp->_lson,i);
            T rval = get_val(exp->_rson,i);
            switch (exp->_otype) {
                case plus_:
                    return exp->_coef*(lval + rval);
                    break;
                case minus_:
                    return exp->_coef*(lval - rval);
                    break;
                case product_:
                    return exp->_coef*(lval*rval);
                    break;
                case div_:
                    return exp->_coef*(lval/rval);
                    break;
                case min_:
                    return exp->_coef*(gravity::min(lval,rval));
                    break;
                case max_:
                    return exp->_coef*(gravity::max(lval,rval));
                    break;
                default:
                    throw invalid_argument("Unsupported binary operator");
                    break;
            }
            
        }
        
        template<class T=type, typename enable_if<is_arithmetic<T>::value>::type* = nullptr>
        inline T  eval_bexpr(bexpr<type>* exp, size_t i){
            if (exp->_lson->is_constant() && !exp->_lson->is_evaluated()) {
                exp->_lson->eval_all();
            }
            if (exp->_rson->is_constant() && !exp->_rson->is_evaluated()) {
                exp->_rson->eval_all();
            }
            if(exp->_otype==product_ && (exp->_lson->is_matrix_indexed() || exp->_rson->is_matrix_indexed()))
            {
                auto dim = exp->_lson->get_dim(i);
                if(exp->_rson->is_matrix_indexed()){
                    dim = exp->_rson->get_dim(i);
                }
                if(dim==0){
                    return 0;
                }
                T res = 0.;
                for (auto idx = 0; idx <dim; idx++) {
                    T lval = get_val(exp->_lson,i,idx);
                    T rval = get_val(exp->_rson,i,idx);
                    res += exp->_coef*(lval*rval);
                }
                return res;
            }
            T lval = get_val(exp->_lson,i);
            T rval = get_val(exp->_rson,i);
            switch (exp->_otype) {
                case plus_:
                    return exp->_coef*(lval + rval);
                    break;
                case minus_:
                    return exp->_coef*(lval - rval);
                    break;
                case product_:
                    return exp->_coef*(lval*rval);
                    break;
                case div_:
                    return exp->_coef*(lval/rval);
                    break;
                case power_:
                    return exp->_coef*(powl(lval,rval));
                    break;
                case min_:
                    return exp->_coef*(gravity::min(lval,rval));
                    break;
                case max_:
                    return exp->_coef*(gravity::max(lval,rval));
                    break;
                default:
                    throw invalid_argument("Unsupported binary operator");
                    break;
            }
            
        }
        
        
        template<class T=type, typename enable_if<is_arithmetic<T>::value>::type* = nullptr>
        T eval_uexpr(const uexpr<type>* exp, size_t i, size_t j) {
            T res = eval(exp->_son,i,j);
            switch (exp->_otype) {
                case cos_:
                    return exp->_coef*std::cos(res);
                    break;
                case sin_:
                    return exp->_coef*std::sin(res);
                    break;
                case acos_:
                    return exp->_coef*std::acos(res);
                    break;
                case asin_:
                    return exp->_coef*std::asin(res);
                    break;
                case sqrt_:
                    return exp->_coef*std::sqrt(res);
                    break;
                case log_:
                    return exp->_coef*std::log(res);
                    break;
                case exp_:
                    return exp->_coef*std::exp(res);
                    break;
                case relu_:{
                    if(res < 0)
                        res = 0;
                    return exp->_coef*res;
                }
                case unit_step_:{
                    if(res <= 0)
                        return 0;
                    return exp->_coef*1;
                }
                    break;
                default:
                    throw invalid_argument("Unsupported unary operator");
                    break;
            }
        }
        template<class T=type, typename enable_if<is_same<T, Cpx>::value>::type* = nullptr>
        Cpx eval_uexpr(const uexpr<type>* exp, size_t i, size_t j) {
            Cpx res = eval(exp->_son,i,j);
            switch (exp->_otype) {
                case cos_:
                    return exp->_coef*std::cos(res);
                    break;
                case sin_:
                    return exp->_coef*std::sin(res);
                    break;
                case acos_:
                    return exp->_coef*std::acos(res);
                    break;
                case asin_:
                    return exp->_coef*std::asin(res);
                    break;
                case sqrt_:
                    return exp->_coef*std::sqrt(res);
                    break;
                case log_:
                    return exp->_coef*std::log(res);
                    break;
                case exp_:
                    return exp->_coef*std::exp(res);
                    break;
                case relu_:{
                    if(res.real() < 0)
                        res.real(0);
                    if(res.imag() < 0)
                        res.imag(0);
                    return exp->_coef*res;
                }
                case unit_step_:{
                    if(res.real() <= 0 || res.imag() <= 0)
                        return Cpx(0,0);
                    return exp->_coef*Cpx(1,0);
                }
                default:
                    throw invalid_argument("Unsupported unary operator");
                    break;
            }
        }
        
        template<class T=type, typename enable_if<is_same<T, Cpx>::value>::type* = nullptr>
        T  eval_bexpr(const bexpr<type>* exp, size_t i, size_t j){
            if (exp->_lson->is_constant() && !exp->_lson->is_evaluated()) {
                exp->_lson->eval_all();
            }
            if (exp->_rson->is_constant() && !exp->_rson->is_evaluated()) {
                exp->_rson->eval_all();
            }
            T lval = eval(exp->_lson,i,j);
            T rval = eval(exp->_rson,i,j);
            switch (exp->_otype) {
                case plus_:
                    return exp->_coef*(lval + rval);
                    break;
                case minus_:
                    return exp->_coef*(lval - rval);
                    break;
                case product_:
                    return exp->_coef*(lval*rval);
                    break;
                case div_:
                    return exp->_coef*(lval/rval);
                    break;
                case min_:
                    return exp->_coef*(gravity::min(lval,rval));
                    break;
                case max_:
                    return exp->_coef*(gravity::max(lval,rval));
                    break;
                default:
                    throw invalid_argument("Unsupported binary operator");
                    break;
            }
            
        }
        
        template<class T=type, typename enable_if<is_arithmetic<T>::value>::type* = nullptr>
        T  eval_bexpr(const bexpr<type>* exp, size_t i, size_t j){
            if (exp->_lson->is_constant() && !exp->_lson->is_evaluated()) {
                exp->_lson->eval_all();
            }
            if (exp->_rson->is_constant() && !exp->_rson->is_evaluated()) {
                exp->_rson->eval_all();
            }
            T lval = eval(exp->_lson,i,j);
            T rval = eval(exp->_rson,i,j);
            switch (exp->_otype) {
                case plus_:
                    return exp->_coef*(lval + rval);
                    break;
                case minus_:
                    return exp->_coef*(lval - rval);
                    break;
                case product_:
                    return exp->_coef*(lval*rval);
                    break;
                case div_:
                    return exp->_coef*(lval/rval);
                    break;
                case power_:
                    return exp->_coef*(powl(lval,rval));
                    break;
                case min_:
                    return exp->_coef*(gravity::min(lval,rval));
                    break;
                case max_:
                    return exp->_coef*(gravity::max(lval,rval));
                    break;
                default:
                    throw invalid_argument("Unsupported binary operator");
                    break;
            }
            
        }
        
        
        type eval(const string& key) {
            return _val->at(_indices->_keys_map->at(key));
        }
        
        type eval(size_t i, size_t j){
            if(is_zero()){
                return 0.;
            }
            //            if (is_constant() && _evaluated) {
            if (_evaluated) {
                if (func_is_number()){
                    return _val->at(0);
                }
                if (is_indexed() && _indices->_ids->size()>1) {
                    if (_indices->_ids->at(i).at(j) >= _val->size()) {
                        throw invalid_argument("eval(i,j): out of range");
                    }
                    return _val->at(_indices->_ids->at(i).at(j));
                }
                if (!is_matrix()) {
                    return eval(j);
                }
                if (_is_transposed) {
                    return _val->at(j*_dim[0]+i);
                }
                return _val->at(i*_dim[1]+j);
            }
            type res = 0;
            if(!_cst->is_zero())
                res += eval_cst(i,j);
            if(!_lterms->empty()){
                for (auto &pair:*_lterms) {
                    if ((pair.second._coef->_is_transposed || pair.second._coef->_is_transposed || pair.second._coef->is_matrix()) && !pair.second._p->is_matrix()) {
                        auto dim = pair.second._p->get_dim(i);
                        if (pair.second._sign) {
                            for (size_t j = 0; j<dim; j++) {
                                res += eval_coef(pair.second._coef,i,j) * eval(pair.second._p,i,j);
                            }
                        }
                        else {
                            for (size_t j = 0; j<dim; j++) {
                                res -= eval_coef(pair.second._coef,i,j) * eval(pair.second._p,i,j);
                            }
                        }
                    }
                    else {
                        if (pair.second._sign) {
                            res += eval_coef(pair.second._coef,i,j) * eval(pair.second._p,i,j);
                        }
                        else {
                            res -= eval_coef(pair.second._coef,i,j) * eval(pair.second._p,i,j);
                        }
                    }
                }
            }
            //                res += eval_lterms(i);
            if(!_qterms->empty()){
                for (auto &pair:*_qterms) {
                    type qval = 0;
                    if (pair.second._coef_p1_tr) { // qterm = (coef*p1)^T*p2
                        assert(pair.second._p->first->_dim[1]==1 && pair.second._coef->_dim[0]==pair.second._p->second->_dim[0]);
                        for (auto i = 0; i<pair.second._p->first->_dim[0]; i++) {
                            for (auto j = 0; j<pair.second._p->first->_dim[0]; j++) {
                                qval += eval_coef(pair.second._coef,i,j) * eval(pair.second._p->first,i) * eval(pair.second._p->second,j);
                            }
                        }
                    }
                    else if (pair.second._p->first->is_matrix() && !pair.second._p->second->is_matrix() && !pair.second._p->second->_is_transposed) {//matrix * vect
                        for (size_t j = 0; j<pair.second._p->second->_dim[0]; j++) {
                            qval += eval(pair.second._p->first,i,j) * eval(pair.second._p->second,j);
                        }
                        qval *= eval_coef(pair.second._coef,i);
                    }
                    else if (!pair.second._p->first->is_matrix() && pair.second._p->first->_is_transposed && pair.second._p->second->is_matrix() ) {//transposed vect * matrix
                        for (size_t j = 0; j<pair.second._p->first->_dim[0]; j++) {
                            qval += eval(pair.second._p->first,j) * eval(pair.second._p->second,j,i);
                        }
                        qval *= eval_coef(pair.second._coef,i);
                    }
                    else if (!pair.second._p->first->is_matrix() && pair.second._p->first->_is_transposed && !pair.second._p->second->is_matrix() && i==0) {//transposed vect * vec, a dot product of two vectors
                        for (size_t j = 0; j<pair.second._p->first->_dim[1]; j++) {
                            qval += eval(pair.second._p->first,j) * eval(pair.second._p->second,j);
                        }
                        qval *= eval_coef(pair.second._coef,i);
                    }
                    else if (!pair.second._coef->is_matrix() && pair.second._coef->_is_transposed && !pair.second._p->first->is_matrix()) {//transposed vect * vec, a dot product of two vectors
                        for (size_t j = 0; j<pair.second._p->first->_dim[0]; j++) {
                            qval += eval_coef(pair.second._coef,j) * eval(pair.second._p->first,j) * eval(pair.second._p->second,j);
                        }
                    }
                    else {
                        qval += eval_coef(pair.second._coef,i,j) * eval(pair.second._p->first,i,j) * eval(pair.second._p->second,i,j);
                    }
                    if (!pair.second._sign) {
                        qval *= -1.;
                    }
                    res += qval;
                }
            }
            //                res += eval_qterms(i);
            if(!_pterms->empty()){
                for (auto &pair:*_pterms) {
                    type pval = unit<type>().eval();
                    for (auto &vpair: *pair.second._l) {
                        pval *= std::pow(eval(vpair.first, i,j), vpair.second);
                    }
                    pval *= eval_coef(pair.second._coef,i,j);
                    if (!pair.second._sign) {
                        pval *= -1.;
                    }
                    res += pval;
                }
                //                res += eval_pterms(i);
            }
            if(_expr)
                res += eval_expr(_expr,i,j);
            if (func_is_number()) {
                _val->at(0) = res;
                _evaluated = true;
            }
            else {
                //                if (is_constant() && i==_val->size()-1) {
                if(is_matrix_indexed()){
                    _val->at(_indices->_ids->at(i).at(j)) = res;
                }
                else{
                    if (_is_transposed) {
                        if(j*_dim[0]+i>=_val->size()){
                            throw invalid_argument("out of range assignment in eval");
                        }
                        _val->at(j*_dim[0]+i) = res;
                        //                    if (j*_dim[0]+i==_val->size()-1) {
                        //                        _evaluated = true;
                        //                    }
                    }
                    else {
                        if(i*_dim[1]+j>=_val->size()){
                            throw invalid_argument("out of range assignment in eval");
                        }
                        _val->at(i*_dim[1]+j) = res;
                        //                    if (i*_dim[1]+j==_val->size()-1) {
                        //                        _evaluated = true;
                        //                    }
                    }
                }
            }
            return res;
        }
        
        //        type eval(size_t i, size_t j) {
        //
        //
        //            if (is_indexed() && _indices->_ids->size()>1) {
        //                if (_indices->_ids->at(i).at(j) >= _val->size()) {
        //                    throw invalid_argument("eval(i,j): out of range");
        //                }
        //                return _val->at(_indices->_ids->at(i).at(j));
        //            }
        //
        //            if (!is_matrix()) {
        //                return eval(j);
        //            }
        //            if (_is_transposed) {
        //                return _val->at(j*_dim[0]+i);
        //            }
        //            return _val->at(i*_dim[1]+j);
        //        }
        
        template<class T=type, typename enable_if<is_arithmetic<T>::value>::type* = nullptr> bool is_unit() const { /**< Returns true if all values of this paramter are 1 **/
            return (!_is_vector && func_is_number() && _range->first == 1 && _range->second == 1);
            //            return (_range->first == 1 && _range->second == 1);
        }
        
        template<class T=type, typename enable_if<is_same<T, Cpx>::value>::type* = nullptr> bool is_unit() const{
            //            return (func_is_number() && _range->first == Cpx(1,0) && _range->second == Cpx(1,0));
            return (!_is_vector && _range->first == Cpx(1,0) && _range->second == Cpx(1,0));
        }
        
        inline bool is_zero() const { return zero_range();};
        
        template<class T=type, typename enable_if<is_same<T, Cpx>::value>::type* = nullptr> bool zero_range() const{
            //            return (func_is_number() && _range->first == Cpx(0,0) && _range->second == Cpx(0,0));
            return (get_dim()==0 || (_range->first == Cpx(0,0) && _range->second == Cpx(0,0)));
        }
        
        template<typename T=type,
        typename enable_if<is_arithmetic<T>::value>::type* = nullptr> inline bool zero_range() const{
            //            return (get_dim()==0 || (_range->first == 0 && _range->second == 0));
            return (get_dim()==0 || (func_is_number() && _range->first == 0 && _range->second == 0));
        }
        
        
        
        
        bool is_non_positive() const { /**< Returns true if all values of this paramter are <= 0 **/
            auto sgn = get_all_sign();
            return (sgn==non_pos_ || sgn==zero_ || sgn==neg_);
        }
        
        bool is_positive() const { /**< Returns true if all values of this paramter are positive **/
            return (get_all_sign()==pos_);
        }
        
        bool is_non_negative() const { /**< Returns true if all values of this paramter are >= 0 **/
            auto sgn = get_all_sign();
            return (sgn==non_neg_ || sgn==zero_ || sgn==pos_);
        }
        
        bool is_negative() const { /**< Returns true if all values of this paramter are positive **/
            return (get_all_sign()==neg_);
        }
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        func& operator*=(const T2 c){
            return *this *= func<type>(c);
        }
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        func& operator*=(const constant<T2>& c){
            return *this *= func<type>(c);
        }
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        func& operator*=(const param<T2>& p){
            return *this *= func<type>(p);
        }
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        func& operator/=(T2 c){
            return *this /= constant<type>(c);
        }
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        func& operator/=(const constant<T2>& c){
            return *this *= 1./c.eval();
        }
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        func& operator/=(const param<T2>& p){
            auto be = bexpr<type>(div_, make_shared<func>(*this), make_shared<param<T2>>(p));
            auto range = get_div_range(_range,p._range);
            *this = func(be);
            _range = range;
            _evaluated = false;
            _all_convexity = undet_;
            _expr->_range->first = _range->first;
            _expr->_range->second = _range->second;
            _expr->_all_convexity = _all_convexity;
            _expr->_all_sign = _all_sign;
            return *this;
        }
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        func& operator/=(const func<T2>& f){
            if(!is_constant() && f.is_constant()){
                if(f.func_is_number()){
                    return *this *= 1./f.get_val();
                }
                return *this *= 1./f;
            }
            auto be = bexpr<type>(div_, make_shared<func>(*this), make_shared<func>(f));
            auto range = get_div_range(_range,f._range);
            *this = func(be);
            _range = range;
            _evaluated = false;
            _all_convexity = undet_;
            _expr->_range->first = _range->first;
            _expr->_range->second = _range->second;
            _expr->_all_convexity = _all_convexity;
            _expr->_all_sign = _all_sign;
            return *this;
        }
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        func& operator*=(const func<T2>& f){
            if (is_zero()) {
                return *this;
            }
            if (f.is_zero()) {
                reset();
                return *this;
            }
            if (is_unit()) {
                *this = func(f);
                return *this;
            }
            if (f.is_unit()) {
                return *this;
            }
            
            /* Case where c is a number */
            //            if (c.is_number()){
            //                return *this *= constant<T2>(c.eval());
            //            }
            /* Case where the current function is not constant and the other operand is */
            if((!is_constant() && f.is_constant()) || f.func_is_number()) {
                bool transp = false;
                func fc = f;
                if(is_linear() && _is_transposed && f._is_vector){// Situation where (*this)^T * f is transformed into (f^T*(*this))^T
                    fc.transpose();
                    this->transpose();
                    transp = true;
                }
                if (!_cst->is_zero()) {
                    _cst = multiply(_cst,fc);
                }
                for (auto &pair:*_lterms) {
                    pair.second._coef = multiply(pair.second._coef, fc);
                }
                for (auto &pair:*_qterms) {
                    pair.second._coef = multiply(pair.second._coef, fc);
                }
                for (auto &pair:*_pterms) {
                    pair.second._coef = multiply(pair.second._coef, fc);
                }
                if (_expr) {
                    if(_expr->is_uexpr()){
                        if(fc.func_is_number()){
                            _expr->_coef *= fc.eval();
                        }
                        else {
                            _expr = make_shared<bexpr<type>>(bexpr<type>(product_, make_shared<func>(*static_pointer_cast<uexpr<type>>(_expr)), make_shared<func>(fc)));
                        }
                    }
                    else {
                        if(fc.func_is_number()){
                            _expr->_coef *= fc.eval();
                        }
                        else {
                            _expr = make_shared<bexpr<type>>(bexpr<type>(product_, make_shared<func>(*static_pointer_cast<bexpr<type>>(_expr)), make_shared<func>(fc)));
                        }
                    }
                    embed(_expr);
                }
                if (fc.get_all_sign()==unknown_ && !_qterms->empty()) {
                    _all_convexity = undet_;
                }
                
                update_sign_multiply(fc);
                if(f.is_non_positive()){
                    reverse_convexity();
                }
                _evaluated = false;
                _range = get_product_range(_range,fc._range);
                if(transp){
                    this->transpose();
                    _range->first = extended_mult(_range->first,(type)_dim[0]);
                    _range->second = extended_mult(_range->second,(type)_dim[0]);
                }
                update_dot_dim(fc);
                return *this;
            }
            /* Case where the current function is constant and the other operand is not. */
            if (func_is_number() || (is_constant() && !f.is_constant())) {
                auto cpy = this->copy();
                func res = f;
                res.update_dot_dim(*this,f);
                update_sign_multiply(f);
                res._all_sign = _all_sign;
                res._range = get_product_range(_range,f._range);
                if(_is_transposed){
                    res._range->first = extended_mult(res._range->first,(type)_dim[0]);
                    res._range->second = extended_mult(res._range->second,(type)_dim[0]);
                }
                if(is_non_positive()){
                    res.reverse_convexity();
                }
                if (!res._cst->is_zero()) {
                    if (res._cst->is_function()) {
                        auto f_cst = *static_pointer_cast<func<type>>(res._cst);
                        res._cst = multiply(cpy,f_cst);
                    }
                    else if(res._cst->is_param()) {
                        auto f_cst = *static_pointer_cast<param<type>>(res._cst);
                        res._cst = multiply(cpy,f_cst);
                    }
                    else if(res._cst->is_number()) {
                        auto f_cst = *static_pointer_cast<constant<type>>(res._cst);
                        res._cst = multiply(cpy,f_cst);
                    }
                    
                }
                for (auto &pair:*res._lterms) {
                    if (pair.second._coef->is_function()) {
                        auto f_cst = *static_pointer_cast<func<type>>(pair.second._coef);
                        pair.second._coef = multiply(cpy,f_cst);
                    }
                    else if(pair.second._coef->is_param()) {
                        auto f_cst = *static_pointer_cast<param<type>>(pair.second._coef);
                        pair.second._coef = multiply(cpy,f_cst);
                    }
                    else if(pair.second._coef->is_number()) {
                        auto f_cst = *static_pointer_cast<constant<type>>(pair.second._coef);
                        pair.second._coef = multiply(cpy,f_cst);
                    }
                }
                for (auto &pair:*res._qterms) {
                    if (pair.second._coef->is_function()) {
                        auto f_cst = *static_pointer_cast<func<type>>(pair.second._coef);
                        pair.second._coef = multiply(cpy,f_cst);
                    }
                    else if(pair.second._coef->is_param()) {
                        auto f_cst = *static_pointer_cast<param<type>>(pair.second._coef);
                        pair.second._coef = multiply(cpy,f_cst);
                    }
                    else if(pair.second._coef->is_number()) {
                        auto f_cst = *static_pointer_cast<constant<type>>(pair.second._coef);
                        pair.second._coef = multiply(cpy,f_cst);
                    }
                }
                for (auto &pair:*res._pterms) {
                    if (pair.second._coef->is_function()) {
                        auto f_cst = *static_pointer_cast<func<type>>(pair.second._coef);
                        pair.second._coef = multiply(cpy,f_cst);
                    }
                    else if(pair.second._coef->is_param()) {
                        auto f_cst = *static_pointer_cast<param<type>>(pair.second._coef);
                        pair.second._coef = multiply(cpy,f_cst);
                    }
                    else if(pair.second._coef->is_number()) {
                        auto f_cst = *static_pointer_cast<constant<type>>(pair.second._coef);
                        pair.second._coef = multiply(cpy,f_cst);
                    }
                }
                if (res._expr) {
                    if(res._expr->is_uexpr()){
                        if(func_is_number()){
                            res._expr->_coef *= this->eval();
                        }
                        else {
                            res._expr = make_shared<bexpr<type>>(bexpr<type>(product_, make_shared<func<type>>(*this), make_shared<func>(*static_pointer_cast<uexpr<type>>(res._expr))));
                        }
                    }
                    else {
                        if(func_is_number()){
                            res._expr->_coef *= this->eval();
                        }
                        else {
                            res._expr = make_shared<bexpr<type>>(bexpr<type>(product_, make_shared<func<type>>(*this), make_shared<func>(*static_pointer_cast<bexpr<type>>(res._expr))));
                        }
                    }
                    res.embed(res._expr);
                }
                *this = res;
                _evaluated = false;
                return *this;
            }
            /* Both functions are either constant or non-constant at this stage */
            if (_expr || (f._expr)) {
                auto be = bexpr<type>(product_, make_shared<func>(*this), make_shared<func>(f));
                auto res = func(be);
                update_dot_dim(f);
                update_sign_multiply(f);
                res._dim[0] = _dim[0];
                res._dim[1] = _dim[1];
                res._is_vector = _is_vector;
                res._is_transposed = _is_transposed;
                res._all_sign = _all_sign;
                res._range = get_product_range(_range,f._range);
                if(_is_transposed){
                    res._range->first = extended_mult(res._range->first,(type)_dim[0]);
                    res._range->second = extended_mult(res._range->second,(type)_dim[0]);
                }
                *this = move(res);
                _evaluated = false;
                _all_convexity = undet_;
                return *this;
            }
            func res;
            /* Case where the multiplication invlolves multiplying variables/parameters together, i.e., they are both parametric or both include variables  */
            //            if(_to_str==f._to_str){
            if(false){
                auto signp = get_all_sign();
                if(signp==neg_ || signp==pos_){
                    res._all_sign = pos_;
                }
                else {
                    res._all_sign = non_neg_;
                }
                res._range->first=zero<type>().eval();
                if(is_positive()|| is_negative()){
                    res._range->first=extended_mult(_range->first,_range->first);
                }
                res._range->second=extended_mult(_range->second,_range->second);
            }
            else {
                res._range = get_product_range(_range,f._range);
            }
            if(_is_transposed){
                res._range->first = extended_mult(res._range->first,(type)_dim[0]);
                res._range->second = extended_mult(res._range->second,(type)_dim[0]);
            }
            for (auto& t1: *_pterms) {
                if (t1.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial function), we cannot factor the coefficients. Just create a binary expression and return it.
                    auto be = bexpr<type>(product_, make_shared<func>(*this), make_shared<func>(f));
                    *this = func(be);
                    _evaluated = false;
                    _range = res._range;
                    _expr->_range->first = _range->first;
                    _expr->_range->second = _range->second;
                    _expr->_all_convexity = _all_convexity;
                    _expr->_all_sign = _all_sign;
                    return *this;
                }
                for (auto& t2: *f._pterms) {
                    if (t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial function), see comment above.
                        auto be = bexpr<type>(product_, make_shared<func>(*this), make_shared<func>(f));
                        *this = func(be);
                        _evaluated = false;
                        _range = res._range;
                        _expr->_range->first = _range->first;
                        _expr->_range->second = _range->second;
                        _expr->_all_convexity = _all_convexity;
                        _expr->_all_sign = _all_sign;
                        return *this;
                    }
                    auto newl(*t1.second._l);
                    for (auto& it: *t2.second._l) {// TODO check if same l
                        newl.push_back(make_pair<>(it.first, it.second));
                    }
                    if (t2.second._coef->is_function()) {
                        auto f_cst = *static_pointer_cast<func<T2>>(t2.second._coef);
                        auto coef = multiply(t1.second._coef, f_cst);
                        res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
                    }
                    else if(t2.second._coef->is_param()) {
                        auto p_cst = *static_pointer_cast<param<T2>>(t2.second._coef);
                        auto coef = multiply(t1.second._coef, p_cst);
                        res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
                    }
                    else if(t2.second._coef->is_number()) {
                        auto p_cst = *static_pointer_cast<constant<T2>>(t2.second._coef);
                        auto coef = multiply(t1.second._coef, p_cst);
                        res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
                    }
                }
                for (auto& t2: *f._qterms) {
                    if (t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial function), see comment above.
                        auto be = bexpr<type>(product_, make_shared<func>(*this), make_shared<func>(f));
                        *this = func(be);
                        _evaluated = false;
                        _range = res._range;
                        _expr->_range->first = _range->first;
                        _expr->_range->second = _range->second;
                        _expr->_all_convexity = _all_convexity;
                        _expr->_all_sign = _all_sign;
                        return *this;
                    }
                    auto newl(*t1.second._l);
                    newl.push_back(make_pair<>((t2.second._p->first), 1));
                    newl.push_back(make_pair<>((t2.second._p->second), 1));
                    if (t2.second._coef->is_function()) {
                        auto f_cst = *static_pointer_cast<func<T2>>(t2.second._coef);
                        auto coef = multiply(t1.second._coef, f_cst);
                        res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
                    }
                    else if(t2.second._coef->is_param()) {
                        auto p_cst = *static_pointer_cast<param<T2>>(t2.second._coef);
                        auto coef = multiply(t1.second._coef, p_cst);
                        res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
                    }
                    else if(t2.second._coef->is_number()) {
                        auto p_cst = *static_pointer_cast<constant<T2>>(t2.second._coef);
                        auto coef = multiply(t1.second._coef, p_cst);
                        res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
                    }
                }
                for (auto& t2: *f._lterms) {
                    if (t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial function)
                        auto be = bexpr<type>(product_, make_shared<func>(*this), make_shared<func>(f));
                        *this = func(be);
                        _evaluated = false;
                        _range = res._range;
                        _expr->_range->first = _range->first;
                        _expr->_range->second = _range->second;
                        _expr->_all_convexity = _all_convexity;
                        _expr->_all_sign = _all_sign;
                        return *this;
                    }
                    auto newl(*t1.second._l);
                    newl.push_back(make_pair<>((t2.second._p), 1));
                    if (t2.second._coef->is_function()) {
                        auto f_cst = *static_pointer_cast<func<T2>>(t2.second._coef);
                        auto coef = multiply(t1.second._coef, f_cst);
                        res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
                    }
                    else if(t2.second._coef->is_param()) {
                        auto p_cst = *static_pointer_cast<param<T2>>(t2.second._coef);
                        auto coef = multiply(t1.second._coef, p_cst);
                        res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
                    }
                    else if(t2.second._coef->is_number()) {
                        auto p_cst = *static_pointer_cast<constant<T2>>(t2.second._coef);
                        auto coef = multiply(t1.second._coef, p_cst);
                        res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
                    }
                }
                if (!f._cst->is_zero()) {
                    auto newl(*t1.second._l);
                    if (f._cst->is_function()) {
                        auto f_cst = *static_pointer_cast<func<T2>>(f._cst);
                        auto coef = multiply(t1.second._coef, f_cst);
                        res.insert(t1.second._sign, *coef, newl);
                    }
                    else if(f._cst->is_param()) {
                        auto p_cst = *static_pointer_cast<param<T2>>(f._cst);
                        auto coef = multiply(t1.second._coef, p_cst);
                        res.insert(t1.second._sign, *coef, newl);
                    }
                    else if(f._cst->is_number()) {
                        auto p_cst = *static_pointer_cast<constant<T2>>(f._cst);
                        auto coef = multiply(t1.second._coef, p_cst);
                        res.insert(t1.second._sign, *coef, newl);
                    }
                }
            }
            for (auto& t1: *_qterms) {
                if (t1.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(Quadratic term)
                    auto be = bexpr<type>(product_, make_shared<func>(*this), make_shared<func>(f));
                    *this = func(be);
                    _evaluated = false;
                    _range = res._range;
                    return *this;
                }
                for (auto& t2: *f._pterms) {
                    if (t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial term)
                        auto range = get_product_range(_range,f._range);
                        auto be = bexpr<type>(product_, make_shared<func>(*this), make_shared<func>(f));
                        *this = func(be);
                        _evaluated = false;
                    }
                    auto newl(*t2.second._l);
                    newl.push_front(make_pair<>(t1.second._p->first, 1));
                    newl.push_front(make_pair<>(t1.second._p->second, 1));
                    if (t2.second._coef->is_function()) {
                        auto f_cst = *static_pointer_cast<func<T2>>(t2.second._coef);
                        auto coef = multiply(t1.second._coef, f_cst);
                        res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
                    }
                    else if(t2.second._coef->is_param()) {
                        auto p_cst = *static_pointer_cast<param<T2>>(t2.second._coef);
                        auto coef = multiply(t1.second._coef, p_cst);
                        res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
                    }
                    else if(t2.second._coef->is_number()) {
                        auto p_cst = *static_pointer_cast<constant<T2>>(t2.second._coef);
                        auto coef = multiply(t1.second._coef, p_cst);
                        res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
                    }
                }
                for (auto& t2: *f._qterms) {
                    if (t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial term)
                        auto range = get_product_range(_range,f._range);
                        auto be = bexpr<type>(product_, make_shared<func>(*this), make_shared<func>(f));
                        *this = func(be);
                        _evaluated = false;
                    }
                    list<pair<shared_ptr<param_>, int>> newl;
                    newl.push_back(make_pair<>(t1.second._p->first, 1));
                    newl.push_back(make_pair<>(t1.second._p->second, 1));
                    newl.push_back(make_pair<>(t2.second._p->first, 1));
                    newl.push_back(make_pair<>(t2.second._p->second, 1));
                    if (t2.second._coef->is_function()) {
                        auto f_cst = *static_pointer_cast<func<T2>>(t2.second._coef);
                        auto coef = multiply(t1.second._coef, f_cst);
                        res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
                    }
                    else if(t2.second._coef->is_param()) {
                        auto p_cst = *static_pointer_cast<param<T2>>(t2.second._coef);
                        auto coef = multiply(t1.second._coef, p_cst);
                        res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
                    }
                    else if(t2.second._coef->is_number()) {
                        auto p_cst = *static_pointer_cast<constant<T2>>(t2.second._coef);
                        auto coef = multiply(t1.second._coef, p_cst);
                        res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
                    }
                }
                for (auto& t2: *f._lterms) {
                    if (t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial term)
                        auto range = get_product_range(_range,f._range);
                        auto be = bexpr<type>(product_, make_shared<func>(*this), make_shared<func>(f));
                        *this = func(be);
                        _evaluated = false;
                    }
                    list<pair<shared_ptr<param_>, int>> newl;
                    newl.push_back(make_pair<>(t1.second._p->first, 1));
                    newl.push_back(make_pair<>(t1.second._p->second, 1));
                    newl.push_back(make_pair<>(t2.second._p, 1));
                    if (t2.second._coef->is_function()) {
                        auto f_cst = *static_pointer_cast<func<T2>>(t2.second._coef);
                        auto coef = multiply(t1.second._coef, f_cst);
                        res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
                    }
                    else if(t2.second._coef->is_param()) {
                        auto p_cst = *static_pointer_cast<param<T2>>(t2.second._coef);
                        auto coef = multiply(t1.second._coef, p_cst);
                        res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
                    }
                    else if(t2.second._coef->is_number()) {
                        auto p_cst = *static_pointer_cast<constant<T2>>(t2.second._coef);
                        auto coef = multiply(t1.second._coef, p_cst);
                        res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
                    }                    }
                if (!f._cst->is_zero()) {
                    if (f._cst->is_function()) {
                        auto f_cst = *static_pointer_cast<func<T2>>(f._cst);
                        auto coef = multiply(t1.second._coef, f_cst);
                        res.insert(t1.second._sign, *coef, *t1.second._p->first, *t1.second._p->second);
                    }
                    else if(f._cst->is_param()) {
                        auto p_cst = *static_pointer_cast<param<T2>>(f._cst);
                        auto coef = multiply(t1.second._coef, p_cst);
                        res.insert(t1.second._sign, *coef, *t1.second._p->first, *t1.second._p->second);
                    }
                    else if(f._cst->is_number()) {
                        auto p_cst = *static_pointer_cast<constant<T2>>(f._cst);
                        auto coef = multiply(t1.second._coef, p_cst);
                        res.insert(t1.second._sign, *coef, *t1.second._p->first, *t1.second._p->second);
                    }
                }
            }
            for (auto& t1: *_lterms) {
                for (auto& t2: *f._pterms) {
                    if (t1.second._coef->_is_transposed || t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial term)
                        auto be = bexpr<type>(product_, make_shared<func>(*this), make_shared<func>(f));
                        *this = func(be);
                        _evaluated = false;
                    }
                    auto newl(*t2.second._l);
                    newl.push_front(make_pair<>((t1.second._p), 1));
                    if (t2.second._coef->is_function()) {
                        auto f_cst = *static_pointer_cast<func<T2>>(t2.second._coef);
                        auto coef = multiply(t1.second._coef, f_cst);
                        res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
                    }
                    else if(t2.second._coef->is_param()) {
                        auto p_cst = *static_pointer_cast<param<T2>>(t2.second._coef);
                        auto coef = multiply(t1.second._coef, p_cst);
                        res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
                    }
                    else if(t2.second._coef->is_number()) {
                        auto p_cst = *static_pointer_cast<constant<T2>>(t2.second._coef);
                        auto coef = multiply(t1.second._coef, p_cst);
                        res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
                    }
                }
                for (auto& t2: *f._qterms) {
                    if (t1.second._coef->_is_transposed || t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial term)
                        auto be = bexpr<type>(product_, make_shared<func>(*this), make_shared<func>(f));
                        *this = func(be);
                        _evaluated = false;
                    }
                    list<pair<shared_ptr<param_>, int>> newl;
                    newl.push_back(make_pair<>(t1.second._p, 1));
                    newl.push_back(make_pair<>(t2.second._p->first, 1));
                    newl.push_back(make_pair<>(t2.second._p->second, 1));
                    if (t2.second._coef->is_function()) {
                        auto f_cst = *static_pointer_cast<func<T2>>(t2.second._coef);
                        auto coef = multiply(t1.second._coef, f_cst);
                        res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
                    }
                    else if(t2.second._coef->is_param()) {
                        auto p_cst = *static_pointer_cast<param<T2>>(t2.second._coef);
                        auto coef = multiply(t1.second._coef, p_cst);
                        res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
                    }
                    else if(t2.second._coef->is_number()) {
                        auto p_cst = *static_pointer_cast<constant<T2>>(t2.second._coef);
                        auto coef = multiply(t1.second._coef, p_cst);
                        res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
                    }
                }
                for (auto& t2: *f._lterms) {
                    if (t2.second._coef->is_function()) {
                        auto f_cst = *static_pointer_cast<func<T2>>(t2.second._coef);
                        auto coef = multiply(t1.second._coef, f_cst);
                        auto p1 = t1.second._p;
                        if(_is_transposed){
                            coef->transpose();
                            p1->transpose();
                        }
                        res.insert(!(t1.second._sign^t2.second._sign), *coef, *p1, *t2.second._p, _is_transposed);
                    }
                    else if(t2.second._coef->is_param()) {
                        auto p_cst = *static_pointer_cast<param<T2>>(t2.second._coef);
                        auto coef = multiply(t1.second._coef, p_cst);
                        auto p1 = t1.second._p;
                        if(_is_transposed){
                            coef->transpose();
                            p1->transpose();
                        }
                        res.insert(!(t1.second._sign^t2.second._sign), *coef, *p1, *t2.second._p, _is_transposed);
                    }
                    else if(t2.second._coef->is_number()) {
                        auto p_cst = *static_pointer_cast<constant<T2>>(t2.second._coef);
                        auto coef = multiply(t1.second._coef, p_cst);
                        auto p1 = t1.second._p;
                        if(_is_transposed){
                            coef->transpose();
                            p1->transpose();
                        }
                        res.insert(!(t1.second._sign^t2.second._sign), *coef, *p1, *t2.second._p, _is_transposed);
                    }
                }
                if (!f._cst->is_zero()) {
                    if (f._cst->is_function()) {
                        auto f_cst = *static_pointer_cast<func<T2>>(f._cst);
                        auto coef = multiply(t1.second._coef, f_cst);
                        res.insert(t1.second._sign, *coef, *t1.second._p);
                    }
                    else if(f._cst->is_param()) {
                        auto p_cst = *static_pointer_cast<param<T2>>(f._cst);
                        auto coef = multiply(t1.second._coef, p_cst);
                        res.insert(t1.second._sign, *coef, *t1.second._p);
                    }
                    else if(f._cst->is_number()) {
                        auto p_cst = *static_pointer_cast<constant<T2>>(f._cst);
                        auto coef = multiply(t1.second._coef, p_cst);
                        res.insert(t1.second._sign, *coef, *t1.second._p);
                    }
                }
            }
            if (!_cst->is_zero()) {
                for (auto& t2: *f._pterms) {
                    if (t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial term)
                        auto be = bexpr<type>(product_, make_shared<func>(*this), make_shared<func>(f));
                        *this = func(be);
                        _evaluated = false;
                        return *this;
                    }
                    if (t2.second._coef->is_function()) {
                        auto f_cst = *static_pointer_cast<func<T2>>(t2.second._coef);
                        auto coef = multiply(_cst, f_cst);
                        res.insert(t2.second._sign, *coef, *t2.second._l);
                    }
                    else if(t2.second._coef->is_param()) {
                        auto p_cst = *static_pointer_cast<param<T2>>(t2.second._coef);
                        auto coef = multiply(_cst, p_cst);
                        res.insert(t2.second._sign, *coef, *t2.second._l);
                    }
                    else if(t2.second._coef->is_number()) {
                        auto p_cst = *static_pointer_cast<constant<T2>>(t2.second._coef);
                        auto coef = multiply(_cst, p_cst);
                        res.insert(t2.second._sign, *coef, *t2.second._l);
                    }
                }
                for (auto& t2: *f._qterms) {
                    if (t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial term)
                        auto be = bexpr<type>(product_, make_shared<func>(*this), make_shared<func>(f));
                        *this = func(be);
                        _evaluated = false;
                        return *this;
                    }
                    if (t2.second._coef->is_function()) {
                        auto f_cst = *static_pointer_cast<func<T2>>(t2.second._coef);
                        auto coef = multiply(_cst, f_cst);
                        res.insert(t2.second._sign, *coef, *t2.second._p->first, *t2.second._p->second);
                    }
                    else if(t2.second._coef->is_param()) {
                        auto p_cst = *static_pointer_cast<param<T2>>(t2.second._coef);
                        auto coef = multiply(_cst, p_cst);
                        res.insert(t2.second._sign, *coef, *t2.second._p->first, *t2.second._p->second);
                    }
                    else if(t2.second._coef->is_number()) {
                        auto p_cst = *static_pointer_cast<constant<T2>>(t2.second._coef);
                        auto coef = multiply(_cst, p_cst);
                        res.insert(t2.second._sign, *coef, *t2.second._p->first, *t2.second._p->second);
                    }
                }
                for (auto& t2: *f._lterms) {
                    if (t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial term)
                        auto be = bexpr<type>(product_, make_shared<func>(*this), make_shared<func>(f));
                        *this = func(be);
                        _evaluated = false;
                        return *this;
                    }
                    if (t2.second._coef->is_function()) {
                        auto f_cst = *static_pointer_cast<func<T2>>(t2.second._coef);
                        auto coef = multiply(_cst, f_cst);
                        res.insert(t2.second._sign, *coef, *t2.second._p);
                    }
                    else if(t2.second._coef->is_param()) {
                        auto p_cst = *static_pointer_cast<param<T2>>(t2.second._coef);
                        auto coef = multiply(_cst, p_cst);
                        res.insert(t2.second._sign, *coef, *t2.second._p);
                    }
                    else if(t2.second._coef->is_number()) {
                        auto p_cst = *static_pointer_cast<constant<T2>>(t2.second._coef);
                        auto coef = multiply(_cst, p_cst);
                        res.insert(t2.second._sign, *coef, *t2.second._p);
                    }
                }
                if (!f._cst->is_zero()) {
                    if (f._cst->is_function()) {
                        auto f_cst = *static_pointer_cast<func<T2>>(f._cst);
                        res._cst = multiply(_cst, f_cst);
                    }
                    else if(f._cst->is_param()) {
                        auto p_cst = *static_pointer_cast<param<T2>>(f._cst);
                        res._cst = multiply(_cst, p_cst);
                    }
                    else if(f._cst->is_number()) {
                        auto p_cst = *static_pointer_cast<constant<T2>>(f._cst);
                        res._cst = multiply(_cst, p_cst);
                    }
                }
            }
            res.update_dot_dim(*this, f);
            if(res.is_quadratic()){res.update_quad_convexity();}
            else {_all_convexity = undet_;}
            res._all_sign = sign_product(_all_sign, f.get_all_sign());
            *this = move(res);
            _evaluated = false;
            return *this;
        }
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        func& operator+=(const T2 c){
            return *this += func<type>(c);
        }
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        func& operator+=(const constant<T2>& c){
            return *this += func<type>(c);
        }
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        func& operator+=(const param<T2>& p){
            return *this += func<type>(p);
        }
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        func& operator+=(const func<T2>& f){
            if (f.is_zero()) {
                return *this;
            }
            if(is_zero()){
                return *this = f;
            }
            _evaluated = false;
            set_max_dim(f);
            if (is_constant() && !f.is_constant()) {
                func res(f);
                res += *this;
                return *this = res;
            }
            if (!is_constant() && f.is_constant()) {
                this->add_cst(f);
                update_sign_add(f);
                _range = get_plus_range(_range, f._range);
                return *this;
            }
            if (!f.get_cst()->is_zero()) {
                if (f.get_cst()->is_number()) {
                    auto f_cst = static_pointer_cast<constant<T2>>(f.get_cst());
                    add_cst(*f_cst);
                }
                else if (f.get_cst()->is_param()) {
                    auto f_cst = static_pointer_cast<param<T2>>(f.get_cst());
                    add_cst(*f_cst);
                }
                else {
                    auto f_cst = static_pointer_cast<func<T2>>(f.get_cst());
                    add_cst(*f_cst);
                }
                if (_cst->is_function()) {
                    embed(*static_pointer_cast<func>(_cst));
                }
            }
            for (auto &pair:*f._lterms) {
                auto term = pair.second;
                if (term._coef->is_function()) {
                    auto coef = *static_pointer_cast<func<T2>>(term._coef);
                    term._coef = func(coef).copy();
                }
                else if(term._coef->is_param()) {
                    auto coef = *static_pointer_cast<param<T2>>(term._coef);
                    term._coef = param<type>(coef).copy();
                }
                else if(term._coef->is_number()) {
                    auto coef = *static_pointer_cast<constant<T2>>(term._coef);
                    term._coef = constant<type>(coef).copy();//TODO if T2==type no need to cast
                }
                this->insert(term);
            }
            for (auto &pair:*f._qterms) {
                auto term = pair.second;
                if (term._coef->is_function()) {
                    auto coef = *static_pointer_cast<func<T2>>(term._coef);
                    term._coef = func(coef).copy();
                }
                else if(term._coef->is_param()) {
                    auto coef = *static_pointer_cast<param<T2>>(term._coef);
                    term._coef = param<type>(coef).copy();
                }
                else if(term._coef->is_number()) {
                    auto coef = *static_pointer_cast<constant<T2>>(term._coef);
                    term._coef = constant<type>(coef).copy();
                }
                this->insert(term);
            }
            for (auto &pair:*f._pterms) {
                auto term = pair.second;
                if (term._coef->is_function()) {
                    auto coef = *static_pointer_cast<func<T2>>(term._coef);
                    term._coef = func(coef).copy();
                }
                else if(term._coef->is_param()) {
                    auto coef = *static_pointer_cast<param<T2>>(term._coef);
                    term._coef = param<type>(coef).copy();
                }
                else if(term._coef->is_number()) {
                    auto coef = *static_pointer_cast<constant<T2>>(term._coef);
                    term._coef = constant<type>(coef).copy();
                }
                this->insert(term);
            }
            if (_expr && f.get_expr()) {
                func f1,f2;
                if (_expr->is_uexpr()) {
                    f1 = func(*static_pointer_cast<uexpr<type>>(_expr));
                }
                else {
                    f1 = func(*static_pointer_cast<bexpr<type>>(_expr));
                }
                if (f.get_expr()->is_uexpr()) {
                    auto ue = *static_pointer_cast<uexpr<T2>>(f.get_expr());
                    if (ue._son->is_function()) {
                        auto son = static_pointer_cast<func<T2>>(ue._son);
                        ue._son = make_shared<func>(*son);
                    }
                    f2 = func(ue);
                }
                else {
                    auto bexp = *static_pointer_cast<bexpr<T2>>(f.get_expr());
                    if (bexp._lson->is_function()) {
                        auto lson = static_pointer_cast<func<T2>>(bexp._lson);
                        bexp._lson = make_shared<func>(*lson);
                    }
                    if (bexp._rson->is_function()) {
                        auto rson = static_pointer_cast<func<T2>>(bexp._rson);
                        bexp._rson = make_shared<func>(*rson);
                    }
                    f2 = func(bexp);
                }
                _expr = make_shared<bexpr<type>>(bexpr<type>(plus_, f1.copy(), f2.copy()));
                embed(_expr);
            }
            else if (!_expr && f.get_expr()) {
                if (f.get_expr()->is_uexpr()) {
                    auto ue = *static_pointer_cast<uexpr<T2>>(f.get_expr());
                    if (ue._son->is_function()) {
                        auto son = static_pointer_cast<func<T2>>(ue._son);
                        ue._son = make_shared<func>(*son);
                    }
                    _expr = make_shared<uexpr<type>>(ue);
                }
                else {
                    auto bexp = make_shared<bexpr<type>>(*static_pointer_cast<bexpr<T2>>(f.get_expr()));
                    if (bexp->_lson->is_function()) {
                        auto son = static_pointer_cast<func<T2>>(bexp->_lson);
                        bexp->_lson = make_shared<func>(*son);
                    }
                    if (bexp->_rson->is_function()) {
                        auto son = static_pointer_cast<func<T2>>(bexp->_rson);
                        bexp->_rson = make_shared<func>(*son);
                    }
                    _expr = bexp;
                }
                embed(_expr);
                if (!_vars->empty()) {
                    _ftype = nlin_;
                }
            }
            update_sign_add(f);
            if(is_quadratic()){
                update_quad_convexity();
            }
            else {
                update_convexity_add(f._all_convexity);
            }
            _range = get_plus_range(_range, f._range);
            if(func_is_number()){
                _ftype = const_;
                set_range(eval(_cst));
            }
            _evaluated = false;
            return *this;
        }
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        func& operator-=(const T2 c){
            return *this -= func<type>(c);
        }
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        func& operator-=(const constant<T2>& c){
            return *this -= func<type>(c);
        }
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        func& operator-=(const param<T2>& p){
            return *this -= func<type>(p);
        }
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        func& operator-=(const func<T2>& f){
            auto res = f;
            res.reverse_sign();
            return *this += res;
        }
        //
        //        func_ tr() const {
        //            auto f = func_(*this);
        //            f.transpose();
        //            return f;
        //        }
        //
        /** Reset all fields to default values */
        void reset(){
            _to_str = "";
            update_range();
            _all_range = nullptr;
            _vars->clear();
            _val->clear();
            _params->clear();
            if(_dfdx){
                _dfdx->clear();
            }
            if(_hess_link){
                _hess_link->clear();
            }
            _convexity = nullptr;
            _sign = nullptr;
            _expr = nullptr;
            _ftype = const_;
            _all_convexity = linear_;
            _all_sign = zero_;
            _is_transposed = false;
            _is_vector = false;
            _evaluated = true;
            _embedded = false;
            _dim[0] = 1;
            _dim[1] = 1;
            this->_val->clear();//TODO all_range?
            _lterms->clear();
            _qterms->clear();
            _pterms->clear();
            _cst = make_shared<constant<type>>();
            _nb_vars = 0;
            _nnz_h = 0;
            _nnz_j = 0;
        };
        
        func tr() const {
            auto f = *this;
            f.transpose();
            return f;
        }
        
        void update_quad_convexity(){
            if(is_unitary()){
                //TODO check second derivative
            }
            if (!_pterms->empty()) {
                _all_convexity = undet_;
                return;
            }
            if (_qterms->empty() && !_expr) {
                _all_convexity = linear_;
                return;
            }
            if (!_qterms->empty() && !_expr) {
                _all_convexity = get_convexity(_qterms->begin()->second);
                for (auto pair_it = next(_qterms->begin()); pair_it != _qterms->end(); pair_it++) {
                    Convexity conv = get_convexity(pair_it->second);
                    if (_all_convexity==undet_ || conv ==undet_ || (_all_convexity==convex_ && conv==concave_) || (_all_convexity==concave_ && conv==convex_)) {
                        _all_convexity = undet_;
                        return;
                    }
                    else {
                        _all_convexity = conv;
                    }
                }
            }
        }
        
        /**
         Returns the convexity of current function if quadratic term q was to be added.
         @param[in] q quadratic term to be added.
         @return convexity of function if q = coef*p1*p2 was to be added.
         */
        Convexity get_convexity(const qterm& q){
            if(q._p->first == q._p->second){
                if (q._sign && (q._coef->is_positive() || q._coef->is_non_negative())) {
                    return convex_;
                }
                if (q._sign && (q._coef->is_negative() || q._coef->is_non_positive())) {
                    return concave_;
                }
                if (!q._sign && (q._coef->is_negative() || q._coef->is_non_positive())) {
                    return convex_;
                }
                if (!q._sign && (q._coef->is_negative() || q._coef->is_non_positive())) {
                    return concave_;
                }
            }
            // At this stage, we know that q._p->first !=q._p->second
            // Checking if the product can be factorized
            auto sqr1 = get_square(q._p->first);
            auto sqr2 = get_square(q._p->second);
            if (sqr1 && sqr2){
                auto c1 = sqr1->_coef;
                auto c2 = sqr2->_coef;
                if ((sqr1->_sign^c1->is_positive())==(sqr2->_sign^c2->is_positive())) {
                    if (c1->func_is_number() && c2->func_is_number() && q._coef->func_is_number()) {
                        if (2.*std::sqrt(eval<type>(c1)*eval<type>(c2)) >= eval<type>(q._coef)) {
                            if (!(sqr1->_sign^c1->is_positive())) {
                                return convex_;
                            }
                            return concave_;
                        }
                    }
                    return undet_;
                }
            }
            return undet_;
        }
        
        bool check_rotated_soc(){
            if (_qterms->empty() || !_pterms->empty() || _expr) {
                return false;
            }
            unsigned nb_bilinear = 0, nb_quad = 0;
            Sign bilinear_sign = unknown_, quadratic_sign = unknown_, var1_sign = unknown_, var2_sign = unknown_;
            for (auto &qt_pair: *_qterms) {
                if (qt_pair.second._p->first!=qt_pair.second._p->second) {
                    bilinear_sign = qt_pair.second.get_all_sign();
                    var1_sign = qt_pair.second._p->first->get_all_sign();
                    var2_sign = qt_pair.second._p->second->get_all_sign();
                    if (bilinear_sign==unknown_ || var1_sign==neg_ || var2_sign==neg_) {
                        return false;
                    }
                    nb_bilinear++;
                    if (nb_bilinear > 1) {
                        return false;
                    }
                }
                else{
                    nb_quad++;
                    auto sign = qt_pair.second.get_all_sign();
                    if (quadratic_sign!=unknown_ && quadratic_sign!=sign) {
                        return false;
                    }
                    if (quadratic_sign!=unknown_ && quadratic_sign==bilinear_sign) {
                        return false;
                    }
                    else {
                        quadratic_sign = sign;
                    }
                }
            }
            if(nb_quad==0){
                return false;
            }
            if (bilinear_sign==pos_) {
                _all_convexity = concave_;
                return true;
            }
            else if(bilinear_sign==neg_) {
                _all_convexity = convex_;
                return true;
            }
            return false;
        }
        
        
        bool check_soc(){
            if (_qterms->empty() || !_pterms->empty() || _expr) {
                return false;
            }
            unsigned nb_neg = 0, nb_pos = 0;
            for (auto &qt_pair: *_qterms) {
                if (qt_pair.second._p->first!=qt_pair.second._p->second) {
                    return false;
                }
                auto sign = qt_pair.second.get_all_sign();
                if (sign==unknown_) {
                    return false;
                }
                if (sign==pos_) {
                    nb_pos++;
                }
                else if(sign==neg_){
                    nb_neg++;
                }
            }
            if (nb_neg==1 && nb_pos>1) {
                _all_convexity = convex_;
                return true;
            }
            else if (nb_pos==1 && nb_neg>1){
                _all_convexity = concave_;
                return true;
            }
            return false;
        }
        
        
        
        /**
         Subfuntion of embed(func_&& f). Merge variables and parameters with f. If a variable x in f exists in the current funtion, x will now point to the same variable appearing in current function.
         @param[in] f function to merge variables and parameters with.
         */
        void merge_vars(func& f){
            for (auto &pair:*f._lterms) {
                auto coef = pair.second._coef;
                if(coef->is_function()){
                    embed(*static_pointer_cast<func>(coef));
                }
                auto p = pair.second._p;
                if (p->is_var()) {
                    auto pname = p->get_name(false,false);
                    auto it = _vars->find(pname);
                    if (it==_vars->end()) {
                        add_var(f.get_var(pname));
                    }
                    else{
                        pair.second._p = it->second.first;
                        it->second.second++;
                    }
                }
                else {
                    auto pname = p->get_name(false,false);
                    auto it = _params->find(pname);
                    if (it==_params->end()) {
                        add_param(f.get_param(pname));
                    }
                    else{
                        pair.second._p = it->second.first;
                        it->second.second++;
                    }
                }
            }
            for (auto &pair:*f._qterms) {
                auto coef = pair.second._coef;
                if (coef->is_function()){
                    embed(*static_pointer_cast<func>(coef));
                }
                auto p1 = pair.second._p->first;
                auto p2 = pair.second._p->second;
                if (p1->is_var()) {
                    auto it1 = _vars->find(p1->get_name(false,false));
                    if (it1==_vars->end()) {
                        add_var(f.get_var(p1->get_name(false,false)));
                    }
                    else{
                        pair.second._p->first = it1->second.first;
                        it1->second.second++;
                    }
                    auto it2 = _vars->find(p2->get_name(false,false));
                    if (it2==_vars->end()) {
                        add_var(f.get_var(p2->get_name(false,false)));
                    }
                    else{
                        pair.second._p->second = it2->second.first;
                        it2->second.second++;
                    }
                }
                else {
                    auto it1 = _params->find(p1->get_name(false,false));
                    if (it1==_params->end()) {
                        add_param(f.get_param(p1->get_name(false,false)));
                    }
                    else{
                        pair.second._p->first = it1->second.first;
                        it1->second.second++;
                    }
                    auto it2 = _params->find(p2->get_name(false,false));
                    if (it2==_params->end()) {
                        add_param(f.get_param(p2->get_name(false,false)));
                    }
                    else{
                        pair.second._p->second = it2->second.first;
                        it2->second.second++;
                    }
                }
            }
            for (auto &pair:*f._pterms) {
                auto coef = pair.second._coef;
                if(coef->is_function()){
                    embed(*static_pointer_cast<func>(coef));
                }
                auto list = pair.second._l;
                for (auto &ppi: *list) {
                    auto p = ppi.first;
                    if (p->is_var()) {
                        auto it = _vars->find(p->get_name(false,false));
                        if (it==_vars->end()) {
                            add_var(f.get_var(p->get_name(false,false)));
                        }
                        else{
                            ppi.first = it->second.first;
                            it->second.second++;
                        }
                    }
                    else {
                        auto it = _params->find(p->get_name(false,false));
                        if (it==_params->end()) {
                            add_param(f.get_param(p->get_name(false,false)));
                        }
                        else{
                            ppi.first = it->second.first;
                            it->second.second++;
                        }
                    }
                }
            }
            if (f._expr) {
                embed(f._expr);
            }
            if(f._cst->is_function()){
                embed(*static_pointer_cast<func>(f._cst));
            }
            
            auto old_vars = *f._vars;
            for (auto &vp: old_vars) {
                auto vv = (*_vars)[vp.first].first;
                auto vv_f = (*f._vars)[vp.first].first;
                if (vv != vv_f) {
                    //                delete vv_f;
                    f._vars->erase(vp.first);
                    f._vars->insert(make_pair<>(vp.first, make_pair<>(vv, 1)));
                }
            }
            auto old_params = *f._params;
            for (auto &pp: old_params) {
                auto p = (*_params)[pp.first].first;
                auto p_f = (*f._params)[pp.first].first;
                if (p != p_f) {
                    //                delete p_f;
                    f._params->erase(pp.first);
                    f._params->insert(make_pair<>(pp.first, make_pair<>(p, 1)));
                }
            }
        }
        
        /**
         Merge variables and parameters with expression e. If a variable x in e exists in the current funtion, x will now point to the same variable appearing in the current function.
         @param[in] e expression to merge variables and parameters with.
         */
        void embed(shared_ptr<expr<type>> e){
            _evaluated = false;
            switch (e->get_type()) {
                case uexp_c:{
                    auto ue = static_pointer_cast<uexpr<type>>(e);
                    if (ue->_son->is_function()) {
                        embed(*static_pointer_cast<func>(ue->_son));
                    }
                    else if(ue->_son->is_expr()){
                        embed(static_pointer_cast<expr<type>>(ue->_son));
                    }
                    else if (ue->_son->is_param() || ue->_son->is_var() ){
                        auto p = static_pointer_cast<param_>(ue->_son);
                        auto name = p->get_name(false,false);
                        if (p->is_var()) {
                            auto pnew = get_var(name);
                            if (!pnew) {
                                pnew = p;
                                add_var(pnew,1);
                            }
                            else {
                                ue->_son = pnew;
                            }
                        }
                        else {
                            auto pnew = get_param(name);
                            if (!pnew) {
                                pnew = p;
                                add_param(pnew);
                            }
                            else {
                                ue->_son = pnew;
                            }
                        }
                    }
                    break;
                }
                case bexp_c:{
                    auto be = static_pointer_cast<bexpr<type>>(e);
                    if (be->_lson->is_function()) {
                        embed(*static_pointer_cast<func>(be->_lson));
                    }
                    else if(be->_lson->is_expr()){
                        embed(static_pointer_cast<expr<type>>(be->_lson));
                    }
                    else if (be->_lson->is_param() || be->_lson->is_var() ){
                        auto p = static_pointer_cast<param_>(be->_lson);
                        auto name = p->get_name(false,false);
                        if (p->is_var()) {
                            auto pnew = get_var(name);
                            if (!pnew) {
                                pnew = p;
                                add_var(pnew,1);
                            }
                            else {
                                be->_lson = pnew;
                            }
                        }
                        else {
                            auto pnew = get_param(name);
                            if (!pnew) {
                                pnew = p;
                                add_param(pnew);
                            }
                            else {
                                be->_lson = pnew;
                            }
                        }
                    }
                    if (be->_rson->is_function()) {
                        embed(*static_pointer_cast<func>(be->_rson));
                    }
                    else if(be->_rson->is_expr()){
                        embed(static_pointer_cast<expr<type>>(be->_rson));
                    }
                    else if (be->_rson->is_param() || be->_rson->is_var() ){
                        auto p = static_pointer_cast<param_>(be->_rson);
                        auto name = p->get_name(false,false);
                        if (p->is_var()) {
                            auto pnew = get_var(name);
                            if (!pnew) {
                                pnew = p;
                                add_var(pnew,1);
                            }
                            else {
                                be->_rson = pnew;
                            }
                        }
                        else {
                            auto pnew = get_param(name);
                            if (!pnew) {
                                pnew = p;
                                add_param(pnew);
                            }
                            else {
                                be->_rson = pnew;
                            }
                        }
                    }
                }
                default:
                    break;
            }
        }
        
        /**
         Mark f as embeded and merge variables and parameters with f (by calling merge_vars(func_&& f). If a variable x in f exists in the current funtion, x will now point to the same variable appearing in current function.
         @param[in] f function to merge variables and parameters with.
         */
        void embed(func& f){
            f._embedded = true;
            merge_vars(f);
        }
        
        /**
         Relax and replace integer variables with continuous ones provided in argument vars.
         @param[in] vars set with continuous variables replacements.
         */
        void relax(const map<size_t, shared_ptr<param_>>& vars){
            auto new_vars = make_shared<map<string, pair<shared_ptr<param_>, unsigned>>>();
            bool has_int = false;
            for (auto &v_p:*_vars) {
                auto old_var = v_p.second.first;
                auto nb_occ = v_p.second.second;
                auto new_var = vars.at(old_var->get_vec_id())->pcopy();
                new_var->shallow_copy(*old_var);
                (*new_vars)[new_var->get_name(false,false)] = make_pair<>(new_var,nb_occ);
                if (old_var->is_binary() || old_var->is_short() || old_var->is_integer()) {
                    has_int = true;
                    new_var->_is_relaxed = true;
                }
            }
            if (!has_int) {
                return;
            }
            
            for (auto &lt:get_lterms()) {
                lt.second._p = new_vars->at(lt.second._p->get_name(false,false)).first;
            }
            for (auto &lt:get_qterms()) {
                lt.second._p->first = new_vars->at(lt.second._p->first->get_name(false,false)).first;
                lt.second._p->second = new_vars->at(lt.second._p->second->get_name(false,false)).first;
            }
            for (auto &lt:get_pterms()) {
                for (auto &v_p:*lt.second._l) {
                    v_p.first = new_vars->at(v_p.first->get_name(false,false)).first;
                }
            }
            if (_expr) {
                if (_expr->is_uexpr()) {
                    auto ue = static_pointer_cast<uexpr<type>>(_expr);
                    ue->_son->relax(vars);
                }
                else {
                    auto be = static_pointer_cast<bexpr<type>>(_expr);
                    be->_lson->relax(vars);
                    be->_rson->relax(vars);
                }
            }
            _vars = new_vars;
        }
        
        void update_vars(){merge_vars(*this);};
        
        bool insert(const constant_& coef, const param_& p1, const param_& p2, bool coef_p1_tr=false){
            return insert(true, coef, p1, p2, coef_p1_tr);
        };
        
        bool insert(bool sign, const constant_& coef, const param_& p1, const param_& p2, bool c_p1_transposed=false){/**< Adds coef*p1*p2 to the function. Returns true if added new term, false if only updated coef of p1*p2 */
            auto ps1 = p1.get_name(false,false);
            auto ps2 = p2.get_name(false,false);
            auto qname = ps1+","+ps2;
            auto pair_it = _qterms->find(qname);
            if (pair_it == _qterms->end()) {
                qname = ps2+","+ps1;
                pair_it = _qterms->find(qname);
            }
            if (pair_it == _qterms->end()) {
                qname = ps1+","+ps2;
            }
            shared_ptr<param_> p_new1;
            shared_ptr<param_> p_new2;
            _evaluated=false;
            if (_ftype <= lin_ && p1.is_var()) {
                _ftype = quad_;
            }
            if (pair_it == _qterms->end()) {
                if (p1.is_var()) {
                    p_new1 = get_var(ps1);
                    if (!p_new1) {
                        p_new1 = p1.pcopy();
                        add_var(p_new1);
                    }
                    else {
                        incr_occ_var(ps1);
                    }
                }
                else {
                    p_new1 = get_param(ps1);
                    if (!p_new1) {
                        p_new1 = p1.pcopy();
                        add_param(p_new1);
                    }
                    else {
                        incr_occ_param(ps1);
                    }
                    
                }
                if (p2.is_var()) {
                    p_new2 = get_var(ps2);
                    if (!p_new2) {
                        p_new2 = p2.pcopy();
                        add_var(p_new2);
                    }
                    else {
                        incr_occ_var(ps2);
                    }
                }
                else {
                    p_new2 = get_param(ps2);
                    if (!p_new2) {
                        p_new2 = p2.pcopy();
                        add_param(p_new2);
                    }
                    else {
                        incr_occ_param(ps2);
                    }
                }
                auto c_new = coef.copy();
                if (c_new->is_function()) {
                    embed(*static_pointer_cast<func>(c_new));
                }
                qterm q(sign, c_new, p_new1, p_new2);
                q._coef_p1_tr = c_p1_transposed;
                _qterms->insert(make_pair<>(qname, move(q)));
                if(p_new1->is_var()){
                    _evaluated = false;
                }
                //            update_convexity();
                return true;
            }
            else {
                if (pair_it->second._sign == sign) {
                    if (coef.is_function()) {
                        auto coef2 = *(func<type>*)(&coef);
                        pair_it->second._coef = add(pair_it->second._coef,coef2);
                    }
                    else if(coef.is_param()) {
                        auto coef2 = *(param<type>*)(&coef);
                        pair_it->second._coef = add(pair_it->second._coef,coef2);
                    }
                    else if(coef.is_number()) {
                        auto coef2 = *(constant<type>*)(&coef);
                        pair_it->second._coef = add(pair_it->second._coef,coef2);
                    }
                }
                else{
                    if (coef.is_function()) {
                        auto coef2 = *(func<type>*)(&coef);
                        pair_it->second._coef = subtract(pair_it->second._coef,coef2);
                    }
                    else if(coef.is_param()) {
                        auto coef2 = *(param<type>*)(&coef);
                        pair_it->second._coef = subtract(pair_it->second._coef,coef2);
                    }
                    else if(coef.is_number()) {
                        auto coef2 = *(constant<type>*)(&coef);
                        pair_it->second._coef = subtract(pair_it->second._coef,coef2);
                    }
                }
                if (pair_it->second._coef->is_zero()) {
                    if (p1.is_var()) {
                        decr_occ_var(ps1);
                    }
                    else {
                        decr_occ_param(ps1);
                    }
                    if (p2.is_var()) {
                        decr_occ_var(ps2);
                    }
                    else {
                        decr_occ_param(ps2);
                    }
                    _qterms->erase(pair_it);
                    if(_qterms->empty()){
                        _ftype = lin_;
                    }
                    if(is_constant()){
                        _ftype = const_;
                        _val->resize(1);
                    }
                    //                update_sign();
                    //                update_convexity();
                }
                //            else {
                //                update_sign(pair_it->second);
                //                update_convexity(pair_it->second);
                //            }
                return false;
            }
        };/**< Adds coef*p1*p2 to the function. Returns true if added new term, false if only updated coef of p1*p2 */
        
        bool insert(bool sign, const constant_& coef, const list<pair<shared_ptr<param_>, int>>& l){/**< Adds polynomial term to the function. Returns true if added new term, false if only updated corresponding coef */
            _all_convexity = undet_;
            string name;
            string s;
            bool newv = true;
            _evaluated=false;
            //        int i = 0;
            for (auto &pair:l) {
                name += pair.first->get_name(false,false);
                name += "^"+to_string(pair.second);
                name += ",";
            }
            auto pair_it = _pterms->find(name);
            auto p = l.begin()->first;
            shared_ptr<param_> pnew;
            if (_ftype <= quad_ && p->is_var()) {
                _ftype = pol_;
            }
            if (pair_it == _pterms->end()) {
                auto newl = make_shared<list<pair<shared_ptr<param_>, int>>>();
                //            i = 1;
                for (auto &pair:l) {
                    p = pair.first;
                    s = p->get_name(false,false);
                    if (p->is_var()) {
                        pnew = get_var(s);
                        if (!pnew) {
                            pnew = p->pcopy();
                            add_var(pnew,pair.second);
                        }
                        else {
                            incr_occ_var(s);
                        }
                    }
                    else {
                        pnew = get_param(s);
                        if (!pnew) {
                            pnew = p->pcopy();
                            add_param(pnew);
                        }
                        else {
                            incr_occ_param(s);
                        }
                    }
                    newv = true;
                    for (auto& p_it:*newl) {
                        if (p_it.first->get_name(false,false)==s) {
                            p_it.second++;
                            newv = false;
                            break;
                        }
                    }
                    if (newv) {
                        newl->push_back(make_pair<>(pnew, pair.second));
                    }
                }
                auto c_new = coef.copy();
                if (c_new->is_function()) {
                    embed(*static_pointer_cast<func>(c_new));
                }
                pterm p(sign, c_new, newl);
                _pterms->insert(make_pair<>(name, move(p)));
                if(pnew->is_var()){
                    _evaluated = false;
                }
                return true;
            }
            else {
                if (pair_it->second._sign == sign) {
                    if (coef.is_function()) {
                        auto coef2 = *(func<type>*)(&coef);
                        pair_it->second._coef = add(pair_it->second._coef,coef2);
                    }
                    else if(coef.is_param()) {
                        auto coef2 = *(param<type>*)(&coef);
                        pair_it->second._coef = add(pair_it->second._coef,coef2);
                    }
                    else if(coef.is_number()) {
                        auto coef2 = *(constant<type>*)(&coef);
                        pair_it->second._coef = add(pair_it->second._coef,coef2);
                    }
                }
                else{
                    if (coef.is_function()) {
                        auto coef2 = *(func<type>*)(&coef);
                        pair_it->second._coef = subtract(pair_it->second._coef,coef2);
                    }
                    else if(coef.is_param()) {
                        auto coef2 = *(param<type>*)(&coef);
                        pair_it->second._coef = subtract(pair_it->second._coef,coef2);
                    }
                    else if(coef.is_number()) {
                        auto coef2 = *(constant<type>*)(&coef);
                        pair_it->second._coef = subtract(pair_it->second._coef,coef2);
                    }
                }
                if (pair_it->second._coef->is_zero()) {
                    for (auto& it:*pair_it->second._l) {
                        p = it.first;
                        s = p->get_name(false,false);
                        if (p->is_var()) {
                            decr_occ_var(s,it.second);
                        }
                        else {
                            decr_occ_param(s,it.second);
                        }
                    }
                    _pterms->erase(pair_it);
                    if(_pterms->empty()){
                        _ftype = quad_;
                    }
                    if(_qterms->empty()){
                        _ftype = lin_;
                    }
                    if(is_constant()){
                        _ftype = const_;
                        _val->resize(1);
                    }
                    //                update_sign();
                    //                update_convexity();
                }
                //            else {
                //                update_sign(pair_it->second);
                //            }
                return false;
            }
            
        };/**< Adds polynomial term to the function. Returns true if added new term, false if only updated corresponding coef */
        
        bool insert(const constant_& coef, const param_& p, int exp){
            return insert(true, coef, p, exp);
        };
        
        bool insert(bool sign, const constant_& coef, const param_& p, int exp){
            list<pair<shared_ptr<param_>, int>> l;
            l.push_back(make_pair<>(p.pcopy(), exp));
            return insert(sign, coef, l);
        };/**< Adds polynomial term to the function. Returns true if added new term, false if only updated corresponding coef */
        
        
        bool insert(const qterm& term){
            return insert(term._sign, *term._coef, *term._p->first, *term._p->second, term._coef_p1_tr);
        };
        
        bool insert(const pterm& term){return insert(term._sign, *term._coef, *term._l);};
        
        
    };
    
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator+(const func<T1>& f1, const func<T2>& f2){
        return func<T1>(f1)+= f2;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator+(const func<T1>& f1, const func<T2>& f2){
        return func<T2>(f1)+= f2;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator-(const func<T1>& f1, const func<T2>& f2){
        return func<T1>(f1)-= f2;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator-(const func<T1>& f1, const func<T2>& f2){
        return func<T2>(f1)-= f2;
    }
    //
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator*(const func<T1>& f1, const func<T2>& f2){
        if((f1._is_transposed || f1.is_matrix()) && !f2._is_vector){
            return func<T1>(f1)*= f2.vec();
        }
        return func<T1>(f1)*= f2;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator*(const func<T1>& f1, const func<T2>& f2){
        if((f1._is_transposed || f1.is_matrix()) && !f2._is_vector){
            return func<T2>(f1)*= f2.vec();
        }
        return func<T2>(f1)*= f2;
    }
    
    
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator*(const param<T1>& p1, const param<T2>& p2){
        func<T1> res;
        if(p1.is_zero() || p2.is_zero()){
            return res;
        }
        if(p1.is_param() && p2.is_var()){
            if(p1._is_transposed && ! p2._is_vector){
                res.insert(true,p1,p2.vec());
            }
            else {
                res.insert(true,p1,p2);
            }
            res.update_dot_dim(p1,p2);
        }
        else if(p2.is_param() && p1.is_var()){
            if(p1._is_transposed && p2.is_row_vector()){/* transform x^T*p to (p^T*x)^T */
                auto new_p2 = p2.tr();
                auto new_p1 = p1.tr();
                res.insert(true,new_p2,new_p1);
                res.update_dot_dim(new_p2,new_p1);
                res.transpose();
            }
            else {
                res.insert(true,p2,p1);
                res.update_dot_dim(p1,p2);
            }
        }
        else {//Both vars or both params
            if(p1._is_transposed && !p2._is_vector){
                res.insert(true,unit<T1>(),p1,p2.vec());
            }
            else {
                res.insert(true,unit<T1>(),p1,p2);
            }
            res.update_dot_dim(p1,p2);
        }
        
        if(res.has_square()){
            auto signp = p1.get_all_sign();
            if(signp==neg_ || signp==pos_){
                res._all_sign = pos_;
            }
            else {
                res._all_sign = non_neg_;
            }
            res._range->first=zero<T1>().eval();
            if(p1.is_positive() || p1.is_negative()){
                res._range->first=extended_mult(p1._range->first,p1._range->first);
            }
            res._range->second=extended_mult(std::max(std::abs(p1._range->second),std::abs(p1._range->first)),std::max(std::abs(p1._range->second),std::abs(p1._range->first)));
        }
        else {
            res._range = get_product_range(p1._range,p2._range);
            res._all_sign = sign_product(p1.get_all_sign(), p2.get_all_sign());
        }
        if(res.is_quadratic()){res.update_quad_convexity();}
        if(p1._is_transposed){
            res._range->first = extended_mult(res._range->first,(T1)p1._dim[0]);
            res._range->second = extended_mult(res._range->second,(T1)p1._dim[0]);
        }
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator*(const param<T1>& p1, const param<T2>& p2){
        func<T2> res;
        if(p1.is_zero() || p2.is_zero()){
            return res;
        }
        if(p1.is_param() && p2.is_var()){
            if(p1._is_transposed && ! p2._is_vector){
                res.insert(true,param<T2>(p1),p2.vec());
            }
            else {
                res.insert(true,param<T2>(p1),p2);
            }
            res.update_dot_dim(p1,p2);
        }
        else if(p2.is_param() && p1.is_var()){
            if(p1._is_transposed && p2.is_row_vector()){/* transform x^T*p to (p^T*x)^T */
                auto new_p2 = p2.tr();
                auto new_p1 = p1.tr();
                res.insert(true,new_p2,new_p1);
                res.update_dot_dim(new_p2,new_p1);
                res.transpose();
            }
            else {
                res.insert(true,p2,p1);
                res.update_dot_dim(p1,p2);
            }
        }
        else {//Both vars or both params
            if(p1._is_transposed && !p2._is_vector){
                res.insert(true,unit<T2>(),p1,p2.vec());
            }
            else {
                res.insert(true,unit<T2>(),p1,p2);
            }
            res.update_dot_dim(p1,p2);
        }
        
        if(res.has_square()){
            auto signp = p1.get_all_sign();
            if(signp==neg_ || signp==pos_){
                res._all_sign = pos_;
            }
            else {
                res._all_sign = non_neg_;
            }
            res._range->first=zero<T2>().eval();
            if(p1.is_positive() || p1.is_negative()){
                res._range->first=extended_mult(p1._range->first,p1._range->first);
            }
            res._range->second=extended_mult(std::max(std::abs(p1._range->second),std::abs(p1._range->first)),std::max(std::abs(p1._range->second),std::abs(p1._range->first)));
        }
        else {
            res._range = get_product_range(p1._range,p2._range);
            res._all_sign = sign_product(p1.get_all_sign(), p2.get_all_sign());
        }
        if(res.is_quadratic()){res.update_quad_convexity();}
        if(p1._is_transposed){
            res._range->first = extended_mult(res._range->first,(T2)p1._dim[0]);
            res._range->second = extended_mult(res._range->second,(T2)p1._dim[0]);
        }
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) <= sizeof(T1)>::type* = nullptr>
    func<T1> operator+(const param<T1>& p1, const param<T2>& p2){
        func<T1> res;
        res.set_max_dim(p1,p2);
        if(p1.is_param() && p2.is_var()){
            res.insert(true,unit<T1>(),p2);
            res.add_cst(p1);
        }
        else if(p2.is_param() && p1.is_var()){
            res.insert(true,unit<T1>(),p1);
            res.add_cst(param<T1>(p2));
        }
        else {//Both vars or both params
            res.insert(true,unit<T1>(),p1);
            res.insert(true,unit<T1>(),p2);
        }
        res._all_sign = sign_add(p1.get_all_sign(), p2.get_all_sign());
        if(res.is_quadratic()){res.update_quad_convexity();}
        res._range = get_plus_range(p1._range,p2._range);
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T1) < sizeof(T2)>::type* = nullptr>
    func<T2> operator+(const param<T1>& p1, const param<T2>& p2){
        func<T2> res;
        res.set_max_dim(p1,p2);
        if(p1.is_param() && p2.is_var()){
            res.insert(true,unit<T2>(),p2);
            res.add_cst(param<T2>(p1));
        }
        else if(p2.is_param() && p1.is_var()){
            res.insert(true,unit<T2>(),p1);
            res.add_cst(p2);
        }
        else {//Both vars or both params
            res.insert(true,unit<T2>(),p1);
            res.insert(true,unit<T2>(),p2);
        }
        res._all_sign = sign_add(p1.get_all_sign(), p2.get_all_sign());
        if(res.is_quadratic()){res.update_quad_convexity();}
        res._range = get_plus_range(p1._range,p2._range);
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator-(const param<T1>& p1, const param<T2>& p2){
        func<T1> res;
        res.set_max_dim(p1,p2);
        if(p1.is_param() && p2.is_var()){
            res.insert(false,unit<T1>(),p2);
            res.add_cst(p1);
        }
        else if(p2.is_param() && p1.is_var()){
            res.insert(true,unit<T1>(),p1);
            func<T1> newp(p2);
            newp.reverse_sign();
            res.add_cst(newp);
        }
        else {//Both vars or both params
            res.insert(true,unit<T1>(),p1);
            res.insert(false,unit<T1>(),p2);
        }
        res._all_sign = sign_add(p1.get_all_sign(), reverse(p2.get_all_sign()));
        if(res.is_quadratic()){res.update_quad_convexity();}
        res._range = get_minus_range(p1._range, p2._range);
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator-(const param<T1>& p1, const param<T2>& p2){
        func<T2> res;
        res.set_max_dim(p1,p2);
        if(p1.is_param() && p2.is_var()){
            res.insert(false,unit<T2>(),p2);
            res.add_cst(param<T2>(p1));
        }
        else if(p2.is_param() && p1.is_var()){
            res.insert(true,unit<T2>(),p1);
            func<T2> newp(p2);
            newp.reverse_sign();
            res.add_cst(newp);
        }
        else {//Both vars or both params
            res.insert(true,unit<T2>(),p1);
            res.insert(false,unit<T2>(),p2);
        }
        res._all_sign = sign_add(p1.get_all_sign(), reverse(p2.get_all_sign()));
        if(res.is_quadratic()){res.update_quad_convexity();}
        res._range = get_minus_range(p1._range,p2._range);
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) <= sizeof(T1)>::type* = nullptr>
    constant<T1> operator+(const constant<T1>& p1, const constant<T2>& p2){
        constant<T1> res(p1);
        res.set_val(res.eval()+p2.eval());
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T1) < sizeof(T2)>::type* = nullptr>
    constant<T2> operator+(const constant<T1>& p1, const constant<T2>& p2){
        constant<T2> res(p2);
        res.set_val(res.eval()+(T2)p1.eval());
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) <= sizeof(T1)>::type* = nullptr>
    constant<T1> operator-(const constant<T1>& p1, const constant<T2>& p2){
        constant<T1> res(p1);
        res.set_val(res.eval()-p2.eval());
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T1) < sizeof(T2)>::type* = nullptr>
    constant<T2> operator-(const constant<T1>& p1, const constant<T2>& p2){
        constant<T2> res(p2);
        res.set_val((T2)p1.eval() - res.eval());
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) <= sizeof(T1)>::type* = nullptr>
    constant<T1> operator*(const constant<T1>& p1, const constant<T2>& p2){
        constant<T1> res(p1);
        res.set_val(res.eval()*p2.eval());
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T1) < sizeof(T2)>::type* = nullptr>
    constant<T2> operator*(const constant<T1>& p1, const constant<T2>& p2){
        constant<T2> res(p2);
        res.set_val((T2)p1.eval() * res.eval());
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) <= sizeof(T1)>::type* = nullptr>
    func<T1> operator/(const func<T1>& f1, const func<T2>& f2){
        return func<T1>(f1)/=f2;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T1) < sizeof(T2)>::type* = nullptr>
    func<T2> operator/(const func<T1>& f1, const func<T2>& f2){
        return func<T2>(f1)/=f2;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) <= sizeof(T1)>::type* = nullptr>
    func<T1> operator/(const param<T1>& f1, const param<T2>& f2){
        return func<T1>(f1)/=f2;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T1) < sizeof(T2)>::type* = nullptr>
    func<T2> operator/(const param<T1>& f1, const param<T2>& f2){
        return func<T2>(f1)/=f2;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) <= sizeof(T1)>::type* = nullptr>
    func<T1> operator/(const func<T1>& f1, const param<T2>& p2){
        return func<T1>(f1)/=p2;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T1) < sizeof(T2)>::type* = nullptr>
    func<T2> operator/(const func<T1>& f1, const param<T2>& p2){
        return func<T2>(f1)/=p2;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) <= sizeof(T1)>::type* = nullptr>
    func<T1> operator/(T1 f1, const param<T2>& p2){
        return func<T1>(f1)/=p2;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T1) < sizeof(T2)>::type* = nullptr>
    func<T2> operator/(T1 f1, const param<T2>& p2){
        return func<T2>(f1)/=p2;
    }
    
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) <= sizeof(T1)>::type* = nullptr>
    func<T1> operator/(const param<T1>& p1, const func<T2>& f2){
        return func<T1>(p1)/=f2;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T1) < sizeof(T2)>::type* = nullptr>
    func<T2> operator/(const param<T1>& p1, const func<T2>& f2){
        return func<T2>(p1)/=f2;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) <= sizeof(T1)>::type* = nullptr>
    func<T1> operator/(const func<T1>& f1, const constant<T2>& p2){
        return func<T1>(f1)/=p2;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T1) < sizeof(T2)>::type* = nullptr>
    func<T2> operator/(const func<T1>& f1, const constant<T2>& p2){
        return func<T2>(f1)/=p2;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) <= sizeof(T1)>::type* = nullptr>
    func<T1> operator/(const constant<T1>& p1, const func<T2>& f2){
        return func<T1>(p1)/=f2;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T1) < sizeof(T2)>::type* = nullptr>
    func<T2> operator/(const constant<T1>& p1, const func<T2>& f2){
        return func<T2>(p1)/=f2;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) <= sizeof(T1)>::type* = nullptr>
    func<T1> operator/(const func<T1>& f1, const T2 p2){
        return func<T1>(f1)/=p2;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T1) < sizeof(T2)>::type* = nullptr>
    func<T2> operator/(const func<T1>& f1, T2 p2){
        return func<T2>(f1)/=p2;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) <= sizeof(T1)>::type* = nullptr>
    func<T1> operator/(const T1 p1, const func<T2>& f2){
        return func<T1>(p1)/=f2;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T1) < sizeof(T2)>::type* = nullptr>
    func<T2> operator/(const T1 p1, const func<T2>& f2){
        return func<T2>(p1)/=f2;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) <= sizeof(T1)>::type* = nullptr>
    func<T1> operator/(const param<T1>& f1, const constant<T2>& p2){
        return func<T1>(f1)/=p2;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T1) < sizeof(T2)>::type* = nullptr>
    func<T2> operator/(const param<T1>& f1, const constant<T2>& p2){
        return func<T2>(f1)/=p2;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) <= sizeof(T1)>::type* = nullptr>
    func<T1> operator/(const constant<T1>& p1, const param<T2>& f2){
        return func<T1>(p1)/=f2;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T1) < sizeof(T2)>::type* = nullptr>
    func<T2> operator/(const constant<T1>& p1, const param<T2>& f2){
        return func<T2>(p1)/=f2;
    }
    
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) <= sizeof(T1)>::type* = nullptr>
    constant<T1> operator/(const constant<T1>& p1, const constant<T2>& p2){
        constant<T1> res(p1);
        res.set_val(res.eval()/p2.eval());
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T1) < sizeof(T2)>::type* = nullptr>
    constant<T2> operator/(const constant<T1>& p1, const constant<T2>& p2){
        constant<T2> res(p2);
        res.set_val((T2)p1.eval()/res.eval());
        return res;
    }
    
    
    template<class T1>
    constant<T1> min(const constant<T1>& p1,const constant<T1>& p2){
        constant<T1> res(p1);
        res.set_val(gravity::min(res.eval(),p2.eval()));
        return res;
    }
    
    template<class T1>
    constant<T1> max(const constant<T1>& p1,const constant<T1>& p2){
        constant<T1> res(p1);
        res.set_val(gravity::max(res.eval(),p2.eval()));
        return res;
    }
    
    template<class T1>
    constant<T1> log(const constant<T1>& p1){
        constant<T1> res(p1);
        res.set_val(log(res.eval()));
        return res;
    }
    
    template<class T1>
    constant<T1> exp(const constant<T1>& p1){
        constant<T1> res(p1);
        res.set_val(exp(res.eval()));
        return res;
    }
    
    template<class T1>
    constant<T1> sqrt(const constant<T1>& p1){
        constant<T1> res(p1);
        res.set_val(sqrt(res.eval()));
        return res;
    }
    
    template<class T1>
    constant<T1> unit_step(const constant<T1>& p1){
        if(p1.is_non_positive()){
            return zero<T1>();
        }
        return unit<T1>();
    }
    
    template<class T1>
    constant<T1> cos(const constant<T1>& p1){
        constant<T1> res(p1);
        res.set_val(cos(res.eval()));
        return res;
    }
    
    template<class T1>
    constant<T1> sin(const constant<T1>& p1){
        constant<T1> res(p1);
        res.set_val(sin(res.eval()));
        return res;
    }
    
    template<class T1>
    constant<T1> tan(const constant<T1>& p1){
        constant<T1> res(p1);
        res.set_val(tan(res.eval()));
        return res;
    }
    
    template<class T1>
    constant<T1> pow(const constant<T1>& p1, int exp){
        constant<T1> res(p1);
        res.set_val(std::pow(res.eval(),exp));
        return res;
    }
    
    /**
     Return the sign and curvature of unitary operator op on given range.
     @param[in] op Mathematical unitary operator, e.g., cos, sin, log, etc...
     @param[in] range, range of values we're interested in.
     @return a pair<Convexity,Sign> charachterizing the sign and the curvature of the operator op in the given range
     */
    template<class T, typename enable_if<is_arithmetic<T>::value>::type* = nullptr>
    pair<Convexity,Sign> get_sign_curvature(OperatorType op, pair<T,T> range){
        pair<Convexity,Sign> res = {undet_,unknown_};
        switch(op){
            case id_:{// f(x) = x
                res.first = linear_;
                if(range.first>=0){
                    res.second = non_neg_;
                }
                if(range.first>0){
                    res.second = pos_;
                }
                if(range.second<=0){
                    res.second = non_pos_;
                }
                if(range.first<0){
                    res.second = neg_;
                }
                return res;
            }
            case cos_:
                return cos_sign_curvature(range);
            case sin_:{
                range.first += pi/2.;
                range.second += pi/2.;
                return cos_sign_curvature(range);
            }
            case log_:{
                res.first = concave_;
                if(range->first>=1){
                    res.second = non_neg_;
                    if(range->first>1){
                        res.second = pos_;
                    }
                }
                if(range->second<=1){
                    res.second = non_pos_;
                    if(range->second<1){
                        res.second = neg_;
                    }
                }
                return res;
            }
            case exp_:{
                return {convex_,pos_};
            }
            case sqrt_:{
                res.first = concave_;
                res.second = non_neg_;
                if(range.first > 0){
                    res.second = pos_;
                }
                return res;
            }
            case tan_:{
                if(range.first>=0){
                    res.first = convex_;
                    res.second = non_neg_;
                }
                if(range.first>0){
                    res.first = convex_;
                    res.second = pos_;
                }
                if(range.second<=0){
                    res.first = concave_;
                    res.second = non_pos_;
                }
                if(range.second<0){
                    res.first = concave_;
                    res.second = neg_;
                }
                return res;
            }
            case relu_:{
                res.first = convex_;
                res.second = non_neg_;
                if(range.first > 0){
                    res.second = pos_;
                }
                return res;
            }
            default:
                break;
        }
        return res;
    }
    
    template<class T>
    func<T> min(const param<T>& p1, const param<T>& p2){
        func<T> res(bexpr<T>(min_, p1.copy(), p2.copy()));
        res._all_sign = std::min(p1.get_all_sign(),p2.get_all_sign());
        res._all_convexity = undet_;
        res._range->first = gravity::min(p1._range->first,p2._range->first);
        res._range->second = gravity::min(p1._range->second,p2._range->second);
        res._expr->_range->first = res._range->first;
        res._expr->_range->second = res._range->second;
        res._expr->_all_convexity = res._all_convexity;
        res._expr->_all_sign = res._all_sign;
        return res;
    }
    
    template<class T>
    func<T> max(const param<T>& p1, const param<T>& p2){
        func<T> res(bexpr<T>(max_, p1.copy(), p2.copy()));
        res._all_sign = std::max(p1.get_all_sign(),p2.get_all_sign());
        res._all_convexity = undet_;
        res._range->first = gravity::max(p1._range->first,p2._range->first);
        res._range->second = gravity::max(p1._range->second,p2._range->second);
        res._expr->_range->first = res._range->first;
        res._expr->_range->second = res._range->second;
        res._expr->_all_convexity = res._all_convexity;
        res._expr->_all_sign = res._all_sign;
        return res;
    }
    
    template<class T>
    func<T> min(const param<T>& p1, const func<T>& p2){
        func<T> res(bexpr<T>(min_, p1.copy(), p2.copy()));
        res._all_sign = std::min(p1.get_all_sign(),p2.get_all_sign());
        res._all_convexity = undet_;
        res._range->first = gravity::min(p1._range->first,p2._range->first);
        res._range->second = gravity::min(p1._range->second,p2._range->second);
        res._expr->_range->first = res._range->first;
        res._expr->_range->second = res._range->second;
        res._expr->_all_convexity = res._all_convexity;
        res._expr->_all_sign = res._all_sign;
        return res;
    }
    
    template<class T>
    func<T> min(const func<T>& p1, const param<T>& p2){
        func<T> res(bexpr<T>(min_, p1.copy(), p2.copy()));
        res._all_sign = std::min(p1.get_all_sign(),p2.get_all_sign());
        res._all_convexity = undet_;
        res._range->first = gravity::min(p1._range->first,p2._range->first);
        res._range->second = gravity::min(p1._range->second,p2._range->second);
        res._expr->_range->first = res._range->first;
        res._expr->_range->second = res._range->second;
        res._expr->_all_convexity = res._all_convexity;
        res._expr->_all_sign = res._all_sign;
        return res;
    }
    
    template<class T>
    func<T> max(const param<T>& p1, const func<T>& p2){
        func<T> res(bexpr<T>(max_, p1.copy(), p2.copy()));
        res._all_sign = std::max(p1.get_all_sign(),p2.get_all_sign());
        res._all_convexity = undet_;
        res._range->first = gravity::max(p1._range->first,p2._range->first);
        res._range->second = gravity::max(p1._range->second,p2._range->second);
        res._expr->_range->first = res._range->first;
        res._expr->_range->second = res._range->second;
        res._expr->_all_convexity = res._all_convexity;
        res._expr->_all_sign = res._all_sign;
        return res;
    }
    
    template<class T>
    func<T> max(const func<T>& p1, const param<T>& p2){
        func<T> res(bexpr<T>(max_, p1.copy(), p2.copy()));
        res._all_sign = std::max(p1.get_all_sign(),p2.get_all_sign());
        res._all_convexity = undet_;
        res._range->first = gravity::max(p1._range->first,p2._range->first);
        res._range->second = gravity::max(p1._range->second,p2._range->second);
        res._expr->_range->first = res._range->first;
        res._expr->_range->second = res._range->second;
        res._expr->_all_convexity = res._all_convexity;
        res._expr->_all_sign = res._all_sign;
        return res;
    }
    
    template<class T>
    func<T> max(const func<T>& p1, const func<T>& p2){
        func<T> res(bexpr<T>(max_, p1.copy(), p2.copy()));
        res._all_sign = std::max(p1.get_all_sign(),p2.get_all_sign());
        res._all_convexity = undet_;
        res._range->first = gravity::max(p1._range->first,p2._range->first);
        res._range->second = gravity::max(p1._range->second,p2._range->second);
        res._expr->_range->first = res._range->first;
        res._expr->_range->second = res._range->second;
        res._expr->_all_convexity = res._all_convexity;
        res._expr->_all_sign = res._all_sign;
        return res;
    }
    
    template<class T>
    func<T> min(const func<T>& p1, const func<T>& p2){
        func<T> res(bexpr<T>(min_, p1.copy(), p2.copy()));
        res._all_sign = std::min(p1.get_all_sign(),p2.get_all_sign());
        res._all_convexity = undet_;
        res._range->first = gravity::min(p1._range->first,p2._range->first);
        res._range->second = gravity::min(p1._range->second,p2._range->second);
        res._expr->_range->first = res._range->first;
        res._expr->_range->second = res._range->second;
        res._expr->_all_convexity = res._all_convexity;
        res._expr->_all_sign = res._all_sign;
        return res;
    }
    
    template<class T, typename enable_if<is_arithmetic<T>::value>::type* = nullptr>
    func<T> log(const param<T>& p1){
        if(!p1.is_positive()){
            throw invalid_argument("Calling log() with a non-positive argument");
        }
        func<T> res(uexpr<T>(log_, p1.copy()));
        if(p1._range->first>=1){
            res._all_sign = non_neg_;
            if(p1._range->first>1){
                res._all_sign = pos_;
            }
        }
        if(p1._range->second<=1){
            res._all_sign = non_pos_;
            if(p1._range->second<1){
                res._all_sign = neg_;
            }
        }
        if (p1.is_var()) {
            res._all_convexity = concave_;
        }
        res._range->first = std::log(p1._range->first);
        res._range->second = std::log(p1._range->second);
        res._expr->_range->first = res._range->first;
        res._expr->_range->second = res._range->second;
        res._expr->_all_convexity = res._all_convexity;
        res._expr->_all_sign = res._all_sign;
        return res;
    }
    
    template<class T1>
    func<T1> exp(const param<T1>& p1){
        func<T1> res(uexpr<T1>(exp_, p1.copy()));
        res._all_sign = pos_;
        if (p1.is_var()) {
            res._all_convexity = convex_;
        }
        res._range->first = std::exp(p1._range->first);
        res._range->second = std::exp(p1._range->second);
        res._expr->_range->first = res._range->first;
        res._expr->_range->second = res._range->second;
        res._expr->_all_convexity = res._all_convexity;
        res._expr->_all_sign = res._all_sign;
        return res;
    }
    
    template<class T1>
    func<T1> sqrt(const param<T1>& p1){
        if(!p1.is_non_negative()){
            throw invalid_argument("Calling sqrt() with a negative argument");
        }
        func<T1> res(uexpr<T1>(sqrt_, p1.copy()));
        res._all_sign = non_neg_;
        if(p1.is_positive()){
            res._all_sign = pos_;
        }
        if (p1.is_var()) {
            res._all_convexity = concave_;
        }
        res._range->first = std::sqrt(p1._range->first);
        res._range->second = std::sqrt(p1._range->second);
        res._expr->_range->first = res._range->first;
        res._expr->_range->second = res._range->second;
        res._expr->_all_convexity = res._all_convexity;
        res._expr->_all_sign = res._all_sign;
        return res;
    }
    
    template<class T, typename enable_if<is_arithmetic<T>::value>::type* = nullptr>
    pair<Convexity,Sign> cos_sign_curvature(const pair<T,T>& range){
        if(range.first==numeric_limits<T>::lowest() || range.second==numeric_limits<T>::max()){
            return {undet_,unknown_};
        }
        pair<Convexity,Sign> res = {undet_,zero_};
        auto lb = fmod(range.first,(2*pi));
        auto ub = fmod(range.second,(2*pi));
        if(ub<= -3*pi/2){
            res.first = concave_;
            res.second = non_neg_;
            if(ub< -3*pi/2){
                res.second = pos_;
            }
        }
        if(lb>=-3*pi/2 && ub<= -pi/2){
            res.first = convex_;
            res.second = non_pos_;
            if(lb>-3*pi/2 && ub<-pi/2){
                res.second = neg_;
            }
        }
        if(lb>=-pi/2 && ub<= pi/2){
            res.first = concave_;
            res.second = non_neg_;
            if(lb>-pi/2 && ub<pi/2){
                res.second = pos_;
            }
        }
        if(lb>=pi/2 && ub<= 3*pi/2){
            res.first = convex_;
            res.second = non_pos_;
            if(lb>pi/2 && ub<3*pi/2){
                res.second = neg_;
            }
        }
        if(lb >= 3*pi/2){
            res.first = concave_;
            res.second = non_neg_;
            if(lb > 3*pi/2){
                res.second = pos_;
            }
        }
        return res;
    }
    
    template<class T1>
    func<T1> unit_step(const param<T1>& p1){
        func<T1> res(uexpr<T1>(unit_step_, p1.copy()));
        if(p1.is_non_positive()){
            res._range->first = zero<T1>().eval();
            res._range->second = zero<T1>().eval();
        }
        else if(p1.is_positive()){
            res._range->first = unit<T1>().eval();
            res._range->second = unit<T1>().eval();
        }
        else {
            res._range->first = zero<T1>().eval();
            res._range->second = unit<T1>().eval();
        }
        res._all_convexity = undet_;
        res._expr->_range->first = res._range->first;
        res._expr->_range->second = res._range->second;
        res._expr->_all_convexity = res._all_convexity;
        res._expr->_all_sign = res._all_sign;
        return res;
    }
    
    template<class T, typename enable_if<is_arithmetic<T>::value>::type* = nullptr>
    func<T> cos(const param<T>& p1){
        func<T> res(uexpr<T>(cos_, p1.copy()));
        auto conv_sign = cos_sign_curvature(*p1._range);
        if (p1.is_var()) {
            res._all_convexity = conv_sign.first;
        }
        res._all_sign = conv_sign.second;
        res._range->first = gravity::min(std::cos(p1._range->first),std::cos(p1._range->second));
        res._range->second = gravity::max(std::cos(p1._range->first),std::cos(p1._range->second));
        if(p1._range->first <0 && p1._range->second >0){
            res._range->second = 1;
        }
        if((p1._range->first <-pi && p1._range->second >-pi) || (p1._range->first <pi && p1._range->second >pi)){
            res._range->first = -1;
        }
        res._expr->_range->first = res._range->first;
        res._expr->_range->second = res._range->second;
        res._expr->_all_convexity = res._all_convexity;
        res._expr->_all_sign = res._all_sign;
        res._indices = p1._indices;
        return res;
    }
    
    template<class T, typename enable_if<is_arithmetic<T>::value>::type* = nullptr>
    func<T> sin(const param<T>& p1){
        func<T> res(uexpr<T>(sin_, p1.copy()));
        auto shifted_range = *p1._range;
        shifted_range.first += pi/2.;
        shifted_range.second += pi/2.;
        auto conv_sign = cos_sign_curvature(shifted_range);
        if (p1.is_var()) {
            res._all_convexity = conv_sign.first;
        }
        res._all_sign = conv_sign.second;
        res._range->first = gravity::min(std::sin(p1._range->first),std::sin(p1._range->second));
        res._range->second = gravity::max(std::sin(p1._range->first),std::sin(p1._range->second));
        if(shifted_range.first <0 && shifted_range.second >0){
            res._range->second = 1;
        }
        if((shifted_range.first <-pi && shifted_range.second >-pi) || (shifted_range.first <pi && shifted_range.second >pi)){
            res._range->first = -1;
        }
        res._expr->_range->first = res._range->first;
        res._expr->_range->second = res._range->second;
        res._expr->_all_convexity = res._all_convexity;
        res._expr->_all_sign = res._all_sign;
        res._indices = p1._indices;
        return res;
    }
    
    template<class T, typename enable_if<is_arithmetic<T>::value>::type* = nullptr>
    func<T> tan(const param<T>& p1){
        auto centered_range = *p1._range;
        centered_range.first= fmod(centered_range.first, 2*pi);
        centered_range.second=fmod(centered_range.second, 2*pi);
        if(centered_range.first<=-pi/2 || centered_range.second>=pi/2){
            throw invalid_argument("Calling tan() with discontinuous domain");
        }
        func<T> res(uexpr<T>(tan_, p1.copy()));
        if(centered_range.first>=0){
            if (p1.is_var()) {
                res._all_convexity = convex_;
            }
            res._all_sign = non_neg_;
            if(centered_range.first>0){
                res._all_sign = pos_;
            }
        }
        if(centered_range.second<=0){
            if (p1.is_var()) {
                res._all_convexity = concave_;
            }
            res._all_sign = non_pos_;
            if(centered_range.first>0){
                res._all_sign = neg_;
            }
        }
        res._range->first = std::tan(p1._range->first);
        res._range->second = std::tan(p1._range->second);
        res._expr->_range->first = res._range->first;
        res._expr->_range->second = res._range->second;
        res._expr->_all_convexity = res._all_convexity;
        res._expr->_all_sign = res._all_sign;
        res._indices = p1._indices;
        return res;
    }
    
    template<class T1>
    func<T1> ReLU(const param<T1>& p1){
        func<T1> res(uexpr<T1>(relu_, p1.copy()));
        if (p1.is_var()) {
            res._all_convexity = convex_;
        }
        res._all_sign = non_neg_;
        res._range->first = zero<T1>().eval();
        if(p1.is_positive()){
            res._all_sign = pos_;
            res._range->first = p1._range->first;
        }
        res._range->second = gravity::max(zero<T1>().eval(),p1._range->second);
        res._expr->_range->first = res._range->first;
        res._expr->_range->second = res._range->second;
        res._expr->_all_convexity = res._all_convexity;
        res._expr->_all_sign = res._all_sign;
        res._indices = p1._indices;
        return res;
    }
    
    template<class T>
    func<T> pow(const param<T>& p1, int exp){
        if(exp<0){
            func<T> res;
            if(!p1.is_negative() && !p1.is_positive()){
                throw invalid_argument("Calling pow() with a negative exponent on an argument that  can be zero");
            }
            res.insert(p1,exp);
            return res;
        }
        if(exp==0){
            return func<T>();
        }
        if(exp==1){
            return func<T>(p1);
        }
        if(exp==2){
            return p1*p1;
        }
        else {
            func<T> res;
            res.insert(p1,exp);
            res.set_max_dim(p1);
            res._range->first = gravity::min(std::pow(p1._range->first,exp),std::pow(p1._range->second,exp));
            res._range->second = gravity::max(std::pow(p1._range->first,exp),std::pow(p1._range->second,exp));
            if(exp%2==0) {
                res._all_sign = non_neg_;
                if(p1.is_positive()){
                    res._all_sign = pos_;
                }
                if(p1._range->first <0 && p1._range->second >0){
                    res._range->first = 0;
                }
            }
            else {
                res._all_sign = p1.get_all_sign();
            }
            if (p1.is_var()) {
                if(exp%2==0) {
                    res._all_convexity = convex_;
                }
                else if(p1.is_non_negative()){
                    res._all_convexity = convex_;
                }
                else if(p1.is_non_positive()){
                    res._all_convexity = concave_;
                }
                else {
                    res._all_convexity = undet_;
                }
            }
            res._indices = p1._indices;
            return res;
        }
    }
    
    
    template<class T1>
    func<T1> log(const func<T1>& f){
        if(!f.is_positive()){
            throw invalid_argument("Calling log() with a potentially negative/zero argument");
        }
        func<T1> res(uexpr<T1>(log_, f.copy()));
        if(f._range->first>=1){
            res._all_sign = non_neg_;
            if(f._range->first>1){
                res._all_sign = pos_;
            }
        }
        if(f._range->second<=1){
            res._all_sign = non_pos_;
            if(f._range->second<1){
                res._all_sign = neg_;
            }
        }
        if (f.is_linear()) {
            res._all_convexity = concave_;
        }
        else if(!f.is_constant()){
            res._all_convexity = undet_;
        }
        res._range->first = log(f._range->first);
        if(f._range->second==numeric_limits<T1>::max()){
            res._range->second = numeric_limits<T1>::max();
        }
        else {
            res._range->second = log(f._range->second);
        }
        res._expr->_range->first = res._range->first;
        res._expr->_range->second = res._range->second;
        res._expr->_all_convexity = res._all_convexity;
        res._expr->_all_sign = res._all_sign;
        res._indices = f._indices;
        return res;
    }
    
    template<class T1>
    func<T1> exp(const func<T1>& f){
        func<T1> res(uexpr<T1>(exp_, f.copy()));
        res._all_sign = pos_;
        if (f.is_linear()) {
            res._all_convexity = convex_;
        }
        else if(!f.is_constant()){
            res._all_convexity = undet_;
        }
        if(f._range->first==numeric_limits<T1>::lowest() || f._range->second==numeric_limits<T1>::max()){
            res._range->first = numeric_limits<T1>::lowest();
            res._range->second = numeric_limits<T1>::max();
        }
        else {
            res._range->first = std::exp(f._range->first);
            res._range->second = std::exp(f._range->second);
        }
        res._expr->_range->first = res._range->first;
        res._expr->_range->second = res._range->second;
        res._expr->_all_convexity = res._all_convexity;
        res._expr->_all_sign = res._all_sign;
        res._indices = f._indices;
        return res;
    }
    
    template<class T1>
    func<T1> sqrt(const func<T1>& f){
        if(!f.is_non_negative()){
            throw invalid_argument("Calling sqrt() with a potentially negative argument");
        }
        func<T1> res(uexpr<T1>(sqrt_, f.copy()));
        res._all_sign = non_neg_;
        if(f.is_positive()){
            res._all_sign = pos_;
        }
        if (f.is_linear()) {
            res._all_convexity = convex_;
        }
        else if(!f.is_constant()){
            res._all_convexity = undet_;
        }
        res._range->first = std::sqrt(f._range->first);
        if(f._range->second==numeric_limits<T1>::max()){
            res._range->second = numeric_limits<T1>::max();
        }
        else {
            res._range->second = std::sqrt(f._range->second);
        }
        res._expr->_range->first = res._range->first;
        res._expr->_range->second = res._range->second;
        res._expr->_all_convexity = res._all_convexity;
        res._expr->_all_sign = res._all_sign;
        res._indices = f._indices;
        return res;
    }
    
    template<class T, typename enable_if<is_arithmetic<T>::value>::type* = nullptr>
    func<T> asin(const func<T>& f){
        if(f._range->first<-1 || f._range->second>1){
            throw invalid_argument("Calling asin(const func<T1>& f) outside [-1,1]");
        }
        func<T> res(uexpr<T>(asin_, f.copy()));
        if (f.is_linear()) {
            if(f.is_non_positive()){
                res._all_convexity = concave_;
            }
            else if(f.is_non_negative()){
                res._all_convexity = convex_;
            }
        }
        else if(!f.is_constant()){
            res._all_convexity = undet_;
        }
        res._all_sign = f._all_sign;
        res._range->first = std::asin(f._range->first);
        res._range->second = std::asin(f._range->second);
        res._expr->_range->first = res._range->first;
        res._expr->_range->second = res._range->second;
        res._expr->_all_convexity = res._all_convexity;
        res._expr->_all_sign = res._all_sign;
        res._indices = f._indices;
        return res;
    }
    
    template<class T>
    func<T> acos(const param<T>& p){
        return acos(func<T>(p));
    }
    
    template<class T>
    func<T> acos(const var<T>& p){
        return acos(func<T>(p));
    }
    
    template<class T>
    func<T> asin(const param<T>& p){
        return acos(func<T>(p));
    }
    
    template<class T>
    func<T> asin(const var<T>& p){
        return acos(func<T>(p));
    }
    
    template<class T, typename enable_if<is_arithmetic<T>::value>::type* = nullptr>
    func<T> acos(const func<T>& f){
        if(f._range->first<-1 || f._range->second>1){
            throw invalid_argument("Calling acos(const func<T1>& f) outside [-1,1]");
        }
        func<T> res(uexpr<T>(acos_, f.copy()));
        if (f.is_linear()) {
            if(f.is_non_positive()){
                res._all_convexity = convex_;
            }
            else if(f.is_non_negative()){
                res._all_convexity = concave_;
            }
        }
        else if(!f.is_constant()){
            res._all_convexity = undet_;
        }
        res._all_sign = non_neg_;
        res._range->first = std::acos(f._range->second);
        res._range->second = std::acos(f._range->first);
        res._expr->_range->first = res._range->first;
        res._expr->_range->second = res._range->second;
        res._expr->_all_convexity = res._all_convexity;
        res._expr->_all_sign = res._all_sign;
        res._indices = f._indices;
        return res;
    }
    
    template<class T, typename enable_if<is_arithmetic<T>::value>::type* = nullptr>
    func<T> cos(const func<T>& f){
        func<T> res(uexpr<T>(cos_, f.copy()));
        auto conv_sign = cos_sign_curvature(*f._range);
        if (f.is_linear()) {
            res._all_convexity = conv_sign.first;
        }
        else if(!f.is_constant()){
            res._all_convexity = undet_;
        }
        res._all_sign = conv_sign.second;
        if(f._range->first==numeric_limits<T>::lowest() || f._range->second==numeric_limits<T>::max()){
            res._range->first = -1;
            res._range->second = 1;
        }
        else {
            res._range->first = gravity::min(std::cos(f._range->first),std::cos(f._range->second));
            res._range->second = gravity::max(std::cos(f._range->first),std::cos(f._range->second));
        }
        if(f._range->first <0 && f._range->second >0){
            res._range->second = 1;
        }
        if((f._range->first <-pi && f._range->second >-pi) || (f._range->first <pi && f._range->second >pi)){
            res._range->first = -1;
        }
        res._expr->_range->first = res._range->first;
        res._expr->_range->second = res._range->second;
        res._expr->_all_convexity = res._all_convexity;
        res._expr->_all_sign = res._all_sign;
        res._indices = f._indices;
        return res;
    }
    
    template<class T, typename enable_if<is_same<T,Cpx>::value>::type* = nullptr>
    func<T> acos(const func<T>& f){
        func<T> res(uexpr<T>(acos_, f.copy()));
        res._all_convexity = undet_;
        res._expr->_range->first = res._range->first;
        res._expr->_range->second = res._range->second;
        res._expr->_all_convexity = res._all_convexity;
        res._expr->_all_sign = res._all_sign;
        res._indices = f._indices;
        return res;
    }
    
    template<class T, typename enable_if<is_same<T,Cpx>::value>::type* = nullptr>
    func<T> asin(const func<T>& f){
        func<T> res(uexpr<T>(asin_, f.copy()));
        res._all_convexity = undet_;
        res._expr->_range->first = res._range->first;
        res._expr->_range->second = res._range->second;
        res._expr->_all_convexity = res._all_convexity;
        res._expr->_all_sign = res._all_sign;
        res._indices = f._indices;
        return res;
    }
    
    template<class T, typename enable_if<is_same<T,Cpx>::value>::type* = nullptr>
    func<T> cos(const func<T>& f){
        func<T> res(uexpr<T>(cos_, f.copy()));
        res._all_convexity = undet_;
        res._expr->_range->first = res._range->first;
        res._expr->_range->second = res._range->second;
        res._expr->_all_convexity = res._all_convexity;
        res._expr->_all_sign = res._all_sign;
        res._indices = f._indices;
        return res;
    }
    
    template<class T, typename enable_if<is_arithmetic<T>::value>::type* = nullptr>
    func<T> sin(const func<T>& f){
        func<T> res(uexpr<T>(sin_, f.copy()));
        auto shifted_range = *f._range;
        auto conv_sign = cos_sign_curvature(shifted_range);
        if (f.is_linear()) {
            res._all_convexity = conv_sign.first;
        }
        else if(!f.is_constant()){
            res._all_convexity = undet_;
        }
        res._all_sign = conv_sign.second;
        if(f._range->first==numeric_limits<T>::lowest() || f._range->second==numeric_limits<T>::max()){
            res._range->first = -1;
            res._range->second = 1;
        }
        else {
            res._range->first = gravity::min(std::sin(f._range->first),std::sin(f._range->second));
            res._range->second = gravity::max(std::sin(f._range->first),std::sin(f._range->second));
            shifted_range.first += pi/2.;
            shifted_range.second += pi/2.;
            if(shifted_range.first <0 && shifted_range.second >0){
                res._range->second = 1;
            }
            if((shifted_range.first <-pi && shifted_range.second >-pi) || (shifted_range.first <pi && shifted_range.second >pi)){
                res._range->first = -1;
            }
        }
        res._expr->_range->first = res._range->first;
        res._expr->_range->second = res._range->second;
        res._expr->_all_convexity = res._all_convexity;
        res._expr->_all_sign = res._all_sign;
        res._indices = f._indices;
        return res;
    }
    
    template<class T, typename enable_if<is_arithmetic<T>::value>::type* = nullptr>
    func<T> atan2(const param<T>& f1, const param<T>& f2){
        return atan2(func<T>(f1), func<T>(f2));
    }
    
    template<class T, typename enable_if<is_arithmetic<T>::value>::type* = nullptr>
    func<T> atan2(const func<T>& f1, const func<T>& f2){
        func<T> res(bexpr<T>(atan2_, f1.copy(), f2.copy()));
        res._all_convexity = undet_;
        res._all_sign = unknown_;
        res._expr->_range->first = res._range->first;
        res._expr->_range->second = res._range->second;
        res._expr->_all_convexity = res._all_convexity;
        res._expr->_all_sign = res._all_sign;
        return res;
    }
    
    template<class T, typename enable_if<is_arithmetic<T>::value>::type* = nullptr>
    func<T> atan(const func<T>& f){
        //        if(f._range->first<=-pi/2 || f._range->second>=pi/2){
        //            throw invalid_argument("Calling atan(const func<T1>& f) outside ]-pi/2,pi/2[");
        //        }
        func<T> res(uexpr<T>(atan_, f.copy()));
        if (f.is_linear()) {
            if(f.is_non_positive()){
                res._all_convexity = convex_;
            }
            else if(f.is_non_negative()){
                res._all_convexity = concave_;
            }
        }
        else if(!f.is_constant()){
            res._all_convexity = undet_;
        }
        res._all_sign = f._all_sign;
        res._range->first = std::atan(f._range->second);
        res._range->second = std::atan(f._range->first);
        res._expr->_range->first = res._range->first;
        res._expr->_range->second = res._range->second;
        res._expr->_all_convexity = res._all_convexity;
        res._expr->_all_sign = res._all_sign;
        res._indices = f._indices;
        return res;
    }
    
    template<class T, typename enable_if<is_same<T,Cpx>::value>::type* = nullptr>
    func<T> atan(const func<T>& f){
        //        if(f._range->first<=-pi/2 || f._range->second>=pi/2){
        //            throw invalid_argument("Calling atan(const func<T1>& f) outside ]-pi/2,pi/2[");
        //        }
        func<T> res(uexpr<T>(atan_, f.copy()));
        if (f.is_linear()) {
            if(f.is_non_positive()){
                res._all_convexity = convex_;
            }
            else if(f.is_non_negative()){
                res._all_convexity = concave_;
            }
        }
        else if(!f.is_constant()){
            res._all_convexity = undet_;
        }
        res._all_sign = f._all_sign;
        res._range->first = std::atan(f._range->second);
        res._range->second = std::atan(f._range->first);
        res._expr->_range->first = res._range->first;
        res._expr->_range->second = res._range->second;
        res._expr->_all_convexity = res._all_convexity;
        res._expr->_all_sign = res._all_sign;
        res._indices = f._indices;
        return res;
    }
    
    template<class T, typename enable_if<is_same<T,Cpx>::value>::type* = nullptr>
    func<T> sin(const func<T>& f){
        func<T> res(uexpr<T>(sin_, f.copy()));
        res._all_convexity = undet_;
        res._expr->_range->first = res._range->first;
        res._expr->_range->second = res._range->second;
        res._expr->_all_convexity = res._all_convexity;
        res._expr->_all_sign = res._all_sign;
        res._indices = f._indices;
        return res;
    }
    
    template<class T, typename enable_if<is_arithmetic<T>::value>::type* = nullptr>
    func<T> tan(const func<T>& f){
        func<T> res(uexpr<T>(tan_, f.copy()));
        auto centered_range = *f._range;
        if(f._range->first==numeric_limits<T>::lowest() || f._range->second==numeric_limits<T>::max()){
            throw invalid_argument("Calling tan(const func<T1>& f) with discontinuous domain");
        }
        centered_range.first %= 2*pi;
        centered_range.second %= 2*pi;
        if(centered_range.first<=-pi/2 || centered_range.second>=pi/2){
            throw invalid_argument("Calling tan(const func<T1>& f) with discontinuous domain");
        }
        if(centered_range.first>=0){
            if (f.is_linear()) {
                res._all_convexity = convex_;
            }
            else if(!f.is_constant()){
                res._all_convexity = undet_;
            }
            res._all_sign = non_neg_;
            if(centered_range.first>0){
                res._all_sign = pos_;
            }
        }
        if(centered_range.second<=0){
            if (f.is_linear()) {
                res._all_convexity = concave_;
            }
            else if(!f.is_constant()){
                res._all_convexity = undet_;
            }
            res._all_sign = non_pos_;
            if(centered_range.first>0){
                res._all_sign = neg_;
            }
        }
        res._range->first = std::tan(f._range->first);
        res._range->second = std::tan(f._range->second);
        res._expr->_range->first = res._range->first;
        res._expr->_range->second = res._range->second;
        res._expr->_all_convexity = res._all_convexity;
        res._expr->_all_sign = res._all_sign;
        res._indices = f._indices;
        return res;
    }
    template<class T1>
    func<T1> unit_step(const func<T1>& f){
        func<T1> res(uexpr<T1>(unit_step_, f.copy()));
        if(f.is_non_positive()){
            res._range->first = zero<T1>().eval();
            res._range->second = zero<T1>().eval();
        }
        else if(f.is_positive()){
            res._range->first = unit<T1>().eval();
            res._range->second = unit<T1>().eval();
        }
        else {
            res._range->first = zero<T1>().eval();
            res._range->second = unit<T1>().eval();
        }
        res._all_convexity = undet_;
        res._expr->_range->first = res._range->first;
        res._expr->_range->second = res._range->second;
        res._expr->_all_convexity = res._all_convexity;
        res._expr->_all_sign = res._all_sign;
        res._indices = f._indices;
        return res;
    }
    
    template<class T1>
    func<T1> ReLU(const func<T1>& f){
        func<T1> res(uexpr<T1>(relu_, f.copy()));
        if (f.is_linear()) {
            res._all_convexity = convex_;
        }
        else if(!f.is_constant()){
            res._all_convexity = undet_;
        }
        res._all_sign = non_neg_;
        res._range->first = zero<T1>().eval();
        if(f.is_positive()){
            res._all_sign = pos_;
            res._range->first = f._range->first;
        }
        res._range->second = gravity::max(zero<T1>().eval(),f._range->second);
        res._expr->_range->first = res._range->first;
        res._expr->_range->second = res._range->second;
        res._expr->_all_convexity = res._all_convexity;
        res._expr->_all_sign = res._all_sign;
        res._indices = f._indices;
        return res;
    }
    
    template<class T, typename enable_if<is_arithmetic<T>::value>::type* = nullptr>
    func<T> pow(const func<T>& f, int exp){
        if(exp<0){
            return func<T>(bexpr<T>(power_, f.copy(), make_shared<constant<int>>(exp)));
        }
        if(exp==0){
            return func<T>();
        }
        if(exp==1){
            return f;
        }
        //        if(exp==2){
        //            return f*f;
        //        }
        else {
            func<T> res(f);
            for (int i = 1; i < exp; i++) {
                res *= f;
            }
            if(f._range->first==numeric_limits<T>::lowest() || f._range->second==numeric_limits<T>::max()){
                res._range->first = numeric_limits<T>::lowest();
                res._range->second = numeric_limits<T>::max();
            }
            else {
                res._range->first = gravity::min(std::pow(f._range->first,exp),std::pow(f._range->second,exp));
                res._range->second = gravity::max(std::pow(f._range->first,exp),std::pow(f._range->second,exp));
            }
            if(exp%2==0) {
                res._all_sign = non_neg_;
                if(f.is_positive()){
                    res._all_sign = pos_;
                }
                if(f._range->first <0 && f._range->second >0){
                    res._range->first = 0;
                }
            }
            else {
                res._all_sign = f.get_all_sign();
            }
            if (f.is_linear()) {
                if(exp%2==0) {
                    res._all_convexity = convex_;
                }
                else if(f.is_non_negative()){
                    res._all_convexity = convex_;
                }
                else if(f.is_non_positive()){
                    res._all_convexity = concave_;
                }
                else {
                    res._all_convexity = undet_;
                }
            }
            else if(!f.is_constant()){
                res._all_convexity = undet_;
            }
            res._indices = f._indices;
            return res;
        }
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator-(T1 p, const param<T2>& v){
        return constant<T1>(p) - v;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator-(T1 p, const param<T2>& v){
        return constant<T2>(p) - v;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator-(const param<T1>& v, T2 p){
        return v - constant<T1>(p);
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator-(const param<T1>& v, T2 p){
        return v - constant<T2>(p);
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator+(T1 p, const param<T2>& v){
        return constant<T1>(p) + v;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator+(T1 p, const param<T2>& v){
        return constant<T2>(p) + v;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator+(const param<T1>& v, T2 p){
        return v + constant<T1>(p);
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator+(const param<T1>& v, T2 p){
        return v + constant<T2>(p);
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator-(const constant<T1>& p, const param<T2>& v){
        func<T1> res(v);
        res.reverse_sign();
        res.add_cst(p);
        res._range = get_minus_range(p.range(),v._range);
        res.update_all_sign();
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator-(const constant<T1>& p, const param<T2>& v){
        func<T2> res(v);
        res.reverse_sign();
        res.add_cst(p);
        res._range = get_minus_range(p.range(),v._range);
        res.update_all_sign();
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator-(const param<T1>& v, const constant<T2>& p){
        func<T1> res(v);
        func<T1> newp = constant<T1>(p);
        newp.reverse_sign();
        res.add_cst(newp);
        res._range = get_minus_range(v._range,p.range());
        res.update_all_sign();
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator-(const param<T1>& v, const constant<T2>& p){
        func<T2> res(v);
        func<T2> newp(p);
        newp.reverse_sign();
        res.add_cst(newp);
        res._range = get_minus_range(v._range,p.range());
        res.update_all_sign();
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator+(const constant<T1>& p, const param<T2>& v){
        func<T1> res(v);
        res.add_cst(p);
        res._range = get_plus_range(v._range,p.range());
        res.update_all_sign();
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator+(const constant<T1>& p, const param<T2>& v){
        func<T2> res(v);
        res.add_cst(p);
        res._range = get_plus_range(v._range,p.range());
        res.update_all_sign();
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator+(const param<T1>& v, const constant<T2>& p){
        func<T1> res(v);
        res.add_cst(p);
        res._range = get_plus_range(v._range,p.range());
        res.update_all_sign();
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator+(const param<T1>& v, const constant<T2>& p){
        func<T2> res(v);
        res.add_cst(p);
        res._range = get_plus_range(v._range,p.range());
        res.update_all_sign();
        return res;
    }
    
    
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator+(const param<T1>& v, const func<T2>& f){
        func<T1> res(v);
        res += f;
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator+(const param<T1>& v, const func<T2>& f){
        func<T2> res(v);
        res += f;
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator+(const func<T1>& f, const param<T2>& v){
        func<T1> res(v);
        res += f;
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator+(const func<T1>& f, const param<T2>& v){
        func<T2> res(v);
        res += f;
        return res;
    }
    
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator+(const constant<T1>& v, const func<T2>& f){
        func<T1> res(v);
        res += f;
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator+(const constant<T1>& v, const func<T2>& f){
        func<T2> res(v);
        res += f;
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator+(const func<T1>& f, const constant<T2>& v){
        func<T1> res(v);
        res += f;
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator+(const func<T1>& f, const constant<T2>& v){
        func<T2> res(v);
        res += f;
        return res;
    }
    
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator+(T1 v, const func<T2>& f){
        func<T1> res(v);
        res += f;
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator+(T1 v, const func<T2>& f){
        func<T2> res(v);
        res += f;
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator+(const func<T1>& f, T2 v){
        func<T1> res(v);
        res += f;
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator+(const func<T1>& f, T2 v){
        func<T2> res(v);
        res += f;
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator-(const constant<T1>& v, const func<T2>& f){
        func<T1> res(v);
        res -= f;
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator-(const constant<T1>& v, const func<T2>& f){
        func<T2> res(v);
        res -= f;
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator-(const func<T1>& f, const constant<T2>& v){
        func<T1> res(f);
        res -= v;
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator-(const func<T1>& f, const constant<T2>& v){
        func<T2> res(f);
        res -= f;
        return res;
    }
    
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator-(T1 v, const func<T2>& f){
        func<T1> res(v);
        res -= f;
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator-(T1 v, const func<T2>& f){
        func<T2> res(v);
        res -= f;
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator-(const func<T1>& f, T2 v){
        func<T1> res(f);
        res -= v;
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator-(const func<T1>& f, T2 v){
        func<T2> res(f);
        res -= v;//TODO check when v = f
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator-(const func<T1>& f, const param<T2>& v){
        func<T1> res(v);
        res.reverse_sign();
        res += f;
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator-(const func<T1>& f, const param<T2>& v){
        func<T2> res(v);
        res.reverse_sign();
        res += f;
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator-(const param<T1>& p, const func<T2>& f){
        func<T1> res(f);
        res.reverse_sign();
        res += p;
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator-(const param<T1>& p, const func<T2>& f){
        func<T2> res(f);
        res.reverse_sign();
        res += p;
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator*(const param<T1>& v, const func<T2>& f){
        func<T1> res(v);
        func<T1> f2(f);
        if((v._is_transposed || v.is_matrix()) && (!f2.func_is_number() && !f2._is_vector)){
            return res *= f2.vec();
        }
        res *= f2;
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator*(const param<T1>& v, const func<T2>& f){
        func<T2> res(v);
        func<T2> f2(f);
        if((v._is_transposed || v.is_matrix()) && (!f2.func_is_number() && !f2._is_vector)){
            return res *= f2.vec();
        }
        res *= f2;
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator*(const func<T1>& f, const param<T2>& v){
        func<T1> res(f);
        func<T1> f2(v);
        if((f._is_transposed || f.is_matrix()) && (!f2.func_is_number() && !f2._is_vector)){
            return res *= f2.vec();
        }
        return res *= f2;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator*(const func<T1>& f, const param<T2>& v){
        func<T2> res(f);
        func<T2> f2(v);
        if((f._is_transposed || f.is_matrix()) && (!f2.func_is_number() && !f2._is_vector)){
            return res *= f2.vec();
        }
        res *= f2;
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator*(const constant<T1>& c, const param<T2>& p){
        func<T1> res;
        auto new_c(c);
        if(c._is_transposed){/* If this is a dot product resize the constant to match p's number of rows */
            new_c._dim[1] = p._dim[0];
        }
        res.update_dot_dim(new_c,p);
        res.insert(true,new_c,p);
        res._range = get_product_range(new_c.range(), p._range);
        res.update_all_sign();
        if(c._is_transposed){
            res._range->first = extended_mult(res._range->first,(T1)p._dim[0]);
            res._range->second = extended_mult(res._range->second,(T1)p._dim[0]);
        }
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator*(const constant<T1>& c, const param<T2>& p){
        func<T2> res;
        constant<T2> new_c(c);
        if(c._is_transposed){/* If this is a dot product resize the constant to match p's number of rows */
            new_c._dim[1] = p._dim[0];
        }
        res.update_dot_dim(new_c,p);
        res.insert(true,new_c,p);
        res._range = get_product_range(new_c.range(), p._range);
        if(c._is_transposed){
            res._range->first = extended_mult(res._range->first,(T2)p._dim[0]);
            res._range->second = extended_mult(res._range->second,(T2)p._dim[0]);
        }
        return res;
    }
    
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator*(const param<T1>& p, const constant<T2>& c){
        func<T1> res;
        res._range = get_product_range(p._range,c.range());
        res.update_all_sign();
        res.update_dot_dim(p,c);
        if(p._is_transposed){
            constant<T1> new_c(c.tr());
            new_c._dim[1] = p._dim[0];
            res.insert(true,new_c,p.tr());
            res._range->first = extended_mult(res._range->first,(T1)p._dim[0]);
            res._range->second = extended_mult(res._range->second,(T1)p._dim[0]);
            res.transpose();
        }
        else {
            res.insert(true,constant<T1>(c),p);
        }
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator*(const param<T1>& p, const constant<T2>& c){
        func<T2> res;
        res._range = get_product_range(p._range,c.range());
        res.update_all_sign();
        res.update_dot_dim(p,c);
        if(p._is_transposed){
            constant<T2> new_c(c.tr());
            new_c._dim[1] = p._dim[0];
            res.insert(true,new_c,p.tr());
            res._range->first = extended_mult(res._range->first,(T2)p._dim[0]);
            res._range->second = extended_mult(res._range->second,(T2)p._dim[0]);
            res.transpose();
        }
        else {
            res.insert(true,c,p);
        }
        return res;
    }
    
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator*(const constant<T1>& c, const func<T2>& p){
        return func<T1>(c) *= p;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator*(const constant<T1>& c, const func<T2>& p){
        return func<T2>(c) *= p;
    }
    
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator*(const func<T1>& p, const constant<T2>& c){
        return func<T1>(p) *= func<T1>(c);
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator*(const func<T1>& p, const constant<T2>& c){
        return func<T2>(p) *= func<T2>(c);
    }
    
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator*(T1 p, const param<T2>& v){
        return constant<T1>(p)*v;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator*(T1 p, const param<T2>& v){
        return constant<T2>(p)*v;
    }
    
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator*(const param<T1>& v, T2 p){
        return v*constant<T1>(p);
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator*(const param<T1>& v, T2 p){
        return v*constant<T2>(p);
    }
    
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator*(T1 p, const func<T2>& f){
        return constant<T1>(p) * f;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator*(T1 p, const func<T2>& f){
        return constant<T2>(p)*f;
    }
    
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator*(const func<T1>& f, T2 p){
        return f * constant<T1>(p);
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator*(const func<T1>& f, T2 p){
        return f * constant<T2>(p);
    }
    
    
    
    template<typename type>
    func<type> sum(const param<type>& p){
        func<type> res;
        if (p.get_dim()==0) {
            return res;
        }
        if(p.is_matrix_indexed()){
            return (unit<type>().tr()*p.vec()).in(range(0,p._indices->size()-1));
        }
        return unit<type>().tr()*p.vec();
    }
    
    /** Create a matrix sum where each row will be indexed based on the entries starting at start_pos and spaning nb_entries.
     Example:
     dv = {
     [1,8] = 0
     [1,9] = 0
     [1,10] = 0
     [1,11] = 0
     [1,12] = 0
     [2,8] = 0
     [2,9] = 0
     [2,10] = 0
     [2,11] = 0
     [2,12] = 0
     [3,8] = 0
     [3,9] = 0
     [3,10] = 0
     [3,11] = 0
     [3,12] = 0
     };
     sum(dv.in_matrix(0,1)) <= 0 gives: dv[1,8] + dv[2,8] + dv[3,8] <= 0;
     sum(dv.in_matrix(1,1)) <= 0 gives: dv[1,8] + dv[1,9] + dv[1,10] + dv[1,11] + dv[1,12] <= 0;
     */
    template<typename type>
    func<type> sum_ith(const param<type>& p, unsigned start_pos, unsigned nb_entries){
        auto matrix_p = p.in_matrix(start_pos,nb_entries);
        auto res = sum(matrix_p);
        return res.in(range(0,matrix_p._indices->get_nb_rows()-1));
    }
    
    /** Create a matrix sum where each row will be indexed based on the entries starting at start_pos and spaning nb_entries.
     Example:
     dv = {
     [1,8] = 0
     [1,9] = 0
     [1,10] = 0
     [1,11] = 0
     [1,12] = 0
     [2,8] = 0
     [2,9] = 0
     [2,10] = 0
     [2,11] = 0
     [2,12] = 0
     [3,8] = 0
     [3,9] = 0
     [3,10] = 0
     [3,11] = 0
     [3,12] = 0
     };
     sum(dv.in_matrix(0,1)) <= 0 gives: dv[1,8] + dv[2,8] + dv[3,8] <= 0;
     sum(dv.in_matrix(1,1)) <= 0 gives: dv[1,8] + dv[1,9] + dv[1,10] + dv[1,11] + dv[1,12] <= 0;
     */
    template<typename type>
    func<type> sum_ith(const var<type>& p, unsigned start_pos, unsigned nb_entries){
        auto matrix_p = p.in_matrix(start_pos,nb_entries);
        auto res = sum(matrix_p);
        return res.in(range(0,matrix_p._indices->get_nb_rows()-1));
    }
    
    template<typename type>
    func<type> sum(const var<type>& p){
        func<type> res;
        if (p.get_dim()==0) {
            return res;
        }
        if(p.is_matrix_indexed()){
            return (unit<type>().tr()*p.vec()).in(range(0,p._indices->size()-1));
        }
        return unit<type>().tr()*p.vec();
    }
    
    template<typename type>
    func<type> sum(const func<type>& p){
        func<type> res;
        if (p.get_dim()==0) {
            return res;
        }
        if(p.is_matrix_indexed()){
            return (unit<type>().tr()*p.vec()).in(range(0,p._indices->size()-1));
        }
        return unit<type>().tr()*p.vec();
    }
    
    template<typename type>
    func<type> sum(const var<type>& p, const indices& ids){
        func<type> res;
        if (p.get_dim()==0) {
            return res;
        }
        if(p.is_matrix_indexed()){
            return (unit<type>().tr()*(p.vec()).in(ids)).in(ids);
        }
        return unit<type>().tr()*(p.vec()).in(ids);
    }
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> product(const func<T1>& f1, const func<T2>& f2){
        if(f1.is_column_vector() && f2.is_column_vector()){/* This is a dot product */
            return f1.tr()*f2.vec();
        }
        if(f1.is_matrix() && f2.is_column_vector())
            return f1*f2.vec();
        return f1*f2;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> product(const func<T1>& f1, const func<T2>& f2){
        if(f1.is_column_vector() && f2.is_column_vector()){/* This is a dot product */
            return f1.tr()*f2.vec();
        }
        if(f1.is_matrix() && f2.is_column_vector())
            return f1*f2.vec();
        return f1*f2;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> product(const param<T1>& f1, const param<T2>& f2){
        if(f1.is_column_vector() && f2.is_column_vector()){/* This is a dot product */
            return f1.tr()*f2.vec();
        }
        if(f1.is_matrix() && f2.is_column_vector())
            return f1*f2.vec();
        return f1*f2;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> product(const var<T1>& f1, const var<T2>& f2){
        if(f1.is_column_vector() && f2.is_column_vector()){/* This is a dot product */
            return f1.tr()*f2.vec();
        }
        if(f1.is_matrix() && f2.is_column_vector())
            return f1*f2.vec();
        return f1*f2;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> product(const var<T1>& f1, const var<T2>& f2){
        if(f1.is_column_vector() && f2.is_column_vector()){/* This is a dot product */
            return f1.tr()*f2.vec();
        }
        if(f1.is_matrix() && f2.is_column_vector())
            return f1*f2.vec();
        return f1*f2;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> product(const var<T1>& f1, const param<T2>& f2){
        if(f1.is_column_vector() && f2.is_column_vector()){/* This is a dot product */
            return f1.tr()*f2.vec();
        }
        if(f1.is_matrix() && f2.is_column_vector())
            return f1*f2.vec();
        return f1*f2;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> product(const var<T1>& f1, const param<T2>& f2){
        if(f1.is_column_vector() && f2.is_column_vector()){/* This is a dot product */
            return f1.tr()*f2.vec();
        }
        if(f1.is_matrix() && f2.is_column_vector())
            return f1*f2.vec();
        return f1*f2;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> product(const param<T1>& f1, const var<T2>& f2){
        if(f1.is_column_vector() && f2.is_column_vector()){/* This is a dot product */
            return f1.tr()*f2.vec();
        }
        if(f1.is_matrix() && f2.is_column_vector())
            return f1*f2.vec();
        return f1*f2;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> product(const param<T1>& f1, const var<T2>& f2){
        if(f1.is_column_vector() && f2.is_column_vector()){/* This is a dot product */
            return f1.tr()*f2.vec();
        }
        if(f1.is_matrix() && f2.is_column_vector())
            return f1*f2.vec();
        return f1*f2;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> product(const param<T1>& f1, const param<T2>& f2){
        if(f1.is_column_vector() && f2.is_column_vector()){/* This is a dot product */
            return f1.tr()*f2.vec();
        }
        if(f1.is_matrix() && f2.is_column_vector())
            return f1*f2.vec();
        return f1*f2;
    }
    
    
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> product(const param<T1>& f1, const func<T2>& f2){
        if(f1.is_column_vector() && f2.is_column_vector()){/* This is a dot product */
            return f1.tr()*f2.vec();
        }
        if(f1.is_matrix() && f2.is_column_vector())
            return f1*f2.vec();
        return f1*f2;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> product(const param<T1>& f1, const func<T2>& f2){
        if(f1.is_column_vector() && f2.is_column_vector()){/* This is a dot product */
            return f1.tr()*f2.vec();
        }
        if(f1.is_matrix() && f2.is_column_vector())
            return f1*f2.vec();
        return f1*f2;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> product(const func<T1>& f1, const param<T2>& f2){
        if(f1.is_column_vector() && f2.is_column_vector()){/* This is a dot product */
            return f1.tr()*f2.vec();
        }
        if(f1.is_matrix() && f2.is_column_vector())
            return f1*f2.vec();
        return f1*f2;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> product(const func<T1>& f1, const param<T2>& f2){
        if(f1.is_column_vector() && f2.is_column_vector()){/* This is a dot product */
            return f1.tr()*f2.vec();
        }
        if(f1.is_matrix() && f2.is_column_vector())
            return f1*f2.vec();
        return f1*f2;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> product(const var<T1>& f1, const func<T2>& f2){
        if(f1.is_column_vector() && f2.is_column_vector()){/* This is a dot product */
            return f1.tr()*f2.vec();
        }
        if(f1.is_matrix() && f2.is_column_vector())
            return f1*f2.vec();
        return f1*f2;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> product(const var<T1>& f1, const func<T2>& f2){
        if(f1.is_column_vector() && f2.is_column_vector()){/* This is a dot product */
            return f1.tr()*f2.vec();
        }
        if(f1.is_matrix() && f2.is_column_vector())
            return f1*f2.vec();
        return f1*f2;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> product(const func<T1>& f1, const var<T2>& f2){
        if(f1.is_column_vector() && f2.is_column_vector()){/* This is a dot product */
            return f1.tr()*f2.vec();
        }
        if(f1.is_matrix() && f2.is_column_vector())
            return f1*f2.vec();
        return f1*f2;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> product(const func<T1>& f1, const var<T2>& f2){
        if(f1.is_column_vector() && f2.is_column_vector()){/* This is a dot product */
            return f1.tr()*f2.vec();
        }
        if(f1.is_matrix() && f2.is_column_vector())
            return f1*f2.vec();
        return f1*f2;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> product(const param<T1>& f1, const constant<T2>& f2){
        constant<T1> new_f2(f2);
        if(f1.is_column_vector()){/* This is a dot product */
            new_f2._dim[0] = f1._dim[1];
            return f1.tr()*new_f2.vec();
        }
        return f1*new_f2;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> product(const param<T1>& f1, const constant<T2>& f2){
        constant<T2> new_f2(f2);
        if(f1.is_column_vector()){/* This is a dot product */
            new_f2._dim[0] = f1._dim[1];
            return f1.tr()*new_f2.vec();
        }
        return f1*new_f2;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> product(const constant<T1>& f1, const param<T2>& f2){
        constant<T1> new_f1(f1);
        if(f2.is_column_vector()){/* Dot product with a constant */
            new_f1.transpose();
            new_f1._dim[1] = f2._dim[0];
            return new_f1*f2.vec();
        }
        return new_f1*f2;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> product(const constant<T1>& f1, const param<T2>& f2){
        constant<T2> new_f1(f1);
        if(f2.is_column_vector()){/* Dot product with a constant */
            new_f1.transpose();
            new_f1._dim[1] = f2._dim[0];
            return new_f1*f2.vec();
        }
        return new_f1*f2;
    }
    
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> product(const var<T1>& f1, const constant<T2>& f2){
        constant<T1> new_f2(f2);
        if(f1.is_column_vector()){/* This is a dot product */
            new_f2._dim[0] = f1._dim[1];
            return f1.tr()*new_f2.vec();
        }
        return f1*new_f2;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> product(const var<T1>& f1, const constant<T2>& f2){
        constant<T2> new_f2(f2);
        if(f1.is_column_vector()){/* This is a dot product */
            new_f2._dim[0] = f1._dim[1];
            return f1.tr()*new_f2.vec();
        }
        return f1*new_f2;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> product(const constant<T1>& f1, const var<T2>& f2){
        constant<T1> new_f1(f1);
        if(f2.is_column_vector()){/* Dot product with a constant */
            new_f1.transpose();
            new_f1._dim[1] = f2._dim[0];
            return new_f1*f2.vec();
        }
        return new_f1*f2;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> product(const constant<T1>& f1, const var<T2>& f2){
        constant<T2> new_f1(f1);
        if(f2.is_column_vector()){/* Dot product with a constant */
            new_f1.transpose();
            new_f1._dim[1] = f2._dim[0];
            return new_f1*f2.vec();
        }
        return new_f1*f2;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> product(const constant<T1>& f1, const func<T2>& f2){
        constant<T1> new_f1(f1);
        if(f2.is_column_vector()){/* Dot product with a constant */
            new_f1.transpose();
            new_f1._dim[1] = f2._dim[0];
            return new_f1*f2.vec();
        }
        return new_f1*f2;
    }
    
    
    template<class T1=double,class T2=double, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> product(const constant<T1>& f1, const func<T2>& f2){
        constant<T2> new_f1(f1);
        if(f2.is_column_vector()){/* Dot product with a constant */
            new_f1.transpose();
            new_f1._dim[1] = f2._dim[0];
            return new_f1*f2.vec();
        }
        return new_f1*f2;
    }
    
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> product(const func<T1>& f1, const constant<T2>& f2){
        constant<T1> new_f2(f2);
        if(f1.is_column_vector()){/* This is a dot product */
            new_f2._dim[0] = f1._dim[1];
            return f1.tr()*new_f2.vec();
        }
        return f1*new_f2;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> product(const func<T1>& f1, const constant<T2>& f2){
        constant<T2> new_f2(f2);
        if(f1.is_column_vector()){/* This is a dot product */
            new_f2._dim[0] = f1._dim[1];
            return f1.tr()*new_f2.vec();
        }
        return f1*new_f2;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> product(T1 f1, const param<T2>& f2){
        return product(constant<T1>(f1), f2);
    }
    
    template<class T1=double,class T2=double, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> product(T1 f1, const param<T2>& f2){
        return product(constant<T2>(f1), f2);
    }
    
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> product(const param<T1>& f1, T2 f2){
        return product(f1, constant<T1>(f2));
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> product(const param<T1>& f1, T2 f2){
        return product(f1, constant<T2>(f2));
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> product(T1 f1, const var<T2>& f2){
        return product(constant<T1>(f1), f2);
    }
    
    template<class T1=double,class T2=double, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> product(T1 f1, const var<T2>& f2){
        return product(constant<T2>(f1), f2);
    }
    
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> product(const var<T1>& f1, T2 f2){
        return product(f1, constant<T1>(f2));
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> product(const var<T1>& f1, T2 f2){
        return product(f1, constant<T2>(f2));
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> product(T1 f1, const func<T2>& f2){
        return product(constant<T1>(f1), f2);
    }
    
    template<class T1=double,class T2=double, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> product(T1 f1, const func<T2>& f2){
        return product(constant<T2>(f1), f2);
    }
    
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> product(const func<T1>& f1, T2 f2){
        return product(f1, constant<T1>(f2));
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> product(const func<T1>& f1, T2 f2){
        return product(f1, constant<T2>(f2));
    }
    
    
    //    template<typename type>
    //    func<type> sum(const func<type>& f, const indices& ids){
    //        auto ff = f;
    //        ff.index_in(ids);
    //        func<type> res;
    //        return unit<type>().tr()*(ff.vec()).in(ids);
    //    }
    
    /** WARNING, only call if the variables appearing in the function are complex or double */
    pair<func<double>,func<double>> get_real_imag(const func<Cpx>& f);
    pair<func<double>,func<double>> get_mag_ang(const func<Cpx>& f);
    func<double> get_real(constant_* c);
    func<double> get_imag(constant_* c);
    func<double> get_mag(constant_* c);
    func<double> get_ang(constant_* c);
    }
    
    
    
#endif /* func_h */
    
