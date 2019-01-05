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
    class func_ : virtual public constant_{
    private:
        shared_ptr<func_> compute_derivative(const param_& v);  /**< Computes and stores the derivative of f with respect to variable v. Returns a pointer to the stored function. */
        
    protected:
        FType                                                             _ftype = const_; /**< Function type, e.g., constant, linear, quadratic... >>**/
        NType                                                             _return_type = double_; /**< Return type, e.g., bool, integer, complex... >>**/

        shared_ptr<map<string, pair<shared_ptr<param_>, unsigned>>>       _params = nullptr;/**< Set of parameters in current function, stored as a map <parameter name, <paramter pointer, number of times it appears in function>>**/
        shared_ptr<map<string, pair<shared_ptr<param_>, unsigned>>>       _vars = nullptr;/**< Set of variables in current function, stored as a map <variable name, <variable pointer, number of times it appears in function>>**/
        
        shared_ptr<constant_>                                             _cst = nullptr;/**< Constant part of the function */
        shared_ptr<map<string, lterm>>                                    _lterms = nullptr; /**< Set of linear terms, stored as a map <string describing term, term>. */
        shared_ptr<map<string, qterm>>                                    _qterms = nullptr; /**< Set of quadratic terms, stored as a map <string describing term, term>.  */
        shared_ptr<map<string, pterm>>                                    _pterms = nullptr; /**< Set of polynomial terms, stored as a map <string describing term, term>.  */
        shared_ptr<expr>                                                  _expr = nullptr; /**< Nonlinear part of the function, this points to the root node in _DAG */
//        map<string, expr*>*                    _DAG = nullptr; /**< Map of experssions stored in the expression tree (a Directed Acyclic Graph) */
//        deque<shared_ptr<expr>>*               _queue = nullptr; /**< A queue storing the expression tree from the leaves to the root (the root is stored at the end of the queue)*/
        Convexity                                                         _all_convexity = linear_; /**< If all instances of this function have the same convexity type, it stores it here, i.e. linear, convex, concave, otherwise it stores unknown. >>**/
        Sign                                                              _all_sign = zero_; /**< If all instances of this function have the same sign, it stores it here, otherwise it stores unknown. >>**/

        shared_ptr<vector<Convexity>>                                     _convexity = nullptr; /**< Vector of convexity types, i.e., linear, convex, concave or unknown. This is a vector since a function can have multiple instances (different constants coefficients, and bounds, but same structure) >>**/
        shared_ptr<vector<Sign>>                                          _sign = nullptr; /**< vector storing the sign of return value if known. >>**/
        shared_ptr<map<size_t, set<size_t>>>                              _hess_link = nullptr; /**< Set of variables linked to one another in the hessian, stored by variable ids  */
        shared_ptr<map<string,shared_ptr<func_>>>                         _dfdx = nullptr;/**< A map storing the derivatives indexed by variables' names */

        bool                                                              _is_constraint = false;
        bool                                                              _is_hessian = false;
        bool                                                              _embedded = false; /**< If the function is embedded in a mathematical model or in another function, this is used for memory management. >>**/
        bool                                                              _evaluated = true;/**< If the function has already been evaluated, useful for constant funcs */
        string                                                            _to_str = "noname";/**< A string representation of the expression */

        size_t                                                            _nb_vars = 0; /**< Number of variables */
        
        size_t                                                            _nnz_j = 0; /**< Number of nonzeros in the Jacobian **/
        size_t                                                            _nnz_h = 0; /**< Number of nonzeros in the Hessian **/
        
    public:

        /** Accessors */
        FType get_ftype() const;
        NType get_return_type() const;

        map<size_t, set<size_t>>& get_hess_link() { return *_hess_link;};
        map<string, pair<shared_ptr<param_>, unsigned>>& get_vars() { return *_vars;};
        map<string, pair<shared_ptr<param_>, unsigned>>& get_params() { return *_params;};

        /** true/false statements */
        bool has_var(const param_& v) const;
        bool has_var(const string& name) const;
        bool is_convex() const;
        bool is_concave() const;
        bool is_convex(size_t idx) const;
        bool is_concave(size_t idx) const;
        bool is_constant() const;
        bool is_linear() const;
        bool is_quadratic() const;
        bool is_polynomial() const;
        bool is_nonlinear() const;
        bool is_complex() const;
        bool is_transposed() const;
        bool is_number() const{
            return (_vars->empty() && _params->empty());
        }
        
        
        /* Virtual functions */
        
        virtual shared_ptr<func_> fcopy() const{return nullptr;};
        
        /* Modifiers */
        
        /**
         Mark f as embeded and merge variables and parameters with f (by calling merge_vars(func_&& f). If a variable x in f exists in the current funtion, x will now point to the same variable appearing in current function.
         @param[in] f function to merge variables and parameters with.
         */
        void embed(func_& f);
        /**
         Merge variables and parameters with expression e. If a variable x in e exists in the current funtion, x will now point to the same variable appearing in the current function.
         @param[in] e expression to merge variables and parameters with.
         */
        void embed(shared_ptr<expr> e);
        
        /**
         Subfuntion of embed(func_&& f). Merge variables and parameters with f. If a variable x in f exists in the current funtion, x will now point to the same variable appearing in current function.
         @param[in] f function to merge variables and parameters with.
         */
        void merge_vars(func_& f);
        
        /**
         Copy and embed derivatives of f.
         @param[in] f function to copy derivatives from.
         */
        void copy_derivatives(const func_& f);

        void set_first_derivative(const param_& v, func_&& f){
            DebugOff(f.to_str()<<endl);
            (*_dfdx)[v._name] = make_shared<func_>(move(f));
        }

        void set_second_derivative(const param_& v1, const param_& v2, func_&& f){
            DebugOff(f.to_str()<<endl);
            (*_dfdx)[v1._name]->_dfdx->insert(make_pair<>(v2._name, make_shared<func_>(move(f))));
        }

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

        shared_ptr<expr> get_expr() const{
            return _expr;
        }

        shared_ptr<map<string,shared_ptr<func_>>> get_dfdx() const{
            return _dfdx;
        };

        shared_ptr<func_> get_stored_derivative(const string& vid) const; /**< Returns the stored derivative with respect to variable v. */

        func_ get_derivative(const param_& v) const; /**< Computes and returns the derivative with respect to variable v. */

        func_ get_dfdx(const param_& v); /**< Computes all derivatives and returns a copy of the derivative with respect to variable v. */


        void compute_derivatives(); /**< Computes and stores the derivative of f with respect to all variables. */

        
        /**
         Returns a vector of monomials of degree d using the variables in the current function
         @param[in] d degree of monomials
         @return a vector of monomials of degree d using the variables in the current function
         */
        vector<pterm> get_monomials(unsigned d);
        
        
        qterm* get_square(shared_ptr<param_> p); /**< Returns the quadratic term containing a square of p or nullptr if none exists. **/
        
        /**
         Returns the convexity of current function if quadratic term q was to be added.
         @param[in] q quadratic term to be added.
         @return convexity of function if q was to be added.
         */
        Convexity get_convexity(const qterm& q);
        
        /**
         Index the function and its variables/parameters using nodes of a graph
         @param[in] vec vector of nodes
         @return current function
         */
        func_& in(const vector<Node*>& vec);
        
        /**
         Index the function and its variables/parameters using the indices in ids
         @param[in] ids indices
         @return current function
         */
        func_& in(const indices& ids);
        
        /**
         Relax and replace integer variables with continuous ones provided in argument vars.
         @param[in] vars set with continuous variables replacements.
         */
        void relax(const map<size_t, shared_ptr<param_>>& vars);
        
        /**
         Returns the number of variables per-instance.
         @param[in] instance number.
         @return number of variables per-instance.
         */
        size_t get_nb_vars(unsigned inst = 0) const{
            size_t n = 0;
            for (auto &vp:*_vars) {
                if(vp.second.first->_is_vector){
                    n += vp.second.first->get_dim(inst);
                }
                else {
                    n += 1;
                }
            }
            return n;
        };
        
        void update_vars();
        
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
        shared_ptr<param_> get_var(string name) const;

        /**
         Returns a pointer to the parameter matching the name provided.
         @param[in] name variable name.
         @return a pointer to the parameter matching the name provided.
         */
        shared_ptr<param_> get_param(string name) const;

        void add_var(shared_ptr<param_> v, int nb = 1);/**< Inserts the variable in this function input list. nb represents the number of occurences v has. WARNING: Assumes that v has not been added previousely!*/

        void add_param(shared_ptr<param_> v, int nb = 1);/**< Inserts the parameter in this function input list. nb represents the number of occurences v has. WARNING: Assumes that v has not been added previousely!*/

        /**
         Reverse the sign of all terms in the function
         */
        void reverse_sign();

        /**
         Reverse the convexity property of the current function
         */
        void reverse_convexity();

    };

//
//        bool insert(bool sign, const constant_& coef, const param_& p);/**< Adds coef*p to the function. Returns true if added new term, false if only updated coef of p */
//        bool insert(bool sign, const constant_& coef, const param_& p1, const param_& p2, bool c_p1_transposed=false);/**< Adds coef*p1*p2 to the function. Returns true if added new term, false if only updated coef of p1*p2 */
//        bool insert(bool sign, const constant_& coef, const list<pair<param_*, int>>& l);/**< Adds polynomial term to the function. Returns true if added new term, false if only updated corresponding coef */
//       
//        
//        void transpose(){
//            this->constant_::transpose();
//            if(_expr){
//                _expr->transpose();
//            }
//            if(_vars->size()==1 && _params->size()==0){ // If function is a variable.
//                auto vars_cpy = *_vars;
//                for (auto &vp:*_vars) {
//                    vp.second.first->transpose();
//                    vars_cpy.erase(vp.first);
//                    vars_cpy[vp.second.first->get_name(false,false)]= make_pair<>(vp.second.first, vp.second.second);
//                }
//                *_vars = move(vars_cpy);
//            }
//            else if(_vars->size()==0 && _params->size()==1){ // If function is a parameter.
//                auto params_cpy = *_params;
//                for (auto &vp:*_params) {
//                    vp.second.first->transpose();
//                    params_cpy.erase(vp.first);
//                    params_cpy[vp.second.first->get_name(false,false)]= make_pair<>(vp.second.first, vp.second.second);
//                }
//                *_params = move(params_cpy);
//
//            }
//        }

//        
//

//        
//
//        void insert(const lterm& term);
//        
//        void insert(const qterm& term);
//        
//        void insert(const pterm& term);
//        
//        void update_to_str(bool input = false);

//        
//        
//        size_t get_nb_instances() const {
//            return _dim[0];
//        }
////            return max((size_t)1,constant_::get_nb_instances());
////        }
//        
//
//        
//        void add_param(shared_ptr<param_> p, int nb = 1);/**< Inserts the parameter in this function input list. WARNING: Assumes that p has not been added previousely!*/
//        
//        
//        
//        void delete_var(const string& vid);
//        
//        void delete_param(const string& vid);
//        

//        void replace(param_* v, func_& f);/**<  Replace v with function f everywhere it appears */
//
//        
//        void reset();
//        
//
////        void untranspose_derivatives(){
////            for (auto &fp:*_dfdx) {
////                auto df = fp.second;
////                df->transpose();
////                df->_nb_instances = max(df->_nb_instances, _nb_instances);
//////                df->_val->resize(max(df->_val->size(),_nb_instances));
////                df->untranspose_derivatives();
////            }
////        }
//        
//        
//
//
//        func_& operator=(func_&& f);
//        
//        bool operator==(const func_& f) const;
//        
//        bool operator!=(const func_& f) const;
//        
//
//        func_& operator+=(const constant_& f);
//        func_& operator-=(const constant_& f);
//        func_& operator*=(const constant_& f);
//        func_& operator/=(const constant_& f);
//        
//        
//        friend func_ cos(const constant_& c);
//        friend func_ cos(constant_&& c);
//        
//        friend func_ sin(const constant_& c);
//        friend func_ sin(constant_&& c);
//        
//        
//        friend func_ sqrt(const constant_& c);
//        friend func_ sqrt(constant_&& c);
//        
//        friend func_ expo(const constant_& c);
//        friend func_ expo(constant_&& c);
//        
//        friend func_ log(const constant_& c);
//        friend func_ log(constant_&& c);
//
//        
//        template<class T, class = typename enable_if<is_arithmetic<T>::value>::type> func_& operator+=(T c){
//            return *this += constant<T>(c);
//        };
//        
//        template<class T, class = typename enable_if<std::is_arithmetic<T>::value>::type> func_& operator-=(T c){
//            return *this -= constant<T>(c);
//        };
//        
//        template<class T, class = typename enable_if<std::is_arithmetic<T>::value>::type> func_& operator*=(T c){
//            return *this *= constant<T>(c);
//        };
//        
//        template<class T, class = typename enable_if<std::is_arithmetic<T>::value>::type> func_& operator/=(T c){
//            return *this /= constant<T>(c);
//        };
//        
//
//      
//        func_ get_outer_app(); /**< Returns an outer-approximation of the function using the current value of the variables **/
//        
//        Sign get_all_sign() const;
//        pair<double, double>* get_all_range() const;
//        Sign get_sign(size_t idx=0) const;
//        Sign get_all_sign(const lterm& l);
//        Sign get_all_sign(const qterm& l);
//        Sign get_all_sign(const pterm& l);
//
//        
//        
//        void update_sign(const constant_& c);
//        
//        void update_sign(const lterm& l);
//        
//        void update_sign(const qterm& q);
//        
//        void update_sign(const pterm& q);
//        
//        void update_convexity(const qterm& q);
//        
//        void update_convexity();
//        
//        void update_dim(const lterm& l);
//        
//        void update_dim(const qterm& q);
//        
//        bool is_soc();
//        bool is_rotated_soc();
//        
//
//        void update_sign();
//        
//        string to_str() const;
//        string to_str(size_t inst) const;
//        void print_symbolic(bool endline = true, bool display_input = true);
//        void print(size_t index);
//        void print();
//    };
//
    template<typename type = double>
    class func: public func_, public var<type>{
    public:
//        pair<type, type>*                      _all_range = nullptr; /**< Range of the return value considering all instances of the current function. >>**/
//        vector<pair<type, type>>*              _range = nullptr; /**< Bounds of the return value per-instance. >>**/
//        
//
//
        func(){
            constant_::set_type(func_c);
            _lterms = make_shared<map<string, lterm>>();
            _qterms = make_shared<map<string, qterm>>();
            _pterms = make_shared<map<string, pterm>>();
            _vars = make_shared<map<string, pair<shared_ptr<param_>, unsigned>>>();
            _params = make_shared<map<string, pair<shared_ptr<param_>, unsigned>>>();
        };
//            constant_::set_type(func_c);
//            _params = new map<string, pair<shared_ptr<param_>, unsigned>>();
//            _vars = new map<string, pair<shared_ptr<param_>, unsigned>>();
//            _cst = new constant<type>();
//            _lterms = new map<string, lterm>();
//            _qterms = new map<string, qterm>();
//            _pterms = new map<string, pterm>();
//            _dfdx = make_shared<map<string,shared_ptr<func_>>>();
//            _DAG = new map<string, expr*>();
//            _queue = new deque<shared_ptr<expr>>();
//            _all_range = new pair<type,type>(numeric_limits<type>::lowest(),numeric_limits<type>::max());
//            _sign = nullptr;
//            _convexity = nullptr;
//            _range = nullptr;
//            _evaluated = true;
//        };
//        
//        func(const type& c){
//            *this = constant<type>(c);
//        };
//
        Sign get_all_sign() const{ /**< If all instances of the current parameter/variable have the same sign, it returns it, otherwise, it returns unknown. **/
            return param<type>::get_all_sign();
        };
        Sign get_sign(size_t idx = 0) const{ /**< returns the sign of one instance of the current parameter/variable. **/
            return param<type>::get_sign(idx);
        }
        
        template<typename... Args>
        func in(const indices& vec1, Args&&... args) {
            func<type> res(*this);
            res.param<type>::operator=(param<type>::in(vec1, forward<Args>(args)...));
            return res;
        }
        
        string to_str() const {
            return to_str(0, 3);
        }
        
        string to_str(size_t index, int prec) const {
            return param<type>::to_str(index, prec);
        }
        
        string to_str(size_t index1, size_t index2, int prec) const {
            return param<type>::to_str(index1, index2, prec);
        }
        

        void propagate_dim(size_t d){
            if (param<type>::is_matrix()) {
                return;
            }
            if(param<type>::_is_transposed){
                param<type>::_dim[1] = d;
            }
            else {
                param<type>::_dim[0] = d;
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
            var<type>::_val->resize(var<type>::get_dim());
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
//
//        func(constant_&& c):func(){
//            if(c.is_function()){
//                auto f = (func_*)&c;
//                switch (f->get_type()) {
//                    case binary_:
//                        *this = move(*(func<bool>*)f);
//                        break;
//                    case short_:
//                        *this = move(*(func<short>*)f);
//                        break;
//                    case integer_:
//                        *this = move(*(func<int>*)f);
//                        break;
//                    case float_:
//                        *this = move(*(func<float>*)f);
//                        break;
//                    case double_:
//                        *this = move(*(func<double>*)f);
//                        break;
//                    case long_:
//                        *this = move(*(func<long double>*)f);
//                        break;
//                    case complex_:
//                        *this = move(*(func<Cpx>*)f);
//                        break;
//                    default:
//                        throw invalid_argument("unsupported numerical type");
//                        break;
//                }
//            }
//            else {
//                _dim[0] = c._dim[0];
//                _dim[1] = c._dim[1];
//                switch (c.get_type()) {
//                    case binary_c: {
//                        _cst = new constant<bool>(move(*(constant<bool>*)(&c)));
//                        _all_sign = ((constant<bool>*)_cst)->get_sign();
//                        auto val = ((constant<bool>*)_cst)->eval();
//                        _all_range = new pair<bool,bool>(val,val);
//                        var<type>::_val = make_shared<vector<bool>>();
//                        var<type>::_val->push_back(val);
//                        _evaluated = true;
//                        break;
//                    }
//                    case short_c: {
//                        _cst = new constant<short>(move(*(constant<short>*)(&c)));
//                        _all_sign = ((constant<short>*)_cst)->get_sign();
//                        auto val = ((constant<short>*)(&c))->eval();
//                        _all_range = new pair<short,short>(val,val);
//                        var<type>::_val = make_shared<vector<short>>();
//                        var<type>::_val->push_back(val);
//                        _evaluated = true;
//                        break;
//                    }
//                    case integer_c: {
//                        _cst = new constant<int>(move(*(constant<int>*)(&c)));
//                        _all_sign = ((constant<int>*)_cst)->get_sign();
//                        auto val = ((constant<int>*)(&c))->eval();
//                        _all_range = new pair<int,int>(val,val);
//                        var<type>::_val = make_shared<vector<int>>();
//                        var<type>::_val->push_back(((constant<int>*)(&c))->eval());
//                        _evaluated = true;
//                        break;
//                    }
//                    case float_c: {
//                        _cst = new constant<float>(move(*(constant<float>*)(&c)));
//                        _all_sign = ((constant<float>*)_cst)->get_sign();
//                        auto val = ((constant<float>*)(&c))->eval();
//                        _all_range = new pair<float,float>(val,val);
//                        var<type>::_val = make_shared<vector<float>>();
//                        var<type>::_val->push_back(val);
//                        _evaluated = true;
//                        break;
//                    }
//                    case double_c: {
//                        _cst = new constant<double>(move(*(constant<double>*)(&c)));
//                        _all_sign = ((constant<double>*)_cst)->get_sign();
//                        auto val = ((constant<double>*)(&c))->eval();
//                        _all_range = new pair<double,double>(val,val);
//                        var<type>::_val = make_shared<vector<double>>();
//                        var<type>::_val->push_back(val);
//                        _evaluated = true;
//                        break;
//                    }
//                    case long_c: {
//                        _cst = new constant<long double>(move(*(constant<long double>*)(&c)));
//                        _all_sign = ((constant<long double>*)_cst)->get_sign();
//                        auto val = ((constant<long double>*)(&c))->eval();
//                        _all_range = new pair<long double,long double>(val,val);
//                        var<type>::_val = make_shared<vector<long double>>();
//                        var<type>::_val->push_back(((constant<long double>*)(&c))->eval());
//                        _evaluated = true;
//                        break;
//                    }
//                    case complex_c: {
//                        _cst = new constant<Cpx>(*(constant<Cpx>*)(&c));
//                        _all_sign = g
//                        _all_range = new pair<double,double>(numeric_limits<double>::lowest(),numeric_limits<double>::max());
//                        auto val = ((constant<Cpx>*)(&c))->eval();
//                        var<type>::_val = make_shared<vector<Cpx>>();
//                        var<type>::_val->push_back(val);
//                        break;
//                    }
//                    case par_c:{
//                        auto p_c2 =     shared_ptr<param_>((param_*)copy(move(c)));
//                        //                    p_c2->untranspose();
//                        _lterms->insert(make_pair<>(p_c2->get_name(false,false), p_c2.get()));
//                        add_param(p_c2);
//                        _cst = new constant<double>(0);
//                        _all_sign = p_c2->get_all_sign();
//                        _all_range = p_c2->get_range();
//                        _is_transposed = c._is_transposed;
//                        _is_vector = c._is_vector;
//                        //                    _is_matrix = c._is_matrix;
//                        _dim[0] = c._dim[0];
//                        _dim[1] = c._dim[1];
//                        //                    _nb_instances = c._dim[0];
//                        //                    _indices = p_c2->_indices;
//                        _evaluated = false;
//                        break;
//                    }
//                    case var_c:{
//                        auto p_c2 = shared_ptr<param_>((param_*)copy(move(c)));
//                        _lterms->insert(make_pair<>(p_c2->get_name(false,false), p_c2.get()));
//                        add_var(p_c2);
//                        _ftype = lin_;
//                        _is_transposed = c._is_transposed;
//                        _is_vector = c._is_vector;
//                        //                    _is_matrix = c._is_matrix;
//                        _dim[0] = c._dim[0];
//                        _dim[1] = c._dim[1];
//                        //                    _nb_instances = c._dim[0];
//                        //                    _val->resize(_nb_instances);
//                        _cst = new constant<double>(0);
//                        _all_sign = p_c2->get_all_sign();
//                        _all_range = p_c2->get_range();
//                        _is_transposed = p_c2->_is_transposed;
//                        _is_vector = p_c2->_is_vector;
//                        //                    _is_matrix = p_c2->_is_matrix;
//                        _dim[0] = p_c2->_dim[0];
//                        _dim[1] = p_c2->_dim[1];
//                        //                    if (p_c2->_indices) {
//                        //                        _indices = p_c2->_indices;
//                        //                    }
//                        
//                        break;
//                    }
//                    case uexp_c: {
//                        auto ue = (uexpr*)(&c);
//                        auto f = ue->_son;
//                        for (auto &pair:*f->_vars) {
//                            auto p = pair.second.first;
//                            auto it = _vars->find(p->get_name(false,false));
//                            if (it==_vars->end()) {
//                                add_var(p,pair.second.second);
//                            }
//                        }
//                        for (auto &pair:*f->_params) {
//                            auto p = pair.second.first;
//                            auto it = _params->find(p->get_name(false,false));
//                            if (it==_params->end()) {
//                                add_param(p,pair.second.second);
//                            }
//                        }
//                        switch (ue->_otype) {
//                            case sin_:
//                                _all_range = new pair<double,double>(-1*ue->_coef,ue->_coef); // TODO UPDATE
//                                _all_sign = unknown_;
//                                _all_convexity = undet_;
//                                //                            _val->resize(max(_val->size(),ue->_son->_nb_instances));
//                                //                            if (ue->_son->is_constant()) {
//                                //                                for (unsigned inst = 0; inst < _val->size(); inst++) {
//                                //                                    _val->at(inst) = ue->_coef*std::sin(ue->_son->eval(inst));
//                                //                                }
//                                //                                _evaluated = true;
//                                //                            }
//                                break;
//                            case cos_:
//                                _all_range = new pair<double,double>(-1*ue->_coef,ue->_coef); // TODO UPDATE
//                                _all_sign = unknown_;
//                                _all_convexity = undet_;
//                                //                            _val->resize(max(_val->size(),ue->_son->_nb_instances));
//                                //                            if (ue->_son->is_constant()) {
//                                //                                for (unsigned inst = 0; inst < _val->size(); inst++) {
//                                //                                    _val->at(inst) = ue->_coef*std::cos(ue->_son->eval(inst));
//                                //                                }
//                                //                            }
//                                //                            _evaluated = true;
//                                break;
//                            case sqrt_:
//                                _all_range = new pair<double,double>(0,numeric_limits<double>::max()); // TO UPDATE
//                                _all_sign = non_neg_;
//                                _all_convexity = undet_;
//                                //                            _val->resize(max(_val->size(),ue->_son->_nb_instances));
//                                //                            if (ue->_son->is_constant()) {
//                                //                                for (unsigned inst = 0; inst < _val->size(); inst++) {
//                                //                                    _val->at(inst) = ue->_coef*std::sqrt(ue->_son->eval(inst));
//                                //                                }
//                                //                            }
//                                //                            _evaluated = true;
//                                break;
//                            case exp_:
//                                _all_range = new pair<double,double>(numeric_limits<double>::lowest(),numeric_limits<double>::max()); // TO UPDATE
//                                _all_sign = pos_;
//                                _all_convexity = undet_;
//                                //                            _val->resize(max(_val->size(),ue->_son->_nb_instances));
//                                //                            if (ue->_son->is_constant()) {
//                                //                                for (unsigned inst = 0; inst < _val->size(); inst++) {
//                                //                                    _val->at(inst) = ue->_coef*std::exp(ue->_son->eval(inst));
//                                //                                }
//                                //                            }
//                                //                            _evaluated = true;
//                                break;
//                            case log_:
//                                _all_range = new pair<double,double>(numeric_limits<double>::lowest(),numeric_limits<double>::max()); // TO
//                                _all_sign = unknown_;
//                                _all_convexity = undet_;
//                                //                            _val->resize(max(_val->size(),ue->_son->_nb_instances));
//                                //                            if (ue->_son->is_constant()) {
//                                //                                for (unsigned inst = 0; inst < _val->size(); inst++) {
//                                //                                    _val->at(inst) = ue->_coef*std::log(ue->_son->eval(inst));
//                                //                                }
//                                //                            }
//                                //                            _evaluated = true;
//                                break;
//                            default:
//                                break;
//                        }
//                        _cst = new constant<double>(0);
//                        _expr = make_shared<uexpr>(move(*ue));
//                        _is_transposed = c._is_transposed;
//                        _is_vector = c._is_vector;
//                        //                    _is_matrix = c._is_matrix;
//                        _dim[0] = c._dim[0];
//                        _dim[1] = c._dim[1];
//                        //                    _nb_instances = c._dim[0];
//                        embed(_expr);
//                        if (!_vars->empty()) {
//                            _ftype = nlin_;
//                        }
//                        //                    _DAG->insert(make_pair<>(_expr->get_str(), _expr));
//                        _queue->push_back(_expr);
//                        //sign and convexity
//                        break;
//                    }
//                    case bexp_c: {
//                        auto be = (bexpr*)&c;
//                        auto f = be->_lson;
//                        for (auto &pair:*f->_vars) {
//                            auto p = pair.second.first;
//                            auto it = _vars->find(p->get_name(false,false));
//                            if (it==_vars->end()) {
//                                add_var(p,pair.second.second);
//                            }
//                        }
//                        for (auto &pair:*f->_params) {
//                            auto p = pair.second.first;
//                            auto it = _params->find(p->get_name(false,false));
//                            if (it==_params->end()) {
//                                add_param(p,pair.second.second);
//                            }
//                        }
//                        f = be->_rson;
//                        for (auto &pair:*f->_vars) {
//                            auto p = pair.second.first;
//                            auto it = _vars->find(p->get_name(false,false));
//                            if (it==_vars->end()) {
//                                add_var(p,pair.second.second);
//                            }
//                        }
//                        for (auto &pair:*f->_params) {
//                            auto p = pair.second.first;
//                            auto it = _params->find(p->get_name(false,false));
//                            if (it==_params->end()) {
//                                add_param(p,pair.second.second);
//                            }
//                        }
//                        _cst = new constant<double>(0);
//                        _all_range = new pair<double,double>(numeric_limits<double>::lowest(),numeric_limits<double>::max()); // TODO update
//                        _all_sign = be->get_all_sign();
//                        _all_convexity = undet_;// TODO update
//                        _is_transposed = c._is_transposed;
//                        _is_vector = c._is_vector;
//                        //                    _is_matrix = c._is_matrix;
//                        _dim[0] = c._dim[0];
//                        _dim[1] = c._dim[1];
//                        //                    _nb_instances = c._dim[0];
//                        _expr = make_shared<bexpr>(move(*be));
//                        embed(_expr);
//                        if (!_vars->empty()) {
//                            _ftype = nlin_;
//                        }
//                        //                    _DAG->insert(make_pair<>(_expr->get_str(), _expr));
//                        _queue->push_back(_expr);
//                        break;
//                    }
//                    default:
//                        break;
//                }
//            }
//        }
//        
//        
//
        func(const constant<type>& c):func(){
            *this = c;
        }
        
        func(const param<type>& c):func(){
            *this = c;
        }
        
        func& operator=(const constant<type>& c){
            _cst = c.copy();
            _all_sign = _cst->get_sign();
            param<type>::operator=(c.eval());
            _evaluated = true;
            return *this;
        }
        
        func& operator=(const param<type>& c){
            param<type>::operator=(c);            
            return *this;
        }
                    
//                case par_c:{
//                    auto p_c2 =     shared_ptr<param_>((param_*)copy(c));
//                    //                p_c2->untranspose();//TODO what is this doing here?
//                    _lterms->insert(make_pair<>(p_c2->get_name(false,false), p_c2.get()));
//                    add_param(p_c2);
//                    _cst = new constant<double>(0);
//                    _all_sign = p_c2->get_all_sign();
//                    _all_range = p_c2->get_range();
//                    _is_transposed = c._is_transposed;
//                    _is_vector = c._is_vector;
//                    //                _is_matrix = c._is_matrix;
//                    _dim[0] = c._dim[0];
//                    _dim[1] = c._dim[1];
//                    //                _nb_instances = c._dim[0];
//                    //                _indices = p_c2->_indices;
//                    //                _val->resize(_nb_instances);
//                    //                for (unsigned inst = 0; inst < _val->size(); inst++) {
//                    //                    _val->at(inst) = eval(p_c2.get(), inst);
//                    //                }
//                    _evaluated = false;
//                    break;
//                }
//                case var_c:{
//                    auto p_c2 = shared_ptr<param_>((param_*)copy(c));
//                    //                p_c2->untranspose();
//                    _lterms->insert(make_pair<>(p_c2->get_name(false,false), p_c2.get()));
//                    add_var(p_c2);
//                    _ftype = lin_;
//                    _is_transposed = c._is_transposed;
//                    _is_vector = c._is_vector;
//                    //                _is_matrix = c._is_matrix;
//                    _dim[0] = c._dim[0];
//                    _dim[1] = c._dim[1];
//                    //                _nb_instances = c._dim[0];
//                    //                _val->resize(_nb_instances);
//                    _cst = new constant<double>(0);
//                    _all_sign = p_c2->get_all_sign();
//                    _all_range = p_c2->get_range();
//                    //                _indices = p_c2->_indices;
//                    break;
//                }
//                case uexp_c: {
//                    _cst = new constant<double>(0);
//                    auto ue = (uexpr*)(&c);
//                    auto f = ue->_son;
//                    for (auto &pair:*f->_vars) {
//                        auto p = pair.second.first;
//                        auto it = _vars->find(p->get_name(false,false));
//                        if (it==_vars->end()) {
//                            add_var(p,pair.second.second);
//                        }
//                    }
//                    for (auto &pair:*f->_params) {
//                        auto p = pair.second.first;
//                        auto it = _params->find(p->get_name(false,false));
//                        if (it==_params->end()) {
//                            add_param(p,pair.second.second);
//                        }
//                    }
//                    switch (ue->_otype) {
//                        case sin_:
//                            _all_range = new pair<double,double>(-1*ue->_coef,ue->_coef); // TO UPDATE
//                            _all_sign = unknown_;
//                            _all_convexity = undet_;// TODO update
//                            //                            _val->resize(max(_val->size(),ue->_son->_nb_instances));
//                            //                            if (ue->_son->is_constant()) {
//                            //                                for (unsigned inst = 0; inst < _val->size(); inst++) {
//                            //                                    _val->at(inst) = ue->_coef*std::sin(ue->_son->eval(inst));
//                            //                                }
//                            //                                _evaluated = true;
//                            //                            }
//                            break;
//                        case cos_:
//                            _all_range = new pair<double,double>(-1*ue->_coef,ue->_coef); // TODO UPDATE
//                            _all_sign = unknown_;
//                            _all_convexity = undet_;// TODO update
//                            //                            _val->resize(max(_val->size(),ue->_son->_nb_instances));
//                            //                            if (ue->_son->is_constant()) {
//                            //                                for (unsigned inst = 0; inst < _val->size(); inst++) {
//                            //                                    _val->at(inst) = ue->_coef*std::cos(ue->_son->eval(inst));
//                            //                                }
//                            //                            }
//                            //                            _evaluated = true;
//                            break;
//                        case sqrt_:
//                            _all_range = new pair<double,double>(0,numeric_limits<double>::max()); // TO UPDATE
//                            _all_sign = non_neg_;
//                            _all_convexity = undet_;// TODO update
//                            //                            _val->resize(max(_val->size(),ue->_son->_nb_instances));
//                            //                            if (ue->_son->is_constant()) {
//                            //                                for (unsigned inst = 0; inst < _val->size(); inst++) {
//                            //                                    _val->at(inst) = ue->_coef*std::sqrt(ue->_son->eval(inst));
//                            //                                }
//                            //                            }
//                            //                            _evaluated = true;
//                            break;
//                        case exp_:
//                            _all_range = new pair<double,double>(numeric_limits<double>::lowest(),numeric_limits<double>::max()); // TO UPDATE
//                            _all_sign = pos_;
//                            _all_convexity = undet_;// TODO update
//                            //                            _val->resize(max(_val->size(),ue->_son->_nb_instances));
//                            //                            if (ue->_son->is_constant()) {
//                            //                                for (unsigned inst = 0; inst < _val->size(); inst++) {
//                            //                                    _val->at(inst) = ue->_coef*std::exp(ue->_son->eval(inst));
//                            //                                }
//                            //                            }
//                            //                            _evaluated = true;
//                            break;
//                        case log_:
//                            _all_range = new pair<double,double>(numeric_limits<double>::lowest(),numeric_limits<double>::max()); // TO
//                            _all_sign = unknown_;
//                            _all_convexity = undet_;// TODO update
//                            //                            _val->resize(max(_val->size(),ue->_son->_nb_instances));
//                            //                            if (ue->_son->is_constant()) {
//                            //                                for (unsigned inst = 0; inst < _val->size(); inst++) {
//                            //                                    _val->at(inst) = ue->_coef*std::log(ue->_son->eval(inst));
//                            //                                }
//                            //                            }
//                            //                            _evaluated = true;
//                            break;
//                        default:
//                            break;
//                    }
//                    _expr = make_shared<uexpr>(*ue);
//                    embed(_expr);
//                    _is_transposed = c._is_transposed;
//                    _is_vector = c._is_vector;
//                    //                _is_matrix = c._is_matrix;
//                    _dim[0] = c._dim[0];
//                    _dim[1] = c._dim[1];
//                    //                _nb_instances = c._dim[0];
//                    if (!_vars->empty()) {
//                        _ftype = nlin_;
//                    }
//                    //                _DAG->insert(make_pair<>(_expr->get_str(), _expr));
//                    _queue->push_back(_expr);
//                    //sign and convexity
//                    break;
//                }
//                case bexp_c: {
//                    auto be = (bexpr*)&c;
//                    auto f = be->_lson;
//                    for (auto &pair:*f->_vars) {
//                        auto p = pair.second.first;
//                        auto it = _vars->find(p->get_name(false,false));
//                        if (it==_vars->end()) {
//                            add_var(p,pair.second.second);
//                        }
//                    }
//                    for (auto &pair:*f->_params) {
//                        auto p = pair.second.first;
//                        auto it = _params->find(p->get_name(false,false));
//                        if (it==_params->end()) {
//                            add_param(p,pair.second.second);
//                        }
//                    }
//                    f = be->_rson;
//                    for (auto &pair:*f->_vars) {
//                        auto p = pair.second.first;
//                        auto it = _vars->find(p->get_name(false,false));
//                        if (it==_vars->end()) {
//                            add_var(p,pair.second.second);
//                        }
//                    }
//                    for (auto &pair:*f->_params) {
//                        auto p = pair.second.first;
//                        auto it = _params->find(p->get_name(false,false));
//                        if (it==_params->end()) {
//                            add_param(p,pair.second.second);
//                        }
//                    }
//                    _cst = new constant<double>(0);
//                    _expr = make_shared<bexpr>(*be);
//                    _all_range = new pair<double,double>(numeric_limits<double>::lowest(),numeric_limits<double>::max()); // TO UPDATE
//                    _all_sign = be->get_all_sign();
//                    _all_convexity = undet_;// TODO update
//                    if (!_val) {
//                        _val = make_shared<vector<double>>();
//                    }
//                    //                if (be->_lson->is_constant() && be->_rson->is_constant()) {
//                    //                    _val->resize(max(_val->size(),max(be->_lson->_nb_instances,be->_rson->_nb_instances)));
//                    //                    switch (be->_otype) {
//                    //                        case plus_:
//                    //                            for (unsigned inst = 0; inst < _val->size(); inst++) {
//                    //                                _val->at(inst) = be->_coef*(be->_lson->eval(inst) + be->_rson->eval(inst));
//                    //                            }
//                    //                            break;
//                    //                        case minus_:
//                    //                            for (unsigned inst = 0; inst < _val->size(); inst++) {
//                    //                                _val->at(inst) = be->_coef*(be->_lson->eval(inst) - be->_rson->eval(inst));
//                    //                            }
//                    //                            break;
//                    //                        case product_:
//                    //                            for (unsigned inst = 0; inst < _val->size(); inst++) {
//                    //                                _val->at(inst) = be->_coef*(be->_lson->eval(inst) * be->_rson->eval(inst));
//                    //                            }
//                    //                            break;
//                    //                        case div_:
//                    //                            for (unsigned inst = 0; inst < _val->size(); inst++) {
//                    //                                _val->at(inst) = be->_coef*(be->_lson->eval(inst) / be->_rson->eval(inst));
//                    //                            }
//                    //                            break;
//                    //                        default:
//                    //                            break;
//                    //                    }
//                    //                }
//                    embed(_expr);
//                    _is_transposed = c._is_transposed;
//                    _is_vector = c._is_vector;
//                    //                _is_matrix = c._is_matrix;
//                    _dim[0] = c._dim[0];
//                    _dim[1] = c._dim[1];
//                    //                _nb_instances = c._dim[0];
//                    if (!_vars->empty()) {
//                        _ftype = nlin_;
//                    }
//                    //                _DAG->insert(make_pair<>(_expr->get_str(), _expr));
//                    _queue->push_back(_expr);
//                    break;
//                }
//                case func_c: {
//                    *this = *(func_*)&c;
//                }
//                default:
//                    break;
//            }
//            _dim[0] = c._dim[0];
//            _dim[1] = c._dim[1];
//        }

        func(func&& f){
            *this = move(f);
        }
        
        func(const func& f){
            *this = f;
        }
        
        shared_ptr<func_> fcopy() const{return make_shared<func>(*this);};
        
        func& operator=(const func& f){
            constant_::_type = f._type;
            _ftype = f._ftype;
            _return_type = f._return_type;
            _to_str = f._to_str;
            _all_convexity = f._all_convexity;
            _all_sign = f._all_sign;
            _cst = f._cst->copy();
            _lterms = make_shared<map<string, lterm>>(*f._lterms);
            _qterms = make_shared<map<string, qterm>>(*f._qterms);
            _pterms = make_shared<map<string, pterm>>(*f._pterms);
            if(f._expr){
                _expr = make_shared<expr>(*f._expr);
            }
            _vars = make_shared<map<string, pair<shared_ptr<param_>, unsigned>>>();
            _params = make_shared<map<string, pair<shared_ptr<param_>, unsigned>>>();
            merge_vars(*this);// Update vars and params
            if(f._indices){
                param_::_indices = make_shared<indices>(*f._indices);
            }
            param<type>::_range = make_shared<pair<type,type>>(*f._range);
            param<type>::_val = make_shared<vector<type>>(*f._val);
            *_convexity = *f._convexity;
            _sign = f._sign;
            constant_::_is_transposed = f._is_transposed;
            constant_::_is_vector = f._is_vector;
            _is_constraint = f._is_constraint;
            _is_hessian = f._is_hessian;
            constant_::_dim[0] = f._dim[0];
            constant_::_dim[1] = f._dim[1];
            _embedded = f._embedded;
            _dfdx = make_shared<map<string,shared_ptr<func_>>>();
            copy_derivatives(f);
            if(f._hess_link){
                _hess_link = make_shared<map<size_t, set<size_t>>>(*f._hess_link);
            }
            _nnz_j = f._nnz_j;
            _nnz_h = f._nnz_h;
            _hess_link = f._hess_link;
            _nb_vars = f._nb_vars;
            return *this;
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
            param_::_indices = move(f._indices);
            param<type>::_range = move(f._range);
            param<type>::_val = move(f._val);
            _convexity = move(f._convexity);
            _sign = f._sign;
            f._sign = nullptr;
            constant_::_is_transposed = f._is_transposed;
            constant_::_is_vector = f._is_vector;
            _is_constraint = f._is_constraint;
            _is_hessian = f._is_hessian;
            constant_::_dim[0] = f._dim[0];
            constant_::_dim[1] = f._dim[1];
            _embedded = f._embedded;
            _dfdx = move(f._dfdx);
            _nnz_j = f._nnz_j;
            _nnz_h = f._nnz_h;
            _hess_link = f._hess_link;
            _nb_vars = f._nb_vars;
            return *this;
        }

//        
//        func_ tr() const {
//            auto f = func_(*this);
//            f.transpose();
//            return f;
//        }
//
    };
//
//
//
//    func_ operator+(const constant_& c1, const constant_& c2);
////    func_ operator+(func_&& f, const constant_& c);
//    //func_ operator+(const constant_& c, func_&& f);
//
//    template<class T, class = typename enable_if<std::is_arithmetic<T>::value>::type> func_ operator+(func_&& f, T c){
//        return f += c;
//    };
//
//    template<class T, class = typename enable_if<std::is_arithmetic<T>::value>::type> func_ operator+(T c, func_&& f){
//        return f += c;
//    };
//
//    template<class T, class = typename enable_if<std::is_arithmetic<T>::value>::type> func_ operator-(func_&& f, T c){
//        return f -= c;
//    };
//
//    template<class T, class = typename enable_if<std::is_arithmetic<T>::value>::type> func_ operator-(T c, func_&& f){
//        return (f *= -1) += c;
//    };
//
//    template<class T, class = typename enable_if<std::is_arithmetic<T>::value>::type> func_ operator*(func_&& f, T c){
//        return f *= c;
//    };
//
//    template<class T, class = typename enable_if<std::is_arithmetic<T>::value>::type> func_ operator*(T c, func_&& f){
//        return f *= c;
//    };
//
//
//    template<class T, class = typename enable_if<std::is_arithmetic<T>::value>::type> func_ operator/(func_&& f, T c){
//        return f /= c;
//    };
//
//    template<class T, class = typename enable_if<std::is_arithmetic<T>::value>::type> func_ operator+(const constant_& c1, T c2){
//        return func_(c1) += c2;
//    };
//
//    template<class T, class = typename enable_if<std::is_arithmetic<T>::value>::type> func_ operator+(T c2, const constant_& c1){
//        return func_(c1) += c2;
//    };
//
//    template<class T, class = typename enable_if<std::is_arithmetic<T>::value>::type> func_ operator-(const constant_& c1, T c2){
//        return func_(c1) -= c2;
//    };
//
//    template<class T, class = typename enable_if<std::is_arithmetic<T>::value>::type> func_ operator-(T c2, const constant_& c1){
//        return (func_(c1) *= -1) += c2;
//    };
//
//    template<class T, class = typename enable_if<std::is_arithmetic<T>::value>::type> func_ operator*(const constant_& c1, T c2){
//        return func_(c1) *= c2;
//    };
//
//    template<class T, class = typename enable_if<std::is_arithmetic<T>::value>::type> func_ operator*(T c2, const constant_& c1){
//        return func_(c1) *= c2;
//    };
//
//    template<class T, class = typename enable_if<std::is_arithmetic<T>::value>::type> func_ operator/(const constant_& c1, T c2){
//        return func_(c1) *= 1/c2;
//    };
//    
//    template<class T, class = typename enable_if<std::is_arithmetic<T>::value>::type> func_ operator/(T c2, const constant_& c1){
//        return func_(c2) /= c1;
//    };
//
//    func_ operator*(const constant_& c1, const constant_& c2);
//
//
//
//    func_ operator-(const constant_& c1, const constant_& c2);
//
//    func_ operator/(const constant_& c1, const constant_& c2);
//
//
//
//        
//        
//    constant_* add(constant_* c1, const constant_& c2); /**< adds c2 to c1, updates its type and returns the result **/
//    //
//    template<class T, class = typename enable_if<is_arithmetic<T>::value>::type> constant_* add(constant_* c1, const constant<T>& c2){ /**< adds c2 to c1, updates its type and returns the result **/
//        switch (c1->get_type()) {
//            case binary_c: {
//                if (c2.is_binary() ) {
//                    *(constant<bool>*)c1 += c2.eval();
//                }
//                else {
//                    bool val = ((constant<bool>*)c1)->eval();
//                    delete c1;
//                    c1 = new constant<T>(c2.eval() + val);
//                }
//                return c1;
//                break;
//            }
//            case short_c: {
//                if (c2.get_type() <= short_c) {
//                    *((constant<short>*)c1) += c2.eval();
//                }
//                else {
//                    auto val = ((constant<short>*)c1)->eval();
//                    delete c1;
//                    c1 = new constant<T>(c2.eval() + val);
//                }
//                break;
//            }
//            case integer_c: {
//                if (c2.get_type() <= integer_c) {
//                    *((constant<int>*)c1) += c2.eval();
//                }
//                else {
//                    auto val = ((constant<int>*)c1)->eval();
//                    delete c1;
//                    c1 = new constant<T>(c2.eval() + val);
//                }
//                break;
//            }
//            case float_c: {
//                if (c2.get_type() <= float_c) {
//                    *((constant<float>*)c1) += c2.eval();
//                }
//                else {
//                    auto val = ((constant<float>*)c1)->eval();
//                    delete c1;
//                    c1 = new constant<T>(c2.eval() + val);
//                }
//                break;
//            }
//            case double_c: {
//                if (c2.get_type() <= double_c) {
//                    *((constant<double>*)c1) += c2.eval();
//                }
//                else {
//                    auto val = ((constant<double>*)c1)->eval();
//                    delete c1;
//                    c1 = new constant<T>(c2.eval() + val);
//                }
//                break;
//            }
//            case long_c: {
//                if (c2.get_type() <= long_c) {
//                    *((constant<long double>*)c1) += c2.eval();
//                }
//                else {
//                    auto val = ((constant<long double>*)c1)->eval();
//                    delete c1;
//                    c1 = new constant<T>(c2.eval() + val);
//                }
//                break;
//            }
//            case complex_c: {
//                auto f = new func_(*c1);
//                *f += c2;
//                c1 =(constant_*)(f);
//                return c1;
//                break;
//            }
//            case par_c:{            
//                auto f = new func_(*c1);
//                *f += c2;
//                c1 =(constant_*)(f);
//                return c1;
//                break;
//            }
//            case var_c:{
//                auto f = new func_(*c1);
//                *f += c2;
//                c1 =(constant_*)(f);
//                return c1;
//                break;
//            }
//    //        case uexp_c: {
//    //            auto res = new bexpr(*(uexpr*)c1 + c2);
//    //            delete c1;
//    //            c1 = (constant_*)res;
//    //            return c1;
//    //            break;
//    //        }
//    //        case bexp_c: {
//    //            auto res = new bexpr(*(bexpr*)c1 + c2);
//    //            delete c1;
//    //            c1 = (constant_*)res;
//    //            return c1;
//    //            break;
//    //        }
//            case func_c: {
//    //            auto res = new func_((*(func_*)c1) + c2);
//    //            delete c1;
//    //            return c1 = (constant_*)res;
//                (*(func_*)c1) += c2;
//                return c1;
//                break;
//            }
//                
//            default:
//                break;
//        }
//        return c1;
//    }
//
//    //constant_* add(constant_* c1, const func_& f);
//
//    template<class T> constant_* add(constant_* c1, const param<T>& c2){ /**< adds c2 to c1, updates its type and returns the result **/
//        switch (c1->get_type()) {
//            case binary_c: {
//                auto val = ((constant<bool>*)c1)->eval();
//                delete c1;
//                auto f = new func_(c2);
//                *f += val;
//                c1 = (constant_*)(f);
//                return c1;
//                break;
//            }
//            case short_c: {
//                auto val = ((constant<short>*)c1)->eval();
//                delete c1;
//                auto f = new func_(c2);
//                *f += val;
//                c1 = (constant_*)(f);
//                return c1;
//                break;
//            }
//            case integer_c: {
//                auto val = ((constant<int>*)c1)->eval();
//                delete c1;
//                auto f = new func_(c2);
//                *f += val;
//                c1 = (constant_*)(f);
//                return c1;
//                break;
//            }
//            case float_c: {
//                auto val = ((constant<float>*)c1)->eval();
//                delete c1;
//                auto f = new func_(c2);
//                *f += val;
//                c1 = (constant_*)(f);
//                return c1;
//                break;
//            }
//            case double_c: {
//                auto val = ((constant<double>*)c1)->eval();
//                delete c1;
//                auto f = new func_(c2);
//                *f += val;
//                c1 = (constant_*)(f);
//                return c1;
//                break;
//            }
//            case long_c: {
//                auto val = ((constant<long double>*)c1)->eval();
//                delete c1;
//                auto f = new func_(c2);
//                *f += val;
//                c1 = (constant_*)(f);
//                return c1;
//                break;
//            }
//            case complex_c: {
//                auto val = *((constant<Cpx>*)c1);
//                delete c1;
//                auto f = new func_(c2);
//                *f += val;
//                c1 = (constant_*)(f);
//                return c1;
//                break;
//            }
//            case par_c:{
//                auto res = new func_(*c1);
//                delete c1;
//                res->insert(true, constant<double>(1), c2);
//                return c1 = res;
//                break;
//            }
//            case var_c:{
//                auto res = new func_(*c1);
//                delete c1;
//                if (c2.is_var()) {
//                    res->insert(true, constant<double>(1), c2);
//                }
//                else {
//                    auto cst = res->get_cst();
//                    cst = add(cst, c2);
//                }
//                return c1 = res;
//                break;
//            }
//
//    //        case uexp_c: {
//    //            auto res = new bexpr(*(uexpr*)c1 + c2);
//    //            delete c1;
//    //            c1 = (constant_*)res;
//    //            return c1;
//    //            break;
//    //        }
//    //        case bexp_c: {
//    //            auto res = new bexpr(*(bexpr*)c1 + c2);
//    //            delete c1;
//    //            c1 = (constant_*)res;
//    //            return c1;
//    //            break;
//    //        }
//            case func_c: {
//    //            auto res = new func_(*c1);
//    //            delete c1;
//    //            *res += c2;
//    //            return c1 = (constant_*)res;
//                (*(func_*)c1) += c2;
//                return c1;
//                break;
//            }
//            default:
//                break;
//        }
//        return c1;
//    }
//
//    constant_* substract(constant_* c1, const constant_& c2);
//
//
//    template<class T, class = typename enable_if<is_arithmetic<T>::value>::type> constant_* substract(constant_* c1, const param<T>& c2){ /**< Substracts c2 from c1, updates its type and returns the result **/
//        switch (c1->get_type()) {
//            case binary_c: {
//                auto val = ((constant<bool>*)c1)->eval();
//                delete c1;
//                auto f = new func_(c2);
//                (*f *= -1) += val;
//                c1 = (constant_*)(f);
//                return c1;
//                break;
//            }
//            case short_c: {
//                auto val = ((constant<short>*)c1)->eval();
//                delete c1;
//                auto f = new func_(c2);
//                (*f *= -1) += val;
//                c1 = (constant_*)(f);
//                return c1;
//                break;
//            }
//            case integer_c: {
//                auto val = ((constant<int>*)c1)->eval();
//                delete c1;
//                auto f = new func_(c2);
//                (*f *= -1) += val;
//                c1 = (constant_*)(f);
//                return c1;
//                break;
//            }
//            case float_c: {
//                auto val = ((constant<float>*)c1)->eval();
//                delete c1;
//                auto f = new func_(c2);
//                (*f *= -1) += val;
//                c1 = (constant_*)(f);
//                return c1;
//                break;
//            }
//            case double_c: {
//                auto val = ((constant<double>*)c1)->eval();
//                delete c1;
//                auto f = new func_(c2);
//                (*f *= -1) += val;
//                c1 = (constant_*)(f);
//                return c1;
//                break;
//            }
//            case long_c: {
//                auto val = ((constant<long double>*)c1)->eval();
//                delete c1;
//                auto f = new func_(c2);
//                (*f *= -1) += val;
//                c1 = (constant_*)(f);
//                return c1;
//                break;
//            }
//            case complex_c:{
//                auto newcst = (*(constant<Cpx>*)c1);
//                delete c1;
//                auto res = new func_(c2);
//                *res *= -1;
//                auto cst = res->get_cst();
//                cst = add(cst, c2);
//                return c1 = res;
//                break;
//            }
//            case par_c:{
//                auto res = new func_(*c1);
//                delete c1;
//                res->insert(false, constant<double>(1), c2);
//                return c1 = res;
//                break;
//            }
//            case var_c:{
//                auto res = new func_(*c1);
//                delete c1;
//                if (c2.is_var()) {
//                    res->insert(false, constant<double>(1), c2);
//                }
//                else {
//                    auto cst = res->get_cst();
//                    cst = substract(cst, c2);
//                }
//                return c1 = res;
//                break;
//            }
//
//            case uexp_c: {
//                auto res = new bexpr(minus_, make_shared<func_>(func_(*c1)), make_shared<func_>(func_(c2)));
//                delete c1;
//                c1 = (constant_*)res;
//                return c1;
//                break;
//            }
//            case bexp_c: {
//                auto res = new bexpr(minus_, make_shared<func_>(func_(*c1)), make_shared<func_>(func_(c2)));
//                delete c1;
//                c1 = (constant_*)res;
//                return c1;
//                break;
//            }
//            case func_c: {
//    //            auto res = new func_(*c1);
//    //            delete c1;
//    //            *res -= c2;
//    //            return c1 = res;
//                (*(func_*)c1) -= c2;
//                return c1;
//                break;
//            }
//            default:
//                break;
//        }
//        return c1;
//    }
//
//
//
//    template<typename T> constant_* multiply(constant_* c1, const param<T>& c2){ /**< adds c2 to c1, updates its type and returns the result **/
//        switch (c1->get_type()) {
//            case binary_c: {
//                auto val = ((constant<bool>*)c1)->eval();
//                if (val==true) {
//                    delete c1;
//                    return c1 = new param<T>(c2);
//                }
////                c1->update_dot_dim(c2);
//                return c1;
//                break;
//            }
//            case short_c: {
//                auto val = ((constant<short>*)c1)->eval();
//                if (val==0) {
//                    return c1;
//                }
//                delete c1;
////                if (val==1) {
////                    return c1 = new param<T>(c2);
////                }
//                auto f = new func_(c2);
//                *f *= val;
//                c1 = (constant_*)(f);
////                c1->update_dot_dim(c2);
//                return c1;
//                break;
//            }
//            case integer_c: {
//                auto val = ((constant<int>*)c1)->eval();
//                if (val==0) {
//                    return c1;
//                }
//                delete c1;
////                if (val==1) {
////                    return c1 = new param<T>(c2);
////                }
//                auto f = new func_(c2);
//                *f *= val;
//                c1 = (constant_*)(f);
////                c1->update_dot_dim(c2);
//                return c1;
//                break;
//            }
//            case float_c: {
//                auto val = ((constant<float>*)c1)->eval();
//                if (val==0) {
//                    return c1;
//                }
//                delete c1;
////                if (val==1) {
////                    return c1 = new param<T>(c2);
////                }
//                auto f = new func_(c2);
//                *f *= val;
//                c1 = (constant_*)(f);
////                c1->update_dot_dim(c2);
//                return c1;
//                break;
//            }
//            case double_c: {
//                auto val = ((constant<double>*)c1)->eval();
//                if (val==0) {
//                    return c1;
//                }
//                delete c1;
////                if (val==1) {
////                    return c1 = new param<T>(c2);
////                }
//                auto f = new func_(c2);
//                *f *= val;
//                c1 = (constant_*)(f);
////                auto c = new param<T>(c2);
////                
////                if (c->_is_matrix) {
////                    for (unsigned i = 0; i<c->_dim[0]; i++) {
////                        for (unsigned j = 0; j<c->_dim[1]; j++) {
////    //                    <#statements#>
////                    }
////                }
////                *f *= val;
////                c1 = c;
////                c1->update_dot_dim(c2);
//                return c1;
//                break;
//            }
//            case long_c: {
//                auto val = ((constant<long double>*)c1)->eval();
//                if (val==0) {
//                    return c1;
//                }
//                delete c1;
////                if (val==1) {
////                    return c1 = new param<T>(c2);
////                }
//                auto f = new func_(c2);
//                *f *= val;
//                c1 = (constant_*)(f);
////                c1->update_dot_dim(c2);
//                return c1;
//                break;
//            }
//            case complex_c:{
//                auto f = new func_(*c1);
//                delete c1;
//                *f *= c2;
//                c1 =(constant_*)(f);
////                c1->update_dot_dim(c2);
//                return c1;
//                break;
//            }
//            case var_c:{
//                auto f = new func_(*c1);
//                delete c1;
//                *f *= c2;
//                c1 =(constant_*)(f);
////                c1->update_dot_dim(c2);
//                return c1;
//                break;
//            }
//            case par_c:{
//                auto f = new func_(*c1);
//                delete c1;
//                *f *= c2;
//                c1 =(constant_*)(f);
////                c1->update_dot_dim(c2);
//                return c1;
//                break;
//            }
//            case uexp_c: {
//                auto res = new bexpr(product_, make_shared<func_>(func_(*c1)), make_shared<func_>(func_(c2)));
//                delete c1;
//                c1 = (constant_*)res;
////                c1->update_dot_dim(c2);
//                return c1;
//                break;
//            }
//            case bexp_c: {
//                auto res = new bexpr(product_, make_shared<func_>(func_(*c1)), make_shared<func_>(func_(c2)));
//                delete c1;
//                c1 = (constant_*)res;
////                c1->update_dot_dim(c2);
//                return c1;
//                break;
//            }
//            case func_c: {
//    //            auto res = new func_(*c1);
//    //            delete c1;
//    //            *res *= c2;
//    //            return c1 = res;
//                (*(func_*)c1) *= c2;
////                c1->update_dot_dim(c2);
//                return c1;
//                break;
//            }
//            default:
//                break;
//        }
//        return c1;
//    }
//
//
//
//
//
//    template<class T, class = typename enable_if<is_arithmetic<T>::value>::type> constant_* substract(constant_* c1, const constant<T>& c2){ /**< adds c2 to c1, updates its type and returns the result **/
//        switch (c1->get_type()) {
//            case binary_c: {
//                if (c2.is_binary() ) {
//                    *(constant<bool>*)c1 -= c2.eval();
//                }
//                else {
//                    auto val = ((constant<bool>*)c1)->eval();
//                    delete c1;
//                    c1 = new constant<T>(val - c2.eval());
//                }
//                return c1;
//                break;
//            }
//            case short_c: {
//                if (c2.get_type() <= short_c) {
//                    *((constant<short>*)c1) -= c2.eval();
//                }
//                else {
//                    auto val = ((constant<short>*)c1)->eval();
//                    delete c1;
//                    c1 = new constant<T>(val - c2.eval());
//                }
//                break;
//            }
//            case integer_c: {
//                if (c2.get_type() <= integer_c) {
//                    *((constant<int>*)c1) -= c2.eval();
//                }
//                else {
//                    auto val = ((constant<int>*)c1)->eval();
//                    delete c1;
//                    c1 = new constant<T>(val - c2.eval());
//                }
//                break;
//            }
//            case float_c: {
//                if (c2.get_type() <= float_c) {
//                    *((constant<float>*)c1) -= c2.eval();
//                }
//                else {
//                    auto val = ((constant<float>*)c1)->eval();
//                    delete c1;
//                    c1 = new constant<T>(val - c2.eval());
//                }
//                break;
//            }
//            case double_c: {
//                if (c2.get_type() <= double_c) {
//                    *((constant<double>*)c1) -= c2.eval();
//                }
//                else {
//                    auto val = ((constant<double>*)c1)->eval();
//                    delete c1;
//                    c1 = new constant<T>(val - c2.eval());
//                }
//                break;
//            }
//            case long_c: {
//                if (c2.get_type() <= long_c) {
//                    *((constant<long double>*)c1) -= c2.eval();
//                }
//                else {
//                    auto val = ((constant<long double>*)c1)->eval();
//                    delete c1;
//                    c1 = new constant<T>(val - c2.eval());
//                }
//                break;
//            }
//            case complex_c: {
//                (*(constant<Cpx>*)c1) -= c2.eval();
//                break;
//            }
//            case par_c:{
//                auto f = new func_(*c1);
//                *f -= c2;
//                c1 =(constant_*)(f);
//                return c1;
//                break;
//            }
//            case var_c:{
//                auto f = new func_(*c1);
//                *f -= c2;
//                c1 =(constant_*)(f);
//                return c1;
//                break;
//            }
//            case uexp_c: {
//                auto res = new bexpr(minus_, make_shared<func_>(func_(*c1)), make_shared<func_>(func_(c2)));
//                delete c1;
//                c1 = (constant_*)res;
//                return c1;
//                break;
//            }
//            case bexp_c: {
//                auto res = new bexpr(minus_, make_shared<func_>(func_(*c1)), make_shared<func_>(func_(c2)));
//                delete c1;
//                c1 = (constant_*)res;
//                return c1;
//                break;
//            }
//            case func_c: {
//    //            auto res = new func_(*(func_*)c1 - c2);
//    //            delete c1;
//    //            return c1 = (constant_*)res;
//                (*(func_*)c1) -= c2;
//                return c1;
//                break;
//            }        default:
//                break;
//        }
//        return c1;
//    }
//
//    constant_* multiply(constant_* c1, const constant_& c2);
//    constant_* divide(constant_* c1, const constant_& c2);
//
//
//    template<class T, class = typename enable_if<is_arithmetic<T>::value>::type> constant_* multiply(constant_* c1, const constant<T>& c2){ /**< adds c2 to c1, updates its type and returns the result **/
//        switch (c1->get_type()) {
//            case binary_c: {
//                if (c2.is_binary() ) {
//                    *(constant<bool>*)c1 *= c2.eval();
//                }
//                else {
//                    auto val = ((constant<bool>*)c1)->eval();
//                    delete c1;
//                    c1 = new constant<T>(c2.eval() * val);
//                }
//                return c1;
//                break;
//            }
//            case short_c: {
//                if (c2.get_type() <= short_c) {
//                    *((constant<short>*)c1) *= c2.eval();
//                }
//                else {
//                    auto val = ((constant<short>*)c1)->eval();
//                    delete c1;
//                    c1 = new constant<T>(c2.eval() * val);
//                }
//                break;
//            }
//            case integer_c: {
//                if (c2.get_type() <= integer_c) {
//                    *((constant<int>*)c1) *= c2.eval();
//                }
//                else {
//                    auto val = ((constant<int>*)c1)->eval();
//                    delete c1;
//                    c1 = new constant<T>(c2.eval() * val);
//                }
//                break;
//            }
//            case float_c: {
//                if (c2.get_type() <= float_c) {
//                    *((constant<float>*)c1) *= c2.eval();
//                }
//                else {
//                    auto val = ((constant<float>*)c1)->eval();
//                    delete c1;
//                    c1 = new constant<T>(c2.eval() * val);
//                }
//                break;
//            }
//            case double_c: {
//                if (c2.get_type() <= double_c) {
//                    *((constant<double>*)c1) *= c2.eval();
//                }
//                else {
//                    auto val = ((constant<double>*)c1)->eval();
//                    delete c1;
//                    c1 = new constant<T>(c2.eval() * val);
//                }
//                break;
//            }
//            case long_c: {
//                if (c2.get_type() <= long_c) {
//                    *((constant<long double>*)c1) *= c2.eval();
//                }
//                else {
//                    auto val = ((constant<long double>*)c1)->eval();
//                    delete c1;
//                    c1 = new constant<T>(c2.eval() * val);
//                }
//                break;
//            }
//            case complex_c: {
//                (*(constant<Cpx>*)c1) *= c2.eval();
//                break;
//            }
//            case par_c:{
//                auto pc1 = (param_*)(c1);
//                auto l = new func_(*pc1);
//                *l *= c2;
//                c1 =(constant_*)(l);
//                return c1;
//                break;
//            }
//            case uexp_c: {
//                auto res = new bexpr(product_, make_shared<func_>(func_(*c1)), make_shared<func_>(func_(c2)));
//                delete c1;
//                c1 = (constant_*)res;
//                return c1;
//                break;
//            }
//            case bexp_c: {
//                auto res = new bexpr(product_, make_shared<func_>(func_(*c1)), make_shared<func_>(func_(c2)));
//                delete c1;
//                c1 = (constant_*)res;
//                return c1;
//                break;
//            }
//            case func_c: {
//                (*(func_*)c1) *= c2;
//                return c1;
//                break;
//            }
//            default:
//                break;
//        }
//        return c1;
//    }
//
//    template<class T, class = typename enable_if<is_arithmetic<T>::value>::type> constant_* divide(constant_* c1, const constant<T>& c2){ /**< adds c2 to c1, updates its type and returns the result **/
//        if (c2.eval()==0) {
//            throw invalid_argument("dividing by zero!\n");
//        }
//        switch (c1->get_type()) {
//            case binary_c: {
//                if (c2.is_binary() ) {
//                    *(constant<bool>*)c1 /= c2.eval();
//                }
//                else {
//                    auto val = ((constant<bool>*)c1)->eval();
//                    delete c1;
//                    c1 = new constant<T>(c2.eval() / val);
//                }
//                return c1;
//                break;
//            }
//            case short_c: {
//                if (c2.get_type() <= short_c) {
//                    *((constant<short>*)c1) /= c2.eval();
//                }
//                else {
//                    auto val = ((constant<short>*)c1)->eval();
//                    delete c1;
//                    c1 = new constant<T>(c2.eval() / val);
//                }
//                break;
//            }
//            case integer_c: {
//                if (c2.get_type() <= integer_c) {
//                    *((constant<int>*)c1) /= c2.eval();
//                }
//                else {
//                    auto val = ((constant<int>*)c1)->eval();
//                    delete c1;
//                    c1 = new constant<T>(c2.eval() / val);
//                }
//                break;
//            }
//            case float_c: {
//                if (c2.get_type() <= float_c) {
//                    *((constant<float>*)c1) /= c2.eval();
//                }
//                else {
//                    auto val = ((constant<float>*)c1)->eval();
//                    delete c1;
//                    c1 = new constant<T>(c2.eval() / val);
//                }
//                break;
//            }
//            case double_c: {
//                if (c2.get_type() <= double_c) {
//                    *((constant<double>*)c1) /= c2.eval();
//                }
//                else {
//                    auto val = ((constant<double>*)c1)->eval();
//                    delete c1;
//                    c1 = new constant<T>(c2.eval() / val);
//                }
//                break;
//            }
//            case long_c: {
//                if (c2.get_type() <= long_c) {
//                    *((constant<long double>*)c1) /= c2.eval();
//                }
//                else {
//                    auto val = ((constant<long double>*)c1)->eval();
//                    delete c1;
//                    c1 = new constant<T>(c2.eval() / val);
//                }
//                break;
//            }
//            case complex_c: {
//                (*(constant<Cpx>*)c1) /= c2.eval();
//                break;
//            }
//            case par_c:{
//                auto pc1 = (param_*)(c1);
//                auto l = new func_(*pc1);
//                *l /= c2;
//                c1 =(constant_*)(l);
//                return c1;
//                break;
//            }
//            case uexp_c: {
//                auto res = new bexpr(div_, make_shared<func_>(func_(*c1)), make_shared<func_>(func_(c2)));
//                delete c1;
//                c1 = (constant_*)res;
//                return c1;
//                break;
//            }
//            case bexp_c: {
//                auto res = new bexpr(div_, make_shared<func_>(func_(*c1)), make_shared<func_>(func_(c2)));
//                delete c1;
//                c1 = (constant_*)res;
//                return c1;
//                break;
//            }
//            case func_c: {
//                switch (((func_*)c1)->get_ftype()) {
//                    case lin_: {
//                        cerr << "Unsupported yet;\n";
//                        //                    auto res = new func_(*(func_*)c1 * 1/c2);
//                        //                    delete c1;
//                        //                    return c1 = res;
//                        break;
//                    }
//                    default:
//                        cerr << "Unsupported yet;\n";
//                        break;
//                }
//            }
//            default:
//                break;
//        }
//        return c1;
//    }
//
//
//    func cos(const constant_& c);
//
//    func sin(const constant_& c);
//
//
//    func sqrt(const constant_& c);
//
//    func expo(const constant_& c);
//
//    func log(const constant_& c);

//
//
//    //template<typename other_type> bexpr operator+(const other_type& c1, const expr& c2){
//    //    bexpr res;
//    //    res._otype = plus_;
//    //    res._lson = copy((constant_*)&c1);
//    //    res._rson =  copy((constant_*)&c2);
//    //    res._to_str = ::to_str(res._lson) + " + " + ::to_str(res._rson);
//    //    return res;
//    //}
//    //
//    //template<typename other_type> bexpr operator+(const expr& c1, const other_type& c2){
//    //    bexpr res;
//    //    res._otype = plus_;
//    //    res._lson = copy((constant_*)&c1);
//    //    res._rson =  copy((constant_*)&c2);
//    //    res._to_str = ::to_str(res._lson) + " + " + ::to_str(res._rson);
//    //    return res;
//    //}
//    //
//    //
//    //template<typename other_type> bexpr operator-(const other_type& c1, const expr& c2){
//    //    bexpr res;
//    //    res._otype = minus_;
//    //    res._lson = copy((constant_*)&c1);
//    //    res._rson =  copy((constant_*)&c2);
//    //    res._to_str = ::to_str(res._lson) + " - " + ::to_str(res._rson);
//    //    return res;
//    //}
//    //
//    //template<typename other_type> bexpr operator-(const expr& c1, const other_type& c2){
//    //    bexpr res;
//    //    res._otype = minus_;
//    //    res._lson = copy((constant_*)&c1);
//    //    res._rson =  copy((constant_*)&c2);
//    //    res._to_str = ::to_str(res._lson) + " - " + ::to_str(res._rson);
//    //    return res;
//    //}
//    //
//    //
//    //template<typename other_type> bexpr operator*(const other_type& c1, const expr& c2){
//    //    bexpr res;
//    //    res._otype = product_;
//    //    res._lson = copy((constant_*)&c1);
//    //    res._rson =  copy((constant_*)&c2);
//    //    res._to_str = ::to_str(res._lson) + " * " + ::to_str(res._rson);
//    //    return res;
//    //}
//    //
//    //template<typename other_type> bexpr operator*(const expr& c1, const other_type& c2){
//    //    bexpr res;
//    //    res._otype = product_;
//    //    res._lson = copy((constant_*)&c1);
//    //    res._rson =  copy((constant_*)&c2);
//    //    res._to_str = ::to_str(res._lson) + " * " + ::to_str(res._rson);
//    //    return res;
//    //}
//    //
//    //
//    //
//    //template<typename other_type> bexpr operator/(const other_type& c1, const expr& c2){
//    //    bexpr res;
//    //    res._otype = div_;
//    //    res._lson = copy((constant_*)&c1);
//    //    res._rson =  copy((constant_*)&c2);
//    //    res._to_str = ::to_str(res._lson) + " / " + ::to_str(res._rson);
//    //    return res;
//    //}
//    //template<typename other_type> bexpr operator/(const expr& c1, const other_type& c2){
//    //    bexpr res;
//    //    res._otype = div_;
//    //    res._lson = copy((constant_*)&c1);
//    //    res._rson =  copy((constant_*)&c2);
//    //    res._to_str = ::to_str(res._lson) + " / " + ::to_str(res._rson);
//    //    return res;
//    //}
//
//    template<typename type>
//    func_ power(const param<type>& v, unsigned p);
//    
//    func_ power(const func_& f, unsigned p);
//
//    template<typename type>
//    func_ sum(const param<type>& p);
//    
//    template<typename type1, typename type2>
//    func_ product(const param<type1>& p, const param<type2>& v);
//    
//    template<typename type>
//    func_ product(const param<type>& p1, const func_& f);
//    
//    func_ product(const func_& f1, const func_& f2);
//
//    template<typename type>
//    func_ innerproduct(const param<type>& p1, const param<type>& p2);
//    
//    func_ get_poly_derivative(constant_* c, const param_ &v); /*< Get the derivative of c with respect to v) */
//    
//    func_ conj(const func_& f);
//    func_ ang(const func_& f);
//    func_ sqrmag(const func_& f);
//    func_ real(const func_& f);
//    func_ imag(const func_& f);
}



#endif /* func_h */

