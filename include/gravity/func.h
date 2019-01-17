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
    class func_ : public constant_{
    private:
        shared_ptr<func_> compute_derivative(const param_& v);  /**< Computes and stores the derivative of f with respect to variable v. Returns a pointer to the stored function. */
        
    public:
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
        string                                                            _to_str = "";/**< A string representation of the expression */

        size_t                                                            _nb_vars = 0; /**< Number of variables */
        
        size_t                                                            _nnz_j = 0; /**< Number of nonzeros in the Jacobian **/
        size_t                                                            _nnz_h = 0; /**< Number of nonzeros in the Hessian **/
        
        string                                         _name = "noname";
        shared_ptr<indices>                            _indices = nullptr; /*< If indexed, point to the indexing set */
        /** Accessors */
        FType get_ftype() const;
        NType get_return_type() const;

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
        bool is_convex() const;
        bool is_concave() const;
        bool is_convex(size_t idx) const;
        bool is_concave(size_t idx) const;
        bool is_linear() const;
        bool is_quadratic() const;
        bool is_polynomial() const;
        bool is_nonlinear() const;
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
        
        void update_vars(){merge_vars(*this);};
        
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
        
        
        
        bool insert(bool sign, const constant_& coef, const param_& p);/**< Adds coef*p to the function. Returns true if added new term, false if only updated coef of p */
        bool insert(bool sign, const constant_& coef, const param_& p1, const param_& p2, bool c_p1_transposed=false);/**< Adds coef*p1*p2 to the function. Returns true if added new term, false if only updated coef of p1*p2 */
        bool insert(bool sign, const constant_& coef, const list<pair<shared_ptr<param_>, int>>& l);/**< Adds polynomial term to the function. Returns true if added new term, false if only updated corresponding coef */
        
        void insert(const lterm& term);

        void insert(const qterm& term);

        void insert(const pterm& term);

        void update_sign(const constant_& c);
        
        void update_convexity();
        
        void update_sign();
        
        string to_str() const {
            string str;
            int ind = 0;
            string sign = " + ";
            for (auto &pair:*_pterms) {
                str += pair.second.to_str();
            }
            if (!_pterms->empty() && (!_qterms->empty() || !_lterms->empty())) {
                str += " + ";
            }
            ind = 0;
            for (auto &pair:*_qterms) {
                str += pair.second.to_str();
            }
            if (!_qterms->empty() && !_lterms->empty()) {
                str += " + ";
            }
            ind = 0;
            for (auto &pair:*_lterms) {
                str += pair.second.to_str();
            }
            if(!_cst->is_zero()){
                if (_cst->is_number()) {
                    auto val = _cst->to_str();
                    if (val.front()=='-') {
                        str += " - " + val.substr(1);
                    }
                    else if (val != "0"){
                        if (!_pterms->empty() || !_qterms->empty() || !_lterms->empty()) {
                            str += " + ";
                        }
                        str += val;
                    }
                }
                else {
                    if (!_pterms->empty() || !_qterms->empty() || !_lterms->empty()) {
                        str += " + ";
                    }
                    str += "(";
                    str += _cst->to_str();
                    str += ")";
                }
            }
            if (_expr && (!_pterms->empty() || !_qterms->empty() || !_lterms->empty() || !_cst->is_zero())) {
                str += " + ";
            }
            if (_expr) {
                str += _expr->to_str();
            }
            if (_is_vector && (is_number() || _vars->size()>1 || _params->size()>1)) {
                str = "[" + str +"]";
            }
            if (_is_transposed) {
                str += "\u1D40";
            }
            if(str.size()==0){
                str = "0";
            }
            if (str.size() > 2 && str.at(1)=='+') {
                str = str.substr(2);
            }
            return str;
        }
        
        void print_symbolic(bool endline = true, bool display_input = true);
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
    
//
//        void update_sign(const lterm& l);
//        
//        void update_sign(const qterm& q);
//        
//        void update_sign(const pterm& q);
//        
//        void update_convexity(const qterm& q);
//        
    
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
    
//        void print(size_t index);
//        void print();
//    };
//
    template<typename type = double>
    class func: public func_{
    public:
        shared_ptr<vector<type>>                _val = nullptr; /**< vector of values **/
        shared_ptr<pair<type,type>>             _range = nullptr; /**< (Min,Max) values in vals **/
        shared_ptr<vector<pair<type,type>>>     _all_range = nullptr; /**< Vector of (Min,Max) values for each instance of this func **/

        template<typename T=type,
        typename std::enable_if<is_arithmetic<T>::value>::type* = nullptr>
        void update_range(){
            _range = make_shared<pair<type,type>>(make_pair<>(numeric_limits<type>::max(), numeric_limits<type>::lowest()));
        }
        
        template<class T=type, class = typename enable_if<is_same<T, Cpx>::value>::type>
        void update_range() {
            _range = make_shared<pair<type,type>>(make_pair<>(Cpx(numeric_limits<double>::max(), numeric_limits<double>::max()), Cpx(numeric_limits<double>::lowest(), numeric_limits<double>::lowest())));
        }
        
        void update_range(type t) {
            _range->first = t;
            _range->second = t;
        }
        
        func(){
            update_type();
            update_range();
            constant_::set_type(func_c);
            _cst = make_shared<constant<type>>();
            _lterms = make_shared<map<string, lterm>>();
            _qterms = make_shared<map<string, qterm>>();
            _pterms = make_shared<map<string, pterm>>();
            _vars = make_shared<map<string, pair<shared_ptr<param_>, unsigned>>>();
            _params = make_shared<map<string, pair<shared_ptr<param_>, unsigned>>>();
            _val = make_shared<vector<type>>();
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
        void update_type() {
            _type = par_c;
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
            auto ids = indices(ids1,args...);
            if(!_indices || _indices->empty()){/**< No need to add each key individually */
                if(!ids._excluded_keys.empty()){
                    ids.remove_excluded();
                }
                _indices = make_shared<indices>(ids);
                auto dim = _indices->size();
                if(ids._type==matrix_){
                    if(_is_transposed){
                        _dim[0] = ids._dim->at(1);
                        _dim[1] = ids._dim->at(0);
                    }
                    else {
                        _dim[1] = ids._dim->at(0);
                        _dim[0] = ids._dim->at(1);
                    }
                }
                else {
                    if(_is_transposed){
                        _dim[1] = dim;
                    }
                    else {
                        _dim[0] = dim;
                    }
                }
                _name += ".in("+ids._name+")";
            }
            else { /**< Add each key in ids individually */
                _indices->_ids = make_shared<vector<vector<size_t>>>();
                _indices->_ids->resize(1);
                if(ids.empty()){
                    DebugOn("In function param.in(const indices& index_set1, Args&&... args), all index sets are empty!\n. Creating and empty variable! Check your sum/product operators.\n");
                    _name += "_EMPTY";
                    return;
                }
                string key, excluded;
                size_t idx = 0;
                /* Used for truncating extra indices */
                auto nb_sep1 = _indices->_dim->size();
                auto nb_sep2 = ids._dim->size();
                
                for(auto key: *ids._keys){
                    if(ids._excluded_keys.count(idx++)!=0){
                        excluded += key + ",";
                        continue;
                    }
                    if(_indices->_type==to_){
                        key = key.substr(key.find_last_of(",")+1,key.size());
                    }
                    else if(_indices->_type==from_){
                        key = key.substr(0, key.find_last_of(","));
                        key = key.substr(key.find_last_of(",")+1,key.size());
                    }
                    if(nb_sep2>nb_sep1){
                        auto pos = nthOccurrence(key, ",", nb_sep2-nb_sep1);
                        key = key.substr(pos+1,key.size()-1);
                    }
                    auto it1 = _indices->_keys_map->find(key);
                    if (it1 == _indices->_keys_map->end()){
                        throw invalid_argument("In function param.in(const vector<Tobj>& vec), vec has unknown key");
                    }
                    _indices->_ids->at(0).push_back(it1->second);
                }
                if(_is_transposed){
                    _dim[1]=_indices->_ids->at(0).size();
                }
                else {
                    _dim[0]=_indices->_ids->at(0).size();
                }
                _name += ".in("+ids._name+")";
                if(!excluded.empty()){
                    excluded = excluded.substr(0,excluded.size()-1); /* remove last comma */
                    _name += "\{" + excluded + "}";
                }
            }
        }
        
        void print(size_t index, int prec = 10) const {
            cout << to_str(index,prec);
        }
        
        void print(size_t i, size_t j, int prec = 10) const {
            cout << to_str(i,j,prec);
        }
        
        string to_str() const{
            return func_::to_str();
        }
        
        string to_str(size_t index, int prec) const {
            if (is_constant()) {
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
            if (_expr) {
                str += " + ";
                str += _expr->to_str(index, prec);
            }
            if (str.size() > 2 && str.at(1)=='+') {
                str = str.substr(2);
            }
            if (_is_vector && (is_number() || _vars->size()>1 || _params->size()>1)) {
                str = "[" + str +"]";
            }
            if (_is_transposed && (is_number() || _vars->size()>1 || _params->size()>1)) {
                str += "\u1D40";
            }
            return str;
        }
        
        void print(){
            string str;
            if (is_constant()) {
                cout << " (Constant";
            }
            else if (is_linear()) {
                cout << " (Linear";
            }
            else if (is_convex()) {
                cout << " (Convex";
            }
            else if (is_concave()){
                cout << " (Concave";
            }
            else {
                cout << " (Unknown";
            }
            if (is_complex()) {
                cout << " Complex) : ";
            }
            else {
                cout << ") : ";
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
            auto nb_inst = _dim[0];
            allocate_mem();
            for (size_t inst = 0; inst<nb_inst; inst++) {
                eval(inst);
                if (inst>0) {
                    str.insert(str.end(), space_size, ' ');
                }
                str += to_str(inst,3);
                str += "\n";
            }
            str += "\n";
            cout << str;
        }
        string to_str(size_t index1, size_t index2, int prec) const {
            return string();
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
//                        _val = make_shared<vector<bool>>();
//                        _val->push_back(val);
//                        _evaluated = true;
//                        break;
//                    }
//                    case short_c: {
//                        _cst = new constant<short>(move(*(constant<short>*)(&c)));
//                        _all_sign = ((constant<short>*)_cst)->get_sign();
//                        auto val = ((constant<short>*)(&c))->eval();
//                        _all_range = new pair<short,short>(val,val);
//                        _val = make_shared<vector<short>>();
//                        _val->push_back(val);
//                        _evaluated = true;
//                        break;
//                    }
//                    case integer_c: {
//                        _cst = new constant<int>(move(*(constant<int>*)(&c)));
//                        _all_sign = ((constant<int>*)_cst)->get_sign();
//                        auto val = ((constant<int>*)(&c))->eval();
//                        _all_range = new pair<int,int>(val,val);
//                        _val = make_shared<vector<int>>();
//                        _val->push_back(((constant<int>*)(&c))->eval());
//                        _evaluated = true;
//                        break;
//                    }
//                    case float_c: {
//                        _cst = new constant<float>(move(*(constant<float>*)(&c)));
//                        _all_sign = ((constant<float>*)_cst)->get_sign();
//                        auto val = ((constant<float>*)(&c))->eval();
//                        _all_range = new pair<float,float>(val,val);
//                        _val = make_shared<vector<float>>();
//                        _val->push_back(val);
//                        _evaluated = true;
//                        break;
//                    }
//                    case double_c: {
//                        _cst = new constant<double>(move(*(constant<double>*)(&c)));
//                        _all_sign = ((constant<double>*)_cst)->get_sign();
//                        auto val = ((constant<double>*)(&c))->eval();
//                        _all_range = new pair<double,double>(val,val);
//                        _val = make_shared<vector<double>>();
//                        _val->push_back(val);
//                        _evaluated = true;
//                        break;
//                    }
//                    case long_c: {
//                        _cst = new constant<long double>(move(*(constant<long double>*)(&c)));
//                        _all_sign = ((constant<long double>*)_cst)->get_sign();
//                        auto val = ((constant<long double>*)(&c))->eval();
//                        _all_range = new pair<long double,long double>(val,val);
//                        _val = make_shared<vector<long double>>();
//                        _val->push_back(((constant<long double>*)(&c))->eval());
//                        _evaluated = true;
//                        break;
//                    }
//                    case complex_c: {
//                        _cst = new constant<Cpx>(*(constant<Cpx>*)(&c));
//                        _all_sign = g
//                        _all_range = new pair<double,double>(numeric_limits<double>::lowest(),numeric_limits<double>::max());
//                        auto val = ((constant<Cpx>*)(&c))->eval();
//                        _val = make_shared<vector<Cpx>>();
//                        _val->push_back(val);
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
        template<class T2, typename std::enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        shared_ptr<constant_> prod_coef(shared_ptr<constant_> coef, const func<T2>& f){
            if (coef->is_function()) {
                auto f_cst = *static_pointer_cast<func<type>>(coef);
                f_cst *= func<type>(f);
                embed(f_cst);
                return make_shared<constant_>(move(f_cst));
            }
            else if(coef->is_param()) {
                auto p_cst = *static_pointer_cast<param<type>>(coef);
                auto new_cst = f * p_cst;
                return make_shared<func<type>>(move(new_cst));
            }
            else if(coef->is_number()) {
                auto p_cst = func<type>(*static_pointer_cast<constant<type>>(coef));
                return make_shared<func<type>>(p_cst *= f);
            }
        }
        
        template<class T2, typename std::enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        void prod_cst(const constant<T2>& f){
            if (_cst->is_function()) {
                auto f_cst = *static_pointer_cast<func<type>>(_cst);
                f_cst *= func<type>(f);
                *_cst = f_cst;
            }
            else if(_cst->is_param()) {
                auto p_cst = *static_pointer_cast<param<type>>(_cst);
                auto new_cst = f * p_cst;
                _cst = make_shared<func<type>>(move(new_cst));
            }
            else if(_cst->is_number()) {
                auto p_cst = *static_pointer_cast<constant<type>>(_cst);
                auto new_cst = f * p_cst;
                *_cst = new_cst;
            }
        }
        
        template<class T2, typename std::enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        void prod_cst(const param<T2>& f){
            if (_cst->is_function()) {
                auto f_cst = *static_pointer_cast<func<type>>(_cst);
                f_cst *= func<T2>(f);
                embed(f_cst);
            }
            else if(_cst->is_param()) {
                auto p_cst = *static_pointer_cast<param<type>>(_cst);
                auto new_cst = f * p_cst;
                _cst = make_shared<func<type>>(move(new_cst));
            }
            else if(_cst->is_number()) {
                auto p_cst = *static_pointer_cast<constant<type>>(_cst);
                auto new_cst = f * p_cst;
                _cst = make_shared<func<type>>(move(new_cst));
            }
        }
        
        template<class T2, typename std::enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        void prod_cst(const func<T2>& f){
            if (_cst->is_function()) {
                auto f_cst = *static_pointer_cast<func<type>>(_cst);
                f_cst *= f;
                embed(f_cst);
            }
            else if(_cst->is_param()) {
                auto p_cst = func<type>(*static_pointer_cast<param<type>>(_cst));
                _cst = make_shared<func<type>>(p_cst *= f);
            }
            else if(_cst->is_number()) {
                auto p_cst = func<type>(*static_pointer_cast<constant<type>>(_cst));
                _cst = make_shared<func<type>>(p_cst *= f);
            }
        }
        
        template<class T2, typename std::enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        void add_cst(const constant<T2>& f){
            if (_cst->is_function()) {
                auto f_cst = *static_pointer_cast<func<type>>(_cst);
                f_cst += func<type>(f);
                embed(f_cst);
            }
            else if(_cst->is_param()) {
                auto p_cst = *static_pointer_cast<param<type>>(_cst);
                auto new_cst = f + p_cst;
                _cst = make_shared<func<type>>(move(new_cst));
            }
            else if(_cst->is_number()) {
                auto p_cst = *static_pointer_cast<constant<type>>(_cst);
                auto new_cst = f + p_cst;
                *_cst = new_cst;
            }
        }
        
        template<class T2, typename std::enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        void add_cst(const param<T2>& f){
            if (_cst->is_function()) {
                auto f_cst = *static_pointer_cast<func<type>>(_cst);
                func<type>(f_cst) += func<T2>(f);
                embed(f_cst);
            }
            else if(_cst->is_param()) {
                auto p_cst = *static_pointer_cast<param<type>>(_cst);
                auto new_cst = f + p_cst;
                _cst = make_shared<func<type>>(move(new_cst));
            }
            else if(_cst->is_number()) {
                auto p_cst = *static_pointer_cast<constant<type>>(_cst);
                auto new_cst = f + p_cst;
                _cst = make_shared<func<type>>(move(new_cst));
            }
        }
        
        template<class T2, typename std::enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        void add_cst(const func<T2>& f){
            if (_cst->is_function()) {
                auto f_cst = *static_pointer_cast<func<type>>(_cst);
                f_cst += f;
                embed(f_cst);
            }
            else if(_cst->is_param()) {
                auto p_cst = *static_pointer_cast<param<type>>(_cst);
                auto new_cst = f + func<type>(p_cst);
                _cst = make_shared<func<type>>(move(new_cst));
            }
            else if(_cst->is_number()) {
                auto p_cst = *static_pointer_cast<constant<type>>(_cst);
                auto new_cst = f + func<type>(p_cst);
                _cst = make_shared<func<type>>(move(new_cst));
            }
        }
        
        
        template<class T2, typename std::enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        func(T2 c):func(){
            *this = constant<T2>(c);
        }
        
        template<class T2, typename std::enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        func(const constant<T2>& c):func(){
            *this = c;
        }
        template<class T2, class = typename std::enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type>
        func(const param<T2>& c):func(){
            *this = c;
        }
        
        template<class T2, typename std::enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        func(const var<T2>& c):func(){
            *this = c;
        }
        
        template<class T2, class = typename std::enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type>
        func(const func<T2>& f): func(){
            *this = f;
        }
        
        template<class T2, typename std::enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        func& operator=(const constant<T2>& c){
            reset();
            static_pointer_cast<constant<type>>(_cst)->set_val(c.eval());
            _all_sign = _cst->get_sign();
            _val->resize(1);
            _val->at(0) = c.eval();
            update_range(_val->at(0));
            _evaluated = true;
            return *this;
        }
        
        template<class T2, typename std::enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        func& operator=(const param<T2>& c){
            reset();
            insert(lterm(c.pcopy()));
            _val->clear();
            _val->resize(c._val->size());
            for(auto i = 0; i< c._val->size(); i++){
                _val->at(i) = c._val->at(i);
            }
            _range->first = c._range->first;
            _range->second = c._range->second;
            _evaluated = true;
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
        
        shared_ptr<constant_> copy()const{return make_shared<func>(*this);};
        
        template<class T2, typename std::enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        void deep_copy(const func<T2>& f){
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
                _expr = make_shared<expr>(*f._expr);
            }
            if(f._indices){
                _indices = make_shared<indices>(*f._indices);
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
        }
        
        func& operator=(const func& f){
            deep_copy(f);
            return  *this;
        }
        
        template<class T2, typename std::enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
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
        
        
        template<class T=type, class = typename enable_if<is_arithmetic<T>::value>::type> void update_range(T val) {
            if (val < _range->first) {
                _range->first = val;
            }
            if (val > _range->second) {
                _range->second = val;
            }
        }
        
        void add_val(size_t i, type val) {
            if(is_matrix()){
                throw invalid_argument("Cannot call func::add_val(type val) on matrix");
            }
            _dim[0] = max(_dim[0],i+1);
            _val->resize(max(_val->size(),i+1));
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
                throw invalid_argument("in Function size_t set_val(const string& key, type val), unknown key");
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
                _val->resize(max(_val->size(),index+1));
                _dim[0] = max(_dim[0],_val->size());
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
            _dim[0] = max(_dim[0],i+1);
            _dim[1] = max(_dim[1],j+1);
            auto index = _dim[1]*i+j;
            _val->resize(max(_val->size(),index+1));
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
        typename std::enable_if<is_arithmetic<T>::value>::type* = nullptr> Sign get_sign(size_t idx) const{
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
        
        template<class T=type, class = typename enable_if<is_same<T, Cpx>::value>::type> Sign get_sign(size_t idx) const{
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
        
        
        //        template<class T2, typename std::enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        //        shared_ptr<constant_> add(shared_ptr<param<T2>> c1){
        //
        //        }
        
        template<class T=type, class = typename enable_if<is_same<T, Cpx>::value>::type> Sign get_all_sign() const{
            if (_range->first == Cpx(0,0) && _range->second == Cpx(0,0)) {
                return zero_;
            }
            if ((_range->second.real() < 0 && _range->second.imag() < 0)) {
                return neg_;
            }
            if ((_range->second.real() > 0 && _range->second.imag() > 0)) {
                return pos_;
            }
            if (_range->second == Cpx(0,0)) {
                return non_pos_;
            }
            if (_range->first == Cpx(0,0)) {
                return non_neg_;
            }
            return unknown_;
        }
        
        template<typename T=type,
        typename std::enable_if<is_arithmetic<T>::value>::type* = nullptr> Sign get_all_sign() const {
            if (_range->first == 0 && _range->second == 0) {
                return zero_;
            }
            if (_range->second < 0  && _range->first < 0) {
                return neg_;
            }
            if (_range->first > 0 && _range->second > 0) {
                return pos_;
            }
            if (_range->second == 0   && _range->first < 0) {
                return non_pos_;
            }
            if (_range->first == 0  && _range->second > 0) {
                return non_neg_;
            }
            return unknown_;
        }
        
        type eval() const {
            if (is_indexed()) {
                return _val->at(_indices->_ids->at(0).back());
            }
            return _val->back();
        }
        
        type eval(size_t i) const {
            if(is_matrix()){
                throw invalid_argument("eval() should be called with double index here\n");
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
        
        
        type eval(const string& key) const{
            return _val->at(_indices->_keys_map->at(key));
        }
        
        type eval(size_t i, size_t j) const {
            
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
        
        template<class T=type, typename enable_if<is_arithmetic<T>::value>::type* = nullptr> bool is_unit() const { /**< Returns true if all values of this paramter are 1 **/
            return (_range->first == 1 && _range->second == 1);
        }
        
        template<class T=type, class = typename enable_if<is_same<T, Cpx>::value>::type> bool is_unit() const{
            return (_range->first == Cpx(1,1) && _range->second == Cpx(1,1));
        }
        
        bool is_zero() const { return zero_range();};
        
        template<class T=type, class = typename enable_if<is_same<T, Cpx>::value>::type> bool zero_range() const{
            return (_range->first == Cpx(0,0) && _range->second == Cpx(0,0));
        }
        
        template<typename T=type,
        typename std::enable_if<is_arithmetic<T>::value>::type* = nullptr> bool zero_range() const{
            return (_range->first == 0 && _range->second == 0);
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
        
        template<class T2, typename std::enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        func& operator*=(const func<T2>& c){
            if (is_zero()) {
                return *this;
            }
            if (c.is_zero()) {
                reset();
                return *this;
            }
            if (is_unit()) {
                *this = func(c);
                return *this;
            }
            if (c.is_unit()) {
                return *this;
            }
    
            /* Case where c is a number */
//            if (c.is_number()){
//                return *this *= constant<T2>(c.eval());
//            }
            /* Case where the current function is not constant and the other operand is */
            if(!is_constant() && c.is_constant()) {
                bool transp = false;
                auto fc = c;
                if(is_linear() && _is_transposed){// Situation where f^T * c is transformed into (c^T*f)^T
                    fc.transpose();
                    this->transpose();
                    transp = true;
                }
                if (!_cst->is_zero()) {
                    fc.prod_cst(c);
                }
                for (auto &pair:*_lterms) {
                    pair.second._coef = prod_coef(pair.second._coef, fc);
                }
                for (auto &pair:*_qterms) {
                    pair.second._coef = prod_coef(pair.second._coef, fc);
                }
                for (auto &pair:*_pterms) {
                    pair.second._coef = prod_coef(pair.second._coef, fc);
                }
                if (_expr) {
                    _expr = make_shared<bexpr>(bexpr(product_, make_shared<expr>((*_expr)), make_shared<func<T2>>(c)));
                    embed(_expr);
                }
                if (c.is_negative()) {
                    reverse_sign();
                }
                if (c.get_all_sign()==unknown_) {
                    _all_sign = unknown_;
                    if (!_qterms->empty()) {
                        _all_convexity = undet_;
                    }
                }
                _evaluated = false;
                if(transp){
                    this->transpose();
                }
                update_dot_dim(c);
                return *this;
            }
            /* Case where the current function is constant and the other operand is not. */
//            if (is_constant() && (c.is_var() || (c.is_function() && !((func_*)&c)->is_constant()))) {
//                func_ f(c);
//                if (!f._cst->is_zero()) {
//                    if (f._cst->is_function()) {
//                        auto fc = (func_*)f._cst;
//                        *fc = (*this)* (*fc);
//                        f.embed(*fc);
//                    }
//                    else {
//                        f._cst = multiply(f._cst, *this);
//                    }
//                }
//                for (auto &pair:*f._lterms) {
//                    if (pair.second._coef->is_function()) {
//                        auto fc = (func_*)pair.second._coef;
//                        *fc = (*this)* (*fc);
//                        f.embed(*fc);
//                    }
//                    else {
//                        pair.second._coef = multiply(pair.second._coef, *this);
//                    }
//                    if (pair.second._coef->_is_transposed) {
//                        pair.second._p->_is_vector = true;
//                        if (!pair.second._coef->is_number() && pair.second._coef->_dim[1]!=pair.second._p->_dim[0]) {
//                            DebugOn("vector dot product with mismatching dimensions, check your param/var dimensions");
//                        }
//                    }
//                    //                f.update_nb_instances(pair.second);
//
//                }
//                for (auto &pair:*f._qterms) {
//                    if (pair.second._coef->is_function()) {
//                        auto fc = (func_*)pair.second._coef;
//                        *fc = (*this)* (*fc);
//                        f.embed(*fc);
//                    }
//                    else {
//                        pair.second._coef = multiply(pair.second._coef, *this);
//                    }
//                    if (pair.second._coef->_is_transposed) {
//                        pair.second._p->first->_is_vector = true;
//                        pair.second._p->second->_is_vector = true;
//                        if (!pair.second._coef->is_number() && pair.second._coef->_dim[1]!=pair.second._p->first->_dim[0]) {
//                            DebugOn("vector dot product with mismatching dimensions, check your param/var dimensions");
//                        }
//
//                    }
//                    //                update_nb_instances(pair.second);
//                }
//                for (auto &pair:*f._pterms) {
//                    if (pair.second._coef->is_function()) {
//                        auto fc = (func_*)pair.second._coef;
//                        *fc = (*this)* (*fc);
//                        f.embed(*fc);
//                    }
//                    else {
//                        pair.second._coef = multiply(pair.second._coef, *this);
//                    }
//                    //                update_nb_instances(pair.second); // TODO
//                }
//                if (f._expr) {
//                    if (this->is_number()) {
//                        f._expr->_coef *= t_eval(this);
//                    }
//                    else {
//                        f._expr = make_shared<bexpr>(bexpr(product_, make_shared<func_>(*this), make_shared<func_>((*f._expr))));
//                        f.embed(f._expr);
//                        //                _DAG->insert(make_pair<>(_expr->get_str(), _expr));
//                        f._queue->push_back(f._expr);
//                    }
//                }
//                if (this->is_negative()) {
//                    f.reverse_sign();
//                }
//                if (this->get_all_sign()==unknown_) {
//                    f._all_sign = unknown_;
//                    if (!f._qterms->empty()) {
//                        f._all_convexity = undet_;
//                    }
//                }
//                f._evaluated = false;
//                if(update_dot(c)){
//                    f._dim[0] = _dim[0];
//                    f._dim[1] = _dim[1];
//                }
//                return *this = move(f);
//            }
//            if (c.is_param() || c.is_var()) {
//    //            if (c.is_matrix() || is_matrix()) {
//    //                *this = func_(bexpr(product_, make_shared<func_>(*this), make_shared<func_>(c)));
//    //            }
//    //            else {
//                if(_is_transposed && !c._is_vector){
//                    auto new_c = copy(c);
//                    auto new_p = (param_*)new_c;
//                    new_p->_is_vector = true;
//                    new_p->_name = "["+new_p->_name+"]";
//                    func_ f(*new_p);
//                    *this *= f;
//                    delete new_c;
//                }
//                else {
//                    func_ f(c);
//                    *this *= f;
//                }
//    //            }
//                _evaluated = false;
//                return *this;
//            }
//            if (_expr || (c.is_function() && ((func_*)&c)->_expr)) {
//                auto be = bexpr(product_, make_shared<func_>(*this), make_shared<func_>(c));
//                *this = func_(be);
//                _evaluated = false;
//                return *this;
//            }
//            /* Case where the multiplication invlolves multiplying variables/parameters together, i.e., they are both parametric or both include variables  */
//            if (c.is_function()) {
//                func_* f = (func_*)&c;
//                constant_* coef;
//                vector<bool>* is_sum = nullptr;
//                func_ res;
//                for (auto& t1: *_pterms) {
//                    if (t1.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial function), we cannot factor the coefficients. Just create a binary expression and return it.
//                        auto be = bexpr(product_, make_shared<func_>(*this), make_shared<func_>(c));
//                        *this = func_(be);
//                        _evaluated = false;
//                        return *this;
//                    }
//                    for (auto& t2: *f->_pterms) {
//                        is_sum = nullptr;
//                        if (t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial function), see comment above.
//                            auto be = bexpr(product_, make_shared<func_>(*this), make_shared<func_>(c));
//                            *this = func_(be);
//                            _evaluated = false;
//                            return *this;
//                        }
//                        auto newl(*t1.second._l);
//                        for (auto& it: *t2.second._l) {
//                            newl.push_back(make_pair<>(it.first, it.second));
//                        }
//                        coef = copy(*t1.second._coef);
//                        coef = multiply(coef, *t2.second._coef);
//                        res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
//                        delete coef;
//                    }
//                    for (auto& t2: *f->_qterms) {
//                        if (t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial function), see comment above.
//                            auto be = bexpr(product_, make_shared<func_>(*this), make_shared<func_>(c));
//                            *this = func_(be);
//                            _evaluated = false;
//                            return *this;
//                        }
//                        auto newl(*t1.second._l);
//                        newl.push_back(make_pair<>((t2.second._p->first), 1));
//                        newl.push_back(make_pair<>((t2.second._p->second), 1));
//                        coef = copy(*t1.second._coef);
//                        coef = multiply(coef, *t2.second._coef);
//                        res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
//                        delete coef;
//                    }
//                    for (auto& t2: *f->_lterms) {
//                        if (t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial function)
//                            auto be = bexpr(product_, make_shared<func_>(*this), make_shared<func_>(c));
//                            *this = func_(be);
//                            _evaluated = false;
//                            return *this;
//                        }
//                        auto newl(*t1.second._l);
//                        newl.push_back(make_pair<>((t2.second._p), 1));
//                        coef = copy(*t1.second._coef);
//                        coef = multiply(coef, *t2.second._coef);
//                        res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
//                        delete coef;
//                    }
//                    if (!f->_cst->is_zero()) {
//                        auto newl(*t1.second._l);
//                        coef = copy(*f->_cst);
//                        coef = multiply(coef, *t1.second._coef);
//                        res.insert(t1.second._sign, *coef, newl);
//                        delete coef;
//                    }
//                }
//
//                for (auto& t1: *_qterms) {
//                    if (t1.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(Quadratic term)
//                        auto be = bexpr(product_, make_shared<func_>(*this), make_shared<func_>(c));
//                        *this = func_(be);
//                        _evaluated = false;
//                        return *this;
//                    }
//                    for (auto& t2: *f->_pterms) {
//                        is_sum = nullptr;
//                        if (t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial term)
//                            auto be = bexpr(product_, make_shared<func_>(*this), make_shared<func_>(c));
//                            *this = func_(be);
//                            _evaluated = false;
//                            return *this;
//                        }
//                        auto newl(*t2.second._l);
//                        newl.push_front(make_pair<>(t1.second._p->first, 1));
//                        newl.push_front(make_pair<>(t1.second._p->second, 1));
//                        coef = copy(*t1.second._coef);
//                        coef = multiply(coef, *t2.second._coef);
//                        res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
//                        delete coef;
//                    }
//                    for (auto& t2: *f->_qterms) {
//                        if (t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial term)
//                            auto be = bexpr(product_, make_shared<func_>(*this), make_shared<func_>(c));
//                            *this = func_(be);
//                            _evaluated = false;
//                            return *this;
//                        }
//                        coef = copy(*t1.second._coef);
//                        coef = multiply(coef, *t2.second._coef);
//                        list<pair<param_*, int>> newl;
//                        newl.push_back(make_pair<>(t1.second._p->first, 1));
//                        newl.push_back(make_pair<>(t1.second._p->second, 1));
//                        newl.push_back(make_pair<>(t2.second._p->first, 1));
//                        newl.push_back(make_pair<>(t2.second._p->second, 1));
//                        res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
//                        delete coef;
//                    }
//                    for (auto& t2: *f->_lterms) {
//                        is_sum = nullptr;
//                        if (t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial term)
//                            auto be = bexpr(product_, make_shared<func_>(*this), make_shared<func_>(c));
//                            *this = func_(be);
//                            _evaluated = false;
//                            return *this;
//                        }
//                        coef = copy(*t1.second._coef);
//                        coef = multiply(coef, *t2.second._coef);
//                        list<pair<param_*, int>> newl;
//                        newl.push_back(make_pair<>(t1.second._p->first, 1));
//                        newl.push_back(make_pair<>(t1.second._p->second, 1));
//                        newl.push_back(make_pair<>(t2.second._p, 1));
//                        res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
//                        delete coef;
//                    }
//                    if (!f->_cst->is_zero()) {
//                        coef = copy(*t1.second._coef);
//                        coef = multiply(coef, *f->_cst);
//                        res.insert(t1.second._sign, *coef, *t1.second._p->first, *t1.second._p->second);
//                        delete coef;
//                    }
//
//                }
//                for (auto& t1: *_lterms) {
//    //                if (t1.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(Quadratic term)
//    //                    auto be = bexpr(product_, make_shared<func_>(*this), make_shared<func_>(c));
//    //                    *this = func_(be);
//    //                    _evaluated = false;
//    //                    return *this;
//    //                }
//                    for (auto& t2: *f->_pterms) {
//                        is_sum = nullptr;
//                        if (t1.second._coef->_is_transposed || t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial term)
//                            auto be = bexpr(product_, make_shared<func_>(*this), make_shared<func_>(c));
//                            *this = func_(be);
//                            _evaluated = false;
//                            return *this;
//                        }
//                        auto newl(*t2.second._l);
//                        newl.push_front(make_pair<>((t1.second._p), 1));
//                        coef = copy(*t1.second._coef);
//                        coef = multiply(coef, *t2.second._coef);
//                        res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
//                        delete coef;
//                    }
//                    for (auto& t2: *f->_qterms) {
//                        if (t1.second._coef->_is_transposed || t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial term)
//                            auto be = bexpr(product_, make_shared<func_>(*this), make_shared<func_>(c));
//                            *this = func_(be);
//                            _evaluated = false;
//                            return *this;
//                        }
//                        coef = copy(*t1.second._coef);
//                        coef = multiply(coef, *t2.second._coef);
//                        list<pair<param_*, int>> newl;
//                        newl.push_back(make_pair<>(t1.second._p, 1));
//                        newl.push_back(make_pair<>(t2.second._p->first, 1));
//                        newl.push_back(make_pair<>(t2.second._p->second, 1));
//                        res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
//                        delete coef;
//                    }
//                    for (auto& t2: *f->_lterms) {
//    //                    if (t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial term)
//    //                        auto be = bexpr(product_, make_shared<func_>(*this), make_shared<func_>(c));
//    //                        *this = func_(be);
//    //                        _evaluated = false;
//    //                        return *this;
//    //                    }
//                        coef = copy(*t1.second._coef);
//                        coef = multiply(coef, *t2.second._coef);
//                        res.insert(!(t1.second._sign^t2.second._sign), *coef, *t1.second._p, *t2.second._p, _is_transposed);
//                        delete coef;
//                    }
//                    if (!f->_cst->is_zero()) {
//                        coef = copy(*t1.second._coef);
//                        coef = multiply(coef, *f->_cst);
//                        res.insert(t1.second._sign, *coef, *t1.second._p);
//                        delete coef;
//                    }
//                }
//                if (!_cst->is_zero()) {
//                    for (auto& t2: *f->_pterms) {
//                        if (t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial term)
//                            auto be = bexpr(product_, make_shared<func_>(*this), make_shared<func_>(c));
//                            *this = func_(be);
//                            _evaluated = false;
//                            return *this;
//                        }
//                        coef = copy(*_cst);
//                        coef = multiply(coef, *t2.second._coef);
//                        res.insert(t2.second._sign, *coef, *t2.second._l);
//                        delete coef;
//                    }
//                    for (auto& t2: *f->_qterms) {
//                        if (t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial term)
//                            auto be = bexpr(product_, make_shared<func_>(*this), make_shared<func_>(c));
//                            *this = func_(be);
//                            _evaluated = false;
//                            return *this;
//                        }
//                        coef = copy(*_cst);
//                        coef = multiply(coef, *t2.second._coef);
//                        res.insert(t2.second._sign, *coef, *t2.second._p->first, *t2.second._p->second);
//                        delete coef;
//                    }
//                    for (auto& t2: *f->_lterms) {
//                        if (t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial term)
//                            auto be = bexpr(product_, make_shared<func_>(*this), make_shared<func_>(c));
//                            *this = func_(be);
//                            _evaluated = false;
//                            return *this;
//                        }
//                        coef = copy(*_cst);
//                        coef = multiply(coef, *t2.second._coef);
//                        res.insert(t2.second._sign, *coef, *t2.second._p);
//                        delete coef;
//                    }
//                    if (!f->_cst->is_zero()) {
//                        coef = copy(*_cst);
//                        coef = multiply(coef, *f->_cst);
//                        delete res._cst;
//                        res._cst = coef;
//                        if (_cst->is_function()) {
//                            embed(*(func_*)_cst);
//                        }
//                    }
//                }
//                res.update_dot_dim(*this, c);
//                *this = move(res);
//            }
            _evaluated = false;
            return *this;
        }
        
        template<class T2, typename std::enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        func& operator+=(const func<T2>& f){
            if (f.is_zero()) {
                return *this;
            }
            _evaluated = false;
            if (is_constant() && !f.is_constant()) {
                func res(f);
                res += *this;
                return *this = res;
            }
            if (!is_constant() && f.is_constant()) {
                this->add_cst(f);
                update_sign(f);
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
                    embed(*static_pointer_cast<func_>(_cst));
                }
            }
            for (auto &pair:f.get_lterms()) {
                this->insert(pair.second);
            }
            for (auto &pair:f.get_qterms()) {
                this->insert(pair.second);
            }
            for (auto &pair:f.get_pterms()) {
                this->insert(pair.second);
            }
            if (_expr && f.get_expr()) {
                shared_ptr<expr> e1,e2;
                if (_expr->is_uexpr()) {
                    e1 = make_shared<uexpr>(*static_pointer_cast<uexpr>(_expr));
                }
                else {
                    e1 = make_shared<bexpr>(*static_pointer_cast<bexpr>(_expr));
                }
                if (f.get_expr()->is_uexpr()) {
                    e2 = make_shared<uexpr>(*static_pointer_cast<uexpr>(f.get_expr()));
                }
                else {
                    e2 = make_shared<bexpr>(*static_pointer_cast<bexpr>(f.get_expr()));
                }
                _expr = make_shared<bexpr>(bexpr(plus_, e1, e2));
                embed(_expr);
            }
            else if (!_expr && f.get_expr()) {
                if (f.get_expr()->is_uexpr()) {
                    _expr = make_shared<uexpr>(*static_pointer_cast<uexpr>(f.get_expr()));
                }
                else {
                    _expr = make_shared<bexpr>(*static_pointer_cast<bexpr>(f.get_expr()));
                }
                embed(_expr);
                if (!_vars->empty()) {
                    _ftype = nlin_;
                }
            }
            update_sign(f);
            update_convexity();
            return *this;
        }
        
        template<class T2, typename std::enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        func& operator-=(const func<T2>& f){
            return *this += -1*f;
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
            *_cst = constant<type>();
            _nb_vars = 0;
            _nnz_h = 0;
            _nnz_j = 0;
        };

        
    };
    
    
    template<class T1,class T2, typename std::enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator+(const func<T1>& f1, const func<T2>& f2){
        return func<T1>(f1)+= f2;
    }

    template<class T1,class T2, typename std::enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator+(const func<T1>& f1, const func<T2>& f2){
        return func<T2>(f1)+= f2;
    }
//
//    template<class T1,class T2, typename std::enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
//    func<T1> operator-(const func<T1>& f1, const func<T2>& f2){
//        return func<T1>(f1)+= f2;
//    }
//
//    template<class T1,class T2, typename std::enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
//    func<T2> operator-(const func<T1>& f1, const func<T2>& f2){
//        return func<T2>(f1)+= f2;
//    }
//
//    template<class T1,class T2, typename std::enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
//    func<T1> operator*(const func<T1>& f1, const func<T2>& f2){
//        return func<T1>(f1)+= f2;
//    }
//
//    template<class T1,class T2, typename std::enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
//    func<T2> operator*(const func<T1>& f1, const func<T2>& f2){
//        return func<T2>(f1)+= f2;
//    }
    
    
    template<class T1,class T2, typename std::enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator*(const param<T1>& p1, const param<T2>& p2){
        if(p1.is_function() || p2.is_function()){
            
        }
        func<T1> res;
        if(p1.is_param() && p2.is_var()){
            res.insert(true,p1,p2);
        }
        else if(p2.is_param() && p1.is_var()){
            res.insert(true,p2,p1);
        }
        else {//Both vars or both params
            res.insert(true,constant<double>(1),p1,p2);
        }
        res.update_dot_dim(p1,p2);
        return res;
    }

    template<class T1,class T2, typename std::enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator*(const param<T1>& p1, const param<T2>& p2){
        func<T2> res;
        if(p1.is_param() && p2.is_var()){
            res.insert(true,p1,p2);
        }
        else if(p2.is_param() && p1.is_var()){
            res.insert(true,p2,p1);
        }
        else {//Both vars or both params
            res.insert(true,constant<double>(1),p1,p2);
        }
        res.update_dot_dim(p1,p2);
        return res;
    }

    template<class T1,class T2, typename std::enable_if<is_convertible<T2, T1>::value && sizeof(T2) <= sizeof(T1)>::type* = nullptr>
    func<T1> operator+(const param<T1>& p1, const param<T2>& p2){
        func<T1> res;
        auto newp1 = p1.pcopy();
        auto newp2 = p2.pcopy();
        res.update_dim(p1,p2);
        if(p1.is_param() && p2.is_var()){
            res.insert(lterm(newp2));
            res.add_cst(p1);
        }
        else if(p2.is_param() && p1.is_var()){
            res.insert(lterm(newp1));
            res.add_cst(p2);
        }
        else {//Both vars or both params
            res.insert(lterm(newp1));
            res.insert(lterm(newp2));
        }
        return res;
    }
        
    template<class T1,class T2, typename std::enable_if<is_convertible<T1, T2>::value && sizeof(T1) < sizeof(T2)>::type* = nullptr>
    func<T2> operator+(const param<T1>& p1, const param<T2>& p2){
        if(p1.is_function() || p2.is_function()){
            return func<T2>(p2) += func<T1>(p1);
        }
        func<T2> res;
        auto newp1 = p1.pcopy();
        auto newp2 = p2.pcopy();
        res.update_dim(p1,p2);
        if(p1.is_param() && p2.is_var()){
            res.insert(lterm(newp2));
            res.add_cst(p1);
        }
        else if(p2.is_param() && p1.is_var()){
            res.insert(lterm(newp1));
            res.add_cst(p2);
        }
        else {//Both vars or both params
            res.insert(lterm(newp1));
            res.insert(lterm(newp2));
        }
        return res;
    }
    
    template<class T1,class T2, typename std::enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator-(const param<T1>& p1, const param<T2>& p2){
        if(p1.is_function() || p2.is_function()){
            return func<T1>(p1) -= p2;
        }
        func<T1> res;
        auto newp1 = p1.pcopy();
        auto newp2 = p2.pcopy();
        res.update_dim(p1,p2);
        if(p1.is_param() && p2.is_var()){
            res.insert(lterm(false,newp2));
            res.add_cst(p1);
        }
        else if(p2.is_param() && p1.is_var()){
            res.insert(lterm(newp1));
            res.add_cst(-1*p2);
        }
        else {//Both vars or both params
            res.insert(lterm(newp1));
            res.insert(lterm(false,newp2));
        }
        return res;
    }
    
    template<class T1,class T2, typename std::enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator-(const param<T1>& p1, const param<T2>& p2){
        func<T2> res;
        auto newp1 = p1.pcopy();
        auto newp2 = p2.pcopy();
        res.update_dim(p1,p2);
        if(p1.is_param() && p2.is_var()){
            res.insert(lterm(false,newp2));
            res.add_cst(p1);
        }
        else if(p2.is_param() && p1.is_var()){
            res.insert(lterm(newp1));
            res.add_cst(-1*p2);
        }
        else {//Both vars or both params
            res.insert(lterm(newp1));
            res.insert(lterm(false,newp2));
        }
        return res;
    }
        
    template<class T1,class T2, typename std::enable_if<is_convertible<T2, T1>::value && sizeof(T2) <= sizeof(T1)>::type* = nullptr>
    constant<T1> operator+(const constant<T1>& p1, const constant<T2>& p2){
        constant<T1> res(p1);
        res.set_val(res.eval()+p2.eval());
        return res;
    }
    
    template<class T1,class T2, typename std::enable_if<is_convertible<T1, T2>::value && sizeof(T1) < sizeof(T2)>::type* = nullptr>
    constant<T2> operator+(const constant<T1>& p1, const constant<T2>& p2){
        constant<T2> res(p2);
        res.set_val(res.eval()+(T2)p1.eval());
        return res;
    }
        
    template<class T1,class T2, typename std::enable_if<is_convertible<T2, T1>::value && sizeof(T2) <= sizeof(T1)>::type* = nullptr>
    constant<T1> operator-(const constant<T1>& p1, const constant<T2>& p2){
        constant<T1> res(p1);
        res.set_val(res.eval()-p2.eval());
        return res;
    }
    
    template<class T1,class T2, typename std::enable_if<is_convertible<T1, T2>::value && sizeof(T1) < sizeof(T2)>::type* = nullptr>
    constant<T2> operator-(const constant<T1>& p1, const constant<T2>& p2){
        constant<T2> res(p2);
        res.set_val((T2)p1.eval() - res.eval());
        return res;
    }
        
    template<class T1,class T2, typename std::enable_if<is_convertible<T2, T1>::value && sizeof(T2) <= sizeof(T1)>::type* = nullptr>
    constant<T1> operator*(const constant<T1>& p1, const constant<T2>& p2){
        constant<T1> res(p1);
        res.set_val(res.eval()*p2.eval());
        return res;
    }
    
    template<class T1,class T2, typename std::enable_if<is_convertible<T1, T2>::value && sizeof(T1) < sizeof(T2)>::type* = nullptr>
    constant<T2> operator*(const constant<T1>& p1, const constant<T2>& p2){
        constant<T2> res(p2);
        res.set_val((T2)p1.eval() * res.eval());
        return res;
    }

    template<class T1,class T2, typename std::enable_if<is_convertible<T2, T1>::value && sizeof(T2) <= sizeof(T1)>::type* = nullptr>
    constant<T1> operator/(const constant<T1>& p1, const constant<T2>& p2){
        constant<T1> res(p1);
        res.set_val(res.eval()/p2.eval());
        return res;
    }
    
    template<class T1,class T2, typename std::enable_if<is_convertible<T1, T2>::value && sizeof(T1) < sizeof(T2)>::type* = nullptr>
    constant<T2> operator/(const constant<T1>& p1, const constant<T2>& p2){
        constant<T2> res(p2);
        res.set_val((T2)p1.eval()/res.eval());
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
        
//    template<class T1,class T2, typename std::enable_if<is_convertible<T2, T1>::value && sizeof(T2) <= sizeof(T1)>::type* = nullptr>
//    func<T1> operator+(const var<T1>& p1, const var<T2>& p2){
//        func<T1> res;
//        auto newp1 = p1.pcopy();
//        auto newp2 = p2.pcopy();
//        res.update_dim(p1,p2);
//        res.insert(lterm(newp1));
//        res.insert(lterm(newp2));
//        return res;
//    }
//
//    template<class T1,class T2, typename std::enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
//    func<T1> operator+(const param<T1>& p, const var<T2>& v){
//        func<T1> res;
//        res.add_cst(p);
//        auto newv = v.pcopy();
//        auto newp = p.pcopy();
//        res.update_dim(p,v);
//        res.insert(lterm(newv));
//        return res;
//    }
//
//    template<class T1,class T2, typename std::enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
//    func<T2> operator+(const param<T1>& p, const var<T2>& v){
//        func<T2> res;
//        auto newv = v.pcopy();
//        auto newp = p.pcopy();
//        res.add_param(newp);
//        res.update_dim(p,v);
//        res.insert(lterm(newv));
//        return res;
//    }
//
    template<class T1,class T2, typename std::enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator+(const constant<T1>& p, const param<T2>& v){
        func<T1> res(v);
        res.add_cst(p);
        return res;
    }
    
    template<class T1,class T2, typename std::enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator+(const constant<T1>& p, const param<T2>& v){
        func<T2> res(v);
        res.add_cst(p);
        return res;
    }

        
    template<class T1,class T2, typename std::enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator+(const param<T1>& v, const constant<T2>& p){
        func<T1> res(v);
        res.add_cst(p);
        return res;
    }
    
    template<class T1,class T2, typename std::enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator+(const param<T1>& v, const constant<T2>& p){
        func<T2> res(v);
        res.add_cst(p);
        return res;
    }
    
    template<class T1,class T2, typename std::enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator+(const param<T1>& v, const func<T2>& f){
        func<T1> res(v);
        res += f;
        return res;
    }
    
    template<class T1,class T2, typename std::enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator+(const param<T1>& v, const func<T2>& f){
        func<T2> res(v);
        res += f;
        return res;
    }
    
    template<class T1,class T2, typename std::enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator+(const func<T1>& f, const param<T2>& v){
        func<T1> res(v);
        res += f;
        return res;
    }
    
    template<class T1,class T2, typename std::enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator+(const func<T1>& f, const param<T2>& v){
        func<T2> res(v);
        res += f;
        return res;
    }
    
    template<class T1,class T2, typename std::enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator*(const param<T1>& v, const func<T2>& f){
        func<T1> res(v);
        res *= f;
        return res;
    }
    
    template<class T1,class T2, typename std::enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator*(const param<T1>& v, const func<T2>& f){
        func<T2> res(v);
        res *= f;
        return res;
    }
    
    template<class T1,class T2, typename std::enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator*(const func<T1>& f, const param<T2>& v){
        func<T1> res(v);
        res *= f;
        return res;
    }
    
    template<class T1,class T2, typename std::enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator*(const func<T1>& f, const param<T2>& v){
        func<T2> res(v);
        res *= f;
        return res;
    }
        
    template<class T1,class T2, typename std::enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator*(const constant<T1>& p, const param<T2>& v){
        func<T1> res(v);
        res.add_cst(p);
        return res;
    }
    
    template<class T1,class T2, typename std::enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator*(const constant<T1>& c, const param<T2>& p){
        func<T2> res(p);
        res.add_cst(c);
        return res;
    }
    
    
    template<class T1,class T2, typename std::enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator*(const param<T1>& p, const constant<T2>& c){
        func<T1> res(p);
        res.add_cst(c);
        return res;
    }
        
    template<class T1,class T2, typename std::enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator*(const param<T1>& p, const constant<T2>& c){
        func<T2> res(p);
        res.add_cst(c);
        return res;
    }
    
        
    template<class T1,class T2, typename std::enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator*(T1 p, const param<T2>& v){
        func<T1> res(v);
        res.add_cst(constant<T1>(p));
        return res;
    }
    
    template<class T1,class T2, typename std::enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator*(T1 p, const param<T2>& v){
        func<T2> res(v);
        res.add_cst(constant<T1>(p));
        return res;
    }
    
    
    template<class T1,class T2, typename std::enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator*(const param<T1>& v, T2 p){
        func<T1> res(v);
        res.add_cst(constant<T2>(p));
        return res;
    }
    
    template<class T1,class T2, typename std::enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator*(const param<T1>& v, T2 p){
        func<T2> res(v);
        res.add_cst(p);
        return res;
    }
        // Transform var<T2> to param<T2> and check type.
    
    
    
    
//
//    template<typename T1,typename T2>
//    func<T1> operator+(const func<T1>& c1, const func<T2>& c2){
//        func<T1> res;
//        return res;
//    }
//
//    template<typename T1,typename T2>
//    func<T1> operator+(const param<T1>& c1, const var<T2>& c2){
//        func<T1> res;
//        return res;
//    }
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
        
//    shared_ptr<constant_> add(shared_ptr<constant_> c1, shared_ptr<constant_> c2); /**< adds c2 to c1, and returns the result **/
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
    template<typename T> shared_ptr<constant_> multiply(shared_ptr<constant_> c1, const func<T>& c2){ /**< multiplies c1 with c2, returns the result **/
        switch (c1->get_type()) {
            case binary_c: {
                auto val = (static_pointer_cast<constant<bool>>(c1))->eval();
                if (val==true) {
                    return make_shared<func<T>>(c2);
                }
                break;
            }
            case short_c: {
                auto val = (static_pointer_cast<constant<short>>(c1))->eval();
                if (val==0) {
                    return c1;
                }
                if (val==1) {
                    return make_shared<func<T>>(c2);
                }
                return make_shared<func<short>>(val * c2);
                break;
            }
            case integer_c: {
                auto val = (static_pointer_cast<constant<int>>(c1))->eval();
                if (val==0) {
                    return c1;
                }
                if (val==1) {
                    return make_shared<func<T>>(c2);
                }
                return make_shared<func<int>>(val * c2);
                break;
            }
            case float_c: {
                auto val = (static_pointer_cast<constant<float>>(c1))->eval();
                if (val==0) {
                    return c1;
                }
                if (val==1) {
                    return make_shared<func<T>>(c2);
                }
                return make_shared<func<float>>(val * c2);
                break;
            }
            case double_c: {
                auto val = (static_pointer_cast<constant<double>>(c1))->eval();
                if (val==0) {
                    return c1;
                }
                if (val==1) {
                    return make_shared<func<T>>(c2);
                }
                return make_shared<func<double>>(val * c2);
                break;
            }
            case long_c: {
                auto val = (static_pointer_cast<constant<long double>>(c1))->eval();
                if (val==0) {
                    return c1;
                }
                if (val==1) {
                    return make_shared<func<T>>(c2);
                }
                return make_shared<func<long double>>(val * c2);
                break;
            }
            case complex_c: {
                auto val = (static_pointer_cast<constant<Cpx>>(c1))->eval();
                if (val==Cpx(0,0)) {
                    return c1;
                }
                if (val==Cpx(1,1)) {
                    return make_shared<func<T>>(c2);
                }
                return make_shared<func<Cpx>>(val * c2);
                break;
            }
            case var_c: {
                auto v = static_pointer_cast<param_>(c1);
                switch (v->get_intype()) {
                    case binary_: {
                        auto vv = static_pointer_cast<var<bool>>(c1);
                        return make_shared<func<bool>>(vv*c2);
                        break;
                    }
                    case short_: {
                        auto vv = static_pointer_cast<var<short>>(c1);
                        return make_shared<func<short>>(vv*c2);
                        break;
                    }
                    case integer_: {
                        auto vv = static_pointer_cast<var<int>>(c1);
                        return make_shared<func<int>>(vv*c2);
                        break;
                    }
                    case float_: {
                        auto vv = static_pointer_cast<var<float>>(c1);
                        return make_shared<func<float>>(vv*c2);
                        break;
                    }
                    case double_: {
                        auto vv = static_pointer_cast<var<double>>(c1);
                        return make_shared<func<double>>(vv*c2);
                        break;
                    }
                    case long_: {
                        auto vv = static_pointer_cast<var<long double>>(c1);
                        return make_shared<func<long double>>(vv*c2);
                        break;
                    }
                    case complex_: {
                        auto vv = static_pointer_cast<var<Cpx>>(c1);
                        return make_shared<func<Cpx>>(vv*c2);
                        break;
                    }
                    default:
                        break;
                }
            }
        }
        return nullptr;
    }
//                    case integer_c: {
//                        auto val = (static_pointer_cast<constant<int>>(c1))->eval();
//                        if (val==0) {
//                            return c1;
//                        }
//                        if (val==1) {
//                            return make_shared<param<T>>(c2);
//                        }
//                        return make_shared<func<int>>(val * c2);
//                        break;
//                    }
//                    case float_c: {
//                        auto val = (static_pointer_cast<constant<float>>(c1))->eval();
//                        if (val==0) {
//                            return c1;
//                        }
//                        if (val==1) {
//                            return make_shared<param<T>>(c2);
//                        }
//                        return make_shared<func<float>>(val * c2);
//                        break;
//                    }
//                    case double_c: {
//                        auto val = (static_pointer_cast<constant<double>>(c1))->eval();
//                        if (val==0) {
//                            return c1;
//                        }
//                        if (val==1) {
//                            return make_shared<param<T>>(c2);
//                        }
//                        return make_shared<func<double>>(val * c2);
//                        break;
//                    }
//                    case long_c: {
//                        auto val = (static_pointer_cast<constant<long double>>(c1))->eval();
//                        if (val==0) {
//                            return c1;
//                        }
//                        if (val==1) {
//                            return make_shared<param<T>>(c2);
//                        }
//                        return make_shared<func<long double>>(val * c2);
//                        break;
//                    }
//                    case complex_c: {
//                        auto val = (static_pointer_cast<constant<Cpx>>(c1))->eval();
//                        if (val==Cpx(0,0)) {
//                            return c1;
//                        }
//                        if (val==Cpx(1,1)) {
//                            return make_shared<param<T>>(c2);
//                        }
//                        return make_shared<func<Cpx>>(val * c2);
//                        break;
//                    }
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

