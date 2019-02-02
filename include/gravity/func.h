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

//        func_ get_derivative(const param_& v) const; /**< Computes and returns the derivative with respect to variable v. */

//        func_ get_dfdx(const param_& v); /**< Computes all derivatives and returns a copy of the derivative with respect to variable v. */


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
        size_t get_nb_vars(unsigned inst) const{
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
        
        
        virtual bool insert(const constant_& coef, const param_& p1, const param_& p2, bool coef_p1_tr=false){
            return insert(true, coef, p1, p2, coef_p1_tr);
        };
        
        virtual bool insert(bool sign, const constant_& coef, const param_& p1, const param_& p2, bool c_p1_transposed=false){/**< Adds coef*p1*p2 to the function. Returns true if added new term, false if only updated coef of p1*p2 */
            auto ps1 = p1.get_name(false,false);
            auto ps2 = p2.get_name(false,false);
            auto qname = ps1+","+ps2;
            auto pair_it = _qterms->find(qname);
            shared_ptr<param_> p_new1;
            shared_ptr<param_> p_new2;
            
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
                    embed(*dynamic_pointer_cast<func_>(c_new));
                }
                _qterms->insert(make_pair<>(qname, qterm(sign, c_new, p_new1, p_new2)));
                if(p_new1->is_var()){
                    _evaluated = false;
                }
                //            update_convexity();
                return true;
            }
            else {
                if (pair_it->second._sign == sign) {
                    //                pair_it->second._coef = add(pair_it->second._coef, coef);
                }
                else{
                    //                pair_it->second._coef = substract(pair_it->second._coef, coef);
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
        
        virtual bool insert(bool sign, const constant_& coef, const list<pair<shared_ptr<param_>, int>>& l){/**< Adds polynomial term to the function. Returns true if added new term, false if only updated corresponding coef */
            _all_convexity = undet_;
            string name;
            string s;
            bool newv = true;
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
                    embed(*dynamic_pointer_cast<func_>(c_new));
                }
                pterm p(sign, c_new, newl);
                //            update_sign(p);
                _dim[0] = std::max(_dim[0], l.begin()->first->_dim[0]);
                _pterms->insert(make_pair<>(name, move(p)));
                if(pnew->is_var()){
                    _evaluated = false;
                }
                return true;
            }
            else {
                if (pair_it->second._sign == sign) {
                    //                pair_it->second._coef = add(pair_it->second._coef, coef);
                }
                else{
                    //                pair_it->second._coef = substract(pair_it->second._coef, coef);
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
                    //                update_sign();
                    //                update_convexity();
                }
                //            else {
                //                update_sign(pair_it->second);
                //            }
                return false;
            }
            
        };/**< Adds polynomial term to the function. Returns true if added new term, false if only updated corresponding coef */
        
        virtual bool insert(const constant_& coef, const param_& p, int exp){
            return insert(true, coef, p, exp);
        };
        
        virtual bool insert(bool sign, const constant_& coef, const param_& p, int exp){
            list<pair<shared_ptr<param_>, int>> l;
            l.push_back(make_pair<>(p.pcopy(), exp));
            return insert(sign, coef, l);
        };/**< Adds polynomial term to the function. Returns true if added new term, false if only updated corresponding coef */
        

        virtual bool insert(const qterm& term){
            return insert(term._sign, *term._coef, *term._p->first, *term._p->second, term._coef_p1_tr);
        };

        virtual bool insert(const pterm& term){return insert(term._sign, *term._coef, *term._l);};

        void update_sign_add(const constant_& c);
        void update_sign_multiply(const constant_& c);
        
        void update_quad_convexity();
        
        
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
                if (_cst->is_number()) {
                    auto val = _cst->to_str();
                    if (val.front()=='-') {
                        str += " - " + val.substr(1);
                    }
                    else if (val != "0"){
                        str += " + ";
                        str += val;
                    }
                }
                else {
                    str += " + ";
                    str += "(";
                    str += _cst->to_str();
                    str += ")";
                }
            }
            if (_expr) {
                str += " + ";
                str += _expr->to_str();
            }
            if (str.size() > 2 && str.at(1)=='+') {
                str = str.substr(3);
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
            return str;
        }
        
        void print_symbolic(bool endline = true, bool display_input = true);
    };

    template<class T1, class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) <= sizeof(T1)>::type* = nullptr>
    shared_ptr<pair<T1,T1>> get_product_range(shared_ptr<pair<T1,T1>> range1, shared_ptr<pair<T2,T2>> range2){
        shared_ptr<pair<T1,T1>> res = make_shared<pair<T1,T1>>();
        shared_ptr<pair<T1,T1>> cast_range1 = make_shared<pair<T1,T1>>(make_pair<>((T1)range1->first,(T1)range1->second));
        auto min1 = min(cast_range1->first*range2->first, cast_range1->first*range2->second);
        auto max1 = max(cast_range1->second*range2->second, cast_range1->second*range2->second);
        auto min2 = min(cast_range1->second*range2->second, cast_range1->first*range2->second);
        auto max2 = max(cast_range1->first*range2->first, cast_range1->first*range2->second);
        res->first = min(min1,min2);
        res->second = max(max1,max2);
        return res;
    }
    
    template<class T1, class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T1) < sizeof(T2)>::type* = nullptr>
    shared_ptr<pair<T2,T2>> get_product_range(shared_ptr<pair<T1,T1>> range1, shared_ptr<pair<T2,T2>> range2){
        shared_ptr<pair<T2,T2>> res = make_shared<pair<T2,T2>>();
        shared_ptr<pair<T2,T2>> cast_range1 = make_shared<pair<T2,T2>>(make_pair<>((T2)range1->first,(T2)range1->second));
        auto min1 = min(cast_range1->first*range2->first, cast_range1->first*range2->second);
        auto max1 = max(cast_range1->second*range2->second, cast_range1->second*range2->second);
        auto min2 = min(cast_range1->second*range2->second, cast_range1->first*range2->second);
        auto max2 = max(cast_range1->first*range2->first, cast_range1->first*range2->second);
        res->first = min(min1,min2);
        res->second = max(max1,max2);
        return res;
    }
    
    template<typename type = double>
    class func: public func_{
    public:
        shared_ptr<vector<type>>                _val = nullptr; /**< vector of values **/
        shared_ptr<pair<type,type>>             _range = nullptr; /**< (Min,Max) values in vals **/
        shared_ptr<vector<pair<type,type>>>     _all_range = nullptr; /**< Vector of (Min,Max) values for each instance of this func **/

        template<typename T=type,
        typename enable_if<is_arithmetic<T>::value>::type* = nullptr>
        void update_range(){
            _range = make_shared<pair<type,type>>(make_pair<>(numeric_limits<type>::max(), numeric_limits<type>::lowest()));
        }
        
        template<class T=type, class = typename enable_if<is_same<T, Cpx>::value>::type>
        void update_range() {
            _range = make_shared<pair<type,type>>(make_pair<>(Cpx(numeric_limits<double>::max(), numeric_limits<double>::max()), Cpx(numeric_limits<double>::lowest(), numeric_limits<double>::lowest())));
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
        
        template<class T, typename enable_if<is_convertible<T, type>::value && sizeof(T) <= sizeof(type)>::type* = nullptr>
        func get_derivative(shared_ptr<constant_> exp, const var<T>& v) const{
            auto name = v.get_name(false,false);
            if(exp->is_var()){
                auto vv = static_pointer_cast<param_>(exp);
                if(vv->get_name(false,false)==name){
                    return unit<type>();
                }
            }
            else if(exp->is_function()){
                auto f = static_pointer_cast<func>(exp);
                return f->get_derivative(v);
            }
            else if(exp->is_uexpr()){
                func son;
                auto uexp = static_pointer_cast<uexpr>(exp);
                if (uexp->_son->is_function()) {
                    auto f = dynamic_pointer_cast<func>(uexp->_son);
                    son = move(*f);
                }
                else if(uexp->_son->is_var()) {
                    auto vv = dynamic_pointer_cast<param_>(uexp->_son);
                    if(vv->get_name(false,false)==name){
                        son = v;
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
                        return uexp->_coef*get_derivative(uexp->_son,v)*son;
                        break;
                    }
                    case log_:
                        return uexp->_coef*get_derivative(uexp->_son,v)/son;
                        break;
                    default:
                        throw invalid_argument("Unsupported unary operation");
                        break;
                }
            }
            return func();
        }
        
        
        template<class T, typename enable_if<is_convertible<T, type>::value && sizeof(T) <= sizeof(type)>::type* = nullptr>
        func get_derivative(const var<T>& v) const{
            func res;
            if(!has_var(v)){
                return res;
            }
            auto name = v.get_name(false,false);
            for (auto &lt: *_lterms) {
                if (lt.second._p->get_name(false,false) == name) {
                    auto coef = lt.second._coef->copy();
                    if ((coef->_is_vector && coef->_is_transposed)) {
                        coef->transpose();
                    }
                    if (coef->is_function()) {
                        auto f_cst = dynamic_pointer_cast<func<type>>(coef);
                        res = move(*f_cst);
                    }
                    else if(coef->is_param()) {
                        auto p_cst = dynamic_pointer_cast<param<type>>(coef);
                        res = move(*p_cst);
                    }
                    else if(coef->is_number()) {
                        auto p_cst = dynamic_pointer_cast<constant<type>>(coef);
                        res = move(*p_cst);
                    }
                    if(!lt.second._sign){
                        res.reverse_sign();
                    }
                    break;
                }
            }
            for (auto &lt: *_qterms) {
                if (lt.second._p->first->get_name(false,false) == name) {
                    auto coef = lt.second._coef->copy();
                    if ((coef->_is_vector && coef->_is_transposed)) {
                        coef->transpose();//TODO is this needed?
                    }
                    if (coef->is_function()) {
                        auto f_cst = *dynamic_pointer_cast<func<type>>(coef);
                        res.insert(lt.second._sign, f_cst, *lt.second._p->second);
                    }
                    else if(coef->is_param()) {
                        auto p_cst = *dynamic_pointer_cast<param<type>>(coef);
                        res.insert(lt.second._sign, p_cst, *lt.second._p->second);
                    }
                    else if(coef->is_number()) {
                        auto p_cst = *dynamic_pointer_cast<constant<type>>(coef);
                        res.insert(lt.second._sign, p_cst, *lt.second._p->second);
                    }
                    if (lt.second._coef->_is_vector || lt.second._coef->is_matrix()) {
                        res._is_vector = true;
                    }
                }
                if (lt.second._p->second->get_name(false,false) == name) {
                    auto coef = lt.second._coef->copy();
                    if ((coef->_is_vector && coef->_is_transposed)) {
                        coef->transpose();//TODO is this needed?
                    }
                    if (coef->is_function()) {
                        auto f_cst = *dynamic_pointer_cast<func<type>>(coef);
                        res.insert(lt.second._sign, f_cst, *lt.second._p->first);
                    }
                    else if(coef->is_param()) {
                        auto p_cst = *dynamic_pointer_cast<param<type>>(coef);
                        res.insert(lt.second._sign, p_cst, *lt.second._p->first);
                    }
                    else if(coef->is_number()) {
                        auto p_cst = *dynamic_pointer_cast<constant<type>>(coef);
                        res.insert(lt.second._sign, p_cst, *lt.second._p->first);
                    }
                    if (lt.second._coef->_is_vector || lt.second._coef->is_matrix()) {
                        res._is_vector = true;
                    }
                }
            }
            int expo = 0;
            for (auto &lt: *_pterms) {
                bool has_v = false;
                auto newl(*lt.second._l);
                auto it = newl.begin();
                while(it!=newl.end()){
                    auto v = it->first;
                    if (v->get_name(false,false) == name) {
                        has_v = true;
                        expo = it->second;
                        it->second--;
                        if(it->second==0){
                            newl.erase(it);
                        }
                    }
                    it++;
                }
                if(has_v){
                    if (lt.second._coef->is_function()) {
                        auto f_cst = *dynamic_pointer_cast<func<type>>(lt.second._coef);
                        res.insert(lt.second._sign, expo*f_cst, newl);
                    }
                    else if(lt.second._coef->is_param()) {
                        auto p_cst = *dynamic_pointer_cast<param<type>>(lt.second._coef);
                        res.insert(lt.second._sign, expo*p_cst, newl);
                    }
                    else if(lt.second._coef->is_number()) {
                        auto p_cst = *dynamic_pointer_cast<constant<type>>(lt.second._coef);
                        res.insert(lt.second._sign, expo*p_cst, newl);
                    }
                }
            }
            if (!_expr) {
                //            res.untranspose();
                return res;
            }
            else { // f is a composition of functions
                res += get_derivative(_expr,v);
                return res;
            }
        }
        
        
        using func_::insert;
        
        bool insert(bool sign, const constant_& coef, const param_& p){/**< Adds coef*p to the linear function. Returns true if added new term, false if only updated coef of p */
            shared_ptr<param_> p_new;
            auto pname = p.get_name(false,false);
            auto pair_it = _lterms->find(pname);
            if (pair_it != _lterms->end() && pair_it->second._p->get_type() != p.get_type()) {
                throw invalid_argument("param and var with same name: " + pname);
            }
            if (_ftype == const_ && p.is_var()) {
                _ftype = lin_;
            }
            
            if (pair_it == _lterms->end()) {
                auto c_new = coef.copy();
                if (c_new->is_function()) {
                    embed(*dynamic_pointer_cast<func_>(c_new));
                }
                if (p.is_var()) {
                    p_new = get_var(pname);
                    if (!p_new) {
                        p_new = p.pcopy();
                        add_var(p_new);
                    }
                    else {
                        incr_occ_var(pname);
                    }
                }
                else {
                    p_new = get_param(pname);
                    if (!p_new) {
                        p_new = p.pcopy();
                        add_param(p_new);
                    }
                    else {
                        incr_occ_param(pname);
                    }
                }
                _lterms->insert(make_pair<>(pname, lterm(sign, c_new, p_new)));
                if(p_new->is_var()){
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
                    if (p.is_var()) {
                        decr_occ_var(pname);
                    }
                    else{
                        decr_occ_param(pname);
                    }
                    _lterms->erase(pair_it);
                    //update_sign();
                }
                //            else {
                //                update_sign(pair_it->second);
                //            }
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
        
        
        template<class T=type, class = typename enable_if<is_same<T, Cpx>::value>::type> void update_all_sign(){
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
//            //                auto vv = dynamic_pointer_cast<var<short>>(c1);
//            //                return make_shared<func<short>>(vv*c2);
//            //                break;
//            //            }
//            //            case integer_: {
//            //                auto vv = dynamic_pointer_cast<var<int>>(c1);
//            //                return make_shared<func<int>>(vv*c2);
//            //                break;
//            //            }
//            //            case float_: {
//            //                auto vv = dynamic_pointer_cast<var<float>>(c1);
//            //                return make_shared<func<float>>(vv*c2);
//            //                break;
//            //            }
//            //            case double_: {
//            //                auto vv = dynamic_pointer_cast<var<double>>(c1);
//            //                return make_shared<func<double>>(vv*c2);
//            //                break;
//            //            }
//            //            case long_: {
//            //                auto vv = dynamic_pointer_cast<var<long double>>(c1);
//            //                return make_shared<func<long double>>(vv*c2);
//            //                break;
//            //            }
//            //            case complex_: {
//            //                auto vv = dynamic_pointer_cast<var<Cpx>>(c1);
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
        
        void print(size_t index, int prec = 10) {
            cout << to_str(index,prec);
        }
        
        void print(size_t i, size_t j, int prec = 10) {
            cout << to_str(i,j,prec);
        }
        
        string to_str() {
            return func_::to_str();
        }
        
        string to_str(size_t index, int prec) {
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
                str = str.substr(3);
            }
            if (_is_vector && (is_number() || _vars->size()>1 || _params->size()>1)) {
                str = "[" + str +"]";
            }
            if (_is_transposed && (is_number() || _vars->size()>1 || _params->size()>1)) {
                str += "\u1D40";
            }
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
        
        void print(){
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
            auto nb_inst = _dim[0];
            allocate_mem();
            if (is_matrix()) {
                auto max_cell_size = get_max_cell_size();
                for (size_t i = 0; i<_dim[0]; i++) {
                    if (i>0) {
                        str.insert(str.end(), space_size, ' ');
                    }
                    str += "|";
                    for (size_t j = 0; j<_dim[1]; j++) {
                        auto cell = to_str(i,j,5);
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
                    str += to_str(inst,5);
                    str += "\n";
                }
            }
            str += "\n";
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
            if (_is_vector && (is_number() || _vars->size()>1 || _params->size()>1)) {
                str = "[" + str +"]";
            }
            if (_is_transposed && (is_number() || _vars->size()>1 || _params->size()>1)) {
                str += "\u1D40";
            }
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

        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        shared_ptr<constant_> multiply(shared_ptr<constant_> coef, const constant<T2>& c){
            if (coef->is_function()) {
                auto f_cst = *dynamic_pointer_cast<func<type>>(coef);
                f_cst *= func<type>(c);
                return f_cst.copy();
            }
            else if(coef->is_param()) {
                auto p_cst = *dynamic_pointer_cast<param<type>>(coef);
                if(c.is_unit()){
                    return p_cst.pcopy();
                }
                auto new_cst = p_cst * c;
                return new_cst.copy();
            }
            else if(coef->is_number()) {
                auto p_cst = *dynamic_pointer_cast<constant<type>>(coef);
                auto new_cst = p_cst * c;
                return new_cst.copy();
            }
            return nullptr;
        }
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        shared_ptr<constant_> multiply(shared_ptr<constant_> coef, const param<T2>& p){
            if (coef->is_function()) {
                auto f_cst = *dynamic_pointer_cast<func<type>>(coef);
                f_cst *= func<type>(p);
                return f_cst.copy();
            }
            else if(coef->is_param()) {
                auto p_cst = *dynamic_pointer_cast<param<type>>(coef);
                auto new_cst = p_cst * p;
                return new_cst.copy();
            }
            else if(coef->is_number()) {
                auto p_cst = *dynamic_pointer_cast<constant<type>>(coef);
                auto new_cst = p_cst * p;
                return new_cst.copy();
            }
            return nullptr;
        }
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        shared_ptr<constant_> multiply(shared_ptr<constant_> coef, const func<T2>& f){
            if (coef->is_function()) {
                auto f_cst = *dynamic_pointer_cast<func<type>>(coef);
                f_cst *= func<type>(f);
                embed(f_cst);
                return f_cst.copy();
            }
            else if(coef->is_param()) {
                auto p_cst = *dynamic_pointer_cast<param<type>>(coef);
                auto new_cst = p_cst * f;
                return new_cst.copy();
            }
            else if(coef->is_number()) {
                auto p_cst = func<type>(*dynamic_pointer_cast<constant<type>>(coef));
                return make_shared<func<type>>(p_cst *= f);
            }
            return nullptr;
        }
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        shared_ptr<constant_> add(shared_ptr<constant_> coef, const constant<T2>& c){
            if (coef->is_function()) {
                auto f_cst = *dynamic_pointer_cast<func<type>>(coef);
                f_cst += func<type>(c);
                return f_cst.copy();
            }
            else if(coef->is_param()) {
                auto p_cst = *dynamic_pointer_cast<param<type>>(coef);
                auto new_cst = p_cst + c;
                return new_cst.copy();
            }
            else if(coef->is_number()) {
                auto p_cst = *dynamic_pointer_cast<constant<type>>(coef);
                auto new_cst = p_cst + c;
                return new_cst.copy();
            }
            return nullptr;
        }
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        shared_ptr<constant_> add(shared_ptr<constant_> coef, const param<T2>& p){
            if (coef->is_function()) {
                auto f_cst = *dynamic_pointer_cast<func<type>>(coef);
                f_cst += func<type>(p);
                return f_cst.copy();
            }
            else if(coef->is_param()) {
                auto p_cst = *dynamic_pointer_cast<param<type>>(coef);
                auto new_cst = p_cst + p;
                return new_cst.copy();
            }
            else if(coef->is_number()) {
                auto p_cst = *dynamic_pointer_cast<constant<type>>(coef);
                auto new_cst = p_cst + p;
                return new_cst.copy();
            }
            return nullptr;
        }
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        shared_ptr<constant_> add(shared_ptr<constant_> coef, const func<T2>& f){
            if (coef->is_function()) {
                auto f_cst = *dynamic_pointer_cast<func<type>>(coef);
                f_cst += func<type>(f);
                embed(f_cst);
                return f_cst.copy();
            }
            else if(coef->is_param()) {
                auto p_cst = *dynamic_pointer_cast<param<type>>(coef);
                auto new_cst = p_cst + f;
                return new_cst.copy();
            }
            else if(coef->is_number()) {
                auto p_cst = func<type>(*dynamic_pointer_cast<constant<type>>(coef));
                return make_shared<func<type>>(p_cst += f);
            }
            return nullptr;
        }
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        shared_ptr<constant_> subtract(shared_ptr<constant_> coef, const constant<T2>& c){
            if (coef->is_function()) {
                auto f_cst = *dynamic_pointer_cast<func<type>>(coef);
                f_cst -= func<type>(c);
                return f_cst.copy();
            }
            else if(coef->is_param()) {
                auto p_cst = *dynamic_pointer_cast<param<type>>(coef);
                if(c.is_unit()){
                    return p_cst.pcopy();
                }
                auto new_cst = p_cst - c;
                return new_cst.copy();
            }
            else if(coef->is_number()) {
                auto p_cst = *dynamic_pointer_cast<constant<type>>(coef);
                auto new_cst = p_cst - c;
                return new_cst.copy();
            }
            return nullptr;
        }
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        shared_ptr<constant_> subtract(shared_ptr<constant_> coef, const param<T2>& p){
            if (coef->is_function()) {
                auto f_cst = *dynamic_pointer_cast<func<type>>(coef);
                f_cst -= func<type>(p);
                return f_cst.copy();
            }
            else if(coef->is_param()) {
                auto p_cst = *dynamic_pointer_cast<param<type>>(coef);
                auto new_cst = p_cst - p;
                return new_cst.copy();
            }
            else if(coef->is_number()) {
                auto p_cst = *dynamic_pointer_cast<constant<type>>(coef);
                auto new_cst = p_cst - p;
                return new_cst.copy();
            }
            return nullptr;
        }
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        shared_ptr<constant_> subtract(shared_ptr<constant_> coef, const func<T2>& f){
            if (coef->is_function()) {
                auto f_cst = *dynamic_pointer_cast<func<type>>(coef);
                f_cst -= func<type>(f);
                embed(f_cst);
                return f_cst.copy();
            }
            else if(coef->is_param()) {
                auto p_cst = *dynamic_pointer_cast<param<type>>(coef);
                auto new_cst = p_cst - f;
                return new_cst.copy();
            }
            else if(coef->is_number()) {
                auto p_cst = func<type>(*dynamic_pointer_cast<constant<type>>(coef));
                return make_shared<func<type>>(p_cst -= f);
            }
            return nullptr;
        }
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        void add_cst(const constant<T2>& f){
            if (_cst->is_function()) {
                auto f_cst = *dynamic_pointer_cast<func<type>>(_cst);
                f_cst += func<type>(f);
                _cst = make_shared<func<type>>(move(f_cst));
            }
            else if(_cst->is_param()) {
                auto p_cst = *dynamic_pointer_cast<param<type>>(_cst);
                auto new_cst = f + p_cst;
                _cst = make_shared<func<type>>(move(new_cst));
            }
            else if(_cst->is_number()) {
                auto p_cst = *dynamic_pointer_cast<constant<type>>(_cst);
                auto new_cst = f + p_cst;
                _cst = make_shared<constant<type>>(new_cst);
            }
        }
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        void add_cst(const param<T2>& f){
            if (_cst->is_function()) {
                auto f_cst = *dynamic_pointer_cast<func<type>>(_cst);
                f_cst += func<type>(f);
                embed(f_cst);
                _cst = make_shared<func<type>>(move(f_cst));                
            }
            else if(_cst->is_param()) {
                auto p_cst = *dynamic_pointer_cast<param<type>>(_cst);
                auto f_cst = f + p_cst;
                embed(f_cst);
                _cst = make_shared<func<type>>(move(f_cst));
            }
            else if(_cst->is_number()) {
                auto p_cst = *dynamic_pointer_cast<constant<type>>(_cst);
                auto f_cst = f + p_cst;
                embed(f_cst);
                _cst = make_shared<func<type>>(move(f_cst));
            }
        }
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        void add_cst(const func<T2>& f){
            if (_cst->is_function()) {
                auto f_cst = *dynamic_pointer_cast<func<type>>(_cst);
                f_cst += f;
                embed(f_cst);
                _cst = make_shared<func<type>>(move(f_cst));
            }
            else if(_cst->is_param()) {
                auto p_cst = *dynamic_pointer_cast<param<type>>(_cst);
                auto f_cst = f + func<type>(p_cst);
                embed(f_cst);
                _cst = make_shared<func<type>>(move(f_cst));
            }
            else if(_cst->is_number()) {
                auto p_cst = *dynamic_pointer_cast<constant<type>>(_cst);
                auto f_cst = f + func<type>(p_cst);
                embed(f_cst);
                _cst = make_shared<func<type>>(move(f_cst));
            }
        }
        
        func(const uexpr& ue):func(){
            _expr = make_shared<uexpr>(ue);
            embed(_expr);
            if (!is_constant()) {
                _ftype = nlin_;
            }
            _dim[0] = ue._dim[0];
            _dim[1] = ue._dim[1];
            _evaluated = false;
        };
        
        func(const bexpr& be):func(){
            _expr = make_shared<bexpr>(be);
            embed(_expr);
            if (!is_constant()) {
                _ftype = nlin_;
            }
            _dim[0] = be._dim[0];
            _dim[1] = be._dim[1];
            _evaluated = false;
        };
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        func(T2 c):func(){
            *this = constant<T2>(c);
        }
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        func(const constant<T2>& c):func(){
            *this = c;
        }
        template<class T2, class = typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type>
        func(const param<T2>& c):func(){
            *this = c;
        }
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        func(const var<T2>& c):func(){
            *this = c;
        }
        
        template<class T2, class = typename enable_if<is_convertible<T2, type>::value && sizeof(T2) < sizeof(type)>::type>
        func(const func<T2>& f): func(){
            *this = f;
        }
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        func& operator=(const constant<T2>& c){
            reset();
            dynamic_pointer_cast<constant<type>>(_cst)->set_val(c.eval());
            _all_sign = _cst->get_sign();
            _val->resize(1);
            _val->at(0) = c.eval();
            update_range(_val->at(0));
            _all_sign = c.get_all_sign();
            _evaluated = true;
            return *this;
        }
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        func& operator=(const param<T2>& c){
            reset();
            insert(true,unit<type>(),c);
            _val->clear();
            _range->first = c._range->first;
            _range->second = c._range->second;
            _all_sign = c.get_all_sign();
            _evaluated = false;
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
                    _expr = make_shared<uexpr>(*static_pointer_cast<uexpr>(f._expr));
                }
                else {
                    _expr = make_shared<bexpr>(*static_pointer_cast<bexpr>(f._expr));
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
                auto coef = *dynamic_pointer_cast<func<T2>>(f._cst);
                _cst = func(coef).copy();
            }
            else if(f._cst->is_param()) {
                auto coef = *dynamic_pointer_cast<param<T2>>(f._cst);
                _cst = param<type>(coef).copy();
            }
            else if(f._cst->is_number()) {
                auto coef = *dynamic_pointer_cast<constant<T2>>(f._cst);
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
                    auto coef = *dynamic_pointer_cast<func<T2>>(term._coef);
                    term._coef = func(coef).copy();
                }
                else if(term._coef->is_param()) {
                    auto coef = *dynamic_pointer_cast<param<T2>>(term._coef);
                    term._coef = param<type>(coef).copy();
                }
                else if(term._coef->is_number()) {
                    auto coef = *dynamic_pointer_cast<constant<T2>>(term._coef);
                    term._coef = constant<type>(coef.eval()).copy();
                }
                this->insert(term);
            }
            for (auto &pair:*f._qterms) {
                auto term = pair.second;
                if (term._coef->is_function()) {
                    auto coef = *dynamic_pointer_cast<func<T2>>(term._coef);
                    term._coef = func(coef).copy();
                }
                else if(term._coef->is_param()) {
                    auto coef = *dynamic_pointer_cast<param<T2>>(term._coef);
                    term._coef = param<type>(coef).copy();
                }
                else if(term._coef->is_number()) {
                    auto coef = *dynamic_pointer_cast<constant<T2>>(term._coef);
                    term._coef = constant<type>(coef.eval()).copy();
                }
                this->insert(term);
            }
            for (auto &pair:*f._pterms) {
                auto term = pair.second;
                if (term._coef->is_function()) {
                    auto coef = *dynamic_pointer_cast<func<T2>>(term._coef);
                    term._coef = func(coef).copy();
                }
                else if(term._coef->is_param()) {
                    auto coef = *dynamic_pointer_cast<param<T2>>(term._coef);
                    term._coef = param<type>(coef).copy();
                }
                else if(term._coef->is_number()) {
                    auto coef = *dynamic_pointer_cast<constant<T2>>(term._coef);
                    term._coef = constant<type>(coef.eval()).copy();
                }
                this->insert(term);
            }
            if(f._expr){
                if (f._expr->is_uexpr()) {
                    auto uexp = make_shared<uexpr>(*static_pointer_cast<uexpr>(f._expr));
                    if (uexp->_son->is_function()) {
                        auto f = dynamic_pointer_cast<func<T2>>(uexp->_son);
                        uexp->_son = make_shared<func>(*f);
                    }
                    _expr = uexp;
                }
                else {
                    auto bexp = make_shared<bexpr>(*static_pointer_cast<bexpr>(f._expr));
                    if (bexp->_lson->is_function()) {
                        auto f = dynamic_pointer_cast<func<T2>>(bexp->_lson);
                        bexp->_lson = make_shared<func>(*f);
                    }
                    if (bexp->_rson->is_function()) {
                        auto f = dynamic_pointer_cast<func<T2>>(bexp->_rson);
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
        
        type eval(size_t i) {
            if(is_matrix()){
                throw invalid_argument("eval() should be called with double index here\n");
            }
            if (is_constant() && _evaluated) {
                if (is_number()){
                    return _val->at(0);
                }
                return _val->at(i);
            }
            type res = zero<type>().eval();
            res += eval_cst(i);
            res += eval_lterms(i);
            res += eval_qterms(i);
            res += eval_pterms(i);
            if(_expr)
                res += eval_expr(_expr,i);
            return res;
        }
        
        type eval_cst(size_t i) {
            return eval_coef(_cst, i);
        }
        
        template<class T=type, typename enable_if<is_arithmetic<T>::value>::type* = nullptr>
        type eval(shared_ptr<constant_> c, size_t i) {
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
                    return static_pointer_cast<constant<double>>(c)->eval();
                    break;
                case long_c:
                    return static_pointer_cast<constant<long double>>(c)->eval();
                    break;
                case func_c:{
                    shared_ptr<func_> f = static_pointer_cast<func_>(c);
                    switch (f->get_return_type()) {
                        case binary_:
                            return static_pointer_cast<func<bool>>(f)->eval(i);
                            break;
                        case short_:
                            return static_pointer_cast<func<short>>(f)->eval(i);
                            break;
                        case integer_:
                            return static_pointer_cast<func<int>>(f)->eval(i);
                            break;
                        case float_:
                            return static_pointer_cast<func<float>>(f)->eval(i);
                            break;
                        case double_:
                            return static_pointer_cast<func<double>>(f)->eval(i);
                            break;
                        case long_:
                            return static_pointer_cast<func<long double>>(f)->eval(i);
                            break;
                        default:
                            break;
                    }
                    break;
                }
                case uexp_c:{
                    return eval_uexpr(static_pointer_cast<uexpr>(c),i);
                    break;
                }
                case bexp_c:{
                    return eval_bexpr(static_pointer_cast<bexpr>(c),i);
                    break;
                }
                default:{
                    shared_ptr<param_> p = static_pointer_cast<param_>(c);
                    switch (p->get_intype()) {
                        case binary_:
                            return static_pointer_cast<param<bool>>(p)->eval(i);
                            break;
                        case short_:
                            return static_pointer_cast<param<short>>(p)->eval(i);
                            break;
                        case integer_:
                            return static_pointer_cast<param<int>>(p)->eval(i);
                            break;
                        case float_:
                            return static_pointer_cast<param<float>>(p)->eval(i);
                            break;
                        case double_:
                            return static_pointer_cast<param<double>>(p)->eval(i);
                            break;
                        case long_:
                            return static_pointer_cast<param<long double>>(p)->eval(i);
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
        type eval(shared_ptr<constant_> c, size_t i, size_t j) {
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
                    return static_pointer_cast<constant<double>>(c)->eval();
                    break;
                case long_c:
                    return static_pointer_cast<constant<long double>>(c)->eval();
                    break;
                case func_c:{
                    shared_ptr<func_> f = static_pointer_cast<func_>(c);
                    switch (f->get_return_type()) {
                        case binary_:
                            return static_pointer_cast<func<bool>>(f)->eval(i,j);
                            break;
                        case short_:
                            return static_pointer_cast<func<short>>(f)->eval(i,j);
                            break;
                        case integer_:
                            return static_pointer_cast<func<int>>(f)->eval(i,j);
                            break;
                        case float_:
                            return static_pointer_cast<func<float>>(f)->eval(i,j);
                            break;
                        case double_:
                            return static_pointer_cast<func<double>>(f)->eval(i,j);
                            break;
                        case long_:
                            return static_pointer_cast<func<long double>>(f)->eval(i,j);
                            break;
                        default:
                            break;
                    }
                    break;
                }
                case uexp_c:{
                    return eval_uexpr(static_pointer_cast<uexpr>(c),i,j);
                    break;
                }
                case bexp_c:{
                    return eval_bexpr(static_pointer_cast<bexpr>(c),i,j);
                    break;
                }
                default:{
                    shared_ptr<param_> p = static_pointer_cast<param_>(c);
                    switch (p->get_intype()) {
                        case binary_:
                            return static_pointer_cast<param<bool>>(p)->eval(i,j);
                            break;
                        case short_:
                            return static_pointer_cast<param<short>>(p)->eval(i,j);
                            break;
                        case integer_:
                            return static_pointer_cast<param<int>>(p)->eval(i,j);
                            break;
                        case float_:
                            return static_pointer_cast<param<float>>(p)->eval(i,j);
                            break;
                        case double_:
                            return static_pointer_cast<param<double>>(p)->eval(i,j);
                            break;
                        case long_:
                            return static_pointer_cast<param<long double>>(p)->eval(i,j);
                            break;
                        default:
                            break;
                    }
                    break;
                }
            }
            throw invalid_argument("Unsupported type");
        }
        
        template<class T=type, class = typename enable_if<is_same<T, Cpx>::value>::type>
        type eval(shared_ptr<constant_> c, size_t i) {
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
                    return static_pointer_cast<constant<double>>(c)->eval();
                    break;
                case long_c:
                    return static_pointer_cast<constant<long double>>(c)->eval();
                    break;
                case complex_c:
                    return static_pointer_cast<constant<Cpx>>(c)->eval();
                    break;
                case func_c:{
                    shared_ptr<func_> f = static_pointer_cast<func_>(c);
                    switch (f->get_return_type()) {
                        case binary_:
                            return static_pointer_cast<func<bool>>(f)->eval(i);
                            break;
                        case short_:
                            return static_pointer_cast<func<short>>(f)->eval(i);
                            break;
                        case integer_:
                            return static_pointer_cast<func<int>>(f)->eval(i);
                            break;
                        case float_:
                            return static_pointer_cast<func<float>>(f)->eval(i);
                            break;
                        case double_:
                            return static_pointer_cast<func<double>>(f)->eval(i);
                            break;
                        case long_:
                            return static_pointer_cast<func<long double>>(f)->eval(i);
                            break;
                        case complex_:
                            return static_pointer_cast<func<Cpx>>(f)->eval(i);
                            break;
                        default:
                            break;
                    }
                    break;
                }
                case uexp_c:{
                    return eval(static_pointer_cast<uexpr>(c),i);
                    break;
                }
                case bexp_c:{
                    return eval(static_pointer_cast<bexpr>(c),i);
                    break;
                }
                default:{
                    shared_ptr<param_> p = static_pointer_cast<param_>(c);
                    switch (p->get_intype()) {
                        case binary_:
                            return static_pointer_cast<param<bool>>(p)->eval(i);
                            break;
                        case short_:
                            return static_pointer_cast<param<short>>(p)->eval(i);
                            break;
                        case integer_:
                            return static_pointer_cast<param<int>>(p)->eval(i);
                            break;
                        case float_:
                            return static_pointer_cast<param<float>>(p)->eval(i);
                            break;
                        case double_:
                            return static_pointer_cast<param<double>>(p)->eval(i);
                            break;
                        case long_:
                            return static_pointer_cast<param<long double>>(p)->eval(i);
                            break;
                        case complex_:
                            return static_pointer_cast<param<Cpx>>(p)->eval(i);
                            break;
                        default:
                            break;
                    }
                    break;
                }
            }
        }
        
        template<class T=type, class = typename enable_if<is_same<T, Cpx>::value>::type>
        type eval(shared_ptr<constant_> c, size_t i, size_t j) {
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
                    return static_pointer_cast<constant<double>>(c)->eval();
                    break;
                case long_c:
                    return static_pointer_cast<constant<long double>>(c)->eval();
                    break;
                case complex_c:
                    return static_pointer_cast<constant<Cpx>>(c)->eval();
                    break;
                case func_c:{
                    shared_ptr<func_> f = static_pointer_cast<func_>(c);
                    switch (f->get_return_type()) {
                        case binary_:
                            return static_pointer_cast<func<bool>>(f)->eval(i,j);
                            break;
                        case short_:
                            return static_pointer_cast<func<short>>(f)->eval(i,j);
                            break;
                        case integer_:
                            return static_pointer_cast<func<int>>(f)->eval(i,j);
                            break;
                        case float_:
                            return static_pointer_cast<func<float>>(f)->eval(i,j);
                            break;
                        case double_:
                            return static_pointer_cast<func<double>>(f)->eval(i,j);
                            break;
                        case long_:
                            return static_pointer_cast<func<long double>>(f)->eval(i,j);
                            break;
                        case complex_:
                            return static_pointer_cast<func<Cpx>>(f)->eval(i,j);
                            break;
                        default:
                            break;
                    }
                    break;
                }
                case uexp_c:{
                    return eval(static_pointer_cast<uexpr>(c),i,j);
                    break;
                }
                case bexp_c:{
                    return eval(static_pointer_cast<bexpr>(c),i,j);
                    break;
                }
                default:{
                    shared_ptr<param_> p = static_pointer_cast<param_>(c);
                    switch (p->get_intype()) {
                        case binary_:
                            return static_pointer_cast<param<bool>>(p)->eval(i,j);
                            break;
                        case short_:
                            return static_pointer_cast<param<short>>(p)->eval(i,j);
                            break;
                        case integer_:
                            return static_pointer_cast<param<int>>(p)->eval(i,j);
                            break;
                        case float_:
                            return static_pointer_cast<param<float>>(p)->eval(i,j);
                            break;
                        case double_:
                            return static_pointer_cast<param<double>>(p)->eval(i,j);
                            break;
                        case long_:
                            return static_pointer_cast<param<long double>>(p)->eval(i,j);
                            break;
                        case complex_:
                            return static_pointer_cast<param<Cpx>>(p)->eval(i,j);
                            break;
                        default:
                            break;
                    }
                    break;
                }
            }
        }
        
        type eval_coef(shared_ptr<constant_> coef, size_t i) {
            if (coef->is_function()) {
                auto f_cst = dynamic_pointer_cast<func<type>>(coef);
                return f_cst->eval(i);
            }
            else if(coef->is_param()) {
                auto p_cst = dynamic_pointer_cast<param<type>>(coef);
                return p_cst->eval(i);
            }
            else if(coef->is_number()) {
                auto p_cst = dynamic_pointer_cast<constant<type>>(coef);
                return p_cst->eval();
            }
            throw invalid_argument("in function eval_coef(shared_ptr<constant_> coef, size_t i), coef should be a constant");
        }
        
        type eval_coef(shared_ptr<constant_> coef, size_t i, size_t j) {
            if (coef->is_function()) {
                auto f_cst = dynamic_pointer_cast<func<type>>(coef);
                return f_cst->eval(i,j);
            }
            else if(coef->is_param()) {
                auto p_cst = dynamic_pointer_cast<param<type>>(coef);
                return p_cst->eval(i,j);
            }
            else if(coef->is_number()) {
                auto p_cst = dynamic_pointer_cast<constant<type>>(coef);
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
                if (!_sign) {
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
            if (!_sign) {
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
            if (!_sign) {
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
            if (!_sign) {
                res *= -1;
            }
            return res;
        }
        
        type eval_pterm(const pterm& pt, size_t i, size_t j){
            type res = zero<type>().eval();
            if (pt._coef->_is_transposed) {
                throw invalid_argument("Unspported operation\n");
            }// TREAT TRANSPOSED VECTORS IN POLYNOMIAL TERMS HERE
            else {
                res += 1;
                for (auto &pair: *pt._l) {
                    res *= pow(eval(pair.first,i,j), pair.second);
                }
                
                res *= eval_coef(pt._coef,i,j);
            }
            if (!_sign) {
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
        
        type eval_expr(shared_ptr<expr> exp, size_t i) {
            if (exp->is_uexpr()) {
                return eval_uexpr(static_pointer_cast<uexpr>(exp),i);
            }
            return eval_bexpr(static_pointer_cast<bexpr>(exp),i);
        }
        
        template<class T=type, typename enable_if<is_arithmetic<T>::value>::type* = nullptr>
        T eval_uexpr(shared_ptr<uexpr> exp, size_t i) {
            T res = eval(exp->_son,i);
            switch (exp->_otype) {
                case cos_:
                    return exp->_coef*std::cos(res);
                    break;
                case sin_:
                    return exp->_coef*std::sin(res);
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
                    break;
                default:
                    throw invalid_argument("Unsupported unary operator");
                    break;
            }
        }
        template<class T=type, class = typename enable_if<is_same<T, Cpx>::value>::type>
        Cpx eval_uexpr(shared_ptr<uexpr> exp, size_t i) {
            Cpx res = eval(exp->_son,i);
            switch (exp->_otype) {
                case cos_:
                    return exp->_coef*std::cos(res);
                    break;
                case sin_:
                    return exp->_coef*std::sin(res);
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
                    break;
                default:
                    throw invalid_argument("Unsupported unary operator");
                    break;
            }
        }
                       
        template<class T=type, class = typename enable_if<is_same<T, Cpx>::value>::type>
        T  eval_bexpr(shared_ptr<bexpr> exp, size_t i){
            T lval = eval(exp->_lson,i);
            T rval = eval(exp->_rson,i);
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
                default:
                    throw invalid_argument("Unsupported binary operator");
                    break;
            }
            
        }
        
        template<class T=type, typename enable_if<is_arithmetic<T>::value>::type* = nullptr>
        T  eval_bexpr(shared_ptr<bexpr> exp, size_t i){
            if (exp->_lson->is_constant() && !exp->_lson->is_evaluated()) {
                for (auto inst = 0; inst < exp->_lson->get_dim(); inst++) {
                    eval(exp->_lson,inst);
                }
                exp->_lson->evaluated(true);
            }
            if (exp->_rson->is_constant() && !exp->_rson->is_evaluated()) {
                for (auto inst = 0; inst < exp->_rson->get_dim(); inst++) {
                    eval(exp->_rson,inst);
                }
                exp->_rson->evaluated(true);
            }
            T lval = eval(exp->_lson,i);
            T rval = eval(exp->_rson,i);
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
                default:
                    throw invalid_argument("Unsupported binary operator");
                    break;
            }
            
        }
        
        
        template<class T=type, typename enable_if<is_arithmetic<T>::value>::type* = nullptr>
        T eval_uexpr(shared_ptr<uexpr> exp, size_t i, size_t j) {
            T res = eval(exp->_son,i,j);
            switch (exp->_otype) {
                case cos_:
                    return exp->_coef*std::cos(res);
                    break;
                case sin_:
                    return exp->_coef*std::sin(res);
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
                    break;
                default:
                    throw invalid_argument("Unsupported unary operator");
                    break;
            }
        }
        template<class T=type, class = typename enable_if<is_same<T, Cpx>::value>::type>
        Cpx eval_uexpr(shared_ptr<uexpr> exp, size_t i, size_t j) {
            Cpx res = eval(exp->_son,i,j);
            switch (exp->_otype) {
                case cos_:
                    return exp->_coef*std::cos(res);
                    break;
                case sin_:
                    return exp->_coef*std::sin(res);
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
                    break;
                default:
                    throw invalid_argument("Unsupported unary operator");
                    break;
            }
        }
        
        template<class T=type, class = typename enable_if<is_same<T, Cpx>::value>::type>
        T  eval_bexpr(shared_ptr<bexpr> exp, size_t i, size_t j){
            if (exp->_lson->is_constant() && !exp->_lson->is_evaluated()) {
                for (auto inst = 0; inst < exp->_lson->get_dim(); inst++) {
                    eval(exp->_lson,inst);
                }
                exp->_lson->evaluated(true);
            }
            if (exp->_rson->is_constant() && !exp->_rson->is_evaluated()) {
                for (auto inst = 0; inst < exp->_rson->get_dim(); inst++) {
                    eval(exp->_rson,inst);
                }
                exp->_rson->evaluated(true);
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
                default:
                    throw invalid_argument("Unsupported binary operator");
                    break;
            }
            
        }
        
        template<class T=type, typename enable_if<is_arithmetic<T>::value>::type* = nullptr>
        T  eval_bexpr(shared_ptr<bexpr> exp, size_t i, size_t j){
            if (exp->_lson->is_constant() && !exp->_lson->is_evaluated()) {
                for (auto inst = 0; inst < exp->_lson->get_dim(); inst++) {
                    eval(exp->_lson,inst);
                }
                exp->_lson->evaluated(true);
            }
            if (exp->_rson->is_constant() && !exp->_rson->is_evaluated()) {
                for (auto inst = 0; inst < exp->_rson->get_dim(); inst++) {
                    eval(exp->_rson,inst);
                }
                exp->_rson->evaluated(true);
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
                default:
                    throw invalid_argument("Unsupported binary operator");
                    break;
            }
            
        }
        
        
        type eval(const string& key) {
            return _val->at(_indices->_keys_map->at(key));
        }
        
        type eval(size_t i, size_t j) {
            
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
            return (_range->first == Cpx(1,0) && _range->second == Cpx(1,0));
        }
        
        bool is_zero() const { return zero_range();};
        
        template<class T=type, class = typename enable_if<is_same<T, Cpx>::value>::type> bool zero_range() const{
            return (_range->first == Cpx(0,0) && _range->second == Cpx(0,0));
        }
        
        template<typename T=type,
        typename enable_if<is_arithmetic<T>::value>::type* = nullptr> bool zero_range() const{
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
            auto be = bexpr(product_, make_shared<func>(*this), make_shared<constant<T2>>(c));
            *this = func(be);
            _evaluated = false;
            _all_convexity = undet_;
            return *this;
        }
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        func& operator/=(const param<T2>& p){
            auto be = bexpr(product_, make_shared<func>(*this), make_shared<param<T2>>(p));
            *this = func(be);
            _evaluated = false;
            _all_convexity = undet_;
            return *this;
        }
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) <= sizeof(type)>::type* = nullptr>
        func& operator/=(const func<T2>& f){
            if(!is_constant() && f.is_constant()){
                return *this *= 1/f;
            }
            auto be = bexpr(product_, make_shared<func>(*this), make_shared<func>(f));
            *this = func(be);
            _evaluated = false;
            _all_convexity = undet_;
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
            if(!is_constant() && f.is_constant()) {
                bool transp = false;
                auto fc = f;
                if(is_linear() && _is_transposed){// Situation where (*this)^T * f is transformed into (f^T*(*this))^T
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
                        _expr = make_shared<bexpr>(bexpr(product_, make_shared<func>(*static_pointer_cast<uexpr>(_expr)), make_shared<func>(fc)));
                    }
                    else {
                        _expr = make_shared<bexpr>(bexpr(product_, make_shared<func>(*static_pointer_cast<bexpr>(_expr)), make_shared<func>(fc)));
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
                _range = get_product_range(_range,f._range);
                if(transp){
                    this->transpose();
                    _range->first *= _dim[0];
                    _range->second *= _dim[0];
                }
                update_dot_dim(fc);
                return *this;
            }
            /* Case where the current function is constant and the other operand is not. */
            if (is_constant() && !f.is_constant()) {
                auto cpy = this->copy();
                update_dot_dim(f);
                update_sign_multiply(f);
                func res = f;
                res._dim[0] = _dim[0];
                res._dim[1] = _dim[1];
                res._all_sign = _all_sign;
                res._range = get_product_range(_range,f._range);
                if(_is_transposed){
                    res._range->first *= _dim[0];
                    res._range->second *= _dim[0];
                }
                if(is_non_positive()){
                    res.reverse_convexity();
                }
                if (!res._cst->is_zero()) {
                    if (res._cst->is_function()) {
                        auto f_cst = *dynamic_pointer_cast<func<T2>>(res._cst);
                        res._cst = multiply(cpy,f_cst);
                    }
                    else if(res._cst->is_param()) {
                        auto f_cst = *dynamic_pointer_cast<param<T2>>(res._cst);
                        res._cst = multiply(cpy,f_cst);
                    }
                    else if(res._cst->is_number()) {
                        auto f_cst = *dynamic_pointer_cast<constant<T2>>(res._cst);
                        res._cst = multiply(cpy,f_cst);
                    }
                    
                }
                for (auto &pair:*res._lterms) {
                    if (pair.second._coef->is_function()) {
                        auto f_cst = *dynamic_pointer_cast<func<T2>>(pair.second._coef);
                        pair.second._coef = multiply(cpy,f_cst);
                    }
                    else if(pair.second._coef->is_param()) {
                        auto f_cst = *dynamic_pointer_cast<param<T2>>(pair.second._coef);
                        pair.second._coef = multiply(cpy,f_cst);
                    }
                    else if(pair.second._coef->is_number()) {
                        auto f_cst = *dynamic_pointer_cast<constant<T2>>(pair.second._coef);
                        pair.second._coef = multiply(cpy,f_cst);
                    }
                }
                for (auto &pair:*res._qterms) {
                    if (pair.second._coef->is_function()) {
                        auto f_cst = *dynamic_pointer_cast<func<T2>>(pair.second._coef);
                        pair.second._coef = multiply(cpy,f_cst);
                    }
                    else if(pair.second._coef->is_param()) {
                        auto f_cst = *dynamic_pointer_cast<param<T2>>(pair.second._coef);
                        pair.second._coef = multiply(cpy,f_cst);
                    }
                    else if(pair.second._coef->is_number()) {
                        auto f_cst = *dynamic_pointer_cast<constant<T2>>(pair.second._coef);
                        pair.second._coef = multiply(cpy,f_cst);
                    }
                }
                for (auto &pair:*res._pterms) {
                    if (pair.second._coef->is_function()) {
                        auto f_cst = *dynamic_pointer_cast<func<T2>>(pair.second._coef);
                        pair.second._coef = multiply(cpy,f_cst);
                    }
                    else if(pair.second._coef->is_param()) {
                        auto f_cst = *dynamic_pointer_cast<param<T2>>(pair.second._coef);
                        pair.second._coef = multiply(cpy,f_cst);
                    }
                    else if(pair.second._coef->is_number()) {
                        auto f_cst = *dynamic_pointer_cast<constant<T2>>(pair.second._coef);
                        pair.second._coef = multiply(cpy,f_cst);
                    }
                }
                if (res._expr) {
                    if(res._expr->is_uexpr()){
                        res._expr = make_shared<bexpr>(bexpr(product_, make_shared<func<type>>(*this), make_shared<func>(*static_pointer_cast<uexpr>(res._expr))));
                    }
                    else {
                        res._expr = make_shared<bexpr>(bexpr(product_, make_shared<func<type>>(*this), make_shared<func>(*static_pointer_cast<bexpr>(res._expr))));
                    }                    
                    res.embed(res._expr);
                }
                *this = res;
                _evaluated = false;
                return *this;
            }
            //Both functions are non-constants at this stage
            if (_expr || (f._expr)) {
                auto be = bexpr(product_, make_shared<func>(*this), make_shared<func>(f));
                *this = func(be);
                _evaluated = false;
                _all_convexity = undet_;
                return *this;
            }

            /* Case where the multiplication invlolves multiplying variables/parameters together, i.e., they are both parametric or both include variables  */
            func res;
            res._range = get_product_range(_range,f._range);
            if(_is_transposed){
                res._range->first *= _dim[0];
                res._range->second *= _dim[0];
            }
            for (auto& t1: *_pterms) {
                if (t1.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial function), we cannot factor the coefficients. Just create a binary expression and return it.
                    auto be = bexpr(product_, make_shared<func>(*this), make_shared<func>(f));
                    *this = func(be);
                    _evaluated = false;
                    return *this;
                }
                for (auto& t2: *f._pterms) {
                    if (t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial function), see comment above.
                        auto be = bexpr(product_, make_shared<func>(*this), make_shared<func>(f));
                        *this = func(be);
                        _evaluated = false;
                        return *this;
                    }
                    auto newl(*t1.second._l);
                    for (auto& it: *t2.second._l) {// TODO check if same l
                        newl.push_back(make_pair<>(it.first, it.second));
                    }
                    if (t2.second._coef->is_function()) {
                        auto f_cst = *dynamic_pointer_cast<func<T2>>(t2.second._coef);
                        auto coef = multiply(t1.second._coef, f_cst);
                        res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
                    }
                    else if(t2.second._coef->is_param()) {
                        auto p_cst = *dynamic_pointer_cast<param<T2>>(t2.second._coef);
                        auto coef = multiply(t1.second._coef, p_cst);
                        res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
                    }
                    else if(t2.second._coef->is_number()) {
                        auto p_cst = *dynamic_pointer_cast<constant<T2>>(t2.second._coef);
                        auto coef = multiply(t1.second._coef, p_cst);
                        res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
                    }
                }
                for (auto& t2: *f._qterms) {
                    if (t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial function), see comment above.
                        auto be = bexpr(product_, make_shared<func>(*this), make_shared<func>(f));
                        *this = func(be);
                        _evaluated = false;
                        return *this;
                    }
                    auto newl(*t1.second._l);
                    newl.push_back(make_pair<>((t2.second._p->first), 1));
                    newl.push_back(make_pair<>((t2.second._p->second), 1));
                    if (t2.second._coef->is_function()) {
                        auto f_cst = *dynamic_pointer_cast<func<T2>>(t2.second._coef);
                        auto coef = multiply(t1.second._coef, f_cst);
                        res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
                    }
                    else if(t2.second._coef->is_param()) {
                        auto p_cst = *dynamic_pointer_cast<param<T2>>(t2.second._coef);
                        auto coef = multiply(t1.second._coef, p_cst);
                        res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
                    }
                    else if(t2.second._coef->is_number()) {
                        auto p_cst = *dynamic_pointer_cast<constant<T2>>(t2.second._coef);
                        auto coef = multiply(t1.second._coef, p_cst);
                        res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
                    }
                }
                for (auto& t2: *f._lterms) {
                    if (t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial function)
                        auto be = bexpr(product_, make_shared<func>(*this), make_shared<func>(f));
                        *this = func(be);
                        _evaluated = false;
                        return *this;
                    }
                    auto newl(*t1.second._l);
                    newl.push_back(make_pair<>((t2.second._p), 1));
                    if (t2.second._coef->is_function()) {
                        auto f_cst = *dynamic_pointer_cast<func<T2>>(t2.second._coef);
                        auto coef = multiply(t1.second._coef, f_cst);
                        res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
                    }
                    else if(t2.second._coef->is_param()) {
                        auto p_cst = *dynamic_pointer_cast<param<T2>>(t2.second._coef);
                        auto coef = multiply(t1.second._coef, p_cst);
                        res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
                    }
                    else if(t2.second._coef->is_number()) {
                        auto p_cst = *dynamic_pointer_cast<constant<T2>>(t2.second._coef);
                        auto coef = multiply(t1.second._coef, p_cst);
                        res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
                    }
                }
                if (!f._cst->is_zero()) {
                    auto newl(*t1.second._l);
                    if (f._cst->is_function()) {
                        auto f_cst = *dynamic_pointer_cast<func<T2>>(f._cst);
                        auto coef = multiply(t1.second._coef, f_cst);
                        res.insert(t1.second._sign, *coef, newl);
                    }
                    else if(f._cst->is_param()) {
                        auto p_cst = *dynamic_pointer_cast<param<T2>>(f._cst);
                        auto coef = multiply(t1.second._coef, p_cst);
                        res.insert(t1.second._sign, *coef, newl);
                    }
                    else if(f._cst->is_number()) {
                        auto p_cst = *dynamic_pointer_cast<constant<T2>>(f._cst);
                        auto coef = multiply(t1.second._coef, p_cst);
                        res.insert(t1.second._sign, *coef, newl);
                    }
                }
            }
            for (auto& t1: *_qterms) {
                if (t1.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(Quadratic term)
                    auto be = bexpr(product_, make_shared<func>(*this), make_shared<func>(f));
                    *this = func(be);
                    _evaluated = false;
                    return *this;
                }
                for (auto& t2: *f._pterms) {
                    if (t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial term)
                        auto be = bexpr(product_, make_shared<func>(*this), make_shared<func>(f));
                        *this = func(be);
                        _evaluated = false;
                    }
                    auto newl(*t2.second._l);
                    newl.push_front(make_pair<>(t1.second._p->first, 1));
                    newl.push_front(make_pair<>(t1.second._p->second, 1));
                    if (t2.second._coef->is_function()) {
                        auto f_cst = *dynamic_pointer_cast<func<T2>>(t2.second._coef);
                        auto coef = multiply(t1.second._coef, f_cst);
                        res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
                    }
                    else if(t2.second._coef->is_param()) {
                        auto p_cst = *dynamic_pointer_cast<param<T2>>(t2.second._coef);
                        auto coef = multiply(t1.second._coef, p_cst);
                        res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
                    }
                    else if(t2.second._coef->is_number()) {
                        auto p_cst = *dynamic_pointer_cast<constant<T2>>(t2.second._coef);
                        auto coef = multiply(t1.second._coef, p_cst);
                        res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
                    }
                }
                for (auto& t2: *f._qterms) {
                    if (t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial term)
                        auto be = bexpr(product_, make_shared<func>(*this), make_shared<func>(f));
                        *this = func(be);
                        _evaluated = false;
                    }
                    list<pair<shared_ptr<param_>, int>> newl;
                    newl.push_back(make_pair<>(t1.second._p->first, 1));
                    newl.push_back(make_pair<>(t1.second._p->second, 1));
                    newl.push_back(make_pair<>(t2.second._p->first, 1));
                    newl.push_back(make_pair<>(t2.second._p->second, 1));
                    if (t2.second._coef->is_function()) {
                        auto f_cst = *dynamic_pointer_cast<func<T2>>(t2.second._coef);
                        auto coef = multiply(t1.second._coef, f_cst);
                        res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
                    }
                    else if(t2.second._coef->is_param()) {
                        auto p_cst = *dynamic_pointer_cast<param<T2>>(t2.second._coef);
                        auto coef = multiply(t1.second._coef, p_cst);
                        res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
                    }
                    else if(t2.second._coef->is_number()) {
                        auto p_cst = *dynamic_pointer_cast<constant<T2>>(t2.second._coef);
                        auto coef = multiply(t1.second._coef, p_cst);
                        res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
                    }
                }
                for (auto& t2: *f._lterms) {
                    if (t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial term)
                        auto be = bexpr(product_, make_shared<func>(*this), make_shared<func>(f));
                        *this = func(be);
                        _evaluated = false;
                    }
                    list<pair<shared_ptr<param_>, int>> newl;
                    newl.push_back(make_pair<>(t1.second._p->first, 1));
                    newl.push_back(make_pair<>(t1.second._p->second, 1));
                    newl.push_back(make_pair<>(t2.second._p, 1));
                    if (t2.second._coef->is_function()) {
                        auto f_cst = *dynamic_pointer_cast<func<T2>>(t2.second._coef);
                        auto coef = multiply(t1.second._coef, f_cst);
                        res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
                    }
                    else if(t2.second._coef->is_param()) {
                        auto p_cst = *dynamic_pointer_cast<param<T2>>(t2.second._coef);
                        auto coef = multiply(t1.second._coef, p_cst);
                        res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
                    }
                    else if(t2.second._coef->is_number()) {
                        auto p_cst = *dynamic_pointer_cast<constant<T2>>(t2.second._coef);
                        auto coef = multiply(t1.second._coef, p_cst);
                        res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
                    }                    }
                if (!f._cst->is_zero()) {
                    if (f._cst->is_function()) {
                        auto f_cst = *dynamic_pointer_cast<func<T2>>(f._cst);
                        auto coef = multiply(t1.second._coef, f_cst);
                        res.insert(t1.second._sign, *coef, *t1.second._p->first, *t1.second._p->second);
                    }
                    else if(f._cst->is_param()) {
                        auto p_cst = *dynamic_pointer_cast<param<T2>>(f._cst);
                        auto coef = multiply(t1.second._coef, p_cst);
                        res.insert(t1.second._sign, *coef, *t1.second._p->first, *t1.second._p->second);
                    }
                    else if(f._cst->is_number()) {
                        auto p_cst = *dynamic_pointer_cast<constant<T2>>(f._cst);
                        auto coef = multiply(t1.second._coef, p_cst);
                        res.insert(t1.second._sign, *coef, *t1.second._p->first, *t1.second._p->second);
                    }
                }
            }
            for (auto& t1: *_lterms) {
                for (auto& t2: *f._pterms) {
                    if (t1.second._coef->_is_transposed || t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial term)
                        auto be = bexpr(product_, make_shared<func>(*this), make_shared<func>(f));
                        *this = func(be);
                        _evaluated = false;
                    }
                    auto newl(*t2.second._l);
                    newl.push_front(make_pair<>((t1.second._p), 1));
                    if (t2.second._coef->is_function()) {
                        auto f_cst = *dynamic_pointer_cast<func<T2>>(t2.second._coef);
                        auto coef = multiply(t1.second._coef, f_cst);
                        res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
                    }
                    else if(t2.second._coef->is_param()) {
                        auto p_cst = *dynamic_pointer_cast<param<T2>>(t2.second._coef);
                        auto coef = multiply(t1.second._coef, p_cst);
                        res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
                    }
                    else if(t2.second._coef->is_number()) {
                        auto p_cst = *dynamic_pointer_cast<constant<T2>>(t2.second._coef);
                        auto coef = multiply(t1.second._coef, p_cst);
                        res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
                    }
                }
                for (auto& t2: *f._qterms) {
                    if (t1.second._coef->_is_transposed || t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial term)
                        auto be = bexpr(product_, make_shared<func>(*this), make_shared<func>(f));
                        *this = func(be);
                        _evaluated = false;
                    }
                    list<pair<shared_ptr<param_>, int>> newl;
                    newl.push_back(make_pair<>(t1.second._p, 1));
                    newl.push_back(make_pair<>(t2.second._p->first, 1));
                    newl.push_back(make_pair<>(t2.second._p->second, 1));
                    if (t2.second._coef->is_function()) {
                        auto f_cst = *dynamic_pointer_cast<func<T2>>(t2.second._coef);
                        auto coef = multiply(t1.second._coef, f_cst);
                        res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
                    }
                    else if(t2.second._coef->is_param()) {
                        auto p_cst = *dynamic_pointer_cast<param<T2>>(t2.second._coef);
                        auto coef = multiply(t1.second._coef, p_cst);
                        res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
                    }
                    else if(t2.second._coef->is_number()) {
                        auto p_cst = *dynamic_pointer_cast<constant<T2>>(t2.second._coef);
                        auto coef = multiply(t1.second._coef, p_cst);
                        res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
                    }
                }
                for (auto& t2: *f._lterms) {
                    if (t2.second._coef->is_function()) {
                        auto f_cst = *dynamic_pointer_cast<func<T2>>(t2.second._coef);
                        auto coef = multiply(t1.second._coef, f_cst);
                        res.insert(!(t1.second._sign^t2.second._sign), *coef, *t1.second._p, *t2.second._p, _is_transposed);
                    }
                    else if(t2.second._coef->is_param()) {
                        auto p_cst = *dynamic_pointer_cast<param<T2>>(t2.second._coef);
                        auto coef = multiply(t1.second._coef, p_cst);
                        res.insert(!(t1.second._sign^t2.second._sign), *coef, *t1.second._p, *t2.second._p, _is_transposed);
                    }
                    else if(t2.second._coef->is_number()) {
                        auto p_cst = *dynamic_pointer_cast<constant<T2>>(t2.second._coef);
                        auto coef = multiply(t1.second._coef, p_cst);
                        res.insert(!(t1.second._sign^t2.second._sign), *coef, *t1.second._p, *t2.second._p, _is_transposed);
                    }
                }
                if (!f._cst->is_zero()) {
                    if (f._cst->is_function()) {
                        auto f_cst = *dynamic_pointer_cast<func<T2>>(f._cst);
                        auto coef = multiply(t1.second._coef, f_cst);
                        res.insert(t1.second._sign, *coef, *t1.second._p);
                    }
                    else if(f._cst->is_param()) {
                        auto p_cst = *dynamic_pointer_cast<param<T2>>(f._cst);
                        auto coef = multiply(t1.second._coef, p_cst);
                        res.insert(t1.second._sign, *coef, *t1.second._p);
                    }
                    else if(f._cst->is_number()) {
                        auto p_cst = *dynamic_pointer_cast<constant<T2>>(f._cst);
                        auto coef = multiply(t1.second._coef, p_cst);
                        res.insert(t1.second._sign, *coef, *t1.second._p);
                    }
                }
            }
            if (!_cst->is_zero()) {
                for (auto& t2: *f._pterms) {
                    if (t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial term)
                        auto be = bexpr(product_, make_shared<func>(*this), make_shared<func>(f));
                        *this = func(be);
                        _evaluated = false;
                        return *this;
                    }
                    if (t2.second._coef->is_function()) {
                        auto f_cst = *dynamic_pointer_cast<func<T2>>(t2.second._coef);
                        auto coef = multiply(_cst, f_cst);
                        res.insert(t2.second._sign, *coef, *t2.second._l);
                    }
                    else if(t2.second._coef->is_param()) {
                        auto p_cst = *dynamic_pointer_cast<param<T2>>(t2.second._coef);
                        auto coef = multiply(_cst, p_cst);
                        res.insert(t2.second._sign, *coef, *t2.second._l);
                    }
                    else if(t2.second._coef->is_number()) {
                        auto p_cst = *dynamic_pointer_cast<constant<T2>>(t2.second._coef);
                        auto coef = multiply(_cst, p_cst);
                        res.insert(t2.second._sign, *coef, *t2.second._l);
                    }
                }
                for (auto& t2: *f._qterms) {
                    if (t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial term)
                        auto be = bexpr(product_, make_shared<func>(*this), make_shared<func>(f));
                        *this = func(be);
                        _evaluated = false;
                        return *this;
                    }
                    if (t2.second._coef->is_function()) {
                        auto f_cst = *dynamic_pointer_cast<func<T2>>(t2.second._coef);
                        auto coef = multiply(_cst, f_cst);
                        res.insert(t2.second._sign, *coef, *t2.second._p->first, *t2.second._p->second);
                    }
                    else if(t2.second._coef->is_param()) {
                        auto p_cst = *dynamic_pointer_cast<param<T2>>(t2.second._coef);
                        auto coef = multiply(_cst, p_cst);
                        res.insert(t2.second._sign, *coef, *t2.second._p->first, *t2.second._p->second);
                    }
                    else if(t2.second._coef->is_number()) {
                        auto p_cst = *dynamic_pointer_cast<constant<T2>>(t2.second._coef);
                        auto coef = multiply(_cst, p_cst);
                        res.insert(t2.second._sign, *coef, *t2.second._p->first, *t2.second._p->second);
                    }
                }
                for (auto& t2: *f._lterms) {
                    if (t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial term)
                        auto be = bexpr(product_, make_shared<func>(*this), make_shared<func>(f));
                        *this = func(be);
                        _evaluated = false;
                        return *this;
                    }
                    if (t2.second._coef->is_function()) {
                        auto f_cst = *dynamic_pointer_cast<func<T2>>(t2.second._coef);
                        auto coef = multiply(_cst, f_cst);
                        res.insert(t2.second._sign, *coef, *t2.second._p);
                    }
                    else if(t2.second._coef->is_param()) {
                        auto p_cst = *dynamic_pointer_cast<param<T2>>(t2.second._coef);
                        auto coef = multiply(_cst, p_cst);
                        res.insert(t2.second._sign, *coef, *t2.second._p);
                    }
                    else if(t2.second._coef->is_number()) {
                        auto p_cst = *dynamic_pointer_cast<constant<T2>>(t2.second._coef);
                        auto coef = multiply(_cst, p_cst);
                        res.insert(t2.second._sign, *coef, *t2.second._p);
                    }
                }
                if (!f._cst->is_zero()) {
                    if (f._cst->is_function()) {
                        auto f_cst = *dynamic_pointer_cast<func<T2>>(f._cst);
                        res._cst = multiply(f._cst, f_cst);
                    }
                    else if(f._cst->is_param()) {
                        auto p_cst = *dynamic_pointer_cast<param<T2>>(f._cst);
                        res._cst = multiply(f._cst, p_cst);
                    }
                    else if(f._cst->is_number()) {
                        auto p_cst = *dynamic_pointer_cast<constant<T2>>(f._cst);
                        res._cst = multiply(f._cst, p_cst);
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
            _evaluated = false;
            update_dim(f);
            if (is_constant() && !f.is_constant()) {
                func res(f);
                res += *this;
                return *this = res;
            }
            if (!is_constant() && f.is_constant()) {
                this->add_cst(f);
                update_sign_add(f);                
                _range->first += f._range->first;
                _range->second += f._range->second;
                return *this;
            }
            if (!f.get_cst()->is_zero()) {
                if (f.get_cst()->is_number()) {
                    auto f_cst = dynamic_pointer_cast<constant<T2>>(f.get_cst());
                    add_cst(*f_cst);
                }
                else if (f.get_cst()->is_param()) {
                    auto f_cst = dynamic_pointer_cast<param<T2>>(f.get_cst());
                    add_cst(*f_cst);
                }
                else {
                    auto f_cst = dynamic_pointer_cast<func<T2>>(f.get_cst());
                    add_cst(*f_cst);
                }
                if (_cst->is_function()) {
                    embed(*dynamic_pointer_cast<func_>(_cst));
                }
            }
            for (auto &pair:*f._lterms) {
                auto term = pair.second;
                if (term._coef->is_function()) {
                    auto coef = *dynamic_pointer_cast<func<T2>>(term._coef);
                    term._coef = func(coef).copy();
                }
                else if(term._coef->is_param()) {
                    auto coef = *dynamic_pointer_cast<param<T2>>(term._coef);
                    term._coef = param<type>(coef).copy();
                }
                else if(term._coef->is_number()) {
                    auto coef = *dynamic_pointer_cast<constant<T2>>(term._coef);
                    term._coef = constant<type>(coef.eval()).copy();//TODO if T2==type no need to cast
                }
                this->insert(term);
            }
            for (auto &pair:*f._qterms) {
                auto term = pair.second;
                if (term._coef->is_function()) {
                    auto coef = *dynamic_pointer_cast<func<T2>>(term._coef);
                    term._coef = func(coef).copy();
                }
                else if(term._coef->is_param()) {
                    auto coef = *dynamic_pointer_cast<param<T2>>(term._coef);
                    term._coef = param<type>(coef).copy();
                }
                else if(term._coef->is_number()) {
                    auto coef = *dynamic_pointer_cast<constant<T2>>(term._coef);
                    term._coef = constant<type>(coef.eval()).copy();
                }
                this->insert(term);
            }
            for (auto &pair:*f._pterms) {
                auto term = pair.second;
                if (term._coef->is_function()) {
                    auto coef = *dynamic_pointer_cast<func<T2>>(term._coef);
                    term._coef = func(coef).copy();
                }
                else if(term._coef->is_param()) {
                    auto coef = *dynamic_pointer_cast<param<T2>>(term._coef);
                    term._coef = param<type>(coef).copy();
                }
                else if(term._coef->is_number()) {
                    auto coef = *dynamic_pointer_cast<constant<T2>>(term._coef);
                    term._coef = constant<type>(coef.eval()).copy();
                }
                this->insert(term);
            }
            if (_expr && f.get_expr()) {
                shared_ptr<expr> e1,e2;
                if (_expr->is_uexpr()) {
                    e1 = make_shared<uexpr>(*dynamic_pointer_cast<uexpr>(_expr));
                }
                else {
                    e1 = make_shared<bexpr>(*dynamic_pointer_cast<bexpr>(_expr));
                }
                if (f.get_expr()->is_uexpr()) {
                    auto ue = make_shared<uexpr>(*dynamic_pointer_cast<uexpr>(f.get_expr()));
                    if (ue->_son->is_function()) {
                        auto f = dynamic_pointer_cast<func<T2>>(ue->_son);
                        ue->_son = make_shared<func>(*f);
                    }
                    e2 = ue;
                }
                else {
                    auto bexp = make_shared<bexpr>(*dynamic_pointer_cast<bexpr>(f.get_expr()));
                    if (bexp->_lson->is_function()) {
                        auto f = dynamic_pointer_cast<func<T2>>(bexp->_lson);
                        bexp->_lson = make_shared<func>(*f);
                    }
                    if (bexp->_rson->is_function()) {
                        auto f = dynamic_pointer_cast<func<T2>>(bexp->_rson);
                        bexp->_rson = make_shared<func>(*f);
                    }
                    e2 = bexp;
                }
                _expr = make_shared<bexpr>(bexpr(plus_, e1, e2));
                embed(_expr);
            }
            else if (!_expr && f.get_expr()) {
                if (f.get_expr()->is_uexpr()) {
                    auto ue = make_shared<uexpr>(*dynamic_pointer_cast<uexpr>(f.get_expr()));
                    if (ue->_son->is_function()) {
                        auto f = dynamic_pointer_cast<func<T2>>(ue->_son);
                        ue->_son = make_shared<func>(*f);
                    }
                    _expr = ue;
                }
                else {
                    auto bexp = make_shared<bexpr>(*dynamic_pointer_cast<bexpr>(f.get_expr()));
                    if (bexp->_lson->is_function()) {
                        auto f = dynamic_pointer_cast<func<T2>>(bexp->_lson);
                        bexp->_lson = make_shared<func>(*f);
                    }
                    if (bexp->_rson->is_function()) {
                        auto f = dynamic_pointer_cast<func<T2>>(bexp->_rson);
                        bexp->_rson = make_shared<func>(*f);
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
            else if(_all_convexity!=f._all_convexity){
                _all_convexity = undet_;
            }
            _range->first += f._range->first;
            _range->second += f._range->second;
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
        return func<T1>(f1)*= f2;
    }

    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator*(const func<T1>& f1, const func<T2>& f2){
        return func<T2>(f1)*= f2;
    }
    
    
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator*(const param<T1>& p1, const param<T2>& p2){
        func<T1> res;
        if(p1.is_param() && p2.is_var()){
            res.insert(true,p1,p2);
        }
        else if(p2.is_param() && p1.is_var()){
            res.insert(true,param<T1>(p2),p1);
        }
        else {//Both vars or both params
            res.insert(true,unit<T1>(),p1,p2);
        }
        res.update_dot_dim(p1,p2);
        if(res.has_square()){
            auto signp = p1.get_all_sign();
            if(signp==neg_ || signp==pos_){
                res._all_sign = pos_;
            }
            else {
                res._all_sign = non_neg_;
            }
            res._range->first=zero<T1>().eval();
            if(p1.is_positive()|| p1.is_negative()){
                res._range->first=p1._range->first*p1._range->first;
            }
            res._range->second=p1._range->second*p1._range->second;
        }
        else {
            res._range = get_product_range(p1._range,p2._range);
            res._all_sign = sign_product(p1.get_all_sign(), p2.get_all_sign());
        }
        if(res.is_quadratic()){res.update_quad_convexity();}
        if(p1._is_transposed){
            res._range->first *= p1._dim[0];
            res._range->second *= p1._dim[0];
        }
        return res;
    }

    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator*(const param<T1>& p1, const param<T2>& p2){
        func<T2> res;
        if(p1.is_param() && p2.is_var()){
            res.insert(true,param<T2>(p1),p2);
        }
        else if(p2.is_param() && p1.is_var()){
            res.insert(true,p2,p1);
        }
        else {//Both vars or both params
            res.insert(true,unit<T2>(),p1,p2);
        }
        res.update_dot_dim(p1,p2);
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
                res._range->first=p1._range->first*p1._range->first;
            }
            res._range->second=p1._range->second*p1._range->second;
        }
        else {
            res._range = get_product_range(p1._range,p2._range);
            res._all_sign = sign_product(p1.get_all_sign(), p2.get_all_sign());
        }
        if(res.is_quadratic()){res.update_quad_convexity();
            
        }
        if(p1._is_transposed){
            res._range->first *= p1._dim[0];
            res._range->second *= p1._dim[0];
        }
        return res;
    }

    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) <= sizeof(T1)>::type* = nullptr>
    func<T1> operator+(const param<T1>& p1, const param<T2>& p2){
        func<T1> res;
        res.update_dim(p1,p2);
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
        res._range->first = p1._range->first+p2._range->first;
        res._range->second = p1._range->second+p2._range->second;
        return res;
    }
        
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T1) < sizeof(T2)>::type* = nullptr>
    func<T2> operator+(const param<T1>& p1, const param<T2>& p2){
        func<T2> res;
        res.update_dim(p1,p2);
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
        res._range->first = p1._range->first+p2._range->first;
        res._range->second = p1._range->second+p2._range->second;
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator-(const param<T1>& p1, const param<T2>& p2){
        func<T1> res;
        res.update_dim(p1,p2);
        if(p1.is_param() && p2.is_var()){
            res.insert(false,unit<T1>(),p2);
            res.add_cst(p1);
        }
        else if(p2.is_param() && p1.is_var()){
            res.insert(true,unit<T1>(),p1);
            res.add_cst(-1*param<T1>(p2));
        }
        else {//Both vars or both params
            res.insert(true,unit<T1>(),p1);
            res.insert(false,unit<T1>(),p2);
        }
        res._all_sign = sign_add(p1.get_all_sign(), reverse(p2.get_all_sign()));
        if(res.is_quadratic()){res.update_quad_convexity();}
        res._range->first = p1._range->first-p2._range->second;
        res._range->second = p1._range->second-p2._range->first;
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator-(const param<T1>& p1, const param<T2>& p2){
        func<T2> res;
        res.update_dim(p1,p2);
        if(p1.is_param() && p2.is_var()){
            res.insert(false,unit<T2>(),p2);
            res.add_cst(param<T2>(p1));
        }
        else if(p2.is_param() && p1.is_var()){
            res.insert(true,unit<T2>(),p1);
            auto newp(p2);
            newp.reverse_sign();
            res.add_cst(newp);
        }
        else {//Both vars or both params
            res.insert(true,unit<T2>(),p1);
            res.insert(false,unit<T2>(),p2);
        }
        res._all_sign = sign_add(p1.get_all_sign(), reverse(p2.get_all_sign()));
        if(res.is_quadratic()){res.update_quad_convexity();}
        res._range->first = p1._range->first-p2._range->second;
        res._range->second = p1._range->second-p2._range->first;
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
    func<T1> operator/(const func<T1>& f1, const param<T2>& p2){
        return func<T1>(f1)/=p2;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T1) < sizeof(T2)>::type* = nullptr>
    func<T2> operator/(const func<T1>& f1, const param<T2>& p2){
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
    
    template<class T, typename enable_if<is_arithmetic<T>::value>::type* = nullptr>
    func<T> log(const param<T>& p1){
        if(!p1.is_positive()){
            throw invalid_argument("Calling log() with a non-positive argument");
        }
        func<T> res(uexpr(log_, p1.copy()));
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
        return res;
    }
    
    template<class T1>
    func<T1> exp(const param<T1>& p1){
        func<T1> res(uexpr(exp_, p1.copy()));
        res._all_sign = pos_;
        if (p1.is_var()) {
            res._all_convexity = convex_;
        }
        res._range->first = std::exp(p1._range->first);
        res._range->second = std::exp(p1._range->second);
        return res;
    }
    
    template<class T1>
    func<T1> sqrt(const param<T1>& p1){
        if(!p1.is_non_negative()){
            throw invalid_argument("Calling sqrt() with a negative argument");
        }
        func<T1> res(uexpr(sqrt_, p1.copy()));
        res._all_sign = non_neg_;
        if(p1.is_positive()){
            res._all_sign = pos_;
        }
        if (p1.is_var()) {
            res._all_convexity = concave_;
        }
        res._range->first = std::sqrt(p1._range->first);
        res._range->second = std::sqrt(p1._range->second);
        return res;
    }
    
    template<class T, typename enable_if<is_arithmetic<T>::value>::type* = nullptr>
    pair<Convexity,Sign> cos_sign_curvature(const pair<T,T>& range){
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
    
    template<class T, typename enable_if<is_arithmetic<T>::value>::type* = nullptr>
    func<T> cos(const param<T>& p1){
        func<T> res(uexpr(cos_, p1.copy()));//TODO update ranges and convexity
        auto conv_sign = cos_sign_curvature(*p1._range);
        if (p1.is_var()) {
            res._all_convexity = conv_sign.first;
        }
        res._all_sign = conv_sign.second;
        res._range->first = min(cos(p1._range->first),cos(p1._range->second));
        res._range->second = max(cos(p1._range->first),cos(p1._range->second));
        if(p1._range->first <0 && p1._range->second >0){
            res._range->second = 1;
        }
        if((p1._range->first <-pi && p1._range->second >-pi) || (p1._range->first <pi && p1._range->second >pi)){
            res._range->first = -1;
        }
        return res;
    }
    
    template<class T, typename enable_if<is_arithmetic<T>::value>::type* = nullptr>
    func<T> sin(const param<T>& p1){
        func<T> res(uexpr(sin_, p1.copy()));
        auto shifted_range = *p1._range;
        shifted_range.first += pi/2.;
        shifted_range.second += pi/2.;
        auto conv_sign = cos_sign_curvature(shifted_range);
        if (p1.is_var()) {
            res._all_convexity = conv_sign.first;
        }
        res._all_sign = conv_sign.second;
        res._range->first = min(sin(p1._range->first),sin(p1._range->second));
        res._range->second = max(sin(p1._range->first),sin(p1._range->second));
        if(shifted_range.first <0 && shifted_range.second >0){
            res._range->second = 1;
        }
        if((shifted_range.first <-pi && shifted_range.second >-pi) || (shifted_range.first <pi && shifted_range.second >pi)){
            res._range->first = -1;
        }
        return res;
    }
    
    template<class T, typename enable_if<is_arithmetic<T>::value>::type* = nullptr>
    func<T> tan(const param<T>& p1){
        auto centered_range = *p1._range;
        centered_range.first %= 2*pi;
        centered_range.second %= 2*pi;
        if(centered_range.first<=-pi/2 || centered_range.second>=pi/2){
            throw invalid_argument("Calling tan() with discontinuous domain");
        }
        func<T> res(uexpr(tan_, p1.copy()));
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
        res._range->first = tan(p1._range->first);
        res._range->second = tan(p1._range->second);
        return res;
    }
    
    template<class T1>
    func<T1> ReLU(const param<T1>& p1){
        func<T1> res(uexpr(relu_, p1.copy()));
        if (p1.is_var()) {
            res._all_convexity = convex_;
        }
        res._all_sign = non_neg_;
        if(p1.is_positive()){
            res._all_sign = pos_;
        }
        res._range->first = max(zero<T1>(),p1._range->first);
        res._range->second = max(zero<T1>(),p1._range->second);
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
            res._range->first = min(std::pow(p1._range->first,exp),std::pow(p1._range->second,exp));
            res._range->second = max(std::pow(p1._range->first,exp),std::pow(p1._range->second,exp));
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
            return res;
        }
    }
    
    
    template<class T1>
    func<T1> log(const func<T1>& f){
        func<T1> res(uexpr(log_, f.copy()));
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
        res._range->second = log(f._range->second);
        return res;
    }
    
    template<class T1>
    func<T1> exp(const func<T1>& f){
        func<T1> res(uexpr(exp_, f.copy()));
        res._all_sign = pos_;
        if (f.is_linear()) {
            res._all_convexity = convex_;
        }
        else if(!f.is_constant()){
            res._all_convexity = undet_;
        }
        res._range->first = std::exp(f._range->first);
        res._range->second = std::exp(f._range->second);
        return res;
    }
    
    template<class T1>
    func<T1> sqrt(const func<T1>& f){
        if(!f.is_non_negative()){
            throw invalid_argument("Calling sqrt() with a potentially negative argument");
        }
        func<T1> res(uexpr(sqrt_, f.copy()));
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
        res._range->second = std::sqrt(f._range->second);
        return res;
    }
    
    template<class T, typename enable_if<is_arithmetic<T>::value>::type* = nullptr>
    func<T> cos(const func<T>& f){
        func<T> res(uexpr(cos_, f.copy()));
        auto conv_sign = cos_sign_curvature(*f._range);
        if (f.is_linear()) {
            res._all_convexity = conv_sign.first;
        }
        else if(!f.is_constant()){
            res._all_convexity = undet_;
        }
        res._all_sign = conv_sign.second;
        res._range->first = min(std::cos(f._range->first),std::cos(f._range->second));
        res._range->second = max(std::cos(f._range->first),std::cos(f._range->second));
        if(f._range->first <0 && f._range->second >0){
            res._range->second = 1;
        }
        if((f._range->first <-pi && f._range->second >-pi) || (f._range->first <pi && f._range->second >pi)){
            res._range->first = -1;
        }
        return res;
    }
    
    template<class T, typename enable_if<is_arithmetic<T>::value>::type* = nullptr>
    func<T> sin(const func<T>& f){
        func<T> res(uexpr(sin_, f.copy()));
        auto shifted_range = *f._range;
        shifted_range.first += pi/2.;
        shifted_range.second += pi/2.;
        auto conv_sign = cos_sign_curvature(shifted_range);
        if (f.is_linear()) {
            res._all_convexity = conv_sign.first;
        }
        else if(!f.is_constant()){
            res._all_convexity = undet_;
        }
        res._all_sign = conv_sign.second;
        res._range->first = min(std::sin(f._range->first),std::sin(f._range->second));
        res._range->second = max(std::sin(f._range->first),std::sin(f._range->second));
        if(shifted_range.first <0 && shifted_range.second >0){
            res._range->second = 1;
        }
        if((shifted_range.first <-pi && shifted_range.second >-pi) || (shifted_range.first <pi && shifted_range.second >pi)){
            res._range->first = -1;
        }
        return res;
    }
    
    template<class T, typename enable_if<is_arithmetic<T>::value>::type* = nullptr>
    func<T> tan(const func<T>& f){
        func<T> res(uexpr(tan_, f.copy()));
        auto centered_range = *f._range;
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
        res._range->first = tan(f._range->first);
        res._range->second = tan(f._range->second);
        return res;
    }
    
    template<class T1>
    func<T1> ReLU(const func<T1>& f){
        func<T1> res(uexpr(relu_, f.copy()));
        if (f.is_linear()) {
            res._all_convexity = convex_;
        }
        else if(!f.is_constant()){
            res._all_convexity = undet_;
        }
        res._all_sign = non_neg_;
        if(f.is_positive()){
            res._all_sign = pos_;
        }
        res._range->first = max(zero<T1>(),f._range->first);
        res._range->second = max(zero<T1>(),f._range->second);
        return res;
    }
    
    template<class T, typename enable_if<is_arithmetic<T>::value>::type* = nullptr>
    func<T> pow(const func<T>& f, int exp){
        if(exp<0){
            return func<T>(bexpr(power_, f.copy(), make_shared<constant<int>>(exp)));
        }
        if(exp==0){
            return func<T>();
        }
        if(exp==1){
            return f;
        }
        if(exp==2){
            return f*f;
        }
        else {
            func<T> res(f);
            for (int i = 1; i < exp; i++) {
                res *= f;
            }
            res._range->first = min(std::pow(f._range->first,exp),std::pow(f._range->second,exp));
            res._range->second = max(std::pow(f._range->first,exp),std::pow(f._range->second,exp));
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
            return res;
        }
    }
    
    

    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator+(const constant<T1>& p, const param<T2>& v){
        func<T1> res(v);
        res.add_cst(p);
        res._range->first = p.eval()+v._range->first;
        res._range->second = p.evla()+v._range->second;
        res.update_all_sign();
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator+(const constant<T1>& p, const param<T2>& v){
        func<T2> res(v);
        res.add_cst(p);
        res._range->first = p.eval()+v._range->first;
        res._range->second = p.eval()+v._range->second;
        res.update_all_sign();
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator-(const constant<T1>& p, const param<T2>& v){
        func<T1> res(v);
        res.reverse_sign();
        res.add_cst(p);
        res._range->first = p.eval()-v._range->second;
        res._range->second = p.eval()-v._range->first;
        res.update_all_sign();
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator-(const constant<T1>& p, const param<T2>& v){
        func<T2> res(v);
        res.reverse_sign();
        res.add_cst(p);
        res._range->first = p.eval()-v._range->second;
        res._range->second = p.eval()-v._range->first;
        res.update_all_sign();
        return res;
    }

        
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator+(const param<T1>& v, const constant<T2>& p){
        func<T1> res(v);
        res.add_cst(p);
        res._range->first = p.eval()+v._range->first;
        res._range->second = p.evla()+v._range->second;
        res.update_all_sign();
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator+(const param<T1>& v, const constant<T2>& p){
        func<T2> res(v);
        res.add_cst(p);
        res._range->first = p.eval()+v._range->first;
        res._range->second = p.eval()+v._range->second;
        res.update_all_sign();
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator-(const param<T1>& v, const constant<T2>& p){
        func<T1> res(v);
        auto newp = constant<T1>(p);
        newp.reverse_sign();
        res.add_cst(newp);
        res._range->first = v._range->first-p.eval();
        res._range->second = v._range->second-p.eval();
        res.update_all_sign();
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator-(const param<T1>& v, const constant<T2>& p){
        func<T2> res(v);
        auto newp(p);
        newp.reverse_sign();
        res.add_cst(newp);
        res._range->first = v._range->first-p.eval();
        res._range->second = v._range->second-p.eval();
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
        res -= f;
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
        res *= f;
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator*(const param<T1>& v, const func<T2>& f){
        func<T2> res(v);
        res *= f;
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator*(const func<T1>& f, const param<T2>& v){
        func<T1> res(f);
        res *= func<T1>(v);
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator*(const func<T1>& f, const param<T2>& v){
        func<T2> res(f);
        res *= func<T2>(v);
        return res;
    }
        
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator*(const constant<T1>& c, const param<T2>& p){
        func<T1> res;
        res.update_dot_dim(c,p);
        res.insert(true,c,p);
        res._range->first = min(c.eval()*p._range->first, c.eval()*p._range->second);
        res._range->second = max(c.eval()*p._range->first, c.eval()*p._range->second);
        res.update_all_sign();
        if(c._is_transposed){
            res._range->first *= p._dim[0];
            res._range->second *= p._dim[0];
        }
        return res;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator*(const constant<T1>& c, const param<T2>& p){
        func<T2> res;
        res.update_dot_dim(c,p);
        res.insert(true,constant<T2>(c),p);
        res._range->first = min(c.eval()*p._range->first, c.eval()*p._range->second);
        res._range->second = max(c.eval()*p._range->first, c.eval()*p._range->second);
        res.update_all_sign();
        if(c._is_transposed){
            res._range->first *= p._dim[0];
            res._range->second *= p._dim[0];
        }
        return res;
    }
    
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator*(const param<T1>& p, const constant<T2>& c){
        func<T1> res;
        res._range->first = min(c.eval()*p._range->first, c.eval()*p._range->second);
        res._range->second = max(c.eval()*p._range->first, c.eval()*p._range->second);
        res.update_all_sign();
        res.update_dot_dim(p,c);
        if(p._is_transposed){
            res.insert(true,constant<T1>(c).tr(),p.tr());
            res._range->first *= p._dim[0];
            res._range->second *= p._dim[0];
        }
        else {
            res.insert(true,constant<T1>(c),p);
        }
        return res;
    }
        
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator*(const param<T1>& p, const constant<T2>& c){
        func<T2> res;
        res._range->first = min(c.eval()*p._range->first, c.eval()*p._range->second);
        res._range->second = max(c.eval()*p._range->first, c.eval()*p._range->second);
        res.update_all_sign();
        res.update_dot_dim(p,c);
        if(p._is_transposed){
            res.insert(true,c.tr(),p.tr());
            res._range->first *= p._dim[0];
            res._range->second *= p._dim[0];
        }
        else {
            res.insert(true,c,p);
        }
        return res;
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
        return constant<T1>(p)*v;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator*(const param<T1>& v, T2 p){
        return constant<T2>(p)*v;
    }
    
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator*(T1 p, const func<T2>& f){
        return func<T1>(p) * f;
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator*(T1 p, const func<T2>& f){
        func<T2> res(p);
        return res*= f;
    }
    
    
    template<class T1,class T2, typename enable_if<is_convertible<T2, T1>::value && sizeof(T2) < sizeof(T1)>::type* = nullptr>
    func<T1> operator*(const func<T1>& f, T2 p){
        return f * func<T1>(p);
    }
    
    template<class T1,class T2, typename enable_if<is_convertible<T1, T2>::value && sizeof(T2) >= sizeof(T1)>::type* = nullptr>
    func<T2> operator*(const func<T1>& f, T2 p){
        return f * func<T2>(p);
    }
    
}



#endif /* func_h */

