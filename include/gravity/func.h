//
//  func.h
//  Gravity
//
//  Created by Hijazi, Hassan on 21/10/16.
//
//

#ifndef func_h
#define func_h

#include <gravity/var.h>
#include <stdio.h>
#include <map>
#include <iterator>
#include <queue>
#include <list>
#include <limits>
#include <set>

using namespace std;

namespace gravity {
    class func_;
    void reverse_sign(constant_* c); /**< Reverses the sign of the constant. */
    constant_* copy(const constant_& c2); /**< Copy c2 into a new constant_* detecting the right class, i.e., constant<>, param<>, var<> or function */
    constant_* copy(constant_&& c2); /**< Copy c2 into a new constant_* detecting the right class, i.e., constant<>, param<>, var<> or function */
    bool equals(const constant_* c1, const constant_* c2);
    void poly_set_val(unsigned i, double val, param_* p);

    template<typename type> type eval(ind i, const constant_* c1);
    template<typename type> type eval(const constant_* c1);
    template<typename type> param<type> vec(type t, int n); /**< Returns a vector of dimension n and storing the value t */


    /** Backbone class for unary and binary expressions. */

    class expr: public constant_{
    protected:
    public:
        double                                 _coef = 1.; /**< coefficient multpying the expression */
        string                                 _to_str; /**< A string representation of the expression */
        
        string get_str();
        double eval(size_t i) const;
        double eval(size_t i, size_t j) const;
        func_ get_derivative(const param_ &v) const;
        void reset_val();
        void untranspose();
        virtual ~expr(){};
    };


    /** Class uexpr (unary expression), stores a unary expression tree. */
    class uexpr: public expr{
        
    public:
        OperatorType                _otype;
        shared_ptr<func_>           _son;
        
        uexpr();
        uexpr(const uexpr& exp);
        uexpr(uexpr&& exp);
        uexpr(OperatorType ot, shared_ptr<func_> son);
        uexpr& operator=(const uexpr& e);
        uexpr& operator=(uexpr&& e);
        
        
        ~uexpr(){};
                
        void reset_val();
        
        void reset(){
            _son = nullptr;
            _otype = id_;
            _to_str = "noname";
            _coef = 1.;
        };
        
        
        OperatorType get_otype() const{
            return _otype;
        };
        
        /** Operators */
        
        void untranspose();
        
        bool operator==(const uexpr &c)const;
        
        bool operator!=(const uexpr& c) const{
            return !(*this==c);
        };
        
        Sign get_all_sign() const{
            return unknown_;// TO UPDATE
        }
        
        double eval(size_t i) const;
        double eval(size_t i, size_t j) const;
        
        double eval() const{
            return eval(0);
        }
        
        string to_str() const;
        void print(bool endline = true) const;
        func_ get_derivative(const param_ &v) const;
        vector<shared_ptr<param_>> get_nl_vars() const;
    };


    class bexpr: public expr{
    private:
        
    public:
        OperatorType    _otype;
        shared_ptr<func_>      _lson;
        shared_ptr<func_>      _rson;
        
        bexpr();
        
        bexpr(OperatorType otype, shared_ptr<func_> lson, shared_ptr<func_> rson);
        
        bexpr(const bexpr& exp);
        
        bexpr(bexpr&& exp);
        
        bexpr& operator=(const bexpr& e);
        
        bexpr& operator=(bexpr&& e);
        
        ~bexpr(){}
        
        void reset_val();
        
        void reset(){
            _otype = id_;
            _to_str = "noname";
            _coef = 1.;            
            _lson = nullptr;
            _rson = nullptr;
        };
        
        void untranspose();
        
        shared_ptr<func_> get_lson() const{
            return _lson;
        };
        
        shared_ptr<func_> get_rson() const{
            return _rson;
        };
        
        void set_lson(shared_ptr<func_> c){
            _lson = c;
        };
        
        void set_rson(shared_ptr<func_> c){
            _rson = c;
        };
        
        OperatorType get_otype() const {
            return _otype;
        };
        
        
        bool operator==(const bexpr &c)const;
        
        bool operator!=(const bexpr& c) const{
            return !(*this==c);
        };
        
        Sign get_all_sign() const{
            return unknown_; // TO UPDATE
        }
        
        
        template<typename other_type> bexpr& operator +=(const other_type& v);
        template<typename other_type> bexpr& operator -=(const other_type& v);
        template<typename other_type> bexpr& operator *=(const other_type& v);
        template<typename other_type> bexpr& operator /=(const other_type& v);
        
        
        string to_str() const;
        
        void print(bool endline = true) const;
        
        void print_tree() const;
        
        double eval(ind i) const;
        double eval(ind i, ind j) const;
        
        func_ get_derivative(const param_ &v) const;
        
        vector<shared_ptr<param_>> get_nl_vars() const;
        
    };



    /** A class to represent a linear term, e.g. 2x. */
    class lterm{
        
    public:
        constant_*              _coef;
        param_*                 _p;
        bool                    _sign = true; /**< True if +, flase if - */    
        
        lterm(){
            _coef = nullptr; // coefficent
            _p = nullptr; // terms. 
        }
        
        lterm(lterm&& t){
            _coef = t._coef;
            t._coef = nullptr;
            _p = t._p;
            t._p = nullptr;
            _sign = t._sign;
        };
        
        
        lterm(param_* p):lterm(true,p){
        };
        
        
        lterm(bool sign, param_* p){
            _coef = new constant<double>(1);
            _p = p;
            _sign = sign;
        };
        
        lterm(constant_* coef, param_* p):lterm(true,coef,p){};
        
        lterm(bool sign, constant_* coef, param_* p);
        
        void reverse_sign() {
            _sign = ! _sign;
        }
        
        double eval(size_t i) const;
        double eval(size_t i, size_t j) const;
        
        ~lterm(){
            delete _coef;
        };
        
        bool operator==(const lterm& l) const;
        bool operator!=(const lterm& l) const;
        
        lterm& operator=(const lterm& l);
        lterm& operator=(lterm&& l);
        
        string to_str(int ind) const;
        void print(int ind) const;
    };


    /** A class to represent a quadratic term, e.g. 2xy or 3x^2. */
    class qterm{
        
    public:
        constant_*                  _coef;
        pair<param_*,param_*>*      _p;
        bool                        _sign = true; /**< True if +, flase if - */
        
        qterm(){
            _coef = nullptr;
            _p = nullptr;
        }
        
        qterm(qterm&& t){
            _coef = t._coef;
            t._coef = nullptr;
            _p = t._p;
            t._p = nullptr;
            _sign = t._sign;
        };
        
        
        qterm(param_* p1, param_* p2):qterm(true, p1, p2){};
        
        
        qterm(constant_* coef, param_* p1, param_* p2):qterm(true, coef, p1, p2){};
        
        qterm(bool sign, param_* p1, param_* p2):qterm(true, new constant<double>(1), p1, p2){};
        
        qterm(bool sign, constant_* coef, param_* p1, param_* p2){
            _coef = coef;
            _p = new pair<param_*, param_*>(make_pair(p1,p2));
            _sign = sign;
            if (coef->_is_transposed && p1->_is_transposed) {
                throw invalid_argument("Check the transpose operator, there seems to be a dimension issue\n");
            }
//            if (p1->_is_transposed) {
//                if (p1->get_dim() != p2->get_dim()) {
//                    throw invalid_argument("Check the transpose operator, there seems to be a dimension issue\n");
//                }
//            }
            if (p2->_is_transposed) {
                throw invalid_argument("Check the transpose operator, there seems to be a dimension issue\n");
            }
        };
        
        void reverse_sign() {
            _sign = ! _sign;
        }
        
        double eval(size_t i) const;
        double eval(size_t i, size_t j) const;
        
        ~qterm(){
            delete _coef;
            if (_p) {
                delete _p;
            }
        };
        
        
        bool operator==(const qterm& l) const;
        bool operator!=(const qterm& l) const;
        
        qterm& operator=(const qterm& l);
        qterm& operator=(qterm&& l);
        
        string to_str(int ind) const;
        void print(int ind) const;
    };


    /** A class to represent a polynomial term, e.g. 2xyz or 3x^2y. */
    class pterm{
        
    public:
        constant_*                      _coef;
        list<pair<param_*, int>>*       _l; /**< A polynomial term is represented as a list of pairs <param_*,int> where the first element points to the parameter and the second indicates the exponent */
        bool                            _sign = true; /**< True if +, flase if - */
        
        pterm(){
            _coef = nullptr;
            _l = nullptr;
        }
        
        
        pterm(bool sign, constant_* coef, param_* p, int exp){
            _coef = coef;
            if (coef->_is_transposed && p->_is_transposed) {
                throw invalid_argument("Check the transpose operator, there seems to be a dimension issue\n");
            }
            if (p->_is_transposed) {
                throw invalid_argument("Check the transpose operator, there seems to be a dimension issue\n");
            }
            _l = new list<pair<param_*, int>>();
            _l->push_back(make_pair<>(p, exp));
            _sign = sign;
        };
        
        
        pterm(pterm&& t){
            _coef = t._coef;
            t._coef = nullptr;
            _l = t._l;
            t._l = nullptr;
            _sign = t._sign;
        };
        
        pterm(bool sign, constant_* coef, list<pair<param_*, int>>* l){
            param_* p1 = nullptr;
            param_* p2 = nullptr;
            for (auto it = l->begin(); it != l->end(); it++){
                p1 = it->first;
                if (p1->_is_transposed && next(it)!=l->end()) {
                    p2 = next(it)->first;
                    if (p1->get_dim() != p2->get_dim()) {
                        throw invalid_argument("Check the transpose operator, there seems to be a dimension issue\n");
                    }
                }
                if (p1->_is_transposed && next(it)==l->end()) {
                    throw invalid_argument("Check the transpose operator, there seems to be a dimension issue\n");
                }
            }
            _coef = coef;
            _l = l;
            _sign = sign;
        };

        
        void reverse_sign() {
            _sign = ! _sign;
        }
        
        double eval(size_t i) const;
        double eval(size_t i, size_t j) const;
        
        ~pterm(){
            delete _coef;
            if (_l) {
                delete _l;
            }
        };
        
        bool operator==(const pterm& l) const;
        bool operator!=(const pterm& l) const;
        
        pterm& operator=(const pterm& l);
        pterm& operator=(pterm&& l);
        
        string to_str(int ind) const;
        void print(int ind) const;
        
    };

    /** Backbone class for function */
    class func_ : public constant_{
    private:
        
        func_* compute_derivative(const param_& v);  /**< Computes and stores the derivative of f with respect to variable v. Returns a pointer to the stored function. */
        
    protected:
        
        FType                                  _ftype = const_; /**< Function type, e.g., constant, linear, quadratic... >>**/
        NType                                  _return_type = integer_; /**< Return type, e.g., bool, integer, complex... >>**/

        map<string, pair<shared_ptr<param_>, int>>*       _params = nullptr;/**< Set of parameters in current function, stored as a map <parameter name, <paramter pointer, number of times it appears in function>>**/
        map<string, pair<shared_ptr<param_>, int>>*       _vars = nullptr;/**< Set of variables in current function, stored as a map <variable name, <variable pointer, number of times it appears in function>>**/
                
        constant_*                             _cst = nullptr;/**< Constant part of the function */
        map<string, lterm>*                    _lterms = nullptr; /**< Set of linear terms, stored as a map <string describing term, term>. */
        map<string, qterm>*                    _qterms = nullptr; /**< Set of quadratic terms, stored as a map <string describing term, term>.  */
        map<string, pterm>*                    _pterms = nullptr; /**< Set of polynomial terms, stored as a map <string describing term, term>.  */
        shared_ptr<expr>                       _expr = nullptr; /**< Nonlinear part of the function, this points to the root node in _DAG */
        map<string, expr*>*                    _DAG = nullptr; /**< Map of experssions stored in the expression tree (a Directed Acyclic Graph) */
        deque<shared_ptr<expr>>*               _queue = nullptr; /**< A queue storing the expression tree from the leaves to the root (the root is stored at the bottom of the queue)*/
        shared_ptr<map<unique_id,shared_ptr<func_>>>      _dfdx;/**< A map storing the first derivatives of f per variable name*/
        Convexity                              _all_convexity; /**< If all instances of this function have the same convexity type, it stores it here, i.e. linear, convex, concave, otherwise it stores unknown. >>**/
        Sign                                   _all_sign; /**< If all instances of this function have the same sign, it stores it here, otherwise it stores unknown. >>**/
        pair<constant_*, constant_*>*          _all_range = nullptr; /**< Range of the return value considering all instances of the current function. >>**/

        vector<Convexity>*                     _convexity = nullptr; /**< Vector of convexity types, i.e., linear, convex, concave or unknown. This is a vector since a function can have multiple instances (different constants coefficients, and bounds, but same structure) >>**/
        vector<Sign>*                          _sign = nullptr; /**< vector storing the sign of return value if known. >>**/
        vector<pair<constant_*, constant_*>>*  _range = nullptr; /**< Bounds of the return value if known. >>**/
        
        map<unsigned, set<unsigned>>            _hess_link; /**< Set of variables linked to one another in the hessian, indexed by variable ids  */
        
        size_t                                 _nb_vars = 0; /**< Number of variables */
                                                                   
        size_t                                 _nnz_j = 0; /**< Number of nonzeros in the Jacobian **/
        size_t                                 _nnz_h = 0; /**< Number of nonzeros in the Jacobian **/
        
    public:
        size_t                                 _nb_instances = 1; /**< Number of different instances this constraint has (different indices, constant coefficients and bounds, but same structure).>>**/

        bool                                   _is_constraint = false;
        bool                                   _embedded = false; /**< If the function is embedded in
                                                                   a mathematical model or in another function, this is used for memory management. >>**/
        bool                                   _evaluated = false;/**< If the function has already been evaluated, useful for constant funcs */
        shared_ptr<vector<double>>             _val;
        shared_ptr<vector<unsigned>>           _ids = nullptr; /*<<A vector storing all the indices this constraint has in the order they were created */
        string                                 _to_str;
        func_();
        
        func_(const constant_& c);
        
        func_(constant_&& c);
        
        func_(const func_& f); /**< Copy constructor */
        
        func_(func_&& f); /**< Move constructor */

        ~func_();

        map<unsigned, set<unsigned>>& get_hess_link() { return _hess_link;};
        map<string, pair<shared_ptr<param_>, int>>& get_vars() { return *_vars;};
        map<string, pair<shared_ptr<param_>, int>>& get_params() { return *_params;};
        
        bool has_var(const param_& v) const;
        bool has_var(const string& name) const;
        
        void reset_val();
        void replace(param_* p, const param_& p_new);
        
        bool insert(bool sign, const constant_& coef, const param_& p);/**< Adds coef*p to the function. Returns true if added new term, false if only updated coef of p */
        bool insert(bool sign, const constant_& coef, const param_& p1, const param_& p2);/**< Adds coef*p1*p2 to the function. Returns true if added new term, false if only updated coef of p1*p2 */
        bool insert(bool sign, const constant_& coef, const list<pair<param_*, int>>& l);/**< Adds polynomial term to the function. Returns true if added new term, false if only updated corresponding coef */
       
        func_ tr() const {
            auto f = func_(*this);
            f.transpose();
            return f;
        }

        void transpose(){
            _is_transposed = !_is_transposed;
            if (!_is_vector) {
                _is_vector = true;
            }
            if (_is_matrix) {
                auto temp = _dim[0];
                _dim[0] = _dim[1];
                _dim[1] = temp;
//                for (auto &p_t:*_params) {
//                    p_t.second.first->transpose();
//                }
            }
        }
        void set_first_derivative(const param_& v, func_&& f){
            DebugOn(f.to_str()<<endl);
            (*_dfdx)[v._unique_id] = make_shared<func_>(move(f));
        }
        
        void set_second_derivative(const param_& v1, const param_& v2, func_&& f){
            DebugOn(f.to_str()<<endl);
            (*_dfdx)[v1._unique_id]->_dfdx->insert(make_pair<>(v2._unique_id, make_shared<func_>(move(f))));
        }
        
        void untranspose(){/**< Untranspose the output of the current function */
            _is_transposed = false;
            _is_vector = false;
            auto vars_cpy = *_vars;
            for (auto &vp:*_vars) {
                vars_cpy.erase(vp.first);
                vp.second.first->untranspose();
                _nb_instances = max(_nb_instances, vp.second.first->get_dim());
                vars_cpy[vp.second.first->get_name()]= make_pair<>(vp.second.first, vp.second.second);
            }
            *_vars = move(vars_cpy);
            auto params_cpy = *_params;
            for (auto &vp:*_params) {
                params_cpy.erase(vp.first);
                vp.second.first->untranspose();
                _nb_instances = max(_nb_instances, vp.second.first->get_dim());
                params_cpy[vp.second.first->get_name()]= make_pair<>(vp.second.first, vp.second.second);
            }
            *_params = move(params_cpy);
            _val->resize(max(_val->size(),_nb_instances));
            if(_expr){
                _expr->untranspose();
//                embed(_expr);
            }
            
        }
        
        template<typename Tobj>
        func_& in(const vector<Tobj*>& vec);
        
        void insert(const lterm& term);
        
        void insert(const qterm& term);
        
        void insert(const pterm& term);
        
//        void insert(expr& e);
        void update_to_str();
        size_t get_nb_vars() const;
        
        constant_* get_cst();
        
        shared_ptr<param_> get_var(string name);
        
        shared_ptr<param_> get_param(string name);
        
        void add_var(shared_ptr<param_> v, int nb = 1);/**< Inserts the variable in this function input list. nb represents the number of occurences v has. WARNING: Assumes that v has not been added previousely!*/
        
        
        void add_param(shared_ptr<param_> p, int nb = 1);/**< Inserts the parameter in this function input list. WARNING: Assumes that p has not been added previousely!*/
        
        
        
        void delete_var(const string& vid);
        
        void delete_param(const string& vid);
        

        void update_nb_instances(const lterm& term){
            if ((term._coef->_is_matrix || term._coef->_is_vector) && (term._p->_is_matrix || term._p->_is_vector)) {
                if (term._coef->_is_matrix) {
                    if (term._coef->_is_transposed) {//TODO check matrix dims match
                        _nb_instances = max(_nb_instances, term._coef->get_dim(1));
                    }
                    else {
                        _nb_instances = max(_nb_instances, term._coef->get_dim(0));
                    }
                    _dim[0] = _nb_instances;
                    if (term._p->_is_matrix) {//Matrix product
                        _dim[1] = _nb_instances;
                    }
                    else {
                        _dim[1] = 1;
                    }
                    _is_vector = true;
                    _is_matrix = false;
                }
                else {//_coef is a transposed vector
                    if (!term._p->_is_transposed && !term._p->_is_matrix) {//_coef is a transposed vector at this stage, if _p is not transposed, we have a scalar product
                        _dim.resize(1);
                        _dim[0] = 1;
                        _nb_instances = 1;
                        _is_vector = false;
                        _is_matrix = false;
                    }
                    else {//_p is either transposed or a matrix at this stage
                        _dim.resize(1);
                        _dim[0] = term._coef->_dim[0];
                        _nb_instances = _dim[0];
                        _is_transposed = true;
                        _is_vector = true;
                        _is_matrix = false;
                    }
                }
                _is_matrix = term._coef->_is_matrix && term._p->_is_matrix;
            }
            else {
                _nb_instances = max(_nb_instances, term._p->get_dim());
            }
            _val->resize(_nb_instances);
        }
        
        void update_nb_instances(const qterm& term){
            if ((term._p->first->_is_matrix || term._p->first->_is_vector) && (term._p->second->_is_matrix || term._p->second->_is_vector)) {
                if (term._p->first->_is_matrix) {
                    if (term._p->first->_is_transposed) {//TODO check matrix dims match
                        _nb_instances = max(_nb_instances, term._coef->get_dim(1));
                    }
                    else {
                        _nb_instances = max(_nb_instances, term._coef->get_dim(0));
                    }
                    _dim[0] = _nb_instances;
                    if (term._p->second->_is_matrix) {//Matrix product
                        _dim[1] = _nb_instances;
                    }
                    else {
                        _dim[1] = 1;
                    }
                    _is_vector = true;
                    _is_matrix = false;
                }
                else {//_coef is a transposed vector
                    if (!term._p->second->_is_transposed && !term._p->second->_is_matrix) {//_coef is a transposed vector at this stage, if _p is not transposed, we have a scalar product
                        _dim.resize(1);
                        _dim[0] = 1;
                        _nb_instances = 1;
                        _is_vector = false;
                        _is_matrix = false;
                    }
                    else {//_p is either transposed or a matrix at this stage
                        _dim.resize(1);
                        _dim[0] = term._coef->_dim[0];
                        _nb_instances = _dim[0];
                        _is_transposed = true;
                        _is_vector = true;
                        _is_matrix = false;
                    }
                }
//                _is_matrix = term._coef->_is_matrix && term._p->_is_matrix;
            }
            else {
                _nb_instances = max(_nb_instances, term._p->second->get_dim());
            }
            _val->resize(_nb_instances);
        }
        
        int nb_occ_var(string name) const;/**< Returns the number of occurences the variable has in this function. */
        
        int nb_occ_param(string name) const;/**< Returns the number of occurences the parameter has in this function. */
        
        void incr_occ_var(string str);/**< Increases the number of occurences the variable has in this function. */
        
        void incr_occ_param(string str);/**< Increases the number of occurences the parameter has in this function. */
        
        void decr_occ_var(string str, int nb=1);/**< Decreases the number of occurences the variable has in this function by nb. */

        void decr_occ_param(string str, int nb=1);/**< Decreases the number of occurences the parameter has in this function by nb. */
        
        
        pair<ind,func_*> operator[](ind i);
        
        bool is_convex(int idx=0) const;
        bool is_concave(int idx=0) const;
        bool is_number() const;
        bool is_constant() const;
        bool is_linear() const;
        bool is_quadratic() const;
        bool is_polynomial() const;
        bool is_nonlinear() const;
        bool is_zero() const;/*<< A function is zero if it is constant and equals zero or if it is a sum of zero valued parameters */
        bool is_unit() const;
        bool is_transposed() const;
        FType get_ftype() const;
        void embed(func_& f);
        void embed(shared_ptr<expr> e);
        
        void reset();
        
        void reverse_sign(); /*<< Reverse the sign of all terms in the function */
        
        void reverse_convexity();
        
        void untranspose_derivatives(){
            for (auto &fp:*_dfdx) {
                auto df = fp.second;
                df->untranspose();
                df->_nb_instances = max(df->_nb_instances, _nb_instances);
                df->_val->resize(max(df->_val->size(),_nb_instances));
                df->untranspose_derivatives();
            }
        }
        
        
        
        func_ vec() const {
            auto f = func_(*this);
            f._is_vector = true;
//            auto vars_cpy = *f._vars;
//            for (auto &vp:*f._vars) {
//                vars_cpy.erase(vp.first);
//                vp.second.first->_is_vector = true;
//                f._nb_instances = max(f._nb_instances, vp.second.first->get_nb_instances());
//                vars_cpy[vp.second.first->get_name()]= make_pair<>(vp.second.first, vp.second.second);
//            }
//            *f._vars = move(vars_cpy);
//            if (f._expr) {
//                f._expr->_is_vector = true;
//                if(f._expr->is_uexpr()){
//                    auto ue = dynamic_pointer_cast<uexpr>(f._expr);
//                    auto s = (func_*)(ue->_son.get());
//                    *s = s->vec();
//                }
//                else {
//                    auto be = dynamic_pointer_cast<bexpr>(f._expr);
//                    auto fl = (func_*)(be->_lson.get());
//                    *fl = fl->vec();
//                    auto fr = (func_*)(be->_rson.get());
//                    *fr = fr->vec();
//                }
//            }
            return f;
        }
        
        func_& operator=(const func_& f);
        
        func_& operator=(func_&& f);
        
        bool operator==(const func_& f) const;
        
        bool operator!=(const func_& f) const;
        

        func_& operator+=(const constant_& f);
        func_& operator-=(const constant_& f);
        func_& operator*=(const constant_& f);
        func_& operator/=(const constant_& f);
        
        
        friend func_ cos(const constant_& c);
        friend func_ cos(constant_&& c);
        
        friend func_ sin(const constant_& c);
        friend func_ sin(constant_&& c);
        
        
        friend func_ sqrt(const constant_& c);
        friend func_ sqrt(constant_&& c);
        
        friend func_ expo(const constant_& c);
        friend func_ expo(constant_&& c);
        
        friend func_ log(const constant_& c);
        friend func_ log(constant_&& c);

        
        template<class T, class = typename enable_if<is_arithmetic<T>::value>::type> func_& operator+=(T c){
            return *this += constant<T>(c);
        };
        
        template<class T, class = typename enable_if<std::is_arithmetic<T>::value>::type> func_& operator-=(T c){
            return *this -= constant<T>(c);
        };
        
        template<class T, class = typename enable_if<std::is_arithmetic<T>::value>::type> func_& operator*=(T c){
            return *this *= constant<T>(c);
        };
        
        template<class T, class = typename enable_if<std::is_arithmetic<T>::value>::type> func_& operator/=(T c){
            return *this /= constant<T>(c);
        };
        
        qterm* get_square(param_* p); /**< Returns the quadratic term containing a square of p or nullptr if none exists. **/
      
        func_ get_outer_app(); /**< Returns an outer-approximation of the function using the current value of the variables **/
        
        Sign get_all_sign() const;
        Sign get_sign(int idx=0) const;
        Sign get_all_sign(const lterm& l);
        Sign get_all_sign(const qterm& l);
        Sign get_all_sign(const pterm& l);

        
        Convexity get_convexity(const qterm& q);
        
        void update_sign(const constant_& c);
        
        void update_sign(const lterm& l);
        
        void update_sign(const qterm& q);
        
        void update_sign(const pterm& q);
        
        void update_convexity(const qterm& q);
        
        void update_convexity();
        
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
        
        shared_ptr<map<unique_id,shared_ptr<func_>>> get_dfdx() const{
            return _dfdx;
        };
        
        shared_ptr<func_> get_stored_derivative(const unique_id& vid) const; /**< Returns the stored derivative with respect to variable v. */
        
        func_ get_derivative(const param_& v) const; /**< Computes and returns the derivative with respect to variable v. */
        
        
        void compute_derivatives(); /**< Computes and stores the derivative of f with respect to all variables. */
        
        void update_sign();
        
        double get_val(size_t inst) const;
        double get_val(size_t i, size_t j) const;
        double eval(size_t i);
        double eval(size_t i, size_t j);
//        double force_eval(size_t i);
        double eval(){ return eval(0);};
        string to_str(bool display_input=false) const;
        void print(bool endline=false, bool display_input=false) const;
        
    };



    bool is_indexed(const constant_* c);

    size_t get_poly_id(const constant_* c);

    size_t get_poly_id_inst(const constant_* c, unsigned inst = 0);

    double poly_eval(const constant_* c, size_t i=0);
    double poly_eval(const constant_* c, size_t i, size_t j);


    void poly_print(const constant_* c);

    string poly_to_str(const constant_* c);


    func_ operator+(const constant_& c1, const constant_& c2);
//    func_ operator+(func_&& f, const constant_& c);
    //func_ operator+(const constant_& c, func_&& f);

    template<class T, class = typename enable_if<std::is_arithmetic<T>::value>::type> func_ operator+(func_&& f, T c){
        return f += c;
    };

    template<class T, class = typename enable_if<std::is_arithmetic<T>::value>::type> func_ operator+(T c, func_&& f){
        return f += c;
    };

    template<class T, class = typename enable_if<std::is_arithmetic<T>::value>::type> func_ operator-(func_&& f, T c){
        return f -= c;
    };

    template<class T, class = typename enable_if<std::is_arithmetic<T>::value>::type> func_ operator-(T c, func_&& f){
        return (f *= -1) += c;
    };

    template<class T, class = typename enable_if<std::is_arithmetic<T>::value>::type> func_ operator*(func_&& f, T c){
        return f *= c;
    };

    template<class T, class = typename enable_if<std::is_arithmetic<T>::value>::type> func_ operator*(T c, func_&& f){
        return f *= c;
    };


    template<class T, class = typename enable_if<std::is_arithmetic<T>::value>::type> func_ operator/(func_&& f, T c){
        return f /= c;
    };

    template<class T, class = typename enable_if<std::is_arithmetic<T>::value>::type> func_ operator+(const constant_& c1, T c2){
        return func_(c1) += c2;
    };

    template<class T, class = typename enable_if<std::is_arithmetic<T>::value>::type> func_ operator+(T c2, const constant_& c1){
        return func_(c1) += c2;
    };

    template<class T, class = typename enable_if<std::is_arithmetic<T>::value>::type> func_ operator-(const constant_& c1, T c2){
        return func_(c1) -= c2;
    };

    template<class T, class = typename enable_if<std::is_arithmetic<T>::value>::type> func_ operator-(T c2, const constant_& c1){
        return (func_(c1) *= -1) += c2;
    };

    template<class T, class = typename enable_if<std::is_arithmetic<T>::value>::type> func_ operator*(const constant_& c1, T c2){
        return func_(c1) *= c2;
    };

    template<class T, class = typename enable_if<std::is_arithmetic<T>::value>::type> func_ operator*(T c2, const constant_& c1){
        return func_(c1) *= c2;
    };

    template<class T, class = typename enable_if<std::is_arithmetic<T>::value>::type> func_ operator/(const constant_& c1, T c2){
        return func_(c1) /= c2;
    };

    func_ operator*(const constant_& c1, const constant_& c2);



    func_ operator-(const constant_& c1, const constant_& c2);

    func_ operator/(const constant_& c1, const constant_& c2);



        
        
    constant_* add(constant_* c1, const constant_& c2); /**< adds c2 to c1, updates its type and returns the result **/
    //
    template<typename T> constant_* add(constant_* c1, const constant<T>& c2){ /**< adds c2 to c1, updates its type and returns the result **/
        switch (c1->get_type()) {
            case binary_c: {
                if (c2.is_binary() ) {
                    *(constant<bool>*)c1 += c2.eval();
                }
                else {
                    bool val = ((constant<bool>*)c1)->eval();
                    delete c1;
                    c1 = new constant<T>(c2.eval() + val);
                }
                return c1;
                break;
            }
            case short_c: {
                if (c2.get_type() <= short_c) {
                    *((constant<short>*)c1) += c2.eval();
                }
                else {
                    auto val = ((constant<short>*)c1)->eval();
                    delete c1;
                    c1 = new constant<T>(c2.eval() + val);
                }
                break;
            }
            case integer_c: {
                if (c2.get_type() <= integer_c) {
                    *((constant<int>*)c1) += c2.eval();
                }
                else {
                    auto val = ((constant<int>*)c1)->eval();
                    delete c1;
                    c1 = new constant<T>(c2.eval() + val);
                }
                break;
            }
            case float_c: {
                if (c2.get_type() <= float_c) {
                    *((constant<float>*)c1) += c2.eval();
                }
                else {
                    auto val = ((constant<float>*)c1)->eval();
                    delete c1;
                    c1 = new constant<T>(c2.eval() + val);
                }
                break;
            }
            case double_c: {
                if (c2.get_type() <= double_c) {
                    *((constant<double>*)c1) += c2.eval();
                }
                else {
                    auto val = ((constant<double>*)c1)->eval();
                    delete c1;
                    c1 = new constant<T>(c2.eval() + val);
                }
                break;
            }
            case long_c: {
                if (c2.get_type() <= long_c) {
                    *((constant<long double>*)c1) += c2.eval();
                }
                else {
                    auto val = ((constant<long double>*)c1)->eval();
                    delete c1;
                    c1 = new constant<T>(c2.eval() + val);
                }
                break;
            }
            case par_c:{            
                auto f = new func_(*c1);
                *f += c2;
                c1 =(constant_*)(f);
                return c1;
                break;
            }
            case var_c:{
                auto f = new func_(*c1);
                *f += c2;
                c1 =(constant_*)(f);
                return c1;
                break;
            }
    //        case uexp_c: {
    //            auto res = new bexpr(*(uexpr*)c1 + c2);
    //            delete c1;
    //            c1 = (constant_*)res;
    //            return c1;
    //            break;
    //        }
    //        case bexp_c: {
    //            auto res = new bexpr(*(bexpr*)c1 + c2);
    //            delete c1;
    //            c1 = (constant_*)res;
    //            return c1;
    //            break;
    //        }
            case func_c: {
    //            auto res = new func_((*(func_*)c1) + c2);
    //            delete c1;
    //            return c1 = (constant_*)res;
                (*(func_*)c1) += c2;
                return c1;
                break;
            }
                
            default:
                break;
        }
        return c1;
    }

    //constant_* add(constant_* c1, const func_& f);

    template<typename T> constant_* add(constant_* c1, const param<T>& c2){ /**< adds c2 to c1, updates its type and returns the result **/
        switch (c1->get_type()) {
            case binary_c: {
                auto val = ((constant<bool>*)c1)->eval();
                delete c1;
                auto f = new func_(c2);
                *f += val;
                c1 = (constant_*)(f);
                return c1;
                break;
            }
            case short_c: {
                auto val = ((constant<short>*)c1)->eval();
                delete c1;
                auto f = new func_(c2);
                *f += val;
                c1 = (constant_*)(f);
                return c1;
                break;
            }
            case integer_c: {
                auto val = ((constant<int>*)c1)->eval();
                delete c1;
                auto f = new func_(c2);
                *f += val;
                c1 = (constant_*)(f);
                return c1;
                break;
            }
            case float_c: {
                auto val = ((constant<float>*)c1)->eval();
                delete c1;
                auto f = new func_(c2);
                *f += val;
                c1 = (constant_*)(f);
                return c1;
                break;
            }
            case double_c: {
                auto val = ((constant<double>*)c1)->eval();
                delete c1;
                auto f = new func_(c2);
                *f += val;
                c1 = (constant_*)(f);
                return c1;
                break;
            }
            case long_c: {
                auto val = ((constant<long double>*)c1)->eval();
                delete c1;
                auto f = new func_(c2);
                *f += val;
                c1 = (constant_*)(f);
                return c1;
                break;
            }
            case par_c:{
                auto res = new func_(*c1);
                delete c1;
                res->insert(true, constant<double>(1), c2);
                return c1 = res;
                break;
            }
            case var_c:{
                auto res = new func_(*c1);
                delete c1;
                if (c2.is_var()) {
                    res->insert(true, constant<double>(1), c2);
                }
                else {
                    auto cst = res->get_cst();
                    cst = add(cst, c2);
                }
                return c1 = res;
                break;
            }

    //        case uexp_c: {
    //            auto res = new bexpr(*(uexpr*)c1 + c2);
    //            delete c1;
    //            c1 = (constant_*)res;
    //            return c1;
    //            break;
    //        }
    //        case bexp_c: {
    //            auto res = new bexpr(*(bexpr*)c1 + c2);
    //            delete c1;
    //            c1 = (constant_*)res;
    //            return c1;
    //            break;
    //        }
            case func_c: {
    //            auto res = new func_(*c1);
    //            delete c1;
    //            *res += c2;
    //            return c1 = (constant_*)res;
                (*(func_*)c1) += c2;
                return c1;
                break;
            }
            default:
                break;
        }
        return c1;
    }

    constant_* substract(constant_* c1, const constant_& c2);


    template<typename T> constant_* substract(constant_* c1, const param<T>& c2){ /**< Substracts c2 from c1, updates its type and returns the result **/
        switch (c1->get_type()) {
            case binary_c: {
                auto val = ((constant<bool>*)c1)->eval();
                delete c1;
                auto f = new func_(c2);
                (*f *= -1) += val;
                c1 = (constant_*)(f);
                return c1;
                break;
            }
            case short_c: {
                auto val = ((constant<short>*)c1)->eval();
                delete c1;
                auto f = new func_(c2);
                (*f *= -1) += val;
                c1 = (constant_*)(f);
                return c1;
                break;
            }
            case integer_c: {
                auto val = ((constant<int>*)c1)->eval();
                delete c1;
                auto f = new func_(c2);
                (*f *= -1) += val;
                c1 = (constant_*)(f);
                return c1;
                break;
            }
            case float_c: {
                auto val = ((constant<float>*)c1)->eval();
                delete c1;
                auto f = new func_(c2);
                (*f *= -1) += val;
                c1 = (constant_*)(f);
                return c1;
                break;
            }
            case double_c: {
                auto val = ((constant<double>*)c1)->eval();
                delete c1;
                auto f = new func_(c2);
                (*f *= -1) += val;
                c1 = (constant_*)(f);
                return c1;
                break;
            }
            case long_c: {
                auto val = ((constant<long double>*)c1)->eval();
                delete c1;
                auto f = new func_(c2);
                (*f *= -1) += val;
                c1 = (constant_*)(f);
                return c1;
                break;
            }
            case par_c:{
                auto res = new func_(*c1);
                delete c1;
                res->insert(false, constant<double>(1), c2);
                return c1 = res;
                break;
            }
            case var_c:{
                auto res = new func_(*c1);
                delete c1;
                if (c2.is_var()) {
                    res->insert(false, constant<double>(1), c2);
                }
                else {
                    auto cst = res->get_cst();
                    cst = add(cst, c2);
                }
                return c1 = res;
                break;
            }

            case uexp_c: {
                auto res = new bexpr(minus_, make_shared<func_>(func_(*c1)), make_shared<func_>(func_(c2)));
                delete c1;
                c1 = (constant_*)res;
                return c1;
                break;
            }
            case bexp_c: {
                auto res = new bexpr(minus_, make_shared<func_>(func_(*c1)), make_shared<func_>(func_(c2)));
                delete c1;
                c1 = (constant_*)res;
                return c1;
                break;
            }
            case func_c: {
    //            auto res = new func_(*c1);
    //            delete c1;
    //            *res -= c2;
    //            return c1 = res;
                (*(func_*)c1) -= c2;
                return c1;
                break;
            }
            default:
                break;
        }
        return c1;
    }



    template<typename T> constant_* multiply(constant_* c1, const param<T>& c2){ /**< adds c2 to c1, updates its type and returns the result **/
        switch (c1->get_type()) {
            case binary_c: {
                auto val = ((constant<bool>*)c1)->eval();
                delete c1;
                auto f = new func_(c2);
                *f *= val;
                c1 = (constant_*)(f);
                return c1;
                break;
            }
            case short_c: {
                auto val = ((constant<short>*)c1)->eval();
                delete c1;
                auto f = new func_(c2);
                *f *= val;
                c1 = (constant_*)(f);
                return c1;
                break;
            }
            case integer_c: {
                auto val = ((constant<int>*)c1)->eval();
                delete c1;
                auto f = new func_(c2);
                *f *= val;
                c1 = (constant_*)(f);
                return c1;
                break;
            }
            case float_c: {
                auto val = ((constant<float>*)c1)->eval();
                delete c1;
                auto f = new func_(c2);
                *f *= val;
                c1 = (constant_*)(f);
                return c1;
                break;
            }
            case double_c: {
                auto val = ((constant<double>*)c1)->eval();
                delete c1;
                auto f = new func_(c2);
                *f *= val;
                c1 = (constant_*)(f);
                return c1;
                break;
            }
            case long_c: {
                auto val = ((constant<long double>*)c1)->eval();
                delete c1;
                auto f = new func_(c2);
                *f *= val;
                c1 = (constant_*)(f);
                return c1;
                break;
            }
            case par_c:{
                auto f = new func_(*c1);
                delete c1;
                *f *= c2;
                c1 =(constant_*)(f);
                return c1;
                break;
            }
            case uexp_c: {
                auto res = new bexpr(product_, make_shared<func_>(func_(*c1)), make_shared<func_>(func_(c2)));
                delete c1;
                c1 = (constant_*)res;
                return c1;
                break;
            }
            case bexp_c: {
                auto res = new bexpr(product_, make_shared<func_>(func_(*c1)), make_shared<func_>(func_(c2)));
                delete c1;
                c1 = (constant_*)res;
                return c1;
                break;
            }
            case func_c: {
    //            auto res = new func_(*c1);
    //            delete c1;
    //            *res *= c2;
    //            return c1 = res;
                (*(func_*)c1) *= c2;
                return c1;
                break;
            }
            default:
                break;
        }
        return c1;
    }





    template<typename T> constant_* substract(constant_* c1, const constant<T>& c2){ /**< adds c2 to c1, updates its type and returns the result **/
        switch (c1->get_type()) {
            case binary_c: {
                if (c2.is_binary() ) {
                    *(constant<bool>*)c1 -= c2.eval();
                }
                else {
                    auto val = ((constant<bool>*)c1)->eval();
                    delete c1;
                    c1 = new constant<T>(val - c2.eval());
                }
                return c1;
                break;
            }
            case short_c: {
                if (c2.get_type() <= short_c) {
                    *((constant<short>*)c1) -= c2.eval();
                }
                else {
                    auto val = ((constant<short>*)c1)->eval();
                    delete c1;
                    c1 = new constant<T>(val - c2.eval());
                }
                break;
            }
            case integer_c: {
                if (c2.get_type() <= integer_c) {
                    *((constant<int>*)c1) -= c2.eval();
                }
                else {
                    auto val = ((constant<int>*)c1)->eval();
                    delete c1;
                    c1 = new constant<T>(val - c2.eval());
                }
                break;
            }
            case float_c: {
                if (c2.get_type() <= float_c) {
                    *((constant<float>*)c1) -= c2.eval();
                }
                else {
                    auto val = ((constant<float>*)c1)->eval();
                    delete c1;
                    c1 = new constant<T>(val - c2.eval());
                }
                break;
            }
            case double_c: {
                if (c2.get_type() <= double_c) {
                    *((constant<double>*)c1) -= c2.eval();
                }
                else {
                    auto val = ((constant<double>*)c1)->eval();
                    delete c1;
                    c1 = new constant<T>(val - c2.eval());
                }
                break;
            }
            case long_c: {
                if (c2.get_type() <= long_c) {
                    *((constant<long double>*)c1) -= c2.eval();
                }
                else {
                    auto val = ((constant<long double>*)c1)->eval();
                    delete c1;
                    c1 = new constant<T>(val - c2.eval());
                }
                break;
            }
            case par_c:{
                auto f = new func_(*c1);
                *f -= c2;
                c1 =(constant_*)(f);
                return c1;
                break;
            }
            case var_c:{
                auto f = new func_(*c1);
                *f -= c2;
                c1 =(constant_*)(f);
                return c1;
                break;
            }
            case uexp_c: {
                auto res = new bexpr(minus_, make_shared<func_>(func_(*c1)), make_shared<func_>(func_(c2)));
                delete c1;
                c1 = (constant_*)res;
                return c1;
                break;
            }
            case bexp_c: {
                auto res = new bexpr(minus_, make_shared<func_>(func_(*c1)), make_shared<func_>(func_(c2)));
                delete c1;
                c1 = (constant_*)res;
                return c1;
                break;
            }
            case func_c: {
    //            auto res = new func_(*(func_*)c1 - c2);
    //            delete c1;
    //            return c1 = (constant_*)res;
                (*(func_*)c1) -= c2;
                return c1;
                break;
            }        default:
                break;
        }
        return c1;
    }

    constant_* multiply(constant_* c1, const constant_& c2);
    constant_* divide(constant_* c1, const constant_& c2);


    template<typename T> constant_* multiply(constant_* c1, const constant<T>& c2){ /**< adds c2 to c1, updates its type and returns the result **/
        switch (c1->get_type()) {
            case binary_c: {
                if (c2.is_binary() ) {
                    *(constant<bool>*)c1 *= c2.eval();
                }
                else {
                    auto val = ((constant<bool>*)c1)->eval();
                    delete c1;
                    c1 = new constant<T>(c2.eval() * val);
                }
                return c1;
                break;
            }
            case short_c: {
                if (c2.get_type() <= short_c) {
                    *((constant<short>*)c1) *= c2.eval();
                }
                else {
                    auto val = ((constant<short>*)c1)->eval();
                    delete c1;
                    c1 = new constant<T>(c2.eval() * val);
                }
                break;
            }
            case integer_c: {
                if (c2.get_type() <= integer_c) {
                    *((constant<int>*)c1) *= c2.eval();
                }
                else {
                    auto val = ((constant<int>*)c1)->eval();
                    delete c1;
                    c1 = new constant<T>(c2.eval() * val);
                }
                break;
            }
            case float_c: {
                if (c2.get_type() <= float_c) {
                    *((constant<float>*)c1) *= c2.eval();
                }
                else {
                    auto val = ((constant<float>*)c1)->eval();
                    delete c1;
                    c1 = new constant<T>(c2.eval() * val);
                }
                break;
            }
            case double_c: {
                if (c2.get_type() <= double_c) {
                    *((constant<double>*)c1) *= c2.eval();
                }
                else {
                    auto val = ((constant<double>*)c1)->eval();
                    delete c1;
                    c1 = new constant<T>(c2.eval() * val);
                }
                break;
            }
            case long_c: {
                if (c2.get_type() <= long_c) {
                    *((constant<long double>*)c1) *= c2.eval();
                }
                else {
                    auto val = ((constant<long double>*)c1)->eval();
                    delete c1;
                    c1 = new constant<T>(c2.eval() * val);
                }
                break;
            }
            case par_c:{
                auto pc1 = (param_*)(c1);
                auto l = new func_(*pc1);
                *l *= c2;
                c1 =(constant_*)(l);
                return c1;
                break;
            }
            case uexp_c: {
                auto res = new bexpr(product_, make_shared<func_>(func_(*c1)), make_shared<func_>(func_(c2)));
                delete c1;
                c1 = (constant_*)res;
                return c1;
                break;
            }
            case bexp_c: {
                auto res = new bexpr(product_, make_shared<func_>(func_(*c1)), make_shared<func_>(func_(c2)));
                delete c1;
                c1 = (constant_*)res;
                return c1;
                break;
            }
            case func_c: {
                (*(func_*)c1) *= c2;
                return c1;
                break;
            }
            default:
                break;
        }
        return c1;
    }

    template<typename T> constant_* divide(constant_* c1, const constant<T>& c2){ /**< adds c2 to c1, updates its type and returns the result **/
        if (c2.eval()==0) {
            throw invalid_argument("dividing by zero!\n");
        }
        switch (c1->get_type()) {
            case binary_c: {
                if (c2.is_binary() ) {
                    *(constant<bool>*)c1 /= c2.eval();
                }
                else {
                    auto val = ((constant<bool>*)c1)->eval();
                    delete c1;
                    c1 = new constant<T>(c2.eval() / val);
                }
                return c1;
                break;
            }
            case short_c: {
                if (c2.get_type() <= short_c) {
                    *((constant<short>*)c1) /= c2.eval();
                }
                else {
                    auto val = ((constant<short>*)c1)->eval();
                    delete c1;
                    c1 = new constant<T>(c2.eval() / val);
                }
                break;
            }
            case integer_c: {
                if (c2.get_type() <= integer_c) {
                    *((constant<int>*)c1) /= c2.eval();
                }
                else {
                    auto val = ((constant<int>*)c1)->eval();
                    delete c1;
                    c1 = new constant<T>(c2.eval() / val);
                }
                break;
            }
            case float_c: {
                if (c2.get_type() <= float_c) {
                    *((constant<float>*)c1) /= c2.eval();
                }
                else {
                    auto val = ((constant<float>*)c1)->eval();
                    delete c1;
                    c1 = new constant<T>(c2.eval() / val);
                }
                break;
            }
            case double_c: {
                if (c2.get_type() <= double_c) {
                    *((constant<double>*)c1) /= c2.eval();
                }
                else {
                    auto val = ((constant<double>*)c1)->eval();
                    delete c1;
                    c1 = new constant<T>(c2.eval() / val);
                }
                break;
            }
            case long_c: {
                if (c2.get_type() <= long_c) {
                    *((constant<long double>*)c1) /= c2.eval();
                }
                else {
                    auto val = ((constant<long double>*)c1)->eval();
                    delete c1;
                    c1 = new constant<T>(c2.eval() / val);
                }
                break;
            }
            case par_c:{
                auto pc1 = (param_*)(c1);
                auto l = new func_(*pc1);
                *l /= c2;
                c1 =(constant_*)(l);
                return c1;
                break;
            }
            case uexp_c: {
                auto res = new bexpr(div_, make_shared<func_>(func_(*c1)), make_shared<func_>(func_(c2)));
                delete c1;
                c1 = (constant_*)res;
                return c1;
                break;
            }
            case bexp_c: {
                auto res = new bexpr(div_, make_shared<func_>(func_(*c1)), make_shared<func_>(func_(c2)));
                delete c1;
                c1 = (constant_*)res;
                return c1;
                break;
            }
            case func_c: {
                switch (((func_*)c1)->get_ftype()) {
                    case lin_: {
                        cerr << "Unsupported yet;\n";
                        //                    auto res = new func_(*(func_*)c1 * 1/c2);
                        //                    delete c1;
                        //                    return c1 = res;
                        break;
                    }
                    default:
                        break;
                }
            }
            default:
                break;
        }
        return c1;
    }


    func_ cos(const constant_& c);

    func_ sin(const constant_& c);


    func_ sqrt(const constant_& c);

    func_ expo(const constant_& c);

    func_ log(const constant_& c);

    func_ log(constant_&& c);


    //template<typename other_type> bexpr operator+(const other_type& c1, const expr& c2){
    //    bexpr res;
    //    res._otype = plus_;
    //    res._lson = copy((constant_*)&c1);
    //    res._rson =  copy((constant_*)&c2);
    //    res._to_str = ::to_str(res._lson) + " + " + ::to_str(res._rson);
    //    return res;
    //}
    //
    //template<typename other_type> bexpr operator+(const expr& c1, const other_type& c2){
    //    bexpr res;
    //    res._otype = plus_;
    //    res._lson = copy((constant_*)&c1);
    //    res._rson =  copy((constant_*)&c2);
    //    res._to_str = ::to_str(res._lson) + " + " + ::to_str(res._rson);
    //    return res;
    //}
    //
    //
    //template<typename other_type> bexpr operator-(const other_type& c1, const expr& c2){
    //    bexpr res;
    //    res._otype = minus_;
    //    res._lson = copy((constant_*)&c1);
    //    res._rson =  copy((constant_*)&c2);
    //    res._to_str = ::to_str(res._lson) + " - " + ::to_str(res._rson);
    //    return res;
    //}
    //
    //template<typename other_type> bexpr operator-(const expr& c1, const other_type& c2){
    //    bexpr res;
    //    res._otype = minus_;
    //    res._lson = copy((constant_*)&c1);
    //    res._rson =  copy((constant_*)&c2);
    //    res._to_str = ::to_str(res._lson) + " - " + ::to_str(res._rson);
    //    return res;
    //}
    //
    //
    //template<typename other_type> bexpr operator*(const other_type& c1, const expr& c2){
    //    bexpr res;
    //    res._otype = product_;
    //    res._lson = copy((constant_*)&c1);
    //    res._rson =  copy((constant_*)&c2);
    //    res._to_str = ::to_str(res._lson) + " * " + ::to_str(res._rson);
    //    return res;
    //}
    //
    //template<typename other_type> bexpr operator*(const expr& c1, const other_type& c2){
    //    bexpr res;
    //    res._otype = product_;
    //    res._lson = copy((constant_*)&c1);
    //    res._rson =  copy((constant_*)&c2);
    //    res._to_str = ::to_str(res._lson) + " * " + ::to_str(res._rson);
    //    return res;
    //}
    //
    //
    //
    //template<typename other_type> bexpr operator/(const other_type& c1, const expr& c2){
    //    bexpr res;
    //    res._otype = div_;
    //    res._lson = copy((constant_*)&c1);
    //    res._rson =  copy((constant_*)&c2);
    //    res._to_str = ::to_str(res._lson) + " / " + ::to_str(res._rson);
    //    return res;
    //}
    //template<typename other_type> bexpr operator/(const expr& c1, const other_type& c2){
    //    bexpr res;
    //    res._otype = div_;
    //    res._lson = copy((constant_*)&c1);
    //    res._rson =  copy((constant_*)&c2);
    //    res._to_str = ::to_str(res._lson) + " / " + ::to_str(res._rson);
    //    return res;
    //}

    template<typename type>
    func_ power(const param<type>& v, unsigned p);
    
    func_ power(const func_& f, unsigned p);

    template<typename type>
    func_ sum(const param<type>& p);
    
    template<typename type1, typename type2>
    func_ product(const param<type1>& p, const param<type2>& v);
    
    template<typename type>
    func_ product(const param<type>& p1, const func_& f);

    template<typename type>
    func_ innerproduct(const param<type>& p1, const param<type>& p2);
    
    func_ get_poly_derivative(constant_* c, const param_ &v); /*< Get the derivative of c with respect to v) */
}

#endif /* func_h */

