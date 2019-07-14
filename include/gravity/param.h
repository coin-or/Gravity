//
//  param.h
//
//
//  Created by Hassan on 13/05/2016.



#ifndef ____param__
#define ____param__
#include <sstream>
#include <iomanip>
#include <stdio.h>
#include <string>
#include <sstream>
#include <vector>
#include <list>
#include <set>
#include <typeinfo>
#include <gravity/constant.h>
#include <gravity/Arc.h>
#include <gravity/Node.h>
#include <limits>
#include <math.h>
#include <random>
#ifdef USE_QPP
    #include "qpp.h"
#endif

using namespace std;





namespace gravity {    
    
    /** Backbone class for parameter */
    class param_: public constant_ {

    protected:

        NType                                          _intype; /**< internal storage type **/

        void set_intype(NType type) {
            _intype = type;
        }

    public:

        string                                         _name = "noname";
        shared_ptr<size_t>                             _id = make_shared<size_t>(0); /**< index of current param/var */
        shared_ptr<size_t>                             _vec_id = make_shared<size_t>(0);; /**< index of the corresponding vector (for Cplex). **/
        shared_ptr<indices>                            _indices = nullptr; /*< If indexed, point to the indexing set */
        bool                                           _is_relaxed = false; /*< Is this an integer parameter/variable that has been relaxed? */
        bool                                           _new = true; /**< Will become false once this parameter/variable is added to a program. Useful for iterative model solving. */

        bool                                           _is_conjugate = false; /**< True if the parameter/variable is a complex number and is conjugated */
        bool                                           _is_sqrmag = false; /**< True if the parameter/variable is the magnitude squared of a complex number */
        bool                                           _is_angle = false; /**< True if the parameter/variable is the angle of a complex number */
        bool                                           _is_real = false; /**< True if the parameter/variable is the real part of a complex number */
        bool                                           _is_imag = false; /**< True if the parameter/variable is the imaginary part of a complex number */
        
        shared_ptr<param_>                                        _real = nullptr; /**< Pointer to the real variable in case this is a complex var */
        shared_ptr<param_>                                        _imag = nullptr; /**< Pointer to the imaginary variable in case this is a complex var */

        shared_ptr<param_>                                        _mag = nullptr; /**< Pointer to the magnitude variable in case this is a complex var */
        shared_ptr<param_>                                        _ang = nullptr; /**< Pointer to the angle variable in case this is a complex var */

        
        /* For Ipopt Use */
        vector<double>                                 _l_dual; /*<<Dual values for lower bounds */
        vector<double>                                 _u_dual; /*<<Dual values for upper bounds */
        
        /**
         A shallow copy of p (ignoring _val and _range)
         @param[in] p param_ to copy from.
         */
        void shallow_copy(const param_& p) {
            _id = p._id;
            _vec_id = p._vec_id;
            _name = p._name;
            _is_transposed = p._is_transposed;
            _is_vector = p._is_vector;
            _is_angle = p._is_angle;
            _is_sqrmag = p._is_sqrmag;
            _is_conjugate = p._is_conjugate;
            _is_real = p._is_real;
            _is_imag = p._is_imag;
            _real = p._real;
            _imag = p._imag; _mag = p._mag; _ang = p._ang;
            _indices = p._indices;
            _dim[0] = p._dim[0];
            _dim[1] = p._dim[1];
        }

        virtual void reset_range(){};
        
        /** let p share the values and indices of current var */
        virtual void share_vals(const shared_ptr<param_>& p){};
        virtual void initialize_uniform(){};
        virtual void initialize_zero(){};
        virtual shared_ptr<param_> pcopy() const{return nullptr;};
        
        virtual void print(){};
        
        virtual void print_vals(int prec){};
        
        virtual void print_symbolic() const{};
        
        virtual void print(int prec){};                

        void set_id(size_t idx) {
            *_id = idx;
        };

        void set_vec_id(size_t idx) {
            *_vec_id = idx;
        };

        size_t get_id() const {
            return *_id;
        };

        size_t get_vec_id() const {
            return *_vec_id;
        };

        indices get_indices() const{return *_indices;};
        
        inline size_t get_id_inst(size_t inst = 0) const {
            if (is_indexed()) {
                if(_indices->_ids->at(0).size() <= inst){
                    throw invalid_argument("param::get_id_inst(size_t inst) inst is out of range");
                }
                return _indices->_ids->at(0).at(inst);
            }
//            auto dim = get_dim();
//            if(inst > dim-1){
//                throw invalid_argument("param::get_id_inst(size_t inst) inst is out of range");
//            }
            if(_dim[0]==1 && _dim[1]==1){/* unidimensional param */
                return 0;
            }
            return inst;
        };

        size_t get_id_inst(size_t inst1, size_t inst2) const {
            if (is_double_indexed()) {
                if (_indices->_ids->size()<=inst1) {
                    throw invalid_argument("get_id_inst(size_t inst1, size_t inst2) inst1 out of range\n");
                }
                if (_indices->_ids->at(inst1).size()<=inst2) {
                    throw invalid_argument("get_id_inst(size_t inst1, size_t inst2) inst2 out of range\n");
                }
                return _indices->_ids->at(inst1).at(inst2);
            }
            return get_id_inst(inst2);
//            throw invalid_argument("Calling get_id_inst(size_t inst1, size_t inst2) on a non-indexed param\n");
        };


        string get_name(bool in_func, bool exclude_indexing) const{
//            return _name;
            string name = _name;
            if(_indices && exclude_indexing){
                name = name.substr(0, name.find_first_of("."));
            }
            if (!in_func && _is_transposed) {
                name += "\u1D40";
            }
            return name;
        };

        void set_name(const string s) {
            _name = s;
        };

        string get_name(size_t inst) const {/*< Get the name of the indexed version of this variable */
            string name = _name;
            name = name.substr(0, name.find_first_of("."));
            if(_is_imag || _is_real || _is_angle || _is_conjugate){
                if (name.find(")")==std::string::npos) {
                    name += ")";
                }
            }
            else if(_is_sqrmag){
                if (name.find("|²")==std::string::npos) {
                    name += "|²";
                }
            }
            if (_is_vector) {//Remove brackets
                if(name.back()==']'){
                    name = name.substr(1, name.size()-2);
                }
                else {
                    name = name.substr(1, name.size()-1);
                }
            }
            if (is_indexed() && name.find("[")!=std::string::npos) {// Name  already has index
                return name;
            }
            if (is_indexed()) {
                size_t rev_idx = _indices->_ids->at(0).at(inst);
                name += "["+_indices->_keys->at(rev_idx)+"]";
            }
            else if(_indices){
                name += "["+_indices->_keys->at(get_id_inst(inst))+"]";
            }
            else {
                name += "["+to_string(inst)+"]";
            }
            return name;
        };

        string get_name(size_t inst1, size_t inst2) const {
            string name = _name;
            name = name.substr(0, name.find_first_of("."));
            if(_is_imag || _is_real || _is_angle || _is_conjugate){
                if (name.find(")")==std::string::npos) {
                    name += ")";
                }
            }
            else if(_is_sqrmag){
                if (name.find("|²")==std::string::npos) {
                    name += "|²";
                }
            }
            if (_is_vector) {//Remove brackets
                if(name.back()==']'){
                    name = name.substr(1, name.size()-2);
                }
                else if(name.front()=='['){
                    name = name.substr(1, name.size()-1);
                }
            }
            if (is_indexed() && name.find("[")!=std::string::npos) {// Name has index already
                return name;
            }
            if (is_indexed()) {
                name += "["+_indices->_keys->at(_indices->_ids->at(inst1).at(inst2))+"]";
            }
            else if(_indices){
                name += "["+_indices->_keys->at(inst1)+","+_indices->_keys->at(inst2)+"]";
            }
            else{
                name += "["+to_string(inst1)+","+to_string(inst2)+"]";
            }
            return name;
        };

        NType get_intype() const {
            return _intype;
        }

        shared_ptr<map<string,size_t>> get_keys_map() const {
            if(_indices){
                return _indices->_keys_map;
            }
            throw invalid_argument("Calling get_keys_map() on a non-indexed param/var.");
        }

        shared_ptr<vector<string>> get_keys() const {
            if(_indices){
                return _indices->_keys;
            }
            throw invalid_argument("Calling get_keys() on a non-indexed param/var.");
        }


        shared_ptr<vector<vector<size_t>>> get_ids() const {
            if(_indices){
                return _indices->_ids;
            }
            throw invalid_argument("Calling get_ids() on a non-indexed param/var.");
        }


        /* Querries */

        inline bool is_indexed() const{
            return (_indices && _indices->_ids);
        }
        
        bool is_double_indexed() const{
            return (_indices && _indices->_ids && _indices->_ids->size()>1);
        }

        bool is_binary() const {
            return (_intype==binary_);
        };

        bool is_integer() const {
            return (_intype==integer_ || _intype==short_);
        };

        bool is_float() const {
            return (_intype==float_);
        };

        bool is_double() const {
            return (_intype==double_);
        };

        bool is_long() const {
            return (_intype==long_);
        };

        bool is_complex() const {
            return (_intype==complex_);
        };

//        Sign get_all_sign() const{return unknown_;}; /**< If all instances of the current parameter/variable have the same sign, it returns it, otherwise, it returns unknown. **/
//        virtual Sign get_sign(size_t idx = 0) const; /**< returns the sign of one instance of the current parameter/variable. **/
        

        /** Operators */
        
        
        
        /**
         Index the current object using incoming edges for nodes stored in vec. This is a double indexing where each row corresponds to a node, and columns correspond to the edge ids.
         @param[in] vec vector of nodes
         */
        void index_in_arcs(const vector<Node*>& vec) {
            _indices->_ids = make_shared<vector<vector<size_t>>>();
            if(vec.empty()){
                DebugOff("In function param.in_arcs(const vector<Node*>& vec), vec is empty!\n. Creating and empty variable! Check your sum/product operators.\n");
                _name += "_EMPTY";
                return;
            }
            string key;
            size_t inst = 0;
            for(auto it = vec.begin(); it!= vec.end(); it++) {
                if(!(*it)->_active) {
                    continue;
                }
                _indices->_ids->push_back(vector<size_t>());
                for (auto &a:(*it)->get_in()) {
                    if (!a->_active) {
                        continue;
                    }
                    key = a->_name;
                    auto it1 = _indices->_keys_map->find(key);
                    if (it1 == _indices->_keys_map->end()){
                        throw invalid_argument("In function param.index_in_arcs(const vector<Node*>& vec), unknown arc key.");
                    }
                    _indices->_ids->at(inst).push_back(it1->second);
                }
                ++inst;
            }
        }
        
        /**
         Index the current object using outgoing edges for nodes stored in vec. This is a double indexing where each row corresponds to a node, and columns correspond to the edge ids.
         @param[in] vec vector of nodes
         */
        void index_out_arcs(const vector<Node*>& vec) {
            _indices->_ids = make_shared<vector<vector<size_t>>>();
            if(vec.empty()){
                DebugOff("In function param.index_out_arcs(const vector<Node*>& vec), vec is empty!\n. Creating and empty variable! Check your sum/product operators.\n");
                _name += "_EMPTY";
                return;
            }
            string key;
            size_t inst = 0;
            for(auto it = vec.begin(); it!= vec.end(); it++) {
                if(!(*it)->_active) {
                    continue;
                }
                _indices->_ids->push_back(vector<size_t>());
                for (auto &a:(*it)->get_out()) {
                    if (!a->_active) {
                        continue;
                    }
                    key = a->_name;
                    auto it1 = _indices->_keys_map->find(key);
                    if (it1 == _indices->_keys_map->end()){
                        throw invalid_argument("In function param.in_arcs(const vector<Node*>& vec), unknown arc key.");
                    }
                    _indices->_ids->at(inst).push_back(it1->second);
                }
                ++inst;
            }
        }
        
        /** Index param/var based on auxiliary objects attached to nodes in vec
         @param[in] vec vector of nodes
         */
        void index_in_aux(const vector<Node*>& vec, const string& aux_type) {
            _indices->_ids = make_shared<vector<vector<size_t>>>();
            if(vec.empty()){
                DebugOff("In function param.index_in_aux(const vector<Node*>& vec), vec is empty!\n. Creating and empty variable! Check your sum/product operators.\n");
                _name += "_EMPTY";
                return;
            }
            string key;
            size_t inst = 0;
            for(auto it = vec.begin(); it!= vec.end(); it++) {
                if(!(*it)->_active) {
                    continue;
                }
                _indices->_ids->push_back(vector<size_t>());
                for (auto &a:(*it)->get_aux(aux_type)) {
                    if (!a->_active) {
                        continue;
                    }
                    key = a->_name;
                    auto it1 = _indices->_keys_map->find(key);
                    if (it1 == _indices->_keys_map->end()){
                        throw invalid_argument("In function param.in_aux(const vector<Node*>& vec), unknown arc key.");
                    }
                    _indices->_ids->at(inst).push_back(it1->second);
                }
                ++inst;
            }
        }
        
        bool operator==(const param_& p) const {
            return (_id==p._id && _type==p._type && _intype==p._intype && get_name(false,false)==p.get_name(false,false));
        }

        size_t get_dim() const{
            return constant_::get_dim();
        }

        size_t get_dim(size_t i) const{
            if(is_double_indexed()){
                if(i>_indices->_ids->size()){
                    throw invalid_argument("get_dim(size_t i) i out of range\n");
                }
                return _indices->_ids->at(i).size();
            }
            if(is_indexed()){
                return _indices->_ids->at(0).size();
            }
            return this->_dim[0];
        }
        
        size_t get_nb_inst() const{
            if(is_double_indexed())
                return _indices->_ids->size();
            if(is_indexed() && !_is_transposed){
                return _indices->_ids->at(0).size();
            }
            return this->_dim[0];
        }
        
        /** set the _evaluated flag to false for the variable's bounds */
        virtual void reset_bounds(){};
        
        /** Fill x with the variable's values */
        virtual void set_double_val(double* x){};
        /** Fill the variable's values from x */
        virtual void get_double_val(const double* x){};
        
        /** Fill x with the variable's values */
        virtual void get_var(const vector<double>& x){};
        /** Fill the variable's values from x */
        virtual void set_var(const vector<double>& x){};
        
        /** Fill the variable's value at pos to x */
        virtual void get_double_val(size_t pos, double x){};
        
        /** Fill x from the variable's value at pos */
        virtual void set_double_val(size_t pos, double& x){};
        
        /** round the value stored at position i to the nearest integer */
        virtual void round_vals(){};

        /** Fill x with the variable's lower bound values */
        virtual void set_double_lb(double* x){};
        /** Fill x with the variable's upper bound values */
        virtual void set_double_ub(double* x){};
        
        /** Return lower bound violation */
        virtual double get_lb_violation(size_t i){return 0;};
        /** Return upper bound violation */
        virtual double get_ub_violation(size_t i){return 0;};
        
        virtual void copy_vals(const shared_ptr<param_>& p){};
        virtual void copy_bounds(const shared_ptr<param_>& p){};
        virtual double get_double_lb(size_t i) const{return 0;};
        virtual double get_double_ub(size_t i) const{return 0;};
        
        virtual int get_num_partns() const{return 0;};
        virtual int get_cur_partn() const{return 0;};
        
        virtual bool get_lift() const{return 0;};

    };


    /** A parameter can be a bool, a short, an int, a float, a double, a long double or a complex<double>. */
    template<typename type = double>
    class param: public param_ {

    public:

        shared_ptr<vector<type>>                _val = nullptr; /**< vector of values **/
        shared_ptr<pair<type,type>>             _range = nullptr; /**< (Min,Max) values in vals **/

        param() {
            update_type();
            init_range();
            _val = make_shared<vector<type>>();
        }
        
        shared_ptr<param_> pcopy() const{return make_shared<param>(*this);};

        shared_ptr<constant_> copy()const{return make_shared<param>(*this);};
        
        ~param(){};
        template<class T2, typename std::enable_if<is_convertible<T2, type>::value && sizeof(T2) < sizeof(type)>::type* = nullptr>
        param (const param<T2>& p) {
            *this = p;            
        }

        param (const param& p):param() {
            *this = p;
        }
        
        param (param&& p):param() {
            *this = move(p);
        }
        
        /** let this share the values, indices and range of p */
        void share_vals_ids(param<type>& p){
            if(p._indices){
                this->_indices = p._indices;
            }
            this->_dim[0] = p._dim[0];
            this->_dim[1] = p._dim[1];
            this->_val = p._val;
            this->_range = p._range;
        }
        
        /** let this share the values of p */
        template<class T2, typename std::enable_if<!is_same<T2, type>::value>::type* = nullptr>
        void share_vals_(param<T2>& p){
            throw invalid_argument("cannot share vals with different typed params/vars");
        }
        
        /** let this share the values of p */
        template<class T2, typename std::enable_if<is_same<T2, type>::value>::type* = nullptr>
        void share_vals_(param<T2>& pp){
            this->_val = pp._val;
        }
        
        /** let this share the values of p */
        void share_vals(const shared_ptr<param_>& p){
            switch (p->get_intype()) {
                case binary_:{
                    auto pp =  static_pointer_cast<param<bool>>(p);
                    share_vals_(*pp);
                }
                    break;
                case short_:{
                    auto pp =  static_pointer_cast<param<short>>(p);
                    share_vals_(*pp);
                }
                    break;
                case integer_:{
                    auto pp =  static_pointer_cast<param<int>>(p);
                    share_vals_(*pp);
                }
                    break;
                case float_:{
                    auto pp =  static_pointer_cast<param<float>>(p);
                    share_vals_(*pp);
                }
                    break;
                case double_:{
                    auto pp =  (param<double>*)(p.get());
                    share_vals_(*pp);
                }
                    break;
                case long_:{
                    auto pp =  static_pointer_cast<param<long double>>(p);
                    share_vals_(*pp);
                }
                    break;
                case complex_:{
                    auto pp =  static_pointer_cast<param<Cpx>>(p);
                    share_vals_(*pp);
                }
                    break;
                default:
                    break;
            }
        }
        
        param deep_copy() const{
            param res;
            res._type = _type;
            res._polar = _polar;
            res._intype = _intype;
//            res._id = make_shared<size_t>();
//            res._vec_id = make_shared<size_t>();
//            res._val = make_shared<vector<type>>(*_val);
            res._val = make_shared<vector<type>>();
            res._val->resize(_val->size());
//            res._val = make_shared<vector<type>>(*_val);
            res._range = make_shared<pair<type,type>>(*_range);
            res._name = _name;
            res._is_transposed = _is_transposed;
            res._is_vector = _is_vector;
            res._new = _new;
            res._is_relaxed = _is_relaxed;
            res._is_angle = _is_angle;
            res._is_sqrmag = _is_sqrmag;
            res._is_conjugate = _is_conjugate;
            res._is_real = _is_real;
            res._is_imag = _is_imag;
            res._real = _real;
            res._imag = _imag; res._mag = res._mag; res._ang = _ang;
            if(_indices){
                res._indices = make_shared<indices>();
                res._indices->shallow_copy(_indices);
            }
            res._dim[0] = _dim[0];
            res._dim[1] = _dim[1];
            return res;
        }

        param& operator=(const param& p) {
            _type = p._type;
            _polar = p._polar;
            _intype = p._intype;
            _id = p._id;
            _vec_id = p._vec_id;
            _val = p._val;
            _range = p._range;
            _name = p._name;
            _is_transposed = p._is_transposed;
            _is_vector = p._is_vector;
            _new = p._new;
            _is_relaxed = p._is_relaxed;
            _is_angle = p._is_angle;
            _is_sqrmag = p._is_sqrmag;
            _is_conjugate = p._is_conjugate;
            _is_real = p._is_real;
            _is_imag = p._is_imag;
            if(p._real)
                _real = p._real->pcopy();
            if(p._imag)
                _imag = p._imag->pcopy();
            if(p._mag)
                _mag = p._mag->pcopy();
            if(p._ang)
                _ang = p._ang->pcopy();
            if(p._indices){
                _indices = make_shared<indices>();
                _indices->shallow_copy(p._indices);
            }
            _dim[0] = p._dim[0];
            _dim[1] = p._dim[1];
            return *this;
        }
        
        template<class T2, typename std::enable_if<is_convertible<T2, type>::value && sizeof(T2) < sizeof(type)>::type* = nullptr>
        param& operator=(const param<T2>& p) {
            update_type();
            _id = p._id;
            _polar = p._polar;
            _vec_id = p._vec_id;
            _val = make_shared<vector<type>>();
            _val->resize(p._val->size());
            for(auto i = 0; i<p._val->size();i++){
                _val->at(i) = p._val->at(i);
            }
            _range = make_shared<pair<type,type>>();
            _range->first = p._range->first;
            _range->second = p._range->second;
            _name = p._name;
            _is_transposed = p._is_transposed;
            _is_vector = p._is_vector;
            _new = p._new;
            _is_relaxed = p._is_relaxed;
            _is_angle = p._is_angle;
            _is_sqrmag = p._is_sqrmag;
            _is_conjugate = p._is_conjugate;
            _is_real = p._is_real;
            _is_imag = p._is_imag;
            if(p._real)
                _real = p._real->pcopy();
            if(p._imag)
                _imag = p._imag->pcopy();
            if(p._mag)
                _mag = p._mag->pcopy();
            if(p._ang)
                _ang = p._ang->pcopy();
            if(p._indices){
                _indices = make_shared<indices>(*p._indices);
            }
            _dim[0] = p._dim[0];
            _dim[1] = p._dim[1];
            return *this;
        }

        param& operator=(param&& p) {
            _type = p._type;
            _polar = p._polar;
            _intype = p._intype;
            _id = p._id;
            _vec_id = p._vec_id;
            _val = move(p._val);
            _range = move(p._range);
            _name = p._name;
            _is_transposed = p._is_transposed;
            _is_vector = p._is_vector;
            _new = p._new;
            _is_relaxed = p._is_relaxed;
            _is_angle = p._is_angle;
            _is_sqrmag = p._is_sqrmag;
            _is_conjugate = p._is_conjugate;
            _is_real = p._is_real;
            _is_imag = p._is_imag;
            _real = p._real;
            _imag = p._imag; _mag = p._mag; _ang = p._ang;
            _indices = move(p._indices);
            _dim[0] = p._dim[0];
            _dim[1] = p._dim[1];
            return *this;
        }

        param tr() const {
            auto p = param(*this);
            p.constant_::transpose();
            p._name = "["+p._name+"]";
            return p;
        }

        param vec() const {
            auto p = param(*this);
            p._is_vector = true;
            p._name = "["+p._name+"]";
            return p;
        }


        shared_ptr<vector<type>> get_vals() const {
            return _val;
        }


        void update_type() {
            _type = par_c;
            if(typeid(type)==typeid(bool)) {
                _intype = binary_;
                return;
            }
            if(typeid(type)==typeid(short)) {
                _intype = short_;
                return;
            }
            if(typeid(type)==typeid(int)) {
                _intype = integer_;
                return;
            }
            if(typeid(type)==typeid(float)) {
                _intype = float_;
                return;
            }
            if(typeid(type)==typeid(double)) {
                _intype = double_;
                return;
            }
            if(typeid(type)==typeid(long double)) {
                _intype = long_;
                return;
            }
            if(typeid(type)==typeid(Cpx)) {
                _intype = complex_;
                return;
            }
            throw invalid_argument("Unsupported numerical parameter type");
        }


        param(const string& s): param(){
            _name = s;
        }

        NType get_intype() const {
            return _intype;
        }

        inline type eval() const {
            if (is_indexed()) {
                return _val->at(_indices->_ids->at(0).back());
            }
            return _val->back();
        }

        inline type eval(size_t i) const {
//            if(is_matrix()){
//                throw invalid_argument("eval() should be called with double index here\n");
//            }
            auto idx = get_id_inst(i);
//            if (is_indexed()) {
//                if (_indices->_ids->size()>1) {
//                    throw invalid_argument("eval() should be called with double index here\n");
//                }
//                if (_val->size()<=idx){
//                    throw invalid_argument("Param eval out of range");
//                }
//                return _val->at(idx);
//            }
//            if (_val->size()<=idx){
//                throw invalid_argument("Param eval out of range");
//            }
            return _val->at(idx);
        }


        type eval(const string& key) const{
            return _val->at(param_::_indices->_keys_map->at(key));
        }

        inline type eval(size_t i, size_t j) const {

            if (is_indexed() && _indices->_ids->size()>1) {
//                if (_indices->_ids->at(i).at(j) >= _val->size()) {
//                    throw invalid_argument("eval(i,j): out of range");
//                }
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


        /* Modifiers */
        
        void mag_ang(const param<>& pmag, const param<>& pang){
            this->_mag = make_shared<param<>>(pmag);
            this->_ang = make_shared<param<>>(pang);
            this->_polar = true;
        }
        
        void real_imag(const param<>& pr, const param<>& pi){
            this->_real = make_shared<param<>>(pr);
            this->_imag = make_shared<param<>>(pi);
        }
        
        void set_real(const param<>& p){
            _real = make_shared<param<>>(p);
        }
        
        void set_imag(const param<>& p){
            _imag = make_shared<param<>>(p);
        }
        

        
        
        template<typename T=type,
        typename std::enable_if<is_arithmetic<T>::value>::type* = nullptr>
        void init_range() {
            _range = make_shared<pair<type,type>>(make_pair<>(numeric_limits<type>::max(), numeric_limits<type>::lowest()));
        }
        
        
        template<class T=type, class = typename enable_if<is_same<T, Cpx>::value>::type>
        void init_range() {
            _range = make_shared<pair<type,type>>(make_pair<>(Cpx(numeric_limits<double>::max(), numeric_limits<double>::max()), Cpx(numeric_limits<double>::lowest(), numeric_limits<double>::lowest())));
        }

        
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
                throw invalid_argument("Cannot call param::add_val(type val) on matrix");
            }
            _val->push_back(val);
            update_range(val);
            _dim[0] = _val->size();
        }
        
        void set_range(type v){
            _range->first = v;
            _range->second = v;
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
                throw invalid_argument("Cannot call param::add_val(type val) on matrix");
            }
            _dim[0] = max(_dim[0],i+1);
            _val->resize(max(_val->size(),i+1));
            _val->at(i) = val;
            update_range(val);
        }

        void set_val(size_t i, size_t j, type val) {
//            if(!is_matrix()){
//                throw invalid_argument("Function set_val(size_t i, size_t j, type val) should be called on a matrix");
//            }
//            if(_dim[0] <= i || _dim[1] <= j){
//                throw invalid_argument("In Function set_val(size_t i, size_t j, type val), i or j are out of bounds");
//            }
            update_range(val);
            if (_is_transposed) {
//                if(_val->at(_dim[0]*j+i)==_range->first ||  _val->at(_dim[0]*j+i)==_range->second || val<_range->first || val>_range->second){
                    _val->at(_dim[0]*j+i) = val;
//                    reset_range();
//                }
//                else {
//                    _val->at(_dim[0]*j+i) = val;
//                }
            }
           else {
//                if(_val->at(_dim[1]*i+j)==_range->first ||  _val->at(_dim[1]*i+j)==_range->second || val<_range->first || val>_range->second){
//                    _val->at(_dim[1]*i+j) = val;
//                    reset_range();
//                }
//                else {
                    _val->at(_dim[1]*i+j) = val;
//                }
           }
        }



        size_t set_val(const string& key, type val) {
            auto it = _indices->_keys_map->find(key);
            if (it == _indices->_keys_map->end()){
                throw invalid_argument("in Function size_t set_val(const string& key, type val), unknown key"+key);
            }
            if(_val->at(it->second)==_range->first ||  _val->at(it->second)==_range->second || val<_range->first || val>_range->second){
                _val->at(it->second) = val;
                reset_range();
            }
            else {
                _val->at(it->second) = val;
            }
            return it->second;
        }

        size_t add_val(const string& key, type val) {
            if(!_indices){
                _indices = make_shared<indices>();
            }
            auto index = param_::_indices->size();
            auto pp = param_::_indices->_keys_map->insert(make_pair<>(key,index));
            if (pp.second) {//new index inserted
                _val->resize(std::max(_val->size(),index+1));
                _dim[0] = std::max(_dim[0],_val->size());
                _indices->_keys->resize(_val->size());
                _indices->_keys->at(index) = key;
                _val->at(index) = val;
                update_range(val);
                return index;
            }
            else {
                Warning("WARNING: calling add_val(const string& key, T val) with an existing key, overriding existing value" << endl);
                set_val(key,val);
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
                if(_val->at(_indices->_ids->at(0).at(i))==_range->first ||  _val->at(_indices->_ids->at(0).at(i))==_range->second || val<_range->first || val>_range->second){
                    reset_range();
                }
            }
            if (_val->size()<=i){
                throw invalid_argument("Param set_val(size_t i, type val) out of range");
            }
            if(_val->at(i)==_range->first ||  _val->at(i)==_range->second || val<_range->first || val>_range->second){
                _val->at(i) = val;
                reset_range();
            }
            else{
                _val->at(i) = val;
            }
        }

        void set_val(type val) {
            if(is_indexed()){
                update_range(val);
                bool reset = false;
                for(auto &idx: _indices->_ids->at(0)){
                    if(_val->at(idx)==_range->first ||  _val->at(idx)==_range->second){
                        reset = true;
                    }
                    _val->at(idx) = val;
                }
                if(reset){
                    reset_range();
                }
            }
            else {
                for (auto i = 0; i<_val->size() ;i++) {
                    _val->at(i) = val;
                }
                _range->first = val;
                _range->second = val;
            }
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
        Sign get_all_sign() const{
            return get_all_sign_();
        }
        
        template<class T=type, class = typename enable_if<is_same<T, Cpx>::value>::type> Sign get_all_sign_() const{
            if (_range->first == Cpx(0,0) && _range->second == Cpx(0,0)) {
                return zero_;
            }
            if ((_range->second.real() < 0 && _range->second.imag() < 0)) {
                return neg_;
            }
            if ((_range->second.real() > 0 && _range->second.imag() > 0)) {
                return pos_;
            }
            if (_range->second <= Cpx(0,0)) {
                return non_pos_;
            }
            if (_range->first >= Cpx(0,0)) {
                return non_neg_;
            }
            return unknown_;
        }

        template<typename T=type,
        typename std::enable_if<is_arithmetic<T>::value>::type* = nullptr> Sign get_all_sign_() const {
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


        bool is_unit() const{
            return is_unit_();
        }
        
        template<class T=type, typename enable_if<is_arithmetic<T>::value>::type* = nullptr> bool is_unit_() const { /**< Returns true if all values of this paramter are 1 **/
            return (_range->first == 1 && _range->second == 1);
        }
        
        template<class T=type, class = typename enable_if<is_same<T, Cpx>::value>::type> bool is_unit_() const{
            return (_range->first == Cpx(1,0) && _range->second == Cpx(1,0));
        }


        bool is_zero() const{
            return is_zero_();
        }
        
        template<class T=type, typename enable_if<is_arithmetic<T>::value>::type* = nullptr> bool is_zero_() const { /**< Returns true if all values of this paramter are 0 **/
//            return (get_dim()==0 || (_range->first == 0 && _range->second == 0));
            return (get_dim()==0);
        }

        template<class T=type, class = typename enable_if<is_same<T, Cpx>::value>::type> bool is_zero_() const{
            return (get_dim()==0 || (_range->first == Cpx(0,0) && _range->second == Cpx(0,0)));
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

        /** Operators */
        bool operator==(const param& p) const {
            if (get_name(false,false)!=p.get_name(false,false) || _type!=p._type || _intype!=p._intype || _dim[0]!=p._dim[0] || _dim[1]!=p._dim[1]) return false;
            if(_indices==p._indices) return true; /* accounts for both being nullptr */
            if((_indices && !p._indices) || (p._indices && !_indices) || (*_indices != *p._indices)) return false;
            return true;
        }
        
        bool operator!=(const param& p) const {
            return !(*this==p);
        }

        template<class T=type, class = typename enable_if<is_arithmetic<T>::value>::type> void initialize_normal(double mean, double dev) {
            if (!is_double()) {
                throw invalid_argument("Function void initialize_normal(double mean, double dev) is only implemented for double typed params/vars");
            }
            std::default_random_engine generator;
            std::normal_distribution<type> distribution(mean,dev);
            for (size_t i = 0; i<_val->size(); i++) {
                _val->at(i) = distribution(generator);
            }
        }




#ifdef USE_QPP
        /* Matrix representation of a Quantum T gate */
        void QuantumT(size_t qubit_pos, size_t nb_qubits, bool transpose=false) {
            using namespace qpp;
            if (transpose) {
                auto U = gt.expandout(adjoint(gt.T), qubit_pos, nb_qubits);
                Debug("T transpose matrix at position " << to_string(qubit_pos) <<" = " << endl);
                Debug(disp(U) << "\n");
                set_vals(U.sparseView());
            }
            else {
                auto U = gt.expandout(gt.T, qubit_pos, nb_qubits);
                Debug("T matrix at position " << to_string(qubit_pos) <<" = " << endl);
                Debug(disp(U) << "\n");
                set_vals(U.sparseView());
            }
        }

        /* Matrix representation of a Quantum Hadamard gate */
        void QuantumH(size_t qubit_pos, size_t nb_qubits) {
            using namespace qpp;
            auto U = gt.expandout(gt.H, qubit_pos, nb_qubits);
            Debug("H matrix at position " << to_string(qubit_pos) <<" = " << endl);
            Debug(disp(U) << "\n");
            set_vals(U.sparseView());
        }

        /* Matrix representation of a Quantum Rx gate */
        void QuantumRx(size_t qubit_pos, size_t nb_qubits) {
            using namespace qpp;
            auto U = gt.expandout(gt.Rn(pi/2, {1,0,0}), qubit_pos, nb_qubits);
            Debug("Rx matrix at position " << to_string(qubit_pos) <<" = " << endl);
            Debug(disp(U) << "\n");
            set_vals(U.sparseView());
        }

        /* Matrix representation of a Quantum S gate */
        void QuantumS(size_t qubit_pos, size_t nb_qubits) {
            using namespace qpp;
            auto U = gt.expandout(gt.S, qubit_pos, nb_qubits);
            Debug("S matrix at position " << to_string(qubit_pos) <<" = " << endl);
            Debug(disp(U) << "\n");
            set_vals(U.sparseView());
            Debug(to_str(true));
        }


        /* Matrix representation of a Quantum S conjugate gate */
        void QuantumSt(size_t qubit_pos, size_t nb_qubits) {
            using namespace qpp;
            auto U = gt.expandout(adjoint(gt.S), qubit_pos, nb_qubits);
            Debug("S* matrix at position " << to_string(qubit_pos) <<" = " << endl);
            Debug(disp(U) << "\n");
            set_vals(U.sparseView());
            Debug(to_str(true));
        }

        /* Matrix representation of a Quantum Cnot gate with qc as control qubit and qt as target one */
        void QuantumCnot(size_t qc, size_t qt, size_t nb_qubits) {
            using namespace qpp;
            auto U = gt.expandout(gt.CTRL(gt.X, {qc}, {qt}, nb_qubits), 0, 1, pow(2,nb_qubits));
            Debug("Cnot matrix from " << to_string(qc) << " to " << to_string(qt) << " = " << endl);
            Debug(disp(U) << "\n");
            set_vals(U.sparseView());
            Debug(to_str(true));
        }
        /* Matrix representation of a Quantum Swap gate with qc and qt as the selected qubits */
        void QuantumSwap(size_t qc, size_t qt, size_t nb_qubits) {
            using namespace qpp;
            auto M1 = gt.expandout(gt.H, qc, nb_qubits);
            auto M2 = gt.expandout(gt.H, qt, nb_qubits);
            auto M3 = gt.expandout(gt.CTRL(gt.X, {qc}, {qt}, nb_qubits), 0, 1, pow(2,nb_qubits));
            auto U = M3*M1*M2*M3*M1*M2*M3;
            Debug("Swap matrix from " << to_string(qc) << " to " << to_string(qt) << " = " << endl);
            Debug(disp(U) << "\n");
            set_vals(U.sparseView());
            Debug(to_str(true));
        }
#endif

        param& operator=(type v) {
            if(_indices){
                set_val(v);
            }
            else {
                add_val(v);
            }
            return *this;
        }

        param& operator=(const initializer_list<type>& l) {
            if(_indices){
                for(auto &v:l)
                    set_val(v);
            }
            else {
                for(auto &v:l)
                    add_val(v);
            }
            return *this;
        }

        param operator()(size_t i, size_t j){
            if(!is_matrix()){
                throw invalid_argument("Current param/var is not a matrix, cannot call: param(i,j).");
            }
            param res(*this);
            if(!res._indices){
                res._indices = make_shared<indices>();
            }
            res._indices->_ids = make_shared<vector<vector<size_t>>>();
            res._indices->_ids->resize(1);
            if (_is_transposed) {
                res._indices->_ids->at(0).push_back(j*_dim[0]+i);
            }
            else{
                res._indices->_ids->at(0) .push_back(i*_dim[1]+j);
            }
            res._name += "["+to_string(i)+","+to_string(j)+"]";
            return res;
        }

        param operator()(size_t idx) {
            if(!_indices){
                throw invalid_argument("Current param/var is not indexed.");
            }
            return (*this)(this->_indices->_keys->at(idx));
        }

        template<bool...> struct bool_pack;
        template<bool... bs>
        using all_true = std::is_same<bool_pack<bs..., true>, bool_pack<true, bs...>>;
        template<class R, class... Ts>
        using are_all_convertible = all_true<std::is_convertible<Ts, R>::value...>;

        template<typename... Args, typename = typename enable_if<are_all_convertible<string, Args...>::value>::type>
        param operator()(string key1, Args&&... args) {
            if(!_indices){
                throw invalid_argument("Current param/var is not indexed.");
            }
            param res(*this);
            auto key = index_(key1, args...);
            auto it1 = _indices->_keys_map->find(key._name);
            if (it1 == _indices->_keys_map->end()){
                throw invalid_argument("In operator()(string key1, Args&&... args), unknown key");
            }
            res._name += "["+key._name+"]";
            res._indices->_ids = make_shared<vector<vector<size_t>>>();
            res._indices->_ids->resize(1);
            res._indices->_ids->at(0).push_back(it1->second);
            res._dim[0] = 1;
            return res;
        }

        /* Use this for dense indexing, prefer indices() for sparse indexing */
        param& in(const space& s){
            set_size(s._dim);
            if(s._dim.size()==1){ /* We can afford to build indices since this is a 1-d set */
                this->_indices = make_shared<indices>(range(0,s._dim[0]-1));
            }
            return *this;
        }

        param in(const node_pairs& np){
            return this->in(np._keys);
        }

        void reverse_sign(){
            throw invalid_argument("Cannot reverse sign of param");
        }
        
        /** Index parameter/variable in ids, look for the keys starting at the ith position
         @param[in] start_position If ids has keys with additional entries, use the substring starting after the start_position comma separator
         @param[in] ids_ index set
         */
        template<typename... Args>
        param in_ith(unsigned start_position, const indices& ids_) {
            auto ids(ids_);
            if(!ids._excluded_keys.empty()){
                ids.remove_excluded();
            }
            if(!_indices){
                throw invalid_argument("unindexed param/var, first call in()");
            }
            param res(*this);
            if(ids.empty()){
                DebugOff("In function param.in_ith(unsigned position, const indices& ids), ids is empty!\n. Creating and empty variable! Check your sum/product operators.\n");
                res._name += "_EMPTY";
                res._dim[0] = 0;
                res._dim[1] = 0;
                //                res.set_range(0);
                res.reset_range();
                return res;
            }
            string key;
            size_t nb_inst=1;
            /** Number of comma separated keys in current variable */
            auto nb_sep_var = count(_indices->_keys->front().begin(), _indices->_keys->front().end(), ',');
            /** Number of comma separated keys in ids */
            auto nb_sep_ids = count(ids._keys->front().begin(), ids._keys->front().end(), ',');
            if(nb_sep_var > nb_sep_ids){/* ids can have more keys if it's a super set of indices */
                throw invalid_argument("Variable " + _name + ": In function in_ith(unsigned start_position, const indices& ids_), ids keys have less entries than current param/var, check indexing");
            }
            res._indices->_ids = make_shared<vector<vector<size_t>>>();
            res._indices->_ids->resize(1);
            nb_inst = ids.size();
            if(ids.is_indexed()){/* If ids has key references, use those */
                for(auto &key_ref: ids._ids->at(0)){
                    key = _indices->_keys->at(key_ref);
                    auto pos = nthOccurrence(key, ",", start_position);
                    if(pos!=0){
                        key = key.substr(pos+1);
                    }
                    pos = nthOccurrence(key, ",", nb_sep_var+1);
                    if(pos!=0){
                        key = key.substr(0,pos);
                    }
                    auto it = _indices->_keys_map->find(key);
                    if (it == _indices->_keys_map->end()){
                        throw invalid_argument("Variable " + _name + ": In function in_ith(unsigned start_position, const indices& ids_), an index set has unrecognized key: " + key);
                    }
                    res._indices->_ids->at(0).push_back(it->second);
                }
                
            }
            else {
                for(auto key: *ids._keys){
                    auto pos = nthOccurrence(key, ",", start_position);
                    if(pos!=0){
                        key = key.substr(pos+1);
                    }
                    pos = nthOccurrence(key, ",", nb_sep_var+1);
                    if(pos!=0){
                        key = key.substr(0,pos);
                    }
                    auto it = _indices->_keys_map->find(key);
                    if (it == _indices->_keys_map->end()){
                        throw invalid_argument("Variable " + _name + ": In function in_ith(unsigned start_position, const indices& ids_), an index set has unrecognized key: " + key);
                    }
                    res._indices->_ids->at(0).push_back(it->second);
                }
            }

            if(res._is_transposed){
                res._dim[1]=res._indices->_ids->at(0).size();
            }
            else {
                res._dim[0]=res._indices->_ids->at(0).size();
            }
            res._name += ".in("+ids.get_name()+")";
            res._indices->set_name(ids.get_name());
            res.reset_range();
            return res;
        }
        
        /** Index parameter/variable in the product of ids1...args
         */
        template<typename... Args>
        param in(const indices& ids1, Args&&... args) {
            auto ids = indices(ids1,args...);
            if(!ids._excluded_keys.empty()){
                ids.remove_excluded();
            }
            if(!_indices || _indices->empty()){/**< No need to add each key individually */
                _indices = make_shared<indices>(ids);
                auto dim = _indices->size();
                _val->resize(dim);
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
                param res(*this);
                res._name += ".in("+ids.get_name()+")";
                return res;
            }
            string key, excluded;
            size_t nb_inst=1;
            param res(*this);
            if(ids._ids && ids._ids->size()>1){/* Double-indexed set */
                nb_inst = ids._ids->size();
                /* Check that the current param has the keys found in ids */
                for(auto key: *ids._keys){
                    auto it1 = _indices->_keys_map->find(key);
                    if (it1 == _indices->_keys_map->end()){
                        throw invalid_argument("Variable " + _name + ": In function param.in(const indices& index_set1, Args&&... args), an index set has unrecognized key: " + key);
                    }
                }
                res._indices = make_shared<indices>(ids);
                if(res._is_transposed){
                    if(res.is_double_indexed()){
                        res._dim[1]=_indices->size();
                    }
                    else {
                        res._dim[1]=res._indices->size();
                    }
                }
                else {
                    if(res.is_double_indexed()){
                        res._dim[0]=_indices->size();
                    }
                    else {
                        res._dim[0]=res._indices->size();
                    }
                }
                res._name += ".in("+ids.get_name()+")";
                res.reset_range();
                return res;
            }
            res._indices->_ids = make_shared<vector<vector<size_t>>>();
            res._indices->_ids->resize(1);
            nb_inst = ids.size();
            if(ids.is_indexed()){/* If ids has key references, use those */
                for(auto &key_ref: ids._ids->at(0)){
                    key = _indices->_keys->at(key_ref);
                    auto it = _indices->_keys_map->find(key);
                    if (it == _indices->_keys_map->end()){
                        throw invalid_argument("Variable " + _name + ": In function param.in(const indices& index_set1, Args&&... args), an index set has unrecognized key: " + key);
                    }
                    res._indices->_ids->at(0).push_back(it->second);
                }
                
            }
            else {
                for(auto key: *ids._keys){
                    if(_indices->_type==to_){/** Assumed to be the last entry in the key */
                        key = key.substr(key.find_last_of(",")+1,key.size());
                    }
                    else if(_indices->_type==from_){/** Assumed to be the one to last entry in the key */
                        key = key.substr(0, key.find_last_of(","));
                        key = key.substr(key.find_last_of(",")+1,key.size());
                    }
                    auto it = _indices->_keys_map->find(key);
                    if (it == _indices->_keys_map->end()){
                        throw invalid_argument("Variable " + _name + ": In function param.in(const indices& index_set1, Args&&... args), an index set has unrecognized key: " + key);
                    }
                    res._indices->_ids->at(0).push_back(it->second);
                }
            }
            if(res._is_transposed){
                res._dim[1]=res._indices->_ids->at(0).size();
            }
            else {
                res._dim[0]=res._indices->_ids->at(0).size();
            }
            res._name += ".in("+ids.get_name()+")";
            res._indices->set_name(ids.get_name());
            if(!excluded.empty()){
                excluded = excluded.substr(0,excluded.size()-1); /* remove last comma */
                res._name += "\{" + excluded + "}";
            }
            res.reset_range();
            return res;
        }

        template<typename... Args>
        void index_in(const indices& ids1, Args&&... args) {
            *this = this->in(ids1,args...);
        }

        param in_pairs(const indices& ids) {
            param res(*this);
            res._indices->_ids = make_shared<vector<vector<size_t>>>();
            res._indices->_ids->resize(1);
            if(ids.empty()){
                DebugOff("In function param.in(const indices& index_set1, Args&&... args), all index sets are empty!\n. Creating and empty variable! Check your sum/product operators.\n");
                res._name += "_EMPTY";
                return res;
            }
            string key, excluded;
            size_t idx = 0;
            for(auto key: *ids._keys){
                key = key.substr(key.find_first_of(",")+1,key.size());
                if(ids._excluded_keys.count(idx++)!=0){
                    excluded += key + ",";
                    continue;
                }
                auto it1 = _indices->_keys_map->find(key);
                if (it1 == _indices->_keys_map->end()){
                    throw invalid_argument("In function param.in_pairs(const indices& ids), ids has unknown key");
                }
                res._indices->_ids->at(0).push_back(it1->second);
            }
            if(res._is_transposed){
                res._dim[1]=res._indices->_ids->at(0).size();
            }
            else {
                res._dim[0]=res._indices->_ids->at(0).size();
            }
            res._name += ".in_pairs("+ids.get_name()+")";
            if(!excluded.empty()){
                excluded = excluded.substr(0,excluded.size()-1); /* remove last comma */
                res._name += "\{" + excluded + "}";
            }
            res.reset_range();
            return res;
        }

        param from(){
            param<type> res(*this);
            res._name += ".from";
            res._indices->_type = from_;
            res._indices->set_name("from("+res._indices->get_name()+")");
            return res;
        }
        
        param from(const indices& ids){
            return this->from().in(ids);
        }
        
        
        param to(){
            param<type> res(*this);
            res._name += ".to";
            res._indices->_type = to_;
            res._indices->set_name("to("+res._indices->get_name()+")");
            return res;
        }
        
        param to(const indices& ids){
            return this->to().in(ids);
        }

        /* Index param/var based on incoming arcs out of nodes in vec */
        param in_arcs(const vector<Node*>& vec) {
            param res(*this);
            res._indices->_ids = make_shared<vector<vector<size_t>>>();
            res._indices->_ids->resize(1);
            if(vec.empty()){
                DebugOff("In function param.in_arcs(const vector<Node*>& vec), vec is empty!\n. Creating and empty variable! Check your sum/product operators.\n");
                res._name += "_EMPTY";
                return res;
            }
            string key;
            size_t inst = 0;
            for(auto it = vec.begin(); it!= vec.end(); it++) {
                if(!(*it)->_active) {
                    continue;
                }
                if (inst>0) {
                    res._indices->_ids->push_back(vector<size_t>());
                }
                for (auto &a:(*it)->get_in()) {
                    if (!a->_active) {
                        continue;
                    }
                    key = a->_name;
                    auto it1 = _indices->_keys_map->find(key);
                    if (it1 == _indices->_keys_map->end()){
                        throw invalid_argument("In function param.in_arcs(const vector<Node*>& vec), unknown arc key.");
                    }
                    res._indices->_ids->at(inst).push_back(it1->second);
                }
                ++inst;
            }
//            res._name += ".in_arcs";
            return res;
        }



        /* Index param/var based on outgoing arcs out of nodes in vec */
        param out_arcs(const vector<Node*>& vec) {
            param res(*this);
            res._indices->_ids = make_shared<vector<vector<size_t>>>();
            res._indices->_ids->resize(1);
            if(vec.empty()){
                DebugOff("In function param.in_arcs(const vector<Node*>& vec), vec is empty!\n. Creating and empty variable! Check your sum/product operators.\n");
                res._name += "_EMPTY";
                return res;
            }
            string key;
            size_t inst = 0;
            for(auto it = vec.begin(); it!= vec.end(); it++) {
                if(!(*it)->_active) {
                    continue;
                }
                if (inst>0) {
                    res._indices->_ids->push_back(vector<size_t>());
                }
                for (auto &a:(*it)->get_out()) {
                    if (!a->_active) {
                        continue;
                    }
                    key = a->_name;
                    auto it1 = _indices->_keys_map->find(key);
                    if (it1 == _indices->_keys_map->end()){
                        throw invalid_argument("In function param.in_arcs(const vector<Node*>& vec), unknown arc key.");
                    }
                    res._indices->_ids->at(inst).push_back(it1->second);
                }
                ++inst;
            }
//            res._name += ".out_arcs";
            return res;
        }


        /* Index param/var based on auxiliary objects attached to nodes in vec */
        param in_aux(const vector<Node*>& vec, const string& aux_type) {
            param res(*this);
            res._indices->_ids = make_shared<vector<vector<size_t>>>();
            if(vec.empty()){
                DebugOff("In function param.in_aux(const vector<Node*>& vec), vec is empty!\n. Creating and empty variable! Check your sum/product operators.\n");
                res._name += "_EMPTY";
                return res;
            }
            string key;
            size_t inst = 0;
            for(auto it = vec.begin(); it!= vec.end(); it++) {
                if(!(*it)->_active) {
                    continue;
                }
                res._indices->_ids->push_back(vector<size_t>());
                for (auto &a:(*it)->get_aux(aux_type)) {
                    if (!a->_active) {
                        continue;
                    }
                    key = a->_name;
                    auto it1 = _indices->_keys_map->find(key);
                    if (it1 == _indices->_keys_map->end()){
                        throw invalid_argument("In function param.in_aux(const vector<Node*>& vec), unknown arc key.");
                    }
                    res._indices->_ids->at(inst).push_back(it1->second);
                }
                ++inst;
            }
//            res._name += "("+aux_type+")";
            return res;
        }


        /** Output */

        string to_str(size_t index1, size_t index2, int prec) {
            if (is_matrix()){
                return to_string_with_precision(eval(index1,index2),prec);
            }
            if (is_indexed()) {
                return to_string_with_precision(eval(index1,index2),prec);
            }
            else {
                return to_string_with_precision(eval(index2),prec);
            }
        }

        string to_str(){
            return get_name(false,false);
        }
        
        string to_str(size_t index, int prec) {
            if (is_indexed()) {
                return to_string_with_precision(eval(index), prec);
            }
            else {
                return to_string_with_precision(eval(index), prec);
            }
        }
        
        
//        inline size_t get_id_inst(unsigned inst = 0) const {
//            if (is_indexed()) {
////                if(_indices->_ids->at(0).size() <= inst){
////                    throw invalid_argument("get_id_inst out of range");
////                }
//                return _indices->_ids->at(0).at(inst);
//            }
//            return inst;
//        };
//
//        size_t get_id_inst(unsigned inst1, unsigned inst2) const {
//            if (is_indexed()) {
//                if (_indices->_ids->size()==1) {
////                    if(_indices->_ids->at(0).size() <= inst2){
////                        throw invalid_argument("get_id_inst out of range");
////                    }
//                    return _indices->_ids->at(0).at(inst2);
//                }
//                return _indices->_ids->at(inst1).at(inst2);
//            }
//            return inst2;
//        };
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
        
        string to_str_vals(bool vals, int prec = 10) {
            string str = get_name(false,true);
            auto name = str.substr(0, str.find_last_of("."));
            str = name;
            if (vals) {
                str += " = { \n";
                auto space_size = str.size();
                if (is_matrix()) {
                    auto max_cell_size = get_max_cell_size();
                    for (size_t i = 0; i<_dim[0]; i++) {
                        str.insert(str.end(), space_size, ' ');
                        str += "|";
                        for (size_t j = 0; j<_dim[1]; j++) {
                            auto cell = to_str(i,j,prec);
                            auto cell_size = cell.size();
                            cell.insert(0, floor((max_cell_size+1 - cell_size)/2.), ' ');
                            cell.append(ceil((max_cell_size+1 - cell_size)/2.), ' ');
                            str += cell;
                            if(j!=_dim[1]-1){
                                str += " ";
                            }
                        }
                        str += "|\n";
                    }
                    str += "}\n";
                    return str;
                }
                if(_indices) {
                    if (is_indexed()) {
                        for (size_t i = 0; i < _dim[0]; i++) {
                            str += "[" + _indices->_keys->at(get_id_inst(i)) + "] = " + to_string_with_precision(eval(i), prec);
                            str += " \n";
                        }
                    }
                    else {
                        for (size_t i = 0; i < _dim[0]; i++) {
                            str += "[" + _indices->_keys->at(i) + "] = " + to_string_with_precision(eval(i), prec);
                            str += " \n";
                        }
                    }
                }
                else {
                    for (size_t idx = 0; idx < _val->size(); idx++) {
                        str += "["+to_string(idx) + "] = " + to_string_with_precision(eval(idx),prec);
                        str += " \n";
                    }
                }
                str += "};\n";
            }
            return str;
        }
        
        
        void print() {
            print(16);
        }

        void print_vals(int prec){
            cout << this->to_str_vals(true, prec);
        }
        
        void print(int prec) {
            print_vals(prec);
        }


        type getvalue() const {
            return _val->back();
        }
        
        void update_range(const Cpx& val);
        /**
         Recompute range based on stored values.
         */
        void reset_range(){
            init_range();
            if(is_double_indexed()){
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
//        void set_vals(const Eigen::SparseMatrix<Cpx,Eigen::RowMajor>& SM);
        
        /**
         Initialize x with model variables values
         @param[out] x values to initialize
         */
        void get_x(const vector<double>& x){set_x_(x);};
        
        template<class T=type, typename enable_if<is_arithmetic<T>::value>::type* = nullptr>
        void get_x_(const vector<double>& x){
            auto vid = get_id();
            for (size_t i = 0; i < get_dim(); i++) {
                x[vid+i] = (double)_val->at(i);
            }
        }
        
        template<class T=type, typename enable_if<is_same<T, Cpx>::value>::type* = nullptr>
        void get_x_(const vector<double>& x){
            throw invalid_argument("Cannot call set_double_val_ with a non-arithmetic type.");
        }
        
        void set_double_val(double* x){set_double_val_(x);};
            
        template<class T=type, typename enable_if<is_arithmetic<T>::value>::type* = nullptr>
        void set_double_val_(double* x){
            auto vid = get_id();
            for (size_t i = 0; i < get_dim(); i++) {
                x[vid+i] = (double)_val->at(i);
            }
        }
        
        template<class T=type, typename enable_if<is_same<T, Cpx>::value>::type* = nullptr>
        void set_double_val_(double* x){
            throw invalid_argument("Cannot call set_double_val_ with a non-arithmetic type.");
        }
        
        void get_double_val(const double* x){get_double_val_(x);};
        
        template<class T=type, typename enable_if<is_arithmetic<T>::value>::type* = nullptr>
        void get_double_val_(const double* x){
            auto vid = get_id();
            for (size_t i = 0; i < get_dim(); i++) {
                _val->at(i) = x[vid+i];
            }
        }
        
        template<class T=type, typename enable_if<is_same<T, Cpx>::value>::type* = nullptr>
        void get_double_val_(const double* x){
            throw invalid_argument("Cannot call get_double_val_ with a non-arithmetic type.");
        }
        
        /**
         Initialize the model variables using values from x
         @param[in] x values to initialize to
         */
        void set_x(const vector<double>& x){set_x_(x);};
        
        template<class T=type, typename enable_if<is_arithmetic<T>::value>::type* = nullptr>
        void set_x_(const vector<double>& x){
            auto vid = get_id();
            for (size_t i = 0; i < get_dim(); i++) {
                _val->at(i) = x[vid+i];
            }
        }
        
        template<class T=type, typename enable_if<is_same<T, Cpx>::value>::type* = nullptr>
        void set_x_(const vector<double>& x){
            throw invalid_argument("Cannot call get_double_val_ with a non-arithmetic type.");
        }
        
        
        void get_double_val(size_t pos, double x){get_double_val_(pos,x);};
        
        template<class T=type, typename enable_if<is_arithmetic<T>::value>::type* = nullptr>
        void get_double_val_(size_t pos, double x){
            _val->at(pos) = x;
        }
        
        template<class T=type, typename enable_if<is_same<T, Cpx>::value>::type* = nullptr>
        void get_double_val_(size_t i, double x){
            throw invalid_argument("Cannot call get_double_val_ with a non-arithmetic type.");
        }
        
        void set_double_val(size_t pos, double& x){set_double_val_(pos,x);};
        
        template<class T=type, typename enable_if<is_arithmetic<T>::value>::type* = nullptr>
        void set_double_val_(size_t pos, double& x){
            x = _val->at(pos);
        }
        
        template<class T=type, typename enable_if<is_same<T, Cpx>::value>::type* = nullptr>
        void set_double_val_(size_t i, double& x){
            throw invalid_argument("Cannot call get_double_val_ with a non-arithmetic type.");
        }
        
        /** round the value stored at position i to the nearest integer */
        void round_vals(){round_vals_();};
        
        template<class T=type, typename enable_if<is_arithmetic<T>::value>::type* = nullptr>
        void round_vals_(){
            for (size_t i = 0; i < get_dim(); i++) {
                _val->at(i) = std::round((T)_val->at(i));
            }
        };
        
        template<class T=type, typename enable_if<is_same<T, Cpx>::value>::type* = nullptr>
        void round_vals_(){
            for (size_t i = 0; i < get_dim(); i++) {
                _val->at(i).real(std::round(_val->at(i).real()));
                _val->at(i).imag(std::round(_val->at(i).imag()));
            }
        };
        
        
        
        template<typename T, typename=enable_if<is_convertible<T,type>::value>>
        void copy_vals(const param<T>& p){
            auto dim = p.get_dim();
            _val->resize(dim);
            for (size_t i = 0; i < dim; i++) {
                _val->at(i) = p._val->at(i);
            }
        }
        
        void copy_vals(const shared_ptr<param_>& p);
                
       


    };

    param<Cpx> conj(const param<Cpx>& p);
    param<Cpx> ang(const param<Cpx>& p);
    param<Cpx> sqrmag(const param<Cpx>& p);
    param<Cpx> real(const param<Cpx>& p);
    param<Cpx> imag(const param<Cpx>& p);

    template<typename type>
    param<type> diag(const param<type>& p){
        param<type> res("diag("+p._name+")");
        if(p.is_matrix()){
            res.set_size(min(p._dim[0], p._dim[1]));
            for (auto i = 0; i<res._dim[0]; i++) {
                res.set_val(i,p.eval(i,i));
            }
        }
        else{
            res.set_size(p._dim[0], p._dim[0]);
            for (auto i = 0; i<res._dim[0]; i++) {
                res.set_val(i,i,p.eval(i));
            }
        }
        return res;
    }

}
#endif /* defined(____param__) */
