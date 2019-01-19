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
#include <Eigen/Sparse>

using namespace std;





namespace gravity {    
    
    /** Backbone class for parameter */
    class param_: public constant_ {

    protected:

        NType                                          _intype; /**< internal storage type **/

        void set_type(NType type) {
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

        
        virtual ~param_(){};
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
            _indices = p._indices;
            _dim[0] = p._dim[0];
            _dim[1] = p._dim[1];
        }

        virtual shared_ptr<param_> pcopy() const{return nullptr;};

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
        
        size_t get_id_inst(size_t inst = 0) const {            
            if (is_indexed()) {
                if(_indices->_ids->at(0).size() <= inst){
                    throw invalid_argument("param::get_id_inst(size_t inst) inst is out of range");
                }
                return _indices->_ids->at(0).at(inst);
            }
            auto dim = get_dim();
            if(inst > dim-1){
                throw invalid_argument("param::get_id_inst(size_t inst) inst is out of range");
            }
            return inst;
        };

        size_t get_id_inst(size_t inst1, size_t inst2) const {
            if (is_indexed()) {
                if (_indices->_ids->size()==1) {
                    return _indices->_ids->at(0).at(inst2);
                }
                return _indices->_ids->at(inst1).at(inst2);
            }
            return get_id_inst(inst2);
//            throw invalid_argument("Calling get_id_inst(size_t inst1, size_t inst2) on a non-indexed param\n");
        };


        string get_name(bool in_func, bool exclude_indexing) const{
            string name = _name;
            if(_indices && exclude_indexing){
                name = name.substr(0, name.find_last_of("."));
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
                name += "["+to_string(get_id_inst(inst))+"]";
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
                name += "["+_indices->_keys->at(inst2)+"]";
            }
            else{
                name += "["+to_string(inst2)+"]";
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

        bool is_indexed() const{
            return (_indices && _indices->_ids);
        }

        bool is_binary() const {
            return (_intype==binary_);
        };

        bool is_integer() const {
            return (_intype==integer_);
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

        Sign get_all_sign() const{return unknown_;}; /**< If all instances of the current parameter/variable have the same sign, it returns it, otherwise, it returns unknown. **/
//        virtual Sign get_sign(size_t idx = 0) const; /**< returns the sign of one instance of the current parameter/variable. **/
        

        /** Operators */
        
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
        
        /**
         Index the current object using incoming edges for nodes stored in vec. This is a double indexing where each row corresponds to a node, and columns correspond to the edge ids.
         @param[in] vec vector of nodes
         */
        void index_in_arcs(const vector<Node*>& vec) {
            _indices->_ids = make_shared<vector<vector<size_t>>>();
            if(vec.empty()){
                DebugOn("In function param.in_arcs(const vector<Node*>& vec), vec is empty!\n. Creating and empty variable! Check your sum/product operators.\n");
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
                DebugOn("In function param.index_out_arcs(const vector<Node*>& vec), vec is empty!\n. Creating and empty variable! Check your sum/product operators.\n");
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
                DebugOn("In function param.index_in_aux(const vector<Node*>& vec), vec is empty!\n. Creating and empty variable! Check your sum/product operators.\n");
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
            if (is_indexed() && _indices->_ids->size()>i) {
                return _indices->_ids->at(i).size();
            }
            return constant_::get_dim(i);
        }

    };


    /** A parameter can be a bool, a short, an int, a float, a double, a long double or a complex<double>. */
    template<typename type = double>
    class param: public param_ {

    public:

        shared_ptr<vector<type>>                _val = nullptr; /**< vector of values **/
        shared_ptr<pair<type,type>>             _range = nullptr; /**< (Min,Max) values in vals **/

        template<typename T=type,
        typename std::enable_if<is_arithmetic<T>::value>::type* = nullptr>
        param() {
            update_type();
            _val = make_shared<vector<type>>();
            _range = make_shared<pair<type,type>>(make_pair<>(numeric_limits<type>::max(), numeric_limits<type>::lowest()));
        }

        template<class T=type, class = typename enable_if<is_same<T, Cpx>::value>::type>
        param(){
            update_type();
            _val = make_shared<vector<type>>();
            _range = make_shared<pair<type,type>>(make_pair<>(Cpx(numeric_limits<double>::max(), numeric_limits<double>::max()), Cpx(numeric_limits<double>::lowest(), numeric_limits<double>::lowest())));
        }
        
        shared_ptr<param_> pcopy() const{return make_shared<param>(*this);};

        shared_ptr<constant_> copy()const{return make_shared<param>(*this);};
        
        ~param(){};
        template<class T2, typename std::enable_if<is_convertible<T2, type>::value && sizeof(T2) < sizeof(type)>::type* = nullptr>
        param (const param<T2>& p) {
            *this = p;            
        }

        param (const param& p) {
            *this = p;
        }
        
        param (param&& p) {
            *this = move(p);
        }

        param& operator=(const param& p) {
            _type = p._type;
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
            if(p._indices){
                _indices = make_shared<indices>(*p._indices);
            }
            _dim[0] = p._dim[0];
            _dim[1] = p._dim[1];
            return *this;
        }
        
        template<class T2, typename std::enable_if<is_convertible<T2, type>::value && sizeof(T2) < sizeof(type)>::type* = nullptr>
        param& operator=(const param<T2>& p) {
            update_type();
            _id = p._id;
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
            if(p._indices){
                _indices = make_shared<indices>(*p._indices);
            }
            _dim[0] = p._dim[0];
            _dim[1] = p._dim[1];
            return *this;
        }

        param& operator=(param&& p) {
            _type = p._type;
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
            return _val->at(param_::_indices->_keys_map->at(key));
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
                throw invalid_argument("Cannot call param::add_val(type val) on matrix");
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
                throw invalid_argument("Cannot call param::add_val(type val) on matrix");
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
            auto index = param_::_indices->size();
            auto pp = param_::_indices->_keys_map->insert(make_pair<>(key,index));
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



        template<class T=type, typename enable_if<is_arithmetic<T>::value>::type* = nullptr> bool is_unit() const { /**< Returns true if all values of this paramter are 1 **/
            return (_range->first == 1 && _range->second == 1);
        }
        
        template<class T=type, class = typename enable_if<is_same<T, Cpx>::value>::type> bool is_unit() const{
            return (_range->first == Cpx(1,1) && _range->second == Cpx(1,1));
        }

        template<class T=type, typename enable_if<is_arithmetic<T>::value>::type* = nullptr> bool is_zero() const { /**< Returns true if all values of this paramter are 0 **/
            return (_range->first == 0 && _range->second == 0);
        }

        template<class T=type, class = typename enable_if<is_same<T, Cpx>::value>::type> bool is_zero() const{
            return (_range->first == Cpx(0,0) && _range->second == Cpx(0,0));
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


        template<class T=type, class = typename enable_if<is_arithmetic<T>::value>::type> void set_vals(const Eigen::SparseMatrix<Cpx,Eigen::RowMajor>& SM){
            for (size_t k=0; k<SM.outerSize(); ++k) {
                for (Eigen::SparseMatrix<Cpx,Eigen::RowMajor>::InnerIterator it(SM,k); it; ++it){
                    set_val(2*it.row(), 2*it.col(), it.value().real());
                    set_val(2*it.row()+1, 2*it.col()+1, it.value().real());
                    set_val(2*it.row()+1, 2*it.col(), it.value().imag());
                    set_val(2*it.row(), 2*it.col()+1, -it.value().imag());
                }
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
                res._indices->_ids->at(0).push_back(i*_dim[1]+j);
            }
            return res;
        }

        param operator()(size_t idx) {
            if(!_indices){
                throw invalid_argument("Current param/var is not indexed.");
            }
            param res(*this);
            res._indices->_ids = make_shared<vector<vector<size_t>>>();
            res._indices->_ids->resize(1);
            res._indices->_ids->at(0).push_back(idx);
            res._dim[0] = 1;
            return res;
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
                this->_indices = make_shared<indices>(indices(0,s._dim[0]-1));
            }
            return *this;
        }

        param in(const node_pairs& np){
            return this->in(np._keys);
        }

        

        template<typename... Args>
        param in(const indices& vec1, Args&&... args) {
            auto ids = indices(vec1,args...);
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
                res._name += ".in("+ids._name+")";
                return res;
            }
            param res(*this);
            res._indices->_ids = make_shared<vector<vector<size_t>>>();
            res._indices->_ids->resize(1);
            if(ids.empty()){
                DebugOn("In function param.in(const indices& index_set1, Args&&... args), all index sets are empty!\n. Creating and empty variable! Check your sum/product operators.\n");
                res._name += "_EMPTY";
                return res;
            }
            string key, excluded;
            size_t idx = 0;
            for(auto key: *ids._keys){
                if(ids._excluded_keys.count(idx++)!=0){
                    excluded += key + ",";
                    continue;
                }
                if(!res._indices->_time_extended && ids._time_extended){/* truncate time indices */
                    auto pos = nthOccurrence(key, ",", ids._time_pos);
                    key = key.substr(0,pos);
                }
                if(_indices){
                    if(_indices->_type==to_){
                        key = key.substr(key.find_last_of(",")+1,key.size());
                    }
                    else if(_indices->_type==from_){
                        key = key.substr(0, key.find_last_of(","));
                        key = key.substr(key.find_last_of(",")+1,key.size());
                    }
                }
                /* Compare indexing and truncate extra indices */
                auto nb_sep1 = count(_indices->_keys->front().begin(), _indices->_keys->front().end(), ',');
                auto nb_sep2 = count(key.begin(), key.end(), ',');
                if(nb_sep2>nb_sep1){
                    auto pos = nthOccurrence(key, ",", nb_sep2-nb_sep1);
                    key = key.substr(pos+1,key.size()-1);
                }
                auto it1 = _indices->_keys_map->find(key);
                if (it1 == _indices->_keys_map->end()){
                    throw invalid_argument("In function param.in(const vector<Tobj>& vec), vec has unknown key");
                }
                res._indices->_ids->at(0).push_back(it1->second);
            }
            if(res._is_transposed){
                res._dim[1]=res._indices->_ids->at(0).size();
            }
            else {
                res._dim[0]=res._indices->_ids->at(0).size();
            }

            res._name += ".in("+ids._name+")";
            if(!excluded.empty()){
                excluded = excluded.substr(0,excluded.size()-1); /* remove last comma */
                res._name += "\{" + excluded + "}";
            }
            return res;
        }


        param in_pairs(const indices& ids) {
            param res(*this);
            res._indices->_ids = make_shared<vector<vector<size_t>>>();
            res._indices->_ids->resize(1);
            if(ids.empty()){
                DebugOn("In function param.in(const indices& index_set1, Args&&... args), all index sets are empty!\n. Creating and empty variable! Check your sum/product operators.\n");
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
            res._name += ".in_pairs("+ids._name+")";
            if(!excluded.empty()){
                excluded = excluded.substr(0,excluded.size()-1); /* remove last comma */
                res._name += "\{" + excluded + "}";
            }
            return res;
        }

        param in_pairs(){
            param<type> res(*this);
            res._name += ".in_pairs";
            res._indices->_type = in_pairs_;
            return res;
        }

        param from(){
            param<type> res(*this);
            res._name += ".from";
            res._indices->_type = from_;
            return res;
        }


        param to(){
            param<type> res(*this);
            res._name += ".to";
            res._indices->_type = to_;
            return res;
        }


        param out_arcs(){
            param<type> res(*this);
            res._name += ".out_arcs";
            res._indices->_type = out_arcs_;
            return res;
        }

        param in_arcs(){
            param<type> res(*this);
            res._name += ".in_arcs";
            res._indices->_type = in_arcs_;
            return res;
        }

        param in_gens(){
            param<type> res(*this);
            res._name += ".in_gens";
            res._indices->_type = in_gens_;
            return res;
        }

        /* Index param/var based on incoming arcs out of nodes in vec */
        param in_arcs(const vector<Node*>& vec) {
            param res(*this);
            res._indices->_ids = make_shared<vector<vector<size_t>>>();
            res._indices->_ids->resize(1);
            if(vec.empty()){
                DebugOn("In function param.in_arcs(const vector<Node*>& vec), vec is empty!\n. Creating and empty variable! Check your sum/product operators.\n");
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
                DebugOn("In function param.in_arcs(const vector<Node*>& vec), vec is empty!\n. Creating and empty variable! Check your sum/product operators.\n");
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
                DebugOn("In function param.in_aux(const vector<Node*>& vec), vec is empty!\n. Creating and empty variable! Check your sum/product operators.\n");
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

        string to_str(size_t index1, size_t index2, int prec) const {
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

        string to_str() const{
            return get_name(false,false);
        }
        
        string to_str(size_t index, int prec) const {
            if (is_indexed()) {
                return to_string_with_precision(eval(index), prec);
            }
            else {
                return to_string_with_precision(eval(index), prec);
            }
        }

        void print(size_t index, int prec = 10) const {
            cout << to_str(index,prec);
        }

        void print(size_t i, size_t j, int prec = 10) const {
            cout << to_str(i,j,prec);
        }
        
        size_t get_id_inst(unsigned inst = 0) const {
            if (is_indexed()) {
                if(_indices->_ids->at(0).size() <= inst){
                    throw invalid_argument("get_id_inst out of range");
                }
                return _indices->_ids->at(0).at(inst);
            }
            return inst;
        };
        
        size_t get_id_inst(unsigned inst1, unsigned inst2) const {
            if (is_indexed()) {
                if (_indices->_ids->size()==1) {
                    if(_indices->_ids->at(0).size() <= inst2){
                        throw invalid_argument("get_id_inst out of range");
                    }
                    return _indices->_ids->at(0).at(inst2);
                }
                return _indices->_ids->at(inst1).at(inst2);
            }
            return inst2;
        };

        string to_str_vals(bool vals, int prec = 10) const {
            string str = get_name(false,true);
            auto name = str.substr(0, str.find_last_of("."));
            str = name;
            if (vals) {
                str += " = { \n";
                if(is_matrix()){
                    for (size_t i = 0; i < _dim[0]; i++) {
                        for (size_t j = 0; j < _dim[1]; j++) {
                            str += to_string_with_precision(this->eval(i,j),prec);
                            str += " ";
                        }
                        str += "\n";
                    }
                    str += "};\n";
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

        void print(bool vals=true, int prec = 10) const {
            cout << this->to_str_vals(vals, prec);
        }


        type getvalue() const {
            if (_indices) {
                return (_val->at(_indices->last()));
            }
            else {
                return _val->at(0);
            }
        }
        
        void update_range(const Cpx& val);
//        void set_vals(const Eigen::SparseMatrix<Cpx,Eigen::RowMajor>& SM);

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
