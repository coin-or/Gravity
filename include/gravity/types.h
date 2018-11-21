//
//  Type.h
//  Gravity
//
//  Created by Hassan on 3 Jan 2016.
//

#ifndef Gravity___Type_h
#define Gravity___Type_h
#include <memory>
#include <list>
#include <map>
#include <set>
#include <assert.h>
#include <string>
#include <iostream>
#include <algorithm>

namespace gravity{
#define EPS 0.00001
#define Cpx complex<double>
//#define Real double
//#define Integer integer  //same name with a boost graph library.
#define Binary bool
#define Debug(x)
#define DebugOn(x) cout << x
#define Warning(x)
#define DebugOff(x)

    typedef unsigned int uint; /* Index type */
    //typedef std::set<ind> indx; /* Set of indices type */

    typedef enum { linear_, convex_, concave_, undet_} Convexity; /* Convexity Type */
    typedef enum { neg_ = -2, non_pos_ = -1, zero_ = 0, non_neg_ = 1, pos_ = 2, unknown_ = 3} Sign; /* Sign Type */
    typedef enum { binary_, short_, integer_, float_, double_, long_, complex_} NType; /* Number Type */
    typedef enum { binary_c, short_c, integer_c, float_c, double_c, long_c, par_c, uexp_c, bexp_c, var_c, func_c, complex_c} CType; /* Constant type, ancestor to parameter, var and function */
    typedef enum { Unknown,
        PrimalAndDualFeasible,
        PrimalFeasible,
        DualFeasible,
        PrimalInfeasible,
        DualInfeasible,
        PrimalAndDualInfeasible,
        IllPosed,
        PrimalInfeasibleOrUnbounded } Outcome;
    typedef enum { geq, leq, eq } ConstraintType;
    typedef enum { const_, lin_, quad_, pol_, nlin_ } FType;  /* Function type in constraint: Constant, Linear, Quadratic, Polynomial or Nonlinear function */
    typedef enum { lin_m, quad_m, pol_m, nlin_m } MType;  /* Model type: Linear, Quadratic, Polynomial or Nonlinear function */
    typedef enum { minimize, maximize } ObjectiveType;
    typedef enum { id_, number_, plus_, minus_, product_, div_, power_, cos_, sin_, sqrt_, exp_, log_} OperatorType;  /* Operation type in the expression tree */

    typedef enum { R_, R_p_, C_} SpaceType;  /* Real, Positive Reals, Complex */

    typedef enum { ordered_pairs_, unordered_ } SetType;
//    typedef enum { vec_=0, in_ordered_pairs_=1, from_ordered_pairs_=2, to_ordered_pairs_=3, in_arcs_=4, from_arcs_=5, to_arcs_=6, in_nodes_=7, in_set_=8, mask_=9, in_bags_=10, time_expand_ = 11, in_set_at_} IndexType;  /* Index type */

    typedef enum { unindexed_, in_, in_pairs_, out_, from_, to_, prev_, in_at_, in_time_, from_time_, to_time_, in_arcs_, out_arcs_, in_gens_, in_pot_gens_, in_bats_, in_pot_bats_,in_wind_, in_pv_, min_time_, excl_} IndexType;  /* Index type */

    using namespace std;


    /** Class for manipulating indices */
    class index_{
    public:
        string _name;
        string _type_name="indices";
        bool   _active = true;
        index_(const string& name, bool active=true):_name(name), _active(active){};
        index_(const index_& idx):_name(idx._name), _active(idx._active){};
        template<typename... Args>
        index_(string t1, Args&&... args) {
            list<string> indices;
            indices = {forward<Args>(args)...};
            indices.push_front(t1);
            auto it = indices.begin();
            for (size_t i= 0; i < indices.size(); i++) {
                _name += *it;
                if (i< indices.size()-1) {
                    _name += ",";
                }
                it++;
            }
        }
    };
    
    class index_pair{
    public:
        string _name;
        string _type_name="index_pairs";
        bool   _active = true;
        index_* _src = nullptr;
        index_* _dest = nullptr;
        index_pair(const index_& src, const index_& dest, bool active = true):_name(src._name+","+dest._name), _active(active), _src(new index_(src)), _dest(new index_(dest)){};
        index_pair(const string& src, const string& dest, bool active = true):_name(src+","+dest), _active(active), _src(new index_(src)), _dest(new index_(dest)){};
        index_pair(const index_pair& p):_name(p._name), _active(p._active), _src(new index_(*p._src)), _dest(new index_(*p._dest)){};
        ~index_pair(){
            delete _src;
            _src = nullptr;
            delete _dest;
            _dest = nullptr;
        }
    };
    
    class ordered_pairs{
        
    public:
        size_t          _first;
        size_t          _last;
        string          _type_name="ordered_pairs";
        std::vector<index_pair*> _keys;
        ordered_pairs(size_t p1 ,size_t p2){
            _first = p1;
            _last = p2;
            auto n = p2 - p1 + 1;
            assert(n >= 0);
            _keys.resize(n*(n-1)/2);
            string key;
            size_t index = 0;
            for (int i = p1-1; i < p2; i++){
                for (int j = i+1; j < p2; j++){
                    _keys[index++] = new index_pair(index_(to_string(i)), index_(to_string(j)));
                }
            }
        }
        ordered_pairs(size_t p1 ,size_t p2, bool rev){
            _first = p1;
            _last = p2;
            auto n = p2 - p1 + 1;
            assert(n >= 0);
            _keys.resize(n*(n-1)/2);
            string key;
            size_t index = 0;
            for (int i = p1-1; i < p2; i++){
                for (int j = i+1; j < p2; j++){
                    _keys[index++] = new index_pair(index_(to_string(j)), index_(to_string(i)));
                }
            }
        }
        ~ordered_pairs(){
            for (auto p: _keys) { delete p;}
        }
    };
    
    class indices{
        
    public:
        
        string                                  _name;/**< index set can be given a name */
        IndexType                               _type = unindexed_;/**< index type */
        bool                                    _time_extended = false;/*<< indices are time extended */
        size_t                                  _time_pos = 0;/*<< number of commas before time extension */
        shared_ptr<vector<string>>              _keys = nullptr; /*<< A vector storing all the keys */
        
        shared_ptr<map<string,size_t>>          _keys_map = nullptr; /*<< A map storing all the indices, the size_t number indicates the right position in the _keys vector */

        set<size_t>                             _excluded_keys; /*<< A set storing all indices that should be excluded */
        shared_ptr<vector<vector<size_t>>>      _ids = nullptr;

        
        indices(string name){
            _name = name;
            _keys_map = make_shared<map<string,size_t>>();
            _keys = make_shared<vector<string>>();
        }
        
        indices(){
            _keys_map = make_shared<map<string,size_t>>();
            _keys = make_shared<vector<string>>();
        }
        
        indices(const ordered_pairs& pairs){
            auto n = pairs._keys.size();
            _keys_map = make_shared<map<string,size_t>>();
            _keys = make_shared<vector<string>>();
            _keys->resize(n);
            size_t index = 0;
            string key;
            for (int i = 0; i < n; i++){
                key = pairs._keys.at(index)->_name;
                (*_keys_map)[key]= index;
                (*_keys)[index++] = key;
            }
        }
        
        
        indices(size_t p1 ,size_t p2){
            auto n = p2 - p1 + 1;
            assert(n >= 0);
            _keys_map = make_shared<map<string,size_t>>();
            _keys = make_shared<vector<string>>();
            _keys->resize(n);
            size_t index = 0;
            for (int i = p1; i <= p2; i++){
                (*_keys_map)[to_string(i)]= index;
                (*_keys)[index++] = to_string(i);
            }
        }
        
        
        
        template<typename... Args>
        indices(string idx1, Args&&... args) {
            list<string> indices;
            indices = {forward<string>(args)...};
            indices.push_front(idx1);
            auto n = indices.size();
            _keys_map = make_shared<map<string,size_t>>();
            _keys = make_shared<vector<string>>();
            _keys->resize(n);
            auto it = indices.begin();
            for (size_t i= 0; i< n; i++) {
                (*_keys_map)[*it]= i;
                (*_keys)[i] = (*it);
                it++;
            }
        }
        
        
        indices& operator=(const indices& cpy){
            _name = cpy._name;
            _type = cpy._type;
            _keys_map = cpy._keys_map;
            _excluded_keys = cpy._excluded_keys;
            _keys = cpy._keys;
            if(cpy._ids){
                _ids = make_shared<vector<vector<size_t>>>(*cpy._ids);
            }
            _time_extended = cpy._time_extended;
            _time_pos = cpy._time_pos;
            return *this;
        }
        
        indices& operator=(indices&& cpy){
            if(!cpy._name.empty())
                _name = cpy._name;
            _type = cpy._type;
            _keys_map = move(cpy._keys_map);
            _excluded_keys = move(cpy._excluded_keys);
            _keys = move(cpy._keys);
            _ids = move(cpy._ids);
            _time_extended = cpy._time_extended;
            _time_pos = cpy._time_pos;
            return *this;
        }
        
        indices(const indices& cpy){
            *this=cpy;
        }
        
        indices(indices&& cpy){
            *this=move(cpy);
        }
        
//        template<typename Tobj>
//        indices(const vector<Tobj*>& vec){
//            _keys_map = make_shared<map<string,size_t>>();
//            _keys = make_shared<vector<string>>();
//            size_t i = 0;
//            for (auto idx:vec) {
//                if(idx->_active){
//                    _keys->push_back(idx->_name);
//                    (*_keys_map)[idx->_name]= i;
//                    i++;
//                }
//            }
//        }
//        
//        template<typename Tobj>
//        indices(const vector<Tobj>& vec){
//            _keys_map = make_shared<map<string,size_t>>();
//            _keys = make_shared<vector<string>>();
//            size_t i = 0;
//            for (auto idx:vec) {
//                if(idx._active){
//                    _keys->push_back(idx._name);
//                    (*_keys_map)[idx._name]= i;                    
//                    i++;
//                }
//            }
//        }

        template<typename Tobj>
        indices(const vector<Tobj*>& vec){
            _keys_map = make_shared<map<string,size_t>>();
            _keys = make_shared<vector<string>>();
            size_t i = 0;
            for (auto idx:vec) {
                if(idx->_active){
                    (*_keys_map)[idx->_name]= i;
                    _keys->push_back(idx->_name);
                    i++;
                }
            }
            if (_keys->size()>0) {
                _name = vec.front()->_type_name;
            }
        }
        
        template<typename Tobj>
        indices(const vector<Tobj>& vec){
            _keys_map = make_shared<map<string,size_t>>();
            _keys = make_shared<vector<string>>();
            size_t i = 0;
            for (auto idx:vec) {
                if(idx._active){
                    (*_keys_map)[idx._name]= i;
                    _keys->push_back(idx._name);
                    i++;
                }
            }
            if (_keys->size()>0) {
                _name = vec.front()->_type_name;
            }
        }
        
        template<typename... Args>
        indices(const indices& vec1, Args&&... args) {
            _keys_map = make_shared<map<string,size_t>>();
            _keys = make_shared<vector<string>>();
            list<indices> vecs;
            vecs = {forward<Args>(args)...};
            vecs.push_front(vec1);
            size_t dim = 1;
            size_t time_pos= 0, nb_ids = 0;
            vector<size_t> dims;
            for(auto &vec: vecs){
                if(vec.empty()){
                    Warning("\n\nWARNING: Defining indices with an empty vector!\n\n");
//                    exit(-1);
                }
                if(vec._time_extended){
                    _time_extended = true;
                    time_pos = vec._time_pos;
                }
                else{
                    nb_ids++;
                }
                dim *= vec.size();
                dims.push_back(vec.size());
                _name += vec._name+",";
            }
            _name = _name.substr(0,_name.size()-1); /* remove last comma */
            if(_time_extended && !vec1._time_extended){
                if(time_pos==0 && !vec1.empty()) {//TODO CHECK
                    _time_pos = std::count(vec1._keys->front().begin(), vec1._keys->front().end(), ',')+1;
                }
                else {
                    _time_pos = time_pos+nb_ids;
                }
            }
            _keys->resize(dim);
            size_t den = 1;
            size_t real_idx = 0;
            bool excluded = false;
            for(size_t idx = 0; idx < dim ; idx++){
                string key;
                den = dim;
                excluded = false;
                for(auto it = vecs.begin(); it!= vecs.end(); it++) {
                    auto vec = &(*it);
                    den /= vec->size();
                    real_idx = (idx/den)%vec->size();
                    if (vec->_excluded_keys.count(real_idx)==1) {
                        excluded = true;
                    }
                    key += vec->_keys->at(real_idx);
                    if(next(it)!=vecs.end()){
                        key += ",";
                    }
                }
                (*_keys)[idx] = key;
                (*_keys_map)[key] = idx;
                if (excluded) {
                    _excluded_keys.insert(idx);
                }
            }
        }
        template<typename... Args>
        void add(string s1, Args&&... args) {
            list<string> indices;
            indices = {forward<Args>(args)...};
            indices.push_front(s1);
            auto it = indices.begin();
            for (size_t i= 0; i < indices.size(); i++) {
                auto idx = _keys->size();
                auto pp = _keys_map->insert(make_pair<>(*it,idx));
                if (pp.second) {//new index inserted
                    _keys->push_back(*it);
                }
                else{
                    if(!_ids){
                        _ids = make_shared<vector<vector<size_t>>>();
                        _ids->resize(1);
                    }
                    _ids->at(0).push_back(pp.first->second);
                }
                it++;
            }
        }
        bool is_indexed() const{
            return (_ids!=nullptr);
        }
        
        void reindex(){
            _ids->at(0).clear();
            for(auto idx = 0; idx<_keys->size();idx++){
                if(_excluded_keys.count(idx)==0){
                    _ids->at(0).push_back(idx);
                }
            }
        }
        
        indices exclude(string key){
            auto res =  *this;
            res._excluded_keys.insert(res._keys_map->at(key));
            if(!is_indexed()){
                res._ids = make_shared<vector<vector<size_t>>>();
                res._ids->resize(1);
            }
            res.reindex();
            return res;
        }
        
        size_t size() const {
            if(is_indexed()){
                return _ids->at(0).size();
            }
            return _keys->size();
        };
        
        size_t nb_active_keys() const {return _keys->size() - _excluded_keys.size();};

        bool empty() const {
            return _keys->size() - _excluded_keys.size() == 0;
        }
        void print() const{
            cout << endl;
            auto i = 0;
            for(auto &key:*_keys){
                if (_excluded_keys.count(i++)==0) {
                    cout << key << " ";
                }
            }
            cout << endl;
        }
        
        string first() const{
            return _keys->front();
        }
        
        string last() const{
            return _keys->back();
        }
    };
        
    
    class node_pairs{
        
    public:
        string _name;
        std::vector<gravity::index_pair*> _keys;
        node_pairs(){
            _keys.resize(0);
        };
        node_pairs(string name):_name(name){_keys.resize(0);};
        ~node_pairs(){
            clear();
        }
        void clear() {
            for (auto p: _keys) { delete p;}
            _keys.clear();
        }
    };
    
    typedef enum { ipopt, gurobi, bonmin, cplex, sdpa, Mosek} SolverType;  /* Solver type */

    // settings of solvers. used by solvers like sdpa.
    typedef enum {unsolved = -1, penalty=0, fast=1, medium=2, stable=3} SolverSettings;

    typedef pair<shared_ptr<size_t>,shared_ptr<indices>> unique_id; /* A unique identifier is defined as a pair<variable id, index address> */

    template <class T>
    std::string type_name(const T& t) {
        return t._type_name;
    }
}

#endif
