//
//  param.h
//
//
//  Created by Hassan on 13/05/2016.
//
//

#ifndef ____param__
#define ____param__

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

using namespace std;

namespace gravity {
/** Backbone class for parameter */
class param_: public constant_ {

protected:
    string                                 _name;
    int                                    _id = -1;
    int                                    _vec_id = -1; /**< index in the vector array (useful for Cplex). **/
    NType                                  _intype;
    shared_ptr<map<string,unsigned>>       _indices = nullptr; /*<< A map storing all the indices this parameter has, the key is represented by a string, while the entry indicates the right position in the values and bounds
                                       vectors */
    unique_ptr<vector<unsigned>>           _ids = nullptr; /*<<A vector storing all the indices this parameter has in the order they were created */

    /* (Guanglei) added this part to record the indices of sdp variables. SDP should be indexed by a pair of integers. This is true for all SDP solvers. */
    shared_ptr<map<string,pair<unsigned, unsigned>>> _sdpindices;

public:

    unique_id                              _unique_id = make_tuple<>(-1,scalar_,0,0,0); /* */

    bool                                   _is_indexed = false;

    virtual ~param_() {};

    void set_id(size_t idx) {
        _id = idx;
        get<0>(_unique_id) = idx;
    };

    void set_vec_id(size_t idx) {
        _vec_id = idx;
    };

    size_t get_id() const {
//            if (_is_indexed) {
//                return _id + _indices->begin()->second;
//            }
        return _id;
    };

    size_t get_vec_id() const {
        return _vec_id;
    };

    size_t get_id_inst(unsigned inst = 0) const {
        if (_is_indexed && inst < _ids->size()) {
            return _ids->at(inst);
        }
//            throw invalid_argument("This is a non-indexed variable!\n");
        return inst;
    };

    // newly added part by guanglei
    pair<size_t,size_t> get_sdpid() const {
        if (_is_indexed) {
            return _sdpindices->begin()->second;
        }
        return make_pair(0, 0);
    };

    string get_name(bool indices=true) const;
    void set_name(const string s) {
        _name = s;
    };
    NType get_intype() const {
        return _intype;
    }
    size_t get_dim() const {
//        if (_is_indexed) {
//            return _ids->size();
//        }
//        return _dim;
        return max(_dim, _indices->size());
    }

    size_t get_nb_instances() const {
        return _ids->size();
    }

    shared_ptr<map<string,unsigned>> get_indices() const {
        return _indices;
    }


    /**  (guanglei) added sdpvar */
    shared_ptr<map<string, pair<unsigned, unsigned>>> get_sdpindices() const {
        return _sdpindices;
    }

    vector<unsigned>& get_ids() const {
        return *_ids;
    }

    void set_type(NType type) {
        _intype = type;
    }

    /** Querries */


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

    Sign get_all_sign() const; /**< If all instances of the current parameter/variable have the same sign, it returns it, otherwise, it returns unknown. **/
    Sign get_sign(int idx = 0) const; /**< returns the sign of one instance of the current parameter/variable. **/
    pair<constant_*, constant_*>* get_range() const;

    /** Operators */
    bool operator==(const param_& p) const {
        return (_id==p._id && *_ids==*p._ids);
        //        return (_id==p._id && _type==p._type && _intype==p._intype && get_name()==p.get_name());
    }
};


/** A parameter can be a bool, a short, an int, a float or a double*/
template<typename type = double>
class param: public param_ {
protected:
    shared_ptr<vector<type>>                _val;

public:
    shared_ptr<pair<type,type>>                         _range; /**< (Min,Max) values in vals **/

    param() {
        _type = par_c;
        _name = "noname";
        //    throw invalid_argument("Please enter a name in the parameter constructor");
        update_type();
        _val = make_shared<vector<type>>();
        _indices = make_shared<map<string,unsigned>>();
        _ids = unique_ptr<vector<unsigned>>(new vector<unsigned>());
        _sdpindices = make_shared<map<string,pair<unsigned, unsigned>>>();
        _range = make_shared<pair<type,type>>(make_pair<>(numeric_limits<type>::max(), numeric_limits<type>::lowest()));
    }

    ~param() {
    }


    param (const param& p) {
        _type = par_c;
        _intype = p._intype;
        _id = p._id;
        _unique_id = p._unique_id;
        _vec_id = p._vec_id;
        _val = p._val;
        _name = p._name;
        _indices = p._indices;
        _ids = unique_ptr<vector<unsigned>>(new vector<unsigned>(*p._ids));
        _sdpindices = p._sdpindices;
        _range = p._range;
        _is_transposed = p._is_transposed;
        _is_vector = p._is_vector;
        _is_indexed = p._is_indexed;
        _dim = p._dim;
    }

    param (param&& p) {
        _type = par_c;
        _intype = p._intype;
        _id = p._id;
        _unique_id = p._unique_id;
        _vec_id = p._vec_id;
        _val = p._val;
        _name = p._name;
        _indices = p._indices;
        _ids = move(p._ids);
        _sdpindices = p._sdpindices;
        _range = p._range;
        _is_transposed = p._is_transposed;
        _is_vector = p._is_vector;
        _is_indexed = p._is_indexed;
        _dim = p._dim;
    }

    param& operator=(const param& p) {
        _type = par_c;
        _intype = p._intype;
        _id = p._id;
        _unique_id = p._unique_id;
        _vec_id = p._vec_id;
        _val = p._val;
        _name = p._name;
        _indices = p._indices;
        _ids = unique_ptr<vector<unsigned>>(new vector<unsigned>(*p._ids));
        _sdpindices = p._sdpindices;
        _range = p._range;
        _is_transposed = p._is_transposed;
        _is_vector = p._is_vector;
        _is_indexed = p._is_indexed;
        _dim = p._dim;
        return *this;
    }

    param& operator=(param&& p) {
        _type = par_c;
        _intype = p._intype;
        _id = p._id;
        _unique_id = p._unique_id;
        _vec_id = p._vec_id;
        _val = p._val;
        _name = p._name;
        _indices = p._indices;
        _ids = move(p._ids);
        _sdpindices = p._sdpindices;
        _range = p._range;
        _is_transposed = p._is_transposed;
        _is_vector = p._is_vector;
        _is_indexed = p._is_indexed;
        _dim = p._dim;
        return *this;
    }

    param tr() const {
        auto p = param(*this);
        p._is_transposed = true;
        return p;
    }

    shared_ptr<vector<type>> get_vals() const {
        return _val;
    }

    void set_type(CType t) {
        _type = t;
    }

    void set_intype(NType t) {
        _intype = t;
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
        throw bad_alloc();
    }


    param(const string& s) {
        _name = s;
        update_type();
        _val = make_shared<vector<type>>();
        _indices = make_shared<map<string,unsigned>>();
        _ids = unique_ptr<vector<unsigned>>(new vector<unsigned>());
        _sdpindices = make_shared<map<string,pair<unsigned, unsigned>>>();
        _range = make_shared<pair<type,type>>(make_pair<>(numeric_limits<type>::max(), numeric_limits<type>::lowest()));
    }

    NType get_intype() const {
        return _intype;
    }

    type eval() const {
        if (_is_indexed) {
            return _val->at(_ids->back());
        }
        return _val->back();
    }

    type eval(int i) const {
        if (_is_indexed) {
            if (i >= _ids->size()) {
                return _val->at(_ids->at(0));
            }
            return _val->at(_ids->at(i));
        }
        return _val->at(i);
    }


    /* Modifiers */
    void   set_size(size_t s, type val = 0) {
        _val->resize(s,val);
        _dim = s;        
    };

    void add_val(type val) {
        _val->push_back(val);
        update_range(val);
        _dim++;
    }

    void update_range(type val) {
        if (val < _range->first) {
            _range->first = val;
        }
        if (val > _range->second) {
            _range->second = val;
        }
    }

    void set_val(unsigned i, type val) {
        if (_val->size() <= i) {
            throw out_of_range("set_val(int i, type val)");
        }
        _val->at(i) = val;
        update_range(val);
    }

    void set_val(type val) {
        for (auto &v: _val) {
            v = val;
        }
        _range->first = val;
        _range->second = val;
    }

    Sign get_sign(int idx = 0) const {
        assert(idx < _val->size());
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



    Sign get_all_sign() const {
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

    bool is_unit() const { /**< Returns true if all values of this paramter are 1 **/
        return (_range->first == 1 && _range->second == 1);
    }

    bool is_zero() const { /**< Returns true if all values of this paramter are 0 **/
        return (_range->first == 0 && _range->second == 0);
    }

    bool is_non_positive() const { /**< Returns true if all values of this paramter are <= 0 **/
        return (_range->second <= 0   && _range->first <= 0);
    }

    bool is_positive() const { /**< Returns true if all values of this paramter are positive **/
        return (_range->first > 0 && _range->second > 0);
    }

    bool is_non_negative() const { /**< Returns true if all values of this paramter are >= 0 **/
        return (_range->first >= 0  && _range->second >= 0);
    }

    bool is_negative() const { /**< Returns true if all values of this paramter are positive **/
        return (_range->second < 0  && _range->first < 0);
    }

    /** Operators */
    bool operator==(const param& p) const {
        return (get_name()==p.get_name() && _type==p._type && _intype==p._intype && _dim==p._dim && _indices==p._indices && _sdpindices==p._sdpindices && _val==p._val);
        //return (get_name()==p.get_name() && _type==p._type && _intype==p._intype && _dim==p._dim && _indices==p._indices && _val==p._val);
    }

    param& operator^(size_t d) {
        set_size(d);
        return *this;
    }

    void initialize_all(type v) {
        for (int i = 0; i<_val->size(); i++) {
            _val->at(i) = v;
        }
    }

    void initialize(size_t i, type v) {
        set_val(i,v);
    }

    param& operator=(type v) {
        if (_is_indexed) {
            _val->at(_ids->back()) = v;
        }
        else {
            _val->push_back(v);
            _dim++;
        }
        update_range(v);
        return *this;
    }

    template<typename... Args>
    param operator()(size_t t1, Args&&... args) {
        param res(this->_name);
        res._id = this->_id;
        res._vec_id = this->_vec_id;
        res._intype = this->_intype;
        res._range = this->_range;
        res._val = this->_val;
        list<size_t> indices;
        indices = {forward<size_t>(args)...};
        indices.push_front(t1);
        string key;
        auto it = indices.begin();
        for (size_t i= 0; i< indices.size(); i++) {
            key += to_string(*it);
            if (i<indices.size()-1) {
                key += ",";
            }
            it++;
        }
        auto pp = param_::_indices->insert(make_pair<>(key,param_::_indices->size()));
        _val->resize(max(_val->size(),param_::_indices->size()));
        if(pp.second) { //new index inserted
            if(res._indices->insert(make_pair<>(key,param_::_indices->size()-1)).second) {
                res._dim++;
                //_dim++;
            };
            res._ids->push_back(param_::_indices->size()-1);
        }
        else {
            if(res._indices->insert(make_pair<>(key,pp.first->second)).second) {
                res._dim++;
            }
            res._ids->push_back(pp.first->second);
        }
        res._name += "["+key+"]";
        res._unique_id = make_tuple<>(res._id,scalar_,typeid(type).hash_code(), res._ids->at(0), res._ids->at(res._ids->size()-1));
        res._is_indexed = true;
        //_is_indexed = true; // Guanglei added this line.
        return res;
    }

    template<typename... Args>
    param operator()(string t1, Args&&... args) {
        param res(this->_name);
        res._id = this->_id;
        res._vec_id = this->_vec_id;
        res._intype = this->_intype;
        res._range = this->_range;
        res._val = this->_val;
        list<string> indices;
        //indices = {forward<size_t>(args)...};
        indices = {forward<Args>(args)...};
        indices.push_front(t1);
        string key;
        auto it = indices.begin();
        for (int i= 0; i < indices.size(); i++) {
            key += *it;
            if (i< indices.size()-1) {
                key += ",";
            }
            it++;
        }
        auto pp = param_::_indices->insert(make_pair<>(key,param_::_indices->size()));
        _val->resize(max(_val->size(),param_::_indices->size()));
        if(pp.second) { //new index inserted
            if(res._indices->insert(make_pair<>(key,param_::_indices->size()-1)).second) {
                res._dim++;
            }
            res._ids->push_back(param_::_indices->size()-1);
        }
        else {
            if(res._indices->insert(make_pair<>(key,pp.first->second)).second) {
                res._dim++;
            }
            res._ids->push_back(pp.first->second);
        }
        //res._dim++;
        res._name += "["+key+"]";
        res._unique_id = make_tuple<>(res._id,scalar_,typeid(type).hash_code(), res._ids->at(0), res._ids->at(res._ids->size()-1));
        res._is_indexed = true;
        return res;
    }



    template<typename Tobj>
    param in(const vector<Tobj*>& vec) {
        param res(this->_name);
        res._id = this->_id;
        res._vec_id = this->_vec_id;
        res._intype = this->_intype;
        res._range = this->_range;
        res._val = this->_val;
        DebugOff(_name << " = ");
        string key;
        for(auto it = vec.begin(); it!= vec.end(); it++) {
            if(!(*it)->_active) {
                continue;
            }
            key = (*it)->_name;
            DebugOff(key<< ", ");
            auto pp = param_::_indices->insert(make_pair<>(key, param_::_indices->size()));
            _val->resize(max(_val->size(),param_::_indices->size()));
            if(pp.second) { //new index inserted
                if(res._indices->insert(make_pair<>(key, param_::_indices->size() - 1)).second) {
                    res._dim++;
                }
                res._ids->push_back(param_::_indices->size() - 1);
            }
            else {
                if(res._indices->insert(make_pair<>(key,pp.first->second)).second) {
                    res._dim++;
                }
                res._ids->push_back(pp.first->second);
            }
        }
        DebugOff(endl);
        res._name += ".in_" + string(typeid(Tobj).name());
        res._unique_id = make_tuple<>(res._id,in_, typeid(Tobj).hash_code(), res.get_id_inst(0),res.get_id_inst(res.get_dim()));
        res._is_indexed = true;
        return res;
    }

    template<typename Tobj>
    param in_pairs(const vector<Tobj*>& vec) {
        param res(this->_name);
        res._id = this->_id;
        res._vec_id = this->_vec_id;
        res._intype = this->_intype;
        res._range = this->_range;
        res._val = this->_val;
        DebugOff(_name << " = ");
        string key;
        for(auto it = vec.begin(); it!= vec.end(); it++) {
            if(!(*it)->_active) {
                continue;
            }
            key = (*it)->_src->_name + "," + (*it)->_dest->_name;
            DebugOff(key<< ", ");
            auto pp = param_::_indices->insert(make_pair<>(key, param_::_indices->size()));
            _val->resize(max(_val->size(),param_::_indices->size()));
            if(pp.second) { //new index inserted
                if(res._indices->insert(make_pair<>(key, param_::_indices->size() - 1)).second) {
                    res._dim++;
                }
                res._ids->push_back(param_::_indices->size() - 1);
            }
            else {
                if(res._indices->insert(make_pair<>(key,pp.first->second)).second) {
                    res._dim++;
                }
                res._ids->push_back(pp.first->second);
            }
        }
        DebugOff(endl);
        res._name += ".in_" + string(typeid(Tobj).name());
        res._unique_id = make_tuple<>(res._id,in_,typeid(Tobj).hash_code(), res.get_id_inst(0),res.get_id_inst(res.get_dim()));
        res._is_indexed = true;
        return res;
    }

    template<typename Tobj>
    param in_pairs(const vector<Tobj*>& vec, unsigned T) {
        assert(T > 0);
        param res(this->_name);
        res._id = this->_id;
        res._vec_id = this->_vec_id;
        res._intype = this->_intype;
        res._range = this->_range;
        res._val = this->_val;
        DebugOff(_name << " = ");
        string key;
        for (unsigned t = 0; t < T; t++) {
            for(auto it = vec.begin(); it!= vec.end(); it++) {
                if(!(*it)->_active) {
                    continue;
                }
                key = (*it)->_src->_name + "," + (*it)->_dest->_name;
                key += ",";
                key += to_string(t);
                Debug(key<< ", " << endl);
                auto pp = param_::_indices->insert(make_pair<>(key, param_::_indices->size()));
                _val->resize(max(_val->size(),param_::_indices->size()));
                if(pp.second) { //new index inserted
                    if(res._indices->insert(make_pair<>(key, param_::_indices->size() - 1)).second) {
                        res._dim++;
                    }
                    res._ids->push_back(param_::_indices->size() - 1);
                }
                else {
                    if(res._indices->insert(make_pair<>(key,pp.first->second)).second) {
                        res._dim++;
                    }
                    res._ids->push_back(pp.first->second);
                }
            }
        }
        DebugOff(endl);
        res._name += ".in_" + string(typeid(Tobj).name()) + "_time_" + to_string(T);
        res._unique_id = make_tuple<>(res._id, in_time_, typeid(Tobj).hash_code(),
                res.get_id_inst(0), res.get_id_inst(res.get_dim()));
        res._is_indexed = true;
        return res;
    }

    template<typename Tobj>
    param from(const vector<Tobj*>& vec) {
        param res(this->_name);
        res._id = this->_id;
        res._vec_id = this->_vec_id;
        res._intype = this->_intype;
        res._range = this->_range;
        res._val = this->_val;
        DebugOff(_name << " = ");
        string key;
        for(auto it = vec.begin(); it!= vec.end(); it++) {
            if(!(*it)->_active) {
                continue;
            }
            key = (*it)->_src->_name;
            DebugOff(key<< ", ");
            auto pp = param_::_indices->insert(make_pair<>(key, param_::_indices->size()));
            _val->resize(max(_val->size(),param_::_indices->size()));
            if(pp.second) { //new index inserted
                if(res._indices->insert(make_pair<>(key, param_::_indices->size() - 1)).second) {
                    res._dim++;
                    // _dim++;
                }
                res._ids->push_back(param_::_indices->size() - 1);
            }
            else {
                if(res._indices->insert(make_pair<>(key,pp.first->second)).second) {
                    res._dim++;
                }
                res._ids->push_back(pp.first->second);
            }
        }
        DebugOff(endl);
        res._name += ".from_" + string(typeid(Tobj).name());
        res._unique_id = make_tuple<>(res._id,from_,typeid(Tobj).hash_code(), res.get_id_inst(0),res.get_id_inst(res.get_dim()));
        res._is_indexed = true;
        return res;
    }

    template<typename Tobj>
    param to(const vector<Tobj*>& vec) {
        param res(this->_name);
        res._id = this->_id;
        res._vec_id = this->_vec_id;
        res._intype = this->_intype;
        res._range = this->_range;
        res._val = this->_val;
        DebugOff(_name << " = ");
        string key;
        for(auto it = vec.begin(); it!= vec.end(); it++) {
            if(!(*it)->_active) {
                continue;
            }
            key = (*it)->_dest->_name;
            DebugOff(key<< ", ");
            auto pp = param_::_indices->insert(make_pair<>(key, param_::_indices->size()));
            _val->resize(max(_val->size(),param_::_indices->size()));
            if(pp.second) { //new index inserted
                if(res._indices->insert(make_pair<>(key, param_::_indices->size() - 1)).second) {
                    res._dim++;
                }
                res._ids->push_back(param_::_indices->size() - 1);
            }
            else {
                if(res._indices->insert(make_pair<>(key,pp.first->second)).second) {
                    res._dim++;
                }
                res._ids->push_back(pp.first->second);
            }
        }
        DebugOff(endl);
        res._name += ".to_" + string(typeid(Tobj).name());
        res._unique_id = make_tuple<>(res._id,to_,typeid(Tobj).hash_code(), res.get_id_inst(0),res.get_id_inst(res.get_dim()));
        res._is_indexed = true;
        return res;
    }

    param in(const string nm, unsigned T) {
        param res(this->_name);
        res._id = this->_id;
        res._vec_id = this->_vec_id;
        res._intype = this->_intype;
        res._range = this->_range;
        res._val = this->_val;
        for (unsigned t = 0; t < T; t++) {
            string key = nm;
            //if (t > 0) {
                key += ",";
                key += to_string(t);
           // }
            //cout << "key: " << key << endl;
            Debug("key: " << key << endl);
            auto pp = param_::_indices->insert(make_pair<>(key, param_::_indices->size()));
            _val->resize(max(_val->size(),param_::_indices->size()));
            if(pp.second) { //new index inserted
                if(res._indices->insert(make_pair<>(key, param_::_indices->size() - 1)).second) {
                    res._dim++;
                    //_dim++;
                }
                res._ids->push_back(param_::_indices->size() - 1);
            }
            else {
                if(res._indices->insert(make_pair<>(key,pp.first->second)).second) {
                    res._dim++;
                }
                res._ids->push_back(pp.first->second);
            }
        }
        res._name += ".in_" +  string(typeid(nm).name()) + "_time_" + to_string(T);
        res._unique_id = make_tuple<>(res._id, in_time_, typeid(nm).hash_code(), res.get_id_inst(0), res.get_id_inst(param_::get_dim()));
        res._is_indexed = true;
        return res;
    }



    template<typename Tobj>
    param in_at(const vector<Tobj*>& nodes, unsigned t) {
        param res(this->_name);
        res._id = this->_id;
        res._vec_id = this->_vec_id;
        res._intype = this->_intype;
        res._range = this->_range;
        res._val = this->_val;
        string key;
        for(auto it = nodes.begin(); it!= nodes.end(); it++) {
            if(!(*it)->_active) {
                continue;
            }
            key = (*it)->_name;
            key += ",";
            key += to_string(t);
            auto pp = param_::_indices->insert(make_pair<>(key, param_::_indices->size()));
            _val->resize(max(_val->size(),param_::_indices->size()));
            if(pp.second) { //new index inserted
                if(res._indices->insert(make_pair<>(key, param_::_indices->size() - 1)).second) {
                    res._dim++;
                }
                res._ids->push_back(param_::_indices->size() - 1);
            }
            else {
                if(res._indices->insert(make_pair<>(key,pp.first->second)).second) {
                    res._dim++;
                }
                res._ids->push_back(pp.first->second);
            }
        }
        res._name += ".in_" + string(typeid(Tobj).name()) + "_at_" + to_string(t);
        res._unique_id = make_tuple<>(res._id, in_at_, typeid(Tobj).hash_code(), res.get_id_inst(0), res.get_id_inst(param_::get_dim()));
        res._is_indexed = true;
        return res;
    }

    template<typename Tobj>
    param in(const vector<Tobj*>& nodes, unsigned T) {
        param res(this->_name);
        res._id = this->_id;
        res._vec_id = this->_vec_id;
        res._intype = this->_intype;
        res._range = this->_range;
        res._val = this->_val;
        string key;
        for (unsigned t = 0; t < T; t++) {
            for(auto it = nodes.begin(); it!= nodes.end(); it++) {
                if(!(*it)->_active) {
                    continue;
                }
                key = (*it)->_name;
                //if (t > 0) {
                    key += ",";
                    key += to_string(t);
                //}
                Debug("_val: " << _val->size() << endl);
                Debug("_indices: " << param_::_indices->size() << endl);
                auto pp = param_::_indices->insert(make_pair<>(key, param_::_indices->size()));
                _val->resize(max(_val->size(),param_::_indices->size()));
                if(pp.second) { //new index inserted
                    if(res._indices->insert(make_pair<>(key, param_::_indices->size() - 1)).second) {
                        res._dim++;
                        //_dim++;
                    }
                    res._ids->push_back(param_::_indices->size() - 1);
                }
                else {
                    if(res._indices->insert(make_pair<>(key,pp.first->second)).second) {
                        res._dim++;
                    }
                    res._ids->push_back(pp.first->second);
                }
            }
        }
        res._name += ".in_" + string(typeid(Tobj).name()) + "_time_" + to_string(T);
        res._unique_id = make_tuple<>(res._id,in_time_,typeid(Tobj).hash_code(), res.get_id_inst(0), res.get_id_inst(param_::get_dim()));
        res._is_indexed = true;
        return res;
    }

    template<typename Tobj>
    param from(const vector<Tobj*>& arcs, unsigned T) {
        param res(this->_name);
        res._id = this->_id;
        res._vec_id = this->_vec_id;
        res._intype = this->_intype;
        res._range = this->_range;
        res._val = this->_val;
        DebugOff(_name << " = ");
        string key;
        for (unsigned t = 0; t < T; t++) {
            for(auto it = arcs.begin(); it!= arcs.end(); it++) {
                if(!(*it)->_active){  // || !(*it)->_src->_active || !(*it)->_dest->_active ) {
                    continue;
                }
                key = (*it)->_src->_name;
                //if (t > 0) {
                    key += ",";
                    key += to_string(t);
                //}
                auto pp = param_::_indices->insert(make_pair<>(key,param_::_indices->size()));
                _val->resize(max(_val->size(),param_::_indices->size()));
                if(pp.second) { //new index inserted
                    if(res._indices->insert(make_pair<>(key,param_::_indices->size()-1)).second) {
                        res._dim++;
                        // _dim++;
                    }
                    res._ids->push_back(param_::_indices->size()-1);
                }
                else {
                    if(res._indices->insert(make_pair<>(key,pp.first->second)).second) {
                        res._dim++;
                    }
                    res._ids->push_back(pp.first->second);
                }
            }
        }
        res._name += ".from_" + string(typeid(Tobj).name()) + "_time_" + to_string(T);
        res._unique_id = make_tuple<>(res._id,from_time_,typeid(Tobj).hash_code(), res.get_id_inst(0), res.get_id_inst(param_::get_dim()));
        res._is_indexed = true;
        return res;
    }

    template<typename Tobj>
    param to(const vector<Tobj*>& arcs, unsigned T) {
        param res(this->_name);
        res._id = this->_id;
        res._vec_id = this->_vec_id;
        res._intype = this->_intype;
        res._range = this->_range;
        res._val = this->_val;
        string key;
        for (unsigned t = 0; t < T; t++) {
            for(auto it = arcs.begin(); it!= arcs.end(); it++) {
                if(!(*it)->_active){ // || !(*it)->_src->_active || !(*it)->_dest->_active ) {
                    continue;
                }
                key = (*it)->_dest->_name;
                //if (t > 0) {
                    key += ",";
                    key += to_string(t);
                //}
                auto pp = param_::_indices->insert(make_pair<>(key,param_::_indices->size()));
                _val->resize(max(_val->size(),param_::_indices->size()));
                if(pp.second) { //new index inserted
                    if(res._indices->insert(make_pair<>(key,param_::_indices->size()-1)).second) {
                        res._dim++;
                        //_dim++;
                    }
                    res._ids->push_back(param_::_indices->size() - 1);
                }
                else {
                    // already exists
                    if(res._indices->insert(make_pair<>(key,pp.first->second)).second) {
                        res._dim++;
                    }
                    res._ids->push_back(pp.first->second);
                }

            }
        }
        res._name += ".to_"  + string(typeid(Tobj).name()) + "_time_"+ to_string(T);
        res._unique_id = make_tuple<>(res._id, to_time_, typeid(Tobj).hash_code(), res.get_id_inst(0), res.get_id_inst(param_::get_dim()));
        res._is_indexed = true;
        return res;
    }

    // T copies of a parameter
    void time_expand(unsigned T) {
        assert(T >= 1);
        auto dim = param_::get_dim();
        set_size(dim*T);
        /* update the indices of the old parameter*/
        string key;
        auto map_temp = *param::_indices;
//        for (map<std::string, unsigned>::iterator it= param::_indices->begin(); it != param::_indices->end(); it++){
//            key = it->first;
//            map_temp.insert(make_pair(key, it->second));
//        }
        //CLEAR OLD ENTRIES
        _indices->clear();
        _ids->clear();
        //STORE NEW ENTRIES
        for(unsigned t = 0; t < T; t ++ ) {
            for (auto &entry: map_temp){
                key = entry.first;
                key += ",";
                key += to_string(t);
                _val->at(param_::_indices->size()) = _val->at(entry.second);
                param_::_indices->insert(make_pair<>(key, param_::_indices->size()));
            }
        }
        _name += ".time_expanded";
        _unique_id = make_tuple<>(_id,in_time_,typeid(type).hash_code(),0,dim*T);
    }

    /** Output */
    void print(bool vals=false) const {
        cout << this->to_str(vals);
    }

    string to_str(bool vals=false) const {
        string str = get_name();
        if(vals) {
            str += " = [ ";
            for(int i = 0 ; i < param_::get_dim(); i++) {
                str += std::to_string(_val->at(i));
                str += " ";
            }
            str += "];";
        }
        return str;
    }

    type getvalue() const {
        if (_is_indexed) {
            return (_val->at(_indices->begin()->second));
        }
        else {
            return _val->at(0);
        }
    }
};

template<typename type>
param<type> all(const param<type>& p) {
    auto pp = param<type>(p);
    pp._is_vector = true;
    return pp;
}
}
#endif /* defined(____param__) */
