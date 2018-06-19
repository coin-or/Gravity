//
//  param.h
//
//
//  Created by Hassan on 13/05/2016.
//
//

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

using namespace std;


template <typename T>
std::string to_string_with_precision(const T a_value, const int n)
{
    std::ostringstream out;
    out << std::setprecision(n) << a_value;
    return out.str();
}


namespace gravity {
    /** Backbone class for parameter */
    class param_: public constant_ {
        
    protected:
        
        NType                                  _intype;
        shared_ptr<map<string,unsigned>>       _indices = nullptr; /*<< A map storing all the indices this parameter has, the key is represented by a string, while the entry indicates the right position in the values and bounds vectors */
        shared_ptr<vector<string>>              _rev_indices = nullptr; /*<< A vector storing all the indices this parameter has */
        
        
        
        /* (Guanglei) added this part to record the indices of sdp variables. SDP should be indexed by a pair of integers. This is true for all SDP solvers. */
        shared_ptr<map<string,pair<unsigned, unsigned>>> _sdpindices;
        
    public:
        
        shared_ptr<int>                                _id = make_shared<int>(-1);
        shared_ptr<int>                                _vec_id; /**< index in the vector array (useful for Cplex). **/
        
        string                                         _name;
        shared_ptr<vector<vector<unsigned>>>           _ids = nullptr; /*<<A vector storing all the indices this parameter has in the order they were created */
        unique_id                              _unique_id = make_tuple<>(_id,unindexed_,0,0,0); /* */
        
        bool                                   _is_indexed = false;
        bool                                   _is_relaxed = false;
        bool                                   _time_extended = false;
        bool                                   _new = true; /**< Will become false once this param is added to a program. Can be useful for iterative model solving. */
        vector<double>                         _l_dual; /*<<Dual values for lower bounds */
        vector<double>                         _u_dual; /*<<Dual values for upper bounds */
        
        
        virtual ~param_() {};
        
        void set_id(size_t idx) {
            *_id = idx;
            get<0>(_unique_id) = _id;
        };
        
        void set_vec_id(size_t idx) {
            *_vec_id = idx;
        };
        
        int get_id() const {
            return *_id;
        };
        
        size_t get_vec_id() const {
            return *_vec_id;
        };
        
        size_t get_id_inst(unsigned inst = 0) const {
            if (_is_indexed) {
                if(_ids->at(0).size() <= inst){
                    throw invalid_argument("get_id_inst out of range");
                }
                return _ids->at(0).at(inst);
            }
            return inst;
        };
        
        size_t get_id_inst(unsigned inst1, unsigned inst2) const {
            if (_is_indexed) {
                if (_ids->size()==1) {
                    return _ids->at(0).at(inst2);
                }
                return _ids->at(inst1).at(inst2);
            }
            return inst2;
            //        throw invalid_argument("Calling get_id_inst on a non indexed parpam\n");
        };
        
        void copy(const param_& p) {
            _id = p._id;
            _unique_id = p._unique_id;
            _vec_id = p._vec_id;
            _name = p._name;
            _indices = p.get_indices();
            _rev_indices = p.get_rev_indices();
            _ids =p._ids;
            _is_transposed = p._is_transposed;
            _is_vector = p._is_vector;
            _is_matrix = p._is_matrix;
            _is_indexed = p._is_indexed;
            _dim = p._dim;
            _time_extended = p._time_extended;
        }

        pair<size_t,size_t> get_sdp_inst(unsigned inst = 0) const {
            int idx = _ids->at(0).at(inst);
            auto it = _indices->begin();
            for(it =_indices->begin(); it != _indices->end(); ++it) {
                if(it->second==idx) break;
            }
            string key = it->first;
            DebugOff("\nkey found: " << it->first.at(0) << "," << it->first.at(2));
            pair<size_t,size_t> res;
            res = make_pair(std::stoi(it->first.substr(0,1)),std::stoi(it->first.substr(2,1)));
            return res;
        };
        
        // newly added part by guanglei
        pair<size_t,size_t> get_sdpid() const {
            if (_is_indexed) {
                return _sdpindices->begin()->second;
            }
            return make_pair(0, 0);
        };
        
        string get_name(bool indices=false) const;
        void set_name(const string s) {
            _name = s;
        };
        
        string get_name(size_t inst) const {
            string name = _name;
            name = name.substr(0, name.find(".", 0));
            if (_is_indexed && name.find("[")!=std::string::npos) {// Name has index already
                return name;
            }
            if (_is_indexed) {
                int rev_idx = _ids->at(0).at(inst);
                name += "["+_rev_indices->at(rev_idx)+"]";
            }
            else {
                name += "["+_rev_indices->at(inst)+"]";
            }
            return name;
        };
        
        string get_name(size_t inst1, size_t inst2) const {
            string name = _name;
            name = name.substr(0, name.find(".", 0));
            if (_is_indexed && name.find("[")!=std::string::npos) {// Name has index already
                return name;
            }
            if (_is_indexed) {
                //                if (<#condition#>) {
                //                    <#statements#>
                //                }
                name += "["+_rev_indices->at(_ids->at(inst1).at(inst2))+"]";
            }
            else {
                name += "["+_rev_indices->at(inst2)+"]";
            }
            return name;
        };
        
        NType get_intype() const {
            return _intype;
        }
        
        
        
        
        
        shared_ptr<map<string,unsigned>> get_indices() const {
            return _indices;
        }
        
        shared_ptr<vector<string>> get_rev_indices() const {
            return _rev_indices;
        }
        
        
        /**  (guanglei) added sdpvar */
        shared_ptr<map<string, pair<unsigned, unsigned>>> get_sdpindices() const {
            return _sdpindices;
        }
        
        shared_ptr<vector<vector<unsigned>>> get_ids() const {
            return _ids;
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
        pair<Real,Real>* get_range() const;
        
        /** Operators */
        bool operator==(const param_& p) const {
            return (_unique_id == p._unique_id);
            //        return (_id==p._id && _type==p._type && _intype==p._intype && get_name()==p.get_name());
        }
        size_t get_nb_instances() const {
            if (_is_indexed) {
                if (_ids->size()>1) {
                    throw invalid_argument("get_nb_instances should be called with index here\n");
                }
                return _ids->at(0).size();
            }
            return get_dim();
        }
        
        size_t get_nb_instances(unsigned inst) const {
            if (_is_indexed) {
                return _ids->at(inst).size();
            }
            return get_dim();
        }
    };
    
    
    /** A parameter can be a bool, a short, an int, a float or a double*/
    template<typename type = double>
    class param: public param_ {
    protected:
        
        
    public:
        shared_ptr<vector<type>>                _val;
        shared_ptr<pair<type,type>>             _range; /**< (Min,Max) values in vals **/
        
        param() {
            _type = par_c;
            _name = "noname";
            _vec_id = make_shared<int>(-1);
            //    throw invalid_argument("Please enter a name in the parameter constructor");
            update_type();
            _val = make_shared<vector<type>>();
            _dim.resize(1,0);
            _indices = make_shared<map<string,unsigned>>();
            _rev_indices = make_shared<vector<string>>();
            _ids = make_shared<vector<vector<unsigned>>>();
            _ids->resize(1);
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
            _rev_indices = p._rev_indices;
            _ids = make_shared<vector<vector<unsigned>>>(*p._ids);
            _sdpindices = p._sdpindices;
            _range = p._range;
            _is_transposed = p._is_transposed;
            _is_vector = p._is_vector;
            _is_matrix = p._is_matrix;
            _is_indexed = p._is_indexed;
            _dim = p._dim;
            _time_extended = p._time_extended;
        }
        
        param (param&& p) {
            _type = par_c;
            _intype = p._intype;
            _id = p._id;
            _unique_id = p._unique_id;
            _vec_id = p._vec_id;
            _val = move(p._val);
            _name = p._name;
            _indices = move(p._indices);
            _rev_indices = move(p._rev_indices);
            _ids = move(p._ids);
            _sdpindices = p._sdpindices;
            _range = p._range;
            _is_transposed = p._is_transposed;
            _is_vector = p._is_vector;
            _is_matrix = p._is_matrix;
            _is_indexed = p._is_indexed;
            _dim = p._dim;
            _time_extended = p._time_extended;
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
            _rev_indices = p._rev_indices;
            _ids = make_shared<vector<vector<unsigned>>>(*p._ids);
            _sdpindices = p._sdpindices;
            _range = p._range;
            _is_transposed = p._is_transposed;
            _is_vector = p._is_vector;
            _is_matrix = p._is_matrix;
            _is_indexed = p._is_indexed;
            _dim = p._dim;
            _time_extended = p._time_extended;
            return *this;
        }
        
        param& operator=(param&& p) {
            _type = par_c;
            _intype = p._intype;
            _id = p._id;
            _unique_id = p._unique_id;
            _vec_id = p._vec_id;
            _val = move(p._val);
            _name = p._name;
            _indices = move(p._indices);
            _rev_indices = move(p._rev_indices);
            _ids = move(p._ids);
            _sdpindices = p._sdpindices;
            _range = p._range;
            _is_transposed = p._is_transposed;
            _is_vector = p._is_vector;
            _is_matrix = p._is_matrix;
            _is_indexed = p._is_indexed;
            _dim = p._dim;
            _time_extended = p._time_extended;
            return *this;
        }
        
        void transpose(){
            _is_transposed = !_is_transposed;
            if (!_is_vector) {
                _is_vector = true;
            }
            if (_is_matrix) {
                //            auto new_val = make_shared<vector<type>>();
                //            new_val->resize(get_dim());
                //            for (int i = 0; i<_dim[1]; i++) {
                //                for (int j = 0; j<_dim[0]; j++) {
                //                    new_val->at(i*_dim[0]+j) = eval(j, i);
                //                }
                //            }
                auto temp = _dim[0];
                _dim[0] = _dim[1];
                _dim[1] = temp;
                
                //            _val = new_val;
            }
        }
        
        param tr() const {
            auto p = param(*this);
            p.transpose();
            //        p._is_transposed = true;
            //        p._is_vector = true;
            return p;
        }
        
        param vec() const {
            auto p = param(*this);
            p._is_vector = true;
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
            
            _vec_id = make_shared<int>(-1);
            update_type();
            _val = make_shared<vector<type>>();
            _dim.resize(1,0);
            _indices = make_shared<map<string,unsigned>>();
            _rev_indices = make_shared<vector<string>>();
            _ids = make_shared<vector<vector<unsigned>>>();
            _ids->resize(1);
            _sdpindices = make_shared<map<string,pair<unsigned, unsigned>>>();
            _range = make_shared<pair<type,type>>(make_pair<>(numeric_limits<type>::max(), numeric_limits<type>::lowest()));
        }
        
        NType get_intype() const {
            return _intype;
        }
        
        type eval() const {
            if (_is_indexed) {
                return _val->at(_ids->at(0).back());
            }
            return _val->back();
        }
        
        type eval(unsigned i) const {
            if (_is_indexed) {
//                if (_ids->size()>1) {
//                    throw invalid_argument("eval() should be called with double index here\n");
//                }
                if (_val->size()<=_ids->at(0).at(i)){
                    throw invalid_argument("Param eval out of range");
                }
                return _val->at(_ids->at(0).at(i));
            }
            if (_val->size()<=i){
                throw invalid_argument("Param eval out of range");
            }
            return _val->at(i);
        }
        
        
        type eval(const string& key) const{
            return _val->at(param_::_indices->at(key));
        }
        
        type eval(unsigned i, unsigned j) const {
            
            if (_is_indexed && _ids->size()>1) {
//                if (_ids->at(i).at(j) >= _val->size()) {
//                    throw invalid_argument("eval(i,j): out of range");
//                }
                return _val->at(_ids->at(i).at(j));
            }
//
//            
            if (!_is_matrix) {
                return eval(j);
            }
            //        if (_is_indexed) {
            //            if (i >= _ids->size()) {
            ////                if (i>=_val->size()) {
            //                    throw invalid_argument("error");
            ////                }
            //                return _val->at(_ids->at(0));
            //            }
            //            if (_ids->at(i)>=_val->size()) {
            //                throw invalid_argument("error");
            //            }
            if (_is_transposed) {
                return _val->at(j*_dim[0]+i);//TODO same for vars
            }
            return _val->at(i*_dim[1]+j);
            //        }
            //        if (i>=_val->size()) {
            //            throw invalid_argument("error");
            //        }
            //        throw invalid_argument("cannot ")
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
                _dim = dims;
                //TODO allocate val
            }
        }
        
        void   set_size(size_t s1, size_t s2, type val = 0) {
            _is_matrix = true;
            _dim.resize(2);
            _dim[0] = s1;
            _dim[1] = s2;
            auto index = _dim[1]*s1-1+s2-1;
            _val->resize(index+1);
            _val->at(index) = val;
        };
        
        void   set_size(size_t s, type val = 0) {
            _val->resize(s, val);
            _dim[0] = s;
        };
        
        
        void add_val(type val) {
            _val->push_back(val);
            update_range(val);
            _dim[0]++;
        }
        
        void update_range(type val) {
            if (val < _range->first) {
                _range->first = val;
            }
            if (val > _range->second) {
                _range->second = val;
            }
        }
        
        void set_val(size_t i, size_t j, type val) {
            _is_matrix = true;
            _dim.resize(2);
            _dim[0] = max(_dim[0],i+1);
            _dim[1] = max(_dim[1],j+1);
            auto index = _dim[1]*i+j;
            _val->resize(max(_val->size(),index+1));
            _val->at(index) = val;
            update_range(val);
        }
        
        size_t set_val(const string& key, type val) {
            auto index = param_::_indices->size();
            auto pp = param_::_indices->insert(make_pair<>(key,index));
            if (pp.second) {//new index inserted
                _val->resize(max(_val->size(),index+1));
                _dim[0] = max(_dim[0],_val->size());
                _rev_indices->resize(_val->size());
                _rev_indices->at(index) = key;
                _val->at(index) = val;
                update_range(val);
                return index;
            }
            else {
                _val->at(pp.first->second) = val;
                update_range(val);
                return pp.first->second;
            }
        }
        
        void set_val(size_t i, type val) {
            _dim[0] = max(_dim[0],i+1);
            _val->resize(max(_val->size(),i+1));
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
        
//        param& operator^(size_t d) {
//            set_size(d);
//            return *this;
//        }
        
        void initialize_normal(double mean, double dev) {
            std::default_random_engine generator;
            std::normal_distribution<double> distribution(mean,dev);
            for (int i = 0; i<_val->size(); i++) {
                _val->at(i) = distribution(generator);
            }
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
                _val->at(_ids->at(0).back()) = v;
            }
            else {
                _val->push_back(v);
                _dim[0]++;
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
            res._is_vector = this->_is_vector;
            res._is_matrix = this->_is_matrix;
            res._is_transposed = _is_transposed;
            res._rev_indices = this->_rev_indices; res._indices = this->_indices;
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
            size_t index = 0;
            if (indices.size()==2) {
                _is_matrix = true;
                _dim.resize(2);
                _dim[0] = max(_dim[0],indices.front()+1);
                _dim[1] = max(_dim[1],indices.back()+1);
                index = _dim[1]*indices.front()+indices.back();
            }
            else {
                _dim[0] = max(_dim[0],indices.front()+1);
                index = indices.front();
            }
            auto pp = param_::_indices->insert(make_pair<>(key,index));
            if(pp.second) { //new index inserted
                _val->resize(max(_val->size(),index+1));
                _rev_indices->resize(_val->size());
                _rev_indices->at(index) = key;
                res._ids->at(0).push_back(index);
            }
            else {
                res._ids->at(0).push_back(pp.first->second);
            }
            res._dim[0]=1;
            res._name += "["+key+"]";
            res._unique_id = make_tuple<>(res._id,unindexed_,typeid(type).hash_code(), indices.front(), indices.back());
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
            res._is_vector = this->_is_vector;
            res._is_matrix = this->_is_matrix;
            res._is_transposed = _is_transposed;
            res._rev_indices = this->_rev_indices; res._indices = this->_indices;
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
            if (indices.size()==2) {
                _is_matrix = true;
            }
            auto index = param_::_indices->size();
            auto pp = param_::_indices->insert(make_pair<>(key,index));
            if(pp.second) { //new index inserted
                _val->resize(max(_val->size(),index+1));
                _dim[0] = max(_dim[0],_val->size());
                _rev_indices->resize(_val->size());
                _rev_indices->at(index) = key;
                res._ids->at(0).push_back(param_::_indices->size()-1);
            }
            else {
                res._ids->at(0).push_back(pp.first->second);
            }
            res._dim[0]=1;
            res._name += "["+key+"]";
            res._unique_id = make_tuple<>(res._id,unindexed_,typeid(type).hash_code(), res._ids->at(0).at(0), res._ids->at(0).at(res._ids->at(0).size()-1));
            res._is_indexed = true;
            return res;
        }
        
        
        param in(const node_pairs& np){
            auto p = this->in(np._keys);
            p._name += np._name;
            return p;
        }
                
//        template<typename Tobj> param in(const vector<Tobj>& vec) {
//            return in(get_ptr_vec(vec));
//        }
//        
//        template<typename Tobj> param from(const vector<Tobj>& vec) {
//            return from(get_ptr_vec(vec));
//        }
//        
//        template<typename Tobj> param to(const vector<Tobj>& vec) {
//            return to(get_ptr_vec(vec));
//        }
        
        template<typename Tobj> param in(const vector<Tobj>& vec) {
            param res(this->_name);
            res._id = this->_id;
            res._vec_id = this->_vec_id;
            res._intype = this->_intype;
            res._range = this->_range;
            res._val = this->_val;
            res._is_vector = this->_is_vector;
            res._is_matrix = this->_is_matrix;
            res._is_transposed = _is_transposed;
            res._rev_indices = this->_rev_indices; res._indices = this->_indices;
            if(vec.empty()){
                DebugOn("In function param.in(const vector<Tobj*>& vec), vec is empty!\n. Creating and empty variable! Check your sum/product operators.\n");
                res._name += "EMPTY_VAR";
                res._is_indexed = true;
                return res;
            }
            DebugOff(_name << " = ");
            string key, excluded;
            for(auto it = vec.begin(); it!= vec.end(); it++) {
                if(!(*it)._active) {
                    excluded += (*it)._name + ",";
                    continue;
                }
                key = (*it)._name;
                auto index = _indices->size();
                auto pp = param_::_indices->insert(make_pair<>(key, index));
                if(pp.second) { //new index inserted
                    _val->resize(max(_val->size(),index+1));
                    _dim[0] = max(_dim[0],_val->size());
                    _rev_indices->resize(_val->size());
                    _rev_indices->at(index) = key;
                    res._ids->at(0).push_back(index);
                }
                else {
                    res._ids->at(0).push_back(pp.first->second);
                }
            }
            DebugOff(endl);
            res._dim[0]=res._ids->at(0).size();
            res._name += ".in_" + vec.front()._type_name+ "," + excluded + "," + to_string(res._dim[0]);
            res._unique_id = make_tuple<>(res._id,in_, typeid(Tobj).hash_code(), 0,res._dim[0]);
            res._is_indexed = true;
            return res;
        }
        
        template<typename... Args>
        param prev(const indices& vec1, Args&&... args) {
            param res(this->_name);
            res._id = this->_id;
            res._vec_id = this->_vec_id;
            res._intype = this->_intype;
            res._range = this->_range;
            res._val = this->_val;
            res._is_vector = this->_is_vector;
            res._is_matrix = this->_is_matrix;
            res._is_transposed = _is_transposed;
            res._rev_indices = this->_rev_indices; res._indices = this->_indices;
            list<indices> vecs;
            vecs = {forward<Args>(args)...};
            vecs.push_front(vec1);
            if(vecs.size()==1 && res._indices->size()==0){
                res._indices = vec1._indices_map;
                res._rev_indices = vec1._indices;
                _indices = vec1._indices_map;
                _rev_indices = vec1._indices;
                auto dim = res._indices->size();
                for(unsigned idx = 0; idx < dim;idx++){
                    if(vec1._excluded_indices.count(idx)==0){
                        res._ids->at(0).push_back(idx-1);
                    }
                }
                dim =res._ids->at(0).size();
                _val->resize(dim);
                _dim[0] = dim;
                res._val->resize(dim);
                res._dim[0]= dim;
                res._name += ".prev";
                res._unique_id = make_tuple<>(res._id,prev_, typeid(indices).hash_code(), 0,res._dim[0]);
                res._is_indexed = true;
                return res;
                
            }
            size_t dim = 1;
            vector<size_t> dims;
            for(auto &vec: vecs){
                dim *= vec.size();
                dims.push_back(vec.size());
            }
            unsigned den = 1;
            size_t real_idx = 0;
            for(size_t idx = 0; idx < dim ; idx++){
                bool excluded = false;
                string key;
                den = dim;
                for(auto it = vecs.begin(); it!= vecs.end(); it++) {
                    auto vec = &(*it);
                    den /= vec->size();
                    real_idx = (idx/den)%vec->size();
                    if(vec->_excluded_indices.count(real_idx)!=0){
                        excluded = true;
                        break;
                    }
                    key += vec->_indices->at(real_idx-1);
                    if(next(it)!=vecs.end()){
                        key += ",";
                    }
                }
                if(excluded){
                    continue;
                }
                auto index = _indices->size();
                auto pp = param_::_indices->insert(make_pair<>(key, index));
                if(pp.second) { //new index inserted
                    _val->resize(max(_val->size(),index+1));
                    _dim[0] = max(_dim[0],_val->size());
                    _rev_indices->resize(_val->size());
                    _rev_indices->at(index) = key;
                    res._ids->at(0).push_back(index);
                    Debug(key << "; ");
                }
                else {
                    res._ids->at(0).push_back(pp.first->second);
                }
            }
            Debug(endl);
            res._dim[0]=res._ids->at(0).size();
            res._name += ".prev";
            res._unique_id = make_tuple<>(res._id,prev_, typeid(indices).hash_code(), 0,res._dim[0]);
            res._is_indexed = true;
            return res;
        }
                
        template<typename Tobj> param min_time(const vector<Tobj*>& vec, const indices& ids, const param<int>& time){
            param res(this->_name);
            res._id = this->_id;
            res._vec_id = this->_vec_id;
            res._intype = this->_intype;
            res._range = this->_range;
            res._val = this->_val;
            res._rev_indices = this->_rev_indices; res._indices = this->_indices;
            res._is_vector = _is_vector;
            res._is_matrix = _is_matrix;
            res._is_transposed = _is_transposed;
            if(vec.empty()){
                DebugOn("In function param.out_arcs(const vector<Tobj*>& vec), vec is empty!\n. Creating and empty variable! Check your sum/product operators.\n");
                res._name += "EMPTY_VAR";
                res._is_indexed = true;
                return res;
            }
            DebugOff(_name << " = ");
            string key;
            unsigned inst = 0;
            for(auto it = vec.begin(); it!= vec.end(); it++) {
                if(!(*it)->_active) {
                    continue;
                }
                
                for (auto i=0; i< ids.size();i++) {
                    if (inst>0) {
                        res._ids->push_back(vector<unsigned>());
                        res._dim.push_back(0);
                    }
                    if(ids._excluded_indices.count(i)!=0){
                        continue;
                    }
                    for (auto j=i; j< min(i+time.eval((*it)->_name),(int)ids.size()); j++) {
                        if(ids._excluded_indices.count(j)!=0){
                            continue;
                        }
                        key = (*it)->_name +","+ids._indices->at(j);
                        auto index = _indices->size();
                        auto pp = param_::_indices->insert(make_pair<>(key, index));
                        if(pp.second) { //new index inserted
                            _val->resize(max(_val->size(),index+1));
                            _dim[0] = max(_dim[0],_val->size());
                            _rev_indices->resize(_val->size());
                            _rev_indices->at(index) = key;
                            res._ids->at(inst).push_back(index);
                        }
                        else {
                            res._ids->at(inst).push_back(pp.first->second);
                        }
                    }
                    res._dim[inst]=res._ids->at(inst).size();
                    ++inst;
                }
            }
            res._unique_id = make_tuple<>(res._id,min_time_, typeid(Tobj).hash_code(), 0,res._dim.size());
            res._is_indexed = true;
            return res;
        }
        
        int nthOccurrence(const std::string& str, const std::string& findMe, int nth)
        {
            size_t  pos = 0;
            int     cnt = 0;
            
            while( cnt != nth )
            {
                pos+=1;
                pos = str.find(findMe, pos);
                if ( pos == std::string::npos )
                    return -1;
                cnt++;
            }
            return pos;
        }

        
        template<typename... Args>
        param in(const indices& vec1, Args&&... args) {            
            auto ids = indices(vec1,args...);
            param res(this->_name);
            res._id = this->_id;
            res._vec_id = this->_vec_id;
            res._intype = this->_intype;
            res._range = this->_range;
            res._val = this->_val;
            res._is_vector = this->_is_vector;
            res._is_matrix = this->_is_matrix;
            res._is_transposed = _is_transposed;
            res._time_extended = _time_extended;
            if(_indices->empty()){
                res._time_extended = vec1._time_extended;
                _time_extended = vec1._time_extended;
            }
            res._rev_indices = this->_rev_indices; res._indices = this->_indices;
            if(ids._excluded_indices.empty() && res._time_extended==ids._time_extended && (res._indices->size()==ids.size() || res._indices->size()==0)){
                _indices = ids._indices_map;
                _rev_indices = ids._indices;
                auto dim = _indices->size();
                _val->resize(dim);
                _dim[0] = dim;
                res._indices = _indices;
                res._rev_indices = _rev_indices;
                res._dim = _dim;
                res._unique_id = _unique_id;                
                return res;
            }
            if(res._indices->size()==0){
                res._time_extended = ids._time_extended;
                _indices = ids._indices_map;
                _rev_indices = ids._indices;
                res._indices = _indices;
                res._rev_indices = _rev_indices;
                auto dim = _indices->size();
                _val->resize(dim);
                _dim[0] = dim;
                for(unsigned idx = 0; idx < dim;idx++){
                    if(ids._excluded_indices.count(idx)==0){
                        res._ids->at(0).push_back(idx);
                    }
                }
                dim = res._ids->at(0).size();
                res._dim[0]= dim;
                if(dim>0){
                    res._name += ".in_indices("+vec1._indices->at(0)+","+vec1._indices->at(vec1._indices->size()-1)+")";
                }
                else{
                    res._name += ".empty_index_set";
                }
                res._unique_id = make_tuple<>(res._id,in_, typeid(indices).hash_code(), 0,dim);
                res._is_indexed = true;
                return res;
            }
            size_t idx = 0;
            for(auto key: *ids._indices){
                if(ids._excluded_indices.count(idx++)!=0){
                    continue;
                }
                if(!res._time_extended && ids._time_extended){
                    auto pos = nthOccurrence(key, ",", ids._time_pos);
                    key = key.substr(0,pos);
                }
                if(res._indices->size()==1 && res._rev_indices->at(0).compare(key)!=0){
                    key = res._rev_indices->at(0);
                }
                if(ids._time_extended){//TODO is this if needed here?
                    auto nb_sep1 = count(_rev_indices->front().begin(), _rev_indices->front().end(), ',');
                    auto nb_sep2 = count(key.begin(), key.end(), ',');
                    if(nb_sep2>nb_sep1){//some extra indices to be removed
                        auto pos = nthOccurrence(key, ",", nb_sep2-nb_sep1);
                        key = key.substr(pos+1,key.size()-1);
                    }
                }
                auto index = _indices->size();
                auto pp = param_::_indices->insert(make_pair<>(key, index));
                if(pp.second) { //new index inserted
                    _val->resize(max(_val->size(),index+1));
                    _dim[0] = max(_dim[0],_val->size());
                    _rev_indices->resize(_val->size());
                    _rev_indices->at(index) = key;
                    res._ids->at(0).push_back(index);
                    Debug(key << "; ");
                }
                else {
                    res._ids->at(0).push_back(pp.first->second);
                }
            }
            Debug(endl);
            res._dim[0]=res._ids->at(0).size();
            if(res._dim[0]>0){
                res._name += ".in_indices("+vec1._indices->at(0)+","+vec1._indices->at(vec1._indices->size()-1)+")";
            }
            else{
                res._name += ".empty_index_set";
            }
            res._unique_id = make_tuple<>(res._id,in_, typeid(indices).hash_code(), 0,res._dim[0]*res._dim.size());
            res._is_indexed = true;
            return res;
        }
        
        template<typename Tobj> param in(const vector<Tobj*>& vec) {
            param res(this->_name);
            res._id = this->_id;
            res._vec_id = this->_vec_id;
            res._intype = this->_intype;
            res._range = this->_range;
            res._val = this->_val;
            res._is_vector = this->_is_vector;
            res._is_matrix = this->_is_matrix;
            res._is_transposed = _is_transposed;
            res._rev_indices = this->_rev_indices; res._indices = this->_indices;
            if(vec.empty()){
                DebugOn("In function param.in(const vector<Tobj*>& vec), vec is empty!\n. Creating and empty variable! Check your sum/product operators.\n");
                res._name += "EMPTY_VAR";
                res._is_indexed = true;
                return res;
            }
            DebugOff(_name << " = ");
            string key, excluded;
            for(auto it = vec.begin(); it!= vec.end(); it++) {
                if(!(*it)->_active) {
                    excluded += (*it)->_name + ",";
                    continue;
                }
                key = (*it)->_name;
                auto index = _indices->size();
                auto pp = param_::_indices->insert(make_pair<>(key, index));
                if(pp.second) { //new index inserted
                    _val->resize(max(_val->size(),index+1));
                    _dim[0] = max(_dim[0],_val->size());
                    _rev_indices->resize(_val->size());
                    _rev_indices->at(index) = key;
                    res._ids->at(0).push_back(index);
                }
                else {
                    res._ids->at(0).push_back(pp.first->second);
                }
            }
            DebugOff(endl);
            res._dim[0]=res._ids->at(0).size();
            res._name += ".in_" + vec.front()->_type_name + "," + excluded + "," + to_string(res._dim[0]);
            res._unique_id = make_tuple<>(res._id,in_, typeid(Tobj).hash_code(), 0,res._dim[0]);
            res._is_indexed = true;
            return res;
        }
        
    
        
        
        
        param in_arcs(const vector<Node*>& vec) {
            param res(this->_name);
            res._id = this->_id;
            res._vec_id = this->_vec_id;
            res._intype = this->_intype;
            res._range = this->_range;
            res._val = this->_val;
            res._rev_indices = this->_rev_indices; res._indices = this->_indices;
            res._is_vector = _is_vector;
            res._is_matrix = _is_matrix;
            res._is_transposed = _is_transposed;
            if(vec.empty()){
                DebugOn("In function param.in_arcs(const vector<Tobj*>& vec), vec is empty!\n. Creating and empty variable! Check your sum/product operators.\n");
                res._name += "EMPTY_VAR";
                res._is_indexed = true;
                return res;
            }
            DebugOff(_name << " = ");
            string key;
            unsigned inst = 0;
            for(auto it = vec.begin(); it!= vec.end(); it++) {
                if(!(*it)->_active) {
                    continue;
                }
                if (inst>0) {
                    res._ids->push_back(vector<unsigned>());
                    res._dim.push_back(0);
                }
                for (auto &a:(*it)->get_in()) {
                    if (!a->_active) {
                        continue;
                    }
                    key = a->_name;
                    auto index = _indices->size();
                    auto pp = param_::_indices->insert(make_pair<>(key, index));
                    if(pp.second) { //new index inserted
                        _val->resize(max(_val->size(),index+1));
                        _dim[0] = max(_dim[0],_val->size());
                        _rev_indices->resize(_val->size());
                        _rev_indices->at(index) = key;
                        res._ids->at(inst).push_back(index);
                    }
                    else {
                        res._ids->at(inst).push_back(pp.first->second);
                    }
                }
                res._dim[inst]=res._ids->at(inst).size();
                ++inst;
            }
            res._unique_id = make_tuple<>(res._id,in_arcs_, typeid(Node).hash_code(), 0,res._dim[0]);
            res._is_indexed = true;
            return res;
        }
        
        param in_arcs(const vector<Node*>& vec, const indices& T) {
            param res(this->_name);
            res._id = this->_id;
            res._vec_id = this->_vec_id;
            res._intype = this->_intype;
            res._range = this->_range;
            res._val = this->_val;
            res._rev_indices = this->_rev_indices; res._indices = this->_indices;
            res._is_vector = _is_vector;
            res._is_matrix = _is_matrix;
            res._is_transposed = _is_transposed;
            if(vec.empty()){
                DebugOn("In function param.in_arcs(const vector<Tobj*>& vec), vec is empty!\n. Creating and empty variable! Check your sum/product operators.\n");
                res._name += "EMPTY_VAR";
                res._is_indexed = true;
                return res;
            }
            DebugOff(_name << " = ");
            string key;
            unsigned inst = 0;
            for(auto it = vec.begin(); it!= vec.end(); it++) {
                if(!(*it)->_active) {
                    continue;
                }
                unsigned t=0;
                for (auto &idx: *T._indices) {
                    if (inst>0) {
                        res._ids->push_back(vector<unsigned>());
                        res._dim.push_back(0);
                    }
                    if(T._excluded_indices.count(t++)!=0){
                        continue;
                    }
                    for (auto &a:(*it)->get_in()) {
                        if (!a->_active) {
                            continue;
                        }
                        key = a->_name;
                        key += ","+idx;
                        auto index = _indices->size();
                        auto pp = param_::_indices->insert(make_pair<>(key, index));
                        if(pp.second) { //new index inserted
                            _val->resize(max(_val->size(),index+1));
                            _dim[0] = max(_dim[0],_val->size());
                            _rev_indices->resize(_val->size());
                            _rev_indices->at(index) = key;
                            res._ids->at(inst).push_back(index);
                        }
                        else {
                            res._ids->at(inst).push_back(pp.first->second);
                        }
                    }
                    res._dim[inst]=res._ids->at(inst).size();
                    ++inst;
                }
            }
            res._unique_id = make_tuple<>(res._id,in_arcs_, typeid(Node).hash_code(), 0,res._dim[0]);
            res._is_indexed = true;
            return res;
        }
        
        param out_arcs(const vector<Node*>& vec) {
            param res(this->_name);
            res._id = this->_id;
            res._vec_id = this->_vec_id;
            res._intype = this->_intype;
            res._range = this->_range;
            res._val = this->_val;
            res._rev_indices = this->_rev_indices; res._indices = this->_indices;
            res._is_vector = _is_vector;
            res._is_matrix = _is_matrix;
            res._is_transposed = _is_transposed;
            if(vec.empty()){
                DebugOn("In function param.out_arcs(const vector<Tobj*>& vec), vec is empty!\n. Creating and empty variable! Check your sum/product operators.\n");
                res._name += "EMPTY_VAR";
                res._is_indexed = true;
                return res;
            }
            DebugOff(_name << " = ");
            string key;
            unsigned inst = 0;
            for(auto it = vec.begin(); it!= vec.end(); it++) {
                if(!(*it)->_active) {
                    continue;
                }
                if (inst>0) {
                    res._ids->push_back(vector<unsigned>());
                    res._dim.push_back(0);
                }
                for (auto &a:(*it)->get_out()) {
                    if (!a->_active) {
                        continue;
                    }
                    key = a->_name;
                    auto index = _indices->size();
                    auto pp = param_::_indices->insert(make_pair<>(key, index));
                    if(pp.second) { //new index inserted
                        _val->resize(max(_val->size(),index+1));
                        _dim[0] = max(_dim[0],_val->size());
                        _rev_indices->resize(_val->size());
                        _rev_indices->at(index) = key;
                        res._ids->at(inst).push_back(index);
                    }
                    else {
                        res._ids->at(inst).push_back(pp.first->second);
                    }
                }
                res._dim[inst]=res._ids->at(inst).size();
                ++inst;
            }
            res._unique_id = make_tuple<>(res._id,out_arcs_, typeid(Node).hash_code(), 0,res._dim[0]);
            res._is_indexed = true;
            return res;
        }
        
        param out_arcs(const vector<Node*>& vec, const indices& T) {
            param res(this->_name);
            res._id = this->_id;
            res._vec_id = this->_vec_id;
            res._intype = this->_intype;
            res._range = this->_range;
            res._val = this->_val;
            res._rev_indices = this->_rev_indices; res._indices = this->_indices;
            res._is_vector = _is_vector;
            res._is_matrix = _is_matrix;
            res._is_transposed = _is_transposed;
            if(vec.empty()){
                DebugOn("In function param.out_arcs(const vector<Tobj*>& vec), vec is empty!\n. Creating and empty variable! Check your sum/product operators.\n");
                res._name += "EMPTY_VAR";
                res._is_indexed = true;
                return res;
            }
            DebugOff(_name << " = ");
            string key;
            unsigned inst = 0;
            for(auto it = vec.begin(); it!= vec.end(); it++) {
                if(!(*it)->_active) {
                    continue;
                }
                auto t = 0;
                for (auto &idx: *T._indices) {
                    if (inst>0) {
                        res._ids->push_back(vector<unsigned>());
                        res._dim.push_back(0);
                    }
                    if(T._excluded_indices.count(t++)!=0){
                        continue;
                    }
                    for (auto &a:(*it)->get_out()) {
                        if (!a->_active) {
                            continue;
                        }
                        key = a->_name;
                        key += ","+idx;
                        auto index = _indices->size();
                        auto pp = param_::_indices->insert(make_pair<>(key, index));
                        if(pp.second) { //new index inserted
                            _val->resize(max(_val->size(),index+1));
                            _dim[0] = max(_dim[0],_val->size());
                            _rev_indices->resize(_val->size());
                            _rev_indices->at(index) = key;
                            res._ids->at(inst).push_back(index);
                        }
                        else {
                            res._ids->at(inst).push_back(pp.first->second);
                        }
                    }
                    res._dim[inst]=res._ids->at(inst).size();
                    ++inst;
                }
            }
            res._unique_id = make_tuple<>(res._id,out_arcs_, typeid(Node).hash_code(), 0,res._dim[0]);
            res._is_indexed = true;
            return res;
        }
        
        
        param in_gens(const vector<Node*>& vec) {
            param res(this->_name);
            res._id = this->_id;
            res._vec_id = this->_vec_id;
            res._intype = this->_intype;
            res._range = this->_range;
            res._val = this->_val;
            res._rev_indices = this->_rev_indices; res._indices = this->_indices;
            res._is_vector = _is_vector;
            res._is_matrix = _is_matrix;
            res._is_transposed = _is_transposed;
            if(vec.empty()){
                DebugOn("In function param.in_arcs(const vector<Tobj*>& vec), vec is empty!\n. Creating and empty variable! Check your sum/product operators.\n");
                res._name += "EMPTY_VAR";
                res._is_indexed = true;
                return res;
            }
            DebugOff(_name << " = ");
            string key;
            unsigned inst = 0;
            for(auto it = vec.begin(); it!= vec.end(); it++) {
                if(!(*it)->_active) {
                    continue;
                }
                if (inst>0) {
                    res._ids->push_back(vector<unsigned>());
                    res._dim.push_back(0);
                }
                for (auto &g:(*it)->get_gens()) {
                    if (!g->_active) {
                        continue;
                    }
                    key = g->_name;
                    auto index = _indices->size();
                    auto pp = param_::_indices->insert(make_pair<>(key, index));
                    if(pp.second) { //new index inserted
                        _val->resize(max(_val->size(),index+1));
                        _dim[0] = max(_dim[0],_val->size());
                        _rev_indices->resize(_val->size());
                        _rev_indices->at(index) = key;
                        res._ids->at(inst).push_back(index);
                    }
                    else {
                        res._ids->at(inst).push_back(pp.first->second);
                    }
                }
                res._dim[inst]=res._ids->at(inst).size();
                ++inst;
            }
            res._unique_id = make_tuple<>(res._id,in_gens_, typeid(Node).hash_code(), 0,res._dim[0]);
            res._is_indexed = true;
            return res;
        }
        
        
        param in_pot_gens(const vector<Node*>& vec) {
            param res(this->_name);
            res._id = this->_id;
            res._vec_id = this->_vec_id;
            res._intype = this->_intype;
            res._range = this->_range;
            res._val = this->_val;
            res._rev_indices = this->_rev_indices; res._indices = this->_indices;
            res._is_vector = _is_vector;
            res._is_matrix = _is_matrix;
            res._is_transposed = _is_transposed;
            if(vec.empty()){
                DebugOn("In function param.in_pot_gens(const vector<Tobj*>& vec), vec is empty!\n. Creating and empty variable! Check your sum/product operators.\n");
                res._name += "EMPTY_VAR";
                res._is_indexed = true;
                return res;
            }
            DebugOff(_name << " = ");
            string key;
            unsigned inst = 0;
            for(auto it = vec.begin(); it!= vec.end(); it++) {
                if(!(*it)->_active) {
                    continue;
                }
                if (inst>0) {
                    res._ids->push_back(vector<unsigned>());
                    res._dim.push_back(0);
                }
                for (auto &g:(*it)->get_pot_gens()) {
                    if (!g->_active) {
                        continue;
                    }
                    key = g->_name;
                    auto index = _indices->size();
                    auto pp = param_::_indices->insert(make_pair<>(key, index));
                    if(pp.second) { //new index inserted
                        _val->resize(max(_val->size(),index+1));
                        _dim[0] = max(_dim[0],_val->size());
                        _rev_indices->resize(_val->size());
                        _rev_indices->at(index) = key;
                        res._ids->at(inst).push_back(index);
                    }
                    else {
                        res._ids->at(inst).push_back(pp.first->second);
                    }
                }
                res._dim[inst]=res._ids->at(inst).size();
                ++inst;
            }
            res._unique_id = make_tuple<>(res._id,in_pot_gens_, typeid(Node).hash_code(), 0,res._dim[0]);
            res._is_indexed = true;
            return res;
        }
        
        param in_pot_bats(const vector<Node*>& vec) {
            param res(this->_name);
            res._id = this->_id;
            res._vec_id = this->_vec_id;
            res._intype = this->_intype;
            res._range = this->_range;
            res._val = this->_val;
            res._rev_indices = this->_rev_indices; res._indices = this->_indices;
            res._is_vector = _is_vector;
            res._is_matrix = _is_matrix;
            res._is_transposed = _is_transposed;
            if(vec.empty()){
                DebugOn("In function param.in_pot_bats(const vector<Tobj*>& vec), vec is empty!\n. Creating and empty variable! Check your sum/product operators.\n");
                res._name += "EMPTY_VAR";
                res._is_indexed = true;
                return res;
            }
            DebugOff(_name << " = ");
            string key;
            unsigned inst = 0;
            for(auto it = vec.begin(); it!= vec.end(); it++) {
                if(!(*it)->_active) {
                    continue;
                }
                if (inst>0) {
                    res._ids->push_back(vector<unsigned>());
                    res._dim.push_back(0);
                }
                for (auto &g:(*it)->get_pot_bats()) {
                    if (!g->_active) {
                        continue;
                    }
                    key = g->_name;
                    auto index = _indices->size();
                    auto pp = param_::_indices->insert(make_pair<>(key, index));
                    if(pp.second) { //new index inserted
                        _val->resize(max(_val->size(),index+1));
                        _dim[0] = max(_dim[0],_val->size());
                        _rev_indices->resize(_val->size());
                        _rev_indices->at(index) = key;
                        res._ids->at(inst).push_back(index);
                    }
                    else {
                        res._ids->at(inst).push_back(pp.first->second);
                    }
                }
                res._dim[inst]=res._ids->at(inst).size();
                ++inst;
            }
            res._unique_id = make_tuple<>(res._id,in_pot_bats_, typeid(Node).hash_code(), 0,res._dim[0]);
            res._is_indexed = true;
            return res;
        }
        
        param in_gens(const vector<Node*>& vec, const indices& T) {
            param res(this->_name);
            res._id = this->_id;
            res._vec_id = this->_vec_id;
            res._intype = this->_intype;
            res._range = this->_range;
            res._val = this->_val;
            res._rev_indices = this->_rev_indices; res._indices = this->_indices;
            res._is_vector = _is_vector;
            res._is_matrix = _is_matrix;
            res._is_transposed = _is_transposed;
            if(vec.empty()){
                DebugOn("In function param.in_arcs(const vector<Tobj*>& vec), vec is empty!\n. Creating and empty variable! Check your sum/product operators.\n");
                res._name += "EMPTY_VAR";
                res._is_indexed = true;
                return res;
            }
            DebugOff(_name << " = ");
            string key;
            unsigned inst = 0;
            for(auto it = vec.begin(); it!= vec.end(); it++) {
                if(!(*it)->_active) {
                    continue;
                }
                auto t = 0;
                for (auto &idx: *T._indices) {
                    if (inst>0) {
                        res._ids->push_back(vector<unsigned>());
                        res._dim.push_back(0);
                    }
                    if(T._excluded_indices.count(t++)!=0){
                        continue;
                    }
                    for (auto &g:(*it)->get_gens()) {
                        if (!g->_active) {
                            continue;
                        }
                        key = g->_name;
                        key += ","+idx;
                        auto index = _indices->size();
                        auto pp = param_::_indices->insert(make_pair<>(key, index));
                        if(pp.second) { //new index inserted
                            _val->resize(max(_val->size(),index+1));
                            _dim[0] = max(_dim[0],_val->size());
                            _rev_indices->resize(_val->size());
                            _rev_indices->at(index) = key;
                            res._ids->at(inst).push_back(index);
                        }
                        else {
                            res._ids->at(inst).push_back(pp.first->second);
                        }
                    }
                    res._dim[inst]=res._ids->at(inst).size();
                    ++inst;
                }
            }
            res._unique_id = make_tuple<>(res._id,in_gens_, typeid(Node).hash_code(), 0,res._dim[0]);
            res._is_indexed = true;
            return res;
        }
        
        param in_bats(const vector<Node*>& vec, const indices& T) {//TODO MERGE THIS WITH IN GENS AND OTHER FUNCS
            param res(this->_name);
            res._id = this->_id;
            res._vec_id = this->_vec_id;
            res._intype = this->_intype;
            res._range = this->_range;
            res._val = this->_val;
            res._rev_indices = this->_rev_indices; res._indices = this->_indices;
            res._is_vector = _is_vector;
            res._is_matrix = _is_matrix;
            res._is_transposed = _is_transposed;
            if(vec.empty()){
                DebugOn("In function param.in_bats(const vector<Tobj*>& vec, const indices& T), vec is empty!\n. Creating and empty variable! Check your sum/product operators.\n");
                res._name += "EMPTY_VAR";
                res._is_indexed = true;
                return res;
            }
            DebugOff(_name << " = ");
            string key;
            unsigned inst = 0;
            for(auto it = vec.begin(); it!= vec.end(); it++) {
                if(!(*it)->_active) {
                    continue;
                }
                auto t = 0;
                for (auto &idx: *T._indices) {
                    if (inst>0) {
                        res._ids->push_back(vector<unsigned>());
                        res._dim.push_back(0);
                    }
                    if(T._excluded_indices.count(t++)!=0){
                        continue;
                    }
                    for (auto g:(*it)->get_bats()) {
                        if (!g->_active) {
                            continue;
                        }
                        key = g->_name;
                        key += ","+idx;
                        auto index = _indices->size();
                        auto pp = param_::_indices->insert(make_pair<>(key, index));
                        if(pp.second) { //new index inserted
                            _val->resize(max(_val->size(),index+1));
                            _dim[0] = max(_dim[0],_val->size());
                            _rev_indices->resize(_val->size());
                            _rev_indices->at(index) = key;
                            res._ids->at(inst).push_back(index);
                        }
                        else {
                            res._ids->at(inst).push_back(pp.first->second);
                        }
                    }
                    res._dim[inst]=res._ids->at(inst).size();
                    ++inst;
                }
            }
            res._unique_id = make_tuple<>(res._id,in_bats_, typeid(Node).hash_code(), 0,res._dim[0]);
            res._is_indexed = true;
            return res;
        }
        
        param in_pv(const vector<Node*>& vec, const indices& T) {//TODO MERGE THIS WITH IN GENS AND OTHER FUNCS
            param res(this->_name);
            res._id = this->_id;
            res._vec_id = this->_vec_id;
            res._intype = this->_intype;
            res._range = this->_range;
            res._val = this->_val;
            res._rev_indices = this->_rev_indices; res._indices = this->_indices;
            res._is_vector = _is_vector;
            res._is_matrix = _is_matrix;
            res._is_transposed = _is_transposed;
            if(vec.empty()){
                DebugOn("In function param.in_pv(const vector<Tobj*>& vec, const indices& T), vec is empty!\n. Creating and empty variable! Check your sum/product operators.\n");
                res._name += "EMPTY_VAR";
                res._is_indexed = true;
                return res;
            }
            DebugOff(_name << " = ");
            string key;
            unsigned inst = 0;
            for(auto it = vec.begin(); it!= vec.end(); it++) {
                if(!(*it)->_active) {
                    continue;
                }
                auto t = 0;
                for (auto &idx: *T._indices) {
                    if (inst>0) {
                        res._ids->push_back(vector<unsigned>());
                        res._dim.push_back(0);
                    }
                    if(T._excluded_indices.count(t++)!=0){
                        continue;
                    }
                    for (auto &g:(*it)->get_pv()) {
                        if (!g->_active) {
                            continue;
                        }
                        key = g->_name;
                        key += ","+idx;
                        auto index = _indices->size();
                        auto pp = param_::_indices->insert(make_pair<>(key, index));
                        if(pp.second) { //new index inserted
                            _val->resize(max(_val->size(),index+1));
                            _dim[0] = max(_dim[0],_val->size());
                            _rev_indices->resize(_val->size());
                            _rev_indices->at(index) = key;
                            res._ids->at(inst).push_back(index);
                        }
                        else {
                            res._ids->at(inst).push_back(pp.first->second);
                        }
                    }
                    res._dim[inst]=res._ids->at(inst).size();
                    ++inst;
                }
            }
            res._unique_id = make_tuple<>(res._id,in_pv_, typeid(Node).hash_code(), 0,res._dim[0]);
            res._is_indexed = true;
            return res;
        }
        
        
        param in_wind(const vector<Node*>& vec, const indices& T) {//TODO MERGE THIS WITH IN GENS AND OTHER FUNCS
            param res(this->_name);
            res._id = this->_id;
            res._vec_id = this->_vec_id;
            res._intype = this->_intype;
            res._range = this->_range;
            res._val = this->_val;
            res._rev_indices = this->_rev_indices; res._indices = this->_indices;
            res._is_vector = _is_vector;
            res._is_matrix = _is_matrix;
            res._is_transposed = _is_transposed;
            if(vec.empty()){
                DebugOn("In function param.in_wind(const vector<Tobj*>& vec, const indices& T), vec is empty!\n. Creating and empty variable! Check your sum/product operators.\n");
                res._name += "EMPTY_VAR";
                res._is_indexed = true;
                return res;
            }
            DebugOff(_name << " = ");
            string key;
            unsigned inst = 0;
            for(auto it = vec.begin(); it!= vec.end(); it++) {
                if(!(*it)->_active) {
                    continue;
                }
                auto t = 0;
                for (auto &idx: *T._indices) {
                    if (inst>0) {
                        res._ids->push_back(vector<unsigned>());
                        res._dim.push_back(0);
                    }
                    if(T._excluded_indices.count(t++)!=0){
                        continue;
                    }
                    for (auto &g:(*it)->get_wind()) {
                        if (!g->_active) {
                            continue;
                        }
                        key = g->_name;
                        key += ","+idx;
                        auto index = _indices->size();
                        auto pp = param_::_indices->insert(make_pair<>(key, index));
                        if(pp.second) { //new index inserted
                            _val->resize(max(_val->size(),index+1));
                            _dim[0] = max(_dim[0],_val->size());
                            _rev_indices->resize(_val->size());
                            _rev_indices->at(index) = key;
                            res._ids->at(inst).push_back(index);
                        }
                        else {
                            res._ids->at(inst).push_back(pp.first->second);
                        }
                    }
                    res._dim[inst]=res._ids->at(inst).size();
                    ++inst;
                }
            }
            res._unique_id = make_tuple<>(res._id,in_wind_, typeid(Node).hash_code(), 0,res._dim[0]);
            res._is_indexed = true;
            return res;
        }
        
        template<typename Tobj> param in_pairs(const vector<Tobj*>& vec) {
            param res(this->_name);
            res._id = this->_id;
            res._vec_id = this->_vec_id;
            res._intype = this->_intype;
            res._range = this->_range;
            res._val = this->_val;
            res._is_vector = this->_is_vector;
            res._is_matrix = this->_is_matrix;
            res._is_transposed = _is_transposed;
            res._rev_indices = this->_rev_indices; res._indices = this->_indices;
            if(vec.empty()){
                DebugOff("In function param.in_pairs(const vector<Tobj*>& vec), vec is empty!\n. Creating and empty variable! Check your sum/product operators.\n");
                res._name += "EMPTY_VAR";
                return res;
            }
            DebugOff(_name << " = ");
            string key;
            for(auto it = vec.begin(); it!= vec.end(); it++) {
                if(!(*it)->_active) {
                    continue;
                }
                key = (*it)->_src->_name + "," + (*it)->_dest->_name;
                DebugOff(key<< ", ");
                auto index = _indices->size();
                auto pp = param_::_indices->insert(make_pair<>(key, index));
                if(pp.second) { //new index inserted
                    _val->resize(max(_val->size(),index+1));
                    _dim[0] = max(_dim[0],_val->size());
                    _rev_indices->resize(_val->size());
                    _rev_indices->at(index) = key;
                    res._ids->at(0).push_back(index);
                }
                else {
                    res._ids->at(0).push_back(pp.first->second);
                }
            }
            res._dim[0]=res._ids->at(0).size();
            DebugOff(endl);
            res._name += ".in_" + vec.front()->_type_name;
            res._unique_id = make_tuple<>(res._id,in_,typeid(Tobj).hash_code(), 0,res._dim[0]);
            res._is_indexed = true;
            return res;
        }
        
        
        template<typename Tobj> param in_pairs(const vector<Tobj*>& vec, const indices& ids) {
            param res(this->_name);
            res._id = this->_id;
            res._vec_id = this->_vec_id;
            res._intype = this->_intype;
            res._range = this->_range;
            res._val = this->_val;
            res._is_vector = this->_is_vector;
            res._is_matrix = this->_is_matrix;
            res._is_transposed = _is_transposed;
            res._rev_indices = this->_rev_indices; res._indices = this->_indices;
            if(vec.empty()){
                DebugOff("In function param.in_pairs(const vector<Tobj*>& vec), vec is empty!\n. Creating and empty variable! Check your sum/product operators.\n");
                res._name += "EMPTY_VAR";
                return res;
            }
            DebugOff(_name << " = ");
            string key;
            for(auto it = vec.begin(); it!= vec.end(); it++) {
                if(!(*it)->_active) {
                    continue;
                }
                auto t = 0;
                for (auto &idx: *ids._indices) {
                    if(ids._excluded_indices.count(t++)!=0){
                        continue;
                    }
                    key = (*it)->_src->_name + "," + (*it)->_dest->_name;
                    key += ","+idx;
                    DebugOff(key<< ", ");
                    
                    auto index = _indices->size();
                    auto pp = param_::_indices->insert(make_pair<>(key, index));
                    if(pp.second) { //new index inserted
                        _val->resize(max(_val->size(),index+1));
                        _dim[0] = max(_dim[0],_val->size());
                        _rev_indices->resize(_val->size());
                        _rev_indices->at(index) = key;
                        res._ids->at(0).push_back(index);
                    }
                    else {
                        res._ids->at(0).push_back(pp.first->second);
                    }
                }
            }
            res._dim[0]=res._ids->at(0).size();
            DebugOff(endl);
            res._name += ".in_" + vec.front()->_type_name;
            res._unique_id = make_tuple<>(res._id,in_,typeid(Tobj).hash_code(), 0,res._dim[0]);
            res._is_indexed = true;
            return res;
        }
        
        
        template<typename Tobj> param in_pairs(const vector<Tobj>& vec) {
            param res(this->_name);
            res._id = this->_id;
            res._vec_id = this->_vec_id;
            res._intype = this->_intype;
            res._range = this->_range;
            res._val = this->_val;
            res._is_vector = this->_is_vector;
            res._is_matrix = this->_is_matrix;
            res._is_transposed = _is_transposed;
            res._rev_indices = this->_rev_indices; res._indices = this->_indices;
            if(vec.empty()){
                DebugOff("In function param.in_pairs(const vector<Tobj*>& vec), vec is empty!\n. Creating and empty variable! Check your sum/product operators.\n");
                res._name += "EMPTY_VAR";
                return res;
            }
            DebugOff(_name << " = ");
            string key;
            for(auto it = vec.begin(); it!= vec.end(); it++) {
                if(!(*it)._active) {
                    continue;
                }
                key = (*it)._src->_name + "," + (*it)._dest->_name;
                DebugOff(key<< ", ");
                auto index = _indices->size();
                auto pp = param_::_indices->insert(make_pair<>(key, index));
                if(pp.second) { //new index inserted
                    _val->resize(max(_val->size(), index+1));
                    _dim[0] = max(_dim[0],_val->size());
                    _rev_indices->resize(_val->size());
                    _rev_indices->at(index) = key;
                    res._ids->at(0).push_back(index);
                }
                else {
                    res._ids->at(0).push_back(pp.first->second);
                }
            }
            DebugOff(endl);
            res._dim[0]=res._ids->at(0).size();
            res._name += ".in_" + vec.front()._type_name;
            res._unique_id = make_tuple<>(res._id,in_,typeid(Tobj).hash_code(), 0,res._dim[0]);
            res._is_indexed = true;
            return res;
        }
        
        template<typename Tobj> param in_pairs(const vector<Tobj*>& vec, unsigned T) {
            assert(T > 0);
            param res(this->_name);
            res._id = this->_id;
            res._vec_id = this->_vec_id;
            res._intype = this->_intype;
            res._range = this->_range;
            res._val = this->_val;
            res._is_vector = this->_is_vector;
            res._is_matrix = this->_is_matrix;
            res._is_transposed = _is_transposed;
            res._rev_indices = this->_rev_indices; res._indices = this->_indices;
            if(vec.empty()){
                DebugOff("In function param.in_pairs(const vector<Tobj*>& vec, unsigned T), vec is empty!\n. Creating and empty variable! Check your sum/product operators.\n");
                res._name += "EMPTY_VAR";
                return res;
            }
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
                    auto index = _indices->size();
                    auto pp = param_::_indices->insert(make_pair<>(key, index));
                    if(pp.second) { //new index inserted
                        _val->resize(max(_val->size(), index+1));
                        _dim[0] = max(_dim[0],_val->size());
                        _rev_indices->resize(_val->size());
                        _rev_indices->at(index) = key;
                        res._ids->at(0).push_back(index);
                    }
                    else {
                        res._ids->at(0).push_back(pp.first->second);
                    }
                }
            }
            DebugOff(endl);
            res._dim[0]=res._ids->at(0).size();
            res._name += ".in_" + string(typeid(Tobj).name()) + "_time_" + to_string(T);
            res._unique_id = make_tuple<>(res._id, in_time_, typeid(Tobj).hash_code(),0,res._dim[0]);
            res._is_indexed = true;
            return res;
        }
        
        template<typename Tobj> param from(const vector<Tobj*>& vec) {
            param res(this->_name);
            res._id = this->_id;
            res._vec_id = this->_vec_id;
            res._intype = this->_intype;
            res._range = this->_range;
            res._val = this->_val;
            res._is_vector = this->_is_vector;
            res._is_matrix = this->_is_matrix;
            res._is_transposed = _is_transposed;
            res._rev_indices = this->_rev_indices; res._indices = this->_indices;
            if(vec.empty()){
                DebugOff("In function param.from(const vector<Tobj*>& vec), vec is empty!\n. Creating and empty variable! Check your sum/product operators.\n");
                res._name += "EMPTY_VAR";
                return res;
            }
            DebugOff(_name << " = ");
            string key;
            for(auto it = vec.begin(); it!= vec.end(); it++) {
                if(!(*it)->_active) {
                    continue;
                }
                key = (*it)->_src->_name;
                DebugOff(key<< ", ");
                auto index = _indices->size();
                auto pp = param_::_indices->insert(make_pair<>(key, index));
                if(pp.second) { //new index inserted
                    _val->resize(max(_val->size(), index+1));
                    _dim[0] = max(_dim[0],_val->size());
                    _rev_indices->resize(_val->size());
                    _rev_indices->at(index) = key;
                    res._ids->at(0).push_back(index);
                }
                else {
                    res._ids->at(0).push_back(pp.first->second);
                }
            }
            DebugOff(endl);
            res._dim[0]=res._ids->at(0).size();
            if (get<1>(_unique_id)!=from_) {
                res._name += ".from";
            }
            res._name += "_"+ vec.front()->_type_name;
            res._unique_id = make_tuple<>(res._id,from_,typeid(Tobj).hash_code(), 0,res._dim[0]);
            res._is_indexed = true;
            return res;
        }
        
        template<typename Tobj> param from(const vector<Tobj*>& vec, const indices& ids) {
            param res(this->_name);
            res._id = this->_id;
            res._vec_id = this->_vec_id;
            res._intype = this->_intype;
            res._range = this->_range;
            res._val = this->_val;
            res._is_vector = this->_is_vector;
            res._is_matrix = this->_is_matrix;
            res._is_transposed = _is_transposed;
            res._rev_indices = this->_rev_indices; res._indices = this->_indices;
            if(vec.empty()){
                DebugOff("In function param.from(const vector<Tobj*>& vec), vec is empty!\n. Creating and empty variable! Check your sum/product operators.\n");
                res._name += "EMPTY_VAR";
                return res;
            }
            DebugOff(_name << " = ");
            string key;
            for(auto it = vec.begin(); it!= vec.end(); it++) {
                if(!(*it)->_active) {
                    continue;
                }
                auto t = 0;
                for (auto &idx: *ids._indices) {
                    if(ids._excluded_indices.count(t++)!=0){
                        continue;
                    }
                    key = (*it)->_src->_name;
                    key += ","+idx;
                    DebugOff(key<< ", ");
                    auto index = _indices->size();
                    auto pp = param_::_indices->insert(make_pair<>(key, index));
                    if(pp.second) { //new index inserted
                        _val->resize(max(_val->size(), index+1));
                        _dim[0] = max(_dim[0],_val->size());
                        _rev_indices->resize(_val->size());
                        _rev_indices->at(index) = key;
                        res._ids->at(0).push_back(index);
                    }
                    else {
                        res._ids->at(0).push_back(pp.first->second);
                    }
                }
            }
            DebugOff(endl);
            res._dim[0]=res._ids->at(0).size();
            if (get<1>(_unique_id)!=from_) {
                res._name += ".from";
            }
            res._name += "_"+ vec.front()->_type_name;
            res._unique_id = make_tuple<>(res._id,from_,typeid(Tobj).hash_code(), 0,res._dim[0]);
            res._is_indexed = true;
            return res;
        }
        
        template<typename Tobj> param from(const vector<Tobj>& vec) {
            param res(this->_name);
            res._id = this->_id;
            res._vec_id = this->_vec_id;
            res._intype = this->_intype;
            res._range = this->_range;
            res._val = this->_val;
            res._is_vector = this->_is_vector;
            res._is_matrix = this->_is_matrix;
            res._is_transposed = _is_transposed;
            res._rev_indices = this->_rev_indices; res._indices = this->_indices;
            if(vec.empty()){
                DebugOff("In function param.from(const vector<Tobj*>& vec), vec is empty!\n. Creating and empty variable! Check your sum/product operators.\n");
                res._name += "EMPTY_VAR";
                return res;
            }
            DebugOff(_name << " = ");
            string key;
            for(auto it = vec.begin(); it!= vec.end(); it++) {
                if(!(*it)._active) {
                    continue;
                }
                key = (*it)._src->_name;
                DebugOff(key<< ", ");
                auto index = _indices->size();
                auto pp = param_::_indices->insert(make_pair<>(key, index));
                if(pp.second) { //new index inserted
                    _val->resize(max(_val->size(), index+1));
                    _dim[0] = max(_dim[0],_val->size());
                    _rev_indices->resize(_val->size());
                    _rev_indices->at(index) = key;
                    res._ids->at(0).push_back(index);
                }
                else {
                    res._ids->at(0).push_back(pp.first->second);
                }
            }
            DebugOff(endl);
            res._dim[0]=res._ids->at(0).size();
            if (get<1>(_unique_id)!=from_) {
                res._name += ".from";
            }
            res._name += "_"+ vec.front()._type_name;
            res._unique_id = make_tuple<>(res._id,from_,typeid(Tobj).hash_code(), 0,res._dim[0]);
            res._is_indexed = true;
            return res;
        }

        
        
        
        template<typename Tobj> param to(const vector<Tobj*>& vec) {
            param res(this->_name);
            res._id = this->_id;
            res._vec_id = this->_vec_id;
            res._intype = this->_intype;
            res._range = this->_range;
            res._val = this->_val;
            res._is_vector = this->_is_vector;
            res._is_matrix = this->_is_matrix;
            res._is_transposed = _is_transposed;
            res._rev_indices = this->_rev_indices; res._indices = this->_indices;
            if(vec.empty()){
                DebugOff("In function param.to(const vector<Tobj*>& vec), vec is empty!\n. Creating and empty variable! Check your sum/product operators.\n");
                res._name += "EMPTY_VAR";
                return res;
            }
            DebugOff(_name << " = ");
            string key;
            for(auto it = vec.begin(); it!= vec.end(); it++) {
                if(!(*it)->_active) {
                    continue;
                }
                key = (*it)->_dest->_name;
                DebugOff(key<< ", ");
                auto index = _indices->size();
                auto pp = param_::_indices->insert(make_pair<>(key, index));
                if(pp.second) { //new index inserted
                    _val->resize(max(_val->size(),index+1));
                    _dim[0] = max(_dim[0],_val->size());
                    _rev_indices->resize(_val->size());
                    _rev_indices->at(index) = key;
                    res._ids->at(0).push_back(index);
                }
                else {
                    res._ids->at(0).push_back(pp.first->second);
                }
            }
            res._dim[0]=res._ids->at(0).size();
            if (get<1>(_unique_id)!=to_) {
                res._name += ".to";
            }
            res._name += "_"+ vec.front()->_type_name;
            res._unique_id = make_tuple<>(res._id,to_,typeid(Tobj).hash_code(), 0,res._dim[0]);
            res._is_indexed = true;
            return res;
        }
        
        template<typename Tobj> param to(const vector<Tobj*>& vec, const indices& ids) {
            param res(this->_name);
            res._id = this->_id;
            res._vec_id = this->_vec_id;
            res._intype = this->_intype;
            res._range = this->_range;
            res._val = this->_val;
            res._is_vector = this->_is_vector;
            res._is_matrix = this->_is_matrix;
            res._is_transposed = _is_transposed;
            res._rev_indices = this->_rev_indices; res._indices = this->_indices;
            if(vec.empty()){
                DebugOff("In function param.to(const vector<Tobj*>& vec), vec is empty!\n. Creating and empty variable! Check your sum/product operators.\n");
                res._name += "EMPTY_VAR";
                return res;
            }
            DebugOff(_name << " = ");
            string key;
            for(auto it = vec.begin(); it!= vec.end(); it++) {
                if(!(*it)->_active) {
                    continue;
                }
                auto t = 0;
                for (auto &idx: *ids._indices) {
                    if(ids._excluded_indices.count(t++)!=0){
                        continue;
                    }
                    key = (*it)->_dest->_name;
                    key += ","+idx;
                
                    DebugOff(key<< ", ");
                    auto index = _indices->size();
                    auto pp = param_::_indices->insert(make_pair<>(key, index));
                    if(pp.second) { //new index inserted
                        _val->resize(max(_val->size(),index+1));
                        _dim[0] = max(_dim[0],_val->size());
                        _rev_indices->resize(_val->size());
                        _rev_indices->at(index) = key;
                        res._ids->at(0).push_back(index);
                    }
                    else {
                        res._ids->at(0).push_back(pp.first->second);
                    }
                }
            }
            res._dim[0]=res._ids->at(0).size();
            if (get<1>(_unique_id)!=to_) {
                res._name += ".to";
            }
            res._name += "_"+ vec.front()->_type_name;
            res._unique_id = make_tuple<>(res._id,to_,typeid(Tobj).hash_code(), 0,res._dim[0]);
            res._is_indexed = true;
            return res;
        }
        
        template<typename Tobj> param to(const vector<Tobj>& vec) {
            param res(this->_name);
            res._id = this->_id;
            res._vec_id = this->_vec_id;
            res._intype = this->_intype;
            res._range = this->_range;
            res._val = this->_val;
            res._is_vector = this->_is_vector;
            res._is_matrix = this->_is_matrix;
            res._is_transposed = _is_transposed;
            res._rev_indices = this->_rev_indices; res._indices = this->_indices;
            if(vec.empty()){
                DebugOff("In function param.to(const vector<Tobj*>& vec), vec is empty!\n. Creating and empty variable! Check your sum/product operators.\n");
                res._name += "EMPTY_VAR";
                return res;
            }
            DebugOff(_name << " = ");
            string key;
            for(auto it = vec.begin(); it!= vec.end(); it++) {
                if(!(*it)._active) {
                    continue;
                }
                key = (*it)._dest->_name;
                DebugOff(key<< ", ");
                auto index = _indices->size();
                auto pp = param_::_indices->insert(make_pair<>(key, index));
                if(pp.second) { //new index inserted
                    _val->resize(max(_val->size(),index+1));
                    _dim[0] = max(_dim[0],_val->size());
                    _rev_indices->resize(_val->size());
                    _rev_indices->at(index) = key;
                    res._ids->at(0).push_back(index);
                }
                else {
                    res._ids->at(0).push_back(pp.first->second);
                }
            }
            res._dim[0]=res._ids->at(0).size();
            if (get<1>(_unique_id)!=to_) {
                res._name += ".to";
            }
            res._name += "_"+ vec.front()._type_name;
            res._unique_id = make_tuple<>(res._id,to_,typeid(Tobj).hash_code(), 0,res._dim[0]);
            res._is_indexed = true;
            return res;
        }
        
        
        param to() {
            auto res(*this);
            get<1>(res._unique_id) = to_;
            res._name += ".to";
            return res;
        }
        
        param from() {
            auto res(*this);
            get<1>(res._unique_id) = from_;
            res._name += ".from";
            return res;
        }
        
        param prev() {
            auto res(*this);
            get<1>(res._unique_id) = prev_;
            res._name += ".prev";
            return res;
        }
        
        param min_time() {
            auto res(*this);
            get<1>(res._unique_id) = min_time_;
            res._name += ".min_time";
            return res;
        }
        
        
        param in_arcs() {
            auto res(*this);
            get<1>(res._unique_id) = in_arcs_;
            res._name += ".incoming_arcs";
            return res;
        }
        
        param out_arcs() {
            auto res(*this);
            get<1>(res._unique_id) = out_arcs_;
            res._name += ".outgoing_arcs";
            return res;
        }
        
        param in_gens() {
            auto res(*this);
            get<1>(res._unique_id) = in_gens_;
            res._name += ".in_gens";
            return res;
        }
        
        param in_pot_gens() {
            auto res(*this);
            get<1>(res._unique_id) = in_pot_gens_;
            res._name += ".in_pot_gens";
            return res;
        }
        
        param in_bats() {
            auto res(*this);
            get<1>(res._unique_id) = in_bats_;
            res._name += ".in_bats";
            return res;
        }
        
        param in_pot_bats() {
            auto res(*this);
            get<1>(res._unique_id) = in_pot_bats_;
            res._name += ".in_pot_bats";
            return res;
        }
        
        param in_wind() {
            auto res(*this);
            get<1>(res._unique_id) = in_wind_;
            res._name += ".in_wind";
            return res;
        }
        
        param in_pv() {
            auto res(*this);
            get<1>(res._unique_id) = in_pv_;
            res._name += ".in_pv";
            return res;
        }
        
        param in_pairs() {
            auto res(*this);
            get<1>(res._unique_id) = in_pairs_;
            res._name += ".in_pairs";
            return res;
        }
        
        param excl(unsigned index) {
            param res(this->_name);
            res._id = this->_id;
            res._vec_id = this->_vec_id;
            res._intype = this->_intype;
            res._range = this->_range;
            res._val = this->_val;
            res._is_vector = this->_is_vector;
            res._is_matrix = this->_is_matrix;
            res._is_transposed = _is_transposed;
            res._rev_indices = this->_rev_indices; res._indices = this->_indices; res._indices = this->_indices;
            for (unsigned i = 0; i<get_nb_instances(); i++) {
                if (i!=index) {
                    res._ids->at(0).push_back(i);
                    res._dim[0]++;
                }
            }
            res._name += ".excl(" +  to_string(index) +")"; // _name and _unique_id should be really unique.
            res._unique_id = make_tuple<>(res._id, in_ ,typeid(type).hash_code(), index ,index);
            res._is_indexed = true;
            return res;
        }
        
        template<typename Tobj>
        param in(const Tobj nm, unsigned T) {
            param res(this->_name);
            res._id = this->_id;
            res._vec_id = this->_vec_id;
            res._intype = this->_intype;
            res._range = this->_range;
            res._val = this->_val;
            res._is_vector = this->_is_vector;
            res._is_matrix = this->_is_matrix;
            res._is_transposed = _is_transposed;
            res._rev_indices = this->_rev_indices; res._indices = this->_indices; res._indices = this->_indices;
            string key;
            if (nm->_active) {
                for (unsigned t = 0; t < T; t++) {
                    key = nm->_name;
                    key += ",";
                    key += to_string(t);
                    Debug("key: " << key << endl);
                    auto index = _indices->size();
                    auto pp = param_::_indices->insert(make_pair<>(key, index));
                    if(pp.second) { //new index inserted
                        _val->resize(max(_val->size(), index+1));
                        _dim[0] = max(_dim[0],_val->size());
                        _rev_indices->resize(_val->size());
                        _rev_indices->at(index) = key;
                        res._ids->at(0).push_back(index);
                    }
                    else {
                        res._ids->at(0).push_back(pp.first->second);
                    }
                }
            }
            res._dim[0]=res._ids->at(0).size();
            res._name += ".in_time_" +  nm->_name + "_time_" + to_string(T); // _name and _unique_id should be really unique.
            hash<string> str_hash;
            res._unique_id = make_tuple<>(res._id, in_time_,str_hash(nm->_name), 0,res._dim[0]);
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
            res._is_vector = this->_is_vector;
            res._is_matrix = this->_is_matrix;
            res._is_transposed = _is_transposed;
            res._rev_indices = this->_rev_indices; res._indices = this->_indices; res._indices = this->_indices;
            if(nodes.empty()){
                DebugOff("In function in_at(const vector<Tobj*>& nodes, unsigned t), nodes is empty!\n. Creating and empty variable! Check your sum/product operators.\n");
                res._name += "EMPTY_VAR";
                return res;
            }
            string key;
            for(auto it = nodes.begin(); it!= nodes.end(); it++) {
                if(!(*it)->_active) {
                    continue;
                }
                key = (*it)->_name;
                key += ",";
                key += to_string(t);
                auto index = _indices->size();
                auto pp = param_::_indices->insert(make_pair<>(key, index));
                if(pp.second) { //new index inserted
                    _val->resize(max(_val->size(), index+1));
                    _dim[0] = max(_dim[0],_val->size());
                    _rev_indices->resize(_val->size());
                    _rev_indices->at(index) = key;
                    res._ids->at(0).push_back(index);
                }
                else {
                    res._ids->at(0).push_back(pp.first->second);
                }
            }
            res._dim[0]=res._ids->at(0).size();
            res._name += ".in_" + string(typeid(Tobj).name()) + "_at_" + to_string(t);
            res._unique_id = make_tuple<>(res._id, in_at_, typeid(Tobj).hash_code(), 0,res._dim[0]);
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
            res._is_vector = this->_is_vector;
            res._is_matrix = this->_is_matrix;
            res._is_transposed = _is_transposed;
            res._rev_indices = this->_rev_indices; res._indices = this->_indices; res._indices = this->_indices;
            if(nodes.empty()){
                DebugOff("In function param.in(const vector<Tobj*>& nodes, unsigned T), nodes is empty!\n. Creating and empty variable! Check your sum/product operators.\n");
                res._name += "EMPTY_VAR";
                return res;
            }
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
                    auto index = _indices->size();
                    auto pp = param_::_indices->insert(make_pair<>(key, index));
                    if(pp.second) { //new index inserted
                        _val->resize(max(_val->size(), index+1));
                        _dim[0] = max(_dim[0],_val->size());
                        _rev_indices->resize(_val->size());
                        _rev_indices->at(index) = key;
                        res._ids->at(0).push_back(index);
                    }
                    else {
                        res._ids->at(0).push_back(pp.first->second);
                    }
                }
            }
            res._dim[0]=res._ids->at(0).size();
            res._name += ".in_" + string(typeid(Tobj).name()) + "_time_" + to_string(T);
            res._unique_id = make_tuple<>(res._id,in_time_,typeid(Tobj).hash_code(), 0,res._dim[0]);
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
            res._is_vector = this->_is_vector;
            res._is_matrix = this->_is_matrix;
            res._is_transposed = _is_transposed;
            res._rev_indices = this->_rev_indices; res._indices = this->_indices;
            if(arcs.empty()){
                DebugOff("In function param.from(const vector<Tobj*>& arcs, unsigned T), arcs is empty!\n. Creating and empty variable! Check your sum/product operators.\n");
                res._name += "EMPTY_VAR";
                return res;
            }
            string key;
            for (unsigned t = 0; t < T; t++) {
                for(auto it = arcs.begin(); it!= arcs.end(); it++) {
                    if(!(*it)->_active) { // || !(*it)->_src->_active || !(*it)->_dest->_active ) {
                        continue;
                    }
                    key = (*it)->_src->_name;
                    key += ",";
                    key += to_string(t);
                    auto index = _indices->size();
                    auto pp = param_::_indices->insert(make_pair<>(key,index));
                    if(pp.second) { //new index inserted
                        _val->resize(max(_val->size(), index+1));
                        _dim[0] = max(_dim[0],_val->size());
                        _rev_indices->resize(_val->size());
                        _rev_indices->at(index) = key;
                        res._ids->at(0).push_back(index);
                    }
                    else {
                        res._ids->at(0).push_back(pp.first->second);
                    }
                }
            }
            res._dim[0]=res._ids->at(0).size();
            res._name += ".from_" + string(typeid(Tobj).name()) + "_time_" + to_string(T);
            res._unique_id = make_tuple<>(res._id,from_time_,typeid(Tobj).hash_code(), 0,res._dim[0]);
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
            res._is_vector = this->_is_vector;
            res._is_matrix = this->_is_matrix;
            res._is_transposed = _is_transposed;
            res._rev_indices = this->_rev_indices; res._indices = this->_indices;
            if(arcs.empty()){
                DebugOff("In function param.to(const vector<Tobj*>& arcs, unsigned T), arcs is empty!\n. Creating and empty variable! Check your sum/product operators.\n");
                res._name += "EMPTY_VAR";
                return res;
            }
            string key;
            for (unsigned t = 0; t < T; t++) {
                for(auto it = arcs.begin(); it!= arcs.end(); it++) {
                    if(!(*it)->_active) {
                        continue;
                    }
                    key = (*it)->_dest->_name;
                    key += ",";
                    key += to_string(t);
                    auto index = _indices->size();
                    auto pp = param_::_indices->insert(make_pair<>(key,index));
                    if(pp.second) { //new index inserted
                        _val->resize(max(_val->size(),index+1));
                        _dim[0] = max(_dim[0],_val->size());
                        _rev_indices->resize(_val->size());
                        _rev_indices->at(index) = key;
                        res._ids->at(0).push_back(index);
                    }
                    else {
                        // already exists
                        res._ids->at(0).push_back(pp.first->second);
                    }
                }
            }
            res._dim[0]=res._ids->at(0).size();
            res._name += ".to_"  + string(typeid(Tobj).name()) + "_time_"+ to_string(T);
            res._unique_id = make_tuple<>(res._id, to_time_, typeid(Tobj).hash_code(), 0,res._dim[0]);
            res._is_indexed = true;
            return res;
        }
        
        // T copies of a parameter
        void time_expand(const indices& T) {
            _time_extended = true;
            auto dim = param_::get_dim()*(T.size() - T._excluded_indices.size());
            /* update the indices of the old parameter*/
            string key;
            auto indices = *_rev_indices;
            auto val_temp = *param::_val;
            set_size(dim);
            //        for (map<std::string, unsigned>::iterator it= param::_indices->begin(); it != param::_indices->end(); it++){
            //            key = it->first;
            //            map_temp.insert(make_pair(key, it->second));
            //        }
            //CLEAR OLD ENTRIES
            _indices->clear();
            _rev_indices->clear();
            _rev_indices->resize(dim);
            _ids->at(0).clear();
            //STORE NEW ENTRIES
            unsigned index = 0;
            for (auto i = 0; i<indices.size();i++) {
                for(unsigned t = 0; t < T.size(); t++ ) {
                    if(T._excluded_indices.count(t)!=0){
                        continue;
                    }
                    key = indices[i];
                    key += ",";
                    key += T._indices->at(t);
                    _val->at(index) = val_temp.at(i);
                    param_::_indices->insert(make_pair<>(key, index));
                    _rev_indices->at(index++) = key;
                }
            }
            _name += ".time_expanded";
            _unique_id = make_tuple<>(_id,unindexed_,typeid(type).hash_code(),0,dim*T.size());
        }
        
        param time_collapse(const indices& T) {
            param res(this->_name);
            res._time_extended = false;
            res._id = this->_id;
            res._vec_id = this->_vec_id;
            res._intype = this->_intype;
            res._range = this->_range;
            res._is_vector = this->_is_vector;
            res._is_matrix = this->_is_matrix;
            auto dim = param_::get_dim()/(T.size() - T._excluded_indices.size());
            auto nt = param_::get_dim()/dim;
            /* update the indices of the old parameter*/
            string key;
            auto indices = *_rev_indices;
            auto val_temp = *param::_val;
            res.set_size(dim);
            //        for (map<std::string, unsigned>::iterator it= param::_indices->begin(); it != param::_indices->end(); it++){
            //            key = it->first;
            //            map_temp.insert(make_pair(key, it->second));
            //        }
            //CLEAR OLD ENTRIES
            res._indices->clear();
            res._rev_indices->clear();
            res._rev_indices->resize(dim);
            res._ids->at(0).clear();
            res._val->resize(dim);
            auto tsize = T._indices->at(0).size();
            //STORE NEW ENTRIES
            unsigned index = 0;
            for (auto i = 0; i<dim*nt; i+=nt) {
                    key = indices[i];
                    key = key.substr(0,tsize-1);
                    res._val->at(index) = val_temp.at(i);
                    res.param_::_indices->insert(make_pair<>(key, index));
                    res._rev_indices->at(index++) = key;
            }
            res._name = _name.substr(0,_name.size()-string(".time_expanded").size());
            res._unique_id = make_tuple<>(get_id(),unindexed_,typeid(type).hash_code(),0,dim);
            return res;
        }
        
        void time_expand(unsigned T) {
            assert(T >= 1);
            auto dim = param_::get_dim();
            set_size(dim*T);
            /* update the indices of the old parameter*/
            string key;
            auto map_temp = *param::_indices;
            auto val_temp = *param::_val;
            //        for (map<std::string, unsigned>::iterator it= param::_indices->begin(); it != param::_indices->end(); it++){
            //            key = it->first;
            //            map_temp.insert(make_pair(key, it->second));
            //        }
            //CLEAR OLD ENTRIES
            _indices->clear();
            _ids->at(0).clear();
            //STORE NEW ENTRIES
            for(unsigned t = 0; t < T; t ++ ) {
                for (auto &entry: map_temp) {
                    key = entry.first;
                    key += ",";
                    key += to_string(t);
                    //_val->at(param_::_indices->size()) = _val->at(entry.second);
                    _val->at(param_::_indices->size()) = val_temp.at(entry.second);
                    param_::_indices->insert(make_pair<>(key, param_::_indices->size()));
                }
            }
            _name += ".time_expanded";
            _unique_id = make_tuple<>(_id,in_time_,typeid(type).hash_code(),0,dim*T);
        }
        
        
        
        /** Output */
        
        string to_str(size_t index1, size_t index2) const {
            if (_is_indexed) {
                return to_string(_val->at(_ids->at(index1).at(index2)));
            }
            else {
                return to_string(_val->at(index2));
            }
        }
        
        string to_str(size_t index) const {
            if (_is_indexed) {
                return to_string_with_precision(_val->at(_ids->at(0).at(index)), 10);
            }
            else {
                return to_string_with_precision(_val->at(index), 10);
            }
        }
        
        void print(size_t index) const {
            cout << to_str(index);
        }
        
        string to_str(bool vals=false) const {
            string str = get_name();
            auto name = str.substr(0, str.find("."));
            str = name;
            if (vals) {
                str += " = [ \n";
                if(_rev_indices->size()>0) {
                    for (unsigned idx = 0; idx < _val->size(); idx++) {
//                        if (!_is_relaxed || fabs(roundf(_val->at(idx)) - _val->at(idx)) > 1e-4) {
                            str += name+"(" + _rev_indices->at(idx) + ")=" + to_string_with_precision(_val->at(idx), 10);
                            str += " \n";
//                        }
                    }
                }else {
                    for(auto &indsp: *_indices) {
                        str += name+"("+indsp.first + "=" + to_string_with_precision(_val->at(indsp.second),10) + ")";
                        str += " ";
                    }
                }
                str += "];\n";
            }
            return str;
        }
        
        void print(bool vals=false) const {
            cout << this->to_str(vals);
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
