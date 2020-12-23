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
    #define _USE_MATH_DEFINES
    #include <math.h>
    #include <map>
    #include <set>
    #include <assert.h>
    #include <string>
    #include <iostream>
    #include <algorithm>
    #include <complex>      // std::complex

    inline int nthOccurrence(const std::string& str, const std::string& findMe, int nth)
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

    namespace gravity{
    #define EPS 0.00001
    #define Cpx complex<double>
    //#define Real double
    //#define Integer integer  //same name with a boost graph library.
    #define Binary bool
    #define Debug(x)
    #define DebugOn(x) cout << x
    #define Warning(x) cout << x
    #define WarningOff(x)
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
    typedef enum { id_, plus_, minus_, product_, div_, power_, cos_, sin_, sqrt_, exp_, log_, tan_, acos_, asin_, atan_, atan2_, abs_, df_abs_, relu_, unit_step_, min_, max_} OperatorType;  /* Operation type in the expression tree */

    typedef enum { R_, R_p_, C_} SpaceType;  /* Real, Positive Reals, Complex */

    typedef enum { ordered_pairs_, unordered_ } SetType;
    //    typedef enum { vec_=0, in_ordered_pairs_=1, from_ordered_pairs_=2, to_ordered_pairs_=3, in_arcs_=4, from_arcs_=5, to_arcs_=6, in_nodes_=7, in_set_=8, mask_=9, in_bags_=10, time_expand_ = 11, in_set_at_} IndexType;  /* Index type */

    typedef enum { unindexed_, in_, in_pairs_, out_, from_, to_, prev_, in_at_, in_time_, from_time_, to_time_, in_arcs_, out_arcs_, in_gens_, in_pot_gens_, in_bats_, in_pot_bats_,in_wind_, in_pv_, min_time_, excl_, matrix_} IndexType;  /* Index type */





    using namespace std;

    static double pi = 4.*atan(1.);


    class space{
    public:
        SpaceType       _type;
        vector<size_t>  _dim;
    };

    class R: public space{
    public:
        R(){};
        template<typename... Args>
        R(size_t t1, Args&&... args) {
            _type = R_;
            list<size_t> dims = {forward<size_t>(args)...};
            dims.push_front(t1);
            size_t size = dims.size();
            _dim.resize(size);
            auto it = dims.begin();
            size_t index = 0;
            while (it!=dims.end()) {
                _dim[index++] = *it++;
            }
        }
        
        R operator^(size_t n){return R(n);};
        
    };

    class R_p: public space{
    public:
        R_p(){};
        template<typename... Args>
        R_p(size_t t1, Args&&... args) {
            _type = R_p_;
            list<size_t> dims = {forward<size_t>(args)...};
            dims.push_front(t1);
            size_t size = dims.size();
            _dim.resize(size);
            auto it = dims.begin();
            size_t index = 0;
            while (it!=dims.end()) {
                _dim[index++] = *it++;
            }
        }
        
        R_p operator^(size_t n){return R_p(n);};
        
    };

    class C: public space{
    public:
        C(){};
        template<typename... Args>
        C(size_t t1, Args&&... args) {
            _type = C_;
            list<size_t> dims = {forward<size_t>(args)...};
            dims.push_front(t1);
            size_t size = dims.size();
            _dim.resize(size);
            auto it = dims.begin();
            size_t index = 0;
            while (it!=dims.end()) {
                _dim[index++] = *it++;
            }
        }
        
        C operator^(size_t n){return C(n);};
        /* TODO */
    };



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
    private:
        string                                  _name;/**< index set can be given a name */
        
    public:
        
        
        IndexType                               _type = unindexed_;/**< index type */
        bool                                    _time_extended = false;/*<< indices are time extended */
        size_t                                  _time_pos = 0;/*<< number of commas before time extension */
        shared_ptr<vector<size_t>>              _dim = nullptr;/*<< A vector storing the dimension of sub-indices */
        shared_ptr<vector<string>>              _keys = nullptr; /*<< A vector storing all the keys */
        
        shared_ptr<map<string,size_t>>          _keys_map = nullptr; /*<< A map storing all the indices, the size_t number indicates the right position in the _keys vector */
        
        set<size_t>                             _excluded_keys; /*<< A set storing all indices that should be excluded */
        shared_ptr<vector<vector<size_t>>>      _ids = nullptr;
        
        string get_name() const{
            string name = _name;
            if(_type==from_){
                name = "from("+name+")";
            }
            else if(_type==to_){
                name = "to("+name+")";
            }
            return name;
        }
        
        indices(string name){
            _name = name;
            _keys_map = make_shared<map<string,size_t>>();
            _keys = make_shared<vector<string>>();
            _dim = make_shared<vector<size_t>>();
            _dim->resize(2,0);
        }
        
        void set_name(const string& name){
            _name = name;
        }
        
        
        indices(int name){
            _name = to_string(name);
            _keys_map = make_shared<map<string,size_t>>();
            _keys = make_shared<vector<string>>();
            _dim = make_shared<vector<size_t>>();
            _dim->resize(2,0);
        }
        
        indices(){
            _keys_map = make_shared<map<string,size_t>>();
            _keys = make_shared<vector<string>>();
            _dim = make_shared<vector<size_t>>();
            _dim->resize(2,0);
        }
        
        indices(const ordered_pairs& pairs){
            auto n = pairs._keys.size();
            _keys_map = make_shared<map<string,size_t>>();
            _keys = make_shared<vector<string>>();
            _dim = make_shared<vector<size_t>>();
            _keys->resize(n);
            size_t index = 0;
            string key;
            for (int i = 0; i < n; i++){
                key = pairs._keys.at(index)->_name;
                (*_keys_map)[key]= index;
                (*_keys)[index++] = key;
            }
            _dim->resize(1);
            _dim->at(0) = n;
            _dim->resize(2,0);
        }
        
        
        /** Returns the number of comma-separated fields in each key */
        unsigned get_nb_entries() const{
            if(_keys->empty()){
                return 0;
            }
            return count(_keys->front().begin(), _keys->front().end(), ',') + 1;
        }
        
        /** Retain the entries starting at the ith position
         @param[in] start_position Use the substring starting after the start_position comma separator
         @param[in] nb_entries Number of comma separated entries to keep
         */
        indices from_ith(unsigned start_position, unsigned nb_entries) const{
            indices res(this->get_name()+"from_ith("+to_string(start_position)+")");
            string key;
            res._ids = make_shared<vector<vector<size_t>>>();
            res._ids->resize(1);
            if(_type == matrix_){/* If ids is matrix indexed */
                res._type = matrix_;
                auto nb_rows = this->get_nb_rows();
                res._ids->resize(nb_rows);
                for (size_t i = 0; i<nb_rows; i++) {
                    for (size_t j = 0; j<this->_ids->at(i).size(); j++) {
                        auto key_ref = _ids->at(i).at(j);
                        key = _keys->at(key_ref);
                        auto pos = nthOccurrence(key, ",", start_position);
                        if(pos>0){
                            key = key.substr(pos+1);
                        }
                        pos = nthOccurrence(key, ",", nb_entries);
                        if(pos>0){
                            key = key.substr(0,pos);
                        }
                        auto it = res._keys_map->find(key);
                        if (it == res._keys_map->end()){
                            res.insert(key);
                            res._ids->at(i).push_back(res._keys->size()-1);
                        }
                        else{
                            res._ids->at(i).push_back(it->second);
                        }
                    }
                }
            }
            else if(is_indexed()){/* If ids has key references, use those */
                for(auto &key_ref: _ids->at(0)){
                    key = _keys->at(key_ref);
                    auto pos = nthOccurrence(key, ",", start_position);
                    if(pos>0){
                        key = key.substr(pos+1);
                    }
                    pos = nthOccurrence(key, ",", nb_entries);
                    if(pos>0){
                        key = key.substr(0,pos);
                    }
                    auto it = res._keys_map->find(key);
                    if (it == res._keys_map->end()){
                        res.insert(key);
                        res._ids->at(0).push_back(res._keys->size()-1);
                    }
                    else{
                        res._ids->at(0).push_back(it->second);
                    }
                }
                
            }
            else {
                for(auto key: *_keys){
                    auto pos = nthOccurrence(key, ",", start_position);
                    if(pos>0){
                        key = key.substr(pos+1);
                    }
                    pos = nthOccurrence(key, ",", nb_entries);
                    if(pos>0){
                        key = key.substr(0,pos);
                    }
                    auto it = res._keys_map->find(key);
                    if (it == res._keys_map->end()){
                        res.insert(key);
                        res._ids->at(0).push_back(res._keys->size()-1);
                    }
                    else{
                        res._ids->at(0).push_back(it->second);
                    }
                }
            }
            return res;
        }
        
        /** Remove keys starting at the ith position and spanning nb_entries
         @param[in] start_position
         */
        indices ignore_ith(unsigned start_position, unsigned nb_entries) const{
            indices res(this->get_name()+"ignore_ith("+to_string(start_position)+","+to_string(nb_entries)+")");
            string key;
            res._ids = make_shared<vector<vector<size_t>>>();
            res._ids->resize(1);
            string first_part, last_part;
            if(_type == matrix_){/* If ids is matrix indexed */
                res._type = matrix_;
                auto nb_rows = this->get_nb_rows();
                res._ids->resize(nb_rows);
                for (size_t i = 0; i<nb_rows; i++) {
                    for (size_t j = 0; j<this->_ids->at(i).size(); j++) {
                        auto key_ref = _ids->at(i).at(j);
                        key = _keys->at(key_ref);
                        auto pos = nthOccurrence(key, ",", start_position);
                        first_part = key.substr(0,pos);
                        if(pos>0){
                            key = key.substr(pos+1);
                        }
                        pos = nthOccurrence(key, ",", nb_entries);
                        if(pos>0){
                            last_part = key.substr(pos+1);
                        }
                        if(first_part.size()>0 && last_part.size()>0){ /* stitch them together */
                            key = first_part+","+last_part;
                        }
                        else {
                            key = first_part+last_part;
                        }
                        auto it = res._keys_map->find(key);
                        if (it == res._keys_map->end()){
                            res.insert(key);
                            res._ids->at(i).push_back(res._keys->size()-1);
                        }
                        else{
                            res._ids->at(i).push_back(it->second);
                        }
                    }
                }
            }
            else if(is_indexed()){/* If current index set has key references, use those */
                for(auto &key_ref: _ids->at(0)){
                    key = _keys->at(key_ref);
                    auto pos = nthOccurrence(key, ",", start_position);
                    first_part = key.substr(0,pos);
                    if(pos>0){
                        key = key.substr(pos+1);
                    }
                    pos = nthOccurrence(key, ",", nb_entries);
                    if(pos>0){
                        last_part = key.substr(pos+1);
                    }
                    if(first_part.size()>0 && last_part.size()>0){ /* stitch them together */
                        key = first_part+","+last_part;
                    }
                    else {
                        key = first_part+last_part;
                    }
                    auto it = res._keys_map->find(key);
                    if (it == res._keys_map->end()){
                        res.insert(key);
                        res._ids->at(0).push_back(res._keys->size()-1);
                    }
                    else{
                        res._ids->at(0).push_back(it->second);
                    }
                }
                
            }
            else {
                string first_part, last_part;
                for(auto key: *_keys){
                    auto pos = nthOccurrence(key, ",", start_position);
                    first_part = key.substr(0,pos);
                    if(pos>0){
                        key = key.substr(pos+1);
                    }
                    pos = nthOccurrence(key, ",", nb_entries);
                    if(pos>0){
                        last_part = key.substr(pos+1);
                    }
                    if(first_part.size()>0 && last_part.size()>0){ /* stitch them together */
                        key = first_part+","+last_part;
                    }
                    else {
                        key = first_part+last_part;
                    }
                    auto it = res._keys_map->find(key);
                    if (it == res._keys_map->end()){
                        res.insert(key);
                        res._ids->at(0).push_back(res._keys->size()-1);
                    }
                    else{
                        res._ids->at(0).push_back(it->second);
                    }
                }
            }
            return res;
        }

        
        bool has_unique_ids() const{
            set<int> unique_ids;
            if(_type == matrix_){/* If ids is matrix indexed */
                auto nb_rows = this->get_nb_rows();
                for (size_t i = 0; i<nb_rows; i++) {
                    for (size_t j = 0; j<this->_ids->at(i).size(); j++) {
                        auto idx = _ids->at(i).at(j);
                        if(!unique_ids.insert(idx).second)
                            return false;
                    }
                }
            }
            else
            if(is_indexed()){/* If ids has key references, use those */
                for (auto &idx:_ids->at(0)) {
                    if(!unique_ids.insert(idx).second)
                        return false;
                }
            }
            return true;
        }
        
        /** Returns a vector of bools indicating if the ith reference is in ids and in this. The function iterates over key references in _ids. */
        vector<bool> get_common_refs(const indices& ids) const{
            vector<bool> res;
            // assert(_ids);
            set<size_t> unique_ids;
            if(_type == matrix_){/* If ids is matrix indexed */
                auto nb_rows = this->get_nb_rows();
                for (size_t i = 0; i<nb_rows; i++) {
                    for (size_t j = 0; j<this->_ids->at(i).size(); j++) {
                        auto idx = _ids->at(i).at(j);
                        if(ids.has_key(_keys->at(idx))){
                            res.push_back(true);
                        }
                        else {
                            res.push_back(false);
                        }
                    }
                }
            }
            else if(is_indexed()){/* If ids has key references, use those */
                for (auto &idx:_ids->at(0)) {
                    if(ids.has_key(_keys->at(idx))){
                        res.push_back(true);
                    }
                    else {
                        res.push_back(false);
                    }
                }
            }
            else {
                for (auto &key:*_keys) {
                    if(ids.has_key(key)){
                        res.push_back(true);
                    }
                    else {
                        res.push_back(false);
                    }
                }
            }
            return res;
        }
        
        bool has_different_order(const indices& ids) const{
            if(_type == matrix_){/* If ids is matrix indexed */
                throw invalid_argument("Function has_different_order(ids) cannot be called on a matrix-indexed set");
            }
            if(!is_indexed() || !ids.is_indexed() || size()!=ids.size()){
                throw invalid_argument("In has_different_order(ids), both index sets need to be indexed and of same size!");
            }
            return ids._ids->at(0)!=_ids->at(0);
        }
        
        vector<int> get_ids_order(const indices& ids) const{
            vector<int> res;
            if(_type == matrix_){/* If ids is matrix indexed */
                throw invalid_argument("Function get_ids_order(ids) cannot be called on a matrix-indexed set");
            }
            if(!is_indexed() || !ids.is_indexed() || size()!=ids.size()){
                throw invalid_argument("In get_ids_order(ids), both index sets need to be indexed and of same size!");
            }
            for (size_t idx:ids._ids->at(0)) {
                int position = 0;
                for (size_t idx2:_ids->at(0)) {
                    if(idx2==idx){
                        res.push_back(position);
                        break;
                    }
                    position++;
                }
                
            }
            return res;
        }
        
        void reorder_rows(const vector<int>& order) {
            vector<size_t> res;
            if(!is_indexed() || size()!=order.size()){
                throw invalid_argument("In reorder(ids), ids needs to be indexed and of same size!");
            }
            for (int position:order) {
                    res.push_back(_ids->at(0).at(position));
            }
            _ids->at(0) = res;
            _name = to_str();
        }
        
        
        bool has_key(const string& key) const{
             if(is_indexed()){/* If ids has key references, use those */
                 auto idx0 = _keys_map->at(key);
                 for (auto i= 0; i<_ids->size();i++) {
                     for (auto const idx :_ids->at(i)) {
                         if(idx0==idx)
                             return true;
                     }
                 }
                 return false;
             }
            return _keys_map->count(key)!=0;
        }
        
        /** Returns a vector of bools indicating if the ith reference is in this but not in ids. The function iterates over key references in _ids. */
        vector<bool> get_diff_refs(const indices& ids) const{
            vector<bool> res;
            // assert(_ids);
            set<size_t> unique_ids;
            if(_type == matrix_){/* If ids is matrix indexed */
                auto nb_rows = this->get_nb_rows();
                for (size_t i = 0; i<nb_rows; i++) {
                    for (size_t j = 0; j<this->_ids->at(i).size(); j++) {
                        auto idx = _ids->at(i).at(j);
                        if(!ids.has_key(_keys->at(idx))){
                            res.push_back(true);
                        }
                        else {
                            res.push_back(false);
                        }
                    }
                }
            }
            else if(is_indexed()){/* If ids has key references, use those */
                for (auto &idx:_ids->at(0)) {
                    if(!ids.has_key(_keys->at(idx))){
                        res.push_back(true);
                    }
                    else {
                        res.push_back(false);
                    }
                }
            }
            else {
                for (auto &key:*_keys) {
                    if(!ids.has_key(key)){
                        res.push_back(true);
                    }
                    else {
                        res.push_back(false);
                    }
                }
            }
            return res;
        }
        
        /** Returns an index set based on the references stored in _ids. The function iterates over key references in _ids and keeps only the unique entries */
        indices get_unique_keys() const{
            indices res(_name);
            res._type = unindexed_;
            res._keys = _keys;
            res._keys_map = _keys_map;
            res._excluded_keys = _excluded_keys;
            res._dim = make_shared<vector<size_t>>();
            res._dim->resize(1);
            res._dim->at(0) = res._keys->size();
            return res;
        }
        
        /* Return an index set repeating the key at position pos n times
         @param[in] n, number of time repeating the same key
         @param[in] pos, position of correponsing key
         */
        indices repeat_id(int n, int pos=0) const{
            auto key_id = get_id_inst(pos);
            indices res(*this);
            res.set_name(res.get_name() + "repeat_id(" + to_string(n)+ "," + to_string(pos) + ")");
            res._ids = make_shared<vector<vector<size_t>>>();
            res._ids->resize(1);
            res._ids->at(0).resize(n);
            for(int i = 0; i< n; i++){
                res._ids->at(0).at(i) = key_id;
            }
            return res;
        }
        
        /* Delete rows where keep[i] is false. */
        void filter_rows(const vector<bool>& keep){
            bool all_true = true;
            for (auto b: keep) {
                if(!b){
                    all_true = false;
                    break;
                }
            }
            if(all_true)
                return;
            string excluded = "\\{";
            if(_type == matrix_){/* If ids is matrix indexed */
                if(keep.size()!=get_nb_rows()){
                    WarningOff("in filter_refs(const vector<bool>& keep): keep has a different size than index set, aborting");
                    return;
                }
                shared_ptr<vector<vector<size_t>>> new_ids = make_shared<vector<vector<size_t>>>();
                auto nb_rows = this->get_nb_rows();
                for (size_t i = 0; i<nb_rows; i++) {
                    if(keep[i]==true){
                        new_ids->push_back(_ids->at(i));
                    }
                }
                _ids = new_ids;
            }
            else if(_ids){
                if(keep.size()!=_ids->at(0).size()){
                    WarningOff("in filter_refs(const vector<bool>& keep): keep has a different size than index set");
                    return;
                }
                vector<vector<size_t>> new_ids;
                new_ids.resize(1);
                for (auto idx = 0; idx<keep.size();idx++) {
                    if(keep[idx]){
                        new_ids.at(0).push_back(_ids->at(0).at(idx));
                    }
                    else {
                        excluded += "("+_keys->at(_ids->at(0).at(idx)) +"),";
                    }
                }
                *_ids = new_ids;
                _name = to_str();
            }
            else {
                if(keep.size()!=_keys->size()){
                    WarningOff("in filter_refs(const vector<bool>& keep): keep has a different size than index set");
                    return;
                }
                shared_ptr<vector<vector<size_t>>> new_ids = make_shared<vector<vector<size_t>>>();
                new_ids->resize(1);
                for (auto idx = 0; idx<keep.size();idx++) {
                    if(keep[idx]){
                        new_ids->at(0).push_back(idx);
                    }
                    else {
                        excluded += "("+_keys->at(idx) +"),";
                    }
                }
                _ids = new_ids;
                _name = to_str();
            }
        }
        
            /* Delete row if all correspondong cols in keep[i] are false. */
            void filter_matrix_rows(const vector<bool>& keep){
                string excluded = "\\{";
                if(_type == matrix_){
                    auto nb_rows = _ids->size();
                    auto nb_cols = keep.size()/nb_rows;
                    vector<vector<size_t>> new_ids;
                    new_ids.resize(nb_rows);
                    for (auto row = 0; row<nb_rows;row++) {
                        bool exclude_row = true;/* Excluderow if all cols are false */
                        for (auto col = 0; col<nb_cols;col++) {
                            if(keep[row*nb_cols+col]){
                                exclude_row = false;
                            }
                        }
                        if(exclude_row){
                            excluded += "("+_keys->at(_ids->at(row).at(0)) +"),";
                        }
                        else {
                            new_ids.at(row).push_back(_ids->at(row).at(0));
                        }
                    }
                    *_ids = new_ids;
                    _name += excluded.substr(0,excluded.size()-1) + "}";
                    remove_empty_rows();
                }
                else if(_ids){
                    if(keep.size()<=_ids->at(0).size()){
                        throw invalid_argument("in filter_matrix_rows(const vector<bool>& keep): dimension mismatch");
                    }
                    auto nb_rows = _ids->at(0).size();
                    auto nb_cols = keep.size()/nb_rows;
                    vector<vector<size_t>> new_ids;
                    new_ids.resize(1);
                    for (auto row = 0; row<nb_rows;row++) {
                        bool exclude_row = true;/* Excluderow if all cols are false */
                        for (auto col = 0; col<nb_cols;col++) {
                            if(keep[row*nb_cols+col]){
                                exclude_row = false;
                            }
                        }
                        if(exclude_row){
                            excluded += "("+_keys->at(_ids->at(0).at(row)) +"),";
                        }
                        else {
                            new_ids.at(0).push_back(_ids->at(0).at(row));
                        }
                    }
                    *_ids = new_ids;
                    _name += excluded.substr(0,excluded.size()-1) + "}";
                }
                else {
                    if(keep.size()<=_keys->size()){
                        throw invalid_argument("in filter_matrix_rows(const vector<bool>& keep): dimension mismatch");
                    }
                    auto nb_rows = _keys->size();
                    auto nb_cols = keep.size()/nb_rows;
                    shared_ptr<vector<vector<size_t>>> new_ids = make_shared<vector<vector<size_t>>>();
                    new_ids->resize(1);
                    for (auto row = 0; row<nb_rows;row++) {
                        bool exclude_row = true;/* Excluderow if all cols are false */
                        for (auto col = 0; col<nb_cols;col++) {
                            if(keep[row*nb_cols+col]){
                                exclude_row = false;
                            }
                        }
                        if(exclude_row){
                            excluded += "("+_keys->at(row) +"),";
                        }
                        else {
                            new_ids->at(0).push_back(row);
                        }
                    }
                    _ids = new_ids;
                    _name += excluded.substr(0,excluded.size()-1) + "}";
                }
            }
        
        /* Delete indices where keep[i] is false. */
        void filter_refs(const vector<bool>& keep){
            string excluded = "\\{";
            if(_type == matrix_){/* If ids is matrix indexed */
                if(keep.size()!=nb_keys()){
                    throw invalid_argument("in filter_refs(const vector<bool>& keep): keep has a different size than index set");
                }
                shared_ptr<vector<vector<size_t>>> new_ids = make_shared<vector<vector<size_t>>>();
                auto nb_rows = this->get_nb_rows();
                new_ids->resize(nb_rows);
                size_t idx = 0;
                for (size_t i = 0; i<nb_rows; i++) {
                    for (size_t j = 0; j<this->_ids->at(i).size(); j++) {
                        if(keep[idx++]==true){
                            new_ids->at(i).push_back(_ids->at(i).at(j));
                        }
                    }
                }
    //            shared_ptr<vector<vector<size_t>>> nnz_ids = make_shared<vector<vector<size_t>>>();
    //            for (size_t i = 0; i<nb_rows; i++) {
    //                if(new_ids->at(i).size()>0){
    //                    nnz_ids->push_back(new_ids->at(i));
    //                }
    //            }
                _ids = new_ids;
                remove_empty_rows();
            }
            else if(_ids){
                if(keep.size()!=_ids->at(0).size()){
                    throw invalid_argument("in filter_refs(const vector<bool>& keep): keep has a different size than index set");
                }
                vector<vector<size_t>> new_ids;
                new_ids.resize(1);
                for (auto idx = 0; idx<keep.size();idx++) {
                    if(keep[idx]){
                        new_ids.at(0).push_back(_ids->at(0).at(idx));
                    }
                    else {
                        excluded += "("+_keys->at(_ids->at(0).at(idx)) +"),";
                    }
                }
                *_ids = new_ids;
                _name += excluded.substr(0,excluded.size()-1) + "}";
            }
            else {
                if(keep.size()!=_keys->size()){
                    throw invalid_argument("in filter_refs(const vector<bool>& keep): keep has a different size than index set");
                }
                shared_ptr<vector<vector<size_t>>> new_ids = make_shared<vector<vector<size_t>>>();
                new_ids->resize(1);
                for (auto idx = 0; idx<keep.size();idx++) {
                    if(keep[idx]){
                        new_ids->at(0).push_back(idx);
                    }
                    else {
                        excluded += "("+_keys->at(idx) +"),";
                    }
                }
                _ids = new_ids;
                _name += excluded.substr(0,excluded.size()-1) + "}";
            }
        }
        
        bool is_matrix_indexed() const{
            return _type == matrix_;
        }
        
        void remove_empty_rows() {
            if(_type!=matrix_){
                throw invalid_argument("clean_empty_rows() can only be called on a matrix indexed set");
            }
            shared_ptr<vector<vector<size_t>>> new_ids = make_shared<vector<vector<size_t>>>();
            auto nb_rows = this->get_nb_rows();
            for (size_t i = 0; i<nb_rows; i++) {
                if(_ids->at(i).size()>0){
                    new_ids->push_back(move(_ids->at(i)));
                }
            }
            _ids = new_ids;
        }
        
        /** Returns a vector of bools indicating if a reference is unique so far. The function iterates over key references in _ids. */
        vector<bool> get_unique_refs() const{
            vector<bool> res;
            set<size_t> unique_ids;
            if(_ids){
                for (auto &idx:_ids->at(0)) {
                    if(unique_ids.insert(idx).second){
                        res.push_back(true);
                    }
                    else {
                        res.push_back(false);
                    }
                }
            }
            return res;
        }
        /* transform a matrix indexed set to a vector indexed one */
        void flatten() {
            if(_type != matrix_){
                return;
            }
            shared_ptr<vector<vector<size_t>>> new_ids = make_shared<vector<vector<size_t>>>();
            new_ids->resize(1);
            for (size_t i = 0; i<get_nb_rows(); i++) {
                for (size_t j = 0; j<this->_ids->at(i).size(); j++) {
                    new_ids->at(0).push_back(_ids->at(i).at(j));
                }
            }
            _ids = new_ids;
            _type = in_;
        }
    //    /** The function iterates over the ith key references in _ids and deletes the ones where keep[i] is false. */
    //    void filter_refs(const vector<bool>& keep) const{
    //        if(_ids){
    //            if(keep.size()!=_ids->at(0).size()){
    //                throw invalid_argument("in filter_refs(const vector<bool>& keep): keep has a different size than index set");
    //            }
    //            vector<vector<size_t>> new_ids;
    //            new_ids.resize(1);
    //            for (auto idx = 0; idx<keep.size();idx++) {
    //                if(keep[idx]){
    //                    new_ids.at(0).push_back(_ids->at(0).at(idx));
    //                }
    //            }
    //            *_ids = new_ids;
    //        }
    //    }
        
        /* Returns true if current index set is a strict subset of ids */
        bool is_subset(const indices & ids) const{
            return size()<ids.size() && is_subset_eq(ids);
        }
        
        
        /* Returns true if current index set is a subset of ids. Returns true if they are equal.*/
        bool is_subset_eq(const indices & ids) const{
            /*TODO: check if same set using set name!*/
            if(_type == matrix_){/* If ids is matrix indexed */
                auto nb_rows = this->get_nb_rows();
                for (size_t i = 0; i<nb_rows; i++) {
                    for (size_t j = 0; j<this->_ids->at(i).size(); j++) {
                        auto key_ref = _ids->at(i).at(j);
                        auto key = _keys->at(key_ref);
                        if(!ids.has_key(key)){
                            return false;
                        }
                    }
                }
                if(ids._type==matrix_){
                    if(nb_keys()!=ids.nb_keys())
                       return false;
                }
                else {
                    if(nb_keys()!=ids.size())
                    return false;
                }
            }
            else if(_ids){
                for (auto idx = 0; idx<size();idx++) {
                    if(!ids.has_key(_keys->at(_ids->at(0).at(idx)))){
                        return false;
                    }
                }
            }
            else {
                for (auto idx = 0; idx<size();idx++) {
                    if(!ids.has_key(_keys->at(idx))){
                        return false;
                    }
                }
            }
            if(ids._type==matrix_){
                if(size()!=ids.nb_keys())
                   return false;
            }
            return true;
        }
        
        /* Returns true if current index set is a superset of ids */
        bool is_superset(const indices & ids) const{
            return ids.is_subset(*this);
        }
        
        /* Transform the current index set into a matrix indexed set. If transformed (it was not a matrix already), the resulting set will have unidimentional columns. */
        void to_matrix(){
            if(_type!=matrix_){
                *this = get_matrix_form();
            }
        }
        
        /* Get a matrix version of the current index set. If transformed (it was not a matrix already), the resulting set will have unidimentional columns. */
        indices get_matrix_form() const {
            if(_type==matrix_){
                return *this;
            }
            indices res(_name);
            res._type = matrix_;
            res._keys = _keys;
            res._keys_map= _keys_map;
            auto nb_rows = this->size();
            for (size_t i = 0; i<nb_rows; i++) {
                res.add_empty_row();
                if(_ids){
                    res.add_in_row(i, _keys->at(_ids->at(0).at(i)));
                }
                else {
                    res.add_in_row(i, _keys->at(i));
                }
            }
            return res;
        }
        /** The function iterates over key references in _ids and keeps only the unique entries */
        void keep_unique_keys(){
            indices res(_name);
            set<size_t> unique_ids;
            if(_ids){
                for (auto &idx:_ids->at(0)) {
                    if(unique_ids.insert(idx).second){
                        res.add(_keys->at(idx));
                    }
                }
                *this=res;
            }
        }
        //        indices(const space& s){
        //            list<indices> l;
        //            _dim->resize(l.size());
        //            for (auto i = 0; i<s._dim.size(); i++) {
        //                l.push_back(indices(0,s._dim[i]-1));
        //                _dim->at(i) = s._dim[i];
        //            }
        //            *this = indices(l);
        //        }
        
        
        
        template<typename... Args>
        void init(string idx1, Args&&... args) {
            list<string> indices;
            indices = {forward<string>(args)...};
            indices.push_front(idx1);
            auto n = indices.size();
            _keys_map = make_shared<map<string,size_t>>();
            _keys = make_shared<vector<string>>();
            _keys->resize(n);
            _dim = make_shared<vector<size_t>>();
            _dim->resize(1);
            _dim->at(0) = n;
            auto it = indices.begin();
            for (size_t i= 0; i< n; i++) {
                (*_keys_map)[*it]= i;
                (*_keys)[i] = (*it);
                it++;
            }
        }
        
        
        bool operator==(const indices& cpy) const{
            if (_type != cpy._type || _time_extended!=cpy._time_extended || _time_pos != cpy._time_pos || *_dim!=*cpy._dim || _excluded_keys != cpy._excluded_keys || *_keys_map != *cpy._keys_map) return false;
            if(_ids==cpy._ids) return true; /* accounts for both being nullptr */
            if((_ids && !cpy._ids) || (cpy._ids && !_ids) || (*_ids != *cpy._ids)) return false;
            return true;
        }
        
        bool operator!=(const indices& cpy) const{
            return !(*this==cpy);
        }
        
        void shallow_copy(shared_ptr<indices> cpy){
            _name = cpy->_name;
            _type = cpy->_type;
            _keys_map = cpy->_keys_map;
            _keys = cpy->_keys;
            _dim = cpy->_dim;
            _excluded_keys = cpy->_excluded_keys;
            if(cpy->_ids){
                _ids = make_shared<vector<vector<size_t>>>(*cpy->_ids);
            }
            _time_extended = cpy->_time_extended;
            _time_pos = cpy->_time_pos;
        }
        
        void time_expand(const indices& T) {//Fix this to expand ids size
            _time_extended = true;
            auto dim = this->size()*(T.size() - T._excluded_keys.size());
            /* update the indices of the old parameter*/
            string key;
            auto keys = *_keys;
            //CLEAR OLD ENTRIES
            _keys->clear();
            _keys_map->clear();
            _keys->resize(dim);
            _ids->at(0).clear();
            //STORE NEW ENTRIES
            unsigned index = 0;
            for (auto i = 0; i<keys.size();i++) {
                for(unsigned t = 0; t < T.size(); t++ ) {
                    if(T._excluded_keys.count(t)!=0){
                        continue;
                    }
                    key = keys[i];
                    key += ",";
                    key += T._keys->at(t);
                    _keys_map->insert(make_pair<>(key, index));
                    _keys->at(index++) = key;
                }
            }
            _name += ".time_expanded";
        }
        
        indices& operator=(const indices& cpy){
            if(_name.empty())
                _name = cpy._name;
            _type = cpy._type;
            _keys_map = cpy._keys_map;
            _keys = cpy._keys;
            _dim = cpy._dim;
            _excluded_keys = cpy._excluded_keys;
            if(cpy._ids){
                _ids = make_shared<vector<vector<size_t>>>(*cpy._ids);
            }
            else{
                _ids = nullptr;
            }
            _time_extended = cpy._time_extended;
            _time_pos = cpy._time_pos;
            return *this;
        }
        
        indices deep_copy() const{
            indices cpy;
            cpy._name = _name;
            cpy._type = _type;
            cpy._dim = _dim;
            cpy._excluded_keys = _excluded_keys;
            if(_ids){
                cpy._ids = make_shared<vector<vector<size_t>>>(*_ids);
            }
            if(_keys){
                cpy._keys = make_shared<vector<string>>(*_keys);
            }
            if(_keys_map){
                cpy._keys_map = make_shared<map<string,size_t>>(*_keys_map);
            }
            cpy._time_extended = _time_extended;
            cpy._time_pos = _time_pos;
            return cpy;
        }
        
        indices& operator=(indices&& cpy){
            if(_name.empty())
                _name = cpy._name;
            _type = cpy._type;
            _keys_map = move(cpy._keys_map);
            _excluded_keys = move(cpy._excluded_keys);
            _keys = move(cpy._keys);
            _dim = move(cpy._dim);
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
        indices(const vector<Tobj*>& vec, bool include_inactive = false){
            _keys_map = make_shared<map<string,size_t>>();
            _keys = make_shared<vector<string>>();
            _dim = make_shared<vector<size_t>>();
            _dim->resize(1);
            size_t i = 0;
            for (auto idx:vec) {
                if(include_inactive || idx->_active){
                    (*_keys_map)[idx->_name]= i;
                    _keys->push_back(idx->_name);
                    i++;
                }
            }
            if (_keys->size()>0) {
                _name = vec.front()->_type_name;
            }
            _dim->at(0) = _keys->size();
        }
        
        template<typename Tobj>
        indices(const vector<Tobj>& vec, bool include_inactive = false){
            _keys_map = make_shared<map<string,size_t>>();
            _keys = make_shared<vector<string>>();
            _dim = make_shared<vector<size_t>>();
            _dim->resize(1);
            size_t i = 0;
            for (auto idx:vec) {
                if(include_inactive || idx._active){
                    (*_keys_map)[idx._name]= i;
                    _keys->push_back(idx._name);
                    i++;
                }
            }
            if (_keys->size()>0) {
                _name = vec.front()->_type_name;
            }
            _dim->at(0) = _keys->size();
        }
        
        void insert(const string& key){
            _keys->push_back(key);
            _keys_map->insert(make_pair<>(key,_keys->size()-1));
        }
        
        void extend(const indices& T) {
            if(!_ids){
                *this = indices(*this,T);
                return;
            }
            if (_ids->size()>1) {//double indexed
                auto dim = _ids->size();
                auto new_dim = T.size()*dim;
                _ids->resize(new_dim);
                for (auto i =dim; i<new_dim; i++) {
                    _ids->at(i).resize(_ids->at(i%dim).size());
                    for (auto j =0; j<_ids->at(i).size(); j++) {
                        _ids->at(i).at(j) = _ids->at(i%dim).at(j);
                    }
                }
            }
            else {
                auto dim = _ids->at(0).size();
                auto new_dim = T.size()*dim;
                _ids->at(0).resize(new_dim);
                for (auto i = dim; i<new_dim; i++) {
                    _ids->at(0).at(i) = _ids->at(0).at(i%dim);
                }
            }
        }
        
        indices(const list<indices>& vecs) {
            //            if (vecs.size()==2) {
            //                _type = matrix_;
            //            }
            _keys_map = make_shared<map<string,size_t>>();
            _keys = make_shared<vector<string>>();
            _dim = make_shared<vector<size_t>>();
            _dim->resize(vecs.size());
            size_t dim = 1;
            size_t time_pos= 0, nb_ids = 0, idx = 0;
            vector<size_t> dims;
            for(auto &vec: vecs){
                if(vec.empty()){
                    WarningOff("\n\nWARNING: Defining indices with an empty vector!\n\n");
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
                _dim->at(idx++) = vec.size();
            }
            auto vec1 = vecs.front();
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
                    if(vec->is_indexed()){
                        real_idx = vec->_ids->at(0).at(real_idx);
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
        void insert(const string& s1, Args&&... args) {
            add(s1,args...);
        }
        
        void insert(const vector<string>& all_keys) {
            add(all_keys);
        }
        
        void add_empty_row() {
            _type = matrix_;
            if (!_ids) {
                _ids = make_shared<vector<vector<size_t>>>();
            }
            _ids->resize(_ids->size()+1);
        }
        
        /** Add a key in the specified row (for sparse matrix indexing)*/
        void add_in_row(size_t row_nb, const string& key) {
            _type = matrix_;
            if (!_ids) {
                _ids = make_shared<vector<vector<size_t>>>();
            }
            _ids->resize(row_nb+1);
            auto ref_it = _keys_map->find(key);
            if(ref_it==_keys_map->end()){
                auto idx = _keys->size();
                _keys_map->insert(make_pair<>(key,idx));
                _keys->push_back(key);
                _ids->at(row_nb).push_back(idx);
            }
            else{
                _ids->at(row_nb).push_back(ref_it->second);
            }
        }
        
        template<typename... Args>
        void add(const string& s1, Args&&... args) {
            add(vector<string>({s1,args...}));
        }
        
        
        template<typename... Args>
        indices(const indices& vec1, Args&&... args) {
            list<indices> vecs;
            vecs = {vec1,forward<Args>(args)...};
            *this = indices(vecs);
        }
        
        /** Adds all new keys found in ids
         @return index set of added indices.
         */
        template<typename... Args>
        indices add(const indices& ids) {
            indices added("added");
            if(ids.is_matrix_indexed()){
                if(!is_matrix_indexed()){
                    throw invalid_argument("calling add(ids) with a matrix indexed set while current set is not matrix indexed.");
                }
                auto nb_rows = ids.get_nb_rows();
                for (size_t i = 0; i<nb_rows; i++) {
                    for (size_t j = 0; j<ids._ids->at(i).size(); j++) {
                        auto key = ids._keys->at(ids._ids->at(i).at(j));
                        auto idx = _keys->size();
                        auto pp = _keys_map->insert(make_pair<>(key,idx));
                        if (pp.second) {//new index inserted
                            _keys->push_back(key);
                            added.add(key);
                        }
                    }
                }
            }
            else if(ids.is_indexed()){
                for(auto &key_ref: ids._ids->at(0)){
                    auto key = ids._keys->at(key_ref);
                    auto idx = _keys->size();
                    auto pp = _keys_map->insert(make_pair<>(key,idx));
                    if (pp.second) {//new index inserted
                        _keys->push_back(key);
                        added.add(key);
                    }
                }
            }
            else {
                auto it = ids._keys->begin();
                for (size_t i= 0; i < ids.size(); i++) {
                    auto idx = _keys->size();
                    auto pp = _keys_map->insert(make_pair<>(*it,idx));
                    if (pp.second) {//new index inserted
                        _keys->push_back(*it);
                        added.add(*it);
                    }
                    it++;
                }
            }
            return added;
        }
        
        
        /*
         Add refs to all keys found in ids
         */
        template<typename... Args>
        void add_refs(const indices& ids) {
            if(ids.is_matrix_indexed()){
                if(!is_matrix_indexed()){
                    throw invalid_argument("calling add_refs(ids) with a matrix indexed set while current set is not matrix indexed.");
                }
                auto nb_rows = ids.get_nb_rows();
                for (size_t i = 0; i<nb_rows; i++) {
                    _ids->push_back(ids._ids->at(i));
                }
            }
            else if(ids.is_indexed()){
                for(auto &key_ref: ids._ids->at(0)){
                    add_ref(ids._keys->at(key_ref));
                }
            }
            else {
                for(auto &key: *ids._keys){
                    add_ref(key);
                }
            }
        }
        
        
        
        /** Adds a reference in _ids
         */
        void add_ref(int key_id){
            if(!_ids){
                _ids = make_shared<vector<vector<size_t>>>();
                _ids->resize(1);
            }
            _ids->at(0).push_back(key_id);
        }
        
        
        /** Adds a reference to the key specified as argument, i.e., adds the corresponding index in _ids
         @throw invalid_argument if key is not part of _keys_map
         */
        void add_ref(const string& key){
            if(!_ids){
                _ids = make_shared<vector<vector<size_t>>>();
                _ids->resize(1);
            }
            auto ref_it = _keys_map->find(key);
            if(ref_it==_keys_map->end()){
                throw invalid_argument("in indices::add_ref(string), unknown key: " + key);
            }
            _ids->at(0).push_back(ref_it->second);
        }
        
        void add(const vector<string>& all_keys) {
            for (auto &key: all_keys) {
                auto idx = _keys->size();
                auto pp = _keys_map->insert(make_pair<>(key,idx));
                if (pp.second) {//new key inserted
                    _keys->push_back(key);
                }
                else{
                    throw invalid_argument("in indices::add(string...) cannot add same key twice: " + key);
                }
            }
        }
        
        
        inline size_t get_id_inst(size_t inst = 0) const {
            if (_ids) {
                if(_ids->at(0).size() <= inst){
                    throw invalid_argument("indices::get_id_inst(size_t inst) inst is out of range");
                }
                return _ids->at(0).at(inst);
            }
            return inst;
        };
        
        
        
        
        
        bool is_indexed() const{
            return (_ids!=nullptr);
        }
        
        void remove_excluded(){
            _ids = nullptr;
            map<string,size_t> new_keys_map;
            for(auto &key_id: _excluded_keys){
                auto key = _keys->at(key_id);
                _keys_map->erase(key);
            }
            _keys->clear();
            _keys->resize(_keys_map->size());
            size_t idx = 0;
            for(auto &key_id: *_keys_map){
                _keys->at(idx) = key_id.first;
                new_keys_map[key_id.first] = idx++;
            }
            *_keys_map = new_keys_map;
            _excluded_keys.clear();
            _dim->resize(1);
            _dim->at(0) = _keys->size();
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
            res._keys_map = make_shared<map<string,size_t>>(*_keys_map);
            res._keys = make_shared<vector<string>>(*_keys);
            res._excluded_keys.insert(res._keys_map->at(key));
            //            if(!is_indexed()){
            //                res._ids = make_shared<vector<vector<size_t>>>();
            //                res._ids->resize(1);
            //            }
            //            res.reindex();
            res._name += "\\{" + key+"}";
            res.remove_excluded();
            return res;
        }
        
        size_t get_max_nb_columns() const {
            assert(is_indexed());
            if(_type != matrix_){
                throw invalid_argument("cannot call get_nb_cols() on a non-indexed set");
            }
            auto nb_inst = _ids->size();
            size_t max_dim = 0;
            for(size_t inst = 0; inst<nb_inst; inst++){
                auto nb_idx = _ids->at(inst).size();
                max_dim = max(max_dim,nb_idx);
            }
            return max_dim;
        }
        
        size_t get_nb_rows() const {
            if(_type != matrix_){
                throw invalid_argument("cannot call get_nb_rows() on a non-indexed set");
            }
            return _ids->size();
        };
        
        
        /** @return size of index set, if this is a matrix indexing, it returns the number of rows, if the set is a mask (has _ids) it returns the mask size, otherwize it return the total number of keys. **/
        size_t size() const {
            if(is_indexed()){
                if(_type == matrix_){ //Matrix-indexed
                    return get_nb_rows();
                }
                return _ids->at(0).size();
            }
            return _keys->size();
        };
        
        size_t nb_keys() const {
            if(_type != matrix_){
                throw invalid_argument("cannot call nb_keys() on a non-matrix index set");
            }
            size_t n = 0;
            for(auto &vec: *_ids){
                n+=vec.size();
            }
            return n;
        };
        
        size_t nb_active_keys() const {return _keys->size() - _excluded_keys.size();};
        
        bool empty() const {
            return _keys->size() - _excluded_keys.size() == 0;
        }
        
        string to_str() const{
            string str = "{";
            if(is_matrix_indexed()){
                auto nb_rows = get_nb_rows();
                for (size_t i = 0; i<nb_rows; i++) {
                    for (size_t j = 0; j<_ids->at(i).size(); j++) {
                        str += "(" + _keys->at(_ids->at(i).at(j)) + ")";
                        if(j<_ids->at(i).size()-1)
                            str += " ; ";
                    }
                    str += "\n";
                }
            }
            else if(is_indexed()){
                int i = 0, nb_keys = _ids->at(0).size();
                for(int idx = 0; idx < nb_keys; idx++){
                    str+= "(" + _keys->at(_ids->at(0).at(idx)) + ")";
                    if(idx<nb_keys-1)
                        str += " ; ";
                }
            }
            else {

                int i = 0, nb_keys = _keys->size();
                for(int idx = 0; idx < nb_keys; idx++){
                    if (_excluded_keys.count(i++)==0) {
                        str += "(" + _keys->at(idx) + ")";
                    }
                    if(idx<nb_keys-1)
                        str += " ; ";
                }
            }
            str += "}";
            return str;
        }
        
        void print() const{
            cout << to_str() << endl;
        }
        
        string first() const{
            return _keys->front();
        }
        
        string last() const{
            return _keys->back();
        }
    };

    /** Adds all new keys found in ids
     @return index set of added indices.
     **/
    template<typename... Args>
    indices union_ids(const indices& ids1, Args&&... args) {
        vector<indices> all_ids;
        all_ids = {ids1,forward<Args>(args)...};
        indices res("Union(");
        auto nb_entries = ids1.get_nb_entries();
        if(!ids1.is_indexed()){ //if the index set is not indexed, we assume the rest are not indexed as well
            res._keys_map = make_shared<map<string,size_t>>(*ids1._keys_map);
            res._keys = make_shared<vector<string>>(*ids1._keys);
            for (size_t i= 1; i < all_ids.size(); i++) {
                auto ids = all_ids[i];
                if(nb_entries!=ids.get_nb_entries()){
                    throw invalid_argument("union cannot be applied to index sets with different number of entries");
                }
                res.set_name(res.get_name() + ids.get_name()+",");
                auto it = ids._keys->begin();
                for (size_t i= 0; i < ids.size(); i++) {
                    auto kkey = res._keys->size();
                    auto pp = res._keys_map->insert(make_pair<>(*it,kkey));
                    if (pp.second) {//new index inserted
                        res._keys->push_back(*it);
                    }
                    it++;
                }
            }
        }
        else{  //means the ids1 is indexed so we assume all the rest are indexed as well
            //            res._ids = make_shared<vector<vector<size_t>>>(*ids1._ids);
            for (size_t i= 0; i < all_ids.size(); i++) {
                auto ids = all_ids[i];
                if(nb_entries!=ids.get_nb_entries()){
                    throw invalid_argument("union cannot be applied to index sets with different number of entries");
                }
                res.set_name(res.get_name() + ids.get_name()+",");
                auto it_ids = ids._ids->begin()->begin();
                for (size_t i= 0; i < ids.size(); i++) {
                    auto kkey = res._keys->size();
                    auto pp = res._keys_map->insert(make_pair<>(ids._keys->at(*it_ids),kkey));
                    if (pp.second) {//new index inserted
                        res._keys->push_back(ids._keys->at(*it_ids));
                    }
                    it_ids++;
                }
            }
        }
        auto name = res.get_name();
        res.set_name(name.substr(0,name.size()-1) + ")");
        return res;
    }

    indices operator-(const indices& s1, const indices& s2);

    indices intersect(const indices& s1, const indices& s2);

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

    typedef enum { ipopt, gurobi, bonmin, cplex, sdpa, _mosek, clp} SolverType;  /* Solver type */

    // settings of solvers. used by solvers like sdpa.
    typedef enum {unsolved = -1, penalty=0, fast=1, medium=2, stable=3} SolverSettings;

    typedef pair<shared_ptr<size_t>,shared_ptr<indices>> unique_id; /* A unique identifier is defined as a pair<variable id, index address> */

    template <class T>
    std::string type_name(const T& t) {
        return t._type_name;
    }

    indices range(size_t i, size_t j);


    /** Combine each ith key of each index set with the other ith keys in the other index sets
     \warning all index sets must have the same size
     */
    template<typename... Args>
    indices combine(const indices& ids1, Args&&... args){
        bool matrix_indexed = false;
        indices res;
        string res_name;
        list<indices> all_ids = {ids1,forward<Args>(args)...};
        for (auto &ids: all_ids) {
            if(ids.is_matrix_indexed()){
                matrix_indexed = true;
            }
            res_name += ids.get_name()+",";
        }
        res_name = res_name.substr(0, res_name.size()-1);/* remove last comma */
        res.set_name(res_name);
        if(matrix_indexed){
            auto nb_rows = ids1.get_nb_rows();
            res._ids = make_shared<vector<vector<size_t>>>();
            res._ids->resize(nb_rows);
            for (auto &ids: all_ids) {
                if(!ids.is_matrix_indexed())
                    throw invalid_argument("In combine(ids..) all or none of the index sets should be matrix indexed");
                if(ids.get_nb_rows()!=nb_rows)
                    throw invalid_argument("In combine(ids..) all indices should have the same number of rows");
            }
            for (size_t i = 0; i<nb_rows; i++) {
                for (size_t j = 0; j< ids1._ids->at(i).size(); j++) {
                    string combined_key, key = "";
                    for (auto &ids: all_ids) {
                        if(ids._ids->at(i).size()!=ids1._ids->at(i).size())
                            throw invalid_argument("In combine(ids..) all indices should have the same number of entries per row");
                        auto idx = ids._ids->at(i).at(j);
                        key = ids._keys->at(idx);
                        combined_key += key +",";
                    }
                    combined_key = combined_key.substr(0, combined_key.size()-1);/* remove last comma */
                    res.add_in_row(i,combined_key);
                }
            }
        }
        else {
            auto nb_keys = ids1.size();
            for (auto i = 0; i< nb_keys; i++) {
                string combined_key, key = "", prev_key = "";
                //            bool same_key = true;
                for (auto &ids: all_ids) {
                    auto idx = ids.get_id_inst(i);
                    key = ids._keys->at(idx);
                    //                if(prev_key!="" && key!=prev_key){
                    //                    same_key = false;
                    //                }
                    //                prev_key = key;
                    combined_key += key +",";
                }
                //            if(same_key){
                //                combined_key = key;
                //            }
                //            else {
                combined_key = combined_key.substr(0, combined_key.size()-1);/* remove last comma */
                //            }
                res.add(combined_key);
            }
        }
        return res;
    }


    }



    #endif
