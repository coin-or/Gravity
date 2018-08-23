//
//  var.h
//  Gravity
//
//  Created by Hijazi, Hassan on 21/10/16.
//
//

#ifndef var_h
#define var_h

#include <gravity/param.h>
#include <gravity/Net.h>
#include <stdio.h>
#include <string>
#include <set>
#include <list>
#include <limits>
#include <random>



using namespace std;

namespace gravity {
    
    
class func_;
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
        /* TODO */
    };
    
/** Backbone class for parameter */
class var_ {
public:
    virtual ~var_() {};    
};

/** A variable can be a bool, a short, an int, a float or a double*/
template<typename type = double>
// define variable as a parameter with bounds
class var: public param<type>, public var_ {

public:
    shared_ptr<func_>   _lb; /**< Lower Bound */
    shared_ptr<func_>   _ub; /**< Upper Bound */
    bool _in_q_cone = false;
    bool _psd = false;
    
    /* Constructors */
    //@{
    /** Unbounded variable constructor */
    var();
    ~var() {};
    var(const string& name);
    var(const string& name, Sign s);
    var(const var<type>& v);
    var(var<type>&& v);
    //@}

    //@{
    /** Bounded variable constructor */
    var(const string& name, type lb, type ub);
//    var(const string& name, const param<type>& lb, const param<type>& ub);//TODO move version of bounds
    var(const string& name, const func_& lb, const func_& ub);
    var(const string& name, func_&& lb, func_&& ub);
    var(const string& name, const param<type>& sb);// Constructor with symmetric bound: [-sb, sb]
//    var(const string& name, func_&& sb);// Constructor with symmetric bound: [-sb, sb]
    //@}


    // Retrieve specified indexed variable.
    template<typename... Args>  
    var operator()(size_t t1, Args&&... args) {
        var<type> res(this->_name);
        res.param<type>::operator=(param<type>::operator()(t1, args...));
        res.param<type>::set_type(var_c);
        res._lb = this->_lb;
        res._ub = this->_ub;
        return res;
    }

    template<typename... Args>
    var operator()(string t1, Args&&... args) {
        var<type> res(this->_name);
        res.param<type>::operator=(param<type>::operator()(t1, args...));
        res.param<type>::set_type(var_c);
        res._lb = this->_lb;
        res._ub = this->_ub;
        return res;
    }
    
    
    


    
    var in_pairs(){
        var<type> res(this->_name);
        res.param<type>::operator=(param<type>::in_pairs());
        res.param<type>::set_type(var_c);
        res._lb = this->_lb;
        res._ub = this->_ub;
        return res;
    }
    
    var from(){
        var<type> res(this->_name);
        res.param<type>::operator=(param<type>::from());
        res.param<type>::set_type(var_c);
        res._lb = this->_lb;
        res._ub = this->_ub;
        return res;
    }
    
    var prev(){
        var<type> res(this->_name);
        res.param<type>::operator=(param<type>::prev());
        res.param<type>::set_type(var_c);
        res._lb = this->_lb;
        res._ub = this->_ub;
        return res;
    }
    
    var min_time(){
        var<type> res(this->_name);
        res.param<type>::operator=(param<type>::min_time());
        res.param<type>::set_type(var_c);
        res._lb = this->_lb;
        res._ub = this->_ub;
        return res;
    }
    
    
    var out_arcs(){
        var<type> res(this->_name);
        res.param<type>::operator=(param<type>::out_arcs());
        res.param<type>::set_type(var_c);
        res._lb = this->_lb;
        res._ub = this->_ub;
        return res;
    }
    
    var in_arcs(){
        var<type> res(this->_name);
        res.param<type>::operator=(param<type>::in_arcs());
        res.param<type>::set_type(var_c);
        res._lb = this->_lb;
        res._ub = this->_ub;
        return res;
    }
    
    var in_gens(){
        var<type> res(this->_name);
        res.param<type>::operator=(param<type>::in_gens());
        res.param<type>::set_type(var_c);
        res._lb = this->_lb;
        res._ub = this->_ub;
        return res;
    }
    
    var in_pot_gens(){
        var<type> res(this->_name);
        res.param<type>::operator=(param<type>::in_pot_gens());
        res.param<type>::set_type(var_c);
        res._lb = this->_lb;
        res._ub = this->_ub;
        return res;
    }
    
    var in_bats(){
        var<type> res(this->_name);
        res.param<type>::operator=(param<type>::in_bats());
        res.param<type>::set_type(var_c);
        res._lb = this->_lb;
        res._ub = this->_ub;
        return res;
    }
    
    var in_wind(){
        var<type> res(this->_name);
        res.param<type>::operator=(param<type>::in_wind());
        res.param<type>::set_type(var_c);
        res._lb = this->_lb;
        res._ub = this->_ub;
        return res;
    }
    
    var in_pv(){
        var<type> res(this->_name);
        res.param<type>::operator=(param<type>::in_pv());
        res.param<type>::set_type(var_c);
        res._lb = this->_lb;
        res._ub = this->_ub;
        return res;
    }
    
    var to(){
        var<type> res(this->_name);
        res.param<type>::operator=(param<type>::to());
        res.param<type>::set_type(var_c);
        res._lb = this->_lb;
        res._ub = this->_ub;
        return res;
    }
    
    template<typename Tobj>
    var from(const vector<Tobj*>& vec){
        var<type> res(this->_name);
        res.param<type>::operator=(param<type>::from(vec));
        res.param<type>::set_type(var_c);
        res._lb = this->_lb;
        res._ub = this->_ub;
        return res;
    }
    
    template<typename Tobj>
    var from(const vector<Tobj>& vec){
        var<type> res(this->_name);
        res.param<type>::operator=(param<type>::from(vec));
        res.param<type>::set_type(var_c);
        res._lb = this->_lb;
        res._ub = this->_ub;
        return res;
    }
    
    template<typename Tobj>
    var from(const vector<Tobj*>& vec, const indices& T){
        var<type> res(this->_name);
        res.param<type>::operator=(param<type>::from(vec, T));
        res.param<type>::set_type(var_c);
        res._lb = this->_lb;
        res._ub = this->_ub;
        return res;
    }
    
    template<typename Tobj>
    var to(const vector<Tobj*>& vec, const indices& T){
        var<type> res(this->_name);
        res.param<type>::operator=(param<type>::to(vec, T));
        res.param<type>::set_type(var_c);
        res._lb = this->_lb;
        res._ub = this->_ub;
        return res;
    }
    

    template<typename Tobj>
    var from(const vector<Tobj*>& vec, unsigned T){
        var<type> res(this->_name);
        res.param<type>::operator=(param<type>::from(vec, T));
        res.param<type>::set_type(var_c);
        res._lb = this->_lb;
        res._ub = this->_ub;
        return res;
    }
    
    template<typename Tobj>
    var to(const vector<Tobj*>& vec){
        var<type> res(this->_name);
        res.param<type>::operator=(param<type>::to(vec));
        res.param<type>::set_type(var_c);
        res._lb = this->_lb;
        res._ub = this->_ub;
        return res;
    }
    
    template<typename Tobj>
    var to(const vector<Tobj>& vec){
        var<type> res(this->_name);
        res.param<type>::operator=(param<type>::to(vec));
        res.param<type>::set_type(var_c);
        res._lb = this->_lb;
        res._ub = this->_ub;
        return res;
    }

    template<typename Tobj>
    var to(const vector<Tobj*>& vec, unsigned T){
        var<type> res(this->_name);
        res.param<type>::operator=(param<type>::to(vec, T));
        res.param<type>::set_type(var_c);
        res._lb = this->_lb;
        res._ub = this->_ub;
        return res;
    }

    template<typename... Args>
    var in(const indices& vec1, Args&&... args) {
        var<type> res(this->_name);
        res._lb = this->_lb;
        res._ub = this->_ub;
        if(get<1>(param_::_unique_id)==unindexed_){
            if(this->_rev_indices->size()==0 && (!res._ub->is_number() || !res._lb->is_number())){
                auto ids = indices(vec1,args...);
                if(!res._lb->is_number()){
                    (res._lb->in(ids));
                }
                if(!res._ub->is_number()){
                    (res._ub->in(ids));
                }
            }
        }
        res.param<type>::operator=(param<type>::in(vec1, forward<Args>(args)...));
        res.param<type>::set_type(var_c);
        res._is_relaxed = param_::_is_relaxed;
        return res;
    }
    
    template<typename... Args>
    var prev(const indices& vec1, Args&&... args) {
        var<type> res(this->_name);
        res.param<type>::operator=(param<type>::prev(vec1, forward<Args>(args)...));
        res.param<type>::set_type(var_c);
        res._lb = this->_lb;
        res._ub = this->_ub;
        //        if(!this->_lb->is_number()){
        //            *res._lb = this->_lb->in(vec);
        //        }
        //        if(!this->_ub->is_number()){
        //            *res._ub = this->_ub->in(vec);
        //        }
        return res;
    }
    
    
    template<typename Tobj> var min_time(const vector<Tobj*>& vec, const indices& ids, param<int> time) {
        var<type> res(this->_name);
        res.param<type>::operator=(param<type>::min_time(vec, ids, time));
        res.param<type>::set_type(var_c);
        //        if(!this->_lb->is_number()){
        //            *res._lb = this->_lb->in(vec);
        //        }
        //        if(!this->_ub->is_number()){
        //            *res._ub = this->_ub->in(vec);
        //        }
        return res;
    }
    
    
    template<typename Tobj>
    var in(const vector<Tobj*>& vec) {
        var<type> res(this->_name);
        res.param<type>::operator=(param<type>::in(vec));
        res.param<type>::set_type(var_c);
        if(!this->_lb->is_number()){
            this->_lb->in(vec);
            res._lb = this->_lb;
        }
        else {
            res._lb = _lb;
        }
        if(!this->_ub->is_number()){
            this->_ub->in(vec);
            res._ub = this->_ub;
        }
        else {
            res._ub = _ub;
        }
        return res;
    }
    
    template<typename Tobj>
    var in(const vector<Tobj>& vec) {
        var<type> res(this->_name);
        res.param<type>::operator=(param<type>::in(vec));
        res.param<type>::set_type(var_c);
        if(!this->_lb->is_number()){
            *res._lb = this->_lb->in(vec);
        }
        if(!this->_ub->is_number()){
            *res._ub = this->_ub->in(vec);
        }
        return res;
    }
    
    var in(const node_pairs& np){
        return this->in(np._keys);
    }
    
    
    
    
    var in_arcs(const vector<Node*>& vec) {
        var<type> res(this->_name);
        res.param<type>::operator=(param<type>::in_arcs(vec));
        res.param<type>::set_type(var_c);
        res._lb = this->_lb;
        res._ub = this->_ub;
        return res;
    }
    
    var in_arcs(const vector<Node*>& vec, const indices& T) {
        var<type> res(this->_name);
        res.param<type>::operator=(param<type>::in_arcs(vec,T));
        res.param<type>::set_type(var_c);
        res._lb = this->_lb;
        res._ub = this->_ub;
        return res;
    }
    
    var out_arcs(const vector<Node*>& vec) {
        var<type> res(this->_name);
        res.param<type>::operator=(param<type>::out_arcs(vec));
        res.param<type>::set_type(var_c);
        res._lb = this->_lb;
        res._ub = this->_ub;
        return res;
    }
    
    var out_arcs(const vector<Node*>& vec, const indices& T) {
        var<type> res(this->_name);
        res.param<type>::operator=(param<type>::out_arcs(vec,T));
        res.param<type>::set_type(var_c);
        res._lb = this->_lb;
        res._ub = this->_ub;
        return res;
    }
    
    
    var in_gens(const vector<Node*>& vec) {
        var<type> res(this->_name);
        res.param<type>::operator=(param<type>::in_gens(vec));
        res.param<type>::set_type(var_c);
        res._lb = this->_lb;
        res._ub = this->_ub;
        return res;
    }
    
    var in_gens(const vector<Node*>& vec, const indices& T) {
        var<type> res(this->_name);
        res.param<type>::operator=(param<type>::in_gens(vec,T));
        res.param<type>::set_type(var_c);
        res._lb = this->_lb;
        res._ub = this->_ub;
        return res;
    }
    
    var in_bats(const vector<Node*>& vec, const indices& T) {
        var<type> res(this->_name);
        res.param<type>::operator=(param<type>::in_bats(vec,T));
        res.param<type>::set_type(var_c);
        res._lb = this->_lb;
        res._ub = this->_ub;
        return res;
    }
    
    var in_wind(const vector<Node*>& vec, const indices& T) {
        var<type> res(this->_name);
        res.param<type>::operator=(param<type>::in_wind(vec,T));
        res.param<type>::set_type(var_c);
        res._lb = this->_lb;
        res._ub = this->_ub;
        return res;
    }
    
    var in_pv(const vector<Node*>& vec, const indices& T) {
        var<type> res(this->_name);
        res.param<type>::operator=(param<type>::in_pv(vec,T));
        res.param<type>::set_type(var_c);
        res._lb = this->_lb;
        res._ub = this->_ub;
        return res;
    }
    
    template<typename Tobj>
    var in_pairs(const vector<Tobj*>& vec) {
        var<type> res(this->_name);
        res.param<type>::operator=(param<type>::in_pairs(vec));
        res.param<type>::set_type(var_c);
        res._lb = this->_lb;
        res._ub = this->_ub;
        return res;
    }
    
    
    template<typename Tobj>
    var in_pairs(const vector<Tobj*>& vec, const indices& ids) {
        var<type> res(this->_name);
        res.param<type>::operator=(param<type>::in_pairs(vec,ids));
        res.param<type>::set_type(var_c);
        res._lb = this->_lb;
        res._ub = this->_ub;
        return res;
    }
    
    template<typename Tobj>
    var in_pairs(const vector<Tobj>& vec) {
        var<type> res(this->_name);
        res.param<type>::operator=(param<type>::in_pairs(vec));
        res.param<type>::set_type(var_c);
        res._lb = this->_lb;
        res._ub = this->_ub;
        return res;
    }

    template<typename Tobj>
    var in_pairs(const vector<Tobj*>& vec, unsigned T) {
        var<type> res(this->_name);
        res.param<type>::operator=(param<type>::in_pairs(vec, T));
        res.param<type>::set_type(var_c);
        res._lb = this->_lb;
        res._ub = this->_ub;
        return res;
    }


    template<typename Tobj>
    var in(const vector<Tobj>& vec, unsigned T) {
        var<type> res(this->_name);
        res.param<type>::operator=(param<type>::in(vec, T));
        res.param<type>::set_type(var_c);
        res._lb = this->_lb;
        res._ub = this->_ub;
        return res;
    }

    template<typename Tobj>
    var in(const Tobj nm, unsigned t) {
        var<type> res(this->_name);
        res.param<type>::operator=(param<type>::in(nm, t));
        res.param<type>::set_type(var_c);
        res._lb = this->_lb;
        res._ub = this->_ub;
        return res;
    }

    template<typename Tobj>
    var in_at(const vector<Tobj>& vec, unsigned t) {
        var<type> res(this->_name);
        res.param<type>::operator=(param<type>::in_at(vec, t));
        res.param<type>::set_type(var_c);
        res._lb = this->_lb;
        res._ub = this->_ub;
        return res;
    }

    var excl(unsigned index) {
        var<type> res(this->_name);
        res.param<type>::operator=(param<type>::excl(index));
        res.param<type>::set_type(var_c);
        res._lb = this->_lb;
        res._ub = this->_ub;
        return res;
    }
    
    template<typename Tobj>
    var submat(const vector<Tobj>& vec) {
        var<type> res(this->_name);
        res.param<type>::operator=(param<type>::submat(vec));
        res.param<type>::set_type(var_c);
        if(!this->_lb->is_number()){
            *res._lb = this->_lb->in(vec);
        }
        if(!this->_ub->is_number()){
            *res._ub = this->_ub->in(vec);
        }
        return res;
    }
    
    var from(const ordered_pairs& pairs);
    var to(const ordered_pairs& pairs);
    var in(const ordered_pairs& pairs);
    vector<var> in_bags(const std::vector<std::vector<Node*>>& bags, unsigned size);
    vector<var> pairs_in(const std::vector<std::vector<Node*>>& bags, unsigned size);
    vector<var> pairs_in_directed(Net& net, const std::vector<std::vector<Node*>>& bags, unsigned size);

    /* Querries */

    type    get_lb(size_t i = 0) const;

    type    get_ub(size_t i = 0) const;


    bool is_bounded_above(int i = 0) const;

    bool is_bounded_below(int i = 0) const;

    bool is_constant(int i=0) const;

    Sign get_sign(int idx = 0) const;

    /* Modifiers */
    void    set_size(vector<size_t>);
    void    set_size(size_t s, type val = 0);

    void    add_bounds(type lb, type ub);
    void    add_lb_only(type v); /**< Adds a new lower bound and infinity in the corresponding upper bound*/
    void    add_ub_only(type v); /**< Adds a new upper bound and -infinity in the corresponding lower bound*/

    void    set_lb(int i, type v);
    void    set_ub(int i, type v);

    void in_q_cone(); /**< States that the variable is contained in a quadratic cone, for mosek */

    //void    set_lb(string name, type v); // change lb and ub via names.
    //void    set_ub(string name, type v);

    /* Operators */
    var& operator=(const var& v);
    var& operator=(var&& v);

    var& operator=(type v) {
        param<type>::_val->push_back(v);
        param<type>::update_range(v);
        param<type>::_dim[0]++;
        return *this;
    }
    
    bool operator==(const var& v) const;
    bool operator!=(const var& v) const;
    
    var& in(const space& s){
//        if(s._dim.size()>1){
//            throw invalid_argument("2D spaces unsupported yet");
//        }
        set_size(s._dim);
        param_::_rev_indices->resize(s._dim[0]);
        for(unsigned i = 0 ; i< s._dim[0]; i++){
            auto key = to_string(i);
            param_::_indices->insert(make_pair<>(key,i));
            (*param_::_rev_indices)[i] = key;
        }
        return *this;
    }
    
    void initialize_uniform() {
        std::default_random_engine generator;
        for (int i = 0; i<param<type>::_val->size(); i++) {
            std::uniform_real_distribution<double> distribution(get_lb(i),get_ub(i));
            param<type>::_val->at(i) = distribution(generator);
        }
    }
    
//    var& operator^(size_t d) {
//        set_size(d);
//        return *this;
//    }

    var tr() const {
        auto v = var(*this);
        v._is_transposed = true;
        v._is_vector = true;
        return v;
    }
    
    var vec() const {
        auto v = var(*this);
        v._is_vector = true;
        return v;
    }

    /* Output */
    void print(bool bounds=false) const;

};

template<typename type>
var<type> all(const var<type>& p) {
    auto pp = var<type>(p);
    pp._is_vector = true;
    return pp;
}


// In contrast to a general variable, indices of an SDP variable should be
// recorded using _sdpindices, moreover the size of an SDP variable is d (d x d).
template<typename type = double>
//class sdpvar: public param<type>, public var_{
class sdpvar: public var<type> {

public:
    size_t _symdim=0;
    //@{
    /** Unbounded sdp-variable constructor */
    sdpvar();
    ~sdpvar() {};

    sdpvar(const string& name);
    sdpvar(const sdpvar<type>& v);
    sdpvar(sdpvar<type>&& v);
    //@}

    /** bounded sdp-variable constructor */
    sdpvar(const string& name, type lb, type ub);

    template<typename... Args>
    sdpvar operator()(size_t t1, Args&&... args) {
        sdpvar res(this->_name);
        res._id = this->_id;
        res._vec_id = this->_vec_id;
        res._intype = this->_intype;
        res._val = this->_val;
        res._range = this->_range;
        res._lb = this->_lb;
        res._ub = this->_ub;
        list<size_t> indices;
        indices = {forward<size_t>(args)...};
        indices.push_front(t1);
        string key;
        auto it = indices.begin();

        for (size_t i= 0; i<indices.size(); i++) {
            key += to_string(*it);
            if (i<indices.size()-1) {
                key += ",";
            }
            it++;
        }

        auto it2 = param_::_sdpindices->find(key);
        if (it2 == param_::_sdpindices->end()) {
            auto temp = make_pair<>(t1, (*(++indices.begin())));
            res._sdpindices->insert(make_pair<>(key,temp));
            param_::_sdpindices->insert(make_pair<>(key,temp));
//            res._dim[0] = 1;
//            res._symdim = 1;
        }
        else {
            auto temp = param_::_sdpindices->at(key);
            res._sdpindices->insert(make_pair<>(key,temp));
            res._dim[0] = 1;
            res._symdim = 1;
        }
        res._name += "["+key+"]";
        res._is_indexed = true;
        return res;
    }

    /* Modifiers */
    
    
    
    
    
    //void    set_size(size_t s, type val = 0);
    /* Operators */
    sdpvar& operator = (type v) {
        param<type>::_val->push_back(v);
        param<type>::_dim[0]++;
        return *this;
    }

    bool operator == (const sdpvar& v) const;
    bool operator >= (const sdpvar& v) const;
    bool operator != (const sdpvar& v) const;
    sdpvar& operator^(size_t d) {
        // the upper/lower triangular part.
        param<type>::set_size(d*(d+1)/2);
        _symdim = d;
        return *this;
    }

    sdpvar tr() const {
        auto v = sdpvar(*this);
        v._is_transposed = true;
        return v;
    }
    /* Output */
    void print(bool vals = false) const;
};
}
#endif /* sdpvar_h */
