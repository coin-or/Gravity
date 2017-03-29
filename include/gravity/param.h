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
#include <set>
#include <Gravity/constant.h>

using namespace std;


/** Backbone class for parameter */
class param_: public constant_{
protected:
    string              _name;
    NType               _intype;
    map<string,int>*    _indices; /*<< A map storing all the indices this parameter has, the key is represented by a string, while the entry indicates the right position in the values and bounds vectors */
public:
    
    virtual ~param_(){};
    
    string get_name(bool indices=true) const;
    NType get_intype() const { return _intype;}
    
    map<string,int>* get_indices() const {
        return _indices;
    }
    
    void set_type(NType type){ _intype = type;}
    
    /** Querries */
    bool is_binary() const{
        return (_intype==binary_);
    };
    
    bool is_integer() const{
        return (_intype==integer_);
    };
    
    bool is_float() const{
        return (_intype==float_);
    };
    
    bool is_double() const{
        return (_intype==double_);
    };
    
    bool is_long() const{
        return (_intype==long_);
    };
    
    Sign get_all_sign() const; /**< If all instances of the current parameter/variable have the same sign, it returns it, otherwise, it returns unknown. **/
    Sign get_sign(int idx = 0) const; /**< returns the sign of one instance of the current parameter/variable. **/
    
    /** Operators */
    bool operator==(const param_& p) const {
        return (_type==p._type && _intype==p._intype && get_name()==p.get_name());
    }
};


///** A pair <indices, param> */
//class ind_param: public constant_{
//    
//public:
//    set<ind>*               _indices;
//    param_*                 _p;
//    
//    ind_param(ind_param&& p){
//        _indices = p._indices;
//        p._indices = nullptr;
//        _p = p._p;
//    };
//    
////    ind_param(const param_& p){
////        _indices = new set<ind>();
////        _p = (param_*)copy((constant_*)&p);
////    };
//    
//    ind_param(param_* p){
//        _indices = new set<ind>();
//        _p = p;
//    };
//    
//    void add_index(ind i);
//    
//    bool has_index(ind i) const;
//    
//    ~ind_param(){
//        delete _indices;
//        delete _p;
//    };
//    
//    bool operator==(const ind_param& p) const;
//    
//};




/** A parameter can be a bool, a short, an int, a float or a double*/
template<typename type = int>
class param: public param_{
protected:
    shared_ptr<vector<type>>                _val;
    

public:
    pair<type,type>                         _range; /**< (Min,Max) values in vals **/
    
    param(){
        _type = par_c;
        _name = "noname";
        throw invalid_argument("Please enter a name in the parameter constructor");
    }
    
    ~param(){
        delete _indices;
    }
    

    param (const param& p) {
        _type = par_c;
        _intype = p._intype;
        _val = p._val;
        _name = p._name;
        _indices = new map<string, int>(*p._indices);
        _range = p._range;
    }
    
    param (param&& p) {
        _type = par_c;
        _intype = p._intype;
        _val = p._val;
        _name = p._name;
        _indices = p._indices;
        p._indices = nullptr;
        _range = p._range;
    }


    void set_type(CType t) { _type = t;}
    
    void set_intype(NType t) { _intype = t;}
    
    void update_type() {
        _type = par_c;
        if(typeid(type)==typeid(bool)){
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
    
    
    param(const char* s){
        _name = s;
        update_type();
        _val = make_shared<vector<type>>();
        _indices = new map<string,int>();
        _range.first = numeric_limits<type>::max();
        _range.second = numeric_limits<type>::lowest();
    }

    NType get_intype() const { return _intype;}
    
    type eval() const{
        if (_val->size() == 0) {
            throw "No values stored!";
        }
        return _val->at(_val->size()-1);
    }
    
    type eval(int i) const{
        if (_val->size() <= i) {
            throw out_of_range("get_val(int i)");
        }
        return _val->at(i);
    }
    
    
    /* Modifiers */
    void    set_size(int s){
        _val->reserve(s);
    };
    
    void add_val(type val){
        _val->push_back(val);
        update_range(val);
    }
    
    void update_range(type val){
        if (val < _range.first) {
            _range.first = val;
        }
        if (val > _range.second) {
            _range.second = val;
        }
    }
    
    void set_val(int i, type val){
        if (_val->size() <= i) {
            throw out_of_range("set_val(int i, type val)");
        }
        _val->at(i) = val;
        update_range(val);
    }
    
    void set_val(type val){
        for (auto &v: _val) {
            v = val;
        }
        _range.first = val;
        _range.second = val;
    }
    
    Sign get_sign(int idx = 0) const{
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

    
    
    Sign get_all_sign() const{
        if (_range.first == 0 && _range.second == 0) {
            return zero_;
        }
        if (_range.second < 0  && _range.first < 0) {
            return neg_;
        }
        if (_range.first > 0 && _range.second > 0) {
            return pos_;
        }
        if (_range.second == 0   && _range.first < 0){
            return non_pos_;
        }
        if (_range.first == 0  && _range.second > 0) {
            return non_neg_;
        }
        return unknown_;
    }
    
    bool is_unit() const{ /**< Returns true if all values of this paramter are 1 **/
        return (_range.first == 1 && _range.second == 1);
    }

    bool is_zero() const{ /**< Returns true if all values of this paramter are 0 **/
        return (_range.first == 0 && _range.second == 0);
    }

    bool is_non_positive() const{ /**< Returns true if all values of this paramter are <= 0 **/
        return (_range.second <= 0   && _range.first <= 0);
    }
    
    bool is_positive() const{ /**< Returns true if all values of this paramter are positive **/
        return (_range.first > 0 && _range.second > 0);
    }
    
    bool is_non_negative() const{ /**< Returns true if all values of this paramter are >= 0 **/
        return (_range.first >= 0  && _range.second >= 0);
    }
    
    bool is_negative() const{ /**< Returns true if all values of this paramter are positive **/
        return (_range.second < 0  && _range.first < 0);
    }

    /** Operators */
    bool operator==(const param& p) const {
        return (get_name()==p.get_name() && _type==p._type && _intype==p._intype && *_indices==*p._indices && *_val==*p._val);
    }
    
    param& operator=(type v){
        _val->push_back(v);
        update_range(v);
        return *this;
    }
    

    /** Output */
    void print(bool vals=false) const{
        cout << to_string(vals);
    }
    
    string to_string(bool vals=false) const{
        string str = get_name();
        if(vals){
            str += " = [ ";
            for(auto v: *_val){
                str += v;
                str += " ";
            }
            str += "];\n";
        }
        return str;
    }
    


};
#endif /* defined(____param__) */
