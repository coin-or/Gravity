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
    string          _name;
    NType           _intype;
    vector<int>     _indices;
public:
    
    virtual ~param_(){};
    
    string get_name() const;
    NType get_intype() const { return _intype;}
    
    vector<int> get_indices() const {
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
    
    
    /** Operators */
    bool operator==(const param_& p) const {
        return (_type==p._type && _intype==p._intype && _name==p._name);
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
    
    
    param(){
        _type = par_c;
        _name = "noname";
        throw invalid_argument("Please enter a name in the parameter constructor");
    }
    
    ~param(){}
    

    param (const param& p) {
        _type = par_c;
        _intype = p._intype;
        _val = p._val;
        _name = p._name;
        _indices = p._indices;
    }

//    template<class... Args>
//    param operator()(Args&&... args){
//        auto res(*this);
////        va_list arg_list;
//        va_start(args, sizeof...(args));
//        for (int i = 0; i < sizeof...(args); ++i) {
//            res._indices.insert(va_arg(args, int));
//        }
//        va_end(args);
//        return res;
//    }

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
//        _indices.reserve(3);
    }

    NType get_intype() const { return _intype;}
    
//    void set_intype(NType type){ _intype = type;}
    
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
    
    void add_val(type val){_val->push_back(val);}
    
    void set_val(int i, type val){
        if (_val->size() <= i) {
            throw out_of_range("set_val(int i, type val)");
        }
        _val->at(i) = val;
    }
    void set_val(type val){
        for (auto &v: _val) {
            v = val;
        }        
    }

    /** Operators */
    bool operator==(const param& p) const {
        return (_name==p._name && _type==p._type && _intype==p._intype && *_val==*p._val);
    }
    
    param& operator=(type v){
        _val->push_back(v);
        return *this;
    }
    

    /** Output */
    void print(bool vals=false) const{
        cout << get_name();
        if(vals){
            cout << " = [ ";
            for(auto v: *_val){
                cout << v << " ";
            }
            cout << "];\n";
        }
    }


};
#endif /* defined(____param__) */
