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
#include <Gravity/constant.h>

using namespace std;


/** Backbone class for parameter */
class param_: public constant_, public expr{
protected:
    string       _name;
    ParamType   _intype;
public:
    
    virtual ~param_(){};
    
    string get_name() const { return _name;};
    
    ParamType get_intype() const { return _intype;}
    
    void set_type(ParamType type){ _intype = type;}
    
    /** Querries */
    bool is_binary() const{
        return (_intype==binary_p);
    };
    
    bool is_integer() const{
        return (_intype==integer_p);
    };
    
    bool is_float() const{
        return (_intype==float_p);
    };
    
    bool is_double() const{
        return (_intype==double_p);
    };
    
    bool is_long() const{
        return (_intype==long_p);
    };
    
    
    /** Operators */
    bool operator==(const param_& p) const {
        return (_type==p._type && _intype==p._intype && _name==p._name);
    }
};




/** A parameter can be a bool, a short, an int, a float or a double*/
template<typename type = int>
class param: public param_{
protected:
    shared_ptr<vector<type>>                _val;    

public:
    
    param(){
        _type = parameter;
        _name = "noname";
        throw invalid_argument("Please enter a name in the parameter constructor");
    }
    
    ~param(){}
    

    param (const param& p) {
        _type = parameter;
        _intype = p._intype;
        _val = p._val;
        _name = p._name;
    }
    
    void update_type() {
        _type = parameter;
        if(typeid(type)==typeid(bool)){
            _intype = binary_p;
            return;
        }
        if(typeid(type)==typeid(short)) {
            _intype = short_p;
            return;
        }
        if(typeid(type)==typeid(int)) {
            _intype = integer_p;
            return;
        }
        if(typeid(type)==typeid(float)) {
            _intype = float_p;
            return;
        }
        if(typeid(type)==typeid(double)) {
            _intype = double_p;
            return;
        }
        if(typeid(type)==typeid(long double)) {
            _intype = long_p;
            return;
        }
        throw bad_alloc();
    }
    
    
    param(const char* s){
        _name = s;
        update_type();
        _val = make_shared<vector<type>>();
    }

    ParamType get_intype() const { return _intype;}
    
//    void set_intype(ParamType type){ _intype = type;}
    
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
    
    void set_val(type val){_val->at(0) = val;}
    
    void set_val(int i, type val){
        if (_val->size() <= i) {
            throw out_of_range("set_val(int i, type val)");
        }
        _val->at(i) = val;
    }

    /** Operators */
    bool operator==(const param& p) const {
        return (_name==p._name && _type==p._type && _intype==p._intype && *_val==*p._val);
    }
    
    param& operator=(type v){
        _val->push_back(v);
        return *this;
    }
    
//    param& operator=(int v){
//        if (_type<integer_p) {
//            throw invalid_argument("Cannot assign an integer value in this parameter, check original type.");
//        }
//        _val = v;
//        return *this;
//    }
//    
//    
//    param& operator=(double v){
//        if (_type<double_p){
//            throw invalid_argument("Cannot assign a double value in this parameter, check original type.");
//        }
//        _val = v;
//        return *this;
//    }
//    
//    param& operator=(long double v){
//        if (_type<long_p){
//            throw invalid_argument("Cannot assign a long double value in this parameter, check original type.");
//        }
//        _val = v;
//        return *this;
//    }

    /** Output */
    void print(bool vals=false) const{
        cout << _name;
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
