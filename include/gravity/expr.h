//
//  func.h
//  Gravity
//
//  Created by Hijazi, Hassan on 25 Oct 18.
//
//

#ifndef expr_h
#define expr_h

#include <gravity/poly.h>
#include <stdio.h>
#include <map>
#include <iterator>
#include <queue>
#include <list>
#include <limits>
#include <set>

using namespace std;

string operator_str(gravity::OperatorType ot);

namespace gravity {
    class func_;
    /** Backbone class for unary and binary expressions. */
    template<typename type = double>
    class expr: public constant_{
    protected:
    public:
        type                                   _coef = unit<type>().eval(); /**< coefficient multpying the expression */
        Convexity                              _all_convexity = linear_; /**< If all instances of this expression have the same convexity type, it stores it here, i.e. linear, convex, concave, otherwise it stores unknown. >>**/
        Sign                                   _all_sign = zero_; /**< If all instances of this expression have the same sign, it stores it here, otherwise it stores unknown. >>**/

        shared_ptr<pair<type,type>>            _range = make_shared<pair<type,type>>(); /**< (Min,Max) values of expression **/
        string                                 _to_str = "noname"; /**< A string representation of the expression */
        
        virtual void in(const indices& ids){};
        void propagate_dim(size_t d){
            if(_is_transposed){
                _dim[1] = d;
            }
            else {
                _dim[0] = d;
            }
        }
        
        void reverse_sign(){ _coef *= -1.; };
    };


    /** Class uexpr (unary expression), stores a unary expression tree. */
    template<typename type = double>
    class uexpr: public expr<type>{
        
    public:
        OperatorType                    _otype = id_;
        shared_ptr<constant_>           _son = nullptr;
        
        
        
        void reset(){
            _son = nullptr;
            _otype = id_;
            this->_to_str = "noname";
            this->_coef = 1.;
        };
        
        
        OperatorType get_otype() const{
            return _otype;
        };
        
        /** Operators */
        
        shared_ptr<constant_> copy()const{return make_shared<uexpr>(*this);};
        
        
        bool operator!=(const uexpr& c) const{
            return !(*this==c);
        };
        
        void propagate_dim(size_t d){
            if(this->_is_transposed){
                this->_dim[1] = d;
            }
            else {
                this->_dim[0] = d;
            }
            _son->propagate_dim(d);
        }
        
        void uneval(){
            _son->uneval();
        }
        
        /** allocates memory for current and all sub-functions */
        void allocate_mem(){
            _son->allocate_mem();
        };
        
        void in(const indices& ids){
            if(_son->is_function()){
                auto f = static_pointer_cast<func<type>>(_son);
                f->in(ids);
            }
        };
        
        Sign get_all_sign() const{
            return unknown_;// TO UPDATE
        }
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) < sizeof(type)>::type* = nullptr>
        uexpr(const uexpr<T2>& exp){
            *this = exp;
        }
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) < sizeof(type)>::type* = nullptr>
        uexpr(uexpr<T2>&& exp){
            *this = move(exp);
        }
        
        uexpr(const uexpr& exp){
            *this = exp;
        }
        
        uexpr(uexpr&& exp){
            *this = move(exp);
        }
        
        uexpr(OperatorType ot, shared_ptr<constant_> son){
            _otype = ot;
            _son = son;
            this->_type = uexp_c;
            this->_dim[0] = son->_dim[0];
            this->_dim[1] = son->_dim[1];
            this->_to_str = this->to_str();
            this->_is_vector = son->_is_vector;
            this->_is_transposed = son->_is_transposed;
        }
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) < sizeof(type)>::type* = nullptr>
        uexpr& operator=(uexpr<T2>&& exp){
            this->_type = uexp_c;
            _son = move(exp._son);
            _otype = exp._otype;
            this->_all_convexity = exp._all_convexity;
            this->_all_sign = exp._all_sign;
            if(exp._range){
                this->_range = make_shared<pair<type,type>>();
                this->_range->first = exp._range->first;
                this->_range->second = exp._range->second;
            }
            this->_to_str = exp._to_str;
            this->_coef = exp._coef;
            this->_is_vector = exp._is_vector;
            this->_is_transposed = exp._is_transposed;
            this->_dim[0] = exp._dim[0]; this->_dim[1] = exp._dim[1];
            return *this;
        }
        
        
        uexpr& operator=(uexpr&& exp){
            this->_type = uexp_c;
            _son = move(exp._son);
            _otype = exp._otype;
            this->_all_convexity = exp._all_convexity;
            this->_all_sign = exp._all_sign;
            this->_range = move(exp._range);
            this->_to_str = exp._to_str;
            this->_coef = exp._coef;
            this->_is_vector = exp._is_vector;
            this->_is_transposed = exp._is_transposed;
            this->_dim[0] = exp._dim[0]; this->_dim[1] = exp._dim[1];
            return *this;
        }
        
        bool operator==(const uexpr &c)const{
            //        return (_otype == c._otype && equals(_son,c._son));
            return (this->_to_str.compare(c._to_str)==0);
        }
        
        
       
        
        /* UNARY EXPRESSIONS */
        
        uexpr(){
            this->_type = uexp_c;
        }
        
        void print() {
            cout << this->_to_str << endl;
        }
        
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) < sizeof(type)>::type* = nullptr>
        uexpr& operator=(const uexpr<T2>& exp){
            this->_type = uexp_c;
            _son = exp._son->copy();
            _otype = exp._otype;
            this->_all_convexity = exp._all_convexity;
            this->_all_sign = exp._all_sign;
            if(exp._range){
                this->_range = make_shared<pair<type,type>>();
                this->_range->first = exp._range->first;
                this->_range->second = exp._range->second;
            }
            this->_to_str = exp._to_str;
            this->_coef = exp._coef;
            this->_is_vector = exp._is_vector;
            this->_is_transposed = exp._is_transposed;
            this->_dim[0] = exp._dim[0]; this->_dim[1] = exp._dim[1];
            return *this;
        }
        
        uexpr& operator=(const uexpr& exp){
            this->_type = uexp_c;
            _son = exp._son->copy();
            _otype = exp._otype;
            this->_all_convexity = exp._all_convexity;
            this->_all_sign = exp._all_sign;
            if(exp._range){
                this->_range = make_shared<pair<type,type>>();
                this->_range->first = exp._range->first;
                this->_range->second = exp._range->second;
            }
            this->_to_str = exp._to_str;
            this->_coef = exp._coef;
            this->_is_vector = exp._is_vector;
            this->_is_transposed = exp._is_transposed;
            this->_dim[0] = exp._dim[0]; this->_dim[1] = exp._dim[1];
            return *this;
        }
        
        string to_str(){
            string str;
            if (this->_coef!=unit<type>().eval()) {
                if (this->_coef!=-1.*unit<type>().eval()) {
                    str+= to_string_with_precision(this->_coef,3);
                }
                else {
                    str+= "-";
                }
            }
            str += operator_str(_otype) +"("+_son->to_str()+")";
            return str;
        }
        
        
        string to_str(int prec){
            string str;
            if (this->_coef!=unit<type>().eval()) {
                if (this->_coef!=-1.*unit<type>().eval()) {
                    str+= to_string_with_precision(this->_coef,prec);
                }
                else {
                    str+= "-";
                }
            }
            str += operator_str(_otype) +"("+_son->to_str(prec)+")";
            return str;
        }
        
        string to_str(size_t inst, int prec) {
            string str;
            if (this->_coef!=unit<type>().eval()) {
                if (this->_coef!=-1.*unit<type>().eval()) {
                    str+= to_string_with_precision(this->_coef,prec);
                }
                else {
                    str+= "-";
                }
            }
            str += operator_str(_otype) +"("+_son->to_str(inst,prec)+")";
            return str;
        }
        
        string to_str(size_t inst1, size_t inst2, int prec) {
            string str;
            if (this->_coef!=unit<type>().eval()) {
                if (this->_coef!=-1.*unit<type>().eval()) {
                    str+= to_string_with_precision(this->_coef,prec);
                }
                else {
                    str+= "-";
                }
            }
            str += operator_str(_otype) +"("+_son->to_str(inst1,inst2,prec)+")";
            return str;
        }
   
    };

    template<typename type = double>
    class bexpr: public expr<type>{
    private:
        
    public:
        OperatorType               _otype = id_;
        shared_ptr<constant_>      _lson = nullptr;
        shared_ptr<constant_>      _rson = nullptr;
        
       
        
        bexpr(){
            this->_type = bexp_c;
        }
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) < sizeof(type)>::type* = nullptr>
        bexpr(const bexpr<T2>& exp){ /**< Copy constructor from binary expression tree */
            *this = exp;
        };
        
        bexpr(const bexpr& exp){ /**< Copy constructor from binary expression tree */
            *this = exp;
        };
        
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) < sizeof(type)>::type* = nullptr>
        bexpr(bexpr<T2>&& exp){ /**< Move constructor from binary expression tree */
            *this = move(exp);
        };
        
        bexpr(bexpr&& exp){ /**< Move constructor from binary expression tree */
            *this = move(exp);
        };
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) < sizeof(type)>::type* = nullptr>
        bexpr& operator=(bexpr<T2>&& exp){
            this->_type = bexp_c;
            _lson = move(exp._lson);
            _rson = move(exp._rson);
            _otype = exp._otype;
            this->_all_convexity = exp._all_convexity;
            this->_all_sign = exp._all_sign;
            if(exp._range){
                this->_range = make_shared<pair<type,type>>();
                this->_range->first = exp._range->first;
                this->_range->second = exp._range->second;
            }
            this->_to_str = exp._to_str;
            this->_coef = exp._coef;
            this->_is_vector = exp._is_vector;
            this->_is_transposed = exp._is_transposed;
            this->_dim[0] = exp._dim[0]; this->_dim[1] = exp._dim[1];
            return *this;
        }
        
        bexpr& operator=(bexpr&& exp){
            this->_type = bexp_c;
            _lson = move(exp._lson);
            _rson = move(exp._rson);
            _otype = exp._otype;
            this->_all_convexity = exp._all_convexity;
            this->_all_sign = exp._all_sign;
            this->_range = move(exp._range);
            this->_to_str = exp._to_str;
            this->_coef = exp._coef;
            this->_is_vector = exp._is_vector;
            this->_is_transposed = exp._is_transposed;
            this->_dim[0] = exp._dim[0]; this->_dim[1] = exp._dim[1];
            return *this;
        }
        
        bool operator==(const bexpr &c)const{
            //        return (_otype == c._otype && equals(_lson,c._lson) && equals(_rson,c._rson));
            return (this->_to_str.compare(c._to_str)==0);
        }
        
        void in(const indices& ids){
            if(_lson->is_function()){
                auto f = static_pointer_cast<func<type>>(_lson);
                f->in(ids);
            }
            if(_rson->is_function()){
                auto f = static_pointer_cast<func<type>>(_rson);
                f->in(ids);
            }
        };
        
        void propagate_dim(size_t d){
            if(this->_is_transposed){
                this->_dim[1] = d;
            }
            else {
                this->_dim[0] = d;
            }
            _lson->propagate_dim(d);
            _rson->propagate_dim(d);
        }
        
        void uneval(){
            _lson->uneval();
            _rson->uneval();
        }
        
        /** allocates memory for current and all sub-functions */
        void allocate_mem(){
            _lson->allocate_mem();
            _rson->allocate_mem();
        };
        
        void print() {
            cout << this->_to_str << endl;
        }
        
        string to_str() {
            string str;
            if (this->_coef!=unit<type>().eval()) {
                if (this->_coef!=-1.*unit<type>().eval()) {
                    str+= to_string_with_precision(this->_coef, 3);
                }
                else {
                    str+= "-";
                }
                str+="(";
            }
            if((_otype==product_ || _otype==div_) && (_lson->get_type()==uexp_c || _lson->get_type()==bexp_c)) {
                str += "(";
                str+= _lson->to_str();
                str += ")";
            }
            else
                str+= _lson->to_str();
            
            if (_otype==plus_) {
                str+= " + ";
            }
            if (_otype==minus_) {
                str+= " - ";
            }
            if (_otype==product_) {
                str+= " * ";
            }
            if (_otype==div_) {
                str+= "/";
            }
            
            if (_otype==power_) {
                str+= "^";
            }
            
            if (_otype==plus_ || (_rson->get_type()!=uexp_c && _rson->get_type()!=bexp_c)) {
                str+= _rson->to_str();
            }
            else {
                str+= "(";
                str+= _rson->to_str();
                str+= ")";
            }
            if (this->_coef!=unit<type>().eval()) {
                str += ")";
            }
            return str;
        }
        
        string to_str(int prec) {
            string str;
            if (this->_coef!=unit<type>().eval()) {
                if (this->_coef!=-1.*unit<type>().eval()) {
                    str+= to_string_with_precision(this->_coef, prec);
                }
                else {
                    str+= "-";
                }
                str+="(";
            }
            if((_otype==product_ || _otype==div_) && (_lson->get_type()==uexp_c || _lson->get_type()==bexp_c)) {
                str += "(";
                str+= _lson->to_str(prec);
                str += ")";
            }
            else
                str+= _lson->to_str(prec);
            
            if (_otype==plus_) {
                str+= " + ";
            }
            if (_otype==minus_) {
                str+= " - ";
            }
            if (_otype==product_) {
                str+= " * ";
            }
            if (_otype==div_) {
                str+= "/";
            }
            
            if (_otype==power_) {
                str+= "^";
            }
            
            if (_otype==plus_ || (_rson->get_type()!=uexp_c && _rson->get_type()!=bexp_c)) {
                str+= _rson->to_str(prec);
            }
            else {
                str+= "(";
                str+= _rson->to_str(prec);
                str+= ")";
            }
            if (this->_coef!=unit<type>().eval()) {
                str += ")";
            }
            return str;
        }
        
        string to_str(size_t inst,int prec) {
            string str;
            if (this->_coef!=unit<type>().eval()) {
                if (this->_coef!=-1.*unit<type>().eval()) {
                    str+= to_string_with_precision(this->_coef, prec);
                }
                else {
                    str+= "-";
                }
                str+="(";
            }
            if((_otype==product_ || _otype==div_) && (_lson->get_type()==uexp_c || _lson->get_type()==bexp_c)) {
                str += "(";
                str+= _lson->to_str(inst,prec);
                str += ")";
            }
            else
                str+= _lson->to_str(inst,prec);
            
            if (_otype==plus_) {
                str+= " + ";
            }
            if (_otype==minus_) {
                str+= " - ";
            }
            if (_otype==product_) {
                str+= " * ";
            }
            if (_otype==div_) {
                str+= "/";
            }
            
            if (_otype==power_) {
                str+= "^";
            }
            
            if (_otype==plus_ || (_rson->get_type()!=uexp_c && _rson->get_type()!=bexp_c)) {
                str+= _rson->to_str(inst,prec);
            }
            else {
                str+= "(";
                str+= _rson->to_str(inst,prec);
                str+= ")";
            }
            if (this->_coef!=unit<type>().eval()) {
                str += ")";
            }
            return str;
        }
        
        string to_str(size_t inst1,size_t inst2,int prec) {
            string str;
            if (this->_coef!=unit<type>().eval()) {
                if (this->_coef!=-1.*unit<type>().eval()) {
                    str+= to_string_with_precision(this->_coef, prec);
                }
                else {
                    str+= "-";
                }
                str+="(";
            }
            if((_otype==product_ || _otype==div_) && (_lson->get_type()==uexp_c || _lson->get_type()==bexp_c)) {
                str += "(";
                str+= _lson->to_str(inst1,inst2,prec);
                str += ")";
            }
            else
                str+= _lson->to_str(inst1,inst2,prec);
            
            if (_otype==plus_) {
                str+= " + ";
            }
            if (_otype==minus_) {
                str+= " - ";
            }
            if (_otype==product_) {
                str+= " * ";
            }
            if (_otype==div_) {
                str+= "/";
            }
            
            if (_otype==power_) {
                str+= "^";
            }
            
            if (_otype==plus_ || (_rson->get_type()!=uexp_c && _rson->get_type()!=bexp_c)) {
                str+= _rson->to_str(inst1,inst2,prec);
            }
            else {
                str+= "(";
                str+= _rson->to_str(inst1,inst2,prec);
                str+= ")";
            }
            if (this->_coef!=unit<type>().eval()) {
                str += ")";
            }
            return str;
        }
        
        bool is_inner_product() const{
            return _otype==product_ && (_lson->get_dim(1)==_rson->get_dim(0) || (_lson->_is_transposed && _lson->get_dim(0)==_rson->get_dim(0)));
        }
        
        
        bexpr(OperatorType otype, shared_ptr<constant_> lson, shared_ptr<constant_> rson);
        
        template<class T2, typename enable_if<is_convertible<T2, type>::value && sizeof(T2) < sizeof(type)>::type* = nullptr>
        bexpr& operator=(const bexpr<T2>& exp){
            this->_type = bexp_c;
            _lson = exp._lson->copy();
            _rson = exp._rson->copy();
            _otype = exp._otype;
            this->_all_convexity = exp._all_convexity;
            this->_all_sign = exp._all_sign;
            if(exp._range){
                this->_range = make_shared<pair<type,type>>();
                this->_range->first = exp._range->first;
                this->_range->second = exp._range->second;
            }
            this->_to_str = exp._to_str;
            this->_coef = exp._coef;
            this->_is_vector = exp._is_vector;
            this->_is_transposed = exp._is_transposed;
            this->_dim[0] = exp._dim[0]; this->_dim[1] = exp._dim[1];
            return *this;
        }
        
        bexpr& operator=(const bexpr& exp){
            this->_type = bexp_c;
            _lson = exp._lson->copy();
            _rson = exp._rson->copy();
            _otype = exp._otype;
            this->_all_convexity = exp._all_convexity;
            this->_all_sign = exp._all_sign;
            if(exp._range){
                this->_range = make_shared<pair<type,type>>();
                this->_range->first = exp._range->first;
                this->_range->second = exp._range->second;
            }
            this->_to_str = exp._to_str;
            this->_coef = exp._coef;
            this->_is_vector = exp._is_vector;
            this->_is_transposed = exp._is_transposed;
            this->_dim[0] = exp._dim[0]; this->_dim[1] = exp._dim[1];
            return *this;
        }
        
        void reset(){
            _otype = id_;
            this->_to_str = "noname";
            this->_coef = 1.;
            _lson = nullptr;
            _rson = nullptr;
        };
        
        
        shared_ptr<constant_> get_lson() const{
            return _lson;
        };
        
        shared_ptr<constant_> get_rson() const{
            return _rson;
        };
        
        void set_lson(shared_ptr<constant_> c){
            _lson = c;
        };
        
        void set_rson(shared_ptr<constant_> c){
            _rson = c;
        };
        
        OperatorType get_otype() const {
            return _otype;
        };
        
        shared_ptr<constant_> copy()const{return make_shared<bexpr>(*this);};
        
        
        bool operator!=(const bexpr& c) const{
            return !(*this==c);
        };
        
        Sign get_all_sign() const{
            return unknown_; // TO UPDATE
        }
        
                
        void print(size_t inst) {
            cout << to_str(inst) << endl;
        }
        
        
    };

}
#endif /* expr_h */

