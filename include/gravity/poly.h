//
//  func.h
//  Gravity
//
//  Created by Hijazi, Hassan on 25 Oct 18.
//
//

#ifndef poly_h
#define poly_h

#include <gravity/var.h>
#include <gravity/expr.h>
#include <stdio.h>
#include <map>
#include <iterator>
#include <queue>
#include <list>
#include <limits>
#include <set>

using namespace std;

namespace gravity {
    
    class func_;
    double t_eval(const constant_* c, size_t i=0);
    double t_eval(const constant_* c, size_t i, size_t j);
    
    string poly_to_str(const constant_* c);
    string poly_to_str(const constant_* c, size_t inst);
    string poly_to_str(const constant_* c, size_t inst1, size_t inst2);
    
    
    constant_* copy(const constant_& c2); /**< Copy c2 into a new constant_* detecting the right class, i.e., constant<>, param<>, var<> or function */
    bool equals(const constant_* c1, const constant_* c2);/**< Checks if c2 equals c1 detecting the right class, i.e., constant<>, param<>, uexpr or bexpr. */
    void set_val(param_* p, double val, size_t i = 0); /**< Polymorphic set_val */
    /** A class to represent a linear term, e.g. 2x. */
    class lterm{
        
    public:
        constant_*              _coef = nullptr; // coefficent
        param_*                 _p = nullptr; // terms.
        bool                    _sign = true; /**< True if +, flase if - */    
        
        lterm(){};
        
//        lterm(lterm&& t){
//            *this = move(t);
//        };
//
//        lterm(const lterm& t){
//            *this = t;
//        };
        
        
        lterm(param_* p):lterm(true,p){
        };
        
        
        lterm(bool sign, param_* p){
            _coef = new constant<double>(1);
            _p = p;
            _sign = sign;
        };
        
        lterm(constant_* coef, param_* p):lterm(true,coef,p){};
        
        lterm(bool sign, constant_* coef, param_* p);
        
        void reverse_sign() {
            _sign = ! _sign;
        }
        
        Sign get_all_sign() const{
            auto sign = _coef->get_all_sign();
            if (sign==unknown_) {
                return unknown_;
            }
            if (_sign) {
                return sign;
            }
            if (sign==pos_) {
                return neg_;
            }
            if (sign==neg_) {
                return pos_;
            }
            return sign;
        }
        
        double eval(size_t i) const;
        double eval(size_t i, size_t j) const;
        
        ~lterm(){
            delete _coef;
        };
        
        bool operator==(const lterm& l) const;
        bool operator!=(const lterm& l) const;
        
//        lterm& operator=(const lterm& l);
//        lterm& operator=(lterm&& l);
        
        string to_str(size_t ind) const;
        string to_str(size_t ind, size_t inst) const;
        void print(size_t ind) const;
    };


    /** A class to represent a quadratic term, e.g. 2xy or 3x^2. */
    class qterm{
        
    public:
        constant_*                  _coef;
        pair<param_*,param_*>*      _p;
        bool                        _sign = true; /**< True if +, flase if - */
        bool                        _c_p1_transposed = false; /**< True if the qterm is (coef*p1)^T*p2 */
        
        qterm(){
            _coef = nullptr;
            _p = nullptr;
        }
        
//        qterm(qterm&& t){
//            _coef = t._coef;
//            t._coef = nullptr;
//            _p = t._p;
//            t._p = nullptr;
//            _sign = t._sign;
//            _c_p1_transposed = t._c_p1_transposed;
//        };
        
        
        qterm(param_* p1, param_* p2):qterm(true, p1, p2){};
        
        
        qterm(constant_* coef, param_* p1, param_* p2):qterm(true, coef, p1, p2){};
        
        qterm(bool sign, param_* p1, param_* p2):qterm(true, new constant<double>(1), p1, p2){};
        
        qterm(bool sign, constant_* coef, param_* p1, param_* p2){
            _coef = coef;
            _p = new pair<param_*, param_*>(make_pair(p1,p2));
            _sign = sign;
//            if (coef->_is_transposed){
//                p1->_is_vector=true;
//                p2->_is_vector=true;
//            }
//            if(p1->_is_vector){
//                coef->_is_transposed=true;
//                p2->_is_vector=true;
//            }
//            if(p2->_is_vector){
//                coef->_is_transposed=true;
//                p1->_is_vector=true;
//            }
            if (coef->_is_transposed && p1->_is_transposed) {
                throw invalid_argument("Check the transpose operator, there seems to be a dimension issue\n");
            }
//            if (p1->_is_transposed) {
//                if (p1->get_dim() != p2->get_dim()) {
//                    throw invalid_argument("Check the transpose operator, there seems to be a dimension issue\n");
//                }
//            }
            if (p2->_is_transposed) {
                throw invalid_argument("Check the transpose operator, there seems to be a dimension issue\n");
            }
        };
        
        void reverse_sign() {
            _sign = ! _sign;
        }
        
        Sign get_all_sign() const{
            auto sign = _coef->get_all_sign();
            if (sign==unknown_) {
                return unknown_;
            }
            if (_sign) {
                return sign;
            }
            if (sign==pos_) {
                return neg_;
            }
            if (sign==neg_) {
                return pos_;
            }
            return sign;
        }
        
        double eval(size_t i) const;
        double eval(size_t i, size_t j) const;
        
        ~qterm(){
            delete _coef;
            if (_p) {
                delete _p;
            }
        };
        
        func_ get_lb();
        func_ get_ub();
        
        bool operator==(const qterm& l) const;
        bool operator!=(const qterm& l) const;
        
//        qterm& operator=(const qterm& l);
//        qterm& operator=(qterm&& l);
        
        string to_str(size_t ind) const;
        string to_str(size_t ind, size_t inst) const;
        void print(size_t ind) const;
    };


    /** A class to represent a polynomial term, e.g. 2xyz or 3x^2y. */
    class pterm{
        
    public:
        constant_*                      _coef;
        list<pair<param_*, int>>*       _l; /**< A polynomial term is represented as a list of pairs <param_*,int> where the first element points to the parameter and the second indicates the exponent */
        bool                            _sign = true; /**< True if +, flase if - */
        
        pterm(){
            _coef = nullptr;
            _l = nullptr;
        }
        
        
        pterm(bool sign, constant_* coef, param_* p, int exp){
            _coef = coef;
            if (coef->_is_transposed && p->_is_transposed) {
                throw invalid_argument("Check the transpose operator, there seems to be a dimension issue\n");
            }
            if (p->_is_transposed) {
                throw invalid_argument("Check the transpose operator, there seems to be a dimension issue\n");
            }
            _l = new list<pair<param_*, int>>();
            _l->push_back(make_pair<>(p, exp));
            _sign = sign;
        };
        
        
//        pterm(pterm&& t){
//            _coef = t._coef;
//            t._coef = nullptr;
//            _l = t._l;
//            t._l = nullptr;
//            _sign = t._sign;
//        };
        
        pterm(bool sign, constant_* coef, list<pair<param_*, int>>* l){
            param_* p1 = nullptr;
            param_* p2 = nullptr;
            for (auto it = l->begin(); it != l->end(); it++){
                p1 = it->first;
                if (p1->_is_transposed && next(it)!=l->end()) {
                    p2 = next(it)->first;
                    if (p1->get_dim() != p2->get_dim()) {
                        throw invalid_argument("Check the transpose operator, there seems to be a dimension issue\n");
                    }
                }
                if (p1->_is_transposed && next(it)==l->end()) {
                    throw invalid_argument("Check the transpose operator, there seems to be a dimension issue\n");
                }
            }
            _coef = coef;
            _l = l;
            _sign = sign;
        };

        
        void reverse_sign() {
            _sign = ! _sign;
        }
        
        Sign get_all_sign() const{
            auto sign = _coef->get_all_sign();
            if (sign==unknown_) {
                return unknown_;
            }
            if (_sign) {
                return sign;
            }
            if (sign==pos_) {
                return neg_;
            }
            if (sign==neg_) {
                return pos_;
            }
            return sign;
        }
        
        double eval(size_t i) const;
        double eval(size_t i, size_t j) const;
        
        ~pterm(){
            delete _coef;
            if (_l) {
                delete _l;
            }
        };
        
        bool operator==(const pterm& l) const;
        bool operator!=(const pterm& l) const;
        
//        pterm& operator=(const pterm& l);
//        pterm& operator=(pterm&& l);
        
        string to_str(size_t ind) const;
        string to_str(size_t ind, size_t inst) const;
        void print(size_t ind) const;
        
    };

}
#endif /* poly_h */

