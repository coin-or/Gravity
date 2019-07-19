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
    
    
    class lterm{
        
    public:
        shared_ptr<constant_>               _coef = nullptr; /**< Coefficient */
        shared_ptr<param_>                  _p = nullptr; /**< Variable */
        bool                                _sign = true; /**< True if +, false if - */
        bool                                _in_S = false; /**< True if the term is in S (when creating the on/off constraints), false if not */
        
        lterm(){};
        
        lterm(lterm&& t){
            *this = move(t);
        };
        
        lterm(const lterm& t){
            *this = t;
        };
        
        
//        lterm(shared_ptr<param_> p):lterm(true,p){
//        };
//        
//        
//        lterm(bool sign, shared_ptr<param_> p){
//            _coef = make_shared<constant<double>>(1);
//            _p = p;
//            _sign = sign;
//        };
        
        lterm(shared_ptr<constant_> coef, shared_ptr<param_> p):lterm(true,coef,p){};
        
        lterm(bool sign, shared_ptr<constant_> coef, shared_ptr<param_> p);
        
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
        
        
        bool operator==(const lterm& l) const;
        bool operator!=(const lterm& l) const;
        
        lterm& operator=(const lterm& l);
        lterm& operator=(lterm&& l);
        
        string print_transposed(int prec) const;
        string print_transposed(size_t idx, int prec) const;
        string to_str() const;
        string to_str(size_t ind, int prec) const;
        string to_str(size_t ind, size_t inst, int prec) const;
        void print(size_t ind) const;
    };


    /** A class to represent a quadratic term, e.g. 2xy or 3x^2. */
    class qterm{
        
    public:
        shared_ptr<constant_>                                        _coef = nullptr;
        shared_ptr<pair<shared_ptr<param_>,shared_ptr<param_>>>      _p = nullptr;
        bool                                                         _sign = true; /**< True if +, false if - */
        bool                                                         _coef_p1_tr = false; /**< True if the qterm is (coef*p1)^T*p2 */
        
        qterm(){}
        
        qterm(const qterm& t){
            *this = t;
        };
        
        qterm(qterm&& t){
            *this = move(t);
        };
        
        
//        qterm(shared_ptr<param_> p1, shared_ptr<param_> p2):qterm(true, p1, p2){};
        
        
        qterm(shared_ptr<constant_> coef, shared_ptr<param_> p1, shared_ptr<param_> p2):qterm(true, coef, p1, p2){};
        
//        qterm(bool sign, shared_ptr<param_> p1, shared_ptr<param_> p2):qterm(true, unit<type>(), p1, p2){};
        
        qterm(bool sign, shared_ptr<constant_> coef, shared_ptr<param_> p1, shared_ptr<param_> p2){
            _coef = coef;
            _p = make_shared<pair<shared_ptr<param_>, shared_ptr<param_>>>(make_pair(p1,p2));
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
                _coef_p1_tr = true;
//                throw invalid_argument("Check the transpose operator, there seems to be a dimension issue\n");
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
        
        Convexity get_convexity() const {
            if(_p->first == _p->second){
                if (_sign && _coef->is_non_negative()) {
                    return convex_;
                }
                if (_sign && _coef->is_non_positive()) {
                    return concave_;
                }
                if (!_sign && _coef->is_non_negative()) {
                    return concave_;
                }
                if (!_sign && _coef->is_non_positive()) {
                    return convex_;
                }
            }
            return undet_;
        }

        
        bool operator==(const qterm& l) const;
        bool operator!=(const qterm& l) const;
        
        qterm& operator=(const qterm& l);
        qterm& operator=(qterm&& l);

        string to_str() const;
        string to_str(size_t ind, int prec) const;
        string to_str(size_t ind, size_t inst, int prec) const;
        void print(size_t ind) const;
        string print_transposed(int prec) const;
        string print_transposed(size_t idx, int prec) const;
    };


    /** A class to represent a polynomial term, e.g. 2xyz or 3x^2y. */
    class pterm{
        
    public:
        shared_ptr<constant_>                                 _coef = nullptr;
        shared_ptr<list<pair<shared_ptr<param_>, int>>>       _l = nullptr; /**< A polynomial term is represented as a list of pairs <param_*,int> where the first element points to the parameter and the second indicates the exponent */
        bool                                                  _sign = true; /**< True if +, false if - */
        
        pterm(){}
        
        
        pterm(bool sign, shared_ptr<constant_> coef, shared_ptr<param_> p, int exp){
            _coef = coef;
            if (coef->_is_transposed && p->_is_transposed) {
                throw invalid_argument("Check the transpose operator, there seems to be a dimension issue\n");
            }
            if (p->_is_transposed) {
                throw invalid_argument("Check the transpose operator, there seems to be a dimension issue\n");
            }
            _l = make_shared<list<pair<shared_ptr<param_>, int>>>();
            _l->push_back(make_pair<>(p, exp));
            _sign = sign;
        };
        
        pterm(const pterm& t){
            *this = t;
        };
        
        pterm(pterm&& t){
            *this = move(t);
        };
        
        pterm(bool sign, shared_ptr<constant_> coef, shared_ptr<list<pair<shared_ptr<param_>, int>>> l){
            shared_ptr<param_> p1 = nullptr;
            shared_ptr<param_> p2 = nullptr;
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
                
        
        bool operator==(const pterm& l) const;
        bool operator!=(const pterm& l) const;
        
        pterm& operator=(const pterm& l);
        pterm& operator=(pterm&& l);
        string print_poly_vars() const;
        string print_poly_vars(size_t ind) const;
        string print_poly_vars(size_t ind1, size_t ind2) const;
        string to_str() const;
        string to_str(size_t ind, int prec) const;
        string to_str(size_t ind, size_t inst, int prec) const;
        void print(size_t ind) const;
        
    };

}
#endif /* poly_h */

