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

namespace gravity {
    class func_;
    /** Backbone class for unary and binary expressions. */

    class expr: public constant_{
    protected:
    public:
        double                                 _coef = 1.; /**< coefficient multpying the expression */
        string                                 _to_str = "noname"; /**< A string representation of the expression */
        
//        func_ get_derivative(const param_ &v) const;
        void propagate_dim(size_t);/*<< Propagates number of indices */
        void allocate_mem();
        void reverse_sign(){ _coef *= -1; };
    };


    /** Class uexpr (unary expression), stores a unary expression tree. */
    class uexpr: public expr{
        
    public:
        OperatorType                    _otype = id_;
        shared_ptr<constant_>           _son = nullptr;
        
        uexpr();
        uexpr(const uexpr& exp);
        uexpr(uexpr&& exp);
        uexpr(OperatorType ot, shared_ptr<constant_> son);
        uexpr& operator=(const uexpr& e);
        uexpr& operator=(uexpr&& e);
        
        
        ~uexpr(){};
        
        
        void reset(){
            _son = nullptr;
            _otype = id_;
            _to_str = "noname";
            _coef = 1.;
        };
        
        
        OperatorType get_otype() const{
            return _otype;
        };
        
        /** Operators */
        
        shared_ptr<constant_> copy()const{return make_shared<uexpr>(*this);};
        
        bool operator==(const uexpr &c)const;
        
        bool operator!=(const uexpr& c) const{
            return !(*this==c);
        };
        
        Sign get_all_sign() const{
            return unknown_;// TO UPDATE
        }
        
        string to_str() ;
        string to_str(int prec) ;
        string to_str(size_t, int prec) ;
        string to_str(size_t, size_t, int prec) ;
        void print(bool endline = true) ;
//        func_ get_derivative(const param_ &v) const;
        vector<shared_ptr<param_>> get_nl_vars() const;
    };


    class bexpr: public expr{
    private:
        
    public:
        OperatorType               _otype = id_;
        shared_ptr<constant_>      _lson = nullptr;
        shared_ptr<constant_>      _rson = nullptr;
        
        bexpr();
        
        bexpr(OperatorType otype, shared_ptr<constant_> lson, shared_ptr<constant_> rson);
        
        bexpr(const bexpr& exp);
        
        bexpr(bexpr&& exp);
        
        bexpr& operator=(const bexpr& e);
        
        bexpr& operator=(bexpr&& e);
        
        ~bexpr(){}
        
        
        void reset(){
            _otype = id_;
            _to_str = "noname";
            _coef = 1.;
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
        
        bool operator==(const bexpr &c)const;
        
        bool operator!=(const bexpr& c) const{
            return !(*this==c);
        };
        
        Sign get_all_sign() const{
            return unknown_; // TO UPDATE
        }
        
        bool is_inner_product() const;
                
        string to_str() ;
        string to_str(int prec) ;
        string to_str(size_t, int prec) ;
        string to_str(size_t, size_t, int prec) ;
        void print(size_t inst) {
            cout << to_str(inst) << endl;
        }
        
        void print() ;
        
        void print_tree() const;
        
//        func_ get_derivative(const param_ &v) const;
        
        vector<shared_ptr<param_>> get_nl_vars() const;
        
    };

}
#endif /* expr_h */

