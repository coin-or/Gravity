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
        string                                 _to_str; /**< A string representation of the expression */
        
        string get_str();
        string to_str(size_t inst);
        double eval(size_t i) const;
        double eval(size_t i, size_t j) const;
        func_ get_derivative(const param_ &v) const;
        void propagate_dim(size_t);/*<< Propagates number of indices */
        void allocate_mem();
        void reset_val();
        virtual ~expr(){};
    };


    /** Class uexpr (unary expression), stores a unary expression tree. */
    class uexpr: public expr{
        
    public:
        OperatorType                _otype;
        shared_ptr<func_>           _son;
        
        uexpr();
        uexpr(const uexpr& exp);
        uexpr(uexpr&& exp);
        uexpr(OperatorType ot, shared_ptr<func_> son);
        uexpr& operator=(const uexpr& e);
        uexpr& operator=(uexpr&& e);
        
        
        ~uexpr(){};
                
        void reset_val();
        
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
        
        
        bool operator==(const uexpr &c)const;
        
        bool operator!=(const uexpr& c) const{
            return !(*this==c);
        };
        
        Sign get_all_sign() const{
            return unknown_;// TO UPDATE
        }
        
        double eval(size_t i) const;
        double eval(size_t i, size_t j) const;
        
        double eval() const{
            return eval(0);
        }
        
        string to_str() const;
        string to_str(size_t) const;
        void print(bool endline = true) const;
        func_ get_derivative(const param_ &v) const;
        vector<shared_ptr<param_>> get_nl_vars() const;
    };


    class bexpr: public expr{
    private:
        
    public:
        OperatorType    _otype;
        shared_ptr<func_>      _lson;
        shared_ptr<func_>      _rson;
        
        bexpr();
        
        bexpr(OperatorType otype, shared_ptr<func_> lson, shared_ptr<func_> rson);
        
        bexpr(const bexpr& exp);
        
        bexpr(bexpr&& exp);
        
        bexpr& operator=(const bexpr& e);
        
        bexpr& operator=(bexpr&& e);
        
        ~bexpr(){}
        
        void reset_val();
        
        void reset(){
            _otype = id_;
            _to_str = "noname";
            _coef = 1.;            
            _lson = nullptr;
            _rson = nullptr;
        };
        
        
        shared_ptr<func_> get_lson() const{
            return _lson;
        };
        
        shared_ptr<func_> get_rson() const{
            return _rson;
        };
        
        void set_lson(shared_ptr<func_> c){
            _lson = c;
        };
        
        void set_rson(shared_ptr<func_> c){
            _rson = c;
        };
        
        OperatorType get_otype() const {
            return _otype;
        };
        
        
        bool operator==(const bexpr &c)const;
        
        bool operator!=(const bexpr& c) const{
            return !(*this==c);
        };
        
        Sign get_all_sign() const{
            return unknown_; // TO UPDATE
        }
        
        bool is_inner_product() const;
        
        template<typename other_type> bexpr& operator +=(const other_type& v);
        template<typename other_type> bexpr& operator -=(const other_type& v);
        template<typename other_type> bexpr& operator *=(const other_type& v);
        template<typename other_type> bexpr& operator /=(const other_type& v);
        
        
        string to_str() const;
        string to_str(size_t) const;
        void print(size_t inst) const{
            cout << to_str(inst) << endl;
        }
        
        void print() const;
        
        void print_tree() const;
        
        double eval(size_t i) const;
        double eval(size_t i, size_t j) const;
        
        func_ get_derivative(const param_ &v) const;
        
        vector<shared_ptr<param_>> get_nl_vars() const;
        
    };

}
#endif /* expr_h */

