//
//  Type.h
//  Gravity
//
//  Created by Hassan on 3 Jan 2016.
//

#ifndef Gravity___Type_h
#define Gravity___Type_h

#include <list>
#include <assert.h>
#include <string>

namespace gravity{
#define EPS 0.00001
#define Real double
//#define Integer integer  //same name with a boost graph library.
#define Binary bool
#define Debug(x)
#define DebugOn(x) cout << x
#define DebugOff(x)

    typedef unsigned int ind; /* Index type */
    //typedef std::set<ind> indx; /* Set of indices type */

    typedef enum { linear_, convex_, concave_, undet_} Convexity; /* Convexity Type */
    typedef enum { neg_ = -2, non_pos_ = -1, zero_ = 0, non_neg_ = 1, pos_ = 2, unknown_ = 3} Sign; /* Sign Type */
    typedef enum { binary_, short_, integer_, float_, double_, long_} NType; /* Number Type */
    typedef enum { binary_c, short_c, integer_c, float_c, double_c, long_c, par_c, uexp_c, bexp_c, var_c, func_c, sdpvar_c} CType; /* Constant type, ancestor to parameter, var and function */
    typedef enum { infeasible, optimal, suboptimal, unbounded, error} Outcome;
    typedef enum { geq, leq, eq } ConstraintType;
    typedef enum { const_, lin_, quad_, pol_, nlin_ } FType;  /* Function type in constraint: Constant, Linear, Quadratic, Polynomial or Nonlinear function */
    typedef enum { lin_m, quad_m, pol_m, nlin_m } MType;  /* Model type: Linear, Quadratic, Polynomial or Nonlinear function */
    typedef enum { minimize, maximize } ObjectiveType;
    typedef enum { id_, number_, plus_, minus_, product_, div_, power_, cos_, sin_, sqrt_, exp_, log_} OperatorType;  /* Operation type in the expression tree */

    typedef enum { R_, R_p_, C_, C_p_} SpaceType;  /* Operation type in the expression tree */

    typedef enum { ordered_pairs_, unordered_ } SetType;
//    typedef enum { vec_=0, in_ordered_pairs_=1, from_ordered_pairs_=2, to_ordered_pairs_=3, in_arcs_=4, from_arcs_=5, to_arcs_=6, in_nodes_=7, in_set_=8, mask_=9, in_bags_=10, time_expand_ = 11, in_set_at_} IndexType;  /* Index type */

    typedef enum { unindexed_, in_, in_pairs_, out_, from_, to_, in_at_, in_time_, from_time_, to_time_, in_arcs_, out_arcs_, in_gens_, in_pairs_time_, out_arcs_time_, in_arcs_time_, in_gens_time_} IndexType;  /* Index type */

    using namespace std;


    /** Class for manipulating indices */
    class index_{
    public:
        string _name;
        string _type_name="index";
        bool   _active = true;
        index_(const string& name, bool active=true):_name(name), _active(active){};
        index_(const index_& idx):_name(idx._name), _active(idx._active){};
    };
    
    class index_pair{
    public:
        string _name;
        string _type_name="index_pair";
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
        unsigned          _first;
        unsigned          _last;
        std::vector<index_pair*> _keys;
        ordered_pairs(unsigned p1 ,unsigned p2){
            _first = p1;
            _last = p2;
            auto n = p2 - p1 + 1;
            assert(n >= 0);
            _keys.resize(n*(n-1)/2);
            string key;
            unsigned index = 0;
            for (int i = p1-1; i < p2; i++){
                for (int j = i+1; j < p2; j++){
                    _keys[index++] = new index_pair(index_(to_string(i)), index_(to_string(j)));
                }
            }
        }
        ~ordered_pairs(){
            for (auto p: _keys) { delete p;}
        }
    };

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
    
    typedef enum { ipopt, gurobi, bonmin, cplex, sdpa, Mosek} SolverType;  /* Solver type */

    // settings of solvers. used by solvers like sdpa.
    typedef enum {unsolved = -1, penalty=0, fast=1, medium=2, stable=3} SolverSettings;

    typedef tuple<unsigned,IndexType,size_t,unsigned,unsigned> unique_id; /* A unique identifier is defined as a tuple<variable index, index type, variable type, first_index, last_index */

    template <class T>
    std::string type_name(const T& t) {
        return t._type_name;
    }
}

#endif
