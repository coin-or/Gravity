//
//  Type.h
//  Gravity
//
//  Created by Hassan on 3 Jan 2016.
//

#ifndef Gravity___Type_h
#define Gravity___Type_h

#include <list>
#include <string>
namespace gravity{
    typedef unsigned int ind; /* Index type */
    //typedef std::set<ind> indx; /* Set of indices type */


    typedef enum { linear_, convex_, concave_, undet_} Convexity; /* Convexity Type */
    typedef enum { neg_ = -2, non_pos_ = -1, zero_ = 0, non_neg_ = 1, pos_ = 2, unknown_ = 3} Sign; /* Sign Type */
    typedef enum { binary_, short_, integer_, float_, double_, long_} NType; /* Number Type */
    typedef enum { binary_c, short_c, integer_c, float_c, double_c, long_c, par_c, uexp_c, bexp_c, var_c, func_c, sdpvar_c} CType; /* Constant type, ancestor to parameter, var and function */
    typedef enum { infeasible, optimal, suboptimal, unbounded, error} Outcome;
    typedef enum { geq, leq, eq } ConstraintType;
    typedef enum { const_, lin_, quad_, pol_, nlin_ } FType;  /* Function type in constraint: Constant, Linear, Quadratic, Polynomial or Nonlinear function */
    typedef enum { minimize, maximize } ObjectiveType;
    typedef enum { id_, number_, plus_, minus_, product_, div_, power_, cos_, sin_, sqrt_, exp_, log_} OperatorType;  /* Operation type in the expression tree */

    typedef enum { ordered_pairs_, unordered_ } SetType;
    typedef enum { vec_=0, in_ordered_pairs_=1, from_ordered_pairs_=2, to_ordered_pairs_=3, in_arcs_=4, from_arcs_=5, to_arcs_=6, in_nodes=7, in_set=8} IndexType;  /* Index type */

    using namespace std;


    /** Class for manipulating indices */
    class ordered_pairs{
        
    public:
        unsigned          _first;
        unsigned          _last;
        std::list<string> _keys;
        std::list<string> _from;
        std::list<string> _to;
        ordered_pairs(unsigned p1 ,unsigned p2){
            _first = p1;
            _last = p2;
            string key;
        for (int i = p1-1; i < p2; i++){
            for (int j = i+1; j < p2; j++){
                _keys.push_back(to_string(i) + "," + to_string(j));
                _from.push_back(to_string(i));
                _to.push_back(to_string(j));
            }
        }

        }
    };

    typedef enum { ipopt, gurobi, bonmin, cplex, sdpa_, mosek_} SolverType;  /* Solver type */

    // settings of solvers. used by solvers like sdpa.
    typedef enum {unsolved = -1, penalty=0, fast=1, medium=2, stable=3} SolverSettings;

    typedef tuple<unsigned,IndexType,unsigned,unsigned> unique_id; /* A unique identifier is defined as a tuple<variable index, index_type, first_index, last_index */
}

#endif
