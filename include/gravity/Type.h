//
//  Type.h
//  Gravity
//
//  Created by Hassan on 3 Jan 2016.
//

#ifndef Gravity___Type_h
#define Gravity___Type_h
typedef unsigned int ind; /* Index type */
//typedef std::set<ind> indx; /* Set of indices type */

typedef enum { linear_, convex_, concave_, indet_} Convexity; /* Number Type */

typedef enum { binary_, short_, integer_, float_, double_, long_} NType; /* Number Type */
typedef enum { binary_c, short_c, integer_c, float_c, double_c, long_c, par_c, uexp_c, bexp_c, var_c, func_c} CType; /* Constant type, ancestor to parameter, var and function */
typedef enum { infeasible, optimal, suboptimal, unbounded, error } Outcome;
typedef enum { geq, leq, eq } ConstraintType;
typedef enum { const_, lin_, quad_, pol_, nlin_ } FType;  /* Function type in constraint: Constant, Linear, Quadratic, Polynomial or Nonlinear function */
typedef enum { minimize, maximize } ObjectiveType;
typedef enum { id_, number_, plus_, minus_, product_, div_, power_, cos_, sin_, sqrt_, exp_, log_} OperatorType;  /* Operation type in the expression tree */
typedef enum { ipopt, gurobi, bonmin } SolverType;  /* Solver type */
#endif
