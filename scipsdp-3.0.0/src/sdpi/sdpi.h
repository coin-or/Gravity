/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/* This file is part of SCIPSDP - a solving framework for mixed-integer      */
/* semidefinite programs based on SCIP.                                      */
/*                                                                           */
/* Copyright (C) 2011-2013 Discrete Optimization, TU Darmstadt               */
/*                         EDOM, FAU Erlangen-Nürnberg                       */
/*               2014-2017 Discrete Optimization, TU Darmstadt               */
/*                                                                           */
/*                                                                           */
/* This program is free software; you can redistribute it and/or             */
/* modify it under the terms of the GNU Lesser General Public License        */
/* as published by the Free Software Foundation; either version 3            */
/* of the License, or (at your option) any later version.                    */
/*                                                                           */
/* This program is distributed in the hope that it will be useful,           */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU Lesser General Public License  */
/* along with this program; if not, write to the Free Software               */
/* Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.*/
/*                                                                           */
/*                                                                           */
/* Based on SCIP - Solving Constraint Integer Programs                       */
/* Copyright (C) 2002-2017 Zuse Institute Berlin                             */
/* SCIP is distributed under the terms of the SCIP Academic Licence,         */
/* see file COPYING in the SCIP distribution.                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   sdpi.h
 * @brief  General interface methods for SDP-preprocessing (mainly fixing variables and removing empty rows/cols)
 * @author Tristan Gally
 *
 * This file specifies a generic SDP-solver interface used by SCIP to create, modify, and solve semidefinite programs of
 * the (dual) form
 *
 *   \f{eqnarray*}{
 *   	  \min & & b^T y \\
 *      \mbox{s.t.} & & \sum_{j=1}^n A_j^i y_j - A_0^i \succeq 0 \quad \forall i \leq m \\
 *      & & Dy \geq d \\
 *      & & \ell \leq y \leq u
 *   \f}
 * for symmetric matrices \f$ A_j^i \in S_{k_i} \f$ and a matrix \f$ D \in \mathbb{R}^{k_0 \times n} \f$ and query information about the solution.
 *
 * All indexing (rows, columns, blocks and variables) starts at 0.
 *
 * Although it includes a few SCIP header files, e.g., because it uses SCIP's return codes, it can be used independently of
 * any SCIP instance.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SDPI_H__
#define __SCIP_SDPI_H__


#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"
#include "sdpi/type_sdpi.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Miscellaneous Methods
 */

/**@name Miscellaneous Methods */
/**@{ */

/** gets name and potentially version of SDP-solver */
EXTERN
const char* SCIPsdpiGetSolverName(
   void
   );

/** gets description of SDP-solver (developer, webpage, ...) */
EXTERN
const char* SCIPsdpiGetSolverDesc(
   void
   );

/** gets pointer for SDP-solver - use only with great care
 *
 *  The behavior of this function depends on the solver and its use is
 *  therefore only recommended if you really know what you are
 *  doing. In general, it returns a pointer to the SDP-solver object.
 */
EXTERN
void* SCIPsdpiGetSolverPointer(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   );

/** gets default feasibility tolerance for SDP-solver in SCIP-SDP */
EXTERN
SCIP_Real SCIPsdpiGetDefaultSdpiSolverFeastol(
   void
   );

/** gets default number of increases of penalty parameter for SDP-solver in SCIP-SDP */
EXTERN
int SCIPsdpiGetDefaultSdpiSolverNpenaltyIncreases(
   void
   );

/**@} */




/*
 * SDPI Creation and Destruction Methods
 */

/**@name SDPI Creation and Destruction Methods */
/**@{ */

/** creates an sdpi object */
EXTERN
SCIP_RETCODE SCIPsdpiCreate(
   SCIP_SDPI**           sdpi,               /**< pointer to an SDP-interface structure */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler to use for printing messages, or NULL */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   BMS_BUFMEM*           bufmem              /**< buffer memory */
   );

/** deletes an sdpi object */
EXTERN
SCIP_RETCODE SCIPsdpiFree(
   SCIP_SDPI**           sdpi                /**< pointer to an SDP-interface structure */
   );

/** cloning method of the general SDP-Interface
 *
 *  @note The solver specific interface is created anew and not copied. */
EXTERN
SCIP_RETCODE SCIPsdpiClone(
   SCIP_SDPI*            oldsdpi,            /**< pointer to the SDP-interface structure that should be cloned */
   SCIP_SDPI*            newsdpi             /**< pointer to an SDP-interface structure to clone into */
   );

/**@} */




/*
 * Modification Methods
 */

/**@name Modification Methods */
/**@{ */

/** copies SDP data into SDP-solver
 *
 *  @note As the SDP-constraint-matrices are symmetric, only the upper triangular part of them must be specified.
 *  @note There must be at least one variable, the SDP- and/or LP-part may be empty.
 */
EXTERN
SCIP_RETCODE SCIPsdpiLoadSDP(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int                   nvars,              /**< number of variables */
   SCIP_Real*            obj,                /**< objective function values of variables */
   SCIP_Real*            lb,                 /**< lower bounds of variables */
   SCIP_Real*            ub,                 /**< upper bounds of variables */
   int                   nsdpblocks,         /**< number of SDP-blocks */
   int*                  sdpblocksizes,      /**< sizes of the SDP-blocks (may be NULL if nsdpblocks = sdpconstnnonz = sdpnnonz = 0) */
   int*                  sdpnblockvars,      /**< number of variables in each SDP-block (may be NULL if nsdpblocks = sdpconstnnonz = sdpnnonz = 0) */
   int                   sdpconstnnonz,      /**< number of nonzero elements in the constant matrices of the SDP-blocks */
   int*                  sdpconstnblocknonz, /**< number of nonzeros for each variable in the constant part, also the i-th entry gives the
                                               *  number of entries  of sdpconst row/col/val [i] */
   int**                 sdpconstrow,        /**< pointer to row-indices of constant matrix for each block (may be NULL if sdpconstnnonz = 0) */
   int**                 sdpconstcol,        /**< pointer to column-indices of constant matrix for each block (may be NULL if sdpconstnnonz = 0) */
   SCIP_Real**           sdpconstval,        /**< pointer to values of entries of constant matrix for each block (may be NULL if sdpconstnnonz = 0) */
   int                   sdpnnonz,           /**< number of nonzero elements in the SDP-constraint-matrices */
   int**                 sdpnblockvarnonz,   /**< sdpnblockvarnonz[i][j] gives the number of nonzeros for the j-th variable (not necessarly
                                               *  variable j) in the i-th block, this is also the length of row/col/val[i][j] */
   int**                 sdpvar,             /**< sdpvar[i][j] gives the global index of the j-th variable (according to the sorting for row/col/val)
                                               *  in the i-th block */
   int***                sdprow,             /**< pointer to the row-indices for each block and variable in this block, so row[i][j][k] gives
                                               *  the k-th nonzero of the j-th variable (not necessarly variable j) in the i-th block
                                               *  (may be NULL if sdpnnonz = 0) */
   int***                sdpcol,             /**< pointer to the column-indices for each block and variable in this block (may be NULL if sdptnnonz = 0) */
   SCIP_Real***          sdpval,             /**< pointer to the values of the nonzeros for each block and variable in this block (may be NULL if sdpnnonz = 0) */
   int                   nlpcons,            /**< number of LP-constraints */
   SCIP_Real*            lplhs,              /**< left-hand sides of LP rows (may be NULL if nlpcons = 0) */
   SCIP_Real*            lprhs,              /**< right-hand sides of LP rows (may be NULL if nlpcons = 0) */
   int                   lpnnonz,            /**< number of nonzero elements in the LP-constraint-matrix */
   int*                  lprow,              /**< row-index for each entry in lpval-array (may be NULL if lpnnonz = 0) */
   int*                  lpcol,              /**< column-index for each entry in lpval-array (may be NULL if lpnnonz = 0) */
   SCIP_Real*            lpval               /**< values of LP-constraint matrix entries (may be NULL if lpnnonz = 0) */
   );

/** adds rows to the LP-Block
 *
 *  @note Arrays are not checked for duplicates, problems may appear if indices are added more than once.
 */
EXTERN
SCIP_RETCODE SCIPsdpiAddLPRows(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int                   nrows,              /**< number of rows to be added */
   const SCIP_Real*      lhs,                /**< left-hand sides of new rows */
   const SCIP_Real*      rhs,                /**< right-hand sides of new rows */
   int                   nnonz,              /**< number of nonzero elements to be added to the LP constraint matrix */
   const int*            row,                /**< row-indices of constraint-matrix entries, going from 0 to nrows - 1, these will be changed
                                                *  to nlpcons + i */
   const int*            col,                /**< column-indices of constraint-matrix entries */
   const SCIP_Real*      val                 /**< values of constraint-matrix entries */
   );

/** deletes all rows in the given range from the LP-Block */
EXTERN
SCIP_RETCODE SCIPsdpiDelLPRows(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int                   firstrow,           /**< first row to be deleted */
   int                   lastrow             /**< last row to be deleted */
   );

/** deletes LP-rows from SDP-interface */
EXTERN
SCIP_RETCODE SCIPsdpiDelLPRowset(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int*                  dstat               /**< deletion status of LP rows <br>
                                              *   input:  1 if row should be deleted, 0 otherwise <br>
                                              *   output: new position of row, -1 if row was deleted */
   );

/** clears the whole SDP */
EXTERN
SCIP_RETCODE SCIPsdpiClear(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   );

/** changes objective coefficients of variables */
EXTERN
SCIP_RETCODE SCIPsdpiChgObj(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int                   nvars,              /**< number of variables to change objective coefficients for */
   const int*            ind,                /**< variables indices */
   const SCIP_Real*      obj                 /**< new objective coefficients */
   );

/** changes lower and upper bounds of variables */
EXTERN
SCIP_RETCODE SCIPsdpiChgBounds(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int                   nvars,              /**< number of variables to change bounds for */
   const int*            ind,                /**< variables indices */
   const SCIP_Real*      lb,                 /**< values for the new lower bounds */
   const SCIP_Real*      ub                  /**< values for the new upper bounds */
   );

/** changes left- and right-hand sides of LP rows */
SCIP_RETCODE SCIPsdpiChgLPLhRhSides(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int                   nrows,              /**< number of LP rows to change right hand sides for */
   const int*            ind,                /**< row indices between 1 and nlpcons */
   const SCIP_Real*      lhs,                /**< new values for left-hand sides */
   const SCIP_Real*      rhs                 /**< new values for right-hand sides */
   );


/*
 * Data Accessing Methods
 */

/**@name Data Accessing Methods */
/**@{ */

/** gets the number of LP-rows in the SDP */
EXTERN
SCIP_RETCODE SCIPsdpiGetNLPRows(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int*                  nlprows             /**< pointer to store the number of rows */
   );

/** gets the number of SDP-Blocks in the SDP */
EXTERN
SCIP_RETCODE SCIPsdpiGetNSDPBlocks(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int*                  nsdpblocks          /**< pointer to store the number of blocks */
   );

/** gets the number of variables in the SDP */
EXTERN
SCIP_RETCODE SCIPsdpiGetNVars(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int*                  nvars               /**< pointer to store the number of variables */
   );

/** gets the number of nonzero elements in the SDP-constraint-matrices */
EXTERN
SCIP_RETCODE SCIPsdpiGetSDPNNonz(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int*                  nnonz               /**< pointer to store the number of nonzeros in the SDP-constraint-matrices */
   );

/** gets the number of nonzero elements in the constant matrices of the SDP-Blocks */
EXTERN
SCIP_RETCODE SCIPsdpiGetConstNNonz(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int*                  nnonz               /**< pointer to store the number of nonzeros in the constant matrices of the SDP-Blocks */
   );

/** gets the number of nonzero elements in the LP-Matrix */
EXTERN
SCIP_RETCODE SCIPsdpiGetLPNNonz(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int*                  nnonz               /**< pointer to store the number of nonzeros in the LP Matrix */
   );

/** gets objective coefficients from SDP-interface */
EXTERN
SCIP_RETCODE SCIPsdpiGetObj(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int                   firstvar,           /**< first variable to get objective coefficient for */
   int                   lastvar,            /**< last variable to get objective coefficient for */
   SCIP_Real*            vals                /**< pointer to store objective coefficients (memory of size lastvar - firstvar + 1 needs to be allocated) */
   );

/** gets current variable lower and/or upper bounds from SDP-interface */
EXTERN
SCIP_RETCODE SCIPsdpiGetBounds(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int                   firstvar,           /**< first variable to get bounds for */
   int                   lastvar,            /**< last variable to get bounds for */
   SCIP_Real*            lbs,                /**< pointer to store lower bound values (memory of size lastvar - firstvar + 1 needs to be allocated), or NULL */
   SCIP_Real*            ubs                 /**< pointer to store upper bound values (memory of size lastvar - firstvar + 1 needs to be allocated), or NULL */
   );

/** gets current left-hand sides from SDP-interface */
EXTERN
SCIP_RETCODE SCIPsdpiGetLhSides(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int                   firstrow,           /**< first row to get sides for */
   int                   lastrow,            /**< last row to get sides for */
   SCIP_Real*            lhss                /**< pointer to store left-hand side values (memory of size lastvar - firstvar + 1 needs to be allocated) */
   );

/** gets current right-hand sides from SDP-interface */
EXTERN
SCIP_RETCODE SCIPsdpiGetRhSides(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int                   firstrow,           /**< first row to get sides for */
   int                   lastrow,            /**< last row to get sides for */
   SCIP_Real*            rhss                /**< pointer to store right-hand side values (memory of size lastvar - firstvar + 1 needs to be allocated) */
   );


/**@} */




/*
 * Solving Methods
 */

/**@name Solving Methods */
/**@{ */

/** solves the SDP, as start optionally a starting point for the solver may be given, if it is NULL, the solver will start from scratch */
EXTERN
SCIP_RETCODE SCIPsdpiSolve(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   SCIP_Real*            start,              /**< NULL or a starting point for the solver, this should have length nvars */
   SCIP_SDPSOLVERSETTING startsettings,      /**< settings used to start with in SDPA, currently not used for DSDP or MOSEK, set this to
                                               *  SCIP_SDPSOLVERSETTING_UNSOLVED to ignore it and start from scratch */
   SCIP_Bool             enforceslatercheck, /**< always check for Slater condition in case the problem could not be solved and printf the solution
                                                  of this check */
   SCIP_Real             timelimit           /**< after this many seconds solving will be aborted (currently only implemented for DSDP and MOSEK) */
   );

/**@} */




/*
 * Solution Information Methods
 */

/**@name Solution Information Methods */
/**@{ */

/** returns whether a solve method was successfully called after the last modification of the SDP */
EXTERN
SCIP_Bool SCIPsdpiWasSolved(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   );

/** returns whether the original problem was solved, if SCIPsdpiWasSolved = true and SCIPsdpiSolvedOrig = false, then a penalty formulation was solved */
EXTERN
SCIP_Bool SCIPsdpiSolvedOrig(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   );

/** returns true if the solver could determine whether the problem is feasible, so it returns true if the
 *  solver knows that the problem is feasible/infeasible/unbounded, it returns false if the solver does not know
 *  anything about the feasibility status and thus the functions IsPrimalFeasible etc. should not be used */
EXTERN
SCIP_Bool SCIPsdpiFeasibilityKnown(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   );

/** gets information about primal and dual feasibility of the current SDP-solution */
EXTERN
SCIP_RETCODE SCIPsdpiGetSolFeasibility(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   SCIP_Bool*            primalfeasible,     /**< pointer to store the primal feasibility status */
   SCIP_Bool*            dualfeasible        /**< pointer to store the dual feasibility status */
   );

/** returns TRUE iff SDP is proven to be primal unbounded;
 *  returns FALSE with a debug-message if the solver could not determine feasibility */
EXTERN
SCIP_Bool SCIPsdpiIsPrimalUnbounded(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   );

/** returns TRUE iff SDP is proven to be primal infeasible;
 *  returns FALSE with a debug-message if the solver could not determine feasibility */
EXTERN
SCIP_Bool SCIPsdpiIsPrimalInfeasible(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   );

/** returns TRUE iff SDP is proven to be primal feasible;
 *  returns FALSE with a debug-message if the solver could not determine feasibility */
EXTERN
SCIP_Bool SCIPsdpiIsPrimalFeasible(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   );

/** returns TRUE iff SDP is proven to be dual unbounded;
 *  returns FALSE with a debug-message if the solver could not determine feasibility */
EXTERN
SCIP_Bool SCIPsdpiIsDualUnbounded(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   );

/** returns TRUE iff SDP is proven to be dual infeasible;
 *  returns FALSE with a debug-message if the solver could not determine feasibility */
EXTERN
SCIP_Bool SCIPsdpiIsDualInfeasible(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   );

/** returns TRUE iff SDP is proven to be dual feasible;
 *  returns FALSE with a debug-message if the solver could not determine feasibility */
EXTERN
SCIP_Bool SCIPsdpiIsDualFeasible(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   );

/** returns TRUE iff the solver converged */
EXTERN
SCIP_Bool SCIPsdpiIsConverged(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   );

/** returns TRUE iff the objective limit was reached */
EXTERN
SCIP_Bool SCIPsdpiIsObjlimExc(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   );

/** returns TRUE iff the iteration limit was reached */
EXTERN
SCIP_Bool SCIPsdpiIsIterlimExc(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   );

/** returns TRUE iff the time limit was reached */
EXTERN
SCIP_Bool SCIPsdpiIsTimelimExc(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   );

/** returns the internal solution status of the solver, which has the following meaning:<br>
 *  -1: solver was not started<br>
 *  0: converged<br>
 *  1: infeasible start<br>
 *  2: numerical problems<br>
 *  3: objective limit reached<br>
 *  4: iteration limit reached<br>
 *  5: time limit reached<br>
 *  6: user termination<br>
 *  7: other */
EXTERN
int SCIPsdpiGetInternalStatus(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   );

/** returns TRUE iff SDP was solved to optimality, meaning the solver converged and returned primal and dual feasible solutions */
EXTERN
SCIP_Bool SCIPsdpiIsOptimal(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   );

/** returns TRUE iff SDP was solved to optimality or some other status was reached
 * that is still acceptable inside a Branch & Bound framework */
EXTERN
SCIP_Bool SCIPsdpiIsAcceptable(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   );

/** gets objective value of solution */
EXTERN
SCIP_RETCODE SCIPsdpiGetObjval(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   SCIP_Real*            objval              /**< pointer to store the objective value */
   );

/** gets the best lower bound on the objective (this is equal to objval, if the problem was solved successfully, but can also give a bound
 *  if we did not get a feasible solution using the penalty approach) */
EXTERN
SCIP_RETCODE SCIPsdpiGetLowerObjbound(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   SCIP_Real*            objlb               /**< pointer to store the lower bound on the objective value */
   );

/** gets dual solution vector for feasible SDPs, if dualsollength isn't equal to the number of variables this will return the needed length and
 *  a debug message */
EXTERN
SCIP_RETCODE SCIPsdpiGetSol(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   SCIP_Real*            objval,             /**< pointer to store the objective value, may be NULL if not needed */
   SCIP_Real*            dualsol,            /**< pointer to store the dual solution vector, may be NULL if not needed */
   int*                  dualsollength       /**< length of the dualsol vector, must be 0 if dualsol is NULL, if this is less than the number
                                               *  of variables in the SDP, a debug-message will be thrown and this is set to the needed value */
   );

/** gets the primal variables corresponding to the lower and upper variable-bounds in the dual problem, the last input should specify the length
 *  of the arrays, if this is less than the number of variables, the needed length will be returned and a debug-message thrown
 *
 *  @note If a variable is either fixed or unbounded in the dual problem, a zero will be returned for the non-existent primal variable. */
EXTERN
SCIP_RETCODE SCIPsdpiGetPrimalBoundVars(
   SCIP_SDPI*            sdpi,               /**< pointer to an SDP-interface structure */
   SCIP_Real*            lbvars,             /**< pointer to store the optimal values of the variables corresponding to lower bounds in the dual problems */
   SCIP_Real*            ubvars,             /**< pointer to store the optimal values of the variables corresponding to upper bounds in the dual problems */
   int*                  arraylength         /**< input: length of lbvars and ubvars<br>
                                               *  output: number of elements inserted into lbvars/ubvars (or needed length if it was not sufficient) */
   );

/** gets the number of SDP-iterations of the last solve call */
EXTERN
SCIP_RETCODE SCIPsdpiGetIterations(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int*                  iterations          /**< pointer to store the number of iterations of the last solve call */
   );

/** gets the number of calls to the SDP-solver for the last solve call */
EXTERN
SCIP_RETCODE SCIPsdpiGetSdpCalls(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int*                  calls               /**< pointer to store the number of calls to the SDP-solver for the last solve call */
   );

/** returns which settings the SDP-solver used in the last solve call */
EXTERN
SCIP_RETCODE SCIPsdpiSettingsUsed(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   SCIP_SDPSOLVERSETTING* usedsetting        /**< the setting used by the SDP-solver */
   );

/** returns which settings the SDP-solver used in the last solve call and whether primal and dual Slater condition were fullfilled */
EXTERN
SCIP_RETCODE SCIPsdpiSlaterSettings(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   SCIP_SDPSLATERSETTING* slatersetting      /**< the combination of Slater conditions and successfull settings */
   );

/** returns whether primal and dual Slater condition held for last solved SDP */
EXTERN
SCIP_RETCODE SCIPsdpiSlater(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   SCIP_SDPSLATER*       primalslater,       /**< pointer to save whether primal Slater condition held */
   SCIP_SDPSLATER*       dualslater          /**< pointer to save whether dual Slater condition held */
   );

/**@} */




/*
 * SDPi State Methods
 */



/*
 * Numerical Methods
 */

/**@name Numerical Methods */
/**@{ */

/** returns value treated as infinity in the SDP-solver */
EXTERN
SCIP_Real SCIPsdpiInfinity(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   );

/** checks if given value is treated as (plus or minus) infinity in the SDP-solver */
EXTERN
SCIP_Bool SCIPsdpiIsInfinity(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   SCIP_Real             val                 /**< value to be checked for infinity */
   );

/** gets floating point parameter of SDP-interface */
EXTERN
SCIP_RETCODE SCIPsdpiGetRealpar(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   SCIP_SDPPARAM         type,               /**< parameter number */
   SCIP_Real*            dval                /**< pointer to store the parameter value */
   );

/** sets floating point parameter of SDP-interface */
EXTERN
SCIP_RETCODE SCIPsdpiSetRealpar(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   SCIP_SDPPARAM         type,               /**< parameter number */
   SCIP_Real             dval                /**< parameter value */
   );

/** gets integer parameter of SDP-interface */
EXTERN
SCIP_RETCODE SCIPsdpiGetIntpar(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   SCIP_SDPPARAM         type,               /**< parameter number */
   int*                  ival                /**< pointer to store the parameter value */
   );

/** sets integer parameter of SDP-interface */
EXTERN
SCIP_RETCODE SCIPsdpiSetIntpar(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   SCIP_SDPPARAM         type,               /**< parameter number */
   int                   ival                /**< parameter value */
   );

/** compute and set lambdastar (only used for SDPA) */
EXTERN
SCIP_RETCODE SCIPsdpiComputeLambdastar(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   SCIP_Real             maxguess            /**< maximum guess for lambda star of all SDP-constraints */
   );

/** compute and set the penalty parameter */
EXTERN
SCIP_RETCODE SCIPsdpiComputePenaltyparam(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   SCIP_Real             maxcoeff,           /**< maximum objective coefficient */
   SCIP_Real*            penaltyparam        /**< the computed penalty parameter */
   );

/** compute and set the maximum penalty parameter */
EXTERN
SCIP_RETCODE SCIPsdpiComputeMaxPenaltyparam(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   SCIP_Real             penaltyparam,       /**< the initial penalty parameter */
   SCIP_Real*            maxpenaltyparam     /**< the computed maximum penalty parameter */
   );

/**@} */




/*
 * File Interface Methods
 */

/**@name File Interface Methods */
/**@{ */

/** reads SDP from a file */
EXTERN
SCIP_RETCODE SCIPsdpiReadSDP(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   const char*           fname               /**< file name */
   );

/** writes SDP to a file */
EXTERN
SCIP_RETCODE SCIPsdpiWriteSDP(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   const char*           fname               /**< file name */
   );

/**@} */

#ifdef __cplusplus
}
#endif

#endif
