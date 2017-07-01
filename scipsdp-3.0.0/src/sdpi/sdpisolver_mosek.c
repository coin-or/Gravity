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

/*#define SCIP_DEBUG*/
/*#define SCIP_MORE_DEBUG*/
/*#define SCIP_DEBUG_PRINTTOFILE  *//* prints each problem inserted into MOSEK to the file mosek.task */

/**@file   sdpisolver_mosek.c
 * @brief  interface for MOSEK
 * @author Tristan Gally
 *
 * As MOSEK solves the primal instead of the dual problem, for solving the problem
 *
 *   \f{eqnarray*}{
 *      \min & & b^T y \\
 *      \mbox{s.t.} & & \sum_{i \in I} A_i^{(k)} y_i - A_0^{(k)} \succeq 0 \quad \forall \ k \in K, \\
 *      & & \sum_{i \in I} d_{ij} y_i \geq c_j \quad \forall \ j \in J, \\
 *      & & \ell_i \leq y_i \leq u_i \quad \forall \ i \in I
 *   \f}
 *
 * we insert the problem
 *
 *   \f{eqnarray*}{
 *      \max & & \sum_{k \in K} A_0^{(k)} \bullet X^{(k)} + \sum_{j \in J} c_j x_j - \sum_{i \in I_u} u_i v_i + \sum_{i \in I_\ell} \ell_i w_i \\
 *      \mbox{s.t.} & & \sum_{k \in K} A_i^{(k)} \bullet X^{(k)} + \sum_{j \in J} d_{ij} x_j - 1_{\{u_i < \infty\}} v_i + 1_{\{\ell_i > -\infty\}} w_i = b_i \quad \forall \ i \in I,\\
 *      & & X^{(k)} \succeq 0 \quad \forall \ k \in K, \\
 *      & & x_j \geq 0 \quad \forall \ j \in J,\\
 *      & & v_i \geq 0 \quad \forall \ i \in I_u,\\
 *      & & w_i \geq 0 \quad \forall \ i \in I_\ell,
 *   \f}
 *
 * where \f$I_\ell := \{i \in I: \ell_i > -\infty\}\f$ and \f$I_u := \{i \in I: u < \infty\}\f$.
 */

#include <assert.h>
#include <sys/time.h>

#include "sdpi/sdpisolver.h"

#include "blockmemshell/memory.h"            /* for memory allocation */
#include "scip/def.h"                        /* for SCIP_Real, _Bool, ... */
#include "scip/pub_misc.h"                   /* for sorting */
#include "mosek.h"                           /* for MOSEK routines */
#include "sdpi/sdpsolchecker.h"              /* to check solution with regards to feasibility tolerance */

/* TODO: use  MSK_putexitfunc to catch errors
 * TODO: Think about what to do with near optimality etc. (If MOSEK cannot compute a solution that has the prescribed accuracy, then it will
 * multiply the termination tolerances with MSK_DPAR_INTPNT_CO_TOL_NEAR_REL. If the solution then satisfies the termination criteria, then
 * the solution is denoted near optimal, near feasible and so forth.) */

#define MIN_PENALTYPARAM            1e5      /**< if the penalty parameter is to be computed, this is the minimum value it will take */
#define MAX_PENALTYPARAM            1e10     /**< if the penalty parameter is to be computed, this is the maximum value it will take */
#define PENALTYPARAM_FACTOR         1e6      /**< if the penalty parameter is to be computed, the maximal objective coefficient will be multiplied by this */
#define MAX_MAXPENALTYPARAM         1e15     /**< if the maximum penaltyparameter is to be computed, this is the maximum value it will take */
#define MAXPENALTYPARAM_FACTOR      1e6      /**< if the maximum penaltyparameter is to be computed, it will be set to penaltyparam * this */
#define TOLERANCE_FACTOR            0.1      /**< all tolerances will be multiplied by this factor since MOSEK does not adhere to its own tolerances */
#define PENALTYBOUNDTOL             1E-3     /**< if the relative gap between Tr(X) and penaltyparam for a primal solution of the penaltyformulation
                                              *   is bigger than this value, it will be reported to the sdpi */
#define INFEASFEASTOLCHANGE         0.1      /**< change feastol by this factor if the solution was found to be infeasible with regards to feastol */
#define INFEASMINFEASTOL            1E-9     /**< minimum value for feasibility tolerance when encountering problems with regards to tolerance */
#define CONVERT_ABSOLUTE_TOLERANCES TRUE     /**< should absolute tolerances be converted to relative tolerances for MOSEK */

/** data used for SDP interface */
struct SCIP_SDPiSolver
{
   SCIP_MESSAGEHDLR*     messagehdlr;        /**< messagehandler for printing messages, or NULL */
   BMS_BLKMEM*           blkmem;             /**< block memory */
   BMS_BUFMEM*           bufmem;             /**< buffer memory */
   MSKenv_t              mskenv;             /**< MOSEK environement */
   MSKtask_t             msktask;            /**< MOSEK task */
   int                   nvars;              /**< number of input variables */
   int                   nactivevars;        /**< number of variables present in the dual problem in MOSEK (nvars minus the number of variables with lb = ub) */
   int*                  inputtomosekmapper; /**< entry i gives the index of input variable i in MOSEK (starting from 0) or
                                               *  -j (j=1, 2, ..., nvars - nactivevars) if the variable is fixed, the value and objective value of
                                               *  this fixed variable can be found in entry j-1 of fixedval/obj */
   int*                  mosektoinputmapper; /**< entry i gives the original index of the i-th variable in MOSEK (indices go from 0 to nactivevars-1) */
   SCIP_Real*            fixedvarsval;       /**< entry i gives the lower and upper bound of the i-th fixed variable */
   SCIP_Real             fixedvarsobjcontr;  /**< total contribution to the objective of all fixed variables, computed as sum obj * val */
   SCIP_Real*            objcoefs;           /**< objective coefficients of all active variables */
   int                   nvarbounds;         /**< number of variable bounds given to MOSEK, length of varboundpos */
   int*                  varboundpos;        /**< maps position of primal variable corresponding to variable bound to the positions
                                               *  of the corresponding variables, -n means lower bound of variable n, +n means upper bound;
                                               *  entry i gives variable bound corresponding to the primal variable in the i-th position
                                               *  of the boundblock */
   SCIP_Bool             solved;             /**< Was the SDP solved since the problem was last changed? */
   int                   sdpcounter;         /**< used for debug messages */
   SCIP_Real             epsilon;            /**< tolerance used for absolute checks */
   SCIP_Real             gaptol;             /**< this is used for checking if primal and dual objective are equal */
   SCIP_Real             feastol;            /**< feasibility tolerance that should be achieved */
   SCIP_Real             sdpsolverfeastol;   /**< feasibility tolerance for the SDP-solver */
   SCIP_Real             objlimit;           /**< objective limit for SDP solver */
   SCIP_Bool             sdpinfo;            /**< Should the SDP solver output information to the screen? */
   SCIP_Bool             penalty;            /**< was the problem last solved using a penalty formulation */
   SCIP_Bool             feasorig;           /**< was the last problem solved with a penalty formulation and with original objective coefficents
                                               *  and the solution was feasible for the original problem? */
   SCIP_Bool             rbound;             /**< was the penalty parameter bounded during the last solve call */
   MSKrescodee           terminationcode;    /**< reason for termination of the last call to the MOSEK-optimizer */
   SCIP_Bool             timelimit;          /**< was the solver stopped because of the time limit? */
   SCIP_Bool             timelimitinitial;   /**< was the problem not even given to the solver because of the time limit? */
   int                   nthreads;           /**< number of threads the SDP solver should use, currently only supported for MOSEK (-1 = number of cores) */
   int                   niterations;        /**< number of SDP-iterations since the last solve call */
   int                   nsdpcalls;          /**< number of SDP-calls since the last solve call */
};

/*
 * Local Methods
 */

/** Calls a MOSEK function and transforms the return-code to a SCIP_LPERROR if needed. */
#define MOSEK_CALL(x)  do                                                                                    \
                      {                                                                                      \
                         int _mosekerrorcode_;                                                               \
                         if ( (_mosekerrorcode_ = (x)) != MSK_RES_OK )                                       \
                         {                                                                                   \
                            SCIPerrorMessage("MOSEK-Error <%d> in function call.\n", _mosekerrorcode_);      \
                            return SCIP_LPERROR;                                                             \
                         }                                                                                   \
                      }                                                                                      \
                      while( FALSE )

/** Same as MOSEK_CALL, but used for functions returning a boolean. */
#define MOSEK_CALL_BOOL(x)  do                                                                               \
                      {                                                                                      \
                         int _mosekerrorcode_;                                                               \
                         if ( (_mosekerrorcode_ = (x)) != MSK_RES_OK )                                       \
                         {                                                                                   \
                            SCIPerrorMessage("MOSEK-Error <%d> in function call.\n", _mosekerrorcode_);      \
                            return FALSE;                                                                    \
                         }                                                                                   \
                      }                                                                                      \
                      while( FALSE )

/** Same as MOSEK_CALL, but this will be used for initialization methods with memory allocation and return a SCIP_NOMEMORY if an error is produced. */
#define MOSEK_CALLM(x) do                                                                                    \
                      {                                                                                      \
                         int _mosekerrorcode_;                                                               \
                         if ( (_mosekerrorcode_ = (x)) != MSK_RES_OK )                                       \
                         {                                                                                   \
                            SCIPerrorMessage("MOSEK-Error <%d> in function call.\n", _mosekerrorcode_);      \
                            return SCIP_NOMEMORY;                                                            \
                         }                                                                                   \
                      }                                                                                      \
                      while( FALSE )

/** Checks if a BMSallocMemory-call was successfull, otherwise returns SCIP_NOMEMORY. */
#define BMS_CALL(x)   do                                                                                     \
                      {                                                                                      \
                         if( NULL == (x) )                                                                   \
                         {                                                                                   \
                            SCIPerrorMessage("No memory in function call.\n");                               \
                            return SCIP_NOMEMORY;                                                            \
                         }                                                                                   \
                      }                                                                                      \
                      while( FALSE )

/** Calls a DSDP-Function and transforms the return-code to a SCIP_LPERROR if needed. */
#define TIMEOFDAY_CALL(x)  do                                                                                \
                      {                                                                                      \
                         int _errorcode_;                                                                    \
                         if ( (_errorcode_ = (x)) != 0 )                                                     \
                         {                                                                                   \
                            SCIPerrorMessage("Error in gettimeofday! \n");                                   \
                            return SCIP_ERROR;                                                               \
                         }                                                                                   \
                      }                                                                                      \
                      while( FALSE )

/** This will be called in all functions that want to access solution information to check if the problem was solved since the last change of the problem. */
#define CHECK_IF_SOLVED(sdpisolver)  do                                                                      \
                      {                                                                                      \
                         if (!(sdpisolver->solved))                                                          \
                         {                                                                                   \
                            SCIPerrorMessage("Tried to access solution information for SDP %d ahead of solving!\n", sdpisolver->sdpcounter);  \
                            return SCIP_LPERROR;                                                             \
                         }                                                                                   \
                      }                                                                                      \
                      while( FALSE )

/** This is the same as CHECK_IF_SOLVED, but will be called for methods returning a bool instead of a SCIP_RETURNCODE */
#define CHECK_IF_SOLVED_BOOL(sdpisolver)  do                                                                      \
                      {                                                                                      \
                         if (!(sdpisolver->solved))                                                          \
                         {                                                                                   \
                            SCIPerrorMessage("Tried to access solution information for SDP %d ahead of solving!\n", sdpisolver->sdpcounter);  \
                            assert( 0 );                                                                     \
                            return FALSE;                                                                    \
                         }                                                                                   \
                      }                                                                                      \
                      while( FALSE )

/** prints MOSEK output to the console */
static
void MSKAPI printstr(
   void*                 handle,             /**< A user-defined handle which is passed to the user-defined function */
   MSKCONST              char str[]          /**< String to print */
   )
{/*lint --e{715,818}*/
  printf("%s",str);
}

#ifndef NDEBUG
/** Test if a lower bound lb is not smaller than an upper bound ub, meaning that lb > ub - epsilon */
static
SCIP_Bool isFixed(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   SCIP_Real             lb,                 /**< lower bound */
   SCIP_Real             ub                  /**< upper bound */
   )
{/*lint --e{818}*/
   assert( sdpisolver != NULL );
   assert( lb < ub + sdpisolver->feastol );

   return (ub-lb <= sdpisolver->epsilon);
}
#else
#define isFixed(sdpisolver,lb,ub) (ub-lb <= sdpisolver->epsilon)
#endif

/*
 * Miscellaneous Methods
 */

/**@name Miscellaneous Methods */
/**@{ */

char name[SCIP_MAXSTRLEN];

/** gets name and version (if available) of SDP-solver */
const char* SCIPsdpiSolverGetSolverName(
   void
   )
{
   int major = 0;/*lint !e123*/
   int minor = 0;/*lint !e123*/
   int build = 0;
   int revision = 0;
   MSKrescodee rescodee;
#ifndef NDEBUG
   int snprintfreturn; /* used to check the return code of snprintf */
#endif

   rescodee = MSK_getversion(&major, &minor, &build, &revision);/*lint !e123*/

   if ( rescodee != MSK_RES_OK )
      return "MOSEK";

#ifndef NDEBUG
   snprintfreturn = SCIPsnprintf(name, SCIP_MAXSTRLEN, "MOSEK %d.%d.%d.%d", major, minor, build, revision);/*lint !e123*/
   assert( snprintfreturn < SCIP_MAXSTRLEN ); /* check whether the name fits into the string */
#else
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "MOSEK %d.%d.%d.%d", major, minor, build, revision);
#endif

   return name;
}

/** gets description of SDP-solver (developer, webpage, ...) */
const char* SCIPsdpiSolverGetSolverDesc(
   void
   )
{
   return "Homogeneous and self-dual interior-point solver for semidefinite programming developed by MOSEK ApS"
         "(http://www.mosek.com)";
}

/** gets pointer to SDP-solver - use only with great care
 *
 *  The behavior of this function depends on the solver and its use is
 *  therefore only recommended if you really know what you are
 *  doing. In general, it returns a pointer to the SDP-solver object.
 */
void* SCIPsdpiSolverGetSolverPointer(
   SCIP_SDPISOLVER*      sdpisolver           /**< pointer to an SDP interface solver structure */
   )
{/*lint --e{818}*/
   assert( sdpisolver != NULL );
   return (void*) NULL;
}

/** gets default feasibility tolerance for SDP-solver in SCIP-SDP */
SCIP_Real SCIPsdpiSolverGetDefaultSdpiSolverFeastol(
   void
   )
{
   return 1E-7;
}

/** gets default number of increases of penalty parameter for SDP-solver in SCIP-SDP */
int SCIPsdpiSolverGetDefaultSdpiSolverNpenaltyIncreases(
   void
   )
{
   return 8;
}

/**@} */


/*
 * SDPI Creation and Destruction Methods
 */

/**@name SDPI Creation and Destruction Methods */
/**@{ */

/** creates an SDP solver interface */
SCIP_RETCODE SCIPsdpiSolverCreate(
   SCIP_SDPISOLVER**     sdpisolver,         /**< pointer to an SDP interface solver structure */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler to use for printing messages, or NULL */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   BMS_BUFMEM*           bufmem              /**< buffer memory */
   )
{
   assert( sdpisolver != NULL );
   assert( blkmem != NULL );
   assert( bufmem != NULL );

   SCIPdebugMessage("Calling SCIPsdpiCreate \n");

   BMS_CALL( BMSallocBlockMemory(blkmem, sdpisolver) );

   (*sdpisolver)->messagehdlr = messagehdlr;
   (*sdpisolver)->blkmem = blkmem;
   (*sdpisolver)->bufmem = bufmem;

   MOSEK_CALLM( MSK_makeenv(&((*sdpisolver)->mskenv), NULL) );/*lint !e641*/ /* the NULL-argument is a debug file, but setting this will spam the whole folder */

   /* this will be properly initialized then calling solve */
   (*sdpisolver)->msktask = NULL;

   (*sdpisolver)->nvars = 0;
   (*sdpisolver)->nactivevars = 0;
   (*sdpisolver)->inputtomosekmapper = NULL;
   (*sdpisolver)->mosektoinputmapper = NULL;
   (*sdpisolver)->fixedvarsval = NULL;
   (*sdpisolver)->fixedvarsobjcontr = 0.0;
   (*sdpisolver)->objcoefs = NULL;
   (*sdpisolver)->nvarbounds = 0;
   (*sdpisolver)->varboundpos = NULL;
   (*sdpisolver)->solved = FALSE;
   (*sdpisolver)->sdpcounter = 0;

   (*sdpisolver)->epsilon = 1e-9;
   (*sdpisolver)->gaptol = 1e-4;
   (*sdpisolver)->feastol = 1e-6;
   (*sdpisolver)->sdpsolverfeastol = 1e-6;
   (*sdpisolver)->objlimit = SCIPsdpiSolverInfinity(*sdpisolver);
   (*sdpisolver)->sdpinfo = FALSE;
   (*sdpisolver)->nthreads = -1;
   (*sdpisolver)->timelimit = FALSE;
   (*sdpisolver)->timelimitinitial = FALSE;

   return SCIP_OKAY;
}

/** deletes an SDP solver interface */
SCIP_RETCODE SCIPsdpiSolverFree(
   SCIP_SDPISOLVER**     sdpisolver          /**< pointer to an SDP interface solver structure */
   )
{
   assert( sdpisolver != NULL );
   assert( *sdpisolver != NULL );

   SCIPdebugMessage("Freeing SDPISolver\n");

   if ( ((*sdpisolver)->msktask) != NULL )
   {
      MOSEK_CALL( MSK_deletetask(&((*sdpisolver)->msktask)) );/*lint !e641*/
   }

   if ( ((*sdpisolver)->mskenv) != NULL )
   {
      MOSEK_CALL( MSK_deleteenv(&((*sdpisolver)->mskenv)) );/*lint !e641*/
   }

   BMSfreeBlockMemoryArrayNull((*sdpisolver)->blkmem, &(*sdpisolver)->varboundpos, 2 * (*sdpisolver)->nactivevars); /*lint !e647*/
   BMSfreeBlockMemoryArrayNull((*sdpisolver)->blkmem, &(*sdpisolver)->inputtomosekmapper, (*sdpisolver)->nvars);/*lint !e737*/
   BMSfreeBlockMemoryArrayNull((*sdpisolver)->blkmem, &(*sdpisolver)->mosektoinputmapper, (*sdpisolver)->nactivevars);/*lint !e737*/
   BMSfreeBlockMemoryArrayNull((*sdpisolver)->blkmem, &(*sdpisolver)->fixedvarsval, (*sdpisolver)->nvars - (*sdpisolver)->nactivevars); /*lint !e776*/
   BMSfreeBlockMemoryArrayNull((*sdpisolver)->blkmem, &(*sdpisolver)->objcoefs, (*sdpisolver)->nactivevars); /*lint !e776*/

   BMSfreeBlockMemory((*sdpisolver)->blkmem, sdpisolver);

   return SCIP_OKAY;
}

/** increases the SDP-Counter */
SCIP_RETCODE SCIPsdpiSolverIncreaseCounter(
   SCIP_SDPISOLVER*      sdpisolver          /**< SDP interface solver structure */
   )
{
   assert( sdpisolver != NULL );

   sdpisolver->sdpcounter++;

   return SCIP_OKAY;
}

/** reset the SDP-Counter to zero */
SCIP_RETCODE SCIPsdpiSolverResetCounter(
   SCIP_SDPISOLVER*      sdpisolver          /**< SDP interface solver structure */
   )
{
   assert( sdpisolver != NULL );

   SCIPdebugMessage("Resetting counter of SDP-Interface from %d to 0.\n", sdpisolver->sdpcounter);
   sdpisolver->sdpcounter = 0;

   return SCIP_OKAY;
}

/**@} */


/*
 * Solving Methods
 */

/**@name Solving Methods */
/**@{ */

/** loads and solves an SDP
 *
 *  For the non-constant SDP- and the LP-part, the original arrays before fixings should be given, for the constant
 *  SDP-part the arrays AFTER fixings should be given. In addition, an array needs to be given, that for every block and
 *  every row/col index within that block either has value -1, meaning that this index should be deleted, or a
 *  non-negative integer stating the number of indices before it that are to be deleated, meaning that this index will
 *  be decreased by that number, in addition to that the total number of deleted indices for each block should be given.
 *  Optionally an array start may be given with a starting point for the solver (if this is NULL then the solver should
 *  start from scratch).
 *
 *  @warning Depending on the solver, the given lp arrays might get sorted in their original position.
 */
SCIP_RETCODE SCIPsdpiSolverLoadAndSolve(
   SCIP_SDPISOLVER*      sdpisolver,         /**< SDP interface solver structure */
   int                   nvars,              /**< number of variables */
   SCIP_Real*            obj,                /**< objective function values of variables */
   SCIP_Real*            lb,                 /**< lower bounds of variables */
   SCIP_Real*            ub,                 /**< upper bounds of variables */
   int                   nsdpblocks,         /**< number of SDP-blocks */
   int*                  sdpblocksizes,      /**< sizes of the SDP-blocks (may be NULL if nsdpblocks = sdpconstnnonz = sdpnnonz = 0) */
   int*                  sdpnblockvars,      /**< number of variables that exist in each block */
   int                   sdpconstnnonz,      /**< number of nonzero elements in the constant matrices of the SDP-blocks AFTER FIXINGS */
   int*                  sdpconstnblocknonz, /**< number of nonzeros for each variable in the constant part, also the i-th entry gives the
                                              *   number of entries  of sdpconst row/col/val [i] AFTER FIXINGS */
   int**                 sdpconstrow,        /**< pointers to row-indices for each block AFTER FIXINGS*/
   int**                 sdpconstcol,        /**< pointers to column-indices for each block AFTER FIXINGS */
   SCIP_Real**           sdpconstval,        /**< pointers to the values of the nonzeros for each block AFTER FIXINGS */
   int                   sdpnnonz,           /**< number of nonzero elements in the SDP-constraint matrix */
   int**                 sdpnblockvarnonz,   /**< entry [i][j] gives the number of nonzeros for block i and variable j, this is exactly
                                              *   the number of entries of sdp row/col/val [i][j] */
   int**                 sdpvar,             /**< sdpvar[i][j] gives the sdp-index of the j-th variable (according to the sorting for row/col/val)
                                              *   in the i-th block */
   int***                sdprow,             /**< pointer to the row-indices for each block and variable */
   int***                sdpcol,             /**< pointer to the column-indices for each block and variable */
   SCIP_Real***          sdpval,             /**< values of SDP-constraint matrix entries (may be NULL if sdpnnonz = 0) */
   int**                 indchanges,         /**< changes needed to be done to the indices, if indchanges[block][nonz]=-1, then the index can
                                              *   be removed, otherwise it gives the number of indices removed before this */
   int*                  nremovedinds,       /**< the number of rows/cols to be fixed for each block */
   int*                  blockindchanges,    /**< block indizes will be modified by these, see indchanges */
   int                   nremovedblocks,     /**< number of empty blocks that should be removed */
   int                   nlpcons,            /**< number of active (at least two nonzeros) LP-constraints */
   int                   noldlpcons,         /**< number of LP-constraints including those with less than two active nonzeros */
   SCIP_Real*            lplhs,              /**< left hand sides of active LP rows after fixings (may be NULL if nlpcons = 0) */
   SCIP_Real*            lprhs,              /**< right hand sides of active LP rows after fixings (may be NULL if nlpcons = 0) */
   int*                  rownactivevars,     /**< number of active variables for each LP constraint */
   int                   lpnnonz,            /**< number of nonzero elements in the LP-constraint matrix */
   int*                  lprow,              /**< row-index for each entry in lpval-array, might get sorted (may be NULL if lpnnonz = 0) */
   int*                  lpcol,              /**< column-index for each entry in lpval-array, might get sorted (may be NULL if lpnnonz = 0) */
   SCIP_Real*            lpval,              /**< values of LP-constraint matrix entries, might get sorted (may be NULL if lpnnonz = 0) */
   SCIP_Real*            start,              /**< NULL or a starting point for the solver, this should have length nvars */
   SCIP_SDPSOLVERSETTING startsettings,      /**< settings used to start with in SDPA, currently not used for DSDP and MOSEK, set this to
                                               *  SCIP_SDPSOLVERSETTING_UNSOLVED to ignore it and start from scratch */
   SCIP_Real             timelimit           /**< after this many seconds solving will be aborted (currently only implemented for DSDP and MOSEK) */
   )
{
   return SCIPsdpiSolverLoadAndSolveWithPenalty(sdpisolver, 0.0, TRUE, FALSE, nvars, obj, lb, ub, nsdpblocks, sdpblocksizes, sdpnblockvars, sdpconstnnonz,
               sdpconstnblocknonz, sdpconstrow, sdpconstcol, sdpconstval, sdpnnonz, sdpnblockvarnonz, sdpvar, sdprow, sdpcol, sdpval, indchanges,
               nremovedinds, blockindchanges, nremovedblocks, nlpcons, noldlpcons, lplhs, lprhs, rownactivevars, lpnnonz, lprow, lpcol, lpval, start,
               startsettings, timelimit, NULL, NULL);
}

/** loads and solves an SDP using a penalty formulation
 *
 *  The penalty formulation of the SDP is:
 *      \f{eqnarray*}{
 *      \min & & b^T y + \Gamma r \\
 *      \mbox{s.t.} & & \sum_{j=1}^n A_j^i y_j - A_0^i + r \cdot \mathbb{I} \succeq 0 \quad \forall i \leq m \\
 *      & & Dy + r \cdot \mathbb{I} \geq d \\
 *      & & l \leq y \leq u \\
 *      & & r \geq 0.\f}
 *  Alternatively withobj can be set to false to set b to 0 and only check for feasibility (if the optimal objective value is
 *  bigger than 0 the problem is infeasible, otherwise it's feasible), and rbound can be set to false to remove the non-negativity condition on r.
 *  For the non-constant SDP- and the LP-part the original arrays before fixings should be given, for the constant SDP-part the arrays AFTER fixings
 *  should be given. In addition, an array needs to be given, that for every block and every row/col index within that block either has value
 *  -1, meaning that this index should be deleted, or a non-negative integer stating the number of indices before it that are to be deleated,
 *  meaning that this index will be decreased by that number. Moreover, the total number of deleted indices for each block should be given.
 *  An optional starting point for the solver may be given; if it is NULL, the solver will start from scratch.
 *
 *  @warning Depending on the solver, the given lp arrays might get sorted in their original position.
 */
SCIP_RETCODE SCIPsdpiSolverLoadAndSolveWithPenalty(
   SCIP_SDPISOLVER*      sdpisolver,         /**< SDP interface solver structure */
   SCIP_Real             penaltyparam,       /**< the Gamma above, needs to be >= 0 */
   SCIP_Bool             withobj,            /**< if this is false the objective is set to 0 */
   SCIP_Bool             rbound,             /**< should r be non-negative ? */
   int                   nvars,              /**< number of variables */
   SCIP_Real*            obj,                /**< objective function values of variables */
   SCIP_Real*            lb,                 /**< lower bounds of variables */
   SCIP_Real*            ub,                 /**< upper bounds of variables */
   int                   nsdpblocks,         /**< number of SDP-blocks */
   int*                  sdpblocksizes,      /**< sizes of the SDP-blocks (may be NULL if nsdpblocks = sdpconstnnonz = sdpnnonz = 0) */
   int*                  sdpnblockvars,      /**< number of variables that exist in each block */
   int                   sdpconstnnonz,      /**< number of nonzero elements in the constant matrices of the SDP-Blocks AFTER FIXINGS */
   int*                  sdpconstnblocknonz, /**< number of nonzeros for each variable in the constant part, also the i-th entry gives the
                                              *   number of entries  of sdpconst row/col/val [i] AFTER FIXINGS */
   int**                 sdpconstrow,        /**< pointers to row-indices for each block AFTER FIXINGS */
   int**                 sdpconstcol,        /**< pointers to column-indices for each block AFTER FIXINGS */
   SCIP_Real**           sdpconstval,        /**< pointers to the values of the nonzeros for each block AFTER FIXINGS */
   int                   sdpnnonz,           /**< number of nonzero elements in the SDP-constraint matrix */
   int**                 sdpnblockvarnonz,   /**< entry [i][j] gives the number of nonzeros for block i and variable j, this is exactly
                                              *   the number of entries of sdp row/col/val [i][j] */
   int**                 sdpvar,             /**< sdpvar[i][j] gives the sdp-index of the j-th variable (according to the sorting for row/col/val)
                                              *   in the i-th block */
   int***                sdprow,             /**< pointer to the row-indices for each block and variable */
   int***                sdpcol,             /**< pointer to the column-indices for each block and variable */
   SCIP_Real***          sdpval,             /**< values of SDP-constraint matrix entries (may be NULL if sdpnnonz = 0) */
   int**                 indchanges,         /**< changes needed to be done to the indices, if indchanges[block][nonz]=-1, then
                                              *   the index can be removed, otherwise it gives the number of indices removed before this */
   int*                  nremovedinds,       /**< the number of rows/cols to be fixed for each block */
   int*                  blockindchanges,    /**< block indizes will be modified by these, see indchanges */
   int                   nremovedblocks,     /**< number of empty blocks that should be removed */
   int                   nlpcons,            /**< number of active (at least two nonzeros) LP-constraints */
   int                   noldlpcons,         /**< number of LP-constraints including those with less than two active nonzeros */
   SCIP_Real*            lplhs,              /**< left hand sides of active LP rows after fixings (may be NULL if nlpcons = 0) */
   SCIP_Real*            lprhs,              /**< right hand sides of active LP rows after fixings (may be NULL if nlpcons = 0) */
   int*                  rownactivevars,     /**< number of active variables for each LP constraint */
   int                   lpnnonz,            /**< number of nonzero elements in the LP-constraint matrix */
   int*                  lprow,              /**< row-index for each entry in lpval-array, might get sorted (may be NULL if lpnnonz = 0) */
   int*                  lpcol,              /**< column-index for each entry in lpval-array, might get sorted (may be NULL if lpnnonz = 0) */
   SCIP_Real*            lpval,              /**< values of LP-constraint matrix entries, might get sorted (may be NULL if lpnnonz = 0) */
   SCIP_Real*            start,              /**< NULL or a starting point for the solver, this should have length nvars */
   SCIP_SDPSOLVERSETTING startsettings,      /**< settings used to start with in SDPA, currently not used for DSDP and MOSEK, set this to
                                               *  SCIP_SDPSOLVERSETTING_UNSOLVED to ignore it and start from scratch */
   SCIP_Real             timelimit,          /**< after this many seconds solving will be aborted (currently only implemented for DSDP and MOSEK) */
   SCIP_Bool*            feasorig,           /**< pointer to store if the solution to the penalty-formulation is feasible for the original problem
                                               *  (may be NULL if penaltyparam = 0) */
   SCIP_Bool*            penaltybound        /**< pointer to store if the primal solution reached the bound Tr(X) <= penaltyparam in the primal problem,
                                               *  this is also an indication of the penalty parameter being to small (may be NULL if not needed) */
)
{/*lint --e{715,818}*/
   int b;
   int i;
   int j;
   int blockvar;
   int v;
   int k;
   long long ind;
   int mosekind = 0;
   SCIP_Real* mosekvarbounds;
   int nfixedvars;
   int oldnactivevars;
   int* vartorowmapper; /* maps the lpvars to the corresponding left- and right-hand-sides of the LP constraints */
   int* vartolhsrhsmapper; /* maps the lpvars to the corresponding entries in lplhs and lprhs */
   int nlpvars;
   int pos;
   int newpos;
   int* mosekblocksizes;
   SCIP_Real one; /* MOSEK always wants a pointer to factors for a sum of matrices, we always use a single matrix with factor one */
   int* mosekrow;
   int* mosekcol;
   SCIP_Real* mosekval;
   int row;
   SCIP_Real val;
   struct timeval starttime;
   struct timeval currenttime;
   SCIP_Real startseconds;
   SCIP_Real currentseconds;
   SCIP_Real elapsedtime;
#ifdef SCIP_MORE_DEBUG
   int nmosekconss;
   int nmosekvars;
   int nmosekcones;
   char varname[SCIP_MAXSTRLEN];
#endif
#if CONVERT_ABSOLUTE_TOLERANCES
   SCIP_Real maxrhscoef; /* MOSEK uses a relative feasibility tolerance, the largest rhs-coefficient is needed for converting the absolute tolerance */
#endif

   assert( sdpisolver != NULL );
   assert( sdpisolver->mskenv != NULL );
   assert( penaltyparam > -1 * sdpisolver->epsilon );
   assert( penaltyparam < sdpisolver->epsilon || ( feasorig != NULL ) );
   assert( nvars > 0 );
   assert( obj != NULL );
   assert( lb != NULL );
   assert( ub != NULL );
   assert( nsdpblocks >= 0 );
   assert( nsdpblocks == 0 || sdpblocksizes != NULL );
   assert( nsdpblocks == 0 || sdpnblockvars != NULL );
   assert( sdpconstnnonz >= 0 );
   assert( nsdpblocks == 0 || sdpconstnnonz == 0 || sdpconstnblocknonz != NULL );
   assert( nsdpblocks == 0 || sdpconstnnonz == 0 || sdpconstrow != NULL );
   assert( nsdpblocks == 0 || sdpconstnnonz == 0 || sdpconstcol != NULL );
   assert( nsdpblocks == 0 || sdpconstnnonz == 0 || sdpconstval != NULL );
   assert( sdpnnonz >= 0 );
   assert( nsdpblocks == 0 || sdpnblockvarnonz != NULL );
   assert( nsdpblocks == 0 || sdpvar != NULL );
   assert( nsdpblocks == 0 || sdprow != NULL );
   assert( nsdpblocks == 0 || sdpcol != NULL );
   assert( nsdpblocks == 0 || sdpval != NULL );
   assert( nsdpblocks == 0 || indchanges != NULL );
   assert( nsdpblocks == 0 || nremovedinds != NULL );
   assert( nsdpblocks == 0 || blockindchanges != NULL );
   assert( 0 <= nremovedblocks && nremovedblocks <= nsdpblocks );
   assert( nlpcons >= 0 );
   assert( noldlpcons >= nlpcons );
   assert( nlpcons == 0 || lplhs != NULL );
   assert( nlpcons == 0 || lprhs != NULL );
   assert( nlpcons == 0 || rownactivevars != NULL );
   assert( lpnnonz >= 0 );
   assert( nlpcons == 0 || lprow != NULL );
   assert( nlpcons == 0 || lpcol != NULL );
   assert( nlpcons == 0 || lpval != NULL );

   /* check the timelimit */
   if ( timelimit <= 0.0 )
   {
      sdpisolver->timelimit = TRUE;
      sdpisolver->timelimitinitial = TRUE;
      sdpisolver->solved = FALSE;
      return SCIP_OKAY;
   }
   sdpisolver->timelimit = FALSE;
   sdpisolver->timelimitinitial = FALSE;
   sdpisolver->feasorig = FALSE;

   /* start the timing */
   TIMEOFDAY_CALL( gettimeofday(&starttime, NULL) );/*lint !e438, !e550, !e641 */

   one = 1.0;
#if CONVERT_ABSOLUTE_TOLERANCES
   maxrhscoef = 0.0;
#endif

   /* create an empty task (second and third argument are guesses for maximum number of constraints and variables), if there already is one, delete it */
   if ((sdpisolver->msktask) != NULL)
   {
      MOSEK_CALL( MSK_deletetask(&(sdpisolver->msktask)) );/*lint !e641*/
   }
   if ( penaltyparam < sdpisolver->epsilon )
   {
      MOSEK_CALLM( MSK_maketask(sdpisolver->mskenv, nvars, nsdpblocks - nremovedblocks + nlpcons + 2 * nvars, &(sdpisolver->msktask)) );/*lint !e641*/
   }
   else
   {
      MOSEK_CALLM( MSK_maketask(sdpisolver->mskenv, nvars + 1, nsdpblocks - nremovedblocks + nlpcons + 2 * nvars, &(sdpisolver->msktask)) );/*lint !e641*/
   }

#ifdef SCIP_MORE_DEBUG
   MOSEK_CALL( MSK_linkfunctotaskstream (sdpisolver->msktask, MSK_STREAM_LOG, NULL, printstr) ); /* output to console */
#else
   /* if sdpinfo is true, redirect output to console */
   if ( sdpisolver->sdpinfo )
   {
      MOSEK_CALL( MSK_linkfunctotaskstream (sdpisolver->msktask, MSK_STREAM_LOG, NULL, printstr) );/*lint !e641*/
   }
#endif

   /* set number of threads */
   if ( sdpisolver->nthreads > 0 )
   {
      MOSEK_CALL( MSK_putintparam(sdpisolver->msktask, MSK_IPAR_NUM_THREADS, sdpisolver->nthreads) );/*lint !e641*/
   }

   /* only increase the counter if we don't use the penalty formulation to stay in line with the numbers in the general interface (where this is still the
    * same SDP) */
   if ( penaltyparam < sdpisolver->epsilon )
      SCIPdebugMessage("Inserting Data into MOSEK for SDP (%d) \n", ++sdpisolver->sdpcounter);
   else
      SCIPdebugMessage("Inserting Data again into MOSEK for SDP (%d) \n", sdpisolver->sdpcounter);

   /* set the penalty and rbound flags accordingly */
   sdpisolver->penalty = (penaltyparam < sdpisolver->epsilon) ? FALSE : TRUE;
   sdpisolver->rbound = rbound;

   /* allocate memory for inputtomosekmapper, mosektoinputmapper and the fixed and active variable information, for the latter this will
    * later be shrinked if the needed size is known */
   BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &(sdpisolver->inputtomosekmapper), sdpisolver->nvars, nvars) );
   BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &(sdpisolver->mosektoinputmapper), sdpisolver->nactivevars, nvars) );
   BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &(sdpisolver->fixedvarsval), sdpisolver->nvars - sdpisolver->nactivevars, nvars) ); /*lint !e776*/
   BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &(sdpisolver->objcoefs), sdpisolver->nactivevars, nvars) ); /*lint !e776*/

   oldnactivevars = sdpisolver->nactivevars; /* we need to save this to realloc the varboundpos-array if needed */
   sdpisolver->nvars = nvars;
   sdpisolver->nactivevars = 0;
   nfixedvars = 0;

   /* find the fixed variables */
   sdpisolver->fixedvarsobjcontr = 0.0;
   for (i = 0; i < nvars; i++)
   {
      if ( isFixed(sdpisolver, lb[i], ub[i]) )
      {
         sdpisolver->fixedvarsobjcontr += obj[i] * lb[i]; /* this is the value this fixed variable contributes to the objective */
         sdpisolver->fixedvarsval[nfixedvars] = lb[i]; /* if lb=ub, then this is the value the variable will have in every solution */
         nfixedvars++;
         sdpisolver->inputtomosekmapper[i] = -nfixedvars;
         SCIPdebugMessage("Fixing variable %d locally to %f for SDP %d in MOSEK\n", i, lb[i], sdpisolver->sdpcounter);
      }
      else
      {
         sdpisolver->mosektoinputmapper[sdpisolver->nactivevars] = i;
         sdpisolver->inputtomosekmapper[i] = sdpisolver->nactivevars;
         sdpisolver->objcoefs[sdpisolver->nactivevars] = obj[i];
         sdpisolver->nactivevars++;
#ifdef SCIP_MORE_DEBUG
         SCIPdebugMessage("Variable %d becomes variable %d for SDP %d in MOSEK\n", i, sdpisolver->inputtomosekmapper[i], sdpisolver->sdpcounter);
#endif
      }
   }
   assert( sdpisolver->nactivevars + nfixedvars == sdpisolver->nvars );

   /* if we want to solve without objective, we reset fixedvarsobjcontr */
   if ( ! withobj )
      sdpisolver->fixedvarsobjcontr = 0.0;

   /* shrink the fixedvars, objcoefs and mosektoinputmapper arrays to the right size */
   BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &(sdpisolver->objcoefs), nvars, sdpisolver->nactivevars) );
   BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &(sdpisolver->fixedvarsval), nvars, nfixedvars) );
   BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &(sdpisolver->mosektoinputmapper), nvars, sdpisolver->nactivevars) );

   /* compute number of variable bounds and save them in mosekvarbounds */
   sdpisolver->nvarbounds = 0;
   BMS_CALL( BMSallocBufferMemoryArray(sdpisolver->bufmem, &mosekvarbounds, 2 * sdpisolver->nactivevars) ); /*lint !e647*/

   if ( sdpisolver->nactivevars != oldnactivevars )
   {
      if ( sdpisolver->varboundpos == NULL )
      {
         BMS_CALL( BMSallocBlockMemoryArray(sdpisolver->blkmem, &(sdpisolver->varboundpos), 2 * sdpisolver->nactivevars) ); /*lint !e647*/
      }
      else
      {
         BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &(sdpisolver->varboundpos), 2 * oldnactivevars, 2 * sdpisolver->nactivevars) ); /*lint !e647*/
      }
   }
   assert( sdpisolver->varboundpos != NULL );

   for (i = 0; i < sdpisolver->nactivevars; i++)
   {
      assert( 0 <= sdpisolver->mosektoinputmapper[i] && sdpisolver->mosektoinputmapper[i] < nvars );
      if ( ! SCIPsdpiSolverIsInfinity(sdpisolver, lb[sdpisolver->mosektoinputmapper[i]]) )
      {
         mosekvarbounds[sdpisolver->nvarbounds] = lb[sdpisolver->mosektoinputmapper[i]];
         sdpisolver->varboundpos[sdpisolver->nvarbounds] = -(i + 1); /* negative sign means lower bound */
         (sdpisolver->nvarbounds)++;
      }
      if ( ! SCIPsdpiSolverIsInfinity(sdpisolver, ub[sdpisolver->mosektoinputmapper[i]]) )
      {
         mosekvarbounds[sdpisolver->nvarbounds] = -1 * ub[sdpisolver->mosektoinputmapper[i]]; /* we give the upper bounds a negative sign for the objective */
         sdpisolver->varboundpos[sdpisolver->nvarbounds] = +(i + 1); /* positive sign means upper bound */
         (sdpisolver->nvarbounds)++;
      }
   }

   if ( nlpcons > 0 )
   {
      /* allocate memory to save which lpvariable corresponds to which lp constraint, negative signs correspond to left-hand-sides of lp constraints,
       * entry i or -i corresponds to the constraint in position |i|-1, as we have to add +1 to make the entries strictly positive or strictly negative */
      BMS_CALL( BMSallocBufferMemoryArray(sdpisolver->bufmem, &vartorowmapper, 2*noldlpcons) ); /*lint !e647*/ /*lint !e530*/
      /* allocate memory to save at which indices the corresponding lhss and rhss of the lpvars are saved */
      BMS_CALL( BMSallocBufferMemoryArray(sdpisolver->bufmem, &vartolhsrhsmapper, 2*noldlpcons) ); /*lint !e647*/ /*lint !e530*/

      /* compute the number of LP constraints after splitting the ranged rows and compute the rowtovarmapper */
      pos = 0;
      newpos = 0; /* the position in the lhs and rhs arrays */
      for (i = 0; i < noldlpcons; i++)
      {
         assert( newpos <= nlpcons );
         if ( rownactivevars[i] >= 2 )
         {
            if ( lplhs[newpos] > - SCIPsdpiSolverInfinity(sdpisolver) )
            {
               vartorowmapper[pos] = -(i+1);
               vartolhsrhsmapper[pos] = newpos;
               pos++;

#if CONVERT_ABSOLUTE_TOLERANCES
               /* update largest rhs-entry */
               if ( REALABS(lplhs[newpos]) > maxrhscoef )
                  maxrhscoef = REALABS(lplhs[newpos]);
#endif

            }
            if ( lprhs[newpos] < SCIPsdpiSolverInfinity(sdpisolver) )
            {
               vartorowmapper[pos] = i+1;
               vartolhsrhsmapper[pos] = newpos;
               pos++;

#if CONVERT_ABSOLUTE_TOLERANCES
               /* update largest rhs-entry */
               if ( REALABS(lprhs[newpos]) > maxrhscoef )
                  maxrhscoef = REALABS(lprhs[newpos]);
#endif
            }
            newpos++;
         }
      }
      nlpvars = pos;
      assert( nlpvars <= 2*nlpcons ); /* *2 comes from left- and right-hand-sides */
   }
   else
      nlpvars = 0;

   /* create matrix variables */
   BMS_CALL( BMSallocBufferMemoryArray(sdpisolver->bufmem, &mosekblocksizes, nsdpblocks - nremovedblocks) ); /*lint !e679 !e776*/

   for (b = 0; b < nsdpblocks; b++)
   {
      if ( blockindchanges[b] > -1 )
      {
         assert( 0 <= blockindchanges[b] && blockindchanges[b] <= b && (b - blockindchanges[b]) <= (nsdpblocks - nremovedblocks) );
         mosekblocksizes[b - blockindchanges[b]] = sdpblocksizes[b] - nremovedinds[b];/*lint !e679*/
      }
   }
   MOSEK_CALLM( MSK_appendbarvars(sdpisolver->msktask, nsdpblocks - nremovedblocks, mosekblocksizes) );/*lint !e641*/

   /* create scalar variables (since we solve the primal problem, these are not the active variables but the dual variables for the
    * lp constraints and variable bounds) */
   MOSEK_CALLM( MSK_appendvars(sdpisolver->msktask, nlpvars + sdpisolver->nvarbounds) );/*lint !e641*/

   /* all of these are non-negative */
   for (v = 0; v < nlpvars + sdpisolver->nvarbounds; v++)
   {
      MOSEK_CALL( MSK_putvarbound(sdpisolver->msktask, v, MSK_BK_LO, 0.0, (double) MSK_DPAR_DATA_TOL_BOUND_INF) );/*lint !e641*/
   }

   /* append empty constraints (since we solve the primal problem, we get one constraint for each active variable) */
   MOSEK_CALLM( MSK_appendcons(sdpisolver->msktask, (penaltyparam < sdpisolver->epsilon) ? sdpisolver->nactivevars : sdpisolver->nactivevars + 1) );/*lint !e641*/

   /* set objective values for the matrix variables */
   i = 0;

   if ( sdpconstnnonz > 0 )
   {
      for (b = 0; b < nsdpblocks; b++)
      {
         if ( blockindchanges[b] > -1 )
         {
            /* if some indices in the block were removed, we need to change indices accordingly */
            if ( nremovedinds[b] > 0 )
            {
               int* moseksdpconstrow;
               int* moseksdpconstcol;

               BMS_CALL( BMSallocBufferMemoryArray(sdpisolver->bufmem, &moseksdpconstrow, sdpconstnblocknonz[b]) );
               BMS_CALL( BMSallocBufferMemoryArray(sdpisolver->bufmem, &moseksdpconstcol, sdpconstnblocknonz[b]) );

               for (k = 0; k < sdpconstnblocknonz[b]; k++)
               {
                  /* rows and cols with active nonzeros should not be removed */
                  assert( -1 < indchanges[b][sdpconstrow[b][k]] && indchanges[b][sdpconstrow[b][k]] <= sdpconstrow[b][k] );
                  assert( -1 < indchanges[b][sdpconstcol[b][k]] && indchanges[b][sdpconstcol[b][k]] <= sdpconstcol[b][k] );

                  assert( 0 <= sdpconstrow[b][k] && sdpconstrow[b][k] <= sdpblocksizes[b] );
                  assert( 0 <= sdpconstcol[b][k] && sdpconstcol[b][k] <= sdpblocksizes[b] );

                  moseksdpconstrow[k] = sdpconstrow[b][k] - indchanges[b][sdpconstrow[b][k]];
                  moseksdpconstcol[k] = sdpconstcol[b][k] - indchanges[b][sdpconstcol[b][k]];

#if CONVERT_ABSOLUTE_TOLERANCES
                  /* update largest rhs-entry */
                  if ( REALABS(sdpconstval[b][k]) > maxrhscoef )
                     maxrhscoef = REALABS(sdpconstval[b][k]);
#endif
               }

               MOSEK_CALL( MSK_appendsparsesymmat(sdpisolver->msktask, mosekblocksizes[b - blockindchanges[b]], sdpconstnblocknonz[b],
                     moseksdpconstrow, moseksdpconstcol, sdpconstval[b], &ind) );/*lint !e641, !e679, !e747*/

               BMSfreeBufferMemoryArray(sdpisolver->bufmem, &moseksdpconstcol);
               BMSfreeBufferMemoryArray(sdpisolver->bufmem, &moseksdpconstrow);
            }
            else
            {
               MOSEK_CALL( MSK_appendsparsesymmat(sdpisolver->msktask, mosekblocksizes[b - blockindchanges[b]], sdpconstnblocknonz[b],
                     sdpconstrow[b], sdpconstcol[b], sdpconstval[b], &ind) );/*lint !e641, !e679, !e747*/
            }
            MOSEK_CALL( MSK_putbarcj(sdpisolver->msktask, i, 1, &ind, &one) );/*lint !e641, !e747*/
            i++;
         }
      }
   }

   /* set objective values for the scalar variables */
   /* first for those corresponding to LP constraints in the dual */
   for (i = 0; i < nlpvars; i++)
   {
      if ( vartorowmapper[i] > 0 )/*lint !e644*/ /* right-hand side */
      {
         MOSEK_CALL( MSK_putcj(sdpisolver->msktask, i, -1 * lprhs[vartolhsrhsmapper[i]]) );/*lint !e641, !e644*/
#ifdef SCIP_MORE_DEBUG
         /* give the variable a meaningful name for debug output */
         (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "rhs-%d", vartorowmapper[i] - 1);
         MOSEK_CALL( MSK_putvarname ( sdpisolver->msktask, i, varname) );
#endif
      }
      else /* left-hand side */
      {
         assert( vartorowmapper[i] < 0 ); /* we should not have value 0 so that we can clearly differentiate between positive and negative */
         MOSEK_CALL( MSK_putcj(sdpisolver->msktask, i, lplhs[vartolhsrhsmapper[i]]) );/*lint !e641*/
#ifdef SCIP_MORE_DEBUG
         /* give the variable a meaningful name for debug output */
         (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "lhs-%d", -1 * vartorowmapper[i] - 1);
         MOSEK_CALL( MSK_putvarname ( sdpisolver->msktask, i, varname) );
#endif
      }
   }

   /* finally for those corresponding to variable bounds in the dual */
   for (i = 0; i < sdpisolver->nvarbounds; i++)
   {
      MOSEK_CALL( MSK_putcj(sdpisolver->msktask, nlpvars + i, mosekvarbounds[i]) );/*lint !e641*/ /* for the ub's we already added a negative sign in mosekvarbounds*/
#ifdef SCIP_MORE_DEBUG
      if ( sdpisolver->varboundpos[i] < 0 ) /* lower bound */
      {
         /* give the variable a meaningful name for debug output */
         (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "lb-%d", sdpisolver->mosektoinputmapper[-1 * sdpisolver->varboundpos[i] - 1]);
         MOSEK_CALL( MSK_putvarname ( sdpisolver->msktask, nlpvars + i, varname) );
      }
      else /* upper bound */
      {
         assert( sdpisolver->varboundpos[i] > 0 ); /* we should not have value 0 so that we can clearly differentiate between positive and negative */
         /* give the variable a meaningful name for debug output */
         (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "ub-%d", sdpisolver->mosektoinputmapper[sdpisolver->varboundpos[i] - 1]);
         MOSEK_CALL( MSK_putvarname ( sdpisolver->msktask, nlpvars + i, varname) );
      }
#endif
   }

   BMSfreeBufferMemoryArray(sdpisolver->bufmem, &mosekvarbounds);

   /* set objective sense (since we want to minimize in the dual, we maximize in the primal) */
   MOSEK_CALL( MSK_putobjsense(sdpisolver->msktask, MSK_OBJECTIVE_SENSE_MAXIMIZE) );/*lint !e641*/

   /* start inserting the constraints */

   /* first add the matrices A_i */
   for (b = 0; b < nsdpblocks; b++)
   {
      if ( blockindchanges[b] > -1 )
      {
         for (blockvar = 0; blockvar < sdpnblockvars[b]; blockvar++)
         {
            v = sdpisolver->inputtomosekmapper[sdpvar[b][blockvar]];

            /* check if the variable is active */
            if ( v > -1 )
            {
               assert( v < sdpisolver->nactivevars );
               /* if there are removed indices, we have to adjust the column and row indices accordingly */
               if ( nremovedinds[b] > 0 )
               {
                  BMS_CALL( BMSallocBufferMemoryArray(sdpisolver->bufmem, &mosekrow, sdpnblockvarnonz[b][blockvar]) );
                  BMS_CALL( BMSallocBufferMemoryArray(sdpisolver->bufmem, &mosekcol, sdpnblockvarnonz[b][blockvar]) );

                  for (k = 0; k < sdpnblockvarnonz[b][blockvar]; k++)
                  {
                     /* rows and cols with active nonzeros should not be removed */
                     assert( -1 < indchanges[b][sdprow[b][blockvar][k]] && indchanges[b][sdprow[b][blockvar][k]] <= sdprow[b][blockvar][k] );
                     assert( -1 < indchanges[b][sdpcol[b][blockvar][k]] && indchanges[b][sdpcol[b][blockvar][k]] <= sdpcol[b][blockvar][k] );

                     assert( 0 <= sdprow[b][blockvar][k] && sdprow[b][blockvar][k] < sdpblocksizes[b] );
                     assert( 0 <= sdpcol[b][blockvar][k] && sdpcol[b][blockvar][k] < sdpblocksizes[b] );

                     mosekrow[k] = sdprow[b][blockvar][k] - indchanges[b][sdprow[b][blockvar][k]];
                     mosekcol[k] = sdpcol[b][blockvar][k] - indchanges[b][sdpcol[b][blockvar][k]];
                  }

                  MOSEK_CALL( MSK_appendsparsesymmat(sdpisolver->msktask, mosekblocksizes[b - blockindchanges[b]], (long long) sdpnblockvarnonz[b][blockvar],
                        mosekrow, mosekcol, sdpval[b][blockvar], &ind) );/*lint !e641, !e679*/

                  BMSfreeBufferMemoryArray(sdpisolver->bufmem, &mosekcol);
                  BMSfreeBufferMemoryArray(sdpisolver->bufmem, &mosekrow);
               }
               else
               {
                  MOSEK_CALL( MSK_appendsparsesymmat(sdpisolver->msktask, mosekblocksizes[b - blockindchanges[b]], (long long) sdpnblockvarnonz[b][blockvar],
                        sdprow[b][blockvar], sdpcol[b][blockvar], sdpval[b][blockvar], &ind) );/*lint !e641, !e679*/
               }

               MOSEK_CALL( MSK_putbaraij(sdpisolver->msktask, v, b - blockindchanges[b], (long long) 1, &ind, &one) );/*lint !e641*/
            }
         }
      }
   }

   /* add the identity matrix corresponding to the penalty variable */
   if ( penaltyparam >= sdpisolver->epsilon )
   {
      int* identityindices;
      SCIP_Real* identityvalues;

      for (b = 0; b < nsdpblocks; b++)
      {
         if ( blockindchanges[b] > -1 )
         {
            BMS_CALL( BMSallocBufferMemoryArray(sdpisolver->bufmem, &identityindices, mosekblocksizes[b - blockindchanges[b]]) );
            BMS_CALL( BMSallocBufferMemoryArray(sdpisolver->bufmem, &identityvalues, mosekblocksizes[b - blockindchanges[b]]) );

            for (i = 0; i < mosekblocksizes[b - blockindchanges[b]]; i++)
            {
               identityindices[i] = i;
               identityvalues[i] = 1.0;
            }
            MOSEK_CALL( MSK_appendsparsesymmat(sdpisolver->msktask, mosekblocksizes[b - blockindchanges[b]], (long long) mosekblocksizes[b - blockindchanges[b]],
                                    identityindices, identityindices, identityvalues, &ind) );/*lint !e641, !e679*/
            MOSEK_CALL( MSK_putbaraij(sdpisolver->msktask, sdpisolver->nactivevars, b - blockindchanges[b], (long long) 1, &ind, &one) );/*lint !e641, !e679*/

            BMSfreeBufferMemoryArray(sdpisolver->bufmem, &identityvalues);
            BMSfreeBufferMemoryArray(sdpisolver->bufmem, &identityindices);
         }
      }
   }

   /* now add the entries corresponding to the lp-constraints in the dual problem */
   if ( penaltyparam < sdpisolver->epsilon )
   {
      BMS_CALL( BMSallocBufferMemoryArray(sdpisolver->bufmem, &mosekrow, lpnnonz) );
      BMS_CALL( BMSallocBufferMemoryArray(sdpisolver->bufmem, &mosekval, lpnnonz) );
   }
   else
   {
      /* one extra entry for the penalty-constraint Trace = Gamma */
      BMS_CALL( BMSallocBufferMemoryArray(sdpisolver->bufmem, &mosekrow, lpnnonz + 1) );/*lint !e776*/
      BMS_CALL( BMSallocBufferMemoryArray(sdpisolver->bufmem, &mosekval, lpnnonz + 1) );/*lint !e776*/
   }

   ind = 0;
   for (i = 0; i < nlpvars; i++)
   {
      if ( vartorowmapper[i] < 0 ) /* left-hand side */
      {
         mosekind = 0;
         /* find the first lp-entry belonging to this variable (those in between have to belong to constraints with less than two active variables and
          * will therefore not be used) */
         while (ind < lpnnonz && lprow[ind] < -1 * vartorowmapper[i] - 1)
            ind++;
         /* iterate over all nonzeros to input them to the array given to MOSEK */
         while (ind < lpnnonz && lprow[ind] == -1 * vartorowmapper[i] - 1) /* they should already be sorted by rows in the sdpi */
         {
            v = sdpisolver->inputtomosekmapper[lpcol[ind]];
            if ( v > -1 )
            {
               mosekrow[mosekind] = v;
               mosekval[mosekind] = lpval[ind];
               mosekind++;
            }
            ind++;
         }
         assert( mosekind <= lpnnonz );
      }
      else /* right-hand side */
      {
         assert( vartorowmapper[i] > 0 ); /* we should not have value 0 so that we can clearly differentiate between positive and negative */

         if ( i > 0 && vartorowmapper[i] == -1 * vartorowmapper[i - 1] )
         {
            /* we already iterated over this constraint in the lp-arrays, so we keep the current ind position and as we currenlty have
             * the corresponding entries in the mosek-arrays, we iterate over them again just changing the signs (except for the penalty-entry) */
            for (j = 0; j < (penaltyparam < sdpisolver->epsilon ? mosekind : mosekind - 1); j++)/*lint !e644*/
               mosekval[j] *= -1;
         }
         else
         {
            /* no left hand side for this constraint exists, so we didnot yet iterate over this constraint in the lp arrays */
            mosekind = 0;
            /* find the first lp-entry belonging to this variable (those in between have to belong to constraints with less than two active variables and
             * will therefore not be used) */
            while (lprow[ind] < vartorowmapper[i] - 1)
               ind++;
            while (ind < lpnnonz && lprow[ind] == vartorowmapper[i] - 1)
            {
               v = sdpisolver->inputtomosekmapper[lpcol[ind]];
               if ( v > -1 )
               {
                  mosekrow[mosekind] = v;
                  mosekval[mosekind] = -1 * lpval[ind]; /* because we need to change the <= to a >= constraint */
                  mosekind++;
               }
               ind++;
            }
            assert( mosekind <= lpnnonz );
         }
      }

      /* add the additional entry for the penalty constraint Trace = Gamma */
      if ( penaltyparam >= sdpisolver->epsilon )
      {
         /* check if we reset mosekind, in case we did not and just "copied" the lhs-entries for the rhs, we do not need to reset the extra entry,
          * since it is already there */
         if ( ! (i > 0 && vartorowmapper[i] == -1 * vartorowmapper[i - 1] ))
         {
            mosekrow[mosekind] = sdpisolver->nactivevars;
            mosekval[mosekind] = 1.0;
            mosekind++;
         }
         assert( mosekind <= lpnnonz + 1 );
      }

      MOSEK_CALL( MSK_putacol(sdpisolver->msktask, i, mosekind, mosekrow, mosekval) );/*lint !e641*/
   }

   BMSfreeBufferMemoryArrayNull(sdpisolver->bufmem, &mosekval);
   BMSfreeBufferMemoryArrayNull(sdpisolver->bufmem, &mosekrow);

   /* finally add the entries corresponding to varbounds in the dual problem, we get exactly one entry per variable */
   for (i = 0; i < sdpisolver->nvarbounds; i++)
   {
      if ( sdpisolver->varboundpos[i] < 0 )
      {
         /* lower bound */
         row =-1 * sdpisolver->varboundpos[i] - 1; /* minus one because we added one to get strictly positive/negative values */
         val = 1.0;
      }
      else
      {
         /* upper bound */
         assert( sdpisolver->varboundpos[i] > 0 ); /* we should not have a zero as we wanted a clear differentiation between positive and negative */
         row = sdpisolver->varboundpos[i] - 1; /* minus one because we added one to get strictly positive/negative values */
         val = -1.0;
      }
      MOSEK_CALL( MSK_putacol(sdpisolver->msktask, nlpvars + i, 1, &row, &val) );/*lint !e641*/
   }

   /* make all constraints equality constraints with right-hand side b_i (or 0 if we solve without objective) */
   for (i = 0; i < sdpisolver->nactivevars; i++)
   {
      if ( withobj )
      {
         MOSEK_CALL( MSK_putconbound(sdpisolver->msktask, i, MSK_BK_FX, obj[sdpisolver->mosektoinputmapper[i]], obj[sdpisolver->mosektoinputmapper[i]]) );/*lint !e641*/
      }
      else
      {
         MOSEK_CALL( MSK_putconbound(sdpisolver->msktask, i, MSK_BK_FX, 0.0, 0.0) );/*lint !e641*/
      }
#ifdef SCIP_MORE_DEBUG
         /* give the constraint a meaningful name for debug output */
         (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "var-%d", sdpisolver->mosektoinputmapper[i]);
         MOSEK_CALL( MSK_putconname ( sdpisolver->msktask, i, varname) );
#endif
   }

   /* the penalty constraint has right-hand side Gamma, it is a <=-inequality if r was bounded and an equality constraint otherwise */
   if ( penaltyparam >= sdpisolver->epsilon )
   {
      if ( rbound )
      {
         MOSEK_CALL( MSK_putconbound(sdpisolver->msktask, sdpisolver->nactivevars, MSK_BK_UP, (double) -1 * MSK_DPAR_DATA_TOL_BOUND_INF, penaltyparam) );/*lint !e641*/
      }
      else
      {
         MOSEK_CALL( MSK_putconbound(sdpisolver->msktask, sdpisolver->nactivevars, MSK_BK_FX, penaltyparam, penaltyparam) );/*lint !e641*/
      }
   }

   /* write to file if asked to */
#ifdef SCIP_DEBUG_PRINTTOFILE
   SCIP_CALL( SCIPsdpiSolverWriteSDP(sdpisolver, "mosek.task") );
#endif

   /* print whole problem and parameters if asked to */
#ifdef SCIP_MORE_DEBUG
   MOSEK_CALL( MSK_getnumcon (sdpisolver->msktask, &nmosekconss) );
   MOSEK_CALL( MSK_getnumvar (sdpisolver->msktask, &nmosekvars) );
   MOSEK_CALL( MSK_getnumcone (sdpisolver->msktask, &nmosekcones) );

   MOSEK_CALL( MSK_printdata (sdpisolver->msktask, MSK_STREAM_LOG, 0, nmosekconss, 0, nmosekvars, 0, nmosekcones, 1, 1, 1, 1, 1, 1, 1, 1) );
#ifdef SCIP_PRINT_PARAMETERS
   MOSEK_CALL( MSK_printparam (sdpisolver->msktask) );
#endif
#endif

   /* set the time limit */
   startseconds = (SCIP_Real) starttime.tv_sec + (SCIP_Real) starttime.tv_usec / 1e6;

   TIMEOFDAY_CALL( gettimeofday(&currenttime, NULL) );/*lint !e438, !e550, !e641 */
   currentseconds = (SCIP_Real) currenttime.tv_sec + (SCIP_Real) currenttime.tv_usec / 1e6;

   elapsedtime = currentseconds - startseconds;

   if ( timelimit <= elapsedtime )
   {
      sdpisolver->timelimit = TRUE;
      sdpisolver->solved = FALSE;
   }
   else
   {
      SCIP_Real feastol;

      /* set epsilon and feasibility tolerance (note that the dual in MOSEK is the problem we are interested in, so this is where we use feastol,
       * since MOSEK works with relative tolerance, we adjust our absolute tolerance accordingly, so that any solution satisfying the relative
       * tolerance in MOSEK satisfies our absolute tolerance) */
      MOSEK_CALL( MSK_putdouparam(sdpisolver->msktask, MSK_DPAR_INTPNT_CO_TOL_PFEAS, sdpisolver->gaptol) );/*lint !e641*/
#if CONVERT_ABSOLUTE_TOLERANCES
      MOSEK_CALL( MSK_putdouparam(sdpisolver->msktask, MSK_DPAR_INTPNT_CO_TOL_DFEAS, sdpisolver->sdpsolverfeastol / (1 + maxrhscoef)) );/*lint !e641*/
      MOSEK_CALL( MSK_putdouparam(sdpisolver->msktask, MSK_DPAR_INTPNT_CO_TOL_INFEAS, sdpisolver->sdpsolverfeastol / (1 + maxrhscoef)) );/*lint !e641*/
      SCIPdebugMessage("Setting relative feasibility tolerance for MOSEK to %.10f / %f = %.12f\n", sdpisolver->sdpsolverfeastol,
            1+maxrhscoef, sdpisolver->sdpsolverfeastol / (1 + maxrhscoef));
#else
      MOSEK_CALL( MSK_putdouparam(sdpisolver->msktask, MSK_DPAR_INTPNT_CO_TOL_DFEAS, sdpisolver->sdpsolverfeastol) );/*lint !e641*/
      MOSEK_CALL( MSK_putdouparam(sdpisolver->msktask, MSK_DPAR_INTPNT_CO_TOL_INFEAS, sdpisolver->sdpsolverfeastol) );/*lint !e641*/
#endif
      MOSEK_CALL( MSK_putdouparam(sdpisolver->msktask, MSK_DPAR_INTPNT_CO_TOL_MU_RED, sdpisolver->gaptol) );/*lint !e641*/
      MOSEK_CALL( MSK_putdouparam(sdpisolver->msktask, MSK_DPAR_INTPNT_CO_TOL_REL_GAP, sdpisolver->gaptol) );/*lint !e641*/

      if ( ! SCIPsdpiSolverIsInfinity(sdpisolver, timelimit - elapsedtime) )
      {
         MOSEK_CALL( MSK_putdouparam(sdpisolver->msktask, MSK_DPAR_OPTIMIZER_MAX_TIME, timelimit - elapsedtime) );/*lint !e641*/
      }

      /* set objective cutoff */
      if ( ! SCIPsdpiSolverIsInfinity(sdpisolver, sdpisolver->objlimit) )
      {
         MOSEK_CALL( MSK_putdouparam(sdpisolver->msktask, MSK_DPAR_UPPER_OBJ_CUT, sdpisolver->objlimit) );/*lint !e641*/
      }

      /* solve the problem */
      MOSEK_CALL( MSK_optimizetrm(sdpisolver->msktask, &(sdpisolver->terminationcode)) );/*lint !e641*/

      if ( sdpisolver->sdpinfo )
      {
         MOSEK_CALL( MSK_optimizersummary(sdpisolver->msktask, MSK_STREAM_LOG) );/*lint !e641*/
         MOSEK_CALL( MSK_analyzesolution(sdpisolver->msktask, MSK_STREAM_LOG, MSK_SOL_ITR) );/*lint !e641*/
      }

      SCIPdebugMessage("Solving problem using MOSEK, return code %d\n", sdpisolver->terminationcode);

      sdpisolver->solved = TRUE;

      sdpisolver->nsdpcalls = 1;
      MOSEK_CALL( MSK_getnaintinf(sdpisolver->msktask, "MSK_IINF_INTPNT_ITER", &(sdpisolver->niterations)) );/*lint !e641*/

      /* if the problem has been stably solved but did not reach the required feasibility tolerance, even though the solver
       * reports feasibility, resolve it with adjusted tolerance */
#if CONVERT_ABSOLUTE_TOLERANCES
      feastol = sdpisolver->sdpsolverfeastol / (1 + maxrhscoef);
#else
      feastol = sdpisolver->sdpsolverfeastol;
#endif

      while ( SCIPsdpiSolverIsAcceptable(sdpisolver) && SCIPsdpiSolverIsDualFeasible(sdpisolver) && penaltyparam < sdpisolver->epsilon && feastol >= INFEASMINFEASTOL )
      {
         SCIP_Real* solvector;
         int nvarspointer;
         SCIP_Bool infeasible;
         int newiterations;

         /* get current solution */
         BMS_CALL( BMSallocBufferMemoryArray(sdpisolver->bufmem, &solvector, nvars) );
         nvarspointer = nvars;
         SCIP_CALL( SCIPsdpiSolverGetSol(sdpisolver, NULL, solvector, &nvarspointer) );
         assert( nvarspointer == nvars );

         /* check the solution for feasibility with regards to our tolerance */
         SCIP_CALL( SCIPsdpSolcheckerCheck(sdpisolver->bufmem, nvars, lb, ub, nsdpblocks, sdpblocksizes, sdpnblockvars, sdpconstnnonz,
               sdpconstnblocknonz, sdpconstrow, sdpconstcol, sdpconstval, sdpnnonz, sdpnblockvarnonz, sdpvar, sdprow, sdpcol, sdpval,
               indchanges, nremovedinds, blockindchanges, nlpcons, noldlpcons, lplhs, lprhs, rownactivevars, lpnnonz, lprow, lpcol, lpval,
               solvector, sdpisolver->feastol, sdpisolver->epsilon, &infeasible) );

         BMSfreeBufferMemoryArray(sdpisolver->bufmem, &solvector);

         if ( infeasible )
         {
            SCIPdebugMessage("Solution feasible for MOSEK but outside feasibility tolerance, changing MOSEK feasibility tolerance from %f to %f\n",
                  feastol, feastol * INFEASFEASTOLCHANGE);
            feastol *= INFEASFEASTOLCHANGE;

            if ( feastol >= INFEASMINFEASTOL )
            {
               /* update settings */
               MOSEK_CALL( MSK_putdouparam(sdpisolver->msktask, MSK_DPAR_INTPNT_CO_TOL_DFEAS, feastol) );/*lint !e641*/
               MOSEK_CALL( MSK_putdouparam(sdpisolver->msktask, MSK_DPAR_INTPNT_CO_TOL_INFEAS, feastol) );/*lint !e641*/

               /* set the time limit */
               startseconds = (SCIP_Real) starttime.tv_sec + (SCIP_Real) starttime.tv_usec / 1e6;

               TIMEOFDAY_CALL( gettimeofday(&currenttime, NULL) );/*lint !e438, !e550, !e641 */
               currentseconds = (SCIP_Real) currenttime.tv_sec + (SCIP_Real) currenttime.tv_usec / 1e6;

               elapsedtime = currentseconds - startseconds;

               if ( timelimit <= elapsedtime )
               {
                  sdpisolver->timelimit = TRUE;
                  sdpisolver->solved = FALSE;
               }

               /* solve the problem */
               MOSEK_CALL( MSK_optimizetrm(sdpisolver->msktask, &(sdpisolver->terminationcode)) );/*lint !e641*/

               if ( sdpisolver->sdpinfo )
               {
                  MOSEK_CALL( MSK_optimizersummary(sdpisolver->msktask, MSK_STREAM_LOG) );/*lint !e641*/
                  MOSEK_CALL( MSK_analyzesolution(sdpisolver->msktask, MSK_STREAM_LOG, MSK_SOL_ITR) );/*lint !e641*/
               }

               /* update number of SDP-iterations and -calls */
               sdpisolver->nsdpcalls++;
               MOSEK_CALL( MSK_getnaintinf(sdpisolver->msktask, "MSK_IINF_INTPNT_ITER", &newiterations) );/*lint !e641*/
               sdpisolver->niterations += newiterations;
            }
            else
            {
               sdpisolver->solved = FALSE;
               SCIPmessagePrintInfo(sdpisolver->messagehdlr, "MOSEK failed to reach required feasibility tolerance! \n");
            }
         }
         else
            break;
      }
   }


   /* if using a penalty formulation, check if the solution is feasible for the original problem */
   if ( penaltyparam >= sdpisolver->epsilon && ( ! sdpisolver->timelimit ) && ( sdpisolver->terminationcode != MSK_RES_TRM_MAX_TIME ) )
   {
      SCIP_Real* moseksol;
      SCIP_Real trace = 0.0;
      SCIP_Real* x;

      assert( feasorig != NULL );

      /* get the r variable in the dual problem */
      BMSallocBufferMemoryArray(sdpisolver->bufmem, &moseksol, sdpisolver->nactivevars + 1);/*lint !e776*/

      MOSEK_CALL( MSK_gety(sdpisolver->msktask, MSK_SOL_ITR, moseksol) );/*lint !e641*/

      *feasorig = (moseksol[sdpisolver->nactivevars] < sdpisolver->feastol); /*lint !e413*/

      /* only set sdpisolver->feasorig to true if we solved with objective, because only in this case we want to compute
       * the objective value by hand since it is numerically more stable then the result returned by DSDP */
      if ( withobj )
         sdpisolver->feasorig = *feasorig;

      /* if r > 0 also check the primal bound */
      if ( ! *feasorig && penaltybound != NULL )
      {

         SCIPdebugMessage("Solution not feasible in original problem, r = %f\n", moseksol[sdpisolver->nactivevars]);

         /* compute Tr(X) */

         /* start with the diagonal entries of the primal semidefinite variables */
         for (b = 0; b < nsdpblocks; b++)
         {
            if ( blockindchanges[b] > -1 )
            {
               SCIP_Real* X; /* the upper triangular entries of matrix X */
               int size;

               size = sdpblocksizes[b] - nremovedinds[b];

               BMS_CALL( BMSallocBufferMemoryArray(sdpisolver->bufmem, &X, 0.5 * size * (size + 1)) );
               MOSEK_CALL( MSK_getbarxj(sdpisolver->msktask, MSK_SOL_ITR, b - blockindchanges[b], X) );/*lint !e641*/

               /* iterate over all diagonal entries */
               for (i = 0; i < size; i++)
               {
                  /* get index in the lower triangular part */
                  ind = 0.5 * i * (i + 3);/*lint !e776*/ /*  i*(i+1)/2 + i  */
                  assert( ind < 0.5 * size * (size + 1) );
                  trace += X[ind];
               }

               BMSfreeBufferMemoryArray(sdpisolver->bufmem, &X);
            }
         }

         /* add primal lp-variables */
         BMS_CALL( BMSallocBufferMemoryArray(sdpisolver->bufmem, &x, nlpvars + sdpisolver->nvarbounds) );

         MOSEK_CALL( MSK_getxx(sdpisolver->msktask, MSK_SOL_ITR, x) );/*lint !e641*/

         for (i = 0; i < nlpvars; i++)
            trace += x[i];

         BMSfreeBufferMemoryArrayNull(sdpisolver->bufmem, &x);

         /* if the relative gap is smaller than the tolerance, we return equality */
         if ( (penaltyparam - trace) / penaltyparam < PENALTYBOUNDTOL )/*lint !e414*/
         {
            if ( penaltybound != NULL )
               *penaltybound = TRUE;
            SCIPdebugMessage("Tr(X) = %f == %f = Gamma, penalty formulation not exact, Gamma should be increased or problem is infeasible\n",
                  trace, penaltyparam);
         }
         else if ( penaltybound != NULL )
            *penaltybound = FALSE;
      }
      BMSfreeBufferMemoryArray(sdpisolver->bufmem, &moseksol);
   }

   /* free memory */
   BMSfreeBufferMemoryArrayNull(sdpisolver->bufmem, &mosekblocksizes);
   if ( nlpcons > 0 )
   {
      BMSfreeBufferMemoryArrayNull(sdpisolver->bufmem, &vartolhsrhsmapper);
      BMSfreeBufferMemoryArrayNull(sdpisolver->bufmem, &vartorowmapper);
   }

   return SCIP_OKAY;
}
/**@} */




/*
 * Solution Information Methods
 */

/**@name Solution Information Methods */
/**@{ */

/** returns whether a solve method was called after the last modification of the SDP */
SCIP_Bool SCIPsdpiSolverWasSolved(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{/*lint --e{818}*/
   assert( sdpisolver != NULL );

   return sdpisolver->solved;
}

/** returns true if the solver could determine whether the problem is feasible
 *
 *  So it returns true if the solver knows that the problem is feasible/infeasible/unbounded, it returns false if the
 *  solver does not know anything about the feasibility status and thus the functions IsPrimalFeasible etc. should not be
 *  used.
 */
SCIP_Bool SCIPsdpiSolverFeasibilityKnown(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{/*lint --e{818}*/
   MSKsolstae solstat;

   assert( sdpisolver != NULL );
   CHECK_IF_SOLVED_BOOL( sdpisolver );

   MOSEK_CALL_BOOL( MSK_getsolsta(sdpisolver->msktask, MSK_SOL_ITR, &solstat) );/*lint !e641*/

   switch ( solstat )
   {
   case MSK_SOL_STA_UNKNOWN:
   case MSK_SOL_STA_PRIM_FEAS:
   case MSK_SOL_STA_DUAL_FEAS:
   case MSK_SOL_STA_NEAR_OPTIMAL:
   case MSK_SOL_STA_NEAR_PRIM_FEAS:
   case MSK_SOL_STA_NEAR_DUAL_FEAS:
   case MSK_SOL_STA_NEAR_PRIM_AND_DUAL_FEAS:
   case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:
   case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER:
      return FALSE;
   case MSK_SOL_STA_OPTIMAL:
   case MSK_SOL_STA_PRIM_AND_DUAL_FEAS:
   case MSK_SOL_STA_PRIM_INFEAS_CER:
   case MSK_SOL_STA_DUAL_INFEAS_CER:
      return TRUE;
   default:
      SCIPdebugMessage("Unknown return code in SCIPsdpiSolverFeasibilityKnown\n");
      return FALSE;
   }/*lint !e788*/
}

/** gets information about primal and dual feasibility of the current SDP solution */
SCIP_RETCODE SCIPsdpiSolverGetSolFeasibility(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   SCIP_Bool*            primalfeasible,     /**< stores primal feasibility status */
   SCIP_Bool*            dualfeasible        /**< stores dual feasibility status */
   )
{/*lint --e{818}*/
   MSKsolstae solstat;

   assert( sdpisolver != NULL );
   assert( primalfeasible != NULL );
   assert( dualfeasible != NULL );
   CHECK_IF_SOLVED( sdpisolver );

   MOSEK_CALL( MSK_getsolsta(sdpisolver->msktask, MSK_SOL_ITR, &solstat) );/*lint !e641*/

   switch ( solstat )
   {
   case MSK_SOL_STA_OPTIMAL:
   case MSK_SOL_STA_PRIM_AND_DUAL_FEAS:
      *primalfeasible = TRUE;
      *dualfeasible = TRUE;
      break;
   case MSK_SOL_STA_PRIM_INFEAS_CER:
      *primalfeasible = FALSE;
      *dualfeasible = TRUE;
      break;
   case MSK_SOL_STA_DUAL_INFEAS_CER:
      *primalfeasible = TRUE;
      *dualfeasible = FALSE;
      break;
   default:
      SCIPdebugMessage("MOSEK does not know about feasibility of solutions\n");
      return SCIP_LPERROR;
   }/*lint !e788*/

   return SCIP_OKAY;
}

/** returns TRUE iff SDP is proven to be primal unbounded,
 *  returns FALSE with a debug-message if the solver could not determine feasibility */
SCIP_Bool SCIPsdpiSolverIsPrimalUnbounded(
   SCIP_SDPISOLVER*      sdpisolver          /**< SDP interface solver structure */
   )
{/*lint --e{818}*/
   MSKsolstae solstat;

   assert( sdpisolver != NULL );
   CHECK_IF_SOLVED_BOOL( sdpisolver );

   MOSEK_CALL_BOOL( MSK_getsolsta(sdpisolver->msktask, MSK_SOL_ITR, &solstat) );/*lint !e641*/

   switch ( solstat )
   {
   case MSK_SOL_STA_DUAL_INFEAS_CER:
      return TRUE;
   case MSK_SOL_STA_OPTIMAL:
   case MSK_SOL_STA_PRIM_AND_DUAL_FEAS:
   case MSK_SOL_STA_PRIM_INFEAS_CER:
      break;
   default:
      SCIPdebugMessage("MOSEK does not know about feasibility of solutions\n");
      break;
   }/*lint !e788*/
   return FALSE;
}

/** returns TRUE iff SDP is proven to be primal infeasible,
 *  returns FALSE with a debug-message if the solver could not determine feasibility */
SCIP_Bool SCIPsdpiSolverIsPrimalInfeasible(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{/*lint --e{818}*/
   MSKsolstae solstat;

   assert( sdpisolver != NULL );
   CHECK_IF_SOLVED_BOOL( sdpisolver );

   MOSEK_CALL_BOOL( MSK_getsolsta(sdpisolver->msktask, MSK_SOL_ITR, &solstat) );/*lint !e641*/

   switch ( solstat )
   {
   case MSK_SOL_STA_PRIM_INFEAS_CER:
      return TRUE;
   case MSK_SOL_STA_OPTIMAL:
   case MSK_SOL_STA_PRIM_AND_DUAL_FEAS:
   case MSK_SOL_STA_DUAL_INFEAS_CER:
      break;
   default:
      SCIPdebugMessage("MOSEK does not know about feasibility of solutions\n");
      break;
   }/*lint !e788*/
   return FALSE;
}

/** returns TRUE iff SDP is proven to be primal feasible,
 *  returns FALSE with a debug-message if the solver could not determine feasibility */
SCIP_Bool SCIPsdpiSolverIsPrimalFeasible(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{/*lint --e{818}*/
   MSKsolstae solstat;

   assert( sdpisolver != NULL );
   CHECK_IF_SOLVED_BOOL( sdpisolver );

   MOSEK_CALL_BOOL( MSK_getsolsta(sdpisolver->msktask, MSK_SOL_ITR, &solstat) );/*lint !e641*/

   switch ( solstat )
   {
   case MSK_SOL_STA_OPTIMAL:
   case MSK_SOL_STA_PRIM_AND_DUAL_FEAS:
   case MSK_SOL_STA_DUAL_INFEAS_CER:
      return TRUE;
   case MSK_SOL_STA_PRIM_INFEAS_CER:
      break;
   default:
      SCIPdebugMessage("MOSEK does not know about feasibility of solutions\n");
      break;
   }/*lint !e788*/
   return FALSE;
}

/** returns TRUE iff SDP is proven to be dual unbounded,
 *  returns FALSE with a debug-message if the solver could not determine feasibility */
SCIP_Bool SCIPsdpiSolverIsDualUnbounded(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{/*lint --e{818}*/
   MSKsolstae solstat;

   assert( sdpisolver != NULL );
   CHECK_IF_SOLVED_BOOL( sdpisolver );

   MOSEK_CALL_BOOL( MSK_getsolsta(sdpisolver->msktask, MSK_SOL_ITR, &solstat) );/*lint !e641*/

   switch ( solstat )
   {
   case MSK_SOL_STA_PRIM_INFEAS_CER:
      return TRUE;
   case MSK_SOL_STA_OPTIMAL:
   case MSK_SOL_STA_PRIM_AND_DUAL_FEAS:
   case MSK_SOL_STA_DUAL_INFEAS_CER:
      break;
   default:
      SCIPdebugMessage("MOSEK does not know about feasibility of solutions\n");
      break;
   }/*lint !e788*/
   return FALSE;
}

/** returns TRUE iff SDP is proven to be dual infeasible,
 *  returns FALSE with a debug-message if the solver could not determine feasibility */
SCIP_Bool SCIPsdpiSolverIsDualInfeasible(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{/*lint --e{818}*/
   MSKsolstae solstat;

   assert( sdpisolver != NULL );
   CHECK_IF_SOLVED_BOOL( sdpisolver );

   MOSEK_CALL_BOOL( MSK_getsolsta(sdpisolver->msktask, MSK_SOL_ITR, &solstat) );/*lint !e641*/

   switch ( solstat )
   {
   case MSK_SOL_STA_DUAL_INFEAS_CER:
      return TRUE;
   case MSK_SOL_STA_OPTIMAL:
   case MSK_SOL_STA_PRIM_AND_DUAL_FEAS:
   case MSK_SOL_STA_PRIM_INFEAS_CER:
      break;
   default:
      SCIPdebugMessage("MOSEK does not know about feasibility of solutions\n");
      break;
   }/*lint !e788*/
   return FALSE;
}

/** returns TRUE iff SDP is proven to be dual feasible,
 *  returns FALSE with a debug-message if the solver could not determine feasibility */
SCIP_Bool SCIPsdpiSolverIsDualFeasible(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{/*lint --e{818}*/
   MSKsolstae solstat;

   assert( sdpisolver != NULL );
   CHECK_IF_SOLVED_BOOL( sdpisolver );

   MOSEK_CALL_BOOL( MSK_getsolsta(sdpisolver->msktask, MSK_SOL_ITR, &solstat) );/*lint !e641*/

   switch ( solstat )
   {
   case MSK_SOL_STA_OPTIMAL:
   case MSK_SOL_STA_PRIM_AND_DUAL_FEAS:
   case MSK_SOL_STA_PRIM_INFEAS_CER:
      return TRUE;
   case MSK_SOL_STA_DUAL_INFEAS_CER:
      break;
   default:
      SCIPdebugMessage("MOSEK does not know about feasibility of solutions\n");
      break;
   }/*lint !e788*/
   return FALSE;
}

/** returns TRUE iff the solver converged */
SCIP_Bool SCIPsdpiSolverIsConverged(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{/*lint --e{818}*/
   assert( sdpisolver != NULL );

   if ( sdpisolver->timelimit )
      return FALSE;

   CHECK_IF_SOLVED_BOOL( sdpisolver );

   return sdpisolver->terminationcode == MSK_RES_OK;
}

/** returns TRUE iff the objective limit was reached */
SCIP_Bool SCIPsdpiSolverIsObjlimExc(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{/*lint --e{818}*/
   assert( sdpisolver != NULL );
   CHECK_IF_SOLVED_BOOL( sdpisolver );

   return sdpisolver->terminationcode == MSK_RES_TRM_OBJECTIVE_RANGE;
}

/** returns TRUE iff the iteration limit was reached */
SCIP_Bool SCIPsdpiSolverIsIterlimExc(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{/*lint --e{818}*/
   assert( sdpisolver != NULL );
   CHECK_IF_SOLVED_BOOL( sdpisolver );

   return sdpisolver->terminationcode == MSK_RES_TRM_MAX_ITERATIONS;
}

/** returns TRUE iff the time limit was reached */
SCIP_Bool SCIPsdpiSolverIsTimelimExc(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{/*lint --e{818}*/
   assert( sdpisolver != NULL );

   if ( sdpisolver->timelimit )
      return TRUE;

   if ( ! sdpisolver->solved )
      return FALSE;

   return sdpisolver->terminationcode == MSK_RES_TRM_MAX_TIME;
}

/** returns the internal solution status of the solver, which has the following meaning:<br>
 * -1: solver was not started<br>
 *  0: converged<br>
 *  1: infeasible start<br>
 *  2: numerical problems<br>
 *  3: objective limit reached<br>
 *  4: iteration limit reached<br>
 *  5: time limit reached<br>
 *  6: user termination<br>
 *  7: other */
int SCIPsdpiSolverGetInternalStatus(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{/*lint --e{818}*/
   assert( sdpisolver != NULL );

   if ( ! sdpisolver->solved )
      return -1;

   if ( sdpisolver->timelimit )
      return 5;

   switch ( sdpisolver->terminationcode )
   {
   case MSK_RES_OK:
      return 0;
   case MSK_RES_TRM_MAX_NUM_SETBACKS:
   case MSK_RES_TRM_NUMERICAL_PROBLEM:
   case MSK_RES_TRM_STALL:
      return 2;
   case MSK_RES_TRM_OBJECTIVE_RANGE:
      return 3;
   case MSK_RES_TRM_MAX_ITERATIONS:
      return 4;
   case MSK_RES_TRM_MAX_TIME:
      return 5;
   default:
      return 7;
   }/*lint !e788*/
}

/** returns TRUE iff SDP was solved to optimality, meaning the solver converged and returned primal and dual feasible solutions */
SCIP_Bool SCIPsdpiSolverIsOptimal(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{/*lint --e{818}*/
   MSKsolstae solstat;

   assert( sdpisolver != NULL );

   if ( sdpisolver->timelimit )
      return FALSE;

   CHECK_IF_SOLVED_BOOL( sdpisolver );

   if ( sdpisolver->terminationcode != MSK_RES_OK )
      return FALSE;

   MOSEK_CALL_BOOL( MSK_getsolsta(sdpisolver->msktask, MSK_SOL_ITR, &solstat) );/*lint !e641*/

   if ( solstat != MSK_SOL_STA_OPTIMAL )
      return FALSE;

   return TRUE;
}

/** returns TRUE iff SDP was solved to optimality or some other status was reached
 *  that is still acceptable inside a Branch & Bound framework */
SCIP_Bool SCIPsdpiSolverIsAcceptable(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{/*lint --e{818}*/
   assert( sdpisolver != NULL );

   if ( sdpisolver->timelimit )
      return FALSE;

   if ( ! sdpisolver->solved )
      return FALSE;

   return SCIPsdpiSolverIsConverged(sdpisolver) && SCIPsdpiSolverFeasibilityKnown(sdpisolver);
}

/** tries to reset the internal status of the SDP-solver in order to ignore an instability of the last solving call */
SCIP_RETCODE SCIPsdpiSolverIgnoreInstability(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   SCIP_Bool*            success             /**< pointer to store, whether the instability could be ignored */
   )
{/*lint --e{715,818}*/
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_LPERROR;
}

/** gets objective value of solution */
SCIP_RETCODE SCIPsdpiSolverGetObjval(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   SCIP_Real*            objval              /**< pointer to store the objective value */
   )
{/*lint --e{818}*/
   SCIP_Real* moseksol;

   assert( sdpisolver != NULL );
   CHECK_IF_SOLVED( sdpisolver );
   assert( objval != NULL );

   if ( sdpisolver->penalty && ( ! sdpisolver->feasorig ) )
   {
      /* in this case we cannot really trust the solution given by MOSEK, since changes in the value of r much less than epsilon can
       * cause huge changes in the objective, so using the objective value given by MOSEK is numerically more stable */
      MOSEK_CALL( MSK_getdualobj(sdpisolver->msktask, MSK_SOL_ITR, objval) );
   }
   else
   {
      int v;

      /* since the objective value given by MOSEK sometimes differs slightly from the correct value for the given solution,
       * we get the solution from MOSEK and compute the correct objective value */
      BMSallocBufferMemoryArray(sdpisolver->bufmem, &moseksol, sdpisolver->penalty ? sdpisolver->nactivevars + 1 : sdpisolver->nactivevars);
      MOSEK_CALL( MSK_gety(sdpisolver->msktask, MSK_SOL_ITR, moseksol) );/*lint !e641*/

      *objval = 0.0;
      for (v = 0; v < sdpisolver->nactivevars; v++)
      {
         if ( moseksol[v] > sdpisolver->epsilon )
            *objval += moseksol[v] * sdpisolver->objcoefs[v];
      }
   }

   /* as we didn't add the fixed (lb = ub) variables to MOSEK, we have to add their contributions to the objective as well */
   *objval += sdpisolver->fixedvarsobjcontr;

   if ( ( ! sdpisolver->penalty ) || sdpisolver->feasorig)
   {
      BMSfreeBufferMemoryArray(sdpisolver->bufmem, &moseksol);
   }

   return SCIP_OKAY;
}

/** gets dual solution vector for feasible SDPs
 *
 *  If dualsollength isn't equal to the number of variables this will return the needed length and a debug message is thrown.
 */
SCIP_RETCODE SCIPsdpiSolverGetSol(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   SCIP_Real*            objval,             /**< pointer to store the objective value, may be NULL if not needed */
   SCIP_Real*            dualsol,            /**< pointer to store the dual solution vector, may be NULL if not needed */
   int*                  dualsollength       /**< length of the dual sol vector, must be 0 if dualsol is NULL, if this is less than the number
                                              *   of variables in the SDP, a DebugMessage will be thrown and this is set to the needed value */
   )
{/*lint --e{818}*/
   int v;
   SCIP_Real* moseksol;

   assert( sdpisolver != NULL );
   CHECK_IF_SOLVED( sdpisolver );
   assert( dualsollength != NULL );

   if ( *dualsollength > 0 )
   {
      assert( dualsol != NULL );
      if ( *dualsollength < sdpisolver->nvars )
      {
         SCIPdebugMessage("The given array in SCIPsdpiSolverGetSol only had length %d, but %d was needed", *dualsollength, sdpisolver->nvars);
         *dualsollength = sdpisolver->nvars;

         return SCIP_OKAY;
      }

      BMSallocBufferMemoryArray(sdpisolver->bufmem, &moseksol, sdpisolver->penalty ? sdpisolver->nactivevars + 1 : sdpisolver->nactivevars);

      MOSEK_CALL( MSK_gety(sdpisolver->msktask, MSK_SOL_ITR, moseksol) );/*lint !e641*/

      /* insert the entries into dualsol, for non-fixed vars we copy those from MOSEK, the rest are the saved entries from inserting (they equal lb=ub) */
      for (v = 0; v < sdpisolver->nvars; v++)
      {
         if ( sdpisolver->inputtomosekmapper[v] >= 0 )
            dualsol[v] = moseksol[sdpisolver->inputtomosekmapper[v]];
         else
         {
            /* this is the value that was saved when inserting, as this variable has lb=ub */
            assert( -sdpisolver->inputtomosekmapper[v] <= sdpisolver->nvars - sdpisolver->nactivevars );
            dualsol[v] = sdpisolver->fixedvarsval[(-1 * sdpisolver->inputtomosekmapper[v]) - 1]; /*lint !e679*/ /* -1 because we wanted strictly negative vals */
         }
      }

      /* if both solution and objective should be printed, we can use the solution to compute the objective */
      if ( objval != NULL )
      {
         if ( sdpisolver->penalty && ( ! sdpisolver->feasorig ))
         {
            /* in this case we cannot really trust the solution given by MOSEK, since changes in the value of r much less than epsilon can
             * cause huge changes in the objective, so using the objective value given by MOSEK is numerically more stable */
            MOSEK_CALL( MSK_getdualobj(sdpisolver->msktask, MSK_SOL_ITR, objval) );
         }
         else
         {
            /* since the objective value given by MOSEK sometimes differs slightly from the correct value for the given solution,
             * we get the solution from MOSEK and compute the correct objective value */
            *objval = 0.0;
            for (v = 0; v < sdpisolver->nactivevars; v++)
            {
               if ( REALABS(moseksol[v]) > sdpisolver->epsilon )
                  *objval += moseksol[v] * sdpisolver->objcoefs[v];
            }
         }

         /* as we didn't add the fixed (lb = ub) variables to MOSEK, we have to add their contributions to the objective as well */
         *objval += sdpisolver->fixedvarsobjcontr;
      }

      BMSfreeBufferMemoryArray(sdpisolver->bufmem, &moseksol);
   }
   else if ( objval != NULL )
   {
      SCIP_CALL( SCIPsdpiSolverGetObjval(sdpisolver, objval) );
   }

   return SCIP_OKAY;
}

/** gets the primal variables corresponding to the lower and upper variable-bounds in the dual problem
 *
 *  The last input should specify the length of the arrays. If this is less than the number of variables, the needed
 *  length will be returned and a debug message thrown.
 *
 *  @note If a variable is either fixed or unbounded in the dual problem, a zero will be returned for the non-existent primal variable.
 */
SCIP_RETCODE SCIPsdpiSolverGetPrimalBoundVars(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   SCIP_Real*            lbvars,             /**< pointer to store the values of the variables corresponding to lower bounds in the dual problems */
   SCIP_Real*            ubvars,             /**< pointer to store the values of the variables corresponding to upper bounds in the dual problems */
   int*                  arraylength         /**< input: length of lbvars and ubvars <br>
                                              *   output: number of elements inserted into lbvars/ubvars (or needed length if it wasn't sufficient) */
   )
{/*lint --e{818}*/
   SCIP_Real* primalvars;
   int nprimalvars;
   int i;

   assert( sdpisolver != NULL );
   CHECK_IF_SOLVED( sdpisolver );
   assert( arraylength != NULL );
   assert( lbvars != NULL );
   assert( ubvars != NULL );

   /* check if the arrays are long enough */
   if ( *arraylength < sdpisolver->nvars )
   {
      *arraylength = sdpisolver->nvars;
      SCIPdebugMessage("Insufficient length of array in SCIPsdpiSolverGetPrimalBoundVars (gave %d, needed %d)\n", *arraylength, sdpisolver->nvars);
      return SCIP_OKAY;
   }

   /* initialize the return-arrays with zero */
   for (i = 0; i < sdpisolver->nvars; i++)
   {
      lbvars[i] = 0.0;
      ubvars[i] = 0.0;
   }

   /* get number of primal variables in MOSEK */
   MOSEK_CALL( MSK_getnumvar(sdpisolver->msktask, &nprimalvars) );/*lint !e641*/

   BMSallocBufferMemoryArray(sdpisolver->bufmem, &primalvars, nprimalvars);

   MOSEK_CALL( MSK_getxx(sdpisolver->msktask, MSK_SOL_ITR, primalvars) );/*lint !e641*/

   /* iterate over all variable bounds and insert the corresponding primal variables in the right positions of the return-arrays */
   assert( sdpisolver->nvarbounds <= 2 * sdpisolver->nvars );

   for (i = 0; i < sdpisolver->nvarbounds; i++)
   {
      if ( sdpisolver->varboundpos[i] < 0 )
      {
         /* this is a lower bound */
         /* the last nvarbounds entries correspond to the varbounds */
         lbvars[sdpisolver->mosektoinputmapper[-1 * sdpisolver->varboundpos[i] -1]] = primalvars[nprimalvars - sdpisolver->nvarbounds + i]; /*lint !e679, !e834 */
      }
      else
      {
         /* this is an upper bound */

         assert( sdpisolver->varboundpos[i] > 0 );

         /* the last nvarbounds entries correspond to the varbounds */
         ubvars[sdpisolver->mosektoinputmapper[sdpisolver->varboundpos[i] - 1]] = primalvars[nprimalvars - sdpisolver->nvarbounds + i]; /*lint !e679, !e834 */
      }
   }

   BMSfreeBufferMemoryArray(sdpisolver->bufmem, &primalvars);

   return SCIP_OKAY;
}

/** gets the number of SDP iterations of the last solve call */
SCIP_RETCODE SCIPsdpiSolverGetIterations(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   int*                  iterations          /**< pointer to store the number of iterations of the last solve call */
   )
{/*lint --e{818}*/

   if ( sdpisolver->timelimitinitial )
      *iterations = 0;
   else
   {
      *iterations = sdpisolver->niterations;
   }

   return SCIP_OKAY;
}

/** gets the number of calls to the SDP-solver for the last solve call */
SCIP_RETCODE SCIPsdpiSolverGetSdpCalls(
   SCIP_SDPISOLVER*      sdpisolver,         /**< SDP-solver interface */
   int*                  calls               /**< pointer to store the number of calls to the SDP-solver for the last solve call */
   )
{/*lint --e{715,818,1784}*/
   assert( calls != NULL );

   *calls = sdpisolver->timelimitinitial ? 0 : sdpisolver->nsdpcalls;

   return SCIP_OKAY;
}

/** gets the settings used by the SDP solver for the last solve call */
SCIP_RETCODE SCIPsdpiSolverSettingsUsed(
   SCIP_SDPISOLVER*      sdpisolver,         /**< SDP interface solver structure */
   SCIP_SDPSOLVERSETTING* usedsetting        /**< the setting used by the SDP solver */
   )
{/*lint --e{818}*/
   assert( sdpisolver != NULL );
   assert( usedsetting != NULL );

   if ( ! SCIPsdpiSolverIsAcceptable(sdpisolver) )
      *usedsetting = SCIP_SDPSOLVERSETTING_UNSOLVED;
   else if ( sdpisolver->penalty )
      *usedsetting = SCIP_SDPSOLVERSETTING_PENALTY;
   else
      *usedsetting = SCIP_SDPSOLVERSETTING_FAST;

   return SCIP_OKAY;
}

/**@} */




/*
 * Numerical Methods
 */

/**@name Numerical Methods */
/**@{ */

/** returns value treated as infinity in the SDP-solver */
SCIP_Real SCIPsdpiSolverInfinity(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to an SDP interface solver structure */
   )
{/*lint --e{715,818}*/
   return 1.0e16;
}

/** checks if given value is treated as (plus or minus) infinity in the SDP-solver */
SCIP_Bool SCIPsdpiSolverIsInfinity(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   SCIP_Real             val                 /**< value to be checked for infinity */
   )
{/*lint --e{818}*/
   return ((val <= -SCIPsdpiSolverInfinity(sdpisolver)) || (val >= SCIPsdpiSolverInfinity(sdpisolver)));
}

/** gets floating point parameter of SDP-Solver */
SCIP_RETCODE SCIPsdpiSolverGetRealpar(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   SCIP_SDPPARAM         type,               /**< parameter number */
   SCIP_Real*            dval                /**< buffer to store the parameter value */
   )
{/*lint --e{818}*/
   assert( sdpisolver != NULL );
   assert( dval != NULL );

   switch( type )
   {
   case SCIP_SDPPAR_EPSILON:
      *dval = sdpisolver->epsilon;
      break;
   case SCIP_SDPPAR_GAPTOL:
         *dval = sdpisolver->gaptol;
         break;
   case SCIP_SDPPAR_FEASTOL:
      *dval = sdpisolver->feastol;
      break;
   case SCIP_SDPPAR_SDPSOLVERFEASTOL:
      *dval = sdpisolver->sdpsolverfeastol;
      break;
   case SCIP_SDPPAR_OBJLIMIT:
      *dval = sdpisolver->objlimit;
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }/*lint !e788*/

   return SCIP_OKAY;
}

/** sets floating point parameter of SDP-Solver */
SCIP_RETCODE SCIPsdpiSolverSetRealpar(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   SCIP_SDPPARAM         type,               /**< parameter number */
   SCIP_Real             dval                /**< parameter value */
   )
{
   assert( sdpisolver != NULL );

   switch( type )
   {
   case SCIP_SDPPAR_EPSILON:
      sdpisolver->epsilon = dval;
      SCIPdebugMessage("Setting sdpisolver epsilon to %f.\n", dval);
      break;
   case SCIP_SDPPAR_GAPTOL:
      sdpisolver->gaptol = dval;
      SCIPdebugMessage("Setting sdpisolver gaptol to %f.\n", dval);
      break;
   case SCIP_SDPPAR_FEASTOL:
      sdpisolver->feastol = dval;
      SCIPdebugMessage("Setting sdpisolver feastol to %f.\n", dval);
      break;
   case SCIP_SDPPAR_SDPSOLVERFEASTOL:
      sdpisolver->sdpsolverfeastol = dval;
      SCIPdebugMessage("Setting sdpisolver sdpsolverfeastol to %f.\n", dval);
      break;
   case SCIP_SDPPAR_OBJLIMIT:
      SCIPdebugMessage("Setting sdpisolver objlimit to %f.\n", dval);
      sdpisolver->objlimit = dval;
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }/*lint !e788*/

   return SCIP_OKAY;
}

/** gets integer parameter of SDP-Solver */
SCIP_RETCODE SCIPsdpiSolverGetIntpar(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   SCIP_SDPPARAM         type,               /**< parameter number */
   int*                  ival                /**< parameter value */
   )
{/*lint --e{818}*/
   assert( sdpisolver != NULL );

   switch( type )
   {
   case SCIP_SDPPAR_SDPINFO:
      *ival = (int) sdpisolver->sdpinfo;
      SCIPdebugMessage("Getting sdpisolver information output (%d).\n", *ival);
      break;
   case SCIP_SDPPAR_NTHREADS:
      *ival = sdpisolver->nthreads;
      SCIPdebugMessage("Getting sdpisolver number of threads: %d.\n", *ival);
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }/*lint !e788*/

   return SCIP_OKAY;
}

/** sets integer parameter of SDP-Solver */
SCIP_RETCODE SCIPsdpiSolverSetIntpar(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   SCIP_SDPPARAM         type,               /**< parameter number */
   int                   ival                /**< parameter value */
   )
{
   assert( sdpisolver != NULL );

   switch( type )
   {
   case SCIP_SDPPAR_NTHREADS:
      sdpisolver->nthreads = ival;
      SCIPdebugMessage("Setting sdpisolver number of threads to %d.\n", ival);
      break;
   case SCIP_SDPPAR_SDPINFO:
      assert( 0 <= ival && ival <= 1 );
      sdpisolver->sdpinfo = (SCIP_Bool) ival;
      SCIPdebugMessage("Setting sdpisolver information output (%d).\n", ival);
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }/*lint !e788*/

   return SCIP_OKAY;
}

/** compute and set lambdastar (only used for SDPA) */
SCIP_RETCODE SCIPsdpiSolverComputeLambdastar(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   SCIP_Real             maxguess            /**< maximum guess for lambda star of all SDP-constraints */
   )
{/*lint --e{715,818}*/
   SCIPdebugMessage("Lambdastar parameter not used by MOSEK"); /* this parameter is only used by SDPA */

   return SCIP_OKAY;
}

/** compute and set the penalty parameter */
SCIP_RETCODE SCIPsdpiSolverComputePenaltyparam(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   SCIP_Real             maxcoeff,           /**< maximum objective coefficient */
   SCIP_Real*            penaltyparam        /**< the computed penalty parameter */
   )
{/*lint --e{818,1784}*/
   SCIP_Real compval;

   assert( sdpisolver != NULL );
   assert( penaltyparam != NULL );

   compval = PENALTYPARAM_FACTOR * maxcoeff;

   if ( compval < MIN_PENALTYPARAM )
   {
      SCIPdebugMessage("Setting penaltyparameter to %f.\n", MIN_PENALTYPARAM);
      *penaltyparam = MIN_PENALTYPARAM;
   }
   else if ( compval > MAX_PENALTYPARAM )
   {
      SCIPdebugMessage("Setting penaltyparameter to %f.\n", MAX_PENALTYPARAM);
      *penaltyparam = MAX_PENALTYPARAM;
   }
   else
   {
      SCIPdebugMessage("Setting penaltyparameter to %f.\n", compval);
      *penaltyparam = compval;
   }
   return SCIP_OKAY;
}

/** compute and set the maximum penalty parameter */
SCIP_RETCODE SCIPsdpiSolverComputeMaxPenaltyparam(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   SCIP_Real             penaltyparam,       /**< the initial penalty parameter */
   SCIP_Real*            maxpenaltyparam     /**< the computed maximum penalty parameter */
   )
{/*lint --e{818,1784}*/
   SCIP_Real compval;

   assert( sdpisolver != NULL );
   assert( maxpenaltyparam != NULL );

   compval = penaltyparam * MAXPENALTYPARAM_FACTOR;

   if ( compval < MAX_MAXPENALTYPARAM )
   {
      *maxpenaltyparam = compval;
      SCIPdebugMessage("Setting maximum penaltyparameter to %f.\n", compval);
   }
   else
   {
      *maxpenaltyparam = MAX_MAXPENALTYPARAM;
      SCIPdebugMessage("Setting penaltyparameter to %f.\n", MAX_MAXPENALTYPARAM);
   }

   return SCIP_OKAY;
}

/**@} */




/*
 * File Interface Methods
 */

/**@name File Interface Methods */
/**@{ */

/** reads SDP from a file */
SCIP_RETCODE SCIPsdpiSolverReadSDP(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   const char*           fname               /**< file name */
   )
{/*lint --e{715,818}*/
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_LPERROR;
}

/** writes SDP to a file */
SCIP_RETCODE SCIPsdpiSolverWriteSDP(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   const char*           fname               /**< file name */
   )
{/*lint --e{818}*/
   assert( sdpisolver != NULL );
   assert( fname != NULL );

   MOSEK_CALL( MSK_writedata(sdpisolver->msktask, fname) );/*lint !e641*/

   return SCIP_OKAY;
}

/**@} */
