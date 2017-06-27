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

/**@file   sdpisolver_dsdp.c
 * @brief  interface for DSDP
 * @author Tristan Gally
 */

#include <assert.h>
#include <sys/time.h>

#include "sdpi/sdpisolver.h"

/* turn off warning for DSDSP */
#pragma GCC diagnostic ignored "-Wstrict-prototypes"
#include "dsdp5.h"                           /* for DSDPUsePenalty, etc */
#pragma GCC diagnostic warning "-Wstrict-prototypes"


#include "blockmemshell/memory.h"            /* for memory allocation */
#include "scip/def.h"                        /* for SCIP_Real, _Bool, ... */
#include "scip/pub_misc.h"                   /* for sorting */
#include "sdpi/sdpsolchecker.h"              /* to check solution with regards to feasibility tolerance */

/* turn off lint warnings for whole file: */
/*lint --e{788,818}*/

#define PENALTYBOUNDTOL             1E-3     /**< if the relative gap between Tr(X) and penaltyparam for a primal solution of the penaltyformulation
                                              *   is bigger than this value, it will be reported to the sdpi */

#define MIN_PENALTYPARAM            1e5      /**< if the penalty parameter is to be computed, this is the minimum value it will take */
#define MAX_PENALTYPARAM            1e12     /**< if the penalty parameter is to be computed, this is the maximum value it will take */
#define PENALTYPARAM_FACTOR         1e4      /**< if the penalty parameter is to be computed, the maximal objective coefficient will be multiplied by this */
#define MAX_MAXPENALTYPARAM         1e15     /**< if the maximum penaltyparameter is to be computed, this is the maximum value it will take */
#define MAXPENALTYPARAM_FACTOR      1e6      /**< if the maximum penaltyparameter is to be computed, it will be set to penaltyparam * this */
#define INFEASFEASTOLCHANGE         0.1      /**< change feastol by this factor if the solution was found to be infeasible with regards to feastol */
#define INFEASMINFEASTOL            1E-9     /**< minimum value for feasibility tolerance when encountering problems with regards to tolerance */

/** Calls a DSDP-Function and transforms the return-code to a SCIP_LPERROR if needed. */
#define DSDP_CALL(x)  do                                                                                     \
                      {                                                                                      \
                         int _dsdperrorcode_;                                                                \
                         if ( (_dsdperrorcode_ = (x)) != 0 )                                                 \
                         {                                                                                   \
                            SCIPerrorMessage("DSDP-Error <%d> in function call.\n", _dsdperrorcode_);        \
                            return SCIP_LPERROR;                                                             \
                         }                                                                                   \
                      }                                                                                      \
                      while( FALSE )

/** Same as DSDP_CALL, but used for functions returning a boolean. */
#define DSDP_CALL_BOOL(x)  do                                                                                \
                      {                                                                                      \
                         int _dsdperrorcode_;                                                                \
                         if ( (_dsdperrorcode_ = (x)) != 0 )                                                 \
                         {                                                                                   \
                            SCIPerrorMessage("DSDP-Error <%d> in function call.\n", _dsdperrorcode_);        \
                            return FALSE;                                                                    \
                         }                                                                                   \
                      }                                                                                      \
                      while( FALSE )

/** Same as DSDP_CALL, but this will be used for initialization methods with memory allocation and return a SCIP_NOMEMORY if an error is produced. */
#define DSDP_CALLM(x) do                                                                                     \
                      {                                                                                      \
                         int _dsdperrorcode_;                                                                \
                         if ( (_dsdperrorcode_ = (x)) != 0 )                                                 \
                         {                                                                                   \
                            SCIPerrorMessage("DSDP-Error <%d> in function call.\n", _dsdperrorcode_);        \
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

/** This is the same as CHECK_IF_SOLVED, but will be called for methods returning a bool instead of a SCIP_RETURNCODE. */
#define CHECK_IF_SOLVED_BOOL(sdpisolver)  do                                                                      \
                      {                                                                                      \
                         if (!(sdpisolver->solved))                                                          \
                         {                                                                                   \
                            SCIPerrorMessage("Tried to access solution information for SDP %d ahead of solving!\n", sdpisolver->sdpcounter);  \
                            return FALSE;                                                                    \
                         }                                                                                   \
                      }                                                                                      \
                      while( FALSE )


/** data used for SDP interface */
struct SCIP_SDPiSolver
{
   SCIP_MESSAGEHDLR*     messagehdlr;        /**< messagehandler for printing messages, or NULL */
   BMS_BLKMEM*           blkmem;             /**< block memory */
   BMS_BUFMEM*           bufmem;             /**< buffer memory */
   DSDP                  dsdp;               /**< solver-object */
   SDPCone               sdpcone;            /**< sdpcone-object of DSDP for handling SDP-constraints */
   LPCone                lpcone;             /**< lpcone-object of DSDP for handling LP-constraints */
   BCone                 bcone;              /**< bcone-object of DSDP to add variable bounds to */
   int                   nvars;              /**< number of input variables */
   int                   nactivevars;        /**< number of variables present in DSDP (nvars minus the number of variables with lb = ub) */
   int*                  inputtodsdpmapper;  /**< entry i gives the index of input variable i in dsdp (starting from 1) or
                                               *  -j (j=1, 2, ..., nvars - nactivevars) if the variable is fixed, the value and objective value of
                                               *  this fixed variable can be found in entry j-1 of fixedval/obj */
   int*                  dsdptoinputmapper;  /**< entry i gives the original index of the (i+1)-th variable in dsdp (indices go from 0 to nactivevars-1) */
   SCIP_Real*            fixedvarsval;       /**< entry i gives the lower and upper bound of the i-th fixed variable */
   SCIP_Real             fixedvarsobjcontr;  /**< total contribution to the objective of all fixed variables, computed as sum obj * val */
   SCIP_Real*            objcoefs;           /**< objective coefficients of all active variables */
   SCIP_Bool             solved;             /**< Was the SDP solved since the problem was last changed? */
   int                   sdpcounter;         /**< used for debug messages */
   SCIP_Real             epsilon;            /**< tolerance for absolute checks */
   SCIP_Real             gaptol;             /**< this is used for checking if primal and dual objective are equal */
   SCIP_Real             feastol;            /**< feasibility tolerance that should be achieved */
   SCIP_Real             sdpsolverfeastol;   /**< feasibility tolerance given to the SDP-solver */
   SCIP_Real             penaltyparam;       /**< the penalty parameter Gamma used for the penalty formulation if the SDP-solver didn't converge */
   SCIP_Real             objlimit;           /**< objective limit for SDP-solver */
   SCIP_Bool             sdpinfo;            /**< Should the SDP-solver output information to the screen? */
   SCIP_Bool             penalty;            /**< Did the last solve use a penalty formulation? */
   SCIP_Bool             penaltyworbound;    /**< Was a penalty formulation solved without bounding r? */
   SCIP_Bool             feasorig;           /**< was the last problem solved with a penalty formulation and with original objective coefficents
                                               *  and the solution was feasible for the original problem? */
   SCIP_SDPSOLVERSETTING usedsetting;        /**< setting used to solve the last SDP */
   SCIP_Bool             timelimit;          /**< was the solver stopped because of the time limit? */
   SCIP_Bool             timelimitinitial;   /**< was the problem not even given to the solver because of the time limit? */
   int                   niterations;        /**< number of SDP-iterations since the last solve call */
   int                   nsdpcalls;          /**< number of SDP-calls since the last solve call */
};

typedef struct Timings
{
   struct timeval        starttime;          /**< time when solving started */
   SCIP_Real             timelimit;          /**< timelimit in seconds */
   SCIP_Bool             stopped;            /**< was the solver stopped because of the time limit? */
} Timings;


/*
 * Local Functions
 */

/** for given row and column (i,j) computes the position in the lower triangular part
 *  numbered from 0 to n(n+1)/2 - 1, this needs to be called for i >= j
 */
static
int compLowerTriangPos(
   int                   i,                  /**< row index */
   int                   j                   /**< column index */
   )
{
   assert( j >= 0 );
   assert( i >= j );

   return i*(i+1)/2 + j;
}

#ifndef NDEBUG
/** test if a lower bound lb is not smaller than an upper bound ub, meaning that lb > ub - epsilon */
static
SCIP_Bool isFixed(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   SCIP_Real             lb,                 /**< lower bound */
   SCIP_Real             ub                  /**< upper bound */
   )
{
   assert( lb < ub + sdpisolver->feastol );

   return (ub-lb <= sdpisolver->epsilon);
}
#else
#define isFixed(sdpisolver,lb,ub) (ub-lb <= sdpisolver->epsilon)
#endif

/** sort the given row, col and val arrays first by non-decreasing col-indices, than for those with identical col-indices by non-increasing row-indices */
static
void sortColRow(
   int*                  row,                /**< row indices */
   int*                  col,                /**< column indices */
   SCIP_Real*            val,                /**< values */
   int                   length              /**< length of the given arrays */
   )
{
   int firstentry;
   int nextentry = 0;

   /* first sort by col indices */
   SCIPsortIntIntReal(col, row, val, length);

   /* for those with identical col-indices now sort by non-decreasing row-index, first find all entries with the same col-index */
   while (nextentry < length)
   {
      firstentry = nextentry; /* the next col starts where the last one ended */

      while (nextentry < length && col[nextentry] == col[firstentry]) /* as long as the row still matches, increase nextentry */
         ++nextentry;

      /* now sort all entries between firstentry and nextentry-1 by their row-indices */
      SCIPsortIntReal(row + firstentry, val + firstentry, nextentry - firstentry);
   }
}

/** check the time limit after each iteration in DSDP */
static
int checkTimeLimitDSDP(
   DSDP                  dsdp,               /**< DSDP-pointer */
   void*                 ctx                 /**< pointer to data of iteration monitor */
   )
{
   Timings* timings;
   struct timeval currenttime;
   SCIP_Real startseconds;
   SCIP_Real currentseconds;
   SCIP_Real elapsedtime;

   assert( dsdp != NULL );
   assert( ctx != NULL );

   timings = (Timings*) ctx;

   startseconds = (SCIP_Real) (timings->starttime).tv_sec + (SCIP_Real) (timings->starttime).tv_usec / 1e6;

   TIMEOFDAY_CALL( gettimeofday(&currenttime, NULL) );/*lint !e438, !e550, !e641 */
   currentseconds = (SCIP_Real) currenttime.tv_sec + (SCIP_Real) currenttime.tv_usec / 1e6;

   elapsedtime = currentseconds - startseconds;

   if ( elapsedtime > timings->timelimit )
   {
      DSDP_CALL( DSDPSetConvergenceFlag(dsdp, DSDP_USER_TERMINATION) );/*lint !e641 */
      timings->stopped = TRUE;
      SCIPdebugMessage("Time limit reached! Stopping DSDP.\n");
   }

   return 0;
}


/*
 * Miscellaneous Methods
 */

/**@name Miscellaneous Methods */
/**@{ */


/** gets name and version (if available) of SDP-solver*/
const char* SCIPsdpiSolverGetSolverName(
   void
   )
{
   return "DSDP"; /* getting the version is not supported in DSDP */
}

/** gets description of SDP-solver (developer, webpage, ...) */
const char* SCIPsdpiSolverGetSolverDesc(
   void
   )
{
   return "Dual-Scaling Interior Point SDP-Solver by S. Benson, Y. Ye, and X. Zhang (http://www.mcs.anl.gov/hs/software/DSDP/)";
}

/** gets pointer to SDP-solver - use only with great care
 *
 *  The behavior of this function depends on the solver and its use is
 *  therefore only recommended if you really know what you are
 *  doing. In general, it returns a pointer to the SDP-solver object.
 */
void* SCIPsdpiSolverGetSolverPointer(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to an SDP-solver interface */
   )
{
   assert( sdpisolver != NULL );
   return (void*) sdpisolver->dsdp;
}

/** gets default feasibility tolerance for SDP-solver in SCIP-SDP */
SCIP_Real SCIPsdpiSolverGetDefaultSdpiSolverFeastol(
   void
   )
{
   return 1E-6;
}

/** gets default number of increases of penalty parameter for SDP-solver in SCIP-SDP */
int SCIPsdpiSolverGetDefaultSdpiSolverNpenaltyIncreases(
   void
   )
{
   return 10;
}

/**@} */


/*
 * SDPI Creation and Destruction Methods
 */

/**@name SDPI Creation and Destruction Methods */
/**@{ */

/** creates an SDP solver interface */
SCIP_RETCODE SCIPsdpiSolverCreate(
   SCIP_SDPISOLVER**     sdpisolver,         /**< pointer to an SDP-solver interface */
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

   /* the following four variables will be properly initialized only immediatly prior to solving because DSDP and the
    * SDPCone need information about the number of variables and sdpblocks during creation */
   (*sdpisolver)->dsdp = NULL;
   (*sdpisolver)->sdpcone = NULL;
   (*sdpisolver)->lpcone = NULL;
   (*sdpisolver)->bcone = NULL;

   (*sdpisolver)->nvars = 0;
   (*sdpisolver)->nactivevars = 0;
   (*sdpisolver)->inputtodsdpmapper = NULL;
   (*sdpisolver)->dsdptoinputmapper = NULL;
   (*sdpisolver)->fixedvarsval = NULL;
   (*sdpisolver)->fixedvarsobjcontr = 0.0;
   (*sdpisolver)->objcoefs = NULL;
   (*sdpisolver)->solved = FALSE;
   (*sdpisolver)->timelimit = FALSE;
   (*sdpisolver)->timelimitinitial = FALSE;
   (*sdpisolver)->penalty = FALSE;
   (*sdpisolver)->penaltyworbound = FALSE;
   (*sdpisolver)->feasorig = FALSE;
   (*sdpisolver)->sdpcounter = 0;
   (*sdpisolver)->niterations = 0;
   (*sdpisolver)->nsdpcalls = 0;

   (*sdpisolver)->epsilon = 1e-9;
   (*sdpisolver)->gaptol = 1e-4;
   (*sdpisolver)->feastol = 1e-6;
   (*sdpisolver)->sdpsolverfeastol = 1e-6;
   (*sdpisolver)->penaltyparam = 1e5;
   (*sdpisolver)->objlimit = SCIPsdpiSolverInfinity(*sdpisolver);
   (*sdpisolver)->sdpinfo = FALSE;
   (*sdpisolver)->usedsetting = SCIP_SDPSOLVERSETTING_UNSOLVED;

   return SCIP_OKAY;
}

/** deletes an SDP solver interface */
SCIP_RETCODE SCIPsdpiSolverFree(
   SCIP_SDPISOLVER**     sdpisolver          /**< pointer to an SDP-solver interface */
   )
{
   assert( sdpisolver != NULL );
   assert( *sdpisolver != NULL );

   SCIPdebugMessage("Freeing SDPISolver\n");

   if ( (*sdpisolver)->dsdp != NULL )
   {
      DSDP_CALL( DSDPDestroy((*sdpisolver)->dsdp) );
   }

   if ( (*sdpisolver)->nvars > 0 )
      BMSfreeBlockMemoryArray((*sdpisolver)->blkmem, &(*sdpisolver)->inputtodsdpmapper, (*sdpisolver)->nvars);/*lint !e737 */

   if ( (*sdpisolver)->nactivevars > 0 )
   {
      BMSfreeBlockMemoryArray((*sdpisolver)->blkmem, &(*sdpisolver)->dsdptoinputmapper, (*sdpisolver)->nactivevars);/*lint !e737 */
      BMSfreeBlockMemoryArray((*sdpisolver)->blkmem, &(*sdpisolver)->objcoefs, (*sdpisolver)->nactivevars); /*lint !e776*/
   }

   if ( (*sdpisolver)->nvars >= (*sdpisolver)->nactivevars )
      BMSfreeBlockMemoryArrayNull((*sdpisolver)->blkmem, &(*sdpisolver)->fixedvarsval, (*sdpisolver)->nvars - (*sdpisolver)->nactivevars); /*lint !e776*/

   BMSfreeBlockMemory((*sdpisolver)->blkmem, sdpisolver);

   return SCIP_OKAY;
}

/** increases the SDP-Counter */
SCIP_RETCODE SCIPsdpiSolverIncreaseCounter(
   SCIP_SDPISOLVER*      sdpisolver          /**< SDP-solver interface */
   )
{
   assert( sdpisolver != NULL );

   sdpisolver->sdpcounter++;

   return SCIP_OKAY;
}

/** reset the SDP-Counter to zero */
SCIP_RETCODE SCIPsdpiSolverResetCounter(
   SCIP_SDPISOLVER*      sdpisolver          /**< SDP-solver interface */
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
   SCIP_SDPISOLVER*      sdpisolver,         /**< SDP-solver interface */
   int                   nvars,              /**< number of variables */
   SCIP_Real*            obj,                /**< objective coefficients of variables */
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
   int                   sdpnnonz,           /**< number of nonzero elements in the SDP-constraint-matrix */
   int**                 sdpnblockvarnonz,   /**< entry [i][j] gives the number of nonzeros for block i and variable j, this is exactly
                                              *   the number of entries of sdp row/col/val [i][j] */
   int**                 sdpvar,             /**< sdpvar[i][j] gives the sdp-index of the j-th variable (according to the sorting for row/col/val)
                                              *   in the i-th block */
   int***                sdprow,             /**< pointer to the row-indices for each block and variable */
   int***                sdpcol,             /**< pointer to the column-indices for each block and variable */
   SCIP_Real***          sdpval,             /**< values of SDP-constraint-matrix entries (may be NULL if sdpnnonz = 0) */
   int**                 indchanges,         /**< changes needed to be done to the indices, if indchanges[block][nonz]=-1, then
                                              *   the index can be removed, otherwise it gives the number of indices removed before this */
   int*                  nremovedinds,       /**< the number of rows/cols to be fixed for each block */
   int*                  blockindchanges,    /**< block indizes will be modified by these, see indchanges */
   int                   nremovedblocks,     /**< number of empty blocks that should be removed */
   int                   nlpcons,            /**< number of active (at least two nonzeros) LP-constraints */
   int                   noldlpcons,         /**< number of LP-constraints including those with less than two active nonzeros */
   SCIP_Real*            lplhs,              /**< left-hand sides of active LP-rows after fixings (may be NULL if nlpcons = 0) */
   SCIP_Real*            lprhs,              /**< right-hand sides of active LP-rows after fixings (may be NULL if nlpcons = 0) */
   int*                  rownactivevars,     /**< number of active variables for each LP-constraint */
   int                   lpnnonz,            /**< number of nonzero elements in the LP-constraint-matrix */
   int*                  lprow,              /**< row-index for each entry in lpval-array, might get sorted (may be NULL if lpnnonz = 0) */
   int*                  lpcol,              /**< column-index for each entry in lpval-array, might get sorted (may be NULL if lpnnonz = 0) */
   SCIP_Real*            lpval,              /**< values of LP-constraint-matrix entries, might get sorted (may be NULL if lpnnonz = 0) */
   SCIP_Real*            start,              /**< NULL or a starting point for the solver, this should have length nvars */
   SCIP_SDPSOLVERSETTING startsettings,      /**< settings used to start with in SDPA, currently not used for DSDP ans MOSEK, set this to
                                               *  SCIP_SDPSOLVERSETTING_UNSOLVED to ignore it and start from scratch */
   SCIP_Real             timelimit           /**< after this many seconds solving will be aborted (currently only implemented for DSDP and MOSEK) */
   )
{
   return SCIPsdpiSolverLoadAndSolveWithPenalty(sdpisolver, 0.0, TRUE, TRUE, nvars, obj, lb, ub, nsdpblocks, sdpblocksizes, sdpnblockvars,
           sdpconstnnonz, sdpconstnblocknonz, sdpconstrow, sdpconstcol, sdpconstval, sdpnnonz, sdpnblockvarnonz, sdpvar, sdprow, sdpcol, sdpval,
           indchanges, nremovedinds, blockindchanges, nremovedblocks, nlpcons, noldlpcons, lplhs, lprhs, rownactivevars, lpnnonz, lprow, lpcol,
           lpval, start, startsettings, timelimit, NULL, NULL);
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
   SCIP_SDPISOLVER*      sdpisolver,         /**< SDP-solver interface */
   SCIP_Real             penaltyparam,       /**< the Gamma above, needs to be >= 0 */
   SCIP_Bool             withobj,            /**< if this is false the objective is set to 0 */
   SCIP_Bool             rbound,             /**< should r be non-negative ? */
   int                   nvars,              /**< number of variables */
   SCIP_Real*            obj,                /**< objective coefficients of variables */
   SCIP_Real*            lb,                 /**< lower bounds of variables */
   SCIP_Real*            ub,                 /**< upper bounds of variables */
   int                   nsdpblocks,         /**< number of SDP-blocks */
   int*                  sdpblocksizes,      /**< sizes of the SDP-blocks (may be NULL if nsdpblocks = sdpconstnnonz = sdpnnonz = 0) */
   int*                  sdpnblockvars,      /**< number of variables that exist in each block */
   int                   sdpconstnnonz,      /**< number of nonzero elements in the constant matrices of the SDP-blocks AFTER FIXINGS */
   int*                  sdpconstnblocknonz, /**< number of nonzeros for each variable in the constant part, also the i-th entry gives the
                                              *   number of entries  of sdpconst row/col/val [i] AFTER FIXINGS */
   int**                 sdpconstrow,        /**< pointers to row-indices for each block AFTER FIXINGS */
   int**                 sdpconstcol,        /**< pointers to column-indices for each block AFTER FIXINGS */
   SCIP_Real**           sdpconstval,        /**< pointers to the values of the nonzeros for each block AFTER FIXINGS */
   int                   sdpnnonz,           /**< number of nonzero elements in the SDP-constraint-matrix */
   int**                 sdpnblockvarnonz,   /**< entry [i][j] gives the number of nonzeros for block i and variable j, this is exactly
                                              *   the number of entries of sdp row/col/val [i][j] */
   int**                 sdpvar,             /**< sdpvar[i][j] gives the sdp-index of the j-th variable (according to the sorting for row/col/val)
                                              *   in the i-th block */
   int***                sdprow,             /**< pointer to the row-indices for each block and variable */
   int***                sdpcol,             /**< pointer to the column-indices for each block and variable */
   SCIP_Real***          sdpval,             /**< values of SDP-constraint-matrix entries (may be NULL if sdpnnonz = 0) */
   int**                 indchanges,         /**< changes needed to be done to the indices, if indchanges[block][nonz]=-1, then
                                              *   the index can be removed, otherwise it gives the number of indices removed before this */
   int*                  nremovedinds,       /**< the number of rows/cols to be fixed for each block */
   int*                  blockindchanges,    /**< block indizes will be modified by these, see indchanges */
   int                   nremovedblocks,     /**< number of empty blocks that should be removed */
   int                   nlpcons,            /**< number of active (at least two nonzeros) LP-constraints */
   int                   noldlpcons,         /**< number of LP-constraints including those with less than two active nonzeros */
   SCIP_Real*            lplhs,              /**< left-hand sides of active LP-rows after fixings (may be NULL if nlpcons = 0) */
   SCIP_Real*            lprhs,              /**< right-hand sides of active LP-rows after fixings (may be NULL if nlpcons = 0) */
   int*                  rownactivevars,     /**< number of active variables for each LP-constraint */
   int                   lpnnonz,            /**< number of nonzero elements in the LP-constraint-matrix */
   int*                  lprow,              /**< row-index for each entry in lpval-array, might get sorted (may be NULL if lpnnonz = 0) */
   int*                  lpcol,              /**< column-index for each entry in lpval-array, might get sorted (may be NULL if lpnnonz = 0) */
   SCIP_Real*            lpval,              /**< values of LP-constraint-matrix entries, might get sorted (may be NULL if lpnnonz = 0) */
   SCIP_Real*            start,              /**< NULL or a starting point for the solver, this should have length nvars */
   SCIP_SDPSOLVERSETTING startsettings,      /**< settings used to start with in SDPA, currently not used for DSDP and MOSEK, set this to
                                               *  SCIP_SDPSOLVERSETTING_UNSOLVED to ignore it and start from scratch */
   SCIP_Real             timelimit,          /**< after this many seconds solving will be aborted (currently only implemented for DSDP and MOSEK) */
   SCIP_Bool*            feasorig,           /**< pointer to store if the solution to the penalty-formulation is feasible for the original problem
                                               *  (may be NULL if penaltyparam = 0) */
   SCIP_Bool*            penaltybound        /**< pointer to store if the primal solution reached the bound Tr(X) <= penaltyparam in the primal problem,
                                               *  this is also an indication of the penalty parameter being to small (may be NULL if not needed) */
   )
{/*lint --e{413}*/
   int* dsdpconstind = NULL;  /* indices for constant SDP-constraint-matrices, needs to be stored for DSDP during solving and be freed only afterwards */
   SCIP_Real* dsdpconstval = NULL; /* non-zero values for constant SDP-constraint-matrices, needs to be stored for DSDP during solving and be freed only afterwards */
   int* dsdpind = NULL;       /* indices for SDP-constraint-matrices, needs to be stored for DSDP during solving and be freed only afterwards */
   SCIP_Real* dsdpval = NULL;    /* non-zero values for SDP-constraint-matrices, needs to be stored for DSDP during solving and be freed only afterwards */
   int* dsdplpbegcol = NULL;  /* starting-indices for all columns in LP, needs to be stored for DSDP during solving and be freed only afterwards */
   int* dsdplprow = NULL;     /* row indices in LP, needs to be stored for DSDP during solving and be freed only afterwards */
   SCIP_Real* dsdplpval = NULL;  /* nonzeroes in LP, needs to be stored for DSDP during solving and be freed only afterwards */
   int i;
   int j;
   int ind;
   int block;
   int startind;
   int nfixedvars;
   int dsdpnlpnonz = 0;
   int nrnonz = 0;
   SCIP_Real feastol;
   Timings timings;

#ifdef SCIP_DEBUG
   DSDPTerminationReason reason; /* this will later be used to check if DSDP converged */
#endif

   assert( sdpisolver != NULL );
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

   sdpisolver->penalty = penaltyparam > sdpisolver->epsilon;

   if ( timelimit <= 0.0 )
   {
      sdpisolver->timelimit = TRUE;
      sdpisolver->timelimitinitial = TRUE;
      sdpisolver->solved = FALSE;
      return SCIP_OKAY;
   }
   else
      sdpisolver->timelimitinitial = FALSE;

   sdpisolver->feasorig = FALSE;

   /* start the timing */
   TIMEOFDAY_CALL( gettimeofday(&(timings.starttime), NULL) );/*lint !e438, !e550, !e641 */
   timings.timelimit = timelimit;
   timings.stopped = FALSE;

   /* only increase the counter if we don't use the penalty formulation to stay in line with the numbers in the general interface (where this is still the
    * same SDP), also remember settings for statistics */
   if ( penaltyparam < sdpisolver->epsilon )
   {
      SCIPdebugMessage("Inserting Data into DSDP for SDP (%d) \n", ++sdpisolver->sdpcounter);
      sdpisolver->usedsetting = SCIP_SDPSOLVERSETTING_FAST;
   }
   else
   {
      SCIPdebugMessage("Inserting Data again into DSDP for SDP (%d) \n", sdpisolver->sdpcounter);
      sdpisolver->usedsetting = SCIP_SDPSOLVERSETTING_PENALTY;
   }

   /* allocate memory for inputtomosekmapper, mosektoinputmapper and the fixed and active variable information, for the latter this will
    * later be shrinked if the needed size is known */
   BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &(sdpisolver->inputtodsdpmapper), sdpisolver->nvars, nvars) );
   BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &(sdpisolver->dsdptoinputmapper), sdpisolver->nactivevars, nvars) );
   BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &(sdpisolver->fixedvarsval), sdpisolver->nvars - sdpisolver->nactivevars, nvars) ); /*lint !e776*/
   BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &(sdpisolver->objcoefs), sdpisolver->nactivevars, nvars) ); /*lint !e776*/

   sdpisolver->nvars = nvars;
   sdpisolver->nactivevars = 0;
   nfixedvars = 0;
   sdpisolver->niterations = 0;
   sdpisolver->nsdpcalls = 0;

   /* find the fixed variables */
   sdpisolver->fixedvarsobjcontr = 0.0;
   for (i = 0; i < nvars; i++)
   {
      if ( isFixed(sdpisolver, lb[i], ub[i]) )
      {
         nfixedvars++;
         sdpisolver->inputtodsdpmapper[i] = -nfixedvars;
         sdpisolver->fixedvarsobjcontr += obj[i] * lb[i]; /* this is the value this variable contributes to the objective */
         sdpisolver->fixedvarsval[nfixedvars - 1] = lb[i]; /* if lb=ub, than this is the value the variable will have in every solution */
         SCIPdebugMessage("Fixing variable %d locally to %f for SDP %d in DSDP\n", i, lb[i], sdpisolver->sdpcounter);
      }
      else
      {
         sdpisolver->dsdptoinputmapper[sdpisolver->nactivevars] = i;
         sdpisolver->objcoefs[sdpisolver->nactivevars] = obj[i];
         sdpisolver->nactivevars++;
         sdpisolver->inputtodsdpmapper[i] = sdpisolver->nactivevars; /* dsdp starts counting at 1, so we do this after increasing nactivevars */
         SCIPdebugMessage("Variable %d becomes variable %d for SDP %d in DSDP\n", i, sdpisolver->inputtodsdpmapper[i], sdpisolver->sdpcounter);
      }
   }
   assert( sdpisolver->nactivevars + nfixedvars == sdpisolver->nvars );
   if ( penaltyparam > sdpisolver->epsilon && (! rbound) )
   {
      SCIPdebugMessage("Variable %d is the slack variable for the explicit penalty formulation\n", sdpisolver->nactivevars + 1);
   }

   /* if we want to solve without objective, we reset fixedvarsobjcontr */
   if ( ! withobj )
      sdpisolver->fixedvarsobjcontr = 0.0;

   /* shrink the fixedvars, objcoefs and mosektoinputmapper arrays to the right size */
   BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &(sdpisolver->objcoefs), nvars, sdpisolver->nactivevars) );
   BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &(sdpisolver->fixedvarsval), nvars, nfixedvars) );
   BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &(sdpisolver->dsdptoinputmapper), nvars, sdpisolver->nactivevars) );

   /* insert data */

   if ( sdpisolver->dsdp != NULL )
   {
      DSDP_CALL( DSDPDestroy(sdpisolver->dsdp) ); /* if there already exists a DSDP-instance, destroy the old one */
   }

   /* in case we don't want to bound r, we can't use the penalty formulation in DSDP and have to give r explicitly */
   if ( penaltyparam > sdpisolver->epsilon && (! rbound) )
   {
      DSDP_CALLM( DSDPCreate(sdpisolver->nactivevars + 1, &(sdpisolver->dsdp)) );
      sdpisolver->penaltyworbound = TRUE;
   }
   else
   {
      DSDP_CALLM( DSDPCreate(sdpisolver->nactivevars, &(sdpisolver->dsdp)) );
      sdpisolver->penaltyworbound = FALSE;
   }
   DSDP_CALLM( DSDPCreateSDPCone(sdpisolver->dsdp, nsdpblocks - nremovedblocks, &(sdpisolver->sdpcone)) );
   DSDP_CALLM( DSDPCreateLPCone(sdpisolver->dsdp, &(sdpisolver->lpcone)) );
   DSDP_CALLM( DSDPCreateBCone(sdpisolver->dsdp, &(sdpisolver->bcone)) );

#ifdef SCIP_MORE_DEBUG
   SCIPmessagePrintInfo(sdpisolver->messagehdlr, "setting objective values for SDP %d:\n", sdpisolver->sdpcounter);
#endif

   for (i = 0; i < sdpisolver->nactivevars; i++)
   {
      if ( withobj )
      {
         /* insert objective value, DSDP counts from 1 to n instead of 0 to n-1, *(-1) because DSDP maximizes instead of minimizing */
         DSDP_CALL( DSDPSetDualObjective(sdpisolver->dsdp, i+1, -1.0 * obj[sdpisolver->dsdptoinputmapper[i]]) );
#ifdef SCIP_MORE_DEBUG
         SCIPmessagePrintInfo(sdpisolver->messagehdlr, "var %d (was var %d): %f, ", i+1, sdpisolver->dsdptoinputmapper[i], obj[sdpisolver->dsdptoinputmapper[i]]);
#endif
      }
      else
      {
         DSDP_CALL( DSDPSetDualObjective(sdpisolver->dsdp, i+1, 0.0) );
      }

      if ( ! SCIPsdpiSolverIsInfinity(sdpisolver, lb[sdpisolver->dsdptoinputmapper[i]]) )
      {
         /* insert lower bound, DSDP counts from 1 to n instead of 0 to n-1 */
         DSDP_CALL( BConeSetLowerBound(sdpisolver->bcone, i+1, lb[sdpisolver->dsdptoinputmapper[i]]) );
      }

      if ( ! SCIPsdpiSolverIsInfinity(sdpisolver, ub[sdpisolver->dsdptoinputmapper[i]]) )
      {
         /* insert upper bound, DSDP counts from 1 to n instead of 0 to n-1 */
         DSDP_CALL(BConeSetUpperBound(sdpisolver->bcone, i+1, ub[sdpisolver->dsdptoinputmapper[i]]));
      }
   }

   /* insert the objective value for r if solving without rbound, it is variable nactivevars + 1 and the objective is multiplied by -1 as we maximize */
   if ( penaltyparam > sdpisolver->epsilon && (! rbound) )
   {
      DSDP_CALL( DSDPSetDualObjective(sdpisolver->dsdp, sdpisolver->nactivevars + 1, -1.0 * penaltyparam) );
#ifdef SCIP_MORE_DEBUG
      SCIPmessagePrintInfo(sdpisolver->messagehdlr, "slack variable r: %f, ", penaltyparam);
#endif
   }

#ifdef SCIP_MORE_DEBUG
   SCIPmessagePrintInfo(sdpisolver->messagehdlr, "\n");
   SCIPdebugMessage("ATTENTION: BConeView shows the WRONG sign for the lower bound!\n");
   BConeView(sdpisolver->bcone);
#endif

   /* set blocksizes */
   for (block = 0; block < nsdpblocks; ++block)
   {
      /* only insert blocksizes for the blocks we didn't remove */
      if ( blockindchanges[block] > -1 )
      {
         /* (blocks are counted from 0 to m-1) */
         DSDP_CALL( SDPConeSetBlockSize(sdpisolver->sdpcone, block- blockindchanges[block], sdpblocksizes[block] - nremovedinds[block]) );
      }
   }

   /* start inserting the non-constant SDP-Constraint-Matrices */
   if ( sdpnnonz > 0 )
   {
      int v;
      int k;
      int blockvar;

      nrnonz = 0;

      /* allocate memory */
      /* This needs to be one long array, because DSDP uses it for solving so all nonzeros have to be in it and it may not be freed before the
       * problem is solved. The distinct blocks/variables (for the i,j-parts) are then given by dsdpind + startind, which gives a pointer to the
       * first array-element belonging to this block and then the number of elements in this block is given to DSDP for iterating over it. */

      if ( penaltyparam > sdpisolver->epsilon && (! rbound) )
      {
         /* we need to compute the total number of nonzeros for the slack variable r, which equals the total number of diagonal entries */
         for (block = 0; block < nsdpblocks; block++)
            nrnonz += sdpblocksizes[block] - nremovedinds[block];
         assert( nrnonz >= 0 );

         /* indices given to DSDP, for this the elements in the lower triangular part of the matrix are labeled from 0 to n*(n+1)/2 -1 */
         BMS_CALL( BMSallocBlockMemoryArray(sdpisolver->blkmem, &dsdpind, sdpnnonz + nrnonz) ); /*lint !e776*/
         /* values given to DSDP, these will be multiplied by -1 because in DSDP -1* (sum A_i^j y_i - A_0) should be positive semidefinite */
         BMS_CALL( BMSallocBlockMemoryArray(sdpisolver->blkmem, &dsdpval, sdpnnonz + nrnonz) ); /*lint !e776*/
      }
      else
      {
         /* indices given to DSDP, for this the elements in the lower triangular part of the matrix are labeled from 0 to n*(n+1)/2 -1 */
         BMS_CALL( BMSallocBlockMemoryArray(sdpisolver->blkmem, &dsdpind, sdpnnonz) );
         /* values given to DSDP, these will be multiplied by -1 because in DSDP -1* (sum A_i^j y_i - A_0) should be positive semidefinite */
         BMS_CALL( BMSallocBlockMemoryArray(sdpisolver->blkmem, &dsdpval, sdpnnonz) );
      }

      ind = 0; /* this will be used for iterating over the nonzeroes */

      for (block = 0; block < nsdpblocks; block++)
      {
         for (i = 0; i < sdpisolver->nactivevars; i++)
         {
            /* we iterate over all non-fixed variables, so add them to the dsdp arrays for this block/var combination */
            v = sdpisolver->dsdptoinputmapper[i];

            /* find the position of variable v in this block */
            blockvar = -1;
            for (k = 0; k < sdpnblockvars[block]; k++)
            {
               if ( v == sdpvar[block][k] )
               {
                  blockvar = k;
                  break;
               }
            }

            startind = ind;

            if ( blockvar > -1 ) /* the variable exists in this block */
            {
               for (k = 0; k < sdpnblockvarnonz[block][blockvar]; k++)
               {
                  /* rows and cols with active nonzeros should not be removed */
                  assert( indchanges[block][sdprow[block][blockvar][k]] > -1 && indchanges[block][sdpcol[block][blockvar][k]] > -1 );

                  /* substract the number of removed indices before the row and col to get the indices after fixings */
                  dsdpind[ind] = compLowerTriangPos(sdprow[block][blockvar][k] - indchanges[block][sdprow[block][blockvar][k]],
                                                    sdpcol[block][blockvar][k] - indchanges[block][sdpcol[block][blockvar][k]]);
                  dsdpval[ind] = -1.0 * sdpval[block][blockvar][k];  /* *(-1) because in DSDP -1* (sum A_i^j y_i - A_0) should be positive semidefinite */
                  ind++;
               }

               /* sort the arrays for this matrix (by non decreasing indices) as this might help the solving time of DSDP */
               SCIPsortIntReal(dsdpind + startind, dsdpval + startind, sdpnblockvarnonz[block][blockvar]);

               assert( blockindchanges[block] > -1 ); /* we shouldn't insert into blocks we removed */

               /* i + 1 because DSDP starts counting the variables at 1, adding startind shifts the arrays to the first
                * nonzero belonging to this block and this variable */
               DSDP_CALL( SDPConeSetASparseVecMat(sdpisolver->sdpcone, block - blockindchanges[block], i + 1, sdpblocksizes[block] - nremovedinds[block],
                     1.0, 0, dsdpind + startind,dsdpval + startind, sdpnblockvarnonz[block][blockvar]));
            }
         }
      }

      if ( penaltyparam > sdpisolver->epsilon && (! rbound) )
      {
         startind = ind;
         /* add r * Identity for each block */
         for (block = 0; block < nsdpblocks; block++)
         {
            if ( blockindchanges[block] > -1 )
            {
               for (i = 0; i < sdpblocksizes[block] - nremovedinds[block]; i++)
               {
                  dsdpind[ind] = compLowerTriangPos(i, i);
                  dsdpval[ind] = -1.0; /* *(-1) because in DSDP -1* (sum A_i^j y_i - A_0 + r*I) should be positive semidefinite */
                  ind++;
               }
               DSDP_CALL( SDPConeSetASparseVecMat(sdpisolver->sdpcone, block - blockindchanges[block], sdpisolver->nactivevars + 1,
                     sdpblocksizes[block] - nremovedinds[block], 1.0, 0, dsdpind + ind - (sdpblocksizes[block] - nremovedinds[block]) ,
                     dsdpval + ind - (sdpblocksizes[block] - nremovedinds[block]), sdpblocksizes[block] - nremovedinds[block]) ); /*lint !e679*/
            }
         }
         assert( ind - startind == nrnonz );
      }
   }

   /* start inserting the constant matrix */
   if ( sdpconstnnonz > 0 )
   {
      assert( nsdpblocks > 0 );
      assert( sdpconstnblocknonz!= NULL );
      assert( sdpconstcol != NULL );
      assert( sdpconstrow != NULL );
      assert( sdpconstval != NULL );

      /* allocate memory */

      /* DSDP uses these for solving, so they may not be freed before the problem is solved. */

      /* indices given to DSDP, for this the elements in the lower triangular part of the matrix are labeled from 0 to n*(n+1)/2 -1 */
      BMS_CALL( BMSallocBlockMemoryArray(sdpisolver->blkmem, &dsdpconstind, sdpconstnnonz) );
      /* values given to DSDP, for this the original values are mutliplied by -1 because in DSDP -1* (sum A_i^j y_i - A_0) should be positive semidefinite */
      BMS_CALL( BMSallocBlockMemoryArray(sdpisolver->blkmem, &dsdpconstval, sdpconstnnonz) );

      ind = 0;

      for (block = 0; block < nsdpblocks; block++)
      {
         startind = ind; /* starting index of this block in the dsdpconst arrays */

         if ( sdpconstnblocknonz[block] > 0 )
         {
            /* insert the constant-nonzeros */
            for (i = 0; i < sdpconstnblocknonz[block]; i++)
            {
               /* rows and cols with nonzeros should not be removed */
               assert( indchanges[block][sdpconstrow[block][i]] > -1 && indchanges[block][sdpconstcol[block][i]] > -1 );

               /* substract the number of deleted indices before this to get the index after variable fixings */
               dsdpconstind[ind] = compLowerTriangPos(sdpconstrow[block][i] - indchanges[block][sdpconstrow[block][i]],
                                                      sdpconstcol[block][i] - indchanges[block][sdpconstcol[block][i]]);
               dsdpconstval[ind] = -1 * sdpconstval[block][i]; /* *(-1) because in DSDP -1* (sum A_i^j y_i - A_0^j) should be positive semidefinite */
               ind++;
            }

            /* sort the arrays for this Matrix (by non decreasing indices) as this might help the solving time of DSDP */
            SCIPsortIntReal(dsdpconstind + startind, dsdpconstval + startind, sdpconstnblocknonz[block]);

            assert( blockindchanges[block] > -1 ); /* we shouldn't insert into a block we removed */

            /* constant matrix is given as variable 0, the arrays are shifted to the first element of this block by adding
             * startind, ind - startind gives the number of elements for this block */
            DSDP_CALL( SDPConeSetASparseVecMat(sdpisolver->sdpcone, block - blockindchanges[block], 0, sdpblocksizes[block] - nremovedinds[block],
                                               1.0, 0, dsdpconstind + startind, dsdpconstval + startind, ind - startind));
         }
      }
   }

#ifdef SCIP_MORE_DEBUG
   SDPConeView2(sdpisolver->sdpcone);
#endif

   /* start inserting the LP constraints */
   if ( nlpcons > 0 || lpnnonz > 0 || ! SCIPsdpiSolverIsInfinity(sdpisolver, sdpisolver->objlimit) )
   {
      int nextcol;
      int* rowmapper; /* maps the lhs- and rhs-inequalities of the old LP-cons to their constraint numbers in DSDP */
      int pos;
      int newpos;
      int nlpineqs;

      assert( noldlpcons > 0 );
      assert( lprhs != NULL );
      assert( lpcol != NULL );
      assert( lprow != NULL );
      assert( lpval != NULL );

      /* allocate memory to save which lpconstraints are mapped to which index, entry 2i corresponds to the left hand side of row i, 2i+1 to the rhs */
      BMS_CALL( BMSallocBufferMemoryArray(sdpisolver->bufmem, &rowmapper, 2*noldlpcons) ); /*lint !e647*/

      /* compute the rowmapper and the number of inequalities we have to add to DSDP (as we have to split the ranged rows) */
      pos = 0;
      newpos = 0; /* the position in the lhs and rhs arrays */
      for (i = 0; i < noldlpcons; i++)
      {
        if ( rownactivevars[i] >= 2 )
        {
           if ( lplhs[newpos] > - SCIPsdpiSolverInfinity(sdpisolver) )
           {
              rowmapper[2*i] = pos; /*lint !e679*/
              pos++;
           }
           else
              rowmapper[2*i] = -1; /*lint !e679*/

           if ( lprhs[newpos] < SCIPsdpiSolverInfinity(sdpisolver) )
           {
              rowmapper[2*i + 1] = pos; /*lint !e679*/
              pos++;
           }
           else
              rowmapper[2*i + 1] = -1; /*lint !e679*/

           newpos++;
        }
        else
        {
           rowmapper[2*i] = -1; /*lint !e679*/
           rowmapper[2*i + 1] = -1; /*lint !e679*/
        }
      }
      nlpineqs = pos;
      assert( nlpineqs <= 2*nlpcons ); /* *2 comes from left- and right-hand-sides */

      /* memory allocation */

      /* these arrays are needed in DSDP during solving, so they may only be freed afterwards */
      /* dsdplpbegcol[i] gives the number of nonzeros in column 0 (right hand side) till i-1 (i going from 1 till m, with extra entries 0 (always 0)
       * and m+1 (always lpcons + lpnnonz)) */
      if ( penaltyparam > sdpisolver->epsilon && (! rbound) )
      {
         BMS_CALL( BMSallocBlockMemoryArray(sdpisolver->blkmem, &dsdplpbegcol, sdpisolver->nactivevars + 3) ); /* extra entry for r */ /*lint !e776*/
      }
      else
      {
         BMS_CALL( BMSallocBlockMemoryArray(sdpisolver->blkmem, &dsdplpbegcol, sdpisolver->nactivevars + 2) ); /*lint !e776*/
      }

      /* dsdplprow saves the row indices of the LP for DSDP */
      /* worst-case length is 2*lpnnonz + nlpineqs, because left- and right-hand-sides are also included in the vectorand we might have to duplicate the
       * non-zeros when splitting the ranged rows, this will be shortened after the exact length after fixings is known, in case we have an objective limit,
       * this is increased by one entry for the right-hand-side and at most nvars entries for the nonzeros, in case of rbound = FALSE, where we have to add
       * the entries for r ourselves, we have to add another nlpineqs for one entry for r for each active lp-constraint */
      if ( SCIPsdpiSolverIsInfinity(sdpisolver, sdpisolver->objlimit) )
      {
         if ( penaltyparam > sdpisolver->epsilon && (! rbound) )
         {
            BMS_CALL( BMSallocBlockMemoryArray(sdpisolver->blkmem, &dsdplprow, 2 * nlpineqs + 2*lpnnonz) ); /*lint !e647*/
         }
         else
         {
            BMS_CALL( BMSallocBlockMemoryArray(sdpisolver->blkmem, &dsdplprow, nlpineqs + 2*lpnnonz) ); /*lint !e647*/
         }
      }
      else
      {
         if ( penaltyparam > sdpisolver->epsilon && (! rbound) )
         {
            BMS_CALL( BMSallocBlockMemoryArray(sdpisolver->blkmem, &dsdplprow, (nlpineqs + 1) + 2*lpnnonz + nvars + nlpineqs) ); /*lint !e647*/
         }
         else
         {
            BMS_CALL( BMSallocBlockMemoryArray(sdpisolver->blkmem, &dsdplprow, (nlpineqs + 1) + 2*lpnnonz + nvars) ); /*lint !e647*/
         }
      }

      /* values given to DSDP */
      /* dsdplprow saves the row indices of the LP for DSDP */
      /* worst-case length is 2*lpnnonz + nlpineqs, because left- and right-hand-sides are also included in the vectorand we might have to duplicate the
       * non-zeros when splitting the ranged rows, this will be shortened after the exact length after fixings is known, in case we have an objective limit,
       * this is increased by one entry for the right-hand-side and at most nvars entries for the nonzeros, in case of rbound = FALSE, where we have to add
       * the entries for r ourselves, we have to add another nlpineqs for one entry for r for each active lp-constraint */
      if ( SCIPsdpiSolverIsInfinity(sdpisolver, sdpisolver->objlimit) )
      {
         if ( penaltyparam > sdpisolver->epsilon && (! rbound) )
         {
            BMS_CALL( BMSallocBlockMemoryArray(sdpisolver->blkmem, &dsdplpval, 2 * nlpineqs + 2*lpnnonz) ); /*lint !e647*/
         }
         else
         {
            BMS_CALL( BMSallocBlockMemoryArray(sdpisolver->blkmem, &dsdplpval, nlpineqs + 2*lpnnonz) ); /*lint !e647*/
         }
      }
      else
      {
         if ( penaltyparam > sdpisolver->epsilon && (! rbound) )
         {
            BMS_CALL( BMSallocBlockMemoryArray(sdpisolver->blkmem, &dsdplpval, (nlpineqs + 1) + 2*lpnnonz + nvars + nlpineqs) ); /*lint !e647*/
         }
         else
         {
            BMS_CALL( BMSallocBlockMemoryArray(sdpisolver->blkmem, &dsdplpval, (nlpineqs + 1) + 2*lpnnonz + nvars) ); /*lint !e647*/
         }
      }

      /* add all left- and right-hand-sides that are greater than zero (if their corresponding inequalities exist), the pos counter is increased for every
       * active row, to get the correct row numbers, the nonz-counter only if the lhs/rhs is unequal to zero and added to DSDP */
      dsdpnlpnonz = 0;
      pos = 0;
      for (i = 0; i < nlpcons; i++)
      {
         if ( lplhs[i] > - SCIPsdpiSolverInfinity(sdpisolver) )
         {
    	      if ( REALABS(lplhs[i]) > sdpisolver->epsilon )
    	      {
    		      dsdplprow[dsdpnlpnonz] = pos;
    		      dsdplpval[dsdpnlpnonz] = -lplhs[i]; /* we multiply by -1 because DSDP wants <= instead of >= */
    		      dsdpnlpnonz++;
    	      }
    	      pos++;
         }
         if ( lprhs[i] < SCIPsdpiSolverInfinity(sdpisolver) )
         {
            if ( REALABS(lprhs[i]) > sdpisolver->epsilon )
            {
               dsdplprow[dsdpnlpnonz] = pos;
               dsdplpval[dsdpnlpnonz] = lprhs[i];
               dsdpnlpnonz++;
            }
            pos++;
         }
      }
      assert( pos == nlpineqs );

      /* add the right-hand-side for the objective bound */
      if ( ! SCIPsdpiSolverIsInfinity(sdpisolver, sdpisolver->objlimit) )
      {
         if ( REALABS(sdpisolver->objlimit) > sdpisolver->epsilon )
         {
            dsdplprow[dsdpnlpnonz] = nlpcons; /* this is the last lp-constraint, as DSDP counts from 0 to nlpcons-1, this is number nlpcons */
            dsdplpval[dsdpnlpnonz] = sdpisolver->objlimit; /* as we want <= upper bound, this is the correct type of inequality for DSDP */
            dsdpnlpnonz++;
         }
      }

      /* now add the nonzeros */

      /* for this we have to sort the nonzeros by col first and then by row, as this is the sorting DSDP wants */
      sortColRow(lprow, lpcol, lpval, lpnnonz);

      /* iterate over all nonzeros to add the active ones to the dsdp arrays and compute dsdplpbegcol */
      nextcol = 0;
      dsdplpbegcol[0] = 0;
      for (i = 0; i < lpnnonz; i++)
      {
         /* if a new variable starts, set the corresponding dsdplpbegcol-entry */
         if ( lpcol[i] >= nextcol )
         {
            /* set the dsdplpbegcol entries, as there might be active variables which appear only in the sdp but not the lp-part, we also have to set
             * the starting values for all variables in between to the same value (as we also set the entry for the found variable, this for-queue
             * will always have at least one index in the index set) */
            for (j = nextcol; j <= lpcol[i]; j++)
            {
               if ( sdpisolver->inputtodsdpmapper[j] >= 0 )
               {
                  assert( ! (isFixed(sdpisolver, lb[j], ub[j])) );
                  dsdplpbegcol[sdpisolver->inputtodsdpmapper[j]] = dsdpnlpnonz;

                  /* add the entry to the objlimit-lp-constraint for the last variables */
                  if ( (! SCIPsdpiSolverIsInfinity(sdpisolver, sdpisolver->objlimit)) && (REALABS( obj[j] ) > sdpisolver->epsilon))
                  {
                     dsdplprow[dsdpnlpnonz] = nlpcons;
                     dsdplpval[dsdpnlpnonz] = obj[j];
                     dsdpnlpnonz++;
                  }
               }
            }
            nextcol = j;
         }

         /* add the nonzero, if it isn't fixed and the row isn't to be deleted (because it is only a variable bound) */
         if ( ! isFixed(sdpisolver, lb[lpcol[i]], ub[lpcol[i]]) )
         {
            /* add it to the inequality for the lhs of the ranged row, if it exists */
            if ( rowmapper[2*lprow[i]] > -1 ) /*lint !e679*/
            {
               /* the index is adjusted for deleted lp rows, also rows are numbered 0,...,nlpcons-1 in DSDP, as they are
                * here, nlpcons is added to the index as the first nlpcons entries correspond to the right hand sides */
               dsdplprow[dsdpnlpnonz] = rowmapper[2*lprow[i]]; /*lint !e679*/
               dsdplpval[dsdpnlpnonz] = -lpval[i]; /* - because dsdp wants <= instead of >= constraints */
               dsdpnlpnonz++;
            }
            /* add it to the inequality for the rhs of the ranged row, if it exists */
            if ( rowmapper[2*lprow[i] + 1] > -1 ) /*lint !e679*/
            {
               /* the index is adjusted for deleted lp rows, also rows are numbered 0,...,nlpcons-1 in DSDP, as they are
                * here, nlpcons is added to the index as the first nlpcons entries correspond to the right hand sides */
               dsdplprow[dsdpnlpnonz] = rowmapper[2*lprow[i] + 1]; /*lint !e679*/
               dsdplpval[dsdpnlpnonz] = lpval[i];
               dsdpnlpnonz++;
            }
         }
#ifndef SCIP_NDEBUG
         /* if this is an active nonzero for the row, it should have at least one active var */
         else
            assert( isFixed(sdpisolver, lb[lpcol[i]], ub[lpcol[i]]) || rownactivevars[lprow[i]] == 1 );
#endif
      }

      /* set the begcol array for all remaining variables (that aren't part of the LP-part), and also set the objlimit-constraint-entries */
      for (j = nextcol; j < nvars; j++)
      {
         if ( sdpisolver->inputtodsdpmapper[j] >= 0 )
         {
            assert( ! (isFixed(sdpisolver, lb[j], ub[j])) );
            dsdplpbegcol[sdpisolver->inputtodsdpmapper[j]] = dsdpnlpnonz;
            /* add the entry to the objlimit-lp-constraint for the last variables */
            if ( (! SCIPsdpiSolverIsInfinity(sdpisolver, sdpisolver->objlimit)) && (REALABS( obj[j] ) > sdpisolver->epsilon))
            {
               dsdplprow[dsdpnlpnonz] = nlpcons;
               dsdplpval[dsdpnlpnonz] = obj[j];
               dsdpnlpnonz++;
            }
         }
      }

      dsdplpbegcol[sdpisolver->nactivevars + 1] = dsdpnlpnonz; /*lint !e679*/

      /* add r * Identity if using a penalty formulation without a bound on r */
      if ( penaltyparam > sdpisolver->epsilon && (! rbound) )
      {
         for (i = 0; i < nlpineqs; i++)
         {
            dsdplprow[dsdpnlpnonz] = i;
            dsdplpval[dsdpnlpnonz] = -1.0; /* for >=-inequalities we would add a +1, but then we have to multiply these with -1 for DSDP */
            dsdpnlpnonz++;
         }
         dsdplpbegcol[sdpisolver->nactivevars + 2] = dsdpnlpnonz; /*lint !e679*/
      }

      /* free the memory for the rowshifts-array */
      BMSfreeBufferMemoryArray(sdpisolver->bufmem, &rowmapper); /*lint !e647, !e737*/

      /* shrink the dsdplp-arrays */
      if ( SCIPsdpiSolverIsInfinity(sdpisolver, sdpisolver->objlimit) )
      {
         if ( penaltyparam > sdpisolver->epsilon && (! rbound) )
         {
            BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &dsdplprow, 2*nlpineqs + 2*lpnnonz, dsdpnlpnonz) ); /*lint !e647*/
            BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &dsdplpval, 2*nlpineqs + 2*lpnnonz, dsdpnlpnonz) ); /*lint !e647*/
         }
         else
         {
            BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &dsdplprow, nlpineqs + 2*lpnnonz, dsdpnlpnonz) ); /*lint !e647*/
            BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &dsdplpval, nlpineqs + 2*lpnnonz, dsdpnlpnonz) ); /*lint !e647*/
         }
      }
      else
      {
         if ( penaltyparam > sdpisolver->epsilon && (! rbound) )
         {
            BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &dsdplprow, (nlpineqs + 1) + 2*lpnnonz + nvars + nlpineqs, dsdpnlpnonz) ); /*lint !e647*/
            BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &dsdplpval, (nlpineqs + 1) + 2*lpnnonz + nvars + nlpineqs, dsdpnlpnonz) ); /*lint !e647*/
         }
         else
         {
            BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &dsdplprow, (nlpineqs + 1) + 2*lpnnonz + nvars, dsdpnlpnonz) ); /*lint !e647*/
            BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &dsdplpval, (nlpineqs + 1) + 2*lpnnonz + nvars, dsdpnlpnonz) ); /*lint !e647*/
         }
      }
      /* add the arrays to dsdp (in this case we need no additional if for the penalty version without bounds, as we already added the extra var,
       * so DSDP knows, that there is an additional entry in dsdplpbegcol which then gives the higher number of nonzeros) */
      if ( SCIPsdpiSolverIsInfinity(sdpisolver, sdpisolver->objlimit) )
      {
         DSDP_CALL( LPConeSetData(sdpisolver->lpcone, nlpineqs, dsdplpbegcol, dsdplprow, dsdplpval) );
      }
      else
      {
         DSDP_CALL( LPConeSetData(sdpisolver->lpcone, nlpineqs + 1, dsdplpbegcol, dsdplprow, dsdplpval) );
      }
#ifdef SCIP_MORE_DEBUG
      LPConeView(sdpisolver->lpcone);
#endif
   }

   SCIPdebugMessage("Calling DSDP-Solve for SDP (%d) \n", sdpisolver->sdpcounter);

   DSDP_CALL( DSDPSetGapTolerance(sdpisolver->dsdp, sdpisolver->gaptol) );  /* set DSDP's tolerance for duality gap */
   DSDP_CALL( DSDPSetRTolerance(sdpisolver->dsdp, sdpisolver->sdpsolverfeastol) );    /* set DSDP's tolerance for the SDP-constraints */
   if ( sdpisolver-> sdpinfo )
   {
      DSDP_CALL( DSDPSetStandardMonitor(sdpisolver->dsdp, 1) );   /* output DSDP information after every 1 iteration */
   }

   /* set the penalty parameter (only if rbound = TRUE, otherwise we had to add everything ourselves) */
   if ( penaltyparam >= sdpisolver->epsilon && rbound ) /* in sdpisolverSolve this is called with an exact 0 */
   {
      DSDP_CALL( DSDPSetPenaltyParameter(sdpisolver->dsdp, penaltyparam) );
      DSDP_CALL( DSDPUsePenalty(sdpisolver->dsdp, 1) );
   }
   else
   {
      /* set the penalty parameter to the default value */
      DSDP_CALL( DSDPSetPenaltyParameter(sdpisolver->dsdp, sdpisolver->penaltyparam) );
   }

   /* set the starting solution */
   if ( start != NULL )
   {
      for (i = 1; i <= sdpisolver->nactivevars; i++) /* we iterate over the variables in DSDP */
      {
         DSDP_CALL( DSDPSetY0(sdpisolver->dsdp, i, start[sdpisolver->dsdptoinputmapper[i]]) );
      }
   }

   /* start the solving process */
   DSDP_CALLM( DSDPSetup(sdpisolver->dsdp) );
   if ( ! SCIPsdpiSolverIsInfinity(sdpisolver, timelimit) )
   {
      DSDP_CALLM( DSDPSetMonitor(sdpisolver->dsdp, checkTimeLimitDSDP, (void*) &timings) );
   }
   DSDP_CALL( DSDPSolve(sdpisolver->dsdp) );

   sdpisolver->nsdpcalls++;
   DSDP_CALL( DSDPGetIts(sdpisolver->dsdp, &(sdpisolver->niterations)) );

   /* check if solving was stopped because of the time limit */
   if ( timings.stopped )
   {
      sdpisolver->timelimit = TRUE;
      sdpisolver->solved = FALSE;
   }
   else
   {
      sdpisolver->timelimit = FALSE;
      DSDP_CALL( DSDPComputeX(sdpisolver->dsdp) ); /* computes X and determines feasibility and unboundedness of the solution */
      sdpisolver->solved = TRUE;
   }

   /* if the problem has been stably solved but did not reach the required feasibility tolerance, even though the solver
    * reports feasibility, resolve it with adjusted tolerance */
   feastol = sdpisolver->sdpsolverfeastol;

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
         SCIPdebugMessage("Solution feasible for DSDP but outside feasibility tolerance, changing SDPA feasibility tolerance from %f to %f\n",
               feastol, feastol * INFEASFEASTOLCHANGE);
         feastol *= INFEASFEASTOLCHANGE;

         if ( feastol >= INFEASMINFEASTOL )
         {
            /* update settings */
            DSDP_CALL( DSDPSetRTolerance(sdpisolver->dsdp, feastol) );    /* set DSDP's tolerance for the SDP-constraints */

            DSDP_CALL( DSDPSolve(sdpisolver->dsdp) );

            /* update number of SDP-iterations and -calls */
            sdpisolver->nsdpcalls++;
            DSDP_CALL( DSDPGetIts(sdpisolver->dsdp, &newiterations) );
            sdpisolver->niterations += newiterations;

            /* check if solving was stopped because of the time limit */
            if ( timings.stopped )
            {
               sdpisolver->timelimit = TRUE;
               sdpisolver->solved = FALSE;
            }
            else
            {
               sdpisolver->timelimit = FALSE;
               DSDP_CALL( DSDPComputeX(sdpisolver->dsdp) ); /* computes X and determines feasibility and unboundedness of the solution */
               sdpisolver->solved = TRUE;
            }
         }
         else
         {
            sdpisolver->solved = FALSE;
            SCIPmessagePrintInfo(sdpisolver->messagehdlr, "SDPA failed to reach required feasibility tolerance! \n");
         }
      }
      else
         break;
   }

   /*these arrays were used to give information to DSDP and were needed during solving and for computing X, so they may only be freed now*/
   if ( sdpconstnnonz > 0 )
   {
      BMSfreeBlockMemoryArray(sdpisolver->blkmem, &dsdpconstval, sdpconstnnonz);/*lint !e737 */
      BMSfreeBlockMemoryArray(sdpisolver->blkmem, &dsdpconstind, sdpconstnnonz);/*lint !e737 */
   }

   if ( sdpnnonz > 0 )
   {
      if ( penaltyparam > sdpisolver->epsilon && (! rbound) )
      {
         BMSfreeBlockMemoryArray(sdpisolver->blkmem, &dsdpval, sdpnnonz + nrnonz); /*lint !e737, !e776*/
         BMSfreeBlockMemoryArray(sdpisolver->blkmem, &dsdpind, sdpnnonz + nrnonz); /*lint !e737, !e776*/
      }
      else
      {
         BMSfreeBlockMemoryArray(sdpisolver->blkmem, &dsdpval, sdpnnonz);/*lint !e737 */
         BMSfreeBlockMemoryArray(sdpisolver->blkmem, &dsdpind, sdpnnonz);/*lint !e737 */
      }
   }

   if ( nlpcons > 0 || lpnnonz > 0 )
   {
      BMSfreeBlockMemoryArray(sdpisolver->blkmem, &dsdplpval, dsdpnlpnonz);/*lint !e644, !e737*/
      BMSfreeBlockMemoryArray(sdpisolver->blkmem, &dsdplprow, dsdpnlpnonz);/*lint !e737 */
      if ( penaltyparam > sdpisolver->epsilon && (! rbound) )
      {
         BMSfreeBlockMemoryArray(sdpisolver->blkmem, &dsdplpbegcol, sdpisolver->nactivevars + 3); /*lint !e737, !e776*/
      }
      else
      {
         BMSfreeBlockMemoryArray(sdpisolver->blkmem, &dsdplpbegcol, sdpisolver->nactivevars + 2); /*lint !e737, !e776*/
      }
   }

#ifdef SCIP_DEBUG
   DSDP_CALL( DSDPStopReason(sdpisolver->dsdp, &reason) );

   switch ( reason )
   {
   case DSDP_CONVERGED:
      SCIPdebugMessage("DSDP converged!\n");
      break;

   case DSDP_INFEASIBLE_START:
      SCIPdebugMessage("DSDP started with an infeasible point!\n");
      break;

   case DSDP_SMALL_STEPS:
      SCIPdebugMessage("Short step lengths created by numerical difficulties prevented progress in DSDP!\n");
      break;

   case DSDP_INDEFINITE_SCHUR_MATRIX:
      SCIPdebugMessage("Schur Matrix in DSDP was indefinite but should have been positive semidefinite!\n");
      break;

   case DSDP_MAX_IT:
      SCIPdebugMessage("DSDP reached maximum number of iterations!\n");
      break;

   case DSDP_NUMERICAL_ERROR:
      SCIPdebugMessage("A numerical error occured in DSDP!\n");
      break;

   case DSDP_UPPERBOUND:
      SCIPdebugMessage("Dual objective value in DSDP reached upper bound.\n");
      break;

   case DSDP_USER_TERMINATION:
      SCIPdebugMessage("DSDP didn't stop solving, did you?\n");
      break;

   case CONTINUE_ITERATING:
      SCIPdebugMessage("DSDP wants to continue iterating but somehow was stopped!\n");
      break;

   default:
      SCIPdebugMessage("Unknown stopping reason in DSDP!\n");
      break;
   }
#endif

   if ( penaltyparam >= sdpisolver->epsilon && sdpisolver->solved )
   {
      if ( rbound )
      {
         /* in this case we used the penalty-formulation of DSDP, so we can check their value of r */
         SCIP_Real rval;
         SCIP_Real trace;

         DSDP_CALL( DSDPGetR(sdpisolver->dsdp, &rval) );

         *feasorig = (rval < sdpisolver->feastol );

         /* only set sdpisolver->feasorig to true if we solved with objective, because only in this case we want to compute
          * the objective value by hand since it is numerically more stable then the result returned by DSDP */
         if ( withobj )
            sdpisolver->feasorig = *feasorig;

         /* if r > 0 or we are in debug mode, also check the primal bound */
#ifdef NDEBUG
         if ( ! *feasorig )
         {
#endif
            if ( penaltybound != NULL )
            {
               SCIPdebugMessage("Solution not feasible in original problem, r = %f\n", rval);

               /* get the trace of X to compare it with the penalty parameter */
               DSDP_CALL( DSDPGetTraceX(sdpisolver->dsdp, &trace) );

#if 0 /* DSDP doesn't seem to adhere to its own feasiblity tolerance */
               assert( trace < penaltyparam + sdpisolver->feastol ); /* solution should be primal feasible */
#endif

               /* if the relative gap is smaller than the tolerance, we return equality */
               if ( (penaltyparam - trace) / penaltyparam < PENALTYBOUNDTOL )
               {
                  *penaltybound = TRUE;
                  SCIPdebugMessage("Tr(X) = %f == %f = Gamma, penalty formulation not exact, Gamma should be increased or problem is infeasible\n",
                     trace, penaltyparam);
               }
               else
                  *penaltybound = FALSE;
            }
#ifdef NDEBUG
         }
#endif
      }
      else
      {
         SCIP_Real* dsdpsol;
         SCIP_Real trace;

         BMS_CALL( BMSallocBufferMemoryArray(sdpisolver->bufmem, &dsdpsol, sdpisolver->nactivevars + 1) ); /*lint !e776*/
         /* last entry of DSDPGetY needs to be the number of variables, will return an error otherwise */
         DSDP_CALL( DSDPGetY(sdpisolver->dsdp, dsdpsol, sdpisolver->nactivevars + 1) );

         *feasorig = (dsdpsol[sdpisolver->nactivevars] < sdpisolver->feastol); /* r is the last variable in DSDP, so the last entry gives us the value */
#ifdef NDEBUG
         if ( ! *feasorig )
         {
#endif
            if ( penaltybound != NULL )
            {
               SCIPdebugMessage("Solution not feasible in original problem, r = %f\n", dsdpsol[sdpisolver->nactivevars]);

               /* get the trace of X to compare it with the penalty parameter */
               DSDP_CALL( DSDPGetTraceX(sdpisolver->dsdp, &trace) );

#if 0 /* DSDP doesn't seem to adhere to its own feasiblity tolerance */
               assert( trace < penaltyparam + sdpisolver->feastol ); /* solution should be primal feasible */
#endif

               /* if the relative gap is smaller than the tolerance, we return equality */
               if ( (penaltyparam - trace) / penaltyparam < PENALTYBOUNDTOL )
               {
                  *penaltybound = TRUE;
                  SCIPdebugMessage("Tr(X) = %f == %f = Gamma, penalty formulation not exact, Gamma should be increased or problem is infeasible\n",
                        trace, penaltyparam);
               }
               else
                  *penaltybound = FALSE;
            }
#ifdef NDEBUG
         }
#endif
          BMSfreeBufferMemoryArray(sdpisolver->bufmem, &dsdpsol);
      }
   }

   return SCIP_OKAY;
}/*lint !e715 */
/**@} */




/*
 * Solution Information Methods
 */

/**@name Solution Information Methods */
/**@{ */

/** returns whether a solve method was called after the last modification of the SDP */
SCIP_Bool SCIPsdpiSolverWasSolved(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP-solver interface */
   )
{
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
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP-solver interface */
   )
{
   DSDPSolutionType pdfeasible;

   assert( sdpisolver != NULL );
   CHECK_IF_SOLVED_BOOL( sdpisolver );

   DSDP_CALL_BOOL( DSDPGetSolutionType(sdpisolver->dsdp, &pdfeasible) );

   if ( pdfeasible == DSDP_PDUNKNOWN )
      return FALSE;

   return TRUE;
}

/** gets information about primal and dual feasibility of the current SDP solution */
SCIP_RETCODE SCIPsdpiSolverGetSolFeasibility(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   SCIP_Bool*            primalfeasible,     /**< stores primal feasibility status */
   SCIP_Bool*            dualfeasible        /**< stores dual feasibility status */
   )
{
   DSDPSolutionType pdfeasible;

   assert( sdpisolver != NULL );
   assert( primalfeasible != NULL );
   assert( dualfeasible != NULL );
   CHECK_IF_SOLVED( sdpisolver );

   DSDP_CALL( DSDPGetSolutionType(sdpisolver->dsdp, &pdfeasible) );

   switch ( pdfeasible )
   {
   case DSDP_PDFEASIBLE:
      *primalfeasible = TRUE;
      *dualfeasible = TRUE;
      break;

   case DSDP_UNBOUNDED:
      *primalfeasible = FALSE;
      *dualfeasible = TRUE;
      break;

   case DSDP_INFEASIBLE:
      *primalfeasible = TRUE;
      *dualfeasible = FALSE;
      break;

   default: /* should only include DSDP_PDUNKNOWN */
      SCIPerrorMessage("DSDP doesn't know if primal and dual solutions are feasible\n");
      return SCIP_LPERROR;
   }

   return SCIP_OKAY;
}

/** returns TRUE iff SDP is proven to be primal unbounded,
 *  returns FALSE with a debug-message if the solver could not determine feasibility */
SCIP_Bool SCIPsdpiSolverIsPrimalUnbounded(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP-solver interface */
   )
{
   DSDPSolutionType pdfeasible;

   assert( sdpisolver != NULL );
   CHECK_IF_SOLVED_BOOL( sdpisolver );

   DSDP_CALL_BOOL( DSDPGetSolutionType(sdpisolver->dsdp, &pdfeasible) );
   if ( pdfeasible == DSDP_PDUNKNOWN )
   {
/*      SCIPerrorMessage("DSDP doesn't know if primal and dual solutions are feasible");
      SCIPABORT();
      return SCIP_LPERROR;*/
      SCIPdebugMessage("DSDP doesn't know if primal and dual solutions are feasible.");
      return FALSE;
   }
   else if ( pdfeasible == DSDP_INFEASIBLE )
      return TRUE;

   return FALSE;
}

/** returns TRUE iff SDP is proven to be primal infeasible,
 *  returns FALSE with a debug-message if the solver could not determine feasibility */
SCIP_Bool SCIPsdpiSolverIsPrimalInfeasible(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP-solver interface */
   )
{
   DSDPSolutionType pdfeasible;

   assert( sdpisolver != NULL );
   CHECK_IF_SOLVED_BOOL( sdpisolver );

   DSDP_CALL_BOOL( DSDPGetSolutionType(sdpisolver->dsdp, &pdfeasible) );
   if ( pdfeasible == DSDP_PDUNKNOWN )
   {
/*      SCIPerrorMessage("DSDP doesn't know if primal and dual solutions are feasible");
      SCIPABORT();
      return SCIP_LPERROR;*/
      SCIPdebugMessage("DSDP doesn't know if primal and dual solutions are feasible");
      return FALSE;
   }
   else if ( pdfeasible == DSDP_UNBOUNDED )
      return TRUE;

   return FALSE;
}

/** returns TRUE iff SDP is proven to be primal feasible,
 *  returns FALSE with a debug-message if the solver could not determine feasibility */
SCIP_Bool SCIPsdpiSolverIsPrimalFeasible(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP-solver interface */
   )
{
   DSDPSolutionType pdfeasible;

   assert( sdpisolver != NULL );
   CHECK_IF_SOLVED_BOOL( sdpisolver );

   DSDP_CALL_BOOL( DSDPGetSolutionType(sdpisolver->dsdp, &pdfeasible) );
   if ( pdfeasible == DSDP_PDUNKNOWN )
   {
      SCIPdebugMessage("DSDP doesn't know if primal and dual solutions are feasible");
      return FALSE;
   }
   else if ( pdfeasible == DSDP_UNBOUNDED )
      return FALSE;

   return TRUE;
}

/** returns TRUE iff SDP is proven to be dual unbounded,
 *  returns FALSE with a debug-message if the solver could not determine feasibility */
SCIP_Bool SCIPsdpiSolverIsDualUnbounded(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP-solver interface */
   )
{
   DSDPSolutionType pdfeasible;

   assert( sdpisolver != NULL );
   CHECK_IF_SOLVED_BOOL( sdpisolver );

   DSDP_CALL_BOOL( DSDPGetSolutionType(sdpisolver->dsdp, &pdfeasible) );
   if ( pdfeasible == DSDP_PDUNKNOWN )
   {
      SCIPdebugMessage("DSDP doesn't know if primal and dual solutions are feasible");
      return FALSE;
   }
   else if ( pdfeasible == DSDP_UNBOUNDED )
      return TRUE;

   return FALSE;
}

/** returns TRUE iff SDP is proven to be dual infeasible,
 *  returns FALSE with a debug-message if the solver could not determine feasibility */
SCIP_Bool SCIPsdpiSolverIsDualInfeasible(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP-solver interface */
   )
{
   DSDPSolutionType pdfeasible;

   assert( sdpisolver != NULL );
   CHECK_IF_SOLVED_BOOL( sdpisolver );

   DSDP_CALL_BOOL(DSDPGetSolutionType(sdpisolver->dsdp, &pdfeasible));

   if ( pdfeasible == DSDP_PDUNKNOWN )
   {
      SCIPdebugMessage("DSDP doesn't know if primal and dual solutions are feasible");
      return FALSE;
   }
   else if ( pdfeasible == DSDP_INFEASIBLE )
      return TRUE;

   return FALSE;
}

/** returns TRUE iff SDP is proven to be dual feasible,
 *  returns FALSE with a debug-message if the solver could not determine feasibility */
SCIP_Bool SCIPsdpiSolverIsDualFeasible(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP-solver interface */
   )
{
   DSDPSolutionType pdfeasible;

   assert( sdpisolver != NULL );
   CHECK_IF_SOLVED_BOOL( sdpisolver );

   DSDP_CALL_BOOL( DSDPGetSolutionType(sdpisolver->dsdp, &pdfeasible) );

   if ( pdfeasible == DSDP_PDUNKNOWN )
   {
      SCIPdebugMessage("DSDP doesn't know if primal and dual solutions are feasible");
      return FALSE;
   }
   else if ( pdfeasible == DSDP_INFEASIBLE )
      return FALSE;

   return TRUE;
}

/** returns TRUE iff the solver converged */
SCIP_Bool SCIPsdpiSolverIsConverged(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP-solver interface */
   )
{
   DSDPTerminationReason reason;

   assert( sdpisolver != NULL );

   if ( sdpisolver->timelimit )
      return FALSE;

   if ( ! sdpisolver->solved )
      return FALSE;

   DSDP_CALL_BOOL( DSDPStopReason(sdpisolver->dsdp, &reason) );

   if ( reason == DSDP_CONVERGED )
      return TRUE;

   return FALSE;
}

/** returns TRUE iff the objective limit was reached */
SCIP_Bool SCIPsdpiSolverIsObjlimExc(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP-solver interface */
   )
{/*lint --e{715}*/
   SCIPdebugMessage("Method not implemented for DSDP, as objective limit is given as an ordinary LP-constraint, so in case the objective limit was "
         "exceeded, the problem will be reported as infeasible ! \n");

   return FALSE;
}

/** returns TRUE iff the iteration limit was reached */
SCIP_Bool SCIPsdpiSolverIsIterlimExc(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP-solver interface */
   )
{
   DSDPTerminationReason reason;

   assert( sdpisolver != NULL );
   CHECK_IF_SOLVED_BOOL( sdpisolver );

   DSDP_CALL_BOOL( DSDPStopReason(sdpisolver->dsdp, &reason) );

   if ( reason == DSDP_MAX_IT )
      return TRUE;

   return FALSE;
}

/** returns TRUE iff the time limit was reached */
SCIP_Bool SCIPsdpiSolverIsTimelimExc(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP-solver interface */
   )
{
   assert( sdpisolver != NULL );

   return sdpisolver->timelimit;
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
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP-solver interface */
   )
{
   DSDPTerminationReason reason;
   int dsdpreturn;

   assert( sdpisolver != NULL );

   if ( sdpisolver->dsdp == NULL || (! sdpisolver->solved) )
      return -1;

   if ( sdpisolver->timelimit )
      return 5;

   dsdpreturn = DSDPStopReason(sdpisolver->dsdp, &reason);

   if (dsdpreturn != 0)
   {
      SCIPerrorMessage("DSDP-Error <%d> in function call.\n", dsdpreturn);
      return 7;
   }

   switch ( reason )/*lint --e{788}*/
   {
   case DSDP_CONVERGED:
      return 0;

   case DSDP_INFEASIBLE_START:
      return 1;

   case DSDP_SMALL_STEPS:
      return 2;

   case DSDP_INDEFINITE_SCHUR_MATRIX:
      return 2;

   case DSDP_MAX_IT:
      return 4;

   case DSDP_NUMERICAL_ERROR:
      return 2;

   case DSDP_UPPERBOUND:
      return 3;

   case DSDP_USER_TERMINATION:
      return 5;

   default:
      return 7;
   }
}

/** returns TRUE iff SDP was solved to optimality, meaning the solver converged and returned primal and dual feasible solutions */
SCIP_Bool SCIPsdpiSolverIsOptimal(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP-solver interface */
   )
{
   assert( sdpisolver != NULL );
   return (SCIPsdpiSolverIsConverged(sdpisolver) && SCIPsdpiSolverIsPrimalFeasible(sdpisolver) && SCIPsdpiSolverIsDualFeasible(sdpisolver));
}

/** returns TRUE iff SDP was solved to optimality or some other status was reached
 *  that is still acceptable inside a Branch & Bound framework */
SCIP_Bool SCIPsdpiSolverIsAcceptable(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP-solver interface */
   )
{
   assert( sdpisolver != NULL );

   return SCIPsdpiSolverIsConverged(sdpisolver);
}

/** tries to reset the internal status of the SDP-solver in order to ignore an instability of the last solving call */
SCIP_RETCODE SCIPsdpiSolverIgnoreInstability(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   SCIP_Bool*            success             /**< pointer to store, whether the instability could be ignored */
   )
{/*lint --e{715}*/
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_LPERROR;
}

/** gets objective value of solution */
SCIP_RETCODE SCIPsdpiSolverGetObjval(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   SCIP_Real*            objval              /**< pointer to store the objective value */
   )
{
   SCIP_Real* dsdpsol;
   int dsdpnvars;

   assert( sdpisolver != NULL );
   assert( objval != NULL );
   CHECK_IF_SOLVED( sdpisolver );

   dsdpnvars = sdpisolver->penaltyworbound ? sdpisolver->nactivevars + 1 : sdpisolver->nactivevars; /* in the first case we added r as an explicit var */

   if ( sdpisolver->penalty && ( ! sdpisolver->feasorig ))
   {
      /* in this case we cannot really trust the solution given by DSDP, since changes in the value of r much less than epsilon can
       * cause huge changes in the objective, so using the objective value given by DSDP is numerically more stable */
      DSDP_CALL( DSDPGetDObjective(sdpisolver->dsdp, objval) );
      *objval = -1*(*objval); /*DSDP maximizes instead of minimizing, so the objective values were multiplied by -1 when inserted */
   }
   else
   {
      int v;

      /* since the objective value given by DSDP sometimes differs slightly from the correct value for the given solution,
       * we get the solution from DSDP and compute the correct objective value */
      BMS_CALL( BMSallocBlockMemoryArray(sdpisolver->blkmem, &dsdpsol, dsdpnvars) );
      DSDP_CALL( DSDPGetY(sdpisolver->dsdp, dsdpsol, dsdpnvars) ); /* last entry needs to be the number of variables, will return an error otherwise */

      /* use the solution to compute the correct objective value */
      *objval = 0.0;
      for (v = 0; v < sdpisolver->nactivevars; v++)
      {
         if ( dsdpsol[v] > sdpisolver->epsilon )
            *objval += sdpisolver->objcoefs[v] * dsdpsol[v];
      }
   }

   /* as we didn't add the fixed (lb = ub) variables to dsdp, we have to add their contributions to the objective as well */
   *objval += sdpisolver->fixedvarsobjcontr;

   if ( ( ! sdpisolver->penalty ) || sdpisolver->feasorig )
   {
      BMSfreeBlockMemoryArray(sdpisolver->blkmem, &dsdpsol, dsdpnvars);/*lint !e737 */
   }

   return SCIP_OKAY;
}

/** gets dual solution vector for feasible SDPs
 *
 *  If dualsollength isn't equal to the number of variables this will return the needed length and a debug message is thrown.
 */
SCIP_RETCODE SCIPsdpiSolverGetSol(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   SCIP_Real*            objval,             /**< pointer to store the objective value, may be NULL if not needed */
   SCIP_Real*            dualsol,            /**< pointer to store the dual solution vector, may be NULL if not needed */
   int*                  dualsollength       /**< length of the dual sol vector, must be 0 if dualsol is NULL, if this is less than the number
                                              *   of variables in the SDP, a DebugMessage will be thrown and this is set to the needed value */
   )
{
   SCIP_Real* dsdpsol;
   int v;
   int dsdpnvars;

   assert( sdpisolver != NULL );
   assert( dualsollength != NULL );
   CHECK_IF_SOLVED( sdpisolver );

   dsdpnvars = sdpisolver->penaltyworbound ? sdpisolver->nactivevars + 1 : sdpisolver->nactivevars; /* in the first case we added r as an explicit var */

   if ( *dualsollength > 0 )
   {
      assert( dualsol != NULL );
      if ( *dualsollength < sdpisolver->nvars )
      {
         SCIPdebugMessage("The given array in SCIPsdpiSolverGetSol only had length %d, but %d was needed", *dualsollength, sdpisolver->nvars);
         *dualsollength = sdpisolver->nvars;

         return SCIP_OKAY;
      }

      BMS_CALL( BMSallocBlockMemoryArray(sdpisolver->blkmem, &dsdpsol, dsdpnvars) );
      DSDP_CALL( DSDPGetY(sdpisolver->dsdp, dsdpsol, dsdpnvars) ); /* last entry needs to be the number of variables, will return an error otherwise */

      /* insert the entries into dualsol, for non-fixed vars we copy those from dsdp, the rest are the saved entries from inserting (they equal lb=ub) */
      for (v = 0; v < sdpisolver->nvars; v++)
      {
         if (sdpisolver->inputtodsdpmapper[v] > -1)
         {
            /* minus one because the inputtodsdpmapper gives the dsdp indices which start at one, but the array starts at 0 */
            dualsol[v] = dsdpsol[sdpisolver->inputtodsdpmapper[v] - 1];
         }
         else
         {
            /* this is the value that was saved when inserting, as this variable has lb=ub */
            dualsol[v] = sdpisolver->fixedvarsval[(-1 * sdpisolver->inputtodsdpmapper[v]) - 1]; /*lint !e679*/
         }
      }

      if ( objval != NULL )
      {
         if ( sdpisolver->penalty && ( ! sdpisolver->feasorig ))
         {
            /* in this case we cannot really trust the solution given by DSDP, since changes in the value of r much less than epsilon can
             * cause huge changes in the objective, so using the objective value given by DSDP is numerically more stable */
            DSDP_CALL( DSDPGetDObjective(sdpisolver->dsdp, objval) );
            *objval = -1*(*objval); /*DSDP maximizes instead of minimizing, so the objective values were multiplied by -1 when inserted */
         }
         else
         {
            /* use the solution to compute the correct objective value */
            *objval = 0.0;
            for (v = 0; v < sdpisolver->nactivevars; v++)
            {
               if ( REALABS(dsdpsol[v]) > sdpisolver->epsilon )
                  *objval += sdpisolver->objcoefs[v] * dsdpsol[v];
            }
         }

         /* as we didn't add the fixed (lb = ub) variables to dsdp, we have to add their contributions to the objective as well */
         *objval += sdpisolver->fixedvarsobjcontr;
      }

      BMSfreeBlockMemoryArray(sdpisolver->blkmem, &dsdpsol, dsdpnvars);/*lint !e737 */
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
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   SCIP_Real*            lbvars,             /**< pointer to store the values of the variables corresponding to lower bounds in the dual problems */
   SCIP_Real*            ubvars,             /**< pointer to store the values of the variables corresponding to upper bounds in the dual problems */
   int*                  arraylength         /**< input: length of lbvars and ubvars <br>
                                              *   output: number of elements inserted into lbvars/ubvars (or needed length if it wasn't sufficient) */
   )
{
   SCIP_Real* lbvarsdsdp;
   SCIP_Real* ubvarsdsdp;
   int i;

   assert( sdpisolver != NULL );
   assert( lbvars != NULL );
   assert( ubvars != NULL );
   assert( arraylength != NULL );
   assert( *arraylength >= 0 );
   CHECK_IF_SOLVED( sdpisolver );

   /* check if the arrays are long enough */
   if ( *arraylength < sdpisolver->nvars )
   {
      *arraylength = sdpisolver->nvars;
      SCIPdebugMessage("Insufficient length of array in SCIPsdpiSolverGetPrimalBoundVars (gave %d, needed %d)\n", *arraylength, sdpisolver->nvars);
      return SCIP_OKAY;
   }

   /* allocate memory for the arrays given to DSDP */
   BMS_CALL( BMSallocBlockMemoryArray(sdpisolver->blkmem, &lbvarsdsdp, sdpisolver->nactivevars) );
   BMS_CALL( BMSallocBlockMemoryArray(sdpisolver->blkmem, &ubvarsdsdp, sdpisolver->nactivevars) );

   /* get the values for the active variables from DSDP */
   DSDP_CALL( BConeCopyX(sdpisolver->bcone, lbvarsdsdp, ubvarsdsdp, sdpisolver->nactivevars) );

   /* copy them to the right spots of lbvars & ubvars */
   for (i = 0; i < sdpisolver->nvars; i++)
   {
      if ( sdpisolver->inputtodsdpmapper[i] < 0 )
      {
         /* if the variable was fixed, it didn't exist in the relaxation, so we set the value to 0
          * (as DSDP already uses this value for unbounded vars) */
         lbvars[i] = 0;
         ubvars[i] = 0;
      }
      else
      {
         lbvars[i] = lbvarsdsdp[sdpisolver->inputtodsdpmapper[i] - 1];
         ubvars[i] = ubvarsdsdp[sdpisolver->inputtodsdpmapper[i] - 1];
      }
   }

   /* free allocated memory */
   BMSfreeBlockMemoryArrayNull(sdpisolver->blkmem, &ubvarsdsdp, sdpisolver->nactivevars);
   BMSfreeBlockMemoryArrayNull(sdpisolver->blkmem, &lbvarsdsdp, sdpisolver->nactivevars);

   return SCIP_OKAY;
}

/** gets the number of SDP iterations of the last solve call */
SCIP_RETCODE SCIPsdpiSolverGetIterations(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   int*                  iterations          /**< pointer to store the number of iterations of the last solve call */
   )
{
   assert( sdpisolver != NULL );
   assert( iterations != NULL );

   if ( sdpisolver->timelimitinitial )
      *iterations = 0;
   else
      *iterations = sdpisolver->niterations;

   return SCIP_OKAY;
}

/** gets the number of calls to the SDP-solver for the last solve call */
SCIP_RETCODE SCIPsdpiSolverGetSdpCalls(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   int*                  calls               /**< pointer to store the number of calls to the SDP-solver for the last solve call */
   )
{
   assert( sdpisolver != NULL );
   assert( calls != NULL );

   if ( sdpisolver->timelimitinitial )
      *calls = 0;
   else
      *calls = sdpisolver->nsdpcalls;

   return SCIP_OKAY;
}

/** gets the settings used by the SDP solver for the last solve call */
SCIP_RETCODE SCIPsdpiSolverSettingsUsed(
   SCIP_SDPISOLVER*      sdpisolver,         /**< SDP-solver interface */
   SCIP_SDPSOLVERSETTING* usedsetting        /**< the setting used by the SDP-solver */
   )
{
   assert( sdpisolver != NULL );
   assert( usedsetting != NULL );

   if ( ! SCIPsdpiSolverIsAcceptable(sdpisolver) )
      *usedsetting = SCIP_SDPSOLVERSETTING_UNSOLVED;
   else
      *usedsetting = sdpisolver->usedsetting;

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
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to an SDP-solver interface */
   )
{
   return 1E+20; /* default infinity from SCIP */
}/*lint !e715 */

/** checks if given value is treated as (plus or minus) infinity in the SDP-solver */
SCIP_Bool SCIPsdpiSolverIsInfinity(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   SCIP_Real             val                 /**< value to be checked for infinity */
   )
{
   return ((val <= -SCIPsdpiSolverInfinity(sdpisolver)) || (val >= SCIPsdpiSolverInfinity(sdpisolver)));
}

/** gets floating point parameter of SDP-Solver */
SCIP_RETCODE SCIPsdpiSolverGetRealpar(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   SCIP_SDPPARAM         type,               /**< parameter number */
   SCIP_Real*            dval                /**< buffer to store the parameter value */
   )
{
   assert( sdpisolver != NULL );
   assert( dval != NULL );

   switch( type )/*lint --e{788}*/
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
   case SCIP_SDPPAR_PENALTYPARAM:
      *dval = sdpisolver->penaltyparam;
      break;
   case SCIP_SDPPAR_OBJLIMIT:
      *dval = sdpisolver->objlimit;
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }

   return SCIP_OKAY;
}

/** sets floating point parameter of SDP-Solver */
SCIP_RETCODE SCIPsdpiSolverSetRealpar(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   SCIP_SDPPARAM         type,               /**< parameter number */
   SCIP_Real             dval                /**< parameter value */
   )
{
   assert( sdpisolver != NULL );

   switch( type )/*lint --e{788}*/
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
   case SCIP_SDPPAR_PENALTYPARAM:
      sdpisolver->penaltyparam = dval;
      SCIPdebugMessage("Setting sdpisolver penaltyparameter to %f.\n", dval);
      break;
   case SCIP_SDPPAR_OBJLIMIT:
      SCIPdebugMessage("Setting sdpisolver objlimit to %f.\n", dval);
      sdpisolver->objlimit = dval;
      break;
   case SCIP_SDPPAR_LAMBDASTAR:
      SCIPdebugMessage("Parameter SCIP_SDPPAR_LAMBDASTAR not used by DSDP"); /* this parameter is only used by SDPA */
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }

   return SCIP_OKAY;
}

/** gets integer parameter of SDP-Solver */
SCIP_RETCODE SCIPsdpiSolverGetIntpar(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   SCIP_SDPPARAM         type,               /**< parameter number */
   int*                  ival                /**< parameter value */
   )
{
   assert( sdpisolver != NULL );

   switch( type )/*lint --e{788}*/
   {
   case SCIP_SDPPAR_SDPINFO:
      *ival = (int) sdpisolver->sdpinfo;
      SCIPdebugMessage("Getting sdpisolver information output (%d).\n", *ival);
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }

   return SCIP_OKAY;
}

/** sets integer parameter of SDP-Solver */
SCIP_RETCODE SCIPsdpiSolverSetIntpar(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   SCIP_SDPPARAM         type,               /**< parameter number */
   int                   ival                /**< parameter value */
   )
{
   assert( sdpisolver != NULL );

   switch( type )/*lint --e{788}*/
   {
   case SCIP_SDPPAR_SDPINFO:
      sdpisolver->sdpinfo = (SCIP_Bool) ival;
      SCIPdebugMessage("Setting sdpisolver information output (%d).\n", ival);
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }

   return SCIP_OKAY;
}

/** compute and set lambdastar (only used for SDPA) */
SCIP_RETCODE SCIPsdpiSolverComputeLambdastar(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   SCIP_Real             maxguess            /**< maximum guess for lambda star of all SDP-constraints */
   )
{
   SCIPdebugMessage("Lambdastar parameter not used by DSDP"); /* this parameter is only used by SDPA */

   return SCIP_OKAY;
}/*lint !e715 */

/** compute and set the penalty parameter */
SCIP_RETCODE SCIPsdpiSolverComputePenaltyparam(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   SCIP_Real             maxcoeff,           /**< maximum objective coefficient */
   SCIP_Real*            penaltyparam        /**< the computed penalty parameter */
   )
{
   SCIP_Real compval;

   assert( sdpisolver != NULL );
   assert( penaltyparam != NULL );

   compval = PENALTYPARAM_FACTOR * maxcoeff;

   if ( compval < MIN_PENALTYPARAM )
   {
      SCIPdebugMessage("Setting penaltyparameter to %f.\n", MIN_PENALTYPARAM);
      sdpisolver->penaltyparam = MIN_PENALTYPARAM;
      *penaltyparam = MIN_PENALTYPARAM;
   }
   else if ( compval > MAX_PENALTYPARAM )
   {
      SCIPdebugMessage("Setting penaltyparameter to %f.\n", MAX_PENALTYPARAM);
      sdpisolver->penaltyparam = MAX_PENALTYPARAM;
      *penaltyparam = MAX_PENALTYPARAM;
   }
   else
   {
      SCIPdebugMessage("Setting penaltyparameter to %f.\n", compval);
      sdpisolver->penaltyparam = compval;
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
{
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

   /* if the maximum penalty parameter is smaller than the initial penalty paramater, we decrease the initial one correspondingly */
   if ( sdpisolver->penaltyparam > *maxpenaltyparam )
   {
      SCIPdebugMessage("Decreasing penaltyparameter of %f to maximum penalty paramater of %f.\n", sdpisolver->penaltyparam, *maxpenaltyparam);
      sdpisolver->penaltyparam = *maxpenaltyparam;
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
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   const char*           fname               /**< file name */
   )
{/*lint --e{715}*/
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_LPERROR;
}

/** writes SDP to a file */
SCIP_RETCODE SCIPsdpiSolverWriteSDP(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   const char*           fname               /**< file name */
   )
{/*lint --e{715}*/
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_LPERROR;
}

/**@} */
