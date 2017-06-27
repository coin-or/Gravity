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

/* #define SCIP_DEBUG*/
/* #define SCIP_MORE_DEBUG   *//* shows all added nonzero entries */
/* #define SCIP_DEBUG_PRINTTOFILE  *//* prints each problem inserted into SDPA to the file sdpa.dat-s and the starting point to sdpa.ini-s */

/* #define SDPA_RESETPARAMS */ /* this can be used together with an update to the SDPA source code to prevent memory leaks when using SCIP-SDP with SDPA */

/**@file   sdpisolver_sdpa.cpp
 * @brief  interface for SDPA
 * @author Tristan Gally
 * @author Ambros Gleixner
 */

/* turn off lint warnings for whole file: */
/*lint --e{750,788,818}*/

#include <assert.h>

#include "sdpi/sdpisolver.h"

/* turn off warnings for sdpa (doesn't seem to work) */
#pragma GCC diagnostic ignored "-Wstrict-prototypes"
#include "sdpa_call.h"                       /* SDPA callable library interface */
#pragma GCC diagnostic warning "-Wstrict-prototypes"

#include "blockmemshell/memory.h"            /* for memory allocation */
#include "scip/def.h"                        /* for SCIP_Real, _Bool, ... */
#include "scip/pub_misc.h"                   /* for sorting */
#include "sdpi/sdpsolchecker.h"              /* to check solution with regards to feasibility tolerance */

/* local defines */
#define GAPTOLCHANGE                1        /**< change gaptol by this factor when switching from fast to default and from default to stable settings */
#define FEASTOLCHANGE               1        /**< change feastol by this factor when switching from fast to default and from default to stable settings */
#define PENALTYBOUNDTOL             1E-3     /**< if the relative gap between Tr(X) and penaltyparam for a primal solution of the penaltyformulation
                                              *   is bigger than this value, it will be reported to the sdpi */

#define INFEASFEASTOLCHANGE         0.1      /**< change feastol by this factor if the solution was found to be infeasible with regards to feastol */
#define INFEASMINFEASTOL            1E-9     /**< minimum value for feasibility tolerance when encountering problems with regards to tolerance */

#define MIN_LAMBDASTAR              1e0      /**< if lambda star is to be computed, this is the minimum value it will take */
#define MAX_LAMBDASTAR              1e8      /**< if lambda star is to be computed, this is the maximum value it will take */
#define LAMBDASTAR_FACTOR           1e0      /**< if lambda star is to be computed, the biggest guess of the SDP blocks is multiplied by this value */
#define LAMBDASTAR_TWOPOINTS        TRUE     /**< if lambda star is to be computed, should we use only a low and a high value or instead a continuous interval */
#define LAMBDASTAR_THRESHOLD        1e1      /**< if lambda star is to be computed and LAMBDASTAR_TWOPOINTS=TRUE, then we distinguish between low and high using this */
#define LAMBDASTAR_LOW              1.5      /**< if lambda star is to be computed and LAMBDASTAR_TWOPOINTS=TRUE, then this is the value for below the threshold */
#define LAMBDASTAR_HIGH             1e5      /**< if lambda star is to be computed and LAMBDASTAR_TWOPOINTS=TRUE, then this is the value for above the threshold */

#define MIN_PENALTYPARAM            1e5      /**< if the penalty parameter is to be computed, this is the minimum value it will take */
#define MAX_PENALTYPARAM            1e12     /**< if the penalty parameter is to be computed, this is the maximum value it will take */
#define PENALTYPARAM_FACTOR         1e1      /**< if the penalty parameter is to be computed, the maximal objective coefficient will be multiplied by this */
#define MAX_MAXPENALTYPARAM         1e15     /**< if the maximum penaltyparameter is to be computed, this is the maximum value it will take */
#define MAXPENALTYPARAM_FACTOR      1e6      /**< if the maximum penaltyparameter is to be computed, it will be set to penaltyparam * this */

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
#define CHECK_IF_SOLVED_BOOL(sdpisolver)  do                                                                 \
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
   SDPA*                 sdpa;               /**< solver-object */
   int                   nvars;              /**< number of input variables */
   int                   nactivevars;        /**< number of variables present in SDPA (nvars minus the number of variables with lb = ub) */
   int*                  inputtosdpamapper;  /**< entry i gives the index of input variable i in sdpa (starting from 1) or
                                              *   -j (j=1, 2, ..., nvars - nactivevars) if the variable is fixed, the value and objective value of
                                              *   this fixed variable can be found in entry j-1 of fixedval/obj */
   int*                  sdpatoinputmapper;  /**< entry i gives the original index of the (i+1)-th variable in sdpa (indices go from 0 to nactivevars-1) */
   SCIP_Real*            fixedvarsval;       /**< entry i gives the lower and upper bound of the i-th fixed variable */
   SCIP_Real             fixedvarsobjcontr;  /**< total contribution to the objective of all fixed variables, computed as sum obj * val */
   int                   nvarbounds;         /**< number of variable bounds given to sdpa, length of sdpavarboundpos */
   int*                  varboundpos;        /**< maps position of variable bounds in the variable bound part of the LP-block in sdpa to the sdpa-indices
                                              *   of the corresponding variables, -n means lower bound of variable n, +n means upper bound */
   SCIP_Bool             solved;             /**< Was the SDP solved since the problem was last changed? */
   SCIP_Bool             timelimit;          /**< Was the SDP not given to the solver because the time limit was already reached? */
   int                   sdpcounter;         /**< used for debug messages */
   int                   niterations;        /**< number of SDP-iterations since the last solve call */
   int                   nsdpcalls;          /**< number of SDP-calls since the last solve call */
   SCIP_Real             epsilon;            /**< tolerance for absolute checks */
   SCIP_Real             gaptol;             /**< this is used for checking if primal and dual objective are equal */
   SCIP_Real             feastol;            /**< feasibility tolerance that should be achieved */
   SCIP_Real             sdpsolverfeastol;   /**< feasibility tolerance given to SDP-solver */
   SCIP_Real             objlimit;           /**< objective limit for SDP-solver */
   SCIP_Bool             sdpinfo;            /**< Should the SDP-solver output information to the screen? */
   SCIP_Bool             penalty;            /**< was the problem last solved using a penalty formulation */
   SCIP_Bool             rbound;             /**< was the penalty parameter bounded during the last solve call */
   SCIP_SDPSOLVERSETTING usedsetting;        /**< setting used to solve the last SDP */
   SCIP_Real             lambdastar;         /**< lambda star parameter to give to SDPA for initial point */
};


/*
 * Local Functions
 */

#ifndef NDEBUG
/** Test if a lower bound lb is not smaller than an upper bound ub, meaning that lb > ub - epsilon */
static
SCIP_Bool isFixed(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   SCIP_Real             lb,                 /**< lower bound */
   SCIP_Real             ub                  /**< upper bound */
   )
{
   assert( sdpisolver != NULL );
   assert( lb < ub + sdpisolver->feastol );

   return (ub-lb <= sdpisolver->epsilon);
}
#else
#define isFixed(sdpisolver,lb,ub) (ub-lb <= sdpisolver->epsilon)
#endif

/** If the problem is feasible for SDPA but not within our feasibility tolerance, adjust feasibility tolerance in
 *  SDPA and resolve until feasibility in SDPA and feasibility with regards to our tolerance match */
static
SCIP_RETCODE checkFeastolAndResolve(
   SCIP_SDPISOLVER*      sdpisolver,         /**< SDP-solver interface */
   SCIP_Real             penaltyparam,       /**< penalty parameter Gamma */
   int                   nvars,              /**< number of variables */
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
   SCIP_Real***          sdpval,             /**< values of SDP-constraintmmatrix entries (may be NULL if sdpnnonz = 0) */
   int**                 indchanges,         /**< changes needed to be done to the indices, if indchanges[block][ind]=-1, then the index can
                                              *   be removed, otherwise it gives the number of indices removed before this */
   int*                  nremovedinds,       /**< the number of rows/cols to be fixed for each block */
   int*                  blockindchanges,    /**< block indizes will be modified by these, see indchanges */
   int                   nlpcons,            /**< number of active (at least two nonzeros) LP-constraints */
   int                   noldlpcons,         /**< number of LP-constraints including those with less than two active nonzeros */
   SCIP_Real*            lplhs,              /**< left-hand sides of active LP-rows after fixings (may be NULL if nlpcons = 0) */
   SCIP_Real*            lprhs,              /**< right-hand sides of active LP-rows after fixings (may be NULL if nlpcons = 0) */
   int*                  rownactivevars,     /**< number of active variables for each LP-constraint */
   int                   lpnnonz,            /**< number of nonzero elements in the LP-constraint-matrix */
   int*                  lprow,              /**< row-index for each entry in lpval-array, might get sorted (may be NULL if lpnnonz = 0) */
   int*                  lpcol,              /**< column-index for each entry in lpval-array, might get sorted (may be NULL if lpnnonz = 0) */
   SCIP_Real*            lpval,              /**< values of LP-constraint-matrix entries, might get sorted (may be NULL if lpnnonz = 0) */
   SCIP_Real*            feastol             /**< current feasibility tolerance for SDPA */
   )
{
#ifdef SCIP_DEBUG
   char phase_string[15];
#endif

   assert( feastol != NULL );

   while ( SCIPsdpiSolverIsAcceptable(sdpisolver) && SCIPsdpiSolverIsDualFeasible(sdpisolver) && penaltyparam < sdpisolver->epsilon && *feastol >= INFEASMINFEASTOL )
   {
      SCIP_Real* solvector;
      int nvarspointer;
      SCIP_Bool infeasible;

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
         SCIPdebugMessage("Solution feasible for SDPA but outside feasibility tolerance, changing SDPA feasibility tolerance from %f to %f\n",
               *feastol, *feastol * INFEASFEASTOLCHANGE);
         *feastol *= INFEASFEASTOLCHANGE;

         if ( *feastol >= INFEASMINFEASTOL )
         {
            /* update settings */
            sdpisolver->sdpa->setParameterEpsilonDash(*feastol);

#ifdef SCIP_MORE_DEBUG
            sdpisolver->sdpa->printParameters(stdout);
#endif
            sdpisolver->sdpa->setInitPoint(false);
#ifdef SDPA_RESETPARAMS
            sdpisolver->sdpa->resetParameters();
#else
            sdpisolver->sdpa->initializeSolve();
#endif
            sdpisolver->sdpa->solve();

            /* update number of SDP-iterations and -calls */
            sdpisolver->niterations += (int) sdpisolver->sdpa->getIteration();
            sdpisolver->nsdpcalls += 1;

#ifdef SCIP_DEBUG
            /* print the phase value , i.e. whether solving was successfull */
            sdpisolver->sdpa->getPhaseString((char*)phase_string);
            SCIPdebugMessage("SDPA solving finished with status %s (primal and dual here are switched in contrast to our formulation)\n", phase_string);
#endif
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
   return SCIP_OKAY;
}

/*
 * Miscellaneous Methods
 */

/**@name Miscellaneous Methods */
/**@{ */


/** gets name and version (if available) of SDP-solver */
const char* SCIPsdpiSolverGetSolverName(
   void
   )
{/*lint !e1784*/
   return "SDPA"; /* version number can only be printed to file/stdout */
}

/** gets description of SDP-solver (developer, webpage, ...) */
const char* SCIPsdpiSolverGetSolverDesc(
   void
   )
{/*lint !e1784*/
   return "Primal-dual Interior Point Solver for SDPs developed by K. Fujisawa et al. (sdpa.sourceforge.net)";
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
{/*lint !e1784*/
   assert( sdpisolver != NULL );
   return (void*) sdpisolver->sdpa;
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
   return 2;
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
{/*lint !e1784*/
   assert( sdpisolver != NULL );
   assert( blkmem != NULL );
   assert( bufmem != NULL );

   SCIPdebugMessage("Calling SCIPsdpiCreate \n");

   BMS_CALL( BMSallocBlockMemory(blkmem, sdpisolver) );

   (*sdpisolver)->messagehdlr = messagehdlr;
   (*sdpisolver)->blkmem = blkmem;
   (*sdpisolver)->bufmem = bufmem;

   /* this will be properly initialized then calling solve */
   (*sdpisolver)->sdpa = NULL;

   (*sdpisolver)->nvars = 0;
   (*sdpisolver)->nactivevars = 0;
   (*sdpisolver)->inputtosdpamapper = NULL;
   (*sdpisolver)->sdpatoinputmapper = NULL;
   (*sdpisolver)->fixedvarsval = NULL;
   (*sdpisolver)->fixedvarsobjcontr = 0.0;
   (*sdpisolver)->nvarbounds = 0;
   (*sdpisolver)->varboundpos = NULL;
   (*sdpisolver)->solved = FALSE;
   (*sdpisolver)->timelimit = FALSE;
   (*sdpisolver)->sdpcounter = 0;
   (*sdpisolver)->niterations = 0;
   (*sdpisolver)->nsdpcalls = 0;

   (*sdpisolver)->epsilon = 1e-9;
   (*sdpisolver)->gaptol = 1e-4;
   (*sdpisolver)->feastol = 1e-6;
   (*sdpisolver)->sdpsolverfeastol = 1e-6;
   (*sdpisolver)->objlimit = SCIPsdpiSolverInfinity(*sdpisolver);
   (*sdpisolver)->sdpinfo = FALSE;
   (*sdpisolver)->usedsetting = SCIP_SDPSOLVERSETTING_UNSOLVED;

   return SCIP_OKAY;
}

/** deletes an SDP solver interface */
SCIP_RETCODE SCIPsdpiSolverFree(
   SCIP_SDPISOLVER**     sdpisolver          /**< pointer to an SDP-solver interface */
   )
{/*lint !e1784*/
   assert( sdpisolver != NULL );
   assert( *sdpisolver != NULL );

   SCIPdebugMessage("Freeing SDPISolver\n");

   /* free SDPA object using destructor and free memory via blockmemshell */
   if ( (*sdpisolver)->sdpa != NULL)
      delete (*sdpisolver)->sdpa;

   BMSfreeBlockMemoryArrayNull((*sdpisolver)->blkmem, &(*sdpisolver)->varboundpos, 2 * (*sdpisolver)->nactivevars); /*lint !e647*/

   BMSfreeBlockMemoryArrayNull((*sdpisolver)->blkmem, &(*sdpisolver)->inputtosdpamapper, (*sdpisolver)->nvars);/*lint !e737*/

   BMSfreeBlockMemoryArrayNull((*sdpisolver)->blkmem, &(*sdpisolver)->sdpatoinputmapper, (*sdpisolver)->nactivevars);/*lint !e737*/

   BMSfreeBlockMemoryArrayNull((*sdpisolver)->blkmem, &(*sdpisolver)->fixedvarsval, (*sdpisolver)->nvars - (*sdpisolver)->nactivevars); /*lint !e776*/

   BMSfreeBlockMemory((*sdpisolver)->blkmem, sdpisolver);

   return SCIP_OKAY;
}

/** increases the SDP-Counter */
SCIP_RETCODE SCIPsdpiSolverIncreaseCounter(
   SCIP_SDPISOLVER*      sdpisolver          /**< SDP-solver interface */
   )
{/*lint !e1784*/
   assert( sdpisolver != NULL );

   sdpisolver->sdpcounter++;

   return SCIP_OKAY;
}

/** reset the SDP-Counter to zero */
SCIP_RETCODE SCIPsdpiSolverResetCounter(
   SCIP_SDPISOLVER*      sdpisolver          /**< SDP-solver interface */
   )
{/*lint !e1784*/
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
   int*                  lprownactivevars,   /**< number of active variables for each LP-constraint */
   int                   lpnnonz,            /**< number of nonzero elements in the LP-constraint-matrix */
   int*                  lprow,              /**< row-index for each entry in lpval-array, might get sorted (may be NULL if lpnnonz = 0) */
   int*                  lpcol,              /**< column-index for each entry in lpval-array, might get sorted (may be NULL if lpnnonz = 0) */
   SCIP_Real*            lpval,              /**< values of LP-constraint-matrix entries, might get sorted (may be NULL if lpnnonz = 0) */
   SCIP_Real*            start,              /**< NULL or a starting point for the solver, this should have length nvars */
   SCIP_SDPSOLVERSETTING startsettings,      /**< settings used to start with in SDPA, currently not used for DSDP and MOSEK, set this to
                                              *   SCIP_SDPSOLVERSETTING_UNSOLVED to ignore it and start from scratch */
   SCIP_Real             timelimit           /**< after this many seconds solving will be aborted (currently only implemented for DSDP and MOSEK) */
   )
{/*lint !e1784*/
   return SCIPsdpiSolverLoadAndSolveWithPenalty(sdpisolver, 0.0, TRUE, FALSE, nvars, obj, lb, ub, nsdpblocks, sdpblocksizes, sdpnblockvars, sdpconstnnonz,
               sdpconstnblocknonz, sdpconstrow, sdpconstcol, sdpconstval, sdpnnonz, sdpnblockvarnonz, sdpvar, sdprow, sdpcol, sdpval, indchanges,
               nremovedinds, blockindchanges, nremovedblocks, nlpcons, noldlpcons, lplhs, lprhs, lprownactivevars, lpnnonz, lprow, lpcol, lpval, start,
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
                                              *   SCIP_SDPSOLVERSETTING_UNSOLVED to ignore it and start from scratch */
   SCIP_Real             timelimit,          /**< after this many seconds solving will be aborted (currently only implemented for DSDP and MOSEK) */
   SCIP_Bool*            feasorig,           /**< pointer to store if the solution to the penalty-formulation is feasible for the original problem
                                              *   (may be NULL if penaltyparam = 0) */
   SCIP_Bool*            penaltybound        /**< pointer to store if the primal solution reached the bound Tr(X) <= penaltyparam in the primal problem,
                                              *   this is also an indication of the penalty parameter being to small (may be NULL if not needed) */
   )
{/*lint !e1784*/
   SCIP_Real feastol;
   SCIP_Real* sdpavarbounds;
   int i;
   int k;
   int block;
   int nfixedvars;
   bool checkinput; /* should the input be checked for consistency in SDPA ? */
   int lpconsind;
   int lastrow;
   int* rowmapper; /* maps the lhs- and rhs-inequalities of the old LP-cons to their constraint numbers in DSDP */
   int nlpineqs;
   int pos;
   int newpos;
   int oldnactivevars;
#ifdef SCIP_MORE_DEBUG
   int ind;
#endif

#ifdef SCIP_DEBUG
   char phase_string[15];
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

   sdpisolver->niterations = 0;
   sdpisolver->nsdpcalls = 0;

   /* immediately exit if the time limit is negative */
   if ( timelimit <= 0.0 )
   {
      sdpisolver->solved = FALSE;
      sdpisolver->timelimit = TRUE;
      return SCIP_OKAY;
   }
   else
      sdpisolver->timelimit = FALSE;

#ifndef SCIP_NDEBUG
   checkinput = false;
#else
   checkinput = true;
#endif

   sdpisolver->usedsetting = SCIP_SDPSOLVERSETTING_UNSOLVED;

   /* only increase the counter if we don't use the penalty formulation to stay in line with the numbers in the general interface (where this is still the
    * same SDP) */
   if ( penaltyparam < sdpisolver->epsilon )
      SCIPdebugMessage("Inserting Data into SDPA for SDP (%d) \n", ++sdpisolver->sdpcounter);
   else
      SCIPdebugMessage("Inserting Data again into SDPA for SDP (%d) \n", sdpisolver->sdpcounter);

   /* set the penalty and rbound flags accordingly */
   sdpisolver->penalty = (penaltyparam < sdpisolver->epsilon) ? FALSE : TRUE;
   sdpisolver->rbound = rbound;

   /* allocate memory for inputtosdpamapper, sdpatoinputmapper and the fixed variable information, for the latter this will
    * later be shrinked if the needed size is known */
   BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &(sdpisolver->inputtosdpamapper), sdpisolver->nvars, nvars) );
   BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &(sdpisolver->sdpatoinputmapper), sdpisolver->nactivevars, nvars) );
   BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &(sdpisolver->fixedvarsval), sdpisolver->nvars - sdpisolver->nactivevars, nvars) ); /*lint !e776*/

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
         sdpisolver->fixedvarsval[nfixedvars] = lb[i]; /* if lb=ub, than this is the value the variable will have in every solution */
         nfixedvars++;
         sdpisolver->inputtosdpamapper[i] = -nfixedvars;
         SCIPdebugMessage("Fixing variable %d locally to %f for SDP %d in SDPA\n", i, lb[i], sdpisolver->sdpcounter);
      }
      else
      {
         sdpisolver->sdpatoinputmapper[sdpisolver->nactivevars] = i;
         sdpisolver->nactivevars++;
         sdpisolver->inputtosdpamapper[i] = sdpisolver->nactivevars; /* sdpa starts counting at 1, so we do this after increasing nactivevars */
#ifdef SCIP_MORE_DEBUG
         SCIPdebugMessage("Variable %d becomes variable %d for SDP %d in SDPA\n", i, sdpisolver->inputtosdpamapper[i], sdpisolver->sdpcounter);
#endif
      }
   }
   assert( sdpisolver->nactivevars + nfixedvars == sdpisolver->nvars );

   /* if we want to solve without objective, we reset fixedvarsobjcontr */
   if ( ! withobj )
      sdpisolver->fixedvarsobjcontr = 0.0;

   /* shrink the fixedvars and sdpatoinputmapper arrays to the right size */
   BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &(sdpisolver->fixedvarsval), nvars, nfixedvars) );
   BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &(sdpisolver->sdpatoinputmapper), nvars, sdpisolver->nactivevars) );

   /* insert data */
   if ( sdpisolver->sdpa != 0 )
   {
      /* if the SDPA solver has already been created, clear the current problem instance */
      delete sdpisolver->sdpa;
      sdpisolver->sdpa = new SDPA();
   }
   else
      sdpisolver->sdpa = new SDPA();
   assert( sdpisolver->sdpa != 0 );

   /* initialize settings (this needs to be done before inserting the problem as the initial point depends on the settings) */
   if ( penaltyparam >= sdpisolver->epsilon || startsettings == SCIP_SDPSOLVERSETTING_STABLE || startsettings == SCIP_SDPSOLVERSETTING_PENALTY )
   {
      sdpisolver->sdpa->setParameterType(SDPA::PARAMETER_STABLE_BUT_SLOW); /* if we already had problems with this problem, there is no reason to try fast */
      /* as we want to solve with stable settings, we also update epsilon and the feasibility tolerance, as we skip the default settings, we multpy twice */
      sdpisolver->sdpa->setParameterEpsilonStar(GAPTOLCHANGE * GAPTOLCHANGE * sdpisolver->gaptol);
      sdpisolver->sdpa->setParameterEpsilonDash(FEASTOLCHANGE * FEASTOLCHANGE * sdpisolver->sdpsolverfeastol);
      feastol = FEASTOLCHANGE * FEASTOLCHANGE * sdpisolver->sdpsolverfeastol;
      SCIPdebugMessage("Start solving process with stable settings\n");
   }
   else if ( startsettings == SCIP_SDPSOLVERSETTING_UNSOLVED || startsettings == SCIP_SDPSOLVERSETTING_FAST)
   {
      sdpisolver->sdpa->setParameterType(SDPA::PARAMETER_UNSTABLE_BUT_FAST);
      sdpisolver->sdpa->setParameterEpsilonStar(sdpisolver->gaptol);
      sdpisolver->sdpa->setParameterEpsilonDash(sdpisolver->sdpsolverfeastol);
      feastol = sdpisolver->sdpsolverfeastol;
      SCIPdebugMessage("Start solving process with fast settings\n");
   }
   else if ( startsettings == SCIP_SDPSOLVERSETTING_MEDIUM )
   {
      sdpisolver->sdpa->setParameterType(SDPA::PARAMETER_DEFAULT);
      /* as we want to solve with stable settings, we also update epsilon and the feasibility tolerance, as we skip the default settings, we multpy once */
      sdpisolver->sdpa->setParameterEpsilonStar(GAPTOLCHANGE * sdpisolver->gaptol);
      sdpisolver->sdpa->setParameterEpsilonDash(FEASTOLCHANGE * sdpisolver->sdpsolverfeastol);
      feastol = FEASTOLCHANGE * sdpisolver->sdpsolverfeastol;
      SCIPdebugMessage("Start solving process with medium settings\n");
   }
   else
   {
      SCIPdebugMessage("Unknown setting for start-settings: %d!\n", startsettings);  \
      return SCIP_LPERROR;
   }
   sdpisolver->sdpa->setParameterLowerBound(-1e20);


   /* set the objective limit */
   if ( ! SCIPsdpiSolverIsInfinity(sdpisolver, sdpisolver->objlimit) )
      sdpisolver->sdpa->setParameterUpperBound(sdpisolver->objlimit);
   else
      sdpisolver->sdpa->setParameterUpperBound(1e8);

   /* increase Lambda Star, this seems to help the numerics */
   sdpisolver->sdpa->setParameterLambdaStar(sdpisolver->lambdastar);

#ifdef SCIP_MORE_DEBUG
   sdpisolver->sdpa->printParameters(stdout);
#endif

   /* compute number of variable bounds and save them in sdpavarbounds */
   sdpisolver->nvarbounds = 0;
   BMS_CALL( BMSallocBufferMemoryArray(sdpisolver->bufmem, &sdpavarbounds, 2 * sdpisolver->nactivevars) ); /*lint !e647*//*lint !e530*/

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
      assert( 0 <= sdpisolver->sdpatoinputmapper[i] && sdpisolver->sdpatoinputmapper[i] < nvars );
      if ( ! SCIPsdpiSolverIsInfinity(sdpisolver, lb[sdpisolver->sdpatoinputmapper[i]]))
      {
         sdpavarbounds[sdpisolver->nvarbounds] = lb[sdpisolver->sdpatoinputmapper[i]];
         sdpisolver->varboundpos[sdpisolver->nvarbounds] = -(i + 1); /* negative sign means lower bound, i + 1 because sdpa starts counting from one */
         (sdpisolver->nvarbounds)++;
      }
      if ( ! SCIPsdpiSolverIsInfinity(sdpisolver, ub[sdpisolver->sdpatoinputmapper[i]]) )
      {
         sdpavarbounds[sdpisolver->nvarbounds] = ub[sdpisolver->sdpatoinputmapper[i]];
         sdpisolver->varboundpos[sdpisolver->nvarbounds] = +(i + 1); /* positive sign means upper bound, i + 1 because sdpa starts counting from one */
         (sdpisolver->nvarbounds)++;
      }
   }

   if ( nlpcons > 0 )
   {
      /* allocate memory to save which lpconstraints are mapped to which index, entry 2i corresponds to the left hand side of row i, 2i+1 to the rhs */
      BMS_CALL( BMSallocBufferMemoryArray(sdpisolver->bufmem, &rowmapper, 2*noldlpcons) ); /*lint !e647*//*lint !e530*/

      /* compute the number of LP constraints after splitting the ranged rows and compute the rowmapper */
      pos = 1; /* SDPA starts counting the LP-inequalities at one */
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
      nlpineqs = pos - 1; /* minus one because we started at one as SDPA wants them numbered one to nlpineqs */
      assert( nlpineqs <= 2*nlpcons ); /* *2 comes from left- and right-hand-sides */
   }
   else
      nlpineqs = 0;

   /* if we use a penalty formulation, we need the constraint r >= 0 */
   if ( penaltyparam >= sdpisolver->epsilon && rbound )
      sdpisolver->nvarbounds++;

   if ( sdpisolver->sdpinfo )
      sdpisolver->sdpa->setDisplay(stdout);
   else
      sdpisolver->sdpa->setDisplay(0);

#ifdef SCIP_MORE_DEBUG
   FILE* fpOut = fopen("output.tmp", "w");
   if ( ! fpOut )
      exit(-1);
   sdpisolver->sdpa->setResultFile(fpOut);
#endif

   /* initialize blockstruct */
   if ( penaltyparam < sdpisolver->epsilon ) /* we initialize this with an exact 0.0 in Solve without penalty */
      sdpisolver->sdpa->inputConstraintNumber((long long) sdpisolver->nactivevars);/*lint !e747*/
   else
      sdpisolver->sdpa->inputConstraintNumber((long long) sdpisolver->nactivevars + 1);/*lint !e747*/ /* the additional variable is r for the penalty formulation */

   /* if there are any lp-cons/variable-bounds, we get an extra block for those, lastrow - nshifts is the number of lp constraints added */
   sdpisolver->sdpa->inputBlockNumber((long long) ((nlpineqs + sdpisolver->nvarbounds > 0) ?
         nsdpblocks - nremovedblocks + 1 : nsdpblocks - nremovedblocks)); /*lint !e834 !e747 !e776*/

   /* block+1 because SDPA starts counting at 1 */
   for (block = 0; block < nsdpblocks; block++)
   {
      if ( blockindchanges[block] >= 0 )
      {
         SCIPdebugMessage("adding block %d to SDPA as block %d with size %d\n",
               block, block - blockindchanges[block] + 1, sdpblocksizes[block] - nremovedinds[block]); /*lint !e834*/
         sdpisolver->sdpa->inputBlockSize((long long) block - blockindchanges[block] + 1,/*lint !e747, !e834*/
               (long long) sdpblocksizes[block] - nremovedinds[block]); /*lint !e834, !e776, !e747*/
         sdpisolver->sdpa->inputBlockType((long long) block - blockindchanges[block] + 1, SDPA::SDP); /*lint !e834, !e776, !e747*/
      }
   }
   if ( nlpineqs + sdpisolver->nvarbounds > 0 )
   {
      /* the last block is the lp block, the size has a negative sign */
      sdpisolver->sdpa->inputBlockSize((long long) nsdpblocks - nremovedblocks + 1, (long long) -(nlpineqs + sdpisolver->nvarbounds)); /*lint !e834, !e776, !e747*/
      sdpisolver->sdpa->inputBlockType((long long) nsdpblocks - nremovedblocks + 1, SDPA::LP); /*lint !e834, !e776, !e747*/
      SCIPdebugMessage("adding LP block to SDPA as block %d with size %d\n", nsdpblocks - nremovedblocks + 1,/*lint !e834*/
            -(nlpineqs + sdpisolver->nvarbounds)); /*lint !e834*/
   }

   sdpisolver->sdpa->initializeUpperTriangleSpace();

   /* set objective values */
   for (i = 0; i < sdpisolver->nactivevars; i++)
   {
      if ( withobj )
      {
         /* insert objective value, SDPA counts from 1 to n instead of 0 to n-1 */
         sdpisolver->sdpa->inputCVec((long long) i + 1, obj[sdpisolver->sdpatoinputmapper[i]]);/*lint !e747*/
#ifdef SCIP_MORE_DEBUG
         SCIPdebugMessage("inserting objective %f for variable %d which became variable %d in SDPA\n", obj[sdpisolver->sdpatoinputmapper[i]],
               sdpisolver->sdpatoinputmapper[i], i+1);
#endif
      }
      if ( penaltyparam >= sdpisolver->epsilon )
         sdpisolver->sdpa->inputCVec((long long) sdpisolver->nactivevars + 1, penaltyparam);/*lint !e747*/ /* set the objective of the additional var to penaltyparam */
   }

   /* if we want to use a starting point we have to tell SDPA to allocate memory for it */
   if ( start != NULL )
      sdpisolver->sdpa->setInitPoint(true);
   else
      sdpisolver->sdpa->setInitPoint(false);

   /* start inserting the non-constant SDP-Constraint-Matrices */
   if ( sdpnnonz > 0 )
   {
      int v;
      int blockvar;

      assert( nsdpblocks > 0 );
      assert( sdpnblockvarnonz != NULL );
      assert( sdpnblockvars != NULL );
      assert( sdpcol != NULL );
      assert( sdprow != NULL );
      assert( sdpval != NULL );
      assert( sdpvar != NULL );
      assert( indchanges != NULL );
      assert( nremovedinds != NULL );

      for (block = 0; block < nsdpblocks; block++)
      {
         /* if the block has no entries, we skip it */
         if ( blockindchanges[block] == -1 )
            continue;
#ifdef SCIP_MORE_DEBUG
         SCIPdebugMessage("   -> building block %d, which becomes block %d in SDPA (%d)\n", block, block - blockindchanges[block] + 1,sdpisolver->sdpcounter);
#endif
         /* iterate over all variables in this block */
         for (blockvar = 0; blockvar < sdpnblockvars[block]; blockvar++)
         {
            v = sdpisolver->inputtosdpamapper[sdpvar[block][blockvar]];

#ifdef SCIP_MORE_DEBUG
            SCIPdebugMessage("      -> adding coefficient matrix for variable %d which becomes variable %d in SDPA (%d)\n",
               sdpvar[block][blockvar], v, sdpisolver->sdpcounter);
#endif

            /* check if the variable is active */
            if ( v > -1 )
            {
               for (k = 0; k < sdpnblockvarnonz[block][blockvar]; k++)
               {
                  /* rows and cols with active nonzeros should not be removed */
                  assert( indchanges[block][sdprow[block][blockvar][k]] > -1 && indchanges[block][sdpcol[block][blockvar][k]] > -1 );

                  assert( indchanges[block][sdprow[block][blockvar][k]] <= sdprow[block][blockvar][k]);
                  assert( indchanges[block][sdpcol[block][blockvar][k]] <= sdpcol[block][blockvar][k]);

                  assert( 0 <= sdprow[block][blockvar][k] && sdprow[block][blockvar][k] < sdpblocksizes[block] );
                  assert( 0 <= sdpcol[block][blockvar][k] && sdpcol[block][blockvar][k] < sdpblocksizes[block] );

                  /* rows and columns start with one in SDPA, so we have to add 1 to the indices */
#ifdef SCIP_MORE_DEBUG
                  SCIPdebugMessage("         -> adding nonzero %g at (%d,%d) (%d)\n",
                        sdpval[block][blockvar][k],
                        sdpcol[block][blockvar][k] - indchanges[block][sdpcol[block][blockvar][k]] + 1,
                        sdprow[block][blockvar][k] - indchanges[block][sdprow[block][blockvar][k]] + 1,
                        sdpisolver->sdpcounter);
#endif

                  sdpisolver->sdpa->inputElement((long long) v, (long long) block - blockindchanges[block] + 1, /*lint !e834, !e747*/
                        (long long) sdpcol[block][blockvar][k] - indchanges[block][sdpcol[block][blockvar][k]] + 1, /*lint !e834, !e747*/
                        (long long) sdprow[block][blockvar][k] - indchanges[block][sdprow[block][blockvar][k]] + 1, /*lint !e834, !e747*/
                        sdpval[block][blockvar][k], checkinput); /*lint !e776 !e834*/
               }
            }
         }
         /* insert the identity matrix if we are using a penalty formulation */
         if ( penaltyparam >= sdpisolver->epsilon )
         {
#ifdef SCIP_MORE_DEBUG
            SCIPdebugMessage("      -> adding coefficient matrix for penalty variable r in SDPA (%d)\n", sdpisolver->sdpcounter);
#endif
            for (i = 0; i < sdpblocksizes[block] - nremovedinds[block]; i++)
            {
#ifdef SCIP_MORE_DEBUG
                  SCIPdebugMessage("         -> adding nonzero 1.0 at (%d,%d) (%d)\n", i + 1, i + 1, sdpisolver->sdpcounter);
#endif

                  sdpisolver->sdpa->inputElement((long long) sdpisolver->nactivevars + 1, (long long) block - blockindchanges[block] + 1, /*lint !e747, !e776, !e834*/
                        (long long) i + 1, (long long) i + 1, 1.0, checkinput);/*lint !e747*/
            }
         }
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

      for (block = 0; block < nsdpblocks; block++)
      {
         if ( blockindchanges[block] == -1 )
            continue;
#ifdef SCIP_MORE_DEBUG
         SCIPdebugMessage("   -> building block %d (%d)\n", block + 1, sdpisolver->sdpcounter);
#endif
         for (k = 0; k < sdpconstnblocknonz[block]; k++)
         {
            /* rows and cols with active nonzeros should not be removed */
            assert( indchanges[block][sdpconstrow[block][k]] > -1 && indchanges[block][sdpconstcol[block][k]] > -1 );

            assert( indchanges[block][sdpconstrow[block][k]] <= sdpconstrow[block][k]);
            assert( indchanges[block][sdpconstcol[block][k]] <= sdpconstcol[block][k]);

            assert (0 <= sdpconstrow[block][k] && sdpconstrow[block][k] < sdpblocksizes[block]);
            assert (0 <= sdpconstcol[block][k] && sdpconstcol[block][k] < sdpblocksizes[block]);

            /* rows and columns start with one in SDPA, so we have to add 1 to the indices, the constant matrix is given as variable 0 */
#ifdef SCIP_MORE_DEBUG
            SCIPdebugMessage("         -> adding constant nonzero %g at (%d,%d) (%d)\n", sdpconstval[block][k],
               sdpconstcol[block][k] - indchanges[block][sdpconstcol[block][k]] + 1,
               sdpconstrow[block][k] - indchanges[block][sdpconstrow[block][k]] + 1,
               sdpisolver->sdpcounter);
#endif
            sdpisolver->sdpa->inputElement((long long) 0, (long long) block - blockindchanges[block] + 1, /*lint !e747, !e776, !e834*/
               (long long) sdpconstcol[block][k] - indchanges[block][sdpconstcol[block][k]] + 1, /*lint !e747, !e776, !e834*/
               (long long) sdpconstrow[block][k] - indchanges[block][sdpconstrow[block][k]] + 1, /*lint !e747, !e776, !e834*/
               sdpconstval[block][k], checkinput);
         }
      }
   }

#ifdef SCIP_MORE_DEBUG
         SCIPdebugMessage("   -> building LP-block %d (%d)\n", nsdpblocks - nremovedblocks + 1, sdpisolver->sdpcounter);
#endif
   /* inserting LP nonzeros */
   lastrow = -1;
   for (i = 0; i < lpnnonz; i++)
   {
      assert( 0 <= lprow[i] && lprow[i] < noldlpcons );
      assert( 0 <= lpcol[i] && lpcol[i] < nvars );
      assert( REALABS(lpval[i]) > sdpisolver->epsilon );

      /* if the variable is active and the constraint is more than a bound, we add it */
      if ( sdpisolver->inputtosdpamapper[lpcol[i]] > 0 )
      {
       /* as this is an active variable, there should be at least one in the constraint */
         assert( rownactivevars[lprow[i]] > 0 );
         if ( rownactivevars[lprow[i]] > 1 )
         {
            if ( lprow[i] > lastrow )  /* we update the lpcons-counter */
            {
               lastrow = lprow[i];
               /* if we use a penalty formulation, add the r * Identity entry */
               if ( penaltyparam >= sdpisolver->epsilon )
               {
                  /* check for the lhs-inequality */
                  if ( rowmapper[2*lastrow] > -1 ) /*lint !e679*/
                  {
#ifdef SCIP_MORE_DEBUG
                     SCIPdebugMessage("         -> adding nonzero 1.0 at (%d,%d) for penalty variable r in SDPA (%d)\n",
                        rowmapper[2*lastrow], rowmapper[2*lastrow], sdpisolver->sdpcounter);
#endif
                     /* LP nonzeros are added as diagonal entries of the last block (coming after the last SDP-block, with
                      * blocks starting at 1, as are rows), the r-variable is variable nactivevars + 1 */
                     sdpisolver->sdpa->inputElement((long long) sdpisolver->nactivevars + 1, (long long) nsdpblocks - nremovedblocks + 1, /*lint !e747, !e776, !e834*/
                           (long long) rowmapper[2*lastrow], (long long) rowmapper[2*lastrow], 1.0, checkinput); /*lint !e679, !e747, !e834*/
                  }

                  /* check for the rhs-inequality */
                  if ( rowmapper[2*lastrow + 1] > -1 ) /*lint !e679*/
                  {
#ifdef SCIP_MORE_DEBUG
                     SCIPdebugMessage("         -> adding nonzero 1.0 at (%d,%d) for penalty variable r in SDPA (%d)\n",
                        rowmapper[2*lastrow + 1], rowmapper[2*lastrow + 1], sdpisolver->sdpcounter);
#endif
                     /* LP nonzeros are added as diagonal entries of the last block (coming after the last SDP-block, with
                      * blocks starting at 1, as are rows), the r-variable is variable nactivevars + 1 */
                     sdpisolver->sdpa->inputElement((long long) sdpisolver->nactivevars + 1, (long long) nsdpblocks - nremovedblocks + 1, /*lint !e747, !e776, !e834*/
                           (long long) rowmapper[2*lastrow + 1], (long long) rowmapper[2*lastrow + 1], 1.0, checkinput); /*lint !e679, !e747, !e834*/
                  }
               }
            }
            /* add the lp-nonzero to the lhs-inequality if it exists: */
            if ( rowmapper[2*lastrow] > -1 ) /*lint !e679*/
            {
#ifdef SCIP_MORE_DEBUG
               SCIPdebugMessage("         -> adding nonzero %g at (%d,%d) for variable %d which became variable %d in SDPA (%d)\n",
                  lpval[i], rowmapper[2*lastrow], rowmapper[2*lastrow], lpcol[i], sdpisolver->inputtosdpamapper[lpcol[i]], sdpisolver->sdpcounter);
#endif
               /* LP nonzeros are added as diagonal entries of the last block (coming after the last SDP-block, with blocks starting at 1, as are rows) */
               sdpisolver->sdpa->inputElement((long long) sdpisolver->inputtosdpamapper[lpcol[i]], (long long) nsdpblocks - nremovedblocks + 1, /*lint !e747, !e776, !e834*/
                     (long long) rowmapper[2*lastrow], (long long) rowmapper[2*lastrow], lpval[i], checkinput); /*lint !e679, !e747, !e834*/
            }
            /* add the lp-nonzero to the rhs-inequality if it exists: */
            if ( rowmapper[2*lastrow + 1] > -1 ) /*lint !e679*/
            {
#ifdef SCIP_MORE_DEBUG
               SCIPdebugMessage("         -> adding nonzero %g at (%d,%d) for variable %d which became variable %d in SDPA (%d)\n",
                  -1 * lpval[i], rowmapper[2*lastrow + 1], rowmapper[2*lastrow + 1], lpcol[i], sdpisolver->inputtosdpamapper[lpcol[i]], sdpisolver->sdpcounter);
#endif
               /* LP nonzeros are added as diagonal entries of the last block (coming after the last SDP-block, with blocks starting at 1, as are rows),
                * the -1 comes from the fact that this is a <=-constraint, while SDPA works with >= */
               sdpisolver->sdpa->inputElement((long long) sdpisolver->inputtosdpamapper[lpcol[i]], (long long) nsdpblocks - nremovedblocks + 1, /*lint !e747, !e776, !e834*/
                     (long long) rowmapper[2*lastrow + 1], (long long) rowmapper[2*lastrow + 1], -1 * lpval[i], checkinput); /*lint !e679, !e747, !e834*/
            }
         }
      }
   }

   /* inserting LP left- and right-hand-sides for active constraints */
   lpconsind = 1; /* this is the same order we used when computing the rowmapper, so we insert at the right positions */
   for (i = 0; i < nlpcons; i++)
   {
      /* check for left-hand side */
      if ( lplhs[i] > - SCIPsdpiSolverInfinity(sdpisolver) )
      {
#ifdef SCIP_MORE_DEBUG
         SCIPdebugMessage("         -> adding lhs %g at (%d,%d) (%d)\n", lplhs[i], lpconsind, lpconsind, sdpisolver->sdpcounter);
#endif
         /* LP constraints are added as diagonal entries of the last block, left-hand-side is added as variable zero */
         sdpisolver->sdpa->inputElement((long long) 0, (long long) nsdpblocks - nremovedblocks + 1, (long long) lpconsind, /*lint !e747, !e776, !e834*/
               (long long) lpconsind, lplhs[i], checkinput);/*lint !e747*/
         lpconsind++;
      }

      /* check for right-hand side */
      if ( lprhs[i] < SCIPsdpiSolverInfinity(sdpisolver) )
      {
#ifdef SCIP_MORE_DEBUG
         SCIPdebugMessage("         -> adding lhs (originally rhs) %g at (%d,%d) (%d)\n", -1 * lprhs[i], lpconsind, lpconsind, sdpisolver->sdpcounter);
#endif
         /* LP constraints are added as diagonal entries of the last block, right-hand side is added as variable zero, the -1 comes from
          * the fact, that SDPA uses >=, while the rhs is a <=-constraint */
         sdpisolver->sdpa->inputElement((long long) 0, (long long) nsdpblocks - nremovedblocks + 1, (long long) lpconsind, /*lint !e747, !e776, !e834*/
               (long long) lpconsind, -1 * lprhs[i], checkinput);/*lint !e747*/
         lpconsind++;
      }
   }

   assert( lpconsind == nlpineqs + 1 ); /* plus one because we started at one as SDPA wants them numbered one to nlpineqs */

   /* print each LP-constraint as one formatted constraint in addition to the single entries inserted into SDPA */
#ifdef SCIP_MORE_DEBUG
   lastrow = -1;
   ind = -1; /* this is increased once before the first usage */
   for (i = 0; i < lpnnonz; i++)
   {
      /* if the variable is active and the constraint is more than a bound, we added it */
      if ( sdpisolver->inputtosdpamapper[lpcol[i]] > 0 )
      {
         if ( rownactivevars[lprow[i]] > 1 )
         {
            if ( lprow[i] > lastrow )  /* we finished the old row */
            {
               if ( lastrow >= 0 )
               {
                  if ( lprhs[ind] < SCIPsdpiSolverInfinity(sdpisolver) )
                     SCIPmessagePrintInfo(sdpisolver->messagehdlr, " <= %f\n", lprhs[ind]);
                  else
                     SCIPmessagePrintInfo(sdpisolver->messagehdlr, "\n");
               }
               lastrow = lprow[i];
               ind++;
               if ( lplhs[ind] > - SCIPsdpiSolverInfinity(sdpisolver) )
                  SCIPmessagePrintInfo(sdpisolver->messagehdlr, "%f <= ", lplhs[ind]);
            }
            SCIPmessagePrintInfo(sdpisolver->messagehdlr, "+ %f <x%d> ", lpval[i], lpcol[i]);
         }
      }
   }
   if ( lastrow >= 0 )
   {
      if ( lprhs[ind] < SCIPsdpiSolverInfinity(sdpisolver) )
         SCIPmessagePrintInfo(sdpisolver->messagehdlr, " <= %f\n", lprhs[ind]);
      else
         SCIPmessagePrintInfo(sdpisolver->messagehdlr, "\n");
   }
   assert( ind == nlpcons - 1 ); /* -1 because we start indexing at zero and do not increase after the last row */
#endif

   /* free the memory for the rowmapper */
   if ( nlpcons > 0 )
      BMSfreeBufferMemoryArray(sdpisolver->bufmem, &rowmapper);/*lint !e647, !e737*/

   /* insert variable bounds, these are also added as LP-constraints and therefore diagonal entries of the LP block
    * if we work with the penalty formulation, we get an extra entry for r >= 0, but this we will add afterwards */
   for (i = 0; i < ((penaltyparam < sdpisolver->epsilon) || (! rbound) ? sdpisolver->nvarbounds : sdpisolver->nvarbounds - 1); i++)
   {
      assert( 0 < abs(sdpisolver->varboundpos[i]) && abs(sdpisolver->varboundpos[i] <= sdpisolver->nactivevars) ); /* the indices are already those for SDPA */

      /* for lower bound */
      if ( sdpisolver->varboundpos[i] < 0 )
      {
         /* add it as an lp-constraint for this variable (- because we saved -n for the lower bound), at the position
          * (nactivelpcons + 1) + varbound-index, because we have >= the variable has coefficient +1 */
         sdpisolver->sdpa->inputElement((long long) -sdpisolver->varboundpos[i], (long long) nsdpblocks - nremovedblocks + 1, /*lint !e747, !e776, !e834*/
               (long long) nlpineqs + 1 + i, (long long) nlpineqs + 1 + i, 1.0, checkinput);/*lint !e747*/

         if ( REALABS(sdpavarbounds[i]) > sdpisolver->epsilon )
         {
            /* the bound is added as the rhs and therefore variable zero */
#ifdef SCIP_MORE_DEBUG
            SCIPdebugMessage("         -> adding lower bound %g at (%d,%d) for variable %d which became variable %d in SDPA (%d)\n",
                  sdpavarbounds[i], nlpineqs + 1 + i, nlpineqs + 1 + i, sdpisolver->sdpatoinputmapper[-sdpisolver->varboundpos[i] - 1],
                  -sdpisolver->varboundpos[i], sdpisolver->sdpcounter);
#endif
            sdpisolver->sdpa->inputElement((long long) 0, (long long) nsdpblocks - nremovedblocks + 1, (long long) nlpineqs + 1 + i, /*lint !e747, !e776, !e834*/
                  (long long) nlpineqs + 1 + i, sdpavarbounds[i], checkinput);/*lint !e747*/
         }
         else
         {
            /* as the bound is zero, we don't need to add a right hand side */
#ifdef SCIP_MORE_DEBUG
            SCIPdebugMessage("         -> adding lower bound 0 at (%d,%d) for variable %d which became variable %d in SDPA (%d)\n",
                  nlpineqs + 1 + i, nlpineqs + 1 + i, sdpisolver->sdpatoinputmapper[-sdpisolver->varboundpos[i] - 1],
               -sdpisolver->varboundpos[i], sdpisolver->sdpcounter);
#endif
         }
      }
      else
      {
         /* this is an upper bound */

         /* add it as an lp-constraint for this variable, at the position nactivelpcons + varbound-index, because we have >= but we
          * want <= for the upper bound, we have to multiply by -1 and therefore the variable has coefficient -1 */
         sdpisolver->sdpa->inputElement((long long) sdpisolver->varboundpos[i], (long long) nsdpblocks - nremovedblocks + 1, /*lint !e747, !e776, !e834*/
               (long long) nlpineqs + 1 + i, (long long) nlpineqs + 1 + i, -1.0, checkinput);/*lint !e747*/

         if ( REALABS(sdpavarbounds[i]) > sdpisolver->epsilon )
         {
            /* the bound is added as the rhs and therefore variable zero, we multiply by -1 for <= */
#ifdef SCIP_MORE_DEBUG
            SCIPdebugMessage("         -> adding upper bound %g at (%d,%d) for variable %d which became variable %d in SDPA (%d)\n",
                  sdpavarbounds[i], nlpineqs + 1 + i, nlpineqs + 1 + i, sdpisolver->sdpatoinputmapper[sdpisolver->varboundpos[i] - 1],
                  sdpisolver->varboundpos[i], sdpisolver->sdpcounter);
#endif
            sdpisolver->sdpa->inputElement((long long) 0, (long long) nsdpblocks - nremovedblocks + 1, (long long) nlpineqs + 1 + i, /*lint !e747, !e776, !e834*/
                  (long long) nlpineqs + 1 + i, -sdpavarbounds[i], checkinput);/*lint !e747*/
         }
         else
         {
            /* as the bound is zero, we don't need to add a right hand side */
#ifdef SCIP_MORE_DEBUG
            SCIPdebugMessage("         -> adding upper bound 0 at (%d,%d) for variable %d which became variable %d in SDPA (%d)\n",
                  0, nlpineqs + 1 + i, nlpineqs + 1 + i, sdpisolver->sdpatoinputmapper[sdpisolver->varboundpos[i] - 1],
                  sdpisolver->varboundpos[i]);
#endif
         }
      }
   }

   if ( penaltyparam >= sdpisolver->epsilon && rbound )
   {
      /* we add the variable bound r >= 0 */
      sdpisolver->sdpa->inputElement((long long) sdpisolver->nactivevars + 1, (long long) nsdpblocks - nremovedblocks + 1, /*lint !e747, !e776, !e834*/
            (long long) nlpineqs + 1 + i, (long long) nlpineqs + 1 + i, 1.0, checkinput);/*lint !e747*/
#ifdef SCIP_MORE_DEBUG
      SCIPdebugMessage("         -> adding lower bound r >= 0 at (%d,%d)  in SDPA (%d)\n", nlpineqs + 1 + i, nlpineqs + 1 + i, sdpisolver->sdpcounter);
#endif
   }

   /* free the arrays used for counting and saving variable bounds and LP-right-hand-sides */
   BMSfreeBufferMemoryArray(sdpisolver->bufmem, &sdpavarbounds); /*lint !e647, !e737*/

   /* transform the matrices to a more efficient form */
   sdpisolver->sdpa->initializeUpperTriangle();
   sdpisolver->sdpa->initializeSolve();

   /* set the starting solution */
   if ( start != NULL && penaltyparam < sdpisolver->epsilon )
   {
      SCIPdebugMessage("Starting with a previous solution is not yet tested for the interface, only x-vector is given, not y and Z");
      for (i = 1; i <= sdpisolver->nactivevars; i++) /* we iterate over the variables in sdpa */
         sdpisolver->sdpa->inputInitXVec((long long) i, start[sdpisolver->sdpatoinputmapper[i] - 1]);/*lint !e747*/
   }
   else if ( penaltyparam >= sdpisolver->epsilon )
      SCIPdebugMessage("Skipping insertion of starting point, as this is not yet supported for penalty formulation.\n");

#ifdef SCIP_DEBUG_PRINTTOFILE
   /* if necessary, dump input data and initial point */
   sdpisolver->sdpa->writeInputSparse(const_cast<char*>("sdpa.dat-s"), const_cast<char*>("%+8.3e"));
   sdpisolver->sdpa->writeInitSparse(const_cast<char*>("sdpa.ini-s"), const_cast<char*>("%+8.3e"));
#endif

   SCIPdebugMessage("Calling SDPA solve (SDP: %d)\n", sdpisolver->sdpcounter);
   sdpisolver->sdpa->solve();
   sdpisolver->solved = TRUE;

   /* update number of SDP-iterations and -calls */
   sdpisolver->niterations += (int) sdpisolver->sdpa->getIteration();
   sdpisolver->nsdpcalls += 1;

#ifdef SCIP_DEBUG
   /* print the phase value , i.e. whether solving was successfull */
   sdpisolver->sdpa->getPhaseString((char*)phase_string);
   SCIPdebugMessage("SDPA solving finished with status %s (primal and dual here are switched in contrast to our formulation)\n", phase_string);
#endif

   /* remember settings */
   if ( SCIPsdpiSolverIsAcceptable(sdpisolver) && penaltyparam < sdpisolver->epsilon )
      sdpisolver->usedsetting = SCIP_SDPSOLVERSETTING_FAST;
   else if ( penaltyparam >= sdpisolver->epsilon )
      sdpisolver->usedsetting = SCIP_SDPSOLVERSETTING_PENALTY;

   /* if the problem has been stably solved but did not reach the required feasibility tolerance, even though the solver
    * reports feasibility, resolve it with adjusted tolerance */
   SCIP_CALL( checkFeastolAndResolve(sdpisolver, penaltyparam, nvars, lb, ub, nsdpblocks, sdpblocksizes, sdpnblockvars, sdpconstnnonz,
         sdpconstnblocknonz, sdpconstrow, sdpconstcol, sdpconstval, sdpnnonz, sdpnblockvarnonz, sdpvar, sdprow, sdpcol, sdpval, indchanges,
         nremovedinds, blockindchanges, nlpcons, noldlpcons, lplhs, lprhs, rownactivevars, lpnnonz, lprow, lpcol, lpval, &feastol) );

   /* check whether problem has been stably solved, if it wasn't and we didn't yet use the default parametersettings (for the penalty formulation we do so), try
    * again with more stable parameters */
   if ( ! SCIPsdpiSolverIsAcceptable(sdpisolver) && penaltyparam < sdpisolver->epsilon &&
      (startsettings == SCIP_SDPSOLVERSETTING_UNSOLVED || startsettings == SCIP_SDPSOLVERSETTING_FAST) )
   {
      SCIPdebugMessage("Numerical troubles -- solving SDP %d again ...\n", sdpisolver->sdpcounter);

      /* initialize settings */
      sdpisolver->sdpa->setParameterType(SDPA::PARAMETER_DEFAULT);
      sdpisolver->sdpa->setParameterEpsilonStar(GAPTOLCHANGE * sdpisolver->gaptol);
      sdpisolver->sdpa->setParameterEpsilonDash(FEASTOLCHANGE * sdpisolver->sdpsolverfeastol);
      sdpisolver->sdpa->setParameterLowerBound(-1e20);
      /* set the objective limit */
      if ( ! SCIPsdpiSolverIsInfinity(sdpisolver, sdpisolver->objlimit) )
         sdpisolver->sdpa->setParameterUpperBound(sdpisolver->objlimit);
      else
         sdpisolver->sdpa->setParameterUpperBound(1e8);

      /* increase Lambda Star, this seems to help the numerics */
      sdpisolver->sdpa->setParameterLambdaStar(sdpisolver->lambdastar);

#ifdef SCIP_MORE_DEBUG
   sdpisolver->sdpa->printParameters(stdout);
#endif
      sdpisolver->sdpa->setInitPoint(false);
#ifdef SDPA_RESETPARAMS
      sdpisolver->sdpa->resetParameters();
#else
      sdpisolver->sdpa->initializeSolve();
#endif
      sdpisolver->sdpa->solve();
      sdpisolver->solved = TRUE;

      /* update number of SDP-iterations and -calls */
      sdpisolver->niterations += (int) sdpisolver->sdpa->getIteration();
      sdpisolver->nsdpcalls += 1;

      /* remember setting */
      if ( SCIPsdpiSolverIsAcceptable(sdpisolver) )
         sdpisolver->usedsetting = SCIP_SDPSOLVERSETTING_MEDIUM;

#ifdef SCIP_DEBUG
      /* print the phase value , i.e. whether solving was successfull */
      sdpisolver->sdpa->getPhaseString((char*)phase_string);
      SCIPdebugMessage("SDPA solving finished with status %s (primal and dual here are switched in contrast to our formulation)\n", phase_string);
#endif

      /* if the problem has been stably solved but did not reach the required feasibility tolerance, even though the solver
       * reports feasibility, resolve it with adjusted tolerance */
      SCIP_CALL( checkFeastolAndResolve(sdpisolver, penaltyparam, nvars, lb, ub, nsdpblocks, sdpblocksizes, sdpnblockvars, sdpconstnnonz,
            sdpconstnblocknonz, sdpconstrow, sdpconstcol, sdpconstval, sdpnnonz, sdpnblockvarnonz, sdpvar, sdprow, sdpcol, sdpval, indchanges,
            nremovedinds, blockindchanges, nlpcons, noldlpcons, lplhs, lprhs, rownactivevars, lpnnonz, lprow, lpcol, lpval, &feastol) );
   }

   /* if we still didn't converge, and did not yet use the stable settings, set the parameters even more conservativly */
   if ( (! SCIPsdpiSolverIsAcceptable(sdpisolver)) && penaltyparam < sdpisolver->epsilon &&
      (startsettings == SCIP_SDPSOLVERSETTING_UNSOLVED || startsettings == SCIP_SDPSOLVERSETTING_FAST || startsettings == SCIP_SDPSOLVERSETTING_MEDIUM) )
   {
      SCIPdebugMessage("Numerical troubles -- solving SDP %d again^2 ...\n", sdpisolver->sdpcounter);

      /* initialize settings */
      sdpisolver->sdpa->setParameterType(SDPA::PARAMETER_STABLE_BUT_SLOW);
      sdpisolver->sdpa->setParameterEpsilonStar(GAPTOLCHANGE * GAPTOLCHANGE * sdpisolver->gaptol);
      sdpisolver->sdpa->setParameterEpsilonDash(FEASTOLCHANGE * FEASTOLCHANGE * sdpisolver->sdpsolverfeastol);
      sdpisolver->sdpa->setParameterLowerBound(-1e20);
      /* set the objective limit */
      if ( ! SCIPsdpiSolverIsInfinity(sdpisolver, sdpisolver->objlimit) )
         sdpisolver->sdpa->setParameterUpperBound(sdpisolver->objlimit);
      else
         sdpisolver->sdpa->setParameterUpperBound(1e8);

      /* increase Lambda Star, this seems to help the numerics */
      sdpisolver->sdpa->setParameterLambdaStar(sdpisolver->lambdastar);

#ifdef SCIP_MORE_DEBUG
sdpisolver->sdpa->printParameters(stdout);
#endif
      sdpisolver->sdpa->setInitPoint(false);
#ifdef SDPA_RESETPARAMS
      sdpisolver->sdpa->resetParameters();
#else
      sdpisolver->sdpa->initializeSolve();
#endif
      sdpisolver->sdpa->solve();
      sdpisolver->solved = TRUE;

      /* update number of SDP-iterations and -calls */
      sdpisolver->niterations += (int) sdpisolver->sdpa->getIteration();
      sdpisolver->nsdpcalls += 1;

      /* remember setting */
      if ( SCIPsdpiSolverIsAcceptable(sdpisolver) )
         sdpisolver->usedsetting = SCIP_SDPSOLVERSETTING_STABLE;

#ifdef SCIP_DEBUG
      /* print the phase value , i.e. whether solving was successfull */
      sdpisolver->sdpa->getPhaseString((char*)phase_string);
      SCIPdebugMessage("SDPA solving finished with status %s (primal and dual here are switched in constrast to our formulation)\n", phase_string);
#endif

      /* if the problem has been stably solved but did not reach the required feasibility tolerance, even though the solver
       * reports feasibility, resolve it with adjusted tolerance */
      SCIP_CALL( checkFeastolAndResolve(sdpisolver, penaltyparam, nvars, lb, ub, nsdpblocks, sdpblocksizes, sdpnblockvars, sdpconstnnonz,
            sdpconstnblocknonz, sdpconstrow, sdpconstcol, sdpconstval, sdpnnonz, sdpnblockvarnonz, sdpvar, sdprow, sdpcol, sdpval, indchanges,
            nremovedinds, blockindchanges, nlpcons, noldlpcons, lplhs, lprhs, rownactivevars, lpnnonz, lprow, lpcol, lpval, &feastol) );
   }

#ifdef SCIP_MORE_DEBUG
   (void) fclose(fpOut);
#endif

   /* if we solved a penalty formulation, check if the solution is feasible for the original problem (which is the case iff r < feastol) */
   if ( penaltyparam >= sdpisolver->epsilon )
   {
      SCIP_Real* sdpasol;
      SCIP_Real* X;
      int b;
      int nblockssdpa;
      int nrow;
      SCIP_Real trace = 0.0;

      /* in the second case we have r as an additional variable */
      assert( (sdpisolver->nactivevars + 1 == sdpisolver->sdpa->getConstraintNumber()) );  /*lint !e776*/

      sdpasol = sdpisolver->sdpa->getResultXVec();

      /* we get r as the last variable in SDPA */
      *feasorig = (sdpasol[sdpisolver->nactivevars] < sdpisolver->feastol); /*lint !e413*/

      /* if r > 0 or we are in debug mode, also check the primal bound */
#ifdef NDEBUG
      if ( ! *feasorig && penaltybound != NULL )
      {
#endif

         SCIPdebugMessage("Solution not feasible in original problem, r = %f\n", sdpasol[sdpisolver->nactivevars]);

         /* compute Tr(X) */

         /* iterate over all blocks (SDPA starts counting at one and includes the LP block) */
         nblockssdpa = (int) sdpisolver->sdpa->getBlockNumber();
         for (b = 1; b <= nblockssdpa; b++)
         {
            /* get the block from SDPA */
            X = sdpisolver->sdpa->getResultYMat((long long) b);/*lint !e747*/
            nrow = (int) sdpisolver->sdpa->getBlockSize((long long) b);/*lint !e747*/
            assert( nrow >= 0 );

            /* if it is the LP-block, we omit the variable bounds as the penalty variable is not added to them */
            if ( sdpisolver->sdpa->getBlockType((long long) b) == SDPA::LP )/*lint !e747*/
            {
               /* iterate over all diagonal entries (until we reach the varbound part), adding them to the trace */
               for (i = 0; i < nrow - sdpisolver->nvarbounds; i++)
                  trace += X[i]; /* get entry (i+1,i+1) for the diagonal matrix X */
            }
            else
            {
               /* iterate over all diagonal entries and add them to the trace */
               for (i = 0; i < nrow; i++)
                  trace += X[i + i*nrow];/*lint !e679*/ /* get entry (i+1,i+1) in X */
            }
         }

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

#ifdef NDEBUG
      }
#endif
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
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP-solver interface */
   )
{/*lint !e1784*/
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
{/*lint !e1784*/
   SDPA::PhaseType phasetype;

   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL);
   CHECK_IF_SOLVED_BOOL( sdpisolver );

   phasetype = sdpisolver->sdpa->getPhaseValue();

   if ( phasetype == SDPA::noINFO || phasetype == SDPA::pFEAS || phasetype == SDPA::dFEAS || phasetype == SDPA::pdINF )
      return FALSE;

   return TRUE;
}

/** gets information about primal and dual feasibility of the current SDP solution */
SCIP_RETCODE SCIPsdpiSolverGetSolFeasibility(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   SCIP_Bool*            primalfeasible,     /**< stores primal feasibility status */
   SCIP_Bool*            dualfeasible        /**< stores dual feasibility status */
   )
{/*lint !e1784*/
   SDPA::PhaseType phasetype;

   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL );
   assert( primalfeasible != NULL );
   assert( dualfeasible != NULL );
   CHECK_IF_SOLVED( sdpisolver );

   phasetype = sdpisolver->sdpa->getPhaseValue();

   switch ( phasetype )/*lint --e{788}*/
   {
   case SDPA::pdOPT:
      *primalfeasible = TRUE;
      *dualfeasible = TRUE;
      break;
   case SDPA::pdFEAS:
      *primalfeasible = TRUE;
      *dualfeasible = TRUE;
      break;
   case SDPA::pFEAS_dINF:
      *primalfeasible = TRUE;
      *dualfeasible = FALSE;
      break;
   case SDPA::pINF_dFEAS:
      *primalfeasible = FALSE;
      *dualfeasible = TRUE;
      break;
   case SDPA::pUNBD:
      *primalfeasible = TRUE;
      *dualfeasible = FALSE;
      SCIPdebugMessage("SDPA stopped because dual objective became smaller than lower bound\n");
      break;
   case SDPA::dUNBD:
      *primalfeasible = FALSE;
      *dualfeasible = TRUE;
      SCIPdebugMessage("SDPA stopped because primal objective became bigger than upper bound\n");
      break;
   default: /* contains noInfo, pFeas, dFeas, pdInf */
      SCIPerrorMessage("SDPA doesn't know if primal and dual solutions are feasible\n");
      return SCIP_LPERROR;
   }

   return SCIP_OKAY;
}

/** returns TRUE iff SDP is proven to be primal unbounded,
 *  returns FALSE with a debug-message if the solver could not determine feasibility */
SCIP_Bool SCIPsdpiSolverIsPrimalUnbounded(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP-solver interface */
   )
{/*lint !e1784*/
   SDPA::PhaseType phasetype;

   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL );
   CHECK_IF_SOLVED_BOOL( sdpisolver );

   phasetype = sdpisolver->sdpa->getPhaseValue();

   if ( phasetype == SDPA::noINFO || phasetype == SDPA::pFEAS || phasetype == SDPA::pdINF )
   {
      SCIPdebugMessage("SDPA doesn't know if primal problem is unbounded");
      return FALSE;
   }
   else if ( phasetype ==  SDPA::pFEAS_dINF )
      return TRUE;
   else if ( phasetype == SDPA::pUNBD )
   {
      SCIPdebugMessage("SDPA was stopped because primal objective became bigger than upper bound");
      return TRUE;
   }

   return FALSE;
}

/** returns TRUE iff SDP is proven to be primal infeasible,
 *  returns FALSE with a debug-message if the solver could not determine feasibility */
SCIP_Bool SCIPsdpiSolverIsPrimalInfeasible(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP-solver interface */
   )
{/*lint !e1784*/
   SDPA::PhaseType phasetype;

   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL );
   CHECK_IF_SOLVED_BOOL( sdpisolver );

   phasetype = sdpisolver->sdpa->getPhaseValue();

   if ( phasetype == SDPA::noINFO || phasetype == SDPA::dFEAS || phasetype == SDPA::pdINF )
   {
      SCIPdebugMessage("SDPA doesn't know if primal problem is infeasible");
      return FALSE;
   }
   else if ( phasetype ==  SDPA::pINF_dFEAS )
      return TRUE;
   else if ( phasetype == SDPA::dUNBD )
   {
      SCIPdebugMessage("SDPA was stopped because dual objective became smaller than lower bound");
      return TRUE;
   }

   return FALSE;
}

/** returns TRUE iff SDP is proven to be primal feasible,
 *  returns FALSE with a debug-message if the solver could not determine feasibility */
SCIP_Bool SCIPsdpiSolverIsPrimalFeasible(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP-solver interface */
   )
{/*lint !e1784*/
   SDPA::PhaseType phasetype;

   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL );
   CHECK_IF_SOLVED_BOOL( sdpisolver );

   phasetype = sdpisolver->sdpa->getPhaseValue();

   if ( phasetype == SDPA::noINFO || phasetype == SDPA::dFEAS || phasetype == SDPA::pdINF )
   {
      SCIPdebugMessage("SDPA doesn't know if primal problem is feasible");
      return FALSE;
   }
   else if ( phasetype ==  SDPA::pFEAS_dINF || phasetype == SDPA::pdOPT || phasetype == SDPA::pFEAS  || phasetype == SDPA::pdFEAS )
      return TRUE;
   else if ( phasetype == SDPA::dUNBD )
   {
      SCIPdebugMessage("SDPA was stopped because dual objective became smaller than lower bound");
      return TRUE;
   }

   return FALSE;
}

/** returns TRUE iff SDP is proven to be dual unbounded,
 *  returns FALSE with a debug-message if the solver could not determine feasibility */
SCIP_Bool SCIPsdpiSolverIsDualUnbounded(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP-solver interface */
   )
{/*lint !e1784*/
   SDPA::PhaseType phasetype;

   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL);
   CHECK_IF_SOLVED_BOOL( sdpisolver );

   phasetype = sdpisolver->sdpa->getPhaseValue();

   if ( phasetype == SDPA::noINFO || phasetype == SDPA::dFEAS || phasetype == SDPA::pdINF )
   {
      SCIPdebugMessage("SDPA doesn't know if dual problem is unbounded");
      return FALSE;
   }
   else if ( phasetype ==  SDPA::pINF_dFEAS )
      return TRUE;
   else if ( phasetype == SDPA::dUNBD )
   {
      SCIPdebugMessage("SDPA was stopped because dual objective became smaller than lower bound");
      return TRUE;
   }

   return FALSE;
}

/** returns TRUE iff SDP is proven to be dual infeasible,
 *  returns FALSE with a debug-message if the solver could not determine feasibility */
SCIP_Bool SCIPsdpiSolverIsDualInfeasible(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP-solver interface */
   )
{/*lint !e1784*/
   SDPA::PhaseType phasetype;

   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL);
   CHECK_IF_SOLVED_BOOL( sdpisolver );

   phasetype = sdpisolver->sdpa->getPhaseValue();

   if ( phasetype == SDPA::noINFO || phasetype == SDPA::pFEAS || phasetype == SDPA::pdINF )
   {
      SCIPdebugMessage("SDPA doesn't know if dual problem is infeasible");
      return FALSE;
   }
   else if ( phasetype ==  SDPA::pFEAS_dINF )
      return TRUE;
   else if ( phasetype == SDPA::pUNBD )
   {
      SCIPdebugMessage("SDPA was stopped because primal objective became bigger than upper bound");
      return TRUE;
   }

   return FALSE;
}

/** returns TRUE iff SDP is proven to be dual feasible,
 *  returns FALSE with a debug-message if the solver could not determine feasibility */
SCIP_Bool SCIPsdpiSolverIsDualFeasible(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP-solver interface */
   )
{/*lint !e1784*/
   SDPA::PhaseType phasetype;

   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL);
   CHECK_IF_SOLVED_BOOL( sdpisolver );

   phasetype = sdpisolver->sdpa->getPhaseValue();

   if ( phasetype == SDPA::noINFO || phasetype == SDPA::dFEAS || phasetype == SDPA::pdINF )
   {
      SCIPdebugMessage("SDPA doesn't know if primal problem is feasible");
      return FALSE;
   }
   else if ( phasetype ==  SDPA::pINF_dFEAS || phasetype == SDPA::pdOPT || phasetype == SDPA::dFEAS  || phasetype == SDPA::pdFEAS )
      return TRUE;
   else if ( phasetype == SDPA::dUNBD )
   {
      SCIPdebugMessage("SDPA was stopped because dual objective became smaller than lower bound");
      return TRUE;
   }

   return FALSE;
}

/** returns TRUE iff the solver converged */
SCIP_Bool SCIPsdpiSolverIsConverged(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP-solver interface */
   )
{/*lint !e1784*/
   SDPA::PhaseType phasetype;

   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL );
   CHECK_IF_SOLVED_BOOL( sdpisolver );

   phasetype = sdpisolver->sdpa->getPhaseValue();

   if ( phasetype == SDPA::pdOPT )
      return TRUE;

   return FALSE;
}

/** returns TRUE iff the objective limit was reached */
SCIP_Bool SCIPsdpiSolverIsObjlimExc(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP-solver interface */
   )
{/*lint !e1784*/
   SDPA::PhaseType phasetype;

   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL );
   CHECK_IF_SOLVED_BOOL( sdpisolver );

   phasetype = sdpisolver->sdpa->getPhaseValue();

   if ( phasetype == SDPA::pUNBD )
      return TRUE;

   return FALSE;
}

/** returns TRUE iff the iteration limit was reached */
SCIP_Bool SCIPsdpiSolverIsIterlimExc(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP-solver interface */
   )
{/*lint !e1784*/
   SDPA::PhaseType phasetype;

   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL);
   CHECK_IF_SOLVED_BOOL( sdpisolver );

   phasetype = sdpisolver->sdpa->getPhaseValue();

   if ( phasetype == SDPA::noINFO || phasetype == SDPA::pFEAS || phasetype == SDPA::dFEAS || phasetype == SDPA::pdFEAS )
   {
      if ( sdpisolver->sdpa->getParameterMaxIteration() == sdpisolver->sdpa->getIteration() )
         return TRUE;
   }

   return FALSE;
}

/** returns TRUE iff the time limit was reached */
SCIP_Bool SCIPsdpiSolverIsTimelimExc(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP-solver interface */
   )
{/*lint !e1784*/
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
{/*lint !e1784*/
   SDPA::PhaseType phasetype;

   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL );

   if ( sdpisolver->sdpa == NULL || (! sdpisolver->solved) )
      return -1;

   phasetype = sdpisolver->sdpa->getPhaseValue();

   if ( phasetype == SDPA::pdOPT || phasetype == SDPA::pFEAS_dINF || phasetype == SDPA::pINF_dFEAS )
      return 0;
   if ( phasetype == SDPA::pdINF )
      return 1;
   if ( phasetype == SDPA::pUNBD)
      return 3;
   if ( phasetype == SDPA::noINFO || phasetype == SDPA::pFEAS || phasetype == SDPA::dFEAS || phasetype == SDPA::pdFEAS )
      return 4;
   else /* should include dUNBD */
      return 7;
}

/** returns TRUE iff SDP was solved to optimality, meaning the solver converged and returned primal and dual feasible solutions */
SCIP_Bool SCIPsdpiSolverIsOptimal(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP-solver interface */
   )
{/*lint !e1784*/
   SDPA::PhaseType phasetype;

   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL );

   if ( ! sdpisolver->solved )
      return FALSE;

   phasetype = sdpisolver->sdpa->getPhaseValue();

   if ( phasetype == SDPA::pdOPT )
      return TRUE;

   return FALSE;
}

/** returns TRUE iff SDP was solved to optimality or some other status was reached
 *  that is still acceptable inside a Branch & Bound framework */
SCIP_Bool SCIPsdpiSolverIsAcceptable(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP-solver interface */
   )
{/*lint !e1784*/
   SDPA::PhaseType phasetype;

   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL );

   if ( sdpisolver->timelimit )
      return FALSE;

   if ( ! sdpisolver->solved )
      return FALSE;

   phasetype = sdpisolver->sdpa->getPhaseValue();

   /* we are happy if we converged, or we reached the objective limit (pUNBD) or we could show that our problem is
    * infeasible [except for numerics], or unbounded */
   if ( SCIPsdpiSolverIsConverged(sdpisolver) || phasetype == SDPA::pUNBD || phasetype == SDPA::pINF_dFEAS || phasetype == SDPA::pFEAS_dINF )
      return TRUE;

   return FALSE;
}

/** tries to reset the internal status of the SDP-solver in order to ignore an instability of the last solving call */
SCIP_RETCODE SCIPsdpiSolverIgnoreInstability(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   SCIP_Bool*            success             /**< pointer to store, whether the instability could be ignored */
   )
{/*lint !e1784*/
   SCIPdebugMessage("Not implemented yet\n");

   return SCIP_LPERROR;
}/*lint !e715*/

/** gets objective value of solution */
SCIP_RETCODE SCIPsdpiSolverGetObjval(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   SCIP_Real*            objval              /**< pointer to store the objective value */
   )
{/*lint !e1784*/
   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL );
   assert( objval != NULL );
   CHECK_IF_SOLVED( sdpisolver );

   *objval = sdpisolver->sdpa->getPrimalObj();

#ifndef NDEBUG
   SCIP_Real primalval = sdpisolver->sdpa->getDualObj();
   SCIP_Real gap = (REALABS(*objval - primalval) / (0.5 * (REALABS(primalval) + REALABS(*objval)))); /* duality gap used in SDPA */
   if ( gap > sdpisolver->gaptol )
      SCIPdebugMessage("Attention: got objective value (before adding values of fixed variables) of %f in SCIPsdpiSolverGetObjval, "
            "but primal objective is %f with duality gap %f!\n", *objval, primalval, gap );
#endif

   /* as we didn't add the fixed (lb = ub) variables to sdpa, we have to add their contributions to the objective by hand */
   *objval += sdpisolver->fixedvarsobjcontr;

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
{/*lint !e1784*/
   SCIP_Real* sdpasol;
   int v;

   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL );
   assert( dualsollength != NULL );
   CHECK_IF_SOLVED( sdpisolver );

   if ( objval != NULL )
   {
      *objval = sdpisolver->sdpa->getPrimalObj();

#ifndef NDEBUG
      SCIP_Real primalval = sdpisolver->sdpa->getDualObj();
      SCIP_Real gap = (REALABS(*objval - primalval) / (0.5 * (REALABS(primalval) + REALABS(*objval)))); /* duality gap used in SDPA */
      if ( gap > sdpisolver->gaptol )
      {
         SCIPdebugMessage("Attention: got objective value (before adding values of fixed variables) of %f in SCIPsdpiSolverGetSol, "
            "but primal objective is %f with duality gap %f!\n", *objval, primalval, gap );
      }
#endif

      /* as we didn't add the fixed (lb = ub) variables to sdpa, we have to add their contributions to the objective by hand */
      *objval += sdpisolver->fixedvarsobjcontr;
   }

   if ( *dualsollength > 0 )
   {
      assert( dualsol != NULL );
      if ( *dualsollength < sdpisolver->nvars )
      {
         SCIPdebugMessage("The given array in SCIPsdpiSolverGetSol only had length %d, but %d was needed", *dualsollength, sdpisolver->nvars);
         *dualsollength = sdpisolver->nvars;

         return SCIP_OKAY;
      }

      /* get the solution from sdpa */
      assert( (sdpisolver->penalty && sdpisolver->nactivevars + 1 == sdpisolver->sdpa->getConstraintNumber()) || /*lint !e776*/
               sdpisolver->nactivevars == sdpisolver->sdpa->getConstraintNumber() ); /* in the second case we have r as an additional variable */
      sdpasol = sdpisolver->sdpa->getResultXVec();
      /* insert the entries into dualsol, for non-fixed vars we copy those from sdpa, the rest are the saved entries from inserting (they equal lb=ub) */
      for (v = 0; v < sdpisolver->nvars; v++)
      {
         if ( sdpisolver->inputtosdpamapper[v] > 0 )
         {
            /* minus one because the inputtosdpamapper gives the sdpa indices which start at one, but the array starts at 0 */
            dualsol[v] = sdpasol[sdpisolver->inputtosdpamapper[v] - 1];
         }
         else
         {
            /* this is the value that was saved when inserting, as this variable has lb=ub */
            assert( -sdpisolver->inputtosdpamapper[v] <= sdpisolver->nvars - sdpisolver->nactivevars );
            dualsol[v] = sdpisolver->fixedvarsval[(-1 * sdpisolver->inputtosdpamapper[v]) - 1]; /*lint !e679*/
         }
      }
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
{/*lint !e1784*/
   int i;
   SCIP_Real* X; /* block of primal solution matrix corresponding to the LP-part */
   int lpblockind;
   int nlpcons;

   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL );
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

   /* initialize the return-arrays with zero */
   for (i = 0; i < sdpisolver->nvars; i++)
   {
      lbvars[i] = 0.0;
      ubvars[i] = 0.0;
   }

   /* if no variable bounds were added, we return the zero vector (we do this separately, because in this case there might be no LP-block) */
   if ( sdpisolver->nvarbounds == 0 )
   {
      SCIPdebugMessage("Asked for PrimalBoundVars, but there were no variable bounds in sdpa, returning zero vector !");
      return SCIP_OKAY;
   }

   /* get the block of primal solution matrix corresponding to the LP-part from sdpa */
   lpblockind = (int) sdpisolver->sdpa->getBlockNumber(); /* the LP block is the last one and sdpa counts from one */
   assert( sdpisolver->sdpa->getBlockType((long long) lpblockind) == SDPA::LP );/*lint !e747*/
   nlpcons = (int) sdpisolver->sdpa->getBlockSize((long long) lpblockind);/*lint !e747*/
   assert( nlpcons >= 0 );

   X = sdpisolver->sdpa->getResultYMat((long long) lpblockind);/*lint !e747*/

   /* iterate over all variable bounds and insert the corresponding primal variables in the right positions of the return-arrays */
   assert( sdpisolver->nvarbounds <= 2 * sdpisolver->nvars || (sdpisolver->nvarbounds <= 2 * sdpisolver->nvars + 1 && sdpisolver->penalty ) );
   /* if we solved a penalty formulation, the last variable bound belongs to the penalty variable, which we aren't interested in here */
   for (i = 0; i < ((sdpisolver->penalty) ? sdpisolver->nvarbounds - 1 : sdpisolver->nvarbounds); i++)
   {
      if ( sdpisolver->varboundpos[i] < 0 )
      {
         /* this is a lower bound */
         /* the last nvarbounds entries correspond to the varbounds */
         lbvars[sdpisolver->sdpatoinputmapper[-1 * sdpisolver->varboundpos[i] -1]] = X[nlpcons - sdpisolver->nvarbounds + i]; /*lint !e679, !e834 */
      }
      else
      {
         /* this is an upper bound */
         /* the last nvarbounds entries correspond to the varbounds */
         ubvars[sdpisolver->sdpatoinputmapper[+1 * sdpisolver->varboundpos[i] - 1]] = X[nlpcons - sdpisolver->nvarbounds + i]; /*lint !e679, !e834 */
      }
   }

   return SCIP_OKAY;
}

/** gets the number of SDP iterations of the last solve call */
SCIP_RETCODE SCIPsdpiSolverGetIterations(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   int*                  iterations          /**< pointer to store the number of iterations of the last solve call */
   )
{/*lint !e1784*/
   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL );
   assert( iterations != NULL );

   *iterations = sdpisolver->niterations;

   return SCIP_OKAY;
}

/** gets the number of calls to the SDP-solver for the last solve call */
SCIP_RETCODE SCIPsdpiSolverGetSdpCalls(
   SCIP_SDPISOLVER*      sdpisolver,         /**< SDP-solver interface */
   int*                  calls               /**< pointer to store the number of calls to the SDP-solver for the last solve call */
   )
{/*lint !e1784*/
   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL );
   assert( calls != NULL );

   *calls = sdpisolver->nsdpcalls;

   return SCIP_OKAY;
}

/** gets the settings used by the SDP solver for the last solve call */
SCIP_RETCODE SCIPsdpiSolverSettingsUsed(
   SCIP_SDPISOLVER*      sdpisolver,         /**< SDP-solver interface */
   SCIP_SDPSOLVERSETTING* usedsetting        /**< the setting used by the SDP-solver */
   )
{/*lint !e1784*/
   assert( sdpisolver != NULL );
   assert( usedsetting != NULL );
   CHECK_IF_SOLVED(sdpisolver);

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
{/*lint !e1784*/
   return 1E+20; /* default infinity from SCIP */
}/*lint !e715*/

/** checks if given value is treated as (plus or minus) infinity in the SDP-solver */
SCIP_Bool SCIPsdpiSolverIsInfinity(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   SCIP_Real             val                 /**< value to be checked for infinity */
   )
{/*lint !e1784*/
   return ( val <= -SCIPsdpiSolverInfinity(sdpisolver) || val >= SCIPsdpiSolverInfinity(sdpisolver) );
}

/** gets floating point parameter of SDP-Solver */
SCIP_RETCODE SCIPsdpiSolverGetRealpar(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   SCIP_SDPPARAM         type,               /**< parameter number */
   SCIP_Real*            dval                /**< buffer to store the parameter value */
   )
{/*lint !e1784*/
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
      *dval = 0.0;
      SCIPdebugMessage("Parameter SCIP_SDPPAR_PENALTYPARAM not used by SDPA"); /* this parameter is only used by DSDP */
      break;
   case SCIP_SDPPAR_OBJLIMIT:
      *dval = sdpisolver->objlimit;
      break;
   case SCIP_SDPPAR_LAMBDASTAR:
      *dval = sdpisolver->lambdastar;
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
{/*lint !e1784*/
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
      SCIPdebugMessage("Parameter SCIP_SDPPAR_PENALTYPARAM not used by SDPA"); /* this parameter is only used by DSDP */
      break;
   case SCIP_SDPPAR_OBJLIMIT:
      SCIPdebugMessage("Setting sdpisolver objlimit to %f.\n", dval);
      sdpisolver->objlimit = dval;
      break;
   case SCIP_SDPPAR_LAMBDASTAR:
      SCIPdebugMessage("Setting sdpisolver lambdastar parameter to %f.\n", dval);
      sdpisolver->lambdastar = dval;
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
{/*lint !e1784*/
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
{/*lint !e1784*/
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
{/*lint !e1784*/
   SCIP_Real compval;

   assert( sdpisolver != NULL );

   /* we set the value to min{max{MIN_LAMBDASTAR, LAMBDASTAR_FACTOR * MAX_GUESS}, MAX_LAMBDASTAR}, where MAX_GUESS is the maximum of the guesses
    * of the SDP-Blocks, if the define LAMBDASTAR_TWOPOINTS is set, we instead choose either LAMBDASTAR_LOW or HIGH depending on LAMBASTAR_THRESHOLD */

#ifdef LAMBDASTAR_TWOPOINTS
   if ( maxguess < LAMBDASTAR_THRESHOLD )
      compval = LAMBDASTAR_LOW;
   else
      compval = LAMBDASTAR_HIGH;
#else
      compval = LAMBDASTAR_FACTOR * maxguess;
#endif

   if ( compval < MIN_LAMBDASTAR )
   {
      sdpisolver->lambdastar = MIN_LAMBDASTAR;
      SCIPdebugMessage("Setting lambdastar to %f.\n", MIN_LAMBDASTAR);
   }
   else if ( compval > MAX_LAMBDASTAR )
   {
      sdpisolver->lambdastar = MAX_LAMBDASTAR;
      SCIPdebugMessage("Setting lambdastar to %f.\n", MAX_LAMBDASTAR);
   }
   else
   {
      sdpisolver->lambdastar = compval;
      SCIPdebugMessage("Setting lambdastar to %f.\n", compval);
   }

   return SCIP_OKAY;
}

/** compute and set the penalty parameter */
SCIP_RETCODE SCIPsdpiSolverComputePenaltyparam(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   SCIP_Real             maxcoeff,           /**< maximum objective coefficient */
   SCIP_Real*            penaltyparam        /**< the computed penalty parameter */
   )
{/*lint !e1784*/
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
{/*lint !e1784*/
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
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   const char*           fname               /**< file name */
   )
{/*lint !e1784*/
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_LPERROR;
}/*lint !e715*/

/** writes SDP to a file */
SCIP_RETCODE SCIPsdpiSolverWriteSDP(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   const char*           fname               /**< file name */
   )
{/*lint !e1784*/
   assert( fname != NULL );

   sdpisolver->sdpa->writeInputSparse(const_cast<char*>(fname), const_cast<char*>("%8.3f"));

   return SCIP_OKAY;
}

/**@} */
