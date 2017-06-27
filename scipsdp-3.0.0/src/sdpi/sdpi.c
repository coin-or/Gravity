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
/* #define SCIP_MORE_DEBUG*/

/**@file   sdpi.c
 * @brief  General interface methods for SDP-preprocessing (mainly fixing variables and removing empty rows/cols)
 * @author Tristan Gally
 */
#include <assert.h>
#include <time.h>

#include "sdpi/sdpisolver.h"
#include "sdpi/sdpi.h"
#include "scipsdp/SdpVarfixer.h"
#include "sdpi/lapack.h"                     /* to check feasibility if all variables are fixed during preprocessing */

#include "blockmemshell/memory.h"            /* for memory allocation */
#include "scip/def.h"                        /* for SCIP_Real, _Bool, ... */
#include "scip/pub_misc.h"                   /* for sorting */

/* turn off lint warnings for whole file: */
/*lint --e{788,818}*/


/** Checks if a BMSallocMemory-call was successfull, otherwise returns SCIP_NOMEMORY */
#define BMS_CALL(x)   do                                                                                      \
                      {                                                                                       \
                          if( NULL == (x) )                                                                   \
                          {                                                                                   \
                             SCIPerrorMessage("No memory in function call\n");                                \
                             return SCIP_NOMEMORY;                                                            \
                          }                                                                                   \
                      }                                                                                       \
                      while( FALSE )

/** this will be called in all functions that want to access solution information to check if the problem was solved since the last change of the problem */
#define CHECK_IF_SOLVED(sdpi)  do                                                                             \
                      {                                                                                       \
                         if ( ! (sdpi->solved) )                                                              \
                         {                                                                                    \
                            SCIPerrorMessage("Tried to access solution information ahead of solving! \n");    \
                            return SCIP_LPERROR;                                                              \
                         }                                                                                    \
                      }                                                                                       \
                      while( FALSE )

/** same as CHECK_IF_SOLVED, but this will be used in functions returning a boolean value */
#define CHECK_IF_SOLVED_BOOL(sdpi)  do                                                                        \
                      {                                                                                       \
                         if ( ! (sdpi->solved) )                                                              \
                         {                                                                                    \
                            SCIPerrorMessage("Tried to access solution information ahead of solving! \n");    \
                            return FALSE;                                                                     \
                         }                                                                                    \
                      }                                                                                       \
                      while( FALSE )

/** duplicate an array that might be null (in that case null is returned, otherwise BMSduplicateMemory is called) */
#define DUPLICATE_ARRAY_NULL(blkmem, target, source, size) do                                                 \
                      {                                                                                       \
                         if (size > 0)                                                                        \
                            BMS_CALL( BMSduplicateBlockMemoryArray(blkmem, target, source, size) );           \
                         else                                                                                 \
                            *target = NULL;                                                                   \
                      }                                                                                       \
                      while( FALSE )

/** same as SCIP_CALL, but gives a SCIP_PARAMETERUNKNOWN error if it fails */
#define SCIP_CALL_PARAM(x)   do                                                                               \
                      {                                                                                       \
                         SCIP_RETCODE _restat_;                                                               \
                         if ( (_restat_ = (x)) != SCIP_OKAY )                                                 \
                         {                                                                                    \
                            if ( _restat_ != SCIP_PARAMETERUNKNOWN )                                          \
                            {                                                                                 \
                               SCIPerrorMessage("Error <%d> in function call\n", _restat_);                   \
                               SCIPABORT();                                                                   \
                            }                                                                                 \
                            return _restat_;                                                                  \
                         }                                                                                    \
                      }                                                                                       \
                      while( FALSE )

/** same as SCIP_CALL_PARAM, but ignores SCIP_PARAMETERUNKNOWN */
#define SCIP_CALL_PARAM_IGNORE_UNKNOWN(x)   do                                                                \
                      {                                                                                       \
                         SCIP_RETCODE _restat_;                                                               \
                         if ( (_restat_ = (x)) != SCIP_OKAY )                                                 \
                         {                                                                                    \
                            if ( _restat_ != SCIP_PARAMETERUNKNOWN )                                          \
                            {                                                                                 \
                               SCIPerrorMessage("Error <%d> in function call\n", _restat_);                   \
                               SCIPABORT();                                                                   \
                            }                                                                                 \
                         }                                                                                    \
                      }                                                                                       \
                      while( FALSE )

/* #define PRINTSLATER */
#define MIN_GAPTOL                  1e-10    /**< minimum gaptolerance for SDP-solver if decreasing it for a penalty formulation */

#define DEFAULT_SDPSOLVERGAPTOL     1e-4     /**< the stopping criterion for the duality gap the sdpsolver should use */
#define DEFAULT_FEASTOL             1e-6     /**< used to test for feasibility */
#define DEFAULT_EPSILON             1e-9     /**< used to test whether given values are equal */
#define DEFAULT_PENALTYPARAM        1e+5     /**< the starting penalty parameter Gamma used for the penalty formulation if the SDP-solver didn't converge */
#define DEFAULT_MAXPENALTYPARAM     1e+10    /**< the maximum penalty parameter Gamma used for the penalty formulation if the SDP-solver didn't converge */
#define DEFAULT_NPENALTYINCR        8        /**< maximum number of times the penalty parameter will be increased if penalty formulation failed */

/** data for SDPI */
struct SCIP_SDPi
{
   SCIP_SDPISOLVER*      sdpisolver;         /**< pointer to the interface for the SDP-solver */
   SCIP_MESSAGEHDLR*     messagehdlr;        /**< messagehandler to printing messages, or NULL */
   BMS_BLKMEM*           blkmem;             /**< block memory */
   BMS_BUFMEM*           bufmem;             /**< buffer memory */
   int                   nvars;              /**< number of variables */
   SCIP_Real*            obj;                /**< objective function values of variables */
   SCIP_Real*            lb;                 /**< lower bounds of variables */
   SCIP_Real*            ub;                 /**< upper bounds of variables */
   int                   nsdpblocks;         /**< number of SDP-blocks */
   int*                  sdpblocksizes;      /**< sizes of the SDP-blocks */
   int*                  sdpnblockvars;      /**< number of variables in each SDP-block */

   /* constant SDP data: */
   int                   sdpconstnnonz;      /**< number of nonzero elements in the constant matrices of the SDP-Blocks */
   int*                  sdpconstnblocknonz; /**< number of nonzeros for each variable in the constant part, also the i-th entry gives the
                                              *   number of entries  of sdpconst row/col/val [i] */
   int**                 sdpconstrow;        /**< pointers to row-indices for each block */
   int**                 sdpconstcol;        /**< pointers to column-indices for each block */
   SCIP_Real**           sdpconstval;        /**< pointers to the values of the nonzeros for each block */

   /* non-constant SDP data: */
   int                   sdpnnonz;           /**< number of nonzero elements in the SDP-constraint matrices */
   int**                 sdpnblockvarnonz;   /**< sdpnblockvarnonz[i][j] gives the number of nonzeros for the j-th variable (not necessarly
                                              *   variable j) in the i-th block, this is also the length of row/col/val[i][j] */
   int**                 sdpvar;             /**< sdpvar[i][j] gives the sdp-index of the j-th variable (according to the sorting for row/col/val)
                                              *   in the i-th block */
   int***                sdprow;             /**< pointer to the row-indices for each block and variable in this block, so row[i][j][k] gives
                                              *   the k-th nonzero of the j-th variable (not necessarly variable j) in the i-th block */
   int***                sdpcol;             /**< pointer to the column-indices for each block and variable in this block */
   SCIP_Real***          sdpval;             /**< pointer to the values of the nonzeros for each block and variable in this block */

   /* lp data: */
   int                   nlpcons;            /**< number of LP-constraints */
   SCIP_Real*            lplhs;              /**< left hand sides of LP rows */
   SCIP_Real*            lprhs;              /**< right hand sides of LP rows */
   int                   lpnnonz;            /**< number of nonzero elements in the LP-constraint matrix */
   int*                  lprow;              /**< row-index for each entry in lpval-array */
   int*                  lpcol;              /**< column-index for each entry in lpval-array */
   SCIP_Real*            lpval;              /**< values of LP-constraint matrix entries */

   /* other data */
   int                   slatercheck;        /**< should the Slater condition for the dual problem be checked ahead of each solving process */
   int                   sdpid;              /**< counter for the number of SDPs solved */
   int                   niterations;        /**< number of iterations since the last solve call */
   int                   nsdpcalls;          /**< number of calls to the SDP-Solver since the last solve call */
   SCIP_Bool             solved;             /**< was the problem solved since the last change */
   SCIP_Bool             penalty;            /**< was the last solved problem a penalty formulation */
   SCIP_Bool             infeasible;         /**< was infeasibility detected in presolving? */
   SCIP_Bool             allfixed;           /**< could all variables be fixed during presolving? */
   SCIP_Real             epsilon;            /**< tolerance for absolute checks */
   SCIP_Real             gaptol;             /**< (previous: sdpsolverepsilon) this is used for checking if primal and dual objective are equal */
   SCIP_Real             feastol;            /**< this is used to check if the SDP-Constraint is feasible */
   SCIP_Real             penaltyparam;       /**< the starting penalty parameter Gamma used for the penalty formulation if the SDP-solver didn't converge */
   SCIP_Real             maxpenaltyparam;    /**< the maximum penalty parameter Gamma used for the penalty formulation if the SDP-solver didn't converge */
   int                   npenaltyincr;       /**< maximum number of times the penalty parameter will be increased if penalty formulation failed */
   SCIP_Real             bestbound;          /**< best bound computed with a penalty formulation */
   SCIP_SDPSLATER        primalslater;       /**< did the primal slater condition hold for the last problem */
   SCIP_SDPSLATER        dualslater;         /**< did the dual slater condition hold for the last problem */
};

/*
 * Local Functions
 */

/** For given row and column (i,j) checks if i >= j, so that i and j give a position in the lower
 *  triangular part, otherwise i and j will be switched. This function will be called whenever a position in a symmetric matrix
 *  is given, to prevent problems if position (i,j) is given but later (j,i) should be changed.
 */
static
void ensureLowerTriangular(
   int*                  i,                  /**< row index */
   int*                  j                   /**< column index */
   )
{
   if ( *i < *j )
   {
      int temp;
      temp = *i;
      *i = *j;
      *j = temp;
   }
}

#ifndef NDEBUG
/** tests if for a given variable the lower bound is in an epsilon neighborhood of the upper bound */
static
SCIP_Bool isFixed(
   SCIP_SDPI*            sdpi,               /**< pointer to an SDP-interface structure */
   int                   v                   /**< global index of the variable to check this for */
   )
{
   SCIP_Real lb;
   SCIP_Real ub;

   assert ( sdpi != NULL );
   assert ( v < sdpi->nvars );
   assert ( sdpi->lb != NULL );
   assert ( sdpi->ub != NULL );

   lb = sdpi->lb[v];
   ub = sdpi->ub[v];

   assert( lb < ub + sdpi->feastol || sdpi->infeasible );

   return ( ub-lb <= sdpi->epsilon );
}
#else
#define isFixed(sdpi, v) (sdpi->ub[v] - sdpi->lb[v] <= sdpi->epsilon)
#endif

/** Computes the constant matrix after all variables with lb=ub have been fixed and their nonzeros were moved to the constant part. The five variables
 *  other than sdpi are used to return the matrix.
 *
 *  The size of sdpconstnblocknonz and the first pointers of sdpconst row/col/val should be equal to sdpi->nsdpblocks,
 *  the size of sdpconst row/col/val [i], which is given in sdpconstblocknnonz, needs to be sufficient, otherwise the
 *  needed length will be returned in sdpconstnblocknonz and a debug message will be thrown.
 */
static
SCIP_RETCODE compConstMatAfterFixings(
   SCIP_SDPI*            sdpi,               /**< pointer to an SDP-interface structure */
   int*                  sdpconstnnonz,      /**< pointer to store number of nonzero elements in the constant matrices of the SDP-blocks */
   int*                  sdpconstnblocknonz, /**< pointer to store number of nonzeros for each variable in the constant part, also the i-th entry gives the
                                              *   number of entries  of sdpconst row/col/val [i] */
   int**                 sdpconstrow,        /**< pointer to store row-indices for each block */
   int**                 sdpconstcol,        /**< pointer to store column-indices for each block */
   SCIP_Real**           sdpconstval         /**< pointer to store the values of the nonzeros for each block */
   )
{
   int i;
   int v;
   int block;
   int* nfixednonz;
   int** fixedrows;
   int** fixedcols;
   SCIP_Real** fixedvals;

   assert ( sdpi != NULL );
   assert ( sdpconstnnonz != NULL );
   assert ( sdpconstnblocknonz != NULL );
   assert ( sdpconstrow != NULL );
   assert ( sdpconstcol != NULL );
   assert ( sdpconstval != NULL );
#ifndef NDEBUG
   for (block = 0; block < sdpi->nsdpblocks; block++)
   {
      assert ( sdpconstrow[block] != NULL );
      assert ( sdpconstcol[block] != NULL );
      assert ( sdpconstval[block] != NULL );
   }
#endif

   fixedrows = NULL;
   fixedcols = NULL;
   fixedvals = NULL;

   /* allocate memory for the nonzeros that need to be fixed, as this is only temporarly needed, we allocate as much as theoretically possible */
   BMS_CALL( BMSallocBlockMemoryArray(sdpi->blkmem, &nfixednonz, sdpi->nsdpblocks) );
   BMS_CALL( BMSallocBlockMemoryArray(sdpi->blkmem, &fixedrows, sdpi->nsdpblocks) );
   BMS_CALL( BMSallocBlockMemoryArray(sdpi->blkmem, &fixedcols, sdpi->nsdpblocks) );
   BMS_CALL( BMSallocBlockMemoryArray(sdpi->blkmem, &fixedvals, sdpi->nsdpblocks) );

   for (block = 0; block < sdpi->nsdpblocks; block++)
   {
      /* compute the number of fixed nonzeros in this block */
      nfixednonz[block] = 0;
      for (v = 0; v < sdpi->sdpnblockvars[block]; v++)
      {
         if (isFixed(sdpi, sdpi->sdpvar[block][v]))
            nfixednonz[block] += sdpi->sdpnblockvarnonz[block][v];
      }

      fixedrows[block] = NULL;
      fixedcols[block] = NULL;
      fixedvals[block] = NULL;

      BMS_CALL( BMSallocBlockMemoryArray(sdpi->blkmem, &(fixedrows[block]), nfixednonz[block]) );
      BMS_CALL( BMSallocBlockMemoryArray(sdpi->blkmem, &(fixedcols[block]), nfixednonz[block]) );
      BMS_CALL( BMSallocBlockMemoryArray(sdpi->blkmem, &(fixedvals[block]), nfixednonz[block]) );

      /* set nfixednonz to 0 to use it for indexing later (at the end of the next for-block it will again have the same value) */
      nfixednonz[block] = 0;
   }

   /* iterate over all variables, saving the nonzeros of the fixed ones */
   for (block = 0; block < sdpi->nsdpblocks; block++)
   {
      for (v = 0; v < sdpi->sdpnblockvars[block]; v++)
      {
         if (isFixed(sdpi, sdpi->sdpvar[block][v]))
         {
            for (i = 0; i < sdpi->sdpnblockvarnonz[block][v]; i++)
            {
               fixedrows[block][nfixednonz[block]] = sdpi->sdprow[block][v][i];
               fixedcols[block][nfixednonz[block]] = sdpi->sdpcol[block][v][i];
               /* this is the final value to add, so we no longer have to remember, from which variable this nonzero comes,
                * the -1 comes from +y_iA_i but -A_0 */
               fixedvals[block][nfixednonz[block]] = - sdpi->sdpval[block][v][i] * sdpi->lb[sdpi->sdpvar[block][v]];
               nfixednonz[block]++;
            }
         }
      }
   }

   /* compute the constant matrix */
   *sdpconstnnonz = 0;
   for (block = 0; block < sdpi->nsdpblocks; block++)
   {
      SCIP_CALL( SCIPsdpVarfixerMergeArraysIntoNew(sdpi->blkmem, sdpi->epsilon, sdpi->sdpconstrow[block], sdpi->sdpconstcol[block], sdpi->sdpconstval[block],
                                               sdpi->sdpconstnblocknonz[block], fixedrows[block], fixedcols[block], fixedvals[block], nfixednonz[block],
                                               sdpconstrow[block], sdpconstcol[block], sdpconstval[block], &sdpconstnblocknonz[block]) );
      *sdpconstnnonz += sdpconstnblocknonz[block];
   }

   /* free all memory */
   for (block = 0; block < sdpi->nsdpblocks; block++)
   {
      BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &(fixedvals[block]), nfixednonz[block]);
      BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &(fixedcols[block]), nfixednonz[block]);
      BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &(fixedrows[block]), nfixednonz[block]);
   }
   BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &fixedvals, sdpi->nsdpblocks);
   BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &fixedcols, sdpi->nsdpblocks);
   BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &fixedrows, sdpi->nsdpblocks);
   BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &nfixednonz, sdpi->nsdpblocks);

   return SCIP_OKAY;
}

/** takes the sdpi and the computed constant matrix after fixings as input and checks for empty rows and columns in each block, which should be
 *  removed to not harm the Slater condition. It also removes SDP-blocks with no entries left, these are returned in blockindchanges and nremovedblocks.
 */
static
SCIP_RETCODE findEmptyRowColsSDP(
   SCIP_SDPI*            sdpi,               /**< pointer to an SDP-interface structure */
   int*                  sdpconstnblocknonz, /**< number of nonzeros for each variable in the constant part, also the i-th entry gives the
                                              *   number of entries  of sdpconst row/col/val [i] */
   int**                 sdpconstrow,        /**< pointers to row-indices for each block */
   int**                 sdpconstcol,        /**< pointers to column-indices for each block */
   SCIP_Real**           sdpconstval,        /**< pointers to the values of the nonzeros for each block */
   int**                 indchanges,         /**< pointer to store the changes needed to be done to the indices, if indchange[block][nonz]=-1, then
                                              *   the index can be removed, otherwise it gives the number of indices removed before this, i.e.
                                              *   the value to decrease this index by, this array should have memory allocated in the size
                                              *   sdpi->nsdpblocks times sdpi->sdpblocksizes[block] */
   int*                  nremovedinds,       /**< pointer to store the number of rows/cols to be fixed for each block */
   int*                  blockindchanges,    /**< pointer to store index change for each block, system is the same as for indchanges */
   int*                  nremovedblocks      /**< pointer to store the number of blocks to be removed from the SDP */
   )
{
   int block;
   int v;
   int i;
   int nfoundinds;

   assert( sdpi != NULL );
   assert( sdpconstnblocknonz != NULL );
   assert( sdpconstrow != NULL );
   assert( sdpconstcol != NULL );
   assert( sdpconstval != NULL );
   assert( indchanges != NULL );
   assert( nremovedinds != NULL );
   assert( blockindchanges != NULL );
   assert( nremovedblocks != NULL );

   /* initialize indchanges with -1 */
   for (block = 0; block < sdpi->nsdpblocks; block++)
   {
      for (i = 0; i < sdpi->sdpblocksizes[block]; i++)
         indchanges[block][i] = -1;
   }
   *nremovedblocks = 0;

   /* iterate over all active nonzeros, setting the values of indchange for their row and col to 1 (this is an intermediate value to save that the
    * index is still needed, it will later be set to the number of rows/cols deleted earlier) */
   for (block = 0; block < sdpi->nsdpblocks; block++)
   {
      /* the number of indices already found in this block, saved for prematurely stopping the loops */
      nfoundinds = 0;
      for (v = 0; v < sdpi->sdpnblockvars[block]; v++)
      {
         if ( ! (isFixed(sdpi, sdpi->sdpvar[block][v])) )
         {
            for (i = 0; i < sdpi->sdpnblockvarnonz[block][v]; i++)
            {
               assert ( REALABS(sdpi->sdpval[block][v][i]) > sdpi->epsilon); /* this should really be a nonzero */
               if ( indchanges[block][sdpi->sdprow[block][v][i]] == -1 )
               {
                  indchanges[block][sdpi->sdprow[block][v][i]] = 1;
                  nfoundinds++;
               }
               if ( indchanges[block][sdpi->sdpcol[block][v][i]] == -1 )
               {
                  indchanges[block][sdpi->sdpcol[block][v][i]] = 1;
                  nfoundinds++;
               }
               if ( nfoundinds == sdpi->sdpblocksizes[block] )
                  break;   /* we're done for this block */
            }
         }
         if (nfoundinds == sdpi->sdpblocksizes[block])
            break;   /* we're done for this block */
      }

      if ( nfoundinds < sdpi->sdpblocksizes[block] )
      {
         /* if some indices haven't been found yet, look in the constant part for them */
         for (i = 0; i < sdpconstnblocknonz[block]; i++)
         {
            assert ( REALABS(sdpconstval[block][i]) > sdpi->epsilon); /* this should really be a nonzero */
            if ( indchanges[block][sdpconstrow[block][i]] == -1 )
            {
               indchanges[block][sdpconstrow[block][i]] = 1;
               nfoundinds++;
            }
            if ( indchanges[block][sdpconstcol[block][i]] == -1 )
            {
               indchanges[block][sdpconstcol[block][i]] = 1;
               nfoundinds++;
            }
            if ( nfoundinds == sdpi->sdpblocksizes[block] )
               break;   /* we're done for this block */
         }
      }

      /* now iterate over all indices to compute the final values of indchanges, all 0 are set to -1, all 1 are changed to the number of -1 before it */
      nremovedinds[block] = 0;
      for (i = 0; i < sdpi->sdpblocksizes[block]; i++)
      {
         if ( indchanges[block][i] == -1 )
         {
            SCIPdebugMessage("empty row and col %d were removed from block %d of SDP %d\n", i, block, sdpi->sdpid);
            /* this index wasn't found (indchanges was initialized with 0), so it can be removed */
            nremovedinds[block]++;
         }
         else
         {
            /* this index has been found, so set the value to the number of removed inds before it */
            indchanges[block][i] = nremovedinds[block];
         }
      }

      /* check if the block became empty */
      if ( nremovedinds[block] == sdpi->sdpblocksizes[block] )
      {
         SCIPdebugMessage("empty block %d detected in SDP %d, this will be removed", block, sdpi->sdpid);
         blockindchanges[block] = -1;
         (*nremovedblocks)++;
      }
      else
         blockindchanges[block] = *nremovedblocks;
   }

   return SCIP_OKAY;
}

/** computes the number of active variables for each constraint, thereby detecting constraints that
 *  may be removed, and computes the LP-left- and right-hand-sides after including all locally fixed variables
 *  for all constraints with at least two remaining active variables
 */
static
SCIP_RETCODE computeLpLhsRhsAfterFixings(
   SCIP_SDPI*            sdpi,               /**< pointer to an SDP-interface structure */
   int*                  nactivelpcons,      /**< output: number of active LP-constraints */
   SCIP_Real*            lplhsafterfix,      /**< output: first nlpcons (output) entries give left-hand sides of
                                              *           active lp-constraints after fixing variables, these are
                                              *           in the same relative order as before (with non-active rows
                                              *           removed) */
   SCIP_Real*            lprhsafterfix,      /**< output: first nlpcons (output) entries give right-hand sides of
                                              *  	  active lp-constraints after fixing variables, these are
                                              *  	  in the same relative order as before (with non-active rows
                                              *  	  removed) */
   int*                  rownactivevars,     /**< output: number of active variables for every row */
   SCIP_Bool*            fixingsfound        /**< output: returns true if a variable was fixed during this function call */
   )
{
   int i;
   int c;
   int lastrow = -1;
   int nonzind = -1;
   int nonzcol = -1;
   SCIP_Real nonzval;

   assert( sdpi != NULL );
   assert( nactivelpcons != NULL );
   assert( sdpi->nlpcons == 0 || lplhsafterfix != NULL );
   assert( sdpi->nlpcons == 0 || lprhsafterfix != NULL );
   assert( sdpi->nlpcons == 0 || rownactivevars != NULL );
   assert( sdpi->nlpcons == 0 || fixingsfound != NULL );

   /* if there is no LP-part, there is nothing to do */
   if ( sdpi->nlpcons == 0 || sdpi->lpnnonz == 0 )
   {
      *nactivelpcons = 0;
      return SCIP_OKAY;
   }

   /* initialize rownactivevars */
   for (c = 0; c < sdpi->nlpcons; c++)
      rownactivevars[c] = 0;
   *nactivelpcons = 0;

   for (i = 0; i < sdpi->lpnnonz; i++)
   {
      assert( i == 0 || sdpi->lprow[i-1] <= sdpi->lprow[i] );

      /* we reached a new row */
      if ( sdpi->lprow[i] > lastrow )
      {
         /* if the last row had at least two active variables, we keep the lhs- and rhs-value */
         if ( lastrow >= 0 && rownactivevars[lastrow] > 1 )
            (*nactivelpcons)++;
         else if ( lastrow >= 0 && rownactivevars[lastrow] == 1 )
         {
            assert( 0 <= nonzind && nonzind < sdpi->lpnnonz );

            nonzcol = sdpi->lpcol[nonzind];
            assert( 0 <= nonzcol && nonzcol < sdpi->nvars );

            nonzval = sdpi->lpval[nonzind];
            assert( REALABS(nonzval) > sdpi->epsilon );

            /* we have to check if this is an improvement of the current bound */
            if ( nonzval < 0.0 ) /* we have to compare with the upper bound for lhs and lower bound for rhs */
            {
               /* check for the left-hand-side */
               if ( (lplhsafterfix[*nactivelpcons] > - SCIPsdpiInfinity(sdpi)) &&
                  ( (lplhsafterfix[*nactivelpcons] / nonzval) < sdpi->ub[nonzcol] - sdpi->epsilon) )
               {
                  /* this bound is sharper than the original one */
                  SCIPdebugMessage("empty LP-row %d has been removed from SDP %d, upper bound of variable %d has been sharpened to %f "
                     "(originally %f)\n", lastrow, sdpi->sdpid, nonzcol, lplhsafterfix[*nactivelpcons] / nonzval, sdpi->ub[nonzcol]);
                  sdpi->ub[nonzcol] = lplhsafterfix[*nactivelpcons] / nonzval;

                  /* check if this leads to a fixing of this variable */
                  if ( REALABS(sdpi->lb[nonzcol] - sdpi->ub[nonzcol]) < sdpi->epsilon )
                  {
                     *fixingsfound = TRUE;
                     SCIPdebugMessage("computeLpLhsRhsAfterFixings fixed variable %d to value %f in SDP %d\n",
                        nonzcol, sdpi->lb[nonzcol], sdpi->sdpid);
                  }
                  /* check if this makes the problem infeasible */
                  if (sdpi->ub[nonzcol] < sdpi->lb[nonzcol] - sdpi->feastol)
                  {
                     sdpi->infeasible = TRUE;
                     SCIPdebugMessage("We found an upper bound %f that is lower than the lower bound %f for variable %d, so the problem is infeasible !\n",
                           sdpi->ub[nonzcol], sdpi->lb[nonzcol], nonzcol);
                     return SCIP_OKAY;
                  }
               }
               /* check for the right-hand-side */
               if ( (lprhsafterfix[*nactivelpcons] < SCIPsdpiInfinity(sdpi)) &&
                  ( (lprhsafterfix[*nactivelpcons] / nonzval) > sdpi->lb[nonzcol] + sdpi->epsilon) )
               {
                  /* this bound is sharper than the original one */
                  SCIPdebugMessage("empty LP-row %d has been removed from SDP %d, lower bound of variable %d has been sharpened to %f "
                     "(originally %f)\n", lastrow, sdpi->sdpid, nonzcol, lprhsafterfix[*nactivelpcons] / nonzval, sdpi->lb[nonzcol]);
                  sdpi->lb[nonzcol] = lprhsafterfix[*nactivelpcons] / nonzval;

                  /* check if this leads to a fixing of this variable */
                  if ( REALABS(sdpi->lb[nonzcol] - sdpi->ub[nonzcol]) < sdpi->epsilon )
                  {
                     *fixingsfound = TRUE;
                     SCIPdebugMessage("computeLpLhsRhsAfterFixings fixed variable %d to value %f in SDP %d\n",
                        nonzcol, sdpi->lb[nonzcol], sdpi->sdpid);
                  }

                  /* check if this makes the problem infeasible */
                  if (sdpi->ub[nonzcol] < sdpi->lb[nonzcol] - sdpi->feastol)
                  {
                     sdpi->infeasible = TRUE;
                     SCIPdebugMessage("We found an upper bound %f that is lower than the lower bound %f for variable %d, so the problem is infeasible !\n",
                           sdpi->ub[nonzcol], sdpi->lb[nonzcol], nonzcol);
                     return SCIP_OKAY;
                  }
               }
            }
            else	/* we compare with the lower bound for lhs and upper bound for rhs */
            {
               /* check for the left-hand-side */
               if ( (lplhsafterfix[*nactivelpcons] < SCIPsdpiInfinity(sdpi)) &&
                  ( (lplhsafterfix[*nactivelpcons] / nonzval) > sdpi->lb[nonzcol] + sdpi->epsilon) )
               {
                  /* this bound is sharper than the original one */
                  SCIPdebugMessage("empty LP-row %d has been removed from SDP %d, lower bound of variable %d has been sharpened to %f "
                     "(originally %f)\n", lastrow, sdpi->sdpid, nonzcol, lplhsafterfix[*nactivelpcons] / nonzval, sdpi->lb[nonzcol]);
                  sdpi->lb[nonzcol] = lplhsafterfix[*nactivelpcons] / nonzval;

                  /* check if this leads to a fixing of this variable */
                  if ( REALABS(sdpi->lb[nonzcol] - sdpi->ub[nonzcol]) < sdpi->epsilon )
                  {
                     *fixingsfound = TRUE;
                     SCIPdebugMessage("computeLpLhsRhsAfterFixings fixed variable %d to value %f in SDP %d\n",
                        nonzcol, sdpi->lb[nonzcol], sdpi->sdpid);
                  }

                  /* check if this makes the problem infeasible */
                  if (sdpi->ub[nonzcol] < sdpi->lb[nonzcol] - sdpi->feastol)
                  {
                     sdpi->infeasible = TRUE;
                     SCIPdebugMessage("We found an upper bound %f that is lower than the lower bound %f for variable %d, so the problem is infeasible !\n",
                           sdpi->ub[nonzcol], sdpi->lb[nonzcol], nonzcol);
                     return SCIP_OKAY;
                  }
               }
               /* check for the right-hand-side */
               if ( (lprhsafterfix[*nactivelpcons] > - SCIPsdpiInfinity(sdpi)) &&
                  ( (lprhsafterfix[*nactivelpcons] / nonzval) < sdpi->ub[nonzcol] - sdpi->epsilon) )
               {
                  /* this bound is sharper than the original one */
                  SCIPdebugMessage("empty LP-row %d has been removed from SDP %d, upper bound of variable %d has been sharpened to %f "
                     "(originally %f)\n", lastrow, sdpi->sdpid, nonzcol, lprhsafterfix[*nactivelpcons] / nonzval, sdpi->ub[nonzcol]);
                  sdpi->ub[nonzcol] = lprhsafterfix[*nactivelpcons] / nonzval;

                  /* check if this leads to a fixing of this variable */
                  if ( REALABS(sdpi->lb[nonzcol] - sdpi->ub[nonzcol]) < sdpi->epsilon )
                  {
                     *fixingsfound = TRUE;
                     SCIPdebugMessage("computeLpLhsRhsAfterFixings fixed variable %d to value %f in SDP %d\n",
                        nonzcol, sdpi->lb[nonzcol], sdpi->sdpid);
                  }

                  /* check if this makes the problem infeasible */
                  if (sdpi->ub[nonzcol] < sdpi->lb[nonzcol] - sdpi->feastol)
                  {
                     sdpi->infeasible = TRUE;
                     SCIPdebugMessage("We found an upper bound %f that is lower than the lower bound %f for variable %d, so the problem is infeasible !\n",
                           sdpi->ub[nonzcol], sdpi->lb[nonzcol], nonzcol);
                     return SCIP_OKAY;
                  }
               }
            }
         }
         else if ( lastrow >= 0 ) /* because of earlier ifs we have rownactivevars = 0 */
         {
            assert( lastrow == -1 || rownactivevars[lastrow] == 0 );
            /* we have a constraint lhs <= 0 <= rhs, so lhs should be non-positive and rhs non-negative, otherwise the problem is infeasible */
            if ( lplhsafterfix[*nactivelpcons] > sdpi->feastol || lprhsafterfix[*nactivelpcons] < -sdpi->feastol )
            {
               sdpi->infeasible = TRUE;
               SCIPdebugMessage("We found a constraint which with given fixings reads %f <= 0 <= %f, so the current problem is infeasible !\n",
                     lplhsafterfix[*nactivelpcons], lprhsafterfix[*nactivelpcons] );
               return SCIP_OKAY;
            }
         }

         /* update lastrow for new row */
         lastrow = sdpi->lprow[i];

         /* start the next lhr & rhs with the original value */
         lplhsafterfix[*nactivelpcons] = sdpi->lplhs[lastrow];
         lprhsafterfix[*nactivelpcons] = sdpi->lprhs[lastrow];
      }

      /* if the variable is active, we increase rownactivevars */
      if ( ! isFixed(sdpi, sdpi->lpcol[i]) )
      {
         rownactivevars[lastrow]++;
         nonzind = i;
      }
      else
      {
         /* otherwise we add the value (coefficient * value of fixed variable) to the lhs and rhs, the minus comes from +A_i but -A_0 */
         lplhsafterfix[*nactivelpcons] -= sdpi->lpval[i] * sdpi->lb[sdpi->lpcol[i]];
         lprhsafterfix[*nactivelpcons] -= sdpi->lpval[i] * sdpi->lb[sdpi->lpcol[i]];
      }
   }

   /* for the last row of the lp we have to check if it is active, as in the above for-queue we only do so when the next row start */
   if ( rownactivevars[lastrow] > 1 )
      (*nactivelpcons)++;
   else if ( rownactivevars[lastrow] == 1 )
   {
      assert( 0 <= nonzind && nonzind < sdpi->lpnnonz );

      nonzcol = sdpi->lpcol[nonzind];
      assert( 0 <= nonzcol && nonzcol < sdpi->nvars );

      nonzval = sdpi->lpval[nonzind];
      assert( REALABS(nonzval) > sdpi->epsilon );

      /* we have to check if this is an improvement of the current bound */
      if ( nonzval < 0.0 ) /* we have to compare with the upper bound for lhs and lower bound for rhs */
      {
         /* check for the left-hand-side */
         if ( (lplhsafterfix[*nactivelpcons] > SCIPsdpiInfinity(sdpi)) &&
            ( (lplhsafterfix[*nactivelpcons] / nonzval) < sdpi->ub[nonzcol] - sdpi->epsilon) )
         {
            /* this bound is sharper than the original one */
            SCIPdebugMessage("empty LP-row %d has been removed from SDP %d, upper bound of variable %d has been sharpened to %f "
               "(originally %f)\n", lastrow, sdpi->sdpid, nonzcol, lplhsafterfix[*nactivelpcons] / nonzval, sdpi->ub[nonzcol]);
            sdpi->ub[nonzcol] = lplhsafterfix[*nactivelpcons] / nonzval;

            /* check if this leads to a fixing of this variable */
            if ( REALABS(sdpi->lb[nonzcol] - sdpi->ub[nonzcol]) < sdpi->epsilon )
            {
               *fixingsfound = TRUE;
               SCIPdebugMessage("computeLpLhsRhsAfterFixings fixed variable %d to value %f in SDP %d\n",
                  nonzcol, sdpi->lb[nonzcol], sdpi->sdpid);
            }

            /* check if this makes the problem infeasible */
            if ( sdpi->ub[nonzcol] < sdpi->lb[nonzcol] - sdpi->feastol )
            {
               sdpi->infeasible = TRUE;
               SCIPdebugMessage("We found an upper bound that is lower than the lower bound, so the problem is infeasible !\n");
               return SCIP_OKAY;
            }
         }
         /* check for the right-hand-side */
         if ( (lprhsafterfix[*nactivelpcons] < SCIPsdpiInfinity(sdpi)) &&
            ( (lprhsafterfix[*nactivelpcons] / nonzval) > sdpi->lb[nonzcol] - sdpi->epsilon) )
         {
            /* this bound is sharper than the original one */
            SCIPdebugMessage("empty LP-row %d has been removed from SDP %d, lower bound of variable %d has been sharpened to %f "
               "(originally %f)\n", lastrow, sdpi->sdpid, nonzcol, lprhsafterfix[*nactivelpcons] / nonzval, sdpi->lb[nonzcol]);
            sdpi->lb[nonzcol] = lprhsafterfix[*nactivelpcons] / nonzval;

            /* check if this leads to a fixing of this variable */
            if ( REALABS(sdpi->lb[nonzcol] - sdpi->ub[nonzcol]) < sdpi->epsilon )
            {
               *fixingsfound = TRUE;
               SCIPdebugMessage("computeLpLhsRhsAfterFixings fixed variable %d to value %f in SDP %d\n",
                  nonzcol, sdpi->lb[nonzcol], sdpi->sdpid);
            }

            /* check if this makes the problem infeasible */
            if ( sdpi->ub[nonzcol] < sdpi->lb[nonzcol] - sdpi->feastol )
            {
               sdpi->infeasible = TRUE;
               SCIPdebugMessage("We found an upper bound that is lower than the lower bound, so the problem is infeasible !\n");
               return SCIP_OKAY;
            }
         }
      }
      else  /* we compare with the lower bound for lhs and upper bound for rhs */
      {
         /* check for the left-hand-side */
         if ( (lplhsafterfix[*nactivelpcons] < SCIPsdpiInfinity(sdpi)) &&
            ( (lplhsafterfix[*nactivelpcons] / nonzval) > sdpi->lb[nonzcol] + sdpi->epsilon) )
         {
            /* this bound is sharper than the original one */
            SCIPdebugMessage("empty LP-row %d has been removed from SDP %d, lower bound of variable %d has been sharpened to %f "
               "(originally %f)\n", lastrow, sdpi->sdpid, nonzcol, lplhsafterfix[*nactivelpcons] / nonzval, sdpi->lb[nonzcol]);
            sdpi->lb[nonzcol] = lplhsafterfix[*nactivelpcons] / nonzval;

            /* check if this leads to a fixing of this variable */
            if ( REALABS(sdpi->lb[nonzcol] - sdpi->ub[nonzcol]) < sdpi->epsilon )
            {
               *fixingsfound = TRUE;
               SCIPdebugMessage("computeLpLhsRhsAfterFixings fixed variable %d to value %f in SDP %d\n",
                  nonzcol, sdpi->lb[nonzcol], sdpi->sdpid);
            }

            /* check if this makes the problem infeasible */
            if ( sdpi->ub[nonzcol] < sdpi->lb[nonzcol] - sdpi->feastol )
            {
               sdpi->infeasible = TRUE;
               SCIPdebugMessage("We found a lower bound that is bigger than the upper bound, so the problem is infeasible !\n");
               return SCIP_OKAY;
            }
         }
         /* check for the right-hand-side */
         if ( (lprhsafterfix[*nactivelpcons] > SCIPsdpiInfinity(sdpi)) &&
            ( (lprhsafterfix[*nactivelpcons] / nonzval) < sdpi->ub[nonzcol] - sdpi->epsilon) )
         {
            /* this bound is sharper than the original one */
            SCIPdebugMessage("empty LP-row %d has been removed from SDP %d, upper bound of variable %d has been sharpened to %f "
               "(originally %f)\n", lastrow, sdpi->sdpid, nonzcol, lprhsafterfix[*nactivelpcons] / nonzval, sdpi->ub[nonzcol]);
            sdpi->ub[nonzcol] = lplhsafterfix[*nactivelpcons] / nonzval;

            /* check if this leads to a fixing of this variable */
            if ( REALABS(sdpi->lb[nonzcol] - sdpi->ub[nonzcol]) < sdpi->epsilon )
            {
               *fixingsfound = TRUE;
               SCIPdebugMessage("computeLpLhsRhsAfterFixings fixed variable %d to value %f in SDP %d\n",
                  nonzcol, sdpi->lb[nonzcol], sdpi->sdpid);
            }

            /* check if this makes the problem infeasible */
            if ( sdpi->ub[nonzcol] < sdpi->lb[nonzcol] - sdpi->feastol )
            {
               sdpi->infeasible = TRUE;
               SCIPdebugMessage("We found an upper bound that is lower than the lower bound, so the problem is infeasible !\n");
               return SCIP_OKAY;
            }
         }
      }
   }
   else
   {
      assert( lastrow == -1 || rownactivevars[lastrow] == 0 );
      /* we have a constraint lhs <= 0 <= rhs, so lhs should be non-positive and rhs non-negative, otherwise the problem is infeasible */
      if ( lplhsafterfix[*nactivelpcons] > sdpi->feastol || lprhsafterfix[*nactivelpcons] < -sdpi->feastol )
      {
         sdpi->infeasible = TRUE;
         SCIPdebugMessage("We found a constraint which with given fixings reads %f <= 0 <= %f, so the current problem is infeasible !\n",
               lplhsafterfix[*nactivelpcons], lprhsafterfix[*nactivelpcons] );
         return SCIP_OKAY;
      }
   }

   return SCIP_OKAY;
}

/** checks whether all variables are fixed (lb=ub), in that case changes the sdpi->allfixed pointer accordingly
 */
static
SCIP_RETCODE checkAllFixed(
   SCIP_SDPI*            sdpi                /**< pointer to an SDP-interface structure */
   )
{
   int v;

   /* check all variables until we find an unfixed one */
   for (v = 0; v < sdpi->nvars; v++)
   {
      if ( ! isFixed(sdpi, v) )
      {
         sdpi->allfixed = FALSE;

         return SCIP_OKAY;
      }
   }

   /* we did not find an unfixed variable, so all are fixed */
   SCIPdebugMessage("Detected that all variables in SDP %d are fixed.\n", sdpi->sdpid);
   sdpi->allfixed = TRUE;

   return SCIP_OKAY;
}

/** If all variables are fixed, check whether the remaining solution is feasible for the SDP-constraints (LP-constraints should be checked
 *  already when computing the rhs after fixing)
 */
static
SCIP_RETCODE checkFixedFeasibilitySdp(
   SCIP_SDPI*            sdpi,               /**< pointer to an SDP-interface structure */
   int*                  sdpconstnblocknonz, /**< number of nonzeros for each variable in the constant part, also the i-th entry gives the
                                              *   number of entries  of sdpconst row/col/val [i] */
   int**                 sdpconstrow,        /**< pointers to row-indices for each block */
   int**                 sdpconstcol,        /**< pointers to column-indices for each block */
   SCIP_Real**           sdpconstval,        /**< pointers to the values of the nonzeros for each block */
   int**                 indchanges,         /**< pointer to store the changes needed to be done to the indices, if indchanges[block][nonz]=-1, then
                                              *   the index can be removed, otherwise it gives the number of indices removed before this, i.e.
                                              *   the value to decrease this index by, this array should have memory allocated in the size
                                              *   sdpi->nsdpblocks times sdpi->sdpblocksizes[block] */
   int*                  nremovedinds,       /**< pointer to store the number of rows/cols to be fixed for each block */
   int*                  blockindchanges     /**< pointer to store index change for each block, system is the same as for indchanges */
   )
{
   SCIP_Real* fullmatrix; /* we need to give the full matrix to LAPACK */
   int maxsize; /* as we don't want to allocate memory newly for every SDP-block, we allocate memory according to the size of the largest block */
   SCIP_Real fixedval;
   SCIP_Real eigenvalue;
   int size;
   int b;
   int i;
   int v;

   assert( sdpi->allfixed );

   /* compute the maximum blocksize */
   maxsize = -1;

   for (b = 0; b < sdpi->nsdpblocks; b++)
   {
      if ( sdpi->sdpblocksizes[b] - nremovedinds[b] > maxsize )
         maxsize = sdpi->sdpblocksizes[b] - nremovedinds[b];
   }

   /* allocate memory */
   BMS_CALL( BMSallocBlockMemoryArray(sdpi->blkmem, &fullmatrix, maxsize * maxsize) ); /*lint !e647*/

   /* iterate over all SDP-blocks and check if the smallest eigenvalue is non-negative */
   for (b = 0; b < sdpi->nsdpblocks; b++)
   {
      /* if the block is removed, we don't need to do anything, otherwise build the full matrix */
      if ( blockindchanges[b] == -1 )
         continue;

      size = sdpi->sdpblocksizes[b] - nremovedinds[b];

      /* initialize the matrix with zero */
      for (i = 0; i < size * size; i++)
         fullmatrix[i] = 0.0;

      /* add the constant part (with negative sign) */
      for (i = 0; i < sdpconstnblocknonz[b]; i++)
      {
         assert( 0 <= sdpconstrow[b][i] - indchanges[b][sdpconstrow[b][i]] && sdpconstrow[b][i] - indchanges[b][sdpconstrow[b][i]] < size );
         assert( 0 <= sdpconstcol[b][i] - indchanges[b][sdpconstcol[b][i]] && sdpconstcol[b][i] - indchanges[b][sdpconstcol[b][i]] < size );
         fullmatrix[(sdpconstrow[b][i] - indchanges[b][sdpconstrow[b][i]]) * size
                    + sdpconstcol[b][i] - indchanges[b][sdpconstcol[b][i]]] = -1 * sdpconstval[b][i]; /*lint !e679*/
      }

      /* add the contributions of the fixed variables */
      for (v = 0; v < sdpi->sdpnblockvars[b]; v++)
      {
         fixedval = sdpi->lb[sdpi->sdpvar[b][v]];

         /* if the variable is fixed to zero, we can ignore its contributions */
         if ( REALABS(fixedval) < sdpi->epsilon )
            continue;

         /* iterate over all nonzeros */
         for (i = 0; i < sdpi->sdpnblockvarnonz[b][v]; i++)
         {
            assert( 0 <= sdpi->sdprow[b][v][i] - indchanges[b][sdpi->sdprow[b][v][i]] &&
                         sdpi->sdprow[b][v][i] - indchanges[b][sdpi->sdprow[b][v][i]] < size );
            assert( 0 <= sdpi->sdpcol[b][v][i] - indchanges[b][sdpi->sdpcol[b][v][i]] &&
                         sdpi->sdpcol[b][v][i] - indchanges[b][sdpi->sdpcol[b][v][i]] < size );
            fullmatrix[(sdpi->sdprow[b][v][i] - indchanges[b][sdpi->sdprow[b][v][i]]) * size
                       + sdpi->sdpcol[b][v][i] - indchanges[b][sdpi->sdpcol[b][v][i]]] += fixedval * sdpi->sdpval[b][v][i]; /*lint !e679*/
         }
      }

      /* compute the smallest eigenvalue */
      SCIP_CALL( SCIPlapackComputeIthEigenvalue(sdpi->bufmem, FALSE, size, fullmatrix, 1, &eigenvalue, NULL) );

      /* check if the eigenvalue is negative */
      if ( eigenvalue < -1 * sdpi->feastol )
      {
         sdpi->infeasible = TRUE;
         SCIPdebugMessage("Detected infeasibility for SDP %d with all fixed variables!\n", sdpi->sdpid);
         break;
      }
   }

   /* free memory */
   BMSfreeBlockMemoryArray(sdpi->blkmem, &fullmatrix, maxsize * maxsize);/*lint !e737*//*lint !e647*/

   /* if we didn't find an SDP-block with negative eigenvalue, the solution is feasible */
   sdpi->infeasible = FALSE;
   SCIPdebugMessage("Unique solution for SDP %d with all fixed variables is feasible!\n", sdpi->sdpid);

   return SCIP_OKAY;
}

/** checks primal and dual Slater condition and outputs result depending on Slater settings in sdpi as well as updating
 *  sdpisolver->primalslater and sdpisolver->dualslater
 */
static
SCIP_RETCODE checkSlaterCondition(
   SCIP_SDPI*            sdpi,               /**< pointer to an SDP-interface structure */
   SCIP_Real             timelimit,          /**< after this many seconds solving will be aborted (currently only implemented for DSDP) */
   clock_t               starttime,          /**< currenttime - starttime will be substracted from the timelimit given to the solver */
   int*                  sdpconstnblocknonz, /**< number of nonzeros for each variable in the constant part, also the i-th entry gives the
                                              *   number of entries  of sdpconst row/col/val [i] (may be NULL if sdpconstnnonz = 0) */
   int**                 sdpconstrow,        /**< pointer to row-indices of constant matrix for each block (may be NULL if sdpconstnnonz = 0) */
   int**                 sdpconstcol,        /**< pointer to column-indices of constant matrix for each block (may be NULL if sdpconstnnonz = 0) */
   SCIP_Real**           sdpconstval,        /**< pointer to values of constant matrix for each block (may be NULL if sdpconstnnonz = 0) */
   int**                 indchanges,         /**< index changes for each variable in each block; variable v is removed in block b if indchanges[b][v] = -1,
                                                  otherwise it gives the number of removed variables with smaller indices (may be NULL if sdpi->nsdpblocks = 0)*/
   int*                  nremovedinds,       /**< number of removed variables for each block (may be NULL if sdpi->nsdpblocks = 0) */
   SCIP_Real*            lplhsafterfix,      /**< left-hand sides of LP-constraints after fixing variables (may be NULL if nactivelpcons = 0) */
   SCIP_Real*            lprhsafterfix,      /**< right-hand sides of LP-constraints after fixing variables (may be NULL if nactivelpcons = 0) */
   int*                  rowsnactivevars,    /**< number of active variables for each LP-constraint (may be NULL if sdpi->nlpcons = 0) */
   int*                  blockindchanges,    /**< index changes for SDP-blocks; blockindchanges[b] = -1 if SDP-block b should be removed
                                              *   (may be NULL if sdpi->nsdpblocks = 0) */
   int                   sdpconstnnonz,      /**< total number of nonzeros in the constant SDP part */
   int                   nactivelpcons,      /**< number of active LP-constraints */
   int                   nremovedblocks,     /**< number of removed SDP-blocks */
   SCIP_Bool             rootnodefailed      /**< if TRUE we will output a message that the root node could not be solved and whether this was due
                                              *   to the Slater condition, otherwise we will print depending on sdpi->slatercheck */
   )
{
   SCIP_Real objval;
   SCIP_Bool origfeas = FALSE;
   SCIP_Bool penaltybound = FALSE;
   int* slaterlprow;
   int* slaterlpcol;
   SCIP_Real* slaterlpval;
   SCIP_Real* slaterlplhs;
   SCIP_Real* slaterlprhs;
   int* slaterrowsnactivevars;
   int nremovedslaterlpinds;
   int i;
   int v;
   int b;
   int slaternactivelpcons;
   SCIP_Real* slaterlb;
   SCIP_Real* slaterub;
   int slaternremovedvarbounds;
   SCIP_Real solvertimelimit;
   clock_t currenttime;

   assert( sdpi != NULL );
   assert( sdpconstnnonz == 0 || sdpconstnblocknonz != NULL );
   assert( sdpconstnnonz == 0 || sdpconstrow != NULL );
   assert( sdpconstnnonz == 0 || sdpconstcol != NULL );
   assert( sdpconstnnonz == 0 || sdpconstval != NULL );
   assert( sdpi->nsdpblocks == 0 || indchanges != NULL );
   assert( sdpi->nsdpblocks == 0 || nremovedinds != NULL );
   assert( nactivelpcons == 0 || lplhsafterfix != NULL );
   assert( nactivelpcons == 0 || lprhsafterfix != NULL );
   assert( sdpi->nlpcons == 0 || rowsnactivevars != NULL );
   assert( sdpi->nsdpblocks == 0 || blockindchanges != NULL );

   /* first check the Slater condition for the dual problem */

   /* compute the timit limit to set for the solver */
   solvertimelimit = timelimit;
   if ( ! SCIPsdpiIsInfinity(sdpi, solvertimelimit) )
   {
      currenttime = clock();
      solvertimelimit -= (SCIP_Real)(currenttime - starttime) / (SCIP_Real) CLOCKS_PER_SEC;/*lint !e620*/
   }

   /* we solve the problem with a slack variable times identity added to the constraints and trying to minimize this slack variable r, if we are
    * still feasible for r > feastol, then we have an interior point with smallest eigenvalue > feastol, otherwise the Slater condition is harmed */
   SCIP_CALL( SCIPsdpiSolverLoadAndSolveWithPenalty(sdpi->sdpisolver, 1.0, FALSE, FALSE, sdpi->nvars, sdpi->obj, sdpi->lb, sdpi->ub,
         sdpi->nsdpblocks, sdpi->sdpblocksizes, sdpi->sdpnblockvars, sdpconstnnonz,
         sdpconstnblocknonz, sdpconstrow, sdpconstcol, sdpconstval,
         sdpi->sdpnnonz, sdpi->sdpnblockvarnonz, sdpi->sdpvar, sdpi->sdprow, sdpi->sdpcol,
         sdpi->sdpval, indchanges, nremovedinds, blockindchanges, nremovedblocks, nactivelpcons, sdpi->nlpcons, lplhsafterfix, lprhsafterfix,
         rowsnactivevars, sdpi->lpnnonz, sdpi->lprow, sdpi->lpcol, sdpi->lpval, NULL, SCIP_SDPSOLVERSETTING_UNSOLVED, solvertimelimit,
         &origfeas, &penaltybound) );

   if ( ! SCIPsdpiSolverIsOptimal(sdpi->sdpisolver) && ! SCIPsdpiSolverIsDualUnbounded(sdpi->sdpisolver) && ! SCIPsdpiSolverIsDualInfeasible(sdpi->sdpisolver) )
   {
      if ( rootnodefailed )
      {
         SCIPmessagePrintInfo(sdpi->messagehdlr, "Aborting due to failing to solve the root node relaxation, Slater condition for the dual problem could "
               "not be checked, ");
      }
      else if ( sdpi->slatercheck == 2 )
         SCIPmessagePrintInfo(sdpi->messagehdlr, "Unable to check Slater condition for dual problem.\n");
      sdpi->dualslater = SCIP_SDPSLATER_NOINFO;
   }
   else
   {
      if ( SCIPsdpiSolverIsDualUnbounded(sdpi->sdpisolver) )
      {
         if ( rootnodefailed )
            SCIPmessagePrintInfo(sdpi->messagehdlr, "Aborting due to failing to solve the root node relaxation, Slater condition for the dual problem holds "
                  "as smallest eigenvalue maximization problem is unbounded, ");
         else
         {
            SCIPdebugMessage("Slater condition for dual problem for SDP %d fullfilled, smallest eigenvalue maximization problem unbounded.\n", sdpi->sdpid);/*lint !e687*/
         }
         sdpi->dualslater = SCIP_SDPSLATER_HOLDS;
      }
      else if ( SCIPsdpiSolverIsDualInfeasible(sdpi->sdpisolver) )
      {
         if ( rootnodefailed )
         {
            SCIPmessagePrintInfo(sdpi->messagehdlr, "Aborting due to failing to solve the root node relaxation, Slater condition for the dual problem "
                  "not fullfilled as problem is infeasible, ");
         }
         else if ( sdpi->slatercheck == 2 )
            SCIPmessagePrintInfo(sdpi->messagehdlr, "Slater condition for dual problem for SDP %d not fullfilled, problem infeasible.\n", sdpi->sdpid);
         sdpi->dualslater = SCIP_SDPSLATER_NOT;
      }
      else
      {
         SCIP_CALL( SCIPsdpiSolverGetObjval(sdpi->sdpisolver, &objval) );

         if ( objval < - sdpi->feastol )
         {
            if ( rootnodefailed )
            {
               SCIPmessagePrintInfo(sdpi->messagehdlr, "Aborting due to failing to solve the root node relaxation, Slater condition for the dual problem holds"
                     "with smallest eigenvalue %f, ", -1.0 * objval);
            }
            else
               SCIPdebugMessage("Slater condition for SDP %d is fullfilled for dual problem with smallest eigenvalue %f.\n", sdpi->sdpid, -1.0 * objval);/*lint !e687*/
            sdpi->dualslater = SCIP_SDPSLATER_HOLDS;
         }
         else if ( objval < sdpi->feastol )
         {
            if ( rootnodefailed )
            {
               SCIPmessagePrintInfo(sdpi->messagehdlr, "Aborting due to failing to solve the root node relaxation, Slater condition for the dual problem "
                     "not fullfilled with smallest eigenvalue %f, ", -1.0 * objval);
            }
            else if ( sdpi->slatercheck == 2 )
            {
               SCIPmessagePrintInfo(sdpi->messagehdlr, "Slater condition for SDP %d not fullfilled for dual problem as smallest eigenvalue was %f, expect numerical trouble.\n",
                     sdpi->sdpid, -1.0 * objval);
            }
            sdpi->dualslater = SCIP_SDPSLATER_NOT;
         }
         else
         {
            if ( sdpi->slatercheck == 2 )
            {
               SCIPmessagePrintInfo(sdpi->messagehdlr, "Slater condition for SDP %d not fullfilled for dual problem as smallest eigenvalue was %f, problem is infeasible.\n",
                     sdpi->sdpid, -1.0 * objval);
            }
            sdpi->dualslater = SCIP_SDPSLATER_INF;
         }
      }
   }

   /* check the Slater condition also for the primal problem */

   /* As we do not want to give equality constraints to the solver by reformulating the primal problem as a dual problem, we instead
    * solve the primal dual pair
    *
    * (P) max (0 0) * Y' s.t. (A_i            0    ) * Y' = c_i forall i, Y' psd
    *         (0 1)           ( 0  sum_j [(A_i)_jj])
    *
    * (D) min sum_i [c_i x_i] s.t. sum_i [A_i x_i] psd, sum_i[(sum_j [(A_i)_jj]) x_i] >= 1
    *
    * where we also set all finite lhs/rhs of all lp-constraints and varbounds to zero.
    * If the objective is strictly positive, than we now that there exists some r > 0 such that
    * Y is psd and Y+rI is feasible for the equality constraints in our original primal problem,
    * so Y+rI is also feasible for the original primal problem and is strictly positive definite
    * so the primal Slater condition holds
    */

   /* allocate the LP-arrays, as we have to add the additional LP-constraint, because we want to add extra entries, we cannot use BMSduplicate... */
   BMS_CALL( BMSallocBlockMemoryArray(sdpi->blkmem, &slaterlprow, sdpi->lpnnonz + sdpi->nvars) );/*lint !e776*/
   BMS_CALL( BMSallocBlockMemoryArray(sdpi->blkmem, &slaterlpcol, sdpi->lpnnonz + sdpi->nvars) );/*lint !e776*/
   BMS_CALL( BMSallocBlockMemoryArray(sdpi->blkmem, &slaterlpval, sdpi->lpnnonz + sdpi->nvars) );/*lint !e776*/

   /* copy all old LP-entries */
   for (i = 0; i < sdpi->lpnnonz; i++)
   {
      slaterlprow[i] = sdpi->lprow[i];
      slaterlpcol[i] = sdpi->lpcol[i];
      slaterlpval[i] = sdpi->lpval[i];
   }

   /* add the new entries sum_j [(A_i)_jj], for this we have to iterate over the whole sdp-matrices (for all blocks), adding all diagonal entries */
   for (v = 0; v < sdpi->nvars; v++)
   {
      slaterlprow[sdpi->lpnnonz + v] = sdpi->nlpcons;/*lint !e679*/
      slaterlpcol[sdpi->lpnnonz + v] = v;/*lint !e679*/
      slaterlpval[sdpi->lpnnonz + v] = 0.0;/*lint !e679*/
   }
   for (b = 0; b < sdpi->nsdpblocks; b++)
   {
      for (v = 0; v < sdpi->sdpnblockvars[b]; v++)
      {
         for (i = 0; i < sdpi->sdpnblockvarnonz[b][v]; i++)
         {
            if ( sdpi->sdprow[b][v][i] == sdpi->sdpcol[b][v][i] ) /* it is a diagonal entry */
               slaterlpval[sdpi->lpnnonz + sdpi->sdpvar[b][v]] += sdpi->sdpval[b][v][i];/*lint !e679*/
         }
      }
   }

   /* iterate over all added LP-entries and remove all zeros (by shifting further variables) */
   nremovedslaterlpinds = 0;
   for (v = 0; v < sdpi->nvars; v++)
   {
      if ( REALABS(slaterlpval[sdpi->lpnnonz + v]) <= sdpi->epsilon )/*lint !e679*/
         nremovedslaterlpinds++;
      else
      {
         /* shift the entries */
         slaterlprow[sdpi->lpnnonz + v - nremovedslaterlpinds] = slaterlprow[sdpi->lpnnonz + v];/*lint !e679*/
         slaterlpcol[sdpi->lpnnonz + v - nremovedslaterlpinds] = slaterlpcol[sdpi->lpnnonz + v];/*lint !e679*/
         slaterlpval[sdpi->lpnnonz + v - nremovedslaterlpinds] = slaterlpval[sdpi->lpnnonz + v];/*lint !e679*/
      }
   }

   /* allocate memory for l/r-hs */
   BMS_CALL( BMSallocBlockMemoryArray(sdpi->blkmem, &slaterlplhs, nactivelpcons + 1) );/*lint !e776*/
   BMS_CALL( BMSallocBlockMemoryArray(sdpi->blkmem, &slaterlprhs, nactivelpcons + 1) );/*lint !e776*/

   /* set the old entries to zero (if existing), as A_0 (including the LP-part) is removed because of the changed primal objective */
   for (i = 0; i < nactivelpcons; i++)
   {
      if ( SCIPsdpiSolverIsInfinity(sdpi->sdpisolver, lplhsafterfix[i]) )
         slaterlplhs[i] = lplhsafterfix[i];
      else
         slaterlplhs[i] = 0.0;

      if ( SCIPsdpiSolverIsInfinity(sdpi->sdpisolver, lprhsafterfix[i]) )
         slaterlprhs[i] = lprhsafterfix[i];
      else
         slaterlprhs[i] = 0.0;
   }

   /* add the new ones */
   slaterlplhs[nactivelpcons] = 1.0;
   slaterlprhs[nactivelpcons] = SCIPsdpiSolverInfinity(sdpi->sdpisolver);

   /* allocate memory for rowsnactivevars to update it for the added row */
   BMS_CALL( BMSallocBlockMemoryArray(sdpi->blkmem, &slaterrowsnactivevars, sdpi->nlpcons + 1) );/*lint !e776*/

   /* copy the old entries */
   for (i = 0; i < sdpi->nlpcons; i++)
      slaterrowsnactivevars[i] = rowsnactivevars[i];

   /* add the new entry (this equals the number of active variables) */
   slaterrowsnactivevars[sdpi->nlpcons] = 0;
   for (v = 0; v < sdpi->nvars; v++)
   {
      if ( ! (isFixed(sdpi, v)) )
         slaterrowsnactivevars[sdpi->nlpcons]++;
   }

   slaternactivelpcons = (slaterrowsnactivevars[sdpi->nlpcons] > 1) ? nactivelpcons + 1 : nactivelpcons;

   /* copy the varbound arrays to change all finite varbounds to zero */
   DUPLICATE_ARRAY_NULL(sdpi->blkmem, &slaterlb, sdpi->lb, sdpi->nvars);
   DUPLICATE_ARRAY_NULL(sdpi->blkmem, &slaterub, sdpi->ub, sdpi->nvars);

   /* set all finite varbounds to zero */
   slaternremovedvarbounds = 0;
   for (v = 0; v < sdpi->nvars; v++)
   {
      if ( slaterlb[v] > -1 * SCIPsdpiSolverInfinity(sdpi->sdpisolver) )
      {
         slaterlb[v] = 0.0;
         slaternremovedvarbounds++;
      }
      if ( slaterub[v] < SCIPsdpiSolverInfinity(sdpi->sdpisolver) )
      {
         slaterub[v] = 0.0;
         slaternremovedvarbounds++;
      }
   }

   /* if all variables have finite upper and lower bounds these add variables to every constraint of the
    * primal problem that allow us to make the problem feasible for every primal matrix X, so the primal
    * Slater condition holds */
   if ( slaternremovedvarbounds == 2 * sdpi->nvars )
   {
      if ( rootnodefailed )
      {
         SCIPmessagePrintInfo(sdpi->messagehdlr, "Slater condition for primal problem holds since all variables have finite upper and lower bounds \n");
      }
      else
         SCIPdebugMessage("Slater condition for primal problem for SDP %d fullfilled as all variables have finite upper and lower bounds \n", sdpi->sdpid);/*lint !e687*/
      sdpi->primalslater = SCIP_SDPSLATER_HOLDS;
   }
   else
   {
      /* compute the timit limit to set for the solver */
      currenttime = clock();
      solvertimelimit = timelimit - ((SCIP_Real)(currenttime - starttime) / (SCIP_Real) CLOCKS_PER_SEC); /*lint !e620*/

      /* solve the problem to check Slater condition for primal of original problem */
      SCIP_CALL( SCIPsdpiSolverLoadAndSolve(sdpi->sdpisolver, sdpi->nvars, sdpi->obj, slaterlb, slaterub,
            sdpi->nsdpblocks, sdpi->sdpblocksizes, sdpi->sdpnblockvars, 0, NULL, NULL, NULL, NULL,
            sdpi->sdpnnonz, sdpi->sdpnblockvarnonz, sdpi->sdpvar, sdpi->sdprow, sdpi->sdpcol,
            sdpi->sdpval, indchanges, nremovedinds, blockindchanges, nremovedblocks, slaternactivelpcons, sdpi->nlpcons + 1, slaterlplhs, slaterlprhs,
            slaterrowsnactivevars, sdpi->lpnnonz + sdpi->nvars - nremovedslaterlpinds, slaterlprow, slaterlpcol, slaterlpval, NULL,
            SCIP_SDPSOLVERSETTING_UNSOLVED, solvertimelimit) );

      if ( ! SCIPsdpiSolverIsOptimal(sdpi->sdpisolver) && ! SCIPsdpiSolverIsDualUnbounded(sdpi->sdpisolver) && ! SCIPsdpiSolverIsPrimalUnbounded(sdpi->sdpisolver) )
      {
         if ( rootnodefailed )
         {
            SCIPmessagePrintInfo(sdpi->messagehdlr, "unable to check Slater condition for primal problem \n");
         }
         else if ( sdpi->slatercheck == 2 )
            SCIPmessagePrintInfo(sdpi->messagehdlr, "Unable to check Slater condition for primal problem, could not solve auxilliary problem.\n");
         sdpi->primalslater = SCIP_SDPSLATER_NOINFO;
      }
      else if ( SCIPsdpiSolverIsDualUnbounded(sdpi->sdpisolver) )
      {
         if ( rootnodefailed )
         {
            SCIPmessagePrintInfo(sdpi->messagehdlr, " primal Slater condition shows infeasibility \n");
         }
         else if ( sdpi->slatercheck == 2 )
         {
            SCIPmessagePrintInfo(sdpi->messagehdlr, "Slater condition for primal problem for SDP %d not fullfilled "
                  "smallest eigenvalue has to be negative, so primal problem is infeasible (if the dual slater condition holds,"
                  "this means, that the original (dual) problem is unbounded.\n",sdpi->sdpid);
         }
         sdpi->primalslater = SCIP_SDPSLATER_NOT;
      }
      else if ( SCIPsdpiSolverIsPrimalUnbounded(sdpi->sdpisolver) )
      {
         if ( rootnodefailed )
         {
            SCIPmessagePrintInfo(sdpi->messagehdlr, "Slater condition for primal problems holds sunce smallest eigenvalue maximization problem"
                  "is unbounded \n");
         }
         else
            SCIPdebugMessage("Slater condition for primal problem for SDP %d fullfilled, smallest eigenvalue maximization problem unbounded \n", sdpi->sdpid);/*lint !e687*/
         sdpi->primalslater = SCIP_SDPSLATER_HOLDS;
      }
      else
      {
         SCIP_CALL( SCIPsdpiSolverGetObjval(sdpi->sdpisolver, &objval) );

         if ( objval > - sdpi->feastol)
         {
            if ( rootnodefailed )
            {
               SCIPmessagePrintInfo(sdpi->messagehdlr, "Slater condition for primal problem not fullfilled with smallest eigenvalue %f \n", -1.0 * objval);
            }
            else if ( sdpi->slatercheck == 2 )
            {
               SCIPmessagePrintInfo(sdpi->messagehdlr, "Slater condition for primal problem for SDP %d not fullfilled "
                        "as smallest eigenvalue was %f, expect numerical trouble or infeasible problem.\n",sdpi->sdpid, -1.0 * objval);
            }
            sdpi->primalslater = SCIP_SDPSLATER_NOT;
         }
         else
         {
            if ( rootnodefailed )
            {
               SCIPmessagePrintInfo(sdpi->messagehdlr, "Slater condition for primal problem fullfilled with smallest eigenvalue %f \n", -1.0 * objval);
            }
            else
               SCIPdebugMessage("Slater condition for primal problem of SDP %d is fullfilled with smallest eigenvalue %f.\n", sdpi->sdpid, -1.0 * objval);/*lint !e687*/
            sdpi->primalslater = SCIP_SDPSLATER_HOLDS;
         }
      }
   }

   /* free all memory */
   BMSfreeBlockMemoryArray(sdpi->blkmem, &slaterub, sdpi->nvars);/*lint !e737*/
   BMSfreeBlockMemoryArray(sdpi->blkmem, &slaterlb, sdpi->nvars);/*lint !e737*/
   BMSfreeBlockMemoryArray(sdpi->blkmem, &slaterrowsnactivevars, sdpi->nlpcons + 1);/*lint !e737*//*lint !e776*/
   BMSfreeBlockMemoryArray(sdpi->blkmem, &slaterlprhs, nactivelpcons + 1);/*lint !e737*//*lint !e776*/
   BMSfreeBlockMemoryArray(sdpi->blkmem, &slaterlplhs, nactivelpcons + 1);/*lint !e737*//*lint !e776*/
   BMSfreeBlockMemoryArray(sdpi->blkmem, &slaterlpval, sdpi->lpnnonz + sdpi->nvars);/*lint !e737*//*lint !e776*/
   BMSfreeBlockMemoryArray(sdpi->blkmem, &slaterlpcol, sdpi->lpnnonz + sdpi->nvars);/*lint !e737*//*lint !e776*/
   BMSfreeBlockMemoryArray(sdpi->blkmem, &slaterlprow, sdpi->lpnnonz + sdpi->nvars);/*lint !e737*//*lint !e776*/

   return SCIP_OKAY;
}

/*
 * Miscellaneous Methods
 */

/**@name Miscellaneous Methods */
/**@{ */


/** gets name and potentially version of SDP-solver */
const char* SCIPsdpiGetSolverName(
   void
   )
{
   return SCIPsdpiSolverGetSolverName();
}

/** gets description of SDP-solver (developer, webpage, ...) */
const char* SCIPsdpiGetSolverDesc(
   void
   )
{
   return SCIPsdpiSolverGetSolverDesc();
}

/** gets pointer for SDP-solver - use only with great care
 *
 *  The behavior of this function depends on the solver and its use is
 *  therefore only recommended if you really know what you are
 *  doing. In general, it returns a pointer to the SDP-solver object.
 */
void* SCIPsdpiGetSolverPointer(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   )
{
   return SCIPsdpiSolverGetSolverPointer(sdpi->sdpisolver);
}

/** gets default feasibility tolerance for SDP-solver in SCIP-SDP */
SCIP_Real SCIPsdpiGetDefaultSdpiSolverFeastol(
   void
   )
{
   return SCIPsdpiSolverGetDefaultSdpiSolverFeastol();
}

/** gets default number of increases of penalty parameter for SDP-solver in SCIP-SDP */
int SCIPsdpiGetDefaultSdpiSolverNpenaltyIncreases(
   void
   )
{
   return SCIPsdpiSolverGetDefaultSdpiSolverNpenaltyIncreases();
}

/**@} */


/*
 * SDPI Creation and Destruction Methods
 */

/**@name SDPI Creation and Destruction Methods */
/**@{ */

/** creates an sdpi object */
SCIP_RETCODE SCIPsdpiCreate(
   SCIP_SDPI**           sdpi,               /**< pointer to an SDP-interface structure */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler to use for printing messages, or NULL */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   BMS_BUFMEM*           bufmem              /**< buffer memory */
   )
{
   assert ( sdpi != NULL );
   assert ( blkmem != NULL );

   SCIPdebugMessage("Calling SCIPsdpiCreate\n");

   BMS_CALL( BMSallocBlockMemory(blkmem, sdpi) );

   SCIP_CALL( SCIPsdpiSolverCreate(&((*sdpi)->sdpisolver), messagehdlr, blkmem, bufmem) );

   (*sdpi)->messagehdlr = messagehdlr;
   (*sdpi)->blkmem = blkmem;
   (*sdpi)->bufmem = bufmem;
   (*sdpi)->sdpid = 1;
   (*sdpi)->niterations = 0;
   (*sdpi)->nsdpcalls = 0;
   (*sdpi)->nvars = 0;
   (*sdpi)->nsdpblocks = 0;
   (*sdpi)->sdpconstnnonz = 0;
   (*sdpi)->sdpnnonz = 0;
   (*sdpi)->nlpcons = 0;
   (*sdpi)->lpnnonz = 0;
   (*sdpi)->slatercheck = 0;
   (*sdpi)->solved = FALSE;
   (*sdpi)->penalty = FALSE;
   (*sdpi)->infeasible = FALSE;
   (*sdpi)->allfixed = FALSE;

   (*sdpi)->obj = NULL;
   (*sdpi)->lb = NULL;
   (*sdpi)->ub = NULL;
   (*sdpi)->sdpblocksizes = NULL;
   (*sdpi)->sdpnblockvars = NULL;
   (*sdpi)->sdpconstnblocknonz = NULL;
   (*sdpi)->sdpconstrow = NULL;
   (*sdpi)->sdpconstcol = NULL;
   (*sdpi)->sdpconstval = NULL;
   (*sdpi)->sdpnblockvarnonz = NULL;
   (*sdpi)->sdpvar = NULL;
   (*sdpi)->sdprow = NULL;
   (*sdpi)->sdpcol = NULL;
   (*sdpi)->sdpval = NULL;
   (*sdpi)->lplhs = NULL;
   (*sdpi)->lprhs = NULL;
   (*sdpi)->lprow = NULL;
   (*sdpi)->lpcol = NULL;
   (*sdpi)->lpval = NULL;

   (*sdpi)->epsilon = DEFAULT_EPSILON;
   (*sdpi)->gaptol = DEFAULT_SDPSOLVERGAPTOL;
   (*sdpi)->feastol = DEFAULT_FEASTOL;
   (*sdpi)->penaltyparam = DEFAULT_PENALTYPARAM;
   (*sdpi)->maxpenaltyparam = DEFAULT_MAXPENALTYPARAM;
   (*sdpi)->npenaltyincr = DEFAULT_NPENALTYINCR;
   (*sdpi)->bestbound = -SCIPsdpiSolverInfinity((*sdpi)->sdpisolver);
   (*sdpi)->primalslater = SCIP_SDPSLATER_NOINFO;
   (*sdpi)->dualslater = SCIP_SDPSLATER_NOINFO;

   return SCIP_OKAY;
}

/** deletes an sdpi object */
SCIP_RETCODE SCIPsdpiFree(
   SCIP_SDPI**           sdpi                /**< pointer to an SDP-interface structure */
   )
{
   int i;
   int j;

   SCIPdebugMessage("Calling SCIPsdpiFree \n");
   assert ( sdpi != NULL );
   assert ( *sdpi != NULL );

   /* free the LP part */
   BMSfreeBlockMemoryArrayNull((*sdpi)->blkmem, &((*sdpi)->lpval), (*sdpi)->lpnnonz);
   BMSfreeBlockMemoryArrayNull((*sdpi)->blkmem, &((*sdpi)->lpcol), (*sdpi)->lpnnonz);
   BMSfreeBlockMemoryArrayNull((*sdpi)->blkmem, &((*sdpi)->lprow), (*sdpi)->lpnnonz);
   BMSfreeBlockMemoryArrayNull((*sdpi)->blkmem, &((*sdpi)->lprhs), (*sdpi)->nlpcons);
   BMSfreeBlockMemoryArrayNull((*sdpi)->blkmem, &((*sdpi)->lplhs), (*sdpi)->nlpcons);

   /* free the individual nonzeros */
   for (i = 0; i < (*sdpi)->nsdpblocks; i++)
   {
      for (j = 0; j < (*sdpi)->sdpnblockvars[i]; j++)
      {
         BMSfreeBlockMemoryArrayNull((*sdpi)->blkmem, &((*sdpi)->sdpval[i][j]), (*sdpi)->sdpnblockvarnonz[i][j]);
         BMSfreeBlockMemoryArrayNull((*sdpi)->blkmem, &((*sdpi)->sdprow[i][j]), (*sdpi)->sdpnblockvarnonz[i][j]);
         BMSfreeBlockMemoryArrayNull((*sdpi)->blkmem, &((*sdpi)->sdpcol[i][j]), (*sdpi)->sdpnblockvarnonz[i][j]);
      }
      BMSfreeBlockMemoryArrayNull((*sdpi)->blkmem, &((*sdpi)->sdpval[i]), (*sdpi)->sdpnblockvars[i]);
      BMSfreeBlockMemoryArrayNull((*sdpi)->blkmem, &((*sdpi)->sdprow[i]), (*sdpi)->sdpnblockvars[i]);
      BMSfreeBlockMemoryArrayNull((*sdpi)->blkmem, &((*sdpi)->sdpcol[i]), (*sdpi)->sdpnblockvars[i]);
      BMSfreeBlockMemoryArrayNull((*sdpi)->blkmem, &((*sdpi)->sdpvar[i]), (*sdpi)->sdpnblockvars[i]);
      BMSfreeBlockMemoryArrayNull((*sdpi)->blkmem, &((*sdpi)->sdpnblockvarnonz[i]), (*sdpi)->sdpnblockvars[i]);
      BMSfreeBlockMemoryArrayNull((*sdpi)->blkmem, &((*sdpi)->sdpconstval[i]), (*sdpi)->sdpconstnblocknonz[i]);
      BMSfreeBlockMemoryArrayNull((*sdpi)->blkmem, &((*sdpi)->sdpconstrow[i]), (*sdpi)->sdpconstnblocknonz[i]);
      BMSfreeBlockMemoryArrayNull((*sdpi)->blkmem, &((*sdpi)->sdpconstcol[i]), (*sdpi)->sdpconstnblocknonz[i]);
   }

   /* free the rest */
   BMSfreeBlockMemoryArrayNull((*sdpi)->blkmem, &((*sdpi)->sdpnblockvarnonz), (*sdpi)->nsdpblocks);
   BMSfreeBlockMemoryArrayNull((*sdpi)->blkmem, &((*sdpi)->sdpconstnblocknonz), (*sdpi)->nsdpblocks);
   BMSfreeBlockMemoryArrayNull((*sdpi)->blkmem, &((*sdpi)->sdpval), (*sdpi)->nsdpblocks);
   BMSfreeBlockMemoryArrayNull((*sdpi)->blkmem, &((*sdpi)->sdpcol), (*sdpi)->nsdpblocks);
   BMSfreeBlockMemoryArrayNull((*sdpi)->blkmem, &((*sdpi)->sdprow), (*sdpi)->nsdpblocks);
   BMSfreeBlockMemoryArrayNull((*sdpi)->blkmem, &((*sdpi)->sdpvar), (*sdpi)->nsdpblocks);
   BMSfreeBlockMemoryArrayNull((*sdpi)->blkmem, &((*sdpi)->sdpconstval), (*sdpi)->nsdpblocks);
   BMSfreeBlockMemoryArrayNull((*sdpi)->blkmem, &((*sdpi)->sdpconstcol), (*sdpi)->nsdpblocks);
   BMSfreeBlockMemoryArrayNull((*sdpi)->blkmem, &((*sdpi)->sdpconstrow), (*sdpi)->nsdpblocks);
   BMSfreeBlockMemoryArrayNull((*sdpi)->blkmem, &((*sdpi)->sdpnblockvars), (*sdpi)->nsdpblocks);
   BMSfreeBlockMemoryArrayNull((*sdpi)->blkmem, &((*sdpi)->sdpblocksizes), (*sdpi)->nsdpblocks);
   BMSfreeBlockMemoryArrayNull((*sdpi)->blkmem, &((*sdpi)->ub), (*sdpi)->nvars);/*lint !e737*/
   BMSfreeBlockMemoryArrayNull((*sdpi)->blkmem, &((*sdpi)->lb), (*sdpi)->nvars);/*lint !e737*/
   BMSfreeBlockMemoryArrayNull((*sdpi)->blkmem, &((*sdpi)->obj), (*sdpi)->nvars);/*lint !e737*/

   /* free the solver */
   SCIP_CALL( SCIPsdpiSolverFree(&((*sdpi)->sdpisolver)) );

   BMSfreeBlockMemory((*sdpi)->blkmem, sdpi);

   return SCIP_OKAY;
}

/** cloning method of the general SDP-Interface
 *
 *  @note The solver specific interface is created anew and not copied.
 */
SCIP_RETCODE SCIPsdpiClone(
   SCIP_SDPI*            oldsdpi,            /**< pointer to the SDP-interface structure that should be cloned */
   SCIP_SDPI*            newsdpi             /**< pointer to an SDP-interface structure to clone into */
   )
{
   BMS_BLKMEM* blkmem;
   int nvars;
   int nsdpblocks;
   int lpnnonz;
   int b;
   int v;

   assert( oldsdpi != NULL );

   SCIPdebugMessage("Cloning SDPI %d\n", oldsdpi->sdpid);

   /* general data */
   blkmem = oldsdpi->blkmem;
   nvars = oldsdpi->nvars;
   nsdpblocks = oldsdpi->nsdpblocks;
   lpnnonz = oldsdpi->lpnnonz;

   BMS_CALL( BMSallocBlockMemory(blkmem, &newsdpi) );

   SCIP_CALL( SCIPsdpiSolverCreate(&(newsdpi->sdpisolver), oldsdpi->messagehdlr, oldsdpi->blkmem, oldsdpi->bufmem) ); /* create new SDP-Solver Interface */

   newsdpi->messagehdlr = oldsdpi->messagehdlr;
   newsdpi->blkmem = blkmem;
   newsdpi->nvars = nvars;

   BMS_CALL( BMSduplicateBlockMemoryArray(blkmem, &(newsdpi->obj), oldsdpi->obj, nvars) );
   BMS_CALL( BMSduplicateBlockMemoryArray(blkmem, &(newsdpi->lb), oldsdpi->lb, nvars) );
   BMS_CALL( BMSduplicateBlockMemoryArray(blkmem, &(newsdpi->ub), oldsdpi->ub, nvars) );

   newsdpi->nsdpblocks = nsdpblocks;

   BMS_CALL( BMSduplicateBlockMemoryArray(blkmem, &(newsdpi->sdpblocksizes), oldsdpi->sdpblocksizes, nsdpblocks) );
   BMS_CALL( BMSduplicateBlockMemoryArray(blkmem, &(newsdpi->sdpnblockvars), oldsdpi->sdpnblockvars, nsdpblocks) );

   /* constant SDP data */
   newsdpi->sdpconstnnonz = oldsdpi->sdpconstnnonz;

   BMS_CALL( BMSduplicateBlockMemoryArray(blkmem, &(newsdpi->sdpconstnblocknonz), oldsdpi->sdpconstnblocknonz, nsdpblocks) );
   BMS_CALL( BMSallocBlockMemoryArray(blkmem, &(newsdpi->sdpconstrow), nsdpblocks) );
   BMS_CALL( BMSallocBlockMemoryArray(blkmem, &(newsdpi->sdpconstcol), nsdpblocks) );
   BMS_CALL( BMSallocBlockMemoryArray(blkmem, &(newsdpi->sdpconstval), nsdpblocks) );

   for (b = 0; b < nsdpblocks; b++)
   {
      BMS_CALL( BMSduplicateBlockMemoryArray(blkmem, &(newsdpi->sdpconstrow[b]), oldsdpi->sdpconstrow[b], oldsdpi->sdpconstnblocknonz[b]) );
      BMS_CALL( BMSduplicateBlockMemoryArray(blkmem, &(newsdpi->sdpconstcol[b]), oldsdpi->sdpconstcol[b], oldsdpi->sdpconstnblocknonz[b]) );
      BMS_CALL( BMSduplicateBlockMemoryArray(blkmem, &(newsdpi->sdpconstval[b]), oldsdpi->sdpconstval[b], oldsdpi->sdpconstnblocknonz[b]) );
   }

   /* SDP data */
   newsdpi->sdpnnonz = oldsdpi->sdpnnonz;

   BMS_CALL( BMSallocBlockMemoryArray(blkmem, &(newsdpi->sdpnblockvarnonz), nsdpblocks) );
   BMS_CALL( BMSallocBlockMemoryArray(blkmem, &(newsdpi->sdpvar), nsdpblocks) );
   BMS_CALL( BMSallocBlockMemoryArray(blkmem, &(newsdpi->sdprow), nsdpblocks) );
   BMS_CALL( BMSallocBlockMemoryArray(blkmem, &(newsdpi->sdpcol), nsdpblocks) );
   BMS_CALL( BMSallocBlockMemoryArray(blkmem, &(newsdpi->sdpval), nsdpblocks) );

   for (b = 0; b < nsdpblocks; b++)
   {
      BMS_CALL( BMSduplicateBlockMemoryArray(blkmem, &(newsdpi->sdpnblockvarnonz[b]), oldsdpi->sdpnblockvarnonz[b], oldsdpi->sdpnblockvars[b]) );
      BMS_CALL( BMSduplicateBlockMemoryArray(blkmem, &(newsdpi->sdpvar[b]), oldsdpi->sdpvar[b], oldsdpi->sdpnblockvars[b]) );

      BMS_CALL( BMSallocBlockMemoryArray(blkmem, &(newsdpi->sdprow[b]), oldsdpi->sdpnblockvars[b]) );
      BMS_CALL( BMSallocBlockMemoryArray(blkmem, &(newsdpi->sdpcol[b]), oldsdpi->sdpnblockvars[b]) );
      BMS_CALL( BMSallocBlockMemoryArray(blkmem, &(newsdpi->sdpval[b]), oldsdpi->sdpnblockvars[b]) );

      for (v = 0; v < oldsdpi->sdpnblockvars[b]; v++)
      {
         BMS_CALL( BMSduplicateBlockMemoryArray(blkmem, &(newsdpi->sdprow[b][v]), oldsdpi->sdprow[b][v], oldsdpi->sdpnblockvarnonz[b][v]) );
         BMS_CALL( BMSduplicateBlockMemoryArray(blkmem, &(newsdpi->sdpcol[b][v]), oldsdpi->sdpcol[b][v], oldsdpi->sdpnblockvarnonz[b][v]) );
         BMS_CALL( BMSduplicateBlockMemoryArray(blkmem, &(newsdpi->sdpval[b][v]), oldsdpi->sdpval[b][v], oldsdpi->sdpnblockvarnonz[b][v]) );
      }
   }

   /* LP data */
   newsdpi->nlpcons = oldsdpi->nlpcons;

   BMS_CALL( BMSduplicateBlockMemoryArray(blkmem, &(newsdpi->lplhs), oldsdpi->lplhs, oldsdpi->nlpcons) );
   BMS_CALL( BMSduplicateBlockMemoryArray(blkmem, &(newsdpi->lprhs), oldsdpi->lprhs, oldsdpi->nlpcons) );

   newsdpi->lpnnonz = lpnnonz;

   BMS_CALL( BMSduplicateBlockMemoryArray(blkmem, &(newsdpi->lprow), oldsdpi->lprow, lpnnonz) );
   BMS_CALL( BMSduplicateBlockMemoryArray(blkmem, &(newsdpi->lpcol), oldsdpi->lpcol, lpnnonz) );
   BMS_CALL( BMSduplicateBlockMemoryArray(blkmem, &(newsdpi->lpval), oldsdpi->lpval, lpnnonz) );

   /* other data */
   newsdpi->solved = FALSE; /* as we don't copy the sdpisolver, this needs to be set to false */
   newsdpi->penalty = FALSE; /* all things about SDP-solutions are set to false as well, as we didn't solve the problem */
   newsdpi->infeasible = FALSE;
   newsdpi->allfixed = FALSE;
   newsdpi->sdpid = 1000000 + oldsdpi->sdpid; /* this is only used for debug output, setting it to this value should make it clear, that it is a new sdpi */
   newsdpi->epsilon = oldsdpi->epsilon;
   newsdpi->gaptol = oldsdpi->gaptol;
   newsdpi->feastol = oldsdpi->feastol;

   return SCIP_OKAY;
}

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
                                               *  (may be NULL if sdptnnonz = 0)*/
   int***                sdpcol,             /**< pointer to the column-indices for each block and variable in this block (may be NULL if sdpnnonz = 0)*/
   SCIP_Real***          sdpval,             /**< pointer to the values of the nonzeros for each block and variable in this
                                               *  block (may be NULL if sdpnnonz = 0)*/
   int                   nlpcons,            /**< number of LP-constraints */
   SCIP_Real*            lplhs,              /**< left-hand sides of LP rows (may be NULL if nlpcons = 0) */
   SCIP_Real*            lprhs,              /**< right-hand sides of LP rows (may be NULL if nlpcons = 0) */
   int                   lpnnonz,            /**< number of nonzero elements in the LP-constraint-matrix */
   int*                  lprow,              /**< row-index for each entry in lpval-array (may be NULL if lpnnonz = 0) */
   int*                  lpcol,              /**< column-index for each entry in lpval-array (may be NULL if lpnnonz = 0) */
   SCIP_Real*            lpval               /**< values of LP-constraint matrix entries (may be NULL if lpnnonz = 0) */
   )
{
   int i;
   int v;
   int block;

   SCIPdebugMessage("Calling SCIPsdpiLoadSDP (%d)\n",sdpi->sdpid);

   assert ( sdpi != NULL );
   assert ( nvars > 0 );
   assert ( obj != NULL );
   assert ( lb != NULL );
   assert ( ub != NULL );

#ifdef SCIP_DEBUG
   if (sdpconstnnonz > 0 || sdpnnonz > 0 || nsdpblocks > 0)
   {
      assert ( sdpblocksizes != NULL );
      assert ( sdpnblockvars != NULL );
      assert ( nsdpblocks > 0 );
      assert ( sdpconstnblocknonz != NULL );
      assert ( sdpnblockvarnonz != NULL );

      if (sdpconstnnonz > 0)
      {
         assert ( sdpconstrow != NULL );
         assert ( sdpconstcol != NULL );
         assert ( sdpconstval != NULL );

         for (i = 0; i < nsdpblocks; i++)
         {
            if (sdpconstnblocknonz[i] > 0)
            {
               assert ( sdpconstrow[i] != NULL );
               assert ( sdpconstcol[i] != NULL );
               assert ( sdpconstval[i] != NULL );
            }
         }
      }

      if (sdpnnonz > 0)
      {
         assert ( sdprow != NULL );
         assert ( sdpcol != NULL );
         assert ( sdpval != NULL );

         for ( i = 0; i < nsdpblocks; i++ )
         {
            assert ( sdpcol[i] != NULL );
            assert ( sdprow[i] != NULL );
            assert ( sdpval[i] != NULL );

            for ( v = 0; v < sdpnblockvars[i]; v++)
            {
               if (sdpnblockvarnonz[i][v] > 0)
               {
                  assert ( sdpcol[i][v] != NULL );
                  assert ( sdprow[i][v] != NULL );
                  assert ( sdpval[i][v] != NULL );
               }
            }
         }
      }
   }
#endif

   assert ( nlpcons == 0 || lplhs != NULL );
   assert ( nlpcons == 0 || lprhs != NULL );
   assert ( lpnnonz == 0 || lprow != NULL );
   assert ( lpnnonz == 0 || lpcol != NULL );
   assert ( lpnnonz == 0 || lpval != NULL );

   /* memory allocation */

   /* first free the old arrays */
   for (block = sdpi->nsdpblocks - 1; block >= 0; block--)
   {
      for (v = sdpi->sdpnblockvars[block] - 1; v >= 0; v--)
      {
         BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &(sdpi->sdpval[block][v]), sdpi->sdpnblockvarnonz[block][v]);
         BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &(sdpi->sdprow[block][v]), sdpi->sdpnblockvarnonz[block][v]);
         BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &(sdpi->sdpcol[block][v]), sdpi->sdpnblockvarnonz[block][v]);
      }

      BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &(sdpi->sdpval[block]), sdpi->sdpnblockvars[block]);
      BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &(sdpi->sdprow[block]), sdpi->sdpnblockvars[block]);
      BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &(sdpi->sdpcol[block]), sdpi->sdpnblockvars[block]);
      BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &(sdpi->sdpconstval[block]), sdpi->sdpconstnblocknonz[block]);
      BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &(sdpi->sdpconstrow[block]), sdpi->sdpconstnblocknonz[block]);
      BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &(sdpi->sdpconstcol[block]), sdpi->sdpconstnblocknonz[block]);
      BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &(sdpi->sdpnblockvarnonz[block]), sdpi->sdpnblockvars[block]);
      BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &(sdpi->sdpvar[block]), sdpi->sdpnblockvars[block]);
   }

   BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &(sdpi->ub), sdpi->nvars);
   BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &(sdpi->lb), sdpi->nvars);
   BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &(sdpi->obj), sdpi->nvars);

   BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &(sdpi->sdpblocksizes), sdpi->nsdpblocks);
   BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &(sdpi->sdpnblockvars), sdpi->nsdpblocks);
   BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &(sdpi->sdpconstnblocknonz), sdpi->nsdpblocks);

   /* duplicate some arrays */
   BMS_CALL( BMSduplicateBlockMemoryArray(sdpi->blkmem, &(sdpi->obj), obj, nvars) );
   BMS_CALL( BMSduplicateBlockMemoryArray(sdpi->blkmem, &(sdpi->lb), lb, nvars) );
   BMS_CALL( BMSduplicateBlockMemoryArray(sdpi->blkmem, &(sdpi->ub), ub, nvars) );
   DUPLICATE_ARRAY_NULL(sdpi->blkmem, &(sdpi->sdpblocksizes), sdpblocksizes, nsdpblocks);
   DUPLICATE_ARRAY_NULL(sdpi->blkmem, &(sdpi->sdpnblockvars), sdpnblockvars, nsdpblocks);
   DUPLICATE_ARRAY_NULL(sdpi->blkmem, &(sdpi->sdpconstnblocknonz), sdpconstnblocknonz, nsdpblocks);

   /* allocate memory for the sdp arrays & duplicate them */
   BMS_CALL( BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpnblockvarnonz), sdpi->nsdpblocks, nsdpblocks) );
   BMS_CALL( BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstcol), sdpi->nsdpblocks, nsdpblocks) );
   BMS_CALL( BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstrow), sdpi->nsdpblocks, nsdpblocks) );
   BMS_CALL( BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstval), sdpi->nsdpblocks, nsdpblocks) );
   BMS_CALL( BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpvar), sdpi->nsdpblocks, nsdpblocks) );
   BMS_CALL( BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpcol), sdpi->nsdpblocks, nsdpblocks) );
   BMS_CALL( BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdprow), sdpi->nsdpblocks, nsdpblocks) );
   BMS_CALL( BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpval), sdpi->nsdpblocks, nsdpblocks) );

   for (block = 0; block < nsdpblocks; block++)
   {
      DUPLICATE_ARRAY_NULL(sdpi->blkmem, &(sdpi->sdpnblockvarnonz[block]), sdpnblockvarnonz[block], sdpnblockvars[block]);

      DUPLICATE_ARRAY_NULL(sdpi->blkmem, &(sdpi->sdpconstcol[block]), sdpconstcol[block], sdpconstnblocknonz[block]);
      DUPLICATE_ARRAY_NULL(sdpi->blkmem, &(sdpi->sdpconstrow[block]), sdpconstrow[block], sdpconstnblocknonz[block]);
      DUPLICATE_ARRAY_NULL(sdpi->blkmem, &(sdpi->sdpconstval[block]), sdpconstval[block], sdpconstnblocknonz[block]);

      /* make sure that we have a lower triangular matrix */
      for (i = 0; i < sdpi->sdpconstnblocknonz[block]; ++i)
         ensureLowerTriangular(&(sdpconstrow[block][i]), &(sdpconstcol[block][i]));

      DUPLICATE_ARRAY_NULL(sdpi->blkmem, &(sdpi->sdpvar[block]), sdpvar[block], sdpnblockvars[block]);

      BMS_CALL( BMSallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpcol[block]), sdpnblockvars[block]) );
      BMS_CALL( BMSallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdprow[block]), sdpnblockvars[block]) );
      BMS_CALL( BMSallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpval[block]), sdpnblockvars[block]) );

      for (v = 0; v < sdpi->sdpnblockvars[block]; v++)
      {
         DUPLICATE_ARRAY_NULL(sdpi->blkmem, &(sdpi->sdpcol[block][v]), sdpcol[block][v], sdpnblockvarnonz[block][v]);
         DUPLICATE_ARRAY_NULL(sdpi->blkmem, &(sdpi->sdprow[block][v]), sdprow[block][v], sdpnblockvarnonz[block][v]);
         DUPLICATE_ARRAY_NULL(sdpi->blkmem, &(sdpi->sdpval[block][v]), sdpval[block][v], sdpnblockvarnonz[block][v]);

         /* make sure that we have a lower triangular matrix */
         for (i = 0; i < sdpi->sdpnblockvarnonz[block][v]; ++i)
            ensureLowerTriangular(&(sdprow[block][v][i]), &(sdpcol[block][v][i]));
      }
   }

   /* free old and duplicate new arrays for the LP part */
   BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &(sdpi->lpval), sdpi->lpnnonz);
   BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &(sdpi->lpcol), sdpi->lpnnonz);
   BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &(sdpi->lprow), sdpi->lpnnonz);
   BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &(sdpi->lprhs), sdpi->nlpcons);
   BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &(sdpi->lplhs), sdpi->nlpcons);

   DUPLICATE_ARRAY_NULL(sdpi->blkmem, &(sdpi->lplhs), lplhs, nlpcons);
   DUPLICATE_ARRAY_NULL(sdpi->blkmem, &(sdpi->lprhs), lprhs, nlpcons);
   DUPLICATE_ARRAY_NULL(sdpi->blkmem, &(sdpi->lprow), lprow, lpnnonz);
   DUPLICATE_ARRAY_NULL(sdpi->blkmem, &(sdpi->lpcol), lpcol, lpnnonz);
   DUPLICATE_ARRAY_NULL(sdpi->blkmem, &(sdpi->lpval), lpval, lpnnonz);

   /* set the general information */
   sdpi->nvars = nvars;
   sdpi->nsdpblocks = nsdpblocks;

   sdpi->sdpconstnnonz = sdpconstnnonz;
   sdpi->sdpnnonz = sdpnnonz;

   /* LP part */
   sdpi->lpnnonz = lpnnonz;
   sdpi->nlpcons = nlpcons;

   sdpi->solved = FALSE;
   sdpi->infeasible = FALSE;
   sdpi->allfixed = FALSE;
   sdpi->nsdpcalls = 0;
   sdpi->niterations = 0;

   return SCIP_OKAY;
}

/** adds rows to the LP-Block
 *
 *  @note Arrays are not checked for duplicates, problems may appear if indices are added more than once.
 */
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
   )
{
   int i;

   SCIPdebugMessage("Adding %d LP-Constraints to SDP %d.\n", nrows, sdpi->sdpid);

   assert ( sdpi != NULL );

   if ( nrows == 0 )
      return SCIP_OKAY; /* nothing to do in this case */

   assert ( lhs != NULL );
   assert ( rhs != NULL );
   assert ( nnonz >= 0 );
   assert ( row != NULL );
   assert ( col != NULL );
   assert ( val != NULL );

   BMS_CALL( BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lplhs), sdpi->nlpcons, sdpi->nlpcons + nrows) ); /*lint !e776*/
   BMS_CALL( BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lprhs), sdpi->nlpcons, sdpi->nlpcons + nrows) ); /*lint !e776*/

   for (i = 0; i < nrows; i++)
   {
      sdpi->lplhs[sdpi->nlpcons + i] = lhs[i]; /*lint !e679*/
      sdpi->lprhs[sdpi->nlpcons + i] = rhs[i]; /*lint !e679*/
   }

   BMS_CALL( BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lprow), sdpi->lpnnonz, sdpi->lpnnonz + nnonz) ); /*lint !e776*/
   BMS_CALL( BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lpcol), sdpi->lpnnonz, sdpi->lpnnonz + nnonz) ); /*lint !e776*/
   BMS_CALL( BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lpval), sdpi->lpnnonz, sdpi->lpnnonz + nnonz) ); /*lint !e776*/

   for (i = 0; i < nnonz; i++)
   {
      assert ( 0 <= row[i] && row[i] < nrows );
      /* the new rows are added at the end, so the row indices are increased by the old number of LP-constraints */
      sdpi->lprow[sdpi->lpnnonz + i] = row[i] + sdpi->nlpcons; /*lint !e679*/

      assert ( 0 <= col[i] && col[i] < sdpi->nvars ); /* only existing vars should be added to the LP-constraints */
      sdpi->lpcol[sdpi->lpnnonz + i] = col[i]; /*lint !e679*/

      sdpi->lpval[sdpi->lpnnonz + i] = val[i]; /*lint !e679*/
   }

   sdpi->nlpcons = sdpi->nlpcons + nrows;
   sdpi->lpnnonz = sdpi->lpnnonz + nnonz;

   sdpi->solved = FALSE;
   sdpi->infeasible = FALSE;
   sdpi->nsdpcalls = 0;
   sdpi->niterations = 0;

   return SCIP_OKAY;
}

/** deletes all rows in the given range from the LP-Block */
SCIP_RETCODE SCIPsdpiDelLPRows(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int                   firstrow,           /**< first row to be deleted */
   int                   lastrow             /**< last row to be deleted */
   )
{
   int i;
   int deletedrows;
   int firstrowind;
   int lastrowind;
   int deletednonz;

   SCIPdebugMessage("Deleting rows %d to %d from SDP %d.\n", firstrow, lastrow, sdpi->sdpid);

   assert ( sdpi != NULL );
   assert ( firstrow >= 0 );
   assert ( firstrow <= lastrow );
   assert ( lastrow < sdpi->nlpcons );

   /* shorten the procedure if the whole LP-part is to be deleted */
   if (firstrow == 0 && lastrow == sdpi->nlpcons - 1)
   {
      BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpi->lpval), sdpi->lpnnonz);/*lint !e737*/
      BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpi->lprow), sdpi->lpnnonz);/*lint !e737*/
      BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpi->lpcol), sdpi->lpnnonz);/*lint !e737*/
      BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpi->lprhs), sdpi->nlpcons);/*lint !e737*/
      BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpi->lplhs), sdpi->nlpcons);/*lint !e737*/

      sdpi->lplhs = NULL;
      sdpi->lprhs = NULL;
      sdpi->lpcol = NULL;
      sdpi->lprow = NULL;
      sdpi->lpval = NULL;

      sdpi->nlpcons = 0;
      sdpi->lpnnonz = 0;

      sdpi->solved = FALSE;
      sdpi->infeasible = FALSE;
      sdpi->allfixed = FALSE;
      sdpi->nsdpcalls = 0;
      sdpi->niterations = 0;

      return SCIP_OKAY;
   }

   deletedrows = lastrow - firstrow + 1; /*lint !e834*/
   deletednonz = 0;

   /* first delete the left- and right-hand-sides */
   for (i = lastrow + 1; i < sdpi->nlpcons; i++) /* shift all rhs after the deleted rows */
   {
      sdpi->lplhs[i - deletedrows] = sdpi->lplhs[i]; /*lint !e679*/
      sdpi->lprhs[i - deletedrows] = sdpi->lprhs[i]; /*lint !e679*/
   }
   BMS_CALL( BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lplhs), sdpi->nlpcons, sdpi->nlpcons - deletedrows) ); /*lint !e776*/
   BMS_CALL( BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lprhs), sdpi->nlpcons, sdpi->nlpcons - deletedrows) ); /*lint !e776*/

   /* for deleting and reordering the lpnonzeroes, the arrays first have to be sorted to have the rows to be deleted together */
   SCIPsortIntIntReal(sdpi->lprow, sdpi->lpcol, sdpi->lpval, sdpi->lpnnonz); /* sort all arrays by non-decreasing row indices */

   firstrowind = -1;
   /*iterate over the lprowind array to find the first index belonging to a row that should be deleted */
   for (i = 0; i < sdpi->lpnnonz; i++)
   {
      if (sdpi->lprow[i] >= firstrow && sdpi->lprow[i] <= lastrow) /* the and part makes sure that there actually were some nonzeroes in these rows */
      {
         firstrowind = i;
         lastrowind = i;
         i++;
         break;
      }
   }

   if (firstrowind > -1) /* if this is still 0 there are no nonzeroes for the given rows */
   {
      /* now find the last occurence of one of the rows (as these are sorted all in between also belong to deleted rows and will be removed) */
      while (i < sdpi->lpnnonz && sdpi->lprow[i] <= lastrow)
      {
         lastrowind++; /*lint !e644*/
         i++;
      }
      deletednonz = lastrowind - firstrowind + 1; /*lint !e834*/

      /* finally shift all LP-array-entries after the deleted rows */
      for (i = lastrowind + 1; i < sdpi->lpnnonz; i++)
      {
         sdpi->lpcol[i - deletednonz] = sdpi->lpcol[i]; /*lint !e679*/
         /* all rowindices after the deleted ones have to be lowered to still have ongoing indices from 0 to nlpcons-1 */
         sdpi->lprow[i - deletednonz] = sdpi->lprow[i] - deletedrows;  /*lint !e679*/
         sdpi->lpval[i - deletednonz] = sdpi->lpval[i]; /*lint !e679*/
      }
   }

   BMS_CALL( BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lpcol), sdpi->lpnnonz, sdpi->lpnnonz - deletednonz) ); /*lint !e776*/
   BMS_CALL( BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lprow), sdpi->lpnnonz, sdpi->lpnnonz - deletednonz) ); /*lint !e776*/
   BMS_CALL( BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lpval), sdpi->lpnnonz, sdpi->lpnnonz - deletednonz) ); /*lint !e776*/
   sdpi->nlpcons = sdpi->nlpcons - deletedrows;
   sdpi->lpnnonz = sdpi->lpnnonz - deletednonz;

   sdpi->solved = FALSE;
   sdpi->infeasible = FALSE;
   sdpi->allfixed = FALSE;
   sdpi->nsdpcalls = 0;
   sdpi->niterations = 0;

   return SCIP_OKAY;
}

/** deletes LP-rows from SDP-interface */
SCIP_RETCODE SCIPsdpiDelLPRowset(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int*                  dstat               /**< deletion status of LP rows <br>
                                              *   input:  1 if row should be deleted, 0 otherwise <br>
                                              *   output: new position of row, -1 if row was deleted */
   )
{
   int i;
   int oldnlpcons;
   int deletedrows;

   SCIPdebugMessage("Calling SCIPsdpiDelLPRowset for SDP %d.\n", sdpi->sdpid);

   assert ( sdpi != NULL );
   assert ( dstat != NULL );

   oldnlpcons = sdpi->nlpcons;
   deletedrows = 0;

   for (i = 0; i < oldnlpcons; i++)
   {
      if (dstat[i] == 1)
      {
         /* delete this row, it is shifted by - deletedrows, because in this problem the earlier rows have already been deleted */
         SCIP_CALL( SCIPsdpiDelLPRows(sdpi, i - deletedrows, i - deletedrows) );
         dstat[i] = -1;
         deletedrows++;
      }
      else
         dstat[i] = i - deletedrows;
   }

   sdpi->solved = FALSE;
   sdpi->infeasible = FALSE;
   sdpi->allfixed = FALSE;
   sdpi->nsdpcalls = 0;
   sdpi->niterations = 0;

   return SCIP_OKAY;
}

/** clears the whole SDP */
SCIP_RETCODE SCIPsdpiClear(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   )
{
   assert( sdpi != NULL );

   SCIPdebugMessage("Called SCIPsdpiClear in SDP %d.\n", sdpi->sdpid);

   /* we reset all counters */
   sdpi->sdpid = 1;
   SCIP_CALL( SCIPsdpiSolverResetCounter(sdpi->sdpisolver) );

   return SCIP_OKAY;
}

/** changes objective coefficients of variables */
SCIP_RETCODE SCIPsdpiChgObj(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int                   nvars,              /**< number of variables to change objective coefficients for */
   const int*            ind,                /**< variables indices */
   const SCIP_Real*      obj                 /**< new objective coefficients */
   )
{
   int i;

   SCIPdebugMessage("Changing %d objective coefficients in SDP %d\n", nvars, sdpi->sdpid);

   assert( sdpi != NULL );
   assert( ind != NULL );
   assert( obj != NULL );

   for (i = 0; i < nvars; i++)
   {
      assert( 0 <= ind[i] && ind[i] < sdpi->nvars );
      sdpi->obj[ind[i]] = obj[i];
   }

   sdpi->solved = FALSE;
   sdpi->nsdpcalls = 0;
   sdpi->niterations = 0;

   return SCIP_OKAY;
}

/** changes lower and upper bounds of variables */
SCIP_RETCODE SCIPsdpiChgBounds(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int                   nvars,              /**< number of variables to change bounds for */
   const int*            ind,                /**< variables indices */
   const SCIP_Real*      lb,                 /**< values for the new lower bounds */
   const SCIP_Real*      ub                  /**< values for the new upper bounds */
   )
{
   int i;

   SCIPdebugMessage("Changing %d variable bounds in SDP %d\n", nvars, sdpi->sdpid);

   assert( sdpi != NULL );
   assert( ind != NULL );
   assert( lb != NULL );
   assert( ub != NULL );

   for (i = 0; i < nvars; i++)
   {
      assert( 0 <= ind[i] && ind[i] < sdpi->nvars );
      sdpi->lb[ind[i]] = lb[i];
      sdpi->ub[ind[i]] = ub[i];
   }

   sdpi->solved = FALSE;
   sdpi->infeasible = FALSE;
   sdpi->allfixed = FALSE;
   sdpi->nsdpcalls = 0;
   sdpi->niterations = 0;

   return SCIP_OKAY;
}

/** changes left- and right-hand sides of LP rows */
SCIP_RETCODE SCIPsdpiChgLPLhRhSides(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int                   nrows,              /**< number of LP rows to change right hand sides for */
   const int*            ind,                /**< row indices between 1 and nlpcons */
   const SCIP_Real*      lhs,                /**< new values for left-hand sides */
   const SCIP_Real*      rhs                 /**< new values for right-hand sides */
   )
{
   int i;

   SCIPdebugMessage("Changing %d left and right hand sides of SDP %d\n", nrows, sdpi->sdpid);

   assert( sdpi != NULL );
   assert( 0 <= nrows && nrows <= sdpi->nlpcons );
   assert( ind != NULL );
   assert( lhs != NULL );
   assert( rhs != NULL );

   for (i = 0; i < nrows; i++)
   {
      assert ( ind[i] >= 0 );
      assert ( ind[i] < sdpi->nlpcons );
      sdpi->lplhs[ind[i]] = lhs[i];
      sdpi->lprhs[ind[i]] = rhs[i];
   }

   sdpi->solved = FALSE;
   sdpi->infeasible = FALSE;
   sdpi->allfixed = FALSE;
   sdpi->nsdpcalls = 0;
   sdpi->niterations = 0;

   return SCIP_OKAY;
}


/*
 * Data Accessing Methods
 */

/**@name Data Accessing Methods */
/**@{ */

/** gets the number of LP-rows in the SDP */
SCIP_RETCODE SCIPsdpiGetNLPRows(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int*                  nlprows             /**< pointer to store the number of rows */
   )
{
   assert( sdpi != NULL );
   assert( nlprows != NULL );

   *nlprows = sdpi->nlpcons;

   return SCIP_OKAY;
}

/** gets the number of SDP-Blocks in the SDP */
SCIP_RETCODE SCIPsdpiGetNSDPBlocks(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int*                  nsdpblocks          /**< pointer to store the number of blocks */
   )
{
   assert( sdpi != NULL );
   assert( nsdpblocks != NULL );

   *nsdpblocks = sdpi->nsdpblocks;

   return SCIP_OKAY;
}

/** gets the number of variables in the SDP */
SCIP_RETCODE SCIPsdpiGetNVars(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int*                  nvars               /**< pointer to store the number of variables */
   )
{
   assert( sdpi != NULL );
   assert( nvars != NULL );

   *nvars = sdpi->nvars;

   return SCIP_OKAY;
}

/** gets the number of nonzero elements in the SDP-constraint-matrices */
SCIP_RETCODE SCIPsdpiGetSDPNNonz(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int*                  nnonz               /**< pointer to store the number of nonzeros in the SDP-constraint-matrices */
   )
{
   assert( sdpi != NULL );
   assert( nnonz != NULL );

   *nnonz = sdpi->sdpnnonz;

   return SCIP_OKAY;
}

/** gets the number of nonzero elements in the constant matrices of the SDP-Blocks */
SCIP_RETCODE SCIPsdpiGetConstNNonz(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int*                  nnonz               /**< pointer to store the number of nonzeros in the constant matrices of the SDP-Blocks */
   )
{
   assert( sdpi != NULL );
   assert( nnonz != NULL );

   *nnonz = sdpi->sdpconstnnonz;

   return SCIP_OKAY;
}

/** gets the number of nonzero elements in the LP-Matrix */
SCIP_RETCODE SCIPsdpiGetLPNNonz(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int*                  nnonz               /**< pointer to store the number of nonzeros in the LP Matrix */
   )
{
   assert( sdpi != NULL );
   assert( nnonz != NULL );

   *nnonz = sdpi->lpnnonz;

   return SCIP_OKAY;
}

/** gets objective coefficients from SDP-interface */
SCIP_RETCODE SCIPsdpiGetObj(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int                   firstvar,           /**< first variable to get objective coefficient for */
   int                   lastvar,            /**< last variable to get objective coefficient for */
   SCIP_Real*            vals                /**< pointer to store objective coefficients (memory of size lastvar - firstvar + 1 needs to be allocated) */
   )
{
   int i;

   assert( sdpi != NULL );
   assert( firstvar >= 0 );
   assert( firstvar <= lastvar );
   assert( lastvar < sdpi->nvars);
   assert( vals != NULL );

   for (i = 0; i < lastvar - firstvar + 1; i++) /*lint !e834*/
      vals[i] = sdpi->obj[firstvar + i]; /*lint !e679*/

   return SCIP_OKAY;
}

/** gets current variable lower and/or upper bounds from SDP-interface */
SCIP_RETCODE SCIPsdpiGetBounds(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int                   firstvar,           /**< first variable to get bounds for */
   int                   lastvar,            /**< last variable to get bounds for */
   SCIP_Real*            lbs,                /**< pointer to store lower bound values (memory of size lastvar - firstvar + 1 needs to be allocated), or NULL */
   SCIP_Real*            ubs                 /**< pointer to store upper bound values (memory of size lastvar - firstvar + 1 needs to be allocated), or NULL */
   )
{
   int i;

   assert( sdpi != NULL );
   assert( firstvar >= 0 );
   assert( firstvar <= lastvar );
   assert( lastvar < sdpi->nvars);
   assert( lbs != NULL );
   assert( ubs != NULL );

   for (i = 0; i < lastvar - firstvar + 1; i++) /*lint !e834*/
   {
      if (lbs != NULL)
         lbs[i] = sdpi->lb[firstvar + i]; /*lint !e679*/
      if (ubs != NULL)
         ubs[i] = sdpi->ub[firstvar + i]; /*lint !e679*/
   }
   return SCIP_OKAY;
}

/** gets current left-hand sides from SDP-interface */
SCIP_RETCODE SCIPsdpiGetLhSides(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int                   firstrow,           /**< first row to get sides for */
   int                   lastrow,            /**< last row to get sides for */
   SCIP_Real*            lhss                /**< pointer to store left-hand side values (memory of size lastvar - firstvar + 1 needs to be allocated) */
   )
{
   int i;

   assert( sdpi != NULL );
   assert( firstrow >= 0 );
   assert( firstrow <= lastrow );
   assert( lastrow < sdpi->nlpcons);
   assert( lhss != NULL );

   for (i = 0; i < lastrow - firstrow + 1; i++) /*lint !e834*/
      lhss[firstrow + i] = sdpi->lplhs[i]; /*lint !e679*/

   return SCIP_OKAY;
}

/** gets current right-hand sides from SDP-interface */
SCIP_RETCODE SCIPsdpiGetRhSides(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int                   firstrow,           /**< first row to get sides for */
   int                   lastrow,            /**< last row to get sides for */
   SCIP_Real*            rhss                /**< pointer to store right-hand side values (memory of size lastvar - firstvar + 1 needs to be allocated) */
   )
{
   int i;

   assert( sdpi != NULL );
   assert( firstrow >= 0 );
   assert( firstrow <= lastrow );
   assert( lastrow < sdpi->nlpcons);
   assert( rhss != NULL );

   for (i = 0; i < lastrow - firstrow + 1; i++) /*lint !e834*/
      rhss[firstrow + i] = sdpi->lprhs[i]; /*lint !e679*/

   return SCIP_OKAY;
}


/**@} */



/*
 * Solving Methods
 */

/**@name Solving Methods */
/**@{ */

/** solves the SDP, as start optionally a starting point for the solver may be given, if it is NULL, the solver will start from scratch */
SCIP_RETCODE SCIPsdpiSolve(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   SCIP_Real*            start,              /**< NULL or a starting point for the solver, this should have length nvars */
   SCIP_SDPSOLVERSETTING startsettings,      /**< settings used to start with in SDPA, currently not used for DSDP or MOSEK, set this to
                                              *   SCIP_SDPSOLVERSETTING_UNSOLVED to ignore it and start from scratch */
   SCIP_Bool             enforceslatercheck, /**< always check for Slater condition in case the problem could not be solved and printf the solution
                                              *   of this check */
   SCIP_Real             timelimit           /**< after this many seconds solving will be aborted (currently only implemented for DSDP and MOSEK) */
   )
{
   int* sdpconstnblocknonz = NULL;
   int** sdpconstrow = NULL;
   int** sdpconstcol = NULL;
   SCIP_Real** sdpconstval = NULL;
   int** indchanges = NULL;
   int* nremovedinds = NULL;
   SCIP_Real* lplhsafterfix;
   SCIP_Real* lprhsafterfix;
   SCIP_Real solvertimelimit;
   SCIP_Bool fixingfound;
   clock_t starttime;
   clock_t currenttime;
   int* rowsnactivevars;
   int* blockindchanges;
   int sdpconstnnonz;
   int nactivelpcons;
   int nremovedblocks = 0;
   int block;
   int naddediterations;
   int naddedsdpcalls;

   assert( sdpi != NULL );

   starttime = clock();

   SCIPdebugMessage("Forwarding SDP %d to solver!\n", sdpi->sdpid);

   sdpi->penalty = FALSE;
   sdpi->bestbound = -SCIPsdpiSolverInfinity(sdpi->sdpisolver);
   sdpi->solved = FALSE;
   sdpi->nsdpcalls = 0;
   sdpi->niterations = 0;

   /* allocate memory for computing the constant matrix after fixings and finding empty rows and columns, this is as much as might possibly be
    * needed, this will be shrinked again before solving */
   BMS_CALL( BMSallocBlockMemoryArray(sdpi->blkmem, &sdpconstnblocknonz, sdpi->nsdpblocks) );
   BMS_CALL( BMSallocBlockMemoryArray(sdpi->blkmem, &sdpconstrow, sdpi->nsdpblocks) );
   BMS_CALL( BMSallocBlockMemoryArray(sdpi->blkmem, &sdpconstcol, sdpi->nsdpblocks) );
   BMS_CALL( BMSallocBlockMemoryArray(sdpi->blkmem, &sdpconstval, sdpi->nsdpblocks) );
   BMS_CALL( BMSallocBlockMemoryArray(sdpi->blkmem, &indchanges, sdpi->nsdpblocks) );
   BMS_CALL( BMSallocBlockMemoryArray(sdpi->blkmem, &nremovedinds, sdpi->nsdpblocks) );
   BMS_CALL( BMSallocBlockMemoryArray(sdpi->blkmem, &blockindchanges, sdpi->nsdpblocks) );
   BMS_CALL( BMSallocBlockMemoryArray(sdpi->blkmem, &lplhsafterfix, sdpi->nlpcons) );
   BMS_CALL( BMSallocBlockMemoryArray(sdpi->blkmem, &lprhsafterfix, sdpi->nlpcons) );
   BMS_CALL( BMSallocBlockMemoryArray(sdpi->blkmem, &rowsnactivevars, sdpi->nlpcons) );

   for (block = 0; block < sdpi->nsdpblocks; block++)
   {
      sdpconstrow[block] = NULL;
      sdpconstcol[block] = NULL;
      sdpconstval[block] = NULL;
      indchanges[block] = NULL;
      BMS_CALL( BMSallocBlockMemoryArray(sdpi->blkmem, &(indchanges[block]), sdpi->sdpblocksizes[block]) );
      BMS_CALL( BMSallocBlockMemoryArray(sdpi->blkmem, &(sdpconstrow[block]), sdpi->sdpnnonz + sdpi->sdpconstnnonz) ); /*lint !e776*/
      BMS_CALL( BMSallocBlockMemoryArray(sdpi->blkmem, &(sdpconstcol[block]), sdpi->sdpnnonz + sdpi->sdpconstnnonz) ); /*lint !e776*/
      BMS_CALL( BMSallocBlockMemoryArray(sdpi->blkmem, &(sdpconstval[block]), sdpi->sdpnnonz + sdpi->sdpconstnnonz) ); /*lint !e776*/
   }

   /* compute the lplphss and lprhss, detect empty rows and check for additional variable fixings caused by boundchanges from
    * lp rows with a single active variable */
   do
   {
      fixingfound = FALSE;
      SCIP_CALL( computeLpLhsRhsAfterFixings(sdpi, &nactivelpcons, lplhsafterfix, lprhsafterfix, rowsnactivevars, &fixingfound) );
   }
   while ( fixingfound );

   /* initialize sdpconstnblocknonz */
   for (block = 0; block < sdpi->nsdpblocks; block++)
      sdpconstnblocknonz[block] = sdpi->sdpnnonz + sdpi->sdpconstnnonz;

   SCIP_CALL( compConstMatAfterFixings(sdpi, &sdpconstnnonz, sdpconstnblocknonz, sdpconstrow, sdpconstcol, sdpconstval) );

   /* shrink the constant arrays after the number of fixed nonzeros is known */
   for (block = 0; block < sdpi->nsdpblocks; block++)
   {
      assert ( sdpconstnblocknonz[block] <= sdpi->sdpnnonz + sdpi->sdpconstnnonz ); /* otherwise the memory wasn't sufficient,
                                                                                     * but we allocated more than enough */
      BMS_CALL( BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpconstrow[block]), sdpi->sdpnnonz + sdpi->sdpconstnnonz, sdpconstnblocknonz[block]) ); /*lint !e776*/
      BMS_CALL( BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpconstcol[block]), sdpi->sdpnnonz + sdpi->sdpconstnnonz, sdpconstnblocknonz[block]) ); /*lint !e776*/
      BMS_CALL( BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpconstval[block]), sdpi->sdpnnonz + sdpi->sdpconstnnonz, sdpconstnblocknonz[block]) ); /*lint !e776*/
   }

   SCIP_CALL( findEmptyRowColsSDP(sdpi, sdpconstnblocknonz, sdpconstrow, sdpconstcol, sdpconstval, indchanges, nremovedinds, blockindchanges, &nremovedblocks) );

   /* check if all variables are fixed, if this is the case, check if the remaining solution if feasible (we only need to check the SDP-constraint,
    * the linear constraints were already checked in computeLpLhsRhsAfterFixings) */
   SCIP_CALL( checkAllFixed(sdpi) );
   if ( sdpi->allfixed && ! sdpi->infeasible )
   {
      SCIP_CALL( checkFixedFeasibilitySdp(sdpi, sdpconstnblocknonz, sdpconstrow, sdpconstcol, sdpconstval, indchanges, nremovedinds, blockindchanges) );
   }

   if ( sdpi->infeasible )
   {
      SCIPdebugMessage("SDP %d not given to solver, as infeasibility was detected during presolving!\n", sdpi->sdpid++);
      SCIP_CALL( SCIPsdpiSolverIncreaseCounter(sdpi->sdpisolver) );

      sdpi->solved = TRUE;
      sdpi->dualslater = SCIP_SDPSLATER_NOINFO;
      sdpi->primalslater = SCIP_SDPSLATER_NOINFO;
   }
   else if ( sdpi->allfixed )
   {
      SCIPdebugMessage("SDP %d not given to solver, as all variables were fixed during presolving (the solution was feasible)!\n", sdpi->sdpid++);
      SCIP_CALL( SCIPsdpiSolverIncreaseCounter(sdpi->sdpisolver) );

      sdpi->solved = TRUE;
      sdpi->dualslater = SCIP_SDPSLATER_NOINFO;
      sdpi->primalslater = SCIP_SDPSLATER_NOINFO;
   }
   else
   {
      if ( sdpi->slatercheck )
      {
         SCIP_CALL( checkSlaterCondition(sdpi, timelimit, starttime, sdpconstnblocknonz, sdpconstrow, sdpconstcol, sdpconstval, indchanges,
                  nremovedinds, lplhsafterfix, lprhsafterfix, rowsnactivevars, blockindchanges, sdpconstnnonz, nactivelpcons, nremovedblocks, FALSE) );
      }

      /* compute the timit limit to set for the solver */
      solvertimelimit = timelimit;
      if ( ! SCIPsdpiIsInfinity(sdpi, solvertimelimit) )
      {
         currenttime = clock();
         solvertimelimit -= (SCIP_Real)(currenttime - starttime) / (SCIP_Real) CLOCKS_PER_SEC;/*lint !e620*/
      }

      /* try to solve the problem */
      SCIP_CALL( SCIPsdpiSolverLoadAndSolve(sdpi->sdpisolver, sdpi->nvars, sdpi->obj, sdpi->lb, sdpi->ub,
            sdpi->nsdpblocks, sdpi->sdpblocksizes, sdpi->sdpnblockvars, sdpconstnnonz,
            sdpconstnblocknonz, sdpconstrow, sdpconstcol, sdpconstval,
            sdpi->sdpnnonz, sdpi->sdpnblockvarnonz, sdpi->sdpvar, sdpi->sdprow, sdpi->sdpcol,
            sdpi->sdpval, indchanges, nremovedinds, blockindchanges, nremovedblocks, nactivelpcons, sdpi->nlpcons, lplhsafterfix, lprhsafterfix,
            rowsnactivevars, sdpi->lpnnonz, sdpi->lprow, sdpi->lpcol, sdpi->lpval, start, startsettings, solvertimelimit) );

      sdpi->solved = TRUE;

      /* add iterations and sdpcalls */
      naddediterations = 0;
      SCIP_CALL( SCIPsdpiSolverGetIterations(sdpi->sdpisolver, &naddediterations) );
      sdpi->niterations += naddediterations;
      naddedsdpcalls = 0;
      SCIP_CALL( SCIPsdpiSolverGetSdpCalls(sdpi->sdpisolver, &naddedsdpcalls) );
      sdpi->nsdpcalls += naddedsdpcalls;

      /* if the solver didn't produce a satisfactory result, we have to try with a penalty formulation */
      if ( ! SCIPsdpiSolverIsAcceptable(sdpi->sdpisolver) && ! SCIPsdpiSolverIsTimelimExc(sdpi->sdpisolver) )
      {
         SCIP_Real penaltyparam;
         SCIP_Real penaltyparamfact;
         SCIP_Real gaptol;
         SCIP_Real gaptolfact;
         SCIP_Bool feasorig;
         SCIP_Bool penaltybound;
         SCIP_Real objbound;
         SCIP_Real objval;

         feasorig = FALSE;
         penaltybound = TRUE;

         /* first check feasibility using the penalty approach */

         /* compute the timit limit to set for the solver */
         solvertimelimit = timelimit;
         if ( ! SCIPsdpiIsInfinity(sdpi, solvertimelimit) )
         {
            currenttime = clock();
            solvertimelimit -= (SCIP_Real)(currenttime - starttime) / (SCIP_Real) CLOCKS_PER_SEC;/*lint !e620*/
         }

         /* we solve the problem with a slack variable times identity added to the constraints and trying to minimize this slack variable r, if
          * the optimal objective is bigger than feastol, then we know that the problem is infeasible */
         SCIP_CALL( SCIPsdpiSolverLoadAndSolveWithPenalty(sdpi->sdpisolver, 1.0, FALSE, FALSE, sdpi->nvars, sdpi->obj, sdpi->lb, sdpi->ub,
               sdpi->nsdpblocks, sdpi->sdpblocksizes, sdpi->sdpnblockvars, sdpconstnnonz,
               sdpconstnblocknonz, sdpconstrow, sdpconstcol, sdpconstval,
               sdpi->sdpnnonz, sdpi->sdpnblockvarnonz, sdpi->sdpvar, sdpi->sdprow, sdpi->sdpcol,
               sdpi->sdpval, indchanges, nremovedinds, blockindchanges, nremovedblocks, nactivelpcons, sdpi->nlpcons, lplhsafterfix, lprhsafterfix,
               rowsnactivevars, sdpi->lpnnonz, sdpi->lprow, sdpi->lpcol, sdpi->lpval, start, SCIP_SDPSOLVERSETTING_UNSOLVED, solvertimelimit,
               &feasorig, &penaltybound) );

         /* add iterations and sdpcalls */
         naddediterations = 0;
         SCIP_CALL( SCIPsdpiSolverGetIterations(sdpi->sdpisolver, &naddediterations) );
         sdpi->niterations += naddediterations;
         naddedsdpcalls = 0;
         SCIP_CALL( SCIPsdpiSolverGetSdpCalls(sdpi->sdpisolver, &naddedsdpcalls) );
         sdpi->nsdpcalls += naddedsdpcalls;

         /* get objective value */
         if ( SCIPsdpiSolverWasSolved(sdpi->sdpisolver) )
         {
            SCIP_CALL( SCIPsdpiSolverGetObjval(sdpi->sdpisolver, &objval) );
         }
         else
            objval = -SCIPsdpiInfinity(sdpi);

         /* If the penalty formulation was successfully solved and has a strictly positive objective value, we know that
          * the problem is infeasible. Note that we need to check against the maximum of feastol and gaptol, since this
          * is the objective of an SDP which is only exact up to gaptol, and cutting a feasible node off is an error
          * while continueing with an infeasible problem only takes additional time until we found out again later.
          */
         if ( (SCIPsdpiSolverIsOptimal(sdpi->sdpisolver) && (objval > (sdpi->feastol > sdpi->gaptol ? sdpi->feastol : sdpi->gaptol))) ||
               (SCIPsdpiSolverWasSolved(sdpi->sdpisolver) && SCIPsdpiSolverIsDualInfeasible(sdpi->sdpisolver)) )
         {
            SCIPdebugMessage("SDP %d found infeasible using penalty formulation, maximum of smallest eigenvalue is %f.\n", sdpi->sdpid, -1.0 * objval);
            sdpi->penalty = TRUE;
            sdpi->infeasible = TRUE;
         }
         else
         {
            feasorig = FALSE;
            penaltybound = TRUE;

            penaltyparam = sdpi->penaltyparam;

            /* we compute the factor to increase with as n-th root of the total increase until the maximum, where n is the number of iterations
             * (for npenaltyincr = 0 we make sure that the parameter is too large after the first change)
             */
            penaltyparamfact = sdpi->npenaltyincr > 0 ? pow((sdpi->maxpenaltyparam / sdpi->penaltyparam), 1.0/sdpi->npenaltyincr) :
                  2*sdpi->maxpenaltyparam / sdpi->penaltyparam;
            gaptol = sdpi->gaptol;
            gaptolfact = sdpi->npenaltyincr > 0 ? pow((MIN_GAPTOL / sdpi->gaptol), 1.0/sdpi->npenaltyincr) : 0.5 * MIN_GAPTOL / sdpi->gaptol;

            /* increase penalty-param and decrease feasibility tolerance until we find a feasible solution or reach the final bound for either one of them */
            while ( ( ! SCIPsdpiSolverIsAcceptable(sdpi->sdpisolver) || ! feasorig ) &&
                  ( penaltyparam < sdpi->maxpenaltyparam + sdpi->epsilon ) && ( gaptol > 0.99 * MIN_GAPTOL ) && ( ! SCIPsdpiSolverIsTimelimExc(sdpi->sdpisolver) ))
            {
               SCIPdebugMessage("Solver did not produce an acceptable result, trying SDP %d again with penaltyparameter %f\n", sdpi->sdpid, penaltyparam);

               /* compute the timit limit to set for the solver */
               solvertimelimit = timelimit;
               if ( ! SCIPsdpiIsInfinity(sdpi, solvertimelimit) )
               {
                  currenttime = clock();
                  solvertimelimit -= (SCIP_Real)(currenttime - starttime) / (SCIP_Real) CLOCKS_PER_SEC;/*lint !e620*/

                  if ( solvertimelimit <= 0 )
                     break;
               }

               SCIP_CALL( SCIPsdpiSolverLoadAndSolveWithPenalty(sdpi->sdpisolver, penaltyparam, TRUE, TRUE, sdpi->nvars, sdpi->obj,
                     sdpi->lb, sdpi->ub, sdpi->nsdpblocks, sdpi->sdpblocksizes, sdpi->sdpnblockvars, sdpconstnnonz,
                     sdpconstnblocknonz, sdpconstrow, sdpconstcol, sdpconstval,
                     sdpi->sdpnnonz, sdpi->sdpnblockvarnonz, sdpi->sdpvar, sdpi->sdprow, sdpi->sdpcol,
                     sdpi->sdpval, indchanges, nremovedinds, blockindchanges, nremovedblocks, nactivelpcons, sdpi->nlpcons, lplhsafterfix, lprhsafterfix,
                     rowsnactivevars, sdpi->lpnnonz, sdpi->lprow, sdpi->lpcol, sdpi->lpval, start, startsettings, solvertimelimit, &feasorig, &penaltybound) );

               /* add iterations and sdpcalls */
               naddediterations = 0;
               SCIP_CALL( SCIPsdpiSolverGetIterations(sdpi->sdpisolver, &naddediterations) );
               sdpi->niterations += naddediterations;
               naddedsdpcalls = 0;
               SCIP_CALL( SCIPsdpiSolverGetSdpCalls(sdpi->sdpisolver, &naddedsdpcalls) );
               sdpi->nsdpcalls += naddedsdpcalls;

               /* If the solver did not converge, we increase the penalty parameter */
               if ( ! SCIPsdpiSolverIsAcceptable(sdpi->sdpisolver) )
               {
                  penaltyparam *= penaltyparamfact;
                  SCIPdebugMessage("Solver did not converge even with penalty formulation, increasing penaltyparameter.\n");
                  continue;
               }

               /* if we succeeded to solve the problem, update the bound */
               SCIP_CALL( SCIPsdpiSolverGetObjval(sdpi->sdpisolver, &objbound) );
               if ( objbound > sdpi->bestbound + sdpi->gaptol )
                  sdpi->bestbound = objbound;

               /* If we don't get a feasible solution to our original problem we have to update either Gamma (if the penalty bound was active
                * in the primal problem) or gaptol (otherwise) */
               if ( ! feasorig )
               {
                  if ( penaltybound )
                  {
                     penaltyparam *= penaltyparamfact;
                     SCIPdebugMessage("Penalty formulation produced a result which is infeasible for the original problem, increasing penaltyparameter\n");
                  }
                  else
                  {
                     gaptol *= gaptolfact;
                     SCIP_CALL_PARAM( SCIPsdpiSolverSetRealpar(sdpi->sdpisolver, SCIP_SDPPAR_GAPTOL, gaptol) );
                     SCIPdebugMessage("Penalty formulation produced a result which is infeasible for the original problem, even though primal penalty "
                           "bound was not reached, decreasing tolerance for duality gap in SDP-solver\n");
                  }
               }
            }

            /* reset the tolerance in the SDP-solver */
            if ( gaptol > sdpi->gaptol )
            {
               SCIP_CALL_PARAM( SCIPsdpiSolverSetRealpar(sdpi->sdpisolver, SCIP_SDPPAR_GAPTOL, sdpi->gaptol) );
            }

            /* check if we were able to solve the problem in the end */
            if ( SCIPsdpiSolverIsAcceptable(sdpi->sdpisolver) && feasorig )
            {
               sdpi->penalty = TRUE;
               sdpi->solved = TRUE;
            }
#if 0 /* we don't really know if it is infeasible or just ill-posed (no KKT-point) */
            else if ( SCIPsdpiSolverIsAcceptable(sdpi->sdpisolver) && ! feasorig )
            {
               SCIPdebugMessage("Problem was found to be infeasible using a penalty formulation \n");
               sdpi->infeasible = TRUE;
               sdpi->penalty = TRUE;
               sdpi->solved = TRUE;
            }
#endif
            else
            {
               SCIPdebugMessage("SDP-Solver could not solve the problem even after using a penalty formulation \n");
               sdpi->solved = FALSE;
               sdpi->penalty = TRUE;
            }

            /* if we still didn't succeed and enforceslatercheck was set, we finally test for the Slater condition to give a reason for failure */
            if ( sdpi->solved == FALSE && enforceslatercheck)
            {
               SCIP_CALL( checkSlaterCondition(sdpi, timelimit, starttime, sdpconstnblocknonz, sdpconstrow, sdpconstcol, sdpconstval, indchanges,
                                 nremovedinds, lplhsafterfix, lprhsafterfix, rowsnactivevars, blockindchanges, sdpconstnnonz, nactivelpcons, nremovedblocks, TRUE) );
            }
            else if ( sdpi->solved == FALSE )
            {
#if 0
               SCIPmessagePrintInfo(sdpi->messagehdlr, "Numerical trouble\n");
#else
               SCIPdebugMessage("SDP-Interface was unable to solve SDP %d\n", sdpi->sdpid);/*lint !e687*/
#endif
            }
         }
      }
   }

   /* empty the memory allocated here */
   for (block = 0; block < sdpi->nsdpblocks; block++)
   {
      BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpconstval[block]), sdpconstnblocknonz[block]);/*lint !e737*/
      BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpconstcol[block]), sdpconstnblocknonz[block]);/*lint !e737*/
      BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpconstrow[block]), sdpconstnblocknonz[block]);/*lint !e737*/
      BMSfreeBlockMemoryArray(sdpi->blkmem, &(indchanges[block]), sdpi->sdpblocksizes[block]);/*lint !e737*/
   }
   BMSfreeBlockMemoryArray(sdpi->blkmem, &rowsnactivevars, sdpi->nlpcons);/*lint !e737*/
   BMSfreeBlockMemoryArray(sdpi->blkmem, &lprhsafterfix, sdpi->nlpcons);/*lint !e737*/
   BMSfreeBlockMemoryArray(sdpi->blkmem, &lplhsafterfix, sdpi->nlpcons);/*lint !e737*/
   BMSfreeBlockMemoryArray(sdpi->blkmem, &blockindchanges, sdpi->nsdpblocks);/*lint !e737*/
   BMSfreeBlockMemoryArray(sdpi->blkmem, &nremovedinds, sdpi->nsdpblocks);/*lint !e737*/
   BMSfreeBlockMemoryArray(sdpi->blkmem, &indchanges, sdpi->nsdpblocks);/*lint !e737*/
   BMSfreeBlockMemoryArray(sdpi->blkmem, &sdpconstval, sdpi->nsdpblocks);/*lint !e737*/
   BMSfreeBlockMemoryArray(sdpi->blkmem, &sdpconstcol, sdpi->nsdpblocks);/*lint !e737*/
   BMSfreeBlockMemoryArray(sdpi->blkmem, &sdpconstrow, sdpi->nsdpblocks);/*lint !e737*/
   BMSfreeBlockMemoryArray(sdpi->blkmem, &sdpconstnblocknonz, sdpi->nsdpblocks);/*lint !e737*/

   sdpi->sdpid++;

   return SCIP_OKAY;
}




/*
 * Solution Information Methods
 */

/**@name Solution Information Methods */
/**@{ */

/** returns whether a solve method was successfully called after the last modification of the SDP */
SCIP_Bool SCIPsdpiWasSolved(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   )
{
   assert( sdpi != NULL );

   return ( sdpi->solved && SCIPsdpiSolverWasSolved(sdpi->sdpisolver) );
}

/** returns whether the original problem was solved, if SCIPsdpiWasSolved = true and SCIPsdpiSolvedOrig = false, then a penalty formulation was solved */
SCIP_Bool SCIPsdpiSolvedOrig(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   )
{
   assert( sdpi != NULL );

   return ( SCIPsdpiWasSolved(sdpi) && (! sdpi->penalty) );
}

/** returns true if the solver could determine whether the problem is feasible, so it returns true if the
 *  solver knows that the problem is feasible/infeasible/unbounded, it returns false if the solver does not know
 *  anything about the feasibility status and thus the functions IsPrimalFeasible etc. should not be used */
SCIP_Bool SCIPsdpiFeasibilityKnown(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   )
{
   assert( sdpi != NULL );
   CHECK_IF_SOLVED_BOOL(sdpi);

   if ( sdpi->infeasible || sdpi->allfixed )
      return TRUE;

   return SCIPsdpiSolverFeasibilityKnown(sdpi->sdpisolver);
}

/** gets information about primal and dual feasibility of the current SDP-solution */
SCIP_RETCODE SCIPsdpiGetSolFeasibility(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   SCIP_Bool*            primalfeasible,     /**< pointer to store the primal feasibility status */
   SCIP_Bool*            dualfeasible        /**< pointer to store the dual feasibility status */
   )
{
   assert( sdpi != NULL );
   CHECK_IF_SOLVED(sdpi);

   if ( sdpi->infeasible )
   {
      SCIPdebugMessage("Problem was found infeasible during preprocessing, primal feasibility not available\n");
      *dualfeasible = FALSE;
      return SCIP_OKAY;
   }
   else if ( sdpi->allfixed )
   {
      SCIPdebugMessage("All variables were fixed during preprocessing, dual problem is feasible, primal feasibility not available\n");
      *dualfeasible = TRUE;
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPsdpiSolverGetSolFeasibility(sdpi->sdpisolver, primalfeasible, dualfeasible) );

   return SCIP_OKAY;
}

/** returns TRUE iff SDP is proven to be primal unbounded;
 *  returns FALSE with a debug-message if the solver could not determine feasibility */
SCIP_Bool SCIPsdpiIsPrimalUnbounded(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   )
{
   assert( sdpi != NULL );
   CHECK_IF_SOLVED_BOOL(sdpi);

   if ( sdpi->infeasible )
   {
      SCIPdebugMessage("Problem was found infeasible during preprocessing, primal unboundedness not available\n");
      return FALSE;
   }
   else if ( sdpi->allfixed )
   {
      SCIPdebugMessage("All variables were fixed during preprocessing, primal unboundedness not available\n");
      return FALSE;
   }

   return SCIPsdpiSolverIsPrimalUnbounded(sdpi->sdpisolver);
}

/** returns TRUE iff SDP is proven to be primal infeasible;
 *  returns FALSE with a debug-message if the solver could not determine feasibility */
SCIP_Bool SCIPsdpiIsPrimalInfeasible(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   )
{
   assert( sdpi != NULL );
   CHECK_IF_SOLVED_BOOL(sdpi);

   if ( sdpi->infeasible )
   {
      SCIPdebugMessage("Problem was found infeasible during preprocessing, primal feasibility not available\n");
      return FALSE;
   }
   else if ( sdpi->allfixed )
   {
      SCIPdebugMessage("All variables were fixed during preprocessing, primal feasibility not available\n");
      return FALSE;
   }

   return SCIPsdpiSolverIsPrimalInfeasible(sdpi->sdpisolver);
}

/** returns TRUE iff SDP is proven to be primal feasible;
 *  returns FALSE with a debug-message if the solver could not determine feasibility */
SCIP_Bool SCIPsdpiIsPrimalFeasible(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   )
{
   assert(sdpi != NULL );
   CHECK_IF_SOLVED_BOOL(sdpi);

   if ( sdpi->infeasible )
   {
      SCIPdebugMessage("Problem was found infeasible during preprocessing, primal feasibility not available\n");
      return FALSE;
   }
   else if ( sdpi->allfixed )
   {
      SCIPdebugMessage("All variables fixed during preprocessing, primal feasibility not available\n");
      return FALSE;
   }

   return SCIPsdpiSolverIsPrimalFeasible(sdpi->sdpisolver);
}

/** returns TRUE iff SDP is proven to be dual unbounded;
 *  returns FALSE with a debug-message if the solver could not determine feasibility */
SCIP_Bool SCIPsdpiIsDualUnbounded(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   )
{
   assert( sdpi != NULL );
   CHECK_IF_SOLVED_BOOL(sdpi);

   if ( sdpi->infeasible )
   {
      SCIPdebugMessage("Problem was found infeasible during preprocessing, therefore is not unbounded\n");
      return FALSE;
   }
   else if ( sdpi->allfixed )
   {
      SCIPdebugMessage("All variables were fixed during preprocessing, therefore the problem is not unbounded\n");
      return FALSE;
   }

   return SCIPsdpiSolverIsDualUnbounded(sdpi->sdpisolver);
}

/** returns TRUE iff SDP is proven to be dual infeasible;
 *  returns FALSE with a debug-message if the solver could not determine feasibility */
SCIP_Bool SCIPsdpiIsDualInfeasible(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   )
{
   assert( sdpi != NULL );
   CHECK_IF_SOLVED_BOOL(sdpi);

   if ( sdpi->infeasible )
   {
      SCIPdebugMessage("Problem was found infeasible during preprocessing\n");
      return TRUE;
   }
   else if ( sdpi->allfixed )
   {
      SCIPdebugMessage("All variables were fixed during preprocessing, solution is feasible\n");
      return FALSE;
   }

   return SCIPsdpiSolverIsDualInfeasible(sdpi->sdpisolver);
}

/** returns TRUE iff SDP is proven to be dual feasible;
 *  returns FALSE with a debug-message if the solver could not determine feasibility */
SCIP_Bool SCIPsdpiIsDualFeasible(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   )
{
   assert( sdpi != NULL );
   CHECK_IF_SOLVED_BOOL(sdpi);

   if ( sdpi->infeasible )
   {
      SCIPdebugMessage("Problem was found infeasible during preprocessing\n");
      return FALSE;
   }
   else if ( sdpi->allfixed )
   {
      SCIPdebugMessage("All variables fixed during preprocessing, solution is feasible\n");
      return TRUE;
   }

   return SCIPsdpiSolverIsDualFeasible(sdpi->sdpisolver);
}

/** returns TRUE iff the solver converged */
SCIP_Bool SCIPsdpiIsConverged(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   )
{
   assert( sdpi != NULL );
   CHECK_IF_SOLVED_BOOL(sdpi);

   if ( sdpi->infeasible )
   {
      SCIPdebugMessage("Problem was found infeasible during preprocessing, this counts as converged.\n");
      return TRUE;
   }
   else if ( sdpi->allfixed )
   {
      SCIPdebugMessage("All variables were fixed during preprocessing, this counts as converged.\n");
      return TRUE;
   }

   return SCIPsdpiSolverIsConverged(sdpi->sdpisolver);
}

/** returns TRUE iff the objective limit was reached */
SCIP_Bool SCIPsdpiIsObjlimExc(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   )
{
   assert( sdpi != NULL );
   CHECK_IF_SOLVED_BOOL(sdpi);

   if ( sdpi->infeasible )
   {
      SCIPdebugMessage("Problem was found infeasible during preprocessing, no objective limit available.\n");
      return FALSE;
   }
   else if ( sdpi->allfixed )
   {
      SCIPdebugMessage("All variables were fixed during preprocessing, no objective limit available.\n");
      return FALSE;
   }

   return SCIPsdpiSolverIsObjlimExc(sdpi->sdpisolver);
}

/** returns TRUE iff the iteration limit was reached */
SCIP_Bool SCIPsdpiIsIterlimExc(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   )
{
   assert( sdpi != NULL );
   CHECK_IF_SOLVED_BOOL(sdpi);

   if ( sdpi->infeasible )
   {
      SCIPdebugMessage("Problem was found infeasible during preprocessing, no iteration limit available.\n");
      return FALSE;
   }
   else if ( sdpi->allfixed )
   {
      SCIPdebugMessage("All variables were fixed during preprocessing, no iteration limit available.\n");
      return FALSE;
   }

   return SCIPsdpiSolverIsIterlimExc(sdpi->sdpisolver);
}

/** returns TRUE iff the time limit was reached */
SCIP_Bool SCIPsdpiIsTimelimExc(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   )
{
   assert( sdpi != NULL );

   if ( sdpi->infeasible )
   {
      SCIPdebugMessage("Problem was found infeasible during preprocessing, no time limit available.\n");
      return FALSE;
   }
   else if ( sdpi->allfixed )
   {
      SCIPdebugMessage("All variables were fixed during preprocessing, no time limit available.\n");
      return FALSE;
   }
   else if ( ! sdpi->solved )
   {
      SCIPdebugMessage("Problem was not solved, time limit not exceeded.\n");
      return FALSE;
   }

   return SCIPsdpiSolverIsTimelimExc(sdpi->sdpisolver);
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
int SCIPsdpiGetInternalStatus(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   )
{
   assert( sdpi != NULL );

   if ( ! sdpi->solved )
   {
      SCIPdebugMessage("Problem wasn't solved yet.\n");
      return -1;
   }
   else if ( sdpi->infeasible )
   {
      SCIPdebugMessage("Problem was found infeasible during preprocessing, no internal status available.\n");
      return 0;
   }
   else if ( sdpi->allfixed )
   {
      SCIPdebugMessage("All variables were fixed during preprocessing, no internal status available.\n");
      return 0;
   }

   return SCIPsdpiSolverGetInternalStatus(sdpi->sdpisolver);
}

/** returns TRUE iff SDP was solved to optimality, meaning the solver converged and returned primal and dual feasible solutions */
SCIP_Bool SCIPsdpiIsOptimal(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   )
{
   assert( sdpi != NULL );
   CHECK_IF_SOLVED_BOOL(sdpi);

   if ( sdpi->infeasible )
   {
      SCIPdebugMessage("Problem was found infeasible during preprocessing, therefore there is no optimal solution.\n");
      return FALSE;
   }
   else if ( sdpi->allfixed )
   {
      SCIPdebugMessage("All variables were fixed during preprocessing, therefore there is no optimal solution.\n");
      return FALSE;
   }

   return SCIPsdpiSolverIsOptimal(sdpi->sdpisolver);
}

/** returns TRUE iff SDP was solved to optimality or some other status was reached
 * that is still acceptable inside a Branch & Bound framework */
SCIP_Bool SCIPsdpiIsAcceptable(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   )
{
   assert( sdpi != NULL );

   if ( sdpi->infeasible )
   {
      SCIPdebugMessage("Problem was found infeasible during preprocessing, this is acceptable in a B&B context.\n");
      return TRUE;
   }
   else if ( sdpi->allfixed )
   {
      SCIPdebugMessage("All variables fixed during preprocessing, this is acceptable in a B&B context.\n");
      return TRUE;
   }
   else if ( ! sdpi->solved )
   {
      SCIPdebugMessage("Problem not solved succesfully, this is not acceptable in a B&B context.\n");
      return FALSE;
   }

   return SCIPsdpiSolverIsAcceptable(sdpi->sdpisolver);
}

/** gets objective value of solution */
SCIP_RETCODE SCIPsdpiGetObjval(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   SCIP_Real*            objval              /**< pointer to store the objective value */
   )
{
   assert( sdpi != NULL );
   assert( objval != NULL );
   CHECK_IF_SOLVED(sdpi);

   if ( sdpi->infeasible )
   {
      SCIPdebugMessage("Problem was found infeasible during preprocessing, no objective value available.\n");
      return SCIP_OKAY;
   }

   if ( sdpi->allfixed )
   {
      int v;

      /* As all variables were fixed during preprocessing, we have to compute it ourselves here */
      *objval = 0;

      for (v = 0; v < sdpi->nvars; v++)
         *objval += sdpi->lb[v] * sdpi->obj[v];

      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPsdpiSolverGetObjval(sdpi->sdpisolver, objval) );

   return SCIP_OKAY;
}

/** gets the best lower bound on the objective (this is equal to objval, if the problem was solved successfully, but can also give a bound
 *  if we did not get a feasible solution using the penalty approach) */
SCIP_RETCODE SCIPsdpiGetLowerObjbound(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   SCIP_Real*            objlb               /**< pointer to store the lower bound on the objective value */
   )
{
   assert( sdpi != NULL );
   assert( objlb != NULL );

   /* if we could successfully solve the problem, the best bound is the optimal objective */
   if ( sdpi->solved )
   {
      if ( sdpi->infeasible )
      {
         SCIPdebugMessage("Problem was found infeasible during preprocessing, no objective value available.\n");
         return SCIP_OKAY;
      }

      if ( sdpi->allfixed )
      {
         int v;

         /* As all variables were fixed during preprocessing, we have to compute it ourselves here */
         *objlb = 0;

         for (v = 0; v < sdpi->nvars; v++)
            *objlb += sdpi->lb[v] * sdpi->obj[v];

         return SCIP_OKAY;
      }

      SCIP_CALL( SCIPsdpiSolverGetObjval(sdpi->sdpisolver, objlb) );
      return SCIP_OKAY;
   }

   /* if we could not solve it, but tried the penalty formulation, we take the best bound computed by the penalty approach */
   if ( sdpi->penalty )
   {
      *objlb = sdpi->bestbound;
      return SCIP_OKAY;
   }

   /* if we could not solve it and did not use the penalty formulation (e.g. because the time limit was reached), we have no information */

   *objlb = -SCIPsdpiInfinity(sdpi);
   return SCIP_OKAY;
}

/** gets dual solution vector for feasible SDPs, if dualsollength isn't equal to the number of variables this will return the needed length and
 *  a debug message */
SCIP_RETCODE SCIPsdpiGetSol(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   SCIP_Real*            objval,             /**< pointer to store the objective value, may be NULL if not needed */
   SCIP_Real*            dualsol,            /**< pointer to store the dual solution vector, may be NULL if not needed */
   int*                  dualsollength       /**< length of the dualsol vector, must be 0 if dualsol is NULL, if this is less than the number
                                               *  of variables in the SDP, a debug-message will be thrown and this is set to the needed value */
   )
{
   assert( sdpi != NULL );
   assert( dualsollength != NULL );
   assert( *dualsollength == 0 || dualsol != NULL );
   CHECK_IF_SOLVED(sdpi);

   if ( sdpi->infeasible )
   {
      SCIPdebugMessage("Problem was found infeasible during preprocessing, no solution available.\n");
      return SCIP_OKAY;
   }
   else if ( sdpi->allfixed )
   {
      if ( objval != NULL )
      {
         SCIP_CALL( SCIPsdpiGetObjval(sdpi, objval) );
      }
      if ( *dualsollength > 0 )
      {
         int v;

         assert( dualsol != NULL );
         if ( *dualsollength < sdpi->nvars )
         {
            SCIPdebugMessage("The given array in SCIPsdpiGetSol only had length %d, but %d was needed", *dualsollength, sdpi->nvars);
            *dualsollength = sdpi->nvars;

            return SCIP_OKAY;
         }

         /* we give the fixed values as the solution */
         for (v = 0; v < sdpi->nvars; v++)
            dualsol[v] = sdpi->lb[v];

         return SCIP_OKAY;
      }
   }

   SCIP_CALL( SCIPsdpiSolverGetSol(sdpi->sdpisolver, objval, dualsol, dualsollength) );

   return SCIP_OKAY;
}

/** gets the primal variables corresponding to the lower and upper variable-bounds in the dual problem, the last input should specify the length
 *  of the arrays, if this is less than the number of variables, the needed length will be returned and a debug-message thrown
 *
 *  @note If a variable is either fixed or unbounded in the dual problem, a zero will be returned for the non-existent primal variable. */
SCIP_RETCODE SCIPsdpiGetPrimalBoundVars(
   SCIP_SDPI*            sdpi,               /**< pointer to an SDP-interface structure */
   SCIP_Real*            lbvars,             /**< pointer to store the values of the variables corresponding to lower bounds in the dual problems */
   SCIP_Real*            ubvars,             /**< pointer to store the values of the variables corresponding to upper bounds in the dual problems */
   int*                  arraylength         /**< input: length of lbvars and ubvars<br>
                                                  output: number of elements inserted into lbvars/ubvars (or needed length if it was not sufficient) */
   )
{
   assert( sdpi != NULL );
   assert( lbvars != NULL );
   assert( ubvars != NULL );
   assert( arraylength != NULL );
   assert( *arraylength >= 0 );
   CHECK_IF_SOLVED(sdpi);

   if ( sdpi->infeasible )
   {
      SCIPdebugMessage("Problem was found infeasible during preprocessing, no primal variables available.\n");
      return SCIP_OKAY;
   }
   else if ( sdpi->allfixed )
   {
      SCIPdebugMessage("All variables fixed during preprocessing, no primal variables available.\n");
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPsdpiSolverGetPrimalBoundVars(sdpi->sdpisolver, lbvars, ubvars, arraylength) );

   return SCIP_OKAY;
}

/** gets the number of SDP-iterations of the last solve call */
SCIP_RETCODE SCIPsdpiGetIterations(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int*                  iterations          /**< pointer to store the number of iterations of the last solve call */
   )
{
   assert( sdpi != NULL );
   assert( iterations != NULL );

   *iterations = sdpi->niterations;

   return SCIP_OKAY;
}

/** gets the number of calls to the SDP-solver for the last solve call */
SCIP_RETCODE SCIPsdpiGetSdpCalls(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int*                  calls               /**< pointer to store the number of calls to the SDP-solver for the last solve call */
   )
{
   assert( sdpi != NULL );
   assert( calls != NULL );

   *calls = sdpi->nsdpcalls;

   return SCIP_OKAY;
}

/** returns which settings the SDP-solver used in the last solve call */
SCIP_RETCODE SCIPsdpiSettingsUsed(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   SCIP_SDPSOLVERSETTING* usedsetting        /**< the setting used by the SDP-solver */
   )
{
   assert( sdpi != NULL );
   assert( usedsetting != NULL );

   if ( ! sdpi->solved )
   {
      SCIPdebugMessage("Problem was not solved successfully.\n");
      *usedsetting = SCIP_SDPSOLVERSETTING_UNSOLVED;
      return SCIP_OKAY;
   }
   else if ( sdpi->infeasible && ! sdpi->penalty ) /* if we solved the penalty formulation, we may also set infeasible if it is infeasible for the original problem */
   {
      SCIPdebugMessage("Problem was found infeasible during preprocessing, no settings used.\n");
      *usedsetting = SCIP_SDPSOLVERSETTING_UNSOLVED;
      return SCIP_OKAY;
   }
   else if ( sdpi->allfixed )
   {
      SCIPdebugMessage("All varialbes fixed during preprocessing, no settings used.\n");
      *usedsetting = SCIP_SDPSOLVERSETTING_UNSOLVED;
      return SCIP_OKAY;
   }
   else if ( sdpi->penalty )
   {
      *usedsetting = SCIP_SDPSOLVERSETTING_PENALTY;
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPsdpiSolverSettingsUsed(sdpi->sdpisolver, usedsetting) );
   return SCIP_OKAY;
}

/** returns which settings the SDP-solver used in the last solve call and whether primal and dual Slater condition were fullfilled */
SCIP_RETCODE SCIPsdpiSlaterSettings(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   SCIP_SDPSLATERSETTING* slatersetting      /**< the combination of Slater conditions and successfull settings */
   )
{
   SCIP_SDPSOLVERSETTING usedsetting;

   assert( sdpi != NULL );
   assert( slatersetting != NULL );

   if ( ! sdpi->solved )
   {
      SCIPdebugMessage("Problem was not solved successfully");
      if ( sdpi->bestbound > -SCIPsdpiSolverInfinity(sdpi->sdpisolver) )
      {
         SCIPdebugMessage(", but we could at least compute a lower bound. \n");
         if ( sdpi->dualslater == SCIP_SDPSLATER_INF)
            *slatersetting = SCIP_SDPSLATERSETTING_BOUNDEDINFEASIBLE;
         else
         {
            switch( sdpi->primalslater )/*lint --e{788}*/
            {
               case SCIP_SDPSLATER_NOINFO:
                  if ( sdpi->dualslater == SCIP_SDPSLATER_NOT )
                     *slatersetting = SCIP_SDPSLATERSETTING_BOUNDEDNOSLATER;
                  else
                     *slatersetting = SCIP_SDPSLATERSETTING_NOINFO;
                  break;
               case SCIP_SDPSLATER_NOT:
                  *slatersetting = SCIP_SDPSLATERSETTING_BOUNDEDNOSLATER;
                  break;
               case SCIP_SDPSLATER_HOLDS:
                  switch( sdpi->dualslater )/*lint --e{788}*/
                  {
                     case SCIP_SDPSLATER_NOINFO:
                        *slatersetting = SCIP_SDPSLATERSETTING_NOINFO;
                        break;
                     case SCIP_SDPSLATER_NOT:
                        *slatersetting = SCIP_SDPSLATERSETTING_BOUNDEDNOSLATER;
                        break;
                     case SCIP_SDPSLATER_HOLDS:
                        *slatersetting = SCIP_SDPSLATERSETTING_BOUNDEDWSLATER;
                        break;
                     default:
                        *slatersetting = SCIP_SDPSLATERSETTING_NOINFO;
                        break;
                  }
                  break;
                  default:
                     *slatersetting = SCIP_SDPSLATERSETTING_NOINFO;
                     break;
            }
         }
      }
      else
      {
         SCIPdebugMessage(".\n");
         if ( sdpi->dualslater == SCIP_SDPSLATER_INF)
            *slatersetting = SCIP_SDPSLATERSETTING_UNSOLVEDINFEASIBLE;
         else
         {
            switch( sdpi->primalslater )/*lint --e{788}*/
            {
               case SCIP_SDPSLATER_NOINFO:
                  if ( sdpi->dualslater == SCIP_SDPSLATER_NOT )
                     *slatersetting = SCIP_SDPSLATERSETTING_UNSOLVEDNOSLATER;
                  else
                     *slatersetting = SCIP_SDPSLATERSETTING_NOINFO;
                  break;
               case SCIP_SDPSLATER_NOT:
                  *slatersetting = SCIP_SDPSLATERSETTING_UNSOLVEDNOSLATER;
                  break;
               case SCIP_SDPSLATER_HOLDS:
                  switch( sdpi->dualslater )/*lint --e{788}*/
                  {
                     case SCIP_SDPSLATER_NOINFO:
                        *slatersetting = SCIP_SDPSLATERSETTING_NOINFO;
                        break;
                     case SCIP_SDPSLATER_NOT:
                        *slatersetting = SCIP_SDPSLATERSETTING_UNSOLVEDNOSLATER;
                        break;
                     case SCIP_SDPSLATER_HOLDS:
                        *slatersetting = SCIP_SDPSLATERSETTING_UNSOLVEDWSLATER;
                        break;
                     default:
                        *slatersetting = SCIP_SDPSLATERSETTING_NOINFO;
                        break;
                  }
                  break;
               default:
                  *slatersetting = SCIP_SDPSLATERSETTING_NOINFO;
                  break;
            }
         }
      }
      return SCIP_OKAY;
   }
   else if ( sdpi->infeasible && ( ! sdpi->penalty ) ) /* if we solved the penalty formulation, we may also set infeasible if it is infeasible for the original problem */
   {
      SCIPdebugMessage("Problem was found infeasible during preprocessing, no settings used.\n");
      *slatersetting = SCIP_SDPSLATERSETTING_NOINFO;
      return SCIP_OKAY;
   }
   else if ( sdpi->allfixed )
   {
      SCIPdebugMessage("All varialbes fixed during preprocessing, no settings used.\n");
      *slatersetting = SCIP_SDPSLATERSETTING_NOINFO;
      return SCIP_OKAY;
   }
   else if ( sdpi->penalty )
   {
      switch( sdpi->primalslater )/*lint --e{788}*/
      {
         case SCIP_SDPSLATER_NOINFO:
            if ( sdpi->dualslater == SCIP_SDPSLATER_NOT )
               *slatersetting = SCIP_SDPSLATERSETTING_PENALTYNOSLATER;
            else if ( sdpi->dualslater == SCIP_SDPSLATER_INF )
               *slatersetting = SCIP_SDPSLATERSETTING_PENALTYINFEASIBLE;
            else
               *slatersetting = SCIP_SDPSLATERSETTING_NOINFO;
            break;
         case SCIP_SDPSLATER_NOT:
            if ( sdpi->dualslater == SCIP_SDPSLATER_INF )
               *slatersetting = SCIP_SDPSLATERSETTING_PENALTYINFEASIBLE;
            else
               *slatersetting = SCIP_SDPSLATERSETTING_PENALTYNOSLATER;
            break;
         case SCIP_SDPSLATER_HOLDS:
            switch( sdpi->dualslater )/*lint --e{788}*/
            {
               case SCIP_SDPSLATER_NOINFO:
                  *slatersetting = SCIP_SDPSLATERSETTING_NOINFO;
                  break;
               case SCIP_SDPSLATER_NOT:
                  *slatersetting = SCIP_SDPSLATERSETTING_PENALTYNOSLATER;
                  break;
               case SCIP_SDPSLATER_HOLDS:
                  *slatersetting = SCIP_SDPSLATERSETTING_PENALTYWSLATER;
                  break;
               case SCIP_SDPSLATER_INF:
                  *slatersetting = SCIP_SDPSLATERSETTING_PENALTYINFEASIBLE;
                  break;
               default:
                  *slatersetting = SCIP_SDPSLATERSETTING_NOINFO;
                  break;
            }
            break;
         default:
            *slatersetting = SCIP_SDPSLATERSETTING_NOINFO;
            break;
      }
      return SCIP_OKAY;
   }

   switch( sdpi->primalslater )/*lint --e{788}*/
   {
   case SCIP_SDPSLATER_NOINFO:
      if ( sdpi->dualslater == SCIP_SDPSLATER_NOT )
      {
         usedsetting = SCIP_SDPSOLVERSETTING_UNSOLVED;
         SCIP_CALL( SCIPsdpiSolverSettingsUsed(sdpi->sdpisolver, &usedsetting) );
         switch( usedsetting )/*lint --e{788}*/
         {
            case SCIP_SDPSOLVERSETTING_FAST:
               *slatersetting = SCIP_SDPSLATERSETTING_STABLENOSLATER;
               break;
            case SCIP_SDPSOLVERSETTING_MEDIUM:
            case SCIP_SDPSOLVERSETTING_STABLE:
               *slatersetting = SCIP_SDPSLATERSETTING_UNSTABLENOSLATER;
               break;
            default:
               *slatersetting = SCIP_SDPSLATERSETTING_NOINFO;
               break;
         }
      }
      if ( sdpi->dualslater == SCIP_SDPSLATER_INF )
      {
         usedsetting = SCIP_SDPSOLVERSETTING_UNSOLVED;
         SCIP_CALL( SCIPsdpiSolverSettingsUsed(sdpi->sdpisolver, &usedsetting) );
         switch( usedsetting )/*lint --e{788}*/
         {
            case SCIP_SDPSOLVERSETTING_FAST:
               *slatersetting = SCIP_SDPSLATERSETTING_STABLEINFEASIBLE;
               break;
            case SCIP_SDPSOLVERSETTING_MEDIUM:
            case SCIP_SDPSOLVERSETTING_STABLE:
               *slatersetting = SCIP_SDPSLATERSETTING_UNSTABLEINFEASIBLE;
               break;
            default:
               *slatersetting = SCIP_SDPSLATERSETTING_NOINFO;
               break;
         }
      }
      else
         *slatersetting = SCIP_SDPSLATERSETTING_NOINFO;
      break;
   case SCIP_SDPSLATER_NOT:
      if ( sdpi->dualslater == SCIP_SDPSLATER_INF )
      {
         usedsetting = SCIP_SDPSOLVERSETTING_UNSOLVED;
         SCIP_CALL( SCIPsdpiSolverSettingsUsed(sdpi->sdpisolver, &usedsetting) );
         switch( usedsetting )/*lint --e{788}*/
         {
            case SCIP_SDPSOLVERSETTING_FAST:
               *slatersetting = SCIP_SDPSLATERSETTING_STABLEINFEASIBLE;
               break;
            case SCIP_SDPSOLVERSETTING_MEDIUM:
            case SCIP_SDPSOLVERSETTING_STABLE:
               *slatersetting = SCIP_SDPSLATERSETTING_UNSTABLEINFEASIBLE;
               break;
            default:
               *slatersetting = SCIP_SDPSLATERSETTING_NOINFO;
               break;
         }
      }
      else
      {
         usedsetting = SCIP_SDPSOLVERSETTING_UNSOLVED;
         SCIP_CALL( SCIPsdpiSolverSettingsUsed(sdpi->sdpisolver, &usedsetting) );
         switch( usedsetting )/*lint --e{788}*/
         {
            case SCIP_SDPSOLVERSETTING_FAST:
               *slatersetting = SCIP_SDPSLATERSETTING_STABLENOSLATER;
               break;
            case SCIP_SDPSOLVERSETTING_MEDIUM:
            case SCIP_SDPSOLVERSETTING_STABLE:
               *slatersetting = SCIP_SDPSLATERSETTING_UNSTABLENOSLATER;
               break;
            default:
               *slatersetting = SCIP_SDPSLATERSETTING_NOINFO;
               break;
         }
      }
      break;
   case SCIP_SDPSLATER_HOLDS:
      switch( sdpi->dualslater )/*lint --e{788}*/
      {
      case SCIP_SDPSLATER_NOINFO:
         *slatersetting = SCIP_SDPSLATERSETTING_NOINFO;
         break;
      case SCIP_SDPSLATER_NOT:
         usedsetting = SCIP_SDPSOLVERSETTING_UNSOLVED;
         SCIP_CALL( SCIPsdpiSolverSettingsUsed(sdpi->sdpisolver, &usedsetting) );
         switch( usedsetting )/*lint --e{788}*/
         {
            case SCIP_SDPSOLVERSETTING_FAST:
               *slatersetting = SCIP_SDPSLATERSETTING_STABLENOSLATER;
               break;
            case SCIP_SDPSOLVERSETTING_MEDIUM:
            case SCIP_SDPSOLVERSETTING_STABLE:
               *slatersetting = SCIP_SDPSLATERSETTING_UNSTABLENOSLATER;
               break;
            default:
               *slatersetting = SCIP_SDPSLATERSETTING_NOINFO;
               break;
         }
         break;
         case SCIP_SDPSLATER_INF:
            usedsetting = SCIP_SDPSOLVERSETTING_UNSOLVED;
            SCIP_CALL( SCIPsdpiSolverSettingsUsed(sdpi->sdpisolver, &usedsetting) );
            switch( usedsetting )/*lint --e{788}*/
            {
               case SCIP_SDPSOLVERSETTING_FAST:
                  *slatersetting = SCIP_SDPSLATERSETTING_STABLEINFEASIBLE;
                  break;
               case SCIP_SDPSOLVERSETTING_MEDIUM:
               case SCIP_SDPSOLVERSETTING_STABLE:
                  *slatersetting = SCIP_SDPSLATERSETTING_UNSTABLEINFEASIBLE;
                  break;
               default:
                  *slatersetting = SCIP_SDPSLATERSETTING_NOINFO;
                  break;
            }
            break;
         case SCIP_SDPSLATER_HOLDS:
            usedsetting = SCIP_SDPSOLVERSETTING_UNSOLVED;
            SCIP_CALL( SCIPsdpiSolverSettingsUsed(sdpi->sdpisolver, &usedsetting) );
            switch( usedsetting )/*lint --e{788}*/
            {
               case SCIP_SDPSOLVERSETTING_FAST:
                  *slatersetting = SCIP_SDPSLATERSETTING_STABLEWSLATER;
                  break;
               case SCIP_SDPSOLVERSETTING_MEDIUM:
               case SCIP_SDPSOLVERSETTING_STABLE:
                  *slatersetting = SCIP_SDPSLATERSETTING_UNSTABLEWSLATER;
                  break;
               default:
                  *slatersetting = SCIP_SDPSLATERSETTING_NOINFO;
                  break;
            }
            break;
         default:
            *slatersetting = SCIP_SDPSLATERSETTING_NOINFO;
            break;
      }
      break;
   default:
      *slatersetting = SCIP_SDPSLATERSETTING_NOINFO;
      break;
   }

   return SCIP_OKAY;
}

/** returns whether primal and dual Slater condition held for last solved SDP */
SCIP_RETCODE SCIPsdpiSlater(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   SCIP_SDPSLATER*       primalslater,       /**< pointer to save whether primal Slater condition held */
   SCIP_SDPSLATER*       dualslater          /**< pointer to save whether dual Slater condition held */
   )
{
   assert( sdpi != NULL );
   assert( primalslater != NULL );
   assert( dualslater != NULL );


   if ( sdpi->infeasible )
   {
      *primalslater = SCIP_SDPSLATER_NOINFO;
      *dualslater = sdpi->dualslater;
      return SCIP_OKAY;
   }

   if (sdpi->allfixed )
   {
      *primalslater = SCIP_SDPSLATER_NOINFO;
      *dualslater = SCIP_SDPSLATER_NOINFO;
      return SCIP_OKAY;
   }

   *primalslater = sdpi->primalslater;
   *dualslater = sdpi->dualslater;

   return SCIP_OKAY;
}

/**@} */




/*
 * Numerical Methods
 */

/**@name Numerical Methods */
/**@{ */

/** returns value treated as infinity in the SDP-solver */
SCIP_Real SCIPsdpiInfinity(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   )
{
   assert( sdpi != NULL  );

   return SCIPsdpiSolverInfinity(sdpi->sdpisolver);
}

/** checks if given value is treated as (plus or minus) infinity in the SDP-solver */
SCIP_Bool SCIPsdpiIsInfinity(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   SCIP_Real             val                 /**< value to be checked for infinity */
   )
{
   assert( sdpi != NULL );

   return ((val <= -SCIPsdpiInfinity(sdpi)) || (val >= SCIPsdpiInfinity(sdpi)));
}

/** gets floating point parameter of SDP-interface */
SCIP_RETCODE SCIPsdpiGetRealpar(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   SCIP_SDPPARAM         type,               /**< parameter number */
   SCIP_Real*            dval                /**< pointer to store the parameter value */
   )
{
   assert( sdpi != NULL );
   assert( sdpi->sdpisolver != NULL );
   assert( dval != NULL );

   switch( type )/*lint --e{788}*/
   {
   case SCIP_SDPPAR_EPSILON:
      *dval = sdpi->epsilon;
      break;
   case SCIP_SDPPAR_GAPTOL:
      *dval = sdpi->gaptol;
      break;
   case SCIP_SDPPAR_FEASTOL:
      *dval = sdpi->feastol;
      break;
   case SCIP_SDPPAR_SDPSOLVERFEASTOL:
      SCIP_CALL_PARAM( SCIPsdpiSolverGetRealpar(sdpi->sdpisolver, type, dval) );
      break;
   case SCIP_SDPPAR_OBJLIMIT:
      SCIP_CALL_PARAM( SCIPsdpiSolverGetRealpar(sdpi->sdpisolver, type, dval) );
      break;
   case SCIP_SDPPAR_PENALTYPARAM:
      *dval = sdpi->penaltyparam;
      break;
   case SCIP_SDPPAR_MAXPENALTYPARAM:
      *dval = sdpi->maxpenaltyparam;
      break;
   case SCIP_SDPPAR_LAMBDASTAR:
      SCIP_CALL_PARAM( SCIPsdpiSolverGetRealpar(sdpi->sdpisolver, type, dval) );
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }

   return SCIP_OKAY;
}

/** sets floating point parameter of SDP-interface */
SCIP_RETCODE SCIPsdpiSetRealpar(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   SCIP_SDPPARAM         type,               /**< parameter number */
   SCIP_Real             dval                /**< parameter value */
   )
{
   assert( sdpi != NULL );

   switch( type )/*lint --e{788}*/
   {
   case SCIP_SDPPAR_EPSILON:
      sdpi->epsilon = dval;
      SCIP_CALL_PARAM( SCIPsdpiSolverSetRealpar(sdpi->sdpisolver, type, dval) );
      break;
   case SCIP_SDPPAR_GAPTOL:
      sdpi->gaptol = dval;
      SCIP_CALL_PARAM( SCIPsdpiSolverSetRealpar(sdpi->sdpisolver, type, dval) );
      break;
   case SCIP_SDPPAR_FEASTOL:
      sdpi->feastol = dval;
      SCIP_CALL_PARAM( SCIPsdpiSolverSetRealpar(sdpi->sdpisolver, type, dval) );
      break;
   case SCIP_SDPPAR_SDPSOLVERFEASTOL:
      SCIP_CALL_PARAM( SCIPsdpiSolverSetRealpar(sdpi->sdpisolver, type, dval) );
      break;
   case SCIP_SDPPAR_OBJLIMIT:
      SCIP_CALL_PARAM( SCIPsdpiSolverSetRealpar(sdpi->sdpisolver, type, dval) );
      break;
   case SCIP_SDPPAR_PENALTYPARAM:
      sdpi->penaltyparam = dval;
      SCIP_CALL_PARAM_IGNORE_UNKNOWN( SCIPsdpiSolverSetRealpar(sdpi->sdpisolver, type, dval) );
      break;
   case SCIP_SDPPAR_MAXPENALTYPARAM:
      sdpi->maxpenaltyparam = dval;
      break;
   case SCIP_SDPPAR_LAMBDASTAR:
      SCIP_CALL_PARAM( SCIPsdpiSolverSetRealpar(sdpi->sdpisolver, type, dval) );
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }

   return SCIP_OKAY;
}

/** gets integer parameter of SDP-interface */
SCIP_RETCODE SCIPsdpiGetIntpar(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   SCIP_SDPPARAM         type,               /**< parameter number */
   int*                  ival                /**< pointer to store the parameter value */
   )
{
   assert( sdpi != NULL );
   assert( sdpi->sdpisolver != NULL );
   assert( ival != NULL );

   switch( type )/*lint --e{788}*/
   {
   case SCIP_SDPPAR_SDPINFO:
   case SCIP_SDPPAR_NTHREADS:
      SCIP_CALL_PARAM( SCIPsdpiSolverGetIntpar(sdpi->sdpisolver, type, ival) );
      break;
   case SCIP_SDPPAR_SLATERCHECK:
      *ival = sdpi->slatercheck;
      break;
   case SCIP_SDPPAR_NPENALTYINCR:
      *ival = sdpi->npenaltyincr;
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }

   return SCIP_OKAY;
}

/** sets integer parameter of SDP-interface */
SCIP_RETCODE SCIPsdpiSetIntpar(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   SCIP_SDPPARAM         type,               /**< parameter number */
   int                   ival                /**< parameter value */
   )
{
   assert( sdpi != NULL );
   assert( sdpi->sdpisolver != NULL );

   switch( type )/*lint --e{788}*/
   {
   case SCIP_SDPPAR_SDPINFO:
      assert( ival == 0 || ival == 1 ); /* this is a boolean parameter */
      SCIP_CALL_PARAM( SCIPsdpiSolverSetIntpar(sdpi->sdpisolver, type, ival) );
      break;
   case SCIP_SDPPAR_NTHREADS:
      SCIP_CALL_PARAM( SCIPsdpiSolverSetIntpar(sdpi->sdpisolver, type, ival) );
      break;
   case SCIP_SDPPAR_SLATERCHECK:
      sdpi->slatercheck = ival;
      break;
   case SCIP_SDPPAR_NPENALTYINCR:
      sdpi->npenaltyincr = ival;
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }

   return SCIP_OKAY;
}

/** compute and set lambdastar (only used for SDPA) */
SCIP_RETCODE SCIPsdpiComputeLambdastar(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   SCIP_Real             maxguess            /**< maximum guess for lambda star of all SDP-constraints */
   )
{
   return SCIPsdpiSolverComputeLambdastar(sdpi->sdpisolver, maxguess);
}

/** compute and set the penalty parameter */
SCIP_RETCODE SCIPsdpiComputePenaltyparam(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   SCIP_Real             maxcoeff,           /**< maximum objective coefficient */
   SCIP_Real*            penaltyparam        /**< the computed penalty parameter */
   )
{
   SCIP_CALL( SCIPsdpiSolverComputePenaltyparam(sdpi->sdpisolver, maxcoeff, penaltyparam) );

   sdpi->penaltyparam = *penaltyparam;

   return SCIP_OKAY;
}

/** compute and set the maximum penalty parameter */
SCIP_RETCODE SCIPsdpiComputeMaxPenaltyparam(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   SCIP_Real             penaltyparam,       /**< the initial penalty parameter */
   SCIP_Real*            maxpenaltyparam     /**< the computed maximum penalty parameter */
   )
{
   SCIP_CALL( SCIPsdpiSolverComputeMaxPenaltyparam(sdpi->sdpisolver, penaltyparam, maxpenaltyparam) );

   sdpi->maxpenaltyparam = *maxpenaltyparam;

   /* if the initial penalty parameter is smaller than the maximum one, we decrease the initial correspondingly */
   /* if the maximum penalty parameter is smaller than the initial penalty paramater, we decrease the initial one correspondingly */
   if ( sdpi->penaltyparam > *maxpenaltyparam )
   {
      SCIPdebugMessage("Decreasing penaltyparameter of %f to maximum penalty paramater of %f.\n", sdpi->penaltyparam, *maxpenaltyparam);
      sdpi->penaltyparam = *maxpenaltyparam;
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
SCIP_RETCODE SCIPsdpiReadSDP(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   const char*           fname               /**< file name */
   )
{/*lint --e{715}*/
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_LPERROR;
}

/** writes SDP to a file */
SCIP_RETCODE SCIPsdpiWriteSDP(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   const char*           fname               /**< file name */
   )
{/*lint --e{715}*/
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_LPERROR;
}

/**@} */
