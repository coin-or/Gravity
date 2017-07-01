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

/**@file   sdpsolchecker.c
 * @brief  checks a given SDP solution for feasibility
 * @author Tristan Gally
 */

/*#define SCIP_DEBUG*/

#include "sdpi/sdpsolchecker.h"
#include "sdpi/lapack.h"

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

/** Given a solution, an SDP instance and a feasibility tolerance, checks whether
 *  the smallest eigenvalue is >= -feastol for a given feasibility tolerance.
 */
SCIP_RETCODE SCIPsdpSolcheckerCheck(
   BMS_BUFMEM*           bufmem,             /**< buffer memory */
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
   SCIP_Real*            solvector,          /**< values of all variables (including fixed ones) in the solution that should be checked */
   SCIP_Real             feastol,            /**< feasibility tolerance to check feasibility for */
   SCIP_Real             epsilon,            /**< tolerance used to check for fixed variables */
   SCIP_Bool*            infeasible          /**< pointer to store whether solution is feasible */
)
{/*lint --e{818}*/
   int i;
   int j;
   int b;
   int v;
   int ind;

   assert( bufmem != NULL );
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
   assert( nlpcons >= 0 );
   assert( noldlpcons >= nlpcons );
   assert( nlpcons == 0 || lplhs != NULL );
   assert( nlpcons == 0 || lprhs != NULL );
   assert( nlpcons == 0 || rownactivevars != NULL );
   assert( lpnnonz >= 0 );
   assert( nlpcons == 0 || lprow != NULL );
   assert( nlpcons == 0 || lpcol != NULL );
   assert( nlpcons == 0 || lpval != NULL );
   assert( solvector != NULL );
   assert( feastol >= 0 );
   assert( infeasible != NULL );

   /* check variable bounds */
   for (i = 0; i < nvars; i++)
   {
      if ( solvector[i] < lb[i] - feastol || solvector[i] > ub[i] + feastol )
      {
         SCIPdebugMessage("solution found infeasible (feastol=%f) for variable bounds: x[%d] = %f <|= [%f, %f]\n",
               feastol, i, solvector[i], lb[i], ub[i]);
         *infeasible = TRUE;
         return SCIP_OKAY;
      }
   }

   /* check linear constraints (since DSDP sorts the lp-arrays by cols, we cannot expect them to be sorted) */
   if ( nlpcons > 0 )
   {
      SCIP_Real* lpconsvals;

      BMS_CALL( BMSallocBufferMemoryArray(bufmem, &lpconsvals, noldlpcons) );

      /* initialize all rows with zero */
      for (i = 0; i < noldlpcons; i++)
         lpconsvals[i] = 0;

      /* compute the values of all rows */
      for (i = 0; i < lpnnonz; i++)
      {
         if ( lb[lpcol[i]] < ub[lpcol[i]] - epsilon ) /* fixed variables are already included in lhs/rhs */
            lpconsvals[lprow[i]] += solvector[lpcol[i]] * lpval[i];
      }

      /* check all active constraints for feasibility */
      ind = 0; /* used to iterate over active constraints */
      for (i = 0; i < noldlpcons; i++)
      {
         if ( rownactivevars[i] > 1 )
         {
            if ( lpconsvals[i] < lplhs[ind] - feastol || lpconsvals[i] > lprhs[ind] + feastol)
            {
               SCIPdebugMessage("solution found infeasible (feastol=%f) for lp constraint: LP-%d = %f <|= [%f,%f]\n",
                     feastol, i, lpconsvals[i], lplhs[ind], lprhs[ind]);
               BMSfreeBufferMemoryArray(bufmem, &lpconsvals);
               *infeasible = TRUE;
               return SCIP_OKAY;
            }

            ind++;
         }
      }
      BMSfreeBufferMemoryArray(bufmem, &lpconsvals);
   }

   /* check sdp constraints */
   if ( nsdpblocks > 0 )
   {
      SCIP_Real* fullsdpmatrix;
      SCIP_Real eigenvalue;
      int maxblocksize = 0;

      /* allocate memory */
      if ( nsdpblocks == 1 )
         maxblocksize = sdpblocksizes[0] - nremovedinds[0];
      else
      {
         /* calculate maximum size of any SDP block to not have to reallocate memory in between */
         for (b = 0; b < nsdpblocks; b++)
         {
            maxblocksize = ((sdpblocksizes[b] - nremovedinds[b]) > maxblocksize) ? sdpblocksizes[b] - nremovedinds[b] : maxblocksize;
         }
      }

      BMS_CALL( BMSallocBufferMemoryArray(bufmem, &fullsdpmatrix, maxblocksize * maxblocksize) );/*lint !e647*/

      for (b = 0; b < nsdpblocks; b++)
      {
         if ( blockindchanges[b] > -1 )
         {
            /* initialize lower triangular part of fullsdpmatrix with zero */
            for (i = 0; i < sdpblocksizes[b] - nremovedinds[b]; i++)
            {
               for (j = 0; j <= i; j++)
               {
                  fullsdpmatrix[i * (sdpblocksizes[b] - nremovedinds[b]) + j] = 0.0; /*lint !e679*/
               }
            }

            /* iterate over all non-fixed variables and add the corresponding nonzeros */
            for (v = 0; v < sdpnblockvars[b]; v++)
            {
               if ( lb[sdpvar[b][v]] < ub[sdpvar[b][v]] - epsilon )
               {
                  for (i = 0; i < sdpnblockvarnonz[b][v]; i++)
                  {
                     fullsdpmatrix[((sdprow[b][v][i] - indchanges[b][sdprow[b][v][i]]) * (sdpblocksizes[b] - nremovedinds[b])) +
                                   sdpcol[b][v][i] - indchanges[b][sdpcol[b][v][i]]] += solvector[sdpvar[b][v]] * sdpval[b][v][i];/*lint !e679*/
                  }
               }
            }

            /* add constant matrix */
            for (i = 0; i < sdpconstnblocknonz[b]; i++)
            {
               fullsdpmatrix[((sdpconstrow[b][i] - indchanges[b][sdpconstrow[b][i]]) * (sdpblocksizes[b] - nremovedinds[b])) +
                             sdpconstcol[b][i] - indchanges[b][sdpconstcol[b][i]]] -= sdpconstval[b][i];/*lint !e679*/
            }

            /* extend to full symmetric matrix for LAPACK */
            for (i = 0; i < sdpblocksizes[b] - nremovedinds[b]; i++)
            {
               for (j = 0; j < i; j++)
               {
                  fullsdpmatrix[j * (sdpblocksizes[b] - nremovedinds[b]) + i] = fullsdpmatrix[i * (sdpblocksizes[b] - nremovedinds[b]) + j];/*lint !e679*/
               }
            }

            /* compute smallest eigenvalue using LAPACK */
            SCIP_CALL( SCIPlapackComputeIthEigenvalue(bufmem, FALSE, sdpblocksizes[b] - nremovedinds[b], fullsdpmatrix, 1, &eigenvalue, NULL) );

            if ( eigenvalue < - feastol )
            {
               SCIPdebugMessage("solution found infeasible (feastol=%.10f) for sdp constraint %d, smallest eigenvector %.10f\n",
                     feastol, b, eigenvalue);
               BMSfreeBufferMemoryArray(bufmem, &fullsdpmatrix);
               *infeasible = TRUE;
               return SCIP_OKAY;
            }
         }
      }

      BMSfreeBufferMemoryArray(bufmem, &fullsdpmatrix);
   }

   *infeasible = FALSE;
   return SCIP_OKAY;
}
