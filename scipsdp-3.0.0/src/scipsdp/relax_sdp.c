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

/**@file   relax_sdp.c
 * @ingroup RELAXATORS
 * @brief  SDP-relaxator
 * @author Sonja Mars
 * @author Tristan Gally
 */

/* #define SCIP_DEBUG*/
/* #define SCIP_MORE_DEBUG   *//* displays complete solution for each relaxation */
/* #define SCIP_EVEN_MORE_DEBUG  *//* shows number of deleted empty cols/rows for every relaxation and variable status &
 * bounds as well as all constraints in the beginning */
/* #define SLATERSOLVED_ABSOLUTE *//* uncomment this to return the absolute number of nodes for, e.g., solved fast with slater in addition to percentages */

#include "relax_sdp.h"

#include "assert.h"                     /*lint !e451*/
#include "string.h"                     /* for strcmp */

#include "SdpVarmapper.h"
#include "sdpi/sdpi.h"
#include "scipsdp/cons_sdp.h"
#include "scipsdp/cons_savedsdpsettings.h"

/* turn off lint warnings for whole file: */
/*lint --e{788,818}*/

#define RELAX_NAME                  "SDP"
#define RELAX_DESC                  "SDP-relaxator"
#define RELAX_PRIORITY              1
#define RELAX_FREQ                  1


/* default values for parameters: */
#define DEFAULT_SDPSOLVERGAPTOL     1e-4     /**< the stopping criterion for the duality gap the sdpsolver should use */
#define DEFAULT_PENALTYPARAM        -1.0     /**< the penalty parameter Gamma used for the penalty formulation if the SDP solver didn't converge */
#define DEFAULT_LAMBDASTAR          -1.0     /**< the parameter lambda star used by SDPA to set the initial point */
#define DEFAULT_MAXPENALTYPARAM     -1.0     /**< the penalty parameter Gamma used for the penalty formulation if the SDP solver didn't converge */
#define DEFAULT_SLATERCHECK         0        /**< Should the Slater condition be checked ? */
#define DEFAULT_OBJLIMIT            FALSE    /**< Should an objective limit be given to the SDP-Solver ? */
#define DEFAULT_RESOLVE             TRUE     /**< Are we allowed to solve the relaxation of a single node multiple times in a row (outside of probing) ? */
#define DEFAULT_TIGHTENVB           TRUE     /**< Should Big-Ms in varbound-like constraints be tightened before giving them to the SDP-solver ? */
#define DEFAULT_SDPINFO             FALSE    /**< Should the SDP solver output information to the screen? */
#define DEFAULT_DISPLAYSTAT         TRUE     /**< Should statistics about SDP iterations and solver settings/success be printed after quitting SCIP-SDP ? */
#define DEFAULT_SETTINGSRESETFREQ   -1       /**< frequency for resetting parameters in SDP solver and trying again with fastest settings */
#define DEFAULT_SETTINGSRESETOFS    0        /**< frequency offset for resetting parameters in SDP solver and trying again with fastest settings */
#define DEFAULT_SDPSOLVERTHREADS    -1       /**< number of threads the SDP solver should use, currently only supported for MOSEK (-1 = number of cores) */

/*
 * Data structures
 */

/** relaxator data */
struct SCIP_RelaxData
{
   SCIP_SDPI*            sdpi;               /**< general SDP Interface that is given the data to presolve the SDP and give it so a solver specific interface */
   SdpVarmapper*         varmapper;          /**< maps SCIP variables to their global SDP indices and vice versa */
   SCIP_Real             objval;             /**< objective value of the last SDP-relaxation */
   SCIP_Bool             origsolved;         /**< solved original problem to optimality (not only a penalty or probing formulation) */
   SCIP_Bool             probingsolved;      /**< was the last probing SDP solved successfully? */
   SCIP_Real             sdpsolvergaptol;    /**< the stopping criterion for the duality gap the sdpsolver should use */
   SCIP_Real             sdpsolverfeastol;   /**< the feasibility tolerance the SDP solver should use for the SDP constraints */
   SCIP_Real             penaltyparam;       /**< the starting penalty parameter Gamma used for the penalty formulation if the SDP solver didn't converge */
   SCIP_Real             maxpenaltyparam;    /**< the maximum penalty parameter Gamma used for the penalty formulation if the SDP solver didn't converge */
   int                   npenaltyincr;       /**< maximum number of times the penalty parameter will be increased if penalty formulation failed */
   SCIP_Real             lambdastar;         /**< the parameter lambda star used by SDPA to set the initial point */
   int                   sdpiterations;      /**< saves the total number of sdp-iterations */
   int                   solvedfast;         /**< number of SDPs solved with fast settings */
   int                   solvedmedium;       /**< number of SDPs solved with medium settings */
   int                   solvedstable;       /**< number of SDPs solved with stable settings */
   int                   solvedpenalty;      /**< number of SDPs solved using penalty formulation */
   int                   unsolved;           /**< number of SDPs that could not be solved even using a penalty formulation */
   int                   slatercheck;        /**< Should the Slater condition for the dual problem be check ahead of solving every SDP ? */
   SCIP_Bool             sdpinfo;            /**< Should the SDP solver output information to the screen? */
   SCIP_Bool             displaystat;        /**< Should statistics about SDP iterations and solver settings/success be printed after quitting SCIP-SDP ? */
   SCIP_Bool             objlimit;           /**< Should an objective limit be given to the SDP solver? */
   SCIP_Bool             resolve;            /**< Are we allowed to solve the relaxation of a single node multiple times in a row (outside of probing) ? */
   SCIP_Bool             tightenvb;          /**< Should Big-Ms in varbound-like constraints be tightened before giving them to the SDP-solver ? */
   int                   settingsresetfreq;  /**< frequency for resetting parameters in SDP solver and trying again with fastest settings */
   int                   settingsresetofs;   /**< frequency offset for resetting parameters in SDP solver and trying again with fastest settings */
   int                   sdpsolverthreads;   /**< number of threads the SDP solver should use, currently only supported for MOSEK (-1 = number of cores) */
   int                   sdpcalls;           /**< number of solved SDPs (used to compute average SDP iterations), different settings tried are counted as multiple calls */
   int                   sdpinterfacecalls;  /**< number of times the SDP interfaces was called (used to compute slater statistics) */
   long int              lastsdpnode;        /**< number of the SCIP node the current SDP-solution belongs to */
   SCIP_Bool             feasible;           /**< was the last solved SDP feasible */
   int                   stablewslater;      /**< number of instances solved with fastest settings where primal and dual slater held */
   int                   unstablewslater;    /**< number of instances solved with stable settings where primal and dual slater held */
   int                   penaltywslater;     /**< number of instances solved with penalty formulation where primal and dual slater held */
   int                   boundedwslater;     /**< number of instances we could compute a bound for via the penalty approach where primal and dual slater held */
   int                   unsolvedwslater;    /**< number of instances that could not be solved where primal and dual slater held */
   int                   stablenoslater;     /**< number of instances solved with fastest setting where either primal or dual slater did not hold */
   int                   unstablenoslater;   /**< number of instances solved with stable settings where either primal or dual slater did not hold */
   int                   penaltynoslater;    /**< number of instances solved with penalty formulation where either primal or dual slater did not hold */
   int                   boundednoslater;    /**< number of instances we could compute a bound for via the penalty approach where either primal or dual slater did not hold */
   int                   unsolvednoslater;   /**< number of instances that could not be solved where either primal or dual slater did not hold */
   int                   nslaterholds;       /**< number of SDPs for which primal and dual slater condition held */
   int                   nnoslater;          /**< number of SDPs for which either primal or dual slater condition did not hold (including those where we could not check the other) */
   int                   nslatercheckfailed; /**< number of SDPs for which we failed to check the slater condition (this only includes SDPs where both checks failed or
                                              *   one checked returned slater holds and the other failed but not those where the first check already returned that it does not hold */
   int                   npslaterholds;      /**< number of SDPs for which primal slater condition held */
   int                   npnoslater;         /**< number of SDPs for which primal slater condition did not hold */
   int                   npslatercheckfailed;/**< number of SDPs for which we failed to check the dual slater condition */
   int                   ndslaterholds;      /**< number of SDPs for which dual slater condition held */
   int                   ndnoslater;         /**< number of SDPs for which dual slater condition did not hold */
   int                   ndslatercheckfailed;/**< number of SDPs for which we failed to check the dual slater condition */
   int                   nslaterinfeasible;  /**< number of SDPs for which we detected infeasibility during the Slater check */
   int                   stableinfeasible;   /**< number of instances solved with fastest settings where the dual slater check showed that the problem is infeasible */
   int                   unstableinfeasible; /**< number of instances solved with stable settings where the dual slater check showed that the problem is infeasible */
   int                   penaltyinfeasible;  /**< number of instances solved with penalty formulation where the dual slater check showed that the problem is infeasible */
   int                   boundedinfeasible;  /**< number of instances we could compute a bound for via the penalty approach where the dual slater check showed that the problem is infeasible */
   int                   unsolvedinfeasible; /**< number of instances that could not be solved where the dual slater check showed that the problem is infeasible */
};

/** inserts all the SDP data into the corresponding SDP Interface */
static
SCIP_RETCODE putSdpDataInInterface(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SdpVarmapper*         varmapper           /**< maps SCIP variables to their global SDP indices and vice versa */
   )
{
   SCIP_CONSHDLR* conshdlr;
   const char* conshdlrname;
   SCIP_CONS** conss;
   SCIP_VAR** blockvars;
   SCIP_VAR** vars;
   SCIP_Real*** val;
   SCIP_Real** constval;
   SCIP_Real* obj;
   SCIP_Real* lb;
   SCIP_Real* ub;
   SCIP_Real param;
   int*** row;
   int*** col;
   int** nblockvarnonz;
   int** constrow;
   int** constcol;
   int** sdpvar;
   int* sdpblocksizes;
   int* nblockvars;
   int* nconstblocknonz;
   int constnnonzcounter;
   int blocknnonz;
   int sdpconstnnonz;
   int sdpnnonz;
   int nsdpblocks;
   int constlength;
   int nvars;
   int nconss;
   int ind;
   int i;
   int j;

   SCIP_CALL( SCIPgetRealParam(scip, "relaxing/SDP/sdpsolvergaptol", &param) );

   SCIPdebugMessage("Putting SDP Data in general SDP interface!\n");

   assert( scip != NULL );
   assert( sdpi != NULL );
   assert( varmapper != NULL );

   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   /* prepare arrays of objective values and bounds */
   SCIP_CALL( SCIPallocBufferArray(scip, &obj, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lb, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ub, nvars) );

   for (i = 0; i < nvars; i++)
   {
      obj[i] = SCIPvarGetObj(vars[i]);
      lb[i] = SCIPvarGetLbLocal(vars[i]);
      ub[i] = SCIPvarGetUbLocal(vars[i]);
   }

   nconss = SCIPgetNConss(scip);
   conss = SCIPgetConss(scip);

   /* count the number of sdpblocks and compute the number of nonzeros */
   nsdpblocks = 0;
   sdpnnonz = 0;
   sdpconstnnonz = 0;

   for (i = 0; i < nconss; i++)
   {
      conshdlr = SCIPconsGetHdlr(conss[i]);
      assert( conshdlr != NULL );

      conshdlrname = SCIPconshdlrGetName(conshdlr);

#ifdef SCIP_EVEN_MORE_DEBUG
      SCIP_CALL( SCIPprintCons(scip, conss[i], NULL) );
      SCIPinfoMessage(scip, NULL, "\n");
#endif

      if ( strcmp(conshdlrname, "SDP") == 0 )
      {
         nsdpblocks++;

         SCIP_CALL( SCIPconsSdpGetNNonz(scip, conss[i], &blocknnonz, &constnnonzcounter) );
         sdpnnonz += blocknnonz;
         sdpconstnnonz += constnnonzcounter;
      }
   }

   /* create the sdp- and sdpconst-arrays */
   SCIP_CALL( SCIPallocBufferArray(scip, &sdpblocksizes, nsdpblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nblockvarnonz, nsdpblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nconstblocknonz, nsdpblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &col, nsdpblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &row, nsdpblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &val, nsdpblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &constcol, nsdpblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &constrow, nsdpblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &constval, nsdpblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nblockvars, nsdpblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sdpvar, nsdpblocks) );

   for (i = 0; i < nsdpblocks; i++)
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &(nblockvarnonz[i]), nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &col[i], nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &row[i], nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &val[i], nvars) );
   }

   /* get the SDP-data */
   ind = 0; /* index of the current sdp block in the complete sdp */
   SCIP_CALL( SCIPallocBufferArray(scip, &blockvars, nvars) );

   for (i = 0; i < nconss; i++)
   {
      conshdlr = SCIPconsGetHdlr(conss[i]);
      assert( conshdlr != NULL );

      conshdlrname = SCIPconshdlrGetName(conshdlr);

      if ( strcmp(conshdlrname, "SDP") == 0 )
      {
         assert( ind < nsdpblocks );

         /* allocate memory for the constant nonzeros */
         SCIP_CALL( SCIPconsSdpGetNNonz(scip, conss[i], NULL, &constlength) );
         nconstblocknonz[ind] = constlength;
         SCIP_CALL( SCIPallocBufferArray(scip, &(constcol[ind]), constlength) );
         SCIP_CALL( SCIPallocBufferArray(scip, &(constrow[ind]), constlength) );
         SCIP_CALL( SCIPallocBufferArray(scip, &(constval[ind]), constlength) );

         /* get the data */
         SCIP_CALL( SCIPconsSdpGetData(scip, conss[i], &nblockvars[ind], &blocknnonz, &sdpblocksizes[ind], &nvars, nblockvarnonz[ind], col[ind],
            row[ind], val[ind], blockvars, &nconstblocknonz[ind], constcol[ind], constrow[ind], constval[ind]) );

         /* nvars and nconstblocknonz[ind] would have been overwritten if the space in the given arrays hadn't been sufficient */
         assert( nvars == SCIPgetNVars(scip) );
         assert( nconstblocknonz[ind] <= constlength );

         SCIP_CALL( SCIPallocBufferArray(scip, &(sdpvar[ind]), nblockvars[ind]) );

         /* get global variable indices */
         for (j = 0; j < nblockvars[ind]; j++)
            sdpvar[ind][j] = SCIPsdpVarmapperGetSdpIndex(varmapper, blockvars[j]);

         ind++;
      }
   }

   /* free the memory that is no longer needed */
   SCIPfreeBufferArray(scip, &blockvars);

   /* load data into SDPI */
   SCIP_CALL( SCIPsdpiLoadSDP(sdpi, nvars,  obj, lb, ub, nsdpblocks, sdpblocksizes, nblockvars, sdpconstnnonz, nconstblocknonz, constrow,
                            constcol, constval, sdpnnonz, nblockvarnonz, sdpvar, row, col,  val, 0,
                            NULL, NULL, 0, NULL, NULL, NULL) ); /* insert the SDP part, add an empty LP part */

   /* free the remaining memory */
   for (i = 0; i < nsdpblocks; i++)
   {
      SCIPfreeBufferArrayNull(scip, &(sdpvar[i]));
      SCIPfreeBufferArrayNull(scip, &val[i]);
      SCIPfreeBufferArrayNull(scip, &row[i]);
      SCIPfreeBufferArrayNull(scip, &col[i]);
      SCIPfreeBufferArrayNull(scip, &(nblockvarnonz[i]));
      SCIPfreeBufferArrayNull(scip, &(constval[i]));
      SCIPfreeBufferArrayNull(scip, &(constrow[i]));
      SCIPfreeBufferArrayNull(scip, &(constcol[i]));
   }

   SCIPfreeBufferArrayNull(scip, &sdpvar);
   SCIPfreeBufferArrayNull(scip, &nblockvars);
   SCIPfreeBufferArrayNull(scip, &constval);
   SCIPfreeBufferArrayNull(scip, &constrow);
   SCIPfreeBufferArrayNull(scip, &constcol);
   SCIPfreeBufferArrayNull(scip, &val);
   SCIPfreeBufferArrayNull(scip, &row);
   SCIPfreeBufferArrayNull(scip, &col);
   SCIPfreeBufferArrayNull(scip, &nconstblocknonz);
   SCIPfreeBufferArrayNull(scip, &nblockvarnonz);
   SCIPfreeBufferArrayNull(scip, &sdpblocksizes);
   SCIPfreeBufferArray(scip, &ub);
   SCIPfreeBufferArray(scip, &lb);
   SCIPfreeBufferArray(scip, &obj);

   return SCIP_OKAY;
}

/** inserts all the LP data (including bounds and objective) into the corresponding SDP Interface */
static
SCIP_RETCODE putLpDataInInterface(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SdpVarmapper*         varmapper           /**< maps SCIP variables to their global SDP indices and vice versa */
   )
{
   SCIP_VAR** vars;
   SCIP_COL** rowcols;
   SCIP_ROW** rows;
   SCIP_Bool tightenvb;
   SCIP_Real* rowvals;
   SCIP_Real* lhs;
   SCIP_Real* rhs;
   SCIP_Real* obj;
   SCIP_Real* lb;
   SCIP_Real* ub;
   SCIP_Real* val;
   SCIP_Real sciplhs;
   SCIP_Real sciprhs;
   int* inds;
   int* objinds;
   int* rowind;
   int* colind;
   int nrowssdpi;
   int nrows;
   int rownnonz;
   int nvars;
   int nconss;
   int scipnnonz;
   int nnonz;
   int i;
   int j;

   assert( scip != NULL );
   assert( sdpi != NULL );
   assert( varmapper != NULL );

   nvars = SCIPgetNVars(scip);
   assert( nvars > 0 );

   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );
   SCIP_CALL( SCIPgetBoolParam(scip, "relaxing/SDP/tightenvb", &tightenvb) );

   SCIPdebugMessage("inserting %d LPRows into the interface.\n", nrows);

   /* compute the total number of LP nonzeroes in SCIP */
   scipnnonz = 0;
   for (i = 0; i < nrows; i++)
   {
      assert( rows[i] != NULL );
      scipnnonz += SCIProwGetNNonz(rows[i]);
   }

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &lhs, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rhs, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rowind, scipnnonz) );
   SCIP_CALL( SCIPallocBufferArray(scip, &colind, scipnnonz) );
   SCIP_CALL( SCIPallocBufferArray(scip, &val, scipnnonz) );

   /* insert the nonzeroes */
   nnonz = 0; /* this is recomputed for the sdpi, because of the possible duplication of non-zeroes for lhs and rhs */
   nconss = 0; /* this will be increased for each finite lhs and rhs */

   for (i = 0; i < nrows; i++)
   {
      SCIP_ROW* row;
      SCIP_Bool tightened = FALSE;
      SCIP_Real tightenedval = 0.0;
      SCIP_Bool swapped = FALSE;

      row = rows[i];
      assert( row != 0 );
      rownnonz = SCIProwGetNNonz(row);

      rowvals = SCIProwGetVals(row);
      rowcols = SCIProwGetCols(row);
      sciplhs = SCIProwGetLhs(row) - SCIProwGetConstant(row);
      sciprhs = SCIProwGetRhs(row) - SCIProwGetConstant(row);

      /* check whether we have a variable bound and can strenghten the big-M */
      if ( tightenvb && rownnonz == 2 && (SCIPisZero(scip, sciplhs) || SCIPisZero(scip, sciprhs) ) )
      {
         SCIP_VAR* var1;
         SCIP_VAR* var2;
         SCIP_Real val1;
         SCIP_Real val2;

         val1 = rowvals[0];
         val2 = rowvals[1];

         assert( rowcols[0] != NULL );
         assert( rowcols[1] != NULL );
         var1 = SCIPcolGetVar(rowcols[0]);
         var2 = SCIPcolGetVar(rowcols[1]);
         assert( var1 != NULL );
         assert( var2 != NULL );

         /* check that variables are not locally fixed */
         if ( ! SCIPisEQ(scip, SCIPvarGetLbLocal(var1), SCIPvarGetUbLocal(var1)) && ! SCIPisEQ(scip, SCIPvarGetLbLocal(var2), SCIPvarGetUbLocal(var2)) )
         {
            /* one coefficient must be 1 and the other negative */
            if ( (SCIPisEQ(scip, val1, 1.0) || SCIPisEQ(scip, val2, 1.0)) && ( SCIPisNegative(scip, val1) || SCIPisNegative(scip, val2) ) )
            {
               /* We want x - a z <= 0 or x - a z >= 0, where var1 = x and var2 = z; possibly swap variables otherwise */
               if ( ! SCIPisEQ(scip, val1, 1.0) || ! SCIPisNegative(scip, val2) )
               {
                  SCIPswapPointers((void**) &var1, (void**) &var2);

                  val2 = val1;
                  swapped = TRUE;
               }

               /* var2 needs to be binary */
               if ( SCIPvarIsBinary(var2) )
               {
                  if ( SCIPisZero(scip, sciprhs) )
                  {
                     if ( SCIPisLT(scip, SCIPvarGetUbLocal(var1), REALABS(val2)) )
                     {
                        SCIPdebugMessage("Big-M in %s changed from %f to %f\n", SCIProwGetName(row), REALABS(val2), SCIPvarGetUbLocal(var1));

                        tightened = TRUE;
                        tightenedval = -SCIPvarGetUbLocal(var1); /* negative sign because the coefficient needs to be negative */
                     }
                  }

                  if ( SCIPisZero(scip, sciplhs) )
                  {
                     if ( SCIPisGT(scip, SCIPvarGetLbLocal(var1), REALABS(val2)) )
                     {
                        SCIPdebugMessage("Big-M in %s changed from %f to %f\n", SCIProwGetName(row), REALABS(val2), SCIPvarGetLbLocal(var1));

                        tightened = TRUE;
                        tightenedval = -SCIPvarGetUbLocal(var1); /* negative sign because the coefficient needs to be negative */
                     }
                  }
               }
            }
         }
      }

      for (j = 0; j < rownnonz; j++)
      {
         /* if the Big-M was tightened, we use the new value (the position where this new value is used is dependant on wheter we needed to swap) */
         if ( tightened && ( (swapped && (j == 0)) || ((! swapped) && (j == 1)) ) ) /* use the tightened value */
         {
            if ( SCIPisFeasGT(scip, REALABS(tightenedval), 0.0) )
            {
               assert( SCIPcolGetVar(rowcols[j]) != 0 );
               colind[nnonz] = SCIPsdpVarmapperGetSdpIndex(varmapper, SCIPcolGetVar(rowcols[j]));
               rowind[nnonz] = nconss;
               val[nnonz] = tightenedval;
               nnonz++;
            }
         }
         else if ( SCIPisFeasGT(scip, REALABS(rowvals[j]), 0.0))
         {
            assert( SCIPcolGetVar(rowcols[j]) != 0 );
            colind[nnonz] = SCIPsdpVarmapperGetSdpIndex(varmapper, SCIPcolGetVar(rowcols[j]));
            rowind[nnonz] = nconss;
            val[nnonz] = rowvals[j];
            nnonz++;
         }
      }
      lhs[nconss] = sciplhs;
      rhs[nconss] = sciprhs;
      nconss++;
   }

   /* delete the old LP-block from the sdpi */
   SCIP_CALL( SCIPsdpiGetNLPRows(sdpi, &nrowssdpi) );
   if ( nrowssdpi > 0 )
   {
      SCIP_CALL( SCIPsdpiDelLPRows(sdpi, 0, nrowssdpi - 1) );
   }

   /* add the LP-block to the sdpi */
   SCIP_CALL( SCIPsdpiAddLPRows(sdpi, nconss, lhs, rhs, nnonz, (const int*)rowind, (const int*)colind, val) );

   /* free the remaining arrays */
   SCIPfreeBufferArray(scip, &val);
   SCIPfreeBufferArray(scip, &colind);
   SCIPfreeBufferArray(scip, &rowind);
   SCIPfreeBufferArray(scip, &rhs);
   SCIPfreeBufferArray(scip, &lhs);

   /* update bounds */

   /* get the variables */
   vars = SCIPgetVars(scip);
   assert( vars != NULL );

   /* prepare arrays of bounds */
   SCIP_CALL( SCIPallocBufferArray(scip, &lb, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ub, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &inds, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &obj, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &objinds, nvars) );

   /* get new bounds and objective coefficients */
   for (i = 0; i < nvars; i++)
   {
      assert( vars[i] != NULL );
      lb[i] = SCIPvarGetLbLocal(vars[i]);
      ub[i] = SCIPvarGetUbLocal(vars[i]);
      inds[i] = i; /* we want to change all bounds, so all indices are included in inds */
      obj[i] = SCIPvarGetObj(vars[i]);
      objinds[i] = i;
   }

   /* inform interface */
   SCIP_CALL( SCIPsdpiChgBounds(sdpi, nvars, inds, lb, ub) );
   SCIP_CALL( SCIPsdpiChgObj(sdpi, nvars, objinds, obj) );

   /* free the bounds-arrays */
   SCIPfreeBufferArray(scip, &objinds);
   SCIPfreeBufferArray(scip, &obj);
   SCIPfreeBufferArray(scip, &inds);
   SCIPfreeBufferArray(scip, &ub);
   SCIPfreeBufferArray(scip, &lb);

   return SCIP_OKAY;
}

/** calculate relaxation and process the relaxation results */
static
SCIP_RETCODE calcRelax(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RELAXDATA*       relaxdata,          /**< data of the relaxator */
   SCIP_RESULT*          result,             /**< pointer to store result of relaxation process */
   SCIP_Real*            lowerbound          /**< pointer to store lowerbound */
   )
{
   char saveconsname[SCIP_MAXSTRLEN];
   SCIP_SDPSOLVERSETTING startsetting;
   SCIP_SDPSOLVERSETTING usedsetting;
   SCIP_CONS* savedsetting;
   SCIP_CONS** conss;
   SCIP_VAR** vars;
   SCIP_SDPI* sdpi;
   SCIP_Bool rootnode;
   SCIP_Bool enforceslater;
   SCIP_Real timelimit;
   SCIP_Real objforscip;
   SCIP_Real* solforscip;
   SCIP_SDPSLATERSETTING slatersetting;
   SCIP_SDPSLATER primalslater;
   SCIP_SDPSLATER dualslater;
   int naddediters;
   int naddedsdpcalls;
   int nvars;
   int v;

   SCIPdebugMessage("calcRelax called\n");

   assert( scip != NULL );
   assert( relaxdata != NULL );
   assert( result != NULL );
   assert( lowerbound != NULL );

   nvars = SCIPgetNVars(scip);
   assert( nvars > 0 );
   vars = SCIPgetVars (scip);

   sdpi = relaxdata->sdpi;
   assert( sdpi != NULL );

   if ( relaxdata->objlimit )
   {
      /* set the objective limit */
      assert( SCIPgetUpperbound(scip) > -SCIPsdpiInfinity(sdpi) );
      SCIP_CALL( SCIPsdpiSetRealpar(sdpi, SCIP_SDPPAR_OBJLIMIT, SCIPgetUpperbound(scip)) );
   }
   /* if this is the root node and we cannot solve the problem, we want to check for the Slater condition independent from the SCIP parameter */
   rootnode = ! SCIPnodeGetParent(SCIPgetCurrentNode(scip));

   /* find settings to use for this relaxation */
   if ( rootnode || (SCIPnodeGetDepth(SCIPgetCurrentNode(scip)) == relaxdata->settingsresetofs) ||
      ( relaxdata->settingsresetfreq > 0 && ((SCIPnodeGetDepth(SCIPgetCurrentNode(scip)) - relaxdata->settingsresetofs) % relaxdata->settingsresetfreq == 0)) )
   {
      startsetting = SCIP_SDPSOLVERSETTING_UNSOLVED; /* in the root node we have no information, at each multiple of resetfreq we reset */
   }
   else
   {
      SCIP_CONSHDLR* conshdlr;
      int parentconsind;

      /* get constraint handler */
      conshdlr = SCIPfindConshdlr(scip, "Savedsdpsettings");
      if ( conshdlr == NULL )
      {
         SCIPerrorMessage("Savedsdpsettings constraint handler not found!\n");
         return SCIP_PLUGINNOTFOUND;
      }

      /* get startsettings of parent node, usually it will be the last active constraint of the corresponding constraint handler, so we iterate from
       * the end of the list until we find the correct one */
      conss = SCIPconshdlrGetConss(conshdlr);
      parentconsind = SCIPconshdlrGetNActiveConss(conshdlr) - 1;
      (void) SCIPsnprintf(saveconsname, SCIP_MAXSTRLEN, "savedsettings_node_%d", SCIPnodeGetNumber(SCIPnodeGetParent(SCIPgetCurrentNode(scip))));

      while ( parentconsind >= 0 && strcmp(saveconsname, SCIPconsGetName(conss[parentconsind])) )
         parentconsind--;
      if ( parentconsind >= 0 )
         startsetting = SCIPconsSavedsdpsettingsGetSettings(scip, conss[parentconsind]);
      else
      {
         SCIPdebugMessage("Startsetting from parent node not found, restarting with fastest settings!\n");
         startsetting = SCIP_SDPSOLVERSETTING_UNSOLVED;
      }

   }

   /* set time limit */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
   if ( ! SCIPisInfinity(scip, timelimit) )
   {
      timelimit -= SCIPgetSolvingTime(scip);
      if ( timelimit <= 0.0 )
      {
         *result = SCIP_DIDNOTRUN;
         return SCIP_OKAY;
      }
   }

   /* if no dual bound is known (we are in the root node and not only repropagating), we will have to abort, so we want
    * to check the Slater condition in this case */
   enforceslater = SCIPisInfinity(scip, -1 * SCIPnodeGetLowerbound(SCIPgetCurrentNode(scip)));

   /* solve the problem */
   SCIP_CALL( SCIPsdpiSolve(sdpi, NULL, startsetting, enforceslater, timelimit) );
   relaxdata->lastsdpnode = SCIPnodeGetNumber(SCIPgetCurrentNode(scip));

   /* update calls, iterations and stability numbers (only if the SDP-solver was actually called) */
   relaxdata->sdpinterfacecalls++;
   naddedsdpcalls = 0;
   SCIP_CALL( SCIPsdpiGetSdpCalls(relaxdata->sdpi, &naddedsdpcalls) );
   usedsetting = SCIP_SDPSOLVERSETTING_UNSOLVED;
   if ( naddedsdpcalls )
   {
      relaxdata->sdpcalls += naddedsdpcalls;
      naddediters = 0;
      SCIP_CALL( SCIPsdpiGetIterations(relaxdata->sdpi, &naddediters) );
      relaxdata->sdpiterations += naddediters;

      SCIP_CALL( SCIPsdpiSettingsUsed(relaxdata->sdpi, &usedsetting) );

      switch( usedsetting )/*lint --e{788}*/
      {
      case SCIP_SDPSOLVERSETTING_PENALTY:
         relaxdata->solvedpenalty++;
         break;
      case SCIP_SDPSOLVERSETTING_FAST:
         relaxdata->solvedfast++;
         break;
      case SCIP_SDPSOLVERSETTING_MEDIUM:
         relaxdata->solvedmedium++;
         break;
      case SCIP_SDPSOLVERSETTING_STABLE:
         relaxdata->solvedstable++;
         break;
      case SCIP_SDPSOLVERSETTING_UNSOLVED:
         relaxdata->unsolved++;
         break;
      default:
         break;
      }
      primalslater = SCIP_SDPSLATER_NOINFO;
      dualslater = SCIP_SDPSLATER_NOINFO;
      SCIP_CALL( SCIPsdpiSlater(relaxdata->sdpi, &primalslater, &dualslater) );
      switch( primalslater )/*lint --e{788}*/
      {
         case SCIP_SDPSLATER_NOINFO:
            relaxdata->npslatercheckfailed++;
            switch( dualslater )/*lint --e{788}*/
            {
               case SCIP_SDPSLATER_NOINFO:
                  relaxdata->ndslatercheckfailed++;
                  relaxdata->nslatercheckfailed++;
                  break;
               case SCIP_SDPSLATER_NOT:
                  relaxdata->ndnoslater++;
                  relaxdata->nnoslater++;
                  break;
               case SCIP_SDPSLATER_HOLDS:
                  relaxdata->ndslaterholds++;
                  relaxdata->nslatercheckfailed++;
                  break;
               case SCIP_SDPSLATER_INF:
                  relaxdata->nslaterinfeasible++;
                  break;
               default:
                  relaxdata->ndslatercheckfailed++;
                  relaxdata->nslatercheckfailed++;
                  break;
            }
            break;
         case SCIP_SDPSLATER_NOT:
            relaxdata->npnoslater++;
            switch( dualslater )/*lint --e{788}*/
            {
               case SCIP_SDPSLATER_NOINFO:
                  relaxdata->ndslatercheckfailed++;
                  relaxdata->nnoslater++;
                  break;
               case SCIP_SDPSLATER_NOT:
                  relaxdata->ndnoslater++;
                  relaxdata->nnoslater++;
                  break;
               case SCIP_SDPSLATER_HOLDS:
                  relaxdata->ndslaterholds++;
                  relaxdata->nnoslater++;
                  break;
               case SCIP_SDPSLATER_INF:
                  relaxdata->nslaterinfeasible++;
                  break;
               default:
                  relaxdata->ndslatercheckfailed++;
                  relaxdata->nnoslater++;
                  break;
            }
            break;
         case SCIP_SDPSLATER_HOLDS:
            relaxdata->npslaterholds++;
            switch( dualslater )/*lint --e{788}*/
            {
               case SCIP_SDPSLATER_NOINFO:
                  relaxdata->ndslatercheckfailed++;
                  relaxdata->nslatercheckfailed++;
                  break;
               case SCIP_SDPSLATER_NOT:
                  relaxdata->ndnoslater++;
                  relaxdata->nnoslater++;
                  break;
               case SCIP_SDPSLATER_HOLDS:
                  relaxdata->ndslaterholds++;
                  relaxdata->nslaterholds++;
                  break;
               case SCIP_SDPSLATER_INF:
                  relaxdata->nslaterinfeasible++;
                  break;
               default:
                  relaxdata->ndslatercheckfailed++;
                  relaxdata->nslatercheckfailed++;
                  break;
            }
            break;
            default:
               relaxdata->npslatercheckfailed++;
               relaxdata->ndslatercheckfailed++;
               relaxdata->nslatercheckfailed++;
               break;
      }
      slatersetting = SCIP_SDPSLATERSETTING_NOINFO;
      SCIP_CALL( SCIPsdpiSlaterSettings(relaxdata->sdpi, &slatersetting) );
      switch( slatersetting )/*lint --e{788}*/
      {
         case SCIP_SDPSLATERSETTING_STABLEWSLATER:
            relaxdata->stablewslater++;
            break;
         case SCIP_SDPSLATERSETTING_UNSTABLEWSLATER:
            relaxdata->unstablewslater++;
            break;
         case SCIP_SDPSLATERSETTING_PENALTYWSLATER:
            relaxdata->penaltywslater++;
            break;
         case SCIP_SDPSLATERSETTING_BOUNDEDWSLATER:
            relaxdata->boundedwslater++;
            break;
         case SCIP_SDPSLATERSETTING_UNSOLVEDWSLATER:
            relaxdata->unsolvedwslater++;
            break;
         case SCIP_SDPSLATERSETTING_STABLENOSLATER:
            relaxdata->stablenoslater++;
            break;
         case SCIP_SDPSLATERSETTING_UNSTABLENOSLATER:
            relaxdata->unstablenoslater++;
            break;
         case SCIP_SDPSLATERSETTING_PENALTYNOSLATER:
            relaxdata->penaltynoslater++;
            break;
         case SCIP_SDPSLATERSETTING_BOUNDEDNOSLATER:
            relaxdata->boundednoslater++;
            break;
         case SCIP_SDPSLATERSETTING_UNSOLVEDNOSLATER:
            relaxdata->unsolvednoslater++;
            break;
         case SCIP_SDPSLATERSETTING_STABLEINFEASIBLE:
            relaxdata->stableinfeasible++;
            break;
         case SCIP_SDPSLATERSETTING_UNSTABLEINFEASIBLE:
            relaxdata->unstableinfeasible++;
            break;
         case SCIP_SDPSLATERSETTING_PENALTYINFEASIBLE:
            relaxdata->penaltyinfeasible++;
            break;
         case SCIP_SDPSLATERSETTING_BOUNDEDINFEASIBLE:
            relaxdata->boundedinfeasible++;
            break;
         case SCIP_SDPSLATERSETTING_UNSOLVEDINFEASIBLE:
            relaxdata->unsolvedinfeasible++;
            break;
         default:
            break;
      }
   }

   /* remember settings */
   (void) SCIPsnprintf(saveconsname, SCIP_MAXSTRLEN, "savedsettings_node_%d", SCIPnodeGetNumber(SCIPgetCurrentNode(scip)));
   SCIP_CALL( createConsSavedsdpsettings(scip, &savedsetting, saveconsname, usedsetting) );
   SCIP_CALL( SCIPaddCons(scip, savedsetting) );
   SCIP_CALL( SCIPreleaseCons(scip, &savedsetting) );

   if ( ! SCIPsdpiWasSolved(sdpi) )
      relaxdata->feasible = FALSE;

   if ( SCIPinProbing(scip) )
      relaxdata->probingsolved = SCIPsdpiWasSolved(sdpi);
   else
      relaxdata->origsolved = SCIPsdpiSolvedOrig(sdpi);


   if ( SCIPsdpiIsAcceptable(sdpi) )
   {
#ifdef SCIP_MORE_DEBUG /* print the optimal solution */
      {
         int sollength;
         int i;
         SCIP_CALL( SCIPallocBufferArray(scip, &solforscip, nvars) );
         sollength = nvars;
         SCIP_CALL( SCIPsdpiGetSol(sdpi, &objforscip, solforscip, &sollength) ); /* get both the objective and the solution from the SDP solver */

         assert( sollength == nvars ); /* If this isn't true any longer, the getSol-call was unsuccessfull, because the given array wasn't long enough,
                                        * but this can't happen, because the array has enough space for all SDP variables. */

         if ( SCIPsdpiFeasibilityKnown(sdpi) )
         {
            SCIPdebugMessage("optimal solution: objective = %f, dual feasible: %u, primal feasible: %u.\n",
                  objforscip, SCIPsdpiIsDualFeasible(sdpi), SCIPsdpiIsPrimalFeasible(sdpi));
         }
         else
         {
            SCIPdebugMessage("The solver could not determine feasibility ! ");
         }

         /* output solution */
         for (i = 0; i < nvars; ++i)
         {
            SCIPdebugMessage("<%s> = %f\n", SCIPvarGetName(vars[i]), solforscip[i]);
         }
         SCIPfreeBufferArray(scip, &solforscip);
      }
#endif

      if ( SCIPsdpiIsDualInfeasible(sdpi) )
      {
         SCIPdebugMessage("Node cut off due to infeasibility.\n");
         relaxdata->feasible = FALSE;
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
      else if ( SCIPsdpiIsObjlimExc(sdpi) )
      {
         SCIPdebugMessage("Node cut off due to objective limit.\n");
         relaxdata->feasible = FALSE;
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
      else if ( SCIPsdpiIsDualUnbounded(sdpi) )
      {
         SCIPdebugMessage("Node unbounded.");
         relaxdata->feasible = TRUE;
         *result = SCIP_SUCCESS;
         *lowerbound = -SCIPinfinity(scip);
         return SCIP_OKAY;
      }
      else if ( SCIPsdpiIsPrimalFeasible(sdpi) && SCIPsdpiIsDualFeasible(sdpi) )
      {
         SCIP_SOL* scipsol;
         int slength;

         /* get solution w.r.t. SCIP variables */
         SCIP_CALL( SCIPallocBufferArray(scip, &solforscip, nvars) );
         slength = nvars;
         SCIP_CALL( SCIPsdpiGetSol(sdpi, &objforscip, solforscip, &slength) ); /* get both the objective and the solution from the SDP solver */

         assert( slength == nvars ); /* If this isn't true any longer, the getSol-Call was unsuccessfull, because the given array wasn't long enough,
                                      * but this can't happen, because the array has enough space for all sdp variables. */

         /* create SCIP solution */
         SCIP_CALL( SCIPcreateSol(scip, &scipsol, NULL) );
         for (v = 0; v < nvars; v++)
         {
            SCIP_CALL( SCIPsetSolVal(scip, scipsol, vars[v], solforscip[SCIPsdpVarmapperGetSdpIndex(relaxdata->varmapper, vars[v])]));
         }

         *lowerbound = objforscip;
         relaxdata->objval = objforscip;

         /* copy solution */
         SCIP_CALL( SCIPsetRelaxSolValsSol(scip, scipsol) );

         SCIP_CALL( SCIPmarkRelaxSolValid(scip) );

         relaxdata->feasible = TRUE;
         *result = SCIP_SUCCESS;

         SCIPfreeBufferArray(scip, &solforscip);
         SCIP_CALL( SCIPfreeSol(scip, &scipsol) );
      }
   }
   else
   {
      SCIP_Real objlb;

      if ( SCIPsdpiIsTimelimExc(relaxdata->sdpi) )
      {
         *result = SCIP_DIDNOTRUN;
         return SCIP_OKAY;
      }

      /* if we used the penalty approach, we might have calculated a good lower bound, even if we did not produce a feasible solution, otherwise we
       * keep the current bound, if the current bound is -infty, we abort */
      objlb = -SCIPinfinity(scip);
      SCIP_CALL( SCIPsdpiGetLowerObjbound(relaxdata->sdpi, &objlb) );
      if ( ! SCIPisInfinity(scip, objlb) )
      {
         *lowerbound = objlb;
         SCIPdebugMessage("The relaxation could not be solved, using best computed bound from penalty formulation.\n");
      }
      else if ( ! SCIPisInfinity(scip, -1 * SCIPnodeGetLowerbound(SCIPgetCurrentNode(scip))) )
      {
         *lowerbound = SCIPnodeGetLowerbound(SCIPgetCurrentNode(scip));
         SCIPdebugMessage("The relaxation could not be solved, keeping old bound.\n");
      }
      else
      {
         *result = SCIP_SUSPENDED;
         SCIPerrorMessage("The relaxation of the root node could not be solved, so there is no hope to solve this instance.\n");
         return SCIP_ERROR;
      }

      *result = SCIP_SUCCESS;
      return SCIP_OKAY;
   }

   return SCIP_OKAY;
}

/** checks whether all variables are fixed */
static
SCIP_Bool allVarsFixed(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_VAR** vars;
   int i;

   assert( scip != NULL );

   vars = SCIPgetVars(scip);

   /* try to find a variable that is not fixed */
   for (i = 0; i < SCIPgetNVars(scip); i++)
   {
      if ( SCIPisLT(scip, SCIPvarGetLbLocal(vars[i]), SCIPvarGetUbLocal(vars[i])) )
         return FALSE;
   }

   /* if no variable with lower bound strictly lower than upper bound has been found, all variables are fixed */
   return TRUE;
}

/** execution method of relaxator */
static
SCIP_DECL_RELAXEXEC(relaxExecSdp)
{
   SCIP_RELAXDATA* relaxdata;
   SCIP_VAR** vars;
   SCIP_Real* ubs;
   SCIP_Bool cutoff;
   SCIP_SOL* scipsol;
   int nconss;
   int nvars;
   int i;
#ifdef SCIP_EVEN_MORE_DEBUG
   SCIP_VAR** varsfordebug = SCIPgetVars(scip);
   const int nvarsfordebug = SCIPgetNVars(scip);
#endif

   SCIPdebugMessage("Calling relaxExecSdp.\n");

   relaxdata = SCIPrelaxGetData(relax);
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   /* don't run again if we already solved the current node (except during probing), and we solved the correct problem */
   if ( (relaxdata->lastsdpnode == SCIPnodeGetNumber(SCIPgetCurrentNode(scip)) && ( ! SCIPinProbing(scip) ) ) && relaxdata->origsolved && ! relaxdata->resolve )
   {
      SCIP_COL** cols;
      SCIP_Real objforscip;
      SCIP_Real* solforscip;
      int ncols;
      int slength;

      SCIPdebugMessage("Already solved SDP-relaxation for node %ld, returning with SCIP_SUCCESS so that no other relaxator is called.\n",
            SCIPrelaxGetData(relax)->lastsdpnode);

      if ( SCIPsdpiIsDualUnbounded(relaxdata->sdpi) )
      {
         relaxdata->feasible = TRUE;
         *result = SCIP_SUCCESS;
         *lowerbound = -SCIPinfinity(scip);
         return SCIP_OKAY;
      }

      /* get solution w.r.t. SCIP variables */
      SCIP_CALL( SCIPallocBufferArray(scip, &solforscip, nvars) );
      slength = nvars;
      SCIP_CALL( SCIPsdpiGetSol(relaxdata->sdpi, &objforscip, solforscip, &slength) ); /* get both the objective and the solution from the SDP solver */

      assert( slength == nvars ); /* If this isn't true any longer, the getSol-Call was unsuccessfull, because the given array wasn't long enough,
                                   * but this can't happen, because the array has enough space for all sdp variables. */

      /* create SCIP solution */
      SCIP_CALL( SCIPcreateSol(scip, &scipsol, NULL) );
      SCIP_CALL( SCIPsetSolVals(scip, scipsol, nvars, vars, solforscip) );

      /* Update the lower bound. Note that we cannot use the objective value given by the SDP-solver since this might
       * vary from the value SCIP computes internally because of rounding errors when extracting the solution from the
       * SDP-solver */
      *lowerbound = SCIPgetSolTransObj(scip, scipsol);

      /* copy solution */
      SCIP_CALL( SCIPgetLPColsData(scip, &cols, &ncols) );
      for (i = 0; i < ncols; i++)
      {
         SCIP_CALL( SCIPsetRelaxSolVal(scip, SCIPcolGetVar(cols[i]), SCIPgetSolVal(scip, scipsol, SCIPcolGetVar(cols[i]))) );
      }

      SCIP_CALL( SCIPmarkRelaxSolValid(scip) );
      *result = SCIP_SUCCESS;

      SCIPfreeBufferArray(scip, &solforscip);
      SCIP_CALL( SCIPfreeSol(scip, &scipsol) );

      *result = SCIP_SUCCESS;
      return SCIP_OKAY;
   }

   /* if we are solving a probing SDP, remember that we didn't solve the original problem */
   relaxdata->origsolved = FALSE;

   /* construct the lp and make sure, that everything is where it should be */
   SCIP_CALL( SCIPconstructLP(scip, &cutoff) );

   if ( cutoff )
   {
      relaxdata->feasible = FALSE;
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   /* very important to call flushLP */
   SCIP_CALL( SCIPflushLP(scip) );

   /* get varmapper */
   nconss = SCIPgetNConss(scip);

#ifdef SCIP_EVEN_MORE_DEBUG
   for (i = 0; i < nvarsfordebug; i++)
   {
      SCIPdebugMessage("variable %s: status = %u, integral = %u, bounds = [%f, %f] \n", SCIPvarGetName(varsfordebug[i]), SCIPvarGetStatus(varsfordebug[i]),
         SCIPvarIsIntegral(varsfordebug[i]), SCIPvarGetLbLocal(varsfordebug[i]), SCIPvarGetUbLocal(varsfordebug[i]));
   }
#endif

   if ( nconss == 0 )
   {
      /* if there are no constraints, there is nothing to do */
      relaxdata->feasible = TRUE;
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   if ( allVarsFixed(scip) )
   {
      SCIP_Bool feasible;

      /* if all variables, really all, are fixed, we give this fixed solution to SCIP */

      SCIP_CALL( SCIPallocBufferArray(scip, &ubs, nvars) );

      *lowerbound = 0.0;
      for (i = 0; i < nvars; i++)
      {
         ubs[i] = SCIPvarGetUbLocal(vars[i]);
         *lowerbound += SCIPvarGetObj(vars[i]) * ubs[i];
         assert( SCIPisFeasEQ(scip, SCIPvarGetUbLocal(vars[i]), SCIPvarGetLbLocal(vars[i])));
      }

      SCIPdebugMessage("EVERYTHING IS FIXED, objective value = %f\n", *lowerbound);

      SCIP_CALL( SCIPcreateSol(scip, &scipsol, NULL) );
      SCIP_CALL( SCIPsetSolVals(scip, scipsol, nvars, vars, ubs) );

      /* set the relaxation solution */
      for (i = 0; i < nvars; i++)
      {
         SCIP_CALL( SCIPsetRelaxSolVal(scip, vars[i], SCIPvarGetLbLocal(vars[i])) );
      }
      SCIP_CALL( SCIPmarkRelaxSolValid(scip) );

      /* check if the solution really is feasible */
      SCIP_CALL( SCIPcheckSol(scip, scipsol, FALSE, TRUE, TRUE, TRUE, TRUE, &feasible) );

      relaxdata->feasible = feasible;

      SCIP_CALL( SCIPfreeSol(scip, &scipsol) );

      SCIPfreeBufferArray(scip, &ubs);

      *result = SCIP_SUCCESS;
      return SCIP_OKAY;
   }

   /* update LP Data in Interface */
   SCIP_CALL( putLpDataInInterface(scip, relaxdata->sdpi, relaxdata->varmapper) );

   SCIP_CALL( calcRelax(scip, relaxdata, result, lowerbound));

   return SCIP_OKAY;
}


/** this method is called after presolving is finished, at this point the varmapper is prepared and the SDP Interface is initialized and gets
 *  the SDP information from the constraints */
static
SCIP_DECL_RELAXINITSOL(relaxInitSolSdp)
{
   SCIP_RELAXDATA* relaxdata;
   SCIP_RETCODE retcode;
   SCIP_VAR** vars;
   SCIP_Real gaptol;
   SCIP_Real feastol;
   SCIP_Real penaltyparam;
   SCIP_Real maxpenaltyparam;
   int npenaltyincr;
   SCIP_Bool sdpinfo;
   SCIP_Real givenpenaltyparam;
   int nthreads;
   int slatercheck;
   int nvars;

   assert( relax != NULL );

   relaxdata = SCIPrelaxGetData(relax);
   assert( relaxdata != NULL );

   relaxdata->objval = 0.0;
   relaxdata->origsolved = FALSE;
   relaxdata->probingsolved = FALSE;
   relaxdata->sdpcalls = 0;
   relaxdata->sdpinterfacecalls = 0;
   relaxdata->sdpiterations = 0;
   relaxdata->solvedfast = 0;
   relaxdata->solvedmedium = 0;
   relaxdata->solvedstable = 0;
   relaxdata->solvedpenalty = 0;
   relaxdata->stablewslater = 0;
   relaxdata->unstablewslater = 0;
   relaxdata->boundedwslater = 0;
   relaxdata->unsolvedwslater = 0;
   relaxdata->stablenoslater = 0;
   relaxdata->unsolvednoslater = 0;
   relaxdata->boundednoslater = 0;
   relaxdata->unsolvednoslater = 0;
   relaxdata->nslaterholds = 0;
   relaxdata->nnoslater = 0;
   relaxdata->nslatercheckfailed = 0;
   relaxdata->npslaterholds = 0;
   relaxdata->npnoslater = 0;
   relaxdata->npslatercheckfailed = 0;
   relaxdata->ndslaterholds = 0;
   relaxdata->ndnoslater = 0;
   relaxdata->ndslatercheckfailed = 0;
   relaxdata->nslaterinfeasible = 0;
   relaxdata->stableinfeasible = 0;
   relaxdata->unstableinfeasible = 0;
   relaxdata->penaltyinfeasible = 0;
   relaxdata->boundedinfeasible = 0;
   relaxdata->unsolvedinfeasible = 0;
   relaxdata->unsolved = 0;
   relaxdata->feasible = FALSE;

   nvars = SCIPgetNVars(scip);
   vars = SCIPgetVars(scip);

   /* all SCIPvars will be added to this list, and 3/4 seems like a good load factor (java uses this factor) */
   SCIP_CALL( SCIPsdpVarmapperCreate(scip, &(relaxdata->varmapper), (int) ceil(1.33 * nvars)) );
   SCIP_CALL( SCIPsdpVarmapperAddVars(scip, relaxdata->varmapper, nvars, vars) );

   if ( SCIPgetNVars(scip) > 0 )
   {
      SCIP_CALL( putSdpDataInInterface(scip, relaxdata->sdpi, relaxdata->varmapper) );
   }

   /* set the parameters of the SDP-Solver */
   SCIP_CALL( SCIPgetRealParam(scip, "relaxing/SDP/sdpsolvergaptol", &gaptol) );
   retcode = SCIPsdpiSetRealpar(relaxdata->sdpi, SCIP_SDPPAR_GAPTOL, gaptol);
   if ( retcode == SCIP_PARAMETERUNKNOWN )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL,
         "SDP Solver <%s>: gaptol setting not available -- SCIP parameter has no effect.\n",
         SCIPsdpiGetSolverName());
   }
   else
   {
      SCIP_CALL( retcode );
   }

   SCIP_CALL( SCIPgetRealParam(scip, "relaxing/SDP/sdpsolverfeastol", &feastol) );
   retcode = SCIPsdpiSetRealpar(relaxdata->sdpi, SCIP_SDPPAR_SDPSOLVERFEASTOL, feastol);
   if ( retcode == SCIP_PARAMETERUNKNOWN )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL,
         "SDP Solver <%s>: sdpsolverfeastol setting not available -- SCIP parameter has no effect.\n",
         SCIPsdpiGetSolverName());
   }
   else
   {
      SCIP_CALL( retcode );
   }

   retcode = SCIPsdpiSetRealpar(relaxdata->sdpi, SCIP_SDPPAR_EPSILON, SCIPepsilon(scip));
   if ( retcode == SCIP_PARAMETERUNKNOWN )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL,
         "SDP Solver <%s>: epsilon setting not available -- SCIP parameter has no effect.\n",
         SCIPsdpiGetSolverName());
   }
   else
   {
      SCIP_CALL( retcode );
   }

   retcode = SCIPsdpiSetRealpar(relaxdata->sdpi, SCIP_SDPPAR_FEASTOL, SCIPfeastol(scip));
   if ( retcode == SCIP_PARAMETERUNKNOWN )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL,
         "SDP Solver <%s>: feastol setting not available -- SCIP parameter has no effect.\n",
         SCIPsdpiGetSolverName());
   }
   else
   {
      SCIP_CALL( retcode );
   }

   /* set/compute the starting penalty parameter */
   SCIP_CALL( SCIPgetRealParam(scip, "relaxing/SDP/penaltyparam", &penaltyparam) );
   if ( SCIPisGE(scip, penaltyparam, 0.0) )
   {
      retcode = SCIPsdpiSetRealpar(relaxdata->sdpi, SCIP_SDPPAR_PENALTYPARAM, penaltyparam);
      givenpenaltyparam = penaltyparam;
      if ( retcode == SCIP_PARAMETERUNKNOWN )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL,
            "SDP Solver <%s>: penaltyparam setting not available -- SCIP parameter has no effect\n",
            SCIPsdpiGetSolverName());
      }
      else
      {
         SCIP_CALL( retcode );
      }
   }
   else
   {
      SCIP_Real maxcoeff;
      int v;

      /* compute the maximum coefficient in the objective */
      maxcoeff = 0.0;
      for (v = 0; v < nvars; v++)
      {
         if ( SCIPisGT(scip, REALABS(SCIPvarGetObj(vars[v])), maxcoeff) )
            maxcoeff = REALABS(SCIPvarGetObj(vars[v]));
      }

      SCIP_CALL( SCIPsdpiComputePenaltyparam(relaxdata->sdpi, maxcoeff, &givenpenaltyparam) );
   }

   /* set/compute the maximum penalty parameter */
   SCIP_CALL( SCIPgetRealParam(scip, "relaxing/SDP/maxpenaltyparam", &maxpenaltyparam) );
   if ( SCIPisGE(scip, maxpenaltyparam, 0.0) )
   {
      retcode = SCIPsdpiSetRealpar(relaxdata->sdpi, SCIP_SDPPAR_MAXPENALTYPARAM, maxpenaltyparam);

      if ( retcode == SCIP_PARAMETERUNKNOWN )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL,
            "SDP Solver <%s>: maxpenaltyparam setting not available -- SCIP parameter has no effect.\n",
            SCIPsdpiGetSolverName());
      }
      else
      {
         SCIP_CALL( retcode );
      }

      /* check if the starting value is not bigger than the maximum one, otherwise update it */
      if ( SCIPisLT(scip, givenpenaltyparam, maxpenaltyparam) )
      {
         SCIPdebugMessage("Penalty parameter %f overwritten by maxpenaltyparam %f! \n", givenpenaltyparam, maxpenaltyparam);
         SCIP_CALL( SCIPsdpiSetRealpar(relaxdata->sdpi, SCIP_SDPPAR_PENALTYPARAM, maxpenaltyparam) );
      }
   }
   else
   {
      SCIP_Real givenmaxpenaltyparam;

      SCIP_CALL( SCIPsdpiComputeMaxPenaltyparam(relaxdata->sdpi, givenpenaltyparam, &givenmaxpenaltyparam) );
   }

   /* set maximum number of penalty increasing rounds */
   SCIP_CALL( SCIPgetIntParam(scip, "relaxing/SDP/npenaltyincr", &npenaltyincr) );
   retcode = SCIPsdpiSetIntpar(relaxdata->sdpi, SCIP_SDPPAR_NPENALTYINCR, npenaltyincr);
   if ( retcode == SCIP_PARAMETERUNKNOWN )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL,
         "SDP Solver <%s>: npenaltyincr setting not available -- SCIP parameter has no effect.\n",
         SCIPsdpiGetSolverName());
   }
   else
   {
      SCIP_CALL( retcode );
   }


   /* set/compute lambda star if SDPA is used as the SDP-Solver */
   if ( strcmp(SCIPsdpiGetSolverName(), "SDPA") == 0.0 )
   {
      SCIP_Real lambdastar;

      SCIP_CALL( SCIPgetRealParam(scip, "relaxing/SDP/lambdastar", &lambdastar) );
      if ( SCIPisGE(scip, lambdastar, 0.0) )
      {
         retcode = SCIPsdpiSetRealpar(relaxdata->sdpi, SCIP_SDPPAR_LAMBDASTAR, lambdastar);
      }
      else
      {
         SCIP_Real guess;
         SCIP_Real maxguess;
         SCIP_CONS** conss;
         int nconss;
         int c;

         /* iterate over all SDP-constraints to compute the biggest guess for lambdastar  */
         conss = SCIPgetConss(scip);
         nconss = SCIPgetNConss(scip);
         maxguess = 0.0;

         for (c = 0; c < nconss; c++)
         {
            /* only check the SDP constraints */
            if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "SDP") == 0 )
            {
               SCIP_CALL( SCIPconsSdpGuessInitialPoint(scip, conss[c], &guess) );
               if ( (! SCIPisInfinity(scip, maxguess) ) && SCIPisGT(scip, guess, maxguess) )
                  maxguess = guess;
            }
         }

         SCIP_CALL( SCIPsdpiComputeLambdastar(relaxdata->sdpi, maxguess) );
      }
   }

   if ( retcode == SCIP_PARAMETERUNKNOWN )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL,
         "SDP Solver <%s>: lambdastar setting not available -- SCIP parameter has no effect.\n",
         SCIPsdpiGetSolverName());
   }
   else
   {
      SCIP_CALL( retcode );
   }

   SCIP_CALL( SCIPgetBoolParam(scip, "relaxing/SDP/sdpinfo", &sdpinfo) );
   retcode = SCIPsdpiSetIntpar(relaxdata->sdpi, SCIP_SDPPAR_SDPINFO, (int) sdpinfo);
   if ( retcode == SCIP_PARAMETERUNKNOWN )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL,
         "SDP Solver <%s>: sdpinfo setting not available -- SCIP parameter has no effect.\n",
         SCIPsdpiGetSolverName());
   }
   else
   {
      SCIP_CALL( retcode );
   }

   SCIP_CALL( SCIPgetIntParam(scip, "relaxing/SDP/sdpsolverthreads", &nthreads) );
   /* only try to set nthreads if the value differs from the default to prevent unnecessary warning messages for unknown parameter */
   if ( nthreads != DEFAULT_SDPSOLVERTHREADS )
   {
      retcode = SCIPsdpiSetIntpar(relaxdata->sdpi, SCIP_SDPPAR_NTHREADS, nthreads);
      if ( retcode == SCIP_PARAMETERUNKNOWN )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL,
            "SDP Solver <%s>: nthreads setting not available -- SCIP parameter has no effect.\n",
            SCIPsdpiGetSolverName());
      }
      else
      {
         SCIP_CALL( retcode );
      }
   }

   SCIP_CALL( SCIPgetIntParam(scip, "relaxing/SDP/slatercheck", &slatercheck) );
   retcode = SCIPsdpiSetIntpar(relaxdata->sdpi, SCIP_SDPPAR_SLATERCHECK, slatercheck);
   if ( retcode == SCIP_PARAMETERUNKNOWN )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL,
         "SDP Solver <%s>: slatercheck setting not available -- SCIP parameter has no effect.\n",
         SCIPsdpiGetSolverName());
   }
   else
   {
      SCIP_CALL( retcode );
   }

   /* initialize objective limit in case it was set in an earlier optimize call */
   SCIP_CALL( SCIPsdpiSetRealpar(relaxdata->sdpi, SCIP_SDPPAR_OBJLIMIT, SCIPsdpiInfinity(relaxdata->sdpi)) );

   return SCIP_OKAY;
}

/** copy method for SDP-relaxation handler (called when SCIP copies plugins) */
static
SCIP_DECL_RELAXCOPY(relaxCopySdp)
{
   assert( scip != NULL );
   assert( relax != NULL );
   assert(strcmp(SCIPrelaxGetName(relax), RELAX_NAME) == 0);

   SCIP_CALL( SCIPincludeRelaxSdp(scip) );

   return SCIP_OKAY;
}

/** reset the relaxator's data */
static
SCIP_DECL_RELAXEXIT(relaxExitSdp)
{
   SCIP_RELAXDATA* relaxdata;

   assert( scip != NULL );
   assert( relax != NULL );

   relaxdata = SCIPrelaxGetData(relax);
   assert( relaxdata != NULL );

   SCIPdebugMessage("Exiting Relaxation Handler.\n");

   if ( relaxdata->displaystat && SCIPgetSubscipDepth(scip) == 0 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "\nSDP iterations:\t\t\t\t%6d\n", relaxdata->sdpiterations);
      if ( relaxdata->sdpcalls )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Average SDP-iterations:\t\t\t%6.2f \n", (SCIP_Real) relaxdata->sdpiterations / (SCIP_Real) relaxdata->sdpcalls );
      }
      if ( relaxdata->sdpinterfacecalls )
      {
         if ( strcmp(SCIPsdpiGetSolverName(), "SDPA") == 0 )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Percentage 'fastest settings' solved:\t%6.2f \n", 100.0 * (SCIP_Real) relaxdata->solvedfast / (SCIP_Real) relaxdata->sdpinterfacecalls);
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Percentage 'medium settings' solved:\t%6.2f \n", 100.0 * (SCIP_Real) relaxdata->solvedmedium / (SCIP_Real) relaxdata->sdpinterfacecalls);
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Percentage 'stable settings' solved:\t%6.2f \n", 100.0 * (SCIP_Real) relaxdata->solvedstable / (SCIP_Real) relaxdata->sdpinterfacecalls);
         }
         else
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Percentage 'default formulation' solved:\t%6.2f \n", 100.0 * (SCIP_Real) relaxdata->solvedfast / (SCIP_Real) relaxdata->sdpinterfacecalls);
         }
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Percentage penalty formulation used:\t%6.2f \n", 100.0 * (SCIP_Real) relaxdata->solvedpenalty / (SCIP_Real) relaxdata->sdpinterfacecalls);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Percentage unsolved even with penalty:\t%6.2f \n", 100.0 * (SCIP_Real) relaxdata->unsolved / (SCIP_Real) relaxdata->sdpinterfacecalls);
      }
      if ( relaxdata->slatercheck )
      {
         if ( relaxdata->sdpinterfacecalls )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Percentage primal Slater condition held:\t%6.2f \n", 100.0 * (SCIP_Real) relaxdata->npslaterholds / (SCIP_Real) relaxdata->sdpinterfacecalls);
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Percentage primal Slater condition did not hold:\t%6.2f \n", 100.0 * (SCIP_Real) relaxdata->npnoslater / (SCIP_Real) relaxdata->sdpinterfacecalls);
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Percentage primal Slater check failed:\t%6.2f \n", 100.0 * (SCIP_Real) relaxdata->npslatercheckfailed / (SCIP_Real) relaxdata->sdpinterfacecalls);

            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Percentage dual Slater condition held:\t%6.2f \n", 100.0 * (SCIP_Real) relaxdata->ndslaterholds / (SCIP_Real) relaxdata->sdpinterfacecalls);
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Percentage dual Slater condition did not hold:\t%6.2f \n", 100.0 * (SCIP_Real) relaxdata->ndnoslater / (SCIP_Real) relaxdata->sdpinterfacecalls);
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Percentage dual Slater check failed:\t%6.2f \n", 100.0 * (SCIP_Real) relaxdata->ndslatercheckfailed / (SCIP_Real) relaxdata->sdpinterfacecalls);
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Percentage dual Slater check detected infeasibility:\t%6.2f \n", 100.0 * (SCIP_Real) relaxdata->nslaterinfeasible / (SCIP_Real) relaxdata->sdpinterfacecalls);
         }
         if ( relaxdata->nslaterholds )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Percentage 'fastest settings' with primal and dual slater holding:\t%6.2f \n",
                  100.0 * (SCIP_Real) relaxdata->stablewslater / (SCIP_Real) relaxdata->nslaterholds);
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Percentage 'stable settings' with primal and dual slater holding:\t%6.2f \n",
                  100.0 * (SCIP_Real) relaxdata->unstablewslater / (SCIP_Real) relaxdata->nslaterholds);
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Percentage 'penalty' with primal and dual slater holding:\t%6.2f \n",
                  100.0 * (SCIP_Real) relaxdata->penaltywslater / (SCIP_Real) relaxdata->nslaterholds);
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Percentage 'computed infeasible lower bound' with primal and dual slater holding:\t%6.2f \n",
                  100.0 * (SCIP_Real) relaxdata->boundedwslater / (SCIP_Real) relaxdata->nslaterholds);
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Percentage 'unsolved' with primal and dual slater holding:\t%6.2f \n",
                  100.0 * (SCIP_Real) relaxdata->unsolvedwslater / (SCIP_Real) relaxdata->nslaterholds);
         }
         if ( relaxdata->nnoslater )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Percentage 'fastest settings' with either primal or dual slater not holding:\t%6.2f \n",
                  100.0 * (SCIP_Real) relaxdata->stablenoslater / (SCIP_Real) relaxdata->nnoslater);
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Percentage 'stable settings' with either primal or dual slater not holding:\t%6.2f \n",
                  100.0 * (SCIP_Real) relaxdata->unstablenoslater / (SCIP_Real) relaxdata->nnoslater);
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Percentage 'penalty' with either primal or dual slater not holding:\t%6.2f \n",
                  100.0 * (SCIP_Real) relaxdata->penaltynoslater / (SCIP_Real) relaxdata->nnoslater);
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Percentage 'computed infeasible lower bound' with either primal or dual slater not holding:\t%6.2f \n",
                  100.0 * (SCIP_Real) relaxdata->boundednoslater / (SCIP_Real) relaxdata->nnoslater);
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Percentage 'unsolved' with either primal or dual slater not holding:\t%6.2f \n",
                  100.0 * (SCIP_Real) relaxdata->unsolvednoslater / (SCIP_Real) relaxdata->nnoslater);
         }
         if ( relaxdata->nslaterinfeasible )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Percentage 'fastest settings' with slater check showing infeasibility:\t%6.2f \n",
                  100.0 * (SCIP_Real) relaxdata->stableinfeasible / (SCIP_Real) relaxdata->nslaterinfeasible);
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Percentage 'stable settings' with slater check showing infeasibility:\t%6.2f \n",
                  100.0 * (SCIP_Real) relaxdata->unstableinfeasible / (SCIP_Real) relaxdata->nslaterinfeasible);
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Percentage 'penalty' with slater check showing infeasibility:\t%6.2f \n",
                  100.0 * (SCIP_Real) relaxdata->penaltyinfeasible / (SCIP_Real) relaxdata->nslaterinfeasible);
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Percentage 'computed infeasible lower bound' with slater check showing infeasibility:\t%6.2f \n",
                  100.0 * (SCIP_Real) relaxdata->boundedinfeasible / (SCIP_Real) relaxdata->nslaterinfeasible);
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Percentage 'unsolved' with slater check showing infeasibility:\t%6.2f \n",
                  100.0 * (SCIP_Real) relaxdata->unsolvedinfeasible / (SCIP_Real) relaxdata->nslaterinfeasible);
         }
#ifdef SLATERSOLVED_ABSOLUTE
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Number of nodes with primal and dual slater holding:\t%d \n", relaxdata->nslaterholds);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Number of nodes with 'fastest settings' and primal and dual slater holding:\t%d \n", relaxdata->stablewslater);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Number of nodes with 'stable settings' and primal and dual slater holding:\t%d \n", relaxdata->unstablewslater);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Number of nodes with 'penalty' and primal and dual slater holding:\t%d \n", relaxdata->penaltywslater);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Number of nodes with 'computed infeasible lower bound' and primal and dual slater holding:\t%d \n", relaxdata->boundedwslater);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Number of nodes with 'unsolved' and primal and dual slater holding:\t%d \n", relaxdata->unsolvedwslater);

         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Number of nodes with either primal or dual slater not holding:\t%d \n", relaxdata->nnoslater);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Number of nodes with 'fastest settings' and either primal or dual slater not holding:\t%d \n", relaxdata->stablenoslater);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Number of nodes with 'stable settings' and either primal or dual slater not holding:\t%d \n", relaxdata->unstablenoslater);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Number of nodes with 'penalty' and either primal or dual slater not holding:\t%d \n", relaxdata->penaltynoslater);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Number of nodes with 'computed infeasible lower bound' and either primal or dual slater not holding:\t%d \n", relaxdata->boundednoslater);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Number of nodes with 'unsolved' and either primal or dual slater not holding:\t%d \n", relaxdata->unsolvednoslater);

         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Number of infeasible nodes:\t%d \n", relaxdata->nslaterinfeasible);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Number of infeasible nodes with 'fastest settings':\t%d \n", relaxdata->stableinfeasible);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Number of infeasible nodes with 'stable settings':\t%d \n", relaxdata->unstableinfeasible);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Number of infeasible nodes with 'penalty':\t%d \n", relaxdata->penaltyinfeasible);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Number of infeasible nodes with 'computed infeasible lower bound':\t%d \n", relaxdata->boundedinfeasible);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Number of infeasible nodes with 'unsolved':\t%d \n", relaxdata->unsolvedinfeasible);
#endif
      }
   }

   if ( relaxdata->varmapper != NULL )
   {
      SCIP_CALL( SCIPsdpVarmapperFree(scip, &(relaxdata->varmapper)) );
   }

   relaxdata->objval = 0.0;
   relaxdata->origsolved = FALSE;
   relaxdata->probingsolved = FALSE;
   relaxdata->feasible = FALSE;
   relaxdata->sdpiterations = 0;
   relaxdata->sdpcalls = 0;
   relaxdata->sdpinterfacecalls = 0;
   relaxdata->lastsdpnode = 0;
   SCIP_CALL( SCIPsdpiClear(relaxdata->sdpi) );

   return SCIP_OKAY;
}

/** free the relaxator's data */
static
SCIP_DECL_RELAXFREE(relaxFreeSdp)
{/*lint --e{715}*/
   SCIP_RELAXDATA* relaxdata;

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   if ( relaxdata->sdpi != NULL )
   {
      SCIP_CALL( SCIPsdpiFree(&(relaxdata->sdpi)) );
   }

   SCIPfreeMemory(scip, &relaxdata);

   SCIPrelaxSetData(relax, NULL);

   return SCIP_OKAY;
}

/** creates the SDP-relaxator and includes it in SCIP */
SCIP_RETCODE SCIPincludeRelaxSdp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_RELAXDATA* relaxdata = NULL;
   SCIP_RELAX* relax;
   SCIP_SDPI* sdpi;

   assert( scip != NULL );

   /* create SDP-relaxator data */
   SCIP_CALL( SCIPallocMemory(scip, &relaxdata) );
   SCIP_CALL( SCIPsdpiCreate(&sdpi, SCIPgetMessagehdlr(scip), SCIPblkmem(scip), SCIPbuffer(scip)) );

   relaxdata->sdpi = sdpi;
   relaxdata->lastsdpnode = -1;

   /* include relaxator */
   SCIP_CALL( SCIPincludeRelaxBasic(scip, &relax, RELAX_NAME, RELAX_DESC, RELAX_PRIORITY, RELAX_FREQ, TRUE, relaxExecSdp, relaxdata) );
   assert( relax != NULL );

   /* include additional callbacks */
   SCIP_CALL( SCIPsetRelaxInitsol(scip, relax, relaxInitSolSdp) );
   SCIP_CALL( SCIPsetRelaxExit(scip, relax, relaxExitSdp) );
   SCIP_CALL( SCIPsetRelaxFree(scip, relax, relaxFreeSdp) );
   SCIP_CALL( SCIPsetRelaxCopy(scip, relax, relaxCopySdp) );

   /* add parameters for SDP-solver */
   SCIP_CALL( SCIPaddRealParam(scip, "relaxing/SDP/sdpsolvergaptol",
         "the stopping criterion for the duality gap the sdpsolver should use",
         &(relaxdata->sdpsolvergaptol), TRUE, DEFAULT_SDPSOLVERGAPTOL, 1e-20, 0.001, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "relaxing/SDP/sdpsolverfeastol",
         "the feasibility tolerance for the SDP solver",
         &(relaxdata->sdpsolverfeastol), TRUE, SCIPsdpiGetDefaultSdpiSolverFeastol(), 1e-17, 0.001, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "relaxing/SDP/penaltyparam",
         "the starting value of the penalty parameter Gamma used for the penalty formulation if the "
         "SDP solver didn't converge; set this to a negative value to compute the parameter depending on the given problem", &(relaxdata->penaltyparam),
         TRUE, DEFAULT_PENALTYPARAM, -1.0, 1e+20, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "relaxing/SDP/maxpenaltyparam",
         "the maximum value of the penalty parameter Gamma used for the penalty formulation if the "
         "SDP solver didn't converge; set this to a negative value to compute the parameter depending on the given problem", &(relaxdata->maxpenaltyparam),
         TRUE, DEFAULT_MAXPENALTYPARAM, -1.0, 1e+20, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "relaxing/SDP/npenaltyincr",
         "maximum number of times the penalty parameter will be increased if the penalty formulation failed", &(relaxdata->npenaltyincr), TRUE,
         SCIPsdpiGetDefaultSdpiSolverNpenaltyIncreases(), 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "relaxing/SDP/lambdastar",
         "the parameter lambda star used by SDPA to set the initial point;"
         "set this to a negative value to compute the parameter depending on the given problem", &(relaxdata->lambdastar),
         TRUE, DEFAULT_LAMBDASTAR, -1.0, 1e+20, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "relaxing/SDP/slatercheck",
         "Should the Slater condition for the primal and dual problem be checked ahead of solving each SDP? 0: no, 1: yes and output statistics, 2: yes and print warning for "
         "every problem not satisfying primal and dual Slater condition", &(relaxdata->slatercheck), TRUE, DEFAULT_SLATERCHECK, 0, 2, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "relaxing/SDP/sdpinfo",
         "Should the SDP solver output information to the screen?",
         &(relaxdata->sdpinfo), TRUE, DEFAULT_SDPINFO, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "relaxing/SDP/objlimit",
         "Should an objective limit be given to the SDP-Solver?",
         &(relaxdata->objlimit), TRUE, DEFAULT_OBJLIMIT, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "relaxing/SDP/resolve",
         "Should the relaxation be resolved after bound-tightenings were found during propagation (outside of probing)?",
         &(relaxdata->resolve), TRUE, DEFAULT_RESOLVE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "relaxing/SDP/tightenvb",
         "Should Big-Ms in varbound-like constraints be tightened before giving them to the SDP-solver ?",
         &(relaxdata->tightenvb), TRUE, DEFAULT_TIGHTENVB, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "relaxing/SDP/displaystatistics",
         "Should statistics about SDP iterations and solver settings/success be printed after quitting SCIP-SDP ?",
         &(relaxdata->displaystat), TRUE, DEFAULT_DISPLAYSTAT, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "relaxing/SDP/settingsresetfreq",
         "frequency for resetting parameters in SDP solver and trying again with fastest settings (-1: never, 0: only at depth settingsresetofs);"
         "currently only supported for SDPA",
         &(relaxdata->settingsresetfreq), TRUE, DEFAULT_SETTINGSRESETFREQ, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "relaxing/SDP/settingsresetofs",
         "frequency offset for resetting parameters in SDP solver and trying again with fastest settings; currently only supported for SDPA",
         &(relaxdata->settingsresetofs), TRUE, DEFAULT_SETTINGSRESETOFS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "relaxing/SDP/sdpsolverthreads",
         "number of threads the SDP solver should use (-1 = number of cores); currently only supported for MOSEK",
         &(relaxdata->sdpsolverthreads), TRUE, DEFAULT_SDPSOLVERTHREADS, -1, INT_MAX, NULL, NULL) );


   /* add description of SDP-solver */
   SCIP_CALL( SCIPincludeExternalCodeInformation(scip, SCIPsdpiGetSolverName(), SCIPsdpiGetSolverDesc()) );

   return SCIP_OKAY;
}


/* external functions */

/** gets the primal variables corresponding to the lower and upper variable-bounds in the dual problem
 *
 *  The last input should specify the length of the arrays. If this is less than the number of variables, the needed
 *  length will be returned and a debug message thrown.
 *
 *  @note If a variable is either fixed or unbounded in the dual problem, a zero will be returned for the non-existent
 *  primal variable.
 */
SCIP_RETCODE SCIPrelaxSdpGetPrimalBoundVars(
   SCIP_RELAX*           relax,              /**< SDP-relaxator to get information for */
   SCIP_Real*            lbvars,             /**< pointer to store the values of the variables corresponding to lower bounds in the dual problems */
   SCIP_Real*            ubvars,             /**< pointer to store the values of the variables corresponding to upper bounds in the dual problems */
   int*                  arraylength         /**< input: length of lbvars and ubvars <br>
                                              *   output: number of elements inserted into lbvars/ubvars (or needed length if it wasn't sufficient) */
   )
{
   SCIP_RELAXDATA* relaxdata;

   assert( relax != NULL );
   assert( lbvars != NULL );
   assert( ubvars != NULL );
   assert( arraylength != NULL );
   assert( *arraylength >= 0 );

   relaxdata = SCIPrelaxGetData(relax);
   assert( relaxdata != NULL );

   SCIP_CALL( SCIPsdpiGetPrimalBoundVars(relaxdata->sdpi, lbvars, ubvars, arraylength) );

   return SCIP_OKAY;
}

/** returns optimal objective value of the current SDP-relaxation if the last SDP-relaxation was successfully solved */
SCIP_RETCODE SCIPrelaxSdpRelaxVal(
   SCIP_RELAX*           relax,              /**< SDP-relaxator to get objective value for */
   SCIP_Bool*            success,            /**< pointer to store whether the last SDP-relaxation was solved successfully */
   SCIP_Real*            objval              /**< pointer to store the optimal objective value of the SDP-relaxation */
   )
{
   SCIP_RELAXDATA* relaxdata;

   assert( relax != NULL );
   assert( success != NULL );
   assert( objval != NULL );

   relaxdata = SCIPrelaxGetData(relax);
   assert( relaxdata != NULL );

   *success = relaxdata->origsolved;
   *objval = relaxdata->objval;

   return SCIP_OKAY;
}

/** returns values of all variables in the solution of the current SDP-relaxation if the last SDP-relaxation was successfully solved */
SCIP_RETCODE SCIPrelaxSdpGetRelaxSol(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_RELAX*           relax,              /**< SDP-relaxator to get solution for */
   SCIP_Bool*            success,            /**< pointer to store whether the last SDP-relaxation was solved successfully */
   SCIP_Real*            solarray,           /**< pointer to store the solution, this has to be at least length nvars */
   int*                  sollength           /**< length of the solarray */
   )
{
   SCIP_RELAXDATA* relaxdata;

   assert( relax != NULL );
   assert( success != NULL );
   assert( solarray != NULL );

   relaxdata = SCIPrelaxGetData(relax);
   assert( relaxdata != NULL );

   *success = relaxdata->origsolved;

   if ( *sollength >= SCIPgetNVars(scip) )
   {
      SCIP_CALL( SCIPsdpiGetSol(relaxdata->sdpi, NULL, solarray, sollength) );
   }
   else
   {
      SCIPdebugMessage("Called SCIPrelaxSdpGetRelaxSol with an array that wasn't big enough, needed length %d, given %d!\n", SCIPgetNVars(scip), *sollength);
      *sollength = SCIPgetNVars(scip);
   }

   return SCIP_OKAY;
}

/** get the number of the SCIP-node which the current SDP solution belongs to */
long int SCIPrelaxSdpGetSdpNode(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get solution for */
   )
{
   assert( relax != NULL );
   assert( SCIPrelaxGetData(relax) != NULL );

   return SCIPrelaxGetData(relax)->lastsdpnode;
}

/** Was the original problem solved for the last SDP-node (or a penalty or probing formulation) ? */
SCIP_Bool SCIPrelaxSdpSolvedOrig(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get solution for */
   )
{
   SCIP_RELAXDATA* relaxdata;

   assert( relax != NULL );

   relaxdata = SCIPrelaxGetData(relax);

   assert( relaxdata != NULL );
   assert( relaxdata->sdpi != NULL );

   return relaxdata->origsolved && SCIPsdpiSolvedOrig(relaxdata->sdpi);
}

/** Was the last probing SDP solved successfully ? */
SCIP_Bool SCIPrelaxSdpSolvedProbing(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get solution for */
   )
{
   SCIP_RELAXDATA* relaxdata;

   assert( relax != NULL );

   relaxdata = SCIPrelaxGetData(relax);

   assert( relaxdata != NULL );
   assert( relaxdata->sdpi != NULL );

   return relaxdata->probingsolved;
}

/** returns whether the last solved problem was feasible */
SCIP_Bool SCIPrelaxSdpIsFeasible(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get feasibility for */
   )
{
   assert( relax != NULL );
   assert( SCIPrelaxGetData(relax) != NULL );

   return ( SCIPrelaxGetData(relax)->feasible );
}

/** returns total number of SDP-iterations */
int SCIPrelaxSdpGetNIterations(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get the iterations for */
   )
{
   assert( relax != NULL );
   assert( SCIPrelaxGetData(relax) != NULL );

   return SCIPrelaxGetData(relax)->sdpiterations;
}

/** returns number of SDPs solved by SDP-solver (including multiple calls for penalty formulation etc.) */
int SCIPrelaxSdpGetNSdpCalls(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get the number of calls for */
   )
{
   assert( relax != NULL );
   assert( SCIPrelaxGetData(relax) != NULL );

   return ( SCIPrelaxGetData(relax)->sdpcalls );
}

/** returns number of solved SDP-relaxations */
int SCIPrelaxSdpGetNSdpInterfaceCalls(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get the number of calls for */
   )
{
   assert( relax != NULL );
   assert( SCIPrelaxGetData(relax) != NULL );

   return ( SCIPrelaxGetData(relax)->sdpinterfacecalls );
}

/** returns number of SDP-relaxations solved with fast settings */
int SCIPrelaxSdpGetNSdpFast(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get the number of calls for */
   )
{
   assert( relax != NULL );
   assert( SCIPrelaxGetData(relax) != NULL );

   return ( SCIPrelaxGetData(relax)->solvedfast );
}

/** returns number of SDP-relaxations solved with medium settings */
int SCIPrelaxSdpGetNSdpMedium(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get the number of calls for */
   )
{
   assert( relax != NULL );
   assert( SCIPrelaxGetData(relax) != NULL );

   return ( SCIPrelaxGetData(relax)->solvedmedium );
}

/** returns number of SDP-relaxations solved with stable settings */
int SCIPrelaxSdpGetNSdpStable(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get the number of calls for */
   )
{
   assert( relax != NULL );
   assert( SCIPrelaxGetData(relax) != NULL );

   return ( SCIPrelaxGetData(relax)->solvedstable );
}

/** returns number of SDP-relaxations solved with penalty formulation */
int SCIPrelaxSdpGetNSdpPenalty(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get the number of calls for */
   )
{
   assert( relax != NULL );
   assert( SCIPrelaxGetData(relax) != NULL );

   return ( SCIPrelaxGetData(relax)->solvedpenalty );
}

/** returns number of SDP-relaxations unsolved even when using a penalty formulation */
int SCIPrelaxSdpGetNSdpUnsolved(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get the number of calls for */
   )
{
   assert( relax != NULL );
   assert( SCIPrelaxGetData(relax) != NULL );

   return ( SCIPrelaxGetData(relax)->unsolved );
}
