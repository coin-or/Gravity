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

/**@file   heur_sdprand.c
 * @brief  randomized rounding heuristic for SDPs
 * @author Marc Pfetsch
 * @author Tristan Gally
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/*#define SCIP_DEBUG*/
#include <assert.h>
#include <string.h>

#include "heur_sdprand.h"
#include "relax_sdp.h"

/* turn off lint warnings for whole file: */
/*lint --e{788,818}*/

#define HEUR_NAME             "sdprand"
#define HEUR_DESC             "randomized rounding heuristic for SDPs"
#define HEUR_DISPCHAR         '~'
#define HEUR_PRIORITY         -1001000
#define HEUR_FREQ             1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE
#define HEUR_USESSUBSCIP      FALSE  /* does the heuristic use a secondary SCIP instance? */


/*
 * Default parameter settings
 */

#define DEFAULT_NROUNDS                 2    /**< default number of rounding rounds */
#define DEFAULT_GENERALINTS             FALSE/**< Should randomized rounding also be applied if there are general integer variables and not only binary variables ? */
#define DEFAULT_RANDSEED                211  /**< default random seed */

/* locally defined heuristic data */
struct SCIP_HeurData
{
   SCIP_SOL*             sol;                /**< working solution */
   int                   nrounds;            /**< number of rounding rounds */
   SCIP_RANDNUMGEN*      randnumgen;         /**< random number generator */
   SCIP_Bool             generalints;        /**< Should randomized rounding also be applied if there are general integer variables and not only binary variables ? */
};


/*
 * Callback methods
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopySdprand)
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( heur != NULL );
   assert( strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0 );

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurSdpRand(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeSdprand)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert( heur != NULL );
   assert( strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0 );
   assert( scip != NULL );

   /* free heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}

/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitSdprand)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert( heur != NULL );
   assert( strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0 );

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   /* create working solution and random number generator */
   SCIP_CALL( SCIPcreateSol(scip, &heurdata->sol, heur) );
   SCIP_CALL( SCIPrandomCreate(&(heurdata->randnumgen), SCIPblkmem(scip), SCIPinitializeRandomSeed(scip, DEFAULT_RANDSEED)) );

   return SCIP_OKAY;
}

/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitSdprand)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert( heur != NULL );
   assert( strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0 );

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   /* free working solution and random number generator*/
   SCIP_CALL( SCIPfreeSol(scip, &heurdata->sol) );
   SCIPrandomFree(&(heurdata->randnumgen));

   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecSdprand)
{  /*lint --e{715}*/
   /* the current bugfix branch (3.2.1) does not have SCIPsolveProbingRelax() -> do nothing */
#if ( (SCIP_VERSION > 321 || SCIP_SUBVERSION > 0) )
   SCIP_HEURDATA* heurdata;
   SCIP_RELAX* relaxsdp;
   SCIP_Real* sdpcandssol;
   SCIP_VAR** sdpcands;
   SCIP_VAR** vars;
   int nsdpcands = 0;
   int ncontvars;
   int nvars;
   int iter;
   int v;

   assert( heur != NULL );
   assert( strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0 );
   assert( scip != NULL );
   assert( result != NULL );

   *result = SCIP_DELAYED;

   /* do not run if relaxation solution is not available */
   if ( ! SCIPisRelaxSolValid(scip) )
      return SCIP_OKAY;

   /* do not call heuristic if node was already detected to be infeasible */
   if ( nodeinfeasible )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTRUN;

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   /* only run if there are no general integer variables or the corresponding parameter is set */
   if ( (! heurdata->generalints) && SCIPgetNIntVars(scip) > 0 )
      return SCIP_OKAY;

   /* get relaxator - exit if not found (use LP randomized rounding) */
   relaxsdp = SCIPfindRelax(scip, "SDP");
   if ( relaxsdp == NULL )
      return SCIP_OKAY;

   /* get number of continuous variables */
   ncontvars = SCIPgetNContVars(scip) +  SCIPgetNImplVars(scip);

   /* save old SDP solution */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   SCIP_CALL( SCIPallocBufferArray(scip, &sdpcands, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sdpcandssol, nvars) );
   if ( ncontvars == 0)
   {
      /* in this case we do not need to solve an SDP again, so we only need to save the fractional values */
      for (v = 0; v < nvars; ++v)
      {
         SCIP_Real val;

         val = SCIPgetRelaxSolVal(scip, vars[v]);
         if ( SCIPvarIsIntegral(vars[v]) && ! SCIPisFeasIntegral(scip, val) )
         {
            sdpcands[nsdpcands] = vars[v];
            sdpcandssol[nsdpcands] = val;
            ++nsdpcands;
         }
      }
   }
   else
   {
      /* if there are continuous variables, we need to save all values to later be able to fix all current integer values */
      for (v = 0; v < nvars; ++v)
      {
         SCIP_Real val;

         val = SCIPgetRelaxSolVal(scip, vars[v]);
         sdpcands[v] = vars[v];
         sdpcandssol[v] = val;
         if ( SCIPvarIsIntegral(vars[v]) && ! SCIPisFeasIntegral(scip, val) )
         {
            ++nsdpcands;
         }
      }
   }

   /* do not proceed, if there are no fractional variables */
   if ( nsdpcands == 0 )
   {
      SCIPfreeBufferArray(scip, &sdpcandssol);
      SCIPfreeBufferArray(scip, &sdpcands);
      return SCIP_OKAY;
   }

   *result = SCIP_DIDNOTFIND;

   SCIPdebugMessage("node %"SCIP_LONGINT_FORMAT") executing SDP randomized rounding heuristic: depth=%d, %d fractionals.\n",
      SCIPgetNNodes(scip), SCIPgetDepth(scip), nsdpcands);

   /* perform rounding rounds */
   for (iter = 0; iter < heurdata->nrounds; ++iter)
   {
      SCIP_Bool success;
      SCIP_Bool cutoff;
      SCIP_VAR* var;
      SCIP_Real r;
      int cnt = 0;

      /* if there are no continuous variables, we can simply try the solution */
      if ( ncontvars == 0 )
      {
         SCIP_CALL( SCIPlinkRelaxSol(scip, heurdata->sol) );

         /* set values */
         for (v = 0; v < nsdpcands; ++v)
         {
            var = sdpcands[v];
            assert( SCIPvarIsIntegral(var) );

            /* if the variable is not fixed and its value is fractional */
            if ( SCIPvarGetLbLocal(var) < 0.5 && SCIPvarGetUbLocal(var) > 0.5 && ! SCIPisFeasIntegral(scip, sdpcandssol[v]) )
            {
               r = SCIPrandomGetReal(heurdata->randnumgen, 0.0, 1.0);

               /* depending on random value, set variable to 0 or 1 */
               if ( SCIPfeasFrac(scip, sdpcandssol[v]) <= r )
               {
                  SCIP_CALL( SCIPsetSolVal(scip, heurdata->sol, var, SCIPfeasFloor(scip, sdpcandssol[v])) );
               }
               else
               {
                  SCIP_CALL( SCIPsetSolVal(scip, heurdata->sol, var, SCIPfeasCeil(scip, sdpcandssol[v])) );
               }
            }
         }

         /* try to add solution to SCIP - do not need to check integrality here */
         SCIP_CALL( SCIPtrySol(scip, heurdata->sol, FALSE, FALSE, FALSE, FALSE, TRUE, &success) );

         if ( success )
            SCIPdebugMessage("Iteration %d: found solution for full binary instance.\n", iter);
         else
            SCIPdebugMessage("Iteration %d: solution not feasible for full binary instance.\n", iter);
      }
      else
      {
         /* if there are continuous variables, we need to solve a final SDP */
         SCIP_CALL( SCIPstartProbing(scip) );

         for (v = 0; v < nvars; ++v)
         {
            SCIP_Real val;

            var = vars[v];
            val = sdpcandssol[v];

            /* if the variable is not fixed and its value is fractional */
            if ( SCIPvarGetLbLocal(var) < 0.5 && SCIPvarGetUbLocal(var) > 0.5 && SCIPvarIsIntegral(var) && ! SCIPisFeasIntegral(scip, val) )
            {
               r = SCIPrandomGetReal(heurdata->randnumgen, 0.0, 1.0);

               /* depending on random value, fix variable to 0 or 1 */
               if ( val <= r )
               {
                  SCIP_CALL( SCIPchgVarUbProbing(scip, var, 0.0) );
               }
               else
               {
                  SCIP_CALL( SCIPchgVarLbProbing(scip, var, 1.0) );
               }
               ++cnt;
            }
            else if ( SCIPvarIsIntegral(var) && SCIPisFeasIntegral(scip, val) && SCIPisFeasGT(scip, val, SCIPvarGetLbLocal(var)) )
            {
               /* if an integral variable already attained an integer value, we fix it */
               SCIP_CALL( SCIPchgVarLbProbing(scip, var, val) );
               ++cnt;
            }
            else if ( SCIPvarIsIntegral(var) && SCIPisFeasIntegral(scip, val) && SCIPisFeasLT(scip, val, SCIPvarGetUbLocal(var)) )
            {
               /* if an integral variable already attained an integer value, we fix it */
               SCIP_CALL( SCIPchgVarUbProbing(scip, var, val) );
               ++cnt;
            }
         }

         /* if no value was changed */
         if ( cnt == 0 )
         {
            /* We can exit since there will be no chance to be successful in future runs. */
            SCIPdebugMessage("Iteration %d: All variables were fixed or their values were integral -> exit.\n", iter);
            break;
         }

         /* apply domain propagation (use parameter settings for maximal number of rounds) */
         SCIP_CALL( SCIPpropagateProbing(scip, 0, &cutoff, NULL) );
         if ( ! cutoff )
         {
            /* solve SDP */
            SCIP_CALL( SCIPsolveProbingRelax(scip, &cutoff) );

            /* if solving was successfull */
            if (SCIPrelaxSdpSolvedProbing(relaxsdp) && SCIPrelaxSdpIsFeasible(relaxsdp) )
            {
               /* check solution */
               SCIP_CALL( SCIPlinkRelaxSol(scip, heurdata->sol) );

               /* try to add solution to SCIP: check all constraints, including integrality */
               SCIP_CALL( SCIPtrySol(scip, heurdata->sol, FALSE, TRUE, TRUE, TRUE, TRUE, &success) );
               /* check, if solution was feasible and good enough */
               if ( success )
               {
                  SCIPdebugMessage("Iteration %d: solution was feasible and good enough.\n", iter);
                  *result = SCIP_FOUNDSOL;
               }
               else
                  SCIPdebugMessage("Iteration %d: solution was not feasible.\n", iter);
            }
         }

         /* free local problem */
         SCIP_CALL( SCIPendProbing(scip) );

         if ( SCIPisStopped(scip) )
            break;
      }
   }

   SCIPfreeBufferArray(scip, &sdpcandssol);
   SCIPfreeBufferArray(scip, &sdpcands);

   SCIPdebugMessage("finished randomized rounding heuristic.\n");

   return SCIP_OKAY;

#else
   *result = SCIP_DIDNOTRUN;

   return SCIP_OKAY;
#endif
}


/*
 * heuristic specific interface methods
 */

/** creates the randomized rounding heuristic for SDPs and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurSdpRand(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create Fracdiving primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecSdprand, heurdata) );

   assert( heur != NULL );

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopySdprand) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeSdprand) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitSdprand) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitSdprand) );

   /* fracdiving heuristic parameters */
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/sdprand/nrounds",
         "number of rounding rounds",
         &heurdata->nrounds, FALSE, DEFAULT_NROUNDS, 0, 10000, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "heuristics/sdprand/generalints",
         "Should randomized rounding also be applied if there are general integer variables and not only binary variables ?",
         &heurdata->generalints, FALSE, DEFAULT_GENERALINTS, NULL, NULL) );

   return SCIP_OKAY;
}
