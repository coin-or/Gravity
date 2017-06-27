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

/**@file   heur_sdpfracdiving.c
 * @brief  SDP diving heuristic that chooses fixings w.r.t. the fractionalities
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/* #define SCIP_DEBUG */
/* #define SCIP_MORE_DEBUG */ /* shows all diving decisions */

#include <assert.h>
#include <string.h>

#include "heur_sdpfracdiving.h"
#include "relax_sdp.h"

/* turn off lint warnings for whole file: */
/*lint --e{788,818}*/

#define HEUR_NAME             "sdpfracdiving"
#define HEUR_DESC             "SDP diving heuristic that chooses fixings w.r.t. the fractionalities"
#define HEUR_DISPCHAR         'f'
#define HEUR_PRIORITY         -1003000
#define HEUR_FREQ             -1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE
#define HEUR_USESSUBSCIP      FALSE  /* does the heuristic use a secondary SCIP instance? */


/*
 * Default parameter settings
 */

#define DEFAULT_MINRELDEPTH         0.0 /**< minimal relative depth to start diving */
#define DEFAULT_MAXRELDEPTH         1.0 /**< maximal relative depth to start diving */
#define DEFAULT_MAXDIVEUBQUOT       0.8 /**< maximal quotient (curlowerbound - lowerbound)/(cutoffbound - lowerbound)
                                         *   where diving is performed (0.0: no limit) */
#define DEFAULT_MAXDIVEAVGQUOT      0.0 /**< maximal quotient (curlowerbound - lowerbound)/(avglowerbound - lowerbound)
                                         *   where diving is performed (0.0: no limit) */
#define DEFAULT_MAXDIVEUBQUOTNOSOL  0.1 /**< maximal UBQUOT when no solution was found yet (0.0: no limit) */
#define DEFAULT_MAXDIVEAVGQUOTNOSOL 0.0 /**< maximal AVGQUOT when no solution was found yet (0.0: no limit) */
#define DEFAULT_BACKTRACK          TRUE /**< use one level of backtracking if infeasibility is encountered? */


/* locally defined heuristic data */
struct SCIP_HeurData
{
   SCIP_SOL*             sol;                /**< working solution */
   SCIP_Real             minreldepth;        /**< minimal relative depth to start diving */
   SCIP_Real             maxreldepth;        /**< maximal relative depth to start diving */
   SCIP_Real             maxdiveubquot;      /**< maximal quotient (curlowerbound - lowerbound)/(cutoffbound - lowerbound)
                                              *   where diving is performed (0.0: no limit) */
   SCIP_Real             maxdiveavgquot;     /**< maximal quotient (curlowerbound - lowerbound)/(avglowerbound - lowerbound)
                                              *   where diving is performed (0.0: no limit) */
   SCIP_Real             maxdiveubquotnosol; /**< maximal UBQUOT when no solution was found yet (0.0: no limit) */
   SCIP_Real             maxdiveavgquotnosol;/**< maximal AVGQUOT when no solution was found yet (0.0: no limit) */
   SCIP_Bool             backtrack;          /**< use one level of backtracking if infeasibility is encountered? */
   int                   nsuccess;           /**< number of runs that produced at least one feasible solution */
};

/*
 * Callback methods
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopySdpFracdiving)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurSdpFracdiving(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeSdpFracdiving)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);

   /* free heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}


/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitSdpFracdiving)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* create working solution */
   SCIP_CALL( SCIPcreateSol(scip, &heurdata->sol, heur) );

   /* initialize data */
   heurdata->nsuccess = 0;

   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitSdpFracdiving)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* free working solution */
   SCIP_CALL( SCIPfreeSol(scip, &heurdata->sol) );

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecSdpFracdiving)
{  /*lint --e{715}*/
   /* the current bugfix branch (3.2.1) does not have SCIPsolveProbingRelax() -> do nothing */
#if ( (SCIP_VERSION > 321 || SCIP_SUBVERSION > 0) )
   SCIP_HEURDATA* heurdata;
   SCIP_VAR** vars;
   SCIP_VAR* var;
   SCIP_VAR** sdpcands;
   SCIP_Bool cutoff;
   SCIP_Bool bestcandmayrounddown;
   SCIP_Bool bestcandmayroundup;
   SCIP_Bool bestcandroundup;
   SCIP_Bool mayrounddown;
   SCIP_Bool mayroundup;
   SCIP_Bool roundup;
   SCIP_Bool backtracked;
   SCIP_Bool backtrack;
   SCIP_Real* sdpcandssol;
   SCIP_Real* sdpcandsfrac;
   SCIP_Real searchubbound;
   SCIP_Real searchavgbound;
   SCIP_Real searchbound;
   SCIP_Real objval;
   SCIP_Real oldobjval;
   SCIP_Real obj;
   SCIP_Real objgain;
   SCIP_Real bestobjgain;
   SCIP_Real frac;
   SCIP_Real bestfrac;
   int nvars;
   int nsdpcands;
   int startnsdpcands;
   int depth;
   int maxdepth;
   int maxdivedepth;
   int divedepth;
   int bestcand;
   int c;
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

   /* don't dive two times at the same node */
   if ( SCIPgetLastDivenode(scip) == SCIPgetNNodes(scip) && SCIPgetDepth(scip) > 0 )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTRUN;

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   /* only try to dive, if we are in the correct part of the tree, given by minreldepth and maxreldepth */
   depth = SCIPgetDepth(scip);
   maxdepth = SCIPgetMaxDepth(scip);
   maxdepth = MAX(maxdepth, 30);
   if ( depth < heurdata->minreldepth*maxdepth || depth > heurdata->maxreldepth*maxdepth )
      return SCIP_OKAY;

   /* get fractional variables that should be integral */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   SCIP_CALL( SCIPallocBufferArray(scip, &sdpcands, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sdpcandssol, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sdpcandsfrac, nvars) );

   nsdpcands = 0;
   for (v = 0; v < nvars; ++v)
   {
      SCIP_Real val;

      val = SCIPgetRelaxSolVal(scip, vars[v]);
      frac = SCIPfeasFrac(scip, val);

      if ( SCIPvarIsIntegral(vars[v]) && ( ! SCIPisFeasZero(scip, frac) ) )
      {
         sdpcands[nsdpcands] = vars[v];
         sdpcandssol[nsdpcands] = val;
         sdpcandsfrac[nsdpcands] = frac;
         ++nsdpcands;
      }
   }

   /* don't try to dive, if there are no fractional variables */
   if ( nsdpcands == 0 )
   {
      SCIPfreeBufferArray(scip, &sdpcandsfrac);
      SCIPfreeBufferArray(scip, &sdpcandssol);
      SCIPfreeBufferArray(scip, &sdpcands);
      return SCIP_OKAY;
   }

   /* calculate the objective search bound */
   if ( SCIPgetNSolsFound(scip) == 0 )
   {
      if ( heurdata->maxdiveubquotnosol > 0.0 )
         searchubbound = SCIPgetLowerbound(scip) + heurdata->maxdiveubquotnosol * (SCIPgetCutoffbound(scip) - SCIPgetLowerbound(scip));
      else
         searchubbound = SCIPinfinity(scip);
      if ( heurdata->maxdiveavgquotnosol > 0.0 )
         searchavgbound = SCIPgetLowerbound(scip) + heurdata->maxdiveavgquotnosol * (SCIPgetAvgLowerbound(scip) - SCIPgetLowerbound(scip));
      else
         searchavgbound = SCIPinfinity(scip);
   }
   else
   {
      if ( heurdata->maxdiveubquot > 0.0 )
         searchubbound = SCIPgetLowerbound(scip) + heurdata->maxdiveubquot * (SCIPgetCutoffbound(scip) - SCIPgetLowerbound(scip));
      else
         searchubbound = SCIPinfinity(scip);
      if ( heurdata->maxdiveavgquot > 0.0 )
         searchavgbound = SCIPgetLowerbound(scip) + heurdata->maxdiveavgquot * (SCIPgetAvgLowerbound(scip) - SCIPgetLowerbound(scip));
      else
         searchavgbound = SCIPinfinity(scip);
   }
   searchbound = MIN(searchubbound, searchavgbound);
   if ( SCIPisObjIntegral(scip) )
      searchbound = SCIPceil(scip, searchbound);

   /* calculate the maximal diving depth: 10 * min{number of integer variables, max depth} */
   maxdivedepth = SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip);
   maxdivedepth = MIN(maxdivedepth, maxdepth);
   maxdivedepth *= 10;

   *result = SCIP_DIDNOTFIND;

   /* start diving */
   SCIP_CALL( SCIPstartProbing(scip) );

   /* enables collection of variable statistics during probing */
   SCIPenableVarHistory(scip);

   /* get SDP objective value*/
   objval = SCIPgetRelaxSolObj(scip);

   SCIPdebugMessage("(node %"SCIP_LONGINT_FORMAT") executing SDP fracdiving heuristic: depth=%d, %d fractionals, dualbound=%g, searchbound=%g\n",
      SCIPgetNNodes(scip), SCIPgetDepth(scip), nsdpcands, SCIPgetDualbound(scip), SCIPretransformObj(scip, searchbound));

   /* dive as long we are in the given objective, depth and iteration limits and fractional variables exist, but
    * - if possible, we dive at least with the depth 10
    * - if the number of fractional variables decreased at least with 1 variable per 2 dive depths, we continue diving
    */
   cutoff = FALSE;
   divedepth = 0;
   bestcandmayrounddown = FALSE;
   bestcandmayroundup = FALSE;
   startnsdpcands = nsdpcands;
   roundup = FALSE;

   while ( ! cutoff && nsdpcands > 0
      && ( divedepth < 10
         || nsdpcands <= startnsdpcands - divedepth/2
         || (divedepth < maxdivedepth && objval < searchbound))
      && ! SCIPisStopped(scip) )
   {
      SCIP_CALL( SCIPnewProbingNode(scip) );
      divedepth++;

      /* choose variable fixing:
       * - prefer variables that may not be rounded without destroying SDP feasibility:
       *   - of these variables, round least fractional variable in corresponding direction
       * - if all remaining fractional variables may be rounded without destroying SDP feasibility:
       *   - round variable with least increasing objective value
       */
      bestcand = -1;
      bestobjgain = SCIPinfinity(scip);
      bestfrac = SCIP_INVALID;
      bestcandmayrounddown = TRUE;
      bestcandmayroundup = TRUE;
      bestcandroundup = FALSE;
      for (c = 0; c < nsdpcands; ++c)
      {
         var = sdpcands[c];
         mayrounddown = SCIPvarMayRoundDown(var);
         mayroundup = SCIPvarMayRoundUp(var);
         frac = sdpcandsfrac[c];
         obj = SCIPvarGetObj(var);
         if ( mayrounddown || mayroundup )
         {
            /* the candidate may be rounded: choose this candidate only, if the best candidate may also be rounded */
            if( bestcandmayrounddown || bestcandmayroundup )
            {
               /* choose rounding direction:
                * - if variable may be rounded in both directions, round corresponding to the fractionality
                * - otherwise, round in the infeasible direction, because feasible direction is tried by rounding
                *   the current fractional solution
                */
               if ( mayrounddown && mayroundup )
                  roundup = (frac > 0.5);
               else
                  roundup = mayrounddown;

               if ( roundup )
               {
                  frac = 1.0 - frac;
                  objgain = frac*obj;
               }
               else
                  objgain = -frac*obj;

               /* penalize too small fractions */
               if ( frac < 0.01 )
                  objgain *= 1000.0;

               /* prefer decisions on binary variables */
               if ( !SCIPvarIsBinary(var) )
                  objgain *= 1000.0;

               /* check, if candidate is new best candidate */
               if ( SCIPisLT(scip, objgain, bestobjgain) || (SCIPisEQ(scip, objgain, bestobjgain) && frac < bestfrac) )
               {
                  bestcand = c;
                  bestobjgain = objgain;
                  bestfrac = frac;
                  bestcandmayrounddown = mayrounddown;
                  bestcandmayroundup = mayroundup;
                  bestcandroundup = roundup;
               }
            }
         }
         else
         {
            /* the candidate may not be rounded */
            if ( frac < 0.5 )
               roundup = FALSE;
            else
            {
               roundup = TRUE;
               frac = 1.0 - frac;
            }

            /* penalize too small fractions */
            if ( frac < 0.01 )
               frac += 10.0;

            /* prefer decisions on binary variables */
            if ( !SCIPvarIsBinary(var) )
               frac *= 1000.0;

            /* check, if candidate is new best candidate: prefer unroundable candidates in any case */
            if ( bestcandmayrounddown || bestcandmayroundup || frac < bestfrac )
            {
               bestcand = c;
               bestfrac = frac;
               bestcandmayrounddown = FALSE;
               bestcandmayroundup = FALSE;
               bestcandroundup = roundup;
            }
            assert( bestfrac < SCIP_INVALID );
         }
      }
      assert( bestcand != -1 );

      /* if all candidates are roundable, try to round the solution */
      if ( bestcandmayrounddown || bestcandmayroundup )
      {
         SCIP_Bool success;

         /* create solution from diving SDP and try to round it */
         SCIP_CALL( SCIPlinkRelaxSol(scip, heurdata->sol) );
         SCIP_CALL( SCIProundSol(scip, heurdata->sol, &success) );

         if ( success )
         {
            SCIPdebugMessage("SDP fracdiving found roundable primal solution: obj=%g\n", SCIPgetSolOrigObj(scip, heurdata->sol));

            /* try to add solution to SCIP */
            SCIP_CALL( SCIPtrySol(scip, heurdata->sol, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );

            /* check, if solution was feasible and good enough */
            if ( success )
            {
               SCIPdebugMessage(" -> solution was feasible and good enough\n");
               *result = SCIP_FOUNDSOL;
            }
         }
      }

      var = sdpcands[bestcand];

      backtracked = FALSE;
      do
      {
         backtrack = FALSE;
         /* if the variable is already fixed or if the solution value is outside the domain, numerical troubles may have
          * occured or variable was fixed by propagation while backtracking => Abort diving!
          */
         if ( SCIPvarGetLbLocal(var) >= SCIPvarGetUbLocal(var) - 0.5 )
         {
            SCIPdebugMessage("Selected variable <%s> already fixed to [%g,%g] (solval: %.9f), diving aborted \n",
               SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), sdpcandssol[bestcand]);
            cutoff = TRUE;
            break;
         }
         if ( SCIPisFeasLT(scip, sdpcandssol[bestcand], SCIPvarGetLbLocal(var)) || SCIPisFeasGT(scip, sdpcandssol[bestcand], SCIPvarGetUbLocal(var)) )
         {
            SCIPdebugMessage("selected variable's <%s> solution value is outside the domain [%g,%g] (solval: %.9f), diving aborted\n",
               SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), sdpcandssol[bestcand]);
#if 0
            assert( backtracked ); /* this may happen if we didn't resolve after propagation, in that case we will also abort (or resolve now and start again?) */
#endif
            break;
         }

         /* apply rounding of best candidate */
         if ( bestcandroundup != backtracked )
         {
            /* round variable up */
#ifdef SCIP_MORE_DEBUG
            SCIPdebugMessage("  dive %d/%d: var <%s>, round=%u/%u, sol=%g, oldbounds=[%g,%g], newbounds=[%g,%g]\n",
               divedepth, maxdivedepth, SCIPvarGetName(var), bestcandmayrounddown, bestcandmayroundup,
               sdpcandssol[bestcand], SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var),
               SCIPfeasCeil(scip, sdpcandssol[bestcand]), SCIPvarGetUbLocal(var));
#endif
            SCIP_CALL( SCIPchgVarLbProbing(scip, var, SCIPfeasCeil(scip, sdpcandssol[bestcand])) );
            roundup = TRUE;
         }
         else
         {
            /* round variable down */
#ifdef SCIP_MORE_DEBUG
            SCIPdebugMessage("  dive %d/%d: var <%s>, round=%u/%u, sol=%g, oldbounds=[%g,%g], newbounds=[%g,%g]\n",
               divedepth, maxdivedepth, SCIPvarGetName(var), bestcandmayrounddown, bestcandmayroundup,
               sdpcandssol[bestcand], SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var),
               SCIPvarGetLbLocal(var), SCIPfeasFloor(scip, sdpcandssol[bestcand]));
#endif
            SCIP_CALL( SCIPchgVarUbProbing(scip, var, SCIPfeasFloor(scip, sdpcandssol[bestcand])) );
            roundup = FALSE;
         }

         /* apply domain propagation */
         SCIP_CALL( SCIPpropagateProbing(scip, 0, &cutoff, NULL) );
         if ( !cutoff )
         {
            SCIP_RELAX* relaxsdp;

            /* resolve the diving SDP */
            SCIP_CALL( SCIPsolveProbingRelax(scip, &cutoff) );

            /* as cutoff doesn't work for relax sdp, we have to check ourselves, if we didn't manage to solve successfully, we abort diving */
            relaxsdp = SCIPfindRelax(scip, "SDP");

            /* if solving was unsuccessfull, abort */
            if (! SCIPrelaxSdpSolvedProbing(relaxsdp))
            {
               SCIPdebugMessage("SDP fracdiving heuristic aborted, as we could not solve one of the diving SDPs.\n");

               SCIPfreeBufferArray(scip, &sdpcandsfrac);
               SCIPfreeBufferArray(scip, &sdpcandssol);
               SCIPfreeBufferArray(scip, &sdpcands);

               *result = SCIP_DIDNOTRUN;
               SCIP_CALL( SCIPendProbing(scip) );

               return SCIP_OKAY;
            }

            cutoff = ! SCIPrelaxSdpIsFeasible(relaxsdp);
         }

         /* perform backtracking if a cutoff was detected */
         if ( cutoff && ! backtracked && heurdata->backtrack )
         {
#ifdef SCIP_MORE_DEBUG
            SCIPdebugMessage("  *** cutoff detected at level %d - backtracking\n", SCIPgetProbingDepth(scip));
#endif
            SCIP_CALL( SCIPbacktrackProbing(scip, SCIPgetProbingDepth(scip)-1) );
            SCIP_CALL( SCIPnewProbingNode(scip) );
            backtracked = TRUE;
            backtrack = TRUE;
         }
         else
            backtrack = FALSE;
      }
      while ( backtrack );

      if ( !cutoff )
      {
         /* get new objective value */
         oldobjval = objval;
         objval = SCIPgetRelaxSolObj(scip);

         /* update pseudo cost values */
         if ( SCIPisGT(scip, objval, oldobjval) )
         {
            if( roundup )
            {
               assert(bestcandroundup || backtracked);
               SCIP_CALL( SCIPupdateVarPseudocost(scip, sdpcands[bestcand], 1.0 - sdpcandsfrac[bestcand], objval - oldobjval, 1.0) );
            }
            else
            {
               assert(!bestcandroundup || backtracked);
               SCIP_CALL( SCIPupdateVarPseudocost(scip, sdpcands[bestcand], 0.0 - sdpcandsfrac[bestcand], objval - oldobjval, 1.0) );
            }
         }

         /* get new fractional variables */
         nsdpcands = 0;
         for (v = 0; v < nvars; ++v)
         {
            SCIP_Real val;

            val = SCIPgetRelaxSolVal(scip, vars[v]);
            frac = SCIPfeasFrac(scip, val);

            if ( SCIPvarIsIntegral(vars[v]) && ( ! SCIPisFeasZero(scip, frac) ) )
            {
               sdpcands[nsdpcands] = vars[v];
               sdpcandssol[nsdpcands] = val;
               sdpcandsfrac[nsdpcands] = frac;
               ++nsdpcands;
            }
         }
      }
#ifdef SCIP_MORE_DEBUG
      SCIPdebugMessage("   -> objval=%g/%g, nfrac=%d\n", objval, searchbound, nsdpcands);
#endif
   }

   /* check if a solution has been found */
   if ( nsdpcands == 0 && !cutoff )
   {
      SCIP_Bool success;

      /* create solution from diving SDP */
      SCIP_CALL( SCIPlinkRelaxSol(scip, heurdata->sol) );
      SCIPdebugMessage("SDP fracdiving found primal solution: obj=%g\n", SCIPgetSolOrigObj(scip, heurdata->sol));

      /* try to add solution to SCIP */
      SCIP_CALL( SCIPtrySol(scip, heurdata->sol, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );

      /* check, if solution was feasible and good enough */
      if ( success )
      {
         SCIPdebugMessage(" -> solution was feasible and good enough\n");
         *result = SCIP_FOUNDSOL;
      }
   }

   /* end diving */
   SCIP_CALL( SCIPendProbing(scip) );

   if ( *result == SCIP_FOUNDSOL )
      heurdata->nsuccess++;

   SCIPdebugMessage("SDP fracdiving heuristic finished\n");

   SCIPfreeBufferArray(scip, &sdpcandsfrac);
   SCIPfreeBufferArray(scip, &sdpcandssol);
   SCIPfreeBufferArray(scip, &sdpcands);

   return SCIP_OKAY;

#else
   *result = SCIP_DIDNOTRUN;

   return SCIP_OKAY;
#endif
}


/*
 * heuristic specific interface methods
 */

/** creates the SDP fracdiving heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurSdpFracdiving(
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
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecSdpFracdiving, heurdata) );

   assert( heur != NULL );

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopySdpFracdiving) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeSdpFracdiving) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitSdpFracdiving) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitSdpFracdiving) );

   /* fracdiving heuristic parameters */
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/sdpfracdiving/minreldepth",
         "minimal relative depth to start diving",
         &heurdata->minreldepth, TRUE, DEFAULT_MINRELDEPTH, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/sdpfracdiving/maxreldepth",
         "maximal relative depth to start diving",
         &heurdata->maxreldepth, TRUE, DEFAULT_MAXRELDEPTH, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/sdpfracdiving/maxdiveubquot",
         "maximal quotient (curlowerbound - lowerbound)/(cutoffbound - lowerbound) where diving is performed (0.0: no limit)",
         &heurdata->maxdiveubquot, TRUE, DEFAULT_MAXDIVEUBQUOT, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/sdpfracdiving/maxdiveavgquot",
         "maximal quotient (curlowerbound - lowerbound)/(avglowerbound - lowerbound) where diving is performed (0.0: no limit)",
         &heurdata->maxdiveavgquot, TRUE, DEFAULT_MAXDIVEAVGQUOT, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/sdpfracdiving/maxdiveubquotnosol",
         "maximal UBQUOT when no solution was found yet (0.0: no limit)",
         &heurdata->maxdiveubquotnosol, TRUE, DEFAULT_MAXDIVEUBQUOTNOSOL, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/sdpfracdiving/maxdiveavgquotnosol",
         "maximal AVGQUOT when no solution was found yet (0.0: no limit)",
         &heurdata->maxdiveavgquotnosol, TRUE, DEFAULT_MAXDIVEAVGQUOTNOSOL, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "heuristics/sdpfracdiving/backtrack",
         "use one level of backtracking if infeasibility is encountered?",
         &heurdata->backtrack, FALSE, DEFAULT_BACKTRACK, NULL, NULL) );

   return SCIP_OKAY;
}

