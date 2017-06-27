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

/**@file   branch_sdpinfobjective.c
 * @brief  combined infeasibility and absolute objective branching rule for SCIP-SDP
 * @author Tristan Gally
 *
 * Branch on integral variable with highest product of fractionality/integral-infeasibility and absolute objective value in the SDP.
 *
 * Will do nothing for continuous variables, since these are what the external callbacks of the SCIP branching rules are for.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/* #define SCIP_DEBUG */

#include <assert.h>
#include <string.h>

#include "branch_sdpinfobjective.h"

/* turn off lint warnings for whole file: */
/*lint --e{788,818}*/

#define BRANCHRULE_NAME            "sdpinfobjective"
#define BRANCHRULE_DESC            "branch on variable with highest product of fractionality/integral-infeasibility and absolute objective of the SDP"
#define BRANCHRULE_PRIORITY        2000000
#define BRANCHRULE_MAXDEPTH        -1
#define BRANCHRULE_MAXBOUNDDIST    1.0
#define DEFAULT_COUPLEDVARS        FALSE     /**< if all branching candidates have objective zero, should we use the sum of the absolute objectives of all
                                               *  continuous variables coupled with the candidate through constraints */
#define DEFAULT_SINGLECOUPLEDVARS  FALSE     /**< if all branching candidates have objective zero, should we use the sum of the absolute objectives of all
                                               *  continuous variables coupled with the candidate through constraints in which no other candidate appears */

/** branching rule data */
struct SCIP_BranchruleData
{
   SCIP_Bool             coupledvars;        /**< If all branching candidates have objective zero, should we use the sum of the absolute objectives of all
                                               *  continuous variables coupled with the candidate through linear constraints ? */
   SCIP_Bool             singlecoupledvars;  /**< If all branching candidates have objective zero, should we use the sum of the absolute objectives of all
                                               *  continuous variables coupled with the candidate through linear constraints in which no other candidate appears ? */
};


/*
 * Data structures
 */

/*
 * Local methods
 */

/* put your local methods here, and declare them static */


/*
 * Callback methods of branching rule
 */

/** copy method for branchrule plugins (called when SCIP copies plugins) */
static
SCIP_DECL_BRANCHCOPY(branchCopySdpinfobjective)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);

   /* call inclusion method of branchrule */
   SCIP_CALL( SCIPincludeBranchruleSdpinfobjective(scip) );

   return SCIP_OKAY;
}

/** branching execution method for external candidates */
static
SCIP_DECL_BRANCHEXECEXT(branchExecextSdpinfobjective)
{/*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_VAR** cands = NULL;
   SCIP_Real* candssol;              /* solution values of all candidates */
   SCIP_Real* candsscore;            /* scores of all candidates */
   SCIP_VAR* maxtargetvar = NULL;    /* variable with currently highest target value, i.e., product of integer infeasibility and absolute objective */
   SCIP_Real maxtargettarget = -1.0; /* target value of current candidate with highest target value, i.e., product of integer infeasibility and absolute objective */
   SCIP_Real maxtargetval = 0.0;     /* value of current candidate with highest target value, i.e., product of integer infeasibility and absolute objective */
   SCIP_Real maxtargetscore = 0.0;   /* score of current candidate with highest target value, i.e., product of integer infeasibility and absolute objective */
   SCIP_Real currentfrac;            /* fractionality of the current candidate */
   SCIP_Real currenttarget;          /* target value, i.e., product of integer infeasibility and absolute objective, of the current candidate */
   int ncands;
   int i;

   assert( scip != NULL );
   assert( branchrule != NULL );
   assert( result != NULL );

   SCIPdebugMessage("Executing External Branching method of SDP-integer-infeasibility-objective!\n");

   /* Get the external candidates, as we use the score only as a tiebreaker, we aren't interested in the number of
    * variables of different types with maximal score, so these return values are set to NULL. */
   SCIP_CALL( SCIPgetExternBranchCands(scip, &cands, &candssol, &candsscore, &ncands, NULL, NULL, NULL, NULL) );

   assert( ncands > 0 ); /* branchExecext should only be called if the list of external branching candidates is non-empty */

   SCIPdebugMessage("branching candidates for SDP-objective:\n");

   /* iterate over all candidates and find the one with the highest absolute objective times integral infeasibility, use score as tiebreaker */
   for (i = 0; i < ncands; i++)
   {
      /* we skip all continuous variables, since we first want to branch on integral variables */
      if ( SCIPvarGetType(cands[i]) == SCIP_VARTYPE_CONTINUOUS )
      {
         SCIPdebugMessage("skipping continuous variable %s\n", SCIPvarGetName(cands[i]));
         continue;
      }
      /* compute the infeasibility for the integrality constraint */
      currentfrac = SCIPfeasFrac(scip, candssol[i]);
      currenttarget = (currentfrac <= 0.5) ? (currentfrac * REALABS(SCIPvarGetObj(cands[i]))) : ((1.0 - currentfrac) * REALABS(SCIPvarGetObj(cands[i])));

      SCIPdebugMessage("%s, value = %f, objective = %f, objective * integer infeasibility = %f, score = %f\n",
         SCIPvarGetName(cands[i]), candssol[i], SCIPvarGetObj(cands[i]), currenttarget, candsscore[i]);

      /* a candidate is better than the current one if:
       * - the absolute objective * integer infeasibility is (epsilon-)bigger than before or
       * - the absolute objective * integer infeasibility is (epsilon-)equal and the score is (epsilon-)bigger or
       * - both are (epsilon-)equal and the index is smaller */
      if ( SCIPisGT(scip, currenttarget, maxtargettarget) ||
         ( SCIPisEQ(scip, currenttarget, maxtargettarget) && SCIPisGT(scip, candsscore[i], maxtargetscore) ) ||
         ( SCIPisEQ(scip, currenttarget, maxtargettarget) && SCIPisEQ(scip, candsscore[i], maxtargetscore) &&
            SCIPvarGetIndex(cands[i]) < SCIPvarGetIndex(maxtargetvar) ) )
      {
         maxtargetvar = cands[i];
         maxtargettarget = currenttarget;
         maxtargetval = candssol[i];
         maxtargetscore = candsscore[i];
      }
   }

   /* if all variables were continuous, we return DIDNOTRUN and let one of the SCIP branching rules decide */
   if ( maxtargettarget == -1.0 )
   {
      SCIPdebugMessage("Skipping SDP-infobj branching rule since all branching variables are continuous\n");
      *result = SCIP_DIDNOTFIND;
      return SCIP_OKAY;
   }

   assert( SCIPisFeasGE(scip, maxtargettarget, 0.0) );
   assert( maxtargetvar != NULL );

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert( branchruledata != NULL );

   /* If all candidates have objective zero, we look for other variables that are coupled with the candidates and check
    * their objective values if the coupledvars or singlecoupledvars parameter is set to true. */
   if ( SCIPisEQ(scip, maxtargettarget, 0.0) && (branchruledata->coupledvars || branchruledata->singlecoupledvars) )
   {
      SCIP_VAR** vars;
      int nconss;
      SCIP_CONS** conss;
      int nvarsincons;
      SCIP_VAR** varsincons;
      SCIP_Bool** coupledvars; /* is there a constraint coupling candidate i and variable j ? */
      SCIP_Bool** singlecoupledvars; /* is variable j coupled with candidate i AND with no other candidate */
      int** candsincons; /* candsincons[i] gives a list of all candidates (indexed as in cands) appearing in cons i */
      int* ncandsincons; /* ncandsincons[i] gives the length of candsincons[i] */
      SCIP_Bool success;
      int coupledcand;
      SCIP_Real currentobj;
      int cand;
      int candpos;
      int nvars;
      int j;
      int c;
      int v;

      SCIPdebugMessage("All branching candidates have objective 0.0, combined integral infeasibility and objective branching proceeds to check coupled "
                       "variables, updated values for candidates:\n");

      nvars = SCIPgetNVars(scip);
      vars = SCIPgetVars(scip);
      assert( vars != NULL );
      nconss = SCIPgetNConss(scip);
      conss = SCIPgetConss(scip);
      assert( conss != NULL );

      /* allocate memory to save the coupled variables and initialize the arrays */
      SCIP_CALL( SCIPallocBufferArray(scip, &ncandsincons, nconss) );
      for (i = 0; i < nconss; i++)
         ncandsincons[i] = 0;

      SCIP_CALL( SCIPallocBufferArray(scip, &candsincons, nconss) );
      for (i = 0; i < nconss; i++)
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &(candsincons[i]), ncands) );
      }

      SCIP_CALL( SCIPallocBufferArray(scip, &coupledvars, ncands) );
      for (i = 0; i < ncands; i++)
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &(coupledvars[i]), nvars) );
         for (j = 0; j < nvars; j++)
            coupledvars[i][j] = FALSE;
      }
      SCIP_CALL( SCIPallocBufferArray(scip, &varsincons, nvars) );

      /* find all variables that are coupled to a candidate */
      for (c = 0; c < nconss; c++)
      {
         /* first check which candidates appear in which constraints */
         SCIP_CALL( SCIPgetConsNVars(scip, conss[c], &nvarsincons, &success) );
         if ( ! success )
         {
            SCIPdebugMessage("couldn't get variable information from constraint %s, so ignoring it for computing coupled variables\n", SCIPconsGetName(conss[c]));
            continue; /* if we can't get the variables of this constraint, we can't include variables coupled through this constraint */
         }

         /* nothing to do for this constraint if there are no variables (this can happen if all vars are fixed, as the constraint is non-trivial to check) */
         if ( nvarsincons == 0)
            continue;

         SCIP_CALL( SCIPgetConsVars(scip, conss[c], varsincons, nvarsincons, &success) );
         assert( success ); /* we allocated enough memory */
         assert( varsincons != NULL );

         for (v = 0; v < nvarsincons; v++)
         {
            for (cand = 0; cand < ncands; cand++)
            {
               if ( varsincons[v] == cands[cand] )
               {
                  candsincons[c][ncandsincons[c]] = cand;
                  ncandsincons[c]++;
               }
            }
         }

         /* now save which variables are coupled to each candidate by adding all those that appear in this constraint to
          * all candidates appearing in this constraint */
         for (candpos = 0; candpos < ncandsincons[c]; candpos++)
         {
            for (v = 0; v < nvarsincons; v++)
            {
               /* the coupledvars-index corresponding to a variable is its variable index - nvars, because we work on the transformed variables which
                * have indices nvars to 2*nvars - 1, as their indices start after those of the original variables */
               coupledvars[candsincons[c][candpos]][SCIPvarGetIndex(varsincons[v]) - nvars] = TRUE; /*lint !e679*/
            }
         }
      }

      if ( branchruledata->singlecoupledvars )
      {
         /* allocate memory for singlecoupledvars */
         SCIP_CALL( SCIPallocBufferArray(scip, &singlecoupledvars, ncands) );
         for (i = 0; i < ncands; i++)
         {
            SCIP_CALL( SCIPallocBufferArray(scip, &(singlecoupledvars[i]), nvars) );
            for (j = 0; j < nvars; j++)
               singlecoupledvars[i][j] = FALSE;
         }
         /* finally remove all variables, that are coupled to multiple candidates */
         for (v = 0; v < nvars; v++)
         {
            /* we use coupledcand to save if we have already found a candidate this variable is coupled with (otherwise coupledcand = -1), if we find one,
             * we set coupledcand to that index, to easily set the corresponding entry to TRUE if we don't find another candidate it is coupled with */
            coupledcand = -1;
            for (cand = 0; cand < ncands; cand++)
            {
               if ( coupledvars[cand][v] )
               {
                  /* check if this is the first or the second found candidate for this variable */
                  if ( coupledcand == -1 )
                  {
                     /* this is the first candidate this is coupled with, so it might be the only one and we save it to
                      * potentially later set singlecoupledvars to true */
                     coupledcand = cand;
                  }
                  else
                  {
                     /* we found a second candidate, so this variable won't be taken into account for the branching rule, so we set coupledcand
                      * to -2 to not set the corresponding entry in singlecoupledvars to TRUE and continue with the next variable */
                     coupledcand = -2;
                     break;
                  }
               }
            }
            if ( coupledcand > -1 )
            {
               /* as we found exactly one candidate this variable is coupled with, we set the corresponding singlecoupledvars-entry to TRUE */
               singlecoupledvars[coupledcand][v] = TRUE;
            }
         }
      }

      /* iterate over all variables and compute the total absolute objective multiplied of all coupled variables */
      for (cand = 0; cand < ncands; cand++)
      {
         currentobj = 0.0;
         for (v = 0; v < nvars; v++)
         {
            if ( branchruledata->singlecoupledvars && singlecoupledvars[cand][v] ) /*lint !e644*/
               currentobj += REALABS(SCIPvarGetObj(vars[v]));
            else if ( coupledvars[cand][v] )
               currentobj += REALABS(SCIPvarGetObj(vars[v]));
         }

         /* multiply it with the integral infeasibility of the candidate */
         currentfrac = SCIPfeasFrac(scip, candssol[cand]);
         currenttarget = (currentfrac <= 0.5) ? (currentfrac * currentobj) : ((1 - currentfrac) * currentobj);

         assert( SCIPisGE(scip, currentobj, 0.0) );

#ifdef SCIP_DEBUG
         SCIPdebugMessage("candidate %s, coupled with ", SCIPvarGetName(cands[cand]));
         for (v = 0; v < nvars; v++)
         {
            if (coupledvars[cand][v])
               SCIPdebugMessage("%s, ", SCIPvarGetName(vars[v]));
         }
         SCIPdebugMessage("out of those ");
         for (v = 0; v < nvars; v++)
         {
            if (singlecoupledvars[cand][v])
               SCIPdebugMessage("%s, ", SCIPvarGetName(vars[v]));
         }
         SCIPdebugMessage("are only coupled with this candidate, total objective = %f, integral infeasibility = %f, total objective * candidate's fractionality = %f,"
                  "score = %f\n", currentobj, (currentfrac <= 0.5) ? currentfrac : (1 - currentfrac), currenttarget, candsscore[cand]);
#endif

         /* a candidate is better than the current one if:
          * - the absolute objective * integer infeasibility is (epsilon-)bigger than before or
          * - the absolute objective * integer infeasibility is (epsilon-)equal and the score is (epsilon-)bigger or
          * - both are (epsilon-)equal and the index is smaller */
         if ( SCIPisGT(scip, currenttarget, maxtargettarget) ||
             (SCIPisEQ(scip, currenttarget, maxtargettarget) && SCIPisGT(scip, candsscore[cand], maxtargetscore)) ||
             (SCIPisEQ(scip, currenttarget, maxtargettarget) && SCIPisEQ(scip, candsscore[cand], maxtargetscore) &&
                   SCIPvarGetIndex(cands[cand]) < SCIPvarGetIndex(maxtargetvar)) )
         {
            maxtargetvar = cands[cand];
            maxtargettarget = currenttarget;
            maxtargetval = candssol[cand];
            maxtargetscore = candsscore[cand];
         }
      }

      /* free Memory */
      if ( branchruledata->singlecoupledvars )
      {
         for (i = 0; i < ncands; i++)
         {
            SCIPfreeBufferArray(scip, &(singlecoupledvars[i]));
         }
         SCIPfreeBufferArray(scip, &singlecoupledvars);
      }
      SCIPfreeBufferArray(scip, &varsincons);
      for (i = 0; i < ncands; i++)
      {
         SCIPfreeBufferArray(scip, &(coupledvars[i]));
      }
      SCIPfreeBufferArray(scip, &coupledvars);
      for (i = 0; i < nconss; i++)
      {
         SCIPfreeBufferArray(scip, &(candsincons[i]));
      }
      SCIPfreeBufferArray(scip, &candsincons);
      SCIPfreeBufferArray(scip, &ncandsincons);
   }

   /* if the objective values of all integer variables (and all coupled variables, if this settings was used) is zero, skip this branching rule */
   if ( SCIPisGT(scip, maxtargettarget, 0.0) )
   {
      /* branch */
      SCIPdebugMessage("branching on variable %s with value %f, absolute objective * integer infeasibility = %f and score %f\n",
            SCIPvarGetName(maxtargetvar), maxtargetval, maxtargettarget, maxtargetscore);
      SCIP_CALL( SCIPbranchVarVal(scip, maxtargetvar, maxtargetval, NULL, NULL, NULL) );

      *result = SCIP_BRANCHED;
   }
   else
   {
      /* skip */
      *result = SCIP_DIDNOTRUN;
   }

   return SCIP_OKAY;
}

/** destructor of branching rule to free user data (called when SCIP is exiting) */
static
SCIP_DECL_BRANCHFREE(branchFreeSdpinfobjective)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;

   /* free branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   SCIPfreeMemory(scip, &branchruledata);
   SCIPbranchruleSetData(branchrule, NULL);

   return SCIP_OKAY;
}

/*
 * branching rule specific interface methods
 */

/** creates the SDP combined infeasibility and absolute objective branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleSdpinfobjective(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_BRANCHRULE* branchrule;

   /* create branching rule data */
   SCIP_CALL( SCIPallocMemory(scip, &branchruledata) );

   branchrule = NULL;

   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchruleBasic(scip, &branchrule, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY,
         BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST, branchruledata) );

   assert(branchrule != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetBranchruleCopy(scip, branchrule, branchCopySdpinfobjective) );
   SCIP_CALL( SCIPsetBranchruleExecExt(scip, branchrule, branchExecextSdpinfobjective) );
   SCIP_CALL( SCIPsetBranchruleFree(scip, branchrule, branchFreeSdpinfobjective) );

   /* add parameters for the branching rule */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "branching/sdpinfobjective/coupledvars",
         "if all branching candidates have objective zero, should we use the sum of the absolute objectives of all continuous variables coupled with the "
         "candidate through constraints ?",
         &branchruledata->coupledvars, TRUE, DEFAULT_COUPLEDVARS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "branching/sdpinfobjective/singlecoupledvars",
         "if all branching candidates have objective zero, should we use the sum of the absolute objectives of all continuous variables coupled with the "
         "candidate through constraints in which no other candidate appears ?",
         &branchruledata->singlecoupledvars, TRUE, DEFAULT_SINGLECOUPLEDVARS, NULL, NULL) );

   return SCIP_OKAY;
}
