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

/**@file   branch_sdpmostfrac.c
 * @brief  most fractional branching rule for SCIP-SDP
 * @author Tristan Gally
 *
 * Branch on the most fractional variable in the current SDP-relaxation, i.e. the variable maximizing \f$x-\lfloor x \rfloor \f$.
 *
 * Will do nothing for continuous variables, since these are what the external callbacks of the SCIP branching rules are for.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/* #define SCIP_DEBUG */

#include <assert.h>
#include <string.h>

#include "branch_sdpmostfrac.h"

/* turn off lint warnings for whole file: */
/*lint --e{788,818}*/

#define BRANCHRULE_NAME            "sdpmostfrac"
#define BRANCHRULE_DESC            "branch on the most fractional variable of the SDP"
#define BRANCHRULE_PRIORITY        500000
#define BRANCHRULE_MAXDEPTH        -1
#define BRANCHRULE_MAXBOUNDDIST    1.0


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
SCIP_DECL_BRANCHCOPY(branchCopySdpmostfrac)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);

   /* call inclusion method of branchrule */
   SCIP_CALL( SCIPincludeBranchruleSdpmostfrac(scip) );

   return SCIP_OKAY;
}

/** branching execution method for external candidates */
static
SCIP_DECL_BRANCHEXECEXT(branchExecextSdpmostfrac)
{/*lint --e{715}*/
   int i;
   int ncands;
   SCIP_VAR** cands = NULL;
   SCIP_Real* candssol; /* solution values of all candidates */
   SCIP_Real* candsscore; /* scores of all candidates */
   SCIP_Real mostfracfrac; /* fractionality of the current most fractional variable */
   SCIP_Real mostfracscore; /* score of the current most fractional variable */
   SCIP_Real mostfracobj; /* objective of the current most fractional variable */
   SCIP_Real mostfracval; /* value of the current most fractional variable */
   SCIP_VAR* mostfracvar = NULL; /* variable with the highest current fractionality */


   assert( scip != NULL );
   assert( result != NULL );

   SCIPdebugMessage("Executing External Branching method of SDP-mostfrac!\n");

   /* get the external candidates, as we use the score only as a tiebreaker, we aren't interested in the number of
    * variables of different types with maximal score, so these return values are set to NULL */
   SCIP_CALL( SCIPgetExternBranchCands(scip, &cands, &candssol, &candsscore, &ncands, NULL, NULL, NULL, NULL) );

   assert( ncands > 0 ); /* branchExecext should only be called if the list of external branching candidates is non-empty */

#ifdef SCIP_DEBUG
   SCIPdebugMessage("branching candidates for SDP-mostfrac:\n");
   for (i = 0; i < ncands; i++)
      SCIPdebugMessage("%s, value = %f, score = %f\n", SCIPvarGetName(cands[i]), candssol[i], candsscore[i]);
#endif

   mostfracfrac = -1.0;
   mostfracscore = 0.0;
   mostfracval = 0.0;
   mostfracobj = -1.0;

   /* iterate over all solution candidates to find the one with the highest fractionality */
   for (i = 0; i < ncands; i++)
   {
      /* we skip all continuous variables, since we first want to branch on integral variables */
      if ( SCIPvarGetType(cands[i]) == SCIP_VARTYPE_CONTINUOUS )
      {
         SCIPdebugMessage("skipping continuous variable %s\n", SCIPvarGetName(cands[i]));
         continue;
      }

      /* a candidate is better than the current one if:
       * - the fractionality is (feastol-)bigger than before or
       * - the fractionality is (feastol-)equal and the score is (epsilon-)bigger or
       * - the fractionality and score are (feastol-/epsilon-)equal and the objective is (epsilon) bigger or
       * - all three are (feastol-/epsilon-)bigger and the index is smaller */
      if ( SCIPisFeasGT(scip, SCIPfeasFrac(scip, candssol[i]), mostfracfrac) ||
          (SCIPisFeasEQ(scip, SCIPfeasFrac(scip, candssol[i]), mostfracfrac) && SCIPisGT(scip, candsscore[i], mostfracscore)) ||
          (SCIPisFeasEQ(scip, SCIPfeasFrac(scip, candssol[i]), mostfracfrac) && SCIPisEQ(scip, candsscore[i], mostfracscore)
              && SCIPisGT(scip, SCIPvarGetObj(cands[i]), mostfracobj)) ||
          (SCIPisFeasEQ(scip, SCIPfeasFrac(scip, candssol[i]), mostfracfrac) && SCIPisEQ(scip, candsscore[i], mostfracscore)
              && SCIPisEQ(scip, SCIPvarGetObj(cands[i]), mostfracobj) && SCIPvarGetIndex(cands[i]) < SCIPvarGetIndex(mostfracvar)) )
      {
         /* update the current best candidate */
         mostfracfrac = SCIPfeasFrac(scip, candssol[i]);
         mostfracscore = candsscore[i];
         mostfracobj = REALABS(SCIPvarGetObj(cands[i]));
         mostfracval = candssol[i];
         mostfracvar = cands[i];
      }
   }

   /* if all variables were continuous, we return DIDNOTRUN and let one of the SCIP branching rules decide */
   if ( mostfracfrac == -1.0 )
   {
      SCIPdebugMessage("Skipping SDP-mostfrac branching rule since all branching variables are continuous\n");
      *result = SCIP_DIDNOTFIND;
      return SCIP_OKAY;
   }

   assert( mostfracvar != NULL );
   assert( SCIPisFeasGT(scip, mostfracfrac, 0.0) );

   /* branch */
   SCIPdebugMessage("branching on variable %s with value %f and score %f\n", SCIPvarGetName(mostfracvar), mostfracval, mostfracscore);
   SCIP_CALL( SCIPbranchVarVal(scip, mostfracvar, mostfracval, NULL, NULL, NULL) );

   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}

/*
 * branching rule specific interface methods
 */

/** creates the SDP most fractional branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleSdpmostfrac(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_BRANCHRULE* branchrule;

   /* create empty branching rule data */
   branchruledata = NULL;

   branchrule = NULL;

   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchruleBasic(scip, &branchrule, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY,
         BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST, branchruledata) );

   assert(branchrule != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetBranchruleCopy(scip, branchrule, branchCopySdpmostfrac) );
   SCIP_CALL( SCIPsetBranchruleExecExt(scip, branchrule, branchExecextSdpmostfrac) );

   return SCIP_OKAY;
}
