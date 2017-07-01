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

/**@file   prop_sdpredcost.c
 * @brief  reduced cost / dual fixing for SDPs
 * @author Tristan Gally
 */

/*#define SCIP_DEBUG*/

#include "prop_sdpredcost.h"
#include "scip/def.h"                        /* for SCIP_Real, _Bool, ... */
#include "relax_sdp.h"                       /* to get relaxation value */
#include "sdpi/sdpi.h"                       /* to get values of primal variables */

#include <string.h>
#include <assert.h>                          /*lint !e451*/

/* turn off lint warnings for whole file: */
/*lint --e{788,818}*/

/**@name Propagator properties
 * @{
 */

#define PROP_NAME                   "sdpredcost"
#define PROP_DESC                   "sdp reduced cost strengthening propagator"
#define PROP_TIMING                 SCIP_PROPTIMING_AFTERLPLOOP
#define PROP_PRIORITY               +1000000 /**< propagator priority */
#define PROP_FREQ                   1        /**< propagator frequency */
#define PROP_DELAY                  FALSE    /**< Should propagation method be delayed, if other propagators found reductions? */
#define DEFAULT_SDPRCBIN            TRUE     /**< Should sdp reduced cost fixing be executed for binary variables? */
#define DEFAULT_SDPRCINTCONT        TRUE     /**< Should sdp reduced cost fixing be executed for integer and continuous variables? */

/**@} */

/** propagator data */
struct SCIP_PropData
{
   /* these could also be freshly allocated for each node, but allocating them only once saves time */
   SCIP_Real*            lbvarvals;          /**< array where the current values of the primal variables corresponding to dual lower variable-bounds are saved */
   SCIP_Real*            ubvarvals;          /**< array where the current values of the primal variables corresponding to dual upper variable-bounds are saved */
   int                   nvars;              /**< number of variables and therefore also length of lbvarvals and ubvarvals */
   SCIP_Bool             forbins;            /**< should sdp reduced cost fixing be executed for binary variables? */
   SCIP_Bool             forintconts;        /**< should sdp reduced cost fixing be executed for integer and continuous variables? */
};

/** reduced cost fixing for binary variables
 *
 *  If the corresponding primal variable for the lower bound is bigger than the cutoff bound minus the
 *  current relaxation value, then the variable can be fixed to zero, if the primal variable for the upper bound is bigger than this value, then it
 *  can be fixed to one.
 */
static
SCIP_RETCODE sdpRedcostFixingBinary(
   SCIP*                 scip,               /**< pointer to SCIP data structure */
   SCIP_VAR*             var,                /**< variable to propagate */
   SCIP_Real             primallbval,        /**< value of the primal variable corresponding to the lower bound */
   SCIP_Real             primalubval,        /**< value of the primal variable corresponding to the upper bound */
   SCIP_Real             cutoffbound,        /**< current cutoffbound in SCIP */
   SCIP_Real             relaxval,           /**< optimal objective value of the current relaxation */
   SCIP_RESULT*          result              /**< pointer to return result */
   )
{
   assert( scip != NULL );
   assert( var != NULL );
   assert( result != NULL );

   /* skip binary variable if it is locally fixed */
   if (SCIPvarGetLbLocal(var) > 0.5 || SCIPvarGetUbLocal(var) < 0.5)
   {
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   /* check if variable can be fixed to zero */
   if ( SCIPisGT(scip, primallbval, cutoffbound - relaxval) )
   {
      SCIP_CALL( SCIPchgVarUb(scip, var, 0.0) );
      SCIPdebugMessage("Variable %s fixed to zero by reduced cost fixing ! \n", SCIPvarGetName(var));
      *result = SCIP_REDUCEDDOM;

      /* check if we would also have to fix the variable to one, in that case, we can cut the node off, as there can't be a new optimal solution */
      if ( SCIPisGT(scip, primalubval, cutoffbound - relaxval) )
      {
         *result = SCIP_CUTOFF;
      }

      return SCIP_OKAY;
   }

   /* check if variable can be fixed to one */
   if ( SCIPisGT(scip, primalubval, cutoffbound - relaxval) )
   {
      SCIP_CALL( SCIPchgVarLb(scip, var, 1.0) );
      SCIPdebugMessage("Variable %s fixed to one by reduced cost fixing ! \n", SCIPvarGetName(var));
      *result = SCIP_REDUCEDDOM;
      return SCIP_OKAY;
   }

   *result = SCIP_DIDNOTFIND;
   return SCIP_OKAY;
}

/** reduced cost fixing for non-binary variables
 *
 *  We propagate the new bounds
 *
 *  \f$ y_j \leq \ell_j + \frac{v_{CO} - \bar{v}}{\bar{X}_{n+m+j,n+m+j}} \f$,
 *
 *  \f$ y_j \geq u_j - \frac{v_{CO} - \bar{v}}{\bar{X}_{n+j,n+j}} \f$
 *
 *  where \f$\bar{v}\f$ is the value of the current relaxation, \f$v_{CO}\f$ is the cutoffbound and \f$\bar{X}_{n+m+j,n+m+j}\f$ the value of the
 *  corresponding primal solution
 */
static
SCIP_RETCODE sdpRedcostFixingIntCont(
   SCIP*                 scip,               /**< pointer to SCIP data structure */
   SCIP_VAR*             var,                /**< variable to propagate */
   SCIP_Real             primallbval,        /**< value of the primal variable corresponding to the lower bound */
   SCIP_Real             primalubval,        /**< value of the primal variable corresponding to the upper bound */
   SCIP_Real             cutoffbound,        /**< current cutoffbound in SCIP */
   SCIP_Real             relaxval,           /**< optimal objective value of the current relaxation */
   SCIP_RESULT*          result              /**< pointer to return result */
   )
{
   SCIP_Real lb;
   SCIP_Real ub;

   assert( scip != NULL );
   assert( var != NULL );
   assert( result != NULL );

   *result = SCIP_DIDNOTFIND;

   /* compute the new lower and upper bound, if we divide by zero (checking > 0 is sufficient, as the variabls are non-negative), the bounds are infinity */
   ub = SCIPisGT(scip, primallbval, 0.0) ? SCIPvarGetLbLocal(var) + (cutoffbound - relaxval) / primallbval : SCIPinfinity(scip);
   lb = SCIPisGT(scip, primalubval, 0.0) ? SCIPvarGetUbLocal(var) - (cutoffbound - relaxval) / primalubval : -SCIPinfinity(scip);

   /* if either bound is infinite, we set it to the corresponding SCIP value */
   if ( SCIPisInfinity(scip, ub) )
      ub = SCIPinfinity(scip);
   else if ( SCIPisInfinity(scip, -ub) )
      ub = -SCIPinfinity(scip);
   if ( SCIPisInfinity(scip, lb) )
      lb = SCIPinfinity(scip);
   else if ( SCIPisInfinity(scip, -lb) )
      lb = -SCIPinfinity(scip);


   /* if after propagation the upper bound is less than the lower bound, the current node is infeasible */
   if ( SCIPisLT(scip, ub, lb) || SCIPisLT(scip, ub, SCIPvarGetLbLocal(var)) || SCIPisLT(scip, SCIPvarGetUbLocal(var), lb) )
   {
      SCIPdebugMessage("Infeasibility of current node detected by prop_sdpredcost! Updated bounds for variable %s: lb = %f > %f = ub !\n",
            SCIPvarGetName(var), SCIPisGT(scip, lb, SCIPvarGetLbLocal(var))? lb : SCIPvarGetLbLocal(var),
            SCIPisLT(scip, ub, SCIPvarGetLbLocal(var)) ? ub : SCIPvarGetUbLocal(var) );
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   /* if the new upper bound is an enhancement, update it */
   if ( SCIPisLT(scip, ub, SCIPvarGetUbLocal(var)) )
   {
      SCIPdebugMessage("Changing upper bound of variable %s from %f to %f because of prop_sdpredcost \n",
            SCIPvarGetName(var), SCIPvarGetUbLocal(var), ub);
      SCIP_CALL( SCIPchgVarUb(scip, var, ub) );
      *result =  SCIP_REDUCEDDOM;
   }

   /* if the new lower bound is an enhancement, update it */
   if ( SCIPisGT(scip, lb, SCIPvarGetLbLocal(var)) )
   {
      SCIPdebugMessage("Changing lower bound of variable %s from %f to %f because of prop_sdpredcost \n",
            SCIPvarGetName(var), SCIPvarGetLbLocal(var), lb);
      SCIP_CALL( SCIPchgVarLb(scip, var, lb) );
      *result =  SCIP_REDUCEDDOM;
   }

   return SCIP_OKAY;
}

/** reduced cost propagation method */
static
SCIP_DECL_PROPEXEC(propExecSdpredcost)
{/*lint --e{715}*/
   int v;
   int nvars;
   SCIP_VAR** vars;
   SCIP_RELAX* relax;
   SCIP_RESULT varresult;
   SCIP_Real cutoffbound;
   SCIP_Real relaxval;
   SCIP_Bool sdpsolved;
   SCIP_PROPDATA* propdata;
   int length;

   SCIPdebugMessage("Calling propExecSdpredcost \n");

   assert( scip != NULL );
   assert( prop != NULL );
   assert( result != NULL );

   if ( SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING )
   {
      /* we can't run before the relaxator is properly initialized */
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   relax = SCIPfindRelax(scip, "SDP"); /* get SDP relaxation handler */
   assert( relax != NULL );

   /* we can only propagate for the last node for which the SDP was solved */
   if ( SCIPrelaxSdpGetSdpNode(relax) != SCIPnodeGetNumber(SCIPgetCurrentNode(scip)) )
   {
      SCIPdebugMessage("Stopped propExecRedcost because current SDP-relaxation doesn't belong to the node the propagator was called for!\n");
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   /* we can only propagate if the SDP in the last node was solved in its original formulation */
   if ( ! SCIPrelaxSdpSolvedOrig(relax) )
   {
      SCIPdebugMessage("Stopped propExecRedcost because current SDP-relaxation was solved using a penalty formulation!\n");
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPrelaxSdpRelaxVal(relax, &sdpsolved, &relaxval) );
   if ( ! sdpsolved )
   {
      SCIPdebugMessage("Stopped propExecRedcost because SDP-relaxation wasn't properly solved!\n");
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   propdata = SCIPpropGetData(prop);
   assert( propdata != NULL );

   *result = SCIP_DIDNOTFIND;

   nvars = SCIPgetNVars(scip);
   vars = SCIPgetVars(scip);

   cutoffbound = SCIPgetCutoffbound(scip);

   length = nvars;

   SCIP_CALL( SCIPrelaxSdpGetPrimalBoundVars(relax, propdata->lbvarvals, propdata->ubvarvals, &length) );

   assert( length == nvars ); /* we should get exactly one value for lower and upper bound-variable per variable in scip */

   for (v = 0; v < nvars; v++)
   {
      if ( SCIPvarIsBinary(vars[v]) && propdata->forbins )
      {
         SCIP_CALL( sdpRedcostFixingBinary(scip, vars[v], propdata->lbvarvals[v], propdata->ubvarvals[v], cutoffbound, relaxval, &varresult) );

         if ( varresult == SCIP_REDUCEDDOM )
            *result = SCIP_REDUCEDDOM;
      }
      else if ( (! SCIPvarIsBinary(vars[v])) && propdata->forintconts )
      {
         SCIP_CALL( sdpRedcostFixingIntCont(scip, vars[v], propdata->lbvarvals[v], propdata->ubvarvals[v], cutoffbound, relaxval, &varresult) );

         if ( varresult == SCIP_REDUCEDDOM )
            *result = SCIP_REDUCEDDOM;
      }
   }

   return SCIP_OKAY;
}

/** free the propagator data */
static
SCIP_DECL_PROPFREE(propFreeSdpredcost)
{/*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   propdata = SCIPpropGetData(prop);
   assert( propdata != NULL );
   SCIPfreeMemory(scip, &propdata);

   SCIPpropSetData(prop, NULL);

   return SCIP_OKAY;
}

/** free memory for the primal variable values */
static
SCIP_DECL_PROPEXIT(propExitSdpredcost)
{
   SCIP_PROPDATA* propdata;

   propdata = SCIPpropGetData(prop);
   assert( propdata != NULL );

   SCIPfreeBlockMemoryArrayNull(scip, &(propdata->lbvarvals), propdata->nvars);
   SCIPfreeBlockMemoryArrayNull(scip, &(propdata->ubvarvals), propdata->nvars);

   return SCIP_OKAY;
}

/** allocate memory for the primal variable values */
static
SCIP_DECL_PROPINITSOL(propInitsolSdpredcost)
{
   SCIP_PROPDATA* propdata;
   int nvars;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);
   nvars = SCIPgetNVars(scip);

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(propdata->lbvarvals), nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(propdata->ubvarvals), nvars) );
   propdata->nvars = nvars;

   return SCIP_OKAY;
}

/** copy method for propagator plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PROPCOPY(propCopySdpredcost)
{
   assert( scip != NULL );
   assert( prop != NULL );
   assert( strcmp(SCIPpropGetName(prop), PROP_NAME) == 0 );

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludePropSdpredcost(scip) );

   return SCIP_OKAY;
}

/** creates the Sdpredcost propagator and includes it in SCIP */
SCIP_RETCODE SCIPincludePropSdpredcost(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROPDATA* propdata = NULL;
   SCIP_PROP* prop;

   /* create propagator data */
   SCIP_CALL( SCIPallocMemory(scip, &propdata) );
   propdata->nvars = 0;

   /* include propagator */
   SCIP_CALL( SCIPincludePropBasic(scip, &prop, PROP_NAME, PROP_DESC, PROP_PRIORITY, PROP_FREQ, PROP_DELAY, PROP_TIMING,
         propExecSdpredcost, propdata) );
   assert( prop != NULL );

   /* set optional callbacks via setter functions */
   SCIP_CALL( SCIPsetPropCopy(scip, prop, propCopySdpredcost) );
   SCIP_CALL( SCIPsetPropInitsol(scip, prop, propInitsolSdpredcost) );
   SCIP_CALL( SCIPsetPropExit(scip, prop, propExitSdpredcost) );
   SCIP_CALL( SCIPsetPropFree(scip, prop, propFreeSdpredcost) );

   /* add additional parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "propagating/sdpredcost/forbins", "Should SDP reduced cost fixing be executed for binary variables?",
         &(propdata->forbins), TRUE, DEFAULT_SDPRCBIN, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "propagating/sdpredcost/forintconts", "Should SDP reduced cost fixing be executed for integer and continuous variables?",
         &(propdata->forintconts), TRUE, DEFAULT_SDPRCINTCONT, NULL, NULL) );

   return SCIP_OKAY;
}
