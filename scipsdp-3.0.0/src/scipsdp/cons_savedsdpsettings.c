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

/**@file   cons_savedsdpsettings.c
 * @brief  constraint handler for saving SDP settings
 * @author Tristan Gally
 *
 * A constraint that is always feasible which can be used to save and recover settings used
 * to solve the SDP-relaxation at the current node.
 */

#include "cons_savedsdpsettings.h"
#include "scip/def.h"                   /* for SCIP_Real, _Bool, ... */
#include <string.h>                     /* for NULL, strcmp */

/* turn off lint warnings for whole file: */
/*lint --e{788,818}*/

/* constraint handler properties */
#define CONSHDLR_NAME          "Savedsdpsettings"
#define CONSHDLR_DESC          "constraint handler to store SDP settings for each node"
#define CONSHDLR_ENFOPRIORITY         0 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY        0 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

/** constraint data to store settings used to solve parent node */
struct SCIP_ConsData
{
   SCIP_SDPSOLVERSETTING settings;           /**< pointer to save settings */
};

/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteSavedsdpsettings)
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( cons != NULL );
   assert( consdata != NULL );
   assert( *consdata != NULL );

   SCIPdebugMessage("Deleting store node data constraint: <%s>.\n", SCIPconsGetName(cons));

   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpSavedsdpsettings)
{/*lint --e{715}*/
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   /* do nothing */
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFORELAX(consEnforelaxSavedsdpsettings)
{/*lint --e{715}*/
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   /* do nothing */
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsSavedsdpsettings)
{/*lint --e{715}*/
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   /* do nothing */
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for primal solutions */
static
SCIP_DECL_CONSCHECK(consCheckSavedsdpsettings)
{/*lint --e{715}*/
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   /* do nothing */
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockSavedsdpsettings)
{/*lint --e{715}*/
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   /* do nothing */
   return SCIP_OKAY;
}


/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopySavedsdpsettings)
{
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( valid != NULL );

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrSavedsdpsettings(scip) );

   *valid = TRUE;

   return SCIP_OKAY;
}


/** constraint copying method of constraint handler */
static
SCIP_DECL_CONSCOPY(consCopySavedsdpsettings)
{  /*lint --e{715}*/

   if ( name )
   {
      SCIP_CALL( createConsSavedsdpsettings(scip, cons, name, SCIPconsSavedsdpsettingsGetSettings(sourcescip, sourcecons)) );
   }
   else
   {
      SCIP_CALL( createConsSavedsdpsettings(scip, cons, SCIPconsGetName(sourcecons), SCIPconsSavedsdpsettingsGetSettings(sourcescip, sourcecons)) );
   }

   *valid = TRUE;
   return SCIP_OKAY;
}


/** include Savedsdpsettings constraint handler */
SCIP_RETCODE SCIPincludeConshdlrSavedsdpsettings(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;

   /* include constraint handler */
   conshdlr = NULL;
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpSavedsdpsettings, consEnfopsSavedsdpsettings, consCheckSavedsdpsettings, consLockSavedsdpsettings,
         NULL) );
   assert( conshdlr != NULL );

   /* set additional callbacks */
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteSavedsdpsettings) );
   SCIP_CALL( SCIPsetConshdlrEnforelax(scip, conshdlr, consEnforelaxSavedsdpsettings) );
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopySavedsdpsettings, consCopySavedsdpsettings) );

   return SCIP_OKAY;
}



/*
 * External functions
 */

/** create a savedsdpsettings constraint, i.e. save the current settings for the SDP-relaxation of this node */
SCIP_RETCODE createConsSavedsdpsettings(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_SDPSOLVERSETTING settings            /**< settings to save */
   )
{
   SCIP_CONSDATA* consdata = NULL;
   SCIP_CONSHDLR* conshdlr;

   assert( scip != NULL );
   assert( name != NULL );

   /* find the node data constraint handler */
   conshdlr = SCIPfindConshdlr(scip, "Savedsdpsettings");
   if ( conshdlr == NULL )
   {
      SCIPerrorMessage("savedsdpsettings constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) );
   consdata->settings = settings;

   SCIPdebugMessage("Creating savedsdpsettings constraint <%s>.\n", name);

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, FALSE, FALSE, FALSE, FALSE, FALSE,
         TRUE, FALSE, TRUE, FALSE, TRUE));

   return SCIP_OKAY;
}

/** get the settings used to solve the SDP relaxation in this node */
SCIP_SDPSOLVERSETTING SCIPconsSavedsdpsettingsGetSettings(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint to get starting point for */
   )
{
   SCIP_CONSDATA* consdata;

   assert( scip != NULL );
   assert( cons != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   return consdata->settings;
}
