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

/**@file   cons_savedsdpsettings.h
 * @brief  constraint handler for saving SDP settings
 * @author Tristan Gally
 *
 * A constraint that is always feasible which can be used to save and recover settings used
 * to solve the SDP-relaxation at the current node.
 */

#ifndef __SCIP_CONS_SAVEDSDPSETTINGS_H_
#define __SCIP_CONS_SAVEDSDPSETTINGS_H_

#include "scip/scip.h"
#include "sdpi/type_sdpi.h"

#ifdef __cplusplus
extern "C" {
#endif

/** include Savedsdpsettings constraint handler */
EXTERN
SCIP_RETCODE SCIPincludeConshdlrSavedsdpsettings(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** create a savedsdpsettings constraint, i.e. save the current settings for the SDP-relaxation of this node */
EXTERN
SCIP_RETCODE createConsSavedsdpsettings(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_SDPSOLVERSETTING settings            /**< settings to save */
   );

/** get the settings used to solve the SDP relaxation in this node */
EXTERN
SCIP_SDPSOLVERSETTING SCIPconsSavedsdpsettingsGetSettings(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint to get starting point for */
   );

#ifdef __cplusplus
}
#endif

#endif
