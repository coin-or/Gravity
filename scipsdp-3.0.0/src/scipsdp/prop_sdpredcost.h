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

/**@file   prop_sdpredcost.h
 * @brief  reduced cost / dual fixing for SDPs
 * @author Tristan Gally
 *
 *  Propagates bounds
 *
 *  \f$ y_j \leq \ell_j + \frac{v_{CO} - \bar{v}}{\bar{X}_{n+m+j,n+m+j}} \f$,
 *
 *  \f$ y_j \geq u_j - \frac{v_{CO} - \bar{v}}{\bar{X}_{n+j,n+j}} \f$
 *
 *  where \f$\bar{v}\f$ is the value of the current SDP-relaxation, \f$v_{CO}\f$ is the cutoffbound and \f$\bar{X}_{n+m+j,n+m+j}\f$ the value of the
 *  corresponding primal solution.
 */

#ifndef __SCIP_PROP_SDPREDCOST_H_
#define __SCIP_PROP_SDPREDCOST_H_

#include <scip/scip.h>

#ifdef __cplusplus
extern "C" {
#endif

/** creates the Sdpredcost propagator and includes it in SCIP */
SCIP_RETCODE SCIPincludePropSdpredcost(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
