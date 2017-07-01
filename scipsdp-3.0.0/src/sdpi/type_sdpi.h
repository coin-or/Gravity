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

/**@file   type_sdpi.h
 * @brief  type definitions for specific SDP-solver interfaces
 * @author Tristan Gally
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_SDPI_H__
#define __SCIP_TYPE_SDPI_H__

#ifdef __cplusplus
extern "C" {
#endif

/* for now, we reuse the enums SCIP_OBJSEN, SCIP_PRICING, and SCIP_BASESTAT from the LPI */
#include "lpi/type_lpi.h"

/** SDP-solver parameters */
enum SCIP_SDPParam
{
   SCIP_SDPPAR_EPSILON        = 0,      /**< absolute tolerance */
   SCIP_SDPPAR_GAPTOL         = 1,      /**< convergence tolerance */
   SCIP_SDPPAR_FEASTOL        = 2,      /**< feasibility tolerance */
   SCIP_SDPPAR_SDPSOLVERFEASTOL = 3,    /**< feasibility tolerance for SDP-solver */
   SCIP_SDPPAR_OBJLIMIT       = 4,      /**< objective limit, if the SDP-solver computes a lower bound for the minimzation
                                         *   problem that is bigger than this, it may stop */
   SCIP_SDPPAR_SDPINFO        = 5,      /**< Should the SDP-solver output information to the screen? */
   SCIP_SDPPAR_SLATERCHECK    = 6,      /**< Should the slater condition for the dual problem be checked before solving each SDP ? */
   SCIP_SDPPAR_PENALTYPARAM   = 7,      /**< the starting penalty parameter Gamma used for the penalty formulation if the SDP-solver did not converge */
   SCIP_SDPPAR_MAXPENALTYPARAM= 8,      /**< the maximum penalty parameter Gamma used for the penalty formulation if the SDP-solver did not converge */
   SCIP_SDPPAR_NPENALTYINCR   = 9,      /**< maximum number of times the penalty parameter will be increased if penalty formulation failed */
   SCIP_SDPPAR_LAMBDASTAR     = 10,     /**< the parameter lambda star used by SDPA to set the initial point */
   SCIP_SDPPAR_NTHREADS       = 11      /**< number of threads the SDP solver should use, currently only supported for MOSEK (-1 = number of cores) */
};
typedef enum SCIP_SDPParam SCIP_SDPPARAM;

/** SDP-solver settings used */
enum SCIP_SDPSolverSetting
{
   SCIP_SDPSOLVERSETTING_UNSOLVED= -1,  /**< problem was not solved */
   SCIP_SDPSOLVERSETTING_PENALTY = 0,   /**< penalty formulation */
   SCIP_SDPSOLVERSETTING_FAST    = 1,   /**< fastest settings */
   SCIP_SDPSOLVERSETTING_MEDIUM  = 2,   /**< medium settings */
   SCIP_SDPSOLVERSETTING_STABLE  = 3    /**< most stable settings */
};
typedef enum SCIP_SDPSolverSetting SCIP_SDPSOLVERSETTING;

/** SDP-solver settings and slater */
enum SCIP_SDPSlaterSetting
{
   SCIP_SDPSLATERSETTING_NOINFO           = 0, /**< Slater check failed or problem not given to solver */
   SCIP_SDPSLATERSETTING_STABLEWSLATER    = 1, /**< solved with fastest settings and primal and dual Slater holding */
   SCIP_SDPSLATERSETTING_UNSTABLEWSLATER  = 2, /**< solved with stable settings and primal and dual Slater holding */
   SCIP_SDPSLATERSETTING_PENALTYWSLATER   = 3, /**< solved with penalty formulation and primal and dual Slater holding */
   SCIP_SDPSLATERSETTING_BOUNDEDWSLATER   = 4, /**< bound computed via penalty approach with primal and dual Slater holding */
   SCIP_SDPSLATERSETTING_UNSOLVEDWSLATER  = 5, /**< unsolved with primal and dual Slater holding */
   SCIP_SDPSLATERSETTING_STABLENOSLATER   = 6, /**< solved with fastest setting and either primal or dual Slater not holding */
   SCIP_SDPSLATERSETTING_UNSTABLENOSLATER = 7, /**< solved with stable settings and either primal or dual Slater not holding */
   SCIP_SDPSLATERSETTING_PENALTYNOSLATER  = 8, /**< solved with penalty formulation and either primal or dual Slater not holding */
   SCIP_SDPSLATERSETTING_BOUNDEDNOSLATER  = 9, /**< bound computed via penalty approach with either primal or dual Slater not holding  */
   SCIP_SDPSLATERSETTING_UNSOLVEDNOSLATER = 10,/**< unsolved with either primal or dual Slater not holding */
   SCIP_SDPSLATERSETTING_STABLEINFEASIBLE = 11, /**< solved with fastest setting and dual Slater check showing infeasibility */
   SCIP_SDPSLATERSETTING_UNSTABLEINFEASIBLE= 12, /**< solved with stable settings and dual Slater check showing infeasibility */
   SCIP_SDPSLATERSETTING_PENALTYINFEASIBLE= 13, /**< solved with penalty formulation and dual Slater check showing infeasibility */
   SCIP_SDPSLATERSETTING_BOUNDEDINFEASIBLE= 14, /**< bound computed via penalty approach with dual Slater check showing infeasibility */
   SCIP_SDPSLATERSETTING_UNSOLVEDINFEASIBLE= 15 /**< unsolved with dual Slater check showing infeasibility */
};
typedef enum SCIP_SDPSlaterSetting SCIP_SDPSLATERSETTING;

/** SDP-solver settings used */
enum SCIP_SDPSlater
{
   SCIP_SDPSLATER_INF    = -2, /**< problem is infeasible */
   SCIP_SDPSLATER_NOINFO = -1, /**< check for Slater condition failed */
   SCIP_SDPSLATER_NOT    = 0,  /**< Slater condition does not hold */
   SCIP_SDPSLATER_HOLDS  = 1   /**< Slater condition holds */
};
typedef enum SCIP_SDPSlater SCIP_SDPSLATER;

typedef struct SCIP_SDPi SCIP_SDPI;                 /**< solver independent SDP interface */

#ifdef __cplusplus
}
#endif

#endif
