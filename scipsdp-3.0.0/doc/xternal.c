/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/* This file is part of SCIPSDP - a solving framework for mixed-integer      */
/* semidefinite programms based on SCIP.                                     */
/*                                                                           */
/* Copyright (C) 2011-2013 Discrete Optimization, TU Darmstadt               */
/*                         EDOM, FAU Erlangen-Nürnberg                       */
/*               2014-2016 Discrete Optimization, TU Darmstadt               */
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
/* Copyright (C) 2002-2016 Zuse Institute Berlin                             */
/* SCIP is distributed under the terms of the SCIP Academic Licence,         */
/* see file COPYING in the SCIP distribution.                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   xternal.c
 * @brief  main document page
 * @author Tristan Gally
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/**@mainpage Overview
 *
 * @version 3.0.0
 * @author Tristan Gally, Marc Pfetsch; Sonja Mars, Lars Schewe
 * @date 2011-2017
 *
 * SCIP-SDP is a plugin for SCIP to solve mixed integer semidefinite programs (MISDPs) of the form
 *
 *   \f{equation*}{
 * 	\begin{aligned}
 *      \inf \quad & b^T y && \\
 *      \mbox{s.t.} \quad & \sum_{i = 1}^m A_i y_i - A_0 \succeq 0,&& \\
 *	& y_i \in \mathbb{Z} && \forall \ i \in \mathcal{I}.
 *	\end{aligned}
 *   \f}
 *
 * It combines the branch-and-bound framework of SCIP with interior-point SDP-solvers to solve MISDPs using either a
 * nonlinear branch-and-bound approach or an outer-approximation-based cutting-plane approach. In addition to providing
 * a constraint handler for SDP-constraints and a relaxator to solve continuous SDP-relaxations using interior-point
 * solvers, SCIPSDP adds several heuristics and propagators to SCIP. The MISDPs can be read in using either an extended
 * SDPA-format or the CBF-format. To use the nonlinear branch-and-bound approach one of the following SDP-solvers needs
 * to be installed:
 *
 * - DSDP
 * - SDPA
 * - MOSEK
 *
 * The solution process of interior-point methods for SDPs is highly dependent on the Slater condition. One of the main
 * purposes of the code is handling cases where the Slater condition does not hold using a penalty approach. However, in
 * some cases the SDP-solvers may still fail because of numerical difficulties or even return wrong results, which cannot
 * be compensated. For this purpose there is the possibility to check the Slater condition for the primal and dual problem
 * before the solution of each SDP by setting a SCIP parameter, for details see the parameters tab.
 */

/** @page PARAMETERS Additional Parameters
 * The following important parameters (with these default values) were added:
 *
 * <table>
 * <tr><td>relaxing/SDP/freq = 1</td> <td>set this to -1 and lp/solvefreq to 1 to solve LP relaxations with eigenvector cuts</td></tr>
 * <tr><td>constraints/SDP/threads = 1</td> <td>number of threads used for openblas; only available with OpenBLAS/OpenMP and compile option OMP=TRUE (default for SDPS=sdpa)</td></tr>
 * <tr><td>relaxing/SDP/sdpsolverthreads = -1</td> <td>number of threads the SDP solver should use (-1 = number of cores); currently only supported for MOSEK</td></tr>
 * <tr><td>relaxing/SDP/displaystatistics = TRUE</td> <td>Should statistics about SDP iterations and solver settings/success be printed after quitting SCIP-SDP ?</td></tr>
 * <tr><td>relaxing/SDP/slatercheck = 0</td> <td>Should the Slater condition for the primal and dual problem be checked ahead of solving each SDP? [0: no, 1: yes and output statistics, 2: yes and print warning for every problem not satisfying primal and dual Slater condition]</td></tr>
 * <tr><td>relaxing/SDP/sdpinfo = FALSE</td> <td>Should output of the SDP-Solver be printed to the console?</td></tr>
 * <tr><td>relaxing/SDP/sdpsolvergaptol = 1e-04</td> <td>sets the tolerance for the duality gap in the SDP-Solver</td></tr>
 * <tr><td>relaxing/SDP/sdpsolverfeastol = 1e-06 (DSDP,SDPA) / 1e-07 (MOSEK)</td> <td>feasibility tolerance for the SDP-Solver (should be less or equal to numerics/feastol)</td></tr>
 * <tr><td>relaxing/SDP/penaltyparam = -1</td> <td>the starting value of the penalty parameter Gamma used for the penalty formulation if the SDP solver didn't converge; set this to a negative value to compute the parameter depending on the given problem</td></tr>
 * <tr><td>relaxing/SDP/lambdastar = -1</td> <td>the parameter lambda star used by SDPA to set the initial point; set this to a negative value to compute the parameter depending on the given problem</td></tr>
 * <tr><td>branching/sdpinfobjective/priority = 2e+06</td> <td>priority of combined infeasibility/objective branching rule; branching rule with highest priority is used</td></tr>
 * <tr><td>branching/sdpinfobjective/coupledvars = FALSE</td> <td>If all branching candidates have objective zero, should we use the sum of the absolute objectives of all continuous variables coupled with the candidate through constraints?</td></tr>
 * <tr><td>branching/sdpinfobjective/singlecoupledvars = FALSE</td> <td>If all branching candidates have objective zero, should we use the sum of the absolute objectives of all continuous variables coupled with the candidate through constraints in which no other candidate appears?</td></tr>
 * <tr><td>branching/sdpmostfrac/priority = 5e+05</td> <td>priority of most fractional (largest fractional part) branching rule; branching rule with highest priority is used</td></tr>
 * <tr><td>branching/sdpmostinf/priority = 1e+06</td> <td>priority of most infeasible (fractional part nearest to 0.5) branching rule; branching rule with highest priority is used</td></tr>
 * <tr><td>branching/sdpobjective/priority = 1.5e+06</td> <td>priority of highest absolute objective branching rule; branching rule with highest priority is used</td></tr>
 * <tr><td>branching/sdpobjective/coupledvars = FALSE</td> <td>If all branching candidates have objective zero, should we use the sum of the absolute objectives of all continuous variables coupled with the candidate through constraints?</td></tr>
 * <tr><td>branching/sdpobjective/singlecoupledvars = FALSE</td> <td>If all branching candidates have objective zero, should we use the sum of the absolute objectives of all continuous variables coupled with the candidate through constraints in which no other candidate appears?</td></tr>
 * <tr><td>constraints/SDP/diaggezerocuts = TRUE</td> <td>Should linear cuts enforcing the non-negativity of diagonal entries of SDP-matrices be added?</td></tr>
 * <tr><td>constraints/SDP/diagzeroimplcuts = TRUE</td> <td>Should linear cuts enforcing the implications of diagonal entries of zero in SDP-matrices be added?</td></tr>
 * <tr><td>display/sdpfastsettings/active = 0</td> <td>Should the percentage of SDP-relaxations solved with the fastest setting (SDPA) or the default formulation (DSDP) be displayed in the console? [0: off, 1: auto, 2:on]</td></tr>
 * <tr><td>display/sdppenalty/active = 0</td> <td>Should the percentage of SDP-relaxations solved using a penalty formulation be displayed in the console? [0: off, 1: auto, 2:on]</td></tr>
 * <tr><td>display/sdpunsolved/active = 1</td> <td>Should the percentage of SDP-relaxations that could not be solved be displayed in the console? [0: off, 1: auto, 2:on]</td></tr>
 * <tr><td>heuristics/sdpfracdiving/freq = -1</td> <td>set this to 0 or more to enable a fractional diving heuristic for SDPs</td></tr>
 * <tr><td>heuristics/sdprand/freq = 1</td> <td>set this to -1 to disable the randomized rounding heuristic</td></tr>
 * <tr><td>heuristics/sdprand/generalints = FALSE</td> <td>Should randomized rounding also be applied if there are general integer variables and not only binary variables ?</td></tr>
 * <tr><td>heuristics/sdprand/nrounds = 2</td> <td>number of rounding rounds</td></tr>
 * <tr><td>propagating/sdp-obbt/freq = -1</td> <td>set this to 0 or more to enable optimization-based bound tightening using SDP-relaxations</td></tr>
 * <tr><td>propagating/sdp-obbt/propcont = TRUE</td> <td>Should optimization-based bound tightening be performed for continuous variables ?</td></tr>
 * <tr><td>propagating/sdp-obbt/propbin = FALSE</td> <td>Should optimization-based bound tightening be performed for binary variables ?</td></tr>
 * <tr><td>propagating/sdpredcost/freq = 1</td> <td>set this to -1 to disable reduced cost fixing for SDPs</td></tr>
 * <tr><td>propagating/sdpredcost/forbins = TRUE</td> <td>Should SDP reduced cost fixing be executed for binary variables?</td></tr>
 * <tr><td>propagating/sdpredcost/forintcons = TRUE</td> <td>Should SDP reduced cost fixing be executed for integer and continuous variables?</td></tr>
 * <tr><td>relaxing/SDP/maxpenaltyparam = -1</td> <td>the maximum value of the penalty parameter Gamma used for the penalty formulation if the SDP solver didn't converge; set this to a negative value to compute the parameter depending on the given problem</td></tr>
 * <tr><td>relaxing/SDP/npenaltyincr = -1</td> <td>maximum number of times the penalty parameter will be increased if the penalty formulation failed</td></tr>
 * <tr><td>relaxing/SDP/objlimit = FALSE</td> <td>Should an objective limit be given to the SDP-Solver?</td></tr>
 * <tr><td>relaxing/SDP/resolve = TRUE</td> <td>Should the relaxation be resolved after bound-tightenings were found during propagation (outside of probing)?</td></tr>
 * <tr><td>relaxing/SDP/settingsresetfreq = -1</td> <td>frequency for resetting parameters in SDP solver and trying again with fastest settings [-1: never, 0: only at depth settingsresetofs, n: all nodes with depth a multiple of n]; currently only supported for SDPA</td></tr>
 * <tr><td>relaxing/SDP/settingsresetofs = 0</td> <td>frequency offset for resetting parameters in SDP solver and trying again with fastest settings; currently only supported for SDPA</td></tr>
 * <tr><td>relaxing/SDP/tightenvb = TRUE</td> <td>Should Big-Ms in varbound-like constraints be tightened before giving them to the SDP-solver ?</td></tr>
 * </table>
 */
