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

/**@file   lapack.h
 * @brief  interface methods for eigenvector computation and matrix multiplication using different versions of LAPACK and BLAS
 * @author Tristan Gally
 *
 * This file is used to call the LAPACK routine DSYEVR (double-symmetric-eigenvector computation) and the
 * BLAS routine DGEMV (double-general-matrix-vector multiplication). It is needed because different SDP-solvers
 * need different BLAS/LAPACK-versions with different data types (for example long long int for
 * Openblas/SDPA vs. int for ATLAS/DSDP).
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_LAPACK_H__
#define __SCIP_LAPACK_H__

#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"

#ifdef __cplusplus
extern "C" {
#endif

/** computes the i-th eigenvalue of a symmetric matrix using LAPACK, where 1 is the smallest and n the largest, matrix has to be given with all \f$n^2\f$ entries */
EXTERN
SCIP_RETCODE SCIPlapackComputeIthEigenvalue(
   BMS_BUFMEM*           bufmem,             /**< buffer memory */
   SCIP_Bool             geteigenvectors,    /**< should also the eigenvectors be computed? */
   int                   n,                  /**< size of matrix */
   SCIP_Real*            A,                  /**< matrix for which eigenvalues should be computed */
   int                   i,                  /**< index of eigenvalue to be computed */
   SCIP_Real*            eigenvalue,         /**< pointer to store eigenvalue */
   SCIP_Real*            eigenvector         /**< pointer to store eigenvector */
   );

/** performs matrix-vector-multiplication using BLAS */
EXTERN
SCIP_RETCODE SCIPlapackMatrixVectorMult(
   int                   nrows,              /**< number of rows in matrix */
   int                   ncols,              /**< number of cols in matrix */
   SCIP_Real*            matrix,             /**< the matrix we want to multiply */
   SCIP_Real*            vector,             /**< vector we want to multiply with the matrix */
   SCIP_Real*            result              /**< pointer to store the resulting vector */
   );

#ifdef __cplusplus
}
#endif

#endif
