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

/**@file   cons_sdp.c
 * @brief  Constraint handler for SDP-constraints
 * @author Sonja Mars
 * @author Lars Schewe
 * @author Tristan Gally
 *
 * Constraint handler for semidefinite constraints of the form \f$ \sum_{j=1}^n A_j y_j - A_0 \succeq 0 \f$,
 * where the matrices \f$A_j\f$ and \f$A_0\f$ need to be symmetric. Only the nonzero entries of the matrices
 * are stored.
 */

/* #define SCIP_DEBUG*/
/* #define SCIP_MORE_DEBUG *//* shows all cuts added and prints constraint after parsing */
/* #define PRINT_HUMAN_READABLE *//* change the output of PRINTCONS to a better readable format (dense instead of sparse), WHICH CAN NO LONGER BE PARSED */

#include "cons_sdp.h"

#include <assert.h>                     /*lint !e451*/
#include <string.h>                     /* for NULL, strcmp */
#include <ctype.h>                      /* for isspace */

#include "sdpi/lapack.h"

#include "scipsdp/SdpVarmapper.h"
#include "scipsdp/SdpVarfixer.h"

#include "scip/cons_linear.h"           /* for SCIPcreateConsLinear */
#include "scip/scip.h"                  /* for SCIPallocBufferArray, etc */

#ifdef OMP
#include "omp.h"                        /* for changing the number of threads */
#endif

/* turn off lint warnings for whole file: */
/*lint --e{788,818}*/

#define CONSHDLR_NAME          "SDP"
#define CONSHDLR_DESC          "SDP constraints of the form \\sum_{j} A_j y_j - A_0 psd"
#define CONSHDLR_SEPAPRIORITY  +1000000 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY  -2000000 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY -2000000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ             1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define CONSHDLR_PRESOLTIMING     SCIP_PRESOLTIMING_FAST
#define PARSE_STARTSIZE               1 /**< initial size of the consdata-arrays when parsing a problem */
#define PARSE_SIZEFACTOR             10 /**< size of consdata-arrays is increased by this factor when parsing a problem */
#define DEFAULT_DIAGGEZEROCUTS     TRUE /**< Should linear cuts enforcing the non-negativity of diagonal entries of SDP-matrices be added? */
#define DEFAULT_DIAGZEROIMPLCUTS   TRUE /**< Should linear cuts enforcing the implications of diagonal entries of zero in SDP-matrices be added? */
#ifdef OMP
#define DEFAULT_NTHREADS              1 /**< number of threads used for OpenBLAS */
#endif

/** constraint data for sdp constraints */
struct SCIP_ConsData
{
   int                   nvars;              /**< number of variables in this SDP constraint */
   int                   nnonz;              /**< number of nonzeros in this SDP constraint */
   int                   blocksize;          /**< size of this SDP-block */
   int*                  nvarnonz;           /**< length of the arrays pointed to by col/row/val, number of nonzeros for each variable */
   int**                 col;                /**< pointers to the column indices of the nonzeros for each variable */
   int**                 row;                /**< pointers to the row indices of the nonzeros for each variable */
   SCIP_Real**           val;                /**< pointers to the values of the nonzeros for each variable */
   SCIP_VAR**            vars;               /**< SCIP_VARiables present in this SDP constraint, ordered by their begvar-indices */
   int                   constnnonz;         /**< number of nonzeros in the constant part of this SDP constraint */
   int*                  constcol;           /**< column indices of the constant nonzeros */
   int*                  constrow;           /**< row indices of the constant nonzeros */
   SCIP_Real*            constval;           /**< values of the constant nonzeros */
   SCIP_Real             maxrhsentry;        /**< maximum entry of constant matrix (needed for DIMACS error norm) */
};

/** SDP constraint handler data */
struct SCIP_ConshdlrData
{
   int                   neigveccuts;        /**< this is used to give the eigenvector-cuts distinguishable names */
   SCIP_Bool             diaggezerocuts;     /**< Should linear cuts enforcing the non-negativity of diagonal entries of SDP-matrices be added? */
   int                   ndiaggezerocuts;    /**< this is used to give the diagGEzero-cuts distinguishable names */
   int                   n1x1blocks;         /**< this is used to give the lp constraints resulting from 1x1 sdp-blocks distinguishable names */
   SCIP_Bool             diagzeroimplcuts;   /**< Should linear cuts enforcing the implications of diagonal entries of zero in SDP-matrices be added? */
#ifdef OMP
   int                   nthreads;           /**< number of threads used for OpenBLAS */
#endif
};

#ifndef NDEBUG
/** for given row and column (i,j) computes the position in the lower triangular part, if
 *  these positions are numbered from 0 to n(n+1)/2 - 1, this needs to be called for i >= j
 */
static
int compLowerTriangPos(
   int                   i,                  /**< row index */
   int                   j                   /**< column index */
   )
{
   assert( j >= 0 );
   assert( i >= j );

   return i*(i+1)/2 + j;
}
#else
#define compLowerTriangPos(i, j) (i*(i+1)/2 + j)
#endif

/** takes a 0.5*n*(n+1) array of a symmetric matrix and expands it to an n*n array of the full matrix to input into LAPACK */
static
SCIP_RETCODE expandSymMatrix(
   int                   size,               /**< size of the matrix, named n above */
   SCIP_Real*            symMat,             /**< symmetric matrix indexed via compLowerTriangPos that should be expanded */
   SCIP_Real*            fullMat             /**< pointer to store the n*n matrix, that is the symmetric expansion of symMat */
   )
{
   int i;
   int j;
   int ind = 0;

   assert( size >= 0 );
   assert( symMat != NULL );
   assert( fullMat != NULL );

   /* traverse the lower triangular part in the order of the indices and copy the values to both lower and upper triangular part */
   for (i = 0; i < size; i++)
   {
      for (j = 0; j <= i; j++)
      {
         assert( ind == compLowerTriangPos(i,j) );
         fullMat[i*size + j] = symMat[ind]; /*lint !e679*/
         fullMat[j*size + i] = symMat[ind]; /*lint !e679*/
         ind++;
      }
   }

   return SCIP_OKAY;
}

/** For a given vector \f$ y \f$ computes the (length of y) * (length of y + 1) /2 -long array of the lower-triangular part
 *  of the SDP-Matrix \f$ \sum_{j=1}^m A_j y_j - A_0 \f$ for this SDP block, indexed by compLowerTriangPos.
 */
static
SCIP_RETCODE computeSdpMatrix(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< the constraint for which the Matrix should be assembled */
   SCIP_SOL*             y,                  /**< solution to separate */
   SCIP_Real*            matrix              /**< pointer to store the SDP-Matrix */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real yval;
   int blocksize;
   int nvars;
   int ind;
   int i;

   assert( cons != NULL );
   assert( matrix != NULL );

   consdata = SCIPconsGetData(cons);
   nvars = consdata->nvars;
   blocksize = consdata->blocksize;

   /* initialize the matrix with 0 */
   for (i = 0; i < (blocksize * (blocksize + 1))/2; i++)
      matrix[i] = 0.0;

   /* add the non-constant-part */
   for (i = 0; i < nvars; i++)
   {
      yval = SCIPgetSolVal(scip, y, consdata->vars[i]);
      for (ind = 0; ind < consdata->nvarnonz[i]; ind++)
         matrix[compLowerTriangPos(consdata->row[i][ind], consdata->col[i][ind])] += yval * consdata->val[i][ind];
   }

   /* substract the constant part */
   for (ind = 0; ind < consdata->constnnonz; ind++)
      matrix[compLowerTriangPos(consdata->constrow[ind], consdata->constcol[ind])] -= consdata->constval[ind];

   return SCIP_OKAY;
}

/** For a given variable-index j and a Vector v computes \f$ v^T A_j v \f$. */
static
SCIP_RETCODE multiplyConstraintMatrix(
   SCIP_CONS*            cons,               /**< the SDP constraint that includes the Matrix \f$ A_j \f$ */
   int                   j,                  /**< variable-index of the matrix to multiply with */
   SCIP_Real*            v,                  /**< vector to multiply with */
   SCIP_Real*            vAv                 /**< pointer to store the the resulting scalar \f$ v^T A_j v \f$ */
   )
{
   SCIP_CONSDATA* consdata;
   int i;

   assert( cons != NULL );
   assert( j >= 0 );
   assert( vAv != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   assert( j < consdata->nvars );

   /* initialize the product with 0 */
   *vAv = 0.0;

   for (i = 0; i < consdata->nvarnonz[j]; i++)
   {
      if (consdata->col[j][i] == consdata->row[j][i])
         *vAv += v[consdata->col[j][i]] * consdata->val[j][i] * v[consdata->row[j][i]];
      else
      {
         /* Multiply by 2, because the matrix is symmetric and there is one identical contribution each from lower and upper triangular part. */
         *vAv += 2.0 * v[consdata->col[j][i]] * consdata->val[j][i] * v[consdata->row[j][i]];
      }
   }

   return SCIP_OKAY;
}

/** Set the maximum absolute value of an entry of the constant matrix.
 *  This must be done before presolving, because otherwise this is influenced by variable fixings (which might lead
 *  to solutions being feasible in presolving no longer being feasible afterwards) */
static
SCIP_RETCODE setMaxRhsEntry(
   SCIP_CONS*            cons                /**< the SDP constraint that includes the Matrix \f$ A_j \f$ */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real max;
   int i;

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   /* initialize max with zero (this is used if there is no constant-matrix) */
   max = 0.0;

   /* iterate over the entries of the constant matrix, updating max if a higher absolute value is found */
   for (i = 0; i < consdata->constnnonz; i++)
   {
      if ( REALABS(consdata->constval[i]) > max )
         max = REALABS(consdata->constval[i]);
   }

   consdata->maxrhsentry = max;

   return SCIP_OKAY;
}


/** separate current solution with a cut using the eigenvectors and -values of the solution matrix
 *
 *  This function computes the eigenvectors of the matrix, takes the one corresponding to the smallest eigenvalue and
 *  multiplies the matrix with it such that \f$ coeff[i] = x^TA_ix , lhs = x^TA_0x \f$.
 */
static
SCIP_RETCODE cutUsingEigenvector(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< the constraint to compute the cut for */
   SCIP_SOL*             sol,                /**< solution to separate */
   SCIP_Real*            coeff,              /**< pointer to store the coefficients of the computed cut */
   SCIP_Real*            lhs                 /**< pointer to store the lhs of the computed cut */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real* eigenvalues = NULL;
   SCIP_Real* matrix = NULL;
   SCIP_Real* fullmatrix = NULL;
   SCIP_Real* fullconstmatrix = NULL;
   SCIP_Real* eigenvector = NULL;
   SCIP_Real* output_vector = NULL;
   int blocksize;
   int j;

   assert( cons != NULL );
   assert( lhs != NULL );

   *lhs = 0.0;

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   blocksize = consdata->blocksize;

   SCIP_CALL( SCIPallocBufferArray(scip, &eigenvalues, blocksize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix, (blocksize * (blocksize+1))/2 ) ); /*lint !e647*/
   SCIP_CALL( SCIPallocBufferArray(scip, &fullmatrix, blocksize * blocksize ) ); /*lint !e647*/
   SCIP_CALL( SCIPallocBufferArray(scip, &fullconstmatrix, blocksize * blocksize) ); /*lint !e647*/
   SCIP_CALL( SCIPallocBufferArray(scip, &eigenvector, blocksize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &output_vector, blocksize) );

   /* compute the matrix \f$ \sum_j A_j y_j - A_0 \f$ */
   SCIP_CALL( computeSdpMatrix(scip, cons, sol, matrix) );

   /* expand it because LAPACK wants the full matrix instead of the lower triangular part */
   SCIP_CALL( expandSymMatrix(blocksize, matrix, fullmatrix) );

   SCIP_CALL( SCIPlapackComputeIthEigenvalue(SCIPbuffer(scip), TRUE, blocksize, fullmatrix, 1, eigenvalues, eigenvector) );

   /* get full constant matrix */
   SCIP_CALL( SCIPconsSdpGetFullConstMatrix(scip, cons, fullconstmatrix) );

   /* multiply eigenvector with constant matrix to get lhs (after multiplying again with eigenvector from the left) */
   SCIP_CALL( SCIPlapackMatrixVectorMult(blocksize, blocksize, fullconstmatrix, eigenvector, output_vector) );

   for (j = 0; j < blocksize; ++j)
      *lhs += eigenvector[j] * output_vector[j];

   /* compute \f$ v^T A_j v \f$ for eigenvector v and each matrix \f$ A_j \f$ to get the coefficients of the LP cut */
   for (j = 0; j < consdata->nvars; ++j)
   {
      SCIP_CALL( multiplyConstraintMatrix(cons, j, eigenvector, &coeff[j]) );
   }

   SCIPfreeBufferArray(scip, &output_vector);
   SCIPfreeBufferArray(scip, &eigenvector);
   SCIPfreeBufferArray(scip, &fullconstmatrix);
   SCIPfreeBufferArray(scip, &fullmatrix);
   SCIPfreeBufferArray(scip, &matrix);
   SCIPfreeBufferArray(scip, &eigenvalues);

   return SCIP_OKAY;
}

/** checks feasibility for a single SDP constraint */
SCIP_RETCODE SCIPconsSdpCheckSdpCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< the constraint the solution should be checked for */
   SCIP_SOL*             sol,                /**< the solution to check feasibility for */
   SCIP_Bool             checkintegrality,   /**< has integrality to be checked? */
   SCIP_Bool             checklprows,        /**< have current LP rows to be checked? */
   SCIP_Bool             printreason,        /**< should the reason for the violation be printed? */
   SCIP_RESULT*          result              /**< pointer to store the result of the feasibility checking call */
   )
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int blocksize;
   SCIP_Real check_value;
   SCIP_Real eigenvalue;
   SCIP_Real* matrix = NULL;
   SCIP_Real* fullmatrix = NULL;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( result != NULL );
   assert( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), "SDP") == 0);

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   blocksize = consdata->blocksize;

   SCIP_CALL( SCIPallocBufferArray(scip, &matrix, (blocksize * (blocksize+1)) / 2) ); /*lint !e647*/
   SCIP_CALL( SCIPallocBufferArray(scip, &fullmatrix, blocksize * blocksize) ); /*lint !e647*/
   SCIP_CALL( computeSdpMatrix(scip, cons, sol, matrix) );
   SCIP_CALL( expandSymMatrix(blocksize, matrix, fullmatrix) );

   SCIP_CALL( SCIPlapackComputeIthEigenvalue(SCIPbuffer(scip), FALSE, blocksize, fullmatrix, 1, &eigenvalue, NULL) );

   /* This enables checking the second DIMACS Error Norm: err=max{0, -lambda_min(x)/(1+maximumentry of rhs)}, in that case it also needs
    * to be changed in the sdpi (and implemented there first), when checking feasibility of problems where all variables are fixed */
#ifdef DIMACS
   check_value = (-eigenvalue) / (1.0 + consdata->maxrhsentry);
#else
   check_value = -eigenvalue;
#endif

   if ( SCIPisFeasLE(scip, check_value, 0.0) )
      *result = SCIP_FEASIBLE;
   else
   {
      *result = SCIP_INFEASIBLE;
#ifdef SCIP_DEBUG
      if ( printreason )
      {
         SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
         SCIP_CALL( SCIPprintSol(scip, sol, NULL, FALSE) );
         SCIPinfoMessage(scip, NULL, "non psd matrix (eigenvalue %f, dimacs error norm = %f).\n", eigenvalue, check_value);
      }
#endif
   }

   SCIPfreeBufferArray(scip, &fullmatrix);
   SCIPfreeBufferArray(scip, &matrix);

   return SCIP_OKAY;
}

/** separates the current solution */
static
SCIP_RETCODE separateSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< the constraint handler itself */
   SCIP_CONS*            cons,               /**< constraint to process */
   SCIP_SOL*             sol,                /**< primal solution that should be separated */
   SCIP_RESULT*          result              /**< pointer to store the result of the separation call */
   )
{
   char cutname[SCIP_MAXSTRLEN];
   SCIP_CONSDATA* consdata;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Real lhs = 0.0;
   SCIP_Real* coeff = NULL;
   SCIP_COL** cols;
   SCIP_Real* vals;
   SCIP_ROW* row;
   int nvars;
   int j;
   int len;
#ifndef NDEBUG
   int snprintfreturn; /* this is used to assert that the SCIP string concatenation works */
#endif

   assert( cons != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   nvars = consdata->nvars;
   SCIP_CALL( SCIPallocBufferArray(scip, &coeff, nvars ) );

   SCIP_CALL( cutUsingEigenvector(scip, cons, sol, coeff, &lhs) );

   SCIP_CALL( SCIPallocBufferArray(scip, &cols, nvars ) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vals, nvars ) );

   len = 0;
   for (j = 0; j < nvars; ++j)
   {
      if ( SCIPisZero(scip, coeff[j]) )
         continue;

      cols[len] = SCIPvarGetCol(SCIPvarGetProbvar(consdata->vars[j]));
      vals[len] = coeff[j];
      ++len;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

#ifndef NDEBUG
   snprintfreturn = SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "sepa_eig_sdp_%d", ++(conshdlrdata->neigveccuts));
   assert( snprintfreturn < SCIP_MAXSTRLEN ); /* check whether name fit into string */
#else
   (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "sepa_eig_sdp_%d", ++(conshdlrdata->neigveccuts));
#endif
   SCIP_CALL( SCIPcreateRowCons(scip, &row, conshdlr, cutname , len, cols, vals, lhs, SCIPinfinity(scip), FALSE, FALSE, TRUE) );

   if ( SCIPisCutEfficacious(scip, sol, row) )
   {
      SCIP_Bool infeasible;
#ifdef SCIP_MORE_DEBUG
      SCIPinfoMessage(scip, NULL, "Added cut %s: ", cutname);
      SCIPinfoMessage(scip, NULL, "%f <= ", lhs);
      for (j = 0; j < len; j++)
         SCIPinfoMessage(scip, NULL, "+ (%f)*%s", vals[j], SCIPvarGetName(SCIPcolGetVar(cols[j])));
      SCIPinfoMessage(scip, NULL, "\n");
#endif

      SCIP_CALL( SCIPaddCut(scip, sol, row, FALSE, &infeasible) );
      if ( infeasible )
         *result = SCIP_CUTOFF;
      else
         *result = SCIP_SEPARATED;
      SCIP_CALL( SCIPresetConsAge(scip, cons) );
   }

   SCIP_CALL( SCIPreleaseRow(scip, &row) );

   SCIPfreeBufferArray(scip, &vals);
   SCIPfreeBufferArray(scip, &cols);
   SCIPfreeBufferArray(scip, &coeff);

   return SCIP_OKAY;
}

/** approximates the sdpcone using the fact that every diagonal entry must be non-negative, so it adds the LP-cut
 *  \f$ \sum_{j = 1}^m (A_j)_{kk} y_j - (A_0)_{kk} \geq 0 \quad \forall k \leq n \f$
 */
static
SCIP_RETCODE diagGEzero(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< array of constraints to add cuts for */
   int                   nconss,             /**< number of constraints to add cuts for */
   int*                  naddconss           /**< pointer to store how many constraints were added */
   )
{
   char cutname[SCIP_MAXSTRLEN];
   SCIP_CONSDATA* consdata;
   SCIP_Real* matrix;
   SCIP_Real* cons_array;
   SCIP_Real* lhs_array;
   SCIP_Real rhs;
   int blocksize;
   int nvars;
   int i;
   int j;
   int k;
   int c;
#ifndef NDEBUG
   int snprintfreturn; /* used to check if sdnprintf worked */
#endif

   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONSHDLR* conshdlr;

      conshdlr = SCIPconsGetHdlr(conss[c]);
      assert( conshdlr != NULL );
#ifndef NDEBUG
      assert( strcmp(SCIPconshdlrGetName(conshdlr), "SDP") == 0);
#endif

      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );

      blocksize = consdata->blocksize;
      nvars = consdata->nvars;
      rhs = SCIPinfinity(scip);

      SCIP_CALL( SCIPallocBufferArray(scip, &matrix, (blocksize * (blocksize + 1)) / 2) ); /*lint !e647*/
      SCIP_CALL( SCIPconsSdpGetLowerTriangConstMatrix(scip, conss[c], matrix) );

      SCIP_CALL( SCIPallocBufferArray(scip, &cons_array, blocksize * nvars) ); /*lint !e647*/
      SCIP_CALL( SCIPallocBufferArray(scip, &lhs_array, blocksize) );

      /* the lhs is the (k,k)-entry of the constant matrix */
      for (k = 0; k < blocksize; ++k)
         lhs_array[k] = matrix[compLowerTriangPos(k,k)];

      /* get the (k,k)-entry of every matrix A_j */
      for (j = 0; j < nvars; ++j)
      {
         for (k = 0; k < blocksize; ++k)
         {
            /* initialize these with 0 */
            cons_array[k * nvars + j] = 0.0; /*lint !e679*/
         }

         /* go through the nonzeros of A_j and look for diagonal entries */
         for (i = 0; i < consdata->nvarnonz[j]; i++)
         {
            if ( consdata->col[j][i] == consdata->row[j][i] )
               cons_array[consdata->col[j][i] * nvars + j] = consdata->val[j][i]; /*lint !e679*/
         }
      }

      /* add the LP-cuts to SCIP */
      for (k = 0; k < blocksize; ++k)
      {
         SCIP_CONS* cons;
         SCIP_CONSHDLRDATA* conshdlrdata;

         conshdlrdata = SCIPconshdlrGetData(conshdlr);
#ifndef NDEBUG
         snprintfreturn = SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "diag_ge_zero_%d", ++(conshdlrdata->ndiaggezerocuts));
         assert( snprintfreturn < SCIP_MAXSTRLEN );
#else
         (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "diag_ge_zero_%d", ++(conshdlrdata->ndiaggezerocuts));
#endif

#ifdef SCIP_MORE_DEBUG
         SCIPinfoMessage(scip, NULL, "Added lp-constraint %s: ", cutname);
         SCIPinfoMessage(scip, NULL, "%f <= ", lhs_array[k]);
         for (i = 0; i < consdata->nvars; i++)
         {
            if ( ! SCIPisZero(scip, cons_array[k * consdata->nvars + i]) )
               SCIPinfoMessage(scip, NULL, "+ (%f)*%s", cons_array[k * consdata->nvars + i], SCIPvarGetName(consdata->vars[i]));
         }
         SCIPinfoMessage(scip, NULL, "\n");
#endif

         SCIP_CALL( SCIPcreateConsLinear(scip, &cons, cutname, consdata->nvars, consdata->vars, cons_array + k * consdata->nvars, lhs_array[k], rhs,
               TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) ); /*lint !e679*/

         SCIP_CALL( SCIPaddCons(scip, cons) );
         SCIP_CALL( SCIPreleaseCons(scip, &cons) );
         ++(*naddconss);
      }

      SCIPfreeBufferArray(scip, &lhs_array);
      SCIPfreeBufferArray(scip, &cons_array);
      SCIPfreeBufferArray(scip, &matrix);
   }

   return SCIP_OKAY;
}

/** Presolve-routine that enforces implications of diagonal entries of zero in SDP-matrices, namely that if \f$X_{ij} > 0\f$,
 *  then also \f$X_{ii} > 0\f$ and \f$ X_{jj} > 0\f$.
 *
 *  More precisely, if \f$ (A_0)_{k\ell} \neq 0\f$, \f$ (A_0)_{kk} = 0\f$, \f$ (A_i)_{k\ell} = 0\f$ for all \f$ i \leq m\f$,
 *  \f$ (A_0)_{kk} = 0\f$ for all continuous variables and \f$ \ell_i \geq 0\f$ for all integer variables, we add the cut
 *  \f$ \sum_{\substack{i \in \mathcal{I}:\\ (A_i)_{kk} > 0}} y_i \geq 1.\f$
 */
static
SCIP_RETCODE diagZeroImpl(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< array of constraints */
   int                   nconss,             /**< number of constraints */
   int*                  naddconss           /**< pointer to store how many constraints were added */
   )
{
   char cutname[SCIP_MAXSTRLEN];
   /* if entry k is >= 0, gives the number of non-diagonal nonzero-entries in row k of the constant matrix (and the number
    * of entries in constnonzeroentries), if -1 C_kk =!= 0 (which means we didnot allocate memory for diagvars), if -2
    * either A_jk =!= 0 for all j with C_jk =!= 0 or A_kk =!= 0 for some continuous variable (which means we did allocate
    * memory for diagvars but cannot use the cut */
   int* nconstnonzeroentries;
   int** constnonzeroentries;
   SCIP_CONSDATA* consdata;
   SCIP_CONS* cons;
   SCIP_VAR** vars;
   SCIP_Real* vals;
   int** diagvars;
   int* ndiagvars;
   int blocksize;
   int i;
   int j;
   int nvars;
   int v;
   int k;
   int l;
   SCIP_Bool anycutvalid;
#ifndef NDEBUG
   int snprintfreturn;
#endif

   assert( scip != NULL );
   assert( conss != NULL );
   assert( nconss >= 0 );
   assert( naddconss != NULL );

   for (i = 0; i < nconss; ++i)
   {
      assert( conss[i] != NULL );
      assert( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[i])), "SDP") == 0 );

      consdata = SCIPconsGetData(conss[i]);
      assert( consdata != NULL );

      blocksize = consdata->blocksize;
      nvars = consdata->nvars;
      SCIP_CALL( SCIPallocBufferArray(scip, &nconstnonzeroentries, blocksize) );
      SCIP_CALL( SCIPallocBufferArray(scip, &constnonzeroentries, blocksize) );
      for (j = 0; j < blocksize; j++)
      {
         nconstnonzeroentries[j] = 0;
         SCIP_CALL( SCIPallocBufferArray(scip, &constnonzeroentries[j], blocksize) );
      }

      /* iterate over all nonzeros of the constant matrix and check which diagonal and non-diagonal entries are nonzero */
      for (j = 0; j < consdata->constnnonz; j++)
      {
         /* if it is a nondiagonal-entry we add this row/column to the constnonzeroentries entries unless we already found a
          * diagonal entry for this row/column */
         if ( (consdata->constcol[j] != consdata->constrow[j]) )
         {
            assert( ! SCIPisZero(scip, consdata->constval[j]) );
            if ( nconstnonzeroentries[consdata->constcol[j]] >= 0 )
            {
               constnonzeroentries[consdata->constcol[j]][nconstnonzeroentries[consdata->constcol[j]]] = consdata->constrow[j];
               nconstnonzeroentries[consdata->constcol[j]]++;
            }
            if ( nconstnonzeroentries[consdata->constrow[j]] >= 0 )
            {
               constnonzeroentries[consdata->constrow[j]][nconstnonzeroentries[consdata->constrow[j]]] = consdata->constcol[j];
               nconstnonzeroentries[consdata->constrow[j]]++;
            }
         }
         /* if we find a diagonal entry in the constant matrix, we remember that we cannot add a cut for this index */
         else
         {
            assert( ! SCIPisZero(scip, consdata->constval[j]) );
            nconstnonzeroentries[consdata->constcol[j]] = -1;
         }
      }

      /* diagvars[j] is an array with all variables with a diagonal entry (j,j) in the corresponding matrix, if nconstnonzeroentries[j] =!= -1 or NULL otherwise
       * the outer array goes over all rows to ease the access, but only for those that are really needed memory will be allocated */
      SCIP_CALL( SCIPallocBufferArray(scip, &diagvars, blocksize) );
      SCIP_CALL( SCIPallocBufferArray(scip, &ndiagvars, blocksize) );
      anycutvalid = FALSE;
      for (j = 0; j < blocksize; ++j)
      {
         ndiagvars[j] = 0;
         if ( nconstnonzeroentries[j] > 0 )
         {
            SCIP_CALL( SCIPallocBufferArray(scip, &(diagvars[j]), nvars) );
            anycutvalid = TRUE;
         }
      }

      /* if no cuts are valid for this block, we free all memory and continue with the next block */
      if ( ! anycutvalid )
      {
         SCIPfreeBufferArray(scip, &ndiagvars);
         SCIPfreeBufferArray(scip, &diagvars);
         for (j = blocksize - 1; j >= 0; j--)
         {
            SCIPfreeBufferArray(scip, &constnonzeroentries[j]);
         }
         SCIPfreeBufferArray(scip, &constnonzeroentries);
         SCIPfreeBufferArray(scip, &nconstnonzeroentries);
         continue;
      }

      /* find all variables with corresponding diagonal entries for a row with nonzero non-diagonal constant entry, also check for entries
       * that prevent the cut from being valid */
      for (v = 0; v < nvars; v++)
      {
         for (j = 0; j < consdata->nvarnonz[v]; j++)
         {
            /* if it is a diagonal entry for an index that might have a valid cut, we add the variable to the corresponding array if it
             * is an integer variable and mark the cut invalid otherwise */
            if ( (consdata->col[v][j] == consdata->row[v][j]) && (nconstnonzeroentries[consdata->col[v][j]] > 0) )
            {
               if ( SCIPvarIsIntegral(consdata->vars[v]) )
               {
                  assert( ! SCIPisEQ(scip, consdata->val[v][j], 0.0) );
                  diagvars[consdata->col[v][j]][ndiagvars[consdata->col[v][j]]] = v;
                  ndiagvars[consdata->col[v][j]]++;
               }
               else
               {
                  nconstnonzeroentries[consdata->col[v][j]] = -2;
               }
            }
            /* If it is a non-diagonal entry, we can no longer use this entry for a cut. If the last entry is removed for a column/row,
             * mark this column/row invalid (but we still have to free memory later, so we have to set it to -2 instead of 0) */
            else if ( consdata->col[v][j] != consdata->row[v][j] )
            {
               if ( nconstnonzeroentries[consdata->col[v][j]] > 0 )
               {
                  /* search for the corresponding row-entry in constnonzeroentries */
                  for (k = 0; k < nconstnonzeroentries[consdata->col[v][j]]; k++)
                  {
                     if ( constnonzeroentries[consdata->col[v][j]][k] == consdata->row[v][j] )
                     {
                        /* if there are remaining entries, we shift them back */
                        if ( nconstnonzeroentries[consdata->col[v][j]] > k + 1 )
                        {
                           for (l = k + 1; l < nconstnonzeroentries[consdata->col[v][j]]; l++)
                              constnonzeroentries[consdata->col[v][j]][l - 1] = constnonzeroentries[consdata->col[v][j]][l];
                        }
                        nconstnonzeroentries[consdata->col[v][j]]--;
                        /* if this was the last entry for this index, we mark it invalid */
                        if ( nconstnonzeroentries[consdata->col[v][j]] == 0 )
                           nconstnonzeroentries[consdata->col[v][j]] = -2;
                        break; /* we should not have another entry for this combination of row and column */
                     }
                  }
               }
               /* do the same for the row */
               if ( nconstnonzeroentries[consdata->row[v][j]] > 0 )
               {
                  /* search for the corresponding row-entry in constnonzeroentries */
                  for (k = 0; k < nconstnonzeroentries[consdata->row[v][j]]; k++)
                  {
                     if ( constnonzeroentries[consdata->row[v][j]][k] == consdata->col[v][j] )
                     {
                        /* if there are remaining entries, we shift them back */
                        if ( nconstnonzeroentries[consdata->row[v][j]] > k + 1 )
                        {
                           for (l = k + 1; l < nconstnonzeroentries[consdata->row[v][j]]; l++)
                              constnonzeroentries[consdata->row[v][j]][l - 1] = constnonzeroentries[consdata->row[v][j]][l];
                        }
                        nconstnonzeroentries[consdata->row[v][j]]--;
                        /* if this was the last entry for this index, we mark it invalid */
                        if ( nconstnonzeroentries[consdata->row[v][j]] == 0 )
                           nconstnonzeroentries[consdata->row[v][j]] = -2;
                        break; /* we should not have another entry for this combination of row and column */
                     }
                  }
               }
            }
         }
      }

      for (j = 0; j < blocksize; ++j)
      {
         if ( nconstnonzeroentries[j] > 0 )
         {
            SCIP_CALL( SCIPallocBufferArray(scip, &vals, ndiagvars[j]) );
            SCIP_CALL( SCIPallocBufferArray(scip, &vars, ndiagvars[j]) );

            /* get the corresponding SCIP variables and set all coefficients to 1 */
            for (v = 0; v < ndiagvars[j]; ++v)
            {
               vars[v] = consdata->vars[diagvars[j][v]];
               vals[v] = 1.0;
            }
#ifndef NDEBUG
            snprintfreturn = SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "diag_zero_impl_block_%d_row_%d", i, j);
            assert( snprintfreturn < SCIP_MAXSTRLEN );  /* check whether name fits into string */
#else
            (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "diag_zero_impl_block_%d_row_%d", i, j);
#endif

#ifdef SCIP_MORE_DEBUG
            SCIPinfoMessage(scip, NULL, "Added lp-constraint %s: ", cutname);
            SCIPinfoMessage(scip, NULL, "1 <=");
            for (l = 0; l < ndiagvars[j]; l++)
               SCIPinfoMessage(scip, NULL, " + (%f)*%s", vals[l], SCIPvarGetName(vars[l]));
            SCIPinfoMessage(scip, NULL, "\n");
#endif

            /* add the linear constraint sum_j 1.0 * diagvars[j] >= 1.0 */
            SCIP_CALL(SCIPcreateConsLinear(scip, &cons, cutname , ndiagvars[j], vars, vals, 1.0, SCIPinfinity(scip), TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE));
            SCIP_CALL(SCIPaddCons(scip, cons));
            SCIP_CALL(SCIPreleaseCons(scip, &cons));
            (*naddconss)++;

            SCIPfreeBufferArray(scip, &vars);
            SCIPfreeBufferArray(scip, &vals);
         }
         if ( nconstnonzeroentries[j] == -2 || nconstnonzeroentries[j] > 0 )
         {
            SCIPfreeBufferArray(scip, &diagvars[j]);
         }
      }

      SCIPfreeBufferArray(scip, &ndiagvars);
      SCIPfreeBufferArray(scip, &diagvars);
      for (j = blocksize - 1; j >= 0; j--)
      {
         SCIPfreeBufferArray(scip, &constnonzeroentries[j]);
      }
      SCIPfreeBufferArray(scip, &constnonzeroentries);
      SCIPfreeBufferArray(scip, &nconstnonzeroentries);
   }

   return SCIP_OKAY;
}

/** detects if there are blocks with size one and transforms them to lp-rows */
static
SCIP_RETCODE move_1x1_blocks_to_lp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< array of constraints to check */
   int                   nconss,             /**< number of constraints to check */
   int*                  naddconss,          /**< pointer to store how many constraints were added */
   int*                  ndelconss,          /**< pointer to store how many constraints were deleted */
   SCIP_RESULT*          result              /**< pointer to store if this routine was successfull or if it detected infeasibility */
   )
{
   char cutname[SCIP_MAXSTRLEN];
   SCIP_CONSHDLR* hdlr;
   SCIP_CONSDATA* consdata;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS* cons;
   SCIP_VAR** vars;
   SCIP_Real* coeffs;
   SCIP_Real rhs;
   int nnonz;
   int nvars;
   int i;
   int j;
   int count;
   int var;
#ifndef NDEBUG
   int snprintfreturn; /* used to assert the return code of snprintf */
#endif

   assert( result != NULL );
   *result = SCIP_SUCCESS;

   for (i = 0; i < nconss; ++i)
   {
      hdlr = SCIPconsGetHdlr(conss[i]);
      assert(hdlr != NULL);

#ifndef NDEBUG
      assert( strcmp(SCIPconshdlrGetName(hdlr), "SDP") == 0);
#endif

      consdata = SCIPconsGetData(conss[i]);
      assert( consdata != NULL );

      /* if there is a 1x1 SDP-Block */
      if ( consdata->blocksize == 1 )
      {
         nvars = consdata->nvars;
         nnonz = consdata->nnonz;
         SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );
         SCIP_CALL( SCIPallocBufferArray(scip, &coeffs, nnonz) );

         /* get all lhs-entries */
         count = 0;

         for (var = 0; var < nvars; var++)
         {
            assert( consdata->nvarnonz[var] <= 1 ); /* in a 1x1 block there may be at most one entry per variable */

            for (j = 0; j < consdata->nvarnonz[var]; j++)
            {
               assert( consdata->col[var][j] == 0 && consdata->row[var][j] == 0 ); /* if the block is size one, all entries should have row and col equal to 0 */
               coeffs[count] = consdata->val[var][j];
               vars[count] = consdata->vars[var];
               count++;
            }
         }

         /* get rhs */
         assert( consdata->constnnonz <= 1 ); /* the 1x1 constant matrix may only have one entry */

         rhs = (consdata->constnnonz == 1) ? consdata->constval[0] : 0.0; /* if this one entry is not 0, than this is the rhs, otherwise it's 0 */

         /* if there is more than one left-hand-side-entry, add a linear constraint, otherwise update the variable bound */
         if ( count > 1 )
         {
            /* add new linear cons */
            conshdlrdata = SCIPconshdlrGetData(hdlr);
#ifndef NDEBUG
            snprintfreturn = SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "1x1block_%d", ++(conshdlrdata->n1x1blocks));
            assert( snprintfreturn < SCIP_MAXSTRLEN ); /* check whether name fits into string */
#else
            (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "1x1block_%d", ++(conshdlrdata->n1x1blocks));
#endif

#ifdef SCIP_MORE_DEBUG
            SCIPinfoMessage(scip, NULL, "Added lp-constraint %s: ", cutname);
            for (i = 0; i < consdata->nvars; i++)
               SCIPinfoMessage(scip, NULL, "+ (%f)*%s", coeffs[i], SCIPvarGetName(vars[i]));
            SCIPinfoMessage(scip, NULL, " <= %f \n", rhs);
#endif

            SCIP_CALL( SCIPcreateConsLinear(scip, &cons, cutname, consdata->nvars, vars, coeffs, rhs, SCIPinfinity(scip),
                  TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );

            SCIP_CALL( SCIPaddCons(scip, cons) );
            SCIP_CALL( SCIPreleaseCons(scip, &cons) );

            (*naddconss)++;
         }
         else
         {
            /* we compare this new variable with the current (local) bounds, we don't do an epsilon check here, because 10^(-21)*x >= 10^(-19) for
             * epsilon = 10^(-20) would be infeasible then, even though it only says x >= 100 */
            if ( coeffs[0] > 0.0 )
            {
               /* this gives a lower bound, if it is bigger than the current one, we need to update it */
               if ( SCIPisFeasGT(scip, rhs / coeffs[0], SCIPvarGetLbLocal(vars[0])) )
               {
                  /* check if the changed bound renders the problem infeasible */
                  if( SCIPisFeasGT(scip, rhs / coeffs[0], SCIPvarGetUbLocal(vars[0])) )
                  {
                     SCIPdebugMessage("Problem detected to be infeasible during presolving, 1x1-SDP-constraint %s caused change"
                           "of lower bound for variable %s from %f to %f, which is bigger than upper bound of %f\n",
                           SCIPconsGetName(conss[i]), SCIPvarGetName(vars[0]), SCIPvarGetLbLocal(vars[0]), rhs / coeffs[0],
                           SCIPvarGetUbLocal(vars[0]));

                     *result = SCIP_CUTOFF;

                     /* delete old 1x1 sdpcone */
                     SCIP_CALL( SCIPdelCons(scip, conss[i]) );
                     (*ndelconss)++;

                     SCIPfreeBufferArray(scip, &coeffs);
                     SCIPfreeBufferArray(scip, &vars);

                     return SCIP_OKAY; /* the node is infeasible, we don't care for the other constraints */
                  }

                  SCIPdebugMessage("Changing lower bound of variable %s from %f to %f because of 1x1-SDP-constraint %s!\n",
                        SCIPvarGetName(vars[0]), SCIPvarGetLbLocal(vars[0]), rhs / coeffs[0], SCIPconsGetName(conss[i]));
                  SCIP_CALL( SCIPchgVarLb(scip, vars[0], rhs / coeffs[0]) );
               }
               else
               {
                  SCIPdebugMessage("Deleting 1x1-SDP-constraint %s, new lower bound %f for variable %s no improvement over old bound %f!\n",
                        SCIPconsGetName(conss[i]), rhs / coeffs[0], SCIPvarGetName(vars[0]), SCIPvarGetLbLocal(vars[0]));
               }
            }
            else if ( coeffs[0] < 0.0 )
            {
               /* this gives an upper bound, if it is lower than the current one, we need to update it */
               if (SCIPisFeasLT(scip, rhs / coeffs[0], SCIPvarGetUbLocal(vars[0])))
               {
                  /* check if the changed bound renders the problem infeasible */
                  if( SCIPisFeasLT(scip, rhs / coeffs[0], SCIPvarGetLbLocal(vars[0])) )
                  {
                     SCIPdebugMessage("Problem detected to be infeasible during presolving, 1x1-SDP-constraint %s caused change"
                           "of upper bound for variable %s from %f to %f, which is less than lower bound of %f\n",
                           SCIPconsGetName(conss[i]), SCIPvarGetName(vars[0]), SCIPvarGetUbLocal(vars[0]), rhs / coeffs[0],
                           SCIPvarGetLbLocal(vars[0]));

                     *result = SCIP_CUTOFF;

                     /* delete old 1x1 sdpcone */
                     SCIP_CALL( SCIPdelCons(scip, conss[i]) );
                     (*ndelconss)++;

                     SCIPfreeBufferArray(scip, &coeffs);
                     SCIPfreeBufferArray(scip, &vars);

                     return SCIP_OKAY; /* the node is infeasible, we don't care for the other constraints */
                  }

                  SCIPdebugMessage("Changing upper bound of variable %s from %f to %f because of 1x1-SDP-constraint %s!\n",
                     SCIPvarGetName(vars[0]), SCIPvarGetUbLocal(vars[0]), -rhs / coeffs[0], SCIPconsGetName(conss[i]));
                  SCIP_CALL( SCIPchgVarUb(scip, vars[0], rhs / coeffs[0]) );
               }
               else
               {
                  SCIPdebugMessage("Deleting 1x1-SDP-constraint %s, new upper bound %f for variable %s no improvement over old bound %f!\n",
                     SCIPconsGetName(conss[i]), rhs / coeffs[0], SCIPvarGetName(vars[0]), SCIPvarGetUbLocal(vars[0]));
               }
            }
            else
            {
               SCIPdebugMessage("Detected 1x1 SDP-block without any nonzero coefficients \n");
               if ( SCIPisFeasGT(scip, rhs, 0.0) )
               {
                  SCIPdebugMessage("Detected infeasibility in 1x1 SDP-block without any nonzero coefficients but with strictly positive rhs\n");
                  *result = SCIP_CUTOFF;

                  /* delete old 1x1 sdpcone */
                  SCIP_CALL( SCIPdelCons(scip, conss[i]) );
                  (*ndelconss)++;

                  SCIPfreeBufferArray(scip, &coeffs);
                  SCIPfreeBufferArray(scip, &vars);

                  return SCIP_OKAY; /* the node is infeasible, we don't care for the other constraints */
               }
            }
         }

         /* delete old 1x1 sdpcone */
         SCIP_CALL( SCIPdelCons(scip, conss[i]) );
         (*ndelconss)++;

         SCIPfreeBufferArray(scip, &coeffs);
         SCIPfreeBufferArray(scip, &vars);
      }
   }
   return SCIP_OKAY;
}

/** local function to perform (parts of) multiaggregation of a single variable within fixAndAggrVars */
static
SCIP_RETCODE multiaggrVar(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            cons,               /**< constraint to multiaggregate for */
   int*                  v,                  /**< position of the variable that gets (multi-)aggregated */
   SCIP_VAR**            aggrvars,           /**< variables this has to be (multi-)aggregated to */
   SCIP_Real*            scalars,            /**< scalar parts to multiply with for each variable this is aggregated to */
   int                   naggrvars,          /**< number of variables this is (multi-)aggregated to */
   SCIP_Real             constant,           /**< the constant part for the (multi-)aggregation */
   int*                  savedcol,           /**< array of columns for nonzeros that need to be added to the constant part */
   int*                  savedrow,           /**< array of rows for nonzeros that need to be added to the constant part */
   SCIP_Real*            savedval,           /**< array of values for nonzeros that need to be added to the constant part */
   int*                  nfixednonz,         /**< length of the arrays of saved nonzeros for the constant part */
   int*                  vararraylength      /**< length of the variable array */
   )
{
   int i;
   SCIP_CONSDATA* consdata;
   int startind;
   int aggrind;
   int aggrtargetlength;
   int globalnvars;
   int aggrconsind;
   SCIP_Real epsilon;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( scalars != NULL );
   assert( naggrvars > 0 );
   assert( savedcol != NULL );
   assert( savedrow != NULL );
   assert( savedval != NULL );
   assert( nfixednonz != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   SCIP_CALL( SCIPgetRealParam(scip, "numerics/epsilon", &epsilon) );

   /* save the current nfixednonz-index, all entries starting from here will need to be added to the variables this is aggregated to */
   startind = *nfixednonz;

   if ( SCIPisEQ(scip, constant, 0.0) )
   {
      /* if there is no constant part, we just save the nonzeros to add them to the variables they are aggregated to, we do this to be able to remove
       * the nonzero-arrays for this variable to be able to fill it with a newly inserted variable, as copying all variables, if we created an empty
       * gap in the variable arrays, will be more time consuming then copying all variables (as we will usually have more variables than nonzeros per
       * variable */
      for (i = 0; i < consdata->nvarnonz[*v]; i++)
      {
         savedcol[*nfixednonz] = consdata->col[*v][i];
         savedrow[*nfixednonz] = consdata->row[*v][i];
         savedval[*nfixednonz] = consdata->val[*v][i];
         (*nfixednonz)++;
      }
   }
   else
   {
      for (i = 0; i < consdata->nvarnonz[*v]; i++)
      {
         savedcol[*nfixednonz] = consdata->col[*v][i];
         savedrow[*nfixednonz] = consdata->row[*v][i];
         savedval[*nfixednonz] = consdata->val[*v][i] * constant;
         (*nfixednonz)++;
      }
   }
   assert( *nfixednonz - startind == consdata->nvarnonz[*v] );

   /* sort them by nondecreasing row and then col to make the search for already existing entries easier (this is done here, because it
    * only needs to be done once and not for each variable this is multiaggregated to), we add startind to the pointers to only start where we started
    * inserting, the number of elements added to the saved arrays for this variable is nfixednonz - startind */
   SCIPsdpVarfixerSortRowCol(savedrow + startind, savedcol + startind, savedval + startind, *nfixednonz - startind);

   /* free the memory for the entries of the aggregated variable */
   SCIPfreeBlockMemoryArray(scip, &(consdata->val[*v]), consdata->nvarnonz[*v]);
   SCIPfreeBlockMemoryArray(scip, &(consdata->row[*v]), consdata->nvarnonz[*v]);
   SCIPfreeBlockMemoryArray(scip, &(consdata->col[*v]), consdata->nvarnonz[*v]);

   /* fill the empty spot of the (multi-)aggregated variable with the last variable of this constraint (as they don't have to be sorted) */
   SCIP_CALL( SCIPreleaseVar(scip, &consdata->vars[*v]) );
   consdata->col[*v] = consdata->col[consdata->nvars - 1];
   consdata->row[*v] = consdata->row[consdata->nvars - 1];
   consdata->val[*v] = consdata->val[consdata->nvars - 1];
   consdata->nvarnonz[*v] = consdata->nvarnonz[consdata->nvars - 1];
   consdata->vars[*v] = consdata->vars[consdata->nvars - 1];
   (consdata->nvars)--;
   (*v)--; /* we need to check again if the variable we just shifted to this position also needs to be (multi-)aggregated */

   /* iterate over all variables this was aggregated to and insert the corresponding nonzeros */
   for (aggrind = 0; aggrind < naggrvars; aggrind++)
   {
      /* check if the variable already exists in this block */
      aggrconsind = -1;
      for (i = 0; i < consdata->nvars; i++)
      {
         if ( consdata->vars[i] == aggrvars[aggrind] )
         {
            aggrconsind = i;
            break;
         }
      }

      if ( aggrconsind > -1 )
      {
         /* if the varialbe to aggregate to is already part of this sdp-constraint, just add the nonzeros of the old variable to it */

         /* resize the arrays to the maximally needed length */
         aggrtargetlength = consdata->nvarnonz[aggrconsind] + *nfixednonz - startind;
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->row[aggrconsind]), consdata->nvarnonz[aggrconsind], aggrtargetlength) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->col[aggrconsind]), consdata->nvarnonz[aggrconsind], aggrtargetlength) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->val[aggrconsind]), consdata->nvarnonz[aggrconsind], aggrtargetlength) );

         if ( SCIPisEQ(scip, constant, 0.0) )
         {
            /* in this case we saved the original values in savedval, we add startind to the pointers to only add those from
             * the current variable, the number of entries is the current position minus the position whre we started */
            SCIP_CALL( SCIPsdpVarfixerMergeArrays(SCIPblkmem(scip), epsilon, savedrow + startind, savedcol + startind, savedval + startind,
                        *nfixednonz - startind, TRUE, scalars[aggrind], consdata->row[aggrconsind], consdata->col[aggrconsind],
                        consdata->val[aggrconsind], &(consdata->nvarnonz[aggrconsind]), aggrtargetlength) );
         }
         else
         {
            /* in this case we saved the original values * constant, so we now have to divide by constant, we add startind to the pointers
             * to only add those from the current variable, the number of entries is the current position minus the position whre we started */
            SCIP_CALL( SCIPsdpVarfixerMergeArrays(SCIPblkmem(scip), epsilon, savedrow + startind, savedcol + startind, savedval + startind,
                        *nfixednonz - startind, TRUE, scalars[aggrind] / constant, consdata->row[aggrconsind], consdata->col[aggrconsind],
                        consdata->val[aggrconsind], &(consdata->nvarnonz[aggrconsind]), aggrtargetlength) );
         }

         /* shrink them again if nonzeros could be combined */
         assert( consdata->nvarnonz[aggrconsind] <= aggrtargetlength );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->row[aggrconsind]), aggrtargetlength, consdata->nvarnonz[aggrconsind]) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->col[aggrconsind]), aggrtargetlength, consdata->nvarnonz[aggrconsind]) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->val[aggrconsind]), aggrtargetlength, consdata->nvarnonz[aggrconsind]) );
      }
      else
      {
         /* the variable has to be added to this constraint */

         SCIPdebugMessage("adding variable %s to SDP constraint %s because of (multi-)aggregation\n", SCIPvarGetName(aggrvars[aggrind]), SCIPconsGetName(cons));

         /* check if we have to enlarge the arrays */
         if ( consdata->nvars == *vararraylength )
         {
            globalnvars = SCIPgetNVars(scip);

            /* we don't want to enlarge this by one for every variable added, so we immediately set it to the maximum possible size */
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->col, *vararraylength, globalnvars) );
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->row, *vararraylength, globalnvars) );
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->val, *vararraylength, globalnvars) );
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->nvarnonz, *vararraylength, globalnvars) );
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->vars, *vararraylength, globalnvars) );
            *vararraylength = globalnvars;
         }

         /* we insert this variable at the last position, as the ordering doesn't matter */
         SCIP_CALL( SCIPcaptureVar(scip, aggrvars[aggrind]) );
         consdata->vars[consdata->nvars] = aggrvars[aggrind];

         /* as there were no nonzeros thus far, we can just duplicate the saved arrays to get the nonzeros for the new variable */
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(consdata->col[consdata->nvars]), savedcol + startind, *nfixednonz - startind) ); /*lint !e776*/
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(consdata->row[consdata->nvars]), savedrow + startind, *nfixednonz - startind) ); /*lint !e776*/

         /* if scalars[aggrind] = constant, we would multiply with 1, if constant = 0, we didn't divide by constant, so in these cases, we can just
          * memcopy the array of nonzero-values */
         /* TODO: only checking scalar and constant for feas-equality might lead to big differences, if the nonzeros they are multiplied with are big,
          * e.g. scalar = 0, constant = 10^(-6), nonzero = 10^(10) leads to new nonzero of 10^4 instead of 0 */
         if ( SCIPisEQ(scip, scalars[aggrind], constant) || SCIPisEQ(scip, constant, 0.0) )
         {
            SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(consdata->val[consdata->nvars]), savedval + startind, *nfixednonz - startind) ); /*lint !e776*/
            consdata->nvarnonz[consdata->nvars] = *nfixednonz - startind; /* as there were no nonzeros thus far, the number of nonzeros equals
                                                                           * the number of nonzeros of the aggregated variable */
         }
         else  /* we have to multiply all entries by scalar before inserting them */
         {
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(consdata->val[consdata->nvars]), *nfixednonz - startind) );

            consdata->nvarnonz[consdata->nvars] = 0;

            for (i = 0; i < *nfixednonz - startind; i++)
            {
               /* if both scalar and savedval are small this might become too small */
               if ( (scalars[i] / constant) * savedval[startind + i] >= epsilon )  /*lint !e679*/
               {
                  consdata->val[consdata->nvars][consdata->nvarnonz[consdata->nvars]] = (scalars[i] / constant) * savedval[startind + i]; /*lint !e679*/
                  consdata->nvarnonz[consdata->nvars]++;
               }
            }
         }

         if ( consdata->nvarnonz[consdata->nvars] > 0 ) /* if scalar and all savedvals were to small */
            consdata->nvars++;
      }
   }

   /* if the constant part is zero, we may delete the nonzero-entries from the saved arrays (by resetting the nfixednonz entry to where
    * it started, so that these entries will be overwritten */
   if ( SCIPisEQ(scip, constant, 0.0) )
      *nfixednonz = startind;

   return SCIP_OKAY;
}


/** presolve routine that looks through the data and handles fixed, (multi-)aggregated and negated variables */
static
SCIP_RETCODE fixAndAggrVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< array with constraints to check */
   int                   nconss,             /**< number of constraints to check */
   SCIP_Bool             aggregate           /**< do we want to (mutli-)aggregate variables ? */
   )
{
   SCIP_CONSDATA* consdata;
   int i;
   int nfixednonz;
   int* savedcol;
   int* savedrow;
   SCIP_Real* savedval;
   int c;
   int v;
   int arraylength;
   SCIP_VAR* var;
   SCIP_VAR** aggrvars;
   SCIP_Real scalar;
   SCIP_Real* scalars;
   int naggrvars;
   SCIP_Real constant;
   int requiredsize;
   int globalnvars;
   int vararraylength;
   SCIP_Bool negated;
   SCIP_Real epsilon;


   /* loop over all variables once, add all fixed to savedrow/col/val, for all multiaggregated variables, if constant-scalar =!= 0, add
    * constant-scalar * entry to savedrow/col/val and call mergeArrays for all aggrvars for savedrow[startindex of this var] and scalar/constant-scalar,
    * if constant-scalar = 0, add 1*entry to savedrow/col/val, call mergeArrays for all aggrvars for savedrow[startindex of this var] and scalar and later
    * reduze the saved size of savedrow/col/val by the number of nonzeros of the mutliagrregated variable to not add them to the constant part later */

   assert( scip != NULL );
   assert( conss != NULL );
   assert( nconss >= 0 );

   SCIPdebugMessage("Calling fixAndAggrVars with aggregate = %u\n", aggregate);

   SCIP_CALL( SCIPgetRealParam(scip, "numerics/epsilon", &epsilon) );

   for (c = 0; c < nconss; ++c)
   {
      assert( conss[c] != NULL );
      assert( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "SDP") == 0);

      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );

      /* allocate memory to save nonzeros that need to be fixed */
      SCIP_CALL( SCIPallocBufferArray(scip, &savedcol, consdata->nnonz) );
      SCIP_CALL( SCIPallocBufferArray(scip, &savedrow, consdata->nnonz) );
      SCIP_CALL( SCIPallocBufferArray(scip, &savedval, consdata->nnonz) );

      /* initialize this with zero for each block */
      nfixednonz = 0;

      vararraylength = consdata->nvars;
      globalnvars = SCIPgetNVars(scip);

      for (v = 0; v < consdata->nvars; v++)/*lint --e{850}*/
      {
         negated = FALSE;
         /* if the variable is negated, get the negation var */
         if ( SCIPvarIsBinary(consdata->vars[v]) && SCIPvarIsNegated(consdata->vars[v]) )
         {
            negated = TRUE;
            var = SCIPvarGetNegationVar(consdata->vars[v]);
         }
         else
            var = consdata->vars[v];

         /* check if the variable is fixed in SCIP */
         if ( SCIPvarGetStatus(var) == SCIP_VARSTATUS_FIXED || SCIPisEQ(scip, SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var)))
         {
            assert( SCIPisEQ(scip, SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var)) );

            SCIPdebugMessage("Globally fixing Variable %s to value %f !\n", SCIPvarGetName(var), SCIPvarGetLbGlobal(var));

            if ( ((! negated) && (! SCIPisEQ(scip, SCIPvarGetLbGlobal(var), 0.0))) || (negated && SCIPisEQ(scip, SCIPvarGetLbGlobal(var), 0.0)) )
            {
               /* the nonzeros are saved to later be inserted into the constant part (this is only done after all nonzeros of fixed variables have been
                * assembled, because we need to sort the constant nonzeros and loop over them, which we only want to do once and not once for each fixed
                * variable) */
               for (i = 0; i < consdata->nvarnonz[v]; i++)
               {
                  savedcol[nfixednonz] = consdata->col[v][i];
                  savedrow[nfixednonz] = consdata->row[v][i];
                  /* this is the final value to add, we no longer have to remember from which variable this comes, minus because we have +A_i but -A_0 */
                  if ( ! negated )
                     savedval[nfixednonz] = consdata->val[v][i] * SCIPvarGetLbGlobal(var);
                  else
                     savedval[nfixednonz] = consdata->val[v][i]; /* if it is the negation of a variable fixed to zero, this variable is fixed to one */

                  nfixednonz++;
                  consdata->nnonz--;
               }
            }
            else
            {
               /* if the variable is fixed to zero, the nonzeros will just vanish, so we only reduce the number of nonzeros */
               consdata->nnonz -= consdata->nvarnonz[v];
            }
            /* free the memory of the corresponding entries in col/row/val */
            SCIPfreeBlockMemoryArrayNull(scip, &(consdata->val[v]), consdata->nvarnonz[v]);
            SCIPfreeBlockMemoryArrayNull(scip, &(consdata->row[v]), consdata->nvarnonz[v]);
            SCIPfreeBlockMemoryArrayNull(scip, &(consdata->col[v]), consdata->nvarnonz[v]);

            /* as the variables don't need to be sorted, we just put the last variable into the empty spot and decrease sizes by one (at the end) */
            SCIP_CALL( SCIPreleaseVar(scip, &(consdata->vars[v])) );
            consdata->col[v] = consdata->col[consdata->nvars - 1];
            consdata->row[v] = consdata->row[consdata->nvars - 1];
            consdata->val[v] = consdata->val[consdata->nvars - 1];
            consdata->nvarnonz[v] = consdata->nvarnonz[consdata->nvars - 1];
            consdata->vars[v] = consdata->vars[consdata->nvars - 1];
            consdata->nvars--;
            v--; /* we need to check again if the variable we just shifted to this position also needs to be fixed */
         }
         else if ( (SCIPvarGetStatus(var) == SCIP_VARSTATUS_AGGREGATED ||
                   SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR)
                  && aggregate )
         {
            SCIP_CALL( SCIPallocBufferArray(scip, &aggrvars, globalnvars) );
            SCIP_CALL( SCIPallocBufferArray(scip, &scalars, globalnvars) );

            /* this is how they should be initialized before calling SCIPgetProbvarLinearSum */
            if (! negated)
            {
               aggrvars[0] = consdata->vars[v];
               naggrvars = 1;
               constant = 0.0;
               scalars[0] = 1.0;
            }
            else
            {
               /* if this variable is the negation of var, than we look for a representation of 1-var */
               aggrvars[0] = consdata->vars[v];
               naggrvars = 1;
               constant = 1.0;
               scalars[0] = -1.0;
            }

            /* get the variables this var was aggregated to */
            SCIP_CALL( SCIPgetProbvarLinearSum(scip, aggrvars, scalars, &naggrvars, globalnvars, &constant, &requiredsize, TRUE) );
            assert( requiredsize <= globalnvars ); /* requiredsize is the number of empty spots in aggrvars needed, globalnvars is the number
                                                    * of spots we provided */

            /* Debugmessages for the (multi-)aggregation */
#ifdef SCIP_DEBUG
            if ( SCIPvarGetStatus(consdata->vars[v]) == SCIP_VARSTATUS_AGGREGATED )
               SCIPdebugMessage("aggregating variable %s to ", SCIPvarGetName(var));
            else
               SCIPdebugMessage("multiaggregating variable %s to ", SCIPvarGetName(var));
            for (i = 0; i < naggrvars; i++)
               SCIPdebugMessage("+ (%f2) * %s ", scalars[i], SCIPvarGetName(aggrvars[i]));
            SCIPdebugMessage("+ (%f2) \n", constant);
#endif

            /* add the nonzeros to the saved-arrays for the constant part, remove the nonzeros for the old variables and add them to the variables this variable
             * was (multi-)aggregated to */
            SCIP_CALL( multiaggrVar(scip, conss[c], &v, aggrvars, scalars, naggrvars, constant, savedcol, savedrow, savedval, &nfixednonz, &vararraylength) );

            SCIPfreeBufferArray(scip, &aggrvars);
            SCIPfreeBufferArray(scip, &scalars);
         }
         else if ( negated && (SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN)
                  && aggregate)
         {
             /* if var1 is the negation of var2, then this is equivalent to it being aggregated to -var2 + 1 = 1 - var2 */

            SCIPdebugMessage("Changing variable %s to negation of variable %s !\n", SCIPvarGetName(consdata->vars[v]), SCIPvarGetName(var));

            scalar = -1.0;

            SCIP_CALL( multiaggrVar(scip, conss[c], &v, &var, &scalar, 1, 1.0, savedcol, savedrow, savedval, &nfixednonz, &vararraylength) );
         }
      }

      /* shrink the variable arrays if they were enlarged too much (or more vars were removed than added) */
      assert( consdata->nvars <= vararraylength );
      if ( consdata->nvars < vararraylength )
      {
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->col, vararraylength, consdata->nvars) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->row, vararraylength, consdata->nvars) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->val, vararraylength, consdata->nvars) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->nvarnonz, vararraylength, consdata->nvars) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->vars, vararraylength, consdata->nvars) );
      }

      /* allocate the maximally needed memory for inserting the fixed variables into the constant part */
      arraylength = consdata->constnnonz + nfixednonz;
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->constcol), consdata->constnnonz, arraylength) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->constrow), consdata->constnnonz, arraylength) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->constval), consdata->constnnonz, arraylength) );

      /* insert the fixed variables into the constant arrays, as we have +A_i but -A_0 we mutliply them by -1 */
      SCIP_CALL( SCIPsdpVarfixerMergeArrays(SCIPblkmem(scip), epsilon, savedrow, savedcol, savedval, nfixednonz, FALSE, -1.0, consdata->constrow,
            consdata->constcol, consdata->constval, &(consdata->constnnonz), arraylength) );

      assert( consdata->constnnonz <= arraylength ); /* the allocated memory should always be sufficient */

      /* shrink the arrays if nonzeros could be combined */
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->constcol), arraylength, consdata->constnnonz) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->constrow), arraylength, consdata->constnnonz) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->constval), arraylength, consdata->constnnonz) );

      /* free the saved arrays */
      SCIPfreeBufferArray(scip, &savedval);
      SCIPfreeBufferArray(scip, &savedrow);
      SCIPfreeBufferArray(scip, &savedcol);

      /* recompute sdpnnonz */
      consdata->nnonz = 0;
      for (v = 0; v < consdata->nvars; v++)
         consdata->nnonz += consdata->nvarnonz[v];
   }

   return SCIP_OKAY;
}

/** enforces the SDP constraints for a given solution (may be NULL to use the LP solution) */
static
SCIP_RETCODE EnforceConstraint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to process */
   int                   nconss,             /**< number of constraints */
   SCIP_SOL*             sol,                /**< solution to enforce (NULL for the LP solution) */
   SCIP_RESULT*          result              /**< pointer to store the result of the enforcing call */
   )
{
   char cutname[SCIP_MAXSTRLEN];
   SCIP_CONSDATA* consdata;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool all_feasible = TRUE;
   SCIP_Bool separated = FALSE;
   SCIP_ROW* row;
   SCIP_Bool infeasible;
   SCIP_Real lhs;
   SCIP_Real* coeff;
   SCIP_Real rhs;
   int nvars;
   int i;
   int j;
#ifndef NDEBUG
   int snprintfreturn; /* used to check the return code of snprintf */
#endif

   *result = SCIP_FEASIBLE;

   for (i = 0; i < nconss; ++i)
   {
      consdata = SCIPconsGetData(conss[i]);
      SCIP_CALL( SCIPconsSdpCheckSdpCons(scip, conss[i], sol, 0, 0, 0, result) );
      if ( *result == SCIP_FEASIBLE )
         continue;

      all_feasible = FALSE;

      nvars = consdata->nvars;
      lhs = 0.0;
      coeff = NULL;

      SCIP_CALL( SCIPallocBufferArray(scip, &coeff, nvars) );
      SCIP_CALL( cutUsingEigenvector(scip, conss[i], sol, coeff, &lhs) );

      rhs = SCIPinfinity(scip);
      conshdlrdata = SCIPconshdlrGetData(conshdlr);

#ifndef NDEBUG
      snprintfreturn = SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "sepa_eig_sdp_%d", ++(conshdlrdata->neigveccuts));
      assert( snprintfreturn < SCIP_MAXSTRLEN ); /* check whether the name fits into the string */
#else
      (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "sepa_eig_sdp_%d", ++(conshdlrdata->neigveccuts));
#endif

      SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, conshdlr, cutname , lhs, rhs, FALSE, FALSE, TRUE) );
      SCIP_CALL( SCIPcacheRowExtensions(scip, row) );

      for (j = 0; j < nvars; ++j)
      {
         SCIP_CALL( SCIPaddVarToRow(scip, row, consdata->vars[j], coeff[j]) );
      }

      SCIP_CALL( SCIPflushRowExtensions(scip, row) );

#ifdef SCIP_MORE_DEBUG
      SCIPinfoMessage(scip, NULL, "Added cut %s: ", cutname);
      SCIPinfoMessage(scip, NULL, "%f <= ", lhs);
      for (j = 0; j < nvars; j++)
         SCIPinfoMessage(scip, NULL, "+ (%f)*%s", coeff[j], SCIPvarGetName(consdata->vars[j]));
      SCIPinfoMessage(scip, NULL, "\n");
#endif

      SCIP_CALL( SCIPaddCut(scip, sol, row, FALSE, &infeasible) );

      if ( infeasible )
      {
         *result = SCIP_CUTOFF;

         SCIP_CALL( SCIPreleaseRow(scip, &row) );
         SCIPfreeBufferArray(scip, &coeff);

         return SCIP_OKAY;
      }
      else
      {
         SCIP_CALL( SCIPaddPoolCut(scip, row) );

         SCIP_CALL( SCIPresetConsAge(scip, conss[i]) );
         *result = SCIP_SEPARATED;
         separated = TRUE;
      }
      SCIP_CALL( SCIPreleaseRow(scip, &row) );
      SCIPfreeBufferArray(scip, &coeff);
   }

   if ( all_feasible )
      return SCIP_OKAY;

   if ( separated )
      *result = SCIP_SEPARATED;

   return SCIP_OKAY;
}

/** informs constraint handler that the presolving process is being started */
static
SCIP_DECL_CONSINITPRE(consInitpreSdp)
{/*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   conshdlrdata->neigveccuts = 0; /* this is used to give the eigenvector-cuts distinguishable names */
   conshdlrdata->ndiaggezerocuts = 0; /* this is used to give the diagGEzero-cuts distinguishable names */
   conshdlrdata->n1x1blocks = 0; /* this is used to give the lp constraints resulting from 1x1 sdp-blocks distinguishable names */

   return SCIP_OKAY;
}

/** locks a variable up if the corresponding constraint matrix is not positive semidefinite, locks it down if it is not negative semidefinite */
static
SCIP_DECL_CONSLOCK(consLockSdp)
{/*lint --e{715}*/
   SCIP_Real* Aj;
   SCIP_CONSDATA* consdata;
   int blocksize;
   int var;
   int nvars;
   SCIP_Real eigenvalue;

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   blocksize = consdata->blocksize;
   nvars = consdata->nvars;

   SCIP_CALL( SCIPallocBufferArray(scip, &Aj, blocksize * blocksize) ); /*lint !e647*/

   for (var = 0; var < nvars; var++)
   {
      SCIP_CALL( SCIPconsSdpGetFullAj(scip, cons, var, Aj) );

      /* compute the smallest eigenvalue */
      SCIP_CALL( SCIPlapackComputeIthEigenvalue(SCIPbuffer(scip), FALSE, blocksize, Aj, 1, &eigenvalue, NULL) );
      if ( SCIPisNegative(scip, eigenvalue) )
      {
         /* as the lowest eigenvalue is negative, the matrix is not positive semidefinite, so adding more of it can remove positive
          * semidefiniteness of the SDP-matrix */
         SCIP_CALL( SCIPaddVarLocks(scip, consdata->vars[var], nlocksneg, nlockspos) );
      }

      /* if the smallest eigenvalue is already positive, we don't need to compute the biggest one */
      if ( SCIPisPositive(scip, eigenvalue) )
      {
         /* as an eigenvalue is positive, the matrix is not negative semidefinite, so substracting more of it can remove positive
          * semidefiniteness of the SDP-matrix */
         SCIP_CALL( SCIPaddVarLocks(scip, consdata->vars[var], nlockspos, nlocksneg) );
      }
      else
      {
         /* compute the biggest eigenvalue */
         SCIP_CALL( SCIPlapackComputeIthEigenvalue(SCIPbuffer(scip), FALSE, blocksize, Aj, blocksize, &eigenvalue, NULL) );
         if ( SCIPisPositive(scip, eigenvalue) )
         {
            /* as the biggest eigenvalue is positive, the matrix is not negative semidefinite, so substracting more of it can remove positive
             * semidefiniteness of the SDP-matrix */
            SCIP_CALL( SCIPaddVarLocks(scip, consdata->vars[var], nlockspos, nlocksneg) );
         }
      }
   }

   SCIPfreeBufferArray(scip, &Aj);

   return SCIP_OKAY;
}

/** after presolving variables are fixed and multiaggregated */
static
SCIP_DECL_CONSEXITPRE(consExitpreSdp)
{/*lint --e{715}*/
   assert( scip != NULL );

   if ( conss == NULL )
      return SCIP_OKAY;

   SCIP_CALL( fixAndAggrVars(scip, conss, nconss, TRUE) );

   return SCIP_OKAY;
}

/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolSdp)
{/*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert( conshdlr != NULL );
   assert( result != 0 );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   if ( nrounds == 0 )
   {
      SCIP_CALL( move_1x1_blocks_to_lp(scip, conss, nconss, naddconss, ndelconss, result) );
      if ( conshdlrdata->diaggezerocuts )
      {
         SCIP_CALL( diagGEzero(scip, conss, nconss, naddconss) );
      }
      if ( conshdlrdata->diagzeroimplcuts )
      {
         SCIP_CALL( diagZeroImpl(scip, conss, nconss, naddconss) );
      }
   }

   return SCIP_OKAY;
}

/** creates transformed constraint */
static
SCIP_DECL_CONSTRANS(consTransSdp)
{/*lint --e{715}*/
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* targetdata;
#ifdef OMP
   SCIP_CONSHDLRDATA* conshdlrdata;
#endif
#ifndef NDEBUG
   int snprintfreturn; /* used to check the return code of snprintf */
#endif
   int i;
   char transname[SCIP_MAXSTRLEN];

   sourcedata = SCIPconsGetData(sourcecons);
   assert( sourcedata != NULL );

  SCIPdebugMessage("Transforming constraint <%s>\n", SCIPconsGetName(sourcecons));

#ifdef OMP
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   SCIPdebugMessage("Setting number of threads to %d via OpenMP in Openblas.\n", conshdlrdata->nthreads);
   omp_set_num_threads(conshdlrdata->nthreads);
#endif

   SCIP_CALL( SCIPallocBlockMemory(scip, &targetdata) );

   /* copy some general data */
   targetdata->nvars = sourcedata->nvars;
   targetdata->nnonz = sourcedata->nnonz;
   targetdata->blocksize = sourcedata->blocksize;

   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(targetdata->nvarnonz), sourcedata->nvarnonz, sourcedata->nvars) );

   /* copy the non-constant nonzeros */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(targetdata->col), sourcedata->nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(targetdata->row), sourcedata->nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(targetdata->val), sourcedata->nvars) );

   for (i = 0; i < sourcedata->nvars; i++)
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(targetdata->col[i]), sourcedata->col[i], sourcedata->nvarnonz[i]) );
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(targetdata->row[i]), sourcedata->row[i], sourcedata->nvarnonz[i]) );
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(targetdata->val[i]), sourcedata->val[i], sourcedata->nvarnonz[i]) );
   }
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(targetdata->vars), sourcedata->nvars));

   /* copy & transform the vars array */
   for (i = 0; i < sourcedata->nvars; i++)
   {
      targetdata->vars[i] = SCIPvarGetTransVar(sourcedata->vars[i]);
      SCIP_CALL( SCIPcaptureVar(scip, targetdata->vars[i]) );
   }

   /* copy the constant nonzeros */
   targetdata->constnnonz = sourcedata->constnnonz;

   if ( sourcedata->constnnonz > 0 )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(targetdata->constcol), sourcedata->constcol, sourcedata->constnnonz));
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(targetdata->constrow), sourcedata->constrow, sourcedata->constnnonz));
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(targetdata->constval), sourcedata->constval, sourcedata->constnnonz));
   }
   else
   {
      targetdata->constcol = NULL;
      targetdata->constrow = NULL;
      targetdata->constval = NULL;
   }

   /* copy the maxrhsentry */
   targetdata->maxrhsentry = sourcedata->maxrhsentry;

   /* name the transformed constraint */
#ifndef NDEBUG
      snprintfreturn = SCIPsnprintf(transname, SCIP_MAXSTRLEN, "t_%s", SCIPconsGetName(sourcecons));
      assert( snprintfreturn < SCIP_MAXSTRLEN ); /* check whether the name fits into the string */
#else
      (void) SCIPsnprintf(transname, SCIP_MAXSTRLEN, "t_%s", SCIPconsGetName(sourcecons));
#endif

   /* create target constraint */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, transname, conshdlr, targetdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
         SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),  SCIPconsIsLocal(sourcecons),
         SCIPconsIsModifiable(sourcecons), SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons),
         SCIPconsIsStickingAtNode(sourcecons)) );

   return SCIP_OKAY;
}

/** checks feasiblity of constraint, e.g. the positive semidefiniteness */
static
SCIP_DECL_CONSCHECK(consCheckSdp)
{/*lint --e{715}*/
   int i;

   assert( scip != NULL );
   assert( result != NULL );

   *result = SCIP_FEASIBLE;

   for (i = 0; i < nconss; ++i)
   {
      SCIP_CALL( SCIPconsSdpCheckSdpCons(scip, conss[i], sol, checkintegrality, checklprows, printreason, result) );
      if ( *result == SCIP_INFEASIBLE )
         return SCIP_OKAY;
   }

   return SCIP_OKAY;
}

/** enforce pseudo solution method
 *
 *  Returns didnotrun if objinfeasible, computes feasibility otherwise.
 */
static
SCIP_DECL_CONSENFOPS(consEnfopsSdp)
{/*lint --e{715}*/
   int i;

   assert( scip != NULL );
   assert( result != NULL );
   assert( conss != NULL );

   *result = SCIP_DIDNOTRUN;

   if ( objinfeasible )
   {
      SCIPdebugMessage("-> pseudo solution is objective infeasible, return.\n");
      return SCIP_OKAY;
   }

   for (i = 0; i < nconss; ++i)
   {
      SCIP_CALL( SCIPconsSdpCheckSdpCons(scip, conss[i], NULL, 0, 0, 0, result) );

      if (*result == SCIP_INFEASIBLE)
      {
         /* if it is infeasible for one SDP constraint, it is infeasible for the whole problem */
         SCIPdebugMessage("-> pseudo solution infeasible for SDP-constraint %s, return.\n", SCIPconsGetName(conss[i]));
         return SCIP_OKAY;
      }
   }

   SCIPdebugMessage("-> pseudo solution feasible for all SDP-constraints.\n");

   return SCIP_OKAY;
}


/** Enforce lp solution; if some block is not psd, an eigenvector cut is added.
 */
static
SCIP_DECL_CONSENFOLP(consEnfolpSdp)
{/*lint --e{715}*/
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( conss != NULL );
   assert( result != NULL );

   return EnforceConstraint(scip, conshdlr, conss, nconss, NULL, result);
}

/** Enforce relaxation solution; if some block is not psd, an eigenvector cut is added.
 */
static
SCIP_DECL_CONSENFORELAX(consEnforelaxSdp)
{/*lint --e{715}*/
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( conss != NULL );
   assert( result != NULL );

   return EnforceConstraint(scip, conshdlr, conss, nconss, sol, result);
}

/** separates a solution using constraint specific ideas, gives cuts to SCIP */
static
SCIP_DECL_CONSSEPASOL(consSepasolSdp)
{/*lint --e{715}*/
   int i;

   assert( result != NULL );
   *result = SCIP_DIDNOTFIND;

   for (i = 0; i < nusefulconss; ++i)
   {
      SCIP_CALL( separateSol(scip, conshdlr, conss[i], sol, result) );
   }

   return SCIP_OKAY;
}

/** separation method of constraint handler for LP solution */
static
SCIP_DECL_CONSSEPALP(consSepalpSdp)
{/*lint --e{715}*/
   int i;

   assert( result != NULL );
   *result = SCIP_DIDNOTFIND;

   for (i = 0; i < nusefulconss; ++i)
   {
      SCIP_CALL( separateSol(scip, conshdlr, conss[i], NULL, result) );
   }

   return SCIP_OKAY;
}

/** delete method of SDP constrainthandler */
static
SCIP_DECL_CONSDELETE(consDeleteSdp)
{/*lint --e{715}*/
   int i;

   assert( cons != NULL );
   assert( consdata != NULL );

   SCIPdebugMessage("deleting SDP constraint <%s>.\n", SCIPconsGetName(cons));

   for (i = 0; i < (*consdata)->nvars; i++)
   {
      SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->val[i], (*consdata)->nvarnonz[i]);
      SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->row[i], (*consdata)->nvarnonz[i]);
      SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->col[i], (*consdata)->nvarnonz[i]);
   }

   /* release all variables */
   for (i = 0; i < (*consdata)->nvars; i++)
   {
      SCIP_CALL( SCIPreleaseVar(scip, &((*consdata)->vars[i])) );
   }

   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->vars, (*consdata)->nvars);
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->constval, (*consdata)->constnnonz);
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->constrow, (*consdata)->constnnonz);
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->constcol, (*consdata)->constnnonz);
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->val, (*consdata)->nvars);
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->row, (*consdata)->nvars);
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->col, (*consdata)->nvars);
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->nvarnonz, (*consdata)->nvars);
   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}

/** free method of SDP constrainthandler */
static
SCIP_DECL_CONSFREE(consFreeSdp)
{/*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   SCIPfreeMemory(scip, &conshdlrdata);
   SCIPconshdlrSetData(conshdlr, NULL);

   return SCIP_OKAY;
}

/** copy an SDP constraint handler */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopySdp)
{
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   SCIP_CALL( SCIPincludeConshdlrSdp(scip) );

   *valid = TRUE;

   return SCIP_OKAY;
}

/** copy an SDP constraint*/
static
SCIP_DECL_CONSCOPY(consCopySdp)
{/*lint --e{715}*/
   SCIP_CONSDATA* sourcedata;
   SCIP_Bool success;
   SCIP_VAR** targetvars;
   SCIP_VAR* var;
   int i;

   assert( scip != NULL );
   assert( sourcescip != NULL );
   assert( sourcecons != NULL );
   assert( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(sourcecons)), CONSHDLR_NAME) == 0 );
   assert( valid != NULL );

   SCIPdebugMessage("Copying SDP constraint <%s>\n", SCIPconsGetName(sourcecons));

   *valid = TRUE;

   /* as we can only map active variables, we have to make sure, that the constraint contains no fixed or (multi-)aggregated vars, after
    * exitpresolve (stage 6) this should always be the case, earlier than that we need to call fixAndAggrVars */
   if ( SCIPgetStage(sourcescip)  <= SCIP_STAGE_EXITPRESOLVE )
   {
      SCIP_CALL( fixAndAggrVars(sourcescip, &sourcecons, 1, TRUE) );
   }

   sourcedata = SCIPconsGetData(sourcecons);
   assert( sourcedata != NULL );

   SCIP_CALL( SCIPallocBufferArray(scip, &targetvars, sourcedata->nvars) );

   /* map all variables in the constraint */
   for (i = 0; i < sourcedata->nvars; i++)
   {
      SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourcedata->vars[i], &var, varmap, consmap, global, &success) );
      if ( success )
         targetvars[i] = var;
      else
         *valid = FALSE;
   }

   /* create the new constraint, using an adjusted source name if no new name was given */
   if ( name )
   {
#ifndef NDEBUG
   int snprintfreturn; /* used to check the return code of snprintf */
#endif
   char copyname[SCIP_MAXSTRLEN];

      /* name the copied constraint */
   #ifndef NDEBUG
         snprintfreturn = SCIPsnprintf(copyname, SCIP_MAXSTRLEN, "c_%s", name);
         assert( snprintfreturn < SCIP_MAXSTRLEN ); /* check whether the name fits into the string */
   #else
         (void) SCIPsnprintf(copyname, SCIP_MAXSTRLEN, "c_%s", name);
   #endif
      SCIP_CALL( SCIPcreateConsSdp( scip, cons, copyname, sourcedata->nvars, sourcedata->nnonz, sourcedata->blocksize, sourcedata->nvarnonz,
                 sourcedata->col, sourcedata->row, sourcedata->val, targetvars, sourcedata->constnnonz,
                 sourcedata->constcol, sourcedata->constrow, sourcedata->constval) );
   }
   else
   {
#ifndef NDEBUG
   int snprintfreturn; /* used to check the return code of snprintf */
#endif
   char copyname[SCIP_MAXSTRLEN];

      /* name the copied constraint */
   #ifndef NDEBUG
         snprintfreturn = SCIPsnprintf(copyname, SCIP_MAXSTRLEN, "c_%s", SCIPconsGetName(sourcecons));
         assert( snprintfreturn < SCIP_MAXSTRLEN ); /* check whether the name fits into the string */
   #else
         (void) SCIPsnprintf(copyname, SCIP_MAXSTRLEN, "c_%s", SCIPconsGetName(sourcecons));
   #endif
      SCIP_CALL( SCIPcreateConsSdp( scip, cons, SCIPconsGetName(sourcecons), sourcedata->nvars, sourcedata->nnonz, sourcedata->blocksize,
                 sourcedata->nvarnonz, sourcedata->col, sourcedata->row, sourcedata->val, targetvars, sourcedata->constnnonz,
                 sourcedata->constcol, sourcedata->constrow, sourcedata->constval) );
   }

   SCIPfreeBufferArray(scip, &targetvars);

   return SCIP_OKAY;
}

/** print an SDP constraint */
static
SCIP_DECL_CONSPRINT(consPrintSdp)
{/*lint --e{715}*/
#ifdef PRINT_HUMAN_READABLE
   SCIP_CONSDATA* consdata;
   SCIP_Real* fullmatrix;
   int v;
   int i;
   int j;

   assert( scip != NULL );
   assert( cons != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   SCIP_CALL( SCIPallocBufferArray(scip, &fullmatrix, consdata->blocksize * consdata->blocksize) );

   /* print the non-constant matrices, for this they first have to be assembled in fullmatrix */
   for (v = 0; v < consdata->nvars; v++)
   {
      /* assemble the matrix */

      /* first initialize it with zero */
      for (i = 0; i < consdata->blocksize; i++)
      {
         for (j = 0; j < consdata->blocksize; j++)
            fullmatrix[i * consdata->blocksize + j] = 0.0; /*lint !e679*/
      }

      /* then add the nonzeros */
      for (i = 0; i < consdata->nvarnonz[v]; i++)
      {
         fullmatrix[consdata->row[v][i] * consdata->blocksize + consdata->col[v][i]] = consdata->val[v][i]; /* lower triangular entry */ /*lint !e679*/
         fullmatrix[consdata->col[v][i] * consdata->blocksize + consdata->row[v][i]] = consdata->val[v][i]; /* upper triangular entry */ /*lint !e679*/
      }

      /* print it */
      SCIPinfoMessage(scip, file, "+\n");
      for (i = 0; i < consdata->blocksize; i++)
      {
         SCIPinfoMessage(scip, file, "( ");
         for (j = 0; j < consdata->blocksize; j++)
            SCIPinfoMessage(scip, file, "%f ", fullmatrix[i * consdata->blocksize + j]); /*lint !e679*/
         SCIPinfoMessage(scip, file, ")\n");
      }
      SCIPinfoMessage(scip, file, "* %s\n", SCIPvarGetName(consdata->vars[v]));
   }

   /* print the constant-matrix */

   /* assemble the matrix */

   /* first initialize it with zero */
   for (i = 0; i < consdata->blocksize; i++)
   {
      for (j = 0; j < consdata->blocksize; j++)
         fullmatrix[i * consdata->blocksize + j] = 0.0; /*lint !e679*/
   }

   /* then add the nonzeros */
   for (i = 0; i < consdata->constnnonz; i++)
   {
      fullmatrix[consdata->constrow[i] * consdata->blocksize + consdata->constcol[i]] = consdata->constval[i]; /* lower triangular entry */ /*lint !e679*/
      fullmatrix[consdata->constcol[i] * consdata->blocksize + consdata->constrow[i]] = consdata->constval[i]; /* upper triangular entry */ /*lint !e679*/
   }

   /* print it */
   SCIPinfoMessage(scip, file, "-\n");
   for (i = 0; i < consdata->blocksize; i++)
   {
      SCIPinfoMessage(scip, file, "( ");
      for (j = 0; j < consdata->blocksize; j++)
         SCIPinfoMessage(scip, file, "%f ", fullmatrix[i * consdata->blocksize + j]); /*lint !e679*/
      SCIPinfoMessage(scip, file, ")\n");
   }
   SCIPinfoMessage(scip, file, ">=0\n");

   SCIPfreeBufferArray(scip, &fullmatrix);

   return SCIP_OKAY;
#else
   SCIP_CONSDATA* consdata;
   int i;
   int v;

   assert( scip != NULL );
   assert( cons != NULL );

   consdata = SCIPconsGetData(cons);

   /* print blocksize */
   SCIPinfoMessage(scip, file, "%d\n", consdata->blocksize);

   /* print A_0 if it exists */
   if ( consdata->constnnonz > 0 )
   {
      SCIPinfoMessage(scip, file, "A_0: ");

      for (i = 0; i < consdata->constnnonz; i++)
      {
         SCIPinfoMessage(scip, file, "(%d,%d):%f, ", consdata->constrow[i], consdata->constcol[i], consdata->constval[i]);
      }
      SCIPinfoMessage(scip, file, "\n");
   }

   /* print other matrices */
   for (v = 0; v < consdata->nvars; v++)
   {
      SCIPinfoMessage(scip, file, "<%s>: ", SCIPvarGetName(consdata->vars[v]));
      for (i = 0; i < consdata->nvarnonz[v]; i++)
      {
         SCIPinfoMessage(scip, file, "(%d,%d):%f, ", consdata->row[v][i], consdata->col[v][i], consdata->val[v][i]);
      }
      /* if this is not the last variable, add a newline */
      if (v < consdata->nvars - 1)
      {
         SCIPinfoMessage(scip, file, "\n");
      }
   }

   return SCIP_OKAY;
#endif
}

/** parse an SDP constraint */
static
SCIP_DECL_CONSPARSE(consParseSdp)
{  /*lint --e{715}*/
   SCIP_Bool parsesuccess;
   SCIP_CONSDATA* consdata = NULL;
   char* pos;
   int currentsize;
   int nvars;
   int i;

   assert( scip != NULL );
   assert( str != NULL );

   nvars = SCIPgetNVars(scip);

   assert( success != NULL );
   *success = TRUE;

   /* create constraint data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) );
   consdata->nvars = 0;
   consdata->nnonz = 0;
   consdata->constnnonz = 0;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->nvarnonz, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->col, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->row, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->val, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->vars, nvars));
   consdata->constcol = NULL;
   consdata->constrow = NULL;
   consdata->constval = NULL;

   /* parse the blocksize */
   parsesuccess = SCIPstrToIntValue(str, &(consdata->blocksize), &pos);
   *success = *success && parsesuccess;

   /* skip whitespace */
   while( isspace((unsigned char)*pos) )
      pos++;

   /* check if there is a constant part */
   if ( pos[0] == 'A' && pos[1] == '_' && pos[2] == '0' )
   {
      pos += 5; /* we skip "A_0: " */

      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->constcol, PARSE_STARTSIZE) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->constrow, PARSE_STARTSIZE) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->constval, PARSE_STARTSIZE) );

      currentsize = PARSE_STARTSIZE;

      /* as long as there is another entry for the constant part, parse it */
      while (pos[0] == '(')
      {
         pos++; /* remove the '(' */

         /* check if we need to enlarge the arrays */
         if ( consdata->constnnonz == currentsize )
         {
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->constcol, currentsize, PARSE_SIZEFACTOR * currentsize) );
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->constrow, currentsize, PARSE_SIZEFACTOR * currentsize) );
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->constval, currentsize, PARSE_SIZEFACTOR * currentsize) );
            currentsize *= PARSE_SIZEFACTOR;
         }

         parsesuccess = SCIPstrToIntValue(pos, &(consdata->constrow[consdata->constnnonz]), &pos);
         *success = *success && parsesuccess;
         assert( consdata->constrow[consdata->constnnonz] < consdata->blocksize );
         pos++; /* remove the ',' */
         parsesuccess = SCIPstrToIntValue(pos, &(consdata->constcol[consdata->constnnonz]), &pos);
         *success = *success && parsesuccess;
         assert( consdata->constcol[consdata->constnnonz] < consdata->blocksize );
         pos += 2; /* remove the "):" */
         parsesuccess = SCIPstrToRealValue(pos, &(consdata->constval[consdata->constnnonz]), &pos);
         *success = *success && parsesuccess;
         pos ++; /* remove the "," */

         /* if we got an entry in the upper triangular part, switch the entries for lower triangular */
         if ( consdata->constcol[consdata->constnnonz] > consdata->constrow[consdata->constnnonz] )
         {
            i = consdata->constcol[consdata->constnnonz];
            consdata->constcol[consdata->constnnonz] = consdata->constrow[consdata->constnnonz];
            consdata->constrow[consdata->constnnonz] = i;
         }

         consdata->constnnonz++;

         /* skip whitespace */
         while( isspace((unsigned char)*pos) )
            pos++;
      }

      /* resize the arrays to their final size */
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->constcol, currentsize, consdata->constnnonz) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->constrow, currentsize, consdata->constnnonz) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->constval, currentsize, consdata->constnnonz) );
   }

   /* skip whitespace */
   while( isspace((unsigned char)*pos) )
      pos++;

   /* parse the non-constant part */

   /* while there is another variable, parse it */
   while (pos[0] == '<')
   {
      /* add the variable to consdata->vars and create the corresponding nonzero arrays */
      SCIP_CALL( SCIPparseVarName(scip, pos, &(consdata->vars[consdata->nvars]), &pos) );
      SCIP_CALL( SCIPcaptureVar(scip, consdata->vars[consdata->nvars]) );

      consdata->nvarnonz[consdata->nvars] = 0;
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(consdata->col[consdata->nvars]), PARSE_STARTSIZE));
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(consdata->row[consdata->nvars]), PARSE_STARTSIZE));
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(consdata->val[consdata->nvars]), PARSE_STARTSIZE));
      consdata->nvars++;
      currentsize = PARSE_STARTSIZE;

      pos += 2; /* remove the ": " */

      /* while there is another entry, parse it */
      while (pos[0] == '(')
      {
         pos++; /* remove the '(' */

         /* check if we need to enlarge the arrays */
         if ( consdata->nvarnonz[consdata->nvars - 1] == currentsize )
         {
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->col[consdata->nvars - 1], currentsize, PARSE_SIZEFACTOR * currentsize) );
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->row[consdata->nvars - 1], currentsize, PARSE_SIZEFACTOR * currentsize) );
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->val[consdata->nvars - 1], currentsize, PARSE_SIZEFACTOR * currentsize) );
            currentsize *= PARSE_SIZEFACTOR;
         }

         parsesuccess = SCIPstrToIntValue(pos, &(consdata->row[consdata->nvars - 1][consdata->nvarnonz[consdata->nvars - 1]]), &pos);
         *success = *success && parsesuccess;
         assert( consdata->row[consdata->nvars - 1][consdata->nvarnonz[consdata->nvars - 1]] < consdata->blocksize );
         pos++; /* remove the ',' */
         parsesuccess = SCIPstrToIntValue(pos, &(consdata->col[consdata->nvars - 1][consdata->nvarnonz[consdata->nvars - 1]]), &pos);
         *success = *success && parsesuccess;
         assert( consdata->col[consdata->nvars - 1][consdata->nvarnonz[consdata->nvars - 1]] < consdata->blocksize );
         pos += 2; /* remove the "):" */
         parsesuccess = SCIPstrToRealValue(pos, &(consdata->val[consdata->nvars - 1][consdata->nvarnonz[consdata->nvars - 1]]), &pos);
         *success = *success && parsesuccess;
         pos ++; /* remove the "," */

         /* if we got an entry in the upper triangular part, switch the entries for lower triangular */
         if ( consdata->col[consdata->nvars - 1][consdata->nvarnonz[consdata->nvars - 1]] >
               consdata->row[consdata->nvars - 1][consdata->nvarnonz[consdata->nvars - 1]] )
         {
            i = consdata->col[consdata->nvars - 1][consdata->nvarnonz[consdata->nvars - 1]];
            consdata->col[consdata->nvars - 1][consdata->nvarnonz[consdata->nvars - 1]] =
                  consdata->row[consdata->nvars - 1][consdata->nvarnonz[consdata->nvars - 1]];
            consdata->row[consdata->nvars - 1][consdata->nvarnonz[consdata->nvars - 1]] = i;
         }

         consdata->nvarnonz[consdata->nvars - 1]++;

         /* skip whitespace */
         while( isspace((unsigned char)*pos) )
            pos++;
      }

      /* resize the arrays to their final size */
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->col[consdata->nvars - 1], currentsize, consdata->nvarnonz[consdata->nvars - 1]) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->row[consdata->nvars - 1], currentsize, consdata->nvarnonz[consdata->nvars - 1]) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->val[consdata->nvars - 1], currentsize, consdata->nvarnonz[consdata->nvars - 1]) );

      /* skip whitespace */
      while( isspace((unsigned char)*pos) )
         pos++;
   }

   /* create the constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate, local, modifiable,
         dynamic, removable, stickingatnode) );

   /* compute maximum rhs entry for later use in the DIMACS Error Norm */
   SCIP_CALL( setMaxRhsEntry(*cons) );

#ifdef SCIP_MORE_DEBUG
   SCIP_CALL( SCIPprintCons(scip, *cons, NULL) );
#endif

   return SCIP_OKAY;
}

/** constraint method of constraint handler which returns the variables */
static
SCIP_DECL_CONSGETVARS(consGetVarsSdp)
{/*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int nvars;
   int i;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( vars != NULL );
   assert( success != NULL );
   assert( varssize >= 0 );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   nvars = consdata->nvars;

   if ( nvars > varssize )
   {
      SCIPdebugMessage("consGetVarsIndicator called for array of size %d, needed size %d.\n", varssize, nvars);
      *success = FALSE;
      return SCIP_OKAY;
   }

   for (i = 0; i < nvars; i++)
      vars[i] = consdata->vars[i];

   *success = TRUE;

   return SCIP_OKAY;
}

/** constraint method of constraint handler which returns the number of variables */
static
SCIP_DECL_CONSGETNVARS(consGetNVarsSdp)
{/*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( nvars != NULL );
   assert( success != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   *nvars = consdata->nvars;
   *success = TRUE;

   return SCIP_OKAY;
}

/** creates the handler for SDP constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrSdp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr = NULL;
   SCIP_CONSHDLRDATA* conshdlrdata = NULL;

   assert( scip != 0 );

   /* allocate memory for the conshdlrdata */
   SCIP_CALL( SCIPallocMemory(scip, &conshdlrdata) );

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpSdp, consEnfopsSdp, consCheckSdp, consLockSdp, conshdlrdata) );

   assert( conshdlr != NULL );

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteSdp) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeSdp) );
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopySdp, consCopySdp) );
   SCIP_CALL( SCIPsetConshdlrInitpre(scip, conshdlr,consInitpreSdp) );
   SCIP_CALL( SCIPsetConshdlrExitpre(scip, conshdlr, consExitpreSdp) );
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolSdp, CONSHDLR_MAXPREROUNDS, CONSHDLR_PRESOLTIMING) );
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpSdp, consSepasolSdp, CONSHDLR_SEPAFREQ,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );
   SCIP_CALL( SCIPsetConshdlrEnforelax(scip, conshdlr, consEnforelaxSdp) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransSdp) );
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintSdp) );
   SCIP_CALL( SCIPsetConshdlrParse(scip, conshdlr, consParseSdp) );
   SCIP_CALL( SCIPsetConshdlrGetVars(scip, conshdlr, consGetVarsSdp) );
   SCIP_CALL( SCIPsetConshdlrGetNVars(scip, conshdlr, consGetNVarsSdp) );

   /* add parameter */
#ifdef OMP
   SCIP_CALL( SCIPaddIntParam(scip, "constraints/SDP/threads", "number of threads used for OpenBLAS",
         &(conshdlrdata->nthreads), TRUE, DEFAULT_NTHREADS, 1, INT_MAX, NULL, NULL) );
#endif
   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/SDP/diaggezerocuts",
         "Should linear cuts enforcing the non-negativity of diagonal entries of SDP-matrices be added?",
         &(conshdlrdata->diaggezerocuts), TRUE, DEFAULT_DIAGGEZEROCUTS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/SDP/diagzeroimplcuts",
         "Should linear cuts enforcing the implications of diagonal entries of zero in SDP-matrices be added?",
         &(conshdlrdata->diagzeroimplcuts), TRUE, DEFAULT_DIAGZEROIMPLCUTS, NULL, NULL) );

   return SCIP_OKAY;
}

/** get the data belonging to a single SDP-constraint
 *
 *  In arraylength the length of the nvarnonz, col, row and val arrays has to be given, if it is not sufficient to store all block-pointers that
 *  need to be inserted, a debug message will be thrown and this variable will be set to the needed length.
 *  constnnonz should give the length of the const arrays, if it is too short it will also give the needed number and a debug message is thrown.
 */
SCIP_RETCODE SCIPconsSdpGetData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< SDP constraint to get data of */
   int*                  nvars,              /**< pointer to store the number of variables in this SDP constraint */
   int*                  nnonz,              /**< pointer to store the number of nonzeros in this SDP constraint */
   int*                  blocksize,          /**< pointer to store the size of this SDP-block */
   int*                  arraylength,        /**< length of the given nvarnonz, col, row and val arrays, if this is too short this will return the needed length*/
   int*                  nvarnonz,           /**< pointer to store the number of nonzeros for each variable, also length of the arrays col/row/val are
                                               *  pointing to */
   int**                 col,                /**< pointer to store the column indices of the nonzeros for each variable */
   int**                 row,                /**< pointer to store the row indices of the nonzeros for each variable */
   SCIP_Real**           val,                /**< pointer to store the values of the nonzeros for each variable */
   SCIP_VAR**            vars,               /**< pointer to store the SCIP variables present in this constraint that correspond to the indices in col/row/val */
   int*                  constnnonz,         /**< pointer to store the number of nonzeros in the constant part of this SDP constraint, also length of
                                               *  the const arrays */
   int*                  constcol,           /**< pointer to store the column indices of the constant nonzeros */
   int*                  constrow,           /**< pointer to store the row indices of the constant nonzeros */
   SCIP_Real*            constval            /**< pointer to store the values of the constant nonzeros */
   )
{
   SCIP_CONSDATA* consdata;
   const char* name;
   int i;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( nvars != NULL );
   assert( nnonz != NULL );
   assert( blocksize != NULL );
   assert( arraylength != NULL );
   assert( nvarnonz != NULL );
   assert( col != NULL );
   assert( row != NULL );
   assert( val != NULL );
   assert( vars != NULL );
   assert( constnnonz != NULL );

   consdata = SCIPconsGetData(cons);
   name = SCIPconsGetName(cons);

   assert( consdata->constnnonz == 0 || ( constcol != NULL && constrow != NULL && constval != NULL ) );

   *nvars = consdata->nvars;
   *nnonz = consdata->nnonz;
   *blocksize = consdata->blocksize;

   for (i = 0; i < consdata->nvars; i++)
      vars[i] = consdata->vars[i];

   /* check that the sdp-arrays are long enough to store the information */
   if ( *arraylength < consdata->nvars )
   {
      SCIPdebugMessage("nvarnonz, col, row and val arrays were not long enough to store the information for cons %s, they need to be at least"
                       "size %d, given was only length %d! \n", name, consdata->nvars, *arraylength);
      *arraylength = consdata->nvars;
   }
   else
   {
      for (i = 0; i < consdata->nvars; i++)
      {
         nvarnonz[i] = consdata->nvarnonz[i];
         /* set the pointers for each variable */
         col[i] = consdata->col[i];
         row[i] = consdata->row[i];
         val[i] = consdata->val[i];
      }
   }

   /* set the constant pointers (if a constant part exists) */
   if ( consdata->constnnonz > 0 )
   {
      if ( consdata->constnnonz > *constnnonz )
      {
         SCIPdebugMessage("The constant nonzeros arrays were not long enough to store the information for cons %s, they need to be at least"
                                "size %d, given was only length %d! \n", name, consdata->constnnonz, *constnnonz);
      }
      else
      {
         for (i = 0; i < consdata->constnnonz; i++)
         {
            constcol[i] = consdata->constcol[i];
            constrow[i] = consdata->constrow[i];
            constval[i] = consdata->constval[i];
         }
      }
   }

   *constnnonz = consdata->constnnonz;

   return SCIP_OKAY;
}

/** gets the number of nonzeros and constant nonzeros for this SDP constraint
 *
 *  Either nnonz or constnnonz may be NULL if only the other one is needed.
 */
SCIP_RETCODE SCIPconsSdpGetNNonz(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< SDP constraint to get number of nonzeros for */
   int*                  nnonz,              /**< pointer to store the number of nonzeros in this SDP constraint */
   int*                  constnnonz          /**< pointer to store the number of nonzeros in the constant part of this SDP constraint */
   )
{
   SCIP_CONSDATA* consdata;

   assert( scip != NULL );
   assert( cons != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   if ( nnonz != NULL )
      *nnonz = consdata->nnonz;

   if ( constnnonz != NULL )
      *constnnonz = consdata->constnnonz;

   return SCIP_OKAY;
}

/** gets the blocksize of the SDP constraint */
int SCIPconsSdpGetBlocksize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< SDP constraint to get blocksize for */
   )
{
   SCIP_CONSDATA* consdata;

   assert( scip != NULL );
   assert( cons != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   return consdata->blocksize;
}

/** gets the full constraint Matrix \f$ A_j \f$ for a given variable j */
SCIP_RETCODE SCIPconsSdpGetFullAj(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< SDP constraint to get matrix for */
   int                   j,                  /**< the variable j to get the corresponding matrix \f$ A_j \f$ for */
   SCIP_Real*            Aj                  /**< pointer to store the full matrix \f$ A_j \f$ */
   )
{
   SCIP_CONSDATA* consdata;
   int blocksize;
   int i;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( j >= 0 );
   assert( Aj != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   blocksize = consdata->blocksize;

   assert( j < consdata->nvars );

   for (i = 0; i < blocksize * blocksize; i++)
      Aj[i] = 0;

   for (i = 0; i < consdata->nvarnonz[j]; i++)
   {
      Aj[consdata->col[j][i] * blocksize + consdata->row[j][i]] = consdata->val[j][i]; /*lint !e679*/
      Aj[consdata->row[j][i] * blocksize + consdata->col[j][i]] = consdata->val[j][i]; /*lint !e679*/
   }

   return SCIP_OKAY;
}

/** gives an n*n-long array with the full constant matrix */
SCIP_RETCODE SCIPconsSdpGetFullConstMatrix(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< SDP constraint to get matrix for */
   SCIP_Real*            mat                 /**< pointer to store the full constant matrix */
   )
{
   SCIP_CONSDATA* consdata;
   int blocksize;
   int i;
   int j;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( mat != NULL );

   consdata = SCIPconsGetData(cons);
   blocksize = consdata->blocksize;

   for (i = 0; i < blocksize; i++)
   {
      for (j = 0; j < blocksize; j++)
         mat[i * blocksize + j] = 0.0; /*lint !e679*/
   }

   for (i = 0; i < consdata->constnnonz; i++)
   {
      mat[consdata->constcol[i] * blocksize + consdata->constrow[i]] = consdata->constval[i]; /*lint !e679*/
      mat[consdata->constrow[i] * blocksize + consdata->constcol[i]] = consdata->constval[i]; /*lint !e679*/
   }

   return SCIP_OKAY;
}

/** gives a 0.5*n*(n+1)-long array with the lower triangular part of the constant matrix indexed by compLowerTriangPos */
SCIP_RETCODE SCIPconsSdpGetLowerTriangConstMatrix(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< SDP constraint to get matrix for */
   SCIP_Real*            mat                 /**< pointer to store the lower triangular part of the constant matrix */
   )
{
   SCIP_CONSDATA* consdata;
   int blocksize;
   int i;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( mat != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   blocksize = consdata->blocksize;

   /* initialize the matrix with 0 */
   for (i = 0; i < (blocksize * (blocksize + 1)) / 2; i++)
      mat[i] = 0.0;

   for (i = 0; i < consdata->constnnonz; i++)
      mat[compLowerTriangPos(consdata->constrow[i], consdata->constcol[i])] = consdata->constval[i];

   return SCIP_OKAY;
}

/** Compute a heuristic guess for a good starting solution \f$ \lambda ^* \cdot I \f$.
 *
 *  The solution is computed as
 *  \f[
 *  \lambda^* = \max \Bigg\{S \cdot \max_{i \in [m]} \{|u_i|, |l_i|\} \cdot \max_{i \in [m]} \|A_i\|_\infty + \|C\|_\infty,
 *  \frac{\max_{i \in [m]} b_i}{S \cdot \min_{i \in [m]} \min_{j, \ell \in [n]} (A_i)_{j\ell} } \Bigg\},
 *  \f]
 *  where \f$ S = \frac{ | \text{nonzero-entries of all } A_i | }{0.5 \cdot \text{ blocksize } (\text{ blocksize } + 1)} \f$
 *  measures the sparsity of the matrices.
 */
SCIP_RETCODE SCIPconsSdpGuessInitialPoint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< the constraint for which the initial point should be constructed */
   SCIP_Real*            lambdastar          /**< pointer to store the guess for the initial point */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real sparsity;
   SCIP_Real maxinfnorm;
   SCIP_Real maxconst;
   SCIP_Real mininfnorm;
   SCIP_Real maxobj;
   SCIP_Real maxbound;
   SCIP_Real primalguess;
   SCIP_Real dualguess;
   SCIP_Real compval;
   int blocksize;
   int i;
   int v;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( lambdastar != NULL );
   assert( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0 );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   /* If there are no nonzeros, we cannot use the usual formula, since it divides through the number of nonzeros. In this case,
    * however, we will not solve an SDP anyways but at most an LP (more likely the problem will be solved in local presolving,
    * if all variables are fixed and not only those in the SDP-part), so we just take the default value of SDPA.
    */
   if ( consdata->nnonz == 0 )
   {
      *lambdastar = 100.0;
      return SCIP_OKAY;
   }

   blocksize = consdata->blocksize;

   sparsity = consdata->nnonz / (0.5 * blocksize * (blocksize + 1));

   /* compute the maximum entry of the A_i */
   maxinfnorm = 0.0;
   mininfnorm = SCIPinfinity(scip);
   for (v = 0; v < consdata->nvars; v++)
   {
      for (i = 0; i < consdata->nvarnonz[v]; i++)
      {
         if ( SCIPisGT(scip, REALABS(consdata->val[v][i]), maxinfnorm ) )
            maxinfnorm = REALABS(consdata->val[v][i]);
         if ( SCIPisLT(scip, REALABS(consdata->val[v][i]), mininfnorm) )
            mininfnorm = REALABS(consdata->val[v][i]);
      }
   }
   maxconst = 0.0;
   for (i = 0; i < consdata->constnnonz; i++)
   {
      if ( SCIPisGT(scip, REALABS(consdata->constval[i]), maxconst ) )
         maxconst = REALABS(consdata->constval[i]);
   }

   assert( SCIPisGT(scip, mininfnorm, 0.0) );

   /* compute maximum b_i and bound */
   maxobj = 0.0;
   maxbound = 0.0;
   for (v = 0; v < consdata->nvars; v++)
   {
      if ( SCIPisGT(scip, REALABS(SCIPvarGetObj(consdata->vars[v])), maxobj) )
         maxobj = REALABS(SCIPvarGetObj(consdata->vars[v]));
      compval = SCIPisInfinity(scip, REALABS(SCIPvarGetUbGlobal(consdata->vars[v]))) ? 1e+6 : REALABS(SCIPvarGetUbGlobal(consdata->vars[v]));
      if ( SCIPisGT(scip, compval, maxbound) )
         maxbound = compval;
      compval = SCIPisInfinity(scip, REALABS(SCIPvarGetLbGlobal(consdata->vars[v]))) ? 1e+6 : REALABS(SCIPvarGetUbGlobal(consdata->vars[v]));
      if ( SCIPisGT(scip, compval, maxbound) )
         maxbound = compval;
   }

   /* if all variables were unbounded, we set the value to 10^6 */
   if ( SCIPisEQ(scip, maxbound, 0.0) )
      maxbound = 1E+6;

   /* compute primal and dual guess */
   primalguess = maxobj / (sparsity * mininfnorm);
   dualguess = sparsity * maxinfnorm * maxbound + maxconst;

   if ( SCIPisGT(scip, primalguess, dualguess) )
      *lambdastar = primalguess;
   else
      *lambdastar = dualguess;

   return SCIP_OKAY;
}

/** creates an SDP-constraint */
SCIP_RETCODE SCIPcreateConsSdp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables in this SDP constraint */
   int                   nnonz,              /**< number of nonzeros in this SDP constraint */
   int                   blocksize,          /**< size of this SDP-block */
   int*                  nvarnonz,           /**< number of nonzeros for each variable, also length of the arrays col/row/val point to */
   int**                 col,                /**< pointer to column indices of the nonzeros for each variable */
   int**                 row,                /**< pointer to row indices of the nonzeros for each variable */
   SCIP_Real**           val,                /**< pointer to values of the nonzeros for each variable */
   SCIP_VAR**            vars,               /**< SCIP_VARiables present in this SDP constraint that correspond to the indices in col/row/val */
   int                   constnnonz,         /**< number of nonzeros in the constant part of this SDP constraint */
   int*                  constcol,           /**< column indices of the constant nonzeros */
   int*                  constrow,           /**< row indices of the constant nonzeros */
   SCIP_Real*            constval            /**< values of the constant nonzeros */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata = NULL;
   int i;
   int j;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( name != NULL );
   assert( nvars >= 0 );
   assert( nnonz >= 0 );
   assert( blocksize >= 0 );
   assert( constnnonz >= 0 );
   assert( nvars == 0 || vars != NULL );
   assert( nnonz == 0 || (nvarnonz != NULL && col != NULL && row != NULL && val != NULL ));
   assert( constnnonz == 0 || (constcol != NULL && constrow != NULL && constval != NULL ));

   conshdlr = SCIPfindConshdlr(scip, "SDP");
   if ( conshdlr == NULL )
   {
      SCIPerrorMessage("SDP constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->nvarnonz, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->col, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->row, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->val, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->constcol, constnnonz) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->constrow, constnnonz) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->constval, constnnonz) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->vars, nvars));

   for (i = 0; i < nvars; i++)
   {
      assert( nvarnonz[i] >= 0 );

      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->col[i], nvarnonz[i]));
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->row[i], nvarnonz[i]));
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->val[i], nvarnonz[i]));
   }

   consdata->nvars = nvars;
   consdata->nnonz = nnonz;
   consdata->constnnonz = constnnonz;
   consdata->blocksize = blocksize;

   for (i = 0; i < nvars; i++)
   {
      consdata->nvarnonz[i] = nvarnonz[i];

      for (j = 0; j < nvarnonz[i]; j++)
      {
         assert( col[i][j] >= 0 );
         assert( col[i][j] < blocksize );
         assert( row[i][j] >= 0 );
         assert( row[i][j] < blocksize );

         consdata->col[i][j] = col[i][j];
         consdata->row[i][j] = row[i][j];
         consdata->val[i][j] = val[i][j];
      }
   }

   for (i = 0; i < constnnonz; i++)
   {
      consdata->constcol[i] = constcol[i];
      consdata->constrow[i] = constrow[i];
      consdata->constval[i] = constval[i];
   }

   for (i = 0; i < nvars; i++)
   {
      consdata->vars[i] = vars[i];
      SCIP_CALL( SCIPcaptureVar(scip, consdata->vars[i]) );
   }
   SCIPdebugMessage("creating cons %s\n", name);
   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* compute maximum rhs entry for later use in the DIMACS Error Norm */
   SCIP_CALL( setMaxRhsEntry(*cons) );

   return SCIP_OKAY;
}
