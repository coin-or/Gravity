/*
 * COPYRIGHT (c) 2011 Science and Technology Facilities Council (STFC)
 * Original date 25 Feburary 2011
 * All rights reserved
 *
 * Written by: Jonathan Hogg
 *
 * THIS FILE ONLY may be redistributed under the below modified BSD licence.
 * All other files distributed as part of the HSL_MA86 package
 * require a licence to be obtained from STFC and may NOT be redistributed
 * without permission. Please refer to your licence for HSL_MA86 for full terms
 * and conditions. STFC may be contacted via hsl(at)stfc.ac.uk.
 *
 * Modified BSD licence (this header file only):
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *  * Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *  * Neither the name of STFC nor the names of its contributors may be used
 *    to endorse or promote products derived from this software without
 *    specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL STFC BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
 * OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

#ifndef HSL_MA86D_H
#define HSL_MA86D_H

#ifndef ma86_default_control
#define ma86_control ma86_control_d
#define ma86_info ma86_info_d
#define ma86_default_control ma86_default_control_d
#define ma86_analyse ma86_analyse_d
#define ma86_factor ma86_factor_d
#define ma86_factor_solve ma86_factor_solve_d
#define ma86_solve ma86_solve_d
#define ma86_finalise ma86_finalise_d
#endif

typedef double ma86pkgtype_d_;
typedef double ma86realtype_d_;

/* Data type for user controls */
struct ma86_control_d {
   /* Note: 0 is false, non-zero is true */

   /* C/Fortran interface related controls */
   int f_arrays; /* Treat arrays as 1-based (Fortran) if true or 0-based (C) if
                    false. */

   /* Printing controls */
   int diagnostics_level; /* Controls diagnostic printing.*/
               /* Possible values are:
                   < 0: no printing.
                     0: error and warning messages only.
                     1: as 0 plus basic diagnostic printing.
                     2: as 1 plus some more detailed diagnostic messages.
                     3: as 2 plus all entries of user-supplied arrays.       */
   int unit_diagnostics;   /* unit for diagnostic messages
                              Printing is suppressed if unit_diagnostics < 0. */
   int unit_error;         /* unit for error messages
                              Printing is suppressed if unit_error  <  0.     */
   int unit_warning;       /* unit for warning messages
                              Printing is suppressed if unit_warning  <  0.   */

   /* Controls used by ma86_analyse */
   int nemin;  /* Node amalgamation parameter. A child node is merged with its
                  parent if they both involve fewer than nemin eliminations.*/
   int nb;     /* Controls the size of the blocks used within each node (used to
                  set nb within node_type)*/

   /* Controls used by ma86_factor and ma86_factor_solve */
   int action; /* Keep going even if matrix is singular if true, or abort
                  if false */
   int nbi;    /* Inner block size for use with ma64*/
   int pool_size; /* Size of task pool arrays*/
   ma86realtype_d_ small_; /* Pivots less than small are treated as zero*/
   ma86realtype_d_ static_;/* Control static pivoting*/
   ma86realtype_d_ u;      /* Pivot tolerance*/
   ma86realtype_d_ umin;   /* Minimum pivot tolerance*/
   int scaling;            /* Scaling algorithm to use */
};

/***************************************************/

/* data type for returning information to user.*/
struct ma86_info_d {
   ma86realtype_d_ detlog;       /* Holds logarithm of abs det A (or 0) */
   int detsign;         /* Holds sign of determinant (+/-1 or 0) */
   int flag;            /* Error return flag (0 on success) */
   int matrix_rank;     /* Rank of matrix */
   int maxdepth;        /* Maximum depth of the tree. */
   int num_delay;       /* Number of delayed pivots */
   long num_factor;     /* Number of entries in the factor. */
   long num_flops;      /* Number of flops for factor. */
   int num_neg;         /* Number of negative pivots */
   int num_nodes;       /* Number of nodes */
   int num_nothresh;    /* Number of pivots not satisfying u */
   int num_perturbed;   /* Number of perturbed pivots */
   int num_two;         /* Number of 2x2 pivots */
   int pool_size;       /* Maximum size of task pool used */
   int stat;            /* STAT value on error return -1. */
   ma86realtype_d_ usmall;       /* smallest threshold parameter used */
};

/* Initialise control with default values */
void ma86_default_control_d(struct ma86_control_d *control);
/* Analyse the sparsity pattern and prepare for factorization */
void ma86_analyse_d(const int n, const int ptr[], const int row[], int order[],
      void **keep, const struct ma86_control_d *control,
      struct ma86_info_d *info);
/* To factorize the matrix */
void ma86_factor_d(const int n, const int ptr[], const int row[],
      const ma86pkgtype_d_ val[], const int order[], void **keep,
      const struct ma86_control_d *control, struct ma86_info_d *info,
      const ma86realtype_d_ scale[]);
/* To factorize the matrix AND solve AX = B */
void ma86_factor_solve_d(const int n, const int ptr[], const int row[],
      const ma86pkgtype_d_ val[], const int order[], void **keep,
      const struct ma86_control_d *control, struct ma86_info_d *info,
      const int nrhs, const int ldx, ma86pkgtype_d_ x[],
      const ma86realtype_d_ scale[]);
/* To solve AX = B using the computed factors */
void ma86_solve_d(const int job, const int nrhs, const int ldx,
      ma86pkgtype_d_ *x, const int order[], void **keep,
      const struct ma86_control_d *control, struct ma86_info_d *info,
      const ma86realtype_d_ scale[]);
/* To clean up memory in keep */
void ma86_finalise_d(void **keep, const struct ma86_control_d *control);

#endif
