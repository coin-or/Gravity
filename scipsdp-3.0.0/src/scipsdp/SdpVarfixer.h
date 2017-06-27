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

/**@file   SdpVarfixer.h
 * @brief  adds the main functionality to fix/unfix/(multi-)aggregate variables by merging two three-tuple-arrays of row/col/val together
 * @author Tristan Gally
 */

#ifndef __SDPVARFIXER_H__
#define __SDPVARFIXER_H__

#include "scip/type_misc.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * sort the given row, col and val arrays first by non-decreasing row-indices, then for those with identical
 * row-indices by non-decreasing col-indices
 */
EXTERN
void SCIPsdpVarfixerSortRowCol(
   int*                  row,                /* row indices */
   int*                  col,                /* column indices */
   SCIP_Real*            val,                /* values */
   int                   length              /* length of the given arrays */
   );

/**
 * merges two three-tuple-arrays together
 *
 * The original arrays (which may have multiple entries for the same row and col) will be mulitplied with
 * scalar and then merged into the target arrays (which may not have multiple entries for the same row and col). If there is already an entry for
 * a row/col combination, these two entries will be combined (their values added together), if they cancel each other out the nonzero entry will
 * be removed. If you think of the matrices described by the two arrays, this is a matrix addition (but only working on the nonzeros for efficiency).
 * The target arrays need to be long enough, otherwise targetlength returns the needed amount and a corresponding debug message is thrown.
 */
SCIP_RETCODE SCIPsdpVarfixerMergeArrays(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_Real             epsilon,            /**< only values bigger than this are counted as nonzeros */
   int*                  originrow,          /** original row-index-array that is going to be merged, may be NULL if originlength = 0 */
   int*                  origincol,          /** original column-index-array that is going to be merged, may be NULL if originlength = 0 */
   SCIP_Real*            originval,          /** original nonzero-values-array that is going to be merged, may be NULL if originlength = 0 */
   int                   originlength,       /** length of the original arrays */
   SCIP_Bool             originsorted,       /** are the origin arrays already sorted by non-decreasing row and in case of ties col */
   SCIP_Real             scalar,             /** scalar that the original nonzero-values will be multiplied with before merging */
   int*                  targetrow,          /** row-index-array the original array will be merged into */
   int*                  targetcol,          /** column-index-array the original array will be merged into */
   SCIP_Real*            targetval,          /** nonzero-values-array the original array will be merged into */
   int*                  targetlength,       /** length of the target arrays the original arrays will be merged into, this will be updated to the
                                               * new length after the mergings */
   int                   targetmemory        /** amount of memory allocated for targetrow, -col, -val, if this isn't sufficient targetlength will
                                               * return the needed amount and a corresponding debug message will be thrown */
   );

/**
 * merges two three-tuple-arrays together
 *
 * If there are multiple entries for a row/col combination, these will be combined (their values added
 * together), if they cancel each other out the nonzero entry will be removed. The first arrays are assumed to have unique row/col-combinations, the
 * second arrays may have duplicates of the same row/col-combination. In constrast to MergeArrays, here the combined arrays will be inserted in
 * the new targetarrays, and not overwrite one of the old arrays. targetlength should give the length of the target arrays, if this is not sufficient,
 * the needed length is returned there and a debug message is thrown.
 */
EXTERN
SCIP_RETCODE SCIPsdpVarfixerMergeArraysIntoNew(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_Real             epsilon,            /**< only values bigger than this are counted as nonzeros */
   int*                  firstrow,           /** first row-index-array that is going to be merged, may be NULL if firstlength = 0 */
   int*                  firstcol,           /** first column-index-array that is going to be merged, may be NULL if firstlength = 0 */
   SCIP_Real*            firstval,           /** first nonzero-values-array that is going to be merged, may be NULL if firstlength = 0 */
   int                   firstlength,        /** length of the first arrays */
   int*                  secondrow,          /** second row-index-array that is going to be merged, may be NULL if secondlength = 0 */
   int*                  secondcol,          /** second column-index-array that is going to be merged, may be NULL if secondlength = 0 */
   SCIP_Real*            secondval,          /** second nonzero-values-array that is going to be merged, may be NULL if secondlength = 0 */
   int                   secondlength,       /** length of the second arrays */
   int*                  targetrow,          /** row-index-array the original arrays will be merged into */
   int*                  targetcol,          /** column-index-array the original arrays will be merged into */
   SCIP_Real*            targetval,          /** nonzero-values-array the original arrays will be merged into */
   int*                  targetlength        /** length of the target arrays the original arrays will be merged into, this will be updated to the
                                               * new length after the mergings */
   );

#ifdef __cplusplus
}
#endif

#endif
