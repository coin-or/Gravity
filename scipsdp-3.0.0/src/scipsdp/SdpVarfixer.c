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

/**@file   SdpVarfixer.c
 * @brief  adds the main functionality to fix/unfix/(multi-)aggregate variables by merging two three-tuple-arrays of row/col/val together
 * @author Tristan Gally
 */

#include "scip/type_misc.h"
#include "scip/def.h"
#include "scip/pub_misc.h" /* for sorting */
#include "SdpVarfixer.h"

/* turn off lint warnings for whole file: */
/*lint --e{788,818}*/

/**
 * sort the given row, col and val arrays first by non-decreasing row-indices, then for those with identical
 * row-indices by non-decreasing col-indices
 */
void SCIPsdpVarfixerSortRowCol(
   int*                  row,                /**< row indices */
   int*                  col,                /**< column indices */
   SCIP_Real*            val,                /**< values */
   int                   length              /**< length of the given arrays */
   )
{
   int firstentry;
   int nextentry;

   /* first sort by row indices */
   SCIPsortIntIntReal(row, col, val, length);

   /* for those with identical row-indices now sort by non-decreasing col-index, first find all entries with the same row-index */
   nextentry = 0;
   while (nextentry < length)
   {
      firstentry = nextentry; /* the next row starts where the last one ended*/

      while (nextentry < length && row[nextentry] == row[firstentry]) /* as long as the row still matches, increase nextentry */
      {
         nextentry++;
      }

      /* now sort all entries between firstentry and nextentry-1 by their col-indices */
      SCIPsortIntReal(col + firstentry, val + firstentry, nextentry - firstentry);
   }
}

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
   int*                  originrow,          /**< original row-index-array that is going to be merged, may be NULL if originlength = 0 */
   int*                  origincol,          /**< original column-index-array that is going to be merged, may be NULL if originlength = 0 */
   SCIP_Real*            originval,          /**< original nonzero-values-array that is going to be merged, may be NULL if originlength = 0 */
   int                   originlength,       /**< length of the original arrays */
   SCIP_Bool             originsorted,       /**< are the origin arrays already sorted by non-decreasing row and in case of ties col */
   SCIP_Real             scalar,             /**< scalar that the original nonzero-values will be multiplied with before merging */
   int*                  targetrow,          /**< row-index-array the original array will be merged into */
   int*                  targetcol,          /**< column-index-array the original array will be merged into */
   SCIP_Real*            targetval,          /**< nonzero-values-array the original array will be merged into */
   int*                  targetlength,       /**< length of the target arrays the original arrays will be merged into, this will be updated to the
                                              *   new length after the mergings */
   int                   targetmemory        /**< amount of memory allocated for targetrow, -col, -val, if this isn't sufficient targetlength will
                                              *   return the needed amount and a corresponding debug message will be thrown */
   )
{  /*lint --e{679}*/
   /*lint --e{850}*/
   int ind;
   int i;
   int nleftshifts; /* if some nonzeros of the target arrays get deleted, this saves the number of spots the following entries have to be moved
                     * to the left */
   int naddednonz;  /* this gives the number of nonzeros that were added to the end of the arrays (this does NOT include those that were added in
                     * the middle of the arrays by decreasing the number of leftshifts) */
   int insertionpos;
   SCIP_Bool debugmsg; /* should a debug message about insufficient length be thrown */

   assert ( blkmem != NULL );
   assert ( originlength == 0 || (originlength > 0 && originrow != NULL && origincol != NULL && originval != NULL) );
   assert ( targetrow != NULL );
   assert ( targetcol != NULL );
   assert ( targetval != NULL );
   assert ( targetlength != NULL );
   assert ( *targetlength >= 0 );

   /* sort the target and origin arrays first by row and then by col to make searching for entries easier */
   SCIPsdpVarfixerSortRowCol(targetrow, targetcol, targetval, *targetlength);

   if ( ! (originsorted) )
      SCIPsdpVarfixerSortRowCol(originrow, origincol, originval, originlength);

   ind = 0; /* this will be used to traverse the nonzeros of the target arrays */
   naddednonz = 0;
   nleftshifts = 0;
   debugmsg = FALSE;

   /* iterate over all nonzeroes */
   for (i = 0; i < originlength; i++)
   {
      /* search the target arrays for an entry at this position, as both the origin and the target arrays are sorted, we go on until
       * we find an entry that is not < as the current entry in the origin arrays according to this sorting, if this has equal row/col,
       * we have found the entry we have to edit, if it is >, then we know, that there is no identical entry, and we can just add a new
       * entry for this row and col */
      while (ind < *targetlength && (targetrow[ind] < originrow[i] || (targetrow[ind] == originrow[i] && targetcol[ind] < origincol[i])))
      {
         /* shift the target nonzeros to the left if needed */
         if ( nleftshifts > 0 )
         {
            targetrow[ind - nleftshifts] = targetrow[ind];
            targetcol[ind - nleftshifts] = targetcol[ind];
            targetval[ind - nleftshifts] = targetval[ind];
         }
         ind++;
      }

      if ( ind < *targetlength && (targetrow[ind] == originrow[i] && targetcol[ind] == origincol[i]) )
      {
         /* add to the old entry */

         /* shift the entry to the left if needed and change the value */
         if ( nleftshifts > 0 )
         {
            targetrow[ind - nleftshifts] = targetrow[ind];
            targetcol[ind - nleftshifts] = targetcol[ind];
         }
         targetval[ind - nleftshifts] = targetval[ind] + scalar * originval[i];

         /* there could be multiple entries to add with identical row and col, so look for further ones in the next entries until there
          * are no more */
         while (i + 1 < originlength && originrow[i + 1] == targetrow[ind - nleftshifts] && origincol[i + 1] == targetcol[ind - nleftshifts])
         {
            targetval[ind - nleftshifts] += scalar * originval[i + 1];
            i++;
         }

         if ( REALABS(targetval[ind - nleftshifts]) < epsilon )
         {
            /* the nonzero became zero */
            nleftshifts++;
         }
         ind++; /* as we already added all origin-entries belonging to this row/col and also shifted the entry, we can continue with the next one */
      }
      else  /* create a new entry */
      {
         if ( nleftshifts > 0 )
         {
            /* we can add the nonzero at one of the empty spots */
            insertionpos = ind - nleftshifts;
            nleftshifts--; /* as one empty spot was filled, all remaining nonzeros should be moved one less position to the left */
         }
         else
         {
            /* add it to the end */
            insertionpos = *targetlength + naddednonz;
            naddednonz++;
         }

         if ( insertionpos < targetmemory )
         {
            /* add the nonzero to the computed position */
            targetrow[insertionpos] = originrow[i];
            targetcol[insertionpos] = origincol[i];
            targetval[insertionpos] = scalar * originval[i];

            /* there could be multiple entries to add with identical row and col, so look for further ones in the next entries until there are no more */
            while (i + 1 < originlength && originrow[i + 1] == targetrow[insertionpos] && origincol[i + 1] == targetcol[insertionpos])
            {
               targetval[insertionpos] += scalar * originval[i + 1];
               i++;
            }

            /* if there were indeed multiple entries, check if they did cancel each other out, in that case remove the entry */
            if ( REALABS(targetval[insertionpos]) < epsilon )
            {
               /* depending on where this actually zero nonzero was added, either add another leftshift to overwrite it or decrease the number of addednonz */
               if ( insertionpos < ind )
                  nleftshifts++;
               else
                  naddednonz--;
            }
         }
         else  /* the memory was not sufficient, so we will throw a debug message, we only wait until we know the final needed size */
            debugmsg = TRUE;
      }
   }
   /* shift the remaining entries of the target arrays */
   if ( nleftshifts > 0 )
   {
      while (ind < *targetlength + naddednonz && ind < targetmemory)
      {
         targetrow[ind - nleftshifts] = targetrow[ind];
         targetcol[ind - nleftshifts] = targetcol[ind];
         targetval[ind - nleftshifts] = targetval[ind];
         ind++;
      }
   }

   if ( debugmsg )
      SCIPdebugMessage("insufficient memory given for SCIPsdpVarfixerMergeArrays, targetmemorys had length %d, would have needed up to %d\n",
                                    targetmemory, *targetlength + naddednonz);

   *targetlength = *targetlength + naddednonz - nleftshifts;

   return SCIP_OKAY;
}


/**
 * merges two three-tuple-arrays together
 *
 * If there are multiple entries for a row/col combination, these will be combined (their values added
 * together), if they cancel each other out the nonzero entry will be removed. The first arrays are assumed to have unique row/col-combinations, the
 * second arrays may have duplicates of the same row/col-combination. In constrast to MergeArrays, here the combined arrays will be inserted in
 * the new targetarrays, and not overwrite one of the old arrays. targetlength should give the length of the target arrays, if this is not sufficient,
 * the needed length is returned there and a debug message is thrown.
 */
SCIP_RETCODE SCIPsdpVarfixerMergeArraysIntoNew(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_Real             epsilon,            /**< only values bigger than this are counted as nonzeros */
   int*                  firstrow,           /**< first row-index-array that is going to be merged, may be NULL if firstlength = 0 */
   int*                  firstcol,           /**< first column-index-array that is going to be merged, may be NULL if firstlength = 0 */
   SCIP_Real*            firstval,           /**< first nonzero-values-array that is going to be merged, may be NULL if firstlength = 0 */
   int                   firstlength,        /**< length of the first arrays */
   int*                  secondrow,          /**< second row-index-array that is going to be merged, may be NULL if secondlength = 0 */
   int*                  secondcol,          /**< second column-index-array that is going to be merged, may be NULL if secondlength = 0 */
   SCIP_Real*            secondval,          /**< second nonzero-values-array that is going to be merged, may be NULL if secondlength = 0 */
   int                   secondlength,       /**< length of the second arrays */
   int*                  targetrow,          /**< row-index-array the original arrays will be merged into */
   int*                  targetcol,          /**< column-index-array the original arrays will be merged into */
   SCIP_Real*            targetval,          /**< nonzero-values-array the original arrays will be merged into */
   int*                  targetlength        /**< length of the target arrays the original arrays will be merged into, this will be updated to the
                                              *   new length after the mergings */
   )
{
   int firstind;
   int secondind;
   int targetind;
   SCIP_Bool debugmsg; /* should we throw a debug message about insufficient memory */

   assert ( blkmem != NULL );
   assert ( firstlength == 0 || (firstlength > 0 && firstrow != NULL && firstcol != NULL && firstval != NULL ) );
   assert ( secondlength == 0 || (secondlength > 0 && secondrow != NULL && secondcol != NULL && secondval != NULL ) );
   assert ( targetrow != NULL );
   assert ( targetcol != NULL );
   assert ( targetval != NULL );
   assert ( targetlength != NULL );
   assert ( *targetlength >= 0 );

   debugmsg = FALSE;

   /* sort both arrays by non-decreasing row and then col indices to make comparisons easier */
   SCIPsdpVarfixerSortRowCol(firstrow, firstcol, firstval, firstlength);
   SCIPsdpVarfixerSortRowCol(secondrow, secondcol, secondval, secondlength);

   /* as both arrays are sorted, traverse them simultanously, always adding the current entry with the lower index of either array to the
    * target arrays (if they both have the same index, we have found entries that need to be merged) */
   firstind = 0;
   secondind = 0;
   targetind = 0;

   while (firstind < firstlength && secondind < secondlength)
   {
      /* if the next entry of the first arrays comes before the next entry of the second arrays according to the row then col sorting, then we can
       * insert the next entry of the first arrays, as there can't be an entry in the second arrays for the same row/col-combination */
      if ( firstrow[firstind] < secondrow[secondind] || (firstrow[firstind] == secondrow[secondind] && firstcol[firstind] < secondcol[secondind]) )
      {
         if ( targetind < *targetlength )
         {
            targetrow[targetind] = firstrow[firstind];
            targetcol[targetind] = firstcol[firstind];
            targetval[targetind] = firstval[firstind];
         }
         else
            debugmsg = TRUE;

         targetind++;
         firstind++;
      }
      /* if the next entry of the second array comes first, we insert it */
      else if ( firstrow[firstind] > secondrow[secondind] || (firstrow[firstind] == secondrow[secondind] && firstcol[firstind] > secondcol[secondind]) )
      {
         if ( targetind < *targetlength )
         {
            targetrow[targetind] = secondrow[secondind];
            targetcol[targetind] = secondcol[secondind];
            targetval[targetind] = secondval[secondind];
         }
         else
            debugmsg = TRUE;
         secondind++;

         /* as the second arrays may have duplicate entries, we have to check the next entry, if it has the same row/col combination, if yes, then we
          * add it's value to the created entry in the target entries and continue */
         while (secondind < secondlength && (secondrow[secondind] == targetrow[targetind] && secondcol[secondind] == targetcol[targetind]))
         {
            if ( targetind < *targetlength )
               targetval[targetind] += secondval[secondind];
            secondind++;
         }

         /* if we combined multiple fixed nonzeros, it is possible that they cancelled each other out, in that case, we shouldn't add a nonzero to the
          * target arrays (if the array was too short we didn't compute the entry, but we add it, as we want to get an upper bound on the needed size) */
         if ( targetind >= *targetlength || REALABS(targetval[targetind]) >= epsilon )
            targetind++;
      }
      /* if the next entries of both arrays are equal according to the row then col sorting, then they need to be combined */
      else
      {
         if ( targetind < *targetlength )
         {
            targetrow[targetind] = firstrow[firstind];
            targetcol[targetind] = firstcol[firstind];
            targetval[targetind] = firstval[firstind] + secondval[secondind];
         }
         else
            debugmsg = TRUE;
         firstind++;
         secondind++;

         /* as the second arrays may have duplicate entries, we have to check the next entry, if it has the same row/col combination, if yes, then we
          * add it's value to the created entry in the target entries and continue */
         while (secondind < secondlength && (secondrow[secondind] == targetrow[targetind] && secondcol[secondind] == targetcol[targetind]))
         {
            if ( targetind < *targetlength )
               targetval[targetind] += secondval[secondind];
            secondind++;
         }

         /* as we combined multiple entires, it is possible that they cancelled each other out, in that case, we shouldn't add a nonzero to the
          * target arrays (if the array was too short we didn't compute the entry, but we add it, as we want to get an upper bound on the needed size) */
         if ( targetind >= *targetlength || REALABS(targetval[targetind]) >= epsilon )
            targetind++;
      }
   }
   /* if we reach the end of one of the two arrays, we can just add the rest of the other array to the target arrays (in case of the second still
    * combining duplicate entries of this same array), so at most one of the following two while-queues will be non-empty, the contents of these
    * queues are exactly the same as the corresponding if-case in the above while-queue (+ checking for the length of the target arrays) */
   while (firstind < firstlength)
   {
      if ( targetind < *targetlength )
      {
         targetrow[targetind] = firstrow[firstind];
         targetcol[targetind] = firstcol[firstind];
         targetval[targetind] = firstval[firstind];
      }
      else
         debugmsg = TRUE;
      firstind++;
      targetind++;
   }

   while (secondind < secondlength)
   {
      if ( targetind < *targetlength )
      {
         targetrow[targetind] = secondrow[secondind];
         targetcol[targetind] = secondcol[secondind];
         targetval[targetind] = secondval[secondind];
      }
      else
         debugmsg = TRUE;
      secondind++;

      /* as the second arrays may have duplicate entries, we have to check the next entry, if it has the same row/col combination, if yes, then we
       * add it's value to the created entry in the target entries and continue */
      while (secondind < secondlength && (secondrow[secondind] == targetrow[targetind] && secondcol[secondind] == targetcol[targetind]))
      {
         if ( targetind < *targetlength )
            targetval[targetind] += secondval[secondind];
         secondind++;
      }

      /* if we combined multiple fixed nonzeros, it is possible that they cancelled each other out, in that case, we shouldn't add a nonzero to the
       * target arrays */
      if ( targetind >= *targetlength || REALABS(targetval[targetind]) >= epsilon )
         targetind++;
   }

   if ( debugmsg )
   {
      SCIPdebugMessage("Insufficient arraylength in SCIPsdpVarfixerMergeArraysIntoNew, given targetarray-length was %d, would have needed %d",
                                  *targetlength, targetind);
   }

   *targetlength = targetind;

   return SCIP_OKAY;
}
