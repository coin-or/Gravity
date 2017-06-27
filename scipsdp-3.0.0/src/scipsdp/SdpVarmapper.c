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

/**@file   SdpVarmapper.c
 * @brief  class that maps SCIP variables to SDP indices (the SCIP variables are given SDP indices in the order in which they were inserted)
 * @author Tristan Gally
 */

/*
#ifndef __SDPVARMAPPER_H__
#define __SDPVARMAPPER_H__
*/

#include "scip/scip.h"
#include "scip/type_misc.h" /* for SCIP Hashmap */
#include "SdpVarmapper.h"

/* turn off lint warnings for whole file: */
/*lint --e{788,818}*/

struct Sdpvarmapper
{
   SCIP_VAR**            sdptoscip;          /**< array of SCIP variables indexed by their SDP indices */
   SCIP_HASHMAP*         sciptosdp;          /**< hashmap that maps SCIP variables to their SDP indices */
   int                   nvars;              /**< number of variables saved in this varmapper */
};

/** creates the SDP varmapper */
SCIP_RETCODE SCIPsdpVarmapperCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SdpVarmapper**        varmapper,          /**< Pointer to the varmapper that should be created */
   int                   size                /**< initial size of the sciptosdp-hashmap */
   )
{
   assert ( scip != NULL );
   assert ( varmapper != NULL );
   assert ( size >= 0 );

   SCIP_CALL( SCIPallocBlockMemory(scip, varmapper) );
   (*varmapper)->nvars = 0;
   (*varmapper)->sdptoscip = NULL;

   if ( size == 0 )
   {
      SCIPdebugMessage("SCIPsdpVarmapperCreate called for size 0!\n");

      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPhashmapCreate(&((*varmapper)->sciptosdp), SCIPblkmem(scip), size) );

   return SCIP_OKAY;
}

/** frees the SDP varmapper */
SCIP_RETCODE SCIPsdpVarmapperFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SdpVarmapper**        varmapper           /**< Pointer to the varmapper that should be freed */
   )
{
   int i;

   SCIPdebugMessage("Freeing SdpVarmapper \n");

   assert ( scip != NULL );
   assert ( varmapper != NULL );

   /* release all vars */
   for (i = 0; i < (*varmapper)->nvars; i++)
   {
      SCIP_CALL( SCIPreleaseVar(scip, &((*varmapper)->sdptoscip[i])) );
   }

   if ( (*varmapper)->nvars )
      SCIPhashmapFree(&((*varmapper)->sciptosdp));

   SCIPfreeBlockMemoryArrayNull(scip, &(*varmapper)->sdptoscip, (*varmapper)->nvars);
   SCIPfreeBlockMemory(scip, varmapper);

   return SCIP_OKAY;
}

/** adds the given variables (if not already existent) to the end of the varmapper */
SCIP_RETCODE SCIPsdpVarmapperAddVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SdpVarmapper*         varmapper,          /**< varmapper to add variables to */
   int                   nvars,              /**< number of variables to add to the varmapper */
   SCIP_VAR**            vars                /**< SCIP variables to add to the varmapper */
   )
{  /*lint --e{818}*/
   int i;
   SCIP_Bool reallocneeded; /* we allocate memory to add nvars variables, but if some of them already existed in the varmapper, we don't add them and
                             * should reallocate later */
   int allocsize;

   if ( nvars == 0 )
      return SCIP_OKAY;

   assert ( scip != NULL );
   assert ( varmapper != NULL );
   assert ( nvars >= 0 );
   assert ( vars != NULL );

   allocsize = varmapper->nvars + nvars;
   SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(varmapper->sdptoscip), varmapper->nvars, allocsize) );

   reallocneeded = FALSE;

   for (i = 0; i < nvars; i++)
   {
      if ( ! (SCIPhashmapExists(varmapper->sciptosdp, vars[i])) ) /* make sure, that there are no duplicates in the lists */
      {
         varmapper->sdptoscip[varmapper->nvars] = vars[i];
         SCIP_CALL( SCIPhashmapInsert(varmapper->sciptosdp, (void*) vars[i], (void*) (size_t) varmapper->nvars) );
         varmapper->nvars++;
         SCIP_CALL( SCIPcaptureVar(scip, vars[i]) );
      }
      else
      {
         SCIPdebugMessage("variable %s was not added to the varmapper as it was already part of it \n", SCIPvarGetName(vars[i]));
         reallocneeded = TRUE;
      }
   }

   if ( reallocneeded )
   {
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(varmapper->sdptoscip), allocsize, varmapper->nvars) );
   }

   return SCIP_OKAY;
}

/** adds the given variable (if not already existent) to the varmapper at the given position */
SCIP_RETCODE SCIPsdpVarmapperInsertVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SdpVarmapper*         varmapper,          /**< varmapper to add variable to */
   SCIP_VAR*             var,                /**< SCIP variable to add to the varmapper */
   int                   pos                 /**< position where the variable should be added */
   )
{
   int i;

   assert ( scip != NULL );
   assert ( varmapper != NULL );
   assert ( var != NULL );
   assert ( pos >= 0 );
   assert ( pos <= varmapper->nvars );

   if ( ! SCIPhashmapExists(varmapper->sciptosdp, var) ) /* make sure, that there are no duplicates in the lists */
   {
      if ( pos == varmapper->nvars )   /* add it to the end */
      {
         SCIP_CALL(SCIPsdpVarmapperAddVars(scip, varmapper, 1, &var));
      }
      else
      {
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &varmapper->sdptoscip, varmapper->nvars, varmapper->nvars + 1) );

         /* move all variables after pos one spot to the right to make room for the new one */
         for (i = varmapper->nvars - 1; i >= pos; i--)
         {
            varmapper->sdptoscip[i + 1] = varmapper->sdptoscip[i]; /*lint !e679*/
            SCIP_CALL( SCIPhashmapSetImage(varmapper->sciptosdp, varmapper->sdptoscip[i + 1], (void*) (size_t) (i + 1)) );
         }

         varmapper->sdptoscip[pos] = var;
         SCIP_CALL( SCIPhashmapInsert(varmapper->sciptosdp, var, (void*) (size_t) pos) );
         varmapper->nvars++;
         SCIP_CALL( SCIPcaptureVar(scip, var) );
      }
   }
   else
      SCIPdebugMessage("variable %s was not added to the varmapper as it was already part of it.\n", SCIPvarGetName(var));

   return SCIP_OKAY;
}

/** gets the number of variables */
int SCIPsdpVarmapperGetNVars(
   SdpVarmapper*         varmapper           /**< varmapper to get number of variables for */
   )
{
   assert ( varmapper != NULL );

   return varmapper->nvars;
}

/** Is the given SCIP variable included in the varmapper? */
SCIP_Bool SCIPsdpVarmapperExistsSCIPvar(
   SdpVarmapper*         varmapper,          /**< varmapper to search in */
   SCIP_VAR*             var                 /**< SCIP variable to search for */
   )
{
   assert ( varmapper != NULL );
   assert ( var != NULL );

   return SCIPhashmapExists(varmapper->sciptosdp, var);
}

/** gets the SDP-index for the given SCIP variable */
int SCIPsdpVarmapperGetSdpIndex(
   SdpVarmapper*         varmapper,          /**< varmapper to get variable index for */
   SCIP_VAR*             var                 /**< SCIP variable to get SDP-index for */
   )
{
   assert ( varmapper != NULL );
   assert ( var != NULL );

   return (int) (size_t) SCIPhashmapGetImage(varmapper->sciptosdp, (void*) var);
}

/** gets the corresponding SCIP variable for the given SDP variable-index */
SCIP_VAR* SCIPsdpVarmapperGetSCIPvar(
   SdpVarmapper*         varmapper,          /**< varmapper to extract variable from */
   int                   ind                 /**< index of the SDP-variable */
   )
{
   assert ( varmapper != NULL );
   assert ( 0 <= ind && ind < varmapper->nvars );

   return varmapper->sdptoscip[ind];
}

/** removes the variable for the given SDP-index from the varmapper, decreasing the indices of all later variables by 1 */
SCIP_RETCODE SCIPsdpVarmapperRemoveSdpIndex(
   SCIP*                 scip,               /**< SCIP data structure */
   SdpVarmapper*         varmapper,          /**< varmapper to remove variable from */
   int                   ind                 /**< index of the SDP-variable */
   )
{
   SCIP_VAR* var;
   int i;

   assert ( scip != NULL );
   assert ( varmapper != NULL );
   assert ( 0 <= ind && ind < varmapper->nvars );

   var = varmapper->sdptoscip[ind];

   assert ( SCIPhashmapExists(varmapper->sciptosdp, var) );

   SCIP_CALL( SCIPhashmapRemove(varmapper->sciptosdp, var) );
   SCIP_CALL( SCIPreleaseVar(scip, &(varmapper)->sdptoscip[ind]) );

   /* shift all entries of the sdptoscip-array behind ind one to the left and update their sciptosdp-entries */
   for (i = ind + 1; i < varmapper->nvars; i++)
   {
      varmapper->sdptoscip[i - 1] = varmapper->sdptoscip[i];
      SCIP_CALL( SCIPhashmapSetImage(varmapper->sciptosdp, varmapper->sdptoscip[i - 1], (void*) (size_t) (i - 1)) );
   }

   /* reallocate memory */
   SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &varmapper->sdptoscip, varmapper->nvars, varmapper->nvars - 1) );

   varmapper->nvars--;

   return SCIP_OKAY;
}

/** swaps all SCIP variables for their transformed counterparts */
SCIP_RETCODE SCIPsdpVarmapperTransform(
   SCIP*                 scip,               /**< SCIP data structure */
   SdpVarmapper*         varmapper           /**< pointer to the varmapper that should be transformed */
   )
{
   SCIP_VAR* var;
   int k;

   assert ( scip != NULL );
   assert ( varmapper != NULL );

   for (k = 0; k < varmapper->nvars; ++k)
   {
      SCIP_CALL( SCIPgetTransformedVar(scip, varmapper->sdptoscip[k], &var) );
      SCIP_CALL( SCIPcaptureVar(scip, var) );

      SCIP_CALL( SCIPhashmapRemove(varmapper->sciptosdp, varmapper->sdptoscip[k]) );
      SCIP_CALL( SCIPhashmapInsert(varmapper->sciptosdp, var, (void*) (size_t) k) );

      SCIP_CALL( SCIPreleaseVar(scip, &varmapper->sdptoscip[k]) );

      varmapper->sdptoscip[k] = var;
   }

   return SCIP_OKAY;
}

/** clones the varmapper in the second argument to the varmapper in the third argument */
SCIP_RETCODE SCIPsdpVarmapperClone(
   SCIP*                 scip,               /**< SCIP data structure */
   SdpVarmapper*         oldmapper,          /**< pointer to the varmapper that should be cloned */
   SdpVarmapper*         newmapper           /**< pointer to the varmapper that should become a clone of the other one */
   )
{
   int nvars;
   int i;

   nvars = oldmapper->nvars;

   newmapper->nvars = nvars;

   /* allocate memory */
   SCIP_CALL( SCIPallocBlockMemory(scip, &newmapper) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &newmapper->sdptoscip, nvars) );

   /* copy entries */
   for (i = 0; i < nvars; i++)
   {
      newmapper->sdptoscip[i] = oldmapper->sdptoscip[i];
      SCIP_CALL( SCIPhashmapInsert(newmapper->sciptosdp, oldmapper->sdptoscip[i], (void*) (size_t) i) );
      SCIP_CALL( SCIPcaptureVar(scip, newmapper->sdptoscip[i]) );
   }

   return SCIP_OKAY;
}
