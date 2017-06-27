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

/**@file   reader_cbf.c
 * @brief  file reader for mixed-integer semidefinite programs in CBF format
 * @author Tristan Gally
 * @author Henrik A. Friberg
 */

/*#define SCIP_DEBUG*/
/*#define SCIP_MORE_DEBUG*/

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>                      /* for strcmp */

#include "scipsdp/reader_cbf.h"
#include "scipsdp/SdpVarmapper.h"
#include "scipsdp/cons_sdp.h"
#include "scip/cons_linear.h"


#define READER_NAME             "cbfreader"
#define READER_DESC             "file reader and writer for MISDPs in cbf format"
#define READER_EXTENSION        "cbf"

#define CBF_VERSION_NR         1         /**< version number for CBF format */
#define CBF_CHECK_NONNEG       TRUE      /**< when writing: check linear constraints and move nonnegativity(-positivity)
                                           *  constraints to definition of variables (which are now defined in non-negative
                                           *  orthant) */
                                          /*  TODO: currently doesn't work for ranged rows (which are not created by sdpa
                                           *  reader) */

/* Use CBF_NAME_FORMAT instead of %s when parsing lines, to avoid buffer overflow. */
#define MACRO_STR_EXPAND(tok) #tok
#define MACRO_STR(tok) MACRO_STR_EXPAND(tok)
#define CBF_NAME_FORMAT "%" MACRO_STR(CBF_MAX_NAME) "s"
#define CBF_MAX_LINE  512       /* Last 3 chars reserved for '\r\n\0' */
#define CBF_MAX_NAME  512

char CBF_LINE_BUFFER[CBF_MAX_LINE];
char CBF_NAME_BUFFER[CBF_MAX_NAME];

struct CBF_Data
{
   int                   firstnegorthvar;    /**< first variable in the createdvars-array which should be non-positive
                                              *   (variables 0 to firstnegorthvar-1 should be non-negative) */
   int                   firstfreevar;       /**< first variable in the createdvars-array which should be free
                                              *   (variables firstnegorthvar to firstfreevar-1 should be non-negative) */
   int                   nvars;              /**< number of variables and length of createdvars-array
                                              *   (variables firstfreevar to nvars-1 should be free) */
   SCIP_VAR**            createdvars;        /**< array of variables created by the CBF reader */
   int                   firstleqcons;       /**< first less or equal constraint in the createdconss-array
                                              *   (constraints 0 to firstleqcons-1 should be greater or equal) */
   int                   firsteqcons;        /**< first equality constraint in the createdconss-array
                                              *   (constrains firstleqcons to firsteqcons-1 should be less or equal) */
   int                   nconss;             /**< number of constraints and length of createdconss-array
                                              *   (constraints firsteqcons to nconss-1 should be equalities) */
   SCIP_CONS**           createdconss;       /**< array of constraints created by the CBF reader */
   int                   nsdpblocks;         /**< number of SDP constraints/blocks */
   int*                  sdpblocksizes;      /**< sizes of the SDP blocks */
   int*                  sdpnblocknonz;      /**< number of nonzeros for each SDP block */
   int*                  sdpnblockvars;      /**< number of variables for each SDP block */
   int**                 nvarnonz;           /**< number of nonzeros for each block and each variable */
   SCIP_VAR***           sdpblockvars;       /**< SCIP variables appearing in each block */
   int**                 sdprow;             /**< array of all row indices for each SDP block */
   int**                 sdpcol;             /**< array of all column indices for each SDP block */
   SCIP_Real**           sdpval;             /**< array of all values of SDP nonzeros for each SDP block */
   int***                rowpointer;         /**< array of pointers to first entries in row-array for each block and variable */
   int***                colpointer;         /**< array of pointers to first entries in row-array for each block and variable */
   SCIP_Real***          valpointer;         /**< array of pointers to first entries in row-array for each block and variable */
   int*                  sdpconstnblocknonz; /**< number of nonzeros for each variable in the constant part, also the i-th entry gives the
                                              *   number of entries of sdpconst row/col/val [i] */
   int**                 sdpconstrow;        /**< pointers to row-indices for each block */
   int**                 sdpconstcol;        /**< pointers to column-indices for each block */
   SCIP_Real**           sdpconstval;        /**< pointers to the values of the nonzeros for each block */
};

typedef struct CBF_Data CBF_DATA;


/*
 * Local methods
 */

/** finds first non-commentary line in given file */
static
SCIP_RETCODE CBFfgets(
   SCIP_FILE*            pFile,              /* file to read from */
   long long int *       linecount           /* current linecount */
   )
{
   assert( pFile != NULL );
   assert( linecount != NULL );

  /* Find first non-commentary line */
  while( SCIPfgets(CBF_LINE_BUFFER, (int) sizeof(CBF_LINE_BUFFER), pFile) != NULL )
  {
     ++(*linecount);

     if (CBF_LINE_BUFFER[0] != '#')
        return SCIP_OKAY;
  }

  return SCIP_READERROR;
}

/** reads objective sense from given CBF-file */
static
SCIP_RETCODE CBFreadObjsense(
   SCIP*                 scip,               /* SCIP data structure */
   SCIP_FILE*            pfile,              /* file to read from */
   long long int *       linecount           /* current linecount */
   )
{
   assert( scip != NULL );
   assert( pfile != NULL );
   assert( linecount != NULL );

   SCIP_CALL( CBFfgets(pfile, linecount) );

   if ( sscanf(CBF_LINE_BUFFER, CBF_NAME_FORMAT, CBF_NAME_BUFFER) == 1 )
   {
      if (strcmp(CBF_NAME_BUFFER, "MIN") == 0)
      {
         SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE) );
      }
      else if (strcmp(CBF_NAME_BUFFER, "MAX") == 0)
      {
         SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MAXIMIZE) );
      }
      else
      {
         SCIPerrorMessage("OBJSENSE should be either MIN or MAX.\n");
         SCIPABORT();
      }
   }
   else
   {
      return SCIP_READERROR;
   }

   return SCIP_OKAY;
}

/** reads the number and type of variables from given CBF-file */
static
SCIP_RETCODE CBFreadVar(
   SCIP*                 scip,               /* SCIP data structure */
   SCIP_FILE*            pfile,              /* file to read from */
   long long int *       linecount,          /* current linecount */
   CBF_DATA*             data                /* data pointer to save the results in */
   )
{
   int nvartypes;
   int nvartypevars;
   int nposorthvars;
   int nnegorthvars;
   int t;
   int v;
   SCIP_VAR* var;
   char varname[SCIP_MAXSTRLEN];
#ifndef NDEBUG
   int nfreevars;
   int snprintfreturn;
#endif

   assert( scip != NULL );
   assert( pfile != NULL );
   assert( linecount != NULL );
   assert( data != NULL );

   SCIP_CALL( CBFfgets(pfile, linecount) );

   if ( sscanf(CBF_LINE_BUFFER, "%i %i", &(data->nvars), &nvartypes) != 2 )
      return SCIP_READERROR;
   else
   {
      assert( data->nvars >= 0 );
      assert( data->nvars == 0 || ((1 <= nvartypes) && (nvartypes <= 3)) );
      /* create variables */
      SCIP_CALL( SCIPallocBufferArray(scip, &(data->createdvars), data->nvars) );
      for (v = 0; v < data->nvars; v++)
      {
#ifndef NDEBUG
         snprintfreturn = SCIPsnprintf(varname, SCIP_MAXSTRLEN, "x_%d", v);
         assert( snprintfreturn < SCIP_MAXSTRLEN);
#else
         (void)SCIPsnprintf(varname, SCIP_MAXSTRLEN, "x_%d", v);
#endif

         SCIP_CALL( SCIPcreateVar(scip, &var, varname,  -SCIPinfinity(scip), SCIPinfinity(scip), 0.0,
               SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL));/*lint !e732*//*lint !e747*/

         SCIP_CALL( SCIPaddVar(scip, var) );
         data->createdvars[v] = var;/*lint !e732*//*lint !e747*/

         /* release variable for the reader. */
         SCIP_CALL( SCIPreleaseVar(scip, &var) );
      }
   }
   nposorthvars = 0;
   nnegorthvars = 0;
#ifndef NDEBUG
   nfreevars = 0;
#endif

   for (t = 0; t < nvartypes; t++)
   {
      SCIP_CALL( CBFfgets(pfile, linecount) );

      if ( sscanf(CBF_LINE_BUFFER, CBF_NAME_FORMAT" %i", CBF_NAME_BUFFER, &nvartypevars) == 2 )
      {
         if (strcmp(CBF_NAME_BUFFER, "L+") == 0)
         {
            if ( nvartypevars >= 0 )
               nposorthvars = nvartypevars;
            else
            {
               SCIPerrorMessage("Number of non-negative variables %d should be non-negative!\n", nvartypevars);
               SCIPABORT();
            }
         }
         else if (strcmp(CBF_NAME_BUFFER, "L-") == 0)
         {
            if ( nvartypevars >= 0 )
               nnegorthvars = nvartypevars;
            else
            {
               SCIPerrorMessage("Number of non-positive variables %d should be non-negative!\n", nvartypevars);
               SCIPABORT();
            }
         }
#ifndef NDEBUG
         else if (strcmp(CBF_NAME_BUFFER, "F") == 0)
         {
            if ( nvartypevars >= 0 )
               nfreevars = nvartypevars;
            else
            {
               SCIPerrorMessage("Number of free variables %d should be non-negative!\n", nvartypevars);
               SCIPABORT();
            }
         }
#endif
         else if (strcmp(CBF_NAME_BUFFER, "F") != 0)
         {
            SCIPerrorMessage("CBF-Reader of SCIP-SDP currently only supports non-negative, non-positive and free variables!\n");
            SCIPABORT();
         }
      }
      else
         return SCIP_READERROR;
   }

   /* set lower bound for non-negative variables */
   data->firstnegorthvar = nposorthvars;
   for (v = 0; v < nposorthvars; v++)
   {
      SCIP_CALL( SCIPchgVarLbGlobal(scip, data->createdvars[v], 0.0) );
   }

   /* set upper bound for non-positive variables */
   data->firstfreevar = data->firstnegorthvar + nnegorthvars;
   for (v = data->firstnegorthvar; v < data->firstfreevar; v++)
   {
      SCIP_CALL( SCIPchgVarUbGlobal(scip, data->createdvars[v], 0.0) );
   }

   assert( data->firstfreevar + nfreevars == data->nvars );

   return SCIP_OKAY;
}

/** reads the number and type of constraints from given CBF-file */
static
SCIP_RETCODE CBFreadCon(
   SCIP*                 scip,               /* SCIP data structure */
   SCIP_FILE*            pfile,              /* file to read from */
   long long int *       linecount,          /* current linecount */
   CBF_DATA*             data                /* data pointer to save the results in */
   )
{
   int nconstypes;
   int nconstypeconss;
   int ngeqconss;
   int nleqconss;
   int t;
   int c;
   SCIP_CONS* cons;
   char consname[SCIP_MAXSTRLEN];
#ifndef NDEBUG
   int neqconss;
   int snprintfreturn;
#endif

   assert( scip != NULL );
   assert( pfile != NULL );
   assert( linecount != NULL );
   assert( data != NULL );

   SCIP_CALL( CBFfgets(pfile, linecount) );

   if ( sscanf(CBF_LINE_BUFFER, "%i %i", &(data->nconss), &nconstypes) != 2 )
      return SCIP_READERROR;
   else
   {
      assert( data->nconss >= 0 );
      assert( data->nconss == 0 || ((1 <= nconstypes) && (nconstypes <= 3)) );
      /* create variables */
      SCIP_CALL( SCIPallocBufferArray(scip, &(data->createdconss), data->nconss) );
      for (c = 0; c < data->nconss; c++)
      {
#ifndef NDEBUG
         snprintfreturn = SCIPsnprintf(consname, SCIP_MAXSTRLEN, "LP_%d", c);
         assert( snprintfreturn < SCIP_MAXSTRLEN);
#else
         (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "linear_%d", c);
#endif

         SCIP_CALL( SCIPcreateConsLinear(scip, &cons, consname, 0, NULL, NULL, -SCIPinfinity(scip), SCIPinfinity(scip),
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));

         SCIP_CALL( SCIPaddCons(scip, cons) );
         data->createdconss[c] = cons;/*lint !e732*//*lint !e747*/

         /* release variable for the reader. */
         SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      }
   }

   ngeqconss = 0;
   nleqconss = 0;
#ifndef NDEBUG
   neqconss = 0;
#endif

   for (t = 0; t < nconstypes; t++)
   {
      SCIP_CALL( CBFfgets(pfile, linecount) );

      if ( sscanf(CBF_LINE_BUFFER, CBF_NAME_FORMAT" %i", CBF_NAME_BUFFER, &nconstypeconss) == 2 )
      {
         if (strcmp(CBF_NAME_BUFFER, "L+") == 0)
         {
            if ( nconstypeconss >= 0 )
               ngeqconss = nconstypeconss;
            else
            {
               SCIPerrorMessage("Number of greater or equal constraints %d should be non-negative!\n", nconstypeconss);
               SCIPABORT();
            }
         }
         else if (strcmp(CBF_NAME_BUFFER, "L-") == 0)
         {
            if ( nconstypeconss >= 0 )
               nleqconss = nconstypeconss;
            else
            {
               SCIPerrorMessage("Number of less or equal constraints %d should be non-negative!\n", nconstypeconss);
               SCIPABORT();
            }
         }
#ifndef NDEBUG
         else if (strcmp(CBF_NAME_BUFFER, "L=") == 0)
         {
            if ( nconstypeconss >= 0 )
               neqconss = nconstypeconss;
            else
            {
               SCIPerrorMessage("Number of equality constraints %d should be non-negative!\n", nconstypeconss);
               SCIPABORT();
            }
         }
#endif
         else if (strcmp(CBF_NAME_BUFFER, "L=") != 0)
         {
            SCIPerrorMessage("CBF-Reader of SCIP-SDP currently only supports linear greater or equal, less or equal and"
                  "equality constraints!\n");
            SCIPABORT();
         }
      }
      else
         return SCIP_READERROR;
   }

   /* set left-hand side for greater or equal constraints */
   data->firstleqcons = ngeqconss;
   for (c = 0; c < ngeqconss; c++)
   {
      SCIP_CALL( SCIPchgLhsLinear(scip, data->createdconss[c], 0.0) );
   }

   /* set right-hand side for less or equal constraints */
   data->firsteqcons = data->firstleqcons + nleqconss;
   for (c = data->firstleqcons; c < data->firsteqcons; c++)
   {
      SCIP_CALL( SCIPchgRhsLinear(scip, data->createdconss[c], 0.0) );
   }

   assert( data->firsteqcons + neqconss == data->nconss );

   return SCIP_OKAY;
}

/** reads integrality conditions from given CBF-file */
static
SCIP_RETCODE CBFreadInt(
   SCIP*                 scip,               /* SCIP data structure */
   SCIP_FILE*            pfile,              /* file to read from */
   long long int *       linecount,          /* current linecount */
   CBF_DATA*             data                /* data pointer to save the results in */
   )
{  /*lint --e{818}*/
   int nintvars;
   int i;
   int v;
   SCIP_Bool infeasible;

   assert( scip != NULL );
   assert( pfile != NULL );
   assert( linecount != NULL );
   assert( data != NULL );

   SCIP_CALL( CBFfgets(pfile, linecount) );

   if ( sscanf(CBF_LINE_BUFFER, "%i", &nintvars) == 1 )
   {
      if ( nintvars >= 0 )
      {
         for (i = 0; i < nintvars; i++)
         {
            SCIP_CALL( CBFfgets(pfile, linecount) );
            if ( sscanf(CBF_LINE_BUFFER, "%i", &v) == 1 )
            {
               SCIP_CALL( SCIPchgVarType(scip, data->createdvars[v], SCIP_VARTYPE_INTEGER, &infeasible)   );

               if ( infeasible )
               {
                  SCIPinfoMessage(scip, NULL, "Infeasibility detected because of integrality of variable %s\n",
                        SCIPvarGetName(data->createdvars[v]));
               }
            }
            else
            {
               SCIPerrorMessage("Number of integrality constraints %d should be non-negative!\n", nintvars);
               SCIPABORT();
            }
         }
      }
      else
      {
         SCIPerrorMessage("Number of integrality constraints %d should be non-negative!\n", nintvars);
         SCIPABORT();
      }
   }
   else
      return SCIP_READERROR;

   return SCIP_OKAY;
}

/** reads SDP-constraint sizes from given CBF-file */
static
SCIP_RETCODE CBFreadPsdcon(
   SCIP*                 scip,               /* SCIP data structure */
   SCIP_FILE*            pfile,              /* file to read from */
   long long int *       linecount,          /* current linecount */
   CBF_DATA*             data                /* data pointer to save the results in */
   )
{
   int b;

   assert( scip != NULL );
   assert( pfile != NULL );
   assert( linecount != NULL );
   assert( data != NULL );

   SCIP_CALL( CBFfgets(pfile, linecount) );

   if ( sscanf(CBF_LINE_BUFFER, "%i", &(data->nsdpblocks)) == 1 )
   {
      if ( data->nsdpblocks >= 0 )
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &(data->sdpblocksizes), data->nsdpblocks) );

         for (b = 0; b < data->nsdpblocks; b++)
         {
            SCIP_CALL( CBFfgets(pfile, linecount) );
            if ( sscanf(CBF_LINE_BUFFER, "%i", &(data->sdpblocksizes[b])) == 1 )
            {
               if ( data->sdpblocksizes[b] <= 0 )
               {
                  SCIPerrorMessage("Size %d of SDP-block %d should be positive!\n", data->sdpblocksizes[b], b);
                  SCIPABORT();
               }
            }
            else
               return SCIP_READERROR;
         }
      }
      else
      {
         SCIPerrorMessage("Number of SDP-blocks %d should be non-negative!\n", data->nsdpblocks);
         SCIPABORT();
      }
   }
   else
      return SCIP_READERROR;

   return SCIP_OKAY;
}

/** reads objective values from given CBF-file */
static
SCIP_RETCODE CBFreadObjacoord(
   SCIP*                 scip,               /* SCIP data structure */
   SCIP_FILE*            pfile,              /* file to read from */
   long long int *       linecount,          /* current linecount */
   CBF_DATA*             data                /* data pointer to save the results in */
   )
{  /*lint --e{818}*/
   int nobjcoefs;
   int i;
   int v;
   SCIP_Real val;

   assert( scip != NULL );
   assert( pfile != NULL );
   assert( linecount != NULL );
   assert( data != NULL );

   SCIP_CALL( CBFfgets(pfile, linecount) );

   if ( sscanf(CBF_LINE_BUFFER, "%i", &nobjcoefs) == 1 )
   {
      if ( nobjcoefs >= 0 )
      {
         for (i = 0; i < nobjcoefs; i++)
         {
            SCIP_CALL( CBFfgets(pfile, linecount) );
            if ( sscanf(CBF_LINE_BUFFER, "%i %lf", &v, &val) == 2 )
            {
               if ( v < 0 || v >= data->nvars )
               {
                  SCIPerrorMessage("Given objective coefficient for variable %d which does not exist!\n", v);
                  SCIPABORT();
               }

               if ( SCIPisZero(scip, val) )
               {
                  SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Ignored objective coefficient of variable %d; value"
                        "%f is smaller than epsilon = %f.\n", v, val, SCIPepsilon(scip));
               }
               else
               {
                  SCIP_CALL( SCIPchgVarObj(scip, data->createdvars[v], val) );
               }
            }
            else
               return SCIP_READERROR;
         }
      }
      else
      {
         SCIPerrorMessage("Number of objective coefficients %d should be non-negative!\n", nobjcoefs);
         SCIPABORT();
      }
   }
   else
      return SCIP_READERROR;

   return SCIP_OKAY;
}

/** reads linear coefficients from given CBF-file */
static
SCIP_RETCODE CBFreadAcoord(
   SCIP*                 scip,               /* SCIP data structure */
   SCIP_FILE*            pfile,              /* file to read from */
   long long int *       linecount,          /* current linecount */
   CBF_DATA*             data                /* data pointer to save the results in */
   )
{  /*lint --e{818}*/
   int ncoefs;
   int c;
   int i;
   int v;
   SCIP_Real val;

   assert( scip != NULL );
   assert( pfile != NULL );
   assert( linecount != NULL );
   assert( data != NULL );

   SCIP_CALL( CBFfgets(pfile, linecount) );

   if ( sscanf(CBF_LINE_BUFFER, "%i", &ncoefs) == 1 )
   {
      if ( ncoefs >= 0 )
      {
         for (i = 0; i < ncoefs; i++)
         {
            SCIP_CALL( CBFfgets(pfile, linecount) );
            if ( sscanf(CBF_LINE_BUFFER, "%i %i %lf", &c, &v, &val) == 3 )
            {
               if ( c < 0 || c >= data->nconss )
               {
                  SCIPerrorMessage("Given linear coefficient for constraint %d which does not exist!\n", c);
                  SCIPABORT();
               }
               if ( v < 0 || v >= data->nvars )
               {
                  SCIPerrorMessage("Given linear coefficient for variable %d which does not exist!\n", v);
                  SCIPABORT();
               }
               if ( SCIPisZero(scip, val) )
               {
                  SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Ignored linear coefficient of constraint %d, variable "
                        "%d; value %.9f is smaller than epsilon = %f.\n", c, v, val, SCIPepsilon(scip));
               }
               else
               {
                  SCIP_CALL( SCIPaddCoefLinear(scip, data->createdconss[c], data->createdvars[v], val) );/*lint !e732*//*lint !e747*/
               }
            }
            else
               return SCIP_READERROR;
         }
      }
      else
      {
         SCIPerrorMessage("Number of linear coefficients %d should be non-negative!\n", ncoefs);
         SCIPABORT();
      }
   }
   else
      return SCIP_READERROR;

   return SCIP_OKAY;
}

/** reads left- and right-hand sides from given CBF-file */
static
SCIP_RETCODE CBFreadBcoord(
   SCIP*                 scip,               /* SCIP data structure */
   SCIP_FILE*            pfile,              /* file to read from */
   long long int *       linecount,          /* current linecount */
   CBF_DATA*             data                /* data pointer to save the results in */
   )
{  /*lint --e{818}*/
   int nsides;
   int c;
   int i;
   SCIP_Real val;

   assert( scip != NULL );
   assert( pfile != NULL );
   assert( linecount != NULL );
   assert( data != NULL );

   SCIP_CALL( CBFfgets(pfile, linecount) );

   if ( sscanf(CBF_LINE_BUFFER, "%i", &nsides) == 1 )
   {
      if ( nsides >= 0 )
      {
         for (i = 0; i < nsides; i++)
         {
            SCIP_CALL( CBFfgets(pfile, linecount) );
            if ( sscanf(CBF_LINE_BUFFER, "%i %lf", &c, &val) == 2 )
            {
               if ( c < 0 || c >= data->nconss )
               {
                  SCIPerrorMessage("Given constant part for constraint %d which does not exist!\n", c);
                  SCIPABORT();
               }
               if ( SCIPisZero(scip, val) )
               {
                  SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Ignored constant part of constraint %d; value %.9f is "
                        "smaller than epsilon = %f.\n", c, val, SCIPepsilon(scip));
               }
               else
               {
                  if ( c < data->firstleqcons )
                  {
                     /* greater or equal constraint -> left-hand side (minus since we have Ax + b >= 0) */
                     SCIP_CALL( SCIPchgLhsLinear(scip, data->createdconss[c], -val) );
                  }
                  else if ( c < data->firsteqcons )
                  {
                     /* less or equal constraint -> right-hand side (minus since we have Ax + b <= 0) */
                     SCIP_CALL( SCIPchgRhsLinear(scip, data->createdconss[c], -val) );
                  }
                  else
                  {
                     /* equality constraint -> left- and right-hand side (minus since we have Ax + b = 0) */
                     SCIP_CALL( SCIPchgLhsLinear(scip, data->createdconss[c], -val) );
                     SCIP_CALL( SCIPchgRhsLinear(scip, data->createdconss[c], -val) );
                  }
               }
            }
            else
               return SCIP_READERROR;
         }
      }
      else
      {
         SCIPerrorMessage("Number of left- and right-hand sides %d should be non-negative!\n", nsides);
         SCIPABORT();
      }
   }
   else
      return SCIP_READERROR;

   return SCIP_OKAY;
}

/** reads nonzero coefficients of SDP-constraints from given CBF-file */
static
SCIP_RETCODE CBFreadHcoord(
   SCIP*                 scip,               /* SCIP data structure */
   SCIP_FILE*            pfile,              /* file to read from */
   long long int *       linecount,          /* current linecount */
   CBF_DATA*             data                /* data pointer to save the results in */
   )
{
   int nnonz;
   int i;
   int b;
   int v;
   int row;
   int col;
   SCIP_Real val;
   int** sdpvar;
   int firstindforvar;
   int nextindaftervar;
   SCIP_Bool varused;

   assert( scip != NULL );
   assert( pfile != NULL );
   assert( linecount != NULL );
   assert( data != NULL );

   /* initialize sdpnblocknonz with 0 */
   SCIP_CALL( SCIPallocBufferArray(scip, &(data->sdpnblocknonz), data->nsdpblocks) );
   for (b = 0; b < data->nsdpblocks; b++)
      data->sdpnblocknonz[b] = 0;

   SCIP_CALL( CBFfgets(pfile, linecount) );

   if ( sscanf(CBF_LINE_BUFFER, "%i", &nnonz) == 1 )
   {
      if ( nnonz >= 0 )
      {
         /* allocate memory (nnonz for each block, since we do not yet no the distribution) */
         SCIP_CALL( SCIPallocBufferArray(scip, &sdpvar, data->nsdpblocks) );
         SCIP_CALL( SCIPallocBufferArray(scip, &(data->sdprow), data->nsdpblocks) );
         SCIP_CALL( SCIPallocBufferArray(scip, &(data->sdpcol), data->nsdpblocks) );
         SCIP_CALL( SCIPallocBufferArray(scip, &(data->sdpval), data->nsdpblocks) );
         for (b = 0; b < data->nsdpblocks; b++)
         {
            SCIP_CALL( SCIPallocBufferArray(scip, &(sdpvar[b]), nnonz) );
            SCIP_CALL( SCIPallocBufferArray(scip, &(data->sdprow[b]), nnonz) );
            SCIP_CALL( SCIPallocBufferArray(scip, &(data->sdpcol[b]), nnonz) );
            SCIP_CALL( SCIPallocBufferArray(scip, &(data->sdpval[b]), nnonz) );
         }

         for (i = 0; i < nnonz; i++)
         {
            SCIP_CALL( CBFfgets(pfile, linecount) );
            if ( sscanf(CBF_LINE_BUFFER, "%i %i %i %i %lf", &b, &v, &row, &col, &val) == 5 )
            {
               if ( b < 0 || b >= data->nsdpblocks )
               {
                  SCIPerrorMessage("Given SDP-coefficient for SDP-constraint %d which does not exist!\n", b);
                  SCIPABORT();
               }
               if ( v < 0 || v >= data->nvars )
               {
                  SCIPerrorMessage("Given SDP-coefficient for variable %d which does not exist!\n", v);
                  SCIPABORT();
               }
               if ( row < 0 || row >= data->sdpblocksizes[b] )
               {
                  SCIPerrorMessage("Row index %d of given SDP coefficient is negative or larger than blocksize %d!\n",
                        row, data->sdpblocksizes[b]);
                  SCIPABORT();
               }
               if ( col < 0 || col >= data->sdpblocksizes[b] )
               {
                  SCIPerrorMessage("Column index %d of given SDP coefficient is negative or larger than blocksize %d!\n",
                        col, data->sdpblocksizes[b]);
                  SCIPABORT();
               }
               if ( SCIPisZero(scip, val) )
               {
                  SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Ignored SDP-coefficient of SDP-constraint %d, variable "
                        "%d, row %d, col %d; value %.9f is smaller than epsilon = %f.\n", b, v, row, col, val, SCIPepsilon(scip));
               }
               else
               {
                  sdpvar[b][data->sdpnblocknonz[b]] = v;
                  data->sdprow[b][data->sdpnblocknonz[b]] = row;
                  data->sdpcol[b][data->sdpnblocknonz[b]] = col;
                  data->sdpval[b][data->sdpnblocknonz[b]] = val;
                  data->sdpnblocknonz[b]++;
               }
            }
            else
               return SCIP_READERROR;
         }

         /* construct pointer arrays */
         SCIP_CALL( SCIPallocBufferArray(scip, &(data->sdpnblockvars), data->nsdpblocks) );
         SCIP_CALL( SCIPallocBufferArray(scip, &(data->sdpblockvars), data->nsdpblocks) );
         SCIP_CALL( SCIPallocBufferArray(scip, &(data->nvarnonz), data->nsdpblocks) );
         SCIP_CALL( SCIPallocBufferArray(scip, &(data->rowpointer), data->nsdpblocks) );
         SCIP_CALL( SCIPallocBufferArray(scip, &(data->colpointer), data->nsdpblocks) );
         SCIP_CALL( SCIPallocBufferArray(scip, &(data->valpointer), data->nsdpblocks) );
         for (b = 0; b < data->nsdpblocks; b++)
         {
            /* sort the nonzeroes by non-decreasing variable indices */
            SCIPsortIntIntIntReal(sdpvar[b], data->sdprow[b], data->sdpcol[b], data->sdpval[b], data->sdpnblocknonz[b]);

            /* create the pointer arrays and insert used variables into vars-array */
            nextindaftervar = 0;
            data->sdpnblockvars[b] = 0;
            SCIP_CALL( SCIPallocBufferArray(scip, &(data->sdpblockvars[b]), data->nvars) );
            SCIP_CALL( SCIPallocBufferArray(scip, &(data->nvarnonz[b]), data->nvars) );
            SCIP_CALL( SCIPallocBufferArray(scip, &(data->rowpointer[b]), data->nvars) );
            SCIP_CALL( SCIPallocBufferArray(scip, &(data->colpointer[b]), data->nvars) );
            SCIP_CALL( SCIPallocBufferArray(scip, &(data->valpointer[b]), data->nvars) );

            for (v = 0; v < data->nvars; v++)
            {
               varused = FALSE;
               firstindforvar = nextindaftervar; /* this variable starts where the last one ended */
               data->nvarnonz[b][data->sdpnblockvars[b]] = 0;
               while (nextindaftervar < data->sdpnblocknonz[b] && sdpvar[b][nextindaftervar] == v) /* get the first index that doesn't belong to this variable */
               {
                  nextindaftervar++;
                  varused = TRUE;
                  data->nvarnonz[b][data->sdpnblockvars[b]]++;
               }
               if (varused)
               {
                  data->sdpblockvars[b][data->sdpnblockvars[b]] = data->createdvars[v];/*lint !e732*//*lint !e747*/ /* if the variable is used, add it to the vars array */
                  data->rowpointer[b][data->sdpnblockvars[b]] = &(data->sdprow[b][firstindforvar]); /* save a pointer to the first nonzero belonging to this variable */
                  data->colpointer[b][data->sdpnblockvars[b]] = &(data->sdpcol[b][firstindforvar]);
                  data->valpointer[b][data->sdpnblockvars[b]] = &(data->sdpval[b][firstindforvar]);
                  data->sdpnblockvars[b]++;
               }
            }

            assert (nextindaftervar == data->sdpnblocknonz[b]);
#ifdef SCIP_DISABLED_CODE
            /* resize arrays and free sdpvar array which is no longer needed */
            SCIP_CALL( SCIPreallocBufferArray(scip, &(data->rowpointer[b]), data->sdpnblocknonz[b]) );
            SCIP_CALL( SCIPreallocBufferArray(scip, &(data->colpointer[b]), data->sdpnblocknonz[b]) );
            SCIP_CALL( SCIPreallocBufferArray(scip, &(data->valpointer[b]), data->sdpnblocknonz[b]) );
#endif
            SCIPfreeBufferArray(scip, &(sdpvar[b]));
         }

         /* free SDP-var array which is no longer needed */
         SCIPfreeBufferArray(scip, &sdpvar);
      }
      else
      {
         SCIPerrorMessage("Number of nonzero coefficients of SDP-constraints %d should be non-negative!\n", nnonz);
         SCIPABORT();
      }
   }
   else
      return SCIP_READERROR;

   return SCIP_OKAY;
}

/** reads constant entries of SDP-constraints from given CBF-file */
static
SCIP_RETCODE CBFreadDcoord(
   SCIP*                 scip,               /* SCIP data structure */
   SCIP_FILE*            pfile,              /* file to read from */
   long long int *       linecount,          /* current linecount */
   CBF_DATA*             data                /* data pointer to save the results in */
   )
{
   int constnnonz;
   int b;
   int i;
   int row;
   int col;
   SCIP_Real val;

   assert( scip != NULL );
   assert( pfile != NULL );
   assert( linecount != NULL );
   assert( data != NULL );

   /* initialize sdpconstnblocknonz with 0 */
   SCIP_CALL( SCIPallocBufferArray(scip, &(data->sdpconstnblocknonz), data->nsdpblocks) );
   for (b = 0; b < data->nsdpblocks; b++)
      data->sdpconstnblocknonz[b] = 0;

   SCIP_CALL( CBFfgets(pfile, linecount) );

   if ( sscanf(CBF_LINE_BUFFER, "%i", &constnnonz) == 1 )
   {
      if ( constnnonz >= 0 )
      {
         /* allocate memory (constnnonz for each block, since we do not yet no the distribution) */
         SCIP_CALL( SCIPallocBufferArray(scip, &(data->sdpconstrow), data->nsdpblocks) );
         SCIP_CALL( SCIPallocBufferArray(scip, &(data->sdpconstcol), data->nsdpblocks) );
         SCIP_CALL( SCIPallocBufferArray(scip, &(data->sdpconstval), data->nsdpblocks) );
         for (b = 0; b < data->nsdpblocks; b++)
         {
            SCIP_CALL( SCIPallocBufferArray(scip, &(data->sdpconstrow[b]), constnnonz) );
            SCIP_CALL( SCIPallocBufferArray(scip, &(data->sdpconstcol[b]), constnnonz) );
            SCIP_CALL( SCIPallocBufferArray(scip, &(data->sdpconstval[b]), constnnonz) );
         }

         for (i = 0; i < constnnonz; i++)
         {
            SCIP_CALL( CBFfgets(pfile, linecount) );
            if ( sscanf(CBF_LINE_BUFFER, "%i %i %i %lf", &b, &row, &col, &val) == 4 )
            {
               if ( b < 0 || b >= data->nsdpblocks )
               {
                  SCIPerrorMessage("Given constant entry for SDP-constraint %d which does not exist!\n", b);
                  SCIPABORT();
               }
               if ( row < 0 || row >= data->sdpblocksizes[b] )
               {
                  SCIPerrorMessage("Row index %d of given constant SDP-entry is negative or larger than blocksize %d!\n",
                        row, data->sdpblocksizes[b]);
                  SCIPABORT();
               }
               if ( col < 0 || col >= data->sdpblocksizes[b] )
               {
                  SCIPerrorMessage("Column index %d of given constant SDP-entry is negative or larger than blocksize %d!\n",
                        col, data->sdpblocksizes[b]);
                  SCIPABORT();
               }
               if ( SCIPisZero(scip, val) )
               {
                  SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Ignored constant entry of SDP-constraint %d, row %d,"
                        " col %d; value %.9f is smaller than epsilon = %f.\n", b, row, col, val, SCIPepsilon(scip));
               }
               else
               {
                  data->sdpconstrow[b][data->sdpconstnblocknonz[b]] = row;
                  data->sdpconstcol[b][data->sdpconstnblocknonz[b]] = col;
                  data->sdpconstval[b][data->sdpconstnblocknonz[b]] = -val;
                  data->sdpconstnblocknonz[b]++;
               }
            }
            else
               return SCIP_READERROR;
         }
#ifdef SCIP_DISABLED_CODE
         for (b = 0; b < data->nsdpblocks; b++)
         {
            /* resize arrays */
            SCIP_CALL( SCIPreallocBufferArray(scip, &(data->sdpconstrow[b]), data->sdpconstnblocknonz[b]) );
            SCIP_CALL( SCIPreallocBufferArray(scip, &(data->sdpconstcol[b]), data->sdpconstnblocknonz[b]) );
            SCIP_CALL( SCIPreallocBufferArray(scip, &(data->sdpconstval[b]), data->sdpconstnblocknonz[b]) );
         }
#endif
      }
      else
      {
         SCIPerrorMessage("Number of constant entries of SDP-constraints %d should be non-negative!\n", constnnonz);
         SCIPABORT();
      }
   }
   else
      return SCIP_READERROR;

   return SCIP_OKAY;
}

/** frees all data allocated for the CBF-data-struct */
static
SCIP_RETCODE CBFfreeData(
   SCIP*                 scip,               /* SCIP data structure */
   CBF_DATA*             data                /* data pointer to save the results in */
   )
{
   int b;

   assert( scip != NULL );
   assert( data != NULL );

   for (b = 0; b < data->nsdpblocks; b++)
   {
      SCIPfreeBufferArrayNull(scip, &(data->sdpconstval[b]));
      SCIPfreeBufferArrayNull(scip, &(data->sdpconstcol[b]));
      SCIPfreeBufferArrayNull(scip, &(data->sdpconstrow[b]));
      SCIPfreeBufferArrayNull(scip, &(data->valpointer[b]));
      SCIPfreeBufferArrayNull(scip, &(data->colpointer[b]));
      SCIPfreeBufferArrayNull(scip, &(data->rowpointer[b]));
      SCIPfreeBufferArrayNull(scip, &(data->sdpval[b]));
      SCIPfreeBufferArrayNull(scip, &(data->sdpcol[b]));
      SCIPfreeBufferArrayNull(scip, &(data->sdprow[b]));
      SCIPfreeBufferArrayNull(scip, &(data->sdpblockvars[b]));
      SCIPfreeBufferArrayNull(scip, &(data->nvarnonz[b]));
   }
   SCIPfreeBufferArrayNull(scip, &data->sdpconstval);
   SCIPfreeBufferArrayNull(scip, &data->sdpconstcol);
   SCIPfreeBufferArrayNull(scip, &data->sdpconstrow);
   SCIPfreeBufferArrayNull(scip, &data->sdpconstnblocknonz);
   SCIPfreeBufferArrayNull(scip, &data->valpointer);
   SCIPfreeBufferArrayNull(scip, &data->colpointer);
   SCIPfreeBufferArrayNull(scip, &data->rowpointer);
   SCIPfreeBufferArrayNull(scip, &data->sdpval);
   SCIPfreeBufferArrayNull(scip, &data->sdpcol);
   SCIPfreeBufferArrayNull(scip, &data->sdprow);
   SCIPfreeBufferArrayNull(scip, &data->sdpblockvars);
   SCIPfreeBufferArrayNull(scip, &data->nvarnonz);
   SCIPfreeBufferArrayNull(scip, &data->sdpnblockvars);
   SCIPfreeBufferArrayNull(scip, &data->sdpnblocknonz);
   SCIPfreeBufferArrayNull(scip, &data->sdpblocksizes);
   SCIPfreeBufferArrayNull(scip, &data->createdconss);
   SCIPfreeBufferArrayNull(scip, &data->createdvars);

   return SCIP_OKAY;
}


/*
 * Callback methods of reader
 */


/** copy method for reader plugins (called when SCIP copies plugins) */
static
SCIP_DECL_READERCOPY(readerCopyCbf)
{  /*lint --e{715,818}*/
   assert(scip != NULL);

   SCIP_CALL( SCIPincludeReaderCbf(scip) );

   return SCIP_OKAY;
}

/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadCbf)
{  /*lint --e{715,818}*/
   SCIP_FILE* scipfile;
   long long int linecount;
   SCIP_Bool versionread;
   SCIP_Bool objread;
   CBF_DATA* data;
   int b;

   *result = SCIP_DIDNOTRUN;
   linecount = 0;
   versionread = FALSE;
   objread = FALSE;

   SCIPdebugMsg(scip, "Reading file %s\n", filename);

   scipfile = SCIPfopen(filename, "r");

   if ( ! scipfile )
      return SCIP_READERROR;

   SCIP_CALL( SCIPallocBuffer(scip, &data) );
   data->nsdpblocks = 0;

   /* create empty problem */
   SCIP_CALL( SCIPcreateProb(scip, filename, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

   while( CBFfgets(scipfile, &linecount) == SCIP_OKAY )
   {
      /* Parse keyword on non-empty lines */
      if ( sscanf(CBF_LINE_BUFFER, CBF_NAME_FORMAT, CBF_NAME_BUFFER) == 1 )
      {
         /* first line should be version number */
         if ( ! versionread )
         {
            if (strcmp(CBF_NAME_BUFFER, "VER") == 0)
            {
               int ver;

               SCIP_CALL( CBFfgets(scipfile, &linecount) );

               if ( sscanf(CBF_LINE_BUFFER, "%i", &ver) == 1 )
               {
                  SCIPdebugMsg(scip, "file version %d\n", ver);
                  if ( ver != CBF_VERSION_NR )
                  {
                     SCIPerrorMessage("Only version number %d is supported!\n", CBF_VERSION_NR);
                     SCIPABORT();
                  }
                  else
                     versionread = TRUE;
               }
               else
                  return SCIP_READERROR;
            }
            else
            {
               SCIPerrorMessage("First keyword should be VER.\n");
               SCIPABORT();
            }
         }
         else
         {
            if ( strcmp(CBF_NAME_BUFFER, "OBJSENSE") == 0 )
            {
               SCIPdebugMsg(scip, "Reading OBJSENSE\n");
               SCIP_CALL( CBFreadObjsense(scip, scipfile, &linecount) );
               objread = TRUE;
            }
            else if (strcmp(CBF_NAME_BUFFER, "VAR") == 0)
            {
               SCIPdebugMsg(scip, "Reading VAR\n");
               SCIP_CALL( CBFreadVar(scip, scipfile, &linecount, data) );
            }
            else if (strcmp(CBF_NAME_BUFFER, "CON") == 0)
            {
               SCIPdebugMsg(scip, "Reading CON\n");
               SCIP_CALL( CBFreadCon(scip, scipfile, &linecount, data) );
            }
            else if (strcmp(CBF_NAME_BUFFER, "INT") == 0)
            {
               SCIPdebugMsg(scip, "Reading INT\n");
               SCIP_CALL( CBFreadInt(scip, scipfile, &linecount, data) );
            }
            else if (strcmp(CBF_NAME_BUFFER, "PSDCON") == 0)
            {
               SCIPdebugMsg(scip, "Reading PSDCON\n");
               SCIP_CALL( CBFreadPsdcon(scip, scipfile, &linecount, data) );
            }
            else if (strcmp(CBF_NAME_BUFFER, "PSDVAR") == 0)
            {
               /* TODO: automatically transform to dual form by introducing a new variable for all elements in the upper
                * triangular part and an SDP-constraint enforcing positive semidefiniteness of the PSD-variable
                */
               SCIPerrorMessage("SDPs in primal form currently not supported, please use PSDCON!\n");
               SCIPABORT();
            }
            else if (strcmp(CBF_NAME_BUFFER, "OBJFCOORD") == 0)
            {
               SCIPerrorMessage("SDPs in primal form currently not supported, please use PSDCON!\n");
               SCIPABORT();
            }
            else if (strcmp(CBF_NAME_BUFFER, "OBJACOORD") == 0)
            {
               SCIPdebugMsg(scip, "Reading OBJACOORD\n");
               SCIP_CALL( CBFreadObjacoord(scip, scipfile, &linecount, data) );
            }
            else if (strcmp(CBF_NAME_BUFFER, "OBJBCOORD") == 0)
            {
               SCIPerrorMessage("constant part in objective value not supported by SCIP!\n");
               SCIPABORT();
            }
            else if (strcmp(CBF_NAME_BUFFER, "FCOORD") == 0)
            {
               SCIPerrorMessage("SDPs in primal form currently not supported, please use PSDCON!\n");
               SCIPABORT();
            }
            else if (strcmp(CBF_NAME_BUFFER, "ACOORD") == 0)
            {
               SCIPdebugMsg(scip, "Reading ACOORD\n");
               SCIP_CALL( CBFreadAcoord(scip, scipfile, &linecount, data) );
            }
            else if (strcmp(CBF_NAME_BUFFER, "BCOORD") == 0)
            {
               SCIPdebugMsg(scip, "Reading BCOORD\n");
               SCIP_CALL( CBFreadBcoord(scip, scipfile, &linecount, data) );
            }
            else if (strcmp(CBF_NAME_BUFFER, "HCOORD") == 0)
            {
               SCIPdebugMsg(scip, "Reading HCOORD\n");
               SCIP_CALL( CBFreadHcoord(scip, scipfile, &linecount, data) );
            }
            else if (strcmp(CBF_NAME_BUFFER, "DCOORD") == 0)
            {
               SCIPdebugMsg(scip, "Reading DCOORD\n");
               SCIP_CALL( CBFreadDcoord(scip, scipfile, &linecount, data) );
            }
            else
            {
               SCIPerrorMessage("Keyword %s not recognized!\n", CBF_NAME_BUFFER);
               SCIPABORT();
            }
         }
      }
   }

   if ( !objread )
   {
      SCIPerrorMessage("Keyword OBJSENSE is missing!\n");
      SCIPABORT();
   }

   /* close the file (and make sure SCIPfclose returns 0) */
   if ( SCIPfclose(scipfile) )
      return SCIP_READERROR;

#ifdef SCIP_MORE_DEBUG
   for (b = 0; b < SCIPgetNConss(scip); b++)
   {
      SCIP_CALL( SCIPprintCons(scip, SCIPgetConss(scip)[b], NULL) );
   }
#endif

   /* create SDP-constraints */
   for (b = 0; b < data->nsdpblocks; b++)
   {
      SCIP_CONS* sdpcons;
      char sdpconname[SCIP_MAXSTRLEN];
#ifndef NDEBUG
      int snprintfreturn;
#endif

      assert( data->sdpblocksizes[b] > 0 );
      assert( ((data->sdpnblockvars[b] > 0) && data->sdpnblocknonz[b] > 0) || (data->sdpconstnblocknonz[b] > 0) );

#ifndef NDEBUG
      snprintfreturn = SCIPsnprintf(sdpconname, SCIP_MAXSTRLEN, "SDP_%d", b);
      assert( snprintfreturn < SCIP_MAXSTRLEN);
#else
      (void) SCIPsnprintf(sdpconname, SCIP_MAXSTRLEN, "SDP_%d", b);
#endif
      SCIP_CALL( SCIPcreateConsSdp(scip, &sdpcons, sdpconname, data->sdpnblockvars[b], data->sdpnblocknonz[b],
            data->sdpblocksizes[b], data->nvarnonz[b], data->colpointer[b], data->rowpointer[b], data->valpointer[b],
            data->sdpblockvars[b], data->sdpconstnblocknonz[b], data->sdpconstcol[b], data->sdpconstrow[b],
            data->sdpconstval[b]) );

#ifdef SCIP_MORE_DEBUG
SCIP_CALL( SCIPprintCons(scip, sdpcons, NULL) );
#endif

      SCIP_CALL( SCIPaddCons(scip, sdpcons) );

      SCIP_CALL( SCIPreleaseCons(scip, &sdpcons) );
   }

   SCIP_CALL( CBFfreeData(scip, data) );
   SCIPfreeBufferNull(scip, &data);

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/** problem writing method of reader */
static
SCIP_DECL_READERWRITE(readerWriteCbf)
{  /*lint --e{715,818}*/
   SdpVarmapper* varmapper;
   SCIP_VAR** linvars;
   SCIP_Real* linvals;
   SCIP_VAR** sdpvars;
   int nsdpconss;
   int sdpnvars;
   int sdpnnonz;
   int totalsdpnnonz;
   int sdpblocksize;
   int sdparraylength;
   int* sdpnvarnonz;
   int** sdpcol;
   int** sdprow;
   SCIP_Real** sdpval;
   int totalsdpconstnnonz;
   int sdpconstnnonz;
   int* sdpconstcol;
   int* sdpconstrow;
   SCIP_Real* sdpconstval;
   int c;
   int i;
   int v;
   int consind;
   int neqconss;
   int nleqconss;
   int ngeqconss;
   int nobjnonz;
   int nnonz;
   int nbnonz;
#ifdef CBF_CHECK_NONNEG
   SCIP_Bool* consdisabled;
   int nposorth;
   SCIP_Bool* posorth;
   int nnegorth;
   SCIP_Bool* negorth;
#endif

   assert( scip != NULL );

   SCIPdebugMsg(scip, "Writing problem in CBF format to %s\n", file);
   *result = SCIP_DIDNOTRUN;

   if ( transformed )
   {
      SCIPerrorMessage("CBF reader currently only supports writing original problems!\n");
      SCIPABORT(); /*lint --e{527}*/
   }

   for (c = 0; c < nconss; c++)
   {
      if ( (strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "linear") != 0)
            && (strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "SDP") != 0 ) )
      {
         SCIPerrorMessage("CBF reader currently only supports linear and SDP constraints!\n");
         SCIPABORT(); /*lint --e{527}*/
      }
   }

#ifndef NDEBUG
   for (v = 0; v < nvars; v++)
   {
      assert( SCIPvarGetStatus(vars[v]) == SCIP_VARSTATUS_ORIGINAL );
   }
#endif

   /* write version number */
   SCIPinfoMessage(scip, file, "VER\n%d\n\n", CBF_VERSION_NR);

   /* write objective sense */
   SCIPinfoMessage(scip, file, "OBJSENSE\n%s\n\n", objsense == SCIP_OBJSENSE_MINIMIZE ? "MIN" : "MAX");

#ifdef CBF_CHECK_NONNEG
   /* check for constraints with only a single non-zero and lhs/rhs zero, in that case disable the constraint
    * and force variable to be in positive or negative orthant
    */
   SCIP_CALL( SCIPallocBufferArray(scip, &consdisabled, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &posorth, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &negorth, nvars) );

   for (v = 0; v < nvars; v++)
   {
      posorth[v] = FALSE;
      negorth[v] = FALSE;
   }
   nposorth = 0;
   nnegorth = 0;

   for (c = 0; c < nconss; c++)
   {
      consdisabled[c] = 0;

      /* only check linear constraints */
      if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "linear") != 0 )
         continue;

      /* only check constraints with exactly one nonzero */
      if ( SCIPgetNVarsLinear(scip, conss[c]) != 1 )
         continue;

      /* the nonzero should be a true nonzero */
      assert( ! SCIPisZero(scip, SCIPgetValsLinear(scip, conss[c])[0]) );

      /* the variable should not be fixed */
      assert( SCIPisLT(scip, SCIPgetLhsLinear(scip, conss[c]), SCIPgetRhsLinear(scip, conss[c])) );

      if ( SCIPgetValsLinear(scip, conss[c])[0] > 0.0 )
      {
         if ( SCIPisZero(scip, SCIPgetLhsLinear(scip, conss[c])) )
         {
            consdisabled[c] = 1;
            posorth[v] = TRUE;
            nposorth++;

            if ( ! SCIPisInfinity(scip, SCIPgetRhsLinear(scip, conss[c])) )
            {
               SCIPerrorMessage("Detection of non-negativity constraints currently not supported for ranged rows!\n");
               SCIPABORT(); /*lint --e{527}*/
            }
         }
         else if ( SCIPisZero(scip, SCIPgetRhsLinear(scip, conss[c])) )
         {
            consdisabled[c] = 1;
            negorth[v] = TRUE;
            nnegorth++;
         }
      }
      else
      {
         if ( SCIPisZero(scip, SCIPgetLhsLinear(scip, conss[c])) )
         {
            consdisabled[c] = 1;
            negorth[v] = TRUE;
            nnegorth++;
         }
         else if ( SCIPisZero(scip, SCIPgetRhsLinear(scip, conss[c])) )
         {
            consdisabled[c] = 1;
            posorth[v] = TRUE;
            nposorth++;
         }
      }
   }

   assert( (0 <= nposorth) && (0 <= nnegorth) && (nposorth + nnegorth <= nvars) );

   /* prepare varmapper that maps SCIP variables to indices for CBF format (3/4 is the load factor java uses) */
   SCIP_CALL( SCIPsdpVarmapperCreate(scip, &(varmapper), (int) ceil(1.33 * nvars)) );
   /* in this case we first need to add all variables in the positive orthant, then all in the negative orthant
    * and finally the free variables, since they will be input and counted in CBF this way
    */
   for (v = 0; v < nvars; v++)
   {
      if ( posorth[v] )
      {
         SCIP_CALL( SCIPsdpVarmapperAddVars(scip, varmapper, 1, &(vars[v])) );
      }
   }
   for (v = 0; v < nvars; v++)
   {
      if ( negorth[v] )
      {
         SCIP_CALL( SCIPsdpVarmapperAddVars(scip, varmapper, 1, &(vars[v])) );
      }
   }
   for (v = 0; v < nvars; v++)
   {
      if ( ( ! posorth[v] ) && ( ! negorth[v] ) )
      {
         SCIP_CALL( SCIPsdpVarmapperAddVars(scip, varmapper, 1, &(vars[v])) );
      }
   }

   /* write variables */
   if ( (nposorth == 0) && (nnegorth == 0) )
      SCIPinfoMessage(scip, file, "VAR\n%d 1\nF %d\n\n", nvars, nvars);
   else if ( (nposorth > 0) && (nnegorth == 0) && (nposorth  < nvars) )
      SCIPinfoMessage(scip, file, "VAR\n%d 2\nL+ %d\nF %d\n\n", nvars, nposorth, nvars - nposorth);
   else if ( (nposorth == 0) && (nnegorth > 0) && (nnegorth < nvars) )
      SCIPinfoMessage(scip, file, "VAR\n%d 2\nL- %d\nF %d\n\n", nvars, nposorth, nvars - nnegorth);
   else if ( nposorth + nnegorth < nvars )
      SCIPinfoMessage(scip, file, "VAR\n%d 3\nL+ %d\nL- %d\nF %d\n\n", nvars, nposorth, nnegorth, nvars - nposorth - nnegorth);/*lint !e834*/
   else if ( (nposorth > 0) && (nnegorth == 0) )
      SCIPinfoMessage(scip, file, "VAR\n%d 1\nL+ %d\n\n", nvars, nposorth);
   else if ( (nposorth == 0) && (nnegorth > 0) )
      SCIPinfoMessage(scip, file, "VAR\n%d 1\nL- %d\n\n", nvars, nnegorth);
   else
      SCIPinfoMessage(scip, file, "VAR\n%d 2\nL+ %d\nL- %d\n\n", nvars, nposorth, nnegorth);

#else
   /* prepare varmapper that maps SCIP variables to indices for CBF format (3/4 is the load factor java uses) */
   SCIP_CALL( SCIPsdpVarmapperCreate(scip, &(varmapper), (int) ceil(1.33 * nvars)) );
   SCIP_CALL( SCIPsdpVarmapperAddVars(scip, varmapper, nvars, vars) );

   /* write variables */
   SCIPinfoMessage(scip, file, "VAR\n%d 1\nF %d\n\n", nvars, nvars);
#endif

   /* write integrality constraints */
   if ( nbinvars + nintvars > 0 )
   {
      SCIPinfoMessage(scip, file, "INT\n%d\n", nbinvars + nintvars);

      for (v = 0; v < nbinvars + nintvars; v++)
      {
         assert( SCIPvarIsIntegral(vars[v]) );
         SCIPinfoMessage(scip, file, "%d\n", SCIPsdpVarmapperGetSdpIndex(varmapper, vars[v]));
      }

      SCIPinfoMessage(scip, file, "\n");
   }

   /* count number of equality/geq/leq constraints */
   neqconss = 0;
   ngeqconss = 0;
   nleqconss = 0;
   for (c = 0; c < nconss; c++)
   {
      /* only count linear constraints */
      if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "linear") != 0 )
         continue;

#ifdef CBF_CHECK_NONNEG
      if ( consdisabled[c] )
         continue;
#endif

      if ( ( ! SCIPisInfinity(scip, -1 * SCIPgetLhsLinear(scip, conss[c]))) &&
            ( ! SCIPisInfinity(scip, SCIPgetRhsLinear(scip, conss[c])) ) &&
            SCIPisEQ(scip, SCIPgetLhsLinear(scip, conss[c]), SCIPgetRhsLinear(scip, conss[c])) )
         neqconss++;
      else
      {
         if ( ! SCIPisInfinity(scip, -1 * SCIPgetLhsLinear(scip, conss[c])) )
            ngeqconss++;
         if ( ! SCIPisInfinity(scip, SCIPgetRhsLinear(scip, conss[c])) )
            nleqconss++;
      }
   }

   /* write constraints */
   if ( neqconss + ngeqconss + nleqconss > 0 )
   {
      if ( ((neqconss > 0) && (ngeqconss == 0) && (nleqconss == 0))
            || ((neqconss == 0) && (ngeqconss > 0) && (nleqconss == 0))
            || ((neqconss == 0) && (ngeqconss == 0) && (nleqconss > 0)) )
         SCIPinfoMessage(scip, file, "CON\n%d 1\n", neqconss + ngeqconss + nleqconss);
      else if ( ((neqconss > 0) && (ngeqconss > 0) && (nleqconss == 0))
            || ((neqconss > 0) && (ngeqconss == 0) && (nleqconss > 0))
            || ((neqconss == 0) && (ngeqconss > 0) && (nleqconss > 0)) )
         SCIPinfoMessage(scip, file, "CON\n%d 2\n", neqconss + ngeqconss + nleqconss);
      else
      {
         assert( (neqconss > 0) && (ngeqconss > 0) && (nleqconss > 0) );
         SCIPinfoMessage(scip, file, "CON\n%d 3\n", neqconss + ngeqconss + nleqconss);
      }

      if ( ngeqconss > 0 )
         SCIPinfoMessage(scip, file, "L+ %d", ngeqconss);
      if ( nleqconss > 0 )
         SCIPinfoMessage(scip, file, "L- %d", nleqconss);
      if ( neqconss > 0 )
         SCIPinfoMessage(scip, file, "L= %d", neqconss);

      SCIPinfoMessage(scip, file, "\n\n");
   }

   /* count number of SDP constraints (conshdlrGetNConss doesn't seem to work before transformation) */
   nsdpconss = 0;
   for (c = 0; c < nconss; c++)
   {
      if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "SDP") == 0 )
         nsdpconss++;
   }

   /* write SDP constraints */
   SCIPinfoMessage(scip, file, "PSDCON\n%d\n", nsdpconss);

   for (c = 0; c < nconss; c++)
   {
      if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "SDP") != 0 )
         continue;

      SCIPinfoMessage(scip, file, "%d\n", SCIPconsSdpGetBlocksize(scip, conss[c]));
   }

   SCIPinfoMessage(scip, file, "\n");

   /* count number of nonzero objective coefficients */
   nobjnonz = 0;
   for (v = 0; v < nvars; v++)
   {
      if ( ! SCIPisZero(scip, SCIPvarGetObj(vars[v])) )
         nobjnonz++;
   }

   /* write objective */
   SCIPinfoMessage(scip, file, "OBJACOORD\n%d\n", nobjnonz);

   for (v = 0; v < nvars; v++)
   {
      if ( ! SCIPisZero(scip, SCIPvarGetObj(vars[v])) )
      {
         SCIPinfoMessage(scip, file, "%d %.9f\n", SCIPsdpVarmapperGetSdpIndex(varmapper, vars[v]), SCIPvarGetObj(vars[v]));
      }
   }

   SCIPinfoMessage(scip, file, "\n");

   if ( neqconss + ngeqconss + nleqconss > 0 )
   {
      /* count number of nonzero coefficients in linear constraints */
      nnonz = 0;
      nbnonz = 0;

      /* first iterate over all greater or equal constraints */
      for (c = 0; c < nconss; c++)
      {
         if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "linear") != 0 )
            continue;

         if ( SCIPisInfinity(scip, -1 * SCIPgetLhsLinear(scip, conss[c]))
               || SCIPisEQ(scip, SCIPgetLhsLinear(scip, conss[c]), SCIPgetRhsLinear(scip, conss[c])) )
            continue;
#ifdef CBF_CHECK_NONNEG
         if ( consdisabled[c] )
            continue;
#endif

         nnonz += SCIPgetNVarsLinear(scip, conss[c]);
         nbnonz += ( ! SCIPisZero(scip, SCIPgetLhsLinear(scip, conss[c])));
      }
      /* iterate over all less or equal constraints */
      for (c = 0; c < nconss; c++)
      {
         if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "linear") != 0 )
            continue;

         if ( SCIPisInfinity(scip, SCIPgetRhsLinear(scip, conss[c]))
               || SCIPisEQ(scip, SCIPgetLhsLinear(scip, conss[c]), SCIPgetRhsLinear(scip, conss[c])) )
            continue;
#ifdef CBF_CHECK_NONNEG
         if ( consdisabled[c] )
            continue;
#endif

         nnonz += SCIPgetNVarsLinear(scip, conss[c]);
         nbnonz += ( ! SCIPisZero(scip, SCIPgetRhsLinear(scip, conss[c])));
      }
      /* finally iterate over all equality constraints */
      for (c = 0; c < nconss; c++)
      {
         if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "linear") != 0 )
            continue;

         if ( SCIPisInfinity(scip, -1 * SCIPgetLhsLinear(scip, conss[c]))
               || SCIPisInfinity(scip, SCIPgetRhsLinear(scip, conss[c]))
               || ( ! SCIPisEQ(scip, SCIPgetLhsLinear(scip, conss[c]), SCIPgetRhsLinear(scip, conss[c]))) )
            continue;
#ifdef CBF_CHECK_NONNEG
         if ( consdisabled[c] )
            continue;
#endif

         nnonz += SCIPgetNVarsLinear(scip, conss[c]);
         nbnonz += ( ! SCIPisZero(scip, SCIPgetLhsLinear(scip, conss[c])));
      }

      /* write linear nonzero coefficients */
      SCIPinfoMessage(scip, file, "ACOORD\n%d\n", nnonz);
      consind = 0;

      /* first iterate over all greater or equal constraints */
      for (c = 0; c < nconss; c++)
      {
         if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "linear") != 0 )
            continue;

         if ( SCIPisInfinity(scip, -1 * SCIPgetLhsLinear(scip, conss[c]))
               || SCIPisEQ(scip, SCIPgetLhsLinear(scip, conss[c]), SCIPgetRhsLinear(scip, conss[c])) )
            continue;
#ifdef CBF_CHECK_NONNEG
         if ( consdisabled[c] )
            continue;
#endif

         linvars = SCIPgetVarsLinear(scip, conss[c]);
         linvals = SCIPgetValsLinear(scip, conss[c]);

         for (v = 0; v < SCIPgetNVarsLinear(scip, conss[c]); v++)
         {
            SCIPinfoMessage(scip, file, "%d %d %.9f\n", consind, SCIPsdpVarmapperGetSdpIndex(varmapper, linvars[v]), linvals[v]);
         }
         consind++;
      }
      /* iterate over all less or equal constraints */
      for (c = 0; c < nconss; c++)
      {
         if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "linear") != 0 )
            continue;

         if ( SCIPisInfinity(scip, SCIPgetRhsLinear(scip, conss[c]))
               || SCIPisEQ(scip, SCIPgetLhsLinear(scip, conss[c]), SCIPgetRhsLinear(scip, conss[c])) )
            continue;
#ifdef CBF_CHECK_NONNEG
         if ( consdisabled[c] )
            continue;
#endif

         linvars = SCIPgetVarsLinear(scip, conss[c]);
         linvals = SCIPgetValsLinear(scip, conss[c]);

         for (v = 0; v < SCIPgetNVarsLinear(scip, conss[c]); v++)
         {
            SCIPinfoMessage(scip, file, "%d %d %.9f\n", consind, SCIPsdpVarmapperGetSdpIndex(varmapper, linvars[v]), linvals[v]);
         }
         consind++;
      }
      /* finally iterate over all equality constraints */
      for (c = 0; c < nconss; c++)
      {
         if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "linear") != 0 )
            continue;

         if ( SCIPisInfinity(scip, -1 * SCIPgetLhsLinear(scip, conss[c]))
               || SCIPisInfinity(scip, SCIPgetRhsLinear(scip, conss[c]))
               || ( ! SCIPisEQ(scip, SCIPgetLhsLinear(scip, conss[c]), SCIPgetRhsLinear(scip, conss[c]))) )
            continue;
#ifdef CBF_CHECK_NONNEG
         if ( consdisabled[c] )
            continue;
#endif

         linvars = SCIPgetVarsLinear(scip, conss[c]);
         linvals = SCIPgetValsLinear(scip, conss[c]);

         for (v = 0; v < SCIPgetNVarsLinear(scip, conss[c]); v++)
         {
            SCIPinfoMessage(scip, file, "%d %d %.9f\n", consind, SCIPsdpVarmapperGetSdpIndex(varmapper, linvars[v]), linvals[v]);
         }
         consind++;
      }

      SCIPinfoMessage(scip, file, "\n");

      /* write constant part of linear constraints */
      SCIPinfoMessage(scip, file, "BCOORD\n%d\n", nbnonz);
      consind = 0;

      /* first iterate over all greater or equal constraints */
      for (c = 0; c < nconss; c++)
      {
         if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "linear") != 0 )
            continue;

         if ( SCIPisInfinity(scip, -1 * SCIPgetLhsLinear(scip, conss[c]))
               || SCIPisEQ(scip, SCIPgetLhsLinear(scip, conss[c]), SCIPgetRhsLinear(scip, conss[c])) )
            continue;
#ifdef CBF_CHECK_NONNEG
         if ( consdisabled[c] )
            continue;
#endif

         /* write the entry; *(-1) because we have Ax - lhs >= 0 */
         if ( ! SCIPisZero(scip, SCIPgetLhsLinear(scip, conss[c])) )
         {
            SCIPinfoMessage(scip, file, "%d %.9f\n", consind, -1 * SCIPgetLhsLinear(scip, conss[c]));
         }
         consind++; /* counting the constraint numbers is independent of whether the lhs is nonzero */
      }
      /* iterate over all less or equal constraints */
      for (c = 0; c < nconss; c++)
      {
         if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "linear") != 0 )
            continue;

         if ( SCIPisInfinity(scip, SCIPgetRhsLinear(scip, conss[c]))
               || SCIPisEQ(scip, SCIPgetLhsLinear(scip, conss[c]), SCIPgetRhsLinear(scip, conss[c])) )
            continue;
#ifdef CBF_CHECK_NONNEG
         if ( consdisabled[c] )
            continue;
#endif

         /* write the entry; *(-1) because we have Ax - rhs <= 0 */
         if ( ! SCIPisZero(scip, SCIPgetRhsLinear(scip, conss[c])) )
         {
            SCIPinfoMessage(scip, file, "%d %.9f\n", consind, -1 * SCIPgetRhsLinear(scip, conss[c]));
         }
         consind++; /* counting the constraint numbers is independent of whether the rhs is nonzero */
      }
      /* finally iterate over all equality constraints */
      for (c = 0; c < nconss; c++)
      {
         if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "linear") != 0 )
            continue;

         if ( SCIPisInfinity(scip, -1 * SCIPgetLhsLinear(scip, conss[c]))
               || SCIPisInfinity(scip, SCIPgetRhsLinear(scip, conss[c]))
               || ( ! SCIPisEQ(scip, SCIPgetLhsLinear(scip, conss[c]), SCIPgetRhsLinear(scip, conss[c]))) )
            continue;
#ifdef CBF_CHECK_NONNEG
         if ( consdisabled[c] )
            continue;
#endif

         /* write the entry; *(-1) because we have Ax - lhs = 0 */
         if ( ! SCIPisZero(scip, SCIPgetLhsLinear(scip, conss[c])) )
         {
            SCIPinfoMessage(scip, file, "%d %.9f\n", consind, -1 * SCIPgetLhsLinear(scip, conss[c]));
         }
         consind++;/* counting the constraint numbers is independent of whether the lhs is nonzero */
      }
      SCIPinfoMessage(scip, file, "\n");
   }

   /*count SDP nonzeros */
   totalsdpnnonz = 0;
   totalsdpconstnnonz = 0;
   for (c = 0; c < nconss; c++)
   {
      if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])),"SDP") != 0 )
         continue;

      SCIP_CALL( SCIPconsSdpGetNNonz(scip, conss[c], &sdpnnonz, &sdpconstnnonz) );
      totalsdpnnonz += sdpnnonz;
      totalsdpconstnnonz += sdpconstnnonz;
   }

   /* allocate memory for SDPdata */
   SCIP_CALL( SCIPallocBufferArray(scip, &sdpnvarnonz, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sdpcol, totalsdpnnonz) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sdprow, totalsdpnnonz) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sdpval, totalsdpnnonz) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sdpvars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sdpconstcol, totalsdpconstnnonz) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sdpconstrow, totalsdpconstnnonz) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sdpconstval, totalsdpconstnnonz) );

   sdparraylength = totalsdpnnonz;
   sdpconstnnonz = totalsdpconstnnonz;

   /* write SDP nonzeros */
   SCIPinfoMessage(scip, file, "HCOORD\n%d\n", totalsdpnnonz);
   consind = 0;
   for (c = 0; c < nconss; c++)
   {
      if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])),"SDP") != 0 )
         continue;

      /* initialization for SDPconsSDPGetData-call */
      sdparraylength = totalsdpnnonz;
      sdpconstnnonz = totalsdpconstnnonz;

      SCIP_CALL( SCIPconsSdpGetData(scip, conss[c], &sdpnvars, &sdpnnonz, &sdpblocksize, &sdparraylength, sdpnvarnonz,
            sdpcol, sdprow, sdpval, sdpvars, &sdpconstnnonz, sdpconstcol, sdpconstrow, sdpconstval) );

      assert( sdpconstnnonz <= totalsdpconstnnonz );
      assert( sdparraylength <= totalsdpnnonz);

      for (v = 0; v < sdpnvars; v++)
      {
         for (i = 0; i < sdpnvarnonz[v]; i++)
         {
            SCIPinfoMessage(scip, file, "%d %d %d %d %.9f\n", consind, SCIPsdpVarmapperGetSdpIndex(varmapper, sdpvars[v]),
                  sdprow[v][i], sdpcol[v][i], sdpval[v][i]);
         }
      }
      consind++;
   }
   SCIPinfoMessage(scip, file, "\n");

   /* write nonzeros of constant part of SDP constraint */
   SCIPinfoMessage(scip, file, "DCOORD\n%d\n", totalsdpconstnnonz);
   consind = 0;
   for (c = 0; c < nconss; c++)
   {
      if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "SDP") != 0 )
         continue;

      /* initialization for SDPconsSDPGetData-call */
      sdparraylength = totalsdpnnonz;
      sdpconstnnonz = totalsdpconstnnonz;

      SCIP_CALL( SCIPconsSdpGetData(scip, conss[c], &sdpnvars, &sdpnnonz, &sdpblocksize, &sdparraylength, sdpnvarnonz,
            sdpcol, sdprow, sdpval, sdpvars, &sdpconstnnonz, sdpconstcol, sdpconstrow, sdpconstval) );

      assert( sdpconstnnonz <= totalsdpconstnnonz );
      assert( sdparraylength <= totalsdpnnonz);

      for (i = 0; i < sdpconstnnonz; i++)
      {
         SCIPinfoMessage(scip, file, "%d %d %d %.9f\n", consind, sdpconstrow[i], sdpconstcol[i], -1* sdpconstval[i]);
      }
      consind++;
   }

   SCIPfreeBufferArray(scip, &sdpconstval);
   SCIPfreeBufferArray(scip, &sdpconstrow);
   SCIPfreeBufferArray(scip, &sdpconstcol);
   SCIPfreeBufferArray(scip, &sdpvars);
   SCIPfreeBufferArray(scip, &sdpval);
   SCIPfreeBufferArray(scip, &sdprow);
   SCIPfreeBufferArray(scip, &sdpcol);
   SCIPfreeBufferArray(scip, &sdpnvarnonz);
   SCIP_CALL( SCIPsdpVarmapperFree(scip, &varmapper) );
#ifdef CBF_CHECK_NONNEG
   SCIPfreeBufferArray(scip, &negorth);
   SCIPfreeBufferArray(scip, &posorth);
   SCIPfreeBufferArray(scip, &consdisabled);
#endif

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/*
 * reader specific interface methods
 */

/** includes the CBF file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderCbf(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READERDATA* readerdata;
   SCIP_READER* reader;

   /* create reader data */
   readerdata = NULL;

   /* include reader */
   SCIP_CALL( SCIPincludeReaderBasic(scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION, readerdata) );

   assert(reader != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetReaderCopy(scip, reader, readerCopyCbf) );
   SCIP_CALL( SCIPsetReaderRead(scip, reader, readerReadCbf) );
   SCIP_CALL( SCIPsetReaderWrite(scip, reader, readerWriteCbf) );

   return SCIP_OKAY;
}
