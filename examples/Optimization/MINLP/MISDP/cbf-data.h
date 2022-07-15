#ifndef cbf_data_h
#define cbf_data_h
// Copyright (c) 2012 by Zuse-Institute Berlin and the Technical University of Denmark.
//
// This software is provided 'as-is', without any express or implied
// warranty. In no event will the authors be held liable for any damages
// arising from the use of this software.
//
// Permission is granted to anyone to use this software for any purpose,
// including commercial applications, and to alter it and redistribute it
// freely, subject to the following restrictions:
//
// 1. The origin of this software must not be misrepresented; you must not
//    claim that you wrote the original software. If you use this software
//    in a product, an acknowledgment in the product documentation would be
//    appreciated but is not required.
// 2. Altered source versions must be plainly marked as such, and must not be
//    misrepresented as being the original software.
// 3. This notice may not be removed or altered from any source distribution.

#ifndef CBF_CBF_DATA_H
#define CBF_CBF_DATA_H

#define CBF_VERSION     2
#define CBF_MAX_LINE  512       // Last 3 chars reserved for '\r\n\0'
#define CBF_MAX_NAME  512


typedef enum CBFobjsense_enum {
  CBF_OBJ_BEGIN = 0,
  CBF_OBJ_END = 2,

  CBF_OBJ_MINIMIZE = 0,
  CBF_OBJ_MAXIMIZE = 1
} CBFobjsensee;


typedef enum CBFscalarcone_enum {
  CBF_CONE_BEGIN = 0,
  CBF_CONE_END = 8,

  CBF_CONE_FREE = 0,
  CBF_CONE_POS = 1,
  CBF_CONE_NEG = 2,
  CBF_CONE_ZERO = 3,
  CBF_CONE_QUAD = 4,
  CBF_CONE_RQUAD = 5,
  CBF_CONE_PEXP = 6,
  CBF_CONE_DEXP = 7
} CBFscalarconee;


typedef struct CBFdata_struct {

  //
  // Problem format
  //
  int ver;

  //
  // Problem structure
  //
  CBFobjsensee    objsense;

  long long int   mapnum;
  long long int   mapstacknum;
  long long int  *mapstackdim;
  CBFscalarconee *mapstackdomain;

  long long int   varnum;
  long long int   varstacknum;
  long long int  *varstackdim;
  CBFscalarconee *varstackdomain;

  long long int   intvarnum;
  long long int  *intvar;

  int             psdmapnum;
  int            *psdmapdim;

  int             psdvarnum;
  int            *psdvardim;

  //
  // Coefficients of the objective scalar map
  //
  long long int  objfnnz;
  int           *objfsubj;
  int           *objfsubk;
  int           *objfsubl;
  double        *objfval;

  long long int  objannz;
  long long int *objasubj;
  double        *objaval;

  double         objbval;

  //
  // Coefficients of the scalar maps
  //
  long long int  fnnz;
  long long int *fsubi;
  int           *fsubj;
  int           *fsubk;
  int           *fsubl;
  double        *fval;

  long long int  annz;
  long long int *asubi;
  long long int *asubj;
  double        *aval;

  long long int  bnnz;
  long long int *bsubi;
  double        *bval;

  //
  // Coefficients of positive semidefinite maps
  //
  long long int  hnnz;
  int           *hsubi;
  long long int *hsubj;
  int           *hsubk;
  int           *hsubl;
  double        *hval;

  long long int  dnnz;
  int           *dsubi;
  int           *dsubk;
  int           *dsubl;
  double        *dval;

} CBFdata;

#endif



#endif /* cbf_data_h */
