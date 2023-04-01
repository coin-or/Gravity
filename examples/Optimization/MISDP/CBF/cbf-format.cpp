////
////  cbf-format.cpp
////  misdp
////
////  Created by Smitha on 7/7/22.
////
//
//// Copyright (c) 2012 by Zuse-Institute Berlin and the Technical University of Denmark.
////
//// This software is provided 'as-is', without any express or implied
//// warranty. In no event will the authors be held liable for any damages
//// arising from the use of this software.
////
//// Permission is granted to anyone to use this software for any purpose,
//// including commercial applications, and to alter it and redistribute it
//// freely, subject to the following restrictions:
////
//// 1. The origin of this software must not be misrepresented; you must not
////    claim that you wrote the original software. If you use this software
////    in a product, an acknowledgment in the product documentation would be
////    appreciated but is not required.
//// 2. Altered source versions must be plainly marked as such, and must not be
////    misrepresented as being the original software.
//// 3. This notice may not be removed or altered from any source distribution.
//
//#include "cbf-format.h"
//#include <string.h>
//
//
//// -------------------------------------
//// Global variable
//// -------------------------------------
//
//// File format buffers
//char CBF_LINE_BUFFER[CBF_MAX_LINE];
//char CBF_NAME_BUFFER[CBF_MAX_NAME];
//
//// Names of the scalar cones
//const char * CBF_CONENAM_FREE = "F";
//const char * CBF_CONENAM_ZERO = "L=";
//const char * CBF_CONENAM_POS = "L+";
//const char * CBF_CONENAM_NEG = "L-";
//const char * CBF_CONENAM_QUAD = "Q";
//const char * CBF_CONENAM_RQUAD = "QR";
//const char * CBF_CONENAM_PEXP = "EXP";
//const char * CBF_CONENAM_DEXP = "EXP*";
//
//// Names of the objective senses
//const char * CBF_OBJSENSENAM_MIN = "MIN";
//const char * CBF_OBJSENSENAM_MAX = "MAX";
//
//// -------------------------------------
//// Function definitions
//// -------------------------------------
//
//CBFresponsee CBF_conetostr(CBFscalarconee cone, const char **str)
//{
//  switch (cone) {
//  case CBF_CONE_FREE:
//    *str = CBF_CONENAM_FREE;
//    break;
//  case CBF_CONE_POS:
//    *str = CBF_CONENAM_POS;
//    break;
//  case CBF_CONE_NEG:
//    *str = CBF_CONENAM_NEG;
//    break;
//  case CBF_CONE_ZERO:
//    *str = CBF_CONENAM_ZERO;
//    break;
//  case CBF_CONE_QUAD:
//    *str = CBF_CONENAM_QUAD;
//    break;
//  case CBF_CONE_RQUAD:
//    *str = CBF_CONENAM_RQUAD;
//  case CBF_CONE_PEXP:
//    *str = CBF_CONENAM_PEXP;
//  case CBF_CONE_DEXP:
//    *str = CBF_CONENAM_DEXP;
//    break;
//  default:
//    return CBF_RES_ERR;
//  }
//  return CBF_RES_OK;
//}
//
//CBFresponsee CBF_strtocone(const char *str, CBFscalarconee *cone)
//{
//  if (strcmp(str, CBF_CONENAM_FREE) == 0)
//    *cone = CBF_CONE_FREE;
//  else if (strcmp(str, CBF_CONENAM_POS) == 0)
//    *cone = CBF_CONE_POS;
//  else if (strcmp(str, CBF_CONENAM_NEG) == 0)
//    *cone = CBF_CONE_NEG;
//  else if (strcmp(str, CBF_CONENAM_ZERO) == 0)
//    *cone = CBF_CONE_ZERO;
//  else if (strcmp(str, CBF_CONENAM_QUAD) == 0)
//    *cone = CBF_CONE_QUAD;
//  else if (strcmp(str, CBF_CONENAM_RQUAD) == 0)
//    *cone = CBF_CONE_RQUAD;
//  else if (strcmp(str, CBF_CONENAM_PEXP) == 0)
//    *cone = CBF_CONE_PEXP;
//  else if (strcmp(str, CBF_CONENAM_DEXP) == 0)
//    *cone = CBF_CONE_DEXP;
//  else
//    return CBF_RES_ERR;
//
//  return CBF_RES_OK;
//}
//
//CBFresponsee CBF_objsensetostr(CBFobjsensee objsense, const char **str)
//{
//  switch (objsense) {
//  case CBF_OBJ_MINIMIZE:
//    *str = CBF_OBJSENSENAM_MIN;
//    break;
//  case CBF_OBJ_MAXIMIZE:
//    *str = CBF_OBJSENSENAM_MAX;
//    break;
//  default:
//    return CBF_RES_ERR;
//  }
//  return CBF_RES_OK;
//}
//
//CBFresponsee CBF_strtoobjsense(const char *str, CBFobjsensee *objsense)
//{
//  if (strcmp(str, CBF_OBJSENSENAM_MIN) == 0)
//    *objsense = CBF_OBJ_MINIMIZE;
//  else if (strcmp(str, CBF_OBJSENSENAM_MAX) == 0)
//    *objsense = CBF_OBJ_MAXIMIZE;
//  else
//    return CBF_RES_ERR;
//
//  return CBF_RES_OK;
//}
