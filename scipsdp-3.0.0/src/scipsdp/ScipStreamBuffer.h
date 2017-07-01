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

/**@file   ScipStreamBuffer.h
 * @brief  An std::streambuf that uses SCIP I/O routines (suitable for reading)
 * @author Lars Schewe
 */

#ifndef SCIPSTREAM_H
#define SCIPSTREAM_H

#include <streambuf>                    // for streamsize, streambuf

#include <cstddef>                      // for size_t
#include "scip/scip.h"

class ScipStreamBuffer : public std::streambuf
{
 public:
   ScipStreamBuffer(SCIP* scip, SCIP_FILE* file, bool close_on_exit);

   ~ScipStreamBuffer();

 protected:
   /// the underflow function is responsible for the refilling of the buffer
   virtual int underflow();

   virtual std::streamsize xsgetn(char *dest, std::streamsize request);

   SCIP* scip_;
   SCIP_FILE* file_;
   char* g_buffer_; /// pointer to the get-buffer
   size_t g_buffer_size_; /// size of the get-buffer
   bool close_on_exit_;
};
#endif
