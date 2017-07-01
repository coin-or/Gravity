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

/**@file   ScipStreamBuffer.cpp
 * @brief  An std::streambuf that uses SCIP I/O routines (suitable for reading)
 * @author Lars Schewe
 */

#include "ScipStreamBuffer.h"

#include <algorithm>                    // for min
#include <cstring>                      // for memcpy
#include "scip/scip.h"                  // for SCIPallocBufferArray, etc

#define SCIPSTREAMBUFFERSIZE 256

ScipStreamBuffer::ScipStreamBuffer(SCIP* scip, SCIP_FILE* file, bool close_on_exit) : scip_(scip), file_(file), g_buffer_(NULL),
                                                                                      g_buffer_size_(SCIPSTREAMBUFFERSIZE), close_on_exit_(close_on_exit)
{
   SCIP_CALL_ABORT(SCIPallocBufferArray(scip_, &g_buffer_, g_buffer_size_));
   /// set buffer in empty state
   setg(g_buffer_, g_buffer_ + g_buffer_size_ , g_buffer_ + g_buffer_size_);
}

ScipStreamBuffer::~ScipStreamBuffer()
{
   SCIPfreeBufferArray(scip_, &g_buffer_);
   scip_ = NULL;
   g_buffer_= NULL;
   g_buffer_size_ = 0;

   if (close_on_exit_)
   {
      SCIPfclose(file_);
   }
}

int ScipStreamBuffer::underflow()
{
   if (gptr() < egptr())
   {
      // still unused chars in the buffer
      // just return
      return *gptr();
   }

   // refill the buffer
   int nresult = SCIPfread(g_buffer_, 1, g_buffer_size_, file_);

   if (nresult <= 0) /// we read nothing => bad
   {
      return EOF;
   }

   // set up buffer, we read nresult bytes, that's where our new end is
   setg(g_buffer_, g_buffer_, g_buffer_ + nresult);

   return *gptr();
}

std::streamsize ScipStreamBuffer::xsgetn(char *dest, std::streamsize request)
{
   int ndone = 0;

   while (request != 0)
   {
      if (!in_avail())
      {
         // we have an empty buffer, so refresh
         if (underflow() == EOF)
         {
            // buffer remains empty
            break;
         }
      }

      int available = std::min(in_avail(), request);

      // copy the available bytes
      memcpy(dest + ndone, gptr(), available);
      // shift read pointer
      gbump(available);

      ndone += available;
      request -= available;
   }

   return ndone;
}
