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

/**@file   BlockMemoryAllocator.h
 * @brief  An STL allocator class using SCIP block memory
 * @author Lars Schewe
 */

#ifndef BLOCKMEMORYALLOCATOR_H
#define BLOCKMEMORYALLOCATOR_H

/**
 * direct implementation of the Allocator interface
 *
 */
template <class T> class BlockMemoryAllocator
{
public:
  typedef T                 value_type;
  typedef value_type*       pointer;
  typedef const value_type* const_pointer;
  typedef value_type&       reference;
  typedef const value_type& const_reference;
  typedef std::size_t       size_type;
  typedef std::ptrdiff_t    difference_type;

  template <class U>
  struct rebind { typedef BlockMemoryAllocator<U> other; };

 BlockMemoryAllocator(SCIP* scip) : scip_(scip) {}

 BlockMemoryAllocator(const BlockMemoryAllocator& other) : scip_(other.scip_) {}
  template <class U>
     BlockMemoryAllocator(const BlockMemoryAllocator<U>& other) : scip_(other.scip_) {}
  ~BlockMemoryAllocator()
  {
     scip_ = NULL;
  }

  pointer address(reference x) const { return &x; }

  const_pointer address(const_reference x) const
  {
    return x;
  }

  pointer allocate(size_type n, const_pointer = 0)
  {
     void* p;
     SCIP_CALL_ABORT(SCIPallocBlockMemorySize(scip_, &p, n * sizeof(T)));

     if (!p)
        throw std::bad_alloc();
     return static_cast<pointer>(p);
  }

  void deallocate(pointer p, size_type n)
  {
     SCIPfreeBlockMemorySize(scip_, &p, n * sizeof(T));
  }

  size_type max_size() const {
    return static_cast<size_type>(-1) / sizeof(T);
  }

  void construct(pointer p, const value_type& x) {
     new(p) value_type(x);
  }
  void destroy(pointer p) { p->~value_type(); }

  template<typename S> friend inline bool operator==(const BlockMemoryAllocator<S>& left, const BlockMemoryAllocator<S>& right);
  template<typename S> friend inline bool operator!=(const BlockMemoryAllocator<S>& left, const BlockMemoryAllocator<S>& right);

  void operator=(BlockMemoryAllocator const &b)
  {
     scip_ = b.scip_;
  }

 private:
  SCIP* scip_;
};

template<> class BlockMemoryAllocator<void>
{
  typedef void        value_type;
  typedef void*       pointer;
  typedef const void* const_pointer;

  template <class U>
  struct rebind { typedef BlockMemoryAllocator<U> other; };
};


template <class T>
inline bool operator==(const BlockMemoryAllocator<T>& left,
                       const BlockMemoryAllocator<T>& right) {
   return left.scip_ == right.scip_;
}

template <class T>
inline bool operator!=(const BlockMemoryAllocator<T>& left,
                       const BlockMemoryAllocator<T>& right) {
   return !(left == right);
}

#endif
