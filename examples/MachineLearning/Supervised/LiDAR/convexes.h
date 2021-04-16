/* Author: Samuel Hornus <samuel.hornus@inria.fr>
 * Copyright © Inria, 2017
 * Licence: Creative Commons CC BY-ND 3.0 license available online at
 * http://creativecommons.org/licenses/by-nd/3.0/
 *
 Changes made to convexes, copying only the bits without CGAL and adding a constructor to VerticesOnly smitha gopinath */

#ifndef RAY_CONVEXES_H
#define RAY_CONVEXES_H

#include "vec.h"
#include <algorithm>
#include <vector>
#include <stack>
#include <map>

extern unsigned long planeStatPerPair;

// ---------------------------------------------------------------

struct RandomPointOnSphere {
//  Point_3 point3() const {
//    float z = 2.0f * drand48() - 1.0f;
//    float angle = 2.0f * M_PI * drand48();
//    float rad = ::sqrtf(1.0f - z*z);
//    return Point_3(rad * cosf(angle), rad * sinf(angle), z);
//  }
  Vec3f vec3f() const {
    float z = 2.0f * drand48() - 1.0f;
    float angle = 2.0f * M_PI * drand48();
    float rad = ::sqrtf(1.0f - z*z);
    return Vec3f(rad * cosf(angle), rad * sinf(angle), z);
  }
  RandomPointOnSphere & operator++(int) {
    return (*this);
  }
  RandomPointOnSphere & operator++() {
    return (*this);
  }
};

struct VerticesOnly {
  std::vector<Vec3f> vertices_;
  Vec3f center_;

  public:

  const Vec3f & vertex(const int i) const { return vertices_[i]; }
  const Vec3f & center() const { return center_; }
  size_t size() const { return vertices_.size(); }
  size_t nbVertices() const { return vertices_.size(); }

  //VerticesOnly(const Polyhedron_3 &) { exit(-1); }
  VerticesOnly(int n, float shift) {
    Vec3f s(shift, 0.0f, 0.0f);
    RandomPointOnSphere rps;
    Vec3f center;
    for( int i = 0; i < n; ++i ) {
      vertices_.push_back(rps.vec3f() + s);
      center = center + vertices_.back();
    }
    center_ = (1.0f / n) * center;
  }
  VerticesOnly(std::vector<std::vector<double>> points) {
        Vec3f center;
          for( int i = 0; i < points.size(); ++i ) {
              auto fl_x=float(points[i][0]);
              auto fl_y=float(points[i][1]);
              auto fl_z=float(points[i][2]);
              Vec3f s(fl_x,fl_y,fl_z);
            vertices_.push_back(s);
            center = center + vertices_.back();
          }
          center_ = (1.0f / points.size()) * center;
        }

//  void makeFrustum() {
//    vertices_.clear();
//    vertices_.reserve(8);
//    Frustum f;
//    for( int i = 0; i < 8; ++i )
//      vertices_.push_back(f.vertex(i));
//    center_ = f.center();
//  }

  void minimizeInDirection(const Vec3f & dir, float & mini, int & vA) const {
    mini = vertices_[0] | dir;
    vA = 0;
    const int na = vertices_.size();
    for( int i = 1; i < na; ++i ) {
      float tempA = vertices_[i] | dir;
      if( tempA < mini ) { mini = tempA; vA = i; }
    }
  }

  void maximizeInDirection(const Vec3f & dir, float & maxi, int & vA) const {
    maxi = vertices_[0] | dir;
    vA = 0;
    const int na = vertices_.size();
    for( int i = 1; i < na; ++i ) {
      float tempA = vertices_[i] | dir;
      if( tempA > maxi ) { maxi = tempA; vA = i; }
    }
  }

  void differenceInDir(const VerticesOnly & B, int & vA, int & vB, float & maxOverA, float & minOverB, const Vec3f & dir) const {
    // difference above is: A - B
    ++planeStatPerPair;
    maxOverA = vertices_[0] | dir;
    minOverB = B.vertices_[0] | dir;
    vA = vB = 0;
    const int na = vertices_.size();
    for( int i = 1; i < na; ++i ) {
      float tempA = vertices_[i] | dir;
      if( tempA > maxOverA ) { maxOverA = tempA; vA = i; }
    }
    const int nb = B.vertices_.size();
    for( int i = 1; i < nb; ++i ) {
      float tempB = B.vertices_[i] | dir;
      if( tempB < minOverB ) { minOverB = tempB; vB = i; }
    }
  }

  bool differenceCoversZeroInDir(const VerticesOnly & B, int & vA, int & vB, const Vec3f & dir) const {
    // difference above is: A - B
    ++planeStatPerPair;
    float maxOverA = vertices_[0] | dir;
    float minOverB = B.vertices_[0] | dir;
    vA = vB = 0;
    const int na = vertices_.size();
    for( int i = 1; i < na; ++i ) {
      float tempA = vertices_[i] | dir;
      if( tempA > maxOverA ) { maxOverA = tempA; vA = i; }
    }
    const int nb = B.vertices_.size();
    for( int i = 1; i < nb; ++i ) {
      float tempB = B.vertices_[i] | dir;
      if( tempB < minOverB ) { minOverB = tempB; vB = i; }
    }
    return maxOverA >= minOverB;
  }

  bool differenceCoversZeroInDirLazy(const VerticesOnly & B, int & vA, int & vB, const Vec3f & dir) const {
    // difference above is: A - B
    ++planeStatPerPair;
    float maxOverA = vertices_[0] | dir;
    float minOverB = B.vertices_[0] | dir;
    vA = vB = 0;
    //if( maxOverA >= minOverB ) return true;
    const int na = vertices_.size();
    const int nb = B.vertices_.size();
    int i;
    if( na <= nb ) {
      for( i = 1; i < na; ++i ) {
        float tempA = vertices_[i] | dir;
        float tempB = B.vertices_[i] | dir;
        if( tempA > maxOverA ) { maxOverA = tempA; vA = i; }
        if( tempB < minOverB ) { minOverB = tempB; vB = i; }
        if( maxOverA >= minOverB ) return true;
      }
      for( ; i < nb; ++i ) {
        float tempB = B.vertices_[i] | dir;
        if( tempB < minOverB ) {
          minOverB = tempB; vB = i;
          if( maxOverA >= minOverB ) return true;
        }
      }
    } else {
      for( i = 1; i < nb; ++i ) {
        float tempA = vertices_[i] | dir;
        float tempB = B.vertices_[i] | dir;
        if( tempA > maxOverA ) { maxOverA = tempA; vA = i; }
        if( tempB < minOverB ) { minOverB = tempB; vB = i; }
        if( maxOverA >= minOverB ) return true;
      }
      for( ; i < na; ++i ) {
        float tempA = vertices_[i] | dir;
        if( tempA > maxOverA ) {
          maxOverA = tempA; vA = i;
          if( maxOverA >= minOverB ) return true;
        }
      }
    }
    return false;
  }
};
#endif // RAY_CONVEXES_H
