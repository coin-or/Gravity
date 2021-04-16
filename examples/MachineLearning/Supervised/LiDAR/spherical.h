/* Author: Samuel Hornus <samuel.hornus@inria.fr>
 * Copyright © Inria, 2017
 * Licence: Creative Commons CC BY-ND 3.0 license available online at
 * http://creativecommons.org/licenses/by-nd/3.0/
 */

#ifndef RAY_PHERICAL_H
#define RAY_PHERICAL_H

#include "vec.h"
#include <vector>
#include <iostream>

struct SphericalPolygonElement {
  Vec3f vertex_;
  Vec3f north_;
  Vec3f silVertex_;
  SphericalPolygonElement(){}
  SphericalPolygonElement(const Vec3f & sil)
    : north_(sil.normalized()), silVertex_(sil) {}
  //SphericalPolygonElement(const Vec3f & v, const Vec3f & n)
  //  : vertex_(v), north_(n), silVertex_(n)  {}
  SphericalPolygonElement(const Vec3f & v, const Vec3f & n, const Vec3f & sil)
    : vertex_(v), north_(n), silVertex_(sil) {}
};

struct SphericalPolygon : public std::vector<SphericalPolygonElement> {

  typedef std::vector<SphericalPolygonElement> Base;
  typedef Base::iterator iterator;
  typedef Base::const_iterator const_iterator;
  SphericalPolygon() {
    reserve(16);
  }

  Vec3f averageDirection() const {
    // PRECONDITION : all northes are normalized.
    switch( size() ) {
      case 0 : return Vec3f(0.0f, 0.0f, 0.0f); break;
      case 1 : return begin()->north_; break;
      case 2 : return (*this)[0].north_ + (*this)[1].north_; break;
      default : {
                  Vec3f avg;
                  for( const SphericalPolygonElement & v : *this )
                    avg = avg + v.vertex_;
                  return avg;
                } break;
    }
  }

  void set_to_triangle(const Vec3f pts[3]) {
    clear();
    emplace_back(Vec3f::cross(pts[0], pts[1]).LInfNormalized(), pts[1].LInfNormalized(), pts[1]);
    emplace_back(Vec3f::cross(pts[1], pts[2]).LInfNormalized(), pts[2].LInfNormalized(), pts[2]);
    emplace_back(Vec3f::cross(pts[2], pts[0]).LInfNormalized(), pts[0].LInfNormalized(), pts[0]);
  }

  void clip(const Vec3f & OrigVertex, const Vec3f & silVertex, SphericalPolygon & result, bool doClean=true) const {
    // PRECONDITION : clipNorth, and all northes are normalized.
#define _ray_spherical_eps 1e-6f
    const int n = size();
    result.clear();
    switch( n ) {
      case 0 : break;
      case 1 : {
                 result = (*this);
                 Vec3f clipNorth = silVertex.normalized();
                 float dot = begin()->north_ | clipNorth;
                 if( dot < -0.99984769515 ) { // about one degree
                   // intersection of two almost opposite hemispheres ==> empty
                   result.clear();
                   break;
                 } else if( dot > 0.99984769515 ) {
                   break;
                 }
                 Vec3f v(Vec3f::cross(clipNorth, begin()->north_).LInfNormalized());
                 result.begin()->vertex_ = v;
                 result.emplace_back(v.negated(), clipNorth, OrigVertex);
                 break;
               }
      case 2 : {
                 result = (*this);
                 Vec3f clipNorth = silVertex.LInfNormalized();
                 iterator next = result.begin();
                 iterator cur = next++;
                 float vDot = begin()->vertex_ | clipNorth;
                 if( vDot >= _ray_spherical_eps ) {
                   // we'll get a triangle
                   next->vertex_ = Vec3f::cross(clipNorth, next->north_).LInfNormalized();
                   Vec3f v(Vec3f::cross(cur->north_, clipNorth).LInfNormalized());
                   result.emplace(next, v, clipNorth, OrigVertex);
                 } else if( vDot <= - _ray_spherical_eps ) {
                   // we'll get a triangle
                   cur->vertex_ = Vec3f::cross(clipNorth, cur->north_).LInfNormalized();
                   Vec3f v(Vec3f::cross(next->north_, clipNorth).LInfNormalized());
                   result.emplace_back(v, clipNorth, OrigVertex);
                 } else {
                   // we keep a moon crescent
                   float curTest(clipNorth | Vec3f::cross(cur->north_, cur->vertex_));
                   Vec3f nextTest(Vec3f::cross(next->north_, next->vertex_));
                   if( curTest > 0.0f ) {
                     if( (clipNorth | nextTest) <= 0.0f ) {
                       next->north_ = clipNorth;
                       next->silVertex_ = OrigVertex;
                       cur->vertex_ = Vec3f::cross(next->north_, cur->north_);
                       cur->vertex_.LInfNormalize();
                       next->vertex_ = cur->vertex_;
                       next->vertex_.negate();
                     } else {
                       // the crescent is unchanged
                       //std::cerr << "kept a crescent\n";
                     }
                   } else {
                     if( (clipNorth | nextTest) > 0.0f ) {
                       cur->north_ = clipNorth;
                       cur->silVertex_ = OrigVertex;
                       next->vertex_ = Vec3f::cross(cur->north_, next->north_);
                       next->vertex_.LInfNormalize();
                       cur->vertex_ = next->vertex_;
                       cur->vertex_.negate();
                     } else {
                       //std::cerr << "killed a crescent\n";
                       result.clear();
                     }
                   }
                 }
                 break;
               }
      default : { // n >= 3
                  int nbKept(0);
                  const_iterator cur = begin();
                  Vec3f clipNorth = silVertex.LInfNormalized();
                  float nextDot, curDot = clipNorth | cur->vertex_;
                  while( cur != end() ) {
                    if( cur+1 == end() )
                      nextDot = clipNorth | begin()->vertex_;
                    else
                      nextDot = clipNorth | (cur+1)->vertex_;
                    if( curDot >= _ray_spherical_eps ) { // cur is "IN"
                      ++nbKept;
                      result.push_back(*cur);
                      if( nextDot <= -_ray_spherical_eps ) { // next is "OUT"
                        result.emplace_back(Vec3f::cross(cur->north_, clipNorth).LInfNormalized(), clipNorth, OrigVertex);
                      }
                    } else if( curDot > -_ray_spherical_eps ) { // cur is "ON" the clipping plane
                      ++nbKept;
                      if ( nextDot <= -_ray_spherical_eps ) // next is "OUT"
                        result.emplace_back(cur->vertex_, clipNorth, OrigVertex);
                      else
                        result.push_back(*cur);
                    } else { // cur is "OUT"
                      if ( nextDot >= _ray_spherical_eps ) { // next is "IN"
                        result.emplace_back(Vec3f::cross(clipNorth, cur->north_).LInfNormalized(), cur->north_, cur->silVertex_);
                      }
                    }
                    curDot = nextDot;
                    ++cur;
                  }
                  if( (result.size() < 3/*too small*/) || ((nbKept == n)/*no change*/ && doClean) ) {
                    result.clear();
                  }
                  //if( nbKept == n ) {
                    //std::cerr << "**";
                  //}
                  break;
                }
    }
  }
  void clip(const Vec3f & silVertex, SphericalPolygon & result, bool doClean=true) const {
    clip(silVertex, silVertex, result, doClean);
  }
};

#endif // RAY_PHERICAL_H
