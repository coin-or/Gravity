/* Author: Samuel Hornus <samuel.hornus@inria.fr>
 * Copyright © Inria, 2017
 * Licence: Creative Commons CC BY-ND 3.0 license available online at
 * http://creativecommons.org/licenses/by-nd/3.0/
 */
#ifndef __VECTOR_STUFF__
#define __VECTOR_STUFF__

#include <cassert>
#include <vector>
#include <cmath>
#include <x86intrin.h>

/*
 *
 *  Vec2
 *
 */

template< typename NT >
struct Vec2
{
    typedef Vec2<NT>    Self;
    typedef NT          NumberType;

#define vecset(x,y) data_[0]=(x); data_[1]=(y);
    Vec2() { vecset(NT(0), NT(0)) }
    Vec2(const Self & s) { vecset(s.x(), s.y()) }
    explicit Vec2(const NT & x) { vecset(x,x) }
    Vec2(const NT & x, const NT & y) { vecset(x,y) }
    inline void set(const NT & x, const NT & y) {
      vecset(x,y)
    }
#undef vecset

    inline NT operator()(const int i) const { assert(2 > i && 0 <= i); return data_[i]; }
    inline NT operator[](const int i) const { return data_[i]; }
    inline NT x() const { return data_[0]; }
    inline NT y() const { return data_[1]; }
    inline NT & operator()(const int i) { assert(2 > i && 0 <= i); return data_[i]; }
    inline NT & x() { return data_[0]; }
    inline NT & y() { return data_[1]; }
    const NT * data() const { return data_; }

    inline NT squaredLength() const { return x()*x()+y()*y(); }

    inline Self operator-(const Self & rhs) const { return Self(x()-rhs.x(), y()-rhs.y()); }
    inline Self operator+(const Self & rhs) const { return Self(x()+rhs.x(), y()+rhs.y()); }
    inline NT operator|(const Self & rhs) const { return x()*rhs.x() + y()*rhs.y(); }
    inline void negate() { data_[0] = - data_[0]; data_[1] = - data_[1]; }
    inline void rotate() { NT temp(data_[0]); data_[0] = - data_[1]; data_[1] = temp; }

    NT data_[2];
};

template< typename NT >
Vec2<NT> operator*(const NT & lhs, const Vec2<NT> & rhs)
{
    return Vec2<NT>(lhs*rhs.x(), lhs*rhs.y());
}

template< typename NT >
Vec2<NT> vecmin(const Vec2<NT> & lhs, const Vec2<NT> & rhs)
{
    return Vec2<NT>(std::min(lhs.x(), rhs.x()), std::min(lhs.y(), rhs.y()));
}

template< typename NT >
Vec2<NT> vecmax(const Vec2<NT> & lhs, const Vec2<NT> & rhs)
{
    return Vec2<NT>(std::max(lhs.x(), rhs.x()), std::max(lhs.y(), rhs.y()));
}

template< typename NT >
bool operator==(const Vec2<NT> & lhs, const Vec2<NT> & rhs)
{
    return ( lhs.x() == rhs.x() ) && ( lhs.y() == rhs.y() );
}

template< typename NT >
bool operator!=(const Vec2<NT> & lhs, const Vec2<NT> & rhs)
{
    return ! ( lhs == rhs );
}

template< typename Out, typename NT >
Out & operator<<(Out & out, const Vec2<NT> & v)
{
    out << v.x() << ' ' << v.y();
    return out;
}

template< typename NT >
bool leftTurn(const Vec2<NT> & a, const Vec2<NT> & b, const Vec2<NT> & c)
{
    Vec2<NT> p(a.y()-b.y(), b.x()-a.x());
    return (p | (c-a)) >= NT(0);
}

extern const Vec2<float> vec2_zero;

/*
 *
 *  Vec3
 *
 */

template< typename NT >
struct Vec3
{
    typedef Vec3<NT>    Self;
    typedef NT          NumberType;

#define vecset(x,y,z) data_[0]=(x); data_[1]=(y); data_[2]=(z);
    Vec3() { vecset(NT(0), NT(0), NT(0)) }
    Vec3(const Self & s) { vecset(s.x(), s.y(), s.z()) }
    Vec3(const Vec2<NT> & v2, const NT & z) { vecset(v2.x(),v2.y(),z) }
    Vec3(const NT & x, const NT & y, const NT & z) { vecset(x,y,z) }
    inline void set(const NT & x, const NT & y, const NT & z) {
      vecset(x,y,z)
    }
#undef vecset

    inline NT operator()(const int i) const { assert(3 > i && 0 <= i); return data_[i]; }
    inline NT operator[](const int i) const { return data_[i]; }
    inline NT & operator[](const int i) { return data_[i]; }
    inline NT x() const { return data_[0]; }
    inline NT y() const { return data_[1]; }
    inline NT z() const { return data_[2]; }
    inline NT & operator()(const int i) { assert(3 > i && 0 <= i); return data_[i]; }
    inline NT & x() { return data_[0]; }
    inline NT & y() { return data_[1]; }
    inline NT & z() { return data_[2]; }
    const NT * data() const { return data_; }

    inline NT squaredLength() const { return x()*x()+y()*y()+z()*z(); }
    inline NT l1norm() const { return std::max(std::abs(x()), std::max(std::abs(y()), std::abs(z()))); }
    inline NT norm2() const { return x()*x()+y()*y()+z()*z(); }
    inline NT length() const { return ::sqrt(squaredLength()); }
    inline void normalize() { const NT l = length(); x() /= l; y() /= l; z() /= l; }
    inline Self normalized() const { const NT l = length(); return Self(x()/l, y()/l, z()/l); }
    inline void LInfNormalize() { const NT l = l1norm(); x() /= l; y() /= l; z() /= l; }
    inline Self LInfNormalized() const { const NT l = l1norm(); return Self(x()/l, y()/l, z()/l); }
    inline static Self cross(const Self & lhs, const Self & rhs) {
      return Self(lhs.y() * rhs.z() - lhs.z() * rhs.y(),
          lhs.z() * rhs.x() - lhs.x() * rhs.z(),
          lhs.x() * rhs.y() - lhs.y() * rhs.x());
    }

    inline Self operator-(const Self & rhs) const { return Self(x()-rhs.x(), y()-rhs.y(), z()-rhs.z()); }
    inline Self operator+(const Self & rhs) const { return Self(x()+rhs.x(), y()+rhs.y(), z()+rhs.z()); }
    inline NT operator|(const Self & rhs) const { return x()*rhs.x() + y()*rhs.y() + z()*rhs.z(); }
    inline void negate() { data_[0] = - data_[0]; data_[1] = - data_[1]; data_[2] = - data_[2]; }
    inline Self negated() const { return Self(-x(), -y(), -z()); }

    NT data_[3];
};

extern const Vec3<float> vec3_zero, unit_x, unit_y, unit_z;

template< typename NT >
Vec3<NT> operator/(const Vec3<NT> & lhs, const NT & rhs)
{
    return Vec3<NT>(lhs.x()/rhs, lhs.y()/rhs, lhs.z()/rhs);
}

template< typename NT >
Vec3<NT> operator*(const NT & lhs, const Vec3<NT> & rhs)
{
    return Vec3<NT>(lhs*rhs.x(), lhs*rhs.y(), lhs*rhs.z());
}

template< typename NT >
Vec3<NT> vecmin(const Vec3<NT> & lhs, const Vec3<NT> & rhs)
{
    return Vec3<NT>(std::min(lhs.x(), rhs.x()), std::min(lhs.y(), rhs.y()), std::min(lhs.z(), rhs.z()));
}

template< typename NT >
Vec3<NT> vecmax(const Vec3<NT> & lhs, const Vec3<NT> & rhs)
{
    return Vec3<NT>(std::max(lhs.x(), rhs.x()), std::max(lhs.y(), rhs.y()), std::max(lhs.z(), rhs.z()));
}

template< typename Out, typename NT >
Out & operator<<(Out & out, const Vec3<NT> & v)
{
    out << v.x() << ' ' << v.y() << ' ' << v.z();
    return out;
}

template< typename NT >
bool operator==(const Vec3<NT> & lhs, const Vec3<NT> & rhs)
{
    return ( lhs.x() == rhs.x() ) && ( lhs.y() == rhs.y() ) && ( lhs.z() == rhs.z() );
}

template< typename NT >
Vec3<NT> orthonormalVector(const Vec3<NT> & v) {
  Vec3<NT> r = ( v.x() * v.x() + v.y() * v.y() > 1e-4 ) ?
    Vec3<NT>(v.y(), -v.x(), NT(0))
    :
    Vec3<NT>(NT(0), v.z(), -v.y());
  r.normalize();
  return r;
}

/*
 *
 *  Vec4
 *
 */

template< typename NT >
struct alignas(16) Vec4
{
    typedef Vec4<NT>    Self;
    typedef Vec3<NT>    V3;
    typedef NT          NumberType;

#define vecset(x,y,z,w) data_[0]=(x); data_[1]=(y); data_[2]=(z); data_[3] = w;
    Vec4() { vecset(NT(0), NT(0), NT(0), NT(0)) }
    Vec4(const Self & s) { vecset(s.x(), s.y(), s.z(), s.w()) }
    Vec4(const V3 & s) { vecset(s.x(), s.y(), s.z(), NT(1)) }
    Vec4(const NT & x, const NT & y, const NT & z, const NT & w) { vecset(x,y,z,w) }
    inline void set(const NT & x, const NT & y, const NT & z, const NT & w) {
      vecset(x,y,z,w)
    }
#undef vecset

    inline NT operator()(const int i) const { assert(4 > i && 0 <= i); return data_[i]; }
    inline NT operator[](const int i) const { return data_[i]; }
    inline NT & operator[](const int i) { return data_[i]; }
    inline NT x() const { return data_[0]; }
    inline NT y() const { return data_[1]; }
    inline NT z() const { return data_[2]; }
    inline NT w() const { return data_[3]; }
    inline NT & operator()(const int i) { assert(4 > i && 0 <= i); return data_[i]; }
    inline NT & x() { return data_[0]; }
    inline NT & y() { return data_[1]; }
    inline NT & z() { return data_[2]; }
    inline NT & w() { return data_[3]; }
    const NT * data() const { return data_; }

    V3 toVec3() const { return V3(x(), y(), z()); }

    inline NT squaredLength() const { return x()*x()+y()*y()+z()*z()+w()*w(); }
    inline NT length() const { return ::sqrt(squaredLength()); }
    inline void normalize() { NT l = length(); x() /= l; y() /= l; z() /= l; w() /= l; }

    inline Self operator-(const Self & rhs) const { return Self(x()-rhs.x(), y()-rhs.y(), z()-rhs.z(), w()-rhs.w()); }
    inline Self operator+(const Self & rhs) const { return Self(x()+rhs.x(), y()+rhs.y(), z()+rhs.z(), w()+rhs.w()); }
    inline NT operator|(const Self & rhs) const { return x()*rhs.x() + y()*rhs.y() + z()*rhs.z() + w()*rhs.w(); }
    inline NT operator|(const V3 & rhs) const { return x()*rhs.x() + y()*rhs.y() + z()*rhs.z(); }
    inline NT dotAs3(const Self & rhs) const { return x()*rhs.x() + y()*rhs.y() + z()*rhs.z(); }
    //inline NT operator|(const V3 & rhs) const { return x()*rhs.x() + y()*rhs.y() + z()*rhs.z() + w(); }
    inline NT operator()(const V3 & rhs) const { return x()*rhs.x() + y()*rhs.y() + z()*rhs.z() + w(); }
    inline void negate() { data_[0] = - data_[0]; data_[1] = - data_[1]; data_[2] = - data_[2]; data_[3] = - data_[3]; }

    NT data_[4];
};

template< typename NT >
Vec4<NT> operator*(const NT & lhs, const Vec4<NT> & rhs)
{
    return Vec4<NT>(lhs*rhs.x(), lhs*rhs.y(), lhs*rhs.z(), lhs*rhs.w());
}

template< typename NT >
Vec4<NT> vecmin(const Vec4<NT> & lhs, const Vec4<NT> & rhs)
{
    return Vec4<NT>(std::min(lhs.x(), rhs.x()),
                    std::min(lhs.y(), rhs.y()),
                    std::min(lhs.z(), rhs.z()),
                    std::min(lhs.w(), rhs.w()));
}

template< typename NT >
Vec4<NT> vecmax(const Vec4<NT> & lhs, const Vec4<NT> & rhs)
{
    return Vec4<NT>(std::max(lhs.x(), rhs.x()),
                    std::max(lhs.y(), rhs.y()),
                    std::max(lhs.z(), rhs.z()),
                    std::max(lhs.w(), rhs.w()));
}

template< typename Out, typename NT >
Out & operator<<(Out & out, const Vec4<NT> & v)
{
    out << v.x() << ' ' << v.y() << ' ' << v.z() << ' ' << v.w();
    return out;
}

/*
 *
 *  Common vector types
 *
 */

typedef Vec2<double>         Vec2d;
typedef Vec2<float>          Vec2f;
typedef Vec2<int>            Vec2i;
typedef Vec2<unsigned int>   Vec2ui;
typedef Vec3<double>         Vec3d;
typedef Vec3<float>          Vec3f;
typedef Vec3<int>            Vec3i;
typedef Vec3<unsigned int>   Vec3ui;
typedef Vec4<double>         Vec4d;
typedef Vec4<float>          Vec4f;
typedef Vec4<int>            Vec4i;
typedef Vec4<unsigned int>   Vec4ui;
typedef Vec4<unsigned short> Vec4us;

#endif // ifndef __VECTOR_STUFF__
