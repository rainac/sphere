#ifndef static_array_hh
#define static_array_hh
/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

#include <iostream>
#include <cassert>
#include <algorithm>
#include <math.h>

#define point_size_t unsigned

template<class T, point_size_t N>
struct Point {
  typedef T value_type;

  static point_size_t const m_size = N;
//   static point_size_t const m_sizem1 = m_size-1;
//   static point_size_t const m_esize = sizeof(T);
        
  T data[N];
//   bool m_transposed;

  explicit Point(point_size_t = 0) {
    for (point_size_t i = 0; i < N; ++i) {
      data[i] = T();
    }
  }

  Point(Point const &v) {
    for (point_size_t i = 0; i < N; ++i) {
      data[i] = v.data[i];
    }
  }

  // create Point from Point of same size N but with different base type U
  template<class U>
  Point(Point<U, N> const &v) {
    for (point_size_t i = 0; i < N; ++i) {
      data[i] = T(v.data[i]);
    }
  }

  Point(T const &v, point_size_t) {
    for (point_size_t i = 0; i < N; ++i) {
      data[i] = v;
    }
  }

//   Point(T const &v, T const &w) {
//     data[0] = v;
//     data[1] = w;
//     for (point_size_t i = 2; i < N; ++i) {
//       data[i] = T();
//     }
//   }

//   Point(T const &u, T const &v, T const &w) {
//     data[0] = u;
//     data[1] = v;
//     data[2] = w;
//     for (point_size_t i = 3; i < N; ++i) {
//       data[i] = T();
//     }
//   }

  Point(T const *v, point_size_t const n) {
    for (point_size_t i = 0; i < n; ++i) {
      data[i] = v[i];
    }
    for (point_size_t i = n; i < N; ++i) {
      data[i] = T();
    }
  }

  template<template <class> class V, class U>
  explicit Point(V<U> const &v) {
//     std::valarray<U> tmp(v);
//     assert(tmp.size() >= N);
    for (point_size_t i = 0; i < N; ++i) {
      data[i] = T(v[i]);
    }
  }

//   ~Point() { }

  point_size_t size() const { return m_size; }

  T normSquared() const {
    T res = T();
    for (point_size_t i = 0; i < N; ++i) {
      res += data[i]*data[i];
    }
    return res;
  }

  T norm() const {
    return sqrt(normSquared());
  }
//   T abs() const {
//     return norm();
//   }

  inline T       &operator[](point_size_t const i)       { return data[i]; }
  inline T const &operator[](point_size_t const i) const { return data[i]; }

  template<class U>
  Point<bool, N> operator ==(Point<U, N> const &o) const {
    Point<bool, N> res;
    for (point_size_t i = 0; i < N; ++i) {
      res[i] = data[i] == o.data[i];
    }
    return res;
  }

  template<class U>
  Point<bool, N> operator !=(Point<U, N> const &o) const {
    Point<bool, N> res;
    for (point_size_t i = 0; i < N; ++i) {
      res[i] = data[i] != o.data[i];
    }
    return res;
  }

  template<class U>
  Point<bool, N> operator <(Point<U, N> const &o) const {
    Point<bool, N> res;
    for (point_size_t i = 0; i < N; ++i) {
      res[i] = data[i] < o.data[i];
    }
    return res;
  }

  template<class U>
  Point<bool, N> operator >(Point<U, N> const &o) const {
    Point<bool, N> res;
    for (point_size_t i = 0; i < N; ++i) {
      res[i] = data[i] > o.data[i];
    }
    return res;
  }

  template<class U>
  Point<bool, N> operator >=(Point<U, N> const &o) const {
    Point<bool, N> res;
    for (point_size_t i = 0; i < N; ++i) {
      res[i] = data[i] >= o.data[i];
    }
    return res;
  }

  template<class U>
  Point<bool, N> operator <=(Point<U, N> const &o) const {
    Point<bool, N> res;
    for (point_size_t i = 0; i < N; ++i) {
      res[i] = data[i] <= o.data[i];
    }
    return res;
  }

  template<class U>
  Point<bool, N> operator ==(U const &o) const {
    Point<bool, N> res;
    for (point_size_t i = 0; i < N; ++i) {
      res[i] = data[i] == o;
    }
    return res;
  }

  template<class U>
  Point<bool, N> operator !=(U const &o) const {
    Point<bool, N> res;
    for (point_size_t i = 0; i < N; ++i) {
      res[i] = data[i] != o;
    }
    return res;
  }

  template<class U>
  Point<bool, N> operator <(U const &o) const {
    Point<bool, N> res;
    for (point_size_t i = 0; i < N; ++i) {
      res[i] = data[i] < o;
    }
    return res;
  }

  template<class U>
  Point<bool, N> operator >(U const &o) const {
    Point<bool, N> res;
    for (point_size_t i = 0; i < N; ++i) {
      res[i] = data[i] > o;
    }
    return res;
  }

  template<class U>
  Point<bool, N> operator >=(U const &o) const {
    Point<bool, N> res;
    for (point_size_t i = 0; i < N; ++i) {
      res[i] = data[i] >= o;
    }
    return res;
  }

  template<class U>
  Point<bool, N> operator <=(U const &o) const {
    Point<bool, N> res;
    for (point_size_t i = 0; i < N; ++i) {
      res[i] = data[i] <= o;
    }
    return res;
  }

  Point &operator =(value_type const &o) {
    for (point_size_t i = 0; i < N; ++i) {
      data[i] = o;
    }
    return *this;
  }

  template<class U>
  Point &operator =(Point<U, N> const &o) {
    for (point_size_t i = 0; i < N; ++i) {
      data[i] = T(o.data[i]);
    }
    return *this;
  }

  template<template <class> class V, class U>
  Point &operator =(V<U> const &o) {
    for (point_size_t i = 0; i < N; ++i) {
      data[i] = T(o[i]);
    }
    return *this;
  }

  template<class U>
  Point &operator +=(Point<U, N> const &o) {
    for (point_size_t i = 0; i < N; ++i) {
      data[i] += o.data[i];
    }
    return *this;
  }

  template<class U>
  Point &operator -=(Point<U, N> const &o) {
    for (point_size_t i = 0; i < N; ++i) {
      data[i] -= o.data[i];
    }
    return *this;
  }

  template<class U>
  Point &operator *=(Point<U, N> const &o) {
    for (point_size_t i = 0; i < N; ++i) {
      data[i] *= o.data[i];
    }
    return *this;
  }

  template<class U>
  Point &operator /=(Point<U, N> const &o) {
    for (point_size_t i = 0; i < N; ++i) {
      data[i] /= o.data[i];
    }
    return *this;
  }

  template<class U>
  Point &operator %=(Point<U, N> const &o) {
    for (point_size_t i = 0; i < N; ++i) {
      data[i] %= o.data[i];
    }
    return *this;
  }

  template<class U>
  Point &operator <<=(Point<U, N> const &o) {
    for (point_size_t i = 0; i < N; ++i) {
      data[i] <<= o.data[i];
    }
    return *this;
  }

  template<class U>
  Point &operator >>=(Point<U, N> const &o) {
    for (point_size_t i = 0; i < N; ++i) {
      data[i] >>= o.data[i];
    }
    return *this;
  }

  Point &operator +=(T const &o) {
    for (point_size_t i = 0; i < N; ++i) {
      data[i] += o;
    }
    return *this;
  }

  Point &operator -=(T const &o) {
    for (point_size_t i = 0; i < N; ++i) {
      data[i] -= o;
    }
    return *this;
  }

  Point &operator *=(T const &o) {
    for (point_size_t i = 0; i < N; ++i) {
      data[i] *= o;
    }
    return *this;
  }

  Point &operator /=(T const &o) {
    for (point_size_t i = 0; i < N; ++i) {
      data[i] /= o;
    }
    return *this;
  }

  Point &operator &=(T const &o) {
    for (point_size_t i = 0; i < N; ++i) {
      data[i] &= o;
    }
    return *this;
  }

  Point &operator |=(T const &o) {
    for (point_size_t i = 0; i < N; ++i) {
      data[i] |= o;
    }
    return *this;
  }

  Point &operator ^=(T o) {
    for (point_size_t i = 0; i < N; ++i) {
      data[i] ^= o;
    }
    return *this;
  }

  Point &operator <<=(T const &o) {
    for (point_size_t i = 0; i < N; ++i) {
      data[i] <<= o;
    }
    return *this;
  }

  Point &operator >>=(T const &o) {
    for (point_size_t i = 0; i < N; ++i) {
      data[i] >>= o;
    }
    return *this;
  }

  Point operator +() const {
    Point res;
    for (point_size_t i = 0; i < N; ++i) {
      res.data[i] = + data[i];
    }
    return res;
  }

  Point &powassign(T const v) {
    for (point_size_t i = 0; i < N; ++i) {
      data[i] = pow(data[i], v);
    }
    return *this;
  }

  Point &powassign(Point const &p) {
    for (point_size_t i = 0; i < N; ++i) {
      data[i] = pow(data[i], p.data[i]);
    }
    return *this;
  }

  Point operator -() const {
    Point res;
    for (point_size_t i = 0; i < N; ++i) {
      res.data[i] = - data[i];
    }
    return res;
  }

  Point abs() const {
    Point res;
    for (point_size_t i = 0; i < N; ++i) {
      res.data[i] = (data[i] < 0 ? -data[i] : data[i]);
    }
    return res;
  }

  Point sign() const {
    Point res;
    for (point_size_t i = 0; i < N; ++i) {
      res.data[i] = (data[i] > 0) - (data[i] < 0);
    }
    return res;
  }

#define Point_unary_const_mem(name)   \
  Point name() const {                \
    Point res;                        \
    for (point_size_t i = 0; i < N; ++i) {  \
      res.data[i] = ::name(data[i]);  \
    }                                 \
    return res;                       \
  }

//   Point_unary_const_mem(abs)
  Point_unary_const_mem(ceil)
  Point_unary_const_mem(round)
  Point_unary_const_mem(floor)

  Point_unary_const_mem(sin)
  Point_unary_const_mem(cos)
  Point_unary_const_mem(tan)

  Point_unary_const_mem(asin)
  Point_unary_const_mem(acos)
  Point_unary_const_mem(atan)

  Point_unary_const_mem(sinh)
  Point_unary_const_mem(cosh)
  Point_unary_const_mem(tanh)

  Point_unary_const_mem(asinh)
  Point_unary_const_mem(acosh)
  Point_unary_const_mem(atanh)

  Point_unary_const_mem(exp)
  Point_unary_const_mem(log)
  Point_unary_const_mem(log10)

  T max() const {
    T res = data[0];
    for (point_size_t i = 1; i < N; ++i) {
      res = std::max(res, data[i]);
    }
    return res;
  }

  T min() const {
    T res = data[0];
    for (point_size_t i = 1; i < N; ++i) {
      res = std::min(res, data[i]);
    }
    return res;
  }

  template<class T2>
  Point min(Point<T2, N> const &o) const {
    Point res;
    for (point_size_t i = 0; i < N; ++i) {
      res[i] = min(data[i], o[i]);
    }
    return res;
  }

  template<class T2>
  Point max(Point<T2, N> const &o) const {
    Point res;
    for (point_size_t i = 0; i < N; ++i) {
      res[i] = std::max(data[i], o[i]);
    }
    return res;
  }

  template<class T2>
  T scalar(Point<T2, N> const &o) const {
    T res = T();
    for (point_size_t i = 0; i < N; ++i) {
      res += data[i] * o[i];
    }
    return res;
  }

  T sum() const {
    T res = T();
    for (point_size_t i = 0; i < N; ++i) {
      res += data[i];
    }
    return res;
  }

  T prod(T initial = T(1)) const {
    T res = initial;
    for (point_size_t i = 0; i < N; ++i) {
      res *= data[i];
    }
    return res;
  }

  T binOr() const {
    T res = T();
    for (point_size_t i = 0; i < N; ++i) {
      res |= data[i];
    }
    return res;
  }

//   void transpose() {
//     m_transposed = !m_transposed;
//   }

  void printVert(std::ostream &aus) const {
    for (point_size_t i = 0; i < N; ++i) {
      aus << data[i] << "\n";
    }
  }

  void printHoriz(std::ostream &aus) const {
    for (point_size_t i = 0; i < N; ++i) {
      if (i > 0) aus << ", ";
      aus << data[i];
    }
  }

  void print(std::ostream &aus) const {
//     if (m_transposed) {
//       printVert(aus);
//     } else {
      printHoriz(aus);
//     }
  }

  void prettyPrint(std::ostream &aus) const {
    aus << "(";
    for (point_size_t i = 0; i < N; ++i) {
      if (i > 0) aus << ", ";
      aus << data[i];
    }
    aus << ")";
  }

  void writeXML(std::ostream &aus) const {
    for (point_size_t i = 0; i < N; ++i) {
      aus << "<coord index='" << i << "'>"<< data[i] << "</coord>\n";
    }
  }

  void read(std::istream &ein) {
    std::string comma;
    for (point_size_t i = 0; i < N; ++i) {
      if (i > 0 and ein.peek() == ',') 
        ein >> comma;
      ein >> data[i];
    }
  }
};

template<class T, point_size_t N>
inline std::ostream &operator <<(std::ostream &aus, Point<T, N> const &v) {
  v.print(aus);
  return aus;
}

template<class T, point_size_t N>
inline std::istream &operator >>(std::istream &aus, Point<T, N> &v) {
  v.read(aus);
  return aus;
}

#define Point_unary_function(name)                     \
  template<class T, point_size_t N>                                 \
  inline Point<T,N> name(Point<T,N> const &a) {\
    return a.name();						    \
  }                                                                 

Point_unary_function(sin)
Point_unary_function(cos)
Point_unary_function(tan)

Point_unary_function(abs)
Point_unary_function(ceil)
Point_unary_function(floor)
Point_unary_function(round)
Point_unary_function(sign)

Point_unary_function(asin)
Point_unary_function(acos)
Point_unary_function(atan)

Point_unary_function(sinh)
Point_unary_function(cosh)
Point_unary_function(tanh)

Point_unary_function(asinh)
Point_unary_function(acosh)
Point_unary_function(atanh)

Point_unary_function(exp)
Point_unary_function(log)
Point_unary_function(log10)

template<class T, point_size_t N>						
inline Point<T,N> pow(Point<T,N> const &a, double v) {
  Point<T,N> res(a);
  res.powassign(v);
  return res;
}

template<class T, point_size_t N>						
inline Point<T,N> pow(Point<T,N> const &a, Point<T,N> const &b) {
  Point<T,N> res(a);
  res.powassign(b);
  return res;
}

#define Point_Point_binary_operator(name, nameass)                     \
  template<class T, point_size_t N>                                 \
  inline Point<T,N> operator name(Point<T,N> const &a, Point<T,N> const &b) {\
    Point<T,N> res(a);                                                 \
    res nameass b;                                                  \
    return res;                                                     \
  }                                                                 \
                                                                    \
  template<class T, point_size_t N>                                          \
  inline Point<T,N> operator name(Point<T,N> const &a, T const &b) {      \
    Point<T,N> res(a);                                                 \
    res nameass b;                                                  \
    return res;                                                     \
  }                                                                 \
                                                                    \
  template<class T, point_size_t N>                                          \
  inline Point<T,N> operator name(T const &a, Point<T,N> const &b) {      \
    Point<T,N> res(a, N);                                                \
    res nameass b;                                                  \
    return res;                                                     \
  }                                                                 \

Point_Point_binary_operator(+, +=)
Point_Point_binary_operator(-, -=)
Point_Point_binary_operator(*, *=)
Point_Point_binary_operator(/, /=)

Point_Point_binary_operator(&, &=)
Point_Point_binary_operator(|, |=)
Point_Point_binary_operator(^, ^=)

Point_Point_binary_operator(<<, <<=)
Point_Point_binary_operator(>>, >>=)

template<class T, point_size_t N>			
inline T norm(Point<T,N> const &a) { return a.norm(); }

template<class T, point_size_t N>			
inline Point<T,N> normed(Point<T,N> const &a) { return a / a.norm(); }

template<class T, point_size_t N>			
inline T normSquared(Point<T,N> const &a) { return a.normSquared(); }

template<class T1, class T2, point_size_t N>			
inline T1 scalar(Point<T1,N> const &a, Point<T2,N> const &b) { 
  return a.scalar(b); 
}

#endif
