#ifndef sph_permutation_hh
#define sph_permutation_hh
/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

#include <iostream>
#include <valarray>
#include <assert.h>

/// Template Permutation<T> represents a permutation of elements of
/// type T. The constructor takes the desired length len, allocates an
/// array of T of that length an fills it with the identy permuation,
/// i.e. (0, 1, 2, 3, 4, ...)

template<class T=unsigned long>
struct Permutation {
  typedef T value_type;

  size_t m_n;
  std::valarray<T> m_p;

  Permutation(size_t const _n) :
    m_n(_n),
    m_p(_n)
  {
    fill();
  }

  void resize(size_t const _n) { 
    m_n = _n;
    m_p.resize(_n);
    fill();
  }

  void fill() { 
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for(long i = 0; i < long(m_n); ++i) {
      m_p[i] = i;
    }
  }

  size_t   size()  const { return m_p.size(); }

  T       *begin()       { return &m_p[0]; }
  T const *begin() const { return &m_p[0]; }

  T       *end()         { return &m_p[m_n]; }
  T const *end()   const { return &m_p[m_n]; }

  T &operator[] (T i)       { return m_p[i]; }
  T  operator[] (T i) const { return m_p[i]; }

  bool valid() const {
    std::valarray<bool> seen(false, m_n);
    for (T i = 0; i < m_n; ++i) {
      seen[m_p[i]] = 1;
    }
    return seen.min();
  }

  void print(std::ostream &aus) const {
    aus << "(";
    for (T i = 0; i < m_n; ++i) {
      aus << m_p[i] << ", ";
    }    
    aus << ")";
  }
  
};
template<class T>
inline std::ostream &operator <<(std::ostream &aus, Permutation<T> const &v) {
  v.print(aus);
  return aus;
}

#endif
