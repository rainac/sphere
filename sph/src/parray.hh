#ifndef sph_paartarray_hh
#define sph_paartarray_hh
/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

#ifdef SPH_USE_VALARRAY_FOR_PARTICLE_ARRAYS
#warning if SPH_USE_VALARRAY_FOR_PARTICLE_ARRAYS is defined, this header is NOT needed!
#endif

#include "arrays/simple-array.hh"

template<class PArray>
struct PSliceArray {
  typedef typename PArray::value_type value_type;

  PArray &a;
  std::slice const &sl;

  PSliceArray(PArray &a, std::slice const &sl) :
    a(a),
    sl(sl)
  {}

  template<template <class> class V, class T>
  PSliceArray &operator =(V<T> const &v) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (long cnt = 0; cnt < long(sl.size()); ++cnt) {
      a[sl.start() + cnt * sl.stride()] = v[cnt];
    }
    return *this;
  }
};


/// struct PArray is a very simple wrapper around an array of Particle or ParticleVariable
/// 
template<class T>
struct PArray {
  typedef T value_type;
  typedef PArray my_type;

  SimpleArray<value_type> m_data;

  PArray() {}

  explicit PArray(size_t const sz) :
    m_data(sz) 
  {}
//   PartikelArray() {}
//   explicit PartikelArray(PartikelPtrArray const &v) : 
//     m_data(v.size())
//   {
//     *this = v;
//   }

//   PartikelArray &operator =(PartikelPtrArray const &v) {
//     size_t const n = size();
//     for (size_t i = 0; i < n; ++i) {
//       m_data[i] = *v[i];
//     }
//     return *this;
//   }


  PArray &operator +=(PArray const &v) {
    long const n = size();
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (long i = 0; i < n; ++i) {
      m_data[i] += v[i];
    }
    return *this;
  }

  PArray &operator *=(value_type const v) {
    long const n = size();
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (long i = 0; i < n; ++i) {
      m_data[i] *= v;
    }
    return *this;
  }

  PArray &operator =(value_type const v) {
    long const n = size();
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (long i = 0; i < n; ++i) {
      m_data[i] = v;
    }
    return *this;
  }

  void resize(size_t nsz) { 
    m_data.resize(nsz); 
  }
  size_t size() const { return m_data.size(); }
  value_type &operator[](size_t i) { return m_data[i]; }
  value_type const &operator[](size_t i) const { return m_data[i]; }

  PSliceArray<PArray> operator[](std::slice const &sl) { 
    return PSliceArray<PArray>(*this, sl); 
  }

};

template<class T>
inline PArray<T> operator *(PArray<T> const &v, T const &u) {
  PArray<T> res(v);
  res *= u;
  return res;
}

template<class T>
inline std::ostream &operator <<(std::ostream &aus, PArray<T> const &p) {
  aus << p.m_data;
  return aus;
}

#endif
