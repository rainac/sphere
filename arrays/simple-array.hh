#ifndef sph_simplearray_hh
#define sph_simplearray_hh
/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

#ifdef SPH_USE_VALARRAY_FOR_PARTICLE_ARRAYS
#warning if SPH_USE_VALARRAY_FOR_PARTICLE_ARRAYS is defined, this header is NOT needed!
#endif

/// Helper function that delete[]-s the data (an array) and sets
/// the pointer to zero.
//
/// \param ptr reference to a pointer variable that is zeroed after
/// delete-ing the data pointed to.  
//
/// \tparam T the data type pointed to
template<class T>
inline void deleteAndClearV(T * &ptr) {
  if (ptr) {
    delete[] ptr;
    ptr = 0;
  }
}

#ifdef DEB_SIMPLE_ARRAY
#define DEB_SA(x) std::cerr << x << "\n";
#else
#define DEB_SA(x)
#endif

/// template class SimpleArray is a simplification of the
/// std::valarray template the reason to write this was that
/// std::valarray is implemented slightly differently on GNU and Sun
/// Studio compilers, so we need this class

template<class T>
struct SimpleArray {

  typedef T value_type;

  char *m_data;
  T    *m_base;
  size_t m_size;
  
  SimpleArray() : 
    m_data(),
    m_base(), 
    m_size() 
  {
    DEB_SA("construct SimpleArray<" << typeid(T).name()
	   << "> of size " << m_size);
  }

  explicit SimpleArray(size_t const n) : 
    m_data(),
    m_base(), 
    m_size() 
  {
    DEB_SA("construct SimpleArray<" << typeid(T).name()
	   << "> of size " << m_size);
    resize(n);
  }

  explicit SimpleArray(SimpleArray const &o) : 
    m_data(),
    m_base(), 
    m_size() 
  {
    DEB_SA("copy-construct SimpleArray<" << typeid(T).name()
	   << "> of size " << o.size());
    long const nsz = o.size();
    clear();
    m_data = new char[sizeof(T)*nsz];
    m_base = (T*) m_data;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(long i = 0; i < nsz; ++i) {
      new (m_base + i) T(o[i]);
    }
    m_size = nsz;
  }

  SimpleArray(T *ptr, size_t sz) : 
    m_data(),
    m_base(), 
    m_size() 
  {
    DEB_SA("copy-construct SimpleArray<" << typeid(T).name()
	   << "> of size " << o.size());
    long const nsz = sz;
    clear();
    m_data = new char[sizeof(T)*nsz];
    m_base = (T*) m_data;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(long i = 0; i < nsz; ++i) {
      new (m_base + i) T(*(ptr + i));
    }
    m_size = nsz;
  }

  SimpleArray(T const &val, size_t sz) : 
    m_data(),
    m_base(), 
    m_size() 
  {
    DEB_SA("copy-construct SimpleArray<" << typeid(T).name()
	   << "> of size " << o.size());
    long const nsz = sz;
    clear();
    m_data = new char[sizeof(T)*nsz];
    m_base = (T*) m_data;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(long i = 0; i < nsz; ++i) {
      new (m_base + i) T(val);
    }
    m_size = nsz;
  }

  ~SimpleArray() {
    clear();
  }

  /// Shortcut function to m_base.size().
  size_t size() const { return m_size; }

  /// non-const access operator shortcuts to m_base[].
  T       &operator[](size_t const i)       { return m_base[i]; }
  /// const access operator shortcuts to m_base[].
  T const &operator[](size_t const i) const { return m_base[i]; }

  void clear() {
    deleteAndClearV(m_data);
    m_base = 0;
    m_size = 0;
  }

  /// resize() simply resizes m_base.
  void resize(long const nsz) { 
    DEB_SA("resize SimpleArray<" << typeid(T).name()
	   << "> to size " << nsz);
    clear();
    m_data = new char[sizeof(T)*nsz];
    m_base = (T*) m_data;
    m_size = nsz;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(long i = 0; i < nsz; ++i) {
      new (m_base + i) T();
    }
  }

  T sum() const {
    T result = T();
    long const n = size();
#ifdef _OPENMP
#pragma omp parallel for reduction(+:result)
#endif
    for(long i = 0; i < n; ++i) {
      result += operator[](i);
    }
    return result;
  }

  SimpleArray &operator =(SimpleArray const &o) {
    DEB_SA("assign SimpleArray<" << typeid(T).name()
	   << "> of size " << m_size);
    assert(size() <= o.size());
    long const n = size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (long i = 0; i < n; ++i) {
      m_base[i] = o[i];
    }
    return *this;
  }

  SimpleArray &operator =(value_type const &o) {
    DEB_SA("assign SimpleArray<" << typeid(T).name()
	   << "> of size " << m_size);
    long const n = size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (long i = 0; i < n; ++i) {
      m_base[i] = o;
    }
    return *this;
  }

};

template<class T>
inline std::ostream &operator <<(std::ostream &aus, SimpleArray<T> const &v) {
  for(size_t i = 0; aus and i < v.size(); ++i) {
    if (i) aus << " ";
    aus << v[i];
  }
  return aus;
}

#endif
