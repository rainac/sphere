#ifndef jw_histogram2343_hh
#define jw_histogram2343_hh
/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

#include "getpointee.hh"
#include <valarray>
/*
template<class Value>
struct IDGetter {
  typedef Value arg_type;
  typedef Value value_type;
  value_type operator()(arg_type v) const { return v; }
};
*/

template<class TIterator, class ValueGetter /*=IDGetter<typename GetPointee<TIterator>::value_type>*/ ,
	 template<class> class Vector=std::valarray,
	 class Counter=unsigned>
struct ParHistogram {

  typedef typename GetPointee<TIterator>::value_type it_value_type;
  typedef typename ValueGetter::value_type value_type;
  typedef Vector<Counter> vector_type;

  size_t const m_valueRange;
//   vector_type m_counts;
  ValueGetter m_valueGetter;

  ParHistogram(size_t const _valueRange) : 
    m_valueRange(_valueRange)
  {}

  ParHistogram(size_t const _valueRange, ValueGetter valueGetter) : 
    m_valueRange(_valueRange)
    , m_valueGetter(valueGetter)
  {}

  size_t size() const { return m_valueRange; }

//   Counter       &operator[](size_t i)       { return m_counts[i]; }
//   Counter const &operator[](size_t i) const { return m_counts[i]; }

  void sample(TIterator const beg, TIterator const end, vector_type &_counts) {

    long const N = end - beg;

    vector_type *localCounts = 0;

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
      size_t const myid = omp_get_thread_num();
      size_t const nthreads = omp_get_num_threads();
      vector_type *myLocalCount = 0;

      if (myid == 0) {
	localCounts = new vector_type[nthreads];
      }

#ifdef _OPENMP
#pragma omp barrier
#endif
      myLocalCount = localCounts + myid;
      myLocalCount->resize(m_valueRange);

#ifdef _OPENMP
#pragma omp for
#endif
      for(long i = 0; i < N; ++i) {
	value_type v = m_valueGetter(*(beg + i));
	++(*myLocalCount)[v];
      }

#ifdef _OPENMP
#pragma omp for
#endif
      // histogram2: switched these for loops, seems to help some what
      for(long i = 0; i < long(m_valueRange); ++i) {
	for(size_t j = 0; j < nthreads; ++j) {
	  _counts[i] += localCounts[j][i];
	}
      }

      if (myid == 0) {
	delete[] localCounts;
	localCounts = 0;
      }
    } // end parallel region
  }
};


template<class TIterator, class ValueGetter /*=IDGetter<typename GetPointee<TIterator>::value_type>*/ ,
	 template<class> class Vector=std::valarray,
	 class Counter=unsigned>
struct ParHistogramOld {

  typedef typename GetPointee<TIterator>::value_type it_value_type;
  typedef typename ValueGetter::value_type value_type;
  typedef Vector<Counter> vector_type;

  size_t const m_valueRange;
//   vector_type m_counts;
  ValueGetter m_valueGetter;

  ParHistogramOld(size_t const _valueRange) : 
    m_valueRange(_valueRange)
  {}

  ParHistogramOld(size_t const _valueRange, ValueGetter valueGetter) : 
    m_valueRange(_valueRange)
    , m_valueGetter(valueGetter)
  {}

  size_t size() const { return m_valueRange; }

//   Counter       &operator[](size_t i)       { return m_counts[i]; }
//   Counter const &operator[](size_t i) const { return m_counts[i]; }

  void sample(TIterator const beg, TIterator const end, vector_type &_counts) {

    long const N = end - beg;

    vector_type *localCounts = 0;

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
      size_t const myid = omp_get_thread_num();
      size_t const nthreads = omp_get_num_threads();
      vector_type *myLocalCount = 0;

      if (myid == 0) {
	localCounts = new vector_type[nthreads];
      }

#ifdef _OPENMP
#pragma omp barrier
#endif
      myLocalCount = localCounts + myid;
      myLocalCount->resize(m_valueRange);

#ifdef _OPENMP
#pragma omp for
#endif
      for(long i = 0; i < N; ++i) {
	value_type v = m_valueGetter(*(beg + i));
	++(*myLocalCount)[v];
      }

      for(size_t j = 0; j < nthreads; ++j) {
#ifdef _OPENMP
#pragma omp for
#endif
	for(long i = 0; i < long(m_valueRange); ++i) {
	  _counts[i] += localCounts[j][i];
	}
      }

      if (myid == 0) {
	delete[] localCounts;
	localCounts = 0;
      }
    } // end parallel region
  }
};

template<class TIterator, class ValueGetter /*=IDGetter<typename GetPointee<TIterator>::value_type>*/ ,
	 template<class> class Vector=std::valarray,
	 class Counter=unsigned>
struct SerialHistogram {

  typedef typename GetPointee<TIterator>::value_type it_value_type;
  typedef typename ValueGetter::value_type value_type;
  typedef Vector<Counter> vector_type;

  size_t const m_valueRange;
//   vector_type m_counts;
  ValueGetter m_valueGetter;

  SerialHistogram(size_t const _valueRange) : 
    m_valueRange(_valueRange)
  {}

  SerialHistogram(size_t const _valueRange, ValueGetter valueGetter) : 
    m_valueRange(_valueRange)
    , m_valueGetter(valueGetter)
  {}

  size_t size() const { return m_valueRange; }

//   Counter       &operator[](size_t i)       { return m_counts[i]; }
//   Counter const &operator[](size_t i) const { return m_counts[i]; }

  void sample(TIterator const beg, TIterator const end, vector_type &_counts) {

    long const N = end - beg;

    for(long i = 0; i < N; ++i) {
      value_type v = m_valueGetter(*(beg + i));
      ++_counts[v];
    }

  }
};

#endif
