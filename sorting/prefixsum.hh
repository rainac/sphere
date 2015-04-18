#ifndef jw_prefixsum34561_hh
#define jw_prefixsum34561_hh
/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

#include "histogram.hh"
#include <algorithm>
#include <valarray>
#include <assert.h>

template<class TIterator, class ValueGetter /*=IDGetter<typename GetPointee<TIterator>::value_type>*/ , 
	 template<class> class Vector=std::valarray,
	 class Counter=unsigned>
struct ParPrefixSumsViaHistogram {

  typedef typename GetPointee<TIterator>::value_type value_type;
  typedef Vector<Counter> vector_type;
  typedef Counter counter_type;

  typedef ParHistogram<TIterator, ValueGetter, Vector, Counter> Histogram;

  Histogram m_histogram;
  ValueGetter &m_valueGetter;

  ParPrefixSumsViaHistogram(size_t m_valueRange) : 
    m_histogram(m_valueRange),
    m_valueGetter(m_histogram.m_valueGetter)
  {}

  ParPrefixSumsViaHistogram(size_t m_valueRange, ValueGetter valueGetter) : 
    m_histogram(m_valueRange, valueGetter),
    m_valueGetter(m_histogram.m_valueGetter)
  {}

  size_t size() const { return m_histogram.size(); }

  Counter       &operator[](size_t i)       { return m_histogram[i]; }
  Counter const &operator[](size_t i) const { return m_histogram[i]; }

  void sample(TIterator const beg, TIterator const end, vector_type &_counts) {
    m_histogram.sample(beg, end, _counts);

    vector_type totals;

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
      size_t const myid = omp_get_thread_num();
      size_t const nthreads = omp_get_num_threads();

      size_t const myRangeM = ceil(double(size()) / nthreads);
      size_t const myLowM = std::min(myid * myRangeM, size());
      size_t const myHighM = std::min((myid + 1) * myRangeM, size());

      size_t lastValue = 0, curValue = 0;

      if (myid == 0) {
	totals.resize(nthreads);
	assert(totals.sum() == 0);
      }

#ifdef _OPENMP
#pragma omp barrier
#endif

      for(size_t i = myLowM; i < myHighM; ++i) {
	curValue = _counts[i];
	_counts[i] = lastValue;
	lastValue += curValue;
      }

      totals[myid] = lastValue;

#ifdef _OPENMP
#pragma omp barrier
#endif
      if (myid == 0) {
	for(size_t i = 1; i < nthreads; ++i) {
	  totals[i] += totals[i - 1];
	}
      }

#ifdef _OPENMP
#pragma omp barrier
#endif
      if (myid > 0) {
	size_t const myTotal = totals[myid-1];
	for(size_t i = myLowM; i < myHighM; ++i) {
	  _counts[i] += myTotal;
	}
      }

    } // end parallel region

  }
 
};


template<class DataVector, class ValueGetter /*=IDGetter<typename GetPointee<TIterator>::value_type>*/ , 
	 template<class> class Vector=std::valarray,
	 class Counter=unsigned >
struct ParPrefixSums {

  typedef typename ValueGetter::value_type value_type;
  typedef Counter counter_type;
  typedef Vector<Counter> vector_type;

  size_t const m_valueRange;
//   vector_type m_counts;
  ValueGetter m_valueGetter;

  ParPrefixSums(size_t m_valueRange) : 
    m_valueRange(m_valueRange)
  {}

  ParPrefixSums(size_t m_valueRange, ValueGetter valueGetter) : 
    m_valueRange(m_valueRange) 
    , m_valueGetter(valueGetter)
  {}

  size_t size() const { return m_valueRange; }

//   Counter       &operator[](size_t i)       { return m_counts[i]; }
//   Counter const &operator[](size_t i) const { return m_counts[i]; }

  void sample(DataVector const &data, vector_type &_counts) {

    long const N = data.size();

    vector_type *localCounts = 0;
    vector_type totals;

    _counts = typename vector_type::value_type();

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
      size_t const myid = omp_get_thread_num();
      size_t const nthreads = omp_get_num_threads();
      vector_type *myLocalCount = 0;

      size_t const myRangeN = ceil(double(N) / nthreads);
      size_t const myLowN = std::min(long(myid * myRangeN), N);
      size_t const myHighN = std::min(long((myid + 1) * myRangeN), N);

      size_t const myRangeM = ceil(double(m_valueRange) / nthreads);
      size_t const myLowM = std::min(myid * myRangeM, m_valueRange);
      size_t const myHighM = std::min((myid + 1) * myRangeM, m_valueRange);

      if (myid == 0) {
	localCounts = new vector_type[nthreads];
	totals.resize(nthreads);
	assert(totals.sum() == 0);
      }

#ifdef _OPENMP
#pragma omp barrier
#endif
      myLocalCount = localCounts + myid;
      myLocalCount->resize(m_valueRange);
      assert(myLocalCount->sum() == 0);
#ifdef _OPENMP
#pragma omp barrier
#endif

      for(size_t i = myLowN; i < myHighN; ++i) {
	value_type v = m_valueGetter(data[i]);
	++(*myLocalCount)[v];
      }

#ifdef _OPENMP
#pragma omp barrier
#endif

      size_t lastTotal = 0;
      
      for(size_t j = 0; j < nthreads; ++j) {
	for(size_t i = myLowM + 1; i < myHighM; ++i) {
	  localCounts[j][i] += localCounts[j][i - 1];
	  _counts[i] += localCounts[j][i - 1];
	}
	lastTotal += localCounts[j][myHighM - 1];
      }

      totals[myid] = lastTotal;

#ifdef _OPENMP
#pragma omp barrier
#endif
      if (myid > 0) {
	for(size_t i = myLowM; i < myHighM; ++i) {
	  for(size_t j = 0; j < myid; ++j) {
	    _counts[i] += totals[j];
	  }
	}
      }

    } // end parallel region

    delete[] localCounts;
    localCounts = 0;
  }
 
};


template<class TIterator, class ValueGetter /*=IDGetter<typename GetPointee<TIterator>::value_type>*/ , 
	 template<class> class Vector=std::valarray,
	 class Counter=unsigned>
struct SerialPrefixSumsViaHistogram {

  typedef typename GetPointee<TIterator>::value_type value_type;
  typedef Vector<Counter> vector_type;
  typedef Counter counter_type;

  typedef SerialHistogram<TIterator, ValueGetter, Vector, Counter> Histogram;

  Histogram m_histogram;
  ValueGetter &m_valueGetter;

  SerialPrefixSumsViaHistogram(size_t m_valueRange) : 
    m_histogram(m_valueRange),
    m_valueGetter(m_histogram.m_valueGetter)
  {}
  
  SerialPrefixSumsViaHistogram(size_t m_valueRange, ValueGetter valueGetter) : 
    m_histogram(m_valueRange, valueGetter),
    m_valueGetter(m_histogram.m_valueGetter)
  {}

  size_t size() const { return m_histogram.size(); }

  Counter       &operator[](size_t i)       { return m_histogram[i]; }
  Counter const &operator[](size_t i) const { return m_histogram[i]; }

  void sample(TIterator const beg, TIterator const end, vector_type &_counts) {
    m_histogram.sample(beg, end, _counts);

    vector_type totals;
    size_t const M = size();
    size_t lastValue = 0, curValue = 0;

    for(size_t i = 0; i < M; ++i) {
      curValue = _counts[i];
      _counts[i] = lastValue;
      lastValue += curValue;
    }

  }
 
};

#endif
