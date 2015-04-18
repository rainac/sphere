/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

#include <omp.h>
#include <valarray>
#include <string.h>

template<class It>
struct SPHParCSDefValueGetter {
  typedef typename It::value_type value_type;
  typedef typename It::value_type argument_type;
};

template<class It>
struct SPHParCSDefValueGetter<It*> {
  typedef It value_type;
  typedef It argument_type;
  value_type operator()(argument_type const &v) const { return v; }
};

template<class TIterator, class ValueGetter=SPHParCSDefValueGetter<TIterator> >
struct SPHParCountingSorter {

  typedef typename ValueGetter::value_type value_type;
  typedef typename ValueGetter::argument_type data_type;
//   typedef typename GetPointee<TIterator>::value_type it_value_type;

  typedef std::valarray<value_type> vector_type;

  size_t const m_valueRange;
  size_t const m_nthreads;
  ValueGetter m_valueGetter;

  explicit SPHParCountingSorter(unsigned valueRange, unsigned nthreads) : 
    m_valueRange(valueRange),
    m_nthreads(nthreads)
  {}
  explicit SPHParCountingSorter(unsigned valueRange, unsigned nthreads, ValueGetter &valueGetter) : 
    m_valueRange(valueRange)
    , m_nthreads(nthreads)
    , m_valueGetter(valueGetter)
  {}

  ~SPHParCountingSorter() { }

  template<class PrefIterator>
  void sort(TIterator beg, TIterator end, PrefIterator prefixSumBegin) {

    long const N = end - beg;
    
    value_type **localCounts = 0;
    value_type *totals = 0;

    data_type *tmpData = new data_type[N];

    size_t const myRangeN = ceil(double(N) / m_nthreads);
    size_t const myRangeM = ceil(double(m_valueRange) / m_nthreads);
    
    localCounts = new value_type*[m_nthreads];

    totals = new value_type[m_valueRange];

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
      size_t const myid = omp_get_thread_num();

      size_t const myLowN = std::min(long(myid * myRangeN), N);
      size_t const myHighN = std::min(long((myid + 1) * myRangeN), N);
      
      size_t const myLowM = std::min(myid * myRangeM, m_valueRange);
      size_t const myHighM = std::min((myid + 1) * myRangeM, m_valueRange);

      value_type *&myLocalCounts = *(localCounts + myid);

      myLocalCounts = new value_type[m_valueRange];

      // Schritt 1.
      memset(myLocalCounts, 0, sizeof(value_type) * m_valueRange);

      // Schritt 2.
      for(size_t i = myLowN; i < myHighN; ++i) {
	value_type v = m_valueGetter(*(beg + i));
	++myLocalCounts[v];
      }

      // (Schritt 3.): transpose
#ifdef _OPENMP
#pragma omp barrier
#endif

      // Schritt 4.
      for(size_t i = myLowM; i < myHighM; ++i) {
	for(size_t j = 1; j < m_nthreads; ++j) {
	  localCounts[j][i] += localCounts[j - 1][i];
	}
        totals[i] = localCounts[m_nthreads-1][i];
      }

      // Schritt 5. (entfällt)
      
      // (Schritt 6.): transpose
#ifdef _OPENMP
#pragma omp barrier
#ifdef __INTEL_COMPILER
      ;
#endif
#endif

      // Schritt 7. (entfällt)

      // Schritt 8.
      size_t offset = 0;

      // Schritt 9.
      if (myid == m_nthreads - 1) {
        for(size_t i = 0; i < m_valueRange; ++i) {
          myLocalCounts[i] += offset;
          prefixSumBegin[i] = myLocalCounts[i];
          offset += totals[i];
        }
      } else {
        for(size_t i = 0; i < m_valueRange; ++i) {
          myLocalCounts[i] += offset;
          offset += totals[i];
        }
      }

#ifdef DEBUG_COUNTING_SORT
#pragma omp critical
      {
	std::cout << "myid: " << myid << " localCounts: " << myLocalCounts << std::endl;
      }
#endif

      // Schritt 10 + 11: werte umsortieren.
      for(size_t i = myLowN; i < myHighN; ++i) {
	value_type const v = m_valueGetter(*(beg + i));
	tmpData[myLocalCounts[v] - 1] = *(beg + i);
        --myLocalCounts[v];
      }

#ifdef _OPENMP
#pragma omp barrier
#endif

      // Schritt 11.
      for(size_t i = myLowN; i < myHighN; ++i) {
        *(beg + i) = tmpData[i];
      }

      delete[] myLocalCounts;

    } // end parallel
    
    delete[] localCounts;
    delete[] totals;
    delete[] tmpData;
  }
  
};
