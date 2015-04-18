#ifndef jw_sph_prefix_sums_hh
#define jw_sph_prefix_sums_hh
/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/
#include <assert.h>

template<class T>
struct DefaultValueGetter {
  typedef T value_type;
  T const &operator()(T const &t) { return t; }
};

template<class SortedSPHPartikelArray, 
         class PrefixSumArray=std::valarray<unsigned>, 
         class ValueGetter=DefaultValueGetter<typename SortedSPHPartikelArray::value_type> >
struct PrefixSumsOfHistogram {

  size_t                      const m_valueRange;
  ValueGetter                       m_valueGetter;

  PrefixSumsOfHistogram(size_t const valueRange) : 
    m_valueRange(valueRange) 
  {}

  PrefixSumsOfHistogram(size_t const valueRange, ValueGetter const &_valueGetter) : 
    m_valueRange(valueRange),
    m_valueGetter(_valueGetter)
  {}

  void update(SortedSPHPartikelArray const &partikelListe, PrefixSumArray &m_listOffsets) {
    assert( m_listOffsets.size() == m_valueRange );

    unsigned lastIndex = m_valueGetter(partikelListe[0]);
    // unsigned lasti = 0;
    // unsigned lasti_bound = 0;
    m_listOffsets[lastIndex] = 0;
    for(unsigned i = 0; i < partikelListe.size(); ++i) {
      unsigned const curIndex = m_valueGetter(partikelListe[i]);
      if (curIndex != lastIndex) {
// #ifdef _OPENMP
// #pragma omp parallel for schedule(static)
// #endif
        for(unsigned j = lastIndex + 1; j <= curIndex; ++j) {
          m_listOffsets[j] = i;
        }
        lastIndex = curIndex;
        // lasti_bound = 0;
        // lasti = i;
      }
    }
    for(unsigned j = lastIndex + 1; j < m_listOffsets.size(); ++j) {
      m_listOffsets[j] = partikelListe.size();
    }
  }
};

#endif
