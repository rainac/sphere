#ifndef jw_sph_spindex_hh
#define jw_sph_spindex_hh
/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

#ifdef SPH_COMPILER_GCC
#include <parallel/algorithm>
#else
#include <algorithm>
#endif

#include "util.hh"
#include "prefix-sums.hh"
#include "sorting/prefixsum.hh"

  /// A BinTriplet object is the answer to a query for particles in a
  /// certain raster cell.
  struct BinTriplet {
    /// The offset into the sorted particle list, where the run of
    /// particles in the desired cell begins.
    unsigned offset;
    /// The offset into the sorted particle list which is just one
    /// past the last moving particle in that cell. This value is
    /// always greater or equal to offset.
    unsigned endMoving;
    /// The offset into the sorted particle list which is just one
    /// past the last non-moving particle in that cell. This value is
    /// always greater or equal to endMoving.
    unsigned endBoundary;
  };


/// Template class SortedPartikelIndex encapsulates the indirect
/// sorting procedure used to map raster cells to particles.

template<class SPHPartikelArray>
struct SortedPartikelIndexInterface {
  typedef typename SPHPartikelArray::value_type particle_type;

  virtual ~SortedPartikelIndexInterface() {};

  virtual void update(SPHPartikelArray const &partikelListe) = 0;
  virtual bool getBinInfo(unsigned indexHash, BinTriplet &triplet) const = 0;
  virtual unsigned getPartikelIndex(unsigned const i) const = 0;
  virtual unsigned numParticlesOut() const = 0;

  virtual void writeXMLTimings(std::ostream &aus) const = 0;
};

  /// class SortIndexMaker builds a single integer value out of the three
  /// particle properties isOut() | hashVal() | isBoundary()
  /// such that sorting for this integer will result in the particles being sorted
  /// lexicographically by these values.
template<class particle_type>
struct SortIndexMaker {
  typedef unsigned value_type;
  unsigned const numBitsHashVal;
  SortIndexMaker(unsigned const numBitsHashVal) : 
    numBitsHashVal(numBitsHashVal) 
  {
    assert(numBitsHashVal <= sizeof(value_type) * 8 - 2);
  }
  value_type operator()(particle_type const &p) const {
    value_type result = p.isOut();
    result <<= numBitsHashVal;
    result |= p.hashVal();
    result <<= 1;
    result |= p.isBoundary();
    return result;
  }
};

template<class SPHPartikelArray>
struct SPICountingSort : public SortedPartikelIndexInterface<SPHPartikelArray> {

  typedef typename SortedPartikelIndexInterface<SPHPartikelArray>::particle_type 
  particle_type;
 
private:

  template<class ArrayIt>
  struct PartikelSortIndexGetter {
    typedef typename SortIndexMaker<particle_type>::value_type value_type;
    typedef unsigned argument_type;

    ArrayIt arrayOffset;
    SortIndexMaker<particle_type> const &sortIndexMaker;
    
    PartikelSortIndexGetter(ArrayIt arrayOffset, SortIndexMaker<particle_type> const &sortIndexMaker) : 
      arrayOffset(arrayOffset),
      sortIndexMaker(sortIndexMaker)
    {}
    
    value_type operator()(argument_type const a) const {
      particle_type const &p = *(arrayOffset + a);
      return sortIndexMaker(p);
    }
    
  };

  HashedIndex<CoordType, ndim> const &m_index;
  Quader<Point<double, ndim> > const &gebiet;
  double const m_h;

  Permutation<unsigned>              partikelIndizes;

  typedef SimpleArray<unsigned>      ListOfOffsets;
  ListOfOffsets                      m_listOffsets;

  unsigned m_size, m_numParticlesOut;

public:

  SPICountingSort(HashedIndex<CoordType, ndim> const &_index, 
                       Quader<Point<double, ndim> > const &gebiet,
                       double const h) :
    m_index(_index),
    gebiet(gebiet),
    m_h(h),
    partikelIndizes(23),
#if __SUNPRO_CC_COMPAT == 5
    m_listOffsets(ListOfOffsets::value_type(), (m_index.size << 2) + 1),
#else
    m_listOffsets((m_index.size << 2) + 1),
#endif
    m_size(),
    m_numParticlesOut()
  { 
  }

  unsigned size() const { return m_size; }

  bool getBinInfo(unsigned indexHash, BinTriplet &triplet) const {
    assert(indexHash < m_listOffsets.size());
    indexHash <<= 1;
    triplet.offset = m_listOffsets[indexHash];
    triplet.endMoving = m_listOffsets[indexHash + 1];
    triplet.endBoundary = m_listOffsets[indexHash + 2];
    return triplet.endBoundary - triplet.offset;
  }

  unsigned getPartikelIndex(unsigned const i) const { 
    return partikelIndizes[i]; 
  }

  void writeXMLTimings(std::ostream &) const { }

  void update(SPHPartikelArray const &partikelListe) {

    m_size = unsigned(partikelListe.size());

    resizeIfNeeded(partikelIndizes, m_size);

    size_t const numBitsSortIndex = m_index.hashFunktion.m_shift * ndim + 2;

    SortIndexMaker<particle_type> sortIndexMaker(m_index.hashFunktion.m_shift * ndim);

    PartikelSortIndexGetter<particle_type const *> 
      partikelSortIndexGetter(address(partikelListe, 0), sortIndexMaker);

    {
#ifdef _OPENMP
      size_t curNumThreads = omp_get_max_threads();
#else
      size_t curNumThreads = 1;
#endif

      SPHParCountingSorter<unsigned*, PartikelSortIndexGetter<particle_type const *> > 
        parCountingSorter(1 << numBitsSortIndex, curNumThreads, partikelSortIndexGetter);
      
      parCountingSorter.sort(&partikelIndizes.m_p[0], &partikelIndizes.m_p[m_size],
                             &m_listOffsets[1]);
      
      m_numParticlesOut = m_size - m_listOffsets[(m_index.size << 1)];
      
    } 

  }

  unsigned numParticlesOut() const {
    return m_numParticlesOut;
  }

};


/// Template class SortedPartikelIndex encapsulates the indirect
/// sorting procedure used to map raster cells to particles.

template<class SPHPartikelArray>
struct SPISTDSortSerialPrefixSums : public SortedPartikelIndexInterface<SPHPartikelArray> {

  typedef typename SortedPartikelIndexInterface<SPHPartikelArray>::particle_type 
  particle_type;

private:
  /// class SortComparator compares two particles lexicographically
  /// after their properties isOut() | hashVal() | isBoundary()
  struct SortComparator {
    SPHPartikelArray const &partikelListe;

    SortComparator(SPHPartikelArray const &_partikelListe) :
      partikelListe(_partikelListe)
    {}

    bool operator()(unsigned const i, unsigned const j) const {
      particle_type const &p = partikelListe[i];
      particle_type const &o = partikelListe[j];
      if (p.isOut() == o.isOut()) {
        if (p.hashVal() == o.hashVal()) {
          return p.isBoundary() < o.isBoundary();
        } else {
          return p.hashVal() < o.hashVal();
        }
      } else {
        return p.isOut() < o.isOut();
      }
    }
  };

  /// class SortComparator compares two particles lexicographically
  /// after their properties isOut() | hashVal() | isBoundary()
  struct IndexedParticleList {
    SPHPartikelArray       const &partikelListe;
    Permutation<unsigned>  const &partikelIndizes;

    IndexedParticleList(SPHPartikelArray       const &_partikelListe,
                        Permutation<unsigned>  const &_partikelIndizes) :
      partikelListe(_partikelListe),
      partikelIndizes(_partikelIndizes)
    {}

    particle_type const &operator[](unsigned const i) const {
      return partikelListe[partikelIndizes[i]];
    }

    unsigned size() const { return partikelListe.size(); }
  };

  HashedIndex<CoordType, ndim> const &m_index;
  Quader<Point<double, ndim> > const &gebiet;
  double const m_h;

  Permutation<unsigned>              partikelIndizes;

  size_t const                       m_valueRange;

  typedef SimpleArray<unsigned>      ListOfOffsets;
  ListOfOffsets                      m_listOffsets;

  unsigned m_size, m_numParticlesOut;

  Timer timeSort, timePrefixSums;

public:

  SPISTDSortSerialPrefixSums(HashedIndex<CoordType, ndim> const &_index, 
                             Quader<Point<double, ndim> > const &gebiet,
                             double const h) :
    m_index(_index),
    gebiet(gebiet),
    m_h(h),
    partikelIndizes(23),
    m_valueRange((m_index.size << 2) + 1),
#if __SUNPRO_CC_COMPAT == 5
    m_listOffsets(ListOfOffsets::value_type(), m_valueRange),
#else
    m_listOffsets(m_valueRange),
#endif
    m_size(),
    m_numParticlesOut()
  { 
  }

  particle_type const &operator()(size_t i) const { return m_size; }

  unsigned size() const { return m_size; }

  bool getBinInfo(unsigned indexHash, BinTriplet &triplet) const {
    assert(indexHash < m_listOffsets.size());
    indexHash <<= 1;
    triplet.offset = m_listOffsets[indexHash];
    triplet.endMoving = m_listOffsets[indexHash + 1];
    triplet.endBoundary = m_listOffsets[indexHash + 2];
    return triplet.endBoundary - triplet.offset;
  }

  unsigned getPartikelIndex(unsigned const i) const { 
    return partikelIndizes[i]; 
  }

  void writeXMLTimings(std::ostream &aus) const {
    aus << "<sorting algorithm='std::sort' time='" << timeSort.lastDiff() << "'/>\n";
    aus << "<prefix-sums algorithm='serial' time='" << timePrefixSums.lastDiff() << "'/>\n";
  }

  void updateOffsetList(SPHPartikelArray const &partikelListe) {

    IndexedParticleList indexedParticleList(partikelListe, partikelIndizes);
    SortIndexMaker<particle_type> sortIndexMaker(m_index.hashFunktion.m_shift * ndim);

    PrefixSumsOfHistogram<IndexedParticleList, SimpleArray<unsigned>, SortIndexMaker<particle_type> >
      prefixSumsOfHistogram(m_valueRange, sortIndexMaker);

    prefixSumsOfHistogram.update(indexedParticleList, m_listOffsets);
    m_numParticlesOut = m_size - m_listOffsets[(m_index.size << 1)];
  }

  void update(SPHPartikelArray const &partikelListe) {

    m_size = unsigned(partikelListe.size());
    resizeIfNeeded(partikelIndizes, m_size);

    SortComparator sortComparator(partikelListe);
    timeSort.start();
#ifdef SPH_COMPILER_GCC
    __gnu_parallel::
#else
      std::
#endif
      sort(&partikelIndizes[0], &partikelIndizes[m_size], sortComparator);
    timeSort.stop();

    timePrefixSums.start();
    updateOffsetList(partikelListe);
    timePrefixSums.stop();
  }

  unsigned numParticlesOut() const {
    return m_numParticlesOut;
  }
};


/// Template class SortedPartikelIndex encapsulates the indirect
/// sorting procedure used to map raster cells to particles.

template<class SPHPartikelArray>
struct SPISTDSortParPrefixSums : public SortedPartikelIndexInterface<SPHPartikelArray> {

  typedef typename SortedPartikelIndexInterface<SPHPartikelArray>::particle_type 
  particle_type;

private:
  /// class SortComparator compares two particles lexicographically
  /// after their properties isOut() | hashVal() | isBoundary()
  struct SortComparator {
    SPHPartikelArray const &partikelListe;

    SortComparator(SPHPartikelArray const &_partikelListe) :
      partikelListe(_partikelListe)
    {}

    bool operator()(unsigned const i, unsigned const j) const {
      particle_type const &p = partikelListe[i];
      particle_type const &o = partikelListe[j];
      if (p.isOut() == o.isOut()) {
        if (p.hashVal() == o.hashVal()) {
          return p.isBoundary() < o.isBoundary();
        } else {
          return p.hashVal() < o.hashVal();
        }
      } else {
        return p.isOut() < o.isOut();
      }
    }
  };

  /// class SortComparator compares two particles lexicographically
  /// after their properties isOut() | hashVal() | isBoundary()
  struct IndexedParticleList {
    SPHPartikelArray       const &partikelListe;
    Permutation<unsigned>  const &partikelIndizes;

    IndexedParticleList(SPHPartikelArray       const &_partikelListe,
                        Permutation<unsigned>  const &_partikelIndizes) :
      partikelListe(_partikelListe),
      partikelIndizes(_partikelIndizes)
    {}

    particle_type const &operator[](unsigned const i) const {
      return partikelListe[partikelIndizes[i]];
    }

    unsigned size() const { return partikelListe.size(); }
  };

  HashedIndex<CoordType, ndim> const &m_index;
  Quader<Point<double, ndim> > const &gebiet;
  double const m_h;
  bool const m_timings;

  Permutation<unsigned>              partikelIndizes;

  size_t const                       m_valueRange;

  typedef SimpleArray<unsigned>      ListOfOffsets;
  ListOfOffsets                      m_listOffsets;

  unsigned m_size, m_numParticlesOut;

  Timer timeSort, timePrefixSums;

public:

  SPISTDSortParPrefixSums(HashedIndex<CoordType, ndim> const &_index, 
                       Quader<Point<double, ndim> > const &gebiet,
                       double const h, bool const _timings = false) :
    m_index(_index),
    gebiet(gebiet),
    m_h(h),
    m_timings(_timings),
    partikelIndizes(23),
    m_valueRange((m_index.size << 2) + 1),
#if __SUNPRO_CC_COMPAT == 5
    m_listOffsets(ListOfOffsets::value_type(), m_valueRange),
#else
    m_listOffsets(m_valueRange),
#endif
    m_size(),
    m_numParticlesOut()
  { 
  }

  particle_type const &operator()(size_t i) const { return m_size; }

  unsigned size() const { return m_size; }

  bool getBinInfo(unsigned indexHash, BinTriplet &triplet) const {
    assert(indexHash < m_listOffsets.size());
    indexHash <<= 1;
    triplet.offset = m_listOffsets[indexHash];
    triplet.endMoving = m_listOffsets[indexHash + 1];
    triplet.endBoundary = m_listOffsets[indexHash + 2];
    return triplet.endBoundary - triplet.offset;
  }

  unsigned getPartikelIndex(unsigned const i) const { 
    return partikelIndizes[i]; 
  }

  void writeXMLTimings(std::ostream &aus) const {
    aus << "<sorting algorithm='std::sort' time='" << timeSort.lastDiff() << "'/>\n";
    aus << "<prefix-sums algorithm='serial' time='" << timePrefixSums.lastDiff() << "'/>\n";
  }

  void updateOffsetList(SPHPartikelArray const &partikelListe) {

    IndexedParticleList indexedParticleList(partikelListe, partikelIndizes);
    SortIndexMaker<particle_type> sortIndexMaker(m_index.hashFunktion.m_shift * ndim);

    ParPrefixSums<IndexedParticleList, SortIndexMaker<particle_type>, SimpleArray, unsigned > 
      parPrefixSums(m_valueRange, sortIndexMaker);

    parPrefixSums.sample(indexedParticleList, m_listOffsets);
    m_numParticlesOut = m_size - m_listOffsets[(m_index.size << 1)];
  }

  void update(SPHPartikelArray const &partikelListe) {

    m_size = unsigned(partikelListe.size());
    resizeIfNeeded(partikelIndizes, m_size);

    SortComparator sortComparator(partikelListe);
    timeSort.start();
#ifdef SPH_COMPILER_GCC
    __gnu_parallel::
#else
      std::
#endif
      sort(&partikelIndizes[0], &partikelIndizes[m_size], sortComparator);
    timeSort.stop();

    timePrefixSums.start();
    updateOffsetList(partikelListe);
    timePrefixSums.stop();
  }

  unsigned numParticlesOut() const {
    return m_numParticlesOut;
  }
};

template<class SPHPartikelArray>
struct SortedPartikelIndex : public SortedPartikelIndexInterface<SPHPartikelArray> {
  typedef Partikel<SphADataType> particle_type;

#include "sorted-index-types.ncd.enum.hh"
#include "sorted-index-types.ncd.cc"

  SortedPartikelIndexInterface<SPHPartikelArray> *m_theIndex;

  SortedPartikelIndex(std::string const &name,
                      HashedIndex<CoordType, ndim> const &_index, 
                      Quader<Point<double, ndim> > const &gebiet,
                      double const h) : 
    m_theIndex() {
    switch(getSortedIndexTypesValue(name)) {
    case SORTED_INDEX_CSORT:
      m_theIndex = new SPICountingSort<SPHPartikelArray>(_index, gebiet, h);
      break;
    case SORTED_INDEX_STDSORT_SER_INDEX:
      m_theIndex = new SPISTDSortSerialPrefixSums<SPHPartikelArray>(_index, gebiet, h);
      break;
    case SORTED_INDEX_STDSORT_PAR_INDEX:
      m_theIndex = new SPISTDSortParPrefixSums<SPHPartikelArray>(_index, gebiet, h);
      break;
    default:
      cerr << "no such sorted index named `" << name << "'\n";
      break;
    }
  }

  ~SortedPartikelIndex() {
    if (m_theIndex) {
      delete m_theIndex;
      m_theIndex = 0;
    }
  }

  void update(SPHPartikelArray const &partikelListe) {
    m_theIndex->update(partikelListe);
  }
  bool getBinInfo(unsigned indexHash, BinTriplet &triplet) const {
    return m_theIndex->getBinInfo(indexHash, triplet);
  }
  unsigned getPartikelIndex(unsigned const i) const {
    return m_theIndex->getPartikelIndex(i);
  }
  unsigned numParticlesOut() const {
    return m_theIndex->numParticlesOut();
  }
  void writeXMLTimings(std::ostream &aus) const {
    m_theIndex->writeXMLTimings(aus);
  }

};


#endif
