#include <omp.h>
#include "CoreSuite.h"
#include "TestTest.h"
#include "sortCommon.hh"
#include "../src/prefix-sums.hh"
#include "sorting/prefixsum.hh"

using namespace std;
using CPPUNIT_NS::TestFixture;

class TestSTDSort : public TestFixture {

  CPPUNIT_TEST_SUITE( TestSTDSort );
  CPPUNIT_TEST( testSomething );
  CPPUNIT_TEST_EXCEPTION( testThrow, std::invalid_argument );

  CPPUNIT_TEST( resultIsSorted1 );
  CPPUNIT_TEST( resultIsSorted2 );

  CPPUNIT_TEST( sortIsStable1 );

  CPPUNIT_TEST( histogram );
  CPPUNIT_TEST( histogram2 );

  CPPUNIT_TEST_SUITE_END();

  std::valarray<long> randData, randDataWithZeros;
  std::valarray<unsigned> histoOfData;

  static size_t const numData = 1000;
  static size_t const dataRange = 10000;
#ifdef _OPENMP
  static size_t const numThreads = 4;
#else
  static size_t const numThreads = 1;
#endif

public:
  TestSTDSort()  : TestFixture() {
  }
 
  ~TestSTDSort() {
  }
 
  void setUp() {
    randData.resize(numData);
    for(size_t i = 0; i < numData; ++i) {
      randData[i] = lrand48() % dataRange;
    }

    randDataWithZeros.resize(numData);
    randDataWithZeros = randData;
    for(size_t i = 0; i < numData/10; ++i) {
      randDataWithZeros[lrand48() % numData] = 0;
    }

    histoOfData.resize(dataRange + 1);
    for(size_t i = 0; i < numData; ++i) {
      long const d = randDataWithZeros[i];
      ++histoOfData[d + 1];
    }
    for(size_t i = 1; i < dataRange + 1; ++i) {
      histoOfData[i] += histoOfData[i - 1];
    }
  }

  void tearDown() {
  }

  void testSomething() {
    CPPUNIT_ASSERT_MESSAGE("das sollte stimmen: ", true);
  }

  void testThrow() {
    throw std::invalid_argument("just like this");
  }

  void resultIsSorted1() {
    std::valarray<long> data(randData);
    std::valarray<long> data2(data);

#ifdef _OPENMP
    omp_set_num_threads(numThreads);
#endif

    std::sort(&data[0], &data[data.size()]);
    std::stable_sort(&data2[0], &data2[data.size()]);
    
//     cout << "sorted result: (std::sort)" << data << "\n";
//     cout << "sorted result: (std::stable_sort)" << data2 << "\n";

    CPPUNIT_ASSERT(isSorted(data));
    CPPUNIT_ASSERT(isSorted(data2));
  }

  void resultIsSorted2() {
    std::valarray<long> data(randDataWithZeros);
    std::valarray<long> data2(data);

#ifdef _OPENMP
    omp_set_num_threads(numThreads);
#endif

    std::sort(&data[0], &data[data.size()]);
    std::stable_sort(&data2[0], &data2[data.size()]);
    
//     cout << "sorted result: (std::sort)" << data << "\n";
//     cout << "sorted result: (std::stable_sort)" << data2 << "\n";

    CPPUNIT_ASSERT(isSorted(data));
    CPPUNIT_ASSERT(isSorted(data2));
  }

  struct FirstGetter {
    typedef unsigned value_type;
    typedef Data argument_type;
    unsigned operator()(Data const &v) const {
      return v.first;
    }
  };

  void sortIsStable1() {
    std::valarray<Data> data(numData);

    for(size_t i = 0; i < numData; ++i) {
      data[i].first = randDataWithZeros[i];
      data[i].second = i;
    }

#ifdef _OPENMP
    omp_set_num_threads(numThreads);
#endif

    std::stable_sort(&data[0], &data[data.size()]);
    
    CPPUNIT_ASSERT(isStable(data));
  }


  void histogram() {
    std::valarray<long> data(randDataWithZeros);
    std::valarray<unsigned> prefixSumsArray(histoOfData.size());

    std::stable_sort(&data[0], &data[data.size()]);

    PrefixSumsOfHistogram<std::valarray<long> > prefixSumsOfHistogram(histoOfData.size());

    prefixSumsOfHistogram.update(data, prefixSumsArray);

//     cout << "prefix histogram1: " << histoOfData << "\n";
//     cout << "prefix histogram2: " << prefixSumsArray << "\n";

    CPPUNIT_ASSERT((prefixSumsArray == histoOfData).min() == 1);

  }

  void histogram2() {
    std::valarray<long> data(randDataWithZeros);
    std::valarray<unsigned> prefixSumsArray(histoOfData.size());

    std::stable_sort(&data[0], &data[data.size()]);

    ParPrefixSums<std::valarray<long>, DefaultValueGetter<long> > 
      parPrefixSums(histoOfData.size());

    parPrefixSums.sample(data, prefixSumsArray);

//     cout << "prefix histogram1: " << histoOfData << "\n";
//     cout << "prefix histogram2: " << prefixSumsArray << "\n";

    CPPUNIT_ASSERT((prefixSumsArray == histoOfData).min() == 1);

  }

};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestSTDSort,
                                       coreSuiteName() );

