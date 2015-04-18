#include "CoreSuite.h"
#include "TestTest.h"
#include "sortCommon.hh"
#include "../src/counting-sort.hh"

using namespace std;
using CPPUNIT_NS::TestFixture;

class TestCountingSort : public TestFixture {

  CPPUNIT_TEST_SUITE( TestCountingSort );
  CPPUNIT_TEST( testSomething );
  CPPUNIT_TEST_EXCEPTION( testThrow, std::invalid_argument );

  CPPUNIT_TEST( resultIsSorted1 );
  CPPUNIT_TEST( resultIsSorted2 );

//   CPPUNIT_TEST( sortIsStable1 ); \todo algorithm is not stable
  CPPUNIT_TEST( histogram );

  CPPUNIT_TEST_SUITE_END();

  std::valarray<long> randData, randDataWithZeros;
  std::valarray<Data> randDataWithIds;
  std::valarray<unsigned> histoOfData;

  unsigned short randSeed[3];

  static size_t const numData = 1000;
  static size_t const dataRange = 100000;
#ifdef _OPENMP
  static size_t const numThreads = 4;
#else
  static size_t const numThreads = 1;
#endif

public:
  TestCountingSort()  : TestFixture() {
  }
 
  ~TestCountingSort() {
  }
 
  void setUp() {
    randSeed[0] = 123;
    randSeed[1] = 132;
    randSeed[2] = 213;

    randData.resize(numData);
    for(size_t i = 0; i < numData; ++i) {
      randData[i] = nrand48(randSeed) % dataRange;
    }

    randDataWithZeros.resize(numData);
    randDataWithZeros = randData;
    for(size_t i = 0; i < numData/20; ++i) {
      randDataWithZeros[lrand48() % numData] = 0;
      randDataWithZeros[lrand48() % numData] = dataRange/2;
      randDataWithZeros[lrand48() % numData] = dataRange-1;
    }

    histoOfData.resize(dataRange + 1);
    for(size_t i = 0; i < numData; ++i) {
      long const d = randDataWithZeros[i];
      ++histoOfData[d + 1];
    }
    for(size_t i = 1; i < dataRange + 1; ++i) {
      histoOfData[i] += histoOfData[i - 1];
    }

    randDataWithIds.resize(numData);
    for(size_t i = 0; i < numData; ++i) {
      randDataWithIds[i].first = randDataWithZeros[i];
      randDataWithIds[i].second = i;
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
    std::valarray<unsigned> histoPrefix(dataRange);

#ifdef _OPENMP
    omp_set_num_threads(numThreads);
#endif
    SPHParCountingSorter<long *> countingSorter(dataRange, numThreads);

    countingSorter.sort<unsigned*>(&data[0], &data[data.size()], &histoPrefix[0]);

    std::sort(&data2[0], &data2[data.size()]);
    
    CPPUNIT_ASSERT((data == data2).min() == 1);
    CPPUNIT_ASSERT(isSorted(data));
    CPPUNIT_ASSERT(isSorted(data2));
  }

  void resultIsSorted2() {
    std::valarray<long> data(randDataWithZeros);
    std::valarray<long> data2(data);
    std::valarray<unsigned> histoPrefix(dataRange);

#ifdef _OPENMP
    omp_set_num_threads(numThreads);
#endif
    SPHParCountingSorter<long *> countingSorter(dataRange, numThreads);

    countingSorter.sort<unsigned*>(&data[0], &data[data.size()], &histoPrefix[0]);

    std::sort(&data2[0], &data2[data.size()]);
    
    CPPUNIT_ASSERT((data == data2).min() == 1);
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
    std::valarray<Data> data(randDataWithIds);
    std::valarray<Data> data2(data);

    std::valarray<unsigned> histoPrefix(dataRange);

#ifdef _OPENMP
    omp_set_num_threads(numThreads);
#endif
    SPHParCountingSorter<Data *, FirstGetter> countingSorter(dataRange, numThreads);

    countingSorter.sort<unsigned *>(&data[0], &data[data.size()], &histoPrefix[0]);

    std::stable_sort(&data2[0], &data2[data.size()]);
    
    CPPUNIT_ASSERT(isStable(data2));
    CPPUNIT_ASSERT(isStable(data)); // \todo algorithm is not stable
    CPPUNIT_ASSERT((data == data2).min() == 1);
  }

  void histogram() {
    std::valarray<long> data(randDataWithZeros);
    std::valarray<unsigned> histoPrefix(dataRange + 1);

#ifdef _OPENMP
    omp_set_num_threads(numThreads);
#endif
    SPHParCountingSorter<long *> countingSorter(dataRange, numThreads);
    countingSorter.sort<unsigned *>(&data[0], &data[data.size()], &histoPrefix[1]);

    CPPUNIT_ASSERT(isSorted(data));

    for(size_t i = 0; i < dataRange; ++i) {
//       cout << "if there are " << i
//            << "s, they begin at: " << histoPrefix[i] << "\n";
      if (histoPrefix[i] >= data.size()) {
        cout << "no more data\n";
        break;
      }
      CPPUNIT_ASSERT(data[histoPrefix[i]] >= long(i));
      size_t const numEntries = histoPrefix[i + 1] - histoPrefix[i];
      if (numEntries) {
//         cout << "number of " << i << "s: " << numEntries << "\n";
        for(size_t k = 0; k < numEntries; ++k) {
          CPPUNIT_ASSERT(data[histoPrefix[i]  + k] == long(i));
        }
      }
    }

//     cout << "prefix histogram1: " << histoOfData << "\n";
//     cout << "prefix histogram2: " << histoPrefix << "\n";
//     cout << "prefix histograms coincidence " << (histoPrefix == histoOfData) << "\n";
    CPPUNIT_ASSERT((histoPrefix == histoOfData).min() == 1);

  }

};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestCountingSort,
                                       coreSuiteName() );

