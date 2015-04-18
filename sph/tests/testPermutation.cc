#include <omp.h>
#include "CoreSuite.h"
#include "TestTest.h"
#include "../src/permutation.hh"

using CPPUNIT_NS::TestFixture;

class TestPermutation : public TestFixture {

  CPPUNIT_TEST_SUITE( TestPermutation );

  CPPUNIT_TEST( isInitialized );

  CPPUNIT_TEST_SUITE_END();

  std::valarray<long> randData;

  static size_t const numData = 10000;
  static size_t const dataRange = 100000;
#ifdef _OPENMP
  static size_t const numThreads = 4;
#else
  static size_t const numThreads = 1;
#endif

public:
  TestPermutation()  : TestFixture() { }
 
  ~TestPermutation() { }
 
  void setUp() {
#ifdef _OPENMP
    omp_set_num_threads(numThreads);
#endif

    randData.resize(numData);
    for(size_t i = 0; i < numData; ++i) {
      randData[i] = lrand48() % dataRange;
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

  void isInitialized() {
    CPPUNIT_ASSERT_MESSAGE("das sollte stimmen: ", true);

    Permutation<long> p(numData);

    CPPUNIT_ASSERT_MESSAGE("p.size() is wrong: ", 
                           p.size() == numData
                           );

    CPPUNIT_ASSERT_MESSAGE("p.m_p.size() is wrong: ", 
                           p.m_p.size() == numData
                           );

    CPPUNIT_ASSERT_MESSAGE("p.m_n is wrong: ", 
                           p.m_n == numData
                           );

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(long i = 0; i < long(numData); ++i) {
      CPPUNIT_ASSERT_MESSAGE("i-th entry in identity permuation should be i: ", 
                             p[i] == i
                             );
      
    }

  }


};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestPermutation,
                                       coreSuiteName() );

