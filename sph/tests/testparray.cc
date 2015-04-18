#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestSuite.h>
#include <cppunit/TestPath.h>
#include <assert.h>

#include "random-array.hh"
#include "CoreSuite.h"
#include "../src/parray.hh"

using CPPUNIT_NS::TestFixture;

class TestPArray : public TestFixture {

  CPPUNIT_TEST_SUITE( TestPArray );

  CPPUNIT_TEST( zeroInit );
  CPPUNIT_TEST( assignmentByConstant );
  CPPUNIT_TEST( assignmentByPArray );
//   CPPUNIT_TEST( assignmentByValarray );

  CPPUNIT_TEST( assignmentByLongerPArray );

  CPPUNIT_TEST_SUITE_END();

  RandomArray *randArray;

public:
  TestPArray()  : 
    TestFixture(), 
    randArray() 
  { }
 
  ~TestPArray() { }
 
  void setUp() {
    randArray = new RandomArray(100000, 10000);
  }

  void tearDown() {
    delete randArray;
    randArray = 0;
  }

  void zeroInit() {
    PArray<long> p(randArray->size());

    CPPUNIT_ASSERT(p.m_data.size() == randArray->size());
    CPPUNIT_ASSERT(p.m_data.sum() == 0);
    
  }

  void assignmentByConstant() {
    PArray<long> p(100);
    p = 2;
    CPPUNIT_ASSERT(p.m_data.sum() == 200);
  }

  void assignmentByPArray() {
    PArray<long> p1(100), p2(100);
    p1 = 1;
    p2 = p1;
    CPPUNIT_ASSERT(p2.m_data.sum() == 100);
  }

//   void assignmentByValarray() {
//     std::valarray<long> p1(100);
//     PArray<long> p2(100);
//     p1 = 1;
//     p2.m_data = p1;
//     CPPUNIT_ASSERT(p2.m_data.sum() == 100);
//   }

  void assignmentByLongerPArray() {
    PArray<long> p1(200), p2(100);
    p1 = 1;
    p2 = p1;
    CPPUNIT_ASSERT(p2.m_data.sum() == 100);
  }

};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestPArray,
                                       coreSuiteName() );

