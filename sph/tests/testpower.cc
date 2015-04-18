#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestSuite.h>
#include <cppunit/TestPath.h>

#include "random-array.hh"
#include "CoreSuite.h"
#include "../src/power.hh"

using CPPUNIT_NS::TestFixture;

class TestPower : public TestFixture {

  CPPUNIT_TEST_SUITE( TestPower );

  CPPUNIT_TEST( testDynamicPower );
  CPPUNIT_TEST( testStaticPower );

  CPPUNIT_TEST_SUITE_END();

  RandomArray *randArray;

public:
  TestPower()  : 
    TestFixture(), 
    randArray() 
  { }
 
  ~TestPower() { }
 
  void setUp() {
    randArray = new RandomArray(100, 10000);
  }

  void tearDown() {
    delete randArray;
    randArray = 0;
  }

  void testStaticPower() {
    CPPUNIT_ASSERT_MESSAGE("das sollte stimmen: ", true);
    
    for(size_t i = 0; i < randArray->size(); ++i) {
      double v = (*randArray)[i];
      CPPUNIT_ASSERT_DOUBLES_EQUAL(pow(v, 3), static_power<3>(v), pow(v, 3) * 1e-14);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(pow(v, 7), static_power<7>(v), pow(v, 7) * 1e-14);
    }
  }

  void testDynamicPower() {
    CPPUNIT_ASSERT_MESSAGE("das sollte stimmen: ", true);
    
    for(size_t i = 0; i < randArray->size(); ++i) {
      double v = (*randArray)[i];
      CPPUNIT_ASSERT_DOUBLES_EQUAL(pow(v, 3), dynamic_power(v, 3), pow(v, 3)*1e-14);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(pow(v, 7), dynamic_power(v, 7), pow(v, 7)*1e-14);
    }
  }

};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestPower,
                                       coreSuiteName() );

