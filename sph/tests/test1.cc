#include "CoreSuite.h"
#include "TestTest.h"

using CPPUNIT_NS::TestFixture;

class Test1 : public TestFixture {

  CPPUNIT_TEST_SUITE( Test1 );
  CPPUNIT_TEST( testSomething );
  CPPUNIT_TEST_EXCEPTION( testThrow, std::invalid_argument );
  CPPUNIT_TEST_SUITE_END();

public:
  Test1()  : TestFixture() {
  }
 
  ~Test1() {
  }
 
  void setUp() {
  }

  void tearDown() {
  }

  void testSomething() {
    CPPUNIT_ASSERT_MESSAGE("das sollte stimmen: ", true);
  }

  void testThrow() {
    throw std::invalid_argument("just like this");
  }

};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( Test1,
                                       coreSuiteName() );

