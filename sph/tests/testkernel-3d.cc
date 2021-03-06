#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestSuite.h>
#include <cppunit/TestPath.h>
#include <assert.h>

#include "random-array.hh"
#include "CoreSuite.h"

#define SPH_DIM 3

#include "../src/kernels.hh"

#ifdef SPH_AD_YAFAD
#define YAFAD_POINT_DEFINED
#include "yafad/yafad.hh"
typedef yafad::FO::Static::Adouble Adouble;
#endif

using CPPUNIT_NS::TestFixture;

class TestKernel3D : public TestFixture {

  CPPUNIT_TEST_SUITE( TestKernel3D );

  CPPUNIT_TEST( nonZeroInsideReach );
  CPPUNIT_TEST( zeroOutsideReach );

  CPPUNIT_TEST( combinedFunction );
#ifdef SPH_AD_YAFAD
  CPPUNIT_TEST( checkGradientAD );
#endif
  CPPUNIT_TEST( checkIntegralIsOne );

  CPPUNIT_TEST_SUITE_END();

  typedef Point<double, SPH_DIM> Vector;
  typedef KernelInterface<double, Vector, SPH_DIM> Kernel;

  std::vector<std::string> kernelNames;
  std::vector<Kernel *> kernels;

#ifdef SPH_AD_YAFAD
  typedef Point<Adouble, SPH_DIM> ADVector;
  typedef KernelInterface<Adouble, ADVector, SPH_DIM> ADKernel;
  std::vector<ADKernel *> ADkernels;
#endif

  DRandomArray *randArray;

  double const H;
  size_t const numTestPointsGradient;

public:
  TestKernel3D()  : 
    TestFixture(), 
    randArray(),
    H(0.01),
    numTestPointsGradient(200)
  { }
 
  ~TestKernel3D() { }
 
  void setUp() {
    randArray = new DRandomArray(100000);

    kernelNames.push_back("gauss");
    kernelNames.push_back("lucy");
    kernelNames.push_back("spline");
    kernelNames.push_back("morris4");
    kernelNames.push_back("morris5");
    kernelNames.push_back("johnson");
    kernelNames.push_back("wendland0");
    kernelNames.push_back("wendland1");
    kernelNames.push_back("wendland2");
    kernelNames.push_back("wendland3");

    kernels.resize(kernelNames.size());
#ifdef SPH_AD_YAFAD
    ADkernels.resize(kernelNames.size());
#endif
    for(size_t i = 0; i < kernelNames.size(); ++i) {
      kernels[i] = Kernel::makeKernel(kernelNames[i], H);
#ifdef SPH_AD_YAFAD
      ADkernels[i] = ADKernel::makeKernel(kernelNames[i], H);
#endif
    }
  }

  void tearDown() {
    delete randArray;

    for(size_t i = 0; i < kernels.size(); ++i) {
      delete kernels[i];
    }
  }



  void nonZeroInsideReach() {
    for(size_t i = 0; i < kernelNames.size(); ++i) {
      double const r = sqrt(3.) * (H - (kernelNames[i].compare("johnson") == 0 ? 1e-9 : 1e-12));
      Vector dist(SPH_DIM), grad;
      dist[0] = r/SPH_DIM;
      dist[1] = r/SPH_DIM;
      dist[2] = r/SPH_DIM;
      double const dnorm = norm(dist);
//       std::cout << "non-zero test: r = " << dist << ", |r| = " << dnorm << ", H - |r| = " << H - dnorm << "\n";
      double const w = kernels[i]->w(dnorm);
      kernels[i]->gradw(dist, dnorm, grad);
//       std::cout << kernelNames[i] << " W(r) = " << w << " dW(r) = " << grad << "\n";
      CPPUNIT_ASSERT_MESSAGE(kernelNames[i], w > 0);
      CPPUNIT_ASSERT_MESSAGE(kernelNames[i], grad.min() != 0);
      CPPUNIT_ASSERT_MESSAGE(kernelNames[i], grad.max() != 0);
    }
  }

  void zeroOutsideReach() {
    double const r = sqrt(3.) * (H + 1e-12);
    Vector dist(SPH_DIM), grad;
    dist[0] = r/SPH_DIM;
    dist[1] = r/SPH_DIM;
    dist[2] = r/SPH_DIM;
    double const dnorm = norm(dist);
    for(size_t i = 0; i < kernelNames.size(); ++i) {
      CPPUNIT_ASSERT_MESSAGE(kernelNames[i], dnorm >= H);
//       std::cout << dnorm << " " << kernels[i]->w(dnorm) << "\n";
      CPPUNIT_ASSERT_MESSAGE(kernelNames[i], kernels[i]->w(dnorm) == 0);
      kernels[i]->gradw(dist, dnorm, grad);
//       std::cout << dist << " " << dnorm << " " << grad << "\n";
      CPPUNIT_ASSERT_MESSAGE(kernelNames[i], grad.min() == 0);
      CPPUNIT_ASSERT_MESSAGE(kernelNames[i], grad.max() == 0);
    }
  }

  void combinedFunction() {
    double const r = H;
    Vector dist(SPH_DIM), grad1, grad2;
    double w1, w2;
    dist[0] = r/SPH_DIM;
    dist[1] = r/SPH_DIM;
    dist[2] = r/SPH_DIM;
    double const dnorm = norm(dist);
    for(size_t i = 0; i < kernelNames.size(); ++i) {
      w1 = kernels[i]->w(dnorm);
      kernels[i]->gradw(dist, dnorm, grad1);

      w2 = kernels[i]->wAndGrad(dist, dnorm, grad2);

      CPPUNIT_ASSERT_MESSAGE(kernelNames[i], w1 >= 0);
      CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(kernelNames[i], w1, w2, abs(w1)*1e-15);

      for(int k=0; k< SPH_DIM; ++k) {
        CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(kernelNames[i], grad1[k], grad2[k], abs(grad1[k])*1e-15);
      }
    }
  }

#ifdef SPH_AD_YAFAD
  void checkGradientAD() {
    ADVector dist1;
    Vector dist2, grad1, grad2;
    Adouble w1;
    for(size_t iTestPoint = 0; iTestPoint < numTestPointsGradient; ++iTestPoint) {
      for(size_t i = 0; i < SPH_DIM; ++i) {
        dist1[i] = dist2[i] = H/sqrt(SPH_DIM) * (*randArray)[iTestPoint * SPH_DIM + i];
        dist1[i].diff(i) = 1;
      }
      Adouble const dnorm1 = norm(dist1);
      double const dnorm2 = norm(dist2);
//       std::cout << "gradient-test: point x = " << dist1 << ", |x| = " << dnorm1 << "\n";

      assert(dnorm1 <= H);

      for(size_t i = 0; i < kernelNames.size(); ++i) {
        w1 = ADkernels[i]->w(dnorm1);
        
        for(size_t k = 0; k < SPH_DIM; ++k) {
          grad1[k] = w1.diff(k);
        }
        
        kernels[i]->gradw(dist2, dnorm2, grad2);
        
//         std::cout << "gradients: " << kernelNames[i] << " \tgradAD = " << grad1 << " \tgradNum = " << grad2 << "\n";
        
        for(size_t k = 0; k < SPH_DIM; ++k) {
//           std::cout << k << " " << grad1[k] - grad2[k] << "\n";
          CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(kernelNames[i], 
                                               grad1[k], grad2[k], fabs(grad1[k])*1e-13);
        }
      }
    }
  }
#endif

  void checkIntegralIsOne() {
    size_t const N = 100000;
    double const dx = 1.0 / N;

    for(size_t i = 0; i < kernelNames.size(); ++i) {
      double sum = 0;
      for(size_t k = 0; k < N-1; ++k) {
        double const v1 = k * dx;
        double const v2 = (k+1) * dx;
        
        double const f1 = 4 * M_PI * kernels[i]->w(v1) * v1*v1;
        double const f2 = 4 * M_PI * kernels[i]->w(v2) * v2*v2;

        sum += dx * (f1 + f2) / 2;
      }

//       std::cout << kernelNames[i] << " sum = " << sum << "\n";

      if (kernelNames[i] == "gauss") {
        CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(kernelNames[i], 1.0, sum, 1e-3);
      } else if (kernelNames[i] == "wendland0") {
        CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(kernelNames[i], 1.0, sum, 1e-7);
      } else {
        CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(kernelNames[i], 1.0, sum, 1e-12);
      }
    }
  }

};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestKernel3D,
                                       coreSuiteName() );

