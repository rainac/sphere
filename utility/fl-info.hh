/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

#include <math.h>
#include <fenv.h>

struct FLInfo {
  int fltEvalMethod, fpFastFMA, fp_ilogb0, fp_ilogbNaN;

  double hugeVal;
#ifdef HUGE_VALF
  float hugeValF;
#endif
#ifdef HUGE_VALL
  long double hugeValL;
#endif
  float infinity;
  float nan;

  long double ldblEpsilon;
  double       dblEpsilon;
  float        fltEpsilon;

//   double m_e, m_log2e, m_log10e,
//     m_ln2, m_ln10, 
//     m_pi, m_pi_2, m_pi_4,
//     m_1_pi, m_2_pi,
//     m_2_sqrtpi, m_sqrt2, m_sqrt1_2;

  FLInfo() : 
#ifdef __FLT_EVAL_METHOD__
    fltEvalMethod(__FLT_EVAL_METHOD__),
#elif defined FLT_EVAL_METHOD
    fltEvalMethod(FLT_EVAL_METHOD),
#else
#warning no FLT_EVAL_METHOD
#endif
#ifdef __FP_FAST_FMA__
    fpFastFMA(__FP_FAST_FMA__),
#elif defined FP_FAST_FMA
    fpFastFMA(FP_FAST_FMA),
#else
    fpFastFMA(),
// #warning no FP_FAST_FMA
#endif
    hugeVal(HUGE_VAL),
#ifdef HUGE_VALF
    hugeValF(HUGE_VALF),
#endif
#ifdef HUGE_VALL
    hugeValL(HUGE_VALL),
#endif
#ifdef INFINITY
    infinity(INFINITY),
#else
    infinity(1.0 / 0.0),
#endif
#ifdef NAN
    nan(NAN)
#else
    nan(0.0 / 0.0)
#endif
  {
#ifdef __LDBL_EPSILON__
    ldblEpsilon = __LDBL_EPSILON__;
#else
    ldblEpsilon = nextafterl(1, 2) - 1;
#endif
#ifdef __DBL_EPSILON__
    dblEpsilon = __DBL_EPSILON__;
#else
    dblEpsilon = nextafter(1, 2) - 1;
#endif
#ifdef __FLT_EPSILON__
    fltEpsilon = __FLT_EPSILON__;
#else
    fltEpsilon = nextafterf(1, 2) - 1;
#endif
#if ! defined INFINITY && ! defined NAN
    feclearexcept(FE_ALL_EXCEPT);
#endif
  }

  void writeXML(std::ostream &aus, std::string const &indent = "") const {
    size_t const cprec = aus.precision();
    aus.precision(-ilogb(ldblEpsilon));
    aus << indent << "<fl-info\n"
        << indent << " flt-eval-method='" << fltEvalMethod << "'"
        << indent << " fp-fast-fma='" << fpFastFMA << "'"
#ifdef FP_ILOBG0
        << indent << " fp-ilogb-0='" << FP_ILOBG0 << "'"
#endif
#ifdef FP_ILOBGNAN
        << indent << " fp-ilogb-nan='" << FP_ILOBGNAN << "'"
#endif
        << indent << " huge-val='" << HUGE_VAL << "'"
        << indent << " infinity='" << infinity << "'"
        << indent << " nan='" << nan << "'"
        << indent << " float-epsilon='" << fltEpsilon << "'"
        << indent << " double-epsilon='" << dblEpsilon << "'"
        << indent << " long-double-epsilon='" << ldblEpsilon << "'"
        << indent << " m_e='" << M_E << "'"
        << indent << " m_log2e='" << M_LOG2E << "'"
        << indent << " m_log10e='" << M_LOG10E << "'"
        << indent << " m_ln2='" << M_LN2 << "'"
        << indent << " m_ln10='" << M_LN10 << "'"
        << indent << " m_pi='" << M_PI << "'"
        << indent << " m_pi_2='" << M_PI_2 << "'"
        << indent << " m_pi_4='" << M_PI_4 << "'"
        << indent << " m_1_pi='" << M_1_PI << "'"
        << indent << " m_2_pi='" << M_2_PI << "'"
        << indent << " m_2_sqrtpi='" << M_2_SQRTPI << "'"
        << indent << " m_sqrt2='" << M_SQRT2 << "'"
        << indent << " m_sqrt1_2='" << M_SQRT1_2 << "'"
#ifdef MATH_ERRNO
        << indent << " math-errno='" << MATH_ERRNO << "'"
#endif
#ifdef MATH_ERREXCEPT
        << indent << " math-errexcept='" << MATH_ERREXCEPT << "'"
#endif
#ifdef math_errhandling
        << indent << " math-errhandling='" << math_errhandling << "'"
#endif
      "/>";
    aus.precision(cprec);
  }

  typedef int FExceptions;
#include "fe-codes.ncd.cc"

  static void printExceptions(std::ostream &aus, int const fes) {
    aus << "<floting-point-exeptions>";
    if (fes & FE_INVALID) {
      aus << getFExceptionsName(FE_INVALID) << " ";
    }
    if (fes & FE_INEXACT) {
      aus << getFExceptionsName(FE_INEXACT) << " ";
    }
    if (fes & FE_OVERFLOW) {
      aus << getFExceptionsName(FE_OVERFLOW) << " ";
    }
    if (fes & FE_UNDERFLOW) {
      aus << getFExceptionsName(FE_UNDERFLOW) << " ";
    }
    if (fes & FE_DIVBYZERO) {
      aus << getFExceptionsName(FE_DIVBYZERO) << " ";
    }
    aus << "</floting-point-exeptions>";
  }

  static int testAndShowExceptions(std::ostream &aus, int feMask = FE_ALL_EXCEPT) {
    int const fes = fetestexcept(feMask);
    if (fes) {
      printExceptions(aus, fes);
    }
    return fes;
  }

};
