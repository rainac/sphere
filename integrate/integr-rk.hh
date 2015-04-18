#ifndef integr_intgr_rk_hh
#define integr_intgr_rk_hh
/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

#include <iostream>
#include <sstream>
#include <valarray>
#include <deque>
#include <vector>

#ifdef DEBUG_INTEGR
#define DEBUG_RK
#endif

#ifdef DEBUG_RK
#include "../../spline/spline.hh"
#define DEB_RK(x) std::cerr << "rk: " << x << "\n";
#else
#define DEB_RK(x)
#endif

enum Formula {
  EXPLICIT,
  VPZM,
  HEUNM2,
  HEUNM3,
  SIMPSON,
  BS2, 
  BS3,
  RK4_KLASSISCH,
  RK4_3_8,
  CK4,
  CK5,
  RKF4,
  RKF5,
  DP4,
  DP5
};
 
template<enum Formula F>
struct IntegratorCoeffs {
  typedef Formula FormulaType;
};

template<>
struct IntegratorCoeffs<EXPLICIT> {
  static Formula const formula = EXPLICIT;
  static int const NC = 1;
  static double const a[NC];
  static double const b[NC][NC-1];
  static double const c[NC];
};
double const IntegratorCoeffs<EXPLICIT>::a[NC] = {
  1
};
double const IntegratorCoeffs<EXPLICIT>::b[NC][NC-1] = {
};
double const IntegratorCoeffs<EXPLICIT>::c[NC] = {
  1
};


template<>
struct IntegratorCoeffs<VPZM> {
  static Formula const formula = VPZM;
  static int const NC = 2;
  static double const a[NC];
  static double const b[NC][NC-1];
  static double const c[NC];
};

double const IntegratorCoeffs<VPZM>::a[NC] = {
  0, 0.5
};
double const IntegratorCoeffs<VPZM>::b[NC][NC-1] = {
  {},
  {0.5}
};
double const IntegratorCoeffs<VPZM>::c[NC] = {
  0, 1
};




template<>
struct IntegratorCoeffs<HEUNM2> {
  static Formula const formula = HEUNM2;
  static int const NC = 2;
  static double const a[NC];
  static double const b[NC][NC-1];
  static double const c[NC];
};

double const IntegratorCoeffs<HEUNM2>::a[NC] = {
  0, 1
};
double const IntegratorCoeffs<HEUNM2>::b[NC][NC-1] = {
  {},
  {1}
};
double const IntegratorCoeffs<HEUNM2>::c[NC] = {
  0.5, 0.5
};


template<>
struct IntegratorCoeffs<HEUNM3> {
  static Formula const formula = HEUNM3;
  static int const NC = 3;
  static double const a[NC];
  static double const b[NC][NC-1];
  static double const c[NC];
};

double const IntegratorCoeffs<HEUNM3>::a[NC] = {
  0, (1.0/3.0), (2.0/3.0)
};
double const IntegratorCoeffs<HEUNM3>::b[NC][NC-1] = {
  {},
  {(1.0/3.0)},
  {0, (2.0/3.0)}
};
double const IntegratorCoeffs<HEUNM3>::c[NC] = {
  0.25, 0, 0.75
};


template<>
struct IntegratorCoeffs<SIMPSON> {
  static Formula const formula = SIMPSON;
  static int const NC = 3;
  static double const a[NC];
  static double const b[NC][NC-1];
  static double const c[NC];
};

double const IntegratorCoeffs<SIMPSON>::a[NC] = {
  0, 0.5, 1
};
double const IntegratorCoeffs<SIMPSON>::b[NC][NC-1] = {
  {0, 0},
  {0.5, 0},
  {-1, 2}
};
double const IntegratorCoeffs<SIMPSON>::c[NC] = {
  1.0/6, 4.0/6, 1.0/6
};


template<>
struct IntegratorCoeffs<RK4_KLASSISCH> {
  static Formula const formula = RK4_KLASSISCH;
  static int const NC = 4;
  static double const a[NC];
  static double const b[NC][NC-1];
  static double const c[NC];
};

double const IntegratorCoeffs<RK4_KLASSISCH>::a[NC] = {
  0, 0.5, 0.5, 1
};
double const IntegratorCoeffs<RK4_KLASSISCH>::b[NC][NC-1] = {
  {},
  {0.5},
  {0, 0.5},
  {0, 0, 1}
};
double const IntegratorCoeffs<RK4_KLASSISCH>::c[NC] = {
  1.0/6, 2.0/6, 2.0/6, 1.0/6
};



template<>
struct IntegratorCoeffs<RK4_3_8> {
  static Formula const formula = RK4_3_8;
  static int const NC = 4;
  static double const a[NC];
  static double const b[NC][NC-1];
  static double const c[NC];
};
double const IntegratorCoeffs<RK4_3_8>::a[NC] = {
  0, 1.0/3, 2.0/3, 1
};
double const IntegratorCoeffs<RK4_3_8>::b[NC][NC-1] = {
  {},
  {1.0/3},
  {-1.0/3, 1},
  {1, -1, 1}
};
double const IntegratorCoeffs<RK4_3_8>::c[NC] = {
  1.0/8, 3.0/8, 3.0/8, 1.0/8
};


/// Coefficients for 2nd order Bogacki-Shampine method
/// see http://en.wikipedia.org/wiki/List_of_Runge-Kutta_methods
template<>
struct IntegratorCoeffs<BS2> {
  static Formula const formula = BS2;
  static int const NC = 4;
  static double const a[NC];
  static double const b[NC][NC-1];
  static double const c[NC];
};
double const IntegratorCoeffs<BS2>::a[NC] = {
  0, 0.5, 3.0/4, 1
};
double const IntegratorCoeffs<BS2>::b[NC][NC-1] = {
  {},
  { 0.5 },
  { 0,     3.0/4},
  { 2.0/9, 1.0/3, 4.0/9},
};
double const IntegratorCoeffs<BS2>::c[NC] = {
  7.0/24, 1.0/4, 1.0/3, 1.0/8
};

/// Coefficients for 3rd order Bogacki-Shampine method
/// identical to those of 2nd order version, but c coefficients
/// (bottom line of tableau) are different.
template<>
struct IntegratorCoeffs<BS3> : public IntegratorCoeffs<BS2> {
  static Formula const formula = BS3;
  static int const NC = 4;
  static double const c[NC];
};
double const IntegratorCoeffs<BS3>::c[NC] = {
  2.0/9, 1.0/3, 4.0/9, 0
};


/// Coefficients for 4th order Runge-Kutta-Fehlberg method
/// see http://en.wikipedia.org/wiki/Runge-Kutta-Fehlberg_method
template<>
struct IntegratorCoeffs<RKF4> {
  static Formula const formula = RKF4;
  static int const NC = 6;
  static double const a[NC];
  static double const b[NC][NC-1];
  static double const c[NC];
};
double const IntegratorCoeffs<RKF4>::a[NC] = {
  0, 1.0/4, 3.0/8, 12.0/13, 1, 1.0/2
};
double const IntegratorCoeffs<RKF4>::b[NC][NC-1] = {
  {},
  { 1.0/4},
  { 3.0/32,       9.0/32},
  { 1932.0/2197, -7200.0/2197,  7296.0/2197},
  { 439.0/216,   -8,            3680.0/513,    -845.0/4104 },
  {-8.0/27,      2,             -3544.0/2565,    1859.0/4104,  -11.0/40 }
};
double const IntegratorCoeffs<RKF4>::c[NC] = {
  25.0/216, 0, 1408.0/2565, 2197.0/4104, -1.0/5, 0
};

/// Coefficients for 5th order Runge-Kutta-Fehlberg method are
/// identical to those of 4th order version, but c coefficients
/// (bottom line of tableau) are different.
///
/// see http://en.wikipedia.org/wiki/Runge-Kutta-Fehlberg_method
template<>
struct IntegratorCoeffs<RKF5> : public IntegratorCoeffs<RKF4> {
  static Formula const formula = RKF5;
  static int const NC = 6;
  static double const c[NC];
};
double const IntegratorCoeffs<RKF5>::c[NC] = {
  16.0/135, 0, 6656.0/12825, 28561.0/56430, -9.0/50, 2.0/55
};



/// Coefficients for 4th order Cash-Karp method
/// see http://en.wikipedia.org/wiki/Cash-Karp
/// or the original paper http://dx.doi.org/10.1145/79505.79507
template<>
struct IntegratorCoeffs<CK4> {
  static Formula const formula = CK4;
  static int const NC = 6;
  static double const a[NC];
  static double const b[NC][NC-1];
  static double const c[NC];
};
double const IntegratorCoeffs<CK4>::a[NC] = {
  0, 1.0/5, 3.0/10, 3.0/5, 1, 7.0/8
};
double const IntegratorCoeffs<CK4>::b[NC][NC-1] = {
  {},
  {1.0/5},
  {3.0/40,       9.0/40},
  {3.0/10,      -9.0/10,    6.0/5},
  {-11.0/54,     5.0/2,     -70.0/27,    35.0/27 },
  {1631.0/55296, 175.0/512, 575.0/13824, 44275.0/110592, 253.0/4096 }
};
double const IntegratorCoeffs<CK4>::c[NC] = {
  2825.0/27648, 0, 18575.0/48384, 13525.0/55296, 277.0/14336, 1.0/4
};

/// Coefficients for 5th order Cash-Karp method are identical to those
/// of 4th order version, but c coefficients (bottom line of tableau)
/// are different
///
/// see http://en.wikipedia.org/wiki/Cash-Karp
/// or the original paper http://dx.doi.org/10.1145/79505.79507
template<>
struct IntegratorCoeffs<CK5> : public IntegratorCoeffs<CK4> {
  static Formula const formula = CK5;
  static int const NC = 6;
  static double const c[NC];
};
double const IntegratorCoeffs<CK5>::c[NC] = {
  37.0/378, 0, 250.0/621, 125.0/594, 0, 512.0/1771
};



/// Coefficients for 4th order Dormand-Prince method
/// see http://en.wikipedia.org/wiki/Dormand-Prince
/// or http://dx.doi.org/10.1016/0771-050X(80)90013-3
template<>
struct IntegratorCoeffs<DP4> {
  static Formula const formula = DP4;
  static int const NC = 7;
  static double const a[NC];
  static double const b[NC][NC-1];
  static double const c[NC];
};
double const IntegratorCoeffs<DP4>::a[NC] = {
  0, 1.0/5, 3.0/10, 4.0/5, 8.0/9, 1, 1
};
double const IntegratorCoeffs<DP4>::b[NC][NC-1] = {
  {},
  { 1.0/5},
  { 3.0/40,         9.0/40},
  { 44.0/45,       -56.0/15,        32.0/9},
  { 19372.0/6561,  -25360.0/2187,   64448.0/6561,   -212.0/729 },
  { 9017.0/3168,   -355.0/33,       46732.0/5247,    49.0/176,    -5103.0/18656 },
  { 35.0/384,         0,            500.0/1113,        125.0/192,   -2187.0/6784,       11.0/84 }
};
double const IntegratorCoeffs<DP4>::c[NC] = {
  5179.0/57600, 0, 7571.0/16695,   393.0/640, -92097.0/339200, 187.0/2100, 1.0/40
};

/// Coefficients for 5th order Dormand-Prince method are
/// identical to those of 4th order version, but c coefficients
/// (bottom line of tableau) are different.
template<>
struct IntegratorCoeffs<DP5> : public IntegratorCoeffs<DP4> {
  static Formula const formula = DP5;
  static int const NC = 7;
  static double const c[NC];
};
double const IntegratorCoeffs<DP5>::c[NC] = { 
  35.0/384, 0, 500.0/1113,   125.0/192,  -2187.0/6784,    11.0/84, 0 
};



template<class Fn, enum Formula F=VPZM, class _IntegratorCoeffs=IntegratorCoeffs<F> >
struct RKIntegrator : public IntegratorBase<Fn> {
  typedef Fn ForceFunctor; 
  typedef _IntegratorCoeffs IntegratorCoeffs;
  typedef typename Fn::value_type value_type;
  typedef typename Fn::first_argument_type  first_argument_type;
  typedef typename Fn::second_argument_type second_argument_type;
  typedef typename Fn::result_type result_type;
  typedef typename Fn::result_type::value_type element_type;

  class NotConsistent {};

  static int const NC = IntegratorCoeffs::NC;
  
//   value_type _states[NC-1];
//   value_type _forces[NC-1];

  std::vector<result_type> forces;

  Fn &fn;

  static char const *formulaName(Formula const f) {
    switch(f) {
    case EXPLICIT:
      return "explicit";
    case VPZM:
      return "VPZM";
    case HEUNM2:
      return "HEUNM2";
    case HEUNM3:
      return "HEUNM3";
    case BS2:
      return "BS2";
    case BS3:
      return "BS3";
    case SIMPSON:
      return "SIMPSON";
    case RK4_KLASSISCH:
      return "RK4_KLASSISCH";
    case RK4_3_8:
      return "RK4_3_8";
    case CK4:
      return "CK4";
    case CK5:
      return "CK5";
    case RKF4:
      return "RKF4";
    case RKF5:
      return "RKF5";
    case DP4:
      return "DP4";
    case DP5:
      return "DP5";
    }
    return "invalid formula";
  }

  static char const *formulaDescription(Formula const f) {
    switch(f) {
    case EXPLICIT:
      return "explicit euler method";
    case VPZM:
      return "trapezoidal rule, 2nd order Runge-Kutta method";
    case HEUNM2:
      return "2nd order Heun method";
    case HEUNM3:
      return "3rd order Heun method";
    case BS2:
      return "2nd order Bogacki-Shampine method";
    case BS3:
      return "3rd order Bogacki-Shampine method";
    case SIMPSON:
      return "Simpson rule, 3rd order Runge-Kutta method";
    case RK4_KLASSISCH:
      return "classic 4th order Runge-Kutta method Runge-Kutta";
    case RK4_3_8:
      return "4th order Runge-Kutta 3/8 rule";
    case CK4:
      return "Cash-Karp fourth order method";
    case CK5:
      return "Cash-Karp fifth order method";
    case RKF4:
      return "fourth order Runge-Kutta-Fehlberg method";
    case RKF5:
      return "fifth order Runge-Kutta-Fehlberg method";
    case DP4:
      return "fourth order Dormand-Prince method";
    case DP5:
      return "fifth order Dormand-Prince method";
    }
    return "invalid formula";
  }

//   RKIntegrator() : fn() {  }
  RKIntegrator(Fn &fn) : fn(fn) { 
    if (not consistent()) {
      throw NotConsistent();
    }
  }

  bool nearEqual(double v1, double v2) const {
    double const d = fabs(v2 - v1);
    return d < 1e-13;
  }

  bool consistent() const {
    bool ok = true;
    for (int si = 1; si < NC; ++si) {
      double sumbij = 0;
      for (int sj = 0; sj < si; ++sj) {
        sumbij += IntegratorCoeffs::b[si][sj];
      }
      if (not nearEqual(sumbij, IntegratorCoeffs::a[si])) {
        std::cerr << "error: " << formulaName(IntegratorCoeffs::formula)
                  << " sum of b-row " << si << " != a-coeff " << si
                  << ": " << sumbij << " != " << IntegratorCoeffs::a[si] << "\n";
        ok = false;
      }
    }
    return ok;
  }

  result_type operator()(first_argument_type const dt, second_argument_type const &state) {
    second_argument_type const &_state_0 = state;
    DEB_RK("integrate: NC " << NC << " Type " << IntegratorCoeffs::formula);
    DEB_RK("integrate: state y_k is " << _state_0);
    forces.clear();
    result_type update = fn(dt*IntegratorCoeffs::a[0], _state_0);
    forces.push_back( update );
    DEB_RK("integrate: update = " << IntegratorCoeffs::c[0]
	   << " * " << update);
    update *= element_type(dt * IntegratorCoeffs::c[0]);
    DEB_RK("integrate: update = " << update);
    for (int si = 1; si < NC; ++si) {
      DEB_RK("si " << si);
      value_type states_si = value_type(state);
      DEB_RK("y*[" << si << "] = " << states_si);
      for (int bj = 0; bj < si; ++bj) {
	DEB_RK("bj " << bj);
// 	result_type del = forces[bj];
// 	del *= element_type(dt * IntegratorCoeffs::b[si][bj]);
 	DEB_RK("y*[" << si << "] += dt * " << IntegratorCoeffs::b[si][bj]
	       << " * " << forces[bj]);
	fn.update(states_si, forces[bj] * element_type(dt * IntegratorCoeffs::b[si][bj]));
      }
      DEB_RK("y*[" << si << "] = " << states_si);
      forces.push_back( fn(dt*IntegratorCoeffs::a[si], states_si) );
      DEB_RK("f(t + dt*" << IntegratorCoeffs::a[si] << ", y*[" << si << "]) = "
	     << forces[si]);
//       result_type updateUpdate = forces[si];
//       updateUpdate *= element_type(dt * IntegratorCoeffs::c[si]);
//       update += updateUpdate;
      update += forces[si] * element_type(dt * IntegratorCoeffs::c[si]);
      DEB_RK("integrate: update += " << IntegratorCoeffs::c[si]
	     << " * " << forces[si] << "\n"
	     "integrate: state y_k+1 = y_k + " << update);
    }
    return update;
  }

  void update(second_argument_type &state, result_type const &update) {
    fn.update(state, update);
  }

  std::string name() const {
    std::ostringstream str;
    str << "SingleStep(" << formulaName(IntegratorCoeffs::formula) << ")";
    return str.str();
  }

  std::string description() const {
    std::ostringstream str;
    str << "Single-step Method `" 
	<< formulaName(IntegratorCoeffs::formula)
	<< "' of order " << IntegratorCoeffs::NC 
	<< " (" << formulaDescription(IntegratorCoeffs::formula) << ")"
	<< "\n"
      ;
    str << " coefficients a: [" ;
    for (int i = 0; i < IntegratorCoeffs::NC; ++i) {
      if (i) str << ", ";
      str << IntegratorCoeffs::a[i];
    }
    str << "]\n";
    str << " coefficients b: [\t" ;
    for (int i = 0; i < IntegratorCoeffs::NC; ++i) {
      if (i) str << "\n\t\t";
      for (int j = 0; j < IntegratorCoeffs::NC - 1; ++j) {
	if (j) str << ", ";
	str << IntegratorCoeffs::b[i][j];
      }
    }
    str << "]\n";
    str << " coefficients c: [" ;
    for (int i = 0; i < IntegratorCoeffs::NC; ++i) {
      if (i) str << ", ";
      str << IntegratorCoeffs::c[i];
    }
    str << "]";
    return str.str();
  }

  void printXML(std::ostream &str) const {
    str << "<single-step-method name='" 
	<< formulaName(IntegratorCoeffs::formula)
	<< "' order='" << IntegratorCoeffs::NC 
	<< "' descr='" << formulaDescription(IntegratorCoeffs::formula) << "'>"
	<< "\n"
      ;
    str << "<coefficients>\n\t" ;
    for (int i = 0; i < IntegratorCoeffs::NC; ++i) {
      str << "<a> " << IntegratorCoeffs::a[i] << " </a>";
    }
    for (int i = 1; i < IntegratorCoeffs::NC; ++i) {
      str << "\n\t<b-row>";
      for (int j = 0; j < IntegratorCoeffs::NC - 1; ++j) {
	str << "<b> " << IntegratorCoeffs::b[i][j] << " </b>";
      }
      str << "</b-row>";
    }
    str << "\n\t";
    for (int i = 0; i < IntegratorCoeffs::NC; ++i) {
      str << "<c> " << IntegratorCoeffs::c[i] << " </c>";
    }
    str << "\n</coefficients>";
    str << "\n</single-step-method>";
  }
};

#endif
