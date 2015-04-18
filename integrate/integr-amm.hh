#ifndef integr_integr_amm_hh
#define integr_integr_amm_hh
/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

#include <iostream>
#include <valarray>
#include <deque>
#include <vector>
#include <assert.h>
#include "integr-rk.hh"

#ifdef DEBUG_INTEGR
#include "../../spline/spline.hh"
#define DEB_AMM(x) std::cerr << "amm: " << x << "\n";
#else
#define DEB_AMM(x)
#endif

template<int S>
struct AMMCoeffs {
  static int const s = S;
  static double const b[s + 1];
};

// template<int S> double const AMMCoeffs<S>::b[s + 1] = {};

#include "coeffs-amm.hh"
/*
template<> double const AMMCoeffs<1>::b[s + 1] = {
  0.5, 0.5
};
template<> double const AMMCoeffs<2>::b[s + 1] = {
  5.0/12, 2.0 / 3, (-1.0)/12
};
template<> double const AMMCoeffs<3>::b[s + 1] = {
  3.0/8.0, 	     19.0/24.0, 	(-5.0)/24.0,
  1.0/24.0
};
template<> double const AMMCoeffs<4>::b[s + 1] = {
  251.0/720.0, 	     323.0/360.0, 	(-11.0)/30.0,
  53.0/360.0,        (-19.0)/720
};
template<> double const AMMCoeffs<5>::b[s + 1] = {
  95.0/288.0, 	     1427.0/1440.0, 	(-133.0)/240.0,
  241.0/720.0,       (-173.0)/1440,     3.0/160.0
};
template<> double const AMMCoeffs<6>::b[s + 1] = {
  19087.0/60480.0,  2713.0/2520.0,      (-15487.0)/20160.0,
  586.0/945.0,      (-6737.0)/20160.0,  263.0/2520.0,
  (-863.0)/60480.0
};
template<> double const AMMCoeffs<7>::b[s + 1] = {
  5257.0/17280.0,   139849.0/120960.0,   (-4511.0)/4480.0,
  123133.0/120960.0,(-88547.0)/120960.0, 1537.0/4480.0,
  (-11351.0)/120960.0, 275.0/24192.0
};
*/
template<class Fn, 
	 size_t S=2, 
	 class Coeffs=AMMCoeffs<S>,
	 class BootstrapIntegrator=RKIntegrator<Fn, CK5> >
struct AMMIntegrator {
  typedef Fn ForceFunctor; 

  Fn &fn;
  BootstrapIntegrator bootstrapIntegrator;

  typedef typename Fn::value_type value_type;
  typedef typename Fn::first_argument_type   first_argument_type;
  typedef typename Fn::second_argument_type second_argument_type;
  typedef typename Fn::result_type result_type;
  typedef typename Fn::result_type::value_type element_type;

  static int const s = Coeffs::s;

  std::deque<result_type> forces;

  bool initialized;

//   AMMIntegrator() : fn(), bootstrapIntegrator(fn), forces() {  
//   }
  AMMIntegrator(Fn &fn) : fn(fn), bootstrapIntegrator(fn), forces(), initialized() { }
  ~AMMIntegrator() {
    clear();
  }

  result_type operator()(first_argument_type const dt, second_argument_type const &state, 
			value_type const &future) {
    if (forces.size()+1 < size_t(s)) {
      return init(dt, state);
    }

    return operator()(dt, state, future, fn(dt, state));
  }

  void clear() {
    initialized = false;
    forces.clear();
  }

  void accept() {
    if (initialized) {
      forces.pop_back();
    }
  }

  result_type operator()(first_argument_type const dt, second_argument_type const &state, 
			 value_type const &future, result_type const &stateForce) {

    DEB_AMM("AMMIntegr(" << s << " + 1)(" << dt << ", " << state
	    << ", " << future
	    << (forces.size() < size_t(s-1) ? ", init":"") << ")");

    if (not initialized and forces.size() + 1 < size_t(s)) {
      return init(dt, state);
    }
    initialized = true;

    if (forces.size() < size_t(s)) {
      forces.push_front( stateForce );
    } else {
//       assert(forces.front() == stateForce);
    }

    result_type update = fn(dt, future) * element_type(Coeffs::b[0]);
//     DEB_AMM("update = " << Coeffs::b[0] << " * " << update);
//     update *= Coeffs::b[0];

    for(int i = 0; i < s; ++i) {
      DEB_AMM("update += " << Coeffs::b[i+1] << " * " << forces[i]);
//       result_type updateUpdate = forces[i];
//       updateUpdate *= Coeffs::b[i+1];
//       update += updateUpdate;
      update += forces[i] * element_type(Coeffs::b[i+1]);
    }
    update *= element_type(dt);

    assert(forces.size() == size_t(s));

//     forces.pop_back();
//     assert(forces.size() == size_t(s-1));

    return update;
  }  

  std::string name() const {
    std::ostringstream str;
    str << "AMM_" << S
	<< "(" << bootstrapIntegrator.name() << ")";
    return str.str();
  }

  std::string description() const {
    std::ostringstream str;
    str << "Adams-Moulton-Method (AMM) with " << S << "-level history\n"
      ;
    str << "coefficients: [" ;
    for (size_t i = 0; i <= S; ++i) {
      if (i) str << ", ";
      str << Coeffs::b[i];
    }
    str << "],\n"
	<< "bootstrap integrator: " << bootstrapIntegrator.description();
    return str.str();
  }

  void printXML(std::ostream &str) const {
    str << "<adams-moulton-method level='" << S << "'>\n\t";
    str << "<coefficients>" ;
    for (size_t i = 0; i <= S; ++i) {
      str << "<c> " << Coeffs::b[i] << " </c>";
    }
    str << "</coefficients>\n"
	<< "<bootstrap-integrator>\n"
	<< bootstrapIntegrator
	<< "\n</bootstrap-integrator>\n"
	<< "</adams-moulton-method>";
  }

private:
  result_type init(first_argument_type const dt, second_argument_type const &state0) {
    forces.push_front( fn(dt, state0) );
    DEB_AMM("AMMIntegr<" << s << ">::init() push_front " << forces.front());
    /// \todo pass forces computed also on to bootstrap integrator...
    result_type update = bootstrapIntegrator(dt, state0);
    return update;
  }

};

template<class A, size_t B, class D, class E>
inline std::ostream &operator <<(std::ostream &str, AMMIntegrator<A,B,D,E> const &pk) {
  pk.printXML(str);
  return str;
}

#endif
