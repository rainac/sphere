#ifndef integr_integr_abm_hh
#define integr_integr_abm_hh
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

// #include "../../spline/spline.hh"

#ifdef DEBUG
#define DEB_ABM(x) std::cerr << "abm: " << x << "\n";
#else
#define DEB_ABM(x)
#endif

template<int S>
struct ABMCoeffs {
  static int const s = S;
  static double const b[s + 1];
};

#include "coeffs-abm.hh"
  
template<class Fn, 
	 size_t S=2, 
	 class Coeffs=ABMCoeffs<S>,
	 class BootstrapIntegrator=RKIntegrator<Fn, CK5> >
struct ABMIntegrator : public IntegratorBase<Fn> {
  typedef typename IntegratorBase<Fn>::ForceFunctor ForceFunctor;
  typedef typename IntegratorBase<Fn>::value_type value_type;
  typedef typename IntegratorBase<Fn>::first_argument_type first_argument_type;
  typedef typename IntegratorBase<Fn>::second_argument_type second_argument_type;
  typedef typename IntegratorBase<Fn>::result_type result_type;
  typedef typename IntegratorBase<Fn>::element_type element_type;

  Fn &fn;
  BootstrapIntegrator bootstrapIntegrator;

  static int const s = Coeffs::s;

  std::deque<result_type> forces;

//   ABMIntegrator() : fn(), bootstrapIntegrator(fn), forces() {  
//   }
  ABMIntegrator(Fn &fn) : fn(fn), bootstrapIntegrator(fn), forces() { }

  result_type operator()(first_argument_type const dt, second_argument_type const &state) {

    DEB_ABM("ABMIntegr(" << s << " + 1)(" << dt << ", " << state
	    << (forces.size() < size_t(s) ? ", init":"") << ")");

    if (forces.size() < size_t(s)) {
      return init(dt, state);
    }

    forces.push_front( fn(dt, state) );
    result_type update = forces[0];
    DEB_ABM("update = " << Coeffs::b[0] << " * " << update);
    update *= element_type(Coeffs::b[0]);

    for(int i = 1; i <= s; ++i) {
      DEB_ABM("update += " << Coeffs::b[i] << " * " << forces[i]);
//       result_type updateUpdate = forces[i];
//       updateUpdate *= Coeffs::b[i];
//       update += updateUpdate;
      update += forces[i] * element_type(Coeffs::b[i]);
    }
    update *= element_type(dt);

    assert(forces.size() == s + 1);
    forces.pop_back();
    assert(forces.size() == size_t(s));

    return update;
  }  

  void update(second_argument_type &state, result_type const &update) {
    fn.update(state, update);
  }

  std::string name() const {
    std::ostringstream str;
    str << "ABM_" << S
	<< "(" << bootstrapIntegrator.name() << ")";
    return str.str();
  }

  std::string description() const {
    std::ostringstream str;
    str << "Adams-Bashforth-Method (ABM) with " << S << "-level history\n"
      ;
    str << " coefficients: [" ;
    for (size_t i = 0; i <= S; ++i) {
      if (i) str << ", ";
      str << Coeffs::b[i];
    }
    str << "],\n"
	<< " bootstrap integrator: " << bootstrapIntegrator.description();
    return str.str();
  }

  void printXML(std::ostream &str) const {
    str << "<adams-bashforth-method level='" << S << "'>\n\t";
    str << "<coefficients>" ;
    for (size_t i = 0; i <= S; ++i) {
      str << "<c> " << Coeffs::b[i] << " </c>";
    }
    str << "</coefficients>\n"
	<< "<bootstrap-integrator>\n"
	<< bootstrapIntegrator
	<< "\n</bootstrap-integrator>\n"
	<< "</adams-bashforth-method>";
  }

private:
  result_type init(first_argument_type const dt, second_argument_type const &state0) {
    forces.push_front( fn(dt, state0) );
    DEB_ABM("ABMIntegr<" << s << ">::init() push_front " << forces.front());
    /// \todo pass forces computed also on to bootstrap integrator...
    result_type update = bootstrapIntegrator(dt, state0);
    return update;
  }

};

#endif
