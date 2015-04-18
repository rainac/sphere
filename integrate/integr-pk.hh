#ifndef integr_integr_amb_hh
#define integr_integr_amb_hh
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
#include "integr-abm.hh"
#include "integr-amm.hh"

// #include "../../spline/spline.hh"

template<size_t M=2>
struct AbbruchKriteriumM {
  template<class T>
  bool operator()(size_t m, T const &) {
    return m+1 >= M;
  }
};

template<class Fn, 
	 size_t S=2, 
	 size_t M=1, 
	 class Praediktor=ABMIntegrator<Fn, S>,
	 class Korrektor=AMMIntegrator<Fn, S>,
	 class AbbruchKriterium=AbbruchKriteriumM<M>
	 >
struct PKIntegrator : public IntegratorBase<Fn> {

  typedef typename IntegratorBase<Fn>::ForceFunctor ForceFunctor;
  typedef typename IntegratorBase<Fn>::value_type value_type;
  typedef typename IntegratorBase<Fn>::first_argument_type first_argument_type;
  typedef typename IntegratorBase<Fn>::second_argument_type second_argument_type;
  typedef typename IntegratorBase<Fn>::result_type result_type;
  typedef typename IntegratorBase<Fn>::element_type element_type;

  Praediktor praediktor;
  Korrektor korrektor;
  AbbruchKriterium abbruch;
  Fn &fn;

  static int const s = S;
  static int const m = M;

  PKIntegrator(Fn &fn) : 
    praediktor(fn),
    korrektor(fn),
    fn(fn) 
  {}

  result_type operator()(first_argument_type const dt, second_argument_type const &state) {

//     cerr << "PKIntegr(" << s << " + 1)(" << dt << ", " << state << ")\n";
//     std::cerr << "PKIntegr(" << s << " + 1)(" << dt << "): predict\n";

    result_type update = praediktor(dt, state);
    value_type future = value_type(state);
    fn.update(future, update);

    for(size_t i = 0; ; ++i) {
//       std::cerr << "PKIntegr(" << s << " + 1)(" << dt << "): correct\n";

      // pass forces at state computed by predictor on to corrector
      // so it does not compute them again
      update = korrektor(dt, state, future, praediktor.forces[0]);
      //       update = korrektor(dt, state, future);

      if (abbruch(i, update)) break;
      future = state;
      fn.update(future, update);
      
    }

    korrektor.accept();

    return update;
  }  

  void update(second_argument_type &state, result_type const &update) {
    fn.update(state, update);
  }

  std::string name() const {
    std::ostringstream str;
    str << "PK_" << S << "_" << M
	<< "(" << praediktor.name() << "," << korrektor.name() << ")";
    return str.str();
  }

  std::string description() const {
    std::ostringstream str;
    str << "Predictor-Corrector with " << S << "-level history and "
	<< M << " correction steps"
	<< ",\n predictor: " << praediktor.description()
	<< ",\n corrector: " << korrektor.description()
      ;
    return str.str();
  }

  void printXML(std::ostream &str) const {
    str << "<predictor-corrector"
      " level='" << S << "'"
      " corrections='" << M << "'>\n\t"
	<< "<predictor>\n" << praediktor << "\n</predictor>\n"
	<< "<corrector>\n" << korrektor << "\n</corrector>\n"
	<< "</predictor-corrector>";
  }

};

#endif
