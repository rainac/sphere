#ifndef integr_integrate3_hh
#define integr_integrate3_hh
/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

#include "integr-base.hh"
#include "integr-abm.hh"
#include "integr-amm.hh"
#include "integr-rk.hh"
#include "integr-pk.hh"

template<class Fn>
struct Integrator {

  typedef typename IntegratorBase<Fn>::ForceFunctor ForceFunctor;
  typedef typename IntegratorBase<Fn>::value_type value_type;
  typedef typename IntegratorBase<Fn>::first_argument_type  first_argument_type;
  typedef typename IntegratorBase<Fn>::second_argument_type second_argument_type;
  typedef typename IntegratorBase<Fn>::result_type result_type;
  typedef typename IntegratorBase<Fn>::element_type element_type;

  IntegratorBase<Fn> *m_integrator;

private:
  /// \todo implement Copy-Constructor and assignment using
  /// m_integrator->clone or something
  Integrator & operator =(Integrator const &) { 
    return *this;
  }
  Integrator(Integrator const &) { }

public:
  Integrator(Fn &fn, std::string const &integrName) {
    m_integrator = makeIntegrator(fn, integrName);
  }

  ~Integrator() {
    if (m_integrator) {
      delete m_integrator;
      m_integrator = 0;
    }
  }

  static IntegratorBase<Fn> *makeIntegrator(Fn &fn, std::string const &integrName) {
    IntegratorBase<Fn> *result = 0;
    if (integrName == "explicit") {
      result = new RKIntegrator<Fn, EXPLICIT>(fn);

    } else if (integrName == "vpzm" 
               or integrName == "midpoint" 
               or integrName == "rk2") {
      result = new RKIntegrator<ForceFunctor, VPZM>(fn);
    } else if (integrName == "heun2") {
      result = new RKIntegrator<ForceFunctor, HEUNM2>(fn);

    } else if (integrName == "simpson"
               or integrName == "rk3") {
      result = new RKIntegrator<ForceFunctor, SIMPSON>(fn);
    } else if (integrName == "heun3") {
      result = new RKIntegrator<ForceFunctor, HEUNM3>(fn);

    } else if (integrName == "bs2") {
      result = new RKIntegrator<ForceFunctor, BS2>(fn);
    } else if (integrName == "bs3") {
      result = new RKIntegrator<ForceFunctor, BS3>(fn);

    } else if (integrName == "rk4") {
      result = new RKIntegrator<ForceFunctor, RK4_KLASSISCH>(fn);
    } else if (integrName == "rk4_3_8") {
      result = new RKIntegrator<ForceFunctor, RK4_3_8>(fn);

    } else if (integrName == "rkf4") {
      result = new RKIntegrator<ForceFunctor, RKF4>(fn);
    } else if (integrName == "rkf5") {
      result = new RKIntegrator<ForceFunctor, RKF5>(fn);

    } else if (integrName == "ck4") {
      result = new RKIntegrator<ForceFunctor, CK4>(fn);
    } else if (integrName == "ck5") {
      result = new RKIntegrator<ForceFunctor, CK5>(fn);

    } else if (integrName == "dp4") {
      result = new RKIntegrator<ForceFunctor, DP4>(fn);
    } else if (integrName == "dp5") {
      result = new RKIntegrator<ForceFunctor, DP5>(fn);

    } else if (integrName == "abm1") {
      result = new ABMIntegrator<ForceFunctor, 1>(fn);
    } else if (integrName == "abm2") {
      result = new ABMIntegrator<ForceFunctor, 2>(fn);
    } else if (integrName == "abm3") {
      result = new ABMIntegrator<ForceFunctor, 3>(fn);
    } else if (integrName == "abm4") {
      result = new ABMIntegrator<ForceFunctor, 4>(fn);
    } else if (integrName == "abm5") {
      result = new ABMIntegrator<ForceFunctor, 5>(fn);
    } else if (integrName == "abm6") {
      result = new ABMIntegrator<ForceFunctor, 6>(fn);
    } else if (integrName == "abm7") {
      result = new ABMIntegrator<ForceFunctor, 7>(fn);
    } else if (integrName == "abm8") {
      result = new ABMIntegrator<ForceFunctor, 8>(fn);
    } else if (integrName == "abm9") {
      result = new ABMIntegrator<ForceFunctor, 9>(fn);
    } else if (integrName == "abm10") {
      result = new ABMIntegrator<ForceFunctor, 10>(fn);
    } else if (integrName == "abm11") {
      result = new ABMIntegrator<ForceFunctor, 11>(fn);
    } else if (integrName == "abm12") {
      result = new ABMIntegrator<ForceFunctor, 12>(fn);

    } else if (integrName == "pk1_1") {
      result = new PKIntegrator<ForceFunctor, 1, 1>(fn);
    } else if (integrName == "pk1_2") {
      result = new PKIntegrator<ForceFunctor, 1, 2>(fn);
    } else if (integrName == "_pk1_3") {
      result = new PKIntegrator<ForceFunctor, 1, 3>(fn);
    } else if (integrName == "_pk1_4") {
      result = new PKIntegrator<ForceFunctor, 1, 4>(fn);
    } else if (integrName == "_pk1_5") {
      result = new PKIntegrator<ForceFunctor, 1, 5>(fn);

    } else if (integrName == "pk2_1") {
      result = new PKIntegrator<ForceFunctor, 2, 1>(fn);
    } else if (integrName == "pk2_2") {
      result = new PKIntegrator<ForceFunctor, 2, 2>(fn);
    } else if (integrName == "_pk2_3") {
      result = new PKIntegrator<ForceFunctor, 2, 3>(fn);
    } else if (integrName == "_pk2_4") {
      result = new PKIntegrator<ForceFunctor, 2, 4>(fn);
    } else if (integrName == "_pk2_5") {
      result = new PKIntegrator<ForceFunctor, 2, 5>(fn);

    } else if (integrName == "pk3_1") {
      result = new PKIntegrator<ForceFunctor, 3, 1>(fn);
    } else if (integrName == "pk3_2") {
      result = new PKIntegrator<ForceFunctor, 3, 2>(fn);
    } else if (integrName == "_pk3_3") {
      result = new PKIntegrator<ForceFunctor, 3, 3>(fn);
    } else if (integrName == "_pk3_4") {
      result = new PKIntegrator<ForceFunctor, 3, 4>(fn);
    } else if (integrName == "_pk3_5") {
      result = new PKIntegrator<ForceFunctor, 3, 5>(fn);

    } else if (integrName == "pk4_1") {
      result = new PKIntegrator<ForceFunctor, 4, 1>(fn);
    } else if (integrName == "pk4_2") {
      result = new PKIntegrator<ForceFunctor, 4, 2>(fn);
    } else if (integrName == "_pk4_3") {
      result = new PKIntegrator<ForceFunctor, 4, 3>(fn);
    } else if (integrName == "_pk4_4") {
      result = new PKIntegrator<ForceFunctor, 4, 4>(fn);
    } else if (integrName == "_pk4_5") {
      result = new PKIntegrator<ForceFunctor, 4, 5>(fn);

    } else if (integrName == "pk5_1") {
      result = new PKIntegrator<ForceFunctor, 5, 1>(fn);
    } else if (integrName == "pk5_2") {
      result = new PKIntegrator<ForceFunctor, 5, 2>(fn);
    } else if (integrName == "_pk5_3") {
      result = new PKIntegrator<ForceFunctor, 5, 3>(fn);
    } else if (integrName == "_pk5_4") {
      result = new PKIntegrator<ForceFunctor, 3, 4>(fn);
    } else if (integrName == "_pk5_5") {
      result = new PKIntegrator<ForceFunctor, 3, 5>(fn);

    } else if (integrName == "pk6_1") {
      result = new PKIntegrator<ForceFunctor, 6, 1>(fn);
    } else if (integrName == "pk6_2") {
      result = new PKIntegrator<ForceFunctor, 6, 2>(fn);
    } else if (integrName == "_pk6_3") {
      result = new PKIntegrator<ForceFunctor, 6, 3>(fn);
    } else if (integrName == "_pk6_4") {
      result = new PKIntegrator<ForceFunctor, 6, 4>(fn);
    } else if (integrName == "_pk6_5") {
      result = new PKIntegrator<ForceFunctor, 6, 5>(fn);

    } else if (integrName == "pk7_1") {
      result = new PKIntegrator<ForceFunctor, 7, 1>(fn);
    } else if (integrName == "pk8_1") {
      result = new PKIntegrator<ForceFunctor, 8, 1>(fn);
    } else if (integrName == "pk9_1") {
      result = new PKIntegrator<ForceFunctor, 9, 1>(fn);
    } else if (integrName == "pk10_1") {
      result = new PKIntegrator<ForceFunctor, 10, 1>(fn);
    } else if (integrName == "pk11_1") {
      result = new PKIntegrator<ForceFunctor, 11, 1>(fn);
    } else if (integrName == "pk12_1") {
      result = new PKIntegrator<ForceFunctor, 12, 1>(fn);

    } else {
      std::cerr << "error: no integrator available named `" << integrName << "'\n";
      exit(-4);
    }
    return result;
  }

  result_type operator()(first_argument_type const dt, second_argument_type const &state) {
    return m_integrator->operator()(dt, state);
  }

  void update(second_argument_type &state, result_type const &update) {
    return m_integrator->update(state, update);
  }

  std::string name() const {
    return m_integrator->name();
  }

  std::string description() const {
    return m_integrator->description();
  }

  void printXML(std::ostream &str) const {
    m_integrator->printXML(str);
  }

};

template<class A>
inline std::ostream &operator <<(std::ostream &str, Integrator<A> const &pk) {
  pk.printXML(str);
  return str;
}

#endif
