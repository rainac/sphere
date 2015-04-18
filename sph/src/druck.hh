#include "sph-config.hh"
#include "partikel.hh"
/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

/// Functional struct representing the Gamme law of pressure. It can
/// compute and update the pressure of a particle from its density
/// value and global material properties.
template<class T=double>
struct GammaGleichung {
  typedef T value_type;

  SPHEQSConfig<value_type> const &config;

  GammaGleichung(SPHEQSConfig<value_type> const &config) : config(config) {}

  void updatePressure(Partikel<value_type> &p) const {
    p.druck() = computePressure(p);
  }

  value_type computePressure(Partikel<value_type> const &p) const {
    return config.pressureB *
      (dynamic_power(p.dichte() / config.rho0, config.gamma) - 1);
  }

};
