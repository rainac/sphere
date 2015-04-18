#ifndef jw_sph_lj_hh
#define jw_sph_lj_hh
/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

#include "power.hh"

/// Similiar to SplineKernel2 and GaussianKernel, this computes the
/// the Lennard-Jones potential of a Lennard-Jones boundary particle a distance
// \f$ \bf r \f$.
template<class T=double>
struct LennardJones {
  typedef T value_type;

  // should be initial spacing betweeen particles
  value_type const r0;
  // coefficient: suggested value: 5 g H (H: height)
  value_type const D;

  // exponents: p2 = 2*p1
  long const p1;
  //   double p2; // value implied

  // cutoff: return LJ(cutoff) for all x < cutoff, else 0
  value_type const cutoff;

  LennardJones(value_type const r0, value_type const D, unsigned const p1) : 
    r0(r0)
    , D(D)
    , p1(p1)
//     , p2(2*p1)  fixed to 2*p1
//     , cutoff(r0 * pow(2.0, 1.0/p1))  // cutoff bei maximum
    , cutoff(r0)  // cutoff bei nullstelle
  {
//     std::cerr << "lennart-jones: r0 " << r0
// 	      << " D " << D << " p1 " << p1 << "\n";
  }

  // this returns just the weight, multiply by dist.-vector r yourself
  value_type w(value_type dnorm) const {
    if (dnorm > cutoff) {
      return 0;
    }
    value_type const prop = r0 / dnorm;
    value_type const propep1 = dynamic_power(prop, long(p1));
    value_type const propep2 = propep1*propep1;
    value_type const result = D*(propep2 - propep1);
    return result;
  }
};

#endif
