#ifndef sph_counter_hh
#define sph_counter_hh
/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

/// Template class Counter uses a vector of VectorType and a
/// compile-time constant Incr to implement a ndim-dimensional counter.
/// \tparam Vector Type type of the vector
/// \tparam Incr Integer increment

template<class VectorType, typename VectorType::value_type Incr, int NDIM>
struct Counter {
  VectorType v;
  VectorType top;

  /// The constructor 
  /// \param top a vector which holds the maximum value of each coordinate.
  Counter(VectorType const &top) : v(NDIM), top(top) {}

  Counter &operator ++() {
    for (size_t i = 0; i < NDIM; ++i) {
      v[i] += Incr;
      if (v[i] >= top[i]) {
	v[i] = 0;
      } else {
	break;
      }
    }
    return *this;
  }

  bool operator ==(Counter const &o) const { 
    return (v == o.v).min() == 1; 
  }
};

#endif
