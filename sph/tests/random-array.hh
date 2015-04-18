/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/
#include <valarray>

struct RandomArray {

  std::valarray<double> m_data;
  
  RandomArray(size_t const n, size_t const range) : m_data(n) {
    unsigned short randSeed[3] = {(unsigned short)time(0), (unsigned short)time(0), (unsigned short)time(0)};
    for(size_t i = 0; i < n; ++i) {
      m_data[i] = nrand48(randSeed) % range;
    }
  }

  size_t size() const { return m_data.size(); }
  double operator[](size_t const i) const { return m_data[i]; }

};

struct DRandomArray {
  
  std::valarray<double> m_data;
  
  DRandomArray(size_t const n) : m_data(n) {
    unsigned short randSeed[3] = {(unsigned short)time(0), (unsigned short)time(0), (unsigned short)time(0)};
    for(size_t i = 0; i < n; ++i) {
      m_data[i] = erand48(randSeed);
    }
  }

  size_t size() const { return m_data.size(); }
  double operator[](size_t const i) const { return m_data[i]; }

};
