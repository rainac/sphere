#ifndef su_time_hh
#define su_time_hh
/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

#include <time.h>

struct Time {
  time_t m_time;

  Time() : m_time(time(0)) {}
  Time(time_t _time) : m_time(_time) {}

  Time operator +(double const d) const { return Time(m_time + d); }
  Time operator -(double const d) const { return Time(m_time - d); }
  operator time_t() const { return m_time; }
};

#endif
