/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

struct Timer {
  double m_t, m_start, m_lastDiff;

  Timer() : m_t(), m_start(), m_lastDiff() {}
  Timer(bool) : m_t(), m_start(omp_get_wtime()), m_lastDiff() {}
  
  void stop() {
    m_lastDiff = omp_get_wtime() - m_start;
    m_t += m_lastDiff;
  }

  void start() {
    m_start = omp_get_wtime();
  }

  void restart() { reset(); }

  void reset() {
    m_t = 0;
    start();
  }
  
  void print(std::ostream &aus) const {
    aus << m_t << " s";
  }

  double lastDiff() const { return m_lastDiff; }
  double operator()() const { return m_t; }

};

std::ostream &operator << (std::ostream &aus, Timer const &t) {
  t.print(aus);
  return aus;
}
