/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

#include <valarray>
#include <iostream>
#include <omp.h>
#include <assert.h>
#include "utility/timer.hh"

using namespace std;

struct DistArray {
  typedef std::valarray<double> Vector;

  size_t m_n, m_p;
  Vector *m_data;
  Timer m_timer;

  DistArray(size_t p, size_t n) : 
    m_n(n), 
    m_p(p), 
    m_data(new Vector[m_p]) {
    cout << "allocating " << m_p << "x" << m_n << " distrib. array\n";
#pragma omp parallel
    {
      int const myid = omp_get_thread_num();
      m_data[myid].resize(m_n);
      unsigned short rSeed[3];
      rSeed[0] = myid;
      rSeed[1] = myid;
      rSeed[2] = myid;
      for(size_t i = 0; i < m_n; ++i) {
        m_data[myid][i] = erand48(rSeed);
//         m_data[myid][i] = myid;
      }
    }
    m_timer.stop();
    cout << "distrib. allocation complete: " << m_timer << "\n";
  }
  
  void clear() {
    m_timer.start();
    if (m_data) {
      delete[] m_data;
      m_data = 0;
    }
    m_timer.stop();
    cout << "distrib. deallocation : " << m_timer << "\n";
  }

  Vector &operator[](size_t i) { return m_data[i]; }
  Vector const &operator[](size_t i) const { return m_data[i]; }

};

void test() {
}

int main(int const argc, char *argv[]) {
  size_t n =  1000000, nthreads = 128;

  if (argc > 1) {
    n = atoi(argv[1]);
  }
  if (argc > 2) {
    nthreads = atoi(argv[2]);
  }

  omp_set_num_threads(nthreads);
 
  DistArray distArray(nthreads, n);

  Timer timer;

  double totalSum = 0;

  cout << "read test\n";
  for(size_t i = 0; i < nthreads; ++i) {
    totalSum = 0;
    timer.reset();
#pragma omp parallel reduction(+:totalSum)
    {
      int const myid = omp_get_thread_num();
      int const offs = (myid + i) % nthreads;
      double s = distArray[offs].sum();
      double m = distArray[offs].min();
      double n = distArray[offs].max();
      totalSum += s + m + n;
    }
    
    timer.stop();
    cout << "summation of data of neighbour " << i << ": " << totalSum << ", time: " << timer << "\n";

    totalSum = 0;
    timer.reset();
#pragma omp parallel reduction(+:totalSum)
    {
      int const myid = omp_get_thread_num();
      int const offs = myid;
      double s = distArray[offs].sum();
      double m = distArray[offs].min();
      double n = distArray[offs].max();
      totalSum += s + m + n;
    }
    
    timer.stop();
    cout << "summation of own data " <<  i << ": " << totalSum << ", time: " << timer << "\n";
  }


  for(size_t i = 0; i < nthreads; ++i) {
    timer.reset();
#pragma omp parallel
    {
      int const myid = omp_get_thread_num();
      int const offs = (myid + i) % nthreads;
      for(size_t k = 0; k < n; ++k) {
        distArray[offs][k] = myid;
      }
    }
    
    timer.stop();
    cout << "setting data of neighbour " << i << ", time: " << timer << "\n";

    timer.reset();
#pragma omp parallel
    {
      int const myid = omp_get_thread_num();
      int const offs = myid;
      for(size_t k = 0; k < n; ++k) {
        distArray[offs][k] = myid;
      }
    }
   
    timer.stop();
    cout << "setting of own data " <<  i << ", time: " << timer << "\n";
  }
}
