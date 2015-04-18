/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <valarray>
#include <assert.h>
extern "C" char **environ;
#include "utility/getenv.hh"
#include <omp.h>
#include <algorithm>
// #include <parallel/algorithm>
// // #include <mcstl.h>
// // #include <bits/mcstl_sort.h>

#define SORT std::sort
// #define SORT __gnu_parallel::sort

/// compare two numbers with operator < "less"
struct CMPNum {
  /// the type of values to compare
  typedef int value_type;
  /// the type of the first argument is value_type
  typedef value_type first_argument_type;
  /// the type of the second argument is value_type
  typedef value_type second_argument_type;
  /// the type of the result is bool
  typedef bool result_type;
  /// the comparison function
  /// \param a the first object
  /// \param b the second object
  /// \returns a < b
  bool operator()(int const a, int const b) const { return a < b; }
};

/// compare two numbers with operator > "greater"
struct CMPNumRev {
  typedef int value_type;
  typedef value_type first_argument_type;
  typedef value_type second_argument_type;
  typedef bool result_type;
  bool operator()(int const a, int const b) const { return a > b; }
};

template<class T>
static typename T::value_type *address(T &v, size_t i) {
  typename T::value_type &firstItem = v[0];
  typename T::value_type * baseAdress = &firstItem;
//   std::cerr << "base: " << baseAdress << "\n";
//   std::cerr << "offset: " << i << " * " << sizeof(typename T::value_type)
// 	    << " = " << i * sizeof(typename T::value_type) << "\n";
  return baseAdress + i;
}

template<class T>
static typename T::value_type const *address(T const &v, size_t i) {
  typename T::value_type const &firstItem = v[0];
  typename T::value_type const * baseAdress = &firstItem;
  return baseAdress + i;
}

using namespace std;

void f(std::vector<int> &vec) {
  SORT (vec.begin(), vec.end(), CMPNumRev());
//   for (size_t i = 0; i < 100; ++i) {
//     cout << "v1[" << i << "] = " << *vec[i] << "\n";
//   }

  SORT (vec.begin(), vec.end() - 1, CMPNum());
//   for (size_t i = 0; i < 100; ++i) {
//     cout << "v1[" << i << "] = " << *vec[i] << "\n";
//   }
}

void g(std::valarray<int> &vec) {
  int *firstE = address(vec, 0);
  int *lastE = address(vec, vec.size() );
//   cout << "first at " << firstE << ", last at " << lastE << " diff " << lastE - firstE << "\n";
  SORT(firstE, lastE, CMPNumRev());
//   for (size_t i = 0; i < 100; ++i) {
//     cout << "v1[" << i << "] = " << *vec[i] << "\n";
//   }
  
//   mcstl::parallel_sort(address(vec, 0), address(vec, vec.size()), CMPNum());
  SORT(address(vec, 0), address(vec, vec.size()), CMPNum());
//   for (size_t i = 0; i < 100; ++i) {
//     cout << "v1[" << i << "] = " << *vec[i] << "\n";
//   }
}

int main(int argc, char const *argv[]) {

  // get number of threads
  size_t threads = 8;
  if (argc > 1) {
    threads = atoi(argv[1]);
  } else {
    threads = GetEnv("SPH_PARAM_THREADS", threads);
  }

  // get outfile name
  string ofname = GetEnvString("SPH_PARAM_RESULT_FILE", "test.h5");
  ofname += ".xml";

  // get problem size N
  size_t N = 100000000;
  N = GetEnv("SPH_PARAM_N", N);

  omp_set_num_threads(threads);

#ifdef _OPENMP
#pragma omp parallel 
#endif
  {
    if (omp_get_thread_num() == 0) {
      assert(omp_get_num_threads() == int(threads));
      cerr << "num threads: " << omp_get_num_threads() << "\n";
    }
  }

  cerr << "array length: " << N << "\n";

  vector<int> vec1(N);
  valarray<int> vec2(N);
  for (size_t i = 0; i < N; ++i) {
    vec1[i] = i;
    vec2[i] = i;
  }

  double t0 = omp_get_wtime();
  cerr << "teste vector\n";
  f(vec1);

  double t1 = omp_get_wtime();
  cerr << "teste valarray\n";
  g(vec2);

  double t2 = omp_get_wtime();

  ostringstream ostr;
  ostr << "<tests>\n";
  ostr << "<konfig threads='" << threads << "'/>\n";
  ostr << "<run>\n";
  ostr << "<end time='" << t2 - t0 << "'>\n";
  ostr << "<test type='vector' num-threads='" << threads << "' time='" << t1 - t0 << "'/>\n";
  ostr << "<test type='valarray' num-threads='" << threads << "' time='" << t2 - t1 << "'/>\n";
  ostr << "</end>\n";
  ostr << "</run>\n";
  ostr << "</tests>\n";

  ofstream ofile(ofname.c_str(), ios::binary);
  ofile << ostr.str();
  cout << ostr.str();
}
