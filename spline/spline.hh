#ifndef spline_spline_a4lhgivn_hh
#define spline_spline_a4lhgivn_hh
/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

#include <errno.h>
#include <string.h>
#include <valarray>
#include <iostream>
#include <fstream>
#include "point.hh"

static double const dtau = 1;

template<class T>
inline std::istream &operator >>(std::istream &ein, std::valarray<T> &v) {
  size_t i = 0;
  for( ; ein and i < v.size(); ++i) {
    ein >> v[i];
  }
  assert(ein.fail() or i == v.size());
  return ein;
}

template<class T>
inline std::ostream &operator <<(std::ostream &aus, std::valarray<T> const &v) {
  for(size_t i = 0; aus and i < v.size(); ++i) {
//     if (i) aus << " ";
    aus << v[i] << "\n";
  }
  return aus;
}

template<int Order>
struct BaseFuncNBase {
  static double tau(int const &i) {
    return dtau * i;
  }
};

template<class T, int Order>
struct BaseFuncN : public BaseFuncNBase<Order> {
  typedef BaseFuncNBase<Order> Base;

  static size_t const subOrder = Order-1;
  BaseFuncN<T, Order-1> subBaseFunc;

  T operator()(size_t const i, T const &u) {
    T const r1 = subBaseFunc(i, u);
    T const r2 = subBaseFunc(i+1, u);
    return r1 * (u            - Base::tau(i)) / (Base::tau(i+subOrder) - Base::tau(i))
         + r2 * (Base::tau(i+Order) - u)      / (Base::tau(i+Order)    - Base::tau(i+1));
  }

};

template<class T>
struct BaseFuncN<T, 1> : public BaseFuncNBase<1> {

  T operator()(size_t const i, T const &u) {
    if (u >= tau(i) and u < tau(i+1)) {
      return 1;
    }
    return 0;
  }

};


template<class T, int Dim, int Order>
struct Spline {

  size_t const numPoints;
  static size_t const dim = Dim;     // dimension of space
  static size_t const order = Order; // order of splines

  BaseFuncN<double, order> baseFuncN;     // evaluate N
  
  std::valarray<Point<T, dim> > points;

  Spline(size_t numPoints) :
    numPoints(numPoints),
    points(numPoints)
  {}

  void setPoints(std::valarray<Point<T,dim> > const &v) {
    for (size_t i = 0; i < numPoints; ++i) {
      points[i] = v[i];
    }
//     for (size_t i = 0; i < order; ++i) {
//       points[numPoints - order + i] = v[numPoints - order - 1];
//       points[i] = v[order];
//     }
  }

  void readPoints(std::string const &pointsName) {
    std::ifstream einPunkte(pointsName.c_str(), std::ios::binary);
    std::valarray<Point<T, 2> > rawPoints(numPoints);
    einPunkte >> rawPoints;
    setPoints(rawPoints);
    if (einPunkte.fail()) {
      std::cerr << "error: reading points '" << pointsName
                << "' failed: " << strerror(errno) << "\n";
    }
  }

  Point<T, dim> operator()(double const &u) {
    Point<T, dim> res;
//     cerr << "evalat: " << u << "\n";
    size_t const num = numPoints;
    for (size_t i = 0; i < num; ++i) {
//       cerr << "bf: " << i << ", " << u << ", " << T(baseFuncN(i, u)) << "\n";
      res += points[i] * T(baseFuncN(i, u));
    }
    return res;
  }

  T tau(size_t i) const { return BaseFuncNBase<order>::tau(i); }

  void operator()(std::valarray<Point<T,dim> > &uv) {
    size_t const nSim = uv.size();
    double const begin = BaseFuncNBase<order>::tau(order - 1);
    double const end = BaseFuncNBase<order>::tau(numPoints);
    double const deltatau = (end - begin) / (nSim - 1);
//     cerr << "sample: " << begin << ", " << end << ", " << deltatau << "\n";

    for (size_t i = 0; i < nSim; ++i) {
      double const u = begin + i * deltatau;
      uv[i] = operator()(u);
    }
  }

  void faster(std::valarray<Point<T,dim> > &uv) {
    size_t const nSim = uv.size();
    double const begin = BaseFuncNBase<order>::tau(order - 1);
    double const end = BaseFuncNBase<order>::tau(numPoints);
    double const deltatau = (end - begin) / (nSim - 1);

    size_t ltaui = 0;
    double ltaunext = BaseFuncNBase<order>::tau(ltaui + order + 1);
    for (size_t i = 0; i < nSim; ++i) {
      double const u = begin + i * deltatau;
      if (u >= ltaunext) {
	++ltaui;
	ltaunext = BaseFuncNBase<order>::tau(ltaui + order + 1);
      }
      Point<T, dim> res;
      for (size_t ci = 0; ci < order + 1; ++ci) {
	res += points[ltaui + ci] * T(baseFuncN(ltaui + ci, u));
      }
//       cerr << (res - operator()(u)) << "\n";
//       assert(res == operator()(u));
      uv[i] = res;
    }
  }

};

#endif
