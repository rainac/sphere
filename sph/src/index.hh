#ifndef sph_index_hh
#define sph_index_hh
/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

#include <assert.h>
#include "power.hh"

/// \tparam C the coordinate integer type.
/// \tparam DIM the number of dimensions.
///
/// This provides a one-to-one mapping \f$ \alpha \f$ 
/// from \f$D\f$-dimensional spatial raster 
/// coordinates \f$ {\bf I}=(I_0, I_1, ..., I_{D - 1}) \f$
/// to a single integer array offset \f$ i := \alpha(I) \f$. Since hashing has been
/// introduced this class is only used by class HashedIndex.
///
/// The mapping is the classical many-dims-to-one mapping, i.e.
/// \f[ \alpha(I) = i := \sum_{j=0}^{D} I_j \prod_{k=0}^{j-1} d_k, \f] where 
///   \f$ {\bf d} =  (d_0, d_1, ..., d_{D-1})\f$ are the number of items in each dimension 
///   and \f$\prod_{k=0}^{-1} := 1\f$. 
/// It is also required that \f$ 0 \le I_j < d_j \f$ for \f$ 0 \le j < D \f$.
///
/// The factors \f$ \prod_{k=0}^{j-1} d_k\f$ are precomputed for \f$ 0 \le j < D \f$ and
/// stored in the vector \f${\bf w} = (w_0, w_1, w_2, ..., w_{D-1}) = (1, d_0, d_0 d_1, ..., \frac{\prod {\bf d}}{d_{D-1}})\f$. Then the index \f$ i\f$ can be computed simply as 
/// the scalar product \f$i = {\bf w}\cdot {\bf I} \f$.
///
/// The assignment of indices to spatial raster bins is demonstrated in the following image,
/// showing the raster of the computing domain
/// \f$[0,1] \times [0,0.7] \;{\rm m}^2\f$ with \f$H=0.1 \;{\rm m}\f$. This results in
/// \f$ {\bf d} = (10, 7)\f$ and \f$ {\bf w} = (1, 10) \f$.
/// \image html raster-index.png "The spatial raster in 2D of the computing domain [0,1] x [0,0.7]" width=1cm
///
template<class C, int Dim>
struct Index {

  /// Contains the number of bins in each dimension, \f${\bf d}\f$.
  Point<C, Dim> base;
  /// Contains the weights vector \f${\bf w}\f$.
  Point<C, Dim> weights;
  /// The size of the index, aka the number of possible index
  /// values. This should be the same as \f$ \prod_{k=0}^{D-1} d_k \f$, 
  /// the product over all dimensions.
  size_t size;

  /// The constructor computes the weights vector entries and size.
  Index(Point<C, Dim> const &_base) : 
    base(_base) 
  {
    weights[0] = 1;
    for(size_t i = 1; i < Dim; ++i) {
      weights[i] = weights[i-1] * base[i-1];
    }
    size = weights[Dim-1] * base[Dim-1];
  }

  /// This computes the mapping \f$ \alpha(I) \f$.
  /// \param p an integer raster coordninate \f$ I \f$.
  /// \returns the mapping \f$ i = \alpha(I) \f$.
  /// \tparam U the integer type of the argument.
  template<class U>
  C operator()(Point<U, Dim> const &p) const {
    return (weights * p).sum();
  }

  /// This computes the inverse mapping \f$ {\bf I} = \alpha^{-1}(i) \f$.
  /// The return value will contain \f$ D \f$ non-negative integers 
  /// \f$ {\bf \tilde I} = (\tilde I_0, \tilde I_1, ..., \tilde I_{D-1}) \f$ 
  /// such that \f$ \alpha(\tilde {\bf I}) = i\f$.
  /// \param p a integer less than size.
  /// \returns the inverse mapping \f$ I = \alpha(i) \f$.
  Point<C, Dim> invert(C p) const {
     Point<C, Dim> res;
     for(long i = Dim - 1; i >= 0; --i) {
       res[i] = p / weights[i];
       p -= res[i] * weights[i];
       assert(p <= weights[i]);
     }
     return res;
   }
};


/// The hash function used simply keeps the \f$ m \f$ lowest bits and
/// ignores the others. Thus the hash function is given 
/// by \f$ \#(x) = x \& (2^m - 1)\f$, where \f$ \& \f$ is the
/// bitwise \b and operation. \f$ 2^m - 1 \f$ is the integer where the lowest
/// \f$ m \f$ bits are set and the others are zero. \f$ 0 < m < 32 \f$ must hold.
/// In effect, the hashed value is the same as \f$ \#(x) = x \;{\rm mod}\; 2^m  \f$
struct HashFunktion {

  /// This holds the value \f$ m \f$, just a copy of SPHKonfig::m.
  unsigned const m_shift;
  /// This holds the number \f$ 2^m \f$, the number of possible hash values.
  unsigned long const m_dim;
  /// This holds the number \f$ 2^m - 1 \f$, the mask used in the hashing operation.
  unsigned long const m_mask;

  /// The constructor receives the value of SPHKonfig::m and computes
  /// the other constants from it.
  /// \param m_shift the value of SPHKonfig::m
  HashFunktion(unsigned const m_shift) : 
    m_shift(m_shift),
    m_dim(1 << m_shift),
    m_mask(m_dim - 1)
  {}
  
  /// Computes the hashed value \f$ \#(v) = v \& (2^m - 1)\f$ of parameter v.
  /// \param v the value to be hashed.
  /// \returns the hashed value.
  long operator()(long const v) const { return v & m_mask; }
};



/// This class provides the one-to-one mapping \f$ \beta \f$ 
/// from \f$D\f$-dimensional spatial raster 
/// coordinates \f$ {\bf I}=(I_0, I_1, ..., I_{D - 1}) \f$
/// to a single integer array offset \f$ i := \beta(I) \f$. 
///
/// The mapping first hashes each entry in \f$ {\bf I} \f$ using class HashFunktion.
/// Also, the vector \f$ {\bf d} \f$ is replaced by the vector  \f$ (2^m, 2^m, ..., 2^m) \f$, 
/// since  \f$ 2^m \f$ is the number of possible hashed values.
/// 
/// Then the mapping procedes as described in class Index, but the special circumstances
/// allow for a faster computation, as follows:
/// \f[ \beta(I) = i = \sum_{j=0}^{D} [I_j << (jm) ], \f] where  
///  \f$ << \f$ is the bitwise shift left operation and \f$ \sum \f$
///  is the bitwise or operation which, in this case is equivalent to
///  the arithmetic sum.
///
/// The shift amounts \f$ jm \f$ are precomputed and
/// stored in the vector \f${\bf S} = (s_0, s_1, ..., s_{D-1}) = (0, m, 2m, ...)\f$. Then the
/// index can be simply computed as \f$i = \sum_{j=0}^{D-1} ({\bf I} << {\bf S})_j \f$, where
///  \f$ << \f$ is the elementwise shift left operation as implemented by std::valarray.
///
/// The inverse mapping \f$ \beta^{-1} \f$ can be implemented using
/// only bit-operations in a similar manner.
///
/// The assignment of indices to spatial raster bins is demonstrated in the following image,
/// showing the raster of the computing domain
/// \f$[0,1] \times [0,0.7] \;{\rm m}^2\f$ with \f$H=0.1 \;{\rm m}\f$ and \f$ m=2 \f$. 
/// This results in \f$ {\bf S} = (2, 2)\f$.
/// \image html raster-index-hashed.png "The spatial raster in 2D of the computing domain [0,1] x [0,0.7] with hashed indices" width=1cm 
///
/// \tparam C the coordinate integer type.
/// \tparam DIM the number of dimensions.
///
template<class C, int Dim>
struct HashedIndex {

  /// A shorthand for the non-hashed index-class.
  typedef Index<C, Dim> Base;

  /// The object providing the hash function.
  HashFunktion hashFunktion;

  /// Contains the numbers \f$ jm \f$ for \f$ 0 \le j < D \f$. These
  /// are the number of bits to shift each entry of the spatial
  /// integer index by.
  Point<C, Dim> shifts;

  /// The object providing the non-hashed index.
  Base baseIndex;

  /// The number of possible hash values, \f$ (2^m)^D  = 2^{mD} \f$.
  size_t const size;

  /// The constructor computes the shifts vector entries and size.
  HashedIndex(Point<C, Dim> const &_base, long const m) : 
    hashFunktion(m),
    baseIndex(_base),
    size(dynamic_power(hashFunktion.m_dim, Dim))
  {
    for(size_t i = 0; i < Dim; ++i) {
      shifts[i] = hashFunktion.m_shift * i;
//       std::cerr << "hash shifts: w[" << i << "] = " << shifts[i] << "\n";
    }
//     std::cerr << "hash size: m = " << m << " 2^m^" << ndim
// 	      << " = " << size << "\n";
  }

  /// This computes the mapping \f$ \beta(I) \f$.
  /// \param p an integer raster coordninate \f$ I \f$.
  /// \returns the mapping \f$ i = \beta(I) \f$.
  /// \tparam U the integer type of the argument.
  template<class U>
  C operator()(Point<U, Dim> const &p) const {
    Point<U, Dim> ph(p);
    ph &= hashFunktion.m_mask;
    ph <<= shifts;
    return ph.binOr();
  }

  /// This computes the inverse mapping \f$ {\bf I} = \beta^{-1}(i) \f$.
  /// The return value will contain \f$ D \f$ non-negative integers 
  /// \f$ {\bf \tilde I} = (\tilde I_0, \tilde I_1, ..., \tilde I_{D-1}) \f$ 
  /// such that \f$ \beta(\tilde {\bf I}) = i\f$.
  /// \param p a integer less than size.
  /// \returns the inverse mapping \f$ I = \beta(i) \f$.
  Point<C, Dim> invert(C const p) const {
    Point<C, Dim> res;
    for(size_t i = 0; i < Dim; ++i) {
      res[i] =  p >> shifts[i];
      res[i] &= hashFunktion.m_mask;
    }
    return res;
  }

};

#endif
