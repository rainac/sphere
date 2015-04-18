#ifndef jw_indir_sort_hh
#define jw_indir_sort_hh
/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

// #include <parallel/algorithm>
#include <algorithm>
#include <functional>
#include <vector>
#include "permutation.hh"

/// IndirectSort
struct IndirectSort {
  
  /// constructor
  IndirectSort() {}

  /// This template class constructs a binary predicate out of another
  /// binary predicate inserting an indirection. The constructed
  /// predicate provides operator()(unsigned const a, unsigned const
  /// b).  The wrapped predicate is called with arguments *(begin +
  /// a), *(begin + b), where begin is a RandIter that is given
  /// to the constructor WrappedPredicate.
  
  /// This predicate can be used to sort an integer permutation p such
  /// that v[p[i]] is sorted according to the wrapped predicate, and v
  /// is unchanged.

  /// \tparam RandIter the type of the two arguments to compare
  /// \tparam Pred type of the predicate that is being wrapped

  template<class RandIter, class Pred>
  struct WrappedPredicate {

    /// Reference to the begin iterator (base pointer)
    RandIter &begin;
    /// A copy of the wrapped predicate.
    Pred pred;

    /// Constructor.
    /// \param begin Reference to the begin iterator (base pointer)
    /// \param pred A predicate Pred to wrap.
    WrappedPredicate(RandIter &begin, Pred pred) : 
      begin(begin),
      pred(pred) 
    {}
    
    /// The comparison operator that this class provides.
    /// \param a index of first object
    /// \param b index of second object
    bool operator()(unsigned const a, unsigned const b) const {
//       std::cerr << "cmp " << a << " " <<  b << "\n";
//       std::cerr << "cmp " << *( begin + a ) << " " <<  *( begin + b ) << "\n";
      return pred(*( begin + a ), *( begin + b ));
    }

  };

  template<class RandIter, class P > 
  void operator ()(RandIter from, RandIter to, std::valarray<unsigned> &perm, P pred = P()) {

    assert(to - from == long(perm.size()));

    // sort the permutation such that permuted vector is sorted
    // this will swap around only integers
    std::sort(&perm[0], &perm[to - from], WrappedPredicate<RandIter, P>(from, pred));

  }

};

template<class V, class P >
inline void sortPermutation(V from, V to, std::valarray<unsigned> &perm, P pred = P()) {
  IndirectSort is;
  is.operator()<V,P>(from, to, perm, pred);
}

#endif
