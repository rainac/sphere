#ifndef jw_sph_util23411_hh
#define jw_sph_util23411_hh
/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

/// Helper function that resizes an array if the current size() differs from newsz.
/// Calls resize(newsz) if v.size() != newsz. 
//
/// \param v the array object (e.g. an std::valarray<of-something>)
/// \param newsz the desired new size
/// \tparam T1 the array type to resize
template<class T1>
void resizeIfNeeded(T1 &v, size_t const newsz) {
  if (v.size() != newsz) {
    v.resize(newsz);
  }
}

#endif
