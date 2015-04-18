#ifndef jw_getpointee_hh
#define jw_getpointee_hh
/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

/// Helper template GetPointee provides the member type value_type
/// which is the type an iterator of type T points to. The standard
/// version uses T::value_type as the desired
/// type. Specialization for pointers is below.

template<class T>
struct GetPointee {
  typedef typename T::value_type value_type;
};


/// Specialization for pointers of template class GetPointee for types
/// of form T* sets value_type to T.

template<class T>
struct GetPointee<T*> {
  typedef T value_type;
};

#endif
