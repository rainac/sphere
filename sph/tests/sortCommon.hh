/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/
#include <valarray>

template<class T>
inline std::ostream &operator <<(std::ostream &aus, std::valarray<T> &v) {
  for(size_t i = 0; i < v.size(); ++i) {
    aus << v[i] << " ";
  }
  return aus;
}

struct Data {
  Data() : first(), second() {}
  long first, second;
  
  bool operator< (Data const &v) const {
    return first < v.first;
  }
  bool operator== (Data const &v) const {
    return first == v.first and second == v.second;
  }
};

inline std::ostream &operator <<(std::ostream &aus, Data &v) {
  aus << "(" << v.first << ", " << v.second << ")";
  return aus;
}

inline bool isStable(std::valarray<Data> const &data) {
  for(size_t i = 1; i < data.size(); ++i) {
    CPPUNIT_ASSERT(data[i - 1].first <= data[i].first);
    if (data[i - 1].first == data[i].first) {
      if (data[i - 1].second >= data[i].second) {
        return false;
      }
    }
  }
  return true;
}

template<class T>
inline bool isSorted(std::valarray<T> const &data) {
  for(size_t i = 1; i < data.size(); ++i) {
    if (data[i] < data[i - 1]) {
      return false;
    }
  }
  return true;
}

