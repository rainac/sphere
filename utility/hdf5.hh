/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

struct HDF5 {

  HDF5() {
    open();
  }


  ~HDF5() {
    close();
  }

  void open() {
//     cerr << "opening HDF5 library\n";
    H5open();
    H5check_version(H5_VERS_MAJOR, H5_VERS_MINOR, H5_VERS_RELEASE);
  }

  void close() {
//     cerr << "closing HDF5 library\n";
    H5garbage_collect();
    H5close();
  }

  void writeXML(std::ostream &aus) const {
    aus << "<hdf5"
      " major='" << H5_VERS_MAJOR << "'"
      " minor='" << H5_VERS_MINOR << "'"
      " release='" << H5_VERS_RELEASE << "'"
      "/>\n";
  }
  
};

inline std::ostream &operator <<(std::ostream &aus, HDF5 const &hdf5) {
  hdf5.writeXML(aus);
  return aus;
}
