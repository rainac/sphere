#ifndef jw_utility_uname
#define jw_utility_uname
/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

#include <iostream>
#include <string>
#include <sys/utsname.h>

struct Uname {
  bool m_fail;
  std::string sysname, nodename, release, version, machine;
  std::string domainname;
  Uname() : m_fail() {
    struct utsname buffer;
    m_fail = uname(&buffer);
    sysname = buffer.sysname;
    nodename = buffer.nodename;
    release = buffer.release;
    version = buffer.version;
    machine = buffer.machine;
#ifdef _GNU_SOURCE
    domainname = buffer.domainname;
#endif
  }
  bool fail() const { return m_fail; }
  operator bool () const { return !fail(); }
  void writeXML(std::ostream &aus, std::string const &indent = "") const {
    aus << indent;
    print(aus);
  }
  void print(std::ostream &aus) const {
    aus << "<uname"
      " sysname='" << sysname << "'"
      " nodename='" << nodename << "'"
      " release='" << release << "'"
      " version='" << version << "'"
      " machine='" << machine << "'"
      " domainname='" << 
#ifdef _GNU_SOURCE
      domainname
#else
      "&lt;unknown, not compiled with _GNU_SOURCE&gt;"
#endif
	<< "'"
      "/>";
  }
  std::string all() const {
    return sysname + " " + 
      nodename + " " + 
      release + " " + 
      version + " " + 
      machine + " " + 
      domainname;
  }
};
inline std::ostream &operator <<(std::ostream &aus, Uname const &un) {
  un.print(aus);
  return aus;
}

#endif
