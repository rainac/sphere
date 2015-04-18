#ifndef jw_ut_locale_hh
#define jw_ut_locale_hh
/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

#include <locale.h>

struct GlobalLocale {
  std::string const fullstring;

  GlobalLocale() : 
    fullstring()
  {
    init();
  }

  GlobalLocale(std::string const &_fullstring) : 
    fullstring(_fullstring)
  {
    init();
  }
  
  void init() {
    try {
      std::locale loc(fullstring.c_str());
      std::locale::global(loc);
    } catch (std::exception &e) {
      std::cerr << "error: locale construction locale(\"" << fullstring
		<< "\") failed: " << e.what() << endl;
      std::cerr << "error: check your environment settings (LANG, LC_ALL, etc.)" << endl;
      std::cerr << "error: locale is set to: " << setlocale(LC_ALL, 0) << "\n";
    }
  }

  operator std::string() const {
    return name(); 
  }

  std::string name() const {
    return fullstring;
  }

  void print(std::ostream &aus) const {
    aus << "global locale setter: " << name() << "\n";
  }

};

inline std::ostream &operator <<(std::ostream &aus, GlobalLocale const &loc) {
  loc.print(aus);
  return aus;
}

#endif
