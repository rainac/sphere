#ifndef jw_util_teestream_hh
#define jw_util_teestream_hh
/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

#include <iostream>

struct Teebuf : public std::streambuf {

  std::streambuf *sb1, *sb2;

  Teebuf(std::streambuf *sb1, std::streambuf *sb2) : 
    sb1(sb1), sb2(sb2) 
  {}

  int overflow(int c) {
    if (c == EOF) {
      return !EOF;
    }
    int const r1 = sb1->sputc(c);
    int const r2 = sb2->sputc(c);
    return r1 == EOF || r2 == EOF ? EOF : c;
  }

  int sync() {
    int const r1 = sb1->pubsync();
    int const r2 = sb2->pubsync();
    return r1 == 0 && r2 == 0 ? 0 : -1;
  }   

private:
  Teebuf(Teebuf const &) : std::streambuf(), sb1(), sb2() {}
  Teebuf &operator =(Teebuf const &) {
    return *this;
  }
};

struct Teestream : public std::ostream {
  Teebuf tbuf;

  Teestream(std::ostream & o1, std::ostream & o2) : 
    tbuf(o1.rdbuf(), o2.rdbuf())
  {
    rdbuf(&tbuf);
  }

};

#endif
