#ifndef sph_quader_hh
#define sph_quader_hh
/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

#ifdef SPH_RENDER
#include <GL/freeglut.h>
# include <GL/gl.h>
#endif

#include "vector.hh"

/// Struct Quader simply holds two points, untenLinks and obenRechts
/// and also provides the difference vector v.
template<class Vec>
struct Quader {
  typedef Vec vector_type;
  typedef Vector::value_type value_type;

  vector_type untenLinks, obenRechts, v;

  Quader() {}

  Quader(vector_type const &_untenLinks, vector_type const &_obenRechts) :
    untenLinks(_untenLinks), 
    obenRechts(_obenRechts),
    v(untenLinks.size())
  {
    v = obenRechts - untenLinks;
  }
  
  bool isOut(Vector const &p) const {
    return (p < untenLinks).max() == 1 // one coordinate smaller than P1
      or (p > obenRechts).max() == 1;  // one coordinate larger than P2
  }

  value_type laenge(size_t const dim) const { return v[dim]; }
  value_type offset(size_t const dim) const { return untenLinks[dim]; }

  void writeXML(std::ostream &aus, std::string const &indent = "") const {
    aus << indent << "<gebiet>\n"
        << indent << "\t<ul x='" << untenLinks[0]
	<< "' y='" << (ndim > 1 ? untenLinks[1] : 0)
	<< "' z='" << (ndim > 2 ? untenLinks[2] : 0) << "'/>\n"
        << indent << "\t<or x='" << obenRechts[0]
	<< "' y='" << (ndim > 1 ? obenRechts[1] : 0)
	<< "' z='" << (ndim > 2 ? obenRechts[2] : 0) << "'/>\n";
    aus << indent << "</gebiet>\n";
  }
  
#ifdef SPH_RENDER
  void render() const {

    Point<double,3> gul(untenLinks.data, ndim);
    Point<double,3> gv(v.data, ndim);

    glPushMatrix();

#if SPH_DIM != 1
    glBegin(GL_LINE_STRIP);
    {
      glVertex3dv(&gul[0]);

      Point<double,3> d;

      d[0] = gv[0];
      glVertex3dv(&(gul + d)[0]);

      d = 0;
      d[0] = gv[0];
      d[1] = gv[1];
      glVertex3dv(&(gul + d)[0]);

      d = 0;
      d[1] = gv[1];
      glVertex3dv(&(gul + d)[0]);

      glVertex3dv(&gul[0]);

      d = 0;
      d[2] = gv[2];

      glVertex3dv(&(gul + d)[0]);

      d[0] = gv[0];
      glVertex3dv(&(gul + d)[0]);

      d = 0;
      d[2] = gv[2];
      d[0] = gv[0];
      d[1] = gv[1];
      glVertex3dv(&(gul + d)[0]);

      d = 0;
      d[2] = gv[2];
      d[1] = gv[1];
      glVertex3dv(&(gul + d)[0]);

      d = 0;
      d[2] = gv[2];
      glVertex3dv(&(gul + d)[0]);
    }
    glEnd();

    // drei fehlen noch
    glBegin(GL_LINES);
    {
      Point<double, 3> d;
      d = 0;
      d[1] = gv[1];
      glVertex3dv(&(gul + d)[0]);
      d[2] = gv[2];
      glVertex3dv(&(gul + d)[0]);

      d = 0;
      d[0] = gv[0];
      d[1] = gv[1];
      glVertex3dv(&(gul + d)[0]);
      d[2] = gv[2];
      glVertex3dv(&(gul + d)[0]);

      d = 0;
      d[0] = gv[0];
      glVertex3dv(&(gul + d)[0]);
      d[2] = gv[2];
      glVertex3dv(&(gul + d)[0]);
    }
    glEnd();
#else
    // DIM == 1

    glBegin(GL_POINTS);
    {
      glVertex3dv(&gul[0]);
      glVertex3dv(&(gul + gv)[0]);
    }
    glEnd();

#endif
    glPopMatrix();
  }
#endif
};

#endif
