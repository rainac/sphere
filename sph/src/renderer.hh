#ifndef sph_renderer_hh
#define sph_renderer_hh
/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

/// class GLStringStream provides a method to write a string in the GL
/// scene a certain position

enum DisplayListIds {
  none,
  CoordAxesDisplayListId,
  dlIdSolidUnitSphere,
  dlIdWireUnitSphere,
  dlIdSensors
};

struct GLString {
  
  static void renderString(void *font, std::string const &text, 
			   double scale, Point3D const &_pos) {
    Point3D pos(_pos);
    
    scale *= 1;

//     glMatrixMode(GL_MODELVIEW_MATRIX);
//     glPushMatrix();
//     glMatrixMode(GL_PROJECTION_MATRIX);
    glPushMatrix();

    glTranslated(pos[0], pos[1], pos[2]);
    glScaled(scale, scale, scale);

    for(size_t i = 0; i < text.size(); ++i) {
      char const c = text[i];
      switch(c) {
      case '\n':
	glPopMatrix();
	glPushMatrix();
	pos[1] += -125*scale;
	// 	pos += -125*scale*up;
	glTranslated(pos[0], pos[1], pos[2]);
	glScaled(scale, scale, scale);
	break;
      case '\r':
	glPopMatrix();
	glPushMatrix();
	glTranslated(pos[0], pos[1], pos[2]);
	glScaled(scale, scale, scale);
	break;
#if defined SPH_RENDER_TEXT_TABS_AS_SPACES || 1
      case '\t':
	for(size_t ci = 0; ci < 4; ++ci) {
	  glutStrokeCharacter(font, ' ');
	}
	break;
#endif
      case '\a':
	std::cout << c; // beep on cout
	break;
      case '\v':
	glTranslated(0, -125*scale, 0);
	//  	glTranslated(-up[0]*125*scale, -up[1]*125*scale, -up[2]*125*scale);
	glScaled(scale, scale, scale);
	break;
      default:
	glutStrokeCharacter(font, c);
	break;
      }
    }

//     glMatrixMode(GL_MODELVIEW_MATRIX);
//     glPopMatrix();
//     glMatrixMode(GL_PROJECTION_MATRIX);
    glPopMatrix();
  }

};

struct GLMatrix {
  double m_v[16];
  
  GLMatrix(GLint matrixType = GL_PROJECTION_MATRIX) {
//     double buffer[16];
    glGetDoublev(matrixType, m_v);
    int const ec = glGetError();
    if (ec != GL_NO_ERROR) {
      cerr << "failed to download matrix\n";
    }
//     for(int i = 0; i < 4; ++i) {
//       for(int j = 0; j < 4; ++j) {
//         m_v[i][j] = buffer[j*4 + i];
//       }
//     }
  }

  double  operator[](size_t i) { return m_v[i]; }
  double  operator()(size_t i) { return m_v[i]; }
  double  operator()(size_t i, size_t j) { return m_v[i + 4 * j]; }

  void print(std::ostream &aus) const {
    for(int i = 0; i < 4; ++i) {
      for(int j = 0; j < 4; ++j) {
        aus << m_v[i + 4*j] << " ";
      }
      aus << "\n";
    }
  }

};

std::ostream &operator <<(std::ostream &aus, GLMatrix const &m) {
  m.print(aus);
  return aus;
}

#endif
