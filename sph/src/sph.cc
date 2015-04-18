#ifndef jw_sph_cc 
#define jw_sph_cc // this file is included by others...
/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

#include <iostream>
#include <iomanip>
#include <set>
#include <fstream>
#include <cmath>
#ifdef _OPENMP
# include <omp.h>
#endif
#include <errno.h>
#include <signal.h>
#include <fenv.h>
#include <sys/stat.h>

// #define DEBUG
extern "C" char **environ;
#include "sph.hh"
#include "utility/uname.hh"
#include "utility/pid.hh"
// #include "utility/pmap.hh"
#include "utility/strftime.hh"
#include "utility/hdf5.hh"

// #define DEBUG_RK
#include "integrate/integrate.hh"

using namespace std;

#ifdef SPH_RENDER
/// Helper constant used by the renderer.
/// This constant is only present if compiled with macro SPH_RENDER defined.
/// \todo make sure to use same default values here and in SPKConfig constructor.
double const SPH_config_h = GetEnv("SPH_PARAM_H", 1e-2, 1e-100, 1);


/// Helper constant used by the renderer.
/// This constant is only present if compiled with macro SPH_RENDER defined.
/// \todo make sure to use same default values here and in SPKConfig constructor.
long const SPH_config_m = GetEnv("SPH_PARAM_M", 4, 1, 12);

/// Helper constant used by the renderer.
/// This constant is only present if compiled with macro SPH_RENDER defined.
long const SPH_config_hashmask = (1 << SPH_config_m) - 1;

#endif

#ifdef SPH_AD_DEACTIVATE_KERNEL
typedef KernelInterface<double, DoubleVector, ndim> SPHKernel;
#else
typedef KernelInterface<SphADataType, Vector, ndim> SPHKernel;
#endif
typedef PartikelManager<SPHConfig<SphADataType>, SPHKernel> SPHPartikelManager;
typedef SPHPartikelManager::ForceFunctorWrapper ForceFunctor;
typedef Integrator<ForceFunctor>                  SPHIntegrator;
typedef SPH<SPHPartikelManager, SPHIntegrator> SPHT;


static bool beende = false;

#ifdef SPH_RENDER
struct Quaternion {
  double m_scalar;
  Point3D m_vector;
  Quaternion() : m_scalar(), m_vector() {}
  explicit Quaternion(double const &s) : m_scalar(s), m_vector() {}
  explicit Quaternion(Point3D const &v) : m_scalar(), m_vector(v) {}
  Quaternion(double const &s, Point3D const &v) : m_scalar(s), m_vector(v) {}
  
  Quaternion &operator =(double const &o) {
    m_scalar = o;
    m_vector = 0;
    return *this;
  }

  Quaternion &operator =(Quaternion const &o) {
    m_scalar = o.scalar();
    m_vector = o.vector();
    return *this;
  }

  Quaternion &operator *=(double const &v) {
    m_scalar *= v;
    m_vector *= v;
    return *this;
  }

  Point3D rotate(Point3D const &v) const {
//     Quaternion res(0, v);
//     res *= v;
//     res *= invert();
    return (invert() * Quaternion(0, v) * *this).vector();
  }

  Quaternion invert() const {
    Quaternion res(scalar(), -vector());
    res *= 1 / normSquared();
    return res;
  }

  Quaternion operator *(Quaternion const &v) const {
    Quaternion res(*this);
    res *= v;
    return res;
  }

  Quaternion &operator *=(Quaternion const &o) {
    double const nsc = scalar()*o.scalar() - vector().scalar(o.vector());
    m_vector = scalar()*o.vector() + o.scalar()*vector() 
      + cross(vector(), o.vector());
    m_scalar = nsc;
    return *this;
  }
  
  Quaternion &operator *=(Point3D const &o) {
    double const nsc =  -vector().scalar(o);
    m_vector = scalar()*o + cross(vector(), o);
    m_scalar = nsc;
    return *this;
  }
  
  double       angle()      const { 
    double const a = acos(scalar())*2;
//     return (vector().min() < 0 ? 2*M_PI - a : a);
    return a;
  }

  double       normSquared()       const { 
    return scalar()*scalar() + vector().normSquared();
  }

  double       norm()       const { 
    return sqrt(normSquared());
  }

  double       &scalar()       { return m_scalar; }
  double const &scalar() const { return m_scalar; }

  Point3D       &vector()       { return m_vector; }
  Point3D const &vector() const { return m_vector; }

};

struct Rotation : public Quaternion {
  Rotation(double const &s, Point3D const &v) : 
    Quaternion(cos(s/2), sin(s/2)*v) 
  { }
};

static double deg2rad(double const x) { return x / 360 * 2 * M_PI; }
static double rad2deg(double const x) { return x / (2 * M_PI) * 360; }
/*
static double angle(Point3D const &p, Point3D const &o) { 
  return acos(scalar(p,o) / (norm(p)*norm(o)));
}
static double csysAngle(Point3D const &p, Point3D const &o) { 
  return acos(scalar(p,o));
}
*/
struct CoordSystem {
  Point3D axes[3];
  CoordSystem() {
    axes[0][0] = 1;
    axes[1][1] = 1;
    axes[2][2] = 1;
  }
  void renormalize() {
    axes[0] /= norm(axes[0]);
    axes[1] /= norm(axes[1]);
    axes[2] /= norm(axes[2]);
  }
  void rotate(Quaternion const &q) {
    axes[0] = q.rotate(axes[0]);
    axes[1] = q.rotate(axes[1]);
    axes[2] = q.rotate(axes[2]);
  }
  void print(std::ostream &aus) const {
    aus
      << "x: " << axes[0] << "\n"
      << "y: " << axes[1] << "\n"
      << "z: " << axes[2] << "\n"
      ;
  }
  void render() const {
    glBegin(GL_LINES);
    {
      glColor3f(1, 0, 0);
      glVertex3f(0, 0, 0);
      glVertex3f(1, 0, 0);

      glColor3f(0, 1, 0);
      glVertex3f(0, 0, 0);
      glVertex3f(0, 1, 0);

      glColor3f(0, 0, 1);
      glVertex3f(0, 0, 0);
      glVertex3f(0, 0, 1);
    }
    glEnd();
  }
};
inline std::ostream &operator <<(std::ostream &aus, CoordSystem const &p) {
  p.print(aus);
  return aus;
}

struct Camera {

  static const CoordSystem coordSystem;

private:
  bool m_updateView, m_updateScreen;
  Point3D m_position;
  Quaternion m_rotation;
  double fovAngle;

  static double const fovAngleMin, fovAngleMax;

public:
  Camera() :
    m_updateView(), 
    m_updateScreen(),
    fovAngle()
  {}

  void rotate(double const deg, int const axis) {
    Rotation rot(deg2rad(deg), coordSystem.axes[axis]);
    m_rotation = rot*m_rotation;
    m_updateScreen = 1;
  }

  void translate(double const amount, int const axis) {
    Point3D const rotatedAxis = m_rotation.rotate(coordSystem.axes[axis]);
//     cerr << "c.translate: rotated axis " << axis << ": " << rotatedAxis;
    m_position += amount*rotatedAxis;
    m_updateScreen = 1;
  }

  void position(Point3D const &p) { 
    m_position = p;
  }

  void zoomIn(double amDegree) {
    fovAngle -= amDegree;
    fovAngle = max(fovAngle, fovAngleMin);
    m_updateView = 1;
  }

  void zoomOut(double amDegree) {
    fovAngle += amDegree;
    fovAngle = min(fovAngle, fovAngleMax);
    m_updateView = 1;
  }

  void reset() {
    m_rotation = 1;
    m_position = 0;
    fovAngle = 30;
  }

  void clear() { 
    m_updateView = m_updateScreen = false; 
  }

  bool updateView() const { return m_updateView; }
  bool updateScreen() const { return m_updateScreen; }
  double fov() const { return fovAngle; }
  Point3D const &position() const { return m_position; }
  Quaternion const &rotation() const { return m_rotation; }

};

struct Mouse {
  bool inside;
  bool buttons[5];
  int lastClickX, lastClickY;
  Mouse() :
    inside(),
    lastClickX(),
    lastClickY()
  {
    buttons[0] =
      buttons[1] =
      buttons[2] =
      buttons[3] =
      buttons[4] = false;
  }
};

/// The Renderer class gathers the OpenGL rendering code and (static) variables.
struct Renderer {

/// \todo Partikel als Volumen (Flaechen) rendern 
/// \todo Eigenschaften von Partikeln graphisch anzeigen
/// \todo Mehr interaktivitaet beim rendern
/// \todo Steuerung "richtig" machen: Kamera-Objekt, Quaternionen, Maus, Scollrad
/// \todo Kann man ohne GLUT/FReeGLUT auskommen? (Freeglut auf der Sun?)
/// \todo Lineale/Raster anzeigen?
/// \todo Interaktivitaet: Partikel auswaehlen/ausblenden/Infos anzeigen

  static Camera camera;
  static Mouse mouse;

  static long windowId;
  static size_t frameRendered, screenWidth, screenHeight;
  static SPHT *sphPtr;

  static RenderOptions options;

  static Point3D lcsysPos, textLineP1, textLineP2, initCamPos;
  static CoordSystem const coordSystem;
//   static Quaternion axes[3];

  static Quaternion initCamRotation;

  static Mutex rendererReady;
  static Condition rendererInitialized;

  static double const rUnit; // degree
  static double const mUnit; // meter

  static void updateFixedScreenPos() {

    GLMatrix modelviewMatrix(GL_MODELVIEW_MATRIX);
    GLMatrix projectionMatrix(GL_PROJECTION_MATRIX);
    
    GLint viewport[4] = {0};
    glGetIntegerv(GL_VIEWPORT, viewport);
    
    gluUnProject(20, 20, 0.3,
                 modelviewMatrix.m_v, projectionMatrix.m_v,
                 viewport,
                 &lcsysPos[0], &lcsysPos[1], &lcsysPos[2]);

    gluUnProject(40, 20, 0.3,
                 modelviewMatrix.m_v, projectionMatrix.m_v,
                 viewport,
                 &textLineP1[0], &textLineP1[1], &textLineP1[2]);
    
    gluUnProject(viewport[2] - 20, 20, 0.3, 
                 modelviewMatrix.m_v, projectionMatrix.m_v,
                 viewport,
                 &textLineP2[0], &textLineP2[1], &textLineP2[2]);
    
  }

  static void setupView() {
    Point3D camPos, camTarget, camUpVector;
    camTarget[2] = -1;
//     camPos[2] = 1;
    camUpVector[1] = 1;

    glViewport(0, 0, screenWidth, screenHeight);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(camera.fov(), 
                   double(screenWidth) / screenHeight, 1e-2, 100);
    
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(camPos[0], camPos[1], camPos[2],
	      camTarget[0], camTarget[1], camTarget[2], 
	      camUpVector[0], camUpVector[1], camUpVector[2]);

  }

  static void setupViewOverlay() {
    Point3D camPos, camTarget, camUpVector;
    camTarget[2] = -1;
//     camPos[2] = 1;
    camUpVector[1] = 1;

    glViewport(0, 0, screenWidth, screenHeight);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(30, 
                   double(screenWidth) / screenHeight, 1e-1, 100);
    
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(camPos[0], camPos[1], camPos[2],
	      camTarget[0], camTarget[1], camTarget[2], 
	      camUpVector[0], camUpVector[1], camUpVector[2]);

  }

  static void setupOpenGL() {
    setupViewOverlay();
    updateFixedScreenPos();
//     setupView();
  }


  static void idleFunction() {
//     cerr << "render: idling " << frameRendered << "!\n";
  }
  
  static void visibleFunction(int state) {
    switch(state) {
    case GLUT_VISIBLE:
      becomeVisible();
      break;
    case GLUT_NOT_VISIBLE:
      becomeInvisible();
      break;
    default:
      assert(0);
      break;
    }
  }
  static void becomeInvisible() {
    cerr << "render: event: becoming invisible\n";
    glutDisplayFunc(displayFunctionNull);
    options.visible = false;
  }
  static void becomeVisible() {
    cerr << "render: event: becoming visible\n";
    if (options.visible == false) {
      glutDisplayFunc(displayFunction);
    }
    options.visible = true;
  }

  static void displayFunctionNull() { }

  static void displayOverlayFunction() {

    setupViewOverlay();

    glColor3f(0, 0, 1);
    glBegin(GL_LINES);
    glVertex3dv(&textLineP1[0]);
    glVertex3dv(&textLineP2[0]);
    glEnd();
    
    glPushMatrix();
    {
      glTranslated(lcsysPos[0], lcsysPos[1], lcsysPos[2]);
      //         glutSolidSphere(0.001, 16, 8);
      rotateScene();
      glScalef(0.002, 0.002, 0.002);
      coordSystem.render();
    }
    glPopMatrix();
    
    glPushMatrix(); 
    {
      if (sphPtr->renderInfo.numRendered == 0) {
        glColor3f(0, 1, 0);
      } else {
        glColor3f(1, 0.1, 0.1);
      }
      glTranslated(textLineP2[0], textLineP2[1], textLineP2[2]);
      glutSolidSphere(0.001, 16, 8);
    }
    glPopMatrix();
    
    Point3D tlTmp = textLineP1;
    tlTmp[2] += 0.00001;
    
    RenderInfo const &renderInfo = sphPtr->renderInfo;
    
    static StrFTime const guiETAFormat("%H:%M:%S");

    ostringstream infoText;
    std::right(infoText);
    std::fixed(infoText);
    infoText.precision(4);
//     std::scientific(infoText);
    infoText << "#p " << renderInfo.npartikel - renderInfo.numParticlesOut
             << "/" << renderInfo.npartikel
             << " step " << setw(6) << renderInfo.step << ""
             << " t " << setw(6) << renderInfo.t << " s"
             << " time/step " << renderInfo.timePerStep << " s"
             << " thr " << renderInfo.numThreads << ""
             << " fr " << renderInfo.frame
             << " eta " << guiETAFormat(renderInfo.eta)
             << "\n";
    // set font color
    glColor3f(1, 0.721, 0);
    GLString::renderString(GLUT_STROKE_MONO_ROMAN, infoText.str(), 
                           1.1e-5, tlTmp);
    
    tlTmp[1] -= 0.002;
    infoText.str("");
    Point3D rotAxis = camera.rotation().vector();
    double const raLen = norm(rotAxis);
    if (raLen > 0) {
      rotAxis /= raLen;
    }
    //       CoordSystem local;
    //       local.rotate(camera.rotation());
    //       cerr << "local axes: " << local << "\n";
    infoText 
      << "pos. " << std::setprecision(3) << camera.position() << " m"
//         << " rot: " << rad2deg(camera.rotation().angle()) << "° " << rotAxis
//         << " phi: " << rad2deg(angle(local.axes[0], coordSystem.axes[0]))
//         << " theta: " << rad2deg(angle(local.axes[1], coordSystem.axes[1]))
//         << " roll: " << rad2deg(angle(local.axes[2], coordSystem.axes[2]))
      << ", or. " << std::setprecision(2)
      << rad2deg(camera.rotation().angle()) << "° " << rotAxis
      << ", FOV. " << camera.fov() << "°"
      << ", vec. sc.: 2^" << options.vectorLengthExp << ""
      << ", part. sc.: 2^" << options.particleSizeExp << ""
      << std::setprecision(5)
      << ", values: " << renderInfo.minColor << "-" << renderInfo.maxColor;

    // set font color
    glColor3f(0.65, 0.65, 0.65);
    GLString::renderString(GLUT_STROKE_ROMAN, infoText.str(), 
                           1.1e-5, tlTmp);
    
  }

  static void rotateScene() {
    double const rotAngle = rad2deg(camera.rotation().angle());
    if (rotAngle != 0) {
      Point3D rotAxis = camera.rotation().vector();
      rotAxis /= norm(rotAxis);
      glRotated(rotAngle, rotAxis[0], rotAxis[1], rotAxis[2]);
    }
  }

  static void setupScene() {
    setupView();
    rotateScene();
    glTranslated(-camera.position()[0], -camera.position()[1], -camera.position()[2]);
  }

  static void drawRulerLine(Point3D const &a, Point3D const &b) {
    double const tlength = 0.05;
    glBegin(GL_LINES);
    {
      glVertex3dv(&a[0]);
      glVertex3dv(&b[0]);
      Point3D toCross;
      toCross = 1;
      Point3D const diff = b - a;
      double const length = norm(diff);
      Point3D const diffUnit = diff/length;
      Point3D tick1 = cross(toCross, diff);
      tick1 *= tlength / norm(tick1);
      Point3D tick2 = cross(tick1, diff);
      tick2 *= tlength / norm(tick2);
      Point3D mtick1 = tick1 * 5e-1;
      Point3D mtick2 = tick2 * 5e-1;
      glColor3d(.7, .7, .7);
      for (int tick = 0; tick < length; ++tick) {
        if (tick > 0) {
          glVertex3dv(&(a + diffUnit*double(tick) + tick1)[0]);
          glVertex3dv(&(a + diffUnit*double(tick) - tick1)[0]);
          glVertex3dv(&(a + diffUnit*double(tick) + tick2)[0]);
          glVertex3dv(&(a + diffUnit*double(tick) - tick2)[0]);
        }
        glColor3d(.5, .5, .5);
        for (int mtick = 1; mtick < 10 and tick + mtick*1e-1 < length; ++mtick) {
          glVertex3dv(&(a + diffUnit*(tick + mtick*1e-1) + mtick1)[0]);
          glVertex3dv(&(a + diffUnit*(tick + mtick*1e-1) - mtick1)[0]);
          glVertex3dv(&(a + diffUnit*(tick + mtick*1e-1) + mtick2)[0]);
          glVertex3dv(&(a + diffUnit*(tick + mtick*1e-1) - mtick2)[0]);
        }
      }
    }
    glEnd();
  }

  static void compileDisplayLists() {
    dlCompileCoordAxes();
    dlCompileUnitSpheres();
    dlCompileSensors();
  }

  static void dlCompileUnitSpheres() {
    glNewList(dlIdSolidUnitSphere, GL_COMPILE);
    {
      glutSolidSphere(1, 8, 8);
    }
    glEndList();
    glNewList(dlIdWireUnitSphere, GL_COMPILE);
    {
      glutWireSphere(1, 8, 8);
    }
    glEndList();
  }

  static void dlCompileCoordAxes() {
    glNewList(CoordAxesDisplayListId, GL_COMPILE);
    {
      Point3D null, ende;
      
      ende = 0;
      ende[0] = sphPtr->config.gebiet.obenRechts[0];
      glColor3d(1, 0, 0);
      drawRulerLine(null, ende);
      
      glColor3d(1, 0.721, 0);
      Point3D labelPos;
      double length = sphPtr->config.gebiet.obenRechts[0];
      labelPos[1] = -0.1;
      for (int tick = 1; tick < length; ++tick) {
        labelPos[0] = tick;
        std::ostringstream nstr;
        nstr << tick;
        GLString::renderString(GLUT_STROKE_ROMAN, nstr.str(), .6e-3, labelPos);
      }
      GLString::renderString(GLUT_STROKE_ROMAN, "X", .6e-3, ende);

#if SPH_DIM > 1
      ende = 0;
      ende[1] = sphPtr->config.gebiet.obenRechts[1];
      glColor3d(0, 1, 0);
      drawRulerLine(null, ende);
      
      glColor3d(1, 0.721, 0);
      labelPos = 0;
      labelPos[0] = -0.1;
      length = sphPtr->config.gebiet.obenRechts[1];
      for (int tick = 1; tick < length; ++tick) {
        labelPos[1] = tick;
        std::ostringstream nstr;
        nstr << tick;
        GLString::renderString(GLUT_STROKE_ROMAN, nstr.str(), .6e-3, labelPos);
      }
      GLString::renderString(GLUT_STROKE_ROMAN, "Y", .6e-3, ende);
#endif

#if SPH_DIM > 2
      ende = 0;
      ende[2] = sphPtr->config.gebiet.obenRechts[2];
      glColor3d(0, 0, 1);
      drawRulerLine(null, ende);
      
      glColor3d(1, 0.721, 0);
      labelPos = 0;
      labelPos[0] = -0.1;
      length = sphPtr->config.gebiet.obenRechts[2];
      for (int tick = 1; tick < length; ++tick) {
        labelPos[1] = tick;
        std::ostringstream nstr;
        nstr << tick;
        GLString::renderString(GLUT_STROKE_ROMAN, nstr.str(), .6e-3, labelPos);
      }
      GLString::renderString(GLUT_STROKE_ROMAN, "Z", .6e-3, ende);
#endif

    }
    glEndList();
  }

  static void dlCompileSensors() {
    glNewList(dlIdSensors, GL_COMPILE);
    {
      sphPtr->doRenderSensors();
    }
    glEndList();
  }

  static void displayFunction() {
//     cerr << "render: frame " << frameRendered << "!\n";
    ++frameRendered;
    
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glPushMatrix();
    {

      setupScene();

      // draw global coord system
      if (not options.hideCoordAxes) {
        glCallList(CoordAxesDisplayListId);
      }

      sphPtr->render(options);
    }
    glPopMatrix();
    
    displayOverlayFunction();

    glutSwapBuffers();
  }

  static void showKeyboardHelp() {
    cout << 
      "\n"
      "SPH Keyboard Help: \n"
      " * Camera movements *\n"
      "<right>          rotate camera right\n"
      "<left>           rotate camera left\n"
      "<up>             rotate camera up\n"
      "<down>           rotate camera down\n"
      "<page-up>        roll camera left\n"
      "<page-down>      roll camera right\n"
      "w                move camera forward\n"
      "s                move camera back\n"
      "a                move camera left\n"
      "d                move camera right\n"
      "<scroll-up>      zoom in\n"
      "<scroll-down>    zoom out\n"
      "<backspace>      reset to default view\n"
      " shift: 10 times slower, ctrl: 10 times faster\n"
      "\n"
      " * Display *\n"
      "ctrl-<scroll>    in-/decrease particle size (if representation is sphere)\n"
      "shift-<scroll>   in-/decrease vector lengths (if vectors are shown)\n"
      "c                toggle coloring by color/by hashed raster\n"
      "g                show grid\n"
      "<space bar>      update view\n"
      "<right-click>    open on-screen menu\n"
      "\n"
      " * Threads *\n"
      "+             use one more thread\n"
      "-             use one thread less\n"
      "*             use twice the number of threads\n"
      "/             use half the number of threads\n"
      "\n"
      " * General *\n"
      "q             quit regularly\n"
      "Q             call exit(EXIT_FAILURE)\n"
      "A             call abort()\n"
      "h <F1>        show this help\n"
      "\n"
      ;
  }

  static void toggleCoordSystem() {
    options.hideCoordAxes = !options.hideCoordAxes;
  }
  static void toggleRaster() {
    options.showRaster = !options.showRaster;
  }
  static void toggleSensors() {
    options.showSensors = !options.showSensors;
  }
  static void toggleSimulationArea() {
    options.hideSimulationArea = !options.hideSimulationArea;
  }
  static void toggleBoundaryParticles() {
    options.hideBoundaryParticles = !options.hideBoundaryParticles;
  }
  static void toggleOutParticles() {
    options.showOutParticles = !options.showOutParticles;
  }

  static void updateVectorLength() {
    options.vectorLength = ldexp(1.0, options.vectorLengthExp);
  }
  static void updateParticleSize() {
    options.particleSize = ldexp(1.0, options.particleSizeExp);
  }
  static void decrVectorLength() {
    options.vectorLengthExp -= 1;
    options.vectorLengthExp = max(int(options.vectorLengthExp), -127);
    updateVectorLength();
  }
  static void incrVectorLength() {
    options.vectorLengthExp += 1;
    options.vectorLengthExp = min(int(options.vectorLengthExp), 127);
    updateVectorLength();
  }

  static void decrParticleSize() {
    options.particleSizeExp -= 1;
    options.particleSizeExp = max(int(options.particleSizeExp), -127);
    updateParticleSize();
  }
  static void incrParticleSize() {
    options.particleSizeExp += 1;
    options.particleSizeExp = min(int(options.particleSizeExp), 127);
    updateParticleSize();
  }

  static double getRotMovMult(int const modifiers) {
    double movMult = 1;
//     if (modifiers & GLUT_ACTIVE_ALT) 
//       movMult *= 4;
    if (modifiers & GLUT_ACTIVE_SHIFT) 
      movMult *= 0.1;
    if (modifiers & GLUT_ACTIVE_CTRL) 
      movMult *= 10;
    return movMult;
  }

  static void specialFunction(int key, int /* x */, int /* y */) {

    int const modi = glutGetModifiers();
    double const movMult = getRotMovMult(modi);
    
    switch (key) {
    case GLUT_KEY_F1:
      showKeyboardHelp();
      break;

    case GLUT_KEY_RIGHT: {
      camera.rotate(rUnit*movMult, 1);
    }
      break;
    case GLUT_KEY_LEFT: {
      camera.rotate(-rUnit*movMult, 1);
    }
      break;
      
    case GLUT_KEY_UP: {
      camera.rotate(rUnit*movMult, 0);
    }
      break;
    case GLUT_KEY_DOWN:
      camera.rotate(-rUnit*movMult, 0);
      break;
      
    case GLUT_KEY_PAGE_UP: {
      camera.rotate(rUnit*movMult, 2);
    }
      break;
    case GLUT_KEY_PAGE_DOWN: {
      camera.rotate(-rUnit*movMult, 2);
    }
      break;
      
    default:
      cerr << "special key event: shift-" << int(key) << "\n";
      break;
    }

    if (camera.updateScreen()) {
      redisplay();
    }
  }

  static void keyboardFunction(unsigned char key, int , int ) {
    int const modi = glutGetModifiers();
//     cerr << "event: key event " << int(key) << " " << x << " " << y << "\n";
    bool pressedShift = false, pressedAlt = false, pressedCtrl = false;
    if (modi & GLUT_ACTIVE_ALT) {
      pressedAlt = true;
    }
    if (modi & GLUT_ACTIVE_SHIFT) {
      key += 32;
      pressedShift = true;
    }
    if (modi & GLUT_ACTIVE_CTRL) {
      key += 64 + 32;
      pressedCtrl = true;
    }

    double const movMult = getRotMovMult(modi);
    unsigned newNumThreads = 0;
    bool redisp = false;

    switch (key) {

    case '8': 
      camera.rotate(rUnit*movMult, 0);
      break;
    case '2': 
      camera.rotate(-rUnit*movMult, 0);
      break;
      
    case '4': 
      camera.rotate(rUnit*movMult, 1);
      break;
    case '6': 
      camera.rotate(-rUnit*movMult, 1);
      break;

    case '9': 
      camera.rotate(rUnit*movMult, 2);
      break;
    case '3': 
      camera.rotate(-rUnit*movMult, 2);
      break;

//     case '+':
// //       cerr << "zoom in\n";
//       fovAngle -= 1;
// //       break;
//     case '-':
// //       cerr << "zoom out\n";
//       fovAngle += 1;
// //       break;

//     case 'R':
//     case 'r':
// //       rotationAngle[1] += movMult*1;
//       break;

//     case 'L':
//     case 'l':
// //       rotationAngle[1] -= movMult*1;
//       break;

    case 27: // escape
      renderQuit();
      //       resetWindow();
      break;

    case 8: // backspace
      resetView();
      break;

    case 'c':
      if (pressedShift) {
        toggleCoordSystem();
      } else {
        options.colorRaster = !options.colorRaster;
      }
      break;

    case 'g':
      toggleRaster();
      break;
    case 'f':
      toggleSimulationArea();
      break;

     case 'i':
       if (pressedCtrl) {
         glutIconifyWindow();
       }
       break;


    case 's':
      camera.translate(movMult*mUnit, 2);
      break;
    case 'w':
      camera.translate(-movMult*mUnit, 2);
      break;

    case 'a':
      camera.translate(-movMult*mUnit, 0);
      break;
    case 'd':
      camera.translate(movMult*mUnit, 0);
      break;

    case 'q':
      camera.translate(movMult*mUnit, 1);
      break;
    case 'z':
    case 'y':
      camera.translate(-movMult*mUnit, 1);
      break;

    case 'H':
    case 'h':
      showKeyboardHelp();
      break;

    case ' ':
      redisp = true;
      break;

   
    case 139: // ctrl-+
      if (pressedCtrl) {
        camera.zoomIn(1);
      }
      break;
    case 141: // ctrl--
      if (pressedCtrl) {
        camera.zoomOut(1);
      }
      break;

    case '+':
      newNumThreads = sphPtr->numThreads + 1;
      break;
    case '-':
      newNumThreads = sphPtr->numThreads - 1;
      break;
    case '*':
      newNumThreads = sphPtr->numThreads << 1;
      break;
    case '/':
      newNumThreads = sphPtr->numThreads >> 1;
      break;

    default:
      cerr << "event: key event `";
      if (pressedAlt) {
        cerr << "alt-";
      }
      if (pressedShift) {
        cerr << "shift-";
      }
      if (pressedCtrl) {
        cerr << "ctrl-";
      }
      if (isprint(key)) {
        cerr << key << "'\n";
      } else {
        cerr << "code-" << int(key) << "'\n";
      }
      break;

    }

    if (redisp or camera.updateView() or camera.updateScreen()) {
      redisplay();
    }
    
    if (newNumThreads > 0 and newNumThreads <= sphPtr->config.maxNumThreads) {
      sphPtr->numThreads = newNumThreads;
    }
  }

#include "menu-entries.ncd.enum.hh"
#include "menu-entries.ncd.cc"

  static void createMenu() {

    void (*menuFunction)(int const value) = menuHandlerFunction;

    int submenuIdFile = 
      glutCreateMenu(menuFunction);
    glutAddMenuEntry("Quit", MENU_QUIT);
    glutAddMenuEntry("Exit (Success)", MENU_EXIT_SUCCESS);
    glutAddMenuEntry("Exit (Failure)", MENU_EXIT_FAILURE);
    glutAddMenuEntry("Abort", MENU_ABORT);

    int submenuIdHelp = 
      glutCreateMenu(menuFunction);
    glutAddMenuEntry("Help", MENU_HELP);
    glutAddMenuEntry("About", MENU_ABOUT);

#ifdef SPH_AD
    int submenuIdVecPositionAD = 
      glutCreateMenu(menuFunction);
//     glutAddMenuEntry("Value", MENU_PARTICLE_COLOR_POSITION);
    glutAddMenuEntry("D0", MENU_PARTICLE_VEC_COORD_D0);
    glutAddMenuEntry("D1", MENU_PARTICLE_VEC_COORD_D1);
    glutAddMenuEntry("D2", MENU_PARTICLE_VEC_COORD_D2);

    int submenuIdVecVelocityAD = 
      glutCreateMenu(menuFunction);
    glutAddMenuEntry("Value", MENU_PARTICLE_VEC_VELOCITY);
    glutAddMenuEntry("D0", MENU_PARTICLE_VEC_VELOCITY_D0);
    glutAddMenuEntry("D1", MENU_PARTICLE_VEC_VELOCITY_D1);
    glutAddMenuEntry("D2", MENU_PARTICLE_VEC_VELOCITY_D2);

    int submenuIdVecDensityAD = 
      glutCreateMenu(menuFunction);
    glutAddMenuEntry("Value", MENU_PARTICLE_VEC_DENSITY);
    glutAddMenuEntry("D0", MENU_PARTICLE_VEC_DENSITY_D0);
    glutAddMenuEntry("D1", MENU_PARTICLE_VEC_DENSITY_D1);
    glutAddMenuEntry("D2", MENU_PARTICLE_VEC_DENSITY_D2);

    int submenuIdVecPressureAD = 
      glutCreateMenu(menuFunction);
    glutAddMenuEntry("Value", MENU_PARTICLE_VEC_PRESSURE);
    glutAddMenuEntry("D0", MENU_PARTICLE_VEC_PRESSURE_D0);
    glutAddMenuEntry("D1", MENU_PARTICLE_VEC_PRESSURE_D1);
    glutAddMenuEntry("D2", MENU_PARTICLE_VEC_PRESSURE_D2);

    int submenuIdVecHeatAD = 
      glutCreateMenu(menuFunction);
    glutAddMenuEntry("Value", MENU_PARTICLE_VEC_HEAT);
    glutAddMenuEntry("D0", MENU_PARTICLE_VEC_HEAT_D0);
    glutAddMenuEntry("D1", MENU_PARTICLE_VEC_HEAT_D1);
    glutAddMenuEntry("D2", MENU_PARTICLE_VEC_HEAT_D2);
#endif

#ifdef SPH_AD
    int submenuIdPartVec = 
      glutCreateMenu(menuFunction);
    glutAddSubMenu("Position",      submenuIdVecPositionAD);
    glutAddSubMenu("Velocity",      submenuIdVecVelocityAD);
    glutAddSubMenu("Density",       submenuIdVecDensityAD);
    glutAddSubMenu("Pressure",      submenuIdVecPressureAD);
    glutAddSubMenu("Heat",          submenuIdVecHeatAD);
#else
    int submenuIdPartVec = 
      glutCreateMenu(menuFunction);
    glutAddMenuEntry("Velocity",    MENU_PARTICLE_VEC_VELOCITY);
    glutAddMenuEntry("Density",     MENU_PARTICLE_VEC_DENSITY);
    glutAddMenuEntry("Pressure",    MENU_PARTICLE_VEC_PRESSURE);
    glutAddMenuEntry("Heat",        MENU_PARTICLE_VEC_HEAT);
#endif

    int submenuIdPartSize = 
      glutCreateMenu(menuFunction);
    glutAddMenuEntry("Kernel-Reach (h = H/2)",    MENU_PARTICLE_SIZE_H);
    glutAddMenuEntry("Volume",               MENU_PARTICLE_SIZE_VOL);

    int submenuIdPart = 
      glutCreateMenu(menuFunction);
    glutAddMenuEntry("Solid Sphere", MENU_PARTICLE_AS_SOLID_SPHERE);
    glutAddMenuEntry("Wire Sphere", MENU_PARTICLE_AS_WIRE_SPHERE);
    glutAddMenuEntry("Point", MENU_PARTICLE_AS_POINT);
    glutAddSubMenu("Size", submenuIdPartSize);


    int submenuIdColorMap = 
      glutCreateMenu(menuFunction);
    glutAddMenuEntry("Gradient Red-Blue", MENU_COLOR_MAP_GRADIENT_RB);
    glutAddMenuEntry("Gradient Red-Green", MENU_COLOR_MAP_GRADIENT_RG);
    glutAddMenuEntry("Gradient Green-Blue", MENU_COLOR_MAP_GRADIENT_GB);
    glutAddMenuEntry("User Map", MENU_COLOR_MAP_MAP);

    int submenuIdColorRange = 
      glutCreateMenu(menuFunction);
    glutAddMenuEntry("All", MENU_COLOR_RANGE_ALL);
    glutAddMenuEntry("Moving only", MENU_COLOR_RANGE_MOVING);

#ifdef SPH_AD
    int submenuIdPositionAD = 
      glutCreateMenu(menuFunction);
//     glutAddMenuEntry("Value", MENU_PARTICLE_COLOR_POSITION);
    glutAddMenuEntry("D0", MENU_PARTICLE_COLOR_COORD_D0);
    glutAddMenuEntry("D1", MENU_PARTICLE_COLOR_COORD_D1);
    glutAddMenuEntry("D2", MENU_PARTICLE_COLOR_COORD_D2);

    int submenuIdVelocityAD = 
      glutCreateMenu(menuFunction);
    glutAddMenuEntry("Value", MENU_PARTICLE_COLOR_VELOCITY);
    glutAddMenuEntry("D0", MENU_PARTICLE_COLOR_VELOCITY_D0);
    glutAddMenuEntry("D1", MENU_PARTICLE_COLOR_VELOCITY_D1);
    glutAddMenuEntry("D2", MENU_PARTICLE_COLOR_VELOCITY_D2);

    int submenuIdDensityAD = 
      glutCreateMenu(menuFunction);
    glutAddMenuEntry("Value", MENU_PARTICLE_COLOR_DENSITY);
    glutAddMenuEntry("D0", MENU_PARTICLE_COLOR_DENSITY_D0);
    glutAddMenuEntry("D1", MENU_PARTICLE_COLOR_DENSITY_D1);
    glutAddMenuEntry("D2", MENU_PARTICLE_COLOR_DENSITY_D2);

    int submenuIdPressureAD = 
      glutCreateMenu(menuFunction);
    glutAddMenuEntry("Value", MENU_PARTICLE_COLOR_PRESSURE);
    glutAddMenuEntry("D0", MENU_PARTICLE_COLOR_PRESSURE_D0);
    glutAddMenuEntry("D1", MENU_PARTICLE_COLOR_PRESSURE_D1);
    glutAddMenuEntry("D2", MENU_PARTICLE_COLOR_PRESSURE_D2);

    int submenuIdHeatAD = 
      glutCreateMenu(menuFunction);
    glutAddMenuEntry("Value", MENU_PARTICLE_COLOR_HEAT);
    glutAddMenuEntry("D0", MENU_PARTICLE_COLOR_HEAT_D0);
    glutAddMenuEntry("D1", MENU_PARTICLE_COLOR_HEAT_D1);
    glutAddMenuEntry("D2", MENU_PARTICLE_COLOR_HEAT_D2);
#endif

    int submenuIdColor = 
      glutCreateMenu(menuFunction);
#ifdef SPH_AD
    glutAddSubMenu("Position", submenuIdPositionAD);
    glutAddSubMenu("Velocity", submenuIdVelocityAD);
    glutAddSubMenu("Density", submenuIdDensityAD);
    glutAddSubMenu("Pressure", submenuIdPressureAD);
    glutAddSubMenu("Heat", submenuIdHeatAD);
#else
    glutAddMenuEntry("Velocity", MENU_PARTICLE_COLOR_VELOCITY);
    glutAddMenuEntry("Density", MENU_PARTICLE_COLOR_DENSITY);
    glutAddMenuEntry("Pressure", MENU_PARTICLE_COLOR_PRESSURE);
    glutAddMenuEntry("Heat", MENU_PARTICLE_COLOR_HEAT);
#endif
    glutAddMenuEntry("Color", MENU_PARTICLE_COLOR_COLOR);
    glutAddMenuEntry("Class", MENU_PARTICLE_COLOR_CLASS);
    glutAddMenuEntry("Flags", MENU_PARTICLE_COLOR_FLAGS);
    glutAddMenuEntry("Hash", MENU_PARTICLE_COLOR_HASH);
    glutAddMenuEntry("Id", MENU_PARTICLE_COLOR_ID);
    glutAddMenuEntry("Neighbours", MENU_PARTICLE_COLOR_NUM_NEIGHBOURS);
    glutAddMenuEntry("Index", MENU_PARTICLE_COLOR_INDEX);
    glutAddSubMenu("Color Map", submenuIdColorMap);
    glutAddSubMenu("Color Range", submenuIdColorRange);

    int submenuIdGraphics = 
      glutCreateMenu(menuFunction);
    glutAddMenuEntry("Reset Window", MENU_RESET_WINDOW);
    glutAddMenuEntry("Reset Camera", MENU_RESET_CAMERA);
    glutAddMenuEntry("Fullscreen", MENU_FULLSCREEN);

    int submenuIdHideShow = 
      glutCreateMenu(menuFunction);
    glutAddMenuEntry("Vectors", MENU_SHOW_VECTORS);
    glutAddMenuEntry("Coordinate Axes",   MENU_SHOW_COORD_SYSTEM);
    glutAddMenuEntry("Sensors",   MENU_SHOW_SENSORS);
    glutAddMenuEntry("Simulation Area",   MENU_SHOW_SIMULATION_AREA);
    glutAddMenuEntry("Raster",   MENU_SHOW_RASTER);
    glutAddMenuEntry("Boundary particles", MENU_SHOW_BOUNDARY_PARTICLES);
    glutAddMenuEntry("Out particles", MENU_SHOW_OUT_PARTICLES);

    int submenuIdView = 
      glutCreateMenu(menuFunction);
    glutAddSubMenu("Hide/Show",   submenuIdHideShow);
    glutAddSubMenu("Graphics",   submenuIdGraphics);
    glutAddSubMenu("Vectors",   submenuIdPartVec);
    glutAddSubMenu("Particles", submenuIdPart);
    glutAddSubMenu("Color", submenuIdColor);

    glutCreateMenu(menuFunction);
    glutAddSubMenu("File", submenuIdFile);
    glutAddSubMenu("View", submenuIdView);
    glutAddSubMenu("Help", submenuIdHelp);

    glutAttachMenu(GLUT_RIGHT_BUTTON);
  }

  static void menuHandlerFunction(int const value) {
    switch(value) {
    case MENU_QUIT:
      renderQuit();
      break;
    case MENU_ABORT:
      renderAbort();
      break;
    case MENU_EXIT_SUCCESS:
      renderExitSuccess();
      break;
    case MENU_EXIT_FAILURE:
      renderExitFailure();
      break;

    case MENU_RESET_WINDOW:
      resetWindow();
      break;
    case MENU_RESET_CAMERA:
      resetView();
      break;
    case MENU_FULLSCREEN:
      glutFullScreen();
      break;

//     case MENU_PARTICLE_VEC_SHORTER:
//       decrVectorLength();
//       break;
//     case MENU_PARTICLE_VEC_LONGER:
//       incrVectorLength();
//       break;

    case MENU_SHOW_VECTORS:
      options.showParticleVectors = !options.showParticleVectors;
      break;
    case MENU_SHOW_COORD_SYSTEM:
      toggleCoordSystem();
      break;
    case MENU_SHOW_SENSORS:
      toggleSensors();
      break;
    case MENU_SHOW_RASTER:
      toggleRaster();
      break;
    case MENU_SHOW_SIMULATION_AREA:
      toggleSimulationArea();
      break;
    case MENU_SHOW_BOUNDARY_PARTICLES:
      toggleBoundaryParticles();
      break;
    case MENU_SHOW_OUT_PARTICLES:
      toggleOutParticles();
      break;

    case MENU_PARTICLE_SIZE_H:
      sphPtr->particleSizeFunction = &SPHT::particleRadiusH;
      break;
    case MENU_PARTICLE_SIZE_VOL:
      sphPtr->particleSizeFunction = &SPHT::particleRadiusVolume;
      break;

    case MENU_PARTICLE_VEC_VELOCITY:
      sphPtr->particleVectorFunction = &SPHT::particleVectorVelocity;
      break;

    case MENU_PARTICLE_VEC_DENSITY:
      sphPtr->particleVectorFunction = &SPHT::particleVectorDensity;
      break;
    case MENU_PARTICLE_VEC_PRESSURE:
      sphPtr->particleVectorFunction = &SPHT::particleVectorPressure;
      break;
    case MENU_PARTICLE_VEC_HEAT:
      sphPtr->particleVectorFunction = &SPHT::particleVectorHeat;
      break;

#ifdef SPH_AD
    case MENU_PARTICLE_VEC_COORD_D0:
      sphPtr->particleVectorFunction = &SPHT::particleVectorPositionDiff<0>;
      break;
    case MENU_PARTICLE_VEC_COORD_D1:
      sphPtr->particleVectorFunction = &SPHT::particleVectorPositionDiff<1>;
      break;
    case MENU_PARTICLE_VEC_COORD_D2:
      sphPtr->particleVectorFunction = &SPHT::particleVectorPositionDiff<2>;
      break;

    case MENU_PARTICLE_VEC_VELOCITY_D0:
      sphPtr->particleVectorFunction = &SPHT::particleVectorVelocityDiff<0>;
      break;
    case MENU_PARTICLE_VEC_VELOCITY_D1:
      sphPtr->particleVectorFunction = &SPHT::particleVectorVelocityDiff<1>;
      break;
    case MENU_PARTICLE_VEC_VELOCITY_D2:
      sphPtr->particleVectorFunction = &SPHT::particleVectorVelocityDiff<2>;
      break;

    case MENU_PARTICLE_VEC_DENSITY_D0:
      sphPtr->particleVectorFunction = &SPHT::particleVectorDensityDiff<0>;
      break;
    case MENU_PARTICLE_VEC_DENSITY_D1:
      sphPtr->particleVectorFunction = &SPHT::particleVectorDensityDiff<1>;
      break;
    case MENU_PARTICLE_VEC_DENSITY_D2:
      sphPtr->particleVectorFunction = &SPHT::particleVectorDensityDiff<2>;
      break;

    case MENU_PARTICLE_VEC_PRESSURE_D0:
      sphPtr->particleVectorFunction = &SPHT::particleVectorPressureDiff<0>;
      break;
    case MENU_PARTICLE_VEC_PRESSURE_D1:
      sphPtr->particleVectorFunction = &SPHT::particleVectorPressureDiff<1>;
      break;
    case MENU_PARTICLE_VEC_PRESSURE_D2:
      sphPtr->particleVectorFunction = &SPHT::particleVectorPressureDiff<2>;
      break;

    case MENU_PARTICLE_VEC_HEAT_D0:
      sphPtr->particleVectorFunction = &SPHT::particleVectorHeatDiff<0>;
      break;
    case MENU_PARTICLE_VEC_HEAT_D1:
      sphPtr->particleVectorFunction = &SPHT::particleVectorHeatDiff<1>;
      break;
    case MENU_PARTICLE_VEC_HEAT_D2:
      sphPtr->particleVectorFunction = &SPHT::particleVectorHeatDiff<2>;
      break;
#endif

    case MENU_PARTICLE_COLOR_VELOCITY:
      sphPtr->setColorFunction(&SPHT::particleColorVelocity);
      break;
    case MENU_PARTICLE_COLOR_DENSITY:
      sphPtr->setColorFunction(&SPHT::particleColorDensity);
      break;
    case MENU_PARTICLE_COLOR_PRESSURE:
      sphPtr->setColorFunction(&SPHT::particleColorPressure);
      break;
    case MENU_PARTICLE_COLOR_HEAT:
      sphPtr->setColorFunction(&SPHT::particleColorHeat);
      break;

    case MENU_PARTICLE_COLOR_COLOR:
      sphPtr->setColorFunction(&SPHT::particleColorColor);
      break;
    case MENU_PARTICLE_COLOR_CLASS:
      sphPtr->setColorFunction(&SPHT::particleColorClass);
      break;
    case MENU_PARTICLE_COLOR_NUM_NEIGHBOURS:
      sphPtr->setColorFunction(&SPHT::particleColorNumNeighbours);
      break;
    case MENU_PARTICLE_COLOR_ID:
      sphPtr->setColorFunction(&SPHT::particleColorId);
      break;
    case MENU_PARTICLE_COLOR_FLAGS:
      sphPtr->setColorFunction(&SPHT::particleColorFlags);
      break;
    case MENU_PARTICLE_COLOR_HASH:
      sphPtr->setColorFunction(&SPHT::particleColorHash);
      break;
    case MENU_PARTICLE_COLOR_INDEX:
      sphPtr->setColorFunction(&SPHT::particleColorIndex);
      break;

#ifdef SPH_AD
    case MENU_PARTICLE_COLOR_COORD_D0:
      sphPtr->particleColorFunction = &SPHT::particleColorPositionDiff<0>;
      break;
    case MENU_PARTICLE_COLOR_COORD_D1:
      sphPtr->particleColorFunction = &SPHT::particleColorPositionDiff<1>;
      break;
    case MENU_PARTICLE_COLOR_COORD_D2:
      sphPtr->particleColorFunction = &SPHT::particleColorPositionDiff<2>;
      break;

    case MENU_PARTICLE_COLOR_VELOCITY_D0:
      sphPtr->particleColorFunction = &SPHT::particleColorVelocityDiff<0>;
      break;
    case MENU_PARTICLE_COLOR_VELOCITY_D1:
      sphPtr->particleColorFunction = &SPHT::particleColorVelocityDiff<1>;
      break;
    case MENU_PARTICLE_COLOR_VELOCITY_D2:
      sphPtr->particleColorFunction = &SPHT::particleColorVelocityDiff<2>;
      break;

    case MENU_PARTICLE_COLOR_DENSITY_D0:
      sphPtr->particleColorFunction = &SPHT::particleColorDensityDiff<0>;
      break;
    case MENU_PARTICLE_COLOR_DENSITY_D1:
      sphPtr->particleColorFunction = &SPHT::particleColorDensityDiff<1>;
      break;
    case MENU_PARTICLE_COLOR_DENSITY_D2:
      sphPtr->particleColorFunction = &SPHT::particleColorDensityDiff<2>;
      break;

    case MENU_PARTICLE_COLOR_PRESSURE_D0:
      sphPtr->particleColorFunction = &SPHT::particleColorPressureDiff<0>;
      break;
    case MENU_PARTICLE_COLOR_PRESSURE_D1:
      sphPtr->particleColorFunction = &SPHT::particleColorPressureDiff<1>;
      break;
    case MENU_PARTICLE_COLOR_PRESSURE_D2:
      sphPtr->particleColorFunction = &SPHT::particleColorPressureDiff<2>;
      break;

    case MENU_PARTICLE_COLOR_HEAT_D0:
      sphPtr->particleColorFunction = &SPHT::particleColorHeatDiff<0>;
      break;
    case MENU_PARTICLE_COLOR_HEAT_D1:
      sphPtr->particleColorFunction = &SPHT::particleColorHeatDiff<1>;
      break;
    case MENU_PARTICLE_COLOR_HEAT_D2:
      sphPtr->particleColorFunction = &SPHT::particleColorHeatDiff<2>;
      break;
#endif

    case MENU_COLOR_MAP_MAP:
      sphPtr->colorMapFunction = &SPHT::colorMapMap;
      break;
    case MENU_COLOR_MAP_GRADIENT_RG:
      sphPtr->colorMapFunction = &SPHT::colorMapGradientRG;
      break;
    case MENU_COLOR_MAP_GRADIENT_RB:
      sphPtr->colorMapFunction = &SPHT::colorMapGradient;
      break;
    case MENU_COLOR_MAP_GRADIENT_GB:
      sphPtr->colorMapFunction = &SPHT::colorMapGradientGB;
      break;

    case MENU_COLOR_RANGE_ALL:
      options.colorRangeAll = true;
      sphPtr->clearColorScale();
      break;
    case MENU_COLOR_RANGE_MOVING:
      options.colorRangeAll = false;
      sphPtr->clearColorScale();
      break;

    case MENU_PARTICLE_AS_POINT:
      sphPtr->particleReprFunction = &SPHT::particleReprPoint;
      break;
    case MENU_PARTICLE_AS_SOLID_SPHERE:
      sphPtr->particleReprFunction = &SPHT::particleReprSolidSphere;
      break;
    case MENU_PARTICLE_AS_WIRE_SPHERE:
      sphPtr->particleReprFunction = &SPHT::particleReprWireSphere;
      break;

    case MENU_HELP:
      showKeyboardHelp();
      break;

    default:
      cerr << "menu value " << getMenuEntriesName(value) << " not implemented\n";
      break;
    }
  }

  static void redisplay() {
    glutPostRedisplay();
    camera.clear();
  }
  static void renderQuit() {
    beende = true;
  }
  static void renderExitSuccess() {
    exit(EXIT_SUCCESS);
  }
  static void renderExitFailure() {
    exit(EXIT_FAILURE);
  }
  static void renderAbort() {
    abort();
  }

  static void mouseFunction(int button, int state, int , int ) {
    int const modi = glutGetModifiers();
    bool pressedShift = false, 
      // pressedAlt = false, 
      pressedCtrl = false;
    // if (modi & GLUT_ACTIVE_ALT) {
    //   pressedAlt = true;
    // }
    if (modi & GLUT_ACTIVE_SHIFT) {
        pressedShift = true;
    }
    if (modi & GLUT_ACTIVE_CTRL) {
      pressedCtrl = true;
    }

//     cerr << "event: mouse event button " << button << " state "
// 	 << state << " x " << x << " y " << y << "\n";

    mouse.buttons[button] = state;

    if (state == GLUT_UP) {
      switch(button) {
      case 3:
        if (pressedShift) {
          incrVectorLength();
        } else if (pressedCtrl) {
          incrParticleSize();
        } else {
          camera.zoomIn(1);
        }
        break;
      case 4:
        if (pressedShift) {
          decrVectorLength();
        } else if (pressedCtrl) {
          decrParticleSize();
        } else {
          camera.zoomOut(1);
        }
        break;
      }
    }

    if (camera.updateView()) {
//       setupOpenGL();
      redisplay();
    } else if (camera.updateScreen()) {
      redisplay();
    }
  }

  static void entryFunction(int const state) {
    switch(state) {
    case GLUT_ENTERED:
//       cerr << "event: mouse entered\n";
      mouse.inside = true;
      break;
    case GLUT_LEFT:
//       cerr << "event: mouse left\n";
      mouse.inside = false;
      break;
    default:
      assert(0);
      break;
    }
  }

  static void dialsFunction(int const x, int const y) {
    cerr << "event: dials event x " << x << " y " << y << "\n";
  }

  static void motionFunction(int const x, int const y) {
    cerr << "event: motion event x " << x << " y " << y << "\n";
  }

  static void passiveMotionFunction(int const x, int const y) {
    cerr << "event: passive motion event x " << x << " y " << y << "\n";
  }

  static void tabletMotionFunction(int const x, int const y) {
    cerr << "event: tablet motion event x " << x << " y " << y << "\n";
  }
  static void tabletButtonFunction(int const key, int const state, 
                                   int const x, int const y) {
    cerr << "event: tablet button event key " << key << " state " << state
         << " at x " << x << " y " << y << "\n";
  }

  static void spaceballRotateFunction(int const x, int const y, int const z) {
    cerr << "event: spaceball rotate event x " << x << " y " << y << " z " << z << "\n";
  }
  static void spaceballButtonFunction(int const key, int const state) {
    cerr << "event: spaceball button event key " << key << " state " << state << "\n";
  }

  static void reshapeFunction(int const w, int const h) {
    cerr << "event: reshape event w " << w << " h " << h << "\n";
    screenWidth = w;
    screenHeight = h;
    setupOpenGL();
  }

  static void rendererCleanup(SPHT *) {
    // this only results in error:
    // "Function <glutDestroyWindow> called without first calling 'glutInit'"
    // I hope this indicates that freeglut will cleanup itself when being told to glutLeaveMainLoop
    //     glutDestroyWindow(windowId);
  }

  static void resetWindow() {
//     assert(0);
    screenWidth = 640;
    screenHeight = 480;
    glutReshapeWindow(screenWidth, screenHeight);
//     fullscreen = false;
  }

  static void resetView() {
    camera.reset();
    camera.position(initCamPos);
  }

  static int rendererFunction(SPHT *sph) {
//     int dargc = 1;
//     char *dargv[] = {0};

    int dargc = sph->program_argc;
    char **dargv = sph->program_argv;
    glutInit(&dargc, dargv);
    
    glutInitWindowSize(640, 480);
    glutInitWindowPosition(0, 0);

    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
    
    windowId = glutCreateWindow("sph");
    
    float light0_ambient[] = {0.5, 0.5, 0.5, 1};
    float light0_specular[] = {1, 0., 0., 1};
    float light0_diffuse[] = {0.5, 0.5, 0.5, 1};
    float light0_position[] = {0, 1, -10, 0};
    
    glLightfv(GL_LIGHT2, GL_SPECULAR, light0_specular);
    glLightfv(GL_LIGHT2, GL_AMBIENT,  light0_ambient);
    glLightfv(GL_LIGHT2, GL_DIFFUSE,  light0_diffuse);
    glLightfv(GL_LIGHT2, GL_POSITION, light0_position);
    
    glEnable(GL_ALPHA_TEST);
    glEnable(GL_LIGHTING);
    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_LIGHT2);
    
    glEnable(GL_DEPTH_TEST);
    
    sphPtr = sph;
    
    updateVectorLength();
    updateParticleSize();

    initCamPos = sphPtr->config.initCamPos;

    double const iraRad = sphPtr->config.initRotationAngle*2*M_PI/360.;
    initCamRotation = Quaternion(cos(iraRad), 
                                 sin(iraRad) * sphPtr->config.initRotationAxis);

//     camTarget = sphPtr->config.initCamTarget;
//     camUpVector = sphPtr->config.initCamUpVector;
//     camAxis = camTarget - camPos;


//     cerr << "camPos: " << camPos << "\n";
//     cerr << "camTarget: " << camTarget << "\n";
//     cerr << "camUp: " << camUpVector << "\n";

//     resetWindow();
    resetView();

    setupOpenGL();

    glutDisplayFunc(displayFunction);

    glutIdleFunc(0);
    // glutIdleFunc(idleFunction);

    // setup keyboard functions
    glutKeyboardFunc(keyboardFunction);
    glutSpecialFunc(specialFunction);

    // setup mouse functions
    glutEntryFunc(entryFunction);
    glutMouseFunc(mouseFunction);
    glutMotionFunc(motionFunction);
    //     glutPassiveMotionFunc(passiveMotionFunction);

    glutTabletMotionFunc(tabletMotionFunction);
    glutTabletButtonFunc(tabletButtonFunction);
    glutSpaceballRotateFunc(spaceballRotateFunction);
    glutSpaceballButtonFunc(spaceballButtonFunction);
    glutDialsFunc(dialsFunction);

    glutReshapeFunc(reshapeFunction);

    glutVisibilityFunc(visibleFunction);

    // this is need for glutLeaveMainLoop to work
    glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE,
 		  GLUT_ACTION_CONTINUE_EXECUTION);

    createMenu();

    compileDisplayLists();

#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    
    // signal initialization has completed
    rendererInitialized.signal();

    glutMainLoop();

    pthread_exit((void*)42);

    return 42;
  }
  
};


// static Renderer members
Point3D 
  Renderer::lcsysPos, 
  Renderer::textLineP1, 
  Renderer::textLineP2,
  Renderer::initCamPos;
Camera Renderer::camera;
Mouse Renderer::mouse;
Quaternion Renderer::initCamRotation;
CoordSystem const Renderer::coordSystem;
long Renderer::windowId;
size_t Renderer::frameRendered, Renderer::screenWidth = 640, Renderer::screenHeight = 480;
double const Renderer::mUnit = 0.1, Renderer::rUnit = 1;
RenderOptions Renderer::options;
SPHT *Renderer::sphPtr;
Mutex Renderer::rendererReady;
Condition Renderer::rendererInitialized;

// static Camera members
CoordSystem const Camera::coordSystem;
double const Camera::fovAngleMin = 1, Camera::fovAngleMax = 175;

#endif


/// static global variable holds info if one signal has been catched.
static int signaled;


/// static global variable holds id of thread running main().
static pthread_t mainThreadId;

/// Handler fuction for signals: Set "exit at next step" on first,
/// exit(FAILURE) on second. Do nothing if not main thread (Signals
/// are catched by all threads!). Set beende=1, then, if signalled==1
/// exit(FAILURE), then set signalled=1
extern "C" void signalHandler(int const signal) {
  if (pthread_equal(pthread_self(), mainThreadId)) {
    beende = 1;
    cerr << "signal " << signal << " caught\n";
    if (signaled) {
      cerr << "exiting because of signal\n";
      exit(EXIT_FAILURE);
    }
    signaled = 1;
  }
}


#if defined DEBUG_SPH
TSOStream<8> debugOutput(std::cerr);
#endif

#ifdef SPH_PARTICLE_MAX_VELOCITY
double maxVelocity;
double maxRelativeVelocity;
#ifdef _OPENMP
#pragma omp threadprivate(maxVelocity)
#pragma omp threadprivate(maxRelativeVelocity)
#endif
#endif

#ifdef SPH_COUNT_RELATIONS
size_t numPPrelations;
size_t numPPindextests;
size_t numPPdisttests;
#ifdef SPH_COMPILER_SUN_STUDIO
#ifdef _OPENMP
#pragma omp threadprivate(numPPrelations)
#pragma omp threadprivate(numPPindextests)
#pragma omp threadprivate(numPPdisttests)
#endif
#endif
#endif

void on_exitHandler(int status, void *) {
  cerr << "beende mit status " << status << "\n";
}

void atexitHandler() {
  cerr << "beende\n";
}

int curNumThreads = 0;

#if SPH_AD == 1
struct MyADErrorHandler : public yafad::ErrorHandler {
  yafad::ErrorHandler &subHandler;
  bool abortOnError;
  size_t numErrors[yafad::Error::YAFAD_ERROR_MAX];
  double maxDelta[yafad::Error::YAFAD_ERROR_MAX];
  double const tolerance;
  std::ostream &output;

  MyADErrorHandler(yafad::ErrorHandler &_subHandler, 
                   bool const _abortOnError,
                   double const _tolerance,
                   std::ostream &_output = std::cerr) : 
    subHandler(_subHandler),
    abortOnError(_abortOnError),
    tolerance(_tolerance),
    output(_output)
  {
    for (unsigned i = 0; i < yafad::Error::YAFAD_ERROR_MAX; ++i) {
      numErrors[i] = 0;
      maxDelta[i] = 0;
    }
  }

  void printXML(std::ostream &output) {
    for (unsigned i = 0; i < yafad::Error::YAFAD_ERROR_MAX; ++i) {
      if (numErrors[i] != 0) {
        output << "<sph-ad-error "
          " code='" << i << "'"
          " name='" << yafad::Error::getERRORSName(i) << "'"
          " count='" << numErrors[i] << "'"
          " max-delta='" << maxDelta[i] << "'"
          ">at discontinuity of function max the derivatives of both branches differ</sph-ad-error>\n";
      }
    }
  }
  
  void handle(yafad::Error const &e) {

    switch(e.code) {
    case yafad::Error::YAFAD_ERROR_DISCONT_MAX: {
      yafad::Discontinuity<adouble> const &f = dynamic_cast<yafad::Discontinuity<adouble> const &>(e);
      SphADataType::DerivVector const derivDiff = f.arguments[1]->diff() - f.arguments[0]->diff();
      double const dmax = derivDiff.max();
      if (dmax < tolerance) {
        // all derivatives equal
      } else {
#ifdef _OPENMP
#pragma omp critical
#endif
        {
          output << "discontinuity of max at args "
                 << f.arguments[0]->value() << ", " << f.arguments[1]->value() << "\n";
          output << " deriv. of left branch: " << f.arguments[0]->diff() << "\n";
          output << " deriv. of right branch: " << f.arguments[1]->diff() << "\n";
          output << "at discontinuity of max, branches have different derivatives, difference >= " << 
            maxDelta << "\n";
          ++numErrors[yafad::Error::YAFAD_ERROR_DISCONT_MAX];
          maxDelta[yafad::Error::YAFAD_ERROR_DISCONT_MAX] = 
            max(dmax, maxDelta[yafad::Error::YAFAD_ERROR_DISCONT_MAX]);
          if (abortOnError) {
            output << "aborting\n";
            beende = 1;
          }
        }
      }
      break;
    }
    default:
      cerr << "an unexpected AD error occured, aborting\n";
      exit(14);
      break;
    }

    subHandler.handle(e);
  }
};
#endif

int main(int argc, char *argv[]) {
  tzset();

  cerr.precision(12);
  cout.precision(17);

  // setup global variables
  mainThreadId = pthread_self();

  // setup signal handlers
  {
    struct sigaction sa, sao;
    sa.sa_handler = signalHandler;
    // ignoriere weitere signale waehrend signalHandler laeuft nicht
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = 0;
//     sa.sa_restorer = 0;
    sigaction(SIGINT, &sa, &sao);
    sigaction(SIGHUP, &sa, &sao);
    sigaction(SIGTERM, &sa, &sao);
    sigaction(SIGUSR1, &sa, &sao);
    sigaction(SIGUSR2, &sa, &sao);
    sigaction(SIGFPE, &sa, &sao);
  }

  SPHT::Config sphConfig;

  ofstream sphlogFile(sphConfig.logFile.c_str(), ios::binary);

  Teestream sphout(cout, sphlogFile);
  
  sphout << "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n";
  sphout << "<sphrun>\n";

#ifdef SPH_SAVE_H5PART
  HDF5 hdf5;
#endif

  SPHT sph(sphout, sphConfig);
  bool failed = false;

  // setup floating point environment
  {
#ifdef SPH_HAVE_feenableexcept
    if (sph.config.trapFloatingPointErrors) {
      feenableexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
    }
#endif
  }

#if SPH_AD == 1
  // use a pointer and new so we can force the destruction before end
  // of program (and print messages!). see corresp. delete below
  yafad::SummaryErrorHandler *subErrorHandler = new yafad::SummaryErrorHandler(false);
  MyADErrorHandler *myErrorHandler 
//      = new yafad::DefaultErrorHandler();
//     = new yafad::MesagingErrorHandler();
//     = new yafad::CountingErrorHandler();
//     = new yafad::NaNErrorHandler();
//     = new yafad::SummaryErrorHandler();
    = new MyADErrorHandler(*subErrorHandler, false, 
                           sph.config.discontMaxBranchesMaxDelta, sphout);
//   yafad::ErrorHandler::setHandler(myErrorHandler);
  yafad::ErrorHandler::setHandler(subErrorHandler);
#endif

  sphout << "<sph"
    " compiled='" __DATE__ " " __TIME__ "'"
    "/>\n";
  {
    StrFTime ftime("%Y-%m-%d %H:%M:%S %z (%c)");
    sphout << "<date>" << ftime() << "</date>\n";
  }
  sphout << "<timezone"
    " name='" << tzname[daylight ? 1 : 0] << "'"
    " daylight='" << daylight << "'"
    " timezone='" << timezone << "'/>\n";
  sphout << "<pid>" << PID() << "</pid>\n";
  sphout << Uname() << "\n";
  sph.config.writeXML(sphout);

//   {
//     FLInfo flInfo;
//     flInfo.writeXML(cout);
//   }

  sphout << "<integrator name='" << sph.integrator.name() << "'>\n"
       << sph.integrator
       << "\n</integrator>\n";

  omp_set_num_threads(1);
  omp_set_nested(false);
  assert(omp_get_nested() == false);
  omp_set_dynamic(false);
  assert(omp_get_dynamic() == false);

#ifdef _OPENMP
#if _OPENMP > 200804
  // this version supports setting schedule type at runtime
  string const &st = sph.config.openMPConfig.ompSchedule;
  int const stc = SPHOpenMPConfig::getScheduleTypesValue(st.c_str());
  sphout << "<omp-schedule type='" << st << "' code='" << stc << "' chunk-size='"
         << sph.config.openMPConfig.ompChunkSize << "'/>\n";
  omp_set_schedule(omp_sched_t(stc), sph.config.openMPConfig.ompChunkSize);
#endif
  sphout << "<omp-nested value='" << sph.config.openMPConfig.ompNested << "'/>\n";
  omp_set_nested(sph.config.openMPConfig.ompNested);
  sphout << "<omp-dynamic value='" << sph.config.openMPConfig.ompDynamic << "'/>\n";
  omp_set_dynamic(sph.config.openMPConfig.ompDynamic);
#endif

  sph.init(argc, argv);
//   assert(sph.check());

  // the particles must be set up to start the rendering thread
#ifdef SPH_RENDER
  ThreadBase *renderThread = 0;
  if (sph.config.renderEvery) {
    renderThread = makeThread(&Renderer::rendererFunction, &sph);
    renderThread->run();
  }
#endif

  Timer runTime;
  runTime.start();

  sphout << "<steps>\n";

#ifdef SPH_RENDER
  if (sph.config.renderEvery) {
    // wait for renderer to initialize completely
    Renderer::rendererReady.lock();
    Renderer::rendererInitialized.wait(Renderer::rendererReady);
  }
#endif

  for (sph.simstep = 0; ; ++sph.simstep) {
#ifdef SPH_RENDER
    if (sph.config.renderEvery and sph.simstep % sph.config.renderEvery == 0) {
      Renderer::redisplay();
    }
#endif
    if (sph.config.saveRate > 0
        and (sph.simulationTime * sph.config.saveRate) + 1 >= sph.savedFrame) {
      sph.save();
    }

    if (FLInfo::testAndShowExceptions(sphout, FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW)) {
      if (sph.config.abortOnFExcept) {
        sphout << "a floating point exception occured, aborting the program\n";
        beende = 1;
      }
    }

    if (sph.simulationTime > sph.config.tmax or beende)
      break;

    if (long(sph.partikelManager.numMovingParticles) - long(sph.partikelManager.numParticlesOut())
        < long(sph.config.minParticlesLeft)) {
      sphout << "<no-particles-left"
        " out='" << sph.partikelManager.numParticlesOut() << "'"
        " moving='" << sph.partikelManager.numMovingParticles << "'"
        " step='" << sph.simstep << "'"
        "/>\n";
      failed = true;
      break;
    }

    if (curNumThreads != sph.numThreads) {
      cerr << "set_num_threads: " << sph.numThreads << "\n";
      omp_set_num_threads(sph.numThreads);
      curNumThreads = sph.numThreads;
    }

    sph.doStep(sph.deltat);

    if (sph.config.sleepEvery
	and sph.simstep % sph.config.sleepEvery == 0) {
      sleep(sph.config.sleepInterval);
    }
  }
  beende = true;

  sph.save();

  runTime.stop();

  sphout << "\n</steps>\n";

  cerr << "terminate after " << sph.simstep
       << " steps: " << runTime << " s"
       << " s/steps: " << runTime()/sph.simstep << " s\n";

  for(unsigned i = 0; i < sph.config.numSensors; ++i) {
    Sensor<SphADataType> &sensor = sph.sensors[i];
#ifdef SPH_AD
    sphout << "<sensor-value value='" << yafad_value(sensor.m_sensorValue) << "'>\n";
    for (unsigned i = 0; i < SPH_AD_NDIR; ++i) {
      sphout << "  <diff i='" << i << "' order='1' val='" << 
#if SPH_AD == 1
        sensor.m_sensorValue.diff()[i]
#else
        sensor.m_sensorValue.imag()
#endif
             << "'/>\n";
    }
    sphout << "</sensor-value>\n";
#else
    sphout << "<sensor-value value='" << sensor.m_sensorValue << "'/>\n";
#endif
  }

  sphout << "<force-functor"
    " calls='" << sph.forceFunctor.numCalls << "'"
    " updates='" << sph.forceFunctor.numUpdates << "'"
    "/>\n";

  sphout << "\n<end steps='" << sph.simstep << "'"
    " num-particles='" << sph.partikelManager.numParticlesLive() << "'"
    " time='" << runTime() << "'"
    " time-per-step='" << runTime() / double(sph.simstep) << "'"
    " time-steps='" << sph.timeStep() << "'"
    " time-non-steps='" << sph.timeOutOfStep() << "'"
    " time-init='" << sph.timeInit() << "'"
    " time-load='" << sph.timeLoad() << "'"
    " time-save='" << sph.timeSave() << "'"
    "/>\n";

#ifdef SPH_RENDER
  if (sph.config.renderEvery) {
    // tell rendering thread to stop
    glutLeaveMainLoop();

    // join the rendering thread
    void *returnValue = 0;
    int const jc = renderThread->join(&returnValue);
    if (jc) {
      int const err = errno;
      cerr << "error: join with rendering thread failed: "
 	   << strerror(err) << "\n";
    } else {
      //       cerr << "joined thread's return value: " << long(returnValue) << "\n";
      delete renderThread;
    }

    renderThread = 0;
  }
#endif

#if SPH_AD == 1
  subErrorHandler->printXML(sphout);
  myErrorHandler->printXML(sphout);
  delete myErrorHandler;
  delete subErrorHandler;
#endif

  sphout << "</sphrun>\n";
  return failed;
}

#endif
