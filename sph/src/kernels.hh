#ifndef sph_kernels_hh
#define sph_kernels_hh
/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

// #include "vector.hh"
#include "spline/spline.hh"
#include "power.hh"

/*! Abstract Template Class defining the Kernel Interface
 */
template<class T, class Vector, int DIM>
struct KernelInterface {
  typedef T       value_type;
  typedef Vector vector_type;

  virtual ~KernelInterface() {}
  virtual T wAndGrad(Vector const &distance, T const dnorm, Vector &grad) const = 0;
  virtual void gradw(Vector const &distance, T const dnorm, Vector &grad) const = 0;
  virtual T w(T const dnorm) const = 0;
  static KernelInterface *makeKernel(std::string const &name, double const H);
  virtual std::string name() const = 0;
  virtual double reachLength() const = 0;
};

/*! Template defining the proportional factor of the _DIM-dimensional
  GaussianKernel.  This default version is empty to enforce an error
  if used wrongly.  Only the specializations for _DIM=1, _DIM=2 and
  _DIM=3 are to be used for 1D, 2D and 3D, resp.  \sa GaussianAlpha< 1 >,
  GaussianAlpha< 2 > and GaussianAlpha< 3 >.
 */
template<int _DIM> struct GaussianAlpha { };

/// Define the proportional factor of the 1-dimensional GaussianKernel. 	       
/// The integration over one dimension  yields \f$ \alpha_1 = \frac{\pi^{1/2} h}{2} \f$.
/// Hence the prop. factor is \f$ \alpha_1^{-1} \f$.
template<> struct GaussianAlpha<1> { 
  /// Compute the kernel propotionality factor given the smoothing length \f$ h \f$.
  /// \param h the kernel smoothing length \f$ h \f$.
  /// \returns \f$ \alpha_1^{-1}(h) \f$.
  static double value(double const h) { return 1 / (sqrt(M_PI) * h); } 
};

/// Define the proportional factor of the 2-dimensional GaussianKernel. 	       
/// The integration over two dimension  yields \f$ \alpha_2 = \pi h^2 \f$.
/// Hence the prop. factor is \f$ \alpha_2^{-1} \f$.
template<> struct GaussianAlpha<2> { 
  /// Compute the kernel propotionality factor given the smoothing length \f$ h \f$.
  /// \param h the kernel smoothing length \f$ h \f$.
  /// \returns \f$ \alpha_2^{-1}(h) \f$.
  static double value(double const h) { return 1 / (M_PI*h*h); } 
};

/// Define the proportional factor of the 3-dimensional GaussianKernel. 	       
/// The integration over three dimension  yields \f$ \alpha_3 = \pi^{3/2} h^3 \f$.
/// Hence the prop. factor is \f$ \alpha_3^{-1} \f$.
template<> struct GaussianAlpha<3> {
  /// Compute the kernel propotionality factor given the smoothing length \f$ h \f$.
  /// \param h the kernel smoothing length \f$ h \f$.
  /// \returns \f$ \alpha_3^{-1}(h) \f$.
  static double value(double const h) { return 1 / (M_PI*sqrt(M_PI)*h*h*h); }
};


/*! Defines the proportional factor of the _DIM-dimensional
  SplineKernel.  This default version is empty to enforce an error
  if used wrongly.  Only the specializations for _DIM=1, _DIM=2 and
  _DIM=3 are to be used for 1D, 2D and 3D, resp.  \sa SplineAlpha< 1 >,
  SplineAlpha< 2 > and SplineAlpha< 3 >.
 */
template<int _DIM> struct SplineAlpha { };

/// Define the proportional factor of the 1-dimensional SplineKernel. 	       
/// The integration over three dimension  yields \f$ \alpha_1 = \frac{h}{2} \f$.
/// Hence the prop. factor is \f$ \alpha_1^{-1} \f$.
template<> struct SplineAlpha<1> { 
  /// Compute the kernel propotionality factor given the smoothing length \f$ h \f$.
  /// \param h the kernel smoothing length \f$ h \f$.
  /// \returns \f$ \alpha_1^{-1}(h) \f$.
  static double value(double const h) { return 1 / h; } 
};

/// Define the proportional factor of the 2-dimensional SplineKernel. 	       
/// The integration over three dimension  yields \f$ \alpha_2 = \frac{7}{15} \pi h^2 \f$.
/// Hence the prop. factor is \f$ \alpha_2^{-1} \f$.
template<> struct SplineAlpha<2> { 
  /// Compute the kernel propotionality factor given the smoothing length \f$ h \f$.
  /// \param h the kernel smoothing length \f$ h \f$.
  /// \returns \f$ \alpha_2^{-1}(h) \f$.
  static double value(double const h) { return 15 / (7 * M_PI*h*h); } 
};

/// Define the proportional factor of the 3-dimensional SplineKernel. 	       
/// The integration over three dimension  yields \f$ \alpha_3 = \frac{2}{3} \pi h^3 \f$.
/// Hence the prop. factor is \f$ \alpha_3^{-1} \f$.
template<> struct SplineAlpha<3> { 
  /// Compute the kernel propotionality factor given the smoothing length \f$ h \f$.
  /// \param h the kernel smoothing length \f$ h \f$.
  /// \returns \f$ \alpha_3^{-1}(h) \f$.
  static double value(double const h) { return 3 / (2 * M_PI*h*h*h); }
};



/////////////////////////////////////////////////////////////////////////////////
/// Implement the gaussian kernel function.  This is given by 
/// \f$ W(r) = e ^{- \left( \frac{x}{h}\right)^2 } \f$, where \f$ r = || {\bf r} ||_2\f$ is the distance.
///  The cutoff is at 
/// \f$ 3h \f$, so the class should be given \f$ h = H/3\f$ as a construction
/// parameter. In 3D, at \f$ z = H/3 \f$ and with \f$ H = 0.01 \f$ the kernel
/// looks as follows.
/// \image html gauss.out.f.png "The 3D-Gauss Kernel at z=H/3"
/// The first component of the gradient, i.e. the derivative with respect to the first component
/// of the distance vector, \f$ \frac{dW}{dr_0} = \frac{d W}{d r} \frac{dr }{d r_0} \f$
/// is seen in the following image.
/// \image html gauss.out.gx.png "The first component of the 3D-Gauss Kernel gradient at z=H/3"
/////////////////////////////////////////////////////////////////////////////////
template<class T, class Vector, int DIM=3>
struct GaussianKernel : public KernelInterface<T,Vector,DIM> {
  typedef KernelInterface<T,Vector,DIM> Base;
  typedef typename Base::value_type value_type;
  typedef typename Base::vector_type vector_type;

  /// The class providing the proportionality constant.
  typedef GaussianAlpha<DIM> Alpha;

  ///  The smoothing length h.
  double const h;

  /// The square \f$ h^2\f$ of the smoothing length is precomputed and stored.
  double const h2;

  /// Constant holding the proportionality factor.
  double const propFactor;

  /// Constructor.
  /// \param kernelReachH the reach distance of the kernel. The smoothing 
  /// length h is H/3 since this is the Gaussian Kernel.
  GaussianKernel(double const kernelReachH) : 
    h(kernelReachH/3),
    h2(h*h),
    propFactor(Alpha::value(h)) 
  {}

  /// Compute kernel value and gradient at once.
  /// \param dist the difference vector \f$ {\bf r} \f$
  /// \param dnorm the norm \f$ r \f$ of the difference vector \f$ {\bf r} \f$
  /// \param grad output where gradient \f$ \nabla W \f$ is placed
  /// \returns kernel value \f$ W(r) \f$ 
  T wAndGrad(Vector const &dist, T const dnorm, Vector &grad) const {
    if (dnorm >= h*3) {
      grad = 0.0;
      return 0;
    }
    T const d2 = dnorm*dnorm;
    T const expexp = exp(-d2 / h2);
    T const w = propFactor * expexp;
    grad = dist;
    grad *= -(2 / h2) * w;
    return w;
  }

  /// Compute kernel value.
  /// \param dist the difference vector \f$ {\bf r} \f$
  /// \param dnorm the norm \f$ r \f$ of the difference vector \f$ {\bf r} \f$
  /// \returns kernel value \f$ W(r) \f$ 
  T w(T const dnorm) const {
    if (dnorm >= h*3) return 0;
    return propFactor * exp(-(dnorm*dnorm) / h2);
  }

  /// Compute kernel gradient.
  /// \param dist the difference vector \f$ {\bf r} \f$
  /// \param dnorm the norm \f$ r \f$ of the difference vector \f$ {\bf r} \f$
  /// \param gradw output where gradient \f$ \nabla W \f$ is placed
  /// \returns kernel gradient \f$ \nabla W \f$ 
  void gradw(Vector const &dist, T const dnorm, Vector &gradw) const {
    if (dnorm >= h*3) {
      gradw = 0;
      return;
    }
    T const d2 = dnorm*dnorm;
    gradw = dist;
    gradw *= ((-2 / h2) * (propFactor * exp(-d2 / h2)));
  }

  virtual std::string name() const { return "gauss"; }
  virtual double reachLength() const { return h; }
};

template<class T, class Vector, int DIM=3>
struct SlowGaussianKernel : public GaussianKernel<T,Vector,DIM> {
  typedef GaussianKernel<T,Vector,DIM> Base;
  typedef typename Base::value_type value_type;
  typedef typename Base::vector_type vector_type;

  /// Constructor.
  /// \param kernelReachH the reach distance of the kernel. The smoothing 
  /// length h is H/3 since this is the Gaussian Kernel.
  SlowGaussianKernel(double const kernelReachH) : 
    GaussianKernel<T,Vector,DIM>(kernelReachH)
  {}

  static size_t const M = 600;

  static T Exp(T const &x) {
    T res = 1.0, term = 1.0;
    for(size_t i = 1; i < M; ++i) {
      term *= x / T(i);
      res += term;
    }
    return res;
  }

  /// Compute kernel value and gradient at once.
  /// \param dist the difference vector \f$ {\bf r} \f$
  /// \param dnorm the norm \f$ r \f$ of the difference vector \f$ {\bf r} \f$
  /// \param grad output where gradient \f$ \nabla W \f$ is placed
  /// \returns kernel value \f$ W(r) \f$ 
  T wAndGrad(Vector const &dist, T const dnorm, Vector &grad) const {
    T const d2 = dnorm * dnorm;
    T const expexp = Exp(-d2 / this->h2);
    T const w = this->propFactor * expexp;
    grad = dist;
    grad *= -(2 / this->h2) * w;
    return w;
  }

  /// Compute kernel value.
  /// \param dist the difference vector \f$ {\bf r} \f$
  /// \param dnorm the norm \f$ r \f$ of the difference vector \f$ {\bf r} \f$
  /// \returns kernel value \f$ W(r) \f$ 
  T w(T const dnorm) const {
    return this->propFactor * Exp(-(dnorm*dnorm) / this->h2);
  }

  /// Compute kernel gradient.
  /// \param dist the difference vector \f$ {\bf r} \f$
  /// \param dnorm the norm \f$ r \f$ of the difference vector \f$ {\bf r} \f$
  /// \param gradw output where gradient \f$ \nabla W \f$ is placed
  /// \returns kernel gradient \f$ \nabla W \f$ 
  void gradw(Vector const &dist, T const dnorm, Vector &gradw) const {
    T const d2 = dnorm*dnorm;
    gradw = dist;
    gradw *= ((-2 / this->h2) * (this->propFactor * Exp(-d2 / this->h2)));
  }

  virtual std::string name() const { return "slowgauss"; }
};

/** Implement the spline-based kernel function.  This is given by 
 \f{equation} 
   W(r) = \begin{array}{ll}
      \frac{2}{3} - (\frac{r}{h})^2 + 0.5 * (\frac{r}{h})^3, &  0 \le (\frac{r}{h}) < 1, \\
      \frac{1}{6} (2 - \frac{r}{h})^3, &  1 \le (\frac{r}{h}) < 2,
     \end{array}
 \f} where \f$ r = || {\bf r} ||_2\f$ is the distance.
 The cutoff is at 
 \f$ 2h \f$, so the class should be given \f$ h = H/2\f$ as a construction
 parameter. In 3D, at \f$ z = H/3 \f$ and with \f$ H = 0.01 \f$ the kernel
 looks as follows.
 \image html spline.out.f.png "The 3D-Spline Kernel at z=H/3"
 The first component of the gradient, i.e. the derivative with respect to the first component
 of the distance vector, \f$ \frac{dW}{dr_0} = \frac{d W}{d r} \frac{dr }{d r_0} \f$
 is seen in the following image.
 \image html spline.out.gx.png "The first component of the 3D-Spline Kernel gradient at z=H/3"
*/
template<class T, class Vector, int DIM=3>
struct SplineKernel2 : public KernelInterface<T,Vector,DIM> {
  typedef KernelInterface<T,Vector,DIM> Base;
  typedef typename Base::value_type value_type;
  typedef typename Base::vector_type vector_type;

  /// The class providing the proportionality constant.
  typedef SplineAlpha<DIM> Alpha;

  /// The smoothing length.
  double const h;

  /// The square \f$ h^2\f$ of the smoothing length.
  double const h2;

  /// Constant holding the proportionality factor.
  double const propFactor;

  /// Constant holding a constant factor used in the gradient computation.
  double const gradFactor;

  /// Constructor.
  /// \param kernelReachH the reach distance of the kernel. 
  /// The smoothing length h is H/2 (spline).
  SplineKernel2(double const kernelReachH) : 
    h(kernelReachH/2),
    h2(h*h),
    propFactor(Alpha::value(h)), 
    gradFactor(propFactor / h) 
  {}

  /// Compute kernel value and gradient at once.
  /// \param distance the difference vector \f$ {\bf r} \f$
  /// \param norm the norm \f$ r \f$ of the difference vector \f$ {\bf r} \f$
  /// \param grad output where gradient \f$ \nabla W \f$ is placed
  /// \returns kernel value \f$ W(r) \f$ 
  T wAndGrad(Vector const &distance, T const dnorm, Vector &grad) const {
    //     T const d = norm(distance);
    T const R = dnorm / h;
    T res = 0, gweight = gradFactor;
    if (R < 1) {
      T const R2 = R*R;
      T const R3 = R2*R;
      res = propFactor * ((2.0/3.0) - R2 + 0.5 * R3);
      gweight *= (1.5 * R - 2.0) / h;
      grad = distance;
      grad *= gweight;
    } else if (R < 2) {
      res = propFactor * (1.0/6.0) * StaticPower<3>::value(2.0 - R);
      gweight *= -0.5 * (2.0 - R) * (2.0 - R) / dnorm;
      grad = distance;
      grad *= gweight;
    } else {
      grad = 0;
    }
//     cerr << gweight << " ";
    return res;
  }

  /// Compute kernel value.
  /// \param distance the difference vector \f$ {\bf r} \f$
  /// \returns kernel value \f$ W(r) \f$ 
  T w(T const dnorm) const {
//     T const d = norm(distance);
    T const R = dnorm / h;
    T res = 0;
    if (R < 1) {
      T const R2 = R*R;
      T const R3 = R2*R;
      res = propFactor * ((2.0/3.0) - R2 + 0.5 * R3);
    } else if (R < 2) {
      res = propFactor * (1.0/6.0) * StaticPower<3>::value(2.0 - R);
    }
    return res;
  }

  /// Compute kernel gradient.
  /// \param distance the difference vector \f$ {\bf r} \f$
  /// \returns kernel gradient \f$ \nabla W \f$ 
  void gradw(Vector const &distance, T const dnorm, Vector &grad) const {
    T const R = dnorm / h;
    T res;
    if (R < 1) {
      res = gradFactor * (1.5 * R - 2.0) / h;
    } else if (R < 2) {
      res = gradFactor * -0.5 * (2.0 - R) * (2.0 - R) / dnorm;
    } else {
      res = 0;
    }
    grad = distance;
    grad *= res;
  }

  virtual std::string name() const { return "spline"; }
  virtual double reachLength() const { return h; }
};

/*! The Lucy kernel (Lucy77) uses the same function as wendland1 in 1D,
  but in any D.
 */
template<int _DIM> struct LucyAlpha { };
template<> struct LucyAlpha<1> { 
  static double value(double const h) { return 5 / (4 * h); } 
};
template<> struct LucyAlpha<2> { 
  static double value(double const h) { return 5 / (M_PI * h*h); } 
};
template<> struct LucyAlpha<3> { 
  static double value(double const h) { return 105 / (16 * M_PI * h*h*h); } 
};

template<class T, class Vector, int DIM=3>
struct LucyKernel : public KernelInterface<T,Vector,DIM> {
  typedef KernelInterface<T,Vector,DIM> Base;
  typedef typename Base::value_type value_type;
  typedef typename Base::vector_type vector_type;

  /// The class providing the proportionality constant.
  typedef LucyAlpha<DIM> Alpha;

  double const h;
  double const propFactor;
  double const gradFactor;

  LucyKernel(double const kernelReachH) : 
    h(kernelReachH),
    propFactor(Alpha::value(h)),
    gradFactor(propFactor / h) 
  {}

  T wAndGrad(Vector const &distance, T const dnorm, Vector &grad) const {
    T const R = dnorm / h;
    T res = 0;
    if (R < 1) {
      res = propFactor * static_power<3>(1.0 - R)*(3.0*R + 1.0);
      T dres = gradFactor * (static_power<2>(1.0 - R)*((1.0 - R)*3.0 - (9.0*R + 3.0))) / dnorm;
      grad = distance;
      grad *= dres;
    } else {
      grad = 0;
    }
    return res;
  }

  T w(T const dnorm) const {
    T const R = dnorm / h;
    T res = 0;
    if (R < 1) {
      res = propFactor * static_power<3>(1.0 - R)*(3.0*R + 1.0);
    }
    return res;
  }

  void gradw(Vector const &distance, T const dnorm, Vector &grad) const {
    T const R = dnorm / h;
    if (R < 1) {
      T dres = gradFactor * (static_power<2>(1.0 - R)*((1.0 - R)*3.0 - (9.0*R + 3.0))) / dnorm;
      grad = distance;
      grad *= dres;
    } else {
      grad = 0;
    }
  }

  virtual std::string name() const { return "lucy"; }
  virtual double reachLength() const { return h; }
};


/*! The quartic Morris4 kernel (Morris1994).
 */
template<int _DIM> struct Morris4Alpha { };
template<> struct Morris4Alpha<1> { 
  static double value(double const h) { return 1 / (24 * h); } 
};
template<> struct Morris4Alpha<2> { 
  static double value(double const h) { return 96 / (1199 * M_PI * h*h); } 
};
template<> struct Morris4Alpha<3> { 
  static double value(double const h) { return 1 / (20 * M_PI * h*h*h); } 
};

template<class T, class Vector, int DIM=3>
struct Morris4Kernel : public KernelInterface<T,Vector,DIM> {
  typedef KernelInterface<T,Vector,DIM> Base;
  typedef typename Base::value_type value_type;
  typedef typename Base::vector_type vector_type;

  /// The class providing the proportionality constant.
  typedef Morris4Alpha<DIM> Alpha;

  double const h;
  double const propFactor;
  double const gradFactor;

  Morris4Kernel(double const kernelReachH) : 
    h(kernelReachH/2.5),
    propFactor(Alpha::value(h)),
    gradFactor(propFactor / h) 
  {}

  T w(T const dnorm) const {
    T const R = dnorm / h;
    T res = 0;
    if (R < 0.5) {
      res = propFactor * (static_power<4>(R+2.5) - 5.*static_power<4>(R+1.5) + 10.*static_power<4>(R+0.5));
    } else if (R < 1.5) {
      res = propFactor * (static_power<4>(2.5-R) - 5.*static_power<4>(1.5-R));
    } else if (R < 2.5) {
      res = propFactor * static_power<4>(2.5-R);
    }
    return res;
  }
  
  void gradw(Vector const &distance, T const dnorm, Vector &grad) const {
    T const R = dnorm / h;
    if (R < 2.5) {
      T dres = 0;
      if (R < 0.5) {
        dres = 4.*static_power<3>(R+2.5) - 20.*static_power<3>(R+1.5) + 40.*static_power<3>(R+0.5);
      } else if (R < 1.5) {
        dres = -4.*static_power<3>(2.5-R) + 20.*static_power<3>(1.5-R);
      } else {
        dres = - 4.*static_power<3>(2.5-R);
      }
      grad = distance;
      grad *= gradFactor * dres / dnorm;
    } else {
      grad = 0;
    }
  }

  T wAndGrad(Vector const &distance, T const dnorm, Vector &grad) const {
    T const R = dnorm / h;
    T res = 0;
    if (R < 2.5) {
      T dres = 0;
      if (R < 0.5) {
        res = propFactor * (static_power<4>(R+2.5) - 5.*static_power<4>(R+1.5) + 10.*static_power<4>(R+0.5));
        dres = 4.*static_power<3>(R+2.5) - 20.*static_power<3>(R+1.5) + 40.*static_power<3>(R+0.5);
      } else if (R < 1.5) {
        res = propFactor * (static_power<4>(2.5-R) - 5.*static_power<4>(1.5-R));
        dres = -4.*static_power<3>(2.5-R) + 20.*static_power<3>(1.5-R);
      } else {
        res = propFactor * static_power<4>(2.5-R);
        dres = -4. * static_power<3>(2.5-R);
      }
      grad = distance;
      grad *= gradFactor * dres / dnorm;
    } else {
      grad = 0;
    }
    return res;
  }

  virtual std::string name() const { return "morris4"; }
  virtual double reachLength() const { return h; }
};


/*! The quintic Morris5 kernel (Morris1996).
 */
template<int _DIM> struct Morris5Alpha { };
template<> struct Morris5Alpha<1> { 
  static double value(double const h) { return 1 / (120 * h); } 
};
template<> struct Morris5Alpha<2> { 
  static double value(double const h) { return 7 / (478 * M_PI * h*h); } 
};
template<> struct Morris5Alpha<3> { 
  static double value(double const h) { return 1 / (120 * M_PI * h*h*h); } 
};

template<class T, class Vector, int DIM=3>
struct Morris5Kernel : public KernelInterface<T,Vector,DIM> {
  typedef KernelInterface<T,Vector,DIM> Base;
  typedef typename Base::value_type value_type;
  typedef typename Base::vector_type vector_type;

  /// The class providing the proportionality constant.
  typedef Morris5Alpha<DIM> Alpha;

  double const h;
  double const propFactor;
  double const gradFactor;

  Morris5Kernel(double const kernelReachH) : 
    h(kernelReachH/3),
    propFactor(Alpha::value(h)),
    gradFactor(propFactor / h) 
  {}

  T w(T const dnorm) const {
    T const R = dnorm / h;
    T res = 0;
    if (R < 1) {
      res = propFactor * (static_power<5>(3. - R) - 6.*static_power<5>(2. - R) + 15.*static_power<5>(1. - R));
    } else if (R < 2) {
      res = propFactor * (static_power<5>(3. - R) - 6.*static_power<5>(2. - R));
    } else if (R < 3) {
      res = propFactor * static_power<5>(3. - R);
    }
    return res;
  }
  
  void gradw(Vector const &distance, T const dnorm, Vector &grad) const {
    T const R = dnorm / h;
    if (R < 3) {
      T dres = 0;
      if (R < 1) {
        dres = -5.*static_power<4>(3. - R) + 30.*static_power<4>(2. - R) - 75.*static_power<4>(1. - R);
      } else if (R < 2) {
        dres = -5.*static_power<4>(3. - R) + 30.*static_power<4>(2. - R);
      } else {
        dres = -5.*static_power<4>(3. - R);
      }
      grad = distance;
      grad *= gradFactor * dres / dnorm;
    } else {
      grad = 0;
    }
  }

  T wAndGrad(Vector const &distance, T const dnorm, Vector &grad) const {
    T const R = dnorm / h;
    T res = 0;
    if (R < 3) {
      T dres = 0;
      if (R < 1) {
        res = propFactor * (static_power<5>(3. - R) - 6.*static_power<5>(2. - R) + 15.*static_power<5>(1. - R));
        dres = -5.*static_power<4>(3. - R) + 30.*static_power<4>(2. - R) - 75.*static_power<4>(1. - R);
      } else if (R < 2) {
        res = propFactor * (static_power<5>(3. - R) - 6.*static_power<5>(2. - R));
        dres = -5.*static_power<4>(3. - R) + 30.*static_power<4>(2. - R);
      } else {
        res = propFactor * static_power<5>(3. - R);
        dres = -5.*static_power<4>(3. - R);
      }
      grad = distance;
      grad *= gradFactor * dres / dnorm;
    } else {
      grad = 0;
    }
    return res;
  }

  virtual std::string name() const { return "morris5"; }
  virtual double reachLength() const { return h; }
};


/*! The quadratic Johnson kernel (Johnson1996).
 */
template<int _DIM> struct JohnsonAlpha { };
template<> struct JohnsonAlpha<1> { 
  static double value(double const h) { return 1 / h; } 
};
template<> struct JohnsonAlpha<2> { 
  static double value(double const h) { return 2 / (M_PI * h*h); } 
};
template<> struct JohnsonAlpha<3> { 
  static double value(double const h) { return 5 / (4 * M_PI * h*h*h); } 
};

template<class T, class Vector, int DIM=3>
struct JohnsonKernel : public KernelInterface<T,Vector,DIM> {
  typedef KernelInterface<T,Vector,DIM> Base;
  typedef typename Base::value_type value_type;
  typedef typename Base::vector_type vector_type;

  /// The class providing the proportionality constant.
  typedef JohnsonAlpha<DIM> Alpha;

  double const h;
  double const propFactor;
  double const gradFactor;

  JohnsonKernel(double const kernelReachH) : 
    h(kernelReachH/2),
    propFactor(Alpha::value(h)),
    gradFactor(propFactor / h) 
  {}

  T wAndGrad(Vector const &distance, T const dnorm, Vector &grad) const {
    T const R = dnorm / h;
    T res = 0;
    if (R < 2) {
      res = propFactor * (3. * R * R / 16. - 3. * R / 4. + 3.0 / 4.);
      T dres = gradFactor * (6. * R / 16. - 3.0 / 4.) / dnorm;
      grad = distance;
      grad *= dres;
    } else {
      grad = 0;
    }
    return res;
  }

  T w(T const dnorm) const {
    T const R = dnorm / h;
    T res = 0;
    if (R < 2) {
      res = propFactor * (3. * R * R / 16. - 3. * R / 4. + 3.0 / 4.);
    }
    return res;
  }

  void gradw(Vector const &distance, T const dnorm, Vector &grad) const {
    T const R = dnorm / h;
    if (R < 2) {
      T dres = gradFactor * (6. * R / 16. - 3.0 / 4.) / dnorm;
      grad = distance;
      grad *= dres;
    } else {
      grad = 0;
    }
  }

  virtual std::string name() const { return "johnson"; }
  virtual double reachLength() const { return h; }
};


template<int D, int K>
struct WendlandAlpha {
};


template<>
struct WendlandAlpha<1,0> {
  static double value(double const h) { return 1.0/h; }
};
template<>
struct WendlandAlpha<1,1> {
  static double value(double const h) { return 5.0/(4.0 * h); }
};
template<>
struct WendlandAlpha<1,2> {
  static double value(double const h) { return 3.0/(2*h); }
};

template<>
struct WendlandAlpha<2,0> {
  static double value(double const h) { return 6.0/(M_PI * h*h); }
};
template<>
struct WendlandAlpha<2,1> {
  static double value(double const h) { return 7.0/(M_PI * h*h); }
};
template<>
struct WendlandAlpha<2,2> {
  static double value(double const h) { return 3.0/(M_PI * h*h); }
};
template<>
struct WendlandAlpha<2,3> {
  static double value(double const h) { return 78.0/(7.0 * M_PI * h*h); }
};


template<>
struct WendlandAlpha<3,0> {
  static double value(double const h) { return 15.0/(2.0 * M_PI * h*h*h); }
};
template<>
struct WendlandAlpha<3,1> {
  static double value(double const h) { return 21.0/(2.0 * M_PI * h*h*h); }
};
template<>
struct WendlandAlpha<3,2> {
  static double value(double const h) { return 165.0/(32.0 * M_PI * h*h*h); }
};
template<>
struct WendlandAlpha<3,3> {
  static double value(double const h) { return 1365.0/(64.0 * M_PI * h*h*h); }
};

template<class T, int D, int K>
struct WendlandFunction { };

template<class T>
struct WendlandFunction<T, 1, 0> {
  T operator ()(T const r) const {
    if (r > 1) return 0;
    return (1.0 - r);
  }
  T deriv(T const r) const {
    if (r > 1) return 0;
    return -1;
  }
};

template<class T>
struct WendlandFunction<T, 1, 1> {
  T operator ()(T const r) const {
    if (r > 1) return 0;
    return static_power<3>(1.0 - r)*(3.0*r + 1.0);
  }
  T deriv(T const r) const {
    if (r > 1) return 0;
    return static_power<2>(1.0 - r)*((1.0 - r)*3.0 - (9.0*r + 3.0));
  }
};

template<class T>
struct WendlandFunction<T, 1, 2> {
  T operator ()(T const r) const {
    if (r > 1) return 0;
    return static_power<5>(1.0 - r)*(8.0*r*r + 5.0*r + 1.0);
  }
  T deriv(T const r) const {
    if (r > 1) return 0;
    return static_power<4>(1.0 - r)*((1.0 - r)*(16.0*r + 5.0) - 5.0*(8.0*r*r + 5.0*r + 1.0));
  }
};

template<class T>
struct WendlandFunction<T, 3, 0> {
  T operator ()(T const r) const {
    if (r > 1) return 0;
    return (1.0 - r)*(1.0 - r);
  }
  T deriv(T const r) const {
    if (r > 1) return 0;
    return -2.0*(1.0 - r);
  }
};
template<class T>
struct WendlandFunction<T, 3, 1> {
  T operator ()(T const r) const {
    if (r > 1) return 0;
    return static_power<4>(1.0 - r)*(4.0*r + 1.0);
  }
  T deriv(T const r) const {
    if (r > 1) return 0;
    return static_power<3>(1.0 - r)*(4.0*(1.0 - r) - 4.0*(4.0*r + 1.0));
  }
};
template<class T>
struct WendlandFunction<T, 3, 2> {
  T operator ()(T const r) const {
    if (r > 1) return 0;
    return static_power<6>(1.0 - r)*(35.0*r*r + 18.0*r + 3.0);
  }
  T deriv(T const r) const {
    if (r > 1) return 0;
    return static_power<5>(1.0 - r)
      *( (1.0 - r)*(70.0*r + 18.0) - 6.0*(35.0*r*r + 18.0*r + 3.0) );
  }
};
template<class T>
struct WendlandFunction<T, 3, 3> {
  T operator ()(T const r) const {
    if (r > 1) return 0;
    return static_power<8>(1.0 - r)*(32.0*r*r*r + 25.0*r*r + 8.0*r + 1.0);
  }
  T deriv(T const r) const {
    if (r > 1) return 0;
    return static_power<7>(1.0 - r)
      *( (1.0 - r)*(96.0*r*r + 50.0*r + 8.0) - 8.0*(32.0*r*r*r + 25.0*r*r + 8.0*r + 1.0) );
  }
};

/// WendlandFunction<T, 2, K> falls back on 3D version
template<class T, int K>
struct WendlandFunction<T, 2, K> {
  WendlandFunction<T, 3, K> wf3d;
  T operator ()(T const r) const {
    return wf3d(r);
  }
  T deriv(T const r) const {
    return wf3d.deriv(r);
  }
};


// template<class T>
// T wendlandFunction<1, 1>(T const r) {
//   if (r > 1) return 0;
//   return static_power<3>(1 - r)*(3*r + 1);
// }

// template<class T>
// T wendlandFunction<1, 2>(T const r) {
//   if (r > 1) return 0;
//   return static_power<5>(1 - r)*(8*r*r + 5*r + 1);
// }

template<int WK>
struct WendlandName {};

template<>
struct WendlandName<0> {
  static std::string name() { return "wendland1"; }
};
template<>
struct WendlandName<1> {
  static std::string name() { return "wendland1"; }
};
template<>
struct WendlandName<2> {
  static std::string name() { return "wendland2"; }
};
template<>
struct WendlandName<3> {
  static std::string name() { return "wendland3"; }
};

template<class T, class Vector, int DIM=3, int WK=1>
struct WendlandKernel : public KernelInterface<T,Vector,DIM> {
  typedef KernelInterface<T,Vector,DIM> Base;
  typedef typename Base::value_type value_type;
  typedef typename Base::vector_type vector_type;

  typedef WendlandFunction<T,DIM,WK> WF;
  WF const wendlandFunction;

  /// The smoothing length.
  double const h;

  /// The square \f$ h^2\f$ of the smoothing length.
  double const h2;

  /// Constant holding the proportionality factor.
  double const propFactor;

  /// Constant holding a constant factor used in the gradient computation.
  double const gradFactor;

  /// Constructor.
  /// \param kernelReachH the reach distance of the kernel. 
  /// The smoothing length h is H/2 (spline).
  WendlandKernel(double const kernelReachH) : 
    wendlandFunction(),
    h(kernelReachH),
    h2(h*h),
    propFactor(WendlandAlpha<DIM, WK>::value(h)), 
    gradFactor(propFactor / h)
  {}

  /// Compute kernel value.
  /// \param distance the difference vector \f$ {\bf r} \f$
  /// \returns kernel value \f$ W(r) \f$ 
  T w(T const dnorm) const {
    T const R = dnorm / h;
    T res = propFactor * wendlandFunction(R);
    return res;
  }

  /// Compute kernel gradient.
  /// \param distance the difference vector \f$ {\bf r} \f$
  /// \returns kernel gradient \f$ \nabla W \f$ 
  void gradw(Vector const &distance, T const dnorm, Vector &grad) const {
    T const R = dnorm / h;
    T res = propFactor * wendlandFunction.deriv(R) / (h*dnorm);
    grad = distance;
    grad *= res;
  }

  /// Compute kernel value and gradient at once.
  /// \param distance the difference vector \f$ {\bf r} \f$
  /// \param norm the norm \f$ r \f$ of the difference vector \f$ {\bf r} \f$
  /// \param grad output where gradient \f$ \nabla W \f$ is placed
  /// \returns kernel value \f$ W(r) \f$ 
  T wAndGrad(Vector const &distance, T const dnorm, Vector &grad) const {
    T const R = dnorm / h;
    T res = propFactor * wendlandFunction(R);
    T dres = propFactor * wendlandFunction.deriv(R) / (h*dnorm);
    grad = distance;
    grad *= dres;
    return res;
  }

  std::string name() const { return WendlandName<WK>::name(); }
  double reachLength() const { return h; }
};


template<class T, class Vector, int DIM>
KernelInterface<T, Vector, DIM> *
KernelInterface<T, Vector, DIM>::makeKernel(std::string const &kernelName, double const H) {
  KernelInterface<T, Vector, DIM> *kernel = 0;
  if (kernelName.compare("spline") == 0) {
    kernel = new SplineKernel2<T, Vector, DIM>(H);

  } else if (kernelName.compare("gauss") == 0) {
    kernel = new GaussianKernel<T, Vector, DIM>(H);

  } else if (kernelName.compare("lucy") == 0) {
    kernel = new LucyKernel<T, Vector, DIM>(H);

  } else if (kernelName.compare("morris4") == 0) {
    kernel = new Morris4Kernel<T, Vector, DIM>(H);

  } else if (kernelName.compare("morris5") == 0) {
    kernel = new Morris5Kernel<T, Vector, DIM>(H);

  } else if (kernelName.compare("johnson") == 0) {
    kernel = new JohnsonKernel<T, Vector, DIM>(H);

  } else if (kernelName.compare("slowgauss") == 0) {
    kernel = new SlowGaussianKernel<T, Vector, DIM>(H);

  } else if (kernelName.compare("wendland0") == 0) {
    kernel = new WendlandKernel<T, Vector, DIM, 0>(H);

  } else if (kernelName.compare("wendland1") == 0) {
    kernel = new WendlandKernel<T, Vector, DIM, 1>(H);

  } else if (kernelName.compare("wendland2") == 0) {
    kernel = new WendlandKernel<T, Vector, DIM, 2>(H);

  } else if (kernelName.compare("wendland3") == 0) {
#if !defined SPH_DIM
#error please define SPH_DIM
#endif
#if SPH_DIM == 1
    std::cerr <<
      "error: wendland3 is not available in 1D version, sorry\n";
#else 
    kernel = new WendlandKernel<T, Vector, DIM, 3>(H);
#endif
  } else {
    std::cerr << "error: unknown kernel name `" << kernelName << "'\n";
  }
  return kernel;
}

#endif
