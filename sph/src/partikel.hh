#ifndef sph_partikel_hh
#define sph_partikel_hh
/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/
#include <valarray>
#include <string.h>
#include "vector.hh"

template<class T>
struct PartikelBuffer {
  
};

template<class T> class Partikel;

/// Represents the differential variables of a particle: density, heat, position and velocity.
/// Position and velocity are \f$ D\f$-dim. vector quantities.
template<class T>
class PartikelVariablen {
  typedef T value_type;

public:
#ifdef SPH_PARTICLE_HAS_ID
  unsigned m_id;
#endif
#ifdef SPH_PARTICLE_VAR_HAS_ID_SUM
  unsigned m_sumOfIds;
#endif
#ifdef SPH_PARTICLE_VAR_HAS_NEIGHBOUR_SUM
  unsigned m_numNeighbours;
#endif

private:

  /// The density.
  value_type m_dichte;
  /// The heat.
  value_type m_waerme;
  /// The position vector.
  Vector m_pos;
  /// The velocity vector.
  Vector m_vel;

public:

  /// The constructor ensures initialization to 0 of all members.
  PartikelVariablen() :
#ifdef SPH_PARTICLE_HAS_ID
    m_id(),
#endif
#ifdef SPH_PARTICLE_VAR_HAS_ID_SUM
    m_sumOfIds(),
#endif
#ifdef SPH_PARTICLE_VAR_HAS_NEIGHBOUR_SUM
    m_numNeighbours(),
#endif
    m_dichte(),
    m_waerme(),
    m_pos(ndim),
    m_vel(ndim)
  {}

  /// The constructor initializes all members with a given value.
  /// \param v the value
  explicit PartikelVariablen(value_type const v) :
#ifdef SPH_PARTICLE_HAS_ID
    m_id(),
#endif
#ifdef SPH_PARTICLE_VAR_HAS_ID_SUM
    m_sumOfIds(),
#endif
#ifdef SPH_PARTICLE_VAR_HAS_NEIGHBOUR_SUM
    m_numNeighbours(),
#endif
    m_dichte(v),
    m_waerme(v),
    m_pos(v, ndim),
    m_vel(v, ndim)
  {}

  /// Addition with another object of this type.
  /// \param c another PartikelVariablen object.
  /// \returns a reference to *this.
  PartikelVariablen &operator +=(PartikelVariablen const &c) {
    m_dichte += c.m_dichte;
    m_waerme += c.m_waerme;
    m_pos += c.m_pos;
    m_vel += c.m_vel;
#ifdef SPH_PARTICLE_HAS_ID
    m_id  = max(m_id, c.m_id);
#endif
#ifdef SPH_PARTICLE_VAR_HAS_ID_SUM
    m_sumOfIds += c.m_sumOfIds;
#endif
#ifdef SPH_PARTICLE_VAR_HAS_NEIGHBOUR_SUM
    m_numNeighbours += c.m_numNeighbours;
#endif
    return *this;
  }

  /// Subtraction of another object of this type.
  /// \param c another PartikelVariablen object.
  /// \returns a reference to *this.
  PartikelVariablen &operator -=(PartikelVariablen const &c) {
    m_dichte -= c.m_dichte;
    m_waerme -= c.m_waerme;
    m_pos -= c.m_pos;
    m_vel -= c.m_vel;
    return *this;
  }

  /// Scaling by a real constant.
  /// \param c an value_type constant.
  /// \returns a reference to *this.
  PartikelVariablen &operator *=(value_type const c) {
    m_dichte *= c;
    m_waerme *= c;
    m_pos *= c;
    m_vel *= c;
    return *this;
  }

  /// Assigning a real constant.
  /// \param c an value_type constant.
  /// \returns a reference to *this.
  PartikelVariablen &operator =(value_type const c) {
    m_dichte = c;
    m_waerme = c;
    m_pos = c;
    m_vel = c;
    return *this;
  }

  /// Scaling by a real constant.
  /// \param c an value_type constant.
  /// \returns a reference to *this.
  PartikelVariablen &operator *=(PartikelVariablen const &c) {
    m_dichte *= c.m_dichte;
    m_waerme *= c.m_waerme;
    m_pos *= c.m_pos;
    m_vel *= c.m_vel;
    return *this;
  }

  /// Assignment of a Partikel object will copy over it's PartikelVariablen member.
  /// \deprecated remove this operator for creater clarity at point of invocation.
  /// class Partikel contains a member of type PartikelVariablen, so class Partikel
  /// has to be forward declared here.
  /// \param v a reference to a particle object.
  /// \returns a reference to *this.
  PartikelVariablen &operator =(Partikel<T> const &v);

#ifndef NDEBUG
  /// Test for equality of all members.
  /// Only needed by debugging code. Only present if NDEBUG not defined.
  /// \param o another PartikelVariablen object.
  /// \returns true if all members equal, else false.
  bool operator ==(PartikelVariablen const &o) const {
    return (m_dichte == o.m_dichte) and
      (m_waerme == o.m_waerme) and
      (m_pos == o.m_pos).min() == 1 and
      (m_vel == o.m_vel).min() == 1;
  }
#endif

#ifdef SPH_PARTICLE_HAS_ID
  unsigned        &id()       { return m_id; }
  unsigned  const &id() const { return m_id; }
#endif

#ifdef SPH_PARTICLE_VAR_HAS_NEIGHBOUR_SUM
  unsigned        &numNeighbours()       { return m_numNeighbours; }
  unsigned  const &numNeighbours() const { return m_numNeighbours; }
#endif

#ifdef SPH_PARTICLE_VAR_HAS_ID_SUM
  unsigned        &sumOfNeighbourIds()       { return m_sumOfIds; }
  unsigned  const &sumOfNeighbourIds() const { return m_sumOfIds; }
#endif

  /// Non-const accessor method returns a reference to density variable.
  value_type &dichte()    { return m_dichte; }

  /// Non-const accessor method returns a reference to heat variable.
  value_type &waerme()    { return m_waerme; }

  /// Constant accessor method returns a copy of density variable.
  value_type const &dichte()    const { return m_dichte; }

  /// Constant accessor method returns a copy of heat variable.
  value_type const &waerme()    const { return m_waerme; }


  /// Non-const accessor method returns a reference to position variable.
  Vector       &position()       { return m_pos; }
  /// Constant accessor method returns a const-reference to position variable.
  Vector const &position() const { return m_pos; }

  /// Non-const accessor method returns a reference to velocity variable.
  Vector       &velocity()       { return m_vel; }
  /// Constant accessor method returns a const-reference to velocity variable.
  Vector const &velocity() const { return m_vel; }

  /// Print the variable values to an std::ostream object.
  void print(std::ostream &aus) const {
    aus << "p='" << position()
	<< "' v='" << velocity()
	<< "' d='" << dichte()
	<< "' w='" << waerme() << "'";
  }

  /// Print the variable values to an std::ostream object.
  void writeXML(std::ostream &aus) const {
    aus << "<pv p='" << position()
	<< "' v='" << velocity()
	<< "' d='" << dichte()
	<< "' w='" << waerme() << "'/>";
  }

  bool hasNaN() const {
#ifdef SPH_AD
#if SPH_AD == 1
    return isnan(dichte()) or isnan(waerme()) 
      or isnan(velocity().sum()) or isnan(position().sum())
      or isnan(dichte().diff().sum()) or isnan(waerme().diff().sum()) 
      or isnan(velocity().sum().diff().sum()) or isnan(position().sum().diff().sum());
#else
    return isnan(dichte().real()) or isnan(waerme().real()) 
      or isnan(velocity().sum().real()) or isnan(position().sum().real())
      or isnan(dichte().imag()) or isnan(waerme().imag()) 
      or isnan(velocity().sum().imag()) or isnan(position().sum().imag());
#endif
#else
    return isnan(dichte()) or isnan(waerme()) 
      or isnan(velocity().sum()) or isnan(position().sum());
#endif
  }

  bool isZero() const {
    return dichte() == 0 and waerme() == 0
      and velocity().sum() == 0 and position().sum() == 0;
  }

};

template<class T>
inline PartikelVariablen<T> operator *(PartikelVariablen<T> const &p, PartikelVariablen<T> const &q) {
  PartikelVariablen<T> res(p);
  res *= q;
  return res;
}

/// Print PartikelVariablen object to an std::ostream object.
template<class T>
inline std::ostream &operator <<(std::ostream &aus, PartikelVariablen<T> const &p) {
  p.print(aus);
  return aus;
}

// Implementation of PartikelVariablen object to an std::ostream object.
// inline PartikelVariablen operator *(value_type a, PartikelVariablen const &v) {
//   PartikelVariablen r(v);
//   r *= a;
//   return r;
// }

/// Represents a SPH particle. It uses a PartikelVariablen object to
/// hold the differential quantities and also stores the constant and
/// dependant variables of a particle.
template<class T>
struct Partikel {
  typedef T value_type;

  /// Type used for particle id: unsigned, there will not be more than 4G particles.
  typedef unsigned Id;

#ifdef SPH_RENDER
#include "show-vectors.ncd.enum.hh"  
#include "show-vectors.ncd.cc"
#endif

  /// Enumeration of allowed/implemented particle classes.
#include "particle-classes.ncd.enum.hh"  
#include "particle-classes.ncd.cc"

  /// Enumeration of particle flags. Each flag must be an integer with
  /// exacly one bit set, i.e. a power of 2.
#include "particle-flags.ncd.enum.hh"  
#include "particle-flags.ncd.cc"

private:
  /// Holds the differential quantities.
  PartikelVariablen<value_type> m_werte;

  /// The particle mass is constant per particle.
  double m_masse;

  /// The particle pressure depends on the density via the eiuation of state (EQS).
  value_type m_druck;

  /// The \f$ D\f$-dimensional integer raster index of the
  /// particle. This is updated in function
  /// PartikelManager::updatePartikelIndex.
  // \ todo really necessary to save for each particle?
#ifdef SPH_PARTICLE_SAVES_RASTER_INDEX
  Coord         m_index;
#endif

  /// The one-dimensional integer raster index of the particle. This
  /// is updated in function PartikelManager::updatePartikelIndex.
  unsigned indexHash;

  /// The particle flags a a bitwise \b or of the values Partikel::ParticleFlags.
  unsigned char m_flags;

  /// a bitwise \b or of the values Partikel::ParticleFlags.
  unsigned short m_sensorsTriggered;

  /// The class of a particle is an integer in the range [0,
  /// 127]. Currently only the values in Partikel::ParticleClasses are
  /// allowed and implemented.
  unsigned char m_class;

  /// The color of a particle is an arbitrary user-defined integer.
  unsigned short m_color;

public:
#ifdef SPH_PARTICLE_HAS_ID
  unsigned        &id()       { return m_werte.id(); }
  unsigned  const &id() const { return m_werte.id(); }
#endif

#ifdef SPH_PARTICLE_VAR_HAS_NEIGHBOUR_SUM
  unsigned        &numNeighbours()       { return m_werte.numNeighbours(); }
  unsigned  const &numNeighbours() const { return m_werte.numNeighbours(); }
#endif

#ifdef SPH_PARTICLE_VAR_HAS_ID_SUM
  unsigned        &sumOfNeighbourIds()       { return m_werte.sumOfNeighbourIds(); }
  unsigned  const &sumOfNeighbourIds() const { return m_werte.sumOfNeighbourIds(); }
#endif

#ifdef SPH_PARTICLE_SAVES_RASTER_INDEX
  Coord       &index()       { return m_index; }
  Coord const &index() const { return m_index; }
#endif

  unsigned &hashVal()       { return indexHash; }
  unsigned  const &hashVal() const { return indexHash; }

  unsigned char &_class()       { return m_class; }
  unsigned char  const &_class() const { return m_class; }

  unsigned char &flags()        { return m_flags; }
  unsigned char  const &flags()  const { return m_flags; }

  unsigned short &sensorsTriggered()               { return m_sensorsTriggered; }
  unsigned short  const &sensorsTriggered()  const { return m_sensorsTriggered; }

  unsigned short &color()        { return m_color; }
  unsigned short  const &color()  const { return m_color; }

  bool isOut()      const { return flags() & PARTICLE_FLAG_OUT; }
  bool isBoundary() const { return flags() & PARTICLE_FLAG_BOUNDARY; }
  bool isRenderOn() const { return flags() & PARTICLE_FLAG_RENDER; }
//   bool isSaveOn()   const { return flags() & PARTICLE_FLAG_SAVE; }

  PartikelVariablen<T>       &werte()       { return m_werte; }
  PartikelVariablen<T> const &werte() const { return m_werte; }

  double     &masse()     { return m_masse; }
  value_type &dichte()    { return m_werte.dichte(); }
  value_type &druck()     { return m_druck; }
  value_type &waerme()    { return m_werte.waerme(); }

  double     const &masse()     const { return m_masse; }
  value_type const &dichte()    const { return m_werte.dichte(); }
  value_type const &druck()     const { return m_druck; }
  value_type const &waerme()    const { return m_werte.waerme(); }


  Vector       &position()       { return m_werte.position(); }
  Vector const &position() const { return m_werte.position(); }

  Vector       &velocity()       { return m_werte.velocity(); }
  Vector const &velocity() const { return m_werte.velocity(); }

  bool hasNaN() const {
    return m_werte.hasNaN() or isnan(masse()) or isnan(druck());
  }

  Partikel() : 
    m_werte(),
    m_masse(),
    m_druck(),
    indexHash(),
    m_flags(), 
    m_class(),
    m_color()
  {
    _class() = PARTICLE_TYPE_MOVING_DEFAULT;
  }

  void print(std::ostream &aus) const {
    m_werte.print(aus);
    aus
#ifdef SPH_PARTICLE_HAS_ID
      << " id " << id()
#endif
      << " bin " << hashVal() << " ind " << index() << "\n";
    aus << " m: " << masse()
	<< " p: " << druck()
	<< " cl: " << int(_class())
	<< " col: " << int(color())
	<< " fl: " << std::hex << int(flags()) << std::dec << "\n";
  }

  void writeXML(std::ostream &aus) const {
    aus << "<p x='" << position() << "'"
	<< " v='" << velocity() << "'"
	<< " d='" << dichte() << "'"
	<< " p='" << druck() << "'"
	<< " h='" << waerme() << "'"
	<< " c='" << color() << "'"
	<< " cl='" << unsigned(_class()) << "'"
	<< " fl='" << unsigned(flags()) << "'"
	<< " s='" << sensorsTriggered() << "'"
        << "/>";
  }

};
template<class T>
inline std::ostream &operator <<(std::ostream &aus, Partikel<T> const &p) {
  p.print(aus);
  return aus;
}

template<class T>
inline PartikelVariablen<T> &PartikelVariablen<T>::operator =(Partikel<T> const &v) {
  *this = v.werte();
  return *this;
}

#endif
