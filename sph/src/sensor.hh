#ifndef sph_sensor_hh
#define sph_sensor_hh
/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

#include "sensor-def.hh"

template<class T>
struct Sensor {

  typedef T value_type;
  typedef value_type (Sensor::*SensorFunction)(Partikel<value_type> const &p) const;
  typedef bool (Sensor::*SensorAreaTest)(Partikel<value_type> const &p) const;

  /// Enumeration of particle flags. Each flag must be an integer with
  /// exacly one bit set, i.e. a power of 2.
#include "sensor-flags.ncd.enum.hh"  
#include "sensor-flags.ncd.cc"

  SPHConfig<SphADataType> const &config;
  SensorDef const &m_sensorDef;
  value_type m_sensorValue;
  SensorFunction m_sensorFunction;
  SensorAreaTest m_sensorAreaTestFunction;
  double m_sensorArea;
  double m_sensorTemp;
  double m_currentTime, m_currentTimestep;
  unsigned m_index, m_triggerFlag;
  std::ostream *m_logFile;

  Sensor(SPHConfig<SphADataType> const &config, SensorDef const &sensorDef, 
         unsigned sensorNumber, std::ostream *logFile) :
    config(config),
    m_sensorDef(sensorDef),
    m_sensorValue(),
    m_sensorFunction(),
    m_sensorArea(m_sensorDef.m_gebiet.v.prod()),
    m_sensorTemp(),
    m_index(sensorNumber),
    m_triggerFlag(SENSOR_FLAGS_SENSOR_TRIGGERED << sensorNumber),
    m_logFile(logFile)
  {
    switch(m_sensorDef.m_areaType) {
    case SensorDef::SENSOR_AREA_ORTHOGONAL_RECTANGLE:
      m_sensorAreaTestFunction = &Sensor::areaTestOrtho;
      break;
    case SensorDef::SENSOR_AREA_HYPERPLANE:
      m_sensorAreaTestFunction = &Sensor::areaTestHyperplane;
      m_sensorArea = 1;
      break;
    case SensorDef::SENSOR_AREA_HYPERPLANE_ORIENTED:
      m_sensorAreaTestFunction = &Sensor::areaTestHyperplaneOriented;
      m_sensorArea = 1;
      break;
    case SensorDef::SENSOR_AREA_TWO_HYPERPLANES:
      m_sensorAreaTestFunction = &Sensor::areaTestTwoHyperplanes;
      m_sensorArea = 1;
      break;
    case SensorDef::SENSOR_AREA_SPHERE:
      m_sensorAreaTestFunction = &Sensor::areaTestSphere;
      m_sensorTemp = norm(m_sensorDef.m_gebiet.v);
      if (ndim == 1) {
        m_sensorArea = 2*m_sensorTemp;
      } else if (ndim == 2) {
        m_sensorArea = 2*m_sensorTemp*m_sensorTemp*M_PI;
      } else if (ndim == 2) {
        m_sensorArea = 3/4.*m_sensorTemp*m_sensorTemp*m_sensorTemp*M_PI;
      }
      break;
    }
    switch(m_sensorDef.m_property) {
    case SensorDef::SENSOR_TYPE_COUNT:
      m_sensorFunction = &Sensor::sensorFunctionCount;
      break;
    case SensorDef::SENSOR_TYPE_DENSITY:
      m_sensorFunction = &Sensor::sensorFunctionDensity;
      break;
    case SensorDef::SENSOR_TYPE_ENERGY_HEIGHT:
      m_sensorFunction = &Sensor::sensorFunctionEnergyHeight;
      break;
    case SensorDef::SENSOR_TYPE_IMPULSE:
      m_sensorFunction = &Sensor::sensorFunctionImpulse;
      break;
    case SensorDef::SENSOR_TYPE_KINETIC_ENERGY:
      m_sensorFunction = &Sensor::sensorFunctionKineticEnergy;
      break;
    case SensorDef::SENSOR_TYPE_MASS:
      m_sensorFunction = &Sensor::sensorFunctionMass;
      break;
    case SensorDef::SENSOR_TYPE_PRESSURE:
      m_sensorFunction = &Sensor::sensorFunctionPressure;
      break;
    case SensorDef::SENSOR_TYPE_VELOCITY:
      m_sensorFunction = &Sensor::sensorFunctionVelocity;
      break;
    default:
      cerr << "error: sensor configuration: no such property to measure: "
           << m_sensorDef.m_property << "\n";
      exit(16);
      break;
    }
    if (m_logFile) {
      *m_logFile << "<sensor-log-starts index='" << m_index << "'>\n";
      sensorDef.writeXML(*m_logFile);
      *m_logFile << "</sensor-log-starts>\n";
    }
  }

  ~Sensor() {
    if (m_logFile) {
      *m_logFile << "<sensor-log-stops index='" << m_index << "'/>\n";
    }
  }

  void logParticle(Partikel<value_type> &p, double const currentTime) {
    if (m_logFile) {
// #ifndef _OPENMP
#pragma omp critical
// #endif
      {      
        *m_logFile << "<sensor-hit index='" << m_index << "' t='" << currentTime << "'>\n";
        p.writeXML(*m_logFile);
        *m_logFile << "</sensor-hit>\n";
      }
    }
  }

  void processParticle(Partikel<value_type> &p) {
    if ((this->*m_sensorAreaTestFunction)(p)) {
      // cerr << "particle " << p._class() << " has triggered the sensor\n";
      p.sensorsTriggered() = p.sensorsTriggered() | m_triggerFlag;
      value_type const sval = (this->*m_sensorFunction)(p);
      addToSensorValue(p, sval);
      logParticle(p, m_currentTime);
    } else {
      p.sensorsTriggered() = p.sensorsTriggered() & ~m_triggerFlag;
    }
  }

  bool areaTestOrtho(Partikel<value_type> const &p) const {
    return ((p.position() >= m_sensorDef.m_gebiet.untenLinks).min() == 1
            and (p.position() < m_sensorDef.m_gebiet.obenRechts).min() == 1);
  }

  bool areaTestHyperplaneOriented(Partikel<value_type> const &p) const {
    Vector const relPos = p.position() - m_sensorDef.m_gebiet.untenLinks;
    Vector const relPos2 = p.position() + p.velocity()*m_currentTimestep - m_sensorDef.m_gebiet.untenLinks;
    // vektoren in m_gebiet2 sind normiert
    double const relPosP1 = scalar(relPos, m_sensorDef.m_gebiet2.untenLinks);
    double const relPosP2 = scalar(relPos2, m_sensorDef.m_gebiet2.untenLinks);
    return relPosP1 < 0 and relPosP2 >= 0;
  }

  bool areaTestHyperplane(Partikel<value_type> const &p) const {
    Vector const relPos = p.position() - m_sensorDef.m_gebiet.untenLinks;
    Vector const relPos2 = p.position() + p.velocity()*m_currentTimestep - m_sensorDef.m_gebiet.untenLinks;
    double const relPosP1 = scalar(relPos, m_sensorDef.m_gebiet2.untenLinks);
    double const relPosP2 = scalar(relPos2, m_sensorDef.m_gebiet2.untenLinks);
    return (relPosP1 < 0 and relPosP2 >= 0)
      or (relPosP1 >= 0 and relPosP2 < 0);
  }

  bool areaTestTwoHyperplanes(Partikel<value_type> const &p) const {
    Vector const relPos = p.position() - m_sensorDef.m_gebiet.untenLinks;
    Vector const relPos2 = p.position() - m_sensorDef.m_gebiet.obenRechts;
    double const relPos1Side = scalar(relPos, m_sensorDef.m_gebiet2.untenLinks);
    double const relPos2Side = scalar(relPos2, m_sensorDef.m_gebiet2.untenLinks);
    // cerr << "pos: " << p.position() << "\n";
    // cerr << "rel. pos 1: " << relPos << ": " << relPos1Side << "\n";
    // cerr << "rel. pos 2: " << relPos2 << ": " << relPos2Side << "\n";
    return relPos1Side >= 0 and relPos2Side < 0;
  }

  bool areaTestSphere(Partikel<value_type> const &p) const {
    Vector const relPos = p.position() - m_sensorDef.m_gebiet.untenLinks;
    return norm(relPos) < m_sensorTemp;
  }

  value_type sensorFunctionDensity(Partikel<value_type> const &p) const {
    return p.dichte();
  }

  value_type sensorFunctionEnergyHeight(Partikel<value_type> const &p) const {
    value_type const veloc2 = normSquared(p.velocity());
    return veloc2/(2.0*config.gravity)
      + p.druck()/(config.eqsConfig.rho0*config.gravity)
#if SPH_DIM == 1
      + p.position()[0]
#else
      + p.position()[1]
#endif
      ;
  }

  value_type sensorFunctionImpulse(Partikel<value_type> const &p) const {
    return norm(p.velocity()) * p.masse();
  }

  value_type sensorFunctionKineticEnergy(Partikel<value_type> const &p) const {
    return 0.5 * normSquared(p.velocity()) * p.masse();
  }

  value_type sensorFunctionMass(Partikel<value_type> const &p) const {
    return p.masse();
  }

  value_type sensorFunctionPressure(Partikel<value_type> const &p) const {
    return p.druck();
  }

  value_type sensorFunctionVelocity(Partikel<value_type> const &p) const {
    return norm(p.velocity());
  }

  value_type sensorFunctionCount(Partikel<value_type> const &) const {
    return 1;
  }

  void addToSensorValue(Partikel<value_type> const &p, value_type const &sval) {
    m_sensorValue += sval
      * (p.masse() / yafad_value(p.dichte())) / m_sensorArea
      * m_currentTimestep / config.tmax;
  }
};

#endif
