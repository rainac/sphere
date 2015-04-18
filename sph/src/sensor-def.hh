#ifndef sph_sensordef_hh
#define sph_sensordef_hh
/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

#include <string.h>

struct SensorDef {
  typedef Quader<DoubleVector> Gebiet;
  Gebiet m_gebiet;
  Gebiet m_gebiet2;
  int m_areaType;
  int m_property;
  bool m_log;

  SensorDef() :
    m_gebiet(),
    m_property(),
    m_log()
  {}

  SensorDef(Gebiet const &gebiet, int property, Gebiet const &gebiet2 = Gebiet()) :
    m_gebiet(gebiet),
    m_gebiet2(gebiet2),
    m_areaType(SENSOR_AREA_ORTHOGONAL_RECTANGLE),
    m_property(property),
    m_log()
  {
    m_gebiet2.untenLinks /= norm(m_gebiet2.untenLinks);
    m_gebiet2.obenRechts /= norm(m_gebiet2.obenRechts);
  }

  /// Enumeration of sensor types
#include "sensor-types.ncd.enum.hh"  
#include "sensor-types.ncd.cc"

  /// Enumeration of sensor types
#include "sensor-areas.ncd.enum.hh"  
#include "sensor-areas.ncd.cc"

  void writeXML(std::ostream &aus, std::string const &indent = "") const {
    aus << indent << "<sensor"
      " property='" << getSensorTypesName(m_property) << "'"
      " propertyi='" << m_property << "'"
      " area-type='" << getSensorAreasName(m_areaType) << "'"
      " area-typei='" << m_areaType << "'"
      " log-='" << m_log << "'"
      ">\n";
    m_gebiet.writeXML(aus, indent + "\t");
    m_gebiet2.writeXML(aus, indent + "\t");
    aus << indent << "</sensor>\n";
  }
};

#endif
