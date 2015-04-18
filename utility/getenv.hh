#ifndef sph_getenv_hh
#define sph_getenv_hh
/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

#include <sstream>
#include <string>
#include <limits>
#include <algorithm>
#include <map>

template<int NDIM, class Vec=std::valarray<double> >
struct GetEnvV {
  typedef typename Vec::value_type value_type;
  Vec value;
  value_type const m_min, m_max;
  std::string const &vname;
  bool m_fail;
  std::string comma;

  void parseEnv(char const *envString) {
    std::istringstream instr;
    instr.str(std::string(envString));
    double envValue = 0;
    size_t i = 0;
    for(; i < NDIM and not instr.fail(); ++i) {
      if (i > 0 and instr.peek() == ',') 
        instr >> comma;
      instr >> envValue;
      if (not instr) {
	break;
      }
//       cerr << "read item: " << envValue << "\n";
      value[i] = std::min(std::max(envValue, m_min), m_max);
      if (envValue > m_max) {
	std::cerr << "Error: invalid value of env. variable `" << vname
		  << "': value " << envValue 
		  << " too large, limiting to " << m_max << "\n";
      }
      if (envValue < m_min) {
	std::cerr << "Error: invalid value of env. variable `" << vname
		  << "': value " << envValue
		  << " too small, limiting to " << m_min << "\n";
      }
    }
    m_fail = i != NDIM;
//     std::cerr << "parsed env. var. `" << vname << "' : " << value << "\n";
  }

  GetEnvV(std::string const &_vname, double const def, 
	 value_type const _min = -std::numeric_limits<value_type>::max(), 
	 value_type const _max = std::numeric_limits<value_type>::max()) : 
    value(def, NDIM),
    m_min(_min), 
    m_max(_max),
    vname(_vname),
    m_fail()
  {
    char const *envString = ::getenv(vname.c_str());
    if (envString) {
      parseEnv(envString);
    }
  }

  GetEnvV(std::string const &_vname, Vec const &def, 
	 value_type const _min = -std::numeric_limits<value_type>::max(),
	 value_type const _max = std::numeric_limits<value_type>::max()) : 
    value(def),
    m_min(_min), 
    m_max(_max),
    vname(_vname),
    m_fail()
  {
    char const *envString = ::getenv(vname.c_str());
    if (envString) {
      parseEnv(envString);
    }
  }

  bool fail() const { return m_fail; }
  operator Vec () const { return value; }
};

struct GetEnv : public GetEnvV<1> {
  GetEnv(std::string const &vname, double const def = double(), 
	 double const _min = -std::numeric_limits<double>::max(), 
	 double const _max = std::numeric_limits<double>::max()) : 
    GetEnvV<1>(vname, def, _min, _max)
  {}

  operator double () const { return this->value[0]; }
};

struct GetEnvString {
  std::string value;
  bool present;
  GetEnvString(std::string const &vname, std::string const &def = "") : 
    value(def),
    present()
  {
    char const *envString = ::getenv(vname.c_str());
    if (envString) {
      present = 1;
      value = envString;
//       std::cerr << "parsed env. var. `" << vname << "' : " << value << "\n";
    }
  }

  operator std::string const &() const { return this->value; }

  operator bool() const { return this->good(); }
  bool good() const { return this->present; }
  bool fail() const { return !this->good(); }
};

// environ variable must be declared in the user program
// extern char **environ;

struct AllEnvs {
  char const * const *m_envVarField;
  
  typedef std::map<std::string, std::string> EnvVarMap;
  typedef EnvVarMap::const_iterator CIT;
  typedef EnvVarMap::iterator IT;

  EnvVarMap m_envVars;

  AllEnvs(char **envVarField = environ) : 
    m_envVarField(envVarField) { 
    refresh();
  }

  void refresh() {
    if (not m_envVarField)
      return;
    char const *ptr = m_envVarField[0];
    for(int i = 0; ptr; ++i) {
      std::string f = ptr;
      size_t ePos = f.find_first_of('=');
      if (ePos != std::string::npos) {
        std::string name(ptr, ePos);
        std::string val(ptr + ePos + 1, f.size() - ePos - 1);
        std::cerr << "Param: name " << name << " " << val << "\n";
      }
      ptr = m_envVarField[i+1];
    }
  }

  std::string operator[](std::string const &key) const {
    std::string res;
    CIT it = m_envVars.find(key);
    if (it != m_envVars.end()) {
      res = it->second;
    }
    return res;
  }

  std::string operator[](std::string const &key) {
    std::string res;
    IT it = m_envVars.find(key);
    if (it != m_envVars.end()) {
      res = it->second;
    }
    return res;
  }

  bool exist(std::string const &key) const {
    bool res = 0;
    CIT it = m_envVars.find(key);
    if (it != m_envVars.end()) {
      res = 1;
    }
    return res;
  }

  EnvVarMap const &getMap() { return m_envVars; }

};

#endif
