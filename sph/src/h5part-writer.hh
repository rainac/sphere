#ifndef jw_sph_h5part_writer_3437_hh
#define jw_sph_h5part_writer_3437_hh
/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

# include <H5Part.h>
# include <string.h>
# include <sstream>
# include <errno.h>
# include "partikel.hh"
# include "writer.hh"

#if SPH_AD == 1
using yafad::FO::Static::yafad_value;
#endif

/// This default version tries to generates a compiler warning, prints
/// an error message to cerr and assret()s 0, only the specializations
/// are to be used. 
///
/// \sa writeH5PartVec<double>, writeH5PartVec<h5part_int64_t>
/// Only present if macro SPH_SAVE_H5PART defined.
/// \tparam U the data item type
template<class U>
inline void writeH5PartVec(H5PartFile *file, std::string const &name, U const * tmp) {
  int i;
  std::cerr << "error: not implemented\n";
  assert(0);
}

/// this function writes number-of-particles double items to the
/// H5Part-file file, reading from pointer data and using name as the
/// data item name. number-of-particles is the number set previously by
/// call to H5Part API-function H5PartSetNumParticles().
//
/// \param file open H5Part file pointer
/// \param name data item name
/// \param data field of double items
//
/// \sa writeH5PartVec<h5part_int64_t>
/// Only present if macro SPH_SAVE_H5PART defined.
template<>
inline void writeH5PartVec<double>(H5PartFile *file, std::string const &name, double const *data) {
  H5PartWriteDataFloat64(file, name.c_str(), data);
}

/// this function writes number-of-particles 64-bit integer items to the
/// H5Part-file file, reading from pointer data and using name as the
/// data item name. number-of-particles is the number set previously by
/// call to H5Part API-function H5PartSetNumParticles().
//
/// \param file open H5Part file pointer
/// \param name data item name
/// \param data field of double items
//
/// \sa writeH5PartVec<double>
/// Only present if macro SPH_SAVE_H5PART defined.
template<>
inline void writeH5PartVec<h5part_int64_t>(H5PartFile *file, std::string const &name, h5part_int64_t const *data) {
  H5PartWriteDataInt64(file, name.c_str(), data);
}

/// The template class H5PartWriter saves a SPHPartikelArray to a HDF5
/// file using the H5Part library interface.

template<class PartikelArray>
struct H5PartWriter : public Writer<PartikelArray> {

  typedef typename PartikelArray::value_type particle_type;
  typedef typename particle_type::value_type value_type;

  H5PartFile *file;
  
  std::string const filename;
  std::string const fieldsToSave;

  H5PartWriter(std::string const &_filename, std::string const &_fieldsToSave) :
    file(),
    filename(_filename),
    fieldsToSave(_fieldsToSave)
  {
    file = H5PartOpenFile(filename.c_str(), H5PART_WRITE);
    if (file == 0) {
      std::cerr << "error: H5PartOpenFile failed to open `" << filename
		<< "': " << strerror(errno) << "\n";
      exit(EXIT_FAILURE);
    }
  }
  
  ~H5PartWriter() {
    if (file) {
      closeFile();
    }
  }
  
  void closeFile() {
    long const r = H5PartCloseFile(file);
    if (r) {
      std::cerr << "error: H5PartCloseFile has failed: " << r << "\n";
#ifdef H5_USE_16_API
      H5Eprint(stderr);
#else
      H5Eprint(H5E_DEFAULT, stderr);
#endif
    }
    file = 0;
  }

  static std::string vecPropName(std::string const &name, unsigned const ind) {
    std::ostringstream nstr;
    nstr << name << "_" << ind;
    return nstr.str();
  }

  template<class T1, class T2, class U>
  void writeH5PartPartProperty(T1 (U::*pptr)() const, 
			       std::string const &name, 
			       PartikelArray const &partikel,
			       std::valarray<unsigned> const &pindex,
			       std::valarray<T2> &tmp) const {
    assert(pindex.size() == tmp.size());
    size_t const npleft = tmp.size();
    long const von = 0, bis = npleft;
    for (long i = von; i < bis; ++i) {
      long const pi = pindex[i];
      particle_type const &p = partikel[pi];
      tmp[i] = 
#ifdef SPH_AD
        yafad_value((p.*pptr)())
#else
        (p.*pptr)()
#endif
        ;
    }
    writeH5PartVec(file, name, &tmp[0]);
  }


  template<class T2, class P, class U>
  void writeH5PartPartPropertyVec(U (P::*pptr)() const,
				  std::string const &name,
				  PartikelArray const &partikel,
				  std::valarray<unsigned> const &pindex,
				  std::valarray<T2> &tmp) const {

    for (size_t dim = 0; dim < ndim; ++dim) {
      for (size_t i = 0; i < pindex.size(); ++i) {
	long const pi = pindex[i];
	particle_type const &p = partikel[pi];
	U u = (p.*pptr)();
	tmp[i] = yafad_value(u[dim]);
      }

      writeH5PartVec(file, vecPropName(name, dim), &tmp[0]);
    }
  }

#ifdef SPH_AD
  static std::string adPropName(std::string const &name, unsigned const dir) {
    std::ostringstream nstr;
    nstr << "ad_" << name << "_d" << dir;
    return nstr.str();
  }

  template<class T1, class T2, class U>
  void writeH5PartPartPropertyAD(T1 (U::*pptr)() const, 
                                 std::string const &name, 
                                 PartikelArray const &partikel,
                                 std::valarray<unsigned> const &pindex,
                                 std::valarray<T2> &tmp,
                                 unsigned dirIndex) const {
    for (size_t i = 0; i < pindex.size(); ++i) {
      long const pi = pindex[i];
      particle_type const &p = partikel[pi];
#if SPH_AD == 1
      tmp[i] = (p.*pptr)().diff(dirIndex);
#else
      if (dirIndex == 0) {
        tmp[i] = (p.*pptr)().imag() / cvMethEpsilon;
      }
#endif
    }
    writeH5PartVec(file, adPropName(name, dirIndex), &tmp[0]);
  }

  template<class T1, class T2, class U>
  void writeH5PartPartPropertyAD(T1 (U::*pptr)() const, 
                                 std::string const &name, 
                                 PartikelArray const &partikel,
                                 std::valarray<unsigned> const &pindex,
                                 std::valarray<T2> &tmp) const {
    for (unsigned i = 0; i < SPH_AD_NDIR; ++i) {
      writeH5PartPartPropertyAD<T1, T2, U>(pptr, name, partikel, pindex, tmp, i);
    }
  }

  template<class T1, class T2, class U>
  void writeH5PartPartPropertyVecAD(T1 (U::*pptr)() const, 
                                 std::string const &name, 
                                 PartikelArray const &partikel,
                                 std::valarray<unsigned> const &pindex,
                                 std::valarray<T2> &tmp,
                                 unsigned dirIndex) const {
    for (unsigned j = 0; j < ndim; ++j) {
      for (size_t i = 0; i < pindex.size(); ++i) {
        long const pi = pindex[i];
        particle_type const &p = partikel[pi];
#if SPH_AD == 1
        tmp[i] = (p.*pptr)()[j].diff(dirIndex);
#else
        if (dirIndex == 0) {
          tmp[i] = (p.*pptr)()[j].imag() / cvMethEpsilon;
        }
#endif
      }
      writeH5PartVec(file, vecPropName(adPropName(name, dirIndex), j), &tmp[0]);
    }
  }

  template<class T1, class T2, class U>
  void writeH5PartPartPropertyVecAD(T1 (U::*pptr)() const, 
                                 std::string const &name, 
                                 PartikelArray const &partikel,
                                 std::valarray<unsigned> const &pindex,
                                 std::valarray<T2> &tmp) const {
    for (unsigned i = 0; i < SPH_AD_NDIR; ++i) {
      writeH5PartPartPropertyVecAD<T1, T2, U>(pptr, name, partikel, pindex, tmp, i);
    }
  }
#endif

  void writeH5Part(PartikelArray const &partikel, 
		   size_t const savedFrame, 
                   int const saveFlags) const {
    
    std::valarray<unsigned> pindex2;
    this->findParticlesToSave(partikel, saveFlags, pindex2);

    H5PartSetStep(file, savedFrame);
    H5PartSetNumParticles(file, pindex2.size());

    std::valarray<double> tmp(pindex2.size());
    std::valarray<h5part_int64_t> tmpInt(pindex2.size());

    // the explicit naming of template parameters here is only because of Sun Studio
    // which even in the version Studio Express 2008 Jan. 28 is too stupid to auto-detect them
    if (fieldsToSave.find("masse") != std::string::npos)
      writeH5PartPartProperty<double const &, double,  particle_type>
	(&particle_type::masse,    "masse",   partikel, pindex2,  tmp);
    if (fieldsToSave.find("dichte") != std::string::npos)
      writeH5PartPartProperty<value_type const &, double,  particle_type>
	(&particle_type::dichte,   "dichte",  partikel, pindex2,  tmp);
    if (fieldsToSave.find("waerme") != std::string::npos)
      writeH5PartPartProperty<value_type const &, double,  particle_type>
	(&particle_type::waerme,   "waerme",  partikel, pindex2,  tmp);
    if (fieldsToSave.find("druck") != std::string::npos)
      writeH5PartPartProperty<value_type const &, double,  particle_type>
	(&particle_type::druck,    "druck",   partikel, pindex2,  tmp);
    if (fieldsToSave.find("class") != std::string::npos)
      writeH5PartPartProperty<unsigned char const &,   h5part_int64_t,  particle_type>
	(&particle_type::_class,   "class",   partikel, pindex2,  tmpInt);
    if (fieldsToSave.find("color") != std::string::npos)
      writeH5PartPartProperty<unsigned short const &,  h5part_int64_t,  particle_type>
	(&particle_type::color,    "color",   partikel, pindex2,  tmpInt);
    if (fieldsToSave.find("hash") != std::string::npos)
      writeH5PartPartProperty<unsigned const &,        h5part_int64_t,  particle_type>
	(&particle_type::hashVal,  "hash",    partikel, pindex2,  tmpInt);
    if (fieldsToSave.find("flags") != std::string::npos)
      writeH5PartPartProperty<unsigned char const &,   h5part_int64_t,  particle_type>
	(&particle_type::flags,    "flags",   partikel, pindex2,  tmpInt);
    if (fieldsToSave.find("sensors") != std::string::npos)
      writeH5PartPartProperty<unsigned short const &,   h5part_int64_t,  particle_type>
	(&particle_type::sensorsTriggered,    "sensors",   partikel, pindex2,  tmpInt);

#ifdef SPH_PARTICLE_HAS_ID
    if (fieldsToSave.find("id") != std::string::npos)
      writeH5PartPartProperty<unsigned const &,   h5part_int64_t,  particle_type>
	(&particle_type::id,    "id",   partikel, pindex2,  tmpInt);
#endif
#ifdef SPH_PARTICLE_VAR_HAS_NEIGHBOUR_SUM
    if (fieldsToSave.find("neighbours") != std::string::npos)
      writeH5PartPartProperty<unsigned const &,   h5part_int64_t,  particle_type>
	(&particle_type::numNeighbours,    "neighbours",   partikel, pindex2,  tmpInt);
#endif
#ifdef SPH_PARTICLE_VAR_HAS_ID_SUM
    if (fieldsToSave.find("sumOfNeighbourIds") != std::string::npos)
      writeH5PartPartProperty<unsigned const &,   h5part_int64_t,  particle_type>
	(&particle_type::sumOfNeighbourIds,    "sumOfNeighbourIds",   partikel, pindex2,  tmpInt);
#endif

    if (fieldsToSave.find("coord") != std::string::npos)
      writeH5PartPartPropertyVec<double, particle_type, Vector const &>
	(&particle_type::position, "coord",   partikel, pindex2,  tmp);
    if (fieldsToSave.find("vel") != std::string::npos)
      writeH5PartPartPropertyVec<double, particle_type, Vector const &>
	(&particle_type::velocity, "vel",     partikel, pindex2,  tmp);

#ifdef SPH_AD
    if (fieldsToSave.find("ad_dichte") != std::string::npos)
      writeH5PartPartPropertyAD<value_type const &, double,  particle_type>
	(&particle_type::dichte, "dichte", partikel, pindex2,  tmp);
    if (fieldsToSave.find("ad_druck") != std::string::npos)
      writeH5PartPartPropertyAD<value_type const &, double,  particle_type>
	(&particle_type::druck, "druck", partikel, pindex2,  tmp);
    if (fieldsToSave.find("ad_waerme") != std::string::npos)
      writeH5PartPartPropertyAD<value_type const &, double,  particle_type>
	(&particle_type::waerme, "waerme", partikel, pindex2,  tmp);

    if (fieldsToSave.find("ad_coord") != std::string::npos)
      writeH5PartPartPropertyVecAD<Vector const &, double, particle_type>
	(&particle_type::position, "coord", partikel, pindex2,  tmp);
    if (fieldsToSave.find("ad_vel") != std::string::npos)
      writeH5PartPartPropertyVecAD<Vector const &, double, particle_type>
	(&particle_type::velocity, "vel", partikel, pindex2,  tmp);
#endif

  }

};

#endif
