#ifndef jw_sph_vtk_reader_3434_hh
#define jw_sph_vtk_reader_3434_hh
/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>

#include "partikel.hh"

/// The template class VTKReader reads a VTK file using the VTK
/// library and fills a array of type PartikelArray with the values
/// found in the file. The fields to be read in can be configured.

template<class PartikelArray>
struct VTKReader {

  typedef typename PartikelArray::value_type value_type;

  class PropertyNotFound {};
  class FileNotValid {};
  class FileNotUnstructuredGrid {};

  vtkUnstructuredGridReader *reader;
 
  std::string const vtkfname;
  std::string const fieldsToLoad;

  VTKReader(std::string const &_vtkfname, 
	    std::string const &_fieldsToLoad) :
    reader(vtkUnstructuredGridReader::New()),
    vtkfname(_vtkfname),
    fieldsToLoad(_fieldsToLoad)
  {
    reader->SetFileName(vtkfname.c_str());
    // this is the same test as below... 
    //     if (!reader->IsFileValid("unstructured_grid")) {
    //       cerr << "error: file `" << vtkfname << "' is not a valid VTK file\n";
    //       throw FileNotValid();
    //     }
    if (!reader->IsFileUnstructuredGrid()) {
      cerr << "error: file `" << vtkfname << "' is not an UnstructuredGrid\n";
      throw FileNotUnstructuredGrid();
    }
    reader->ReadAllScalarsOn();
    reader->ReadAllVectorsOn();
  }
  
  ~VTKReader() {
    if (reader) {
      reader->Delete();
    }
    reader = 0;
  }

  template<class P, class U>
  void readVTKProperty(PartikelArray &partikel, 
		       U (P::*pptr)(),
		       std::string const &name,
		       vtkPointData &pointData) const {

    vtkDataArray *propArray = pointData.GetArray(name.c_str());
    if (not propArray) {
      cerr << "error: property " << name << " not found in file " << vtkfname << "\n";
      throw PropertyNotFound();
    }
    size_t const bis = partikel.size();
    for (unsigned i = 0; i < bis; ++i) {
      value_type &p = partikel[i];
      (p.*pptr)() = propArray->GetTuple1(i);
    }

//     cerr << "load property `" << name << "'\n";
  }


#ifdef SPH_AD
  template<class P, class U>
  void readVTKPropertyAD(PartikelArray &partikel, 
                         U (P::*pptr)(),
                         std::string const &name,
                         vtkPointData &pointData,
                         int dirIndex) const {

    vtkDataArray *propArray = pointData.GetArray(name.c_str());
    if (not propArray) {
      cerr << "error: property " << name << " not found in file " << vtkfname << "\n";
      throw PropertyNotFound();
    }
    size_t const bis = partikel.size();
    for (unsigned i = 0; i < bis; ++i) {
      value_type &p = partikel[i];
#if SPH_AD == 1
      (p.*pptr)().diff(dirIndex) = propArray->GetTuple1(i);
#else
      if (dirIndex == 0) {
        (p.*pptr)().imag() = propArray->GetTuple1(i) * cvMethEpsilon;
      }
#endif
    }

//     cerr << "load property `" << name << "'\n";
  }
#endif

  template<class P, class U>
  void readVTKPropertyVector(PartikelArray &partikel, 
			     U (P::*pptr)(),
			     std::string const &name,
			     vtkPointData &pointData) const {
    
    vtkDataArray *propArray = pointData.GetArray(name.c_str());
    if (not propArray) {
      cerr << "error: property " << name << " not found in file " << vtkfname << "\n";
      throw PropertyNotFound();
    }

    size_t const bis = partikel.size();
    for (unsigned i = 0; i < bis; ++i) {
      value_type &p = partikel[i];
      for(unsigned j = 0; j < ndim; ++j) {
	(p.*pptr)()[j] = propArray->GetComponent(i, j);
      }
    }

//     cerr << "load vector property `" << name << "'\n";

  }

#ifdef SPH_AD
  template<class P, class U>
  void readVTKPropertyVectorAD(PartikelArray &partikel, 
                               U (P::*pptr)(),
                               std::string const &name,
                               vtkPointData &pointData, 
                               int const dirIndex) const {

    vtkDataArray *propArray = pointData.GetArray(name.c_str());
    if (not propArray) {
      cerr << "error: property " << name << " not found in file " << vtkfname << "\n";
      throw PropertyNotFound();
    }
    cerr << "loading vtk-property " << name << " into derivative direction " << dirIndex << "\n";
    size_t const bis = partikel.size();
    for (unsigned i = 0; i < bis; ++i) {
      value_type &p = partikel[i];
      for(unsigned j = 0; j < ndim; ++j) {
#if SPH_AD == 1
	(p.*pptr)()[j].diff(dirIndex) = propArray->GetComponent(i, j);
#else
        if (dirIndex == 0) {
          (p.*pptr)()[j].imag() = propArray->GetComponent(i, j) * cvMethEpsilon;
        }
#endif
      }
    }
  }
#endif

  void readVTK(PartikelArray &partikel) const {

    bool switchedLocale = false;
    std::locale oldLocale;
    if (reader->GetFileType() == VTK_ASCII) {
      // file type is ASCII
      // set locale to C, otherwise error when parsing ASCII formatted numbers
      oldLocale = std::locale::global(std::locale("C"));
      switchedLocale = true;
    }

    assert(reader);

    // this actually loads the file
    reader->Update();
    // print interesting debug information to cout
//     cout << "<vtk-reader>\n";
//     reader->Print(cout);
//     cout << "</vtk-reader>\n";
    
    // Extract particle infos.
    vtkUnstructuredGrid *grid = reader->GetOutput();
    vtkPointData *pointData = grid->GetPointData();
    
    vtkIdType const numPoints = grid->GetNumberOfPoints();
    partikel.resize(numPoints);

    double x[3];

    for (vtkIdType i=0; i < numPoints; i++) {
      value_type &p = partikel[i];
      grid->GetPoint(i, x);
      for(size_t j = 0; j < ndim; ++j) {
	p.position()[j] = x[j];
      }
    }

    if (fieldsToLoad.find("dichte") != std::string::npos)
      readVTKProperty(partikel, &value_type::dichte, "dichte",  *pointData);
    if (fieldsToLoad.find("masse") != std::string::npos)
      readVTKProperty(partikel, &value_type::masse, "masse",  *pointData);
    if (fieldsToLoad.find("druck") != std::string::npos)
      readVTKProperty(partikel, &value_type::druck, "druck",  *pointData);
    if (fieldsToLoad.find("color") != std::string::npos)
      readVTKProperty(partikel, &value_type::color, "color",  *pointData);
    if (fieldsToLoad.find("hash") != std::string::npos)
      readVTKProperty(partikel, &value_type::hashVal, "hash",  *pointData);
    if (fieldsToLoad.find("class") != std::string::npos)
      readVTKProperty(partikel, &value_type::_class, "class",  *pointData);
    if (fieldsToLoad.find("flags") != std::string::npos)
      readVTKProperty(partikel, &value_type::flags, "flags",  *pointData);
    if (fieldsToLoad.find("sensors") != std::string::npos)
      readVTKProperty(partikel, &value_type::sensorsTriggered, "sensors",  *pointData);
    if (fieldsToLoad.find("waerme") != std::string::npos)
      readVTKProperty(partikel, &value_type::waerme, "waerme",  *pointData);
#ifdef SPH_PARTICLE_VAR_HAS_NEIGHBOUR_SUM
    if (fieldsToLoad.find("neighbours") != std::string::npos)
      readVTKProperty(partikel, &value_type::numNeighbours, "neighbours",  *pointData);
#endif
    
    if (fieldsToLoad.find("vel") != std::string::npos)
      readVTKPropertyVector(partikel, &value_type::velocity, "vel",  *pointData);

#ifdef SPH_AD
    // load derivatives of density
    if (fieldsToLoad.find("ad_dichte") != std::string::npos) {
      for(unsigned dir = 0; dir < SPH_AD_NDIR; ++dir) {
        std::ostringstream propNameStr;
        propNameStr << "ad_dichte_d" << dir;
        std::string propName = propNameStr.str();
        if (pointData->GetArray(propName.c_str())) {
          readVTKPropertyAD(partikel, &value_type::dichte, propName,  *pointData, dir);
        }
      }
    }
    // load derivatives of location
    if (fieldsToLoad.find("ad_coord") != std::string::npos) {
      for(unsigned dir = 0; dir < SPH_AD_NDIR; ++dir) {
        std::ostringstream propNameStr;
        propNameStr << "ad_coord_d" << dir;
        std::string propName = propNameStr.str();
        if (pointData->GetArray(propName.c_str())) {
          readVTKPropertyVectorAD(partikel, &value_type::position, propName,  *pointData, dir);
        }
      }
    }
    // load derivatives of velocity
    if (fieldsToLoad.find("ad_vel") != std::string::npos) {
      for(unsigned dir = 0; dir < SPH_AD_NDIR; ++dir) {
        std::ostringstream propNameStr;
        propNameStr << "ad_vel_d" << dir;
        std::string const propName = propNameStr.str();
        if (pointData->GetArray(propName.c_str())) {
          readVTKPropertyVectorAD(partikel, &value_type::velocity, propName,  *pointData, dir);
        }
      }
    }
#endif // SPH_AD

    if (switchedLocale) {
      std::locale::global(oldLocale);
    }
  }

};

#endif
