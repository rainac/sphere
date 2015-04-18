#ifndef jw_sph_vtk_writer_3438_hh
#define jw_sph_vtk_writer_3438_hh
/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>

#include "utility/strftime.hh"
#include "utility/uname.hh"

#include "partikel.hh"
#include "writer.hh"

#if SPH_AD == 1
using yafad::FO::Static::yafad_value;
#endif

/// The template class VTKWriter saves a SPHPartikelArray to a VTK
/// file using the VTK library itself.
template<class PartikelArray>
struct VTKWriter : public Writer<PartikelArray> {

  typedef typename PartikelArray::value_type particle_type;
  typedef typename particle_type::value_type value_type;

  vtkUnstructuredGridWriter * const writer;
  
  std::string const vtkfname;
  std::string const fieldsToSave;

  VTKWriter(std::string const &_vtkfname, std::string const &fieldsToSave, bool const binary = true) :
    writer(vtkUnstructuredGridWriter::New()),
    vtkfname(_vtkfname),
    fieldsToSave(fieldsToSave)
  {
    if (writer == 0) {
      std::cerr << "error: failed to New() vtkUnstructuredGridWriter: " << strerror(errno) << "\n";
      exit(EXIT_FAILURE);
    }
    writer->SetFileName(vtkfname.c_str());
    std::string header("Written by VTKWriter");
    header += " at ";
    header += StrFTime("%Y-%m-%d %H:%M:%S %z (%c)")();
    header += " on ";
    header += Uname().all();
    writer->SetHeader(header.c_str());
    if (binary) {
      writer->SetFileTypeToBinary();
    } else {
      writer->SetFileTypeToASCII();
    }
  }
  
  ~VTKWriter() {
    if (writer) {
      writer->Delete();
    }
  }

  void writeVTKArray(std::valarray<DoubleVector> const &values,
                     std::string const &name,
                     vtkPointData &pointData, 
                     unsigned const numComp) const {
    vtkDoubleArray *anArray = vtkDoubleArray::New();
    anArray->SetName(name.c_str());
    anArray->SetNumberOfComponents(numComp);
    anArray->SetNumberOfTuples(values.size());

    for (unsigned j = 0; j < numComp; ++j) {
      for (unsigned i = 0; i < values.size(); ++i) {
        anArray->SetComponent(i, j, values[i][j]);
      }
    }

    pointData.AddArray(anArray);
    anArray->Delete();
  }

  template<class P, class U>
  void writeVTKProperty(PartikelArray const &partikel, 
                        std::valarray<DoubleVector> &values,
			std::valarray<unsigned> const &indizes,
			U (P::*pptr)() const,
			std::string const &name,
			vtkPointData &pointData) const {

    values = DoubleVector(); // zero buffer
    
    for (unsigned i = 0; i < indizes.size(); ++i) {
      unsigned const pi = indizes[i];
      particle_type const &p = partikel[pi];
      values[i][0] = 
#ifdef SPH_AD
                        yafad_value((p.*pptr)())
#else
                        (p.*pptr)()
#endif
        ;
    }
    
    writeVTKArray(values, name, pointData, 1);
  }


  template<class P, class U>
  void writeVTKPropertyVector(PartikelArray const &partikel, 
                              std::valarray<DoubleVector> &values,
			      std::valarray<unsigned> const &indizes,
			      U (P::*pptr)() const,
			      std::string const &name,
			      vtkPointData &pointData) const {
    
    values = DoubleVector(); // zero buffer

    for (unsigned i = 0; i < indizes.size(); ++i) {
      unsigned const pi = indizes[i];
      particle_type const &p = partikel[pi];
      for(size_t j = 0; j < ndim; ++j) {
        values[i][j] = yafad_value((p.*pptr)()[j]);
      }
    }
    
    writeVTKArray(values, name, pointData, ndim);
  }


#ifdef SPH_AD

  static std::string adPropName(std::string const &name, unsigned const dir) {
    std::ostringstream nstr;
    nstr << "ad_" << name << "_d" << dir;
    return nstr.str();
  }

  template<class P, class U>
  void writeVTKPropertyAD(PartikelArray const &partikel, 
                          std::valarray<DoubleVector> &values,
                          std::valarray<unsigned> const &indizes,
                          U (P::*pptr)() const,
                          std::string const &name,
                          vtkPointData &pointData, 
                          unsigned const dirIndex) const {

    values = DoubleVector(); // zero buffer
    
    for (unsigned i = 0; i < indizes.size(); ++i) {
      unsigned const pi = indizes[i];
      particle_type const &p = partikel[pi];
#if SPH_AD == 1
      values[i][0] = (p.*pptr)().diff(dirIndex);
#else
      if (dirIndex == 0) {
        values[i][0] = (p.*pptr)().imag() / cvMethEpsilon;
      }
#endif
    }
    
    writeVTKArray(values, adPropName(name, dirIndex), pointData, 1);
  }

  template<class P, class U>
  void writeVTKPropertyAD(PartikelArray const &partikel, 
                          std::valarray<DoubleVector> &values,
                          std::valarray<unsigned> const &indizes,
                          U (P::*pptr)() const,
                          std::string const &name,
                          vtkPointData &pointData) const {
    for (unsigned i = 0; i < SPH_AD_NDIR; ++i) {
      writeVTKPropertyAD<P, U>
        (partikel, values, indizes, pptr, name,  pointData, i);
    }
  }

  template<class P, class U>
  void writeVTKPropertyVectorAD(PartikelArray const &partikel, 
                                std::valarray<DoubleVector> &values,
                                std::valarray<unsigned> const &indizes,
                                U (P::*pptr)() const,
                                std::string const &name,
                                vtkPointData &pointData, 
                                int const dirIndex) const {
    
    values = DoubleVector(); // zero buffer

    for (unsigned i = 0; i < indizes.size(); ++i) {
      unsigned const pi = indizes[i];
      particle_type const &p = partikel[pi];
      for(size_t j = 0; j < ndim; ++j) {
#if SPH_AD == 1
        values[i][j] = (p.*pptr)()[j].diff(dirIndex);
#else
        if (dirIndex == 0) {
          values[i][j] = (p.*pptr)()[j].imag() / cvMethEpsilon;
        }
#endif
      }
    }
    
    writeVTKArray(values, adPropName(name, dirIndex), pointData, ndim);
  }

  template<class P, class U>
  void writeVTKPropertyVectorAD(PartikelArray const &partikel, 
                                std::valarray<DoubleVector> &values,
                                std::valarray<unsigned> const &indizes,
                                U (P::*pptr)() const,
                                std::string const &name,
                                vtkPointData &pointData) const {
    for (unsigned i = 0; i < SPH_AD_NDIR; ++i) {
      writeVTKPropertyVectorAD<P, U>
        (partikel, values, indizes, pptr, name,  pointData, i);
    }
  }
#endif

  void writeVTK(PartikelArray const &partikel, 
		size_t const, 
                int const saveFlags) const {

    assert(writer);
    
    std::valarray<unsigned> pindex2;
    this->findParticlesToSave(partikel, saveFlags, pindex2);
    assert(pindex2.size() <= partikel.size());

    vtkUnstructuredGrid *grid = vtkUnstructuredGrid::New();
    writer->SetInput(grid);
    
    vtkPoints *points = vtkPoints::New();
    points->SetDataTypeToDouble();

    vtkCellArray *cellArray = vtkCellArray::New();

    size_t const npartikel = pindex2.size();

    double xyz[3] = {0, 0, 0};
    for(size_t i = 0; i < npartikel; ++i) {
      particle_type const &p = partikel[pindex2[i]];
      for (size_t j = 0; j < ndim; ++j) {
	xyz[j] = yafad_value(p.position()[j]);
      }
      points->InsertPoint(i, xyz);
      vtkIdType pointId = i;
      cellArray->InsertNextCell(1, &pointId);
    }

    grid->SetPoints(points);
    grid->SetCells(1, cellArray);

    vtkPointData *pointData = grid->GetPointData();

    std::valarray<DoubleVector> values(pindex2.size());

    if (fieldsToSave.find("dichte") != std::string::npos)
      writeVTKProperty<particle_type, value_type const &>
        (partikel, values, pindex2, &particle_type::dichte, "dichte",  *pointData);
    if (fieldsToSave.find("masse") != std::string::npos)
      writeVTKProperty<particle_type, double const &>
        (partikel, values, pindex2, &particle_type::masse, "masse",  *pointData);
    if (fieldsToSave.find("druck") != std::string::npos)
      writeVTKProperty<particle_type, value_type const &>
        (partikel, values, pindex2, &particle_type::druck, "druck",  *pointData);
    if (fieldsToSave.find("color") != std::string::npos)
      writeVTKProperty<particle_type, unsigned short const &>
        (partikel, values, pindex2, &particle_type::color, "color",  *pointData);
    if (fieldsToSave.find("hash") != std::string::npos)
      writeVTKProperty<particle_type, unsigned const &>
        (partikel, values, pindex2, &particle_type::hashVal, "hash",  *pointData);
    if (fieldsToSave.find("class") != std::string::npos)
      writeVTKProperty<particle_type, unsigned char const &>
        (partikel, values, pindex2, &particle_type::_class, "class",  *pointData);
    if (fieldsToSave.find("flags") != std::string::npos)
      writeVTKProperty<particle_type, unsigned char const &>
        (partikel, values, pindex2, &particle_type::flags, "flags",  *pointData);
    if (fieldsToSave.find("sensors") != std::string::npos)
      writeVTKProperty<particle_type, unsigned short const &>
        (partikel, values, pindex2, &particle_type::sensorsTriggered, "sensors",  *pointData);
    if (fieldsToSave.find("waerme") != std::string::npos)
      writeVTKProperty<particle_type, value_type const &>
        (partikel, values, pindex2, &particle_type::waerme, "waerme",  *pointData);

#ifdef SPH_PARTICLE_HAS_ID
    if (fieldsToSave.find("id") != std::string::npos)
      writeVTKProperty<particle_type, unsigned const &>
        (partikel, values, pindex2, &particle_type::id, "id",  *pointData);
#endif
#ifdef SPH_PARTICLE_VAR_HAS_NEIGHBOUR_SUM
    if (fieldsToSave.find("neighbours") != std::string::npos)
      writeVTKProperty<particle_type, unsigned const &>
        (partikel, values, pindex2, &particle_type::numNeighbours, "neighbours",  *pointData);
#endif
#ifdef SPH_PARTICLE_VAR_HAS_ID_SUM
    if (fieldsToSave.find("sumOfNeighbourIds") != std::string::npos)
      writeVTKProperty<particle_type, unsigned const &>
        (partikel, values, pindex2, &particle_type::sumOfNeighbourIds, "sumOfNeighbourIds",  *pointData);
#endif
    
    if (fieldsToSave.find("vel") != std::string::npos)
      writeVTKPropertyVector<particle_type, Vector const &>
        (partikel, values, pindex2, &particle_type::velocity, "vel",  *pointData);

#ifdef SPH_AD
    if (fieldsToSave.find("ad_dichte") != std::string::npos)
      writeVTKPropertyAD<particle_type, value_type const &>
        (partikel, values, pindex2, &particle_type::dichte, "dichte",  *pointData);
    if (fieldsToSave.find("ad_druck") != std::string::npos)
      writeVTKPropertyAD<particle_type, value_type const &>
        (partikel, values, pindex2, &particle_type::druck, "druck",  *pointData);
    if (fieldsToSave.find("ad_waerme") != std::string::npos)
      writeVTKPropertyAD<particle_type, value_type const &>
        (partikel, values, pindex2, &particle_type::waerme, "waerme",  *pointData);

    if (fieldsToSave.find("ad_coord") != std::string::npos)
      writeVTKPropertyVectorAD<particle_type, Vector const &>
        (partikel, values, pindex2, &particle_type::position, "coord",  *pointData);
    if (fieldsToSave.find("ad_vel") != std::string::npos)
      writeVTKPropertyVectorAD<particle_type, Vector const &>
        (partikel, values, pindex2, &particle_type::velocity, "vel",  *pointData);
#endif

    points->Delete();
    cellArray->Delete();

    writer->Update();

    grid->Delete();

  }

};

#endif
