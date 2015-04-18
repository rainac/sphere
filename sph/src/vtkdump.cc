/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

#include <map>
#include <vector>
#include <valarray>
#include <unistd.h>

#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include "spline/point.hh"

static size_t const ndim = 3;
typedef Point<double, ndim> Vector;

typedef std::valarray<double> DoubleArray;
typedef std::valarray<Vector> VectorArray;

typedef std::vector<std::string> StringList;

struct VTKReader {

  typedef double value_type;

  class PropertyNotFound {};
  class FileNotValid {};
  class FileNotUnstructuredGrid {};

  vtkUnstructuredGridReader *reader;
 
  std::string const vtkfname;

  VTKReader(std::string const &_vtkfname) :
    reader(vtkUnstructuredGridReader::New()),
    vtkfname(_vtkfname)
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
    reader->Update();
  }
  
  ~VTKReader() {
    if (reader) {
      reader->Delete();
    }
    reader = 0;
  }

  template<class P, class U>
  void readVTKProperty(DoubleArray &v,
		       vtkDataArray &propArray) const {

    unsigned const bis = propArray.GetNumberOfTuples();
    for (unsigned i = 0; i < bis; ++i) {
      v[i] = propArray.GetTuple1(i);
    }
  }


  void readVTKPropertyVector(VectorArray &v, 
			     vtkDataArray &propArray,
                             unsigned const ncomp) const {
    
    unsigned const bis = propArray.GetNumberOfTuples();
    for (unsigned i = 0; i < bis; ++i) {
      for(unsigned j = 0; j < ncomp; ++j) {
	v[i][j] = propArray.GetComponent(i, j);
      }
    }

  }

  void dumpVTK(std::ostream &aus) const {

    vtkUnstructuredGrid *grid = reader->GetOutput();
    vtkPointData *pointData = grid->GetPointData();
    
    unsigned const numPoints = grid->GetNumberOfPoints();

    VectorArray coordinates(numPoints);

    double x[3];
    for (unsigned i=0; i < numPoints; i++) {
      grid->GetPoint(i, x);
      for(size_t j = 0; j < ndim; ++j) {
	coordinates[i][j] = x[j];
      }
    }

    unsigned const numArrays = pointData->GetNumberOfArrays();
//     cerr << "found " << numArrays << " additional arrays\n";

    std::valarray<VectorArray> data(numArrays);
    std::valarray<unsigned> numComps(numArrays);

    StringList propList(numArrays);

    std::map<std::string, unsigned> propIndizes;

    for (unsigned i = 0; i < numArrays; ++i) {
      std::string const name = pointData->GetArrayName(i);
      vtkDataArray *propArray = pointData->GetArray(name.c_str());
      if (not propArray) {
        cerr << "error: property " << name << " not found in file " << vtkfname << "\n";
        throw PropertyNotFound();
      }

      propList[i] = name;
      propIndizes.insert(make_pair(name, i));

      unsigned const numTupels = propArray->GetNumberOfTuples();
      unsigned const numComp = propArray->GetNumberOfComponents();
//       cerr << " " << name << ": " << numTupels << "x" << numComp << "\n";

      assert(numTupels == numPoints);
      assert(numComp <= ndim);
      
      data[i].resize(numTupels);
      numComps[i] = numComp;
      
      readVTKPropertyVector(data[i], *propArray, numComp);
    }

    std::map<std::string, unsigned>::const_iterator it;

    aus.precision(17);
    aus << "# x y z";
    for (it = propIndizes.begin(); it != propIndizes.end(); ++it) {
      std::string const &name = it->first;
      unsigned i = it->second;
      if (numComps[i] > 1) {
        for (unsigned j = 0; j < numComps[i]; ++j) {
          aus << " " << name << "[" << j << "]";
        }
      } else {
        aus << " " << name;
      }
    }
    aus << "\n";
    for (unsigned pi = 0; pi < numPoints; ++pi) {
      aus << coordinates[pi][0] << " ";
      aus << coordinates[pi][1] << " ";
      aus << coordinates[pi][2] << " ";
      for (it = propIndizes.begin(); it != propIndizes.end(); ++it) {
        unsigned const i = it->second;
        for (unsigned j = 0; j < numComps[i]; ++j) {
          aus << data[i][pi][j] << " ";
        }
      }
      aus << "\n";
    }

  }
};

static int vtkDump(std::string const &vtkname, StringList const &) {

  VTKReader reader(vtkname);

  reader.dumpVTK(cout);

  return 0;
}

static void showVersion() {
  std::cout << "vtkdump " << ndim << "D: dump data in VTK files, one line per time step" << std::endl;
}

int main(int const argc, char *argv[]) {

  std::string inname, outname, fieldlist;

  char getoptResult = 0;
  char availableOptions[20] = "hvP:";

  StringList propList;

  do {
    getoptResult = getopt(argc, argv, availableOptions);
    switch(getoptResult) {
    case '?':
      return 2;
    case ':':
      cerr << "error: argument missing of option: " << char(optopt) << endl;
      break;
    case 'h':
      showVersion();
      cout << "usage: vtkdump VTK-file" << endl;
      cout << "  -h                       show help and exit" << endl;
      cout << "  -v                       show version and exit" << endl;
      return 0;
    case 'v':
      showVersion();
      return 0;
    case 'P': {
      propList.push_back(optarg);
    }
      break;
    case -1:
      break;
    default:
      cerr << "error: unknown return value from getopt: " << int(getoptResult) << endl;
      break;
    }
  } while (getoptResult != -1);

  inname = argv[optind];

  vtkDump(inname, propList);

  return 0;
}

