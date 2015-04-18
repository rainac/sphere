/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

#define yafad_value(x) (x)
#include "vtk-reader.hh"
#include "h5part-writer.hh"
#include "partikel.hh"
#include "utility/locale.hh"
#include "utility/hdf5.hh"

#include <vector>
#include <unistd.h>

typedef std::valarray<Partikel<double> > PartikelArray;

struct VTKSeries2H5Part {

  std::string fieldList;

  VTKSeries2H5Part(std::string const &_fieldList) : 
    fieldList(_fieldList) 
  {
    cerr << "list of fields: " << fieldList << "\n";
  }

  void operator()(std::vector<std::string> const &infiles, 
                  std::string &outfile) const {

    H5PartWriter<PartikelArray> writer(outfile, fieldList);

    for(size_t i = 0; i < infiles.size(); ++i) {
      std::string const &iname = infiles[i];
      PartikelArray partikel;
      {
        VTKReader<PartikelArray> reader(iname, fieldList);
        reader.readVTK(partikel);
        cerr << "loaded VTK file `" << iname << "'\n";
      }
      writer.writeH5Part(partikel, i, 0xffff);
      cerr << "wrote H5Part step " << i << " to file `" << outfile << "'\n";
    }
  }
};


void showVersion() {
  std::cout << "VTK series to H5Part converter" << std::endl;
}

void showHelp() {
      cout << "usage: vtk2hdf5 [Options] {vtkfiles} | pattern" << endl;
      cout << "  -h                       show help and exit" << endl;
      cout << "  -v                       show version and exit" << endl;
      cout << "  -o <name>                set output filename" << endl;
      cout << "  -f <fields>              set list of fields to copy" << endl;
      cout << "  -n <number>              set number of files to copy (if pattern given)" << endl;
      cout << "  -s <number>              set index of first file to copy (if pattern given)" << endl;
}

int main(int const argc, char *argv[]) {

  HDF5 hdf5;

  std::string outname = "out.h5", 
    fieldlist = "masse dichte waerme coord vel class color";

  unsigned nmin = 0, nmax = 0;

  std::vector<std::string> infiles;

  char getoptResult = 0;
  char availableOptions[20] = "hvo:f:s:n:";
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
      showHelp();
      return 0;
    case 'v':
      showVersion();
      return 0;
    case 's':
      nmin = atoi(optarg);
      break;
    case 'n':
      nmax = atoi(optarg);
      break;
    case 'o':
      outname = optarg;
      break;
    case 'f':
      fieldlist = optarg;
      break;
    case -1:
      break;
    default:
      cerr << "error: unknown return value from getopt: " << int(getoptResult) << endl;
      break;
    }
  } while (getoptResult != -1);

  for(int i = optind; i < argc; ++i) {
    infiles.push_back(argv[i]);
  }
 
  if (infiles.size() == 1) {
    if (infiles[0].find("%") != std::string::npos) {
      if (nmin == 0 and nmax == 0) {
        cerr << "if you use a %d pattern in the infile name,"
          " you must give the number of steps with -n\n";
        exit(3);
      }
      if (nmin != 0) {
        nmax += nmin;
      }
      
      std::string pattern = infiles[0];
      infiles.clear();

      unsigned const bufferSize = pattern.size() - 2 + 10;
      char *buffer = new char[bufferSize];

      for(unsigned i = nmin; i < nmax; ++i) {
        snprintf(buffer, bufferSize, pattern.c_str(), i);
        infiles.push_back(buffer);
        cerr << "generating " << i << "-th infilename `" << infiles.back() << "\n";
      }
      
      delete[] buffer;
    }
  }

  VTKSeries2H5Part converter(fieldlist);

  converter(infiles, outname);

  return 0;
}

