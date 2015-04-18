#
# This file is part of Sphere.
# Copyright © 2015 Johannes Willkomm 
# See the file sph/src/sph.hh for copying conditions.  
#

SUBDIRS = 

project_files = 

EXTRA_DIST = 

instbasedir = $(pkgdatadir)
instdir = $(instbasedir)

inst_DATA = $(project_files)

NCD_FILES = $(basename $(wildcard *.ncd.xml))

generated_files = $(addsuffix .cc, $(NCD_FILES)) $(addsuffix .hh, $(NCD_FILES)) $(addsuffix .enum.hh, $(NCD_FILES))

AM_CXXFLAGS =
AM_CPPFLAGS =
AM_LDFLAGS =

if amcond_PARALLEL
AM_CXXFLAGS += -fopenmp
AM_LDFLAGS += -fopenmp
endif

AM_CXXFLAGS += -fwhole-program
AM_CPPFLAGS += -I$(VTKHOME)/include/vtk-$(VTKVERSION) -I$(H5PART)/include -I$(HDF5_INCDIR) -Wno-deprecated
AM_LDFLAGS += -L$(H5PART)/lib -L$(HDF5_LIBDIR)

bin_PROGRAMS = sph-1d sph-2d sph-3d sph-scene vtkdump

sph_1d_SOURCES = sph-1d.cc
sph_2d_SOURCES = sph-2d.cc
sph_3d_SOURCES = sph-3d.cc
sph_scene_SOURCES = sph-scene.cc
vtkdump_SOURCES = vtkdump.cc

SPH_LIBS = -lrt -lvtkIO -lvtkFiltering -lvtkCommon -lpthread -lGLU -lGL -lglut -lH5Part -lhdf5

sph_1d_LDADD = $(SPH_LIBS)
sph_2d_LDADD = $(SPH_LIBS)
sph_3d_LDADD = $(SPH_LIBS)
sph_scene_LDADD = $(SPH_LIBS)
vtkdump_LDADD = $(SPH_LIBS)

sph-1d.cc sph-2d.cc sph-3d.cc: $(generated_files)

gen-files: $(generated_files)

%.ncd.cc %.ncd.hh %.ncd.enum.hh: %.ncd.xml
	$(GENNC) -c -i $<

clean-local:
	$(RM) -f $(generated_files)

distclean-local:
	$(RM) -f version.txt
	$(RM) -rf autom4te.cache

METASOURCES = AUTO
