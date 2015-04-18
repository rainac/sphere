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

generated_files = 

bin_PROGRAMS = cppunittestmain

test_sources = testkernel-1d.cc testpower.cc testSTDSort.cc		\
testkernel-2d.cc testparray.cc testCountingSort.cc testPermutation.cc	\
test1.cc testkernel-3d.cc TestTest.cc MockTestCase.cc

AM_CXXFLAGS =
AM_LDFLAGS =

if amcond_PARALLEL
AM_CXXFLAGS += -fopenmp
AM_LDFLAGS += -fopenmp
endif

cppunittestmain_SOURCES = cppunittestmain.cc $(test_sources)
cppunittestmain_CPPFLAGS = $(shell $(CPPUNIT)/bin/cppunit-config --cflags)
cppunittestmain_LDFLAGS = $(shell $(CPPUNIT)/bin/cppunit-config --libs)

clean-local:
	$(RM) -f $(generated_files)

distclean-local:
	$(RM) -f version.txt
	$(RM) -rf autom4te.cache
