#
# This file is part of Sphere.
# Copyright © 2015 Johannes Willkomm 
# See the file sph/src/sph.hh for copying conditions.  
#

SUBDIRS =

project_files= 

EXTRA_DIST =

instbasedir = $(pkgdatadir)
instdir = $(instbasedir)

inst_DATA= $(project_files)

NCD_FILES = $(basename $(wildcard *.ncd.xml))

generated_files = $(addsuffix .cc, $(NCD_FILES)) $(addsuffix .hh, $(NCD_FILES)) $(addsuffix .enum.hh, $(NCD_FILES))

all-local: $(generated_files)

%.ncd.cc %.ncd.hh %.ncd.enum.hh: %.ncd.xml
	$(GENNC) -c -i $<

clean-local:
	$(RM) -f $(generated_files)

