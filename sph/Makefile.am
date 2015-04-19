#
# This file is part of Sphere.
# Copyright © 2015 Johannes Willkomm 
# See the file sph/src/sph.hh for copying conditions.  
#

SUBDIRS = sh src tests octave

project_files = 

EXTRA_DIST = 

instbasedir = $(pkgdatadir)
instdir = $(instbasedir)

inst_DATA = $(project_files)

generated_files = 

clean-local:
	$(RM) -f $(generated_files)

distclean-local:
	$(RM) -f version.txt
	$(RM) -rf autom4te.cache
