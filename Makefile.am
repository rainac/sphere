#
# This file is part of Sphere.
# Copyright © 2015 Johannes Willkomm 
# See the file sph/src/sph.hh for copying conditions.  
#

SUBDIRS = utility sph

project_files= README.org LICENSE.md AUTHORS TODO NEWS ChangeLog

EXTRA_DIST = configure.ac

instdir = $(datadir)/doc/$(PACKAGE)

inst_DATA= $(project_files)

install-recursive all-recursive: version.txt

version.txt:
	echo "$(VERSION)" > $@

distclean-local:
	$(RM) -f version.txt
	$(RM) -rf autom4te.cache

METASOURCES = AUTO
