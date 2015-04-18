           __________________________________________________

            SPHERE, A SMOOTHED PARTICLE HYDRODYNAMICS (SPH)
                          SIMULATION SOFTWARE

                           Johannes Willkomm
           __________________________________________________


Table of Contents
_________________

1 Set up Sphere
.. 1.1 Preparations
.. 1.2 Set up H5Part
..... 1.2.1 Obtain H5Part library
..... 1.2.2 Prepare system
..... 1.2.3 Compile and install
.. 1.3 Set up GenNC
.. 1.4 Compile SPH
..... 1.4.1 Build test input files (scenes)
..... 1.4.2 Running Sphere





1 Set up Sphere
===============

1.1 Preparations
~~~~~~~~~~~~~~~~

  Create installation dir and setup PREFIX variable used in this
  document's examples:

  ,----
  | export PREFIX=$HOME/sw/sph-new
  `----

  ,----
  | mkdir $PREFIX
  `----

  Install required system packages:
  ,----
  | sudo apt-get install make g++ libhdf5-dev libvtk5-dev freeglut3-dev octave
  `----

  Optionally, for viewing VTK files we recommend to use Paraview and for
  HDF5 files the `hdf5-tools' package contains some utility programs
  ,----
  | sudo apt-get install paraview hdf5-tools
  `----


1.2 Set up H5Part
~~~~~~~~~~~~~~~~~

1.2.1 Obtain H5Part library
---------------------------

  Download H5Part library
  ,----
  | wget https://codeforge.lbl.gov/frs/download.php/170/H5Part-1.6.1.tar.gz
  `----

  Unpack H5Part library
  ,----
  | tar -xzf H5Part-1.6.1.tar.gz
  `----


1.2.2 Prepare system
--------------------

  On a current Debian Testing system the layout of paths of libhdf5-dev
  is different from what H5Part expects. Therefore we create a directory
  with two symlinks that point to the correct locations:

  ,----
  | pushd /usr/local
  | mkdir hdf5
  | cd hdf5
  | ln -sfT /usr/include/hdf5/serial include
  | ln -sfT /usr/lib/x86_64-linux-gnu/hdf5/serial lib
  | popd
  `----

  Check result:
  ,----
  | ls -l /usr/local/hdf5
  `----


1.2.3 Compile and install
-------------------------

  Now, change directory to `h5part-1.6.1' and run `configure' and `make
  install'
  ,----
  | cd H5Part-1.6.1
  `----

  ,----
  | ./configure --prefix=$PREFIX/h5part --with-hdf5path=/usr/local/hdf5
  `----

  ,----
  | make -j8 install
  `----

  ,----
  | cd ..
  `----


1.3 Set up GenNC
~~~~~~~~~~~~~~~~

  Download the Name Constant Generator tool from
  [https://savannah.nongnu.org/projects/named-constant/]:

  ,----
  | wget http://download.savannah.gnu.org/releases/named-constant/named-constant-0.3.1.tar.gz
  `----

  Unpack:
  ,----
  | tar -xzf named-constant-0.3.1.tar.gz
  `----

  ,----
  | cd named-constant-0.3.1
  | make prefix=$PREFIX install
  `----

  ,----
  | cd ..
  `----

  Now set the GENNC_HOME variable in your shell and also add the /bin
  directory to the $PATH:
  ,----
  | export GENNC_HOME=$PREFIX/gennc
  | export PATH=$PATH:$GENNC_HOME/bin
  `----


1.4 Compile SPH
~~~~~~~~~~~~~~~

  First, you can set up some environment variables to control the
  compiler flags during the build. Recommended settings are:
  ,----
  | export CXXFLAGS="-Wall -Wextra -m64 -O3 -march=native"
  | export CPPFLAGS="-DNDEBUG"
  | export LDFLAGS="-m64"
  `----

  Run the script `bootstrap' to invoke the autotools, which will
  generate the `configure' script, and some other files:
  ,----
  | ./bootstrap
  `----

  Now, run configure, setting the location of the H5Part library that we
  installed and also using the make shift HDF5 directory that we
  created:

  ,----
  | ./configure --prefix=$PREFIX/sph-v1.0.0 --with-h5part=$PREFIX/h5part --with-hdf5=/usr/local/hdf5
  `----

  An alternative is to specify the HDF5 include and library directories
  separately:
  ,----
  | ./configure --prefix=$PREFIX/sph-v1.0.0 --with-h5part=$PREFIX/h5part \
  |  --with-hdf5-incdir=/usr/include/hdf5/serial --with-hdf5-libdir=/usr/lib/x86_64-linux-gnu/hdf5/serial
  `----

  ,----
  | make -j8 install
  `----

  ,----
  | export PATH=$PATH:$PREFIX/sph-v1.0.0/bin
  `----


1.4.1 Build test input files (scenes)
-------------------------------------

  This takes about an hour
  ,----
  | cd sph/vtk
  | make
  `----


1.4.2 Running Sphere
--------------------

  ,----
  | sphere -d1e-4 -r10 vtk/test-cases2/2D/monaghan-ellipsis/mittel/test.vtk
  `----
