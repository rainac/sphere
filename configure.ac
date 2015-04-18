# Process this file with autoconf to produce a configure script.
#
# This file is part of Sphere.
# Copyright © 2015 Johannes Willkomm
# See the file sph/src/sph.hh for copying conditions.  
#

AC_INIT(sphere, 0.5.0, [johannes@johannes-willkomm.de], 
                sphere, [https://github.com/rainac/sphere])

AC_CONFIG_AUX_DIR(config)
AC_CANONICAL_TARGET
AM_CONFIG_HEADER([config.h])
AM_INIT_AUTOMAKE()

AC_PROG_MAKE_SET
AC_PROG_CXX
AC_PROG_CPP

AC_CHECK_PROGS([CP], [cp copy])
AC_CHECK_PROGS([RM], [rm del])
AC_CHECK_PROGS([MV], [mv rename])

AC_CHECK_PROGS([GENNC], [gennc])
if test "x$GENNC" = "x"; then
   AC_MSG_ERROR([mandatory program gennc was not found])
fi

AC_HEADER_STDC
AC_CHECK_HEADERS([unistd.h], , AC_MSG_ERROR([mandatory header file not found!]))
AC_CHECK_HEADERS([omp.h sys/stat.h])

AC_DEFINE_UNQUOTED(install_prefix, "$prefix", [The prefix where Sphere is installed.])

AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
AC_TYPE_PID_T
AC_CHECK_TYPES(long double)
# Find a type that is 32 bit wide
AC_CHECK_TYPES(uint32_t)
AC_CHECK_SIZEOF(unsigned int)

AC_CHECK_FUNCS([strchr strerror strrchr strstr strtol strtod strtoul], ,AC_MSG_ERROR([mandatory function not found!]))
AC_CHECK_FUNCS([atexit setenv getenv strdup stricmp strcasecmp memset])

AC_ARG_ENABLE([debug], [  --enable-debug          Debug features.], [
      build_debug=$enableval
], [ 
      build_debug=no
])
AC_MSG_CHECKING([setting debug])
AC_MSG_RESULT([$build_debug])
if test "$build_debug" = "yes"; then
   AC_DEFINE_UNQUOTED([SPH_DEBUG], [1], [Debug features])
fi
AM_CONDITIONAL(amcond_DEBUG, [test "x$build_debug" = "xyes"])

AC_ARG_ENABLE([parallel], [  --enable-parallel       Build with OpenMP.], [
      build_parallel=$enableval
], [ 
      build_parallel=yes
])
AC_MSG_CHECKING([setting parallel])
AC_MSG_RESULT([$build_parallel])
if test "$build_parallel" = "yes"; then
   AC_DEFINE_UNQUOTED([SPH_PARALLEL], [1], [Parallelize with OpenMP])
fi
AM_CONDITIONAL(amcond_PARALLEL, [test "x$build_parallel" = "xyes"])

AC_ARG_ENABLE([render], [  --enable-render         Build with OpenGL rendering.], [
      build_render=$enableval
], [ 
      build_render=yes
])
AC_MSG_CHECKING([setting render])
AC_MSG_RESULT([$build_render])
if test "$build_render" = "yes"; then
   AC_DEFINE_UNQUOTED([SPH_RENDER], [1], [On-line rendering with OpenGL])
fi
AM_CONDITIONAL(amcond_RENDER, [test "x$build_render" = "xyes"])

AC_ARG_ENABLE([save-h5part], [  --enable-save-h5part    Build with save files in HDF5.], [
      build_save_h5part=$enableval
], [ 
      build_save_h5part=yes
])
AC_MSG_CHECKING([setting save_h5part])
AC_MSG_RESULT([$build_save_h5part])
if test "$build_save_h5part" = "yes"; then
   AC_DEFINE_UNQUOTED([SPH_SAVE_H5PART], [1], [Support for saving results with H5Part])
fi
AM_CONDITIONAL(amcond_SAVE_H5PART, [test "x$build_save_h5part" = "xyes"])

AC_MSG_CHECKING([Set VTK home])
AC_ARG_WITH([vtk-home], [  --with-vtk-home         Set the VTK library home directory.], 
                  [VTKHOME=$withval],
                  [VTKHOME=/usr])
AC_MSG_RESULT([$VTKHOME])
AC_SUBST(VTKHOME)

AC_MSG_CHECKING([Set VTK version])
AC_ARG_WITH([vtk-version], [  --with-vtk-version      Set the VTK library version.], 
                  [VTKVERSION=$withval],
                  [VTKVERSION=5.8])
AC_MSG_RESULT([$VTKVERSION])
AC_SUBST(VTKVERSION)

AC_MSG_CHECKING([Freeglut home])
AC_ARG_WITH([freeglut], [  --with-freeglut=path    Specify the path where the Freeglut library resides.],
      [if test ! "x$withval" = "x" ; then
         FREEGLUT3="$withval"
      fi], [
         FREEGLUT3=/usr
      ]
)
AC_SUBST(FREEGLUT3)
AC_MSG_RESULT([$FREEGLUT3])

AC_MSG_CHECKING([H5Part home])
AC_ARG_WITH([h5part], [  --with-h5part=path      Specify the path where the H5part library resides.],
      [if test ! "x$withval" = "x" ; then
         H5PART="$withval"
      fi], [
         H5PART=/usr
      ]
)
AC_SUBST(H5PART)
AC_MSG_RESULT([$H5PART])

AC_MSG_CHECKING([HDF5 home])
AC_ARG_WITH([hdf5], [  --with-hdf5=path        Specify the path where the HDF5 library resides.],
      [if test ! "x$withval" = "x" ; then
         HDF5="$withval"
      fi], [
         HDF5=/usr
      ]
)
AC_SUBST(HDF5)
AC_MSG_RESULT([$HDF5])

AC_MSG_CHECKING([HDF5 includes])
AC_ARG_WITH([hdf5-incdir], [  --with-hdf5-incdir=path    Specify the path where the HDF5 include headers reside.],
      [if test ! "x$withval" = "x" ; then
         HDF5_INCDIR="$withval"
      fi], [
         HDF5_INCDIR=$HDF5/include
      ]
)
AC_SUBST(HDF5_INCDIR)
AC_MSG_RESULT([$HDF5_INCDIR])

AC_MSG_CHECKING([HDF5 libs])
AC_ARG_WITH([hdf5-libdir], [  --with-hdf5-libdir=path    Specify the path where the HDF5 library resides.],
      [if test ! "x$withval" = "x" ; then
         HDF5_LIBDIR="$withval"
      fi], [
         HDF5_LIBDIR=$HDF5/lib
      ]
)
AC_SUBST(HDF5_LIBDIR)
AC_MSG_RESULT([$HDF5_LIBDIR])

AC_MSG_CHECKING([CPPUnit home])
AC_ARG_WITH([cppunit], [  --with-cppunit=path     Specify the path where the CPPUnit library resides.],
      [if test ! "x$withval" = "x" ; then
         CPPUNIT="$withval"
      fi], [
         CPPUNIT=/usr
      ]
)
AC_SUBST(CPPUNIT)
AC_MSG_RESULT([$CPPUNIT])

AC_CONFIG_FILES([Makefile
                 utility/Makefile
                 sph/Makefile
                 sph/sh/Makefile
                 sph/src/Makefile
                 sph/tests/Makefile
])

AC_OUTPUT

