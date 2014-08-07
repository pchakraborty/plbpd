dnl Process this file with autoconf to produce a configure script.

AC_PREREQ([2.50])
AC_INIT([plbpd], [1.0], [pchakraborty at gmail dot com])
AC_CONFIG_SRCDIR([src/timer.cpp])
AM_INIT_AUTOMAKE
AM_CONFIG_HEADER([config.h])
AC_CONFIG_FILES([Makefile src/Makefile])


dnl Checks for programs:  mpic++/mpiCC/mpicxx etc
AC_PROG_CXX
AC_CHECK_PROGS(MPICXX, mpic++ mpicxx mpiCC hcp mpxlC_r mpxlC mpCC cmpic++)
if test -z $MPICXX; then
   AC_MSG_ERROR(MPI wrapper for C++ not found)
fi
CXX="$MPICXX"

dnl check for libgsl (depends on libm and libgslcblas)
AC_CHECK_LIB([m],[cos],[],[AC_MSG_ERROR(libm not found)])
AC_CHECK_LIB([gslcblas],[cblas_dgemm],[],[AC_MSG_ERROR(libgslcblas not found)])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm],[],[AC_MSG_ERROR(libgsl not found)])

dnl check for libhdf5
AC_CHECK_LIB([hdf5],[H5Pset_dxpl_mpio],[],[AC_MSG_ERROR(libhdf5 not found)])

dnl Checks for header files.

dnl Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

dnl Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_CHECK_FUNCS([sqrt])

AC_OUTPUT
