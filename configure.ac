#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT([Berkeley Range-Energy Calculator],
        [1.6.2],
        [benjamin.weaver@nyu.edu],
        [crange],
        [http://cosmo.nyu.edu/~bw55/dedx/])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_CONFIG_SRCDIR([crange.c])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CC

# Checks for libraries.
AC_CHECK_LIB([m],[cos])
AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm])
AC_CHECK_LIB([iniparser], [iniparser_getboolean])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h unistd.h errno.h time.h iniparser.h gsl/gsl_complex.h gsl/gsl_complex_math.h])

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_CHECK_FUNCS([sqrt])

# Add doxygen configuration
DX_INIT_DOXYGEN(crange, doxygen.cfg, doc)
DX_DOXYGEN_FEATURE(ON)
DX_DOT_FEATURE(OFF)
DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)
DX_MAN_FEATURE(ON)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(ON)
DX_PS_FEATURE(OFF)

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
