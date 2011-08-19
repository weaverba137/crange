#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT([Berkeley Range-Energy Calculator],
        m4_esyscmd([./crange_version]),
        [benjamin.weaver@nyu.edu],
        [crange],
        [http://cosmo.nyu.edu/~bw55/dedx/])
AC_REVISION([$Revision$])
# AC_REVISION(AC_SVN_REVISION)
# AC_REVISION(AC_SVN_GET_VERSION(AC_SVN_HEADURL))
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_CONFIG_SRCDIR([crange.c])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CC

# Checks for libraries.
AC_CHECK_LIB([m], [exp])
AC_CHECK_LIB([iniparser], [iniparser_getboolean])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h unistd.h iniparser.h])

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_CHECK_FUNCS([sqrt])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
