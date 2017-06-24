#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

CFLAGS+=" -ansi -Wall -pedantic"
AC_PREREQ([2.69])
AC_INIT([nwalign], [0.1], [nathaniel.echols AT gmail.com])
AC_CONFIG_SRCDIR([src/nwalign.c])

AM_INIT_AUTOMAKE([-Wall -Werror foreign])

# Checks for programs.
AM_PROG_AR
AC_PROG_CC
AC_PROG_LIBTOOL

# Checks for libraries.
PKG_CHECK_MODULES([CHECK], [check >= 0.9.6])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_FUNC_MALLOC

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile src/Makefile tests/Makefile])
AC_OUTPUT
