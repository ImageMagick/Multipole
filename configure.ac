#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([Multipole], [1.0.0], [http://www.multipole.org], [Multipole])
AC_CONFIG_SRCDIR([multipole.h])
AC_CONFIG_HEADERS([config/config.h])

AC_CANONICAL_TARGET([])

# Ensure that make can run correctly
AM_SANITY_CHECK

AM_INIT_AUTOMAKE([1.12 nostdinc color-tests foreign parallel-tests silent-rules subdir-objects dist-bzip2 dist-lzip dist-xz dist-zip])

# Enable support for silent build rules
m4_ifdef([AM_SILENT_RULES],[AM_SILENT_RULES([yes])])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_RANLIB

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([limits.h malloc.h memory.h stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.

#
# Find math library
#
LIB_MATH=''
AC_CHECK_LIB(m,sqrt,LIB_MATH="-lm",,)
LIBS="$LIB_MATH $LIBS"

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([floor pow sqrt])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
