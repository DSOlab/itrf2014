#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_LANG(C++)

AC_INIT([itrf_tools], [1.00], [xanthos@mail.ntua.gr])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/itrf_tools.hpp])
AC_CONFIG_HEADERS([config.h])

# Libtool
LT_PREREQ([2.2])
LT_INIT([dlopen])

# Checks for programs.
AC_PROG_CXX

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lggdatetime':
AC_CHECK_LIB([ggdatetime], [ngpt::dat])

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_TYPE_SIZE_T

# Checks for library functions.

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 bin/Makefile])
AC_OUTPUT
