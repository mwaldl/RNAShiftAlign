AC_INIT([RNAShiftAlign],[0.1.0],[maria@bioinf.uni-leipzig.de])

AC_CONFIG_SRCDIR([src/RNAShiftAlign.cc])
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([subdir-objects])
AC_CONFIG_HEADERS([config.h])

AC_PROG_CXX
AC_PROG_INSTALL

# check for C++14
AX_CXX_COMPILE_STDCXX_14(,[mandatory])

# configure for use of LocARNA library
PKG_CHECK_MODULES([LOCARNA], [LocARNA-2.0 >= 2.0.0])

# set rpath for LocARNA lib directories
LOCARNA_LIBS=[`echo $LOCARNA_LIBS |sed 's/-L\([^ ]*\)/-L\1 -Wl,-rpath \1 /g'`]

CXXFLAGS="$LOCARNA_CFLAGS $CXXFLAGS"
LIBS="$LIBS $LOCARNA_LIBS"

AC_CONFIG_FILES([
 Makefile
 src/Makefile
])

AC_OUTPUT
