AC_PREREQ(2.59)
AC_INIT([diseaseSims], 0.1.0)
AC_CONFIG_MACRO_DIR([m4])
AC_PROG_CXX
AC_LANG(C++)

dnl AC_CHECK_HEADER(diseaseSims/ccintermediate.hpp,,[AC_MSG_ERROR(diseaseSims headers not found -- please install)])
AC_CHECK_HEADER(boost/interprocess/sync/file_lock.hpp,, [AC_MSG_ERROR([boost/interprocess/sync/file_lock.hpp not found])])
AC_CHECK_HEADER(gsl/gsl_rng.h,,[AC_MSG_ERROR([gsl/gsl_rng.h not found.  Please install the GNU Scientific Library])])
AC_CHECK_HEADER(gsl/gsl_randist.h,,[AC_MSG_ERROR([gsl/gsl_rng.h not found.  Please install the GNU Scientific Library])])
AC_CHECK_HEADER(zlib.h,,[AC_MSG_ERROR([zlib.h not found.  zlib >= 1.2.5 is requred])])

AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
