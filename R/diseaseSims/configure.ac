AC_PREREQ(2.59)
AC_INIT([diseaseSims], 0.1.0)
dnl AC_CONFIG_MACRO_DIR([m4])

CC=`R CMD config CC`
CPP=`R CMD config CPP`
CFLAGS=`R CMD config CFLAGS`
CPPFLAGS="`R CMD config CPPFLAGS` $CPPFLAGS"
AC_PROG_CC
AC_PROG_CPP
AC_PROG_CXX
AC_LANG(C++)
CXXFLAGS="-std=c++11 $CXXFLAGS"
AC_CHECK_HEADER(diseaseSims/ccintermediate.hpp,,[AC_MSG_ERROR(diseaseSims headers not found -- please install)])
AC_CHECK_HEADER(diseaseSims/mutation_with_age.hpp,,[AC_MSG_ERROR(diseaseSims headers not found -- please install)])
AC_CHECK_HEADER(diseaseSims/traitValues.hpp,,[AC_MSG_ERROR(diseaseSims headers not found -- please install)])
AC_CHECK_HEADER(diseaseSims/util.hpp,,[AC_MSG_ERROR(diseaseSims headers not found -- please install)])
AC_CHECK_HEADER(fwdpp/diploid.hh,,[AC_MSG_ERROR(fwdpp headers not found -- please install)])
AC_CHECK_HEADER(boost/interprocess/sync/file_lock.hpp,, [AC_MSG_ERROR([boost/interprocess/sync/file_lock.hpp not found])])
AC_CHECK_HEADER(boost/interprocess/sync/scoped_lock.hpp,, [AC_MSG_ERROR([boost/interprocess/sync/scoped_lock.hpp not found])])
AC_CHECK_HEADER(gsl/gsl_rng.h,,[AC_MSG_ERROR([gsl/gsl_rng.h not found.  Please install the GNU Scientific Library])])
AC_CHECK_HEADER(gsl/gsl_randist.h,,[AC_MSG_ERROR([gsl/gsl_rng.h not found.  Please install the GNU Scientific Library])])
AC_CHECK_HEADER(zlib.h,,[AC_MSG_ERROR([zlib.h not found.  zlib >= 1.2.5 is requred])])

dnl check that libraries are accessible
AC_CHECK_LIB([z],gzungetc,,[AC_MSG_ERROR([zlib run time library not found])])
AC_CHECK_LIB([gslcblas],[cblas_dgemm],,[AC_MSG_ERROR([gslcblas run-time library not found])])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm],,[AC_MSG_ERROR([gsl run-time library not found])])

AC_LANG_SAVE
          AC_LANG_CPLUSPLUS
          AC_REQUIRE_CPP
          AC_CHECK_LIB(sequence, main,,[AC_MSG_ERROR([libsequence runtime library does not seem to be here])])
          AC_CHECK_LIB(boost_system, main,,[AC_MSG_ERROR([boost_system runtime library does not seem to be here])])
	  AC_CHECK_LIB(diseaseSimsHelper,main,,[AC_MSG_ERROR([diseaseSimsHelper runtime library does not seem to be here])])
AC_LANG_RESTORE

AC_CONFIG_FILES([src/Makevars])
XTRA_CPPFLAGS=$CPPFLAGS
XTRA_LDFLAGS=$LDFLAGS
AC_SUBST(XTRA_CPPFLAGS)
AC_SUBST(XTRA_LDFLAGS)
AC_OUTPUT
