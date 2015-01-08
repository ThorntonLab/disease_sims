AC_PREREQ(2.59)

AC_INIT([diseaseSims],[0.1.0],[https://github.com/molpopgen/diseaseSims])
AC_CONFIG_SRCDIR([src/TFL2013_ind.cc])
AM_INIT_AUTOMAKE
AC_CONFIG_HEADERS([config.h])

AC_CONFIG_MACRO_DIR([m4])

AC_PROG_CC
AC_C_CONST
AC_PROG_CXX

AM_MAINTAINER_MODE([disable])

AC_LANG(C++)
AX_CXX_COMPILE_STDCXX_11([noext],[mandatory])
AC_PROG_RANLIB
AC_CONFIG_FILES([Makefile diseaseSims/Makefile src/Makefile])

AC_ARG_ENABLE([standard],AS_HELP_STRING([--enable-standard],[Disable use of boost containers.]))
AM_CONDITIONAL([USE_STANDARD],test x$enable_standard = xyes)

AC_ARG_ENABLE([debug],AS_HELP_STRING([--enable-debug],[Compile examples with run-time assert macro enabled.]))
AM_CONDITIONAL([DEBUG],test x$enable_debug = xyes)

dnl check for things that we might need
AC_CHECK_HEADER(Sequence/SimDataIO.hpp,HAVE_LIBSEQ=1 , [AC_MSG_ERROR([Sequence/SimDataIO.hpp not found. libsequence (www.molpopgen.org) >= 1.7.8 required.])])
AC_MSG_CHECKING(that Sequence::SimData is derived from std::pair< std::vector<double>, std::vector<std::string> >)
AC_TRY_COMPILE(
[
#include <Sequence/SimData.hpp>
using namespace Sequence;
using namespace std;
],
[
SimData d;
d.first = vector<double>();
d.second = vector< string >();
],
[AC_MSG_RESULT(yes)
AC_DEFINE([HAVE_SEQUENCE_COALSIM],[],"Sequence::SimData has direct access to underlying data")],
[
AC_MSG_RESULT(no)
AC_MSG_ERROR([Sequence::SimData does not have direct access to underlying data.  Please install libsequence >= 1.8.4])
])

AC_CHECK_HEADER(fwdpp/diploid.hh,HAVE_FWDPP=1 , [AC_MSG_ERROR([fwdpp/diploid.hh not found. fwdpp >= 0.2.5 required])])
AC_CHECK_HEADER(fwdpp/internal/rec_gamete_updater.hpp,HAVE_FWDPP_OK=1 , [AC_MSG_ERROR([fwdpp/rec_gamete_updater.hpp not found. fwdpp >= 0.2.5 required])])
AC_CHECK_HEADER(boost/container/vector.hpp, BVECTOR=1, [AC_MSG_WARN([echo "boost/container/vector.hpp not found.  Will use std::vector to compile."])])
AC_CHECK_HEADER(boost/container/list.hpp, BLIST=1, [echo "boost/container/list.hpp not found.  Will use std::list to compile."])
AC_CHECK_HEADER(boost/interprocess/sync/file_lock.hpp,, [AC_MSG_ERROR([echo "boost/interprocess/sync/file_lock.hpp not found"])])
AC_CHECK_HEADER(boost/archive/binary_iarchive.hpp,, [AC_MSG_ERROR([echo "boost/archive/binary_iarchive.hpp not found"])])
AC_CHECK_HEADER(boost/archive/binary_oarchive.hpp,, [AC_MSG_ERROR([echo "boost/archive/binary_oarchive.hpp not found"])])
AC_CHECK_HEADER(boost/serialization/vector.hpp,, [AC_MSG_ERROR([echo "boost/serialization/vector.hpp not found"])])
AC_CHECK_HEADER(boost/serialization/string.hpp,, [AC_MSG_ERROR([echo "boost/serialization/string.hpp not found"])])
AC_CHECK_HEADER(gsl/gsl_rng.h,HAVE_GSL=1,[AC_MSG_ERROR([gsl/gsl_rng.h not found.  Please install the GNU Scientific Library])])
AC_CHECK_HEADER(gsl/gsl_randist.h,HAVE_GSL=1,[AC_MSG_ERROR([gsl/gsl_rng.h not found.  Please install the GNU Scientific Library])])
AC_CHECK_HEADER(zlib.h,,[AC_MSG_ERROR([zlib.h not found.  zlib >= 1.2.5 is requred])])
AM_CONDITIONAL([BVEC_PRESENT],test x$BVECTOR$enable_boost = x1)
AM_CONDITIONAL([BLIST_PRESENT],test x$BLIST = x1)
AM_CONDITIONAL([BPOOL_PRESENT],test x$BPOOL = x1)
AM_CONDITIONAL([BSET_PRESENT],test x$BSET = x1)
AM_CONDITIONAL([BHASH_PRESENT],test x$BHASH = x1)

dnl check for C run-time libraries
AC_CHECK_LIB([z],gzungetc,,[AC_MSG_ERROR(["zlib run time library not found";exit 1])])
AC_CHECK_LIB([gslcblas],[cblas_dgemm],,[AC_MSG_ERROR(["gslcblas run-time library not found"])])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm],,[AC_MSG_ERROR(["gsl run-time library not found"])])

dnl check for C++ run-time libraries
AC_LANG_SAVE
          AC_LANG_CPLUSPLUS
          AC_REQUIRE_CPP
          AC_CHECK_LIB(sequence, main,,[AC_MSG_ERROR([eecho "libsequence runtime library does not seem to be here";exit 1])])
          AC_CHECK_LIB(boost_system, main,,[AC_MSG_ERROR([echo "boost_system runtime library does not seem to be here";exit 1])])
          AC_CHECK_LIB(boost_program_options, main,,[AC_MSG_ERROR([eecho "boost_program_options runtime library does not seem to be here";exit 1])])
          AC_CHECK_LIB(boost_serialization, main,,[AC_MSG_ERROR([eecho "boost_serialization runtime library does not seem to be here";exit 1])])
AC_LANG_RESTORE

dnl AM_CONFIG_HEADER(config.h)

AC_OUTPUT
