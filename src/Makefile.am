bin_PROGRAMS=TFL2013 TFL2013_ind make_case_control chksimoutput phenoburden sfs sfsSample

noinst_LTLIBRARIES=libcommon.la
lib_LTLIBRARIES=libdiseaseSimsHelper.la
libcommon_la_SOURCES=ccintermediate.cc
libdiseaseSimsHelper_la_SOURCES=ccintermediate.cc
#atomic_locker/atomic_locker

#TFL2013_ind_SOURCES=
TFL2013_SOURCES=TFL2013.cc
TFL2013_ind_SOURCES=TFL2013_ind.cc
make_case_control_SOURCES=make_case_control.cc simindex.cc readSimOutput.cc
#esm_chisq_zscore_SOURCES=esm_chisq_zscore.cc esm.cc chisq_per_marker.cc esm_filter_sites.cc readCC.cc
chksimoutput_SOURCES=chksimoutput.cc simindex.cc
phenoburden_SOURCES=phenoburden.cc
#CCburden_SOURCES=CCburden.cc readCC.cc
sfs_SOURCES=sfs.cc readSimOutput.cc
sfsSample_SOURCES=sfsSample.cc
# atomic_locker/atomic_locker

#esm_chisq_zscore_LDFLAGS=-pthread
make_case_control_LDADD=libcommon.la
AM_CPPFLAGS=-W -Wall -I.

if DEBUG
else !DEBUG
AM_CPPFLAGS+=-DNDEBUG
endif

TFL2013_LDADD=-lboost_iostreams

