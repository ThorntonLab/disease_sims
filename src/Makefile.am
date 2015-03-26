bin_PROGRAMS=TFL2013 TFL2013_ind make_case_control make_case_controlEW chksimoutput phenoburden sfs riskmatrix sampleSFS genburden

noinst_LTLIBRARIES=libcommon.la
lib_LTLIBRARIES=libdiseaseSimsHelper.la
libcommon_la_SOURCES=ccintermediate.cc util.cc
libdiseaseSimsHelper_la_SOURCES=ccintermediate.cc util.cc
#atomic_locker/atomic_locker

#TFL2013_ind_SOURCES=
TFL2013_SOURCES=TFL2013.cc
TFL2013_ind_SOURCES=TFL2013_ind.cc TFL2013_ind_params.cc
make_case_control_SOURCES=make_case_control.cc simindex.cc readSimOutput.cc
make_case_controlEW_SOURCES=make_case_controlEW.cc simindex.cc readSimOutput.cc
#esm_chisq_zscore_SOURCES=esm_chisq_zscore.cc esm.cc chisq_per_marker.cc esm_filter_sites.cc readCC.cc
chksimoutput_SOURCES=chksimoutput.cc simindex.cc
phenoburden_SOURCES=phenoburden.cc
genburden_SOURCES=genburden.cc
#CCburden_SOURCES=CCburden.cc readCC.cc
sfs_SOURCES=sfs.cc readSimOutput.cc
riskmatrix_SOURCES=riskmatrix.cc
sampleSFS_SOURCES=sampleSFS.cc
# atomic_locker/atomic_locker

#esm_chisq_zscore_LDFLAGS=-pthread
make_case_control_LDADD=libcommon.la
make_case_controlEW_LDADD=libcommon.la
riskmatrix_LDADD=libcommon.la

AM_CPPFLAGS=-W -Wall -I.

if DEBUG
else !DEBUG
AM_CPPFLAGS+=-DNDEBUG
endif

TFL2013_LDADD=-lboost_iostreams

