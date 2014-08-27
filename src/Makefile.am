bin_PROGRAMS=TFL2013 TFL2013_ind make_case_control esm_chisq_zscore chksimoutput phenoburden CCburden sfs sfsSample

#atomic_locker/atomic_locker

#TFL2013_ind_SOURCES=
TFL2013_SOURCES=TFL2013.cc
TFL2013_ind_SOURCES=TFL2013_ind.cc
make_case_control_SOURCES=make_case_control.cc ccintermediate.cc locking_routines.cc simindex.cc
esm_chisq_zscore_SOURCES=esm_chisq_zscore.cc esm.cc chisq_per_marker.cc esm_filter_sites.cc readCC.cc locking_routines.cc
chksimoutput_SOURCES=chksimoutput.cc simindex.cc
phenoburden_SOURCES=phenoburden.cc locking_routines.cc
CCburden_SOURCES=CCburden.cc locking_routines.cc readCC.cc
sfs_SOURCES=sfs.cc readSimOutput.cc
sfsSample_SOURCES=sfsSample.cc
# atomic_locker/atomic_locker



AM_CPPFLAGS=-W -Wall -I.

if DEBUG
else !DEBUG
AM_CPPFLAGS+=-DNDEBUG
endif
LDADD=

if BVEC_PRESENT
if BLIST_PRESENT
if BPOOL_PRESENT
if BHASH_PRESENT
if BSET_PRESENT
#On some systems (or at least recent boost versions, the boost_system is needed to use the boost containers
endif
endif
endif
endif
else
#Boost is missing, so we gotta use the standard headers                                                                   
AM_CPPFLAGS+=-DUSE_STANDARD_CONTAINERS
endif

TFL2013_LDADD=-lboost_iostreams
#This will over-ride the above
if USE_STANDARD
AM_CPPFLAGS+=-DUSE_STANDARD_CONTAINERS
endif