include $(top_srcdir)/config/spats.mk

AM_CPPFLAGS = -I$(top_srcdir)/src 

EXTRA_DIST = $(top_srcdir)/AUTHORS $(top_srcdir)/make_bin.sh

AM_CXXFLAGS =

bin_PROGRAMS = \
	relabel_reads \
    match_read_pairs \
    compute_profiles

dist_bin_SCRIPTS = \
	spats.py 

noinst_HEADERS = \
	assert_helpers.h \
	qual.h \
	tokenize.h \
	common.h \
	reads.h \
	common.h \
    bundles.h \
    hits.h

noinst_LIBRARIES = libspats.a

libspats_a_SOURCES = \
	common.cpp \
	qual.cpp \
	reads.cpp \
	tokenize.cpp \
    bundles.cpp \
    hits.cpp

relabel_reads_SOURCES = relabel_reads.cpp  
relabel_reads_LDADD = libspats.a $(BOOST_THREAD_LIB)
relabel_reads_LDFLAGS =  $(BOOST_LDFLAGS)

match_read_pairs_SOURCES = match_read_pairs.cpp  
match_read_pairs_LDADD = libspats.a $(BOOST_THREAD_LIB)
match_read_pairs_LDFLAGS =  $(BOOST_LDFLAGS)

compute_profiles_SOURCES = compute_profiles.cpp  
compute_profiles_LDADD = libspats.a $(BOOST_THREAD_LIB)
compute_profiles_LDFLAGS =  $(BOOST_LDFLAGS)