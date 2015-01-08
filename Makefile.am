AUTOMAKE_OPTIONS = foreign
bin_PROGRAMS=
LIBTOOL_DEPS = @LIBTOOL_DEPS@	
SUBDIRS=diseaseSims src
ACLOCAL_AMFLAGS = -I m4
EXTRA_DIST=src/*.hpp example_workflow/*.sh example_workflow/*.md doc/*.tex R/*.sh R/*.R R/*.md R/diseaseSims/*

libtool: $(LIBTOOL_DEPS)
	$(SHELL) ./config.status libtool
