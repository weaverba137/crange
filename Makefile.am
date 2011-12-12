## Process this file with automake to produce Makefile.in
ACLOCAL_AMFLAGS = -I m4
bin_PROGRAMS = crange
crange_SOURCES = crange.h crange.c
dist_doc_DATA = COPYRIGHT INSTALL LGPL README switch.ini target.ini doc/crange.pdf
dist_bin_SCRIPTS = crange_version

include aminclude.am

if DX_COND_man

# You'd probably want to post-process man pages and installed the patched
# versions.
dist_man1_MANS = @DX_DOCDIR@/man/man1/doxample.c.1

$(dist_man1_MANS): doxygen-doc

endif

EXTRA_DIST = autogen.sh $(DX_CONFIG) doc/html
MOSTLYCLEANFILES = $(DX_CLEANFILES)

