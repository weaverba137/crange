###############################################################################
#
# This Makefile & all Makefiles in this product are GNU make compliant.
# Please help keep them that way.  See
# http://www.gnu.org/software/make/manual/make.html
#
###############################################################################
#
# Install files into this directory for testing
#
INSTALL_DIR = $(HOME)/Sites/crange
#
# Use this shell to interpret shell commands, & pass its value to sub-make
#
export SHELL = /bin/sh
#
# This is like doing 'make -w' on the command line.  This tells make to
# print the directory it is in.
#
MAKEFLAGS = w
#
# This is a list of subdirectories that make should descend into.  Makefiles
# in these subdirectories should also understand 'make all' & 'make clean'.
# This list can be empty, but should still be defined.
#
SUBDIRS =
#
# Find all *.html files.
#
HTML := $(wildcard *.html)
LAST_MODIFIED := $(shell date +"%Y-%m-%d %H:%M:%S %z (%a, %d %b %Y)")
#
# Be careful with sed
#
ifeq ($(OSTYPE), darwin)
export SED = sed -E
else
export SED = sed -r
endif
#
# This should compile all code prior to it being installed
#
all : dedx.js
	@ for f in $(SUBDIRS); do $(MAKE) -C $$f all ; done
#
# This line helps prevent make from getting confused in the case where you
# have a file named 'all'.
#
.PHONY : all
#
# How to compile to JavaScript.
#
%.js : %.coffee
	coffee -b -c $<
#
# Update modification times in HTML files.
#
last_modified.txt : $(HTML)
	@ for f in $?; do \
		echo $(SED) "s|^Last modified: .*</p>|Last modified: $(LAST_MODIFIED)</p>|" $$f; \
		$(SED) "s|^Last modified: .*</p>|Last modified: $(LAST_MODIFIED)</p>|" $$f > $$f.new; \
		echo /bin/mv -f $$f.new $$f; \
		/bin/mv -f $$f.new $$f; \
		done
	$(RM) $@
	echo "$(LAST_MODIFIED)" > $@
#
# Install things in their proper places in $(INSTALL_DIR)
#
install : all
	@ echo "You will be installing in \$$INSTALL_DIR=$(INSTALL_DIR)."
	@ echo "I'll give you 5 seconds to think about it."
	@ echo ""
	@ sleep 5
	- mkdir -p $(INSTALL_DIR)
	rsync --verbose --recursive --copy-links --times --omit-dir-times \
		--exclude-from=excludes.txt \
		./ $(INSTALL_DIR)/
#
# This line helps prevent make from getting confused in the case where you
# have a file named 'install'.
#
.PHONY : install
#
# Test the install command but do not actually copy files.
#
testinstall : all
	@ echo "You will be installing in \$$INSTALL_DIR=$(INSTALL_DIR)."
	@ echo "This test install will also show what should be deleted."
	@ echo ""
	- mkdir -p $(INSTALL_DIR)
	rsync --verbose --recursive --copy-links --times --omit-dir-times \
		--dry-run --delete \
		--exclude-from=excludes.txt \
		./ $(INSTALL_DIR)/
#
# This line helps prevent make from getting confused in the case where you
# have a file named 'testinstall'.
#
.PHONY : testinstall
#
# GNU make pre-defines $(RM).  The - in front of $(RM) causes make to
# ignore any errors produced by $(RM).
#
clean :
	- $(RM) *~ core
	@ for f in $(SUBDIRS); do $(MAKE) -C $$f clean ; done
#
# This line helps prevent make from getting confused in the case where you
# have a file named 'all'.
#
.PHONY : clean
