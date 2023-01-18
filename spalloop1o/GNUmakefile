# spalloop/GNUmakefile
#
#                         SPALLOOP Makefile
#
#             * * * * * * * * IMPORTANT * * * * * * * * *
#             *  This makefile is designed specifically *
#             *  for GNU Make Version 3.70 or higher.   *
#             *  GNU make is free software obtainable   *
#             *  from the Free Software Foundation via  *
#             *  anonymous ftp to prep.ai.mit.edu.      *
#             * * * * * * * * IMPORTANT * * * * * * * * *
#
#
# documentation is in syngas/doc/development/build.txt
#
# ================================== targets ===================================
#
# to compile spalloop only
#    % gmake spalloop1o
#
# to clean up your file space
#    % gmake clean              -> removes temporary files, leaves the binary,
#                                  library, debugger and cross reference files
#    % gmake cleaner            -> clean + removes the binary, library,
#                                  debugger, cores and cross reference files
#    % gmake cleanest           -> cleaner + emacs ~ files
#
# to tar the spalloop source code for distribution
#    % gmake tar                -> tars source into ../spalloop.tar.gz
#
# To count the number of lines of code and comments in the entire source:
#    % gmake loc
#
# ==============================================================================

# get the spalloop1o version number (returned in VERSION)
include options/spalloop1o.tag

# build the spalloop1o executable 
spalloop1o:
	@cd src/build; $(MAKE) DEPENDS=yes spalloop1o

# count the lines of (spalloop1o) source code
loc:
	@cd src/build; $(MAKE) loc

# remove most of the leftovers
myclean:
	-@cd src/build;  $(MAKE) clean
clean:
	-@cd src/build;  $(MAKE) clean
#	-@if [ -d problems/tests ] ; then cd problems/tests; $(MAKE) clean; fi

# remove more of the leftovers
cleaner:
	-@cd src/build;  $(MAKE) cleaner
#	-@if [ -d problems/tests ] ; then cd problems/tests; $(MAKE) cleaner; fi

# remove even the emacs ~ files
cleanest distclean realclean superclean:
	@echo "scrubbing the directory structure"
	-@cd src/build;  $(MAKE) cleanest
#	-@if [ -d problems/tests ] ; then cd problems/tests; $(MAKE) cleanest; fi
	find . \( -name \*~ -o -name .\*~ -o -name core \) -exec $(RM) $(RMFLAGS) {} \;

# create a tar archive of the current tree, after cleaning
tar: cleanest
	@echo "writing tar file"
	@rm -f ../spalloop1o-$(VERSION) spalloop1o-$(VERSION).tar*
	@ln -s `pwd` ../spalloop1o-$(VERSION)
	@cd ..; gtar --create --dereference --exclude-from=spalloop1o-$(VERSION)/src/options/tar.exclude --file=- spalloop1o-$(VERSION) | gzip > spalloop1o-$(VERSION)/spalloop1o-$(VERSION).tar.gz
	@rm -f ../spalloop1o-$(VERSION)
	@echo "left distribution in `pwd`/spalloop1o-$(VERSION).tar.gz"

# ==============================================================================

# targets which involve running problems
# can be invoked with "compare=yes" to run diff on the output

# run_all: problem1 problem2 problem3


# ==============================================================================

# uncomment as needed
diagnostic:
	echo "version: $(VERSION)"
	cd src/build; $(MAKE) diagnostic

# ==============================================================================

