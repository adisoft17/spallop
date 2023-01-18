# spalloop/src/build/GNUmakefile
#
# requires GNU make, tested with GNU make v3.75
#
#
# Targets:
#	spalloop	test compilers, normally invoked as "make DEPENDS=yes spalloop"
#	loc		count lines of source code
#	clean		remove temporaries
#	cleaner		remove almost everything
#	cleanest	even remove emacs ~ files
#	usage		print a message
#	diagnostic	debug the makefile
#
# documentation is in syngate/doc/development/build.txt

# pathnames - relative to spalloop11out/src/build
bindir = ../../bin
libdir = ../../bin
optdir = ../../options
#very good options above; redefine them for the correct path for gdb
srcdir = ..

# determine code name and version, creates variables:
#    CODE
#    VERSION
#    LIB
#    BUILD_DATE
#    ARCHITECTURE
#    HOST_NAME
include $(optdir)/spalloop1o.tag

# determine environment, sets variables:
#    arch
#    os
#    os_ver
#    compiler
#    host
include $(optdir)/spalloop1o.env

# determine definitions of standard UNIX utilities, returns many, including FC
include $(optdir)/unix_utilities

# determine vendor specific compiler flags - this breaks the gui
include $(optdir)/cf.$(os).$(arch).$(compiler)

# determine the location of the various libraries, sets variables:
#    USER_LIBS
# make sure that USER_LIBS is correctly defined
include $(optdir)/spalloop1o.libs

# include C pre-processor flags specific to SYNGATE
# include it if we include source code in C
# include $(optdir)/spalloop.cpp

# create some lists of useful collections of files
files.F90 = $(wildcard $(srcdir)/spalloop_thermo/*.F90)                     \
            $(wildcard $(srcdir)/main/*.F90)                              \
            $(wildcard $(srcdir)/utility/*.F90)                          \
            $(wildcard $(srcdir)/stress/*.F90) 

names.src = $(notdir $(basename $(files.F90) ))
files.d   = $(addsuffix .d, $(names.src))
files.o   = $(libdir)/lib$(LIB).a($(addsuffix .o,$(names.src)))
#files.h   = $(wildcard include/*.h)

# search paths
vpath %.F90 $(srcdir)/spalloop_thermo:                              \
            $(srcdir)/main:                                       \
            $(srcdir)/utility:                                   \
            $(srcdir)/stress:

vpath %.a   $(libdir)

# places to look for included header files - these change with files.F90 and vpath
#INCLUDES = -I$(srcdir)/include
#INCLUDES = 

# include dependency files, only for building code
ifeq ($(DEPENDS),yes)
   -include $(files.d)
endif

# ========== false targets ==========

# declare false targets and precious files
.PHONY: clean cleaner cleanest usage diagnostic FORCE_REBUILD
.PRECIOUS: %.f90

# ========== default rule ==========

# must be first rule to be the default
usage:
	@echo "usage: make DEPENDS=yes spalloop1o"
	@echo "       make loc | clean | cleaner | cleanest | diagnostic"

# ========== pattern rules ==========

# how to generate dependency files from .F90 files
%.d: %.F90
	$(CP) $< $(notdir $(basename $<)).c
	$(CPP) -MM $(INCLUDES) $(notdir $(basename $<)).c | $(SED) -e 's|$*\.c|$*.F90|' -e 's|$*\.o|$*.f90 $@|' > $@
	$(AWK) -f $(optdir)/f90dep.awk -v outfile=$(notdir $(basename $<)) -v lib='$(libdir)/lib$(LIB).a' $< > $@.tmp
	$(SORT) $(SORTFLAGS) $@.tmp >> $@
	$(RM) $(notdir $(basename $<)).c $@.tmp

# how to get .f90 files from .F90 files
%.f90: %.F90
	$(CP) $< $(notdir $(basename $<)).c
	$(CPP) $(CPPFLAGS) $(INCLUDES) $(notdir $(basename $<)).c | $(GREP) -v "^$$" | $(GREP) -v "^#" > $@
	$(RM) $(RMFLAGS) $(notdir $(basename $<)).c

# how to get .o files from .f90 files
%.o: %.f90
	$(FC) -c $(FFLAGS) $(notdir $<)

# ========== rules that force rebuilds ==========

# FORCE_REBUILD is a funny target - anything that depends on FORCE_REBUILD will be rebuilt
FORCE_REBUILD:

# always rebuild main
main.f90: FORCE_REBUILD

# ========== rules ==========

# how to generate the object archive library
$(libdir)/lib$(LIB).a: $(files.o)
ifneq ($(strip $(RANLIB)),)
	$(RANLIB) $@
endif

$(CODE): lib$(LIB).a
	$(FL) $(FFLAGS) -g -O0 $< $(USER_LIBS) -o $(bindir)/$(CODE)

# ========== maintenance targets ==========

# count lines of code (including comments)
loc:
	@PATH=$(optdir):$(PATH); linesofcode $(files.F90)

# remove most of the leftovers
clean:
	$(RM) $(RMFLAGS) *.L
	$(RM) $(RMFLAGS) *.M
	$(RM) $(RMFLAGS) *.T
	$(RM) $(RMFLAGS) *.c
	$(RM) $(RMFLAGS) *.kmo
	$(RM) $(RMFLAGS) *.l
	$(RM) $(RMFLAGS) *.lst
	$(RM) $(RMFLAGS) *.mod
	$(RM) $(RMFLAGS) *.o
	$(RM) $(RMFLAGS) *.out
	$(RM) $(RMFLAGS) *.tmp

# remove more of the leftovers
cleaner: clean
	$(RM) $(RMFLAGS) *.d
	$(RM) $(RMFLAGS) *.f90
	$(RM) $(RMFLAGS) .*.fxref
	$(RM) $(RMFLAGS) $(libdir)/lib$(LIB).a
	$(RM) $(RMFLAGS) $(bindir)/$(CODE)

cleanest: cleaner

# diagnostic target - report the current item of interest
diagnostic:
	@echo "diagnostic target invoked"
	@echo "CFLAGS: $(CFLAGS)"
	@echo "$(CODE)"
	@echo "$(FL) $(FFLAGS) $(CODE)"
	@echo "diagnostic target completed"
