################################################################################
# Makefile for building a HepMC executable with the gcc compiler
#
# Matt.Dobbs@CERN.CH  1.2000
#
# This makefiles also works to compile the example programs.
# I.E.: syntax for compiling the example_MyPythia.cc example:
#       gmake example_MyPythia.exe
# or simply   gmake all     to compile all examples.
#
################################################################################ Define directory paths 
#    You may have to change GENSERdir and/or other variables
#
  HepMCdir             = /afs/cern.ch/work/a/antoniov/Workspace/HEP-Generators/HepMC-2.06.08-x86_64-slc6-gcc481
#  HepMCdir             = /afs/cern.ch/sw/lcg/external/HepMC/2.06.05/x86_64-slc5-gcc46-opt
#  ROOTdir              = 
  HepMClib             = -L$(HepMCdir)/lib -lHepMC
  HepMCfiolib          = -L$(HepMCdir)/lib -lHepMCfio
  GENSERdir            =  
#  CLHEPdir             = /afs/cern.ch/sw/lcg/external/clhep/2.1.4.1/x86_64-slc5-gcc46-opt
  HepPDTdir            = /afs/cern.ch/work/a/antoniov/Workspace/HEP-Generators/HepPDT-3.04.01-x86_64-slc6-gcc481
  HepPDTlib            = -L$(HepPDTdir)/lib -lHepPDT -lHepPID
#  Pythia_LIB	= -L$(GENSERdir)/lib -Wl,-rpath -Wl,$(GENSERdir)/lib \
                  -lpythia6_403 -lpythia6_403_dumm -lpythia6_403_pdfdumm
#  Herwig_LIB	= -L$(GENSERdir)/lib -Wl,-rpath -Wl,$(GENSERdir)/lib \
                  -lherwig6_510 -lherwig6_510_dumm -lherwig6_510_pdfdumm
#  ROOTLIBS = $(shell root-config --libs) -lMinuit -lEG
  ROOTLIBS    = $(shell root-config --libs)
  FASTJETLIBS = $(shell /afs/cern.ch/work/a/antoniov/Workspace/HEP-Generators/fastjet-3.1.3-x86_64-slc6-gcc481/bin/fastjet-config --libs --plugins)
################################################################################ Compiler options
#
  CXX           = g++
  INCLUDES 	= -I$(HepMCdir)/include -I$(CLHEPdir)/include -I$(HepPDTdir)/include
  CXXFLAGS      =  -ansi -pedantic -Wall -Wno-long-long -g -O2 $(shell root-config --cflags) $(shell /afs/cern.ch/work/a/antoniov/Workspace/HEP-Generators/fastjet-3.1.3-x86_64-slc6-gcc481/bin/fastjet-config --cxxflags) $(INCLUDES)
ifeq "$(CXX)" "g++"
   F77           = g77
   FLAGS 	= $(DFLG) -fno-second-underscore $(INCDIR)

else
  F77           = f77
  FLAGS 	= $(DFLG) $(INCDIR)
endif

  LINK_LIBS     =  
UNAME = $(shell uname)
ifeq "$(UNAME)" "Darwin"
  CLHEP_LIB     = -L$(CLHEPdir)/lib -lCLHEP 
else
  CLHEP_LIB     = -L$(CLHEPdir)/lib -lCLHEP -Wl,-r -Wl,$(CLHEPdir)/lib
  LINK_LIBS     += -Wl,-rpath -Wl,$(HepMCdir)/lib -Wl,-rpath -Wl,$(HepPDTdir)/lib
endif

  CLHEP_LIB     = -L$(CLHEPdir)/lib -lCLHEP -Wl,-rpath -Wl,$(CLHEPdir)/lib
  initpydata_OBJ= initpydata.o
  pythia_OBJ    = initPythia.o initpydata.o
  #SRCS	        = example_BuildEventFromScratch.cc	\
		  example_EventSelection.cc	\
		  example_MyHerwig.cc	\
		  example_MyPythia.cc	\
		  example_MyPythiaOnlyToHepMC.cc	\
		  example_UsingIterators.cc	\
		  example_UsingIterators_teste.cc	\
		  example_PythiaStreamIO.cc	\
		  testPythiaCopies.cc	\
		  testHerwigCopies.cc	\
		  initPythia.cc		\
		  initpydata.f
  HDRS          = $(HepMCdir)/include/HepMC/*.h *.h
  #EXAMPLES	= example_BuildEventFromScratch.exe	\
		  example_EventSelection.exe	\
		  example_MyHerwig.exe	\
		  example_MyPythia.exe	\
		  example_MyPythiaOnlyToHepMC.exe	\
		  example_UsingIterators.exe	\
		  example_UsingIterators_teste.exe	\
		  example_PythiaStreamIO.exe	\
		  testPythiaCopies.exe	\
		  testHerwigCopies.exe

################################################################################ 

PLATFORM=$(shell uname)

ifeq "$(PLATFORM)" "Linux"
#    LINK_LIBS     += -lg2c 
    LINK_LIBS     +=
endif

################################################################################ definition of the compiler options
#	-I location of directory containing include files
#	-L location of directory containing libraries
#       -lname include the library from -L location called libname.a
#	   -lg2c is the library containing info on converting fortran to C
#          -lf   is the library containing the intrinsic for HPUX only.
#	-shared make a shared library as output
#	-fPIC produce position independent code
#        necessary on some platforms (including HPUX) for -shared
# 	-fpic ^^ same(?)
#	-O optimizes
#	-g produces output for the debugger
#       -pg produces output for gprof profiler
#       note: if you want to see all warnings and ensure ansi standard 
#             compatibility, use:
#             -pipe -ansi -pedantic -fnonnull-objects \
#             -W -Wall -Wwrite-strings -Wpointer-arith -Wnested-externs \
#             -Woverloaded-virtual -Wbad-function-cast -fnonnull-objects
#       The proper order for cernlib libraries is:
#       -lpawlib -lgraflib -lgrafX11 -lmathlib -lkernlib -lpacklib -ljetset74
#
# makefile syntax:
#        for target thedir/target.suf from source anotherdir/source.suf2
#        ${*D}  = thedir
#        ${*F}  = target
#        $*     = thedir/target
#        $@     = thedir/target.suf
#        $<     = anotherdir/source.suf2
#  

###############################################################################
#
.SUFFIXES:      .o .cxx .f .exe
#all:	$(EXAMPLES)
all:	analysis-SD-Dijets-HepMC.exe analysis-gamgamZZ-HepMC.exe

example-HepMC-ROOT.exe: example-HepMC-ROOT.o
	@echo "Building $@ ..."
	$(CXX) $(FLAGS) example-HepMC-ROOT.o \
		$(HepMClib) \
		$(HepPDTlib) \
                $(ROOTLIBS) \
	        $(LINK_LIBS) -o $@

analysis-SD-Dijets-HepMC.exe: analysis-SD-Dijets-HepMC.o
	@echo "Building $@ ..."
	$(CXX) $(FLAGS) analysis-SD-Dijets-HepMC.o \
		$(HepMClib) \
		$(HepPDTlib) \
                $(ROOTLIBS) \
	        $(LINK_LIBS) -o $@

analysis-gamgamZZ-HepMC.exe: analysis-gamgamZZ-HepMC.o
	@echo "Building $@ ..."
	$(CXX) $(FLAGS) analysis-gamgamZZ-HepMC.o \
		$(HepMClib) \
		$(HepPDTlib) \
                $(ROOTLIBS) \
                $(FASTJETLIBS) \
	        $(LINK_LIBS) -o $@

example_BuildEventFromScratch.exe: example_BuildEventFromScratch.o
	@echo "Building $@ ..."
	$(CXX) $(FLAGS) example_BuildEventFromScratch.o \
		$(HepMClib) \
		$(CLHEP_LIB) \
	        $(LINK_LIBS) -o $@

example_EventSelection.exe: example_EventSelection.o
	@echo "Building $@ ..."
	$(CXX) $(FLAGS) example_EventSelection.o \
		$(HepMClib) \
	        $(LINK_LIBS) -o $@

example_MyHerwig.exe: example_MyHerwig.o
	@echo "Building $@ ..."
	$(CXX) $(FLAGS) example_MyHerwig.o \
		$(HepMClib) $(HepMCfiolib) \
	        $(Herwig_LIB) $(LINK_LIBS) -o $@

example_MyPythia.exe: $(initpydata_OBJ) $(pythia_OBJ) example_MyPythia.o
	@echo "Building $@ ..."
	$(CXX) $(FLAGS) $(pythia_OBJ) example_MyPythia.o \
		$(HepMClib) $(HepMCfiolib) \
	        $(Pythia_LIB) $(LINK_LIBS) -o $@

example_MyPythiaOnlyToHepMC.exe: $(initpydata_OBJ) $(pythia_OBJ) example_MyPythiaOnlyToHepMC.o
	@echo "Building $@ ..."
	$(CXX) $(FLAGS) $(pythia_OBJ) example_MyPythiaOnlyToHepMC.o \
		$(HepMClib) $(HepMCfiolib) \
	        $(Pythia_LIB) $(LINK_LIBS) -o $@

example_PythiaStreamIO.exe: $(initpydata_OBJ) $(pythia_OBJ) example_PythiaStreamIO.o
	@echo "Building $@ ..."
	$(CXX) $(FLAGS) $(pythia_OBJ) example_PythiaStreamIO.o \
		$(HepMClib) $(HepMCfiolib) \
	        $(Pythia_LIB) $(LINK_LIBS) -o $@

example_UsingIterators.exe: example_UsingIterators.o
	@echo "Building $@ ..."
	$(CXX) $(FLAGS) example_UsingIterators.o \
		$(HepMClib) \
	        $(LINK_LIBS) -o $@

example_UsingIterators_teste.exe: example_UsingIterators_teste.o
	@echo "Building $@ ..."
	$(CXX) $(FLAGS) example_UsingIterators_teste.o \
		$(HepMClib) \
                $(ROOTLIBS) \
	        $(LINK_LIBS) -o $@

testPythiaCopies.exe: $(initpydata_OBJ) $(pythia_OBJ) testPythiaCopies.o
	@echo "Building $@ ..."
	$(CXX) $(FLAGS) $(pythia_OBJ) testPythiaCopies.o \
		$(HepMClib) $(HepMCfiolib) \
	        $(Pythia_LIB) $(LINK_LIBS) -o $@

testHerwigCopies.exe: testHerwigCopies.o
	@echo "Building $@ ..."
	$(CXX) $(FLAGS) testHerwigCopies.o \
		$(HepMClib) $(HepMCfiolib) \
	        $(Herwig_LIB) $(LINK_LIBS) -o $@

###############################################################################
# instructions for building a .o file from a .cxx file
#
.cc.o:         $(HDRS) $<
	@echo "Compiling $< with $(CXX) ..."
	@$(CXX) $(CXXFLAGS) -c $< -o $@

###############################################################################
# instructions for building a .o file from a .f file
#
.f.o:           $<
	@echo "Compiling $< with $(F77) ..."
	@$(F77) $(FLAGS) -c $< -o $@

###############################################################################
# gmake clean       removes all garbage from HepMC directories.
#
clean:
	rm -f *.o

###############################################################################
# gmake distclean       removes all compiled libraries, executables, +garbage
#                       to prepare the package for distribution
distclean: 
	$(MAKE) clean --no-print-directory
	rm -f *.exe
	rm -f *.dat

