#--*- makefile -*--------------------------------------------------------------
#
#   Standard makefile
#
#------------------------------------------------------------------------------

# Path to project directory.
UPDIR = ..
# Path to subdirectories.
SUBDIR=
# Path to directory for binaries:
BINPATH = ../bin


#
# Libraries and paths (which are not defined globally).
#

LIBS = -lNazgul -lBiopool2000 -ltools -lEnergy -lTorsion

LIB_PATH = -L.

INC_PATH = -I. -I$(UPDIR)/Peso2 -I$(UPDIR)/Energy -I../Energy/TorsionPotential/ -I$(UPDIR)/Biopool2000 -I$(UPDIR)/tools
#
# Objects and headers
#

SOURCES = RamachandranData.cc VectorTransformation.cc LoopTableEntry.cc LoopMTest.cc\
   LoopTable.cc LoopModel.cc LoopExtractor.cc ranking_helper.cc \
   ranking_helper2.cc globalStatistic.cc RankAnalyzer.cc  

OBJECTS = RamachandranData.o VectorTransformation.o LoopTableEntry.o LoopMTest.o\
   LoopTable.o LoopModel.o   LoopExtractor.o ranking_helper.o ranking_helper2.o \
   globalStatistic.o RankAnalyzer.o

TARGETS = loboLUT LoopTablePlot loboFull ClusterLoopTable loboIter LoopMTest \
   loop2torsion scatEdit lobo loboAuto createLoopTestset backboneAnalyzer  

EXECS = loboLUT loboFull loboIter scatEdit lobo loboAuto createLoopTestset LoopMTest \
        loboLUT_all loboAuto_all  

LIBRARY = libNazgul.a

#
# Install rule
#

compile: all

install: $(LIBRARY) $(TARGETS)
	mv $(LIBRARY) $(UPDIR)/lib
	cp $(EXECS) $(BINPATH)

#
# Call global Makefile to do the job.
#

include ../tools/Makefile.global
export VICTOR_ROOT=/home/layla/victor/  PATH=$PATH:/home/layla/victor/bin/