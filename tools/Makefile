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

LIBS =  -ltools

LIB_PATH = -L.

INC_PATH = -I.

#
# Objects and headers
#

SOURCES = vector3.cc matrix3.cc vglStd.cc config.cc GetArg.cc \
 String2Number.cc timer.cc IoTools.cc StatTools.cc    
OBJECTS = vector3.o matrix3.o vglStd.o config.o GetArg.o \
 String2Number.o timer.o IoTools.o StatTools.o  
TARGETS =  

LIBRARY = libtools.a

#
# Install rule
#

compile: all

install: $(LIBRARY) $(TARGET) 
	mkdir ../lib
	mv libtools.a ../lib
	mkdir ../bin	
	mkdir ../Biopool/data	
	mkdir ../Lobo/data
	mkdir ../Energy/data
	export PATH=$(VICTOR_ROOT)bin/:$(PATH)
	ln -s $(VICTOR_ROOT)data/AminoAcidHydrogenData.txt $(VICTOR_ROOT)Biopool/data/AminoAcidHydrogenData.txt
	ln -s $(VICTOR_ROOT)data/allaminoacids_new.pdb $(VICTOR_ROOT)Energy/data/allaminoacids_new.pdb
	ln -s $(VICTOR_ROOT)data/frst.model $(VICTOR_ROOT)Energy/data/frst.model
	ln -s $(VICTOR_ROOT)data/frst2.range $(VICTOR_ROOT)Energy/data/frst2.range
	ln -s $(VICTOR_ROOT)data/polar.par $(VICTOR_ROOT)Energy/data/polar.par
	ln -s $(VICTOR_ROOT)data/ram.par $(VICTOR_ROOT)Energy/data/ram.par
	ln -s $(VICTOR_ROOT)data/solv.par $(VICTOR_ROOT)Energy/data/solv.par
	ln -s $(VICTOR_ROOT)data/tor.par $(VICTOR_ROOT)Energy/data/tor.par
	ln -s $(VICTOR_ROOT)data/tor.par $(VICTOR_ROOT)Lobo/data/tor.par
	ln -s $(VICTOR_ROOT)data/allaminoacids_new_.pdb $(VICTOR_ROOT)Energy/data/allaminoacids_new_.pdb

#
# Call global Makefile to do the job.
#
export VICTOR_ROOT=$(dir $(patsubst %/,%, $(shell pwd)))
include ../Makefile.global

