#--*- makefile -*--------------------------------------------------------------
#
#   Standard makefile
#
#------------------------------------------------------------------------------

# Path to project directory
UPDIR = ../..
# Path to subdirectories
SUBDIR =
# Path to directory for binaries
BINPATH = ../../bin

#
# Libraries and paths (which are not defined globally)
#

LIBS = -lAlign2 -lBiopool -ltools

LIB_PATH = -L.

INC_PATH = -I. -I ../../Biopool/ -I../../tools/ -I../../Align2/Sources

#
# Objects and headers
#

SOURCES =  

OBJECTS =  

TARGETS =   subali \
 

EXECS =  subali \
 

LIBRARY = APPSlibAlign2.a

#
# Install rule
#

compile: all

install: $(LIBRARY) $(TARGETS)
	mv $(LIBRARY) $(UPDIR)/lib
	mv $(EXECS) $(BINPATH)

#
# Call global Makefile to do the job
#

include ../../Makefile.global
