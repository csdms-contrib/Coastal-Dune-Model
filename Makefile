###############################################################################
#   $Id: Makefile,v 1.26 2005/04/12 14:25:53 duran Exp $
###############################################################################

###### Include and library path for FFTW

FFTW_INC = /usr/local/include
FFTW_LIB = /usr/local/lib

###### Compiler, tools and options

TESTOPT = -g -ggdb

C++	= g++
CFLAGS	= -ansi -Wall -O3 -pedantic-errors -funroll-loops ${TESTOPT}
INCPATH = -I. -I${FFTW_INC}
MY_DEFS =-DPT_VERBOSITY=20 -DPT_SAFETY=1 -DPDE_3D 

LINK    = g++
LFLAGS	= -lm -lfftw3 -L${FFTW_LIB} ${TESTOPT}
LIBS    = $(LIBS_QT) $(LIBS_X11) $(LIBS_OS) 

DOXFILE = dune.dox

###### objects

#OBJECTS = evolution.o dune_evolution.o flux_stationary.o globals.o main.o \
	initsurf.o rotatematrix.o \
	initsurfwaves.o initsurfgauss.o initsurfcone.o initsurfmatlab.o initsurfbeach.o initsurfparabol.o \
	initsurfalea.o analyze_new.o shear.o rfftw12d.o sepbubble.o avalanche.o \
	shear_hlr.o save.o func.o PTG_Func2dScalar.o PTG_Func2dVec.o \
	 wind.o influx.o shore.o vegetation.o storm.o

OBJECTS = evolution.o dune_evolution.o flux_stationary.o globals.o main.o \
	initsurf.o rotatematrix.o \
	initsurfbeach.o initsurfalea.o\
	analyze_new.o shear.o rfftw12d.o sepbubble.o avalanche.o \
	shear_hlr.o save.o func.o PTG_Func2dScalar.o PTG_Func2dVec.o \
	 wind.o influx.o shore.o vegetation.o storm.o

HEADERS_MAIN = $(OBJECTS:.o=.h)
HEADERS = $(HEADERS_MAIN:main.h=)

###### Rules

Dune:	$(OBJECTS) $(PT_OBJECTS)
	$(LINK) -o $@ $^ $(LFLAGS)
	
	cp Dune /Users/orencio/bin/DuneSLC

%.o:	%.cc
	$(C++) -c $(CFLAGS) $(INCPATH) -o $@ $<

%.d:	%.cc
	@echo -e $@ '\c' > $@
	$(C++) -MM $(CFLAGS) $(INCPATH) $< >> $@ || rm -f $@

clean:
	rm -f Dune $(OBJECTS) $(OBJECTS:.o=.d)

doc:	$(OBJECTS:.o=.cc) $(HEADERS) $(DOXFILE)
	doxygen $(DOXFILE)

###### Include object dependencies

-include $(OBJECTS:.o=.d)

