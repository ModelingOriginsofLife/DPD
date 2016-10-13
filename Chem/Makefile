# ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc #
#                                                                     #
#  dpd-simulation                                                     #
#                                      BioMIP Ruhr-Univ. Bochum       #
#                                                                     #
#  (c) 2004                                             Thomas Maeke  #
# ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc #

# $Id: Makefile,v 1.36 2004/09/29 10:19:44 tmaeke Exp $

# uncomment this if not needed

SIMPLE_CONS = 1
WITH_THREAD = 1
# PROFILING = 1
USE_CPP = 1
#PROBABALISTIC_BONDING = 1

###############################
SAMPLES = *.dat
ADDONS = changelog.txt

###############################

###############################
SOURCES =  SDL_etc.c \
           SDL_prim.c \
           dpd_vars.c \
           dpd_stat.c \
           dpd_file.c \
           dpd_util.c \
           dpd_disp.c \
           dpd_rand.c \
           dpd_elec.c \
           dpd_init.c \
           dpd_calc.c \
           dpd_geom.c \
           dpd_flow.c \
	   dpd_dynamic.c \
           dpd_3d.c

OBJS =  dpd_dynamic.o \
		dpd_3d.o \
		SDL_prim.o \
		SDL_etc.o \
		dpd_calc.o \
		dpd_disp.o \
		dpd_elec.o \
		dpd_file.o \
		dpd_flow.o \
		dpd_geom.o \
		dpd_init.o \
		dpd_rand.o \
		dpd_stat.o \
		dpd_util.o \
		dpd_vars.o


HEADERS =  SDL_etc.h \
           SDL_prim.h \
           dpd_vars.h \
           dpd_stat.h \
           dpd_file.h \
           dpd_util.h \
           dpd_disp.h \
           dpd_elec.h \
           dpd_init.h \
           dpd_calc.h \
           dpd_geom.h \
           dpd_flow.h \
           dpd_rand.h \
	   dpd_dynamic.h

###############################

TARGET = dpd3

DEFINES = 

#CFLAGS = -Wall -O5 -I/usr/include/SDL
CFLAGS = -Wall -O5  `sdl-config --cflags`
#CFLAGS = -Wall -g `sdl-config --cflags`

#LFLAGS = -lSDL
LFLAGS = `sdl-config --libs`

DATE := `date +"%Y%m%d"`

CC = gcc
SHELL = bash

SRC = $(SOURCES) 
HDR = $(HEADERS)

WD := $(shell pwd)
DIR := $(shell basename $(WD))

###############################

ifdef USE_CPP
    CC = g++
endif

ifdef PROFILING
    # for Profiler (gprof, kprof)
    LFLAGS += -pg
endif

ifdef WITH_THREAD
    DEFINES += -DWITH_THREAD
endif

ifdef SIMPLE_CONS
    DEFINES += -DSIMPLE_CONS
endif

ifdef PROBABALISTIC_BONDING
   DEFINES  += -DPROBABALISTIC_BONDING
endif

CFLAGS += $(DEFINES)

###############################

$(TARGET):	$(OBJS) $(HDR) Makefile
	$(CC) $(CFLAGS) $(OBJS) $(LFLAGS) -o $(TARGET)

all:	$(TARGET) dist xdist ci

ci:	
	ci -l -mmake $(SOURCES) $(HEADERS) Makefile $(ADDONS)

dist:	
	tar cfvz dpd-$(DATE).tgz \
                    $(SOURCES) $(HEADERS) \
                    Makefile $(TARGET) $(SAMPLES) $(ADDONS)

xdist:	
	tar cfvz dpd3-$(DATE).tgz \
                    $(SOURCES) $(HEADERS) \
                    Makefile $(TARGET) $(SAMPLES) $(ADDONS)

neat:
	rm -f *~ data*.out vertices.out *.reax .reax *.pout .pout .DS_Store

clean:	
	rm -f *~ *.o $(OBJS) $(TARGET) data*.out vertices.out *.reax *.pout .reax .pout .DS_Store

install:
	name=$(DIR) && mv $(TARGET) /usr/local/bin/$${name/#dpd/dpd3}

###
SDL_etc.o: SDL_etc.h
SDL_prim.o: SDL_prim.h
dpd_3d.o: dpd_vars.h dpd_util.h dpd_rand.h dpd_file.h dpd_stat.h dpd_disp.h dpd_elec.h dpd_init.h dpd_calc.h dpd_geom.h dpd_dynamic.h
dpd_calc.o: dpd_calc.h dpd_vars.h dpd_util.h dpd_rand.h dpd_elec.h dpd_flow.h dpd_init.h dpd_dynamic.h
dpd_disp.o: SDL_etc.h SDL_prim.h dpd_vars.h dpd_util.h dpd_stat.h dpd_file.h dpd_disp.h dpd_elec.h dpd_geom.h
dpd_dynamic.o: dpd_dynamic.h dpd_vars.h dpd_util.h
dpd_elec.o: dpd_vars.h dpd_util.h dpd_geom.h dpd_elec.h dpd_disp.h
dpd_file.o: dpd_vars.h dpd_util.h dpd_init.h dpd_file.h dpd_geom.h dpd_elec.h dpd_flow.h dpd_disp.h dpd_dynamic.h
dpd_flow.o: dpd_vars.h dpd_geom.h dpd_flow.h
dpd_geom.o: dpd_vars.h dpd_geom.h dpd_disp.h dpd_init.h
dpd_init.o: dpd_vars.h dpd_util.h dpd_rand.h dpd_disp.h dpd_stat.h dpd_init.h dpd_dynamic.h
dpd_rand.o: dpd_rand.h
dpd_stat.o: dpd_vars.h dpd_util.h dpd_stat.h
dpd_util.o: dpd_vars.h dpd_util.h dpd_dynamic.h
dpd_vars.o: dpd_vars.h SDL_etc.h
