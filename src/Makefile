.SUFFIXES: .h .cpp
#-----------------------------------------------
# Makefile for ModeCode.
# Please modify the variables properly.
# We recommend to use gcc >= 4.8 with openMPI >= 1.8
#-----------------------------------------------


#-----------------------------------------------

CXX = mpicxx
CXXFLAGS = -std=c++11 -Wno-write-strings

INCLUDE = 

CXXL =
LDFLAGS = 

LAPACK =
LIBS = -lfftw3

#-----------------------------------------------
# General rules 
#-----------------------------------------------

.cpp.o:
	${CXX} ${CXXFLAGS} ${INCLUDE} -c $<

PROG = postproc

CXXSRC= pp.cpp main.cpp mem.cpp powerspec.cpp

OBJS= ${CXXSRC:.cpp=.o}

default: postproc

all: ${PROG}

postproc: ${OBJS}
	${CXX} ${CXXFLAGS} -o $@ ${OBJS} ${LIBS} ${LAMMPSLIB}

clean:
	rm -f ${OBJS}

.PHONY: clean

depend:
	gcc -MM *.cpp > .depend

-include .depend
