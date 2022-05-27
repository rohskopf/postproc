.SUFFIXES: .h .cpp
#-----------------------------------------------
# Makefile for ModeCode.
# Please modify the variables properly.
# We recommend to use gcc >= 4.8 with openMPI >= 1.8
#-----------------------------------------------

# Enter the path to your LAMMPS directory
#LAMMPSDIR = /home/drew/lammps-5Jun19/src
#LAMMPSDIR = /home/drew/lammps/src
#LAMMPSDIR = /home/rohskopf/lammps/src
#USE DEEPMD LAMMPS DIR IF YOU WANNA DO FC ON NNP. 
#LAMMPSDIR = /home/rohskopf/deepmd_mod/lammps/src

# Enter path to the shared library (.so) file
#LAMMPSLIB = /home/drew/lammps-5Jun19/src/liblammps_g++_openmpi.so
#LAMMPSLIB = /home/drew/lammps/src/liblammps_g++_openmpi.so
#LAMMPSLIB = /home/rohskopf/lammps/src/liblammps_g++_openmpi.so
#LAMMPSLIB = /home/rohskopf/deepmd_mod/lammps/src/liblammps_g++_openmpi.so

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
