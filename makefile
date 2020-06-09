# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

# $Id: Makefile.in 726 2006-04-17 04:16:00Z andreasw $

##########################################################################
#    You can modify this example makefile to fit for your own program.   #
#    Usually, you only need to change the five CHANGEME entries below.   #
##########################################################################

# To compile other examples, either changed the following line, or
# add the argument EXNAME=example_name to make
EXNAME = PSCG
#EXNAME = investment

# CHANGEME: This should be the name of your executable
EXE = $(EXNAME)

# CHANGEME: Here is the name of all object files corresponding to the source
#           code that you wrote in order to define the problem statement
SRCDIR = src
SRC = $(SRCDIR)/PSCG.cpp 
OBJDIR = obj
OBJS = $(OBJDIR)/PSCG.o

# CHANGEME: Additional libraries
ADDLIBS = 

# CHANGEME: Additional flags for compilation (e.g., include flags)
ADDINCFLAGS =

# CHANGEME: Directory to the sources for the (example) problem definition
# files
COINORPATH = $(HOME)/COIN-OR
CPLEXDIR = $(HOME)/CPLEX128/cplex

##########################################################################
#  Usually, you don't have to change anything below.  Note that if you   #
#  change certain compiler options, you might have to recompile the      #
#  COIN package.                                                         #
##########################################################################

# C++ Compiler command
CXX = g++

# C++ Compiler options
CXXFLAGS = -O3 -pipe -DNDEBUG -Wparentheses -Wreturn-type -Wcast-qual -Wall -Wpointer-arith -Wwrite-strings -Wconversion -Wno-unknown-pragmas -Wno-long-long   -DSMI_BUILD

# Stochastic data directory
# CXXFLAGS += -DDATASTOCHASTICDIR=/homes/bcdandurand/COIN-OR/Data/Stochastic

# additional C++ Compiler options for linking
CXXLINKFLAGS =  -Wl,--rpath -Wl,${COINORPATH}/lib

INCL = -I$(COINORPATH)/include/coin 
INCL += -I$(COINORPATH)/include/coin-or 
CPLEXINC=$(CPLEXDIR)/include/ilcplex
INCL += -I$(CPLEXINC)
INCL += -I./include
INCL += $(ADDINCFLAGS)

# Linker flags
#LIBS = -L$(COINORPATH)/lib -lSmi -lFlopCpp -lOsiCpx -lOsiClp -lOsi -lClp -lCoinUtils 
LIBS = -L$(COINORPATH)/lib -lSmi -lOsiCpx -lOsiClp -lOsi -lClp -lCoinUtils 
CPLEXLIB=$(CPLEXDIR)/bin/x86-64_linux
LIBS += -L$(CPLEXLIB) -lcplex1280


all: $(EXE)

.SUFFIXES: .cpp .o .obj

$(EXE): $(OBJS)
	$(CXX) $(CXXLINKFLAGS) $(CXXFLAGS) -o $(EXE) $(OBJS) $(LIBS) $(ADDLIBS)

$(OBJS): $(SRC) 
	$(CXX) $(CXXFLAGS) $(INCL) -c $(SRC) -o $(OBJS)

clean: 
	rm -rf $(EXE) $(OBJS)
