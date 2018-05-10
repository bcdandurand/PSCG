# Makefile

#LOCALDIRROOT = /lustre/pRMIT0153/infrastructure/SMI_Stuff
#LOCALDIRROOT = /short/ka3/comp_infrastructure
LOCALDIRROOT = /homes/bcdandurand/comp_infrastructure
#SMIDIR = $(LOCALDIRROOT)/CoinSMI
SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic
CPLEXDIR      = $(LOCALDIRROOT)/CPLEX/cplex
CONCERTDIR    = $(LOCALDIRROOT)/CPLEX/concert
#OPENMPIDIR = /system/dccapps/openmpi/1.6.2
#OPENMPIDIR = /apps/openmpi/1.8.8
#OPENMPIDIR = /soft/openmpi/1.8.8
OPENMPIDIR = /usr/lib/openmpi/lib
#LIBNUMADIR = /usr/lib/x86_64-linux-gnu
#MPIDIRDSP = /nfs2/b216449/mpich-install/lib
#DSPDIR = /homes/bcdandurand/DSP
#DSPDIR = $(LOCALDIRROOT)/DSP
#DSPDIR = $(LOCALDIRROOT)/DSP-stable/DSP
DSPDIR = $(LOCALDIRROOT)/DSP
SRC = ./src
OBJDIR = ./obj
OBJS = $(OBJDIR)/PSCGMain.o $(OBJDIR)/$(ALG).o $(OBJDIR)/$(SOLVER).o  
OBJSLIB = $(OBJDIR)/$(ALG).o $(OBJDIR)/$(SOLVER).o 
#OBJS_serial = $(OBJDIR)/PSCGMain_serial.o $(OBJDIR)/$(ALG)_serial.o $(OBJDIR)/$(SOLVER)_serial.o $(OBJDIR)/ProblemDataBodur_serial.o 
LIBDIR = ./lib
BIN = ./bin
#CPLEXDIR      = /usr/local/ibm/ILOG/CPLEX_Studio125/cplex
#CONCERTDIR    = /usr/local/ibm/ILOG/CPLEX_Studio125/concert


# C++ Compiler command
#CXX = g++
#Need to use gcc/4.6.4 with mpicxx
CXX = mpicxx 
#GCCOPENMPI = pgcpp
#Need to use gcc/4.6.4
#GCCOPENMPI = g++

# C++ Compiler options
#CXXFLAGS = -O3 -g -std=c++11 -fpic -Wall
CXXFLAGS = -O3 -g -std=c++11 -fpic
#CXXFLAGS = -g -O3 -pipe -DNDEBUG -pedantic-errors -Wparentheses -Wreturn-type -Wcast-qual -Wall -Wpointer-arith -Wwrite-strings -Wconversion -Wno-unknown-pragmas -Wno-long-long   -DBCPS_BUILD 
#CXXFLAGS = -O3 

# Stochastic data directory

# additional C++ Compiler options for linking
#CXXLINKFLAGS =  -Wl,--rpath -Wl,$(SMIDIR)/lib -Wl,--rpath -Wl,$(LIBDIR)
#CXXLINKFLAGS =  -Wl,--rpath -Wl,$(DSPDIR)/build/lib -Wl,--rpath -Wl,$(SMIDIR)/lib -Wl,--rpath -Wl,$(LIBDIR)
#CXXLINKFLAGS =  -Wl,-rpath,$(DSPDIR)/build/lib -Wl,-rpath,$(SMIDIR)/lib -Wl,--rpath -Wl,$(LIBDIR)
#CCOPT = -m64 -O -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG -DIL_STD -std=c++0x

# Include directories (we use the CYGPATH_W variables to allow compilation with Windows compilers)
INCLTCLAP = -I $(LOCALDIRROOT) 
INCLDSP = -I$(DSPDIR)/build/include/coin -I$(DSPDIR)/src/Model -I$(DSPDIR)/src -I$(DSPDIR)/src/Utility
#INCLSMI = -I$(SMIDIR)/include/coin 
INCLCPLEX = -I$(CPLEXDIR)/include -I$(CONCERTDIR)/include 
INCLMPI = -I$(OPENMPIDIR)/include

INCL = $(INCLTCLAP) $(INCLDSP) $(INCLCPLEX)
#INCL = $(INCLTCLAP) $(INCLDSP) $(INCLSMI) $(INCLCPLEX)
#INCL = -I$(SMIDIR)/include/coin -I$(CPLEXDIR)/include -I$(CONCERTDIR)/include $(INCLDSP)
#INCL += $(ADDINCFLAGS)

#DSPLIBS = -L$(DSPDIR)/build/lib -lDsp -D_GLIBCXX_USE_CXX11_ABI=0
DSPLIBS = -Wl,-rpath,$(DSPDIR)/build/lib -L$(DSPDIR)/build/lib -lDsp
#DSPLIBS =  

#COINORLIBS = -L$(SMIDIR)/lib -lSmi -lOsiCpx -lClp -lClpSolver -lCoinUtils -lOsi -lOsiClp -lOsiCommonTests  
#ALPSLIBS = -L$(ALPSDIR)/lib -lAlps -lBcps 
#ALPSLIBS = -L$(ALPSDIR)/lib -lBcps 

#CPLEXLIBR = -DIL_STD -DILOSTRICTPOD $(CPPFLAGS) $(LDFLAGS) -L$(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT) \

CPLEXLIBR = -DIL_STD $(CPPFLAGS) $(LDFLAGS) -L$(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT) \
	-lilocplex -lcplex -L$(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT) -lconcert 
PSCGLIBR = -Wl,-rpath,$(LIBDIR) -L$(LIBDIR) -lPSCG
OTHERLIBS = -lstdc++ -lm -lpthread -lz -lbz2

ALG = PSCG
#mpicxx -o $(BIN)/$(ALG) $(OBJS) $(CPLEXLIBR) $(COINORLIBS) $(DSPLIBS) $(OTHERLIBS) $(CXXFLAGS) $(CXXLINKFLAGS) -D USING_MPI 
SOLVER = PSCGScen

SRCFILES = $(SRC)/PSCGMain.cpp $(SRC)/$(ALG).cpp $(SRC)/$(ALG).h $(SRC)/$(SOLVER).cpp $(SRC)/$(SOLVER).h  
	
	
all: all_parallel 
#all_serial	

all_parallel: $(ALG)

#all_serial: $(ALG)_Serial


$(ALG): $(SRCFILES) 

	$(CXX) -c -o $(OBJDIR)/PSCGMain.o -D USING_MPI $(INCLMPI) $(INCL) $(SRC)/PSCGMain.cpp $(DSPLIBS) $(CPLEXLIBR) $(OTHERLIBS) $(CXXFLAGS) 
	
	$(CXX) -c -o $(OBJDIR)/$(ALG).o -D USING_MPI $(INCLMPI) $(INCL) $(SRC)/$(ALG).cpp $(DSPLIBS) $(CPLEXLIBR) $(OTHERLIBS) $(CXXFLAGS) 
	
	$(CXX) -c -o $(OBJDIR)/$(SOLVER).o $(SRC)/$(SOLVER).cpp $(CPLEXLIBR) $(OTHERLIBS) $(CXXFLAGS) $(INCL)
	
	mpicxx -shared -o $(LIBDIR)/lib$(ALG).so $(OBJSLIB) $(DSPLIBS) $(CPLEXLIBR) $(OTHERLIBS) $(CXXFLAGS) $(CXXLINKFLAGS) -D USING_MPI 

	mpicxx -o $(BIN)/$(ALG) $(OBJDIR)/PSCGMain.o $(PSCGLIBR) $(DSPLIBS) $(CPLEXLIBR) $(OTHERLIBS) $(CXXFLAGS) $(CXXLINKFLAGS) -D USING_MPI 
	
#$(ALG)_Serial: $(SRCFILES)
#	$(CXX) -c -o $(OBJDIR)/PSCGMain_serial.o $(SRC)/PSCGMain.cpp $(CPLEXLIBR) $(DSPLIBS) $(OTHERLIBS) $(CXXFLAGS) $(INCL)
#	
#	$(CXX) -c -o $(OBJDIR)/$(ALG)_serial.o $(SRC)/$(ALG).cpp $(CPLEXLIBR) $(DSPLIBS) $(OTHERLIBS) $(CXXFLAGS) $(INCL)
#	
#	$(CXX) -c -o $(OBJDIR)/$(SOLVER)_serial.o $(SRC)/$(SOLVER).cpp $(CPLEXLIBR) $(OTHERLIBS) $(CXXFLAGS) $(INCL)
#	
#	$(CXX) -c -o $(OBJDIR)/ProblemDataBodur_serial.o $(SRC)/ProblemDataBodur.cpp $(CPLEXLIBR) $(OTHERLIBS) $(CXXFLAGS) $(INCL)
#	
#	$(CXX) -o $(BIN)/$(ALG)_serial $(OBJS_serial) $(CPLEXLIBR) $(COINORLIBS) $(DSPLIBS) $(OTHERLIBS) $(CXXFLAGS) $(CXXLINKFLAGS) 

clean:
	rm -f $(BIN)/$(ALG) $(BIN)/$(ALG)_serial $(OBJDIR)/*.o $(LIBDIR)/*
