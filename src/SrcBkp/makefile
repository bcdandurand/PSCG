# Makefile

#LOCALDIRROOT = /lustre/pRMIT0153/infrastructure/SMI_Stuff
#LOCALDIRROOT = /short/ka3/comp_infrastructure
LOCALDIRROOT = /homes/bcdandurand/comp_infrastructure
SMIDIR = $(LOCALDIRROOT)/CoinSMI
SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic
CPLEXDIR      = $(LOCALDIRROOT)/CPLEX/cplex
CONCERTDIR    = $(LOCALDIRROOT)/CPLEX/concert
#OPENMPIDIR = /system/dccapps/openmpi/1.6.2
#OPENMPIDIR = /apps/openmpi/1.8.8
#OPENMPIDIR = /soft/openmpi/1.8.8
OPENMPIDIR = /usr/lib/openmpi/lib
LIBNUMADIR = /usr/lib/x86_64-linux-gnu
MPIDIRDSP = /nfs2/b216449/mpich-install/lib
DSPDIR = /homes/bcdandurand/DSP
SRC = ./src
OBJDIR = ./obj
OBJS = $(OBJDIR)/$(ALG).o $(OBJDIR)/$(SOLVER).o $(OBJDIR)/ProblemDataBodur.o 
OBJS_serial = $(OBJDIR)/$(ALG)_serial.o $(OBJDIR)/$(SOLVER)_serial.o $(OBJDIR)/ProblemDataBodur_serial.o 
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
CXXFLAGS = -O3 -g
#CXXFLAGS = -O3 

# Stochastic data directory

# additional C++ Compiler options for linking
CXXLINKFLAGS =  -Wl,--rpath -Wl,$(SMIDIR)/lib
#CXXLINKFLAGS =  -Wl,--rpath -Wl,$(SMIDIR)/lib -Wl,--rpath -Wl,$(DSPDIR)/src -Wl,--rpath -Wl,$(OPENMPIDIR)
#CXXLINKFLAGS = -Wl,--rpath -Wl,$(SMIDIR)/lib -Wl,--rpath -Wl,$(DSPDIR)/src -Wl,--rpath -Wl,$(OPENMPIDIR) -Wl,--rpath -Wl,$(MPIDIRDSP) -Wl,--rpath -Wl,$(LIBNUMADIR)
#CCOPT = -m64 -O -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG -DIL_STD -std=c++0x

# Include directories (we use the CYGPATH_W variables to allow compilation with Windows compilers)
INCLPSCG = -I ./include
INCLDSP = -I$(DSPDIR)/src/Model -I$(DSPDIR)/src
INCLSMI = -I$(SMIDIR)/include/coin 
INCLCPLEX = -I$(CPLEXDIR)/include -I$(CONCERTDIR)/include 
INCLMPI = -I$(OPENMPIDIR)/include

INCL = $(INCLPSCG) $(INCLDSP) $(INCLSMI) $(INCLCPLEX)
#INCL = -I$(SMIDIR)/include/coin -I$(CPLEXDIR)/include -I$(CONCERTDIR)/include $(INCLDSP)
#INCL += $(ADDINCFLAGS)

#DSPLIBS = -L$(DSPDIR)/build/lib -lDsp -D_GLIBCXX_USE_CXX11_ABI=0
DSPLIBS = -L$(DSPDIR)/build/lib -lDsp 
#DSPLIBS =  

COINORLIBS = -L$(SMIDIR)/lib -lSmi -lOsiCpx -lClp -lClpSolver -lCoinUtils -lOsi -lOsiClp -lOsiCommonTests  

CPLEXLIBR = -DIL_STD -DILOSTRICTPOD $(CPPFLAGS) $(LDFLAGS) -L$(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT) \
	-lilocplex -lcplex -L$(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT) -lconcert 
OTHERLIBS = -lstdc++ -lm -lpthread
#OTHERLIBS = -lstdc++ -lgcc_s -lc -lgcc -lm -lpthread

ALG = ParallelSCG
#SOLVER = CPLEXsolverSCG
SOLVER = CPLEXSCGSolver
	
all: all_parallel all_serial	

all_parallel: $(ALG)

all_serial: $(ALG)_Serial

$(ALG): $(SRC)/$(ALG).cpp $(SRC)/$(ALG).h $(SRC)/$(SOLVER).cpp $(SRC)/$(SOLVER).h $(SRC)/ProblemDataBodur.cpp 
	$(CXX) -c -o $(OBJDIR)/$(ALG).o -D USING_MPI $(INCLMPI) $(INCL) $(SRC)/$(ALG).cpp $(CPLEXLIBR) $(DSPLIBS) $(OTHERLIBS) $(CXXFLAGS) 
	
	$(CXX) -c -o $(OBJDIR)/$(SOLVER).o $(SRC)/$(SOLVER).cpp $(CPLEXLIBR) $(OTHERLIBS) $(CXXFLAGS) $(INCL)
	
	$(CXX) -c -o $(OBJDIR)/ProblemDataBodur.o $(SRC)/ProblemDataBodur.cpp $(CPLEXLIBR) $(OTHERLIBS) $(CXXFLAGS) $(INCL) 
	
	#$(CXX) -c -o $(OBJDIR)/ProblemDataSMPS.o $(SRC)/ProblemDataSMPS.cpp $(CPLEXLIBR) $(OTHERLIBS) $(CXXFLAGS) $(INCL) 
	
	mpicxx -o $(BIN)/$(ALG) $(OBJS) $(CPLEXLIBR) $(COINORLIBS) $(DSPLIBS) $(OTHERLIBS) $(CXXFLAGS) $(CXXLINKFLAGS) -D USING_MPI 
	
$(ALG)_Serial: $(SRC)/$(ALG).cpp $(SRC)/$(ALG).h $(SRC)/$(SOLVER).cpp $(SRC)/ProblemDataBodur.cpp
	$(CXX) -c -o $(OBJDIR)/$(ALG)_serial.o $(SRC)/$(ALG).cpp $(CPLEXLIBR) $(DSPLIBS) $(OTHERLIBS) $(CXXFLAGS) $(INCL)
	
	$(CXX) -c -o $(OBJDIR)/$(SOLVER)_serial.o $(SRC)/$(SOLVER).cpp \
	$(CPLEXLIBR) $(OTHERLIBS) $(CXXFLAGS) $(INCL)
	
	$(CXX) -c -o $(OBJDIR)/ProblemDataBodur_serial.o $(SRC)/ProblemDataBodur.cpp \
	$(CPLEXLIBR) $(OTHERLIBS) $(CXXFLAGS) $(INCL)
	
	#$(CXX) -c -o $(OBJDIR)/ProblemDataSMPS_serial.o $(SRC)/ProblemDataSMPS.cpp \
	#$(CPLEXLIBR) $(OTHERLIBS) $(CXXFLAGS) $(INCL)
	
	$(CXX) -o $(BIN)/$(ALG)_serial $(OBJS_serial) \
	$(CPLEXLIBR) $(COINORLIBS) $(DSPLIBS) $(OTHERLIBS) $(CXXFLAGS) $(CXXLINKFLAGS) 
clean:
	rm -f $(BIN)/$(ALG) $(BIN)/$(ALG)_serial $(OBJDIR)/*.o
