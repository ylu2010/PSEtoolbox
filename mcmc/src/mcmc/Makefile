# the compilers
CXX      =  mpicxx		# sets the Cpp compiler
CC       =  mpicc               # sets the C-compiler
FC	 =  mpif90		# sets the fortran compiler

# compiling options
OPT +=

OPTIMIZE =  -O3 -g        # sets optimization and warning flags

# the executable
EXEC   = mcmc_test

# the objects
OBJS   = main.o \
         run_info.o prior.o proposal.o chains.o \
         start_mcmc.o checkpoint_io.o mpi_rng.o \
         parameter.o \
         mcmc_run.o chain_log_io.o mcmc_log.o \
         mcmc_algorithm_metropolis_hastings.o \
         mcmc_algorithm_differential_evolution.o \
         mcmc_algorithm_tempered_transition.o \
         probability.o

GSL_INCL = -I/usr/include
GSL_LIBS = -L/usr/lib

# libraries
LIBS   =   -lm  $(GSL_LIBS) -lgsl -lgslcblas
INCL   =   $(GSL_INCL)

# auxiliary code
include Makefile_Aux.h

# define your project libraries, compile options and object files here
include Makefile_Project.h

# flags for c++ and c and fortran
CXXFLAGS = $(OPT) $(OPTIMIZE) $(INCL)  # -c is automatically included
CFLAGS =   $(OPT) $(OPTIMIZE) $(INCL) -std=gnu99  # -c is automatically included
FFLAGS =   -O3 -g

%.o: %.f90
	$(FC) $(FFLAGS) -c $*.f90

all: $(EXEC)
# c++ compiler does the link
$(EXEC): $(OBJS)
	$(CXX) $(OBJS) $(LIBS)   -o  $(EXEC)
    
clean:
	rm -f $(OBJS) 
tidy:
	rm -f $(OBJS) $(EXEC)
