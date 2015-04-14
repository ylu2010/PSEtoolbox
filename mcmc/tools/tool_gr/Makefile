# the compilers
CXX      =  g++		# sets the Cpp compiler
CC       =  gcc         # sets the C-compiler

# compiling options
OPT +=

OPTIMIZE =  -O3 -Wall -g        # sets optimization and warning flags

# the executable
EXEC   = GRanalyze

# the objects
OBJS   = GR_analyze.o GR.o grubbs.o mix.o student_t.o

# libraries
LIBS   =   -lm  
INCL   =   

# flags for c++ and c
CXXFLAGS = $(OPT) $(OPTIMIZE) $(INCL)  # -c is automatically included
CFLAGS =   $(OPT) $(OPTIMIZE) $(INCL) -std=gnu99  # -c is automatically included


all: $(EXEC)
# c++ compiler does the link
$(EXEC): $(OBJS)
	$(CXX) $(OBJS) $(LIBS)   -o  $(EXEC)
    
clean:
	rm -f $(OBJS) 
tidy:
	rm -f $(OBJS) $(EXEC)

install:
	cp -r $(EXEC) ../../bin/
