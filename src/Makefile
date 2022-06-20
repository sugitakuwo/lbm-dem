CC=g++
export LANG=en_US.utf8
export LC_ALL=en_US.utf8
export CFLAGS=-Qunused-arguments
export CPPFLAGS=-Qunused-arguments


# Set optimization flags: cc is generic, gcc and icc are for PIII
CFLAGS_cc  = -O
CFLAGS_gcc = -O2 -funroll-loops -fomit-frame-pointer -malign-double -march=i686
CFLAGS_icc = -O -tpp6 -xK

# Set linker flags: i.e. mpe libraries or options like -static
LFLAGS_cc  =
LFLAGS_gcc =
LFLAGS_icc =

# Should not need to modify makefile below here
PCC    = mpicc
CFLAGS = $(CFLAGS_$(CC)) -c -fopenmp
LFLAGS = $(LFLAGS_$(CC)) -lm -fopenmp

# Should not need to modify makefile below here

C_OBJ  = DEM.o elmt.o main.o IO.o LB.o myvector.o node.o utils.o

.cpp.o:
	$(CC)  $(CFLAGS) $*.cpp

serial: $(C_OBJ)
	$(CC)  $(LFLAGS) -o susp_fs $(C_OBJ)


clean:
	rm -f *.o susp_fs

install:
	ln -sf ./main ./bin/main

versions:
	make serial
	make mpich
	make mplite

all:
	make clean
	make versions
	make install
