src     = calc2DIR.cpp
exes    = calc2DIR.exe
CC      = g++
LIBDIR  = -L/${MKLROOT}/lib/intel64
LIBS    = -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl -lfftw3 -fopenmp

all: ${exes}

${exes}: ${src} calc2DIR.h
	$(CC) $(src) -o $(exes) $(LIBDIR) $(LIBS) -std=c++11 -fmax-errors=10 -O3

clean:
	rm calc2DIR.exe
