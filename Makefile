CXX = g++
CXXFLAGS = -std=c++0x -Wall -Wextra -Wshadow -Werror -march=native -fopenmp -O3 -DNDEBUG

#-xHost chk what is equivalent in gcc
INCLUDES =
LDFLAGS =
LIBS =

# blas
#INCLUDES += -I/usr/lib64/atlas/include/
#LDFLAGS += -L/usr/lib64/atlas/
#LIBS += -lcblas -latlas

 #likwid
#CXXFLAGS += -DUSE_LIKWID -pthread
#INCLUDES += -I/usr/local/likwid-3.1.2/include/
#LDFLAGS += -L/usr/local/likwid-3.1.2/lib/
#LIBS += -llikwid

#TARGET = matmult
#OBJS = $(TARGET).o

all: mgsolve

mgsolve: Makefile main.o matrixread.o 
	$(CXX) $(CXXFLAGS) -o mgsolve matrixread.o main.o $(LDFLAGS) $(LIBS)

main.o: Makefile main.cpp Timer.h
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) main.cpp

matrixread.o: Makefile matrixread.cpp Timer.h  matrixread.h Makefile 
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) matrixread.cpp

clean:
	@$(RM) -rf *.o mgsolve *~ solution.txt
