CXX = g++
HOME= /usr/local/include
LIB_HOME = ../../../flowstar-toolbox
LIBS = -lmpfr -lgmp -lgsl -lgslcblas -lm -lglpk -lflowstar
CFLAGS = -I . -I $(HOME) -g -O3 -std=c++11
LINK_FLAGS = -g -L$(LIB_HOME) -L/usr/local/lib

all: p53

p53: p53.o
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)


%.o: %.cc
	$(CXX) -O3 -c $(CFLAGS) -o $@ $<
%.o: %.cpp
	$(CXX) -O3 -c $(CFLAGS) -o $@ $<
%.o: %.c
	$(CXX) -O3 -c $(CFLAGS) -o $@ $<


clean: 
	rm -f *.o p53