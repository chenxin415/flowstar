CXX = g++-8
include makefile.local
LIBS = -lmpfr -lgmp -lgsl -lgslcblas -lm -lglpk
CFLAGS = -I . -I $(GMP_HOME) -g -O3 -std=c++11
LINK_FLAGS = -g -L $(GMP_LIB_HOME)
OBJS = Interval.o Variables.o settings.o Matrix.o Geometry.o Constraints.o Continuous.o

all: flowstar lib

flowstar: $(OBJS) lex.yy.o modelParser.tab.o modelParser.o
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)

lib: $(OBJS) lex.yy.o modelParser.tab.o
	ar rcs libflowstar.a $^

%.o: %.cc
	$(CXX) -O3 -c $(CFLAGS) -o $@ $<
%.o: %.cpp
	$(CXX) -O3 -c $(CFLAGS) -o $@ $<
%.o: %.c
	$(CXX) -O3 -c $(CFLAGS) -o $@ $<

modelParser.tab.c: modelParser.y
	bison -d -v modelParser.y

lex.yy.c: modelLexer.l modelParser.tab.c
	flex modelLexer.l

clean: 
	rm -f flowstar *.o libflowstar.a *~ libflowstar.a *~ modelParser.tab.c modelParser.tab.h modelParser.output lex.yy.c
