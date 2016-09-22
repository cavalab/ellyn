CC=g++
CFLAGS= -c -I /usr/include/eigen3 -I /usr/include/python3.5m -std=c++0x -fopenmp -Ofast -fPIC
LDFLAGS= -std=c++0x -fopenmp -Ofast -fPIC -shared -L/usr/lib/x86_64-linux-gnu -Bstatic -lboost_python-py35 -lpython3.5m
LDFLAGS2= -Bstatic
SOURCES= stdafx.cpp ellenGP.cpp AgeBreed.cpp AgeFitSurvival.cpp AgeFitGenSurvival.cpp ApplyGenetics.cpp BruteforceAlgo.cpp Crossover.cpp Datapoint.cpp DC.cpp EpiHC.cpp EpiMut.cpp Eqn2Line.cpp Fitness.cpp FitnessEstimator.cpp general_fns.cpp Generation.cpp HillClimb.cpp InitPop.cpp LexicaseSelect.cpp Line2Eqn.cpp Mutate.cpp NondominatedsortAlgo.cpp pareto.cpp pareto_fc.cpp ParetoSurvival.cpp Prune.cpp runEllenGP.cpp StablesortAlgo.cpp strdist.cpp Tournament.cpp
#SOURCES = $(wildcard *.cpp)
#OBJECTS=$(SOURCES:.cpp=.o)
OBJECTS = $(patsubst %.cpp, $(BIN)/%.o, $(SOURCES))
BIN=./lib
EXECUTABLE=lib/elgp.so


all: lib_dir program
program: $(SOURCES) $(EXECUTABLE)

debug:	CFLAGS += -ggdb
debug: $(SOURCES) $(EXECUTABLE)

lib_dir:
	mkdir -p $(BIN)


$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) $(LDFLAGS2) -o $@

$(BIN)/%.o: %.cpp
	$(CC) $(CFLAGS) $< -o $@

.PHONY : clean
clean :
	- rm $(OBJECTS)
	- rm $(EXECUTABLE)
