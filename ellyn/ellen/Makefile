CC=g++
CFLAGS= -c -fPIC -I /usr/include/eigen3 -I /usr/include/python3.5m -std=c++0x -fopenmp -Ofast
LDFLAGS= -std=c++0x -fopenmp -Ofast -shared -lpython3.5m -L/usr/include/boost_1_61_0/stage/lib/
#-L/usr/lib/x86_64-linux-gnu/
# LDFLAGS2= /usr/include/boost_1_61_0/stage/lib/libboost_python3.a
#LDFLAGS2= include/libboost_python3.a
LDFLAGS2= -lboost_python3 -Wl,-rpath,'/usr/include/boost_1_61_0/stage/lib/'
#-L/usr/lib/x86_64-linux-gnu
SOURCES= stdafx.cpp AgeBreed.cpp AgeFitSurvival.cpp AgeFitGenSurvival.cpp ApplyGenetics.cpp BruteforceAlgo.cpp Crossover.cpp Datapoint.cpp DC.cpp EpiHC.cpp EpiMut.cpp Eqn2Line.cpp Fitness.cpp FitnessEstimator.cpp general_fns.cpp Generation.cpp HillClimb.cpp InitPop.cpp LexicaseSelect.cpp Line2Eqn.cpp Mutate.cpp NondominatedsortAlgo.cpp pareto.cpp pareto_fc.cpp ParetoSurvival.cpp Prune.cpp runEllenGP.cpp StablesortAlgo.cpp StochasticGradient.cpp strdist.cpp Tournament.cpp
#SOURCES = $(wildcard *.cpp)
#OBJECTS=$(SOURCES:.cpp=.o)
OBJECTS = $(patsubst %.cpp, $(BIN)/%.o, $(SOURCES))
BIN=./lib
EXECUTABLE=lib/elgp.so


all: lib_dir program

debug:	CFLAGS += -ggdb
debug: lib_dir program

lib_dir:
	mkdir -p $(BIN)

program: $(SOURCES) $(EXECUTABLE)

link: cleanexe $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) $(LDFLAGS2) -o $@

$(BIN)/%.o: %.cpp
	$(CC) $(CFLAGS) $< -o $@

.PHONY : clean
clean :
	- rm $(OBJECTS)
	- rm $(EXECUTABLE)

cleanexe :
	- rm $(EXECUTABLE)
