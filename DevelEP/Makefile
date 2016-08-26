CC=g++ 
CFLAGS= -ggdb -c -I /usr/include/eigen3 -std=c++0x -fopenmp -Ofast
LDFLAGS= -std=c++0x -fopenmp -Ofast
LDFLAGS2= -Bstatic 
SOURCES= stdafx.cpp ellenGP.cpp AgeBreed.cpp AgeFitSurvival.cpp AgeFitGenSurvival.cpp ApplyGenetics.cpp BruteforceAlgo.cpp Crossover.cpp Datapoint.cpp DC.cpp EpiHC.cpp EpiMut.cpp Eqn2Line.cpp Fitness.cpp FitnessEstimator.cpp general_fns.cpp Generation.cpp HillClimb.cpp InitPop.cpp LexicaseSelect.cpp Line2Eqn.cpp Mutate.cpp NondominatedsortAlgo.cpp pareto.cpp pareto_fc.cpp ParetoSurvival.cpp Prune.cpp runEllenGP.cpp StablesortAlgo.cpp strdist.cpp Tournament.cpp 
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=ellenGP
BIN=./bin/

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) $(LDFLAGS2) -o $@ 

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
