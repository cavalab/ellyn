CC=icpc
CFLAGS=-c -Wall
LDFLAGS=
SOURCES=AgeBreed.cpp AgeFitSelect.cpp ApplyGenetics.cpp BruteforceAlgo.cc Crossover.cpp Datapoint.cc DC.cpp Develep.cpp EpiHC.cpp Eqn2Line.cpp Fitness.cpp FitnessEstimator.cpp general_fns.cpp Generation.cpp HillClimb.cpp InitPop.cpp LexicaseSelect.cpp Line2Eqn.cpp Mutate.cpp NondominatedsortAlgo.cc pareto.cc runDevelep.cpp StablesortAlgo.cc stdafx.cpp strdist.cpp Tournament.cpp 
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=DevelEP

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
