CC=icpc
CFLAGS= -c -I /home/quarrian/boost/boost_1_56_0 -std=c++11 -openmp -pch-use pch/stdafx.pchi 
LDFLAGS= -I /home/quarrian/boost/boost_1_56_0 -std=c++11 -openmp 
LDFLAGS2= -Bstatic -L/home/quarrian/boost/boost_1_56_0/stage/lib -lboost_regex 
SOURCES= stdafx.cpp Develep.cpp AgeBreed.cpp AgeFitSelect.cpp ApplyGenetics.cpp BruteforceAlgo.cpp Crossover.cpp Datapoint.cpp DC.cpp EpiHC.cpp Eqn2Line.cpp Fitness.cpp FitnessEstimator.cpp general_fns.cpp Generation.cpp HillClimb.cpp InitPop.cpp LexicaseSelect.cpp Line2Eqn.cpp Mutate.cpp NondominatedsortAlgo.cpp pareto.cpp runDevelep.cpp StablesortAlgo.cpp strdist.cpp Tournament.cpp 
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=DevelEP

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) $(LDFLAGS2) -o $@ 

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
