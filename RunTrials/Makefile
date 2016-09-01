CC=g++
CFLAGS = -c -I /usr/include/eigen3 -I ../DevelEP -fopenmp -std=c++0x -Ofast
LDFLAGS= -I ../DevelEP -fopenmp -std=c++0x -Ofast
LDFLAGS2= -Bstatic
SOURCES= stdafxRT.cpp RunTrials.cpp
PATH_ELLEN = ../DevelEP/bin
OBJ_ELLEN = stdafx.o AgeBreed.o AgeFitSurvival.o AgeFitGenSurvival.o ApplyGenetics.o BruteforceAlgo.o Crossover.o Datapoint.o DC.o EpiHC.o EpiMut.o Eqn2Line.o Fitness.o FitnessEstimator.o general_fns.o Generation.o HillClimb.o InitPop.o LexicaseSelect.o Line2Eqn.o Mutate.o NondominatedsortAlgo.o pareto.o pareto_fc.o ParetoSurvival.o Prune.o runEllenGP.o StablesortAlgo.o strdist.o Tournament.o
OBJECTS = $(patsubst %.cpp, $(BIN)/%.o, $(SOURCES))$(patsubst %, $(PATH_ELLEN)/%, $(OBJ_ELLEN))
#OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=bin/RunTrials
BIN=./bin

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) $(LDFLAGS2) -o $@ 


$(BIN)/%.o: %.cpp
	$(CC) $(CFLAGS) $< -o $@
dir:
	mkdir -p $(BIN)


