CC=icpc
CFLAGS= -c -I /home/quarrian/boost/boost_1_56_0 -I /home/quarrian/ellen/DevelEP -std=c++11 -openmp -pch-use pch/stdafxRT.pchi 
LDFLAGS= -I /home/quarrian/boost/boost_1_56_0 -I /home/quarrian/ellen/DevelEP -std=c++11 -openmp 
LDFLAGS2= -Bstatic -L/home/quarrian/boost/boost_1_56_0/stage/lib -lboost_regex 
SOURCES= stdafxRT.cpp RunTrials.cpp ../DevelEP/stdafx.cpp ../DevelEP/AgeBreed.cpp ../DevelEP/AgeFitSelect.cpp ../DevelEP/ApplyGenetics.cpp ../DevelEP/BruteforceAlgo.cpp ../DevelEP/Crossover.cpp ../DevelEP/Datapoint.cpp ../DevelEP/DC.cpp ../DevelEP/EpiHC.cpp ../DevelEP/Eqn2Line.cpp ../DevelEP/Fitness.cpp ../DevelEP/FitnessEstimator.cpp ../DevelEP/general_fns.cpp ../DevelEP/Generation.cpp ../DevelEP/HillClimb.cpp ../DevelEP/InitPop.cpp ../DevelEP/LexicaseSelect.cpp ../DevelEP/Line2Eqn.cpp ../DevelEP/Mutate.cpp ../DevelEP/NondominatedsortAlgo.cpp ../DevelEP/pareto.cpp ../DevelEP/runDevelep.cpp ../DevelEP/StablesortAlgo.cpp ../DevelEP/strdist.cpp ../DevelEP/Tournament.cpp 
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=RunTrials

pchmake: stdafxRT.cpp
	$(CC) -pch-create pch/stdafxRT.pchi stdafxRT.cpp -c -I /home/quarrian/boost/boost_1_56_0 -I /home/quarrian/ellen/DevelEP -std=c++11 -openmp

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) $(LDFLAGS2) -o $@ 

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
