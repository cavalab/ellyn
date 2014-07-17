#include "stdafx.h"
#include "pop.h"
#include "state.h"
#include "rnd.h"
#include "Generationfns.h"
#include "FitnessEstimator.h"
#include "Fitness.h"

void AgeBreed(vector<ind>& pop,params& p,vector<Randclass>& r,data& d,state& s,FitnessEstimator& FE)
{
	//boost::progress_timer timer;
	float choice;
	
	vector <ind> tmppop;
	///if(!p.parallel)
	//{
	unsigned int numits=0;
	std::random_shuffle(pop.begin(),pop.end(),r[omp_get_thread_num()]);
	int counter=0;
	while(tmppop.size()<pop.size() && counter<100000)
	{
		
		choice = r[omp_get_thread_num()].rnd_flt(0,1);

		if (choice<p.rep_wheel[1] && pop.size()>numits+1) // crossover, as long as there are at least two parents left
		{			
			// only breed if parents have different fitness 
			if(pop.at(numits).fitness!=pop.at(numits+1).fitness) 
			{				
				//cross parents to make two kids and push kids into tmppop
				Crossover(pop.at(numits),pop.at(numits+1),tmppop,p,r);
				//update ages
				tmppop.back().age=max(pop.at(numits).age,pop.at(numits+1).age);
				tmppop.at(tmppop.size()-2).age=max(pop.at(numits).age,pop.at(numits+1).age);
				pop.at(numits).age++;
				pop.at(numits+1).age++;
				
				numits+=2;
				
			}
			else
			{
				Mutate(pop.at(numits),tmppop,p,r,d);
				Mutate(pop.at(numits+1),tmppop,p,r,d);
				//update ages
				tmppop.at(tmppop.size()-2).age = pop.at(numits).age;
				tmppop.back().age = pop.at(numits+1).age;
				
				pop.at(numits).age++;
				pop.at(numits+1).age++;

				numits+=2;
			}
		}
		else if (choice < p.rep_wheel[2]) //mutation
		{
			Mutate(pop.at(numits),tmppop,p,r,d);
			//update ages
			tmppop.back().age=pop.at(numits).age;
			pop.at(numits).age++;

			numits++;
		}
		counter++;
	}
	if (counter==100000)
		cout << "stuck in while loop in AgeBreed.cpp\n";

	while(tmppop.size()>pop.size())
		tmppop.erase(tmppop.end()-1);
	//get tmppop fitness 
	Fitness(tmppop,p,d,s,FE);
	// genetic stats
	s.setCrossPct(tmppop);
	//add tmppop to pop
	pop.insert(pop.end(),tmppop.begin(),tmppop.end());
	tmppop.clear();
}