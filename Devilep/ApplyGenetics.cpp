#include "stdafx.h"
#include "pop.h"
#include "params.h"
#include "state.h"
#include "rnd.h"
#include "data.h"
#include "FitnessEstimator.h"
#include "Generationfns.h"
#include "Fitness.h"

void ApplyGenetics(vector<ind>& pop,vector<unsigned int>& parloc,params& p,vector<Randclass>& r,data& d,state& s,FitnessEstimator& FE)
{
	//boost::progress_timer timer;
	float choice;
	
	vector <ind> tmppop;
	///if(!p.parallel)
	//{
	unsigned int numits=0;
	if (p.sel==4) //age-fitness pareto
	{
		for (unsigned int i=0;i<pop.size();i++)
			parloc.push_back(i);
	}
	std::random_shuffle(parloc.begin(),parloc.end(),r[omp_get_thread_num()]);
	while(tmppop.size()<pop.size())
	{
		
		choice = r[omp_get_thread_num()].rnd_flt(0,1);

		if(choice < p.rep_wheel[0]) // reproduction
		{
			//pick parent from parloc 
				
			//push parent into tmppop
			tmppop.push_back(pop.at(parloc[numits]));
			//update ages
			tmppop.back().age=pop.at(parloc[numits]).age;
			pop.at(parloc[numits]).age++;
			numits++;
		}
		else if (choice<p.rep_wheel[1] && parloc.size()>numits+1) // crossover, as long as there are at least two parents left
		{			
			// only breed if parents have different fitness 
			if(pop.at(parloc[numits]).fitness!=pop.at(parloc[numits+1]).fitness) 
			{				
				//cross parents to make two kids
				Crossover(pop.at(parloc[numits]),pop.at(parloc[numits+1]),tmppop,p,r);
				//update ages
				tmppop.back().age=max(pop.at(parloc[numits]).age,pop.at(parloc[numits+1]).age);
				tmppop.at(tmppop.size()-2).age=max(pop.at(parloc[numits]).age,pop.at(parloc[numits+1]).age);
				pop.at(parloc[numits]).age++;
				pop.at(parloc[numits+1]).age++;
				numits+=2;
				//push kids into tmppop
			}
			else
			{
				Mutate(pop.at(parloc[numits]),tmppop,p,r,d);
				Mutate(pop.at(parloc[numits+1]),tmppop,p,r,d);
				//update ages
				tmppop.at(tmppop.size()-2).age = pop.at(parloc[numits]).age;
				tmppop.back().age = pop.at(parloc[numits+1]).age;
				
				pop.at(parloc[numits]).age++;
				pop.at(parloc[numits+1]).age++;

				numits+=2;
			}
		}
		else if (choice < p.rep_wheel[2]) //mutation
		{
			Mutate(pop.at(parloc[numits]),tmppop,p,r,d);
			//update ages
			tmppop.back().age=pop.at(parloc[numits]).age;
			pop.at(parloc[numits]).age++;
			numits++;
		}

		
	}

	while(tmppop.size()>pop.size())
		tmppop.erase(tmppop.end()-1);

	if (p.sel==4 || p.lexage) // insert new pop into old pop
	{
		//get tmppop fitness 
		Fitness(tmppop,p,d,s,FE);
		// genetic stats
		s.setCrossPct(tmppop);
		//add tmppop to pop
		pop.insert(pop.end(),tmppop.begin(),tmppop.end());
	}
	else //replace pop with tmppop
		pop.swap(tmppop);

}