#include "stdafx.h"
#include "pop.h"
#include "state.h"
#include "Generationfns.h"

void ApplyGenetics(vector<ind>& pop,vector<unsigned int>& parloc,params& p,vector<Randclass>& r)
{
	//boost::progress_timer timer;
	float choice;
	
	vector <ind> tmppop;
	///if(!p.parallel)
	//{
	unsigned int numits=0;
	std::random_shuffle(parloc.begin(),parloc.end(),r[0]);
	while(tmppop.size()<pop.size())
	{
		
		choice = r[omp_get_thread_num()].rnd_flt(0,1);

		if(choice < p.rep_wheel[0]) // reproduction
		{
			//pick parent from parloc 
				
			//push parent into tmppop
			tmppop.push_back(pop.at(parloc[numits]));
			numits++;
		}
		else if (choice<p.rep_wheel[1] && parloc.size()>numits+1) // crossover, as long as there are at least two parents left
		{			
			// only breed if parents have different fitness 
			if(pop.at(parloc[numits]).fitness!=pop.at(parloc[numits+1]).fitness) 
			{				
				//cross parents to make two kids
				Crossover(pop.at(parloc[numits]),pop.at(parloc[numits+1]),tmppop,p,r);
				numits+=2;
				//push kids into tmppop
			}
		}
		else if (choice < p.rep_wheel[2]) //mutation
		{
			Mutate(pop.at(parloc[numits]),tmppop,p,r);
			numits++;
		}

		
	}

	while(tmppop.size()>pop.size())
		tmppop.erase(tmppop.end()-1);
	//replace pop with tmppop
	pop.swap(tmppop);
	//clear tmppop

}