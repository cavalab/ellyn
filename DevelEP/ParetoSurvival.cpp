#include "stdafx.h"
#include "pop.h"
#include "params.h"
#include "rnd.h"
#include "pareto.h"
#include "state.h"

void ParetoSurvival(vector<ind>& pop,params& p,vector<Randclass>& r,state& s)
{
	//boost::progress_timer timer;
	vector<float> data;

	if(p.PS_sel==1)
		data.reserve(4);
	else if(p.PS_sel==2)
		data.reserve(6);
	
	vector<int> fitindex(2);
	//float minfit=p.max_fit;
	//float minage=p.g; // maximum age = total max generations
	int loser;
	int counter =0;
	bool draw=true;
	int popsize = (int)floor((float)pop.size()/2);
	//vector<ind> tourn_pop; 
	//tourn_pop.reserve(2);
	// remove all frontier pareto individuals from the tournament and give them a free pass:
	//vector<int> tournament; for(int i=0; i<pop.size(); ++i){tournament.push_back(i);}

	while (pop.size()>popsize && counter<p.popsize*10)//pow(float(p.popsize),2))
	{
		
		for (int j=0;j<2;++j)
		{
			fitindex.at(j)=r[omp_get_thread_num()].rnd_int(0,pop.size()-1);
			//tourn_pop.push_back(pop.at(fitindex.at(j)));
			data.push_back(pop.at(fitindex.at(j)).fitness);
			data.push_back(pop.at(fitindex.at(j)).age);
			if (p.PS_sel==2) data.push_back(pop.at(fitindex.at(j)).genty);
		}
		pareto2(data,pop.at(fitindex.at(0)).rank,pop.at(fitindex.at(1)).rank);
		//pareto(tourn_pop,data);
		//float genty1 = abs(tourn_pop[0].fitness-tourn_pop[0].fitness_v)/tourn_pop[0].fitness;
		//float genty2 = abs(tourn_pop[1].fitness-tourn_pop[1].fitness_v)/tourn_pop[1].fitness;
		if (pop.at(fitindex.at(0)).rank != pop.at(fitindex.at(1)).rank)
		{
			if (pop.at(fitindex.at(0)).rank==1)
				loser = fitindex[1];
			else
				loser = fitindex[0];
			draw=false;
		}
		
		// delete losers from population
		if(!draw){
			pop.at(loser).swap(pop.back());
			//std::swap(pop.at(loser),pop.back());
			pop.pop_back();
			//pop.erase(pop.begin()+loser);
		}

		draw=true;
		counter++;
		//tourn_pop.resize(0);
		data.resize(0);
	}
	if (pop.size()>popsize){
		//s.out << "warning: pareto survival tournament exceeded maximum comparisons.\n";
		if (p.PS_sel==1){
			sort(pop.begin(),pop.end(),SortAge());
			stable_sort(pop.begin(),pop.end(),SortFit());
		}
		else if(p.PS_sel==2){
			sort(pop.begin(),pop.end(),SortAge());
			stable_sort(pop.begin(),pop.end(),SortGenty());
			stable_sort(pop.begin(),pop.end(),SortFit());
		}
	}
	while(pop.size()>popsize){
		//pop.erase(pop.end()-1);
		pop.pop_back();
	}
	

}