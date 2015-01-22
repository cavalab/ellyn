#include "stdafx.h"
#include "pop.h"
#include "params.h"
#include "rnd.h"

void AgeFitSurvival(vector<ind>& pop,params& p,vector<Randclass>& r)
{
	//boost::progress_timer timer;
	vector<float> fit(p.tourn_size);
	vector<int> age(p.tourn_size);
	vector<int> fitindex(p.tourn_size);
	//float minfit=p.max_fit;
	//float minage=p.g; // maximum age = total max generations
	int loser;
	int counter =0;
	bool draw=true;
	int popsize = (int)floor((float)pop.size()/2);
	while (pop.size()>popsize && counter<p.popsize*10)
	{
		for (int j=0;j<2;++j)
		{
			fitindex.at(j)=r[omp_get_thread_num()].rnd_int(0,pop.size()-1);
			fit.at(j) = pop.at(fitindex.at(j)).fitness;
			age.at(j) = pop.at(fitindex.at(j)).age;
		}
		if(fit[0]<=fit[1] && age[0]<=age[1])
		{
			if(age[0]<age[1] || fit[0]<fit[1])
			{
				draw=false;
				loser = fitindex[1];
			}
		}
		else if(fit[1]<=fit[0] && age[1]<=age[0]) //fit0 > fit1 or age0 >age1
		{
			if(age[1]<age[0] || fit[1]<fit[0])
			{
				draw=false;
				loser = fitindex[0];
			}
		}

		// delete losers from population
		if(!draw){
			pop.at(loser).swap(pop.back());
			//std::swap(pop.at(loser),pop.back());
			pop.pop_back();
			//pop.erase(pop.begin()+loser);
		}

		//minfit=p.max_fit;
		//minage=p.g;
		draw=true;
		counter++;
	}
	if (pop.size()>popsize)	{
		sort(pop.begin(),pop.end(),SortAge());
		stable_sort(pop.begin(),pop.end(),SortFit());
	}
	while(pop.size()>popsize){
		//pop.erase(pop.end()-1);
		pop.pop_back();
	}
	

}