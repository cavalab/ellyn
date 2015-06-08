#include "stdafx.h"
#include "pop.h"
#include "params.h"
#include "state.h"
#include "rnd.h"
#include "data.h"
#include "FitnessEstimator.h"
#include "Generationfns.h"
#include "Fitness.h"

void LexicaseSelect(vector<ind>& pop,vector<unsigned int>& parloc,params& p,vector<Randclass>& r,data& d)
{
	//boost::progress_timer timer;
	//vector<float> fitcompare;
	vector<int> pool; // pool from which to choose parent. normally it is set to the whole population.
	int numcases;
	if (p.EstimateFitness)
		numcases = p.FE_ind_size*p.train_pct;
	else
		numcases = d.vals.size()*p.train_pct;
	
	if (p.lexpool==1){
		for (int i=0;i<pop.size();++i){
			if(pop[i].error.size() == numcases)
				pool.push_back(i);
		}
	}
	vector <int> starting_pool = pool;
	vector <int> case_order;
	for(int i=0;i<pop[pool[0]].error.size();++i)
		case_order.push_back(i);
	assert(case_order.size()>0);
	float minfit=p.max_fit;
	vector<int> winner;
	bool draw=true;
	bool pass;
	int h;
	
	int tmp;
	for (int i=0;i<pop.size();++i)
	{
		//shuffle test cases
		std::random_shuffle(case_order.begin(),case_order.end(),r[omp_get_thread_num()]);
		
		if (p.lexpool!=1){
			for (int j=0;j<p.lexpool*pop.size();++j){
				tmp = r[omp_get_thread_num()].rnd_int(0,pop.size()-1);
				while(pop[tmp].error.empty())
					tmp = r[omp_get_thread_num()].rnd_int(0,pop.size()-1);
				pool.push_back(tmp);
			}
		}
		else 
			pool = starting_pool;

		pass=true;
		h=0;
		
		while ( pass && h<case_order.size())
		{
			//reset winner and minfit for next case
			winner.resize(0);
			minfit=p.max_fit;

			for (int j=0;j<pool.size();++j)
			{
				
				if (pop[pool[j]].error[case_order[h]]<minfit)
				{
					minfit=pop[pool[j]].error[case_order[h]];
					winner.resize(0);
					winner.push_back(pool[j]);
				}
				else if (pop[pool[j]].error[case_order[h]]==minfit)
				{
					winner.push_back(pool[j]);
				}
				//else{ // individual is worse than minimum
				//	swap(pop[pool[j]],pop.back());
				//	pop.pop_back();
				//}

			}
			if(winner.size()>1 && h<case_order.size()-1) 
			{
				pass=true;
				++h;
				//reduce pool to elite individuals on case case_oder[h]
				pool = winner; 
				
				/*for (int i=0; i<fitindex.size();++i)
					fitcompare.push_back(pop.at(fitindex[i]).error[case_order[h]]);*/
			}
			else 
				pass=false;
		}
		//get lowest fitness, return pop.at(fitindex.at(lowest)).id
		if(winner.size()>1)
			parloc[i]=winner[r[omp_get_thread_num()].rnd_int(0,winner.size()-1)];
		else if (winner.size()==1)
			parloc[i]=winner[0];
		/*else if (winner.empty())
			parloc[i]=pool[r[omp_get_thread_num()].rnd_int(0,pool.size()-1)];*/
		else
			cout << "??";

		minfit=p.max_fit;
	}//for (int i=0;i<pop.size();++i)

}