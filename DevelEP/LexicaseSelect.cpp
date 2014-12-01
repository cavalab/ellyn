#include "stdafx.h"
#include "pop.h"
#include "params.h"
#include "state.h"
#include "rnd.h"
#include "data.h"
#include "FitnessEstimator.h"
#include "Generationfns.h"
#include "Fitness.h"

void LexicaseSelect(vector<ind>& pop,vector<unsigned int>& parloc,params& p,vector<Randclass>& r)
{
	//boost::progress_timer timer;
	vector<float> fitcompare;
	vector<int> fitindex;
	vector <int> fitorder;
	for(int i=0;i<p.numcases;++i){	fitorder.push_back(i);}
	float minfit=p.max_fit;
	vector<int> winner;
	bool draw=true;
	bool pass;
	int h,lexpool;
	for (int i=0;i<pop.size();++i)
	{
		std::random_shuffle(fitorder.begin(),fitorder.end(),r[omp_get_thread_num()]);
		
		for (int j=0;j<p.lexpool;++j)
		{
			fitindex.push_back(r[omp_get_thread_num()].rnd_int(0,pop.size()-1));
			fitcompare.push_back(pop.at(fitindex.at(j)).fitlex[fitorder[0]]);
		}
		pass=true;
		h=0;
		lexpool=p.lexpool;
		while ( pass && h<fitorder.size())
		{
			for (int j=0;j<lexpool;++j)
			{
				
				if (fitcompare.at(j)<minfit)
				{
					minfit=fitcompare.at(j);
					winner.clear();
					winner.push_back(fitindex.at(j));
				}
				else if (fitcompare.at(j)==minfit)
				{
					winner.push_back(fitindex.at(j));
				}
			}
			if(winner.size()>1 && h<fitorder.size()-1) 
			{
				pass=true;
				++h;
				fitindex.clear();
				fitcompare.clear();
				fitindex=winner;
				lexpool = fitindex.size();
				for (int i=0; i<fitindex.size();++i)
					fitcompare.push_back(pop.at(fitindex[i]).fitlex[fitorder[h]]);
			}
			else pass=false;
		}
		//get lowest fitness, return pop.at(fitindex.at(lowest)).id
		if(winner.size()>1)
			parloc.at(i)=winner.at(r[omp_get_thread_num()].rnd_int(0,winner.size()-1));
		else
			parloc.at(i)=winner[0];

		minfit=p.max_fit;
		fitcompare.clear();
		fitindex.clear();
	}//for (int i=0;i<pop.size();++i)

}