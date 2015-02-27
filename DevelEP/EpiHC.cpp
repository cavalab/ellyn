#include "stdafx.h"
#include "pop.h"
#include "params.h"
#include "rnd.h"
#include "state.h"
#include "data.h"
#include "FitnessEstimator.h"
#include "Fitness.h"
#include "Gen2Phen.h"
#include "general_fns.h"

void EpiHC(ind& oldind,params& p,vector<Randclass>& r,data& d,state& s,FitnessEstimator& FE)
{
	////#pragma omp parallel for
	//for (int i=0; i<pop.size(); ++i) // for each individual
	//{
		
		vector<ind> tmp_ind;
		//tmp_ind[0] = oldind;
		//makenew(tmp_ind[0]);
		//tmp_ind[0].clrPhen();
		bool updated = false;
		int numchanges = 0;
		bool pass = true;
		int outstart = 0;
		int linestart = 0;
		vector<float> init_stack;
		for (int j=0;j<p.eHC_its; ++j) // for number of specified iterations
		{
		
			/* if (updated)
			{
				//tmp_ind.clear();
				tmp_ind.push_back(oldind); 
				//makenew(tmp_ind[0]);
				tmp_ind[0].clrPhen(); // clear phenotype
			} */
			// THIS NEEDS TO BE TESTED
			tmp_ind.resize(0);
			tmp_ind.push_back(oldind);
			tmp_ind[0].clrPhen();
			
			for(unsigned int h = 0;h<tmp_ind[0].line.size();++h)
			{
				if(r[omp_get_thread_num()].rnd_flt(0,1)<=p.eHC_prob){
					tmp_ind[0].line.at(h).on = !tmp_ind[0].line.at(h).on;
					numchanges++;
					if (numchanges==1){
						linestart = h;
						for (int x = 0; x<h;x++)
							outstart+=(int)(tmp_ind[0].line.at(x).on);
					}
				}
			}	
			if (numchanges==0){
				int tmp = r[omp_get_thread_num()].rnd_int(0,tmp_ind[0].line.size()-1);
				tmp_ind[0].line.at(tmp).on = !tmp_ind[0].line.at(tmp).on;
				linestart=tmp;
				for (int x = 0; x<tmp;x++)
					outstart+=(int)(tmp_ind[0].line.at(x).on);
			}

			//Gen2Phen(tmp_ind,p);
			//get fitness 
			if (p.eHC_slim){
				tmp_ind[0].outstack.resize(outstart*d.vals.size());
				tmp_ind[0].outstacklen.resize(outstart);
				pass = SlimFitness(tmp_ind[0],p,d,s,FE,linestart,oldind.fitness);
			}
			else
				Fitness(tmp_ind,p,d,s,FE); 
			if(pass){
				if ( tmp_ind[0].fitness < oldind.fitness) // if fitness is better, replace individual
				{
					oldind = tmp_ind[0];
					updated = true;
					//tmp_ind.resize(0);
					s.eHC_updates[omp_get_thread_num()]++;
				}
				else if (tmp_ind[0].fitness == oldind.fitness && tmp_ind[0].complexity < oldind.complexity) // if fitness is same but equation is smaller, replace individual
				{
					oldind = tmp_ind[0];
					updated = true;
					//tmp_ind.resize(0);
					s.eHC_updates[omp_get_thread_num()]++;
				}
				else
					updated = false;
			}
			else 
				pass = 1;

			numchanges=0;
			outstart=0;
		}
		//tmp_ind.clear();
	//}

}