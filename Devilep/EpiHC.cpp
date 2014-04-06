#include "stdafx.h"
#include "pop.h"
#include "params.h"
#include "state.h"
#include "Fitness.h"
#include "Gen2Phen.h"

void EpiHC(ind& oldind,params& p,vector<Randclass>& r,data& d,state& s)
{
	////#pragma omp parallel for
	//for (int i=0; i<pop.size(); i++) // for each individual
	//{
		
		vector<ind> tmp_ind(1);
		tmp_ind[0]=oldind; 
		tmp_ind[0].clrPhen(p.sim_nom_mod);
		bool updated = false;
		for (int j=0;j<p.eHC_its; j++) // for number of specified iterations
		{
			if (updated)
			{
				tmp_ind[0] = oldind; 
				tmp_ind[0].clrPhen(p.sim_nom_mod); // clear phenotype
			}

			for(unsigned int h = 0;h<tmp_ind[0].line.size();h++)
			{
				if(r[omp_get_thread_num()].rnd_flt(0,1)<=p.eHC_prob)
					tmp_ind[0].line.at(h)->on = !tmp_ind[0].line.at(h)->on;
			}	
			//Gen2Phen(tmp_ind,p);
			Fitness(tmp_ind,p,d,s); //get fitness 
			if ( tmp_ind[0].fitness < oldind.fitness) // if fitness is better, replace individual
			{
				oldind = tmp_ind[0];
				updated = true;
				s.eHC_updates[omp_get_thread_num()]++;
			}
			else if (tmp_ind[0].fitness == oldind.fitness && tmp_ind[0].eqn.size() < oldind.eqn.size()) // if fitness is same but equation is smaller, replace individual
			{
				oldind = tmp_ind[0];
				updated = true;
				s.eHC_updates[omp_get_thread_num()]++;
			}
			else
				updated = false;
		}
	//}

}