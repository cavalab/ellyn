#include "stdafx.h"
#include "state.h"
#include "pop.h"
#include "params.h"
#include "Fitness.h"
#include "Gen2Phen.h"

void HillClimb(ind& oldind,params& p,vector<Randclass>& r,data& d,state& s)
{
	//#pragma omp parallel for
	//for (int i=0; i<pop.size(); i++) // for each individual
	//	{
			vector<ind> tmp_ind(1,oldind); 
			tmp_ind[0].clrPhen(p.sim_nom_mod); // clear phenotype
			bool updated=0;

			for (int j=0;j<p.pHC_its; j++) // for number of specified iterations
			{
				if (updated)
				{
					tmp_ind[0] = oldind; 
					tmp_ind[0].clrPhen(p.sim_nom_mod); // clear phenotype
				}
				for (int h= 0; h<tmp_ind[0].line.size();h++) // for length of genotype
				{
					if(tmp_ind[0].line.at(h)>=100 && tmp_ind[0].epiline.at(h)) // insert function with true epiline value
					{
						if (is_number(oldind.args.at(tmp_ind[0].line.at(h)-100))) //if constant, perturb with gaussian noise
						{
							float num = stof(oldind.args.at(tmp_ind[0].line.at(h)-100),0);
							num = num + r[omp_get_thread_num()].gasdev()*num/10; // 10% gaussian noise perturbation
							if (p.ERCints)
								tmp_ind[0].args.at(tmp_ind[0].line.at(h)-100) = to_string(static_cast<long long>(num));
							else
								tmp_ind[0].args.at(tmp_ind[0].line.at(h)-100) = to_string(static_cast<long double>(num));
						}
					}
				}
				Gen2Phen(tmp_ind,p);
				Fitness(tmp_ind,p,d,s); //get fitness 
				if ( tmp_ind[0].fitness < oldind.fitness) // if fitness is better, replace individual
				{
					oldind = tmp_ind[0];
					updated = true;
					s.pHC_updates[omp_get_thread_num()]++;
				}
				//else if (tmp_ind[0].fitness == oldind.fitness && tmp_ind[0].eqn.size() < oldind.eqn.size()) // if fitness is same but equation is smaller, replace individual
				//{
				//	oldind = tmp_ind[0];
				//	updated = true;
				//	s.pHC_updates[omp_get_thread_num()]++;
				//}
				else
					updated = false;
			}
		//}
}