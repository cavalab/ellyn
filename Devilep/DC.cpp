#include "stdafx.h"
#include "pop.h"
#include "params.h"
#include "state.h"
#include "Generationfns.h"
#include "strdist.h"
#include "Gen2Phen.h"
#include "Fitness.h"

// steady-state deterministic crowding algorithm.
// update_locs are the indices of the update locations for passage to hill climbing. 
void DC(vector<ind>& pop,params& p,vector<Randclass>& r,data& d,state& s)
{
	
	vector <ind> tmppop;
	float choice = r[omp_get_thread_num()].rnd_flt(0,1);
	int p1=0;
	int p2=0;

	if (choice<p.rep_wheel[1]) // crossover, as long as there are at least two parents left
	{			
		p1 = r[omp_get_thread_num()].rnd_int(0,pop.size()-1);
		p2 = r[omp_get_thread_num()].rnd_int(0,pop.size()-1);
		// only breed if parents have different fitness 
		if(pop.at(p1).fitness!=pop.at(p2).fitness) 
		{				
			//cross parents to make two kids
			Crossover(pop.at(p1),pop.at(p2),tmppop,p,r);
		}
		else
		{
			Mutate(pop.at(p1),tmppop,p,r);
			Mutate(pop.at(p2),tmppop,p,r);
		}
	}
	else if (choice < p.rep_wheel[2]) //mutation
	{
		p1 = r[omp_get_thread_num()].rnd_int(0,pop.size()-1);
		Mutate(pop.at(p1),tmppop,p,r);
	}
	//if (p.loud ) cout << "  Gen2Phen...\n";
	Gen2Phen(tmppop,p);
	//if (p.loud ) cout << "  Fitness...\n";
	Fitness(tmppop,p,d,s);

	if(tmppop.size()==1)
	{
		if ( tmppop.at(0).fitness <= pop.at(p1).fitness )
		{
			if( tmppop.at(0).fitness < pop.at(p1).fitness )
				s.good_cross[omp_get_thread_num()]=s.good_cross[omp_get_thread_num()]+1;
			else
				s.neut_cross[omp_get_thread_num()]=s.neut_cross[omp_get_thread_num()]+1;

			swap(tmppop.at(0),pop.at(p1));

		}
		else
			s.bad_cross[omp_get_thread_num()]=s.bad_cross[omp_get_thread_num()]+1;

	}
	else if (tmppop.size()==2)
	{
		if (strdist(tmppop.at(0).eqn_form,pop.at(p1).eqn_form) + strdist(tmppop.at(1).eqn_form,pop.at(p2).eqn_form) > strdist(tmppop.at(1).eqn_form,pop.at(p1).eqn_form) + strdist(tmppop.at(0).eqn_form,pop.at(p2).eqn_form))
		{	
			swap(tmppop.at(0),tmppop.at(1));
			tmppop.at(0).parentfitness = pop.at(p1).fitness;
			tmppop.at(1).parentfitness = pop.at(p2).fitness;
		}

		if ( tmppop.at(0).fitness <= pop.at(p1).fitness )
		{
			if( tmppop.at(0).fitness < pop.at(p1).fitness )
				s.good_cross[omp_get_thread_num()]=s.good_cross[omp_get_thread_num()]+1;
			else
				s.neut_cross[omp_get_thread_num()]=s.neut_cross[omp_get_thread_num()]+1;

			swap(tmppop.at(0),pop.at(p1));
		}
		else
			s.bad_cross[omp_get_thread_num()]=s.bad_cross[omp_get_thread_num()]+1;

		if ( tmppop.at(1).fitness <= pop.at(p2).fitness )
		{
			if( tmppop.at(1).fitness < pop.at(p2).fitness )
				s.good_cross[omp_get_thread_num()]=s.good_cross[omp_get_thread_num()]+1;
			else
				s.neut_cross[omp_get_thread_num()]=s.neut_cross[omp_get_thread_num()]+1;

			swap(tmppop.at(1),pop.at(p2));
		}
		else
			s.bad_cross[omp_get_thread_num()]=s.bad_cross[omp_get_thread_num()]+1;
	}
	else
	{
		//cout << "Problem in DC. tmppop size: " << tmppop.size() << endl;
	}
		
}