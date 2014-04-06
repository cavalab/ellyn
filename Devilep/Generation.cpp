#include "stdafx.h"
#include "pop.h"
#include "params.h"
#include "state.h"
#include "Fitness.h"
#include "Generationfns.h"
#include "InitPop.h"

void Generation(vector<ind>& pop,params& p,vector<Randclass>& r,data& d,state& s)
{
	
	//ind (*parloc)[p.popsize-1] = &pop;

	switch (p.sel)
	{
	case 1: // tournament selection
		{
		//if (p.loud) boost::progress_timer timer;
		// return pop ids for parents 
		vector<unsigned int> parloc(pop.size());
		//if (p.loud ) fcout << "     Tournament...";
		Tournament(pop,parloc,p,r);		
		//if (p.loud ) fcout << "     Apply Genetics...";
		try
		{
			ApplyGenetics(pop,parloc,p,r,d);
		}
		catch(...)
				{
					cout<<"40\n";
					throw;
				}
		//if (p.loud ) fcout << "     Gen 2 Phen...";
		//if (p.loud ) fcout << "     Fitness...";
		Fitness(pop,p,d,s);
		//get mutation/crossover stats
		s.setCrossPct(pop);

		/*if (p.pHC_on && p.ERC)
		{
				for(int i=0; i<pop.size(); i++)
					HillClimb(pop.at(i),p,r,d,s);
		}
		if (p.eHC_on) 
		{
				for(int i=0; i<pop.size(); i++)
					EpiHC(pop.at(i),p,r,d,s);
		}*/
		break;
		}
	case 2: // deterministic crowding
		{	
		DC(pop,p,r,d,s);
		break;
		}
	case 3: // lexicase
		{
		break;
		}
	case 4: //age-fitness pareto front
		{//scheme:
			// produce a new population equal in size to the old
			// pool all individuals from both populations
			
			AgeBreed(pop,p,r,d,s);
			// add one new individual
			vector<ind> tmppop(1);
			tmppop[0].age=0;
			InitPop(tmppop,p,r,d);
			
			Fitness(tmppop,p,d,s);
			pop.push_back(tmppop[0]);
			// select new population with tournament size 2, based on pareto age-fitness
			AgeFitSelect(pop,p,r);
		break;

		}
	default:
		cout << "Bad p.sel parameter. " << endl;
		break;
	}
	//if (p.loud ) fcout << "  Gentime...";
}



