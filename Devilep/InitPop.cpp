#include "stdafx.h"
#include "pop.h"
#include "params.h"

//extern params p;
//extern vector<Randclass> r;
using namespace std;
/*Initialize population
		make genotypes
		genotype to phenotype
		calculate fitness
*/
void makeline(ind&,params& p,vector<Randclass>& r);

void InitPop(vector<ind> &pop,params& p, vector<Randclass>& r)
{
	//boost::progress_timer timer;
	#pragma omp parallel for
	for(int i=0;i<pop.size();i++)
	{
		makeline(pop.at(i),p,r);
		pop.at(i).origin = 'i';
		if (p.eHC_on)
		{
			for (int j=0;j<pop.at(i).line.size();j++)
			{
				float tmp = r[omp_get_thread_num()].rnd_flt(0,1);
				//cout << "tmp: " << tmp << " p.eHC: " << p.eHC_init; 
				if (tmp <= p.eHC_init)
				{
					pop.at(i).epiline.push_back(true);
				//	cout << " epi true \n";
				}
				else
				{
					pop.at(i).epiline.push_back(false);
				//	cout << " epi false \n";
				}
			}
		}
		else
			pop.at(i).epiline.assign(pop.at(i).line.size(),1); // set all values to true 
	}
	//cout <<"\nInit Pop time: ";
}

void makeline(ind& newind,params& p,vector<Randclass>& r)
{
	// construct line 
	// obtain a seed from the system clock:
	
	// random number generator
	//mt19937_64 engine(p.seed);
	
	int linelen = r[omp_get_thread_num()].rnd_int(p.min_len,p.max_len);
	
	vector <string> load_choices = p.allblocks;
	
	
	//uniform_int_distribution<int> dist(0,21);
	if (p.ERC)
	{
		float ERC;
		for (int j=0;j<p.numERC;j++)
		{
			ERC = r[omp_get_thread_num()].rnd_flt((float)p.minERC,(float)p.maxERC);
			if(p.ERCints)
				load_choices.push_back(to_string(static_cast<long long>(ERC)));
			else
				load_choices.push_back(to_string(static_cast<long double>(ERC)));
		}
	}
	vector<float> wheel(p.op_weight.size());
	//wheel.resize(p.op_weight.size());
	if (p.weight_ops_on) //fns are weighted
	{
		partial_sum(p.op_weight.begin(), p.op_weight.end(), wheel.begin());
	}
	else
	{}
		
		

	int choice;
	for (int x = 0; x<linelen; x++)
	{
		if (p.weight_ops_on) //fns are weighted
		{
			float tmp = r[omp_get_thread_num()].rnd_flt(0,1);
			if (tmp < wheel.at(0))
				choice=0;
			else
			{
				for (unsigned int k=1;k<wheel.size();k++)
				{
					if(tmp<wheel.at(k) && tmp>=wheel.at(k-1))
						choice = k;
				}
			}
		}
		else 
			choice = r[omp_get_thread_num()].rnd_int(0,p.op_choice.size()-1);

		
		if (choice==0) // store number reference to argument from all blocks in next position in line
		{
			newind.args.push_back(load_choices.at(r[omp_get_thread_num()].rnd_int(0,load_choices.size()-1)));
			newind.line.push_back(100+newind.args.size()-1);
		}
		else 
			newind.line.push_back(p.op_choice.at(choice));
	}
}

