#include "stdafx.h"
#include "pop.h"
#include "general_fns.h"

void Mutate(ind& par,vector<ind>& tmppop,params& p,vector<Randclass>& r,data& d)
{
	vector<unsigned int> ichange;
	vector<ind> kid(1,par);
	makenewcopy(kid[0]);
	kid[0].origin='m';
	kid[0].parentfitness=par.fitness;
	kid[0].clrPhen();
		

	for(unsigned int i = 0;i<kid[0].line.size();i++)
	{
		if(r[omp_get_thread_num()].rnd_flt(0,1)<=p.mut_ar)
			ichange.push_back(i);
	}
	for(unsigned int j=0;j<ichange.size();j++)
	{
		if(kid[0].line.at(ichange.at(j))->type=='n') // insert function
		{
			if (p.ERC)
			{
				    float num = static_pointer_cast<n_num>(kid[0].line.at(ichange.at(j)))->value;
					num = num/2 + r[omp_get_thread_num()].gasdev()*num/2;
					static_pointer_cast<n_num>(kid[0].line.at(ichange.at(j)))->value = num;
			}
		}
		else if (kid[0].line.at(ichange.at(j))->type=='v'){
				//else if variable, pick random variable replacement
			static_pointer_cast<n_sym>(kid[0].line.at(ichange.at(j)))->varname = d.label.at(r[omp_get_thread_num()].rnd_int(0,d.label.size()-1));
		}
		else
		{
			
			NewInstruction(kid[0],ichange.at(j),p,r,d);
		}
	}
	
	tmppop.push_back(kid[0]);
	kid.clear();
		
}
