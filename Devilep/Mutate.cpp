#include "stdafx.h"
#include "pop.h"
#include "general_fns.h"

void Mutate(ind& par,vector<ind>& tmppop,params& p,vector<Randclass>& r,data& d)
{
	vector<unsigned int> ichange;
	//ind kid;

	for(unsigned int i = 0;i<par.line.size();i++)
	{
		if(r[omp_get_thread_num()].rnd_flt(0,1)<=p.mut_ar)
			ichange.push_back(i);
	}
	for(unsigned int j=0;j<ichange.size();j++)
	{
		if(par.line.at(ichange.at(j))->type=='n') // insert function
		{
			if (p.ERC)
			{
				    float num = static_pointer_cast<n_num>(par.line.at(ichange.at(j)))->value;
					num = num/2 + r[omp_get_thread_num()].gasdev()*num/2;
					static_pointer_cast<n_num>(par.line.at(ichange.at(j)))->value = num;
			}
		}
		else if (par.line.at(ichange.at(j))->type=='v'){
				//else if variable, pick random variable replacement
			static_pointer_cast<n_sym>(par.line.at(ichange.at(j)))->varname = d.label.at(r[omp_get_thread_num()].rnd_int(0,d.label.size()-1));
			static_pointer_cast<n_sym>(par.line.at(ichange.at(j)))->valpt = d.datatable.at(static_pointer_cast<n_sym>(par.line.at(ichange.at(j)))->varname);
		}
		else
		{
			
			NewInstruction(par,ichange.at(j),p,r,d);
		}
	}
		par.origin='m';
		par.parentfitness=par.fitness;
		tmppop.push_back(par);
}
