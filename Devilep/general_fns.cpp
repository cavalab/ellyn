#include "stdafx.h"
#include "pop.h"
#include "params.h"

bool is_number(const std::string& s)
{
    std::string::const_iterator it = s.begin();
    while (it != s.end() && (std::isdigit(*it) || (*it=='-') || (*it=='.'))) ++it;
    return !s.empty() && it == s.end();
}


void NewInstruction(ind& newind,int loc,params& p,vector<Randclass>& r)
{
	
	vector <string> load_choices = p.allblocks;
	//cout<<"iterator definition\n";
	std::vector<int>::iterator it;
	it = newind.line.begin();
	//cout<<"end def \n";

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
		choice = r[omp_get_thread_num()].rnd_int(0,p.op_list.size()-1);

		
	if (choice==0) // store number reference to argument from all blocks in next position in line
	{
		
		newind.args.push_back(load_choices.at(r[omp_get_thread_num()].rnd_int(0,load_choices.size()-1)));
		/*if(it+loc>newind.line.size())
		{
			cout<< "it+ loc: " << *it+loc << "line size: " << newind.line.size() << ". too big" <<endl;
			newind.line.push_back(100+newind.args.size()-1);
		}
		else*/
		//cout<< "newind insert\n";
			/*newind.line.insert(it+loc,100+newind.args.size()-1);*/
		newind.line.at(loc) = 100+newind.args.size()-1;
		//cout<<"end newind insert\n";
	}
	else 
	{
		/*if(*it+loc>newind.line.size())
		{
			cout<< "it+ loc: " << it+loc << "line size: " << newind.line.size() <<". too big!\n";
			
			newind.line.push_back(p.op_choice.at(choice));
		}
		else*/
		//cout<< "newind insert\n";
			//newind.line.insert(it+loc,p.op_choice.at(choice));
		newind.line.at(loc) = p.op_choice.at(choice);
			//cout<<"end newind insert\n";
	}

}