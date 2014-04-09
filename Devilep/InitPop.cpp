#include "stdafx.h"
#include "pop.h"
#include "params.h"
//#include "RPN_class.h"
//extern params p;
//extern vector<Randclass> r;
using namespace std;
/*Initialize population
		make genotypes
		genotype to phenotype
		calculate fitness
*/
void makeline(ind&,params& p,vector<Randclass>& r,data& d);

void InitPop(vector<ind> &pop,params& p, vector<Randclass>& r,data& d)
{
	//boost::progress_timer timer;
	//#pragma omp parallel for

	for(int i=0;i<pop.size();i++)
	{
		makeline(pop.at(i),p,r,d);
		//remove dangling numbers and variables from end
		/*while((pop.at(i).line.back()->type=='n' || pop.at(i).line.back()->type=='v') && pop.at(i).line.size()>1)
			pop.at(i).line.pop_back();*/
			
		pop.at(i).origin = 'i';
		if (p.eHC_on)
		{
			for (int j=0;j<pop.at(i).line.size();j++)
			{
				float tmp = r[omp_get_thread_num()].rnd_flt(0,1);
				//cout << "tmp: " << tmp << " p.eHC: " << p.eHC_init; 
				if (tmp > p.eHC_init)
				{
					pop.at(i).line[j]->on=!pop.at(i).line[j]->on;
				}
			}
		}


	}
	//cout <<"\nInit Pop time: ";
}

void makeline(ind& newind,params& p,vector<Randclass>& r,data& d)
{
	// construct line 
	// obtain a seed from the system clock:
	
	// random number generator
	//mt19937_64 engine(p.seed);
	
	int linelen = r[omp_get_thread_num()].rnd_int(p.min_len,p.max_len);
	
	int choice=0;
	
	//uniform_int_distribution<int> dist(0,21);
	vector<float> wheel(p.op_weight.size());
	//wheel.resize(p.op_weight.size());
	if (p.weight_ops_on) //fns are weighted
	{
		partial_sum(p.op_weight.begin(), p.op_weight.end(), wheel.begin());
	}

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
					if(tmp<wheel.at(k) && tmp>=wheel.at(k-1)){
						choice = k;
						break;
					}
				}
			}
		}
		else 
			choice = r[omp_get_thread_num()].rnd_int(0,p.op_choice.size()-1);

		string varchoice; 

		switch (p.op_choice.at(choice))
			{
			case 0: //load number
				if(p.ERCints)
					newind.line.push_back(shared_ptr<node>(new n_num((float)r[omp_get_thread_num()].rnd_int(p.minERC,p.maxERC))));
				else
					newind.line.push_back(shared_ptr<node>(new n_num(r[omp_get_thread_num()].rnd_flt(p.minERC,p.maxERC))));
				break;
			case 1: //load variable
				varchoice = d.label.at(r[omp_get_thread_num()].rnd_int(0,d.label.size()-1));
				newind.line.push_back(shared_ptr<node>(new n_sym(varchoice)));
				break;
			case 2: // +
				newind.line.push_back(shared_ptr<node>(new n_add()));
				break;
			case 3: // -
				newind.line.push_back(shared_ptr<node>(new n_sub()));
				break;
			case 4: // *
				newind.line.push_back(shared_ptr<node>(new n_mul()));
				break;
			case 5: // /
				newind.line.push_back(shared_ptr<node>(new n_div()));
				break;
			case 6: // sin
				newind.line.push_back(shared_ptr<node>(new n_sin()));
				break;
			case 7: // cos
				newind.line.push_back(shared_ptr<node>(new n_cos()));
				break;
			case 8: // exp
				newind.line.push_back(shared_ptr<node>(new n_exp()));
				break;
			case 9: // log
				newind.line.push_back(shared_ptr<node>(new n_log()));
				break;
			}
	}
}

//void makestack(ind& newind,params& p,vector<Randclass>& r)
//{
//	// construct line 
//	// obtain a seed from the system clock:
//	
//	// random number generator
//	//mt19937_64 engine(p.seed);
//	
//	int linelen = r[omp_get_thread_num()].rnd_int(p.min_len,p.max_len);
//	
//	vector <string> load_choices(p.allblocks);
//
//	int choice=0;
//	
//	//uniform_int_distribution<int> dist(0,21);
//	if (p.ERC)
//	{
//		float ERC;
//		for (int j=0;j<p.numERC;j++)
//		{
//			ERC = r[omp_get_thread_num()].rnd_flt((float)p.minERC,(float)p.maxERC);
//			if(p.ERCints)
//				load_choices.push_back(to_string(static_cast<long long>(ERC)));
//			else
//				load_choices.push_back(to_string(static_cast<long double>(ERC)));
//		}
//	}
//	vector<float> wheel(p.op_weight.size());
//	//wheel.resize(p.op_weight.size());
//	if (p.weight_ops_on) //fns are weighted
//	{
//		partial_sum(p.op_weight.begin(), p.op_weight.end(), wheel.begin());
//	}
//
//	for (int x = 0; x<linelen; x++)
//	{
//		if (p.weight_ops_on) //fns are weighted
//		{
//			float tmp = r[omp_get_thread_num()].rnd_flt(0,1);
//			if (tmp < wheel.at(0))
//				choice=0;
//			else
//			{
//				for (unsigned int k=1;k<wheel.size();k++)
//				{
//					if(tmp<wheel.at(k) && tmp>=wheel.at(k-1))
//						choice = k;
//				}
//			}
//		}
//		else 
//			choice = r[omp_get_thread_num()].rnd_int(0,p.op_choice.size()-1);
//
//		
//		if (choice==0) // number reference to argument from all blocks in next position in line
//		{
//			if(p.ERCints)
//				newind.line.push_back(ops(p.op_choice.at(choice),r[omp_get_thread_num()].rnd_int((float)p.minERC,(float)p.maxERC)));
//			else
//				newind.line.push_back(ops(p.op_choice.at(choice),r[omp_get_thread_num()].rnd_flt((float)p.minERC,(float)p.maxERC)));
//		}
//		else if (choice==1) // pick pointer values from mapped variables in data struct
//		{
//		}
//		else 
//			newind.line.push_back(ops(p.op_choice.at(choice)));
//	}
//}