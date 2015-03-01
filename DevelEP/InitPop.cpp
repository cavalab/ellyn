#include "stdafx.h"
#include "pop.h"
#include "params.h"
#include "rnd.h"
#include "data.h"
#include "general_fns.h"
//#include "RPN_class.h"
//extern params p;
//extern vector<Randclass> r;
using namespace std;
/*Initialize population
		make genotypes
		genotype to phenotype
		calculate fitness
*/
void makeline(ind&,params& p,vector<Randclass>& r);
void makeline_rec(ind&,params& p,vector<Randclass>& r,int linelen);
int maketree(vector <node>& line, int level, bool exactlevel, int lastnode,params& p,vector<Randclass>& r);
float round(float d)
{
  return floor(d + 0.5);
}
void InitPop(vector<ind> &pop,params& p, vector<Randclass>& r)
{
	//boost::progress_timer timer;
	//#pragma omp parallel for

	for(int i=0;i<pop.size();++i)
	{
		if (!p.init_trees){
			makeline(pop.at(i),p,r);
			if (p.eHC_on)
			{
				for (int j=0;j<pop.at(i).line.size();++j)
				{
					float tmp = r[omp_get_thread_num()].rnd_flt(0,1);
					//cout << "tmp: " << tmp << " p.eHC: " << p.eHC_init; 
					if (tmp > p.eHC_init)
					{
						pop.at(i).line[j].on=false;
					}
				}
			}
		}
		else{
			int linelen = r[omp_get_thread_num()].rnd_int(p.min_len,p.max_len);
			if (p.eHC_on){
				int onlen = int(linelen*p.eHC_init);
				makeline_rec(pop.at(i),p,r,onlen);
				int offlen = linelen-onlen;
				for (int j=0;j<offlen;++j){
					int loc = r[omp_get_thread_num()].rnd_int(0,pop.at(i).line.size()-1);
					 InsInstruction(pop.at(i),loc,p,r);
					pop.at(i).line[loc].on=false;
				}
			}
			else
				makeline_rec(pop.at(i),p,r,linelen);
		}
		//remove dangling numbers and variables from end
		/*while((pop.at(i).line.back()->type=='n' || pop.at(i).line.back()->type=='v') && pop.at(i).line.size()>1)
			pop.at(i).line.pop_back();*/
			
		pop.at(i).origin = 'i';
		
		assert(!pop.at(i).line.empty());


	}
	//cout <<"\nInit Pop time: ";
}

void makeline(ind& newind,params& p,vector<Randclass>& r)
{
	// construct line 
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
				for (unsigned int k=1;k<wheel.size();++k)
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
		int seedchoice;
		vector<node> tmpstack;

		switch (p.op_choice.at(choice))
			{
			case 0: //load number
				if(p.ERC){ // if ephemeral random constants are on
					if (!p.cvals.empty()){ // if there are constants defined by the user
						if (r[omp_get_thread_num()].rnd_flt(0,1)<.5)
							//newind.line.push_back(shared_ptr<node>(new n_num(p.cvals.at(r[omp_get_thread_num()].rnd_int(0,p.cvals.size()-1)))));
							newind.line.push_back(node(p.cvals.at(r[omp_get_thread_num()].rnd_int(0,p.cvals.size()-1))));
						else{	
							if(p.ERCints)
								//newind.line.push_back(shared_ptr<node>(new n_num((float)r[omp_get_thread_num()].rnd_int(p.minERC,p.maxERC))));
								newind.line.push_back(node((float)r[omp_get_thread_num()].rnd_int(p.minERC,p.maxERC)));
							else
								//newind.line.push_back(shared_ptr<node>(new n_num(r[omp_get_thread_num()].rnd_flt(p.minERC,p.maxERC))));
								newind.line.push_back(node(r[omp_get_thread_num()].rnd_flt(p.minERC,p.maxERC)));
						}
					}
					else{
						if(p.ERCints)
								//newind.line.push_back(shared_ptr<node>(new n_num((float)r[omp_get_thread_num()].rnd_int(p.minERC,p.maxERC))));
						newind.line.push_back(node((float)r[omp_get_thread_num()].rnd_int(p.minERC,p.maxERC)));
							else
								//newind.line.push_back(shared_ptr<node>(new n_num(r[omp_get_thread_num()].rnd_flt(p.minERC,p.maxERC))));
						newind.line.push_back(node(r[omp_get_thread_num()].rnd_flt(p.minERC,p.maxERC)));
					}
				}
				else if (!p.cvals.empty()) // if ERCs are off, but there are constants defined by the user
				{
					//newind.line.push_back(shared_ptr<node>(new n_num(p.cvals.at(r[omp_get_thread_num()].rnd_int(0,p.cvals.size()-1)))));
					newind.line.push_back(node(p.cvals.at(r[omp_get_thread_num()].rnd_int(0,p.cvals.size()-1))));
				}
				break;
			case 1: //load variable
				varchoice = p.allvars.at(r[omp_get_thread_num()].rnd_int(0,p.allvars.size()-1));
				//varchoice = d.label.at(r[omp_get_thread_num()].rnd_int(0,d.label.size()-1));
				newind.line.push_back(node(varchoice));
				break;
			case 2: // +
				//newind.line.push_back(shared_ptr<node>(new n_add()));
				newind.line.push_back(node('+'));
				break;
			case 3: // -
				//newind.line.push_back(shared_ptr<node>(new n_sub()));
				newind.line.push_back(node('-'));
				break;
			case 4: // *
				//newind.line.push_back(shared_ptr<node>(new n_mul()));
				newind.line.push_back(node('*'));
				break;
			case 5: // /
				//newind.line.push_back(shared_ptr<node>(new n_div()));
				newind.line.push_back(node('/'));
				break;
			case 6: // sin
				//newind.line.push_back(shared_ptr<node>(new n_sin()));
				newind.line.push_back(node('s'));
				break;
			case 7: // cos
				//newind.line.push_back(shared_ptr<node>(new n_cos()));
				newind.line.push_back(node('c'));
				break;
			case 8: // exp
				//newind.line.push_back(shared_ptr<node>(new n_exp()));
				newind.line.push_back(node('e'));
				break;
			case 9: // log
				//newind.line.push_back(shared_ptr<node>(new n_log()));
				newind.line.push_back(node('l'));
				break;
			case 10: // seed
				seedchoice = r[omp_get_thread_num()].rnd_int(0,p.seedstacks.size()-1);
				//copystack(p.seedstacks.at(seedchoice),tmpstack);
				tmpstack = p.seedstacks.at(seedchoice);
				
				for(int i=0;i<tmpstack.size(); ++i)
				{
					if (x<p.max_len){
						newind.line.push_back(tmpstack[i]);
						x++;
					}
				}
				tmpstack.clear();
				break;
			}
	}
	/*while (newind.line.back()->type=='n'||newind.line.back()->type=='v')
		newind.line.pop_back();*/
}
void makeline_rec(ind& newind,params& p,vector<Randclass>& r, int linelen)
{
	// recursive version of makeline that creates lines that are complete with respect to stack operations
	// in other words the entire line is guaranteed to be a valid syntax tree
	// construct line 
	
	
	
	int choice=0;
	
	int tmp = maketree(newind.line,linelen,1,0,p,r);

}
void getChoice(int& choice, int min_arity, int max_arity, params& p,vector<Randclass>& r)
{
	vector<int> choices;
	vector<float> op_weight;
	int tmpchoice;
	for (int i=0;i<p.op_arity.size();++i)
	{
		if (p.op_arity[i] >= min_arity && p.op_arity[i] <= max_arity)
		{
			choices.push_back(i);
			if(p.weight_ops_on)
				op_weight.push_back(p.op_weight[i]);
		}
	}
	if(!choices.empty()){
		vector<float> wheel(choices.size());
		//wheel.resize(p.op_weight.size());
		if (p.weight_ops_on) //fns are weighted
			partial_sum(op_weight.begin(), op_weight.end(), wheel.begin());

		if (p.weight_ops_on) //fns are weighted
		{
			float tmp = r[omp_get_thread_num()].rnd_flt(0,*std::max_element(wheel.begin(),wheel.end()));
			if (tmp < wheel.at(0))
				tmpchoice=0;
			else
			{
				for (unsigned int k=1;k<wheel.size();++k)
				{
					if(tmp<wheel.at(k) && tmp>=wheel.at(k-1)){
						tmpchoice = k;
						break;
					}
				}
			}
		}
		else 
			tmpchoice =r[omp_get_thread_num()].rnd_int(0,choices.size()-1);

		choice = choices[tmpchoice];
	}
	else
		choice = -1;
}

void push_back_node(vector <node>& line, int choice, params& p,vector<Randclass>& r)
{
	string varchoice; 

	switch (p.op_choice.at(choice))
	{
		case 0: //load number
			if(p.ERC){ // if ephemeral random constants are on
				if (!p.cvals.empty()){ // if there are constants defined by the user
					if (r[omp_get_thread_num()].rnd_flt(0,1)<.5)
						//line.push_back(shared_ptr<node>(new n_num(p.cvals.at(r[omp_get_thread_num()].rnd_int(0,p.cvals.size()-1)))));
						line.push_back(node(p.cvals.at(r[omp_get_thread_num()].rnd_int(0,p.cvals.size()-1))));
					else{	
						if(p.ERCints)
							/*line.push_back(shared_ptr<node>(new n_num((float)r[omp_get_thread_num()].rnd_int(p.minERC,p.maxERC))));*/
							line.push_back(node((float)r[omp_get_thread_num()].rnd_int(p.minERC,p.maxERC)));
						else
							//line.push_back(shared_ptr<node>(new n_num(r[omp_get_thread_num()].rnd_flt(p.minERC,p.maxERC))));
							line.push_back(node(r[omp_get_thread_num()].rnd_flt(p.minERC,p.maxERC)));
					}
				}
				else{
					if(p.ERCints)
							//line.push_back(shared_ptr<node>(new n_num((float)r[omp_get_thread_num()].rnd_int(p.minERC,p.maxERC))));
							line.push_back(node((float)r[omp_get_thread_num()].rnd_int(p.minERC,p.maxERC)));
						else
							//line.push_back(shared_ptr<node>(new n_num(r[omp_get_thread_num()].rnd_flt(p.minERC,p.maxERC))));
							line.push_back(node(r[omp_get_thread_num()].rnd_flt(p.minERC,p.maxERC)));
				}
			}
			else if (!p.cvals.empty()) // if ERCs are off, but there are constants defined by the user
			{
				//line.push_back(shared_ptr<node>(new n_num(p.cvals.at(r[omp_get_thread_num()].rnd_int(0,p.cvals.size()-1)))));
				line.push_back(node(p.cvals.at(r[omp_get_thread_num()].rnd_int(0,p.cvals.size()-1))));
			}
			break;
		case 1: //load variable
			varchoice = p.allvars.at(r[omp_get_thread_num()].rnd_int(0,p.allvars.size()-1));
			//varchoice = d.label.at(r[omp_get_thread_num()].rnd_int(0,d.label.size()-1));
			/*line.push_back(shared_ptr<node>(new n_sym(varchoice)));*/
			line.push_back(node(varchoice));
			break;
		case 2: // +
			//line.push_back(shared_ptr<node>(new n_add()));
			line.push_back(node('+'));
			break;
		case 3: // -
			//line.push_back(shared_ptr<node>(new n_sub()));
			line.push_back(node('-'));
			break;
		case 4: // *
			//line.push_back(shared_ptr<node>(new n_mul()));
			line.push_back(node('*'));
			break;
		case 5: // /
			//line.push_back(shared_ptr<node>(new n_div()));
			line.push_back(node('/'));
			break;
		case 6: // sin
			//line.push_back(shared_ptr<node>(new n_sin()));
			line.push_back(node('s'));
			break;
		case 7: // cos
			//line.push_back(shared_ptr<node>(new n_cos()));
			line.push_back(node('c'));
			break;
		case 8: // exp
			//line.push_back(shared_ptr<node>(new n_exp()));
			line.push_back(node('e'));
			break;
		case 9: // log
			//line.push_back(shared_ptr<node>(new n_log()));
			line.push_back(node('l'));
			break;
	}
}
void push_front_node(vector <node>& line, int choice, params& p,vector<Randclass>& r)
{
	string varchoice; 

	switch (p.op_choice.at(choice))
	{
		case 0: //load number
			if(p.ERC){ // if ephemeral random constants are on
				if (!p.cvals.empty()){ // if there are constants defined by the user
					if (r[omp_get_thread_num()].rnd_flt(0,1)<.5)
						//line.insert(line.begin(),shared_ptr<node>(new n_num(p.cvals.at(r[omp_get_thread_num()].rnd_int(0,p.cvals.size()-1)))));
						line.insert(line.begin(),node(p.cvals.at(r[omp_get_thread_num()].rnd_int(0,p.cvals.size()-1))));
					else{	
						if(p.ERCints)
							/*line.insert(line.begin(),shared_ptr<node>(new n_num((float)r[omp_get_thread_num()].rnd_int(p.minERC,p.maxERC))));*/
							line.insert(line.begin(),node((float)r[omp_get_thread_num()].rnd_int(p.minERC,p.maxERC)));
						else
							//line.insert(line.begin(),shared_ptr<node>(new n_num(r[omp_get_thread_num()].rnd_flt(p.minERC,p.maxERC))));
							line.insert(line.begin(),node(r[omp_get_thread_num()].rnd_flt(p.minERC,p.maxERC)));
					}
				}
				else{
					if(p.ERCints)
							//line.insert(line.begin(),shared_ptr<node>(new n_num((float)r[omp_get_thread_num()].rnd_int(p.minERC,p.maxERC))));
							line.insert(line.begin(),node((float)r[omp_get_thread_num()].rnd_int(p.minERC,p.maxERC)));
						else
							//line.insert(line.begin(),shared_ptr<node>(new n_num(r[omp_get_thread_num()].rnd_flt(p.minERC,p.maxERC))));
							line.insert(line.begin(),node(r[omp_get_thread_num()].rnd_flt(p.minERC,p.maxERC)));
				}
			}
			else if (!p.cvals.empty()) // if ERCs are off, but there are constants defined by the user
			{
				//line.insert(line.begin(),shared_ptr<node>(new n_num(p.cvals.at(r[omp_get_thread_num()].rnd_int(0,p.cvals.size()-1)))));
				line.insert(line.begin(),node(p.cvals.at(r[omp_get_thread_num()].rnd_int(0,p.cvals.size()-1))));
			}
			break;
		case 1: //load variable
			varchoice = p.allvars.at(r[omp_get_thread_num()].rnd_int(0,p.allvars.size()-1));
			//varchoice = d.label.at(r[omp_get_thread_num()].rnd_int(0,d.label.size()-1));
			/*line.insert(line.begin(),shared_ptr<node>(new n_sym(varchoice)));*/
			line.insert(line.begin(),node(varchoice));
			break;
		case 2: // +
			//line.insert(line.begin(),shared_ptr<node>(new n_add()));
			line.insert(line.begin(),node('+'));
			break;
		case 3: // -
			//line.insert(line.begin(),shared_ptr<node>(new n_sub()));
			line.insert(line.begin(),node('-'));
			break;
		case 4: // *
			//line.insert(line.begin(),shared_ptr<node>(new n_mul()));
			line.insert(line.begin(),node('*'));
			break;
		case 5: // /
			//line.insert(line.begin(),shared_ptr<node>(new n_div()));
			line.insert(line.begin(),node('/'));
			break;
		case 6: // sin
			//line.insert(line.begin(),shared_ptr<node>(new n_sin()));
			line.insert(line.begin(),node('s'));
			break;
		case 7: // cos
			//line.insert(line.begin(),shared_ptr<node>(new n_cos()));
			line.insert(line.begin(),node('c'));
			break;
		case 8: // exp
			//line.insert(line.begin(),shared_ptr<node>(new n_exp()));
			line.insert(line.begin(),node('e'));
			break;
		case 9: // log
			//line.insert(line.begin(),shared_ptr<node>(new n_log()));
			line.insert(line.begin(),node('l'));
			break;
	}
}
int maketree(vector<node>& line, int level, bool exactlevel, int lastnode,params& p,vector<Randclass>& r)
{
	int choice; 
//	int splitnodes; 
	int thisnode = lastnode+1;
	int startsize = line.size();
	if (level==1) // choose a terminal because of the level limitation
		getChoice(choice,0,0,p,r); 	
	else if (exactlevel){
		getChoice(choice,1,level-1,p,r);
		if (choice==-1)
			getChoice(choice,0,0,p,r);
	}
	else
		getChoice(choice,0,level-1,p,r);

	// insert choice into line
	push_front_node(line,choice,p,r);
	int a = p.op_arity[choice];
	int newlevel;
	int nodes=0;
	if (a !=0){
		level = level-1; // change splitnodes so that all trees aren't symmetric
		//splitnodes = int(round(float(level)/float(a)));
		//splitnodes = r[omp_get_thread_num()].rnd_int(1,level-1);
	}
	for (int i=1;i<=a;++i)
	{
		if (i==a)
			//newlevel = level-(splitnodes*(a-1));
			//newlevel = r[omp_get_thread_num()].rnd_int(1,level-nodes);
			newlevel = level-nodes;
		else
			//newlevel = splitnodes;
			newlevel = r[omp_get_thread_num()].rnd_int(1,level-1);
		nodes += maketree(line,newlevel,exactlevel,thisnode,p,r);
	}
	return line.size()-startsize;
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
//		for (int j=0;j<p.numERC;++j)
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
//				for (unsigned int k=1;k<wheel.size();++k)
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
