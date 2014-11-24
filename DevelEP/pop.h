// header file for ind struct
#pragma once
#ifndef POP_H
#define POP_H

//#include "params.h"
//#include "data.h"

//#include "RPN_class.h"
#include "op_node.h"
//#include "general_fns.h"
//#include "pareto.h"

struct ind {
	//unsigned int id;
	/*vector <std::shared_ptr<node> > line;*/
	vector <node> line;
	//std::vector <std::string> args;
	std::string eqn;
	std::string eqn_form; // equation form for string distance comparison to other forms
	int complexity;
	//std::vector <int> ptr;
	char origin; // x: crossover, m: mutation, i: initialization

	std::vector<float> output;
	std::vector<float> output_v;
	std::vector<vector<float>> outstack; //optional node-by-node stack tracing

	float abserror;
	float abserror_v;

	float corr;
	float corr_v;
	float VAF;
	float VAF_v;

	float fitness;
	float fitness_v; 
	vector<float> fitlex; // fitnesses for lexicase selection
	float FEvar; //variance in fitness estimates (for sorting purposes)

	float parentfitness;
	int eff_size;
	int age;
	int rank;
	//typedef exprtk::expression<float> expression_t;
	//expression_t expression;
	/*ind()
	{
		eqn=nominal_model;
	}*/
	ind()
	{
		abserror = 0;
		corr = 0;
		age=1;
	}
	~ind() {
		//cout << "ind destructor\n";
		//if(!line.empty())
		//{
		//	for(vector<node*>::iterator it= line.begin();it!=line.end();it++)
		//		delete(*it);
		//	line.clear();
		//	//cout << "ind destructor deleted line nodes\n";
		//}
		//
	}
	//void init(string& nom_mod)
	//{
	//	eqn = nom_mod;
	//	ptr.push_back(1);
	//	ptr.push_back(nom_mod.size()-2);
	//	//nominal_model=nom_mod;
	//	//expression.register_symbol_table(d.symbol_table);
	//}

	void clrPhen()
	{
		abserror = 0;
		abserror_v=0;
		corr = 0;
		corr_v=0;
		fitness=0;
		fitness_v=0;
		eqn = "";
		eqn_form="";
		output.clear();
		output_v.clear();
		outstack.clear();
		// nominal model must be encased in set of parenthesis. the pointer points to that which is encased.
		//ptr[0]= 1;
		//ptr[1] = nom_mod.size()-2;
	}
//private:
//	string& nominal_model;
};

struct SortFit{
	bool operator() (const ind& i,const ind& j) { return (i.fitness<j.fitness);} 
};

struct SortRank{
	bool operator() (const ind& i,const ind& j) { return (i.rank<j.rank);} 
};

struct revSortRank{
	bool operator() (ind& i,ind& j) { return (i.rank>j.rank);} 
};
struct SortEqnSize{
	bool operator() (const ind& i,const ind& j) { return (i.eqn.size()<j.eqn.size());} 
};
struct SortFEVar{
	bool operator() (const ind& i,const ind& j) { return (i.FEvar>j.FEvar);} 
};
struct SortComplexity{bool operator() (const ind& i,const ind& j) { return (i.complexity<j.complexity);} 
};
struct SortFit_v{
	bool operator() (const ind& i,const ind& j) { return (i.fitness_v<j.fitness_v);}
};
struct SortSize{
	bool operator() (const ind& i,const ind& j) { return (i.line.size()<j.line.size());}
};

struct sameEqn{
	bool operator() (ind& i,ind& j) { return i.eqn==j.eqn;} 
};

struct sameEqnSize{
	bool operator() (ind& i,ind& j) { return i.eqn.size()==j.eqn.size();} 
};

struct sameSizeFit{
	bool operator() (ind& i,ind& j) { return (i.fitness==j.fitness && i.eqn.size()==j.eqn.size());} 
};

struct sameFit{
	bool operator() (ind& i,ind& j) { return (i.fitness==j.fitness);} 
};

struct sameComplexity{bool operator() (const ind& i,const ind& j) { return (i.complexity==j.complexity);} 
};

struct sameFitComplexity{bool operator() (const ind& i,const ind& j) { return (i.fitness==j.fitness && i.complexity==j.complexity);} 
};


struct tribe{ 

	vector <ind> pop; // population
	float best;
	float worst;

	tribe(int size,float& max_fit,float& min_fit)	
	{
		pop.resize(size);
		best=max_fit;
		worst=min_fit;
		/*for(unsigned int i = 1;i<pop.size();i++)
			pop.at(i).init(nom_mod);*/
		maxf=max_fit;
		minf=min_fit;
	}
	~tribe() {}

	float bestFit() // returns best fitness value
	{

		/*#pragma omp parallel
		{
		   float localbest = maxf;

		   #pragma omp for schedule(static)
		   for(int i = 0; i < pop.size(); i++)
			   localbest = min(localbest, pop.at(i).fitness);

		   #pragma omp critical
		   {
			  best = min(localbest, best);
		   }
		}*/
		best = maxf;
		for(int i = 0; i < pop.size(); i++)
			   best = min(best, pop.at(i).fitness);
		return best;
	}
	float bestFit_v() // returns best fitness value
	{

		/*#pragma omp parallel
		{
		   float localbest = maxf;

		   #pragma omp for schedule(static)
		   for(int i = 0; i < pop.size(); i++)
			   localbest = min(localbest, pop.at(i).fitness);

		   #pragma omp critical
		   {
			  best = min(localbest, best);
		   }
		}*/
		best = maxf;
		for(int i = 0; i < pop.size(); i++)
			   best = min(best, pop.at(i).fitness_v);
		return best;
	}
	float worstFit() //worst fitness
	{
		worst = minf;
		/*#pragma omp parallel
		{
		   float localworst = minf;

		   #pragma omp for schedule(static)
		   for(int i = 0; i < pop.size(); i++)
			   localworst = max(localworst, pop.at(i).fitness);

		   #pragma omp critical
		   {
			  worst = max(localworst, worst);
		   }
		}*/
		 for(int i = 0; i < pop.size(); i++)
			 worst = max(worst, pop.at(i).fitness);
		return worst;
	}	
	float medFit() //median fitness
	{
		//vector<ind> tmppop = pop;
		sort(pop.begin(),pop.end(),SortFit());
		int index = (int)floor((float)pop.size()/2);
		return pop.at(index).fitness; 

	} 
	float medFit_v() //median fitness
	{
		//vector<ind> tmppop = pop;
		sort(pop.begin(),pop.end(),SortFit_v());
		int index = (int)floor((float)pop.size()/2);
		return pop.at(index).fitness_v; 

	}
	float meanFit() // mean fitness
	{
		float answer=0;
		//#pragma omp parallel for reduction(+ : answer)
		for(int i=0; i<pop.size(); i++)
		{
			answer+=pop.at(i).fitness;
		}
		return (float)answer/pop.size();
	}
			
	float meanSize() // mean line length
	{
		float answer=0;
		//#pragma omp parallel for reduction(+ : answer)
		for(int i=0; i<pop.size(); i++)
		{
			answer+=pop.at(i).line.size();
		}
		return (float)answer/pop.size();
	}
	float meanEffSize()
	{
		float answer=0;
		//#pragma omp parallel for reduction(+ : answer)
		for(int i=0; i<pop.size(); i++)
		{
			answer+=pop.at(i).eff_size;
		}
		return (float)answer/pop.size();
	}
	int medSize() // median line length
	{		
		//vector<ind> tmppop = pop;
		sort(pop.begin(),pop.end(),SortSize());
		int index = (int)floor((float)pop.size()/2);
		return pop.at(index).line.size(); 
	}

	void topTen(vector <ind>& eqns) //returns address to vector of equation strings
	{
		vector<ind> tmppop = pop;
		sort(tmppop.begin(),tmppop.end(),SortFit());
		unique(tmppop.begin(),tmppop.end(),sameFit());
		
		for (int i=0;i<10;i++)
			eqns.push_back(tmppop.at(i));

		/*vector <float> fitnesses;
		int i=0;
		bool pass=true;
		while(eqns.size()<10 && i<pop.size())
		{
			fitnesses.push_back(pop.at(i).fitness);
			for(unsigned int j=0;j<fitnesses.size()-1;j++)		 
			{
				if(fitnesses.at(j)==fitnesses.back())
				{
					fitnesses.pop_back();
					pass=0;
					break;
				}
			}

			if (pass)
				eqns.push_back(pop.at(i));
			else
				pass=1;
			++i;
		}*/
	}
	void getbestind(ind& bestind)
	{
		//vector<ind> tmppop = pop;
		sort(pop.begin(),pop.end(),SortFit());
		bestind = pop.front();
	}// address of best individual
	void sortpop()
	{
		sort(pop.begin(),pop.end(),SortFit());
	}
	/*
private:
	
	bool fitlow (ind& i,ind& j) { return (i.fitness<j.fitness); }
	bool eqncomp(ind& i,ind& j) { return (i.eqn_form.compare(j.eqn_form)==0); }
	bool fitcomp (ind& i,ind& j) { return (i.fitness==j.fitness); }
	*/

private:
		float maxf;
		float minf;


};


#endif
