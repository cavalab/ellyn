// header file for ind struct
#pragma once
#ifndef POP_H
#define POP_H
#include <iostream>
#include<string>
#include<vector>
#include "params.h"
#include "data.h"


struct ind {
	//unsigned int id;
	std::vector <int> line;
	std::vector <bool> epiline;
	std::vector <std::string> args;
	std::string eqn;
	std::string eqn_form; // equation form for string distance comparison to other forms

	std::vector <int> ptr;
	char origin; // x: crossover, m: mutation, i: initialization

	std::vector<float> output;
	float abserror;
	float corr;
	float fitness;

	float parentfitness;
	int eff_size;
	int age;
	//typedef exprtk::expression<float> expression_t;
	//expression_t expression;
	/*ind()
	{
		eqn=nominal_model;
	}*/
	ind()
	{
		age=1;
	}
	//void init(string& nom_mod)
	//{
	//	eqn = nom_mod;
	//	ptr.push_back(1);
	//	ptr.push_back(nom_mod.size()-2);
	//	//nominal_model=nom_mod;
	//	//expression.register_symbol_table(d.symbol_table);
	//}
	void clrPhen(string& nom_mod)
	{
		eqn = nom_mod;
		output.clear();
		// nominal model must be encased in set of parenthesis. the pointer points to that which is encased.
		//ptr[0]= 1;
		//ptr[1] = nom_mod.size()-2;
	}
//private:
//	string& nominal_model;
};

struct SortFit{
	bool operator() (ind& i,ind& j) { return (i.fitness<j.fitness);} 
};
struct SortSize{
	bool operator() (ind& i,ind& j) { return (i.line.size()<j.line.size());} 
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
		//vector<ind> tmppop = pop;
		sort(pop.begin(),pop.end(),SortFit());
		vector <float> fitnesses;
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
		}
	}
	void getbestind(ind& bestind)
	{
		//vector<ind> tmppop = pop;
		sort(pop.begin(),pop.end(),SortFit());
		bestind = pop.front();
	}// address of best individual
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