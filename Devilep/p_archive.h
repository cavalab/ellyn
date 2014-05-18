#pragma once
#include "stdafx.h"
//#include "pop.h"

struct SortRank{
	bool operator() (ind& i,ind& j) { return (i.rank<j.rank);} 
};
struct revSortRank{
	bool operator() (ind& i,ind& j) { return (i.rank>j.rank);} 
};
struct sameEqn{
	bool operator() (ind& i,ind& j) { return i.eqn.compare(j.eqn)==0;} 
};
struct sameSizeFit{
	bool operator() (ind& i,ind& j) { return (i.fitness==j.fitness && i.complexity==j.complexity);} 
};
//struct SortFit{
//	bool operator() (ind& i,ind& j) { return (i.fitness<j.fitness);} 
//};
struct paretoarchive{
	vector <ind> pop; // population
	int archsize;
	paretoarchive(int size){archsize=size;}
	~paretoarchive(){}

	void update(vector<ind>& newpop)
	{
		vector <ind> tmppop;
		
		
		pop.insert(pop.end(),newpop.begin(),newpop.end());
		sort(pop.begin(),pop.end(),SortFit());
		vector<ind>::iterator it;
		it = std::unique (pop.begin(), pop.end(), sameSizeFit());
		pop.resize(distance(pop.begin(),it));

		int r = 0;
		while (tmppop.size()<archsize && !pop.empty())
		{
			pareto(pop);
			sort(pop.begin(),pop.end(),SortRank());
			while(!pop.empty() && pop.back().rank==1)
			{
				pop.back().rank+=r;
				tmppop.push_back(pop.back());
				pop.pop_back();
			}
			r++;
		}

		pop = tmppop;
		sort(pop.begin(),pop.end(),SortFit());
	}

};