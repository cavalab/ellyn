#pragma once
#ifndef STATE_H
#define STATE_H
#include "logger.h"
//#include "pop.h"

struct state{

	vector<long long> ptevals; //number of instruction evals
	vector<long long> numevals; //number of evals on each thread
	vector<int> genevals; //total evals for each generation
	vector <int> fit_best;
	vector <float> fit_mean;
	vector <float> fit_med;
	vector <float> fit_std;
	vector <float> size_mean;
	vector <int> size_med;
	vector <float> size_std;
	vector <float> eff_size;
	vector<int> eHC_updates;
	int total_eHC_updates;
	float current_eHC_updates;
	vector<int> pHC_updates;
	float current_pHC_updates;
	int total_pHC_updates;
	vector<int> good_cross;
	vector<int> bad_cross;
	vector<int> neut_cross;
	float good_cross_pct;
	float neut_cross_pct;
	float bad_cross_pct; 

	logger out; 

	state()
	{
		ptevals.resize(omp_get_max_threads());
		numevals.resize(omp_get_max_threads());
		genevals.push_back(0);
		eHC_updates.assign(omp_get_max_threads(),0);
		pHC_updates.assign(omp_get_max_threads(),0);
		good_cross.assign(omp_get_max_threads(),0);
		bad_cross.assign(omp_get_max_threads(),0);
		neut_cross.assign(omp_get_max_threads(),0);
		total_eHC_updates=0;
		total_pHC_updates=0;
		good_cross_pct=0;
		neut_cross_pct=0;
	}
	~state() {}

	int getgenevals()
	{
		return genevals.back();
	}
	void setgenevals()
	{
		int gentmp=0;
		for(unsigned int i = 0;i<numevals.size();i++)
			gentmp+=numevals.at(i);

		genevals.push_back(gentmp - totalevals());

	}
	long long totalevals()
	{
		long long te=0;
		for(unsigned int i = 0;i<genevals.size();i++)
			te+=genevals.at(i);
		return te;
	}
	long long totalptevals()
	{
		long long te=0;
		for(unsigned int i = 0;i<ptevals.size();i++)
			te+=ptevals.at(i);
		return te;
	}
	int setPHCupdates()
	{
		int updates=0;
		for(unsigned int i =0;i<pHC_updates.size();i++)
			updates+=pHC_updates[i];

		int val = (updates-total_pHC_updates);
		total_pHC_updates+=val;
		current_pHC_updates = float(val);
		return val; 

	
	}
	int setEHCupdates()
	{
		int updates=0;
		for(unsigned int i =0;i<eHC_updates.size();i++)
			updates+=eHC_updates[i];

		int val = (updates-total_eHC_updates);
		total_eHC_updates+=val;
		current_eHC_updates = float(val);
		return val; 
	}

	float getGoodCrossPct()
	{
		
		float total_good=0;
		float total=0;

		for(unsigned int i=0;i<good_cross.size();i++)
		{
			total_good+=good_cross.at(i);
			total += good_cross.at(i)+bad_cross.at(i)+neut_cross.at(i);
		}
		//good_cross.assign(omp_get_max_threads(),0);
		if (total==0){
			good_cross_pct=0;
			return 0;
		}
		else
		{
			good_cross_pct = total_good/float(total)*100;	
			return good_cross_pct;
		}
		
	}
	float getNeutCrossPct()
	{
	
		float total_neut=0;
		float total=0;

		for(unsigned int i=0;i<neut_cross.size();i++)
		{
			total_neut+=neut_cross.at(i);
			total += neut_cross.at(i)+bad_cross.at(i)+good_cross.at(i);
		}
		//neut_cross.assign(omp_get_max_threads(),0);
		if (total==0){
			neut_cross_pct = 0;
			return 0;
		}
		else
		{
			neut_cross_pct = total_neut/total*100;			
			return neut_cross_pct;
		}

	}
	float getBadCrossPct()
	{
	
		float total_bad=0;
		float total=0;

		for(unsigned int i=0;i<neut_cross.size();i++)
		{
			total_bad+=bad_cross.at(i);
			total += neut_cross.at(i)+bad_cross.at(i)+good_cross.at(i);
		}
		clearCross();
		if (total==0){
			bad_cross_pct = 0;
			return 0;
		}
		else
		{
			bad_cross_pct = total_bad/total*100;
			
			return bad_cross_pct;
		}

	}
	void clearCross()
	{
		good_cross.assign(omp_get_max_threads(),0);
		bad_cross.assign(omp_get_max_threads(),0);
		neut_cross.assign(omp_get_max_threads(),0);

	}
	void setCrossPct(vector<ind> pop)
	{
		
		for(int i=0;i<pop.size();i++)
		{
			if (pop.at(i).parentfitness > pop.at(i).fitness)
				good_cross[omp_get_thread_num()]=good_cross[omp_get_thread_num()]+1;
			else if(pop.at(i).parentfitness == pop.at(i).fitness)
				neut_cross[omp_get_thread_num()]=neut_cross[omp_get_thread_num()]+1;
			else
				bad_cross[omp_get_thread_num()]=bad_cross[omp_get_thread_num()]+1;

		}
	}
	void clear()
	{
		ptevals.clear();
		numevals.clear();
		genevals.clear();
		fit_best.clear();
		fit_mean.clear();
		fit_med.clear();
		fit_std.clear();
		size_mean.clear();
		size_med.clear();
		size_std.clear();
		ptevals.resize(omp_get_max_threads());
		numevals.resize(omp_get_max_threads());
		genevals.push_back(0);
	}
};

#endif