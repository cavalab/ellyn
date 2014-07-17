#include "stdafx.h"
#include "pop.h"
#include "params.h"
#include "rnd.h"
#include "data.h"
#include "Fitness.h"
#include "state.h"
#include "FitnessEstimator.h"
#include "general_fns.h"

//class FitnessEstimator{
//public:
//	vector <int> FEpts; // points of fitness estimation (subset from data vals)
//	vector <float> TrainerFit; // fitness on trainer population
//	float fitness; // self fitness across trainers
//
//	FitnessEstimator(int length,vector<Randclass>& r,data& d)
//	{
//		for (int i=0;i<length;i++){
//			FEpts.push_back(r[omp_get_thread_num()].rnd_int(0,d.vals.size()));
//		}
//	}
//
//	
//};
struct SortFE{
	bool operator() (const FitnessEstimator& i,const FitnessEstimator& j) { return (i.fitness<j.fitness);} 
};

void setFEvals(vector<vector<float>>& FEvals, vector<float>& FEtarget,FitnessEstimator& FE, data& d)
{
	FEvals.clear();
	FEtarget.clear();
	for (int i=0;i<FE.FEpts.size();i++)
	{
		FEvals.push_back(d.vals[FE.FEpts[i]]);
		FEtarget.push_back(d.target[FE.FEpts[i]]);
	}
};
void FitnessFE(vector <FitnessEstimator>& FE, vector <ind>& trainers,params& p,data& d,state& s)
{
	// calculate fitness of each Fitness Estimator on each trainer. 
	/*vector<float> exactfit(trainers.size());
		for (int j=0;j<trainers.size();j++)
		exactfit[j]=trainers[j].fitness;*/

	//get trainer exact fitness
	p.EstimateFitness = 0; // turn fitness estimation off
	Fitness(trainers,p,d,s,FE[0]); //note: FE[0] does not get used here
	p.EstimateFitness = 1; // turn fitness estimation back on

	vector<vector<float>> FEvals;
	vector<float> FEtarget;
	int ndata_t;
	for (int i=0;i<FE.size();i++){
		
		setFEvals(FEvals,FEtarget,FE[i],d);

		if(p.train)
			ndata_t = FEvals.size()/2;
		else
			ndata_t = FEvals.size();

		for (int j=0;j<trainers.size();j++){

			float abserror=0;
			float meantarget = 0;
			float meanout = 0;
			float corr;
			float FEness;
			float target_std;
			vector<float> tmpoutput;
			for(unsigned int sim=0;sim<ndata_t;sim++)
			{	
				abserror += abs(FEtarget[sim]-trainers[j].output[FE[i].FEpts[sim]]);
				meantarget += FEtarget.at(sim);
				meanout += trainers[j].output[FE[i].FEpts[sim]];
				tmpoutput.push_back(trainers[j].output[FE[i].FEpts[sim]]);
			}
				
				// mean absolute error
				abserror = abserror/ndata_t;
				meantarget = meantarget/ndata_t;
				meanout = meanout/ndata_t;
				//calculate correlation coefficient
				corr = getCorr(tmpoutput,FEtarget,meanout,meantarget,0,target_std);
				tmpoutput.clear();

				if(corr < p.min_fit)
					corr=p.min_fit;

				if(trainers[j].output.empty()) 
					FEness=p.max_fit;
				/*else if (*std::max_element(pop.at(count).output.begin(),pop.at(count).output.end())==*std::min_element(pop.at(count).output.begin(),pop.at(count).output.end()))
					pop.at(count).fitness=p.max_fit;*/
				else if ( boost::math::isnan(abserror) || boost::math::isinf(abserror) || boost::math::isnan(corr) || boost::math::isinf(corr))
					FEness=p.max_fit;
				else{
					if (p.fit_type==1)
						FEness = abserror;
					else if (p.fit_type==2)
						FEness = 1-corr;
					else if (p.fit_type==3)
						FEness = abserror/corr;
					if (p.norm_error)
						FEness = FEness/target_std;
				}

				if(FEness>p.max_fit)
					FEness=p.max_fit;
				else if(FEness<p.min_fit)
					(FEness=p.min_fit); 

		//Fitness(trainers,p,d,s,FE[i]); //note: FE[0] does not get used here
		
			FE[i].fitness+=abs(trainers[j].fitness-FEness);
		}
		FE[i].fitness=FE[i].fitness/trainers.size();
	}
	
};
void InitPopFE(vector <FitnessEstimator>& FE,vector<ind> &pop,vector<ind>& trainers,params p,vector<Randclass>& r,data& d,state& s)
{
	//vector <FitnessEstimator> FE; //fitness estimator population 
	//initialize estimator population
	FE.resize(p.FE_pop_size);
	for (int i=0;i<p.FE_pop_size;i++){
		FE[i]= FitnessEstimator(p.FE_ind_size,r,d,p.train);			
	}
	
	//initialize trainer population
	trainers.resize(p.FE_train_size);
	for (int i=0;i<p.FE_train_size;i++){
		trainers[i] = pop[r[omp_get_thread_num()].rnd_int(0,pop.size()-1)];
	}
	//p.EstimateFitness = 0; // turn fitness estimation off (all data points used)
	//Fitness(trainers,p,d,s,FE[0]); //note: FE[0] does not get used here
	//p.EstimateFitness = 1; // turn fitness estimation back on
	//evaluate fitness of estimators
	FitnessFE(FE,trainers,p,d,s);
	//perform selection
	sort(FE.begin(),FE.end(),SortFE());
};


void PickTrainers(vector<ind> pop, vector <FitnessEstimator>& FE,vector <ind>& trainers,params& p,data& d,state& s)
{
	vector <vector<float>> FEfits; //rows: predictors, cols: population
	vector<float> meanfits(pop.size());
	vector <float> varfits(pop.size());
	//vector <ind> tmppop(1);
	// get FE fitness on each individual in population
	vector<ind> tmppop(pop.size());

	for (int j=0;j<pop.size();j++){
		makenewcopy(pop[j],tmppop[j]);
		tmppop[j].clrPhen();
	}
	//i: fitness predictors
	//j: solution population 
	for (int i=0;i<FE.size();i++){

		Fitness(tmppop,p,d,s,FE[i]);	
		FEfits.push_back(vector<float>(tmppop.size()));
		for (int j=0;j<tmppop.size();j++){
				FEfits[i][j]=tmppop[j].fitness;	
				tmppop[j].clrPhen();
			}
		
	}

	// calculate variance in fitness estimates
	for (int j=0;j<tmppop.size();j++){
		for (int i=0;i<FE.size();i++)
			meanfits[j]+=FEfits[i][j];

		meanfits[j]=meanfits[j]/FE.size(); // mean fitness over predictors

		for (int i=0;i<FE.size();i++)
			varfits[j]+=pow(FEfits[i][j]-meanfits[j],2);

		varfits[j]=varfits[j]/FE.size(); // variance in fitness over predictors
		tmppop[j].FEvar=varfits[j];
	}
	
	sort(tmppop.begin(),tmppop.end(),SortFEVar());

	//vector<ind>::iterator it = pop.begin();
	trainers.clear();
	trainers.assign(tmppop.begin(),tmppop.begin()+p.FE_train_size);

	//evaluate fitness of estimators
	FitnessFE(FE,trainers,p,d,s);
	sort(FE.begin(),FE.end(),SortFE());
};
void crossFE(vector <FitnessEstimator>& FE,vector<Randclass>& r)
{
	// select parents
	int nt = 2;
	vector<int> parents;
	vector<int> choices(2);
	vector <FitnessEstimator> newFE;
	float bestfit; 
	int bestchoice;
	for (int i=0; i < FE.size(); i++){
		
		for (int j=0; j<nt; j++){
			
			choices[j] = r[omp_get_thread_num()].rnd_int(0,FE.size()-1);
			if (j==0) {
				bestfit = FE[choices[0]].fitness;
				bestchoice = choices[0];
			}
			else if(FE[choices[j]].fitness<bestfit)
			{
				bestfit = FE[choices[j]].fitness;
				bestchoice = choices[j];
			}
		}
		parents.push_back(bestchoice);		
	}

	//// new pop
	newFE.resize(FE.size());
	for (int i=0;i<parents.size();i++){
		//newFE[i]= FE[parents[i]];	
		newFE[i].FEpts.resize(FE[i].FEpts.size());
	}
	//cross parents
	int p1,p2;
	int crosspt;
	for (int i=0; i < FE.size()-1; i++){
		int point1 = r[omp_get_thread_num()].rnd_int(0,FE[0].FEpts.size()-1);
		
		newFE[i].FEpts.assign(FE[parents[i]].FEpts.begin(),FE[parents[i]].FEpts.begin()+point1);
		newFE[i].FEpts.insert(newFE[i].FEpts.end(),FE[parents[i+1]].FEpts.begin()+point1,FE[parents[i+1]].FEpts.end());

		newFE[i+1].FEpts.assign(FE[parents[i+1]].FEpts.begin(),FE[parents[i+1]].FEpts.begin()+point1);
		newFE[i+1].FEpts.insert(newFE[i+1].FEpts.end(),FE[parents[i]].FEpts.begin()+point1,FE[parents[i]].FEpts.end());

		i++;
	}
	FE = newFE;
	// keep strictly the best predictors	
	/*FE.insert(FE.end(),newFE.begin(),newFE.end());
	sort(FE.begin(),FE.end(),SortFE());
	FE.erase(FE.begin()+newFE.size(),FE.end());*/
};
void mutateFE(vector <FitnessEstimator>& FE,params p,data& d,vector<Randclass>& r)
{
	int pt;
	int lastpt;
	if (p.train) lastpt = d.vals.size()/2-1;
	else lastpt = d.vals.size()-1;

	for (int i=0; i < FE.size(); i++){
		pt = r[omp_get_thread_num()].rnd_int(0,FE[i].FEpts.size()-1);
		FE[i].FEpts[pt] = r[omp_get_thread_num()].rnd_int(0,lastpt);
	}
};
void EvolveFE(vector<ind> &pop, vector <FitnessEstimator>& FE,vector <ind>& trainers,params p,data& d,state& s,vector<Randclass>& r) 
{
	vector <float> FEfitness(FE.size()); //fitness of Fitness Estimators
	
	vector <FitnessEstimator> newFE = FE;
	//cross estimators
	crossFE(newFE,r);

	//mutate estimators
	mutateFE(newFE,p,d,r);
	
	//evaluate fitness of estimators
	FitnessFE(newFE,trainers,p,d,s);

	// keep strictly the best predictors	
	/*FE.insert(FE.end(),newFE.begin(),newFE.end());
	sort(FE.begin(),FE.end(),SortFE());
	FE.erase(FE.begin()+newFE.size(),FE.end());*/

	//keep the elite predictor and the new predictors
	newFE.insert(newFE.end(),FE[0]);
	sort(newFE.begin(),newFE.end(),SortFE());
	FE.assign(newFE.begin(),newFE.end()-1);
	////perform selection
	//sort(FE.begin(),FE.end(),SortFE());
	
	/*//if time to add new fitness trainer
	PickTrainers(pop,FE,trainers,p,d,s);
	*/
};