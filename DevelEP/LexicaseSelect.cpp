#include "stdafx.h"
#include "pop.h"
#include "params.h"
#include "state.h"
#include "rnd.h"
#include "data.h"
#include "FitnessEstimator.h"
#include "Generationfns.h"
#include "Fitness.h"

float std_dev(vector<float>& x) {

	// get mean of x
	float mean = std::accumulate(x.begin(), x.end(), 0.0)/x.size();
	//calculate variance
	float var = 0;
	for (vector<float>::iterator it = x.begin(); it != x.end(); ++it)
		var += pow(*it - mean,2);

	var /= (x.size()-1);
	return sqrt(var);

}
float median(vector<float> x) {
	sort(x.begin(), x.end());

	if (x.size() % 2 == 0) {
		float tmp1 = x[(x.size() / 2) - 1];
		float tmp2 = x[x.size() / 2];
		return (tmp1 + tmp2) / 2;
	}
	else
		return x[x.size() / 2];
}
float median(vector<int> x) {
	sort(x.begin(), x.end());

	if (x.size() % 2 == 0) {
		int tmp1 = x[(x.size() / 2) - 1];
		int tmp2 = x[x.size() / 2];
		return float((tmp1 + tmp2) / 2);
	}
	else
		return float(x[x.size() / 2]);
}
float mad(vector<float>& x) {
	// returns median absolute deviation (MAD)
	// get median of x
	float x_median = median(x);

	//calculate absolute deviation from median
	vector<float> dev;
	for (vector<float>::iterator it = x.begin(); it != x.end(); ++it)
		dev.push_back(abs(*it - x_median));

	return median(dev);

}
void LexicaseSelect(vector<ind>& pop,vector<unsigned int>& parloc,params& p,vector<Randclass>& r,Data& d,state& s)
{
	// add metacases if needed
	for (unsigned i = 0; i<p.lex_metacases.size(); ++i)
	{
		if (p.lex_metacases[i].compare("age")==0){
			for (unsigned j=0;j<pop.size(); ++j)
				pop[j].error.push_back(pop[j].age);
		}
		else if (p.lex_metacases[i].compare("complexity")==0){
			for (unsigned j=0;j<pop.size(); ++j)
				pop[j].error.push_back(pop[j].complexity);
		}
		else if (p.lex_metacases[i].compare("dimensionality") == 0) {
			for (unsigned j = 0; j<pop.size(); ++j)
				pop[j].error.push_back(pop[j].dim);
		}
	}


	//boost::progress_timer timer;
	//vector<float> fitcompare;
	vector<int> pool; // pool from which to choose parent. normally it is set to the whole population.
	int numcases;
	if (p.lex_class)
		numcases = p.number_of_classes + p.lex_metacases.size();
	else if (p.EstimateFitness)
		numcases = p.FE_ind_size*p.train_pct + p.lex_metacases.size();
	else
		numcases = d.vals.size()*p.train_pct + p.lex_metacases.size();

	if (p.lexpool==1){
		for (int i=0;i<pop.size();++i){
			if(pop[i].error.size() == numcases)
				pool.push_back(i);
		}
	}

	vector <int> starting_pool = pool;
	vector <int> case_order;
	for(int i=0;i<pop[pool[0]].error.size();++i)
		case_order.push_back(i);
	assert(case_order.size()>0);
	float minfit=p.max_fit;
	vector<int> winner;
	bool draw=true;
	bool pass;
	int h;
	vector<int> num_passes(numcases-p.lex_metacases.size(),0);
	int tmp;
  vector<float> epsilon; // epsilon values for each case

	// epsilon lexicase implementations
	if (p.lex_eps_error ) // errors within episilon of the global best error are pass, otherwise fail
	{
		if(p.lex_eps_global){
		// get minimum error on each case
			vector<float> min_error(numcases - p.lex_metacases.size(),p.max_fit);

			for (size_t i = 0; i < numcases - p.lex_metacases.size(); ++i) {
				for (size_t j = 0; j < pool.size(); ++j) {
					if (pop[pool[j]].error[i] < min_error[i])
						min_error[i] = pop[pool[j]].error[i];
				}
			}
			for (size_t i = 0; i < pool.size(); ++i) {
				// check if error is within epsilon
				for (size_t j = 0; j < numcases - p.lex_metacases.size(); ++j) {
					if (pop[pool[i]].error[j] <= (1+p.lex_epsilon)*min_error[j]){
						pop[pool[i]].error[j] = 0;
						++num_passes[j];
					}
					else
						pop[pool[i]].error[j] = 1;
				}
			}
		}
		else{// for non-global version, just save the epsilon values for each case
			for (size_t i = 0; i < numcases - p.lex_metacases.size(); ++i)
			 	epsilon.push_back(p.lex_epsilon);
		}
	}
	else if (p.lex_eps_target) // errors within epsilon of the target are pass, otherwise fail
	{
		if (p.lex_eps_global){
			for (size_t i = 0; i < pool.size(); ++i) {
				// check if error is within epsilon
				for (size_t j = 0; j < numcases - p.lex_metacases.size(); ++j) {
					if (pop[pool[i]].error[j] <= p.lex_epsilon)
						pop[pool[i]].error[j] = 0;
					else
						pop[pool[i]].error[j] = 1;
				}
			}
		}
		else{// for non-global version, just save the epsilon values for each case
			for (size_t i = 0; i < numcases - p.lex_metacases.size(); ++i)
			 	epsilon.push_back(p.lex_epsilon);
		}
	}

	else if (p.lex_eps_std) // errors in a standard dev of the best are pass, otherwise fail
	{
		// get minimum error on each case
		vector<float> min_error(numcases - p.lex_metacases.size(), p.max_fit);
		vector<float> std_error(numcases - p.lex_metacases.size());

		for (size_t i = 0; i < numcases - p.lex_metacases.size(); ++i) {
			vector<float> case_error(pool.size());
			for (size_t j = 0; j < pool.size(); ++j) {
				if (pop[pool[j]].error[i] < min_error[i])
					min_error[i] = pop[pool[j]].error[i];

				case_error[j] = pop[pool[j]].error[i];
			}
			//std_error[i] = 0.5*std_dev(case_error)/sqrt(float(pool.size()));
			std_error[i] = std_dev(case_error);

		}
		if (p.lex_eps_global){ // if pass conditions defined relative to whole population (rather than selection pool):

			for (size_t i = 0; i < pool.size(); ++i) {
				// check if error is within epsilon
				for (size_t j = 0; j < numcases - p.lex_metacases.size(); ++j) {
					if (pop[pool[i]].error[j] <= min_error[j]+std_error[j])
						pop[pool[i]].error[j] = 0;
					else
						pop[pool[i]].error[j] = 1;
				}

			}
		}
		else{// for non-global version, just save the epsilon values for each case
				for (size_t i = 0; i < numcases - p.lex_metacases.size(); ++i)
				 	epsilon.push_back(std_error[i]);
		}
	}
	else if (p.lex_eps_target_mad) { // errors within mad of the target pass
									 // get minimum error on each case
		vector<float> mad_error(numcases - p.lex_metacases.size());

		for (size_t i = 0; i < numcases - p.lex_metacases.size(); ++i) {
			vector<float> case_error(pool.size());
			for (size_t j = 0; j < pool.size(); ++j)
				case_error[j] = pop[pool[j]].error[i];

			mad_error[i] = mad(case_error);
		}
		if (p.lex_eps_global){ // if pass conditions defined relative to whole population (rather than selection pool):
			for (size_t i = 0; i < pool.size(); ++i) {
				// check if error is within epsilon
				for (size_t j = 0; j < numcases - p.lex_metacases.size(); ++j) {
					if (pop[pool[i]].error[j] <= mad_error[j])
						pop[pool[i]].error[j] = 0;
					else
						pop[pool[i]].error[j] = 1;
				}
			}
		}
		else{// for non-global version, just save the epsilon values for each case
				for (size_t i = 0; i < numcases - p.lex_metacases.size(); ++i)
					epsilon.push_back(mad_error[i]);
		}
	}
	else if (p.lex_eps_error_mad) { // errors within mad of the best error pass
									 // get minimum error on each case
		vector<float> min_error(numcases - p.lex_metacases.size(), p.max_fit);
		vector<float> mad_error(numcases - p.lex_metacases.size());

		for (size_t i = 0; i < numcases - p.lex_metacases.size(); ++i) {
			vector<float> case_error(pool.size());
			for (size_t j = 0; j < pool.size(); ++j) {
				if (pop[pool[j]].error[i] < min_error[i])
					min_error[i] = pop[pool[j]].error[i];

				case_error[j] = pop[pool[j]].error[i];
			}
			mad_error[i] = mad(case_error);
		}
		if (p.lex_eps_global){ // if pass conditions defined relative to whole population (rather than selection pool):
				for (size_t i = 0; i < pool.size(); ++i) {
					// check if error is within epsilon
					for (size_t j = 0; j < numcases - p.lex_metacases.size(); ++j) {
						if (pop[pool[i]].error[j] <= min_error[j]+mad_error[j])
							pop[pool[i]].error[j] = 0;
						else
							pop[pool[i]].error[j] = 1;
					}
				}
		}
		else{ // for non-global version, just save the epsilon values for each case
				for (size_t i = 0; i < numcases - p.lex_metacases.size(); ++i)
					epsilon.push_back(mad_error[i]);
		}
	}
	// measure median number of cases used
	vector<float> cases_used;
	// measure median pool size at selection
	vector<float> sel_size;
	// for each selection event:
	for (int i=0;i<parloc.size();++i)
	{
		//shuffle test cases
		std::random_shuffle(case_order.begin(),case_order.end(),r[omp_get_thread_num()]);

		// select a subset of the population if lexpool is being used
		if (p.lexpool!=1){
			for (int j=0;j<p.lexpool*pop.size();++j){
				tmp = r[omp_get_thread_num()].rnd_int(0,pop.size()-1);
				int n = 0;
				while (pop[tmp].error.empty() && n < pop.size()) {
					tmp = r[omp_get_thread_num()].rnd_int(0, pop.size() - 1);
					n += 1;
				}
				pool.push_back(tmp);
			}
		}
		else //otherwise use the whole population for selection
			pool = starting_pool;

		pass=true;
		h=0; // count the number of cases that are used
		while ( pass && h<case_order.size()) //while there are remaining cases and more than 1 unique choice
		{
			//reset winner and minfit for next case
			winner.resize(0);
			minfit=p.max_fit;

			if (p.lex_eps_global){
				// loop through the individuals and pick out the elites
				for (int j=0;j<pool.size();++j)
				{

					if (pop[pool[j]].error[case_order[h]]<minfit)
					{
						minfit=pop[pool[j]].error[case_order[h]];
						winner.resize(0);
						winner.push_back(pool[j]);
					}
					else if (pop[pool[j]].error[case_order[h]]==minfit)
					{
						winner.push_back(pool[j]);
					}
				}
			}
			else{// epsilon lexicase with local pool pass conditions
				// get best fitness in selection pool
				for (int j=0;j<pool.size();++j){
					if (pop[pool[j]].error[case_order[h]]<minfit)
						minfit=pop[pool[j]].error[case_order[h]];
				}
				// winners are within epsilon of the local pool's minimum fitness
				for (int j=0;j<pool.size();++j){
					if (pop[pool[j]].error[case_order[h]]<=minfit+epsilon[case_order[h]])
						winner.push_back(pool[j]);
				}
			}
			// if there is more than one elite individual and still more cases to consider
			if(winner.size()>1 && h<case_order.size()-1)
			{
				pass=true;
				++h;
				//reduce pool to elite individuals on case case_oder[h]
				pool = winner;
			} // otherwise, a parent has been chosen or has to be chosen randomly from the remaining pool
			else
				pass=false;


		}
		//if more than one winner, pick randomly
		if(winner.size()>1)
			parloc[i]=winner[r[omp_get_thread_num()].rnd_int(0,winner.size()-1)];
		else if (winner.size()==1) // otherwise make the winner a parent
			parloc[i]=winner[0];
		else{ // otherwise throw an ??
			cout << "??";
			cout << "winner is not correct size\n";
			cout << "minfit: " << minfit << "\n";
			cout << "epsilon:\n";
			for (size_t i = 0; i< epsilon.size(); ++i)
				cout << epsilon[i] << ",";
		}
		assert(winner.size()>0);
		// reset minfit

		minfit=p.max_fit;
		cases_used.push_back(h+1);
		sel_size.push_back(winner.size());
	}//for (int i=0;i<parloc.size();++i)

	std::sort(cases_used.begin(), cases_used.end());
	std::sort(sel_size.begin(), sel_size.end());

	// store median lex cases used, normalized by the number of cases
	s.median_lex_cases[omp_get_thread_num()] = median(cases_used)/case_order.size();

	// store median number of individuals that pass each case
	s.median_passes_per_case[omp_get_thread_num()] = median(num_passes) / pop.size();

	// store median lex pool size at selection, normalized by the population size
	s.median_lex_pool[omp_get_thread_num()] = median(sel_size)/pop.size();

	//cout << "mean num cases used: " << hmean << "\n";

}
