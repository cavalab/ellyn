// header file for ind struct
#pragma once
#ifndef PARAMS_H
#define PARAMS_H
//#include <iostream>
//#include <string>
//#include <vector>
//#include <random>
//#include <array>
#include "op_node.h"
using namespace std;

struct params {
	
	int g; // number of generations (limited by default)
	int popsize; //population size
	bool limit_evals; // limit evals instead of generations
	long long max_evals; // maximum number of evals before termination (only active if limit_evals is true)

	// Generation Settings 
	int sel; // 1: tournament 2: deterministic crowding 3: lexicase selection 4: age-fitness pareto algorithm
	int tourn_size;
	float rt_rep; //probability of reproduction
	float rt_cross; 
	float rt_mut;
	vector<float> rep_wheel;
	float cross_ar; //crossover alternation rate
	float mut_ar;
	int cross; // 1: ultra; 2: one point
	int mutate; // 1: uniform point mutation; 2: sub-tree mutation 
	bool align_dev; // on or off - adds alignment deviation via gaussian random variable
	bool elitism; // on or off - if on, saves best individual each generation
	//float stoperror; // stop condition / convergence condition
	
	// Data settings
	bool init_validate_on; // initial fitness validation of individuals
	bool train; // choice to turn on training for splitting up the data set
	float train_pct; // percent of data to use for training (validation pct = 1-train_pct)
	bool shuffle_data; // shuffle the data

	bool pop_restart; // restart from previous population
	string pop_restart_path; // restart population file path 
	// Results
	string resultspath;
	//bool loud;

	// use fitness estimator coevolution
	bool EstimateFitness; 
	int FE_pop_size;
	int FE_ind_size;
	int FE_train_size;
	int FE_train_gens;
	bool FE_rank;

	bool estimate_generality;
	int G_sel; 
	bool G_shuffle;


	
	// Problem information
	vector <string> intvars; // internal variables
	//vector <string> extvars; // external variables (external forces)
	vector <string> cons;
	vector <float> cvals;
	vector <string> seeds;
	vector <vector <node> > seedstacks;
	bool AR; //auto-regressive output
	int AR_n; //order of auto-regressive output
	bool AR_lookahead; 
	vector <string> allvars;// = intvars.insert(intvars.end(), extvars.begin(), extvars.end());
	vector <string> allblocks;// = allvars.insert(allvars.end(),consvals.begin(),convals.end());
	//allblocks = allblocks.insert(allblocks.end(),seeds.begin(),seeds.end());

	bool init_trees; // do tree-like recursive initialization of genotypes

	bool ERC; // ephemeral random constants
	bool ERCints;
	int maxERC;
	int minERC;
	int numERC;

	//vector <float> target;
	// Fitness Settings 

	int fit_type; // 1: error, 2: corr, 3: combo
	bool norm_error; // normalize fitness by the standard deviation of the target output
	float max_fit;
	float min_fit;
	bool weight_error; // weight error vector by predefined weights from data file
	vector<float> error_weight; // vector of error weights
	// Operator Settings
	vector <string> op_list;
	vector <int> op_choice; // map op list to pointer location in makeline() pointer function
	vector <float> op_weight;
	vector <int> op_arity; // arity list corresponding op_choice for recursive makeline
	vector <char> return_type; // return type of operators for tree construction
	bool weight_ops_on;

	int min_len;
	int max_len;
	int max_len_init; // initial max len

	int complex_measure; // currently not used


	// Hill Climbing Settings

		// generic line hill climber (Bongard)
	bool lineHC_on;
	int lineHC_its;

		// parameter Hill Climber
	bool pHC_on;
	bool pHC_delay_on;
	float pHC_delay;
	int pHC_size;
	int pHC_its;
	float pHC_gauss;

		// epigenetic Hill Climber
	bool eHC_on;
	int eHC_its;
	float eHC_prob;
	float eHC_init;
	bool eHC_mut; // epigenetic mutation rather than hill climbing
	bool eHC_slim; // use SlimFitness
	// Pareto settings

	bool prto_arch_on;
	int prto_arch_size;
	bool prto_sel_on;

	//island model
	bool islands;
	int island_gens;
	int nt; // number of threads

	//int seed;

	// lexicase selection
	float lexpool; // percent of population to use for lexicase selection events
	bool lexage;// currently not used; up for deletion
	bool lex_class; // use class-based fitness objectives instead of raw error 
	vector<string> lex_metacases; // extra cases to be used for lexicase selection
	bool lex_eps_error; // errors within episilon of the best error are pass, otherwise fail
	bool lex_eps_target; // errors within epsilon of the target are pass, otherwise fail
	bool lex_eps_std; // errors in a standard dev of the best are pass, otherwise fail
	float lex_epsilon; 

	// ==== Printing Options === //
	//print initial population
	bool print_init_pop; 
	bool print_homology; // print homology
	bool print_log; //print log
	bool print_every_pop; // print pop every generation
	bool print_genome; // print genome every generation
	bool print_epigenome; // print epigenome every generation
	bool print_novelty; // print novelty in data file

	// number of log points to print (with eval limitation)
	int num_log_pts;
	//pareto survival setting
	int PS_sel;
	
	//classification
	bool classification;
	bool class_bool; // use binary or multiclass
	bool class_m3gp; // use m3gp fitness
	bool class_prune; // prune dimensions of best individual each generation
	int number_of_classes; // number of unique classes

	// stop condition
	bool stop_condition; 
	//print protected operators
	bool print_protected_operators;
	// number of threads
	params(){ //default values
		g=100; // number of generations (limited by default)
		popsize=500; //population size
		limit_evals=false; // limit evals instead of generations
		max_evals=0; // maximum number of evals before termination (only active if limit_evals is true)
		init_trees=1;
		// Generation Settings 
		sel=1; // 1: tournament 2: deterministic crowding 3: lexicase selection 4: age-fitness pareto algorithm
		tourn_size=2;
		rt_rep=0; //probability of reproduction
		rt_cross=0.8; 
		rt_mut=0.2;
		cross_ar=0.025; //crossover alternation rate
		mut_ar=0.025;
		cross=1; // 1: ultra; 2: one point; 3: subtree
		mutate=1; // 1: one point; 2: subtree
		align_dev = 0;
		elitism = 0;

		// ===============   Data settings
		init_validate_on=0; // initial fitness validation of individuals
		train=0; // choice to turn on training for splitting up the data set
		train_pct=0.5; // default split of data is 50/50
		shuffle_data=0; // shuffle the data
		pop_restart = 0; // restart from previous population
		pop_restart_path=""; // restart population file path
		AR = 0;
		AR_n = 1;
		AR_lookahead = 0;
		// ================ Results and printing
		resultspath="";
		//print every population
		print_every_pop=0;
		//print initial population
		print_init_pop = 0;
		// print homology 
		print_homology = 0;
		//print log
		print_log = 1;
		// number of log points to print (with eval limitation)
		num_log_pts = 0;
		// print csv files of genome each print cycle
		print_genome = 0;
		// print csv files of epigenome each print cycle
		print_epigenome = 0;
		// print number of unique output vectors
		print_novelty = 0;

		// ============ Fitness settings
		fit_type = 1; // 1: error, 2: corr, 3: combo
		norm_error = 0 ; // normalize fitness by the standard deviation of the target output
		weight_error = 0; // weight error vector by predefined weights from data file
		max_fit = 1.0E20;
		min_fit = 0.00000000000000000001;
		
		// Fitness estimators
		EstimateFitness=0; 
		FE_pop_size=0;
		FE_ind_size=0;
		FE_train_size=0;
		FE_train_gens=0;
		FE_rank=0;
		estimate_generality=0;
		G_sel=1; 
		G_shuffle=0;
		// Computer Settings
		//bool parallel;
		//int numcores;

		
		// =========== Program Settings
		ERC = 1; // ephemeral random constants
		ERCints =0 ;
		maxERC = 1;
		minERC = -1;
		numERC = 1;

		min_len = 3;
		max_len = 20;
		max_len_init = 0;

		complex_measure=1; // 1: genotype size; 2: symbolic size; 3: effective genotype size

		weight_ops_on=0;

		// Hill Climbing Settings

			// generic line hill climber (Bongard)
		lineHC_on = 0;
		lineHC_its = 0;

			// parameter Hill Climber
		pHC_on = 0;
		//pHC_size;
		pHC_its = 1;
		pHC_gauss = 0;

			// epigenetic Hill Climber
		eHC_on = 0;
		eHC_its = 1;
		eHC_prob = 0.1;
		eHC_init = 0.5;
		eHC_mut = 0; // epigenetic mutation rather than hill climbing
		eHC_slim = 0; // use SlimFitness
		
		// Pareto settings

		prto_arch_on = 0;
		prto_arch_size = 20;
		prto_sel_on = 0;

		//island model
		islands = 0;
		island_gens = 100;
		nt = 1;

		//seed = 0;

		// lexicase selection
		lexpool = 1;
		lexage = 0;
		lex_class = 0;
		lex_eps_error = false; // errors within episilon of the best error are pass, otherwise fail
		lex_eps_target = false; // errors within epsilon of the target are pass, otherwise fail
		lex_eps_std = false; // errors in a standard dev of the best are pass, otherwise fail
		lex_epsilon = 0.1;
		//pareto survival setting
		PS_sel=1;

		// classification
		classification = 0;
		class_bool = 0;
		class_m3gp = 0;
		class_prune = 0;
		//class_multiclass=0; // use multiclass 
		number_of_classes=1; //for use with multiclass

		stop_condition=1;

		print_protected_operators = 0;
	}
	~params(){}

	
	//void clear()
	//{
	//	
	//	rep_wheel.clear();
	//
	//	// Problem information
	//	intvars.clear(); // internal variables
	//	extvars.clear(); // external variables (external forces)
	//	cons.clear();
	//	cvals.clear();
	//	seeds.clear();
	//
	//	allvars.clear();// = intvars.insert(intvars.end(), extvars.begin(), extvars.end());
	//	allblocks.clear();// = allvars.insert(allvars.end(),consvals.begin(),convals.end());
	//	//allblocks = allblocks.insert(allblocks.end(),seeds.begin(),seeds.end());

	//	op_list.clear();
	//	op_choice.clear(); // map op list to pointer location in makeline() pointer function
	//	op_weight.clear();
	//
	//}

};
#endif
