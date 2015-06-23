ellenGP
=======

ellenGP is a genetic programming tool for symbolic regression and multi-class classification that incorporates epigenetic learning and uses a stack-based, linear representation.
It is very much under development!

There are some boost library dependencies, including regex. 

The files have been built in Visual Studio C++ 2010 and in linux with g++ and the intel c++ compiler. 

About
=====
ellenGP uses a stack-based, syntax-free, linear genome for constructing candidate equations. 

It is built to include different evolutionary methods for system identification adapted from literature. The options include  normal tournament selection, deterministic crowding, and age-pareto fitness selection. All algorithm choices are mangaged by one parameter file. 

How to Build
============
I've built the project in Visual Studio 2010 professional as well as C++ Express (which is free from Microsoft), and in linux with g++ and the intel c++ compiler using the make files. If you use VS 2010 Express, the OpenMP files (which were removed from VS 2010) need to be added to the VS path. 

There are a couple of external library dependencies: 

- [boost libraries](http://www.boost.org) - a set of multi-purpose c++ libraries
- [eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) - a c++ template library for linear algebra

How to run ellenGP
==================
Run ellenGP like this:
```
ellenGP sampleparams.txt sampledata.txt
```
As you can see, ellenGP takes two arguments: a parameter file and a data file. The parameter file includes all of the run-time settings for your search. The data file includes all your experimental data. See the sampleparams.txt and sampledata.txt files to see how formatting works.   

How to run RunTrials
====================
RunTrials will run ellenGP for many trials. It uses OpenMP to parallelize the trials. Here is the syntax:
```
RunTrials sampletrials.txt
```
RunTrials takes one input file (sampletrials.txt). The trials input file contains three columns:
```
[#trials] [parameterfile] [datafile]
```
for example,
...
100 ../in/sampleparams.txt ../in/sampledata.txt

These are the simple instructions for running ellenGP. 

RunTrialsMPI
====================
RunTrialsMPI is the same as RunTrials except it is written to be compiled on the clusters (the TACC cluster Stampede as well as the Umass HPCC cluster). MakefileTACC and MakefileUMG has the compilation notes. It has been built using intel icpc and the MPI compiler mvapich2 from OSU, as well as g++ with mpicxx.

Settings
===========
Here is a comprehensive list of all of the options that you can include in the parameter file. 

Setting | Default | Description
------- | ------- | -----------
g	|	100	|	 number of generations
popsize	|	500	|	population size
limit_evals	|	0	|	 limit evals instead of number of generations
max_evals	|	0	|	 max point evaluations
	|		|	
	|		|	
 Generation Settings 	|		|	
sel	|	1	|	 1: tournament 2: deterministic crowding 3: lexicase selection 4: age-fitness pareto algorithm
lexage	|	0	|	 use age survival with lexicase selection
PS_sel	|	1	|	objectives for pareto survival. 1: age + fitness; 2: age+fitness+generality; 3: age+fitness+complexity; 4: class fitnesses (classification ONLY); 5: class fitnesses+ age (classification ONLY)
tourn_size	|	2	|	number  of individuals in each tournament
rt_rep	|	0	|	rate of reproduction
rt_cross	|	0.8	|	rate of crossover
rt_mut	|	0.2	|	rate of mutation
cross	|	2	|	 1: ultra 2: one point1 3: sub-tree
cross_ar	|	0.025	|	 crossover alternation rate (ultra only)
mut_ar	|	0.025	|	mutation alternation rate 
align_dev	|	0	|	 on or off; adds gaussian alignment deviation to crossover 
elitism	|	0	|	 save best individual each generation
	|		|	
init_validate_on	|	0	|	 initial fitness validation of starting population
train	|	0	|	 split data into training and validation sets
train_pct	|	0.5	|	  percent of data to be used in training
shuffle_data	|	0	|	 shuffle the data 
pop_restart	|	0	|	restart run from previous population specified by pop_restart_path
pop_restart_path	|	""	|	filename of restart population with path
	|		|	
 Results and Printing Options	|		|	
resultspath	|	""	|	path where results are saved
print_every_pop	|	0	|	  save printout of population at every generation
print_genome	|	0	|	prints genome for visualization in paraview
num_log_pts	|	0	|	number of log points to print (0 means print each generation)
	|		|	
classification	|	0	|	defines a classification, rather than regression, problem
class_bool	|	0	|	interpret class labels as bit-string conversion of boolean stack output
class_m3gp	|	0	|	use mahalanobis distance classification fitness
 Problem information	|		|	
intvars	|	none	|	variables in data file to use in programs
cons     	|	none	|	constants to use
cvals	|	none	|	values of constants defined in cons
seeds	|	none	|	seed partial solutions
AR	|	0	|	include auto-regressive output variables
AR_n 	|	1	|	  order of auto-regression (number of time-steps back)
AR_lookahead	|	0	|	 just predict one output ahead
	|		|	
ERC	|	1	|	 ephemeral random constants
ERCints	|	0	|	
maxERC	|	1	|	
minERC	|	-1	|	
numERC	|	1	|	
	|		|	
fitness settings	|		|	
fit_type	|	1	|	 1: mean absolute error, 2: corr, 3: combo, 4: VAF
norm_error	|	0	|	  normalize error by the standard deviation of the target data being u
max_fit	|	1.00E+20	|	maximum fitness possible
min_fit	|	1.00E-20	|	minimum fitness possible
estimate_fitness 	|	0	|	 coevolve fitness estimators
FE_pop_size	|	0	|	 fitness estimator population size
FE_ind_size	|	0	|	 number of fitness cases for FE to use
FE_train_size	|	0	|	 trainer population size
FE_train_gens	|	0	|	number of generations between trainer selections
FE_rank	|	0	|	 use rank for FE fitness rather than error  
estimate_generality	|	0	|	estimate how well the solutions generalize using the validation portion of the fitness estimator
G_sel	|	0	|	which fit_type to use to test generality
G_shuffle	|	0	|	shuffles data each generation
	|		|	
op_list	|	n,v,+,-,*,/	|	available operators: n
op_weight	|	empty	|	4
weight_ops_on	|	0	|	weight the operators differently
min_len	|	3	|	minimum program length
max_len	|	20	|	maximum length a program is allowed to be
max_len_init	|	max_len	|	option to specify different max length for initial population
init_trees	|		|	 initialize genotypes as syntactically valid trees rather than randomized stacks
	|		|	
complex_measure	|	2	|	 1: genotype size 2: symbolic size 3: effective genotype size
	|		|	
	|		|	
 Hill Climbing Settings	|		|	
	|		|	
	|		|	
 parameter Hill Climber	|		|	
pHC_on	|	0	|	parameter hill climbing each generation
pHC_its	|	1	|	number of iterations
	|		|	
 epigenetic Hill Climber	|		|	
eHC_on	|	0	|	epigenetic hill climbing
eHC_its	|	1	|	number of iterations
eHC_prob	|	0.1	|	probability of a gene being switched
eHC_init	|	0.5	|	 percent of expressed genes in initial genotypes
eHC_slim	|	0	|	  minimize point evaluations as much as possible
eHC_mut	|	0	|	do mutation rather than hill climbing
	|		|	
 Pareto settings	|		|	
prto_arch_on	|	0	|	
prto_arch_size	|	20	|	
	|		|	
 island model	|		|	
islands	|	0	|	use multiple island populations, one for each core. 
island_gens	|	100	|	number of generations until island populations are shuffled

 
FYI
===
ellenGP
Copyright (C) 2014  William La Cava


This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License (License.txt) for more details.


 

