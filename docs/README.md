Documentation
===

ellyn is fast because it uses a c++ library to do most of the computation. However, once you have it installed, you can use it just like you would any other scikit-learn estimator, which makes it easy to do cross validation, ensemble learning, or to build any other kind of ML pipeline design. Follow the installation guide to get it up and running. 

Installation
============
These instructions are written for an [anaconda3](https://www.continuum.io/downloads) default python installation, but you can easily modify the paths to point to your installation. 

```bash

git clone http://github.com/EpistasisLab/ellyn

cd ellyn

conda env create environment.yml

conda activate ellyn-env

python setup.py install
```

`environment.yml` lists the package dependencies for ellyn, if you'd like to install them yourself.

Usage
===
In a python script, import ellyn:

```python
from ellyn import ellyn
```

ellyn uses the same nomenclature as [sklearn](http://scikit-learn.org/) supervised learning modules. You can initialize a few learner in python as:

```python
learner = ellyn()
```

or specify the generations, population size and selection algorithm as:

```python
learner = ellyn(g = 100, popsize = 25, selection = 'lexicase')
```

Given a set of data with variables X and target Y, fit ellyn using the ```fit()``` method:

```python
learner.fit(X,Y)
```

You have now learned a model for your data. Predict your model's response on a new set of variables as

```python
y_pred = learner.predict(X_unseen)
```

Call ellyn from the terminal as

```bash
python -m ellyn.ellyn data_file_name -g 100 -p 50 -sel lexicase
```

try `python -m ellyn.ellyn --help` to see options.

GP options
=====
ellyn uses a stack-based, syntax-free, linear genome for constructing candidate equations. 

*Selection/Survival options*

 - tournament
 - deterministic crowding
 - lexicase selection
 - age-fitness pareto optimization 
 - SPEA2
 - random

*Variation*

 - subtree, uniform, and point mutation
 - subtree, unifrom, and point crossover

*Parameter learning*

 - stochastic hill climbing

Cite
====

ellyn has been used in several publications. 
Cite the one that best represents your use case, or you can cite my dissertation if you're not sure. 

2018

 - Orzechowski, P., La Cava, W., & Moore, J. H. (2018). 
Where are we now? A large benchmark study of recent symbolic regression methods. 
GECCO 2018. [DOI](https://doi.org/10.1145/3205455.3205539), [Preprint](https://www.researchgate.net/profile/Patryk_Orzechowski/publication/324769381_Where_are_we_now_A_large_benchmark_study_of_recent_symbolic_regression_methods/links/5ae779b70f7e9b837d392dc9/Where-are-we-now-A-large-benchmark-study-of-recent-symbolic-regression-methods.pdf)

2017 

 - La Cava, W., Silva, S., Vanneschi, L., Spector, L., and Moore, J. (2017). "Genetic programming representations for multi-dimensional feature learning in biomedical classification." Evo Applications, EvoStar 2017, Amsterdam, Netherlands. [preprint](http://williamlacava.com/pubs/evobio_m4gp_lacava.pdf)

2016

 - La Cava, William G., "Automatic Development and Adaptation of Concise Nonlinear Models for System Identification" (2016). Doctoral Dissertations May 2014 - current. 731. [link](http://scholarworks.umass.edu/dissertations_2/731/)

 - La Cava, W., Danai, K., Spector, L., (2016). "Inference of Compact Nonlinear Dynamic Models by Epigenetic Local Search." Engineering Applications of Artificial Intelligence. [doi:10.1016/j.engappai.2016.07.004 ](http://authors.elsevier.com/a/1TVk33OWJ8hFJk)

 -  La Cava, W., Spector, L., Danai, K. (2016). "epsilon-Lexicase selection for regression." Proceedings of the Genetic and Evolutionary Computation Conference (GECCO). ACM, Denver, CO. [preprint](http://williamlacava.com/pubs/GECCO_lex_reg-corrected.pdf)

2015 

 - La Cava, W., Danai, K., Spector, L., Fleming, P., Wright, A., Lackner, M. (2015). "Automatic identification of wind turbine models using evolutionary multi-objective optimization." Renewable Energy. [doi:10.1016/j.renene.2015.09.068](http://www.sciencedirect.com/science/article/pii/S0960148115303475)
