ellyn
=======

ellyn is a Python-wrapped version of [ellenGP](http://www.github.com/lacava/ellen) that allows ellenGP to play nice with sci-kitlearn. ellyn's parameter settings are totally accessible from the commandline, whereas ellenGP relies on a parameter file. This can make batch jobs less tedious. 

ellyn is a genetic programming tool for symbolic regression and multi-class classification that incorporates epigenetic learning and uses a stack-based, linear representation.

ellyn also inherits the `BaseEstimator` class used by [sklearn's](http://scikit-learn.org/) supervised learning modules. That means that ellyn can be used with that environment of tools, including the hyperparameter optimization tools and cross validation tools.

ellyn is also very fast due to its c++ underpinning. As a consequence, there are two library dependencies: [Boost](http://www.boost.org) and [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page), that need to be installed in order to make the ellyn library for python usage.

GP options
=====
ellyn uses a stack-based, syntax-free, linear genome for constructing candidate equations. 

It is built to include different evolutionary methods for system identification adapted from literature. The options include  normal tournament selection, deterministic crowding, lexicase selection, and pareto optimization. Variation operators include tree-, uniform- and point-based mutation and crossover operators. There is a built-in stochastic hill climber for parameter learning as well. 


How to Build
============
After downloading the dependencies, go to the download location in terminal. Then type

```bash
cd ellyn/ellen
make
```

In addition to downloading those packages, the paths to them need to be modified in the Makefiles for ellyn and RunTrials. 

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
learner = ellyn(generations = 100, popsize = 25, selection = 'lexicase')
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

FYI
===
ellyn
Copyright (C) 2016  William La Cava


This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License (License.txt) for more details.


 

