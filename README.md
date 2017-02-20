ellyn
=======

ellyn is a Python-wrapped version of [ellenGP](http://www.github.com/lacava/ellen) that allows ellenGP to play nice with sci-kitlearn. ellyn's parameter settings are totally accessible from the commandline, whereas ellenGP relies on a parameter file. This can make batch jobs less tedious. 

ellyn is a genetic programming tool for symbolic regression and multi-class classification that incorporates epigenetic learning and uses a stack-based, linear representation.

ellyn also inherits the `BaseEstimator` class used by [sklearn's](http://scikit-learn.org/) supervised learning modules. That means that ellyn can be used with that environment of tools, including the hyperparameter optimization tools and cross validation tools.

ellyn is also very fast due to its c++ underpinning. As a consequence, there are two library dependencies: [Boost](http://www.boost.org) and [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page), that need to be installed in order to make the ellyn library for python usage.

Head over to the [docs](http://www.github.com/lacava/ellyn/docs) page for more info. 

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


 

