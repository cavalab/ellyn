Develep
=======
Develep is a genetic programming tool for symbolic regression that incorporates epigenetic learning.
It is very much under development! Or you could say, under-develeped...

The files are built in Visual Studio C++ Express 2010 with some boost library dependencies. 

About
=====
Develep uses a developmental, syntax-free, linear genome for constructing candidate equations. 

It is built to include different evolutionary methods for system identification adapted from literature. The options include  normal tournament selection, deterministic crowding, and age-pareto fitness selection. All algorithm choices are mangaged by one parameter file. 

How to Build
============
I've built the project in Visual Studio 2010 professional as well as C++ Express (which is free from Microsoft), but hope to have a build with the GNU compiler in the future. The boost libraries (also free) need to be installed. If you use VS 2010 Express, the OpenMP files (which were removed from VS 2010) need to be added to the VS path. 

How to Run Develep
==================
Run Develep from the command line. Here is the syntax:
```
Develep sampleparams.txt sampledata.txt
```
As you can see, Develep takes two arguments: a parameter file and a data file. The parameter file includes all of the run-time settings for your search. The data file includes all your experimental data. See the sampleparams.txt and sampledata.txt filesto see how formatting works.   

How to run RunTrials
====================
RunTrials will run Develep for many trials. Here is the syntax:
```
RunTrials sampletrials.txt
```
RunTrials takes one input file. The trials input file contains three columns:
```
#trials parameterfile datafile
```
These are the simple instructions for running Develep. 

Develep
Copyright (C) 2014  William La Cava


This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License (License.txt) for more details.


 

