ellenGP
=======

ellenGP is a genetic programming tool for symbolic regression that incorporates epigenetic learning and uses a stack-based, linear representation.
It is very much under development!

There are some boost library dependencies, including regex. 

The files have been built in Visual Studio C++ 2010 and in linux with g++ and the intel c++ compiler. 

About
=====
ellenGP uses a stack-based, syntax-free, linear genome for constructing candidate equations. 

It is built to include different evolutionary methods for system identification adapted from literature. The options include  normal tournament selection, deterministic crowding, and age-pareto fitness selection. All algorithm choices are mangaged by one parameter file. 

How to Build
============
I've built the project in Visual Studio 2010 professional as well as C++ Express (which is free from Microsoft), and in linux with g++ and the intel c++ compiler using the make files. The boost libraries (also free) need to be installed. If you use VS 2010 Express, the OpenMP files (which were removed from VS 2010) need to be added to the VS path. 

How to Run ellenGP
==================
Run ellenGP from the command line. Here is the syntax:
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
RunTrialsMPI is the same as RunTrials except it is written to be compiled on the TACC cluster Stampede. MakefileTACC has the compilation notes. It has been built using intel icpc and the MPI compiler mvapich2 from OSU. 
 
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


 

