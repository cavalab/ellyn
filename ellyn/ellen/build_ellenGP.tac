#!/bin/bash

# load modules
module load intel/14.0.1.106
module load boost/1.51.0

# clean build directory
rm *.o
# run MakefileTACC
make -f MakefileTACC

