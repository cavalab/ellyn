#!/bin/bash

# load modules
module load intel/13.1.1.163
module load boost/1.51.0

# clean build directory
rm *.o
# run MakefileTACC
make -f MakefileTACC

