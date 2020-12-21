from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import subprocess
import os
import sys
from glob import glob

# SOURCES = sorted(list(glob('ellyn/ellen/*.cpp')))
SOURCES = [
           'stdafx.cpp',
           'AgeBreed.cpp',
           'AgeFitSurvival.cpp',
           'AgeFitGenSurvival.cpp',
           'ApplyGenetics.cpp',
           'BruteforceAlgo.cpp',
           'Crossover.cpp',
           'Datapoint.cpp',
           'DC.cpp',
           'EpiHC.cpp',
           'EpiMut.cpp',
           'Eqn2Line.cpp',
           'Fitness.cpp',
           'FitnessEstimator.cpp',
           'general_fns.cpp',
           'Generation.cpp',
           'HillClimb.cpp',
           'InitPop.cpp',
           'LexicaseSelect.cpp',
           'Line2Eqn.cpp',
           'Mutate.cpp',
           'NondominatedsortAlgo.cpp',
           'pareto.cpp', 
           'pareto_fc.cpp', 
           'ParetoSurvival.cpp', 
           'Prune.cpp', 
           'runEllenGP.cpp',
           'StablesortAlgo.cpp', 
           'StochasticGradient.cpp', 
           'strdist.cpp',
           'Tournament.cpp']
SOURCES = ['ellyn/ellen/'+s for s in SOURCES]

CONDA_PATH = os.environ['CONDA_PREFIX']

INCLUDE = [CONDA_PATH + '/include/eigen3',
           CONDA_PATH + '/include/', 
           # '/usr/local/lib/python3.5/dist-packages/numpy/core/include/'
          ]
LIB = [CONDA_PATH + '/lib','-lpython3']
COMPILE_ARGS = ['-std=c++0x','-fopenmp', 
                '-Wno-sign-compare',
                '-Wno-unused-variable',
                '-Wno-unused-value',
               ]
LINK_ARGS = ['-lboost_python39','-lpython3','-fopenmp',
             "-Wl,-rpath,'" + CONDA_PATH +"/lib/'"]
       
# LINK_ARGS=[]

#avoid a gcc warning below:
# cc1plus: warning: command line option ‘-Wstrict-prototypes’ is valid
# for C/ObjC but not for C++
class BuildExt(build_ext):
    def build_extensions(self):
        if '-Wstrict-prototypes' in self.compiler.compiler_so:
            self.compiler.compiler_so.remove('-Wstrict-prototypes')
        super().build_extensions()

setup(name='ellyn',
      version='0.1',
      description='A Python-wrapped Genetic Programming System',
      author='William La Cava',
      author_email='williamlacava@gmail.com',
      url='https://github.com/EpistasisLab/ellyn',
      py_modules = ['ellyn.ellyn'],
      ext_modules = [Extension('ellyn.elgp',
                              SOURCES,
                              include_dirs = INCLUDE,
                              library_dirs = LIB,
                              extra_compile_args = COMPILE_ARGS,
                              extra_link_args = LINK_ARGS 
                             )
                    ],
    cmdclass={'build_ext': BuildExt}
)


