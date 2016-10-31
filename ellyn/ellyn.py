# -*- coding: utf-8 -*-
"""
Copyright 2016 William La Cava

license: GNU/GPLv3

"""

import argparse
# from ._version import __version__

from sklearn.base import BaseEstimator
from sklearn.cross_validation import train_test_split
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error, explained_variance_score, accuracy_score
import numpy as np
import pandas as pd
import warnings
import copy
import itertools as it
import pdb
import ellen.lib.elgp as elgp
from update_checker import update_check
from DistanceClassifier import DistanceClassifier
# from numpy.ctypeslib import ndpointer
# import ctypes
# from joblib import Parallel, delayed

# import multiprocessing as mp
# NUM_THREADS = mp.cpu_count()


class ellyn(BaseEstimator):
    """ellyn uses GP to build its models.

    It uses a stack-based, syntax-free, linear genome for constructing candidate equations.
    It is built to include different evolutionary methods for system identification
    adapted from literature. The options include normal tournament selection,
    deterministic crowding, (epsilon) lexicase selection, and age-pareto fitness selection.

    All algorithm choices are accessible from the command line.

    """
    update_checked = False

    def __init__(self, g=100, popsize=500, limit_evals=False, max_evals = 0,
                 selection='tournament',sel=1,classification=False,islands=False,
                 fit_type=None,verbosity=0,random_state=0,
                 class_m4gp=False,scoring_function=mean_squared_error,print_log=False,
                 class_bool=False,max_len=None,island_gens=50, numERC=None,
                 print_every_pop=None, op_weights=None,
                 FE_pop_size=None, rt_cross=None, ERC_ints=None,
                 train_pct=None, eHC_on=None, PS_sel=None,
                 lex_eps_global=None, lex_pool=None, AR_nka=None,
                 print_homology=None, max_len_init=None, prto_arch_size=None,
                 cvals=None, stop_condition=None, lex_metacases=None,
                 FE_rank=None, EHC_its=None, lex_eps_error_mad=True,
                 ERC=None, erc_ints=None, AR_na=None, rt_mut=None,
                 pop_restart=None, seeds=None, tourn_size=None, prto_arch_on=None,
                 FE_ind_size=None, lex_eps_target_mad=None, maxERC=None,
                 resultspath=None, AR=None, rt_rep=None, estimate_fitness=None,
                 pHC_on=None, INPUT_FILE=None, FE_train_size=None,
                 DISABLE_UPDATE_CHECK=False, AR_lookahead=None, VERBOSITY=0, pop_restart_path=None,
                 INPUT_SEPARATOR=None,
                 AR_nkb=None, num_log_pts=None, train=None,
                 FE_train_gens=None, AR_nb=None, init_trees=None,
                 print_novelty=None, EHC_slim=None, elitism=None,
                 print_genome=None, pHC_its=None, shuffle_data=None, class_prune=None,
                 EHC_init=None, init_validate_on=None, minERC=None,
                 op_list=None, EHC_mut=None, min_len=None):
                # sets up GP.

        # Save params to be recalled later by get_params()
        self.params = locals()  # Must be placed before any local variable definitions
        self.params.pop('self')
        # convert selection method to number
        try:
            self.params['sel']={'tournament': 1,'dc':2,'lexicase': 3,'afp': 4,'rand': 5}[selection]
        except Exception:
            print("sel:",selection)
            raise(Exception)


        if fit_type:
            self.params['scoring_function'] = {'MAE':mean_absolute_error,'MSE':mean_squared_error,
                                    'R2':r2_score,'VAF':explained_variance_score,
                                    'combo':mean_absolute_error}[fit_type]
        elif classification:
            # self.params['fit_type'] = 'MAE'
            self.params['scoring_function'] = accuracy_score
        else:
            self.params['scoring_function'] = mean_squared_error
        #
        # # set verbosity
        if not print_log:
            self.params['print_log'] = self.params['verbosity'] > 1

        self._best_estimator = []
        #convert m4gp argument to m3gp used in ellenGP
        if classification and not (self.params['class_m4gp'] or self.params['class_bool']):
            #default to M4GP in the case no classifier specified
            self.params['class_m4gp'] = True

        # # delete items that have a default value of None so that the elgp defaults will be used
        for k in list(self.params.keys()):
            if self.params[k] is None:
                del self.params[k]

    def fit(self, features, labels):
        """Fit model to data"""

        # pdb.set_trace()
        np.random.seed(self.params['random_state'])
        # run ellenGP
        elgp.runEllenGP(self.params,np.asarray(features,dtype=np.float32,order='C'),
                        np.asarray(labels,dtype=np.float32,order='C'),self._best_estimator)
        # print("best program:",self._best_estimator)

        # if M4GP is used, call Distance Classifier
        # pdb.set_trace()
        if self.params['class_m4gp']:
            if self.params['verbosity'] > 0: print("Storing DistanceClassifier...")
            self.DC = DistanceClassifier()
            self.DC.fit(self._out(self._best_estimator,features),labels)

        ####
        # initial model


    def predict(self, testing_features):
        """predict on a holdout data set."""
        # print("best_inds:",self._best_inds)
        # print("best estimator size:",self._best_estimator.coef_.shape)
        # tmp = self._out(self._best_estimator,testing_features)
        if 'class_m4gp' in self.params:
            if self.params['class_m4gp']:
                return self.DC.predict(self._out(self._best_estimator,testing_features))
        else:
            return self._out(self._best_estimator,testing_features)

    def fit_predict(self, features, labels):
        """Convenience function that fits a pipeline then predicts on the provided features

        Parameters
        ----------
        features: array-like {n_samples, n_features}
            Feature matrix
        labels: array-like {n_samples}
            List of class labels for prediction

        Returns
        ----------
        array-like: {n_samples}
            Predicted labels for the provided features

        """
        self.fit(features, labels)
        return self.predict(features)

    def score(self, testing_features, testing_labels):
        """estimates accuracy on testing set"""
        # print("test features shape:",testing_features.shape)
        # print("testing labels shape:",testing_labels.shape)
        yhat = self.predict(testing_features)
        # pdb.set_trace()
        return self.params['scoring_function'](testing_labels,yhat)

    def export(self, output_file_name):
        """does nothing currently"""

    def get_params(self, deep=None):
        """Get parameters for this estimator

        This function is necessary for ellyn to work as a drop-in feature constructor in,
        e.g., sklearn.cross_validation.cross_val_score

        Parameters
        ----------
        deep: unused
            Only implemented to maintain interface for sklearn

        Returns
        -------
        params: mapping of string to any
            Parameter names mapped to their values
        """
        return self.params

    def _eval(self,n, features, stack_float, stack_bool):
        """evaluation function for best estimator"""
        np.seterr(all='ignore')
        if len(stack_float) >= n[1]:
            stack_float.append(eval_dict[n[0]](n,features,stack_float,stack_bool))
            if any(np.isnan(stack_float[-1])) or any(np.isinf(stack_float[-1])):
                print("problem operator:",n)

    def _out(self,I,features):
        """computes the output for individual I"""
        stack_float = []
        stack_bool = []
        # print("stack:",I.stack)
        # evaulate stack over rows of features,labels
        for n in I:
            self._eval(n,features,stack_float,stack_bool)
            # print("stack_float:",stack_float)
        if 'class_m4gp' in self.params:
            if self.params['class_m4gp']:
                return np.asarray(stack_float).transpose()
        else:
            return stack_float[-1]

eval_dict = {
# float operations
    '+': lambda n,features,stack_float,stack_bool: stack_float.pop() + stack_float.pop(),
    '-': lambda n,features,stack_float,stack_bool: stack_float.pop() - stack_float.pop(),
    '*': lambda n,features,stack_float,stack_bool: stack_float.pop() * stack_float.pop(),
    '/': lambda n,features,stack_float,stack_bool: divs(stack_float.pop(),stack_float.pop()),
    'sin': lambda n,features,stack_float,stack_bool: np.sin(stack_float.pop()),
    'cos': lambda n,features,stack_float,stack_bool: np.cos(stack_float.pop()),
    'exp': lambda n,features,stack_float,stack_bool: np.exp(stack_float.pop()),
    'log': lambda n,features,stack_float,stack_bool: logs(stack_float.pop()),#np.log(np.abs(stack_float.pop())),
    'x':  lambda n,features,stack_float,stack_bool: features[:,n[2]],
    'k': lambda n,features,stack_float,stack_bool: np.ones(features.shape[0])*n[2],
    '^2': lambda n,features,stack_float,stack_bool: stack_float.pop()**2,
    '^3': lambda n,features,stack_float,stack_bool: stack_float.pop()**3,
    'sqrt': lambda n,features,stack_float,stack_bool: np.sqrt(np.abs(stack_float.pop())),
    # 'rbf': lambda n,features,stack_float,stack_bool: np.exp(-(np.norm(stack_float.pop()-stack_float.pop())**2)/2)
# bool operations
    '!': lambda n,features,stack_float,stack_bool: not stack_bool.pop(),
    '&': lambda n,features,stack_float,stack_bool: stack_bool.pop() and stack_bool.pop(),
    '|': lambda n,features,stack_float,stack_bool: stack_bool.pop() or stack_bool.pop(),
    '==': lambda n,features,stack_float,stack_bool: stack_bool.pop() == stack_bool.pop(),
    '>': lambda n,features,stack_float,stack_bool: stack_float.pop() > stack_float.pop(),
    '<': lambda n,features,stack_float,stack_bool: stack_float.pop() < stack_float.pop(),
    '}': lambda n,features,stack_float,stack_bool: stack_float.pop() >= stack_float.pop(),
    '{': lambda n,features,stack_float,stack_bool: stack_float.pop() <= stack_float.pop(),
    # '>_b': lambda n,features,stack_float,stack_bool: stack_bool.pop() > stack_bool.pop(),
    # '<_b': lambda n,features,stack_float,stack_bool: stack_bool.pop() < stack_bool.pop(),
    # '>=_b': lambda n,features,stack_float,stack_bool: stack_bool.pop() >= stack_bool.pop(),
    # '<=_b': lambda n,features,stack_float,stack_bool: stack_bool.pop() <= stack_bool.pop(),
}
def divs(x,y):
    """safe division"""
    # pdb.set_trace()
    tmp = np.ones(y.shape)
    nonzero_y = np.abs(y) >= 0.000001
    # print("nonzero_y.sum:", np.sum(nonzero_y))
    tmp[nonzero_y] = x[nonzero_y]/y[nonzero_y]
    return tmp

def logs(x):
    """safe log"""
    tmp = np.zeros(x.shape)
    nonzero_x = np.abs(x) >= 0.000001
    tmp[nonzero_x] = np.log(np.abs(x[nonzero_x]))
    return tmp

def positive_integer(value):
    """Ensures that the provided value is a positive integer; throws an exception otherwise

    Parameters
    ----------
    value: int
        The number to evaluate

    Returns
    -------
    value: int
        Returns a positive integer
    """
    try:
        value = int(value)
    except Exception:
        raise argparse.ArgumentTypeError('Invalid int value: \'{}\''.format(value))
    if value < 0:
        raise argparse.ArgumentTypeError('Invalid positive int value: \'{}\''.format(value))
    return value

def float_range(value):
    """Ensures that the provided value is a float integer in the range (0., 1.); throws an exception otherwise

    Parameters
    ----------
    value: float
        The number to evaluate

    Returns
    -------
    value: float
        Returns a float in the range (0., 1.)
    """
    try:
        value = float(value)
    except:
        raise argparse.ArgumentTypeError('Invalid float value: \'{}\''.format(value))
    if value < 0.0 or value > 1.0:
        raise argparse.ArgumentTypeError('Invalid float value: \'{}\''.format(value))
    return value

# main functions
def main():
    """Main function that is called when ellyn is run on the command line"""
    parser = argparse.ArgumentParser(description='A genetic programming '
                                                 'system for regression and classification.',
                                     add_help=False)

    parser.add_argument('INPUT_FILE', type=str, help='Data file to run ellyn on; ensure that the target/label column is labeled as "label".')

    parser.add_argument('-h', '--help', action='help', help='Show this help message and exit.')

    parser.add_argument('-is', action='store', dest='INPUT_SEPARATOR', default=None,
                        type=str, help='Character used to separate columns in the input file.')

# GP Runtime Options
    parser.add_argument('-g', action='store', dest='g', default=None,
                        type=positive_integer, help='Number of generations to run ellyn.')

    parser.add_argument('-p', action='store', dest='popsize', default=None,
                        type=positive_integer, help='Number of individuals in the GP population.')

    parser.add_argument('--limit_evals', action='store_true', dest='limit_evals', default=None,
                    help='Limit evaluations instead of generations.')

    parser.add_argument('-me', action='store', dest='max_evals', default=None,
                        type=float_range, help='Max point evaluations.')

    parser.add_argument('-sel', action='store', dest='sel', default='tournament', choices = ['tournament','dc','lexicase','afp','rand'],
                        type=str, help='Selection method (Default: tournament)')

    parser.add_argument('-PS_sel', action='store', dest='PS_sel', default=None, choices = [1,2,3,4,5],
                        type=str, help='objectives for pareto survival. 1: age + fitness; 2: age+fitness+generality; 3: age+fitness+complexity; 4: class fitnesses (classification ONLY); 5: class fitnesses+ age (classification ONLY)')

    parser.add_argument('-tourn_size', action='store', dest='tourn_size', default=None,
                        type=positive_integer, help='Tournament size for tournament selection (Default: 2)')

    parser.add_argument('-mr', action='store', dest='rt_mut', default=None,
                        type=float_range, help='GP mutation rate in the range [0.0, 1.0].')

    parser.add_argument('-xr', action='store', dest='rt_cross', default=None,
                        type=float_range, help='GP crossover rate in the range [0.0, 1.0].')

    parser.add_argument('-rr', action='store', dest='rt_rep', default=None,
                        type=float_range, help='GP reproduction rate in the range [0.0, 1.0].')

    parser.add_argument('--elitism', action='store_true', dest='elitism', default=None,
                help='Flag to force survival of best individual in GP population.')

    parser.add_argument('--init_validate', action='store_true', dest='init_validate_on', default=None,
                help='Flag to guarantee initial population outputs are valid (no NaNs or Infs).')

    parser.add_argument('--no_stop', action='store_false', dest='stop_condition', default=None,
                    help='Flag to keep running even though fitness < 1e-6 has been reached.')
# Data Options
    parser.add_argument('--validate', action='store_true', dest='train', default=None,
                help='Flag to split data into training and validation sets.')

    parser.add_argument('-train_split', action='store', dest='train_pct', default=None, type = float_range,
                help='Fraction of data to use for training in [0,1]. Only used when --validate flag is present.')

    parser.add_argument('--shuffle', action='store_true', dest='shuffle_data', default=None,
                help='Flag to shuffle data samples.')

    parser.add_argument('--restart', action='store_true', dest='pop_restart_path', default=None,
                help='Flag to restart from previous population.')

    parser.add_argument('-pop', action='store_true', dest='pop_restart', default=None,
                help='Population to use in restart. Only used when --restart flag is present.')

    parser.add_argument('--AR', action='store_true', dest='AR', default=None,
                    help='Flag to use auto-regressive variables.')

    parser.add_argument('--AR_lookahead', action='store_true', dest='AR_lookahead', default=None,
                    help='Flag to only estimate on step ahead of data when evaluating candidate models.')

    parser.add_argument('-AR_na', action='store', dest='AR_na', default=None,
                    type=positive_integer, help='Order of auto-regressive output (y).')
    parser.add_argument('-AR_nb', action='store', dest='AR_nb', default=None,
                    type=positive_integer, help='Order of auto-regressive inputs (x).')
    parser.add_argument('-AR_nka', action='store', dest='AR_nka', default=None,
                    type=positive_integer, help='AR state (output) delay.')
    parser.add_argument('-AR_nkb', action='store', dest='AR_nkb', default=None,
                    type=positive_integer, help='AR input delay.')

# Results and Printing Options
    parser.add_argument('-o', action='store', dest='resultspath', default=None,
                        type=str, help='Path where results will be saved.')
    parser.add_argument('--print_log', action='store_true', dest='print_log', default=False,
                    help='Flag to print log to terminal.')

    parser.add_argument('--print_every_pop', action='store_true', dest='print_every_pop', default=None,
                    help='Flag to seed initial GP population with components of the ML model.')

    parser.add_argument('--print_genome', action='store_true', dest='print_genome', default=None,
                    help='Flag to prints genome for visualization.')

    parser.add_argument('--print_novelty', action='store_true', dest='print_novelty', default=None,
                    help='Flag to print number of unique output vectors.')

    parser.add_argument('--print_homology', action='store_true', dest='print_homology', default=None,
                    help='Flag to print genetic homology in programs.')

    parser.add_argument('-num_log_pts', action='store', dest='num_log_pts', default=None,
                        type=positive_integer, help='number of log points to print (0 means print each generation)')

# Classification Options

    parser.add_argument('--class', action='store_true', dest='classification', default=None,
                    help='Flag to define a classification problem instead of regression (the default).')

    parser.add_argument('--class_bool', action='store_true', dest='class_bool', default=None,
                    help='Flag to interpret class labels as bit-string conversion of boolean stack output.')

    parser.add_argument('--class_m4gp', action='store_true', dest='class_m4gp', default=False,
                    help='Flag to use the M4GP algorithm.')

    parser.add_argument('--class_prune', action='store_true', dest='class_prune', default=None,
                        help='Flag to prune the dimensions of the best individual each generation.')

# Terminal and Operator Options
    parser.add_argument('-op_list', nargs = '*', action='store', dest='op_list', default=None,
                    type=str, help='Operator list. Default: +,-,*,/,n,v. available operators: n v + - * / sin cos log exp sqrt = ! < <= > >= if-then if-then-else & |')

    parser.add_argument('-op_weight', nargs = '*',action='store', dest='op_weights', default=None,
                    help='Operator weights for each element in operator list. If not specified all operators have the same weights.')

    parser.add_argument('-constants', nargs = '*', action='store', dest='cvals', default=None,
                        type=float, help='Seed GP initialization with constant values.')

    parser.add_argument('-seeds', nargs = '*', action='store', dest='seeds', default=None,
                        type=str, help='Seed GP initialization with partial solutions, e.g. (x+y). Each partial solution must be enclosed in parentheses.')

    parser.add_argument('-min_len', action='store', dest='min_len', default=None,
                        type=positive_integer, help='Minimum length of GP programs.')

    parser.add_argument('-max_len', action='store', dest='max_len', default=None,
                        type=positive_integer, help='Maximum number of nodes in GP programs.')

    parser.add_argument('-max_len_init', action='store', dest='max_len_init', default=None,
                        type=positive_integer, help='Maximum number of nodes in initialized GP programs.')

    parser.add_argument('--erc_ints', action='store_true', dest='erc_ints', default=None,
                help='Flag to use integer instead of floating point ERCs.')

    parser.add_argument('--no_erc', action='store_false', dest='ERC', default=None,
                help='Flag to turn of ERCs. Useful if you are specifying constant values and don''t want random ones.')

    parser.add_argument('-min_erc', action='store', dest='minERC', default=None,
                    help='Minimum ERC value.')

    parser.add_argument('-max_erc', action='store', dest='maxERC', default=None, type = int,
                    help='Maximum ERC value.')

    parser.add_argument('-num_erc', action='store', dest='numERC', default=None, type = int,
                    help='Number of ERCs to use.')

    parser.add_argument('--trees', action='store_true', dest='init_trees', default=None,
                    help='Flag to initialize genotypes as syntactically valid trees rather than randomized stacks.')
# Fitness Options
    parser.add_argument('-fit', action='store', dest='fit_type', default=None, choices = ['MSE','MAE','R2','VAF','combo'],
                    type=str, help='Fitness metric (Default: mse). combo is mae/r2')

    parser.add_argument('--norm_error', action='store_true', dest='ERC_ints', default=None,
                help='Flag to normalize error by the standard deviation of the target data being used.')

    parser.add_argument('--fe', action='store_true', dest='estimate_fitness', default=None,
            help='Flag to coevolve fitness estimators.')

    parser.add_argument('-fe_p', action='store', dest='FE_pop_size', default=None, type = positive_integer,
                    help='fitness estimator population size.')

    parser.add_argument('-fe_i', action='store', dest='FE_ind_size', default=None, type = positive_integer,
                    help='number of fitness cases for FE to use.')

    parser.add_argument('-fe_t', action='store', dest='FE_train_size', default=None, type = positive_integer,
                    help='trainer population size. Trainers are evaluated on the entire data set.')

    parser.add_argument('-fe_g', action='store', dest='FE_train_gens', default=None, type = positive_integer,
                    help='Number of generations between trainer selections.')

    parser.add_argument('--fe_rank', action='store_true', dest='FE_rank', default=None,
            help='User rank for FE fitness rather than error.')

# Parameter hillclimbing
    parser.add_argument('--phc', action='store_true', dest='pHC_on', default=None,
                        help='Flag to use parameter hillclimbing.')

    parser.add_argument('-phc_its', action='store', dest='pHC_its', default=None, type = positive_integer,
                        help='Number of iterations of parameter hill climbing each generation.')
# Epigenetic hillclimbing
    parser.add_argument('--ehc', action='store_true', dest='eHC_on', default=None,
                        help='Flag to use epigenetic hillclimbing.')

    parser.add_argument('-ehc_its', action='store', dest='EHC_its', default=None, type = positive_integer,
                    help='Number of iterations of epigenetic hill climbing each generation.')

    parser.add_argument('-ehc_init', action='store', dest='EHC_init', default=None, type = positive_integer,
                        help='Fraction of initial population''s genes that are silenced.')

    parser.add_argument('--emut', action='store_true', dest='EHC_mut', default=None,
                        help='Flag to use epigenetic mutation. Only works if ehc flag is present.')

    parser.add_argument('--e_slim', action='store_true', dest='EHC_slim', default=None,
                        help='Flag to store partial program outputs such that point evaluations are minimized during EHC.')
# Pareto Archive
    parser.add_argument('--archive', action='store_true', dest='prto_arch_on', default=None,
                        help='Flag to save the Pareto front of equations (fitness and complexity) during the run.')

    parser.add_argument('-arch_size', action='store', dest='prto_arch_size', default=None,
                        help='Minimum size of the Pareto archive. By default, ellyn will save the entire front, but more individuals will be added if the front is less than arch_size.')
# island model
    parser.add_argument('--islands', action='store_true', dest='islands', default=None,
                    help='Flag to use island populations. Parallel execution across codes on a single node.')

    parser.add_argument('-island_gens', action='store', dest='island_gens', default=None, type = positive_integer,
                    help='Number of generations between synchronized shuffling of island populations.')
# lexicase selection options
    parser.add_argument('-lex_pool', action='store', dest='lex_pool', default=None,
                        type=float_range, help='Fraction of population to use in lexicase selection events [0.0, 1.0].')

    parser.add_argument('-lex_metacases', nargs = '*',action='store', dest='lex_metacases', default=None,
                    help='Specify extra cases for selection. Options: age, complexity.')

    parser.add_argument('--lex_eps_error_mad', action='store_true', dest='lex_eps_error_mad', default=None,
                    help='Flag to use epsilon lexicase with median absolute deviation, error-based epsilons.')

    parser.add_argument('--lex_eps_target_mad', action='store_true', dest='lex_eps_target_mad', default=None,
                    help='Flag to use epsilon lexicase with median absolute deviation, target-based epsilons.')

    parser.add_argument('--lex_eps_dynamic', action='store_false', dest='lex_eps_global', default=None,
                    help='Flag to use dynamic epsilon lexicase selection.')

    parser.add_argument('-s', action='store', dest='random_state', default=np.random.randint(4294967295),
                        type=int, help='Random number generator seed for reproducibility. Note that using multi-threading may '
                                       'make exacts results impossible to reproduce.')

    parser.add_argument('-v', action='store', dest='VERBOSITY', default=0, choices=[0, 1, 2, 3],
                        type=int, help='How much information ellyn communicates while it is running: 0 = none, 1 = minimal, 2 = lots, 3 = all.')

    parser.add_argument('--no-update-check', action='store_true', dest='DISABLE_UPDATE_CHECK', default=False,
                        help='Flag indicating whether the ellyn version checker should be disabled.')

    # parser.add_argument('--version', action='version', version='ellyn {version}'.format(version=__version__),
    #                     help='Show ellyn\'s version number and exit.')

    args = parser.parse_args()
    #
    # if args.VERBOSITY >= 2:
    #     print('\nellyn settings:')
    #     for arg in sorted(args.__dict__):
    #         if arg == 'DISABLE_UPDATE_CHECK':
    #             continue
    #         print('{}\t=\t{}'.format(arg, args.__dict__[arg]))
    #     print('')

    # load data from csv file
    if args.INPUT_SEPARATOR is None:
        input_data = pd.read_csv(args.INPUT_FILE, sep=args.INPUT_SEPARATOR,engine='python')
    else: # use c engine for read_csv is separator is specified
        input_data = pd.read_csv(args.INPUT_FILE, sep=args.INPUT_SEPARATOR)

    if 'Label' in input_data.columns.values:
        input_data.rename(columns={'Label': 'label'}, inplace=True)

    random_state = args.random_state if args.random_state > 0 else None

    train_i, test_i = train_test_split(input_data.index,
                                                        stratify = None,#  stratify=input_data['label'].values,
                                                         train_size=0.75,
                                                         test_size=0.25,
                                                         random_state=random_state)

    training_features = np.array(input_data.loc[train_i].drop('label', axis=1).values,dtype=np.float32, order='C')
    training_labels = np.array(input_data.loc[train_i, 'label'].values,dtype=np.float32,order='C')
    # training_features = np.array(input_data.drop('label', axis=1).values,dtype=np.float32, order='C')
    # training_labels = np.array(input_data['label'].values,dtype=np.float32,order='C')
    testing_features = input_data.loc[test_i].drop('label', axis=1).values
    testing_labels = input_data.loc[test_i, 'label'].values
    # pdb.set_trace()
    learner = ellyn(**args.__dict__)
    # learner = ellyn(generations=args.GENERATIONS, population_size=args.POPULATION_SIZE,
    #             mutation_rate=args.MUTATION_RATE, crossover_rate=args.CROSSOVER_RATE,
    #             machine_learner = args.MACHINE_LEARNER, min_depth = args.MIN_DEPTH,
    #             max_depth = args.MAX_DEPTH, sel = args.SEL, tourn_size = args.TOURN_SIZE,
    #             seed_with_ml = args.SEED_WITH_ML, op_weight = args.OP_WEIGHT,
    #             erc = args.ERC, random_state=args.random_state, verbosity=args.VERBOSITY,
    #             disable_update_check=args.DISABLE_UPDATE_CHECK,fit_choice = args.FIT_CHOICE)

    learner.fit(training_features, training_labels)

    if args.VERBOSITY >= 1:
        print('\nTraining accuracy: {}'.format(learner.score(training_features, training_labels)))
        print('Holdout accuracy: {}'.format(learner.score(testing_features, testing_labels)))

    # if args.OUTPUT_FILE != '':
    #     learner.export(args.OUTPUT_FILE)


if __name__ == '__main__':
    main()
