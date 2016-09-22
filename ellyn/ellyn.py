# -*- coding: utf-8 -*-
"""
Copyright 2016 William La Cava

license: GNU/GPLv3

"""

import argparse
# from ._version import __version__

from sklearn.base import BaseEstimator

import numpy as np
import pandas as pd
import warnings
import copy
import itertools as it
import pdb
import elgp
from update_checker import update_check
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

    def __init__(self, **kwargs):
                # sets up GP.

        # Save params to be recalled later by get_params()
        self.params = locals()  # Must be placed before any local variable definitions
        self.params.pop('self')

        # Do not prompt the user to update during this session if they ever disabled the update check
        if disable_update_check:
            ellyn.update_checked = True

        # Prompt the user if their version is out of date
        if not disable_update_check and not ellyn.update_checked:
            update_check('ellyn', __version__)
            ellyn.update_checked = True

        self.param_dict = kwargs.__dict__
        # self.param_file = ''
        # for arg in sorted(kwargs.__dict__):
        #     self.param_file += ('{}\t{}'.format(lower(arg), kwargs.__dict__[arg]))
        pdb.set_trace()
        self._best_estimator = None
    def fit(self, features, labels):
        """Fit model to data"""
        np.random.seed(self.random_state)
        # setup data
        # setup call to ellenGP library
        # curdir = os.path.split(os.path.abspath( __file__))[0]
        # curdirList = curdir.split(os.path.sep)
        # reporootdir = os.path.sep.join(curdirList[:-2])
        # lib = ctypes.cdll.LoadLibrary(
        # os.path.join(reporootdir,'lib/ellengp.so'))
        # lib.runEllenGP.argtypes = [
        #     ndpointer(ctypes.c_float),
        #     ndpointer(ctypes.c_float),
        #     ctypes.c_int,
        #     ctypes.c_int]
        # run GP code
        elgp.runEllenGP(self.param_dict,features,labels,0,0)
        ####

        # initial model


    def predict(self, testing_features):
        """predict on a holdout data set."""
        # print("best_inds:",self._best_inds)
        # print("best estimator size:",self._best_estimator.coef_.shape)
        if self._best_inds is None:
            return self._best_estimator.predict(testing_features)
        else:
            X_transform = (np.asarray(list(map(lambda I: out(I,testing_features), self._best_inds))))
            return self._best_estimator.predict(X_transform[self.valid_loc(self._best_inds),:].transpose())

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
        return self.scoring_function(testing_labels,yhat)

    def export(self, output_file_name):
        """exports engineered features

        Parameters
        ----------
        output_file_name: string
            String containing the path and file name of the desired output file

        Returns
        -------
        None

        """
        if self._best_estimator is None:
            raise ValueError('A model has not been optimized. Please call fit() first.')

        # Have the exported code import all of the necessary modules and functions
        # write model form from coefficients and features
        model = ''
        sym_model = ''
        x = 0
        for c,i in zip(self._best_estimator.coef_,stacks_2_eqns(self._best_inds)):
           if c != 0:
               if model:
                   model+= '+'
                   sym_model += '+'
               model+= str(c) + '*' + str(i)
               sym_model += "k_" + str(x) + '*' + str(i)
               x += 1

        print_text = "exact_model: " + model
        print_text += "\nsymbolic_model: " + sym_model
        print_text += "\ncoefficients: " + str([c for c in self._best_estimator.coef_ if c != 0])
        print_text += "\nfeatures: " + str([s for s,c in zip(stacks_2_eqns(self._best_inds),self._best_estimator.coef_) if c!=0])

        with open(output_file_name, 'w') as output_file:
            output_file.write(print_text)

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
    parser.add_argument('-g', action='store', dest='GENERATIONS', default=100,
                        type=positive_integer, help='Number of generations to run ellyn.')

    parser.add_argument('-p', action='store', dest='POPULATION_SIZE', default=500,
                        type=positive_integer, help='Number of individuals in the GP population.')

    parser.add_argument('--limit_evals', action='store_true', dest='LIMIT_EVALS', default=False,
                    help='Limit evaluations instead of generations.')

    parser.add_argument('-me', action='store', dest='MAX_EVALS', default=0,
                        type=float_range, help='Max point evaluations.')

    parser.add_argument('-sel', action='store', dest='SEL', default='tournament', choices = ['tournament','lexicase','epsilon_lexicase'],
                        type=str, help='Selection method (Default: tournament)')

    parser.add_argument('-PS_sel', action='store', dest='PS_SEL', default=1, choices = [1,2,3,4,5,6],
                        type=str, help='objectives for pareto survival. 1: age + fitness; 2: age+fitness+generality; 3: age+fitness+complexity; 4: class fitnesses (classification ONLY); 5: class fitnesses+ age (classification ONLY)')

    parser.add_argument('-tourn_size', action='store', dest='TOURN_SIZE', default=2,
                        type=positive_integer, help='Tournament size for tournament selection (Default: 2)')

    parser.add_argument('-mr', action='store', dest='MUTATION_RATE', default=0.2,
                        type=float_range, help='GP mutation rate in the range [0.0, 1.0].')

    parser.add_argument('-xr', action='store', dest='CROSSOVER_RATE', default=0.8,
                        type=float_range, help='GP crossover rate in the range [0.0, 1.0].')

    parser.add_argument('-rr', action='store', dest='REPRODUCTION_RATE', default=0.2,
                        type=float_range, help='GP reproduction rate in the range [0.0, 1.0].')

    parser.add_argument('--elitism', action='store_true', dest='ELITISM', default=False,
                help='Flag to force survival of best individual in GP population.')

    parser.add_argument('--init_validate', action='store_true', dest='INIT_VALIDATE_ON', default=False,
                help='Flag to guarantee initial population outputs are valid (no NaNs or Infs).')

    parser.add_argument('--no_stop', action='store_false', dest='STOP_CONDITION', default=True,
                    help='Flag to keep running even though fitness < 1e-6 has been reached.')
# Data Options
    parser.add_argument('--validate', action='store_true', dest='TRAIN', default=False,
                help='Flag to split data into training and validation sets.')

    parser.add_argument('-train_split', action='store', dest='TRAIN_PCT', default=0.5, type = float_range,
                help='Fraction of data to use for training in [0,1]. Only used when --validate flag is present.')

    parser.add_argument('--shuffle', action='store_true', dest='SHUFFLE_DATA', default=False,
                help='Flag to shuffle data samples.')

    parser.add_argument('--restart', action='store_true', dest='POP_RESTART_PATH', default=False,
                help='Flag to restart from previous population.')

    parser.add_argument('-pop', action='store_true', dest='POP_RESTART', default=False,
                help='Population to use in restart. Only used when --restart flag is present.')

    parser.add_argument('--AR', action='store_true', dest='AR', default=False,
                    help='Flag to use auto-regressive variables.')

    parser.add_argument('--AR_lookahead', action='store_true', dest='AR_LOOKAHEAD', default=False,
                    help='Flag to only estimate on step ahead of data when evaluating candidate models.')

    parser.add_argument('-AR_na', action='store', dest='AR_NA', default=1,
                    type=positive_integer, help='Order of auto-regressive output (y).')
    parser.add_argument('-AR_nb', action='store', dest='AR_NB', default=1,
                    type=positive_integer, help='Order of auto-regressive inputs (x).')
    parser.add_argument('-AR_nka', action='store', dest='AR_NKA', default=1,
                    type=positive_integer, help='AR state (output) delay.')
    parser.add_argument('-AR_nkb', action='store', dest='AR_NKB', default=1,
                    type=positive_integer, help='AR input delay.')

# Results and Printing Options
    parser.add_argument('-o', action='store', dest='RESULTS_PATH', default='',
                        type=str, help='Path where results will be saved.')

    parser.add_argument('--print_every_pop', action='store_true', dest='PRINT_EVERY_POP', default=False,
                    help='Flag to seed initial GP population with components of the ML model.')

    parser.add_argument('--print_genome', action='store_true', dest='PRINT_GENOME', default=False,
                    help='Flag to prints genome for visualization.')

    parser.add_argument('--print_novelty', action='store_true', dest='PRINT_NOVELTY', default=False,
                    help='Flag to print number of unique output vectors.')

    parser.add_argument('--print_homology', action='store_true', dest='PRINT_HOMOLOGY', default=False,
                    help='Flag to print genetic homology in programs.')

    parser.add_argument('-num_log_pts', action='store', dest='num_log_pts', default=0,
                        type=positive_integer, help='number of log points to print (0 means print each generation)')

# Classification Options

    parser.add_argument('--class', action='store_true', dest='CLASSIFICATION', default=False,
                    help='Flag to define a classification problem instead of regression (the default).')

    parser.add_argument('--class_bool', action='store_true', dest='CLASS_BOOL', default=False,
                    help='Flag to interpret class labels as bit-string conversion of boolean stack output.')

    parser.add_argument('--class_m4gp', action='store_true', dest='CLASS_M3GP', default=False,
                    help='Flag to use the M4GP algorithm.')

    parser.add_argument('--class_prune', action='store_true', dest='CLASS_PRUNE', default=False,
                        help='Flag to prune the dimensions of the best individual each generation.')

# Terminal and Operator Options
    parser.add_argument('-op_list', nargs = '*', action='store', dest='OP_LIST', default=None,
                    type=float, help='Operator list. Default: +,-,*,/,n,v. available operators: n v + - * / sin cos log exp sqrt = ! < <= > >= if-then if-then-else & |')

    parser.add_argument('-op_weight', nargs = '*',action='store', dest='OP_WEIGHTS', default=None,
                    help='Operator weights for each element in operator list. If not specified all operators have the same weights.')

    parser.add_argument('-constants', nargs = '*', action='store', dest='CVALS', default=None,
                        type=float, help='Seed GP initialization with constant values.')

    parser.add_argument('-seeds', nargs = '*', action='store', dest='SEEDS', default=None,
                        type=float, help='Seed GP initialization with partial solutions, e.g. (x+y). Each partial solution must be enclosed in parentheses.')

    parser.add_argument('-min_len', action='store', dest='MIN_LEN', default=1,
                        type=positive_integer, help='Minimum length of GP programs.')

    parser.add_argument('-max_len', action='store', dest='MAX_LEN', default=2,
                        type=positive_integer, help='Maximum number of nodes in GP programs.')

    parser.add_argument('-max_len_init', action='store', dest='MAX_LEN_INIT', default=1,
                        type=positive_integer, help='Maximum number of nodes in initialized GP programs.')

    parser.add_argument('--erc_ints', action='store_true', dest='ERC_INTS', default=False,
                help='Flag to use integer instead of floating point ERCs.')

    parser.add_argument('--no_erc', action='store_false', dest='ERC', default=True,
                help='Flag to turn of ERCs. Useful if you are specifying constant values and don''t want random ones.')

    parser.add_argument('-min_erc', action='store', dest='MINERC', default=-1, type = int,
                    help='Minimum ERC value.')

    parser.add_argument('-max_erc', action='store', dest='MAXERC', default=1, type = int,
                    help='Maximum ERC value.')

    parser.add_argument('-num_erc', action='store', dest='NUMERC', default=1, type = int,
                    help='Number of ERCs to use.')

    parser.add_argument('--trees', action='store_true', dest='INIT_TREES', default=False,
                    help='Flag to initialize genotypes as syntactically valid trees rather than randomized stacks.')
# Fitness Options
    parser.add_argument('-fit', action='store', dest='FIT_TYPE', default='mse', choices = ['mse','mae','r2','vaf','combo'],
                    type=str, help='Fitness metric (Default: mse). combo is mae/r2')

    parser.add_argument('--norm_error', action='store_true', dest='ERC_INTS', default=False,
                help='Flag to normalize error by the standard deviation of the target data being used.')

    parser.add_argument('--fe', action='store_true', dest='ESTIMATE_FITNESS', default=False,
            help='Flag to coevolve fitness estimators.')

    parser.add_argument('-fe_p', action='store', dest='FE_POP_SIZE', default=1, type = positive_integer,
                    help='fitness estimator population size.')

    parser.add_argument('-fe_i', action='store', dest='FE_IND_SIZE', default=1, type = positive_integer,
                    help='number of fitness cases for FE to use.')

    parser.add_argument('-fe_t', action='store', dest='FE_TRAIN_SIZE', default=1, type = positive_integer,
                    help='trainer population size. Trainers are evaluated on the entire data set.')

    parser.add_argument('-fe_g', action='store', dest='FE_TRAIN_GENS', default=1, type = positive_integer,
                    help='Number of generations between trainer selections.')

    parser.add_argument('--fe_rank', action='store_true', dest='FE_RANK', default=False,
            help='User rank for FE fitness rather than error.')

# Parameter hillclimbing
    parser.add_argument('--phc', action='store_true', dest='PHC_ON', default=False,
                        help='Flag to use parameter hillclimbing.')

    parser.add_argument('-phc_its', action='store', dest='pHC_its', default=1, type = positive_integer,
                        help='Number of iterations of parameter hill climbing each generation.')
# Epigenetic hillclimbing
    parser.add_argument('--ehc', action='store_true', dest='EHC_ON', default=False,
                        help='Flag to use epigenetic hillclimbing.')

    parser.add_argument('-ehc_its', action='store', dest='EHC_its', default=1, type = positive_integer,
                    help='Number of iterations of epigenetic hill climbing each generation.')

    parser.add_argument('-ehc_init', action='store', dest='EHC_INIT', default=0.5, type = positive_integer,
                        help='Fraction of initial population''s genes that are silenced.')

    parser.add_argument('--emut', action='store_true', dest='EHC_MUT', default=False,
                        help='Flag to use epigenetic mutation. Only works if ehc flag is present.')

    parser.add_argument('--e_slim', action='store_true', dest='EHC_SLIM', default=False,
                        help='Flag to store partial program outputs such that point evaluations are minimized during EHC.')
# Pareto Archive
    parser.add_argument('--archive', action='store_true', dest='PRTO_ARCH_ON', default=False,
                        help='Flag to save the Pareto front of equations (fitness and complexity) during the run.')

    parser.add_argument('-arch_size', action='store', dest='PRTO_ARCH_SIZE', default=20,
                        help='Minimum size of the Pareto archive. By default, ellyn will save the entire front, but more individuals will be added if the front is less than arch_size.')
# island model
    parser.add_argument('--islands', action='store_true', dest='ISLANDS', default=False,
                    help='Flag to use island populations. Parallel execution across codes on a single node.')

    parser.add_argument('-island_gens', action='store', dest='ISLAND_GENS', default=100, type = positive_integer,
                    help='Number of generations between synchronized shuffling of island populations.')
# lexicase selection options
    parser.add_argument('-lex_pool', action='store', dest='LEX_POOL', default=1,
                        type=float_range, help='Fraction of population to use in lexicase selection events [0.0, 1.0].')

    parser.add_argument('-lex_metacases', nargs = '*',action='store', dest='LEX_METACASES', default=None,
                    help='Specify extra cases for selection. Options: age, complexity.')

    parser.add_argument('--lex_eps_error_mad', action='store_true', dest='LEX_EPS_ERROR_MAD', default=False,
                    help='Flag to use epsilon lexicase with median absolute deviation, error-based epsilons.')

    parser.add_argument('-s', action='store', dest='RANDOM_STATE', default=np.random.randint(4294967295),
                        type=int, help='Random number generator seed for reproducibility. Note that using multi-threading may '
                                       'make exacts results impossible to reproduce.')

    parser.add_argument('-v', action='store', dest='VERBOSITY', default=1, choices=[0, 1, 2, 3],
                        type=int, help='How much information ellyn communicates while it is running: 0 = none, 1 = minimal, 2 = lots, 3 = all.')

    parser.add_argument('--no-update-check', action='store_true', dest='DISABLE_UPDATE_CHECK', default=False,
                        help='Flag indicating whether the ellyn version checker should be disabled.')

    # parser.add_argument('--version', action='version', version='ellyn {version}'.format(version=__version__),
    #                     help='Show ellyn\'s version number and exit.')

    args = parser.parse_args()

    if args.VERBOSITY >= 2:
        print('\nellyn settings:')
        for arg in sorted(args.__dict__):
            if arg == 'DISABLE_UPDATE_CHECK':
                continue
            print('{}\t=\t{}'.format(arg, args.__dict__[arg]))
        print('')

    # load data from csv file
    if args.INPUT_SEPARATOR is None:
        input_data = pd.read_csv(args.INPUT_FILE, sep=args.INPUT_SEPARATOR,engine='python')
    else: # use c engine for read_csv is separator is specified
        input_data = pd.read_csv(args.INPUT_FILE, sep=args.INPUT_SEPARATOR)

    if 'Label' in input_data.columns.values:
        input_data.rename(columns={'Label': 'label'}, inplace=True)

    RANDOM_STATE = args.RANDOM_STATE if args.RANDOM_STATE > 0 else None

    train_i, test_i = train_test_split(input_data.index,
                                                        stratify = None,#  stratify=input_data['label'].values,
                                                         train_size=0.75,
                                                         test_size=0.25,
                                                         random_state=RANDOM_STATE)

    training_features = input_data.loc[train_i].drop('label', axis=1).values
    training_labels = input_data.loc[train_i, 'label'].values

    testing_features = input_data.loc[test_i].drop('label', axis=1).values
    testing_labels = input_data.loc[test_i, 'label'].values

    learner = ellyn(args.__dict__)
    # learner = ellyn(generations=args.GENERATIONS, population_size=args.POPULATION_SIZE,
    #             mutation_rate=args.MUTATION_RATE, crossover_rate=args.CROSSOVER_RATE,
    #             machine_learner = args.MACHINE_LEARNER, min_depth = args.MIN_DEPTH,
    #             max_depth = args.MAX_DEPTH, sel = args.SEL, tourn_size = args.TOURN_SIZE,
    #             seed_with_ml = args.SEED_WITH_ML, op_weight = args.OP_WEIGHT,
    #             erc = args.ERC, random_state=args.RANDOM_STATE, verbosity=args.VERBOSITY,
    #             disable_update_check=args.DISABLE_UPDATE_CHECK,fit_choice = args.FIT_CHOICE)

    learner.fit(training_features, training_labels)

    if args.VERBOSITY >= 1:
        print('\nTraining accuracy: {}'.format(learner.score(training_features, training_labels)))
        print('Holdout accuracy: {}'.format(learner.score(testing_features, testing_labels)))

    if args.OUTPUT_FILE != '':
        learner.export(args.OUTPUT_FILE)


if __name__ == '__main__':
    main()
