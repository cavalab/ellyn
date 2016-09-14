# -*- coding: utf-8 -*-
"""
Copyright 2016 William La Cava

license: GNU/GPLv3

"""

import argparse
from ._version import __version__
from .evaluation import out, calc_fitness, f
from .population import *
from .variation import *
from .selection import *

from sklearn.base import BaseEstimator
from sklearn.linear_model import LassoLarsCV
from sklearn.svm import LinearSVR
from sklearn.cross_validation import train_test_split
from sklearn.metrics import r2_score
import numpy as np
import pandas as pd
import warnings
import copy
import itertools as it
import pdb
from update_checker import update_check
from joblib import Parallel, delayed

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

    def __init__(self, population_size=100, generations=100,
                 mutation_rate=0.2, crossover_rate=0.8,
                 machine_learner = 'lasso', min_depth = 1, max_depth = 5, max_depth_init = 5,
                 sel = 'tournament', tourn_size = 2, fit_choice = 'mse', op_weight = False,
                 seed_with_ml = False, erc = False, random_state=np.random.randint(4294967295), verbosity=0, scoring_function=r2_score,
                 disable_update_check=False,elitism=False):
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

        self._best_estimator = None
        self._training_features = None
        self._training_labels = None
        self._best_inds = None

        self.population_size = population_size
        self.generations = generations
        self.mutation_rate = mutation_rate
        self.crossover_rate = crossover_rate
        self.machine_learner = machine_learner
        self.min_depth = min_depth
        self.max_depth = max_depth
        self.max_depth_init = max_depth_init
        self.sel = sel
        self.tourn_size = tourn_size
        self.fit_choice = fit_choice
        self.op_weight = op_weight
        self.seed_with_ml = seed_with_ml
        self.erc = erc
        self.random_state = random_state
        self.verbosity = verbosity
        self.scoring_function = scoring_function
        self.gp_generation = 0
        self.elitism = elitism
        self.max_fit = 99999999.666
        # self.op_weight = op_weight

        if "lexicase" in self.sel and ("_vec" not in self.fit_choice or "_rel" not in self.fit_choice):
            self.fit_choice += "_vec"
        # pdb.set_trace()
        # instantiate sklearn estimator according to specified machine learner
        if (self.machine_learner.lower() == "lasso"):
            self.ml = LassoLarsCV(n_jobs=-1)
        elif (self.machine_learner.lower() == 'svr'):
            self.ml = LinearSVR()
        # elif (self.machine_learner.lower() == "distance"):
        #     self.ml = DistanceClassifier()
        else:
            self.ml = LassoLarsCV(n_jobs=-1)

        # Columns to always ignore when in an operator
        self.non_feature_columns = ['label', 'group', 'guess']

        # function set
        self.func_set = [('+',2),('-',2),('*',2),('/',2),
                         ('sin',1),('cos',1),('exp',1),('log',1),
                         ('^2',1),('^3',1),('sqrt',1)]
        # terminal set
        self.term_set = []

    def fit(self, features, labels):
        """Fit model to data"""
        np.random.seed(self.random_state)
        # setup data

        # Train-test split routine for internal validation
        ####
        train_val_data = pd.DataFrame(data=features)
        train_val_data['labels'] = labels
        # print("train val data:",train_val_data[::10])
        new_col_names = {}
        for column in train_val_data.columns.values:
            if type(column) != str:
                new_col_names[column] = str(column).zfill(10)
        train_val_data.rename(columns=new_col_names, inplace=True)
        # internal training/validation split
        train_i, val_i = train_test_split(train_val_data.index,
                                                             stratify=None,
                                                             train_size=0.75,
                                                             test_size=0.25)

        x_t = train_val_data.loc[train_i].drop('labels',axis=1).values
        x_v = train_val_data.loc[val_i].drop('labels',axis=1).values
        y_t = train_val_data.loc[train_i, 'labels'].values
        y_v = train_val_data.loc[val_i, 'labels'].values

        # Store the training features and classes for later use
        self._training_features = x_t
        self._training_labels = y_t

        ####

        # initial model
        self._best_estimator = copy.deepcopy(self.ml.fit(x_t,y_t))
        self._best_score = self._best_estimator.score(x_v,y_v)
        if self.verbosity > 2: print("initial estimator size:",self._best_estimator.coef_.shape)
        if self.verbosity > 1: print("initial score:",self._best_score)
        # create terminal set
        for i in np.arange(x_t.shape[1]):
            # (.,.,.): node type, arity, feature column index or value
            self.term_set.append(('x',0,i)) # features
            # add ephemeral random constants if flag
            if self.erc:
                self.term_set.append(('k',0,np.random.rand())) # ephemeral random constants

        # Create initial population
        pop = self.init_pop()

        # Evaluate the entire population
        # X represents a matrix of the population outputs (number os samples x population size)
        # pop.X = np.asarray(list(map(lambda I: out(I,x_t,labels), pop.individuals)))
        # pop.X = self.transform(x_t,pop.individuals,y_t)
        pop.X = np.asarray(Parallel(n_jobs=-1,backend="threading")(delayed(out)(I,x_t,y_t) for I in pop.individuals), order = 'F')
        # pdb.set_trace()
        # calculate fitness of individuals
        # fitnesses = list(map(lambda I: fitness(I,y_t,self.machine_learner),pop.X))
        fitnesses = calc_fitness(pop.X,y_t,self.fit_choice)

        # print("fitnesses:",fitnesses)
        # Assign fitnesses to inidividuals in population
        for ind, fit in zip(pop.individuals, fitnesses):
            if isinstance(fit,(list,np.ndarray)): # calc_fitness returned raw fitness values
                fit[fit < 0] = self.max_fit
                fit[np.isnan(fit)] = self.max_fit
                fit[np.isinf(fit)] = self.max_fit
                ind.fitness_vec = fit
                ind.fitness = np.mean(ind.fitness_vec)
            else:
                ind.fitness = np.nanmin([fit,self.max_fit])

        #with Parallel(n_jobs=-1,backend="threading") as parallel:
        ####################
        ### Main GP loop
        # for each generation g
        for g in np.arange(self.generations):
            if self.verbosity > 0: print(".",end='')
            if self.verbosity > 1: print(str(g)+".)",end='')
            # if self.verbosity > 1: print("population:",stacks_2_eqns(pop.individuals))
            if self.verbosity > 2: print("pop fitnesses:", ["%0.2f" % x.fitness for x in pop.individuals])
            if self.verbosity > 1: print("median fitness pop: %0.2f" % np.median([x.fitness for x in pop.individuals]))
            if self.verbosity > 1: print("best fitness pop: %0.2f" % np.min([x.fitness for x in pop.individuals]))

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                try:
                    self.ml.fit(pop.X[self.valid_loc(pop.individuals),:].transpose(),y_t)
                except ValueError as detail:
                    print("warning: ValueError in ml fit. X.shape:",pop.X[self.valid_loc(pop.individuals),:].transpose().shape,"y_t shape:",y_t.shape)
                    print("First ten entries X:",pop.X[self.valid_loc(pop.individuals),:].transpose()[:10])
                    print("First ten entries y_t:",y_t[:10])
                    print("equations:",stacks_2_eqns(pop.individuals))
                    print("ellyn parameters:",self.get_params())
                    if self.verbosity > 1: print("---\ndetailed error message:",detail)
                    raise(ValueError)

            if self.verbosity > 1: print("number of non-zero regressors:",self.ml.coef_.shape[0])
            # keep best model
            # try:
            tmp = self.ml.score(self.transform(x_v,pop.individuals)[self.valid_loc(pop.individuals),:].transpose(),y_v)
            if self.verbosity > 1: print("current ml validation score:",tmp)
            # except Exception:
            #     tmp = 0

            if tmp > self._best_score:
                self._best_estimator = copy.deepcopy(self.ml)
                self._best_score = tmp
                self._best_inds = pop.individuals[:]
                if self.verbosity > 1: print("updated best internal validation score:",self._best_score)

            offspring = []

            # clone individuals for offspring creation
            if self.sel == 'lasso':
                # for lasso, filter individuals with 0 coefficients
                offspring = copy.deepcopy(list(x for i,x in zip(self.ml.coef_, self.valid(pop.individuals)) if  i != 0))
            else:
                offspring = copy.deepcopy(self.valid(pop.individuals))

            if self.elitism: # keep a copy of the elite individual
                elite_index = np.argmin([x.fitness for x in pop.individuals])
                elite = copy.deepcopy(pop.individuals[elite_index])

            # Apply crossover and mutation on the offspring
            for child1, child2 in zip(offspring[::2], offspring[1::2]):
                if np.random.rand() < self.crossover_rate:
                    cross(child1.stack, child2.stack, self.max_depth)
                else:
                    mutate(child1.stack,self.func_set,self.term_set)
                    mutate(child2.stack,self.func_set,self.term_set)
                child1.fitness = -1
                child2.fitness = -1
            while len(offspring) < self.population_size:
                #make new offspring to replace the invalid ones
                offspring.append(Ind())
                make_program(offspring[-1].stack,self.func_set,self.term_set,np.random.randint(self.min_depth,self.max_depth+1))
                offspring[-1].stack = list(reversed(offspring[-1].stack))

            # print("offspring:",stacks_2_eqns(offspring))
            # X_offspring = self.transform(x_t,offspring)
            X_offspring = np.asarray(Parallel(n_jobs=-1,backend="threading")(delayed(out)(O,x_t,y_t) for O in offspring), order = 'F')
            # X_offspring = np.asarray([out(O,x_t,y_t) for O in offspring], order = 'F')

            F_offspring = calc_fitness(X_offspring,y_t,self.fit_choice)
            # F_offspring = parallel(delayed(f[self.fit_choice])(y_t,yhat) for yhat in X_offspring)
            # print("fitnesses:",fitnesses)
            # Assign fitnesses to inidividuals in population
            for ind, fit in zip(offspring, F_offspring):
                if isinstance(fit,(list,np.ndarray)): # calc_fitness returned raw fitness values
                    fit[fit < 0] = self.max_fit
                    fit[np.isnan(fit)] = self.max_fit
                    fit[np.isinf(fit)] = self.max_fit
                    ind.fitness_vec = fit
                    ind.fitness = np.mean(ind.fitness_vec)
                else:
                    # print("fit.shape:",fit.shape)
                    ind.fitness = np.nanmin([fit,self.max_fit])
            # if self.verbosity > 0: print("median fitness offspring: %0.2f" % np.median([x.fitness for x in offspring]))

            # Survival the next generation individuals
            if self.sel == 'tournament':
                survivors, survivor_index = tournament(pop.individuals + offspring, self.tourn_size, num_selections = len(pop.individuals))
            elif self.sel == 'lexicase':
                survivors, survivor_index = lexicase(pop.individuals + offspring, num_selections = len(pop.individuals), survival = True)
            elif self.sel == 'epsilon_lexicase':
                survivors, survivor_index = epsilon_lexicase(pop.individuals + offspring, num_selections = len(pop.individuals), survival = True)

            if self.elitism and min([x.fitness for x in survivors]) > elite.fitness:
                # if the elite individual did not survive and elitism is on, replace worst individual with elite
                rep_index = np.argmax([x.fitness for x in survivors])
                survivors[rep_index] = elite
                survivor_index[rep_index] = elite_index
            # print("current population:",stacks_2_eqns(pop.individuals))
            # print("current pop.X:",pop.X[:,:4])
            # print("offspring:",stacks_2_eqns(offspring))
            # print("current X_offspring:",X_offspring[:,:4])
            # print("survivor index:",survivor_index)
            # print("survivors:",stacks_2_eqns(survivors))
            pop.individuals[:] = survivors
            pop.X = np.vstack((pop.X, X_offspring))[survivor_index,:]
            if pop.X.shape[0] != self.population_size:
                pdb.set_trace()
            # print("new pop.X:",pop.X[:,:4])
            # pdb.set_trace()
            # pop.X = pop.X[survivor_index,:]
            #[[s for s in survivor_index if s<len(pop.individuals)],:],
                                    #  X_offspring[[s-len(pop.individuals) for s in survivor_index if s>=len(pop.individuals)],:]))
            if self.verbosity > 2: print("median fitness survivors: %0.2f" % np.median([x.fitness for x in pop.individuals]))
        # end of main GP loop
            ####################
        if self.verbosity > 0: print("finished. best internal val score:",self._best_score)
        if self.verbosity > 2: print("features:",stacks_2_eqns(self._best_inds))
        return self.score(features,labels)

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
    parser = argparse.ArgumentParser(description='A feature engineering wrapper for '
                                                 'machine learning algorithms using genetic programming.',
                                     add_help=False)

    parser.add_argument('INPUT_FILE', type=str, help='Data file to run ellyn on; ensure that the target/label column is labeled as "label".')

    parser.add_argument('-h', '--help', action='help', help='Show this help message and exit.')

    parser.add_argument('-is', action='store', dest='INPUT_SEPARATOR', default=None,
                        type=str, help='Character used to separate columns in the input file.')

    parser.add_argument('-o', action='store', dest='OUTPUT_FILE', default='',
                        type=str, help='File to export the final model.')

    parser.add_argument('-g', action='store', dest='GENERATIONS', default=100,
                        type=positive_integer, help='Number of generations to run ellyn.')

    parser.add_argument('-p', action='store', dest='POPULATION_SIZE', default=100,
                        type=positive_integer, help='Number of individuals in the GP population.')

    parser.add_argument('-mr', action='store', dest='MUTATION_RATE', default=0.8,
                        type=float_range, help='GP mutation rate in the range [0.0, 1.0].')

    parser.add_argument('-xr', action='store', dest='CROSSOVER_RATE', default=0.2,
                        type=float_range, help='GP crossover rate in the range [0.0, 1.0].')

    parser.add_argument('-ml', action='store', dest='MACHINE_LEARNER', default='lasso', choices = ['lasso','svr'],
                        type=str, help='ML algorithm to pair with features. Default: Lasso')

    parser.add_argument('-min_depth', action='store', dest='MIN_DEPTH', default=1,
                        type=positive_integer, help='Minimum length of GP programs.')

    parser.add_argument('-max_depth', action='store', dest='MAX_DEPTH', default=2,
                        type=positive_integer, help='Maximum number of nodes in GP programs.')

    parser.add_argument('-max_depth_init', action='store', dest='MAX_DEPTH_INIT', default=1,
                        type=positive_integer, help='Maximum number of nodes in initialized GP programs.')

    parser.add_argument('-op_weight', action='store', dest='OP_WEIGHT', default=1,
                        type=bool, help='Weight variables for inclusion in synthesized features based on ML scores. Default: off')

    parser.add_argument('-sel', action='store', dest='SEL', default='tournament', choices = ['tournament','lexicase','epsilon_lexicase'],
                        type=str, help='Selection method (Default: tournament)')

    parser.add_argument('-tourn_size', action='store', dest='TOURN_SIZE', default=2,
                        type=positive_integer, help='Tournament size for tournament selection (Default: 2)')

    parser.add_argument('-fit', action='store', dest='FIT_CHOICE', default='mse', choices = ['mse','mae','mdae','r2','vaf',
                        'mse_rel','mae_rel','mdae_re','r2_rel','vaf_rel'],
                        type=str, help='Fitness metric (Default: mse)')

    parser.add_argument('--seed_with_ml', action='store_true', dest='SEED_WITH_ML', default=False,
                    help='Flag to seed initial GP population with components of the ML model.')

    parser.add_argument('--elitism', action='store_true', dest='ELITISM', default=False,
                    help='Flag to force survival of best individual in GP population.')

    parser.add_argument('--erc', action='store_true', dest='ERC', default=False,
                    help='Flag to use ephemeral random constants in GP feature construction.')

    parser.add_argument('-s', action='store', dest='RANDOM_STATE', default=np.random.randint(4294967295),
                        type=int, help='Random number generator seed for reproducibility. Note that using multi-threading may '
                                       'make exacts results impossible to reproduce.')

    parser.add_argument('-v', action='store', dest='VERBOSITY', default=1, choices=[0, 1, 2, 3],
                        type=int, help='How much information ellyn communicates while it is running: 0 = none, 1 = minimal, 2 = lots, 3 = all.')

    parser.add_argument('--no-update-check', action='store_true', dest='DISABLE_UPDATE_CHECK', default=False,
                        help='Flag indicating whether the ellyn version checker should be disabled.')

    parser.add_argument('--version', action='version', version='ellyn {version}'.format(version=__version__),
                        help='Show ellyn\'s version number and exit.')

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

    learner = ellyn(generations=args.GENERATIONS, population_size=args.POPULATION_SIZE,
                mutation_rate=args.MUTATION_RATE, crossover_rate=args.CROSSOVER_RATE,
                machine_learner = args.MACHINE_LEARNER, min_depth = args.MIN_DEPTH,
                max_depth = args.MAX_DEPTH, sel = args.SEL, tourn_size = args.TOURN_SIZE,
                seed_with_ml = args.SEED_WITH_ML, op_weight = args.OP_WEIGHT,
                erc = args.ERC, random_state=args.RANDOM_STATE, verbosity=args.VERBOSITY,
                disable_update_check=args.DISABLE_UPDATE_CHECK,fit_choice = args.FIT_CHOICE)

    learner.fit(training_features, training_labels)

    if args.VERBOSITY >= 1:
        print('\nTraining accuracy: {}'.format(learner.score(training_features, training_labels)))
        print('Holdout accuracy: {}'.format(learner.score(testing_features, testing_labels)))

    if args.OUTPUT_FILE != '':
        learner.export(args.OUTPUT_FILE)


if __name__ == '__main__':
    main()
