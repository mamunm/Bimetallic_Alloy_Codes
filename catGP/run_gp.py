#!/usr/bin/env python
# -*-coding: utf-8 -*-

#online_learning.py
#Osman Mamun
#DATE CREATED: 09-13-2018

import numpy as np
from sklearn.metrics import mean_absolute_error, mean_squared_error
from importlib import import_module
import matplotlib.pyplot as plt
from time import time

class OMGP():
    '''An object to oversee a gp run and store all the information pertinent
    to that run as a self contained source.'''

    def __init__(self,
                X_train=None,
                X_test=None,
                y_train=None,
                y_test=None,
                kernel_recipe=None,
                kernel=None,
                estimator_param=None,
                scaling=None,
                scaling_params=None,
                scaling_y_train=None,
                scaling_y_test=None
                ):

        self.X_train = X_train
        self.X_test = X_test
        self.y_train = y_train
        self.y_test = y_test
        self.kernel_recipe = kernel_recipe
        self.kernel = kernel
        self.estimator_param = estimator_param
        if  scaling == True:
            self.scaling = scaling
            self.scaling_y_train = scaling_y_train
            self.scaling_y_test = scaling_y_test
            self.alpha = scaling_params['alpha']
            self.gamma = scaling_params['gamma']
        else:
            self.scaling = False

    def run_GP(self):

        est = import_module('sklearn.gaussian_process')
        estimator = getattr(est, 'GaussianProcessRegressor')

        if self.kernel is None:
            kernel = OMGP._cook_kernel(self.kernel_recipe)
        else:
            kernel = self.kernel

        if not self.estimator_param:
            estimator = estimator(kernel=kernel,
                                  n_restarts_optimizer=8,
                                  alpha=0)
        else:
            if 'kernel' not in self.estimator_param:
                 estimator = estimator(kernel=kernel,
                                       **self.estimator_param)
            else:
                 estimator = estimator(**self.estimator_param)
        if self.scaling:
            t0 = time()
            y_train = self.y_train - self.alpha * \
                      self.scaling_y_train - self.gamma
            estimator.fit(self.X_train, y_train)
            self.training_time = time() - t0
            self.delta_pred_test, self.y_pred_test_uncertainty = \
                 estimator.predict(self.X_test, return_std=True)
            self.y_pred_test = self.delta_pred_test + \
                               self.alpha * self.scaling_y_test + self.gamma
            self.testing_set_prediction_time = time() - self.training_time
            self.delta_pred_train, self.y_pred_train_uncertainty = \
                 estimator.predict(self.X_train, return_std=True)
            self.y_pred_train = self.delta_pred_train + \
                               self.alpha * self.scaling_y_train + self.gamma
            self.training_set_prediction_time = \
                time() - self.testing_set_prediction_time
        else:
            t0 = time()
            estimator.fit(self.X_train, self.y_train)
            self.training_time = time() - t0
            self.y_pred_test, self.y_pred_test_uncertainty = \
                 estimator.predict(self.X_test, return_std=True)
            self.testing_set_prediction_time = time() - self.training_time
            self.y_pred_train, self.y_pred_train_uncertainty = \
                 estimator.predict(self.X_train, return_std=True)
            self.training_set_prediction_time = \
                time() - self.testing_set_prediction_time
        self.MAE_train = mean_absolute_error(self.y_train, self.y_pred_train)
        self.MAE_test = mean_absolute_error(self.y_test, self.y_pred_test)
        self.MSE_train = mean_squared_error(self.y_train, self.y_pred_train)
        self.MSE_test = mean_squared_error(self.y_test, self.y_pred_test)
        self.final_kernel = estimator.kernel_
        self.log_marginal_likelihood = estimator.log_marginal_likelihood()

    @staticmethod
    def _cook_kernel(recipe):
        K = None
        for key, values in recipe.items():
            kern_lib = import_module('sklearn.gaussian_process.kernels')
            kern = getattr(kern_lib, key)
            if isinstance(values, list):
                if isinstance(values[0], dict):
                    if K is None:
                        K = OMGP._cook_kernel(values[0]) * kern(**values[1])
                    else:
                        K += OMGP._cook_kernel(values[0]) * kern(**values[1])
                else:
                    if K is None:
                        K = values[0] * kern(**values[1])
                    else:
                        K += values[0] * kern(**values[1])
            else:
                if K is None:
                    K = kern(**values)
                else:
                    K += kern(**values)
        return K

    def print_sample_recipe(self):
        print('''
              To cook a kernel, one can follow the guidelines listed below
              to prepare a delicious recipe for a medium rare kernel.
              1.
              Kernel : 100 * RBF(length_scale=100)
              recipe : {'RBF' : [100, {'length_scale' : 100}]}
              2.
              kernel: 100 * RBF(length_scale=100) +
                      White_kernel(noise_laevel=1)
              recipe : {'RBF' : [100, {'length_scale' : 100}],
                        'White_kernel' : {'noise_level' : 1}}
              3.
              kernel: 35**2 * RBF(length_scale=41) *
                      ExpSineSquared(length_scale=2, periodicity=1) +
                      White_kernel(noise_laevel=1)
              recipe : {'ExpSineSquared' : [{'RBF' : [35**2,
                                          {'length_scale' : 41}]},
                                          {'length_scale' : 2,
                                           'Periodicity' : 1},
                        'White_kernel' : {'noise_level' : 1}}

              ''')

    def plot_learning_curve(self,
                            n_run=10,
                            title="Learning curve"):
        """A utility function to plot the learning curve

        Parameters
        ----------
        n_run : number of datapoints to evaluate to generate the
            leraning curve.
        title : title for the generated plot.
        """
        n_data = len(self.X_train)
        n_samp = np.linspace(n_data//n_run, n_run*(n_data//n_run),
                             n_run, dtype=int)
        est = import_module('sklearn.gaussian_process')
        estimator = getattr(est, 'GaussianProcessRegressor')
        kernel = OMGP._cook_kernel(self.kernel_recipe)
        if not self.estimator_param:
            estimator = estimator(kernel=kernel,
                                  n_restarts_optimizer=4,
                                  alpha=0)
        else:
            if 'kernel' not in self.estimator_param:
                 estimator = estimator(kernel=kernel,
                                       **self.estimator_param)
            else:
                 estimator = estimator(**self.estimator_param)

        _MAE_train = []
        _MAE_test = []
        for i, n in enumerate(n_samp):
            xtrain = self.X_train[:n]
            if self.scaling:
                scaled_part_train = self.alpha * self.scaling_y_train[:n] \
                    + self.gamma
                scaled_part_test = self.alpha * self.scaling_y_test \
                    + self.gamma
                ytrain = self.y_train[:n] - scaled_part_train
                estimator.fit(xtrain, ytrain)
                _MAE_train += [mean_absolute_error(self.y_train[:n],
                            estimator.predict(xtrain) + scaled_part_train)]
                _MAE_test += [mean_absolute_error(self.y_test,
                            estimator.predict(self.X_test) + scaled_part_test)]
            else:
                ytrain = self.y_train[:n]
                estimator.fit(xtrain, ytrain)
                _MAE_train += [mean_absolute_error(ytrain,
                                               estimator.predict(xtrain))]
                _MAE_test += [mean_absolute_error(self.y_test,
                                              estimator.predict(self.X_test))]
        plt.figure()
        plt.title(title)
        plt.grid(color='b', linestyle='-', linewidth=0.5)
        plt.xlabel("# of training samples")
        plt.ylabel("MAE [eV]")
        plt.plot(n_samp, _MAE_train, 'o-', color='r',
                 label="MAE for training set")
        plt.plot(n_samp, _MAE_test, 'o-', color='g',
                 label="MAE for testing set")
        plt.legend(loc='best')

        return plt

    def parity_plot(self, data='train', err_bar=False):
        """A utility function to plot the parity plot along with uncertainties
        of prediction.

        Parameters
        ----------
        data : train | test
        """
        if data not in ['train', 'test']:
            print('data must be either train or test')
            return None
        if data == 'train':
            x = self.y_train
            y = self.y_pred_train
            dy = self.y_pred_train_uncertainty
        if data == 'test':
            x = self.y_test
            y = self.y_pred_test
            dy = self.y_pred_test_uncertainty
        plt.figure()
        plt.title('Parity plot for {}ing data'.format(data))
        plt.grid(color='b', linestyle='-', linewidth=0.5)
        plt.xlabel("DFT energy [eV]")
        plt.ylabel("Predicted energy [eV]")
        if err_bar:
            plt.errorbar(x, y, yerr=dy, fmt='.k', color='g', alpha=0.2)
        plt.scatter(x, y, marker='o', color='g', alpha=0.7,
                 label="ML predicted energy")
        plt.plot(x, x, '-', color='black',
                 label="parity line")
        plt.legend(loc='best')

        return plt

