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


def plot_learning_curve(X_train=None,
                        y_train=None,
                        X_test=None,
                        y_test=None,
                        kernel=None,
                        n_run=10,
                        title="Learning curve"):
    """A utility function to plot the learning curve

    Parameters
    ----------
    X_train : feature vector for all the training data.
    y_train : target value for all the training data.
    X_train : feature vector for all the testing data.
    X_train : target value for all the testing data.
    kernel : scikit learn compatible kernel.
    n_run : number of datapoints to evaluate to generate the
          leraning curve.
    title : title for the generated plot.
    """
    n_data = len(self.X_train)
    n_samp = np.linspace(n_data//n_run, n_run*(n_data//n_run),
                         n_run, dtype=int)
    est = import_module('sklearn.gaussian_process')
    estimator = getattr(est, 'GaussianProcessRegressor')


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
            estimator = estimator(kernel=kernel,
                                  n_restarts_optimizer=4,
                                  alpha=0)
            estimator.fit(xtrain, ytrain)
            _MAE_train += [mean_absolute_error(self.y_train[:n],
                           estimator.predict(xtrain) + scaled_part_train)]
            _MAE_test += [mean_absolute_error(self.y_test,
                          estimator.predict(self.X_test) + scaled_part_test)]
            kernel = estimator.kernel_
        else:
            ytrain = self.y_train[:n]
            estimator = estimator(kernel=kernel,
                                  n_restarts_optimizer=4,
                                  alpha=0)
            estimator.fit(xtrain, ytrain)
            _MAE_train += [mean_absolute_error(ytrain,
                                               estimator.predict(xtrain))]
            _MAE_test += [mean_absolute_error(self.y_test,
                                              estimator.predict(self.X_test))]
            kernel = estimator.kernel_

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

