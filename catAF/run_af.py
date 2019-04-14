#!/usr/bin/env python
# -*-coding: utf-8 -*-

#run_af.py
#Osman Mamun
#DATE CREATED: 10-24-2018

import numpy as np
from importlib import import_module
import matplotlib.pyplot as plt
from time import time
from catGP import OMGP
from scipy.stats import norm
import os
import random

class AFGP:
    '''An object to oversee an incremental learning process using acquisition
    function to accelerate the learning process.'''

    def __init__(self,
                 X=None,
                 y=None,
                 AF=None,
                 kernel_recipe=None,
                 estimator_param=None,
                 scaling=None,
                 scaling_params=None,
                 scaling_y=None,
                 n_steps=None,
                 train_ind=None,
                 af_kwargs=None,
                 save_loc=None,
                 chunk=False,
                 chunk_dict=None
                 ):

        self.X = X
        self.y = y
        self.AF = AF
        self.kernel_recipe = kernel_recipe
        self.estimator_param = estimator_param
        self.train_ind = train_ind
        if  scaling == True:
            self.scaling = scaling
            self.scaling_y = scaling_y
            self.scaling_params = scaling_params
        else:
            self.scaling = False

        self.af_kwargs = af_kwargs
        self.n_iter = 0
        self.n_steps = n_steps
        self.history = {}
        self.save_loc = save_loc
        self.chunk = chunk
        if self.chunk:
            self.chunk_dict = chunk_dict
            self.chunk_ind = {i: [] for i in self.chunk_dict}

        if not os.path.exists(self.save_loc):
            os.mkdir(self.save_loc)

        if self.n_iter == 0:
            if train_ind is not None:
                self.train_ind = train_ind
            else:
                if not self.chunk_dict:
                    self.train_ind = np.array([True if i < len(self.y) // 5
                                      else False for i in range(len(self.y))])
                    np.random.shuffle(self.train_ind)
                else:
                    self.train_ind_nom = []
                    for i in self.chunk_dict:
                        self.train_ind_nom += [random.sample(self.chunk[i],
                                                    len(self.chunk[i] // 5))]
                    self.train_ind = [True if i in train_ind_nom else False
                                      for i in range(len(self.X))]

    @staticmethod
    def _get_af_index(AF, mean, std, ybest, **af_kwargs):
        if AF == 'OPT':
            eval_imp = mean + std - ybest
            return np.argmin(eval_imp)

        if AF == 'PI':
            eval_imp = -((mean - ybest) / (std))
            eval_imp = norm.cdf(eval_imp)
            return np.argmin(eval_imp)

        if AF == 'EI':
            eval_imp = -((mean - ybest) / (std))
            eval_imp = -(mean - ybest) * norm.cdf(eval_imp) \
                - std * norm.cdf(eval_imp)
            return np.argmin(eval_imp)

        if AF == 'ES':
            raise NotImplementedError

        if AF == 'LCB':
            if not isinstance(af_kwargs['kappa'], list):
                eval_imp = mean - af_kwargs['kappa'] * std
                return np.argmin(eval_imp)
            else:
                ret_val = []
                for i in af_kwargs['kappa']:
                    eval_imp = mean - i * std
                    ret_val += [np.argmin(eval_imp)]
                return ret_val

        if AF == 'UCB':
            if not isinstance(af_kwargs['kappa'], list):
                eval_imp = mean + af_kwargs['kappa'] * std
                return np.argmin(eval_imp)
            else:
                ret_val = []
                for i in af_kwargs['kappa']:
                    eval_imp = mean + i * std
                    ret_val += [np.argmax(eval_imp)]
                return ret_val

    def run_AF(self):
        print('Working on iteration: {}'.format(self.n_iter))
        if self.n_iter != 0:
            self.update_train_index()
        X_tr = self.X[self.train_ind]
        X_ts = self.X[~self.train_ind]
        y_tr = self.y[self.train_ind]
        y_ts = self.y[~self.train_ind]

        if not self.scaling:
            MLGP = OMGP(X_train=X_tr,
                        X_test=X_ts,
                        y_train=y_tr,
                        y_test=y_ts,
                        kernel_recipe=self.kernel_recipe)
        else:
            scaling_y_train = self.scaling_y[self.train_ind]
            scaling_y_test = self.scaling_y[~self.train_ind]
            MLGP = OMGP(X_train=X_tr,
                        X_test=X_ts,
                        y_train=y_tr,
                        y_test=y_ts,
                        kernel_recipe=self.kernel_recipe,
                        scaling=True,
                        scaling_params=self.scaling_params,
                        scaling_y_train=scaling_y_train,
                        scaling_y_test=scaling_y_test)

        MLGP.run_GP()
        PP = MLGP.parity_plot(data='train', err_bar=True)
        PP.savefig(self.save_loc + '/parity_plot_train_{}.png'.format(
            self.n_iter))
        PP.close()
        PP = MLGP.parity_plot(data='test', err_bar=True)
        PP.savefig(self.save_loc + '/parity_plot_test_{}.png'.format(
            self.n_iter))
        PP.close()
        self.history[self.n_iter] = MLGP.__dict__
        self.history[self.n_iter]['train_ind'] = self.train_ind.copy()
        self.n_iter += 1


    def update_train_index(self):
        histroy_pre = self.history[self.n_iter-1]
        mean_pred = histroy_pre['y_pred_test']
        std_pred = histroy_pre['y_pred_test_uncertainty']
        y_best = min(histroy_pre['y_train'])
        if self.af_kwargs is not None:
            af_index = AFGP._get_af_index(self.AF,
                                          mean_pred,
                                          std_pred,
                                          y_best,
                                          **self.af_kwargs)
        else:
            if not self.chunk:
                af_index = AFGP._get_af_index(self.AF,
                                              mean_pred,
                                              std_pred,
                                              y_best)
            else:
                '''
                Work on that later
                for key in self.chunk:
                    mean_chunk =
                    std_chunk =
                    af_index = AFGP._get_af_index(self.AF,
                                                  mean_chunk,
                                                  std_chunk,
                                                  y_best)

                '''
                raise NotImplementedError('This is a work in progress.')
        false_ind = [i for i, j in enumerate(self.train_ind)
                     if j == False]
        if isinstance(af_index, list):
            for i in af_index:
                self.train_ind[false_ind[i]] = True
        else:
            self.train_ind[false_ind[af_index]] = True

    def main_run(self):
        while self.n_iter < self.n_steps:
            self.run_AF()
        np.save(self.save_loc + '/af_data.npy', self.history)
        self.plot_learning_curve()

    def plot_learning_curve(self):

        n_samp = []
        _MAE_train = []
        _MAE_test = []
        for key in self.history:
            n_samp += [len(self.history[key]['y_train'])]
            _MAE_train += [self.history[key]['MAE_train']]
            _MAE_test += [self.history[key]['MAE_test']]

        plt.figure()
        plt.title('Learning curve')
        plt.grid(color='b', linestyle='-', linewidth=0.5)
        plt.xlabel("# of training samples")
        plt.ylabel("MAE [eV]")
        plt.plot(n_samp, _MAE_train, 'o-', color='r',
                 label="MAE for training set")
        plt.plot(n_samp, _MAE_test, 'o-', color='g',
                 label="MAE for testing set")
        plt.legend(loc='best')
        plt.savefig(self.save_loc + '/Learning_Curve.png')
        plt.close()


