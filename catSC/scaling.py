#!/usr/bin/env python
# -*-coding: utf-8 -*-

#scaling.py
#Osman Mamun
#DATE CREATED: 10-09-2018

from scipy.stats import linregress
import matplotlib.pyplot as plt

class Linear_Scaling():
    '''An object to get the scaling coefficients'''

    def __init__(self,
                 x=None,
                 y=None,
                 x_name=None,
                 y_name=None,
                 e_unit='eV'):
        self.x = x
        self.y = y
        self.x_name = x_name
        self.y_name = y_name
        self.e_unit = e_unit

    def get_coeff(self):

        lr = linregress(self.x, self.y)

        self.slope = lr[0]
        self.intercept = lr[1]
        self.r_value = lr[2]
        self.p_value = lr[3]
        self.std_err = lr[4]
        self.r_squared = lr[2]**2
        self.y_pred = self.slope * self.x + self.intercept
        self.MAE = sum([abs(i - j) for i, j in zip(self.y_pred, self.y)])/len(self.y)
    def plot_scaling(self):
        '''A utility function to plot the linear scaling trend.'''

        plt.figure()
        title = "Scaling relation between {0}* and {1}*".format(self.y_name,
                                                                self.x_name)
        plt.title(title)
        plt.grid(color='b', linestyle='-', linewidth=0.5)
        plt.xlabel(r"$E_{0}\;[{1}]$".format('{' + self.x_name + '}',
                                            self.e_unit))
        plt.ylabel(r"$E_{0}\;[{1}]$".format('{' + self.y_name + '}',
                                            self.e_unit))
        plt.scatter(self.x, self.y, marker='o', color='g', label='data points')
        plt.plot(self.x, self.y_pred, '-', color='black', alpha=0.7,
                 label="Scaling line")
        text = "slope = {0:0.2f}\nintercept = {1:0.2f}\nr^2 = {2:0.2f}\nMAE={3:0.2f}"
        text = text.format(self.slope, self.intercept, self.r_squared, self.MAE)
        x = max(self.x)/1.2
        y = max(self.y)/5
        plt.text(x, y, text,
                 ha='center', va='center', fontsize=10,
                 bbox=dict(facecolor='w',
                           edgecolor='black',
                           boxstyle='round',
                           pad=1))
        plt.legend(loc='best')

        return plt


