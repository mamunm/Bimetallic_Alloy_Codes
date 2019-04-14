#!/usr/bin/env python
# -*-coding: utf-8 -*-

#data_cleanup.py
#Osman Mamun
#LAST UPDATED: 09-11-2018

import numpy as np
from sklearn.preprocessing import Imputer
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import MinMaxScaler
import sklearn.preprocessing as skpre
from sklearn.decomposition import PCA

class preprocess_data():
    """ An object to hold the data for preprocessing and later it will return
    the processsed data when needed

    TODO: Metadata modification
    """

    def __init__(self,
                 X = None,
                 y = None):
        self.X = X
        self.y = y
        #self.metadata=metadata

    def get_data(self):
        """Returns the data at any point it was called"""

        return self.X, self.y

    def remove_instance(self, null_count=0.5):
        """function to remove any data instance with more than
        null_count * 100% null values"""

        mask = [i for i, XX in enumerate(self.X)
                if np.isnan(XX).mean() < null_count]
        self.X = self.X[mask]
        self.y = self.y[mask]
        #self.metadata = self.metadata[mask]

    def remove_null_features(self, null_count=0.5):
        """function to remove any feature with more than
        null_count * 100% null values"""

        mask = [i for i in range(self.X.shape[1])
                if np.isnan(self.X[:, i]).mean() < null_count]

        self.X = self.X[:, mask]

    def remove_low_variation(self, var_threshold=0.05):
        """function to remove data with little information as
        characterize by their variance"""

        mask = [i for i in range(self.X.shape[1])
                if self.X[:, i].var() > var_threshold]

        self.X = self.X[:, mask]

    def impute_data(self, strategy='mean'):
        """Uses sklearn Imputer to impute null values with mean,
        median or most_frequent depending on the user input"""

        imp = Imputer(strategy=strategy)
        self.X = imp.fit_transform(self.X)

    def scale_data(self, strategy='MinMaxScaler'):
        """Uses skleran StandardScaler or MinMaxScaler to scale data"""

        if isinstance(strategy, str):
            strategy = getattr(skpre, strategy)
        scale = strategy()
        self.X = scale.fit_transform(self.X)

    def clean_data(self):
        """Performs all the operations available."""
        self.remove_instance()
        self.remove_null_features()
        self.remove_low_variation()
        self.impute_data()
        self.scale_data()

    def get_PCA(self, n_pc=10):
        self.X_PCA = PCA(n_components=n_pc).fit_transform(self.X)

    def remove_highly_correlated_features():
        pass
