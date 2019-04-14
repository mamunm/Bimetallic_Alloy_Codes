# Bimetallic_Alloy_Codes

In this repository, I updated all the codes I used to generate the bimetallic alloy datasets and to create machine learning model based on the data.

## Optimization

Bulk structures are generated using ASE. Here, I used a base template for all the bulk material and then used another python script to change the portions relevant for particular metal. All the codes used to generate and optimize bulk can be found [here](bulk_optimization). Optimization of the surface structure and analysis of the calculation output can be found [here](slab_optimization). Optimization of the adsorbate-surface structure and analysis of the calculation output can be found [here](ads_slab_optimization).

## Analysis

Code to 1. identify all the unique sites, 2. to identify surface reconstruction, and 3. adsorbate reorientation can be found [here](Analysis)

## Fireworks

Code to submit jobs to the workflow and manage them later will be found [here](Fireworks)

## catGP

A python wrapper to handle sklearn Gaussian process training and testing and plotting the data.

## catSC

A python code to get the scaling coefficients and to plot the scaling line.

## catAF

A python code to perform active learning with sklearn Gaussian process.

## pawprint

A python code to generate convolution fingerprint for machine learning model generation.

## get_db_fp

Generate all the fingerprints using pawprint and store them in a sqlite db.

## GP

Gaussian process training of the dataset.

## scaling

Scaling relation trianing of the dataset.

## scaling_GP

Delta learning training of the dataset.

## Getting the data from the Catalysis-Hub.org

The whole bimetallic alloy dataset can be fetched from the Catalysis-Hub.org. Code to get the data is avialable [here](get_data/get_data.py). Also, all the log files can be downloaded using the [get_log_file](get_data/get_log_file.py) script.

## Generating the fingerprints

Then I used `pawprint` and `catlearn` to generate the fingerprint to train my machine learning models. Code used to generate the fingerprints are available [here](get_db_fp).
