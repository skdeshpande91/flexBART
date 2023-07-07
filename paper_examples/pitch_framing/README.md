## Pitch framing

This directory contains the data and code used to perform the analysis of the pitch framing data from Section 4.1 of the [paper](https://arxiv.org/abs/2211.04459).

The directory `data` contains an `.RData` file for each season (2013--2019). 
Each individual file contains the design matrices and an object called `folds`, which is a list of 10 lists (one for each CV fold).
Each element of `folds` records the observation indices for the training and testing data in the corresponding fold.


The scripts `BART-default.R`, `BART-alt.R`, `flexBART.R`, and `flexBART-location.R` run a different implementation of BART on a single CV fold for a single season.
You can specify the season and fold number in each script (around lines 5--10).
