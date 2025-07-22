## Update (July 2025)

This directory contains code, data, and instructions for replicating the analyses reported in Sections 4 and 5 of the original [flexBART paper](https://doi.org/10.1080/10618600.2024.2431072). 
These analyses were all performed using earlier versions (<= 1.2.0) of the package.
The code will eventually be updated to showcase the new **flexBART** interface.

## Overview
This directory contains the code, data, and instructions for running the analyses reported in Sections 4 and 5 of the paper.
Each subdirectory corresponds to a different experiment:

  * `grouping_advantage` contains code for the synthetic data experiments reported in Section 4.1
  * `network_constant` contains code for the network-linked regression experiment with a piecewise constant function (Section 4.2)
  * `network_fn` contains code for the network-linked regression experiment with a function that smoothly varies between two extremes (Section 4.2)
  * `benchmark`: contains code and data for the real-world benchmark datasets
  * `pitchFraming`: contains code and data for the baseball data analysis
  * `philly_crime`: contains code and data for the Philadelpha crime data analysis
