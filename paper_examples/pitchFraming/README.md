## Pitch framing

This directory contains the data and code used to perform the analysis of the pitch framing data from Section 5.2 of the [paper](https://arxiv.org/abs/2211.04459).
The subdirectory `data` contains an `.RData` file for each season (2013--2019). 

The file `sim_settings.RData` contains a data frame `sim_settings` that tabulates every combination of method, year, and simulation number.
The script `pitchFraming.R` takes as input a row index (`job_id`) and runs the corresponding method on a training-testing split for the corresponding year.

To run the experiment for a single row `job_id` locally, you should comment out the lines
```
args <- commandArgs(TRUE)
job_id <- as.numeric(args[1])
```
and instead manually set `job_id` (around line 5 of the script)

After running the experiment for all rows, you can tabulate hte results by running the script `compile_results.R`


