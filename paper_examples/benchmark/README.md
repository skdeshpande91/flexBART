## Benchmark datasets

This directory contains the data and code used to perform the analysis of the benchmark datasets from Section 5.1 of the [paper](https://arxiv.org/abs/2211.04459).
The sub-directories `raw_data` and `preprocessing_scripts` contain the raw data files and scripts used to preprocess the data, respectively.
The sub-directory `data` contains `.RData` files for each benchmark datasets.

The file `sim_settings.RData` contains a data frame `sim_settings` that tabulates every combination of dataset and simulation number.
The script `benchmark.R` takes as input a block of rows of `sim_settings` and the runs every method on the corresponding combinations of datasets and simulation replications.

To run the experiment for rows `start_ix` to `end_ix` locally, you should comment out the lines
```
args <- commandArgs(TRUE)
job_id <- as.numeric(args[1])
```
and change the range of the for loop from `job_id in block_starts[block_id]:block_ends[block_id])` to `job_id in start_ix:end_ix`

After running the experiment for all rows, you can tabulate hte results by running the script `compile_results.R`


