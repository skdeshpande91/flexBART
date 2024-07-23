## Philadelphia crime data

This directory contains code for the piecewise constant network-linked regression experiment.
The file `sim_settings.RData` contains a data frame `sim_settings` that tabulates all combinations of methods and simulation number.
The file `philly_obs_cv_simulation.R` loops over a block of rows of `sim_settings`; generates the relevant dataset; and runs the method corresponding to that row.
It was written to be run on a high-throughput computing cluster where different compute nodes looped over different blocks of rows, specified by a command-line argument.

To run the experiment for rows `start_ix` to `end_ix` locally, you should comment out the lines
```
args <- commandArgs(TRUE)
block_id <- as.numeric(args[1])
```
and change the bound of the for loop from `job_id in block_starts[block_id]:block_ends[block_id])` to `job_id in start_ix:end_ix)`.

After running the full experiment, you can tabulate the results by running the script `compile_results.R`

The remaining files are used to generate the data and implement the target-encoded BART and oracle procedures


## Network-splitting strategies

The **flexBART** package implements a few more network-splitting strategies than what is reported in the paper.The strategies reported in the paper correspond to the following specifications of the `graph_split_option` argument of `networkBART`:
  * `gs1` corresponds to `graph_split_option=1`
  * `gs2` corresponds to `graph_split_option=3`
  * `gs3` corresponds to `graph_split_option=5`
  * `gs4` corresponds to `graph_split_option=6`

