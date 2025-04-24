source("../flexBART/R/parse_formula.R")
source("generate_nested_example.R")

frmla <- 
  Y ~ bart(Blkgrp + Tract) + bart(Classroom + School + District) + bart(.)
  

tmp <- parse_formula(frmla, colnames(train_data))
tmp$cov_ensm
# make sure X1, X2 are first