source("preprocessing_scripts/preprocess.R")

raw_data <- read.table("raw_data/machine.data", sep = ",", header = FALSE)
colnames(raw_data) <- c("vendor", "model", "myct", "mmin", "mmax",
                        "cach", "chmin", "chmax", "prp", "erp")

out_name <- c("prp")
cont_names <- c("myct", "mmin", "mmax", "cach", "chmin", "chmax")
cat_names <- c("vendor")
use_cut_names <- cont_names

cpu_data <- preprocess(raw_data, cont_names, cat_names, out_name, use_cut_names)
save(cpu_data, file = "data/cpu.RData")
