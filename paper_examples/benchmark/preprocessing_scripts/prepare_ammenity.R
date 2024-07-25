# Ammenity data
# Source: http://qed.econ.queensu.ca/jae/2002-v17.6/chattopadhyay/readme.ch.txt
source("preprocessing_scripts/preprocess.R")
tmp1 <- read.table(file = "raw_data/ammenity/headdata1.dat", header = FALSE)
tmp2 <- read.table(file = "raw_data/ammenity/headdata2.dat", header = FALSE)


colnames(tmp1) <- c("SPRICE", "NROOMS", "LVAREA", "HAGE", "LSIZE", "AIRCON",
                    "NBATH", "GARAGE", "PTAXES", "COOK", "SSPEND", "MSPEND")
colnames(tmp2) <- c("PCTWHT", "MEDINC", "DFCL", "DFNI", "PMD", "SOD",
                    "PARTICLE", "SULFUR", "RACE", "CHILDREN", "PURINC", "OHARE",
                    "MARSTAT", "RATE")

raw_data <- cbind(tmp1, tmp2)

out_name <- c("SPRICE")
cat_names <- c("AIRCON", "GARAGE")
cont_names <- colnames(raw_data)[!colnames(raw_data) %in% c(out_name, cat_names)]
use_cut_names <- c("NROOMS", "HAGE", "NBATH", "COOK", "RACE", "OHARE", "MARSTAT")

ammenity_data <- preprocess(raw_data, cont_names, cat_names, out_name, use_cut_names, log_y = TRUE)
save(ammenity_data, file = "data/ammenity.RData")
