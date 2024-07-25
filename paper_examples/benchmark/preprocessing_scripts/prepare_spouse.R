source("preprocessing_scripts/preprocess.R")
orig_data <- read.table(file = "raw_data/datafile.co")
colnames(orig_data) <- c("whrswk", "lfp", "hhi", "whi",
                         "hhi2", "ed9_12", "edhs", "smcol",
                         "college", "gradsch", "black", "othrace",
                         "hispanic", "exp", "exp2", "exp3", "kids1t6",
                         "kids618", "husby", "nc", "south", "west", "wght")

n <- nrow(orig_data)
raw_data <- data.frame(whrswk = orig_data[,"whrswk"],
                       lfp = orig_data[,"lfp"],
                       hhi = orig_data[,"hhi"],
                       whi = orig_data[,"whi"],
                       hhi2 = orig_data[,"hhi2"],
                       educ = rep("lt_hs", times = n),
                       race = rep("white", times = n),
                       hispanic = orig_data[,"hispanic"],
                       exp = orig_data[,"exp"],
                       exp2 = orig_data[,"exp2"],
                       exp3 = orig_data[,"exp3"],
                       kids1t6 = orig_data[,"kids1t6"],
                       kids618 = orig_data[,"kids618"],
                       husby = orig_data[,"husby"],
                       region = rep("other", times = n))

raw_data[which(orig_data[,"ed9_12"] == 1),"educ"] <- "some_hs"
raw_data[which(orig_data[,"edhs"] == 1), "educ"] <- "hs"
raw_data[which(orig_data[,"smcol"] == 1), "educ"] <- "some_college"
raw_data[which(orig_data[,"college"] == 1), "educ"] <- "college"
raw_data[which(orig_data[,"gradsch"] == 1), "educ"] <- "grad_school"

raw_data[which(orig_data[,"black"] == 1), "race"] <- "black"
raw_data[which(orig_data[,"othrace"] == 1), "race"] <- "other"

raw_data[which(orig_data[,"nc"] == 1), "region"] <- "north_central"
raw_data[which(orig_data[,"south"] == 1), "region"] <- "south"
raw_data[which(orig_data[,"west"] == 1), "region"] <- "west"

raw_data[,"educ"] <- factor(raw_data[,"educ"])
raw_data[,"race"] <- factor(raw_data[,"race"])
raw_data[,"region"] <- factor(raw_data[,"region"])

out_name <- c("whrswk")
cat_names <- c("educ", "race", "region")
cont_names <- c("lfp", "hhi", "whi", "hhi2",
                "hispanic", "exp", "exp2", "exp3",
                "kids1t6", "kids618", "husby")
use_cut_names <- c("lfp", "hhi", "whi", "hhi2",
                   "hispanic", "kids1t6", "kids618")

spouse_data <- preprocess(raw_data, cont_names, cat_names, out_name, use_cut_names)
save(spouse_data, file = "data/spouse.RData")
