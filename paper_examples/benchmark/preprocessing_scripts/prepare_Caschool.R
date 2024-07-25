source("preprocessing_scripts/preprocess.R")
library(Ecdat)
data("Caschool")

raw_data <- Caschool
raw_data[,"grspan"] <- as.numeric(factor(raw_data[,"grspan"])) - 1 # grspan is binary
out_name <- c("testscr") # average of the reading and math scores
cat_names <- c("county")
cont_names <- c("grspan", "enrltot", "teachers", "calwpct", "mealpct", "computer",
               "testscr", "compstu", "expnstu", "str", "avginc", "elpct")
use_cut_names <- c()

Caschool_data <- preprocess(raw_data, cont_names, cat_names, out_name, use_cut_names)
save(Caschool_data, file = "data/Caschool.RData")

