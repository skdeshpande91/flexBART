source("preprocessing_scripts/preprocess.R")
raw_data <- read.table(file = "raw_data/servo.data", sep = ",",header = FALSE)

colnames(raw_data) <- c("motor", "screw", "pgain", "vgain", "class")
out_name <- c("class")
cont_names <- c("pgain", "vgain")
cat_names <- c("motor", "screw")
use_cut_names <- c("pgain", "vgain")

servo_data <- preprocess(raw_data, cont_names, cat_names, out_name, use_cut_names)
save(servo_data, file = "data/servo.RData")
