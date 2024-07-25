library(GLMsData)
source("preprocessing_scripts/preprocess.R")

data(motorins)

raw_data <- as.data.frame(motorins)
raw_data[,"Zone"] <- factor(motorins$Zone, levels = 1:7)
raw_data[,"Make"] <- factor(motorins$Make, levels = 1:9)
raw_data[,"Payment"] <- log(motorins$Payment+1)
out_name <- c("Payment")
cont_names <- c("Kilometres", "Bonus", "Insured", "Claims")
cat_names <- c("Zone", "Make")
use_cut_names <- c("Kilometres", "Bonus", "Claims")

Insur_data <- preprocess(raw_data = raw_data,
                         cont_names = cont_names,
                         cat_names = cat_names,
                         out_name = out_name,
                         use_cut_names = use_cut_names,
                         log_y = FALSE)

Insur_data$cutpoints_list[["Claims"]] <- min(motorins$Claims):max(motorins$Claims)
save(Insur_data, file = "data/Insur.RData")
