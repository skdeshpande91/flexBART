source("preprocessing_scripts/preprocess.R")

raw_data <- read.table("raw_data/strikes.data", skip = 19, header = FALSE)
colnames(raw_data) <- c("country", "year", "strike_volume",
                        "unemployment", "inflation", "socdem_lbr","union")

raw_data[,"country"] <- factor(raw_data[,"country"], levels = 1:18, labels = 1:18)

out_name <- "strike_volume"
cont_names <- c("year", "unemployment", "inflation", "socdem_lbr", "union")
cat_names <- c("country")
use_cut_names <- c("year")


strike_data <- preprocess(raw_data, cont_names, cat_names, out_name, use_cut_names, log_y = TRUE)
save(strike_data, file = "data/strike.RData")
