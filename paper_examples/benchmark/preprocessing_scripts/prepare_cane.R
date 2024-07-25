source("preprocessing_scripts/preprocess.R")
raw_data <- read.table(file = "raw_data/cane.txt", header = TRUE, stringsAsFactors = TRUE)

# District, DistrictGroup, DistrictPosition, SoilName are categorical
# variety is categorical, ratoon is as well
# outcome is Tonn.Hect
# 

out_name <- c("Tonn.Hect")
cat_names <- c("District", "DistrictGroup", "DistrictPosition", "SoilName", "Variety",
               "HarvestMonth")
cont_names <- c("HarvestDuration", "Age","Fibre",
  "Sugar","Jul.96", "Aug.96", "Sep.96", "Oct.96", "Nov.96", "Dec.96",
                "Jan.97", "Feb.97", "Mar.97", "Apr.97", "May.97", "Jun.97",
                "Jul.97", "Aug.97", "Sep.97", "Oct.97", "Nov.97", "Dec.97")
use_cut_names <- c("Age", "HarvestDuration")

cane_data <- preprocess(raw_data, cont_names, cat_names, out_name, use_cut_names)
save(cane_data, file = "data/cane.RData")
