source("preprocessing_scripts/preprocess.R")
library(UsingR)
library(forcats)
data("MLBattend")

# Brewers moved from AL to NL but didn't physically move
# combine MILA and MILN categories
# Also CAL and ANA are the same team

MLBattend$franchise <- 
  fct_collapse(MLBattend$franchise, CAL=c("CAL", "ANA"),
               MIL=c("MILA", "MILN"))
MLBattend$year[which(MLBattend$year == 0)] <- 100
MLBattend$year <- MLBattend$year + 1900

out_name <- c("attendance")
cont_names <- c("year", "runs.scored", "runs.allowed", "wins", "losses", "games.behind")
cat_names <- c("franchise", "league", "division")
use_cut_names <- c("year", "runs.scored", "runs.allowed", "wins", "losses", "games.behind")

attend_data <- preprocess(raw_data = MLBattend,
                          cont_names = cont_names,
                          cat_names = cat_names,
                          out_name = out_name,
                          use_cut_names = use_cut_names,
                          log_y = TRUE)
attend_data$cutpoints_list[["year"]] <- 1969:2000
attend_data$cutpoints_list[["runs.scored"]] <- 325:1125
attend_data$cutpoints_list[["runs.allowed"]] <- 325:1125
attend_data$cutpoints_list[["wins"]] <- 0:200
attend_data$cutpoints_list[["losses"]] <- 0:60

save(attend_data, file = "data/attend.RData")
