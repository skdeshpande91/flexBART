# Data downloaded from https://www.fueleconomy.gov/feg/download.shtml
source("preprocessing_scripts/preprocess.R")
library(dplyr)

raw_data <- 
  readr::read_csv(file = "raw_data/vehicles.csv")
fuel_econ_data <-
  raw_data %>%
  mutate(volume_hatch = hlv + hpv,
         volume_2 = lv2 + pv2,
         volume_4 = lv4 + pv4) %>%
  filter(!(volume_hatch == 0 & volume_2 == 0 & volume_4 == 0)) %>%
  filter((volume_hatch!=0) + (volume_2!=0) + (volume_4!=0) == 1) %>%
  mutate(volume = volume_hatch + volume_2 + volume_4) %>%
  select(comb08, cylinders, displ, drive,
         fuelCost08, fuelType1, volume, make, year, trany, VClass) %>%
  filter(!fuelType1 %in% c("Electricity", "Hydrogen")) %>%
  filter(!is.na(drive))
  
fuel_econ_df <- as.data.frame(fuel_econ_data)


out_name <- c("comb08")
cont_names <- c("cylinders", "displ", "fuelCost08", "volume", "year")
use_cut_names <- c("cylinders", "year")
cat_names <- c("drive", "fuelType1", "make", "trany", "VClass")

fuelEcon_data <- 
  preprocess(raw_data = fuel_econ_df,
             cont_names = cont_names,
             cat_names = cat_names,
             out_name = out_name,
             use_cut_names = use_cut_names)

save(fuelEcon_data, file = "data/fuelEcon.RData")
