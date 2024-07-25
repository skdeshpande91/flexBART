data_list <- c("abalone", "ais", "Alcohol",
               "ammenity", "attend", "cane",
               "Caschool", "cpu", "Insur","Medicare", 
               "mpg", "servo", 
               "spouse", "strike",
               "engel", "fuelEcon")
sim_number <- 1:50

sim_settings <- 
  expand.grid(sim_number = sim_number,
              dataset = data_list,
              stringsAsFactors = FALSE)

block_starts1 <- seq(from = 1, to = 700, by = 10)
block_ends1 <- c(block_starts1[-1]-1,700)
block_starts2 <- seq(from = 701, to = 800, by = 2)
block_ends2 <- c(block_starts2[-1]-1, 800)

block_starts <- c(block_starts1, block_starts2)
block_ends <- c(block_ends1, block_ends2)

save(sim_settings, block_starts, block_ends, data_list, file = "sim_settings.RData")
