load("philly_crime_data.RData")

method_list <-
  c(paste0("networkBART_", 1:7),
    "flexBART",
    "BART", "DART",
    "BART_latlon", 
    paste0("BART_ase", c(1,3,5)),
    "targetBART")


sim_number <- 1:100

sim_settings <- expand.grid(method = method_list, sim_number = sim_number,
                            stringsAsFactors = FALSE)

block_size <- 2
block_starts <- seq(1, nrow(sim_settings), by = block_size)
block_ends <- c(block_starts[-1]-1, nrow(sim_settings))


save(sim_settings, method_list, 
     block_size, block_starts, block_ends,file = "sim_settings.RData")
