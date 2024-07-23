sim_number <- 1:50
cat_unif <- c(TRUE, FALSE)
dgp <- 1:4
n_all <- c(1000, 5000, 10000)

sim_settings <-
  expand.grid(n_all = n_all,
              dgp = dgp,
              cat_unif = cat_unif,
              sim_number = sim_number)
block_size <- 12
block_starts <- seq(1, nrow(sim_settings), by = block_size)
block_ends <- c(block_starts[-1]-1, nrow(sim_settings))
save(sim_settings, block_size, block_starts, block_ends, 
     file = "sim_settings.RData")