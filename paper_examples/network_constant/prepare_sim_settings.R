networkBART_grid <-
  expand.grid(alg = "networkBART",
              graph_cut_type = 1:7,
              sparse = FALSE,
              embed = NA,
              embed_dim = NA,
              sim_number = 1:50,
              n_obs = c(10, 100))
flexBART_grid <-
  expand.grid(alg = "flexBART",
              graph_cut_type = NA,
              sparse = FALSE,
              embed = NA,
              embed_dim = NA,
              sim_number = 1:50,
              n_obs = c(10, 100))
BART_grid <-
  expand.grid(alg = "BART",
              graph_cut_type = NA,
              sparse = c(TRUE, FALSE),
              embed = c(NA, "ase"),
              embed_dim = c(1, 3, 5),
              sim_number = 1:50,
              n_obs = c(10, 100))

exclude_ix <- 
  which(
    (is.na(BART_grid$embed) & BART_grid$embed_dim %in% c(3,5)) |
      (!is.na(BART_grid$embed) & BART_grid$embed_dim == 1 & BART_grid$sparse == TRUE))

exclude_ix <- unique(exclude_ix)
BART_grid <- BART_grid[-exclude_ix,]
BART_grid[which(is.na(BART_grid$embed)),"embed_dim"] <- NA

sim_settings <-
  rbind(networkBART_grid,
        flexBART_grid,
        BART_grid)

block_size <- 10
block_starts <- seq(1, nrow(sim_settings), by = block_size)
block_ends <- c(block_starts[-1]-1, nrow(sim_settings))

save(sim_settings, block_starts, block_ends, file = "sim_settings.RData")