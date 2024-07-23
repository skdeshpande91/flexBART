method_list <- c("flexBART",
                 "BART",
                 "flexBART_loc",
                 "targetBART")
year <- 2013:2019
sim_number <- 1:10

sim_settings <- 
  expand.grid(method = method_list,
              year = year,
              sim_number = sim_number,
              stringsAsFactors = FALSE)

compute_misclass <- function(phat, Y){
  return( mean(Y != (phat >= 0.5)) )
}

# Compute the brief score
compute_brier <- function(phat, Y){
  return(mean( (Y-phat)^2 ))
}

# Compute the log-loss
compute_logloss <- function(phat,Y, trunc = 1e-60){
  lo_index <- which(phat <= trunc)
  hi_index <- which(phat >= 1-trunc)
  if(length(lo_index) > 0)   phat[lo_index] <- trunc
  if(length(hi_index)) phat[hi_index] <- 1.0 - trunc
  return(-1 * mean( (Y * log(phat) + (1 - Y)*log(1-phat))))
}

save(sim_settings,
     compute_misclass, compute_brier, compute_logloss,
     file = "sim_settings.RData")
