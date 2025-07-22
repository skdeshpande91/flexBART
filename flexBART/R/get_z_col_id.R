get_z_col_id <- function(Z_train, tol = 1e-12){
  # id_z_map has length R; values are from 0 to num_unik_z
  R <- ncol(Z_train)
  z_counter <- 0
  id_z_map <- rep(NA, times = R)
  id_z_map[1] <- 0 # for intercept 
  if(R > 1){
    for(r in 2:R){
      flag <- FALSE
      for(rr in 1:(r-1)){
        if(max(abs(Z_train[,r] - Z_train[,rr])) < tol){
          # Z_train[,r] is identical to an earlier column!
          id_z_map[r] <- id_z_map[rr]
          flag <- TRUE
          break
        }
      }
      if(!flag){
        # Z_train[,r] is unik
        z_counter <- z_counter+1 # increment counter
        id_z_map[r] <- z_counter
      }
    }
  }
  return(id_z_map)
}
