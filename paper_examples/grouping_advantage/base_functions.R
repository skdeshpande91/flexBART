friedman_part1 <- function(X_cont){
  stopifnot(ncol(X_cont) == 10)
  stopifnot(all(abs(X_cont) <= 1))
  tmp_X <- (X_cont+1)/2 # rescale to [0,1]
  return(10 * sin(pi*tmp_X[,1]*tmp_X[,2])) 
}
friedman_part2 <- function(X_cont){
  stopifnot(ncol(X_cont) == 10)
  stopifnot(all(abs(X_cont) <= 1))
  tmp_X <- (X_cont+1)/2 # rescale to [0,1]
  return(10 * (tmp_X[,3] - 0.5)^2)
}
friedman_part3 <- function(X_cont){
  stopifnot(ncol(X_cont) == 10)
  stopifnot(all(abs(X_cont) <= 1))
  tmp_X <- (X_cont+1)/2 # rescale to [0,1]
  return(10 * (tmp_X[,3] - 0.5)^2 + 
           10*tmp_X[,4] + 5 * tmp_X[,5])
}
friedman_full <- function(X_cont, offset = 0){
  f1 <- friedman_part1(X_cont)
  f2 <- friedman_part2(X_cont)
  f3 <- friedman_part3(X_cont)
  return(f1 + f2 + f3 + offset)
}
test_fun <- function(X_cont, offset  = 0, scale = 1){
  stopifnot(ncol(X_cont) == 10)
  stopifnot(all(abs(X_cont) <= 1))
  tmp_X <- (X_cont + 1)/2
  
  tmp <- 3 * tmp_X[,1] + 
    (2 - 5 * (tmp_X[,2] > 0.5)) * sin(pi * tmp_X[,1]) - 
    2 * (tmp_X[,2] > 0.5)
  
  return(offset + scale * tmp)
}