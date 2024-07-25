####################################
# 10 continuous and 10 categorical variables
# f0: 10 * sin(pi * tmp_X[,1] * tmp_X[2])
# f1: 10 ( (tmp_X[,3] - 0.5)^2)
# f2: 10 + 10 * tmp_X[,4] + 5 * tmp_X[,5]
# f3:  12 * tmp_X[,1] + 8 - 20 * (tmp_X[,2] > 0.5)) * sin(pi * tmp_X[,1]) - 8 * (tmp_X[,2] > 0.5))

# DGP1: (0,2,4,8): gets f0+f1+f2 (full friedman function), (1,3,5,6,7,9) get f3
# DGP2: 0 gets f3 and 1--9 get friedman
# DGP3: (0,..., 6) get part of friedman, 7--9 get f3
# DGP4: level l gets (l+1)/10 * f0 + (9-l)/10 * f1

## Update: dgp1 should have offsets and scaling
# DGP3: make 7 be the positive part of the test function, 8 the negative part, and 9 the other one
  # DGP4: should also get the offsets of DGP1

dgp1 <- function(X_cont, X_cat, 
                 offset_fried = 0, 
                 offset_test = 0, scale_test = 1){
  n <- nrow(X_cont)
  ix0 <- which(X_cat[,1] %in% c(0, 2, 4, 8))
  ix1 <- which(X_cat[,1] %in% c(1,3,5,6,7, 9))
  
  mu <- rep(NA, times = n)
  if(length(ix0) > 0){
    mu[ix0] <- 
      friedman_full(X_cont[ix0,], 
                    offset = offset_fried)
  }
  if(length(ix1) > 0){
    mu[ix1] <- test_fun(X_cont[ix1,],
                        offset = offset_test,
                        scale = scale_test)
  }
  return(mu)
}

dgp2 <- function(X_cont, X_cat, 
                 offset_fried = 0,
                 offset_test = 0,
                 scale_test = 1){
  n <- nrow(X_cont)
  ix0 <- which(X_cat[,1] == 0)
  ix1 <- which(X_cat[,1] != 0)
  mu <- rep(NA, times = n)
  if(length(ix0) > 0){
    mu[ix0] <- 
      test_fun(X_cont[ix0,],
               offset = offset_test,
               scale = scale_test)
  }
  if(length(ix1) > 0){
    mu[ix1] <- 
      friedman_full(X_cont[ix1,],
                    offset = offset_fried)
  }
  return(mu)
}

dgp3 <- function(X_cont, X_cat,
                 offset_test = 0,
                 scale_test = 0){
  n <- nrow(X_cont)
  ix0 <- which(X_cat[,1] == 0)
  ix1 <- which(X_cat[,1] == 1)
  ix2 <- which(X_cat[,1] == 2)
  ix3 <- which(X_cat[,1] == 3)
  ix4 <- which(X_cat[,1] == 4)
  ix5 <- which(X_cat[,1] == 5)
  ix6 <- which(X_cat[,1] == 6)
  ix_rem <- which(X_cat[,1] %in% c(7,8,9))
  
  mu <- rep(NA, times = n)
  
  if(length(ix0) > 0){
    mu[ix0] <- 
      friedman_part1(X_cont[ix0,])
  }
  if(length(ix1) > 0){
    mu[ix1] <-
      friedman_part2(X_cont[ix1,])
  }
  if(length(ix2) > 0){
    mu[ix2] <-
      friedman_part3(X_cont[ix2,])
  }
  if(length(ix3) > 0){
    mu[ix3] <-
      friedman_part1(X_cont[ix3,]) + 
      friedman_part2(X_cont[ix3,])
  }
  if(length(ix4) > 0){
    mu[ix4] <-
      friedman_part1(X_cont[ix4,]) + 
      friedman_part3(X_cont[ix4,])
  }
  if(length(ix5) > 0){
    mu[ix5] <-
      friedman_part2(X_cont[ix5,]) + 
      friedman_part3(X_cont[ix5,])
  }
  if(length(ix6) > 0){
    mu[ix6] <-
      friedman_full(X_cont[ix6,])
  }
  if(length(ix_rem) > 0){
    mu[ix_rem] <-
      test_fun(X_cont[ix_rem,],
               offset = offset_test,
               scale = scale_test)
  }
  return(mu)
}
dgp4 <- function(X_cont, X_cat,
                 offset_fried = 0, 
                 offset_test = 0, scale_test = 1){
  tmp_f0 <- friedman_full(X_cont, offset = offset_fried)
  tmp_f1 <- test_fun(X_cont, offset = offset_test, scale = scale_test)
  
  mu <- 
    (X_cat[,1] + 1)/10 * tmp_f0 + 
    (9 - X_cat[,1])/10 * tmp_f1
  return(mu)
}

