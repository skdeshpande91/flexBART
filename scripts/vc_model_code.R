

#X_cont <- matrix(runif(n*2, min = -1, max = 1), nrow = n, ncol = 2)
X_cont <- matrix(runif(N*2, min = -1, max = 1), nrow = N, ncol = 2)
X_cat <- matrix(rep( (1:n) - 1, each = T), nrow = N, ncol = 1)

beta0_true <- function(X_cont){
  scaled_X_cont <- (X_cont + 1)/2 # moves it to [0,1]
  z1 <- scaled_X_cont[,1]
  #z2 <- 1*(X_cont[,2] >= 0)
  z2 <- X_cont[,2]
  return(3 * z1 + (2 - 5 * z2) * sin(pi * z1) - 2 * z2)
}
beta1_true <- function(X_cont){
  scaled_X_cont <- (X_cont + 1)/2
  return( (3 - 3*cos(6*pi*scaled_X_cont[,1]) * scaled_X_cont[,1]^2) * (scaled_X_cont[,1] > 0.6) - (10 * sqrt(scaled_X_cont[,1])) * (scaled_X_cont[,1] < 0.25) )
}
beta2_true <- function(X_cont){
  return(X_cont[,1]^2)
}

beta0 <- beta0_true(X_cont)
beta1 <- beta1_true(X_cont)
beta2 <- beta2_true(X_cont)

mu <- 3*sqrt(n) * L_eig$vector[,n-1]
#mu <- sqrt(n) * (L_eig$vector[,n-1] * beta0 + L_eig$vector[,n-2] * beta1 + L_eig$vector[n-3] * beta2 )
#mu <- beta0 # this works well
mu <- beta0 + sqrt(n) * L_eig$vector[,n-1] # this also seems to work
#mu <- sqrt(n) * beta0 + sqrt(n) * L_eig$vector[,n-1] # also works
# mu <- sqrt(n) * beta0 * L_eig$vector[,n-1] # doesn't work well
#mu <- beta0 + sqrt(n) * X_cont[,1]^2 * L_eig$vectors[,n-1] # worked!
#mu <- beta0 + sqrt(n) * beta1 * L_eig$vectors[,n-1] # did not work that well.

#save(g, A, n, beta0, beta1, beta2, L_eig, fixed_layout, file = "data/regression_data.RData")