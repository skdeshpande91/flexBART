library(Rcpp)
library(igraph)
library(RcppArmadillo)

sourceCpp("flexBART/src/rflexBART.cpp")

p_cont <- 10
p_cat <- 10

p <- p_cont + p_cat

n <- 10000
set.seed(129)

# First five of the continuous variables we will supply cutpoints
# second half will be drawn uniformly (so second half needs to be scaled to [-1,1])
X_cont <- matrix(nrow = n, ncol = p_cont)
cutpoints_list <- list()
unif_cuts <- c(rep(FALSE, times = 5), rep(TRUE, times = p_cont - 5))
for(j in 1:5){
  X_cont[,j] <- rnorm(n, mean = 0, sd = 5)
  sd_x_j <- sd(X_cont[,j])
  cutpoints_list[[j]] <- seq(min(X_cont[,j]) - sd_x_j, max(X_cont[,j]) + sd_x_j, length = 1000)
}
for(j in 6:p_cont){
  X_cont[,j] <- runif(n, min = -1, max = 1)
  cutpoints_list[[j]] <- c(0)
}

# first 3 categorical variables will have some graphical structure
g1 <- erdos.renyi.game(10, 0.75, type = "gnp")
g2 <- erdos.renyi.game(10, 0.5, type = "gnp")
g3 <- erdos.renyi.game(10, 0.25, type = "gnp")





X_cat <- matrix(NA, nrow = n, ncol = p_cat)
cat_levels_list <- list()
adj_support_list <- list()
graph_split <- rep(FALSE, times = p_cat)
for(j in 1:p_cat){
  X_cat[,j] <- as.integer(sample(0:9), size = n, replace = TRUE)
  cat_levels_list[[j]] <- 0:9
  if(j <= 3){
    graph_split[j] <- TRUE
    g <- get(paste0("g",j))
    A <- as_adjacency_matrix(g, sparse = FALSE, names = FALSE)
    A_lower <- A
    A_lower[upper.tri(A_lower)] <- 0
    adj_support_list[[j]] <- which(A_lower != 0) - 1 # C++ is 1-indexed
  } else{
    adj_support_list[[j]] <- c(0) # this corresponds to an invalid adjacency matrix
  }
}


test <- .draw_tree(tX_cont = t(X_cont),
                   tX_cat = t(X_cat),
                   unif_cuts = unif_cuts,
                   cutpoints_list = cutpoints_list,
                   cat_levels_list = cat_levels_list,
                   graph_split = graph_split, graph_cut_type = 1,
                   adj_support_list = adj_support_list, 
                   rc_split = FALSE, prob_rc = 0.0,
                   alpha = 0.95, beta = 2.0,
                   mu0 = 0.0, tau = 1.0)
                   

# let's use fixed 