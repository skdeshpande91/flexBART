library(igraph)
n_side <- 10
n <- 100

g <- make_lattice(length = n_side, dim = 2)
A <- as_adjacency_matrix(g, type = "both", sparse = FALSE)
colnames(A) <- 0:(n-1)
rownames(A) <- 0:(n-1)

col_list <- colorBlindness::Blue2DarkRed18Steps
# make up 5 clusters
cluster1 <- 81:100
cluster1 <- cluster1[!cluster1 %in% c(89,90, 100)]
cluster2 <- c(89, 90, 100, 79, 80, 69, 70, 59, 60, 50, 40)
cluster3 <- rep(1:4, times = 4) + 10 * rep(0:3, each = 4) # lower 4x4 grid
cluster4 <- c(52, 53,63,54, 64, 55,65,56, 66,46, 
              47, 57, 67, 37, 38, 48, 28, 29, 27, 
              17, 16, 18, 6, 7, 8,30, 20, 10, 9, 19, 39, 49)
cluster5 <- which(!1:100 %in% c(cluster1, cluster2, cluster3, cluster4))

cluster_means <- c(7.5, 3, -2, -5, 0)
mu <- rep(NA, times = n)
mu[cluster1] <- cluster_means[1]
mu[cluster2] <- cluster_means[2]
mu[cluster3] <- cluster_means[3]
mu[cluster4] <- cluster_means[4]
mu[cluster5] <- cluster_means[5]

sigma <- 1
t <- 10
set.seed(624)

full_data <- 
  data.frame(X1 = factor(rep(0:(n-1), each = t), levels = 0:(n-1)))
full_data[,"Y"] <- rep(mu, each = t) + sigma * rnorm(n*t, mean = 0, sd = 1)

set.seed(129)
vertex_id_all <- rep(1:n, each = t)
test_vertices <- sample(1:n, size = floor(0.1 * n), replace = FALSE)
train_vertices <- (1:n)[-test_vertices]

train_index <- which(vertex_id_all %in% train_vertices)

train_data <- full_data[train_index,]

test_data <- data.frame(X1 = factor(0:(n-1), levels = 0:(n-1)))
adjacency_list <- list(X1 = A)



