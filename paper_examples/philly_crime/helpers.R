summarize_post_pred <- function(mu_samples, sigma_samples){
  nd <- nrow(mu_samples)
  if(length(sigma_samples) != nd) stop("mu_samples & sigma_samples must be the same length")
  
  n_obs <- ncol(mu_samples)
  post_summary <- 
    data.frame(MEAN = colMeans(mu_samples),
               L95 = rep(NA, times = n_obs),
               U95 = rep(NA, times = n_obs))
  for(i in 1:n_obs){
    ystar <- mu_samples[,i] + sigma_samples * rnorm(n = nd, mean = 0, sd = 1)
    post_summary[i,c("L95", "U95")] <- quantile(ystar, prob = c(0.025, 0.975))
  }
  return(post_summary)
}

interval_score <- function(truth, lower, upper, alpha = 0.05){
  tmp_score <- (upper - lower) + 
    2/alpha * (lower - truth) * (truth < lower) +
    2/alpha * (truth - upper) * (truth > upper)
  return(mean(tmp_score))
}

coverage <- function(truth, lower, upper){
  return(mean( (truth >= lower) & (truth <= upper)))
}
