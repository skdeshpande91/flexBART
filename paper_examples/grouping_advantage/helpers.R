summarize <- function(mu_samples, 
                      sigma_samples = NULL){
  nd <- nrow(mu_samples)
  if(!is.null(sigma_samples)){
    if(length(sigma_samples) != nd) stop("mu_samples & sigma_samples must be the same length!")
    post_summary <-
      flexBART::summarize_post_pred(mu_samples, sigma_samples)
  } else{
    post_summary <-
      flexBART::summarize_post_pred(mu_samples, rep(0, times = nd))
  }
  colnames(post_summary) <- c("MEAN", "L95", "U95")
  return(post_summary)
}
square_error <- function(truth, pred){
  return(mean( (truth - pred)^2 ))
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
