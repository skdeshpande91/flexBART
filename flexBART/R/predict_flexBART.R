predict_flexBART <- function(fit, 
                             X_cont = matrix(0, nrow = 1, ncol = 1), 
                             X_cat = matrix(0, nrow = 1, ncol = 1), 
                             verbose = FALSE, print_every = 50)
{
  tmp <- .predict_flexBART(tree_draws = fit$trees,
                           tX_cont = t(X_cont),
                           tX_cat = t(X_cat),
                           verbose = verbose,
                           print_every = print_every)
  output <- fit$y_mean + fit$y_sd * tmp
  return(output)
}