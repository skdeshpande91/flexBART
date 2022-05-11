predict_flexBART <- function(fit, X_cont, X_cat, verbose, print_every)
{
  tmp <- .predict_flexBART(tree_draws = fit$trees,
                           tX_cont = t(X_cont),
                           tX_cat = t(X_cat),
                           verbose = verbose,
                           print_every = print_every)
  output <- fit$y_mean + fit$y_sd * tmp
  return(output)
}