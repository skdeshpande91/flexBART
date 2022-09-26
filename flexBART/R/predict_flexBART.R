predict_flexBART <- function(fit, 
                             X_cont = matrix(0, nrow = 1, ncol = 1), 
                             X_cat = matrix(0, nrow = 1, ncol = 1), 
                             verbose = FALSE, print_every = 50)
{
  # flexBART does not include an is.probit by default so we'll add it
  if(is.null(fit[["is.probit"]])) fit[["is.probit"]] <- FALSE

  tmp <- .predict_flexBART(tree_draws = fit$trees,
                           tX_cont = t(X_cont),
                           tX_cat = t(X_cat),
                           probit = fit[["is.probit"]],
                           verbose = verbose,
                           print_every = print_every)
  if(!fit[["is.probit"]]) output <- fit$y_mean + fit$y_sd * tmp
  else output <- tmp
  return(output)
}