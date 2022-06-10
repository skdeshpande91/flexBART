\name{predict_flexBART}
\alias{predict_flexBART}
\title{
Predicting new observations with previously fitted flexBART model
}
\description{
\code{predict_flexBART} can take the output of \code{flexBART} and use it to make predictions at new inputs.
}
\usage{
predict_flexBART(fit, 
                 X_cont = matrix(0, nrow = 1, ncol = 1), 
                 X_cat = matrix(0, nrow = 1, ncol = 1), 
                 verbose = FALSE, print_every = 50)
}

\arguments{
  \item{fit}{Object returned by \code{flexBART}}
  \item{X_cont}{Matrix of continuous predictors. Note all predictors must be re-scaled to lie in the interval [-1,1]. Default is a 1x1 matrix, which signals that no continuos predictors are available.}
  \item{X_cat}{Integer matrix of categorical predictors for training data. Note categorical levels should be 0-indexed. That is, if a categorical predictor has 10 levels, the values should run from 0 to 9. Default is a 1x1 matrix, which signals that no categorical predictors.}
  \item{verbose}{Logical, indicating whether or not to print message predictions are being made. Default is \code{FALSE}.}
  \item{print_every}{As the function loops over the MCMC samples, a message is printed to the console every \code{print_every} iterations. Default is 50.}
%%     ~~Describe \code{x} here~~
}
\details{
Make predictions at new inputs based on the output of \code{flexBART}. Useful when the testing dataset is quite large.
If \code{fit} were produced by \code{predict_flexBART} then the function outputs draws of the fitted probabilities.
}
\value{
A matrix containing posterior samples of the regression function evaluated at the supplied inputs.
The rows of the matrix correspond to MCMC iterations and the columns correspond to the observations in the supplied data (i.e. rows of \code{X_cont} and/or \code{X_cat}).
}