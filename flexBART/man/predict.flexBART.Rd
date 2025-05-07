\name{predict.flexBART}
\alias{predict.flexBART}
\title{
Predicting new observations with previously fitted flexBART model
}
\description{
\code{predict.flexBART} can take the output of \code{flexBART} and \code{probit_flexBART} and use it to make predictions at new inputs.
}
\usage{
\method{predict}{flexBART}(object, ...)
}

\arguments{
  \item{object}{object of class inhereting from "flexBART."}
  \item{\dots}{additional arguments affecting the predictions produced.}
  %\item{newdata}{Data frame in whch to look for variables with which to predict. Cannot be omitted.}
  %\item{verbose}{Logical, indicating whether or not to print message predictions are being made. Default is \code{FALSE}.}
}
\details{
Make predictions at new inputs based on the output of \code{flexBART}. Useful when the testing dataset is quite large.
If \code{fit} were produced by \code{predict_flexBART} then the function outputs draws of the fitted probabilities.
}
\value{
When there is only one ensemble, a matrix containing posterior samples of the regression function evaluated at the supplied inputs.
For models with multiple ensembles, a list with three elements: (i) a matrix containing posterior samples of the regression function evaluation and 
(ii) an array containing evaluations of the identified slopes; and (iii) an array containing evaluations of all slopes (see details)
}