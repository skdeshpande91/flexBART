\name{predict.flexBART}
\alias{predict.flexBART}
\title{
Predicting new observations with previously fitted flexBART model
}
\description{
\code{predict.flexBART()} can take the output of \code{\link{flexBART()}} and \code{\link{probit_flexBART()}} and use it to make predictions at new inputs.
}
\usage{
\method{predict}{flexBART}(object, newdata, ...)
}

\arguments{
  \item{object}{object of class inhereting from \dQuote{flexBART}.}
  \item{newdata}{Data frame in whch to look for variables with which to predict. Cannot be omitted.}
  \item{\dots}{Additional optional arguments governing whether to print progress (i.e., \code{verbose} and \code{print_every}).}
}
\details{
Make predictions at new inputs based on the output of \code{\link{flexBART()}}. 
Useful when the testing dataset is quite large.
If \code{fit} were produced by \code{\link{probit_flexBART()}}, then the function outputs draws of the fitted probabilities.
}
\value{
When there is only one ensemble, a matrix containing posterior samples of the regression function evaluated at the supplied inputs.
For models with multiple ensembles, a list with three elements: (i) a matrix containing posterior samples of the regression function evaluation and 
(ii) an array containing evaluations of the identified slopes; and (iii) an array containing evaluations of all slopes on the standardized scale.
}