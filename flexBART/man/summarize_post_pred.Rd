\name{summarize_post_pred}
\alias{summarize_post_pred}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarize the posterior predictive distribution from a flexBART fit.
}
\description{
Compute the posterior predictive mean and 95\% interval for new observations (i.e. f(x) + noise).
}
\usage{
summarize_post_pred(fit_samples, sigma_samples)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit_samples}{Matrix of posterior samples of evaluations of the regression function (e.g. the \code{y_train} object returned by \code{flexBART} or the matrix returned by \code{predict_flexBART}.)}
  \item{sigma_samples}{Vector of posterior samples of the residual standard deviation sigma. These samples are contained in the \code{sigma} object returned by \code{flexBART}.}
}
\details{

}
\value{
A matrix with three columns, one each for the posterior predictive mean, 2.5\% quantile, and 97.5\% quantiles. The rows correspond to the rows of \code{fit_samples}.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
