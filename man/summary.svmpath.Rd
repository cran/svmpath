\name{summary.svmpath}
\alias{summary.svmpath}
\title{produce a summary of an svmpath object}
\description{
printing an svmpath object can produce a lot of lines. The summary
methods gives a more concise description by picking out a subset of the
steps
}
\usage{
\method{summary}{svmpath}(object, nsteps = 5, digits = 6, ...)
}
\arguments{
  \item{object}{the \code{svmpath} object}
  \item{nsteps}{usually omitted, but can be changed to get longer
    summaries}
  \item{digits}{number of significant digits}
  \item{...}{additional arguments to the generic summary function}
}
\details{
Uses the \code{pretty} function to extract the approximately the desired
number of steps. Always includes the first and last step.
}
\value{
returns a dataframe with the steps, value of lambda, training error,
size of elbow, number of support points, and the sum of the overlaps}
\references{The paper \url{http://www-stat.stanford.edu/~hastie/Papers/svmpath.pdf}, as well
as the talk \url{http://www-stat.stanford.edu/~hastie/TALKS/svmpathtalk.pdf}.
}
\author{Trevor Hastie}
\seealso{coef.svmpath, svmpath, predict.svmpath, print.svmpath}
\examples{
data(svmpath)
attach(balanced.overlap)
fit <- svmpath(x,y,trace=TRUE,plot=TRUE)
summary(fit)
detach(2)
}
\keyword{regression}
