\name{print.svmpath}
\alias{print.svmpath}
\title{Print a summary of the SVM path}
\description{
print a summary of the fitted svmpath object
}
\usage{
print.svmpath(x, digits, ...)
}
\arguments{
  \item{x}{object to be printed}
  \item{digits}{number of significant digits (default 6)}
  \item{...}{additional arguments to the generic print function}
  }
\value{
  For each step taken by the algorithm, one or more lines are printed.
  The step is described in terms of the observation number involved, a
  coded version of what happened, such as "L->E" meaning "from the Left
  set" to the "Elbow". Initially all the sets are empty. It gives the
  margin (sum of the xi), the size of the elbow, and the training error.}
\references{The paper \url{http://www-stat.stanford.edu/~hastie/Papers/svmpath.pdf}, as well
as the talk \url{http://www-stat.stanford.edu/~hastie/TALKS/svmpathtalk.pdf}.
}
\author{Trevor Hastie}

\seealso{coef.svmpath, svmpath, predict.svmpath}
\examples{
data(svmpath)
attach(balanced.overlap)
fit <- svmpath(x,y,trace=TRUE,plot=TRUE)
print(fit)
detach(2)
}
\keyword{regression}
