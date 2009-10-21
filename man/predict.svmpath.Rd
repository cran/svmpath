\name{predict.svmpath}
\alias{predict.svmpath}
\title{Make predictions from a "svmpath" object}
\description{
  Provide a value for \code{lambda}, and produce the fitted lagrange \code{alpha}
  values. Provide values for \code{x}, and get fitted function values or
  class labels.
}
\usage{
\method{predict}{svmpath}(object, newx, lambda, type = c("function", "class",
"alpha", "margin"),...)
}
\arguments{
  \item{object}{fitted \code{svmpath} object}
  \item{newx}{values of \code{x} at which prediction are wanted. This is a
    matrix with observations per row}
  \item{lambda}{the value of the regularization parameter. Note that
    \code{lambda} is equivalent to 1/C for the usual parametrization of
    a SVM}
  \item{type}{type of prediction, with default \code{"function"}. For
    \code{type="alpha"} or \code{type="margin"} the  \code{newx} argument is not required}
  \item{...}{Generic compatibility}
}
\details{
This implementation of the SVM uses a parameterization that is slightly
different but equivalent to the usual (Vapnik) SVM. Here
\eqn{\lambda=1/C}{lambda=1/C}.
The Lagrange multipliers are related via
\eqn{\alpha^*_i=\alpha_i/\lambda}{\alphastar = alpha/lambda}, where
\eqn{\alpha^*_i}{alphastar} is the usual multiplier, and
\eqn{\alpha_i}{alpha} our multiplier. Note that if \code{alpha=0}, that
observation is right of the elbow; \code{alpha=1}, left of the elbow;
\code{0<alpha<1} on the elbow. The latter two cases are all support
points.
}
\value{
  In each case, the desired prediction.
}
\references{The paper \url{http://www-stat.stanford.edu/~hastie/Papers/svmpath.pdf}, as well
as the talk \url{http://www-stat.stanford.edu/~hastie/TALKS/svmpathtalk.pdf}.
}
\author{Trevor Hastie}

\seealso{coef.svmpath, svmpath}
\examples{
data(svmpath)
attach(balanced.overlap)
fit <- svmpath(x,y,trace=TRUE,plot=TRUE)
predict(fit, lambda=1,type="alpha")
predict(fit, x, lambda=.9)
detach(2)
}
\keyword{regression}
