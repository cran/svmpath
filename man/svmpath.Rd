\name{svmpath}
\alias{svmpath}
\title{Fit the entire regularization path for a 2-class SVM}
\description{
  The SVM has a regularization or cost parameter C, which controls the
amount by which points overlap their soft margins. Typically either a
default large value for C is chosen (allowing minimal overlap),
or else a few values are compared
using a valication set. This algorithm computes the entire
regularization path (i.e. for all possible values of C for which the
solution changes), with a cost a small (~3)  multiple of the cost of fitting a
single model.
}
\usage{
svmpath(x, y, K, kernel.function = poly.kernel, param.kernel = 1, trace, plot.it, eps = 1e-10, Nmoves = 3 * n, digits = 6, lambda.min = 1e-04, ...)
}
\arguments{
  \item{x}{the data matrix (n x p) with n rows (observations) on p
    variables (columns)}
  \item{y}{The "-1,+1" valued response variable.}
  \item{K}{a n x n kernel matrix, with default value \code{K= kernel.function(x, x)}}
  \item{kernel.function}{This is a user-defined function. Provided are
    \code{poly.kernel} (the default, with parameter set to default to a
    linear kernel) and \code{radial.kernel}}
  \item{param.kernel}{parameter(s) of the kernels}
  \item{trace}{if \code{TRUE}, a progress report is printed as the
    algorithm runs; default is \code{FALSE}}
  \item{plot.it}{a flag indicating whether a plot should be produced
    (default \code{FALSE}; only usable with \code{p=2}}
  \item{eps}{a small machine number which is used to identify minimal
    step sizes}
  \item{Nmoves}{the maximum number of moves}
  \item{digits}{the number of digits in the printout}
  \item{lambda.min}{The smallest value of \code{lambda = 1/C}; default
    is \code{lambda=10e-4}, or \code{C=10000}}
  \item{...}{additional arguments to some of the functions called by
    svmpath. 
  }
}
\details{
  The algorithm used in \code{svmpath()} is described in detail in
  "The Entire Regularization Path for the Support Vector Machine" by
  Hastie, Rosset, Tibshirani and Zhu (2004). It exploits the fact that
  the "hinge" loss-function is piecewise linear, and the penalty term is
  quadratic. This means that in the dual space, the lagrange multipliers
  will be pieceise linear (c.f. \code{lars}).
  }
\value{
a "svmpath" object is returned, for which there are print, summary, coef
and predict methods.}
\references{The paper \url{http://www-stat.stanford.edu/~hastie/Papers/svmpath.pdf}, as well
as the talk \url{http://www-stat.stanford.edu/~hastie/TALKS/svmpathtalk.pdf}.
}
\author{Trevor Hastie}
\note{This implementation of the algorithm does not use updating to solve the "elbow"
  linear equations. This is possible, since the elbow changes by a small
  number of points at a time. Future version of the software will do
  this. The author has encountered numerical problems with early
  attempts at this.}
\section{Warning}{Currently the algorithm can get into machine errors if
  \code{epsilon} is too small, or if \code{lambda.min} is too
  small. Increasing either from their defaults should make the problems
  go away, by terminating the algorithm slightly early.}
\seealso{print, coef, summary, predict, and FilmPath}
\examples{
data(svmpath)
attach(unbalanced.separated)
svmpath(x,y,trace=TRUE,plot=TRUE)
detach(2)
\dontrun{svmpath(x,y,kernel=radial.kernel,param.kernel=.8)}
}
\keyword{regression}
