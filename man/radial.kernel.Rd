\name{radial.kernel}
\alias{radial.kernel}
\alias{poly.kernel}
\title{compute the kernel matrix for svmpath}
\description{
compute the kernel matrix for svmpath
}
\usage{
radial.kernel(x, y=x, param.kernel = 1)
poly.kernel(x, y=x, param.kernel = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an n x p matrix of features}
  \item{y}{an m x p matrix of features (if omitted, it defaults to \code{x})}
  \item{param.kernel}{the parameter(s) for the kernel. For this radial
    kernel, the parameter is known in the fields as "gamma". For the
    polynomial kernel, it is the "degree"}
}
\details{
For the radial kernel, this computes the function \eqn{\exp(-\gamma||x-y||^2)}{exp(-g||x-y||^2)} for each pair of rows x,y
from the input matrices. Here g is param.kernel.
For the polynomial kernel, it computes \eqn{(xy^T+1)^d}{(xy'+1)^d}, where d is \code{param.kernel}.
}
\value{
An n x m matrix.}
\references{The paper \url{http://www-stat.stanford.edu/~hastie/Papers/svmpath.pdf}, as well
as the talk \url{http://www-stat.stanford.edu/~hastie/TALKS/svmpathtalk.pdf}.
}
\author{Trevor Hastie}

\seealso{svmpath}
\examples{
data(svmpath)
attach(balanced.overlap)
fit<-svmpath(x,y,kernel=radial.kernel)
detach(2)
}
\keyword{regression}
