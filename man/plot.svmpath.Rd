\name{plot.svmpath}
\alias{plot.svmpath}
\title{plot the svm solution at a step along the path}
\description{
produces a plot of the svm solution along the path, and optinally
indicates support points}
\usage{
\method{plot}{svmpath}(x, step, Size = 60, elbow.show = TRUE, support.show = TRUE, ...)
}
\arguments{
  \item{x}{the \code{svmpath} object}
  \item{step}{which step to plot; default is the last step. Use
    \code{summary} to see how many steps}
  \item{Size}{If the solution is non-linear, this is the gridsize for \code{countour}}
  \item{elbow.show}{Should the points on the elbow be indicated}
  \item{support.show}{Should the support points be indicated}
  \item{\dots}{additional arguments to plot, allowing one to change, for
    example, "main", "xlab" etc}
}
\details{
A two-dimensional plot is produced of the SVM solution. Makes sense only
if X is two-dimensional. If not, the first two dimensions will be used}
\value{
A list is returned silently, with the ingredients of the plot
}
\references{The paper \url{http://www-stat.stanford.edu/~hastie/Papers/svmpath.pdf}, as well
as the talk \url{http://www-stat.stanford.edu/~hastie/TALKS/svmpathtalk.pdf}.
}
\author{Trevor Hastie}
\seealso{coef.svmpath, svmpath, predict.svmpath, print.svmpath,summary.svmpath}
\examples{
data(svmpath)
attach(balanced.overlap)
fit <- svmpath(x,y,trace=TRUE,plot=FALSE)
plot(fit,step=2)
detach(2)
}
\keyword{regression}
