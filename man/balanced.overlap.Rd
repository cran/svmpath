\name{balanced.overlap}
\alias{balanced.overlap}
\alias{balanced.separated}
\alias{unbalanced.separated}
\alias{mixture.data}
\docType{data}
\title{simple examples for svmpath}
\description{
Datasets for illustrating the svmpath function, that can be plotted
while its running
}
\usage{data(svmpath)}
\format{
In each case a list with a component \code{x} (t column matrix) and a
component \code{y} (vector of +1/-1 values)
"Balanced" refers to whether the number of +1s is the same as the -1s.
"Overlap" indicates whether the classes are linearly separable.
\code{mixture.data} is a balanced dataset with 100 observations in
each class. The others are smaller with between 10-12 obs total.}
\references{The paper \url{http://www-stat.stanford.edu/~hastie/Papers/svmpath.pdf}, as well
as the talk \url{http://www-stat.stanford.edu/~hastie/TALKS/svmpathtalk.pdf}.
}
\examples{
data(svmpath)
attach(balanced.overlap)
svmpath(x,y,trace=TRUE,plot=TRUE)
detach(2)
}
\keyword{datasets}
