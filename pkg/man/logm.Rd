\name{logm}
\alias{logm}
\title{Matrix Logarithm}
\description{
  This function computes the exponential of a square matrix
  \eqn{A} is \eqn{L} such that \eqn{A=exp(L)} where the matrix exponential is 
  defined as the sum from \eqn{r=0} to infinity of
  \eqn{L^r/r!}.
}
\usage{
logm(x, method = c("Eigen"),
		 order = 8,
                 trySym = TRUE, tol = .Machine$double.eps)

}
\arguments{
  \item{x}{a square matrix.}
  \item{method}{
    \code{"Eigen"} tries to diagonalise the matrix \code{x}, if not possible, 
    \code{"Eigen"} raises an error. 
    }
  \item{order}{TODO
  }
  \item{trySym}{TODO}
  \item{tol}{a given tolerance used to check if \code{x} is
    computationally singular when \code{method = "Eigen"}.}
}
\details{
  The exponential of a matrix is defined as the infinite Taylor series
  \deqn{e^M = \sum_{k = 1}^\infty \frac{M^k}{k!}.}{%
    exp(M) = I + M + M^2/2! + M^3/3! + \dots}. The logarithm matrix of \eqn{A} is 
    the matrix \eqn{M} such that \eqn{exp(M) = A}.

}
\value{
  The matrix logarithm of \code{x}.
}
\seealso{
  The package vignette for details on the algorithms and calling the
  function from external packages.
}
\author{
  The \code{"Eigen"} method by Christophe Dutang.
}
\references{
}
\examples{
m <- diag(2)
logm(m)
expm(logm(m))
}
\keyword{algebra}
\keyword{math}
\encoding{UTF-8}
