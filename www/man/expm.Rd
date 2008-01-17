\name{expm}
\alias{expm}
\title{Matrix Exponential}
\description{
  Compute the exponential of a matrix.
}
\usage{
expm(x)
}
\arguments{
  \item{x}{a square matrix.}
}
\details{
  The exponential of a matrix is defined as the infinite Taylor series
  \deqn{e^M = \sum_{k = 1}^\infty \frac{M^k}{k!}.}{%
    exp(M) = I + M + M^2/2! + M^3/3! + \dots}

  The function uses Ward's diagonal \enc{Padé}{Pade} approximation with
  three step preconditioning.
}
\value{
  The matrix exponential of \code{x}.
}
\seealso{
  The package vignette for details on the algorithms and calling the
  function from external packages.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca}, and Christophe
  Dutang, based on code translated by D. Bates and M. Maechler from the
  implementation of the corresponding Octave function contributed by A.
  Scottedward Hodel \email{A.S.Hodel@eng.auburn.edu}.
}
\examples{
x <- matrix(c(-49, -64, 24, 31), 2, 2)
expm(x)
}
\keyword{algebra}
\keyword{math}
\encoding{UTF-8}