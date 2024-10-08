\name{expm.Higham08}
\Rdversion{1.1}
\alias{expm.Higham08}
\title{Matrix Exponential [Higham 2008]}
\description{
  Calculation of matrix exponential \eqn{e^A} with the \sQuote{Scaling &
    Squaring} method with balancing.

  Implementation of Higham's Algorithm from his book (see references),
  Chapter 10, Algorithm 10.20.

  The balancing option is an extra from Michael Stadelmann's Masters thesis.
}
\usage{
expm.Higham08(A, balancing = TRUE)
}
\arguments{
  \item{A}{square matrix, may be a \code{"\link[Matrix:sparseMatrix-class]{sparseMatrix}"},
    currently only if \code{balancing} is false.}
  \item{balancing}{logical indicating if balancing should happen (before
    and after scaling and squaring).}
}
\details{
  The algorithm comprises the following steps
  \describe{
    \item{0.}{Balancing}
    \item{1.}{Scaling}
    \item{2.}{Padé-Approximation}
    \item{3.}{Squaring}
    \item{4.}{Reverse Balancing}
  }
}
\value{
  a matrix of the same dimension as \code{A}, the matrix exponential of \code{A}.
}
\references{
  Higham, Nicholas J. (2008).
  \emph{Functions of Matrices: Theory and Computation};
  SIAM (Society for Industrial and Applied Mathematics), Philadelphia, USA;
  \doi{10.1137/1.9780898717778}

  Michael Stadelmann (2009).
  \emph{Matrixfunktionen; Analyse und Implementierung}.
  [in German] Master's thesis and Research Report 2009-12, SAM, ETH Zurich;
  \url{https://math.ethz.ch/sam/research/reports.html?year=2009}, or the
  pdf directly at
  \url{https://www.sam.math.ethz.ch/sam_reports/reports_final/reports2009/2009-12.pdf}.
}
\author{
  Michael Stadelmann (final polish by Martin Maechler).
}
\note{
  \code{expm.Higham8()} no longer needs to be called directly; rather
  \code{expm(A, "Higham8b")} and \code{expm(A, "Higham8")} correspond to
  the two options of \code{balancing = TRUE || FALSE}.
}
\seealso{
  The other algorithms \code{\link{expm}(x, method = *)}.

  \code{\link{expmCond}}, to compute the exponential-\emph{condition} number.
}
\examples{
## The *same* examples as in ../expm.Rd  {FIXME} --
x <- matrix(c(-49, -64, 24, 31), 2, 2)
expm.Higham08(x)

## ----------------------------
## Test case 1 from Ward (1977)
## ----------------------------
test1 <- t(matrix(c(
    4, 2, 0,
    1, 4, 1,
    1, 1, 4), 3, 3))
expm.Higham08(test1)
##                    [,1]               [,2]               [,3]
## [1,] 147.86662244637000 183.76513864636857  71.79703239999643
## [2,] 127.78108552318250 183.76513864636877  91.88256932318409
## [3,] 127.78108552318204 163.67960172318047 111.96810624637124
## -- these agree with ward (1977, p608)


## ----------------------------
## Test case 2 from Ward (1977)
## ----------------------------
test2 <- t(matrix(c(
    29.87942128909879, .7815750847907159, -2.289519314033932,
    .7815750847907159, 25.72656945571064,  8.680737820540137,
   -2.289519314033932, 8.680737820540137,  34.39400925519054),
           3, 3))
expm.Higham08(test2)
expm.Higham08(test2, balancing = FALSE)
##                   [,1]               [,2]               [,3]
##[1,]   5496313853692405 -18231880972009100 -30475770808580196
##[2,] -18231880972009160  60605228702221760 101291842930249376
##[3,] -30475770808580244 101291842930249200 169294411240850880
## -- in this case a very similar degree of accuracy.

## ----------------------------
## Test case 3 from Ward (1977)
## ----------------------------
test3 <- t(matrix(c(
    -131, 19, 18,
    -390, 56, 54,
    -387, 57, 52), 3, 3))
expm.Higham08(test3)
expm.Higham08(test3, balancing = FALSE)
##                    [,1]                [,2]                [,3]
##[1,] -1.5096441587713636 0.36787943910439874 0.13533528117301735
##[2,] -5.6325707997970271 1.47151775847745725 0.40600584351567010
##[3,] -4.9349383260294299 1.10363831731417195 0.54134112675653534
## -- agrees to 10dp with Ward (1977), p608. ??? (FIXME)

## ----------------------------
## Test case 4 from Ward (1977)
## ----------------------------
test4 <-
    structure(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1e-10,
                1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 1, 0),
    .Dim = c(10, 10))

E4 <- expm.Higham08(test4)
Matrix(zapsmall(E4))

S4 <- as(test4, "sparseMatrix") # some R based expm() methods work for sparse:
ES4 <- expm.Higham08(S4, bal=FALSE)
stopifnot(all.equal(E4, unname(as.matrix(ES4))))
## NOTE: Need much larger sparse matrices for sparse arith to be faster!

##
## example of computationally singular matrix
##
m <- matrix(c(0,1,0,0), 2,2)
eS <- expm.Higham08(m)  # "works"  (hmm ...)
}
\keyword{algebra}
\keyword{math}
\encoding{UTF-8}
