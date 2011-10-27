\name{expAtv}
\alias{expAtv}
\title{Compute Matrix Exponential  exp(A t) * v  directly}
\description{Performs exp(A t) \%*\% v  directly  w/o  evaluating exp(A)
  Originally by  Roger B. Sidje
  EXPOKIT: Software Package for Computing Matrix Exponentials.
  ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
}
\usage{
expAtv(A, v, t = 1,
       tol = 1e-07, btol = 1e-07, m.max = 30, mxrej = 10,
       verbose = getOption("verbose"))
}
\arguments{
  \item{A}{n x n matrix}
  \item{v}{n - vector}
  \item{t}{number (scalar) ___ FIXME allow vector ? ___}
  \item{tol}{}
  \item{btol}{}
  \item{m.max, mxrej}{integer constants you should only change if you know what you're doing}
%  \item{mxrej}{}
  \item{verbose}{flag indicating if the algorithm should be verbose..}
}
\value{a list with components
  \item{eAtv}{}
}
\details{

}
\references{
............
}
\author{Ravi Varadhan, Johns Hopkins University;
  Martin Maechler (cosmetic, generalization to sparse matrices).
}
% \note{
% }

\seealso{
  \code{\link{expm}},....
}
\examples{
source(system.file("demo", "exact-fn.R", package = "expm"))
set.seed(1)
(s5 <- Matrix(m5 <- rnilMat(5))); v <- c(1,6:9)
(em5 <- expm(m5))
r5 <- expAtv(m5, v)
stopifnot(all.equal(em5 \%*\% v, r5$eAtv))

v <- 10:1
with(xct10, all.equal(expm(m), expm))
all.equal(c(xct10$expm \%*\% v),
          expAtv(xct10$m, v)$eAtv)

n <- 500
A <- bandSparse(n,n, -1:1, diag = list(-(2:n), -5*(1:n), 1:(n-1)))
v <- 100*(500:1)
t.v <- system.time(rr  <- expAtv(A, v=v))
## this is an order of magnitude slower :
t.A <- system.time(w. <- (eA <- expm(A)) \%*\% v)
stopifnot(all.equal(rr$eAtv, as.numeric(w.)))
mean((t.A / t.v)[c(1,3)]) ## 23.57  {nb-mm3}
}
\keyword{algebra}
\keyword{math}
