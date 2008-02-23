#### Will be sourced by several R scripts in ../tests/

paste0 <- function(...) paste(..., sep = '')

identical3 <- function(x,y,z)	identical(x,y) && identical (y,z)
identical4 <- function(a,b,c,d) identical(a,b) && identical3(b,c,d)

## Make sure errors are signaled
assertError <- function(expr) {
    d.expr <- deparse(substitute(expr))
    t.res <- try(expr, silent = TRUE)
    if(!inherits(t.res, "try-error"))
	stop(d.expr, "\n\t did not give an error", call. = FALSE)
    invisible(t.res)
}

is.all.equal3 <- function(x,y,z, tol = .Machine$double.eps^0.5)
    isTRUE(all.equal(x,y, tol=tol)) && isTRUE(all.equal(y,z, tol=tol))

is.all.equal4 <- function(x,y,z,u, tol = .Machine$double.eps^0.5)
    is.all.equal3(x,y,z, tol=tol) && isTRUE(all.equal(z,u, tol=tol))


## The relative error typically returned by all.equal:
relErr <- function(target, current)
    mean(abs(target - current)) / mean(abs(target))

expm.t.identity <- function(x, method,
                            tol = .Machine$double.eps^0.5,
                            check.attributes = FALSE,
                            ...)
{
  ## Purpose: Test the identity   expm(A') = (expm(A))'
  ## ----------------------------------------------------------------------
  ## Arguments: method, ... :          arguments to  expm()
  ##            tol, check.attributes: arguments to  all.equal()
  ## ----------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 23 Feb 2008, 17:26
    ex <- expm::expm(x   , method=method, ...)
    et <- expm::expm(t(x), method=method, ...)
    all.equal(t(ex), et, tol = tol, check.attributes = check.attributes)
}

### This is similar to Matrix'  example(spMatrix) :
##
rSpMatrix <- function(nrow, ncol = nrow, nnz,
		      rand.x = function(n) round(100 * rnorm(nnz)))
{
    ## Purpose: random sparse matrix
    ## --------------------------------------------------------------
    ## Arguments: (nrow,ncol): dimension
    ##		nnz  :	number of non-zero entries
    ##	       rand.x:	random number generator for 'x' slot
    ## --------------------------------------------------------------
    ## Author: Martin Maechler, Date: 14.-16. May 2007
    stopifnot((nnz <- as.integer(nnz)) >= 0,
	      nrow >= 0, ncol >= 0,
	      nnz <= nrow * ncol)
##     spMatrix(nrow, ncol,
##		i = sample(nrow, nnz, replace = TRUE),
##		j = sample(ncol, nnz, replace = TRUE),
##		x = rand.x(nnz))
    m <- matrix(0, nrow, ncol)
    m[cbind(i = sample(nrow, nnz, replace = TRUE),
	    j = sample(ncol, nnz, replace = TRUE))] <- rand.x(nnz)
    m
}


## This needs "Matrix" package
rMat <- function(n, R_FUN = rnorm,
                 rcondMin = 1.4 * n ^ -1.6226,
                 iterMax = 100)
{
    ## Purpose: random square matrix "not close to singular"
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## NOTE: needs  Matrix::rcond()
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 19 Jan 2008
    ##
    ##--> /u/maechler/R/MM/Pkg-ex/Matrix/rcondition-numb.R  researches rcond( <random mat>)
    ## Result :
    ##   -log[rcond] = log(Kappa) = 1.051 + 1.6226 * log(n)
    ##   ==================================================
    ##   1/rcond = Kappa = exp(1.051 + 1.6226 * log(n))
    ##                   = 2.8605 * n ^ 1.6226
    ##   ==================================================

    ## since we *search* a bit, take a factor ~ 4  higher rcond:
    ##  4 / 2.8605 ~ 1.4 --> default of rcondMin  above

    it <- 1
    rcOpt <- 0
    repeat {
        M <- matrix(R_FUN(n^2), n,n)
        if((rc <- Matrix::rcond(M)) >= rcondMin) break
        if(rc > rcOpt) {
            rcOpt <- rc
            M.Opt <- M
        }
        if((it <- it+1) > iterMax) {
            warning("No Matrix found with rcond() >= ",format(rcondMin),
                    "\n Achieved rcond() = ", format(rcOpt),"\n")
            M <- M.Opt
            break
        }
    }
    M
}
