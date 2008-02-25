###  M ^ k  for a matrix  M   and  non-negative integer 'k'

## "%^%" <- function(x, k) if(k == 1) x else x %*% (x %^% (k-1))

"%^%" <- function(x, k)
{
    if(!is.matrix(x))
	stop("x must be a matrix")
    d <- dim(x)
    if((n <- d[1]) != d[2])
	stop("x must be a square matrix")
    stopifnot(k == round(k))
    ## TODO: .Call(do_matpow, x, n, k)
    .Fortran(matpowf,
	     x,
	     as.integer(n),
	     as.integer(k),
	     z = matrix(0, n,n))$z
}
