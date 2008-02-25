###  M ^ k  for a matrix  M   and  non-negative integer 'k'

## "%^%" <- function(x, k) if(k == 1) x else x %*% (x %^% (k-1))

"%^%" <- function(x, k)
{
    stopifnot(k == round(k))
    .Call(R_matpow, x, as.integer(k))
}
