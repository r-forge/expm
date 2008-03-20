###  M^k  for a matrix  M and non-negative integer 'k'

matpow <- function(x, k)
    .Call(R_matpow, x, as.integer(k))

"%^%" <- matpow
