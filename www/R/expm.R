### ===== File part of R package expm =====
###
### Function to compute the matrix exponential
###
###    exp(M) = sum(n = 0:Inf; M^n / n!),
###
### where M is an (n x n) matrix.
###
### AUTHORS: Christophe Dutang, Vincent Goulet <vincent.goulet@act.ulaval.ca>

expm <- function(x)
    .Call("do_expm", x)
