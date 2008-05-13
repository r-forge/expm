### ===== File part of R package expm =====
###
### Function to compute the matrix logarithm
###
###    exp(M) = sum(n = 0:Inf; M^n / n!),
###
### where M is an (n x n) matrix.
###




logm <- function(x, method = c("Eigen"),
		 order = 8,
                 trySym = TRUE, tol = .Machine$double.eps)
{
    if (!is.matrix(x))
        stop("invalid (non-matrix) argument")
    method <- match.arg(method)
    if (method == "Eigen") {
    ## AUTHOR: Christophe Dutang
    ## matrix exponential using eigenvalues / spectral decomposition and
    ## Ward(1977) algorithm if x is numerically non diagonalisable
        .Call("do_logm_eigen", x, tol)
    }

}

