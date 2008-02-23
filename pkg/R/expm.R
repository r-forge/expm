### ===== File part of R package expm =====
###
### Function to compute the matrix exponential
###
###    exp(M) = sum(n = 0:Inf; M^n / n!),
###
### where M is an (n x n) matrix.
###

expm <- function(x, method = c("Ward77", "Pade", "Taylor", "PadeO", "TaylorO", "Eigen"),
		 order = 8, trySym = TRUE)
## NOTA BENE:  Matlab uses order = 6  !!!
{
    if (!is.matrix(x))
        stop("invalid (non-matrix) argument")
    method <- match.arg(method)
    if(method == "Ward77")
	## AUTHORS: Christophe Dutang, Vincent Goulet <vincent.goulet@act.ulaval.ca>
	.Call("do_expm", x)
    else if (method == "Eigen") {
	## matrix exponential using eigenvalues / spectral decomposition :
        ## ==  Dubious Way 'Method 14' : is
        ## good for 'symmetric' or 'orthogonal' (or other 'normal' : A'A = AA' ):

        ## MM: improved from mexp2() with 'trySym' and isSymmetric()
        isSym <- if(trySym) isSymmetric.matrix(x) else FALSE
        z <- eigen(x, sym = isSym)
        V <- z$vectors
        Vi <- if(isSym) t(V) else solve(V)
        Re(V %*% (    exp(z$values)   *  Vi)) ## ==
        ##(V %*% diag(exp(z$values)) %*% Vi)
    }
    else {

	## AUTHORS: Marina Shapira and David Firth --------------

	dx <- dim(x)
	if (dx[1] != dx[2])
	    stop("matrix not square")
	if (!is.numeric(order) || order != as.integer(order) || order < 0)
	    stop("order must be an integer number >= 0")

	storage.mode(x) <- "double"
	order <- as.integer(order)
	## MM:	a "silly"  way to code the method / order
	ntaylor <- npade <- as.integer(0)
	if (substr(method,1,4) == "Pade")
	    npade <- order else ntaylor <- order
	res <- .Fortran(if(identical(grep("O$", method), 1L))
			"matrexpO" else "matrexp",
			X = x,
			size = dx[1],
			ntaylor,
			npade,
			accuracy = double(1),
			PACKAGE = "expm")[c("X", "accuracy")]
	structure(res$X, accuracy = res$accuracy)
    }

}

### For quick "compatibility" --- FIXME later
mexp <- function(a,
		 order = 8,
		 ## different default method!
		 method = c("Pade", "Taylor", "PadeO", "TaylorO"))
{
    expm(a, method, order)
}
