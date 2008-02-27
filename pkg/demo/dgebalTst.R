dgebalTst <- function(A) {

    ## Purpose: Consistency checking of	 dgebal()
    ## ----------------------------------------------------------------------
    ## Arguments: a square matrix
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, 20 Feb 2008 and on

    n <- dim(A)[1]
    ## do *the* three calls and look at result
    P <- dgebal(A, "P")
    didPerm <- ((leftP	<- (i1 <- P$i1) != 1L) |
		(rightP <- (i2 <- P$i2) != n))
    if(didPerm) {## *had* permutation -- now check my idea about it
	pp <- as.integer(P$scale)
	## Permute A to become P$z :
	A. <- A
	if(rightP) {## The upper part
	    for(i in n:(i2+1)) { # 'p2' in *reverse* order
		## swap	 i <-> pp[i]   both rows and columns
		tt <- A.[,i]; A.[,i] <- A.[,pp[i]]; A.[,pp[i]] <- tt
		tt <- A.[i,]; A.[i,] <- A.[pp[i],]; A.[pp[i],] <- tt
	    }
	}
	if(leftP) {## The lower part
	    for(i in 1:(i1-1)) { # 'p1' in *forward* order
		tt <- A.[,i]; A.[,i] <- A.[,pp[i]]; A.[,pp[i]] <- tt
		tt <- A.[i,]; A.[i,] <- A.[pp[i],]; A.[pp[i],] <- tt
	    }
	}
	stopifnot(isTRUE(all.equal(A., P$z, tol = 1e-15)))

	## Now the reverse: Use pp[] and permute  A. "back to A":
	if(leftP) {## The lower part
	    for(i in (i1-1):1) { # 'p1' in *reverse* order
		tt <- A.[,i]; A.[,i] <- A.[,pp[i]]; A.[,pp[i]] <- tt
		tt <- A.[i,]; A.[i,] <- A.[pp[i],]; A.[pp[i],] <- tt
	    }
	}
	if(rightP) {## The upper part
	    for(i in (i2+1):n) { # 'p2' in *forward* order
		## swap	 i <-> pp[i]   both rows and columns
		tt <- A.[,i]; A.[,i] <- A.[,pp[i]]; A.[,pp[i]] <- tt
		tt <- A.[i,]; A.[i,] <- A.[pp[i],]; A.[pp[i],] <- tt
	    }
	}
	stopifnot(isTRUE(all.equal(A., A, tol = 1e-15)))

    }
    S <- dgebal(P$z, "S")# "S" starting from result of "P"
    stopifnot(S$i1 == 1, S$i2 == n)

    ## Now check the scaling
    d <- S$scale
    ## A.scaled <- diag(1/d, n) \%*\% P$z \%*\% diag(d, n)
    ## more efficiently:
    A.scaled <- P$z * (rep(d, each = n) / d)
    stopifnot(isTRUE(all.equal(S$z, A.scaled, tol = 1e-15)))
    ## Check the reverse:
    S.rescaled <- S$z * (d * rep(1/d, each = n))
    stopifnot(isTRUE(all.equal(P$z, S.rescaled, tol = 1e-15)))

    B <- dgebal(A, "B")# "B" : B[oth]
    stopifnot(P$i1 == B$i1, P$i2 == B$i2)
    list(P = P, S = S, B = B, Sz.eq.Bz = isTRUE(all.equal(S$z, B$z)))
}

m4. <- rbind(c(-1,-2, 0, 0),
             c( 0, 0,10,11),
             c( 0, 0,12, 0),
             c( 0,13, 0, 0))
str(b4. <- dgebalTst(m4.))

## better (?) example
(m <- matrix(c(0,-1,0,-2,10, rep(0,11)), 4,4))
str(ba <- dgebalTst(m))
## Hmm: here S$z  *differs*  from B$z
## ---  but at least, the scale[] and z[] returned seem ok


## a non-empty ``less-balanced'' example  ---

m4 <- matrix(outer(2^(0:7),c(-1,1)), 4,4)
m4[lower.tri(m4)] <- 0 #--> upper triangular ==> will have many permutations
## now permute it; so dgebal() will find the permutation
p <- c(4,2:1,3); m4 <- m4[p,p]
m4

str(dm4 <- dgebalTst(m4)) # much permutation!  i1 = i2 = 1 !
