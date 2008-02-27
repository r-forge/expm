#### Examples where we know the result "exactly"

library(expm)

source(system.file("test-tools.R", package= "expm"))## -> assertError()...
## rMat(), ..

### For nilpotent matrices A, exp(A) is polynomial in A
###  Mathworld gives the example of  the general  3 x 3  upper triangle
nilA3 <- function(x,y,z) {
    ## Purpose: simple nilpotent matrix 3x3  A (with A^n = 0 for n >= 3)
    ##          and its exact matrix exponential
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 19 Jan 2008
    stopifnot((n <- length(x)) == 1, length(y) == 1, length(z) == 1,
              is.numeric(x), is.numeric(y), is.numeric(z))
    list(A = cbind(0, rbind(matrix(c(x,0,z,y), 2,2), 0)),
         expA = cbind(c(1,0,0), c(x,1,0), c(z + x*y/2, y, 1)))
}



re.nilA3 <- function(xyz, EXPM)
{
### TODO:  allow several kinds of EXPM() simultaneously (for same 'r')
    EXPM <- match.fun(EXPM)
    r <- do.call(nilA3, as.list(xyz))
    relErr(r$expA, EXPM(r$A))
}

set.seed(321)
re <- replicate(1000,
                c(re.nilA3(rlnorm(3),function(x)expm(x,"Pade")),
                  re.nilA3(rnorm(3), function(x)expm(x,"Pade"))))

summary(t(re))
stopifnot(rowMeans(re) < 1e-15,
          apply(re, 1, quantile, 0.80) < 1e-16,
          apply(re, 1, quantile, 0.90) < 2e-15,
          apply(re, 1, max) < c(4e-14, 6e-15))

cat('Time elapsed: ', (p1 <- proc.time()),'\n') # for ``statistical reasons''


## Check *many* random nilpotent matrices:
set.seed(321)
RE <- replicate(1000,
                c(re.nilA3(rlnorm(3), function(x) expm(x, "Ward77")),
                  re.nilA3(rnorm(3),  function(x) expm(x, "Ward77"))))
stopifnot(rowMeans(RE) < 1e-15,
          apply(RE, 1, quantile, 0.80) < 1e-16,
          apply(RE, 1, quantile, 0.90) < 2e-15,
          apply(RE, 1, max) < c(4e-14, 6e-15))

print(summary(t(RE)))
cat("relErr(expm(.,Pade)) - relErr(expm(.,Ward77))  in Machine_eps units:\n")
print(summary(c(re - RE)) / .Machine$double.eps)
##       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
## -0.6183442  0.0000000  0.0000000  1.3650410  0.1399719 94.9809161



cat('Time elapsed: ',(p2 <- proc.time())-p1,'\n') # for ``statistical reasons''

###--- A second group --- where we know the diagonalization of A ---

if(!require("Matrix"))
    q('no')
##  ------  the rest really uses 'Matrix'
##---> now use  expm::expm()  since Matrix has its own may mask the expm one
##              ^^^^^^^^^^^^

## rMat() relies on Matrix::rcond():

## Now with the change default rcondMin, this "works"
system.time(R40 <- rMat(40))
system.time(R80 <- rMat(80))


if(FALSE) # unfinished idea:  Construct  M * diag() * solve(M)
facMat <- function(n, R_FUN) {
    R_FUN <- match.arg(R_FUN)


}


###---
rnilMat <- function(n, R_FUN = function(n) rpois(n, lambda=5))
{
    ## random upper triangular (zero-diagonal) nilpotent  n x n matrix
    m <- matrix(0, n,n)
    ut <- upper.tri(m)
    R_FUN <- match.fun(R_FUN)
    m[ut] <- R_FUN(sum(ut))
    m
}

set.seed(17)
m <- rnilMat(10)
as(m, "sparseMatrix")# for nicer printing
E.m <- expm::expm(m, method="Pade")
as(E.m, "sparseMatrix")

(dN <- 9*7*320) # 20160
stopifnot(abs(round(E.m * dN)  -  (E.m * dN)) < 9e-6)
EmN <- matrix(c(dN, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                3*dN, dN, 0, 0, 0, 0, 0, 0, 0, 0,
                352800, 5*dN, dN, 0, 0, 0, 0, 0, 0, 0,
                1018080, 332640, 5*dN, dN, 0, 0, 0, 0, 0, 0,
                2235240, 786240, 292320, 3*dN, dN, 0, 0, 0, 0, 0,
                9368520, 3483480, 1582560, 413280, 181440, dN, 0, 0, 0, 0,
                24676176, 9598680, 5073600, 1562400, 826560, 161280, dN, 0,0,0,
                43730160, 17451000, 10051440, 3430560, 1955520, 504000,
                5*dN, dN, 0, 0,
                68438436, 27747480, 16853760, 6036240, 3638880, 1038240,
                252000, 3*dN, dN, 0,
                119725855, 49165892, 31046760, 11652480, 7198800, 2264640,
                614880, 191520, 3*dN, dN),
              10, 10)

Em.xct <- EmN / dN

stopifnot(all.equal(E.m, Em.xct,
                    check.attributes = FALSE, tol= 1e-13))
expmList <-
    list(function(x) expm::expm(x,method="Ward77"),
         function(x) expm::expm(x,method="Pade"),
         function(x) expm::expm(x,method="PadeO"),
         function(x) expm::expm(x,method="Taylor"),
         function(x) expm::expm(x,method="TaylorO"))
sapply(expmList, function(EXPM) relErr(Em.xct, EXPM(m)))
## 1.079e-16 4.5054e-14 4.5025e-14 3.716219e-17 7.078512e-18
## "Ward77" expm() is again  much more accurate than  s+Pade+s,
## *but* s+Taylor+s is even more accurate
