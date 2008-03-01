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
epsC <- .Machine$double.eps
cat("relErr(expm(.,Pade)) - relErr(expm(.,Ward77))  in Machine_eps units:\n")
print(summary(c(re - RE)) / epsC)
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

facMat <- function(n, R_FUN, ev = R_FUN(n), M = rMat(n, R_FUN = R_FUN))
{
    ## Purpose: Construct random matrix x of which we "know" expm(x)
    ## because we set  x :=  M %*% diag(ev)  %*% solve(M)
    ## ----------------------------------------------------------------------
    ## Arguments: n:     dimension of matrices
    ##            R_FUN: random number generator function (n)
    ##            ev:    numeric length-n vector of eigenvalues
    ##            M:     n x n matrix. Note that the default,
    ##                   rMat() will give matrices ``not close to singular''
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: Feb 2008
    R_FUN <- match.fun(R_FUN)
    stopifnot(n > 0, is.numeric(ev), length(ev) == n,
              dim(M) == c(n,n), is.numeric(M))

    iM <- solve(M)
    ## D <- diag(ev); A = M %*% D %*% iM
    list(A    = M %*% (ev      * iM),
         expA = M %*% (exp(ev) * iM))
}

re.facMat <- function(n, EXPMlist, rFUN = rnorm, ...)
{
    stopifnot(is.list(EXPMlist))
    r <- facMat(n, rFUN, ...)
    sapply(EXPMlist, function(EXPM) relErr(r$expA, EXPM(r$A)))
}

expm.safe.Eigen <- function(x, silent = FALSE) {
    r <- try(expm::expm(x, "R_Eigen"), silent = silent)
    if(inherits(r, "try-error")) NA else r
}

expmList <-
    list(Ward  = function(x) expm::expm(x, "Ward77"),
         s.P.s = function(x) expm::expm(x, "Pade"),
         s.P.sO= function(x) expm::expm(x, "PadeO"),
         s.T.s = function(x) expm::expm(x, "Taylor"),
         s.T.sO= function(x) expm::expm(x, "TaylorO"),
         Eigen = expm.safe.Eigen,
         hybrid= function(x) expm::expm(x, "hybrid")
         )
expmL.wo.E <- expmList[names(expmList) != "R_Eigen"]


set.seed(12)
re.facMat(20, expmList)
fRE <- replicate(100, re.facMat(20, expmList))

## Now look at that:
boxplot(t(fRE), log="y", notch=TRUE,
        main = "relative errors for 'random' eigen-ok 20 x 20 matrix")

## Now  try an example with badly conditioned "random" M matrix...
## ...
## ... (not yet)


### --- The 2x2 example with bad condition , see A3 in ./ex2.R
m2ex3 <- function(eps = 0) {
    stopifnot(is.numeric(eps), length(eps) == 1)
    A <- rbind(c(-1,     1),
               c(eps^2, -1))
    I.e <- 1 - eps^2 / 2
    V <- I.e* rbind(    c(-1, 1),
                    eps*c( 1, 1))
    D <- c(-1-eps, -1+eps)
    iV <- ## solve(V)
        rbind(c(-1, 1/eps),
              c( 1, 1/eps)) / (2 * I.e)
    ## NOTE:  kappa(V) = condition_number(V) == 1/eps exactly
    useTol <- 2e-16 / eps
    stopifnot(all.equal(diag(2), V %*% iV,       tol=useTol),
              all.equal(A, V %*% diag(D) %*% iV, tol=useTol) )
    ch.e <- cosh(eps)
    sh.e <- sinh(eps)
    list(A = A,
         expA = exp(-1) *
         rbind(c( ch.e,  sh.e/eps),
               c(sh.e*eps, ch.e  )))
}
re.m2ex3 <- function(eps, EXPMlist)
{
    stopifnot(is.list(EXPMlist))
    r <- m2ex3(eps)
    sapply(EXPMlist, function(EXPM) relErr(r$expA, EXPM(r$A)))
}

re.m2ex3(1e-8, expmList)# Ward is best; s.P.s fair; Eigen (and hybrid): ~1e-9

eps <- 10^-(1:18)
t.m2 <- t(sapply(eps, re.m2ex3, EXPMlist = expmList))
## --> 5 error messages from try(... "R_Eigen" ...)
t.m2[, c(1:2,4,6:7)] ## very impressive !!  Ward is best!

## library(RColorBrewer)
## Bcol <- brewer.pal(7,"Dark2")
Bcol <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A",
          "#66A61E", "#E6AB02", "#A6761D")
matplot(eps, t.m2, type = "b", log = "xy", col=Bcol, lty = 1:7, pch=1:7,
        axes=FALSE, frame = TRUE,
        xlab = expression(epsilon), ylab = "relative error",
        main = expression(expm(A, method == "*") *"  relative errors for  " *
            A == bgroup("[", atop({-1} *"  "* 1, {epsilon^2} *"  "*{-1}), "]")))
legend("bottomright",colnames(t.m2),       col=Bcol, lty = 1:7, pch=1:7,
       inset = 0.02)
if(require("sfsmisc")) {
    sfsmisc::eaxis(1, labels=FALSE)
    sfsmisc::eaxis(1, at = eps[c(TRUE,FALSE)])
    sfsmisc::eaxis(2, labels=FALSE)
    op <- par(las=2)
    sfsmisc::eaxis(2, at = axTicks(2,log=TRUE)[c(TRUE,FALSE,FALSE)])
    par(op)
} else {
    axis(1)
    axis(2)
}

## typical case:
ep <- 1e-10
(me <- m2ex3(ep))
me$expA * exp(1) ## the correct value ; numerically identical to simple matrix:
stopifnot(identical(me$expA * exp(1),
                    rbind(c(  1,  1),
                          c(ep^2, 1))))
## The relative error (matrices):
lapply(expmList, function(EXPM) 1 - EXPM(me$A)/me$expA)

## Average number of correct digits [less "extreme" than plot above]
nDig <- sapply(expmList, function(EXPM) -log10(mean(abs(1 - EXPM(me$A)/me$expA))))
round(nDig, 2)
##   Ward  s.P.s s.P.sO  s.T.s s.T.sO  Eigen hybrid
##  16.26  14.65  14.65  14.65  14.65   6.20   6.39  [AMD Opteron 64-bit]
##    Inf  14.65  14.65  14.65  14.65   6.74   6.33  [Pentium-M (32-bit)]

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
sapply(expmL.wo.E, function(EXPM) relErr(Em.xct, EXPM(m)))
## result depends quite a bit on platform:

## Pentium-M 32-bit ubuntu gave
##      Ward     s.P.s    s.P.sO     s.T.s    s.T.sO    hybrid
## 1.079e-16 4.505e-14 4.503e-14 3.716e-17 7.079e-18 1.079e-16

## "Ward77" expm() is again more accurate than  s+Pade+s,
## *but* s+Taylor+s is even more accurate

## but on 64-bit AMD Opterons
##         Ward        s.P.s       s.P.sO        s.T.s       s.T.sO       hybrid
## 4.424070e-17 3.989036e-17 3.989036e-17 8.435227e-17 8.435227e-17 4.424070e-17


cat('Time elapsed: ',(p3 <- proc.time())-p2,'\n') # for ``statistical reasons''
