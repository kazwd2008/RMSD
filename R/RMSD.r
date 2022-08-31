#----------------------------------------------------------
#' Modified Stahel-Donoho Estimators (Single core version)
#'
#' This function is for multivariate outlier detection.
#'    Ver.1.6 2009/07/14   Published at http://www.stat.go.jp/training/2kenkyu/pdf/ihou/67/wada1.pdf (in Japanese)
#'    Ver.1.7 2018/10/19   Modify gso function to stop warning messages
#'    Ver.2   2021/09/10   Added the outlier detection step
#'
#' @param inp imput data (a numeric matrix)
#' @param nb  number of basis
#' @param sd  seed (for reproducibility)
#' @param pt  threshold for outlier detection (probability)
#'
#' @return a list of the following information
#' \itemize{
#'   \item u  final mean vector
#'   \item V  final covariance matrix
#'   \item wt final weights
#'   \item mah Squared mahalanobis distance of each observation
#'   \item FF F test statistics
#'   \item cf threshold to detect outliers (percentile point)
#'   \item ot 1:normal observation, 2:outlier
#' }
#' @importFrom stats mad
#' @importFrom stats mahalanobis
#' @importFrom stats runif
#' @importFrom stats qchisq
#' @importFrom stats median
#' @importFrom stats qf
#' @export

RMSD <- function(inp, nb=0, sd=0, pt=0.999) {

inp_d <- ncol(inp)            # number of variables
inp_n <- nrow(inp)            # number of observations

##############################
# create orthogonal bases
##############################

if (sd != 0) set.seed(sd)

## number of orthogonal bases
bb_n <- ifelse(nb==0, trunc(exp(2.1328+0.8023*inp_d) / inp_d), nb)
rn <- bb_n * inp_d^2           # number of random variables for creating the bases
basis <- array(runif(rn), c(inp_d, inp_d, bb_n))

## orthonormalization by function "gso"
basis <- apply(basis, 3, gso)
basis <- array(basis, c(inp_d, inp_d, bb_n))

##############################
# projection and residual computation
##############################

prj <- array(0, c(inp_n, inp_d, bb_n))  # for projection
res <- array(0, c(inp_n, inp_d, bb_n))  # resudal
wt <- array(0, c(inp_n, inp_d, bb_n))   # weight by observations x variable x bases
wts <- array(0, c(inp_n, bb_n))         # weight by observations x bases
bwt <- rep(0, inp_n)                    # the smallest weight by observations
kijun <- qchisq(0.95, inp_d)            # reference for trimming

Fprj <- function(pj) t(pj %*% t(inp))   # projection
prj <- apply(basis, 3, Fprj)
prj <- array(prj, c(inp_n, inp_d, bb_n))

medi <- apply(prj, c(2, 3), median)     # median
madx <- apply(prj, c(2, 3), mad)        # median absolute deviation (MAD) / 0.674
for (i in 1:bb_n) {                     # robust standardization of residuals
    res[,,i] <- t(abs(t(prj[,,i]) - medi[,i]) / madx[,i])
}

### trimming weight by a Huber-like weight function
k0 <- which(res <= sqrt(kijun))
k1 <- which(res > sqrt(kijun))
wt[k0] <- 1
wt[k1] <- kijun / (res[k1]^2)
wts <- apply(wt, c(1,3), prod)
bwt <- apply(wts, 1, min)               # selecting the smallest weight

### initial robust covariance matrix
u1 <- apply(inp * bwt, 2, sum) / sum(bwt)
V1 <- t(t(t(inp) - u1) * bwt) %*% (t(t(inp) - u1) * bwt) / sum(bwt^2)

### avoiding NaN error
u1 <- ifelse(is.nan(u1), 0, u1)
V1 <- ifelse(is.nan(V1), 0, V1)

### robust PCA (LAPACK)
eg <- eigen(V1, symmetric=TRUE)
ctb <- eg$value / sum(eg$value)         # contribution ratio

##############################
# projection pursuit (PP)
##############################
res2 <- array(0, c(inp_n, inp_d))       # residuals
wt2  <- array(0, c(inp_n, inp_d))       # weight by observations x variables
wts2 <- array(0, inp_n)                 # final weight by observations

prj2 <- t(eg$vector %*% (t(inp) - u1))  # projection
medi2 <- apply(prj2, 2, median)         # median and
madx2 <- apply(prj2, 2, mad)            # MAD for standardization
res2 <- t(abs(t(prj2) - medi2) / madx2) # standardized residuals

# trimming by the Huber-like weight function
k0 <- which(res2 <= sqrt(kijun))
k1 <- which(res2 > sqrt(kijun))
wt2[k0] <- 1
wt2[k1] <- kijun / (res2[k1]^2)

wts2 <- apply(wt2, 1, prod)
wts2 <- pmin(wts2, bwt)

##############################
# final mean vector and covariance matrix
##############################

u2 <- apply(inp * wts2, 2, sum) / sum(wts2)
V2 <- t(t(t(inp) - u2) * wts2) %*% (t(t(inp) - u2) * wts2) / sum(wts2^2)

##############################
# outlier detection
##############################

# compute mahalanobis distance of each observation for outlier detection
mah <- mahalanobis(inp, u2, V2)
FF <- mah * ( inp_n - inp_d )* inp_n /(( inp_n^2 - 1 )* inp_d)	# test statistics
cf <- qf(pt, inp_d, inp_n - inp_d) 	# 99.9 percentile point of F distribution
ot <- rep(1, inp_n)			# outlier flag (1: OK;   2: outlier)
ot[which(FF > cf)] <- 2 		# outliers are those exceed the threshold

return(list(u=u2, V=V2, wt=wts2, mah=mah, FF=FF, cf=cf, ot=ot))
}

###########################################################
# gso: Gram-Schmidt Orthonormalization for function "msd"
###########################################################

gso <- function(basis) {
    bd <- ncol(basis)
    basis[1,] <- as.vector(basis[1,]) / sqrt(as.vector(t(basis[1,]) %*% basis[1,]))
    for (i in 2 : bd ) {
        wk1 <- as.vector(basis[i,])
        for (j in 1:(i-1)) {
            wk2 <- as.vector(basis[j,])
            basis[i,] <- wk1 - (as.vector(t(wk1) %*% wk2) * wk2)
        }
        basis[i,] <- as.vector(basis[i,]) / sqrt(as.vector(t(basis[i,]) %*% basis[i,]))
    }
    return(basis)
}

################################################################################
