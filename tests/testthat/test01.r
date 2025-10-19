
# create test dataset following normal distribution with correlation

# source("http://aoki2.si.gunma-u.ac.jp/R/src/gendat2.R", encoding="euc-jp")  
gendat2 <- function(	nc,  # sample size
			r)   # correlation coefficient
{
	z <- matrix(rnorm(2*nc), ncol=2)
	r2 <- cor(z)
	res <- eigen(r2)
	coeff <-  solve(r2) %*% t(sqrt(res$values)*t(res$vectors))
	z <- scale(z) %*% coeff
	return(z %*% chol(matrix(c(1, r, r, 1), ncol=2)))
} # gendat2 is provided by Prof. Aoki of Gunma Univ.

nc <- 100  # number of records
z <- gendat2(nc, 0.8)
ot1 <- RMSD(z)

cqd99   <- qchisq(0.99, df=2)	# 99 percentile point
cqd999   <- qchisq(0.999, df=2)	# 99.9 percentile point

eg <- eigen(ot1$V)		# calculation of eigen value and vector 
P <- eg$vectors			
D <- matrix(rep(0, 2*2), ncol=2)	
diag(D) <- eg$values
PP <- solve(P)			

ax1 <- sqrt(cqd99 * eg$values)	
ax2 <- sqrt(cqd999 * eg$values)	

nb  <- 0:200			　
dw <- 2 * pi / 200	
w  <- dw * nb
XA1 <- ax1[1] * cos(w)		# 95%点の場合
XA2 <- ax1[2] * sin(w)
XB1 <- ax2[1] * cos(w)		# 99%点の場合
XB2 <- ax2[2] * sin(w)
XX1 <- t(t(PP) %*% t(cbind(XA1, XA2)) )  
XX2 <- t(t(PP) %*% t(cbind(XB1, XB2)))  

plot(z, col=ot1$ot, pch=19, main="Scatterplot (outliers are shown in red)")
lines(XX1, col=2, lwd=2, lty=3)		# 楕円描画　95%点の場合
lines(XX2, col=4, lwd=2, lty=3)  	# 楕円描画　99%点の場合
