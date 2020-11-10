'***********************************
Author:  Jack Young, November 2020 *
************************************
The following R code is based upon MATLAB code created by Shireen Elhabian and Aly A. Farag,
"A Tutorial on Data Reduction:  Linear Discriminant  Analysis (LDA)," University of Louisville, 
CVIP Lab, September 2009, http://www.sci.utah.edu/~shireen/pdfs/tutorials/Elhabian_LDA09.pdf 
***********************************'
# Example usage:
# l <- processthreeclasseslda()
# Projecting points onto eigenvector
# y1_w1 <- as.vector(l$c1 %*% l$eigenvectors[,1]); y1_w1
# y2_w1 <- as.vector(l$c2 %*% l$eigenvectors[,1]); y2_w1
# y3_w1 <- as.vector(l$c3 %*% l$eigenvectors[,1]); y3_w1

processthreeclasseslda2 <- function()
{
  set.seed(262)
  size <- 500
  # Generate random data from uniform distribution
  u1 <- cbind(runif(size), runif(size))
  u2 <- cbind(runif(size), runif(size))
  
  # Using u1 and u2, use Box-Muller method to populate the feature matrix
  # Feature matrix will follow a standard normal distribution
  x1 <- matrix(nrow=size, ncol=2)
  for (i in 1:size)
    {
    x1[i,1] = sqrt(-2*log(u1[i,1]))*cos(2*pi*u2[i,1])
    x1[i,2] = sqrt(-2*log(u1[i,2]))*cos(2*pi*u2[i,2])
  }
  # Base data mean
  mu <- c(5,5)
  # c1 data mean
  mu1 <- c(mu[1]-3,mu[2]+7)
  # c2 data mean
  mu2 <- c(mu[1]-2.5,mu[2]-3.5)
  # c3 data mean
  mu3 <- c(mu[1]+7,mu[2]+5)
  # Create input covariance matrices to transform baseline x1 data
  cov1 <- rbind(c(5,-1),c(-3,3))
  cov2 <- rbind(c(4,0),c(0,4))
  cov3 <- rbind(c(3.5,1),c(3,2.5))
  # Store input eigenvalues and eigenvectors to transform baseline x1 data
  eigen1 <- eigen(cov1)
  eigen2 <- eigen(cov2)
  eigen3 <- eigen(cov3)
  c1 <- x1
  c2 <- x1
  c3 <- x1

  # Modify c1 matrix data with target variance
 for (i in 1:size)
 {
     c1[i,] <- eigen1$vectors %*% sqrt(diag(eigen1$values)) %*%x1[i,]
 }
  # Modify c2 matrix data with target variance
  for (i in 1:size)
  {
    c2[i,] <- eigen2$vectors %*% sqrt(diag(eigen2$values)) %*%x1[i,]
  }
  # Modify c3 matrix data with target variance
  for (i in 1:size)
  {
    c3[i,] <- eigen3$vectors %*% sqrt(diag(eigen3$values)) %*%x1[i,]
  }
  
  # Modify c1 matrix data with target mean
 c1[,1] <- c1[,1] + mu1[1]
 c1[,2] <- c1[,2] + mu1[2]
 
 # Modify c2 matrix data with target mean
 c2[,1] <- c2[,1] + mu2[1]
 c2[,2] <- c2[,2] + mu2[2]
 
 # Modify c3 matrix data with target mean
 c3[,1] <- c3[,1] + mu3[1]
 c3[,2] <- c3[,2] + mu3[2]
 
 # Calculate means from c1 matrix
 mean1 <- apply(c1, 2,mean)
 # Calculate means from c2 matrix
 mean2 <- apply(c2, 2,mean)
 # Calculate means from c3 matrix
 mean3 <- apply(c3, 2,mean)
 # Calculate overall means from c1, c2 & c3
 overallmean <- (mean1 + mean2 + mean3)/3
 
 # Calculate s1, s2, s3 & sw
 s1 <- cov(c1)
 s2 <- cov(c2)
 s3 <- cov(c3)
 sw <- s1 + s2 + s3
 
 # Calculate n1, n2 & n3 
 n1 <- nrow(c1)
 n2 <- nrow(c2)
 n3 <- nrow(c3)
 
 # Calculate sb1, sb2, sb3 & sb
 sb1 <- n1*(mean1 - overallmean)%*%t(mean1 - overallmean)
 sb2 <- n2*(mean2 - overallmean)%*%t(mean2 - overallmean)
 sb3 <- n3*(mean3 - overallmean)%*%t(mean3 - overallmean)
 sb <- sb1 + sb2 + sb3
 # Calculate projection vectors
 ldaprojection <- eigen(solve(sw)%*%sb)
 eigenvals <- ldaprojection$values
 eigenvectors <- ldaprojection$vectors
 par(mfrow=c(2,2),mar=c(4,4,2,2))
 
 # Visualize results

 plot(c1, xlim=c(-15,25), ylim=c(-10,25), col="red", axes = FALSE, panel.first = grid(), 
      xlab="first feature", ylab="second feature", main="LDA (three classes)")
 abline(0, eigenvectors[2,1]/eigenvectors[1,1], col="black", lwd=3)
 abline(0, eigenvectors[2,2]/eigenvectors[1,2], col="red", lwd=3)
 axis(side = 1, at = seq(-15, 25,5))
 axis(side = 2, at = seq(-10, 25,5))
 points(c2, col="darkgreen")
 points(c3, col="darkblue")
 points(mu1[1], mu1[2], pch=19, cex=1.5, col="darkblue")
 points(mu2[1], mu2[2], pch=19, cex=1.5, col="red")
 points(mu3[1], mu3[2], pch=19, cex=1.5, col="yellow")
 legend("topright", legend=c("c1", "c2", "c3"), pch=1, col=c("red","darkgreen", "darkblue"))
 text(18,4, "eigenvector 1", cex=0.6, col="black")
 text(-9,4, "eigenvector 2", cex=0.6, col="red")
 
 y1 <- mvtnorm::dmvnorm(c1,mean1,s1)
 plot(c1 %*% eigenvectors[,1],y1, xlim=c(-5, 25), ylim=c(0, .08), col="red", 
      xlab = "projected values", ylab = "",
      main="Multivariate Normal PDF for Eigenvector 1", cex.main=0.9)
 legend("topright", legend=c("c1", "c2", "c3"), pch=1, col=c("red","darkgreen", "darkblue"))
 
 y2 <- mvtnorm::dmvnorm(c2,mean2,s2)
 points(c2 %*% eigenvectors[,1],y2, col="darkgreen")
 
 y3 <- mvtnorm::dmvnorm(c3,mean3,s3)
 points(c3 %*% eigenvectors[,1],y3, col="darkblue")
 
 y1 <- mvtnorm::dmvnorm(c1,mean1,s1)
 plot(c1 %*% eigenvectors[,2],y1, xlim=c(-5, 25), ylim=c(0, .08), col="red",
      xlab = "projected values", ylab = "",
      main="Multivariate Normal PDF for Eigenvector 2", cex.main=0.9)
 legend("topright", legend=c("c1", "c2", "c3"), pch=1, col=c("red","darkgreen", "darkblue"))
 
 y2 <- mvtnorm::dmvnorm(c2,mean2,s2)
 points(c2 %*% eigenvectors[,2],y2, col="darkgreen")
 
 y3 <- mvtnorm::dmvnorm(c3,mean3,s3)
 points(c3 %*% eigenvectors[,2],y3, col="darkblue")
 
 return(list(eigenvectors=eigenvectors,c1=c1,c2=c2,c3=c3,s1=s1,s2=s2,s3=s3,mean1=mean1,mean2=mean2,mean3=mean3))

}