'***********************************
Author:  Jack Young, November 2020 *
************************************
The following R code is based upon MATLAB code created by Shireen Elhabian and Aly A. Farag,
"A Tutorial on Data Reduction:  Linear Discriminant  Analysis (LDA)," University of Louisville, 
CVIP Lab, September 2009, http://www.sci.utah.edu/~shireen/pdfs/tutorials/Elhabian_LDA09.pdf 
***********************************'
# Example usage:
# l <- processtwoclasseslda()
# Projecting points onto eigenvector
# y1_w1 <- as.vector(l$x1 %*% l$eigenvectors[,1]); y1_w1
# y2_w1 <- as.vector(l$x2 %*% l$eigenvectors[,1]); y2_w1

processtwoclasseslda <- function()
{
  # Create dataset for class 1
  x1 <- rbind(c(4,2), c(2,4),c(2,3),c(3,6),c(4,4))
  # Create dataset for class 2
  x2 <- rbind(c(9,10), c(6,8),c(9,5),c(8,7),c(10,8))
  # Calculate means for each class
  mu1 <- apply(x1, 2,mean)
  mu2 <- apply(x2, 2,mean)
  # Generate covariance matrix (s1) from x1
  x1_11 <- sum((x1[,1] - mu1[1]) * (x1[,1] - mu1[1]))/(nrow(x1) -1)
  x1_12 <- sum((x1[,1] - mu1[1]) * (x1[,2] - mu1[2]))/(nrow(x1) -1)
  x1_21 <- x1_12
  x1_22 <- sum((x1[,2] - mu1[2]) * (x1[,2] - mu1[2]))/(nrow(x1) -1)
  s1 <- rbind(c(x1_11, x1_12),c(x1_21, x1_22))
  
  # Generate covariance matrix (s2) from x2
  x2_11 <- sum((x2[,1] - mu2[1]) * (x2[,1] - mu2[1]))/(nrow(x2) -1)
  x2_12 <- sum((x2[,1] - mu2[1]) * (x2[,2] - mu2[2]))/(nrow(x2) -1)
  x2_21 <- x2_12
  x2_22 <- sum((x2[,2] - mu2[2]) * (x2[,2] - mu2[2]))/(nrow(x2) -1)
  s2 <- rbind(c(x2_11, x2_12),c(x2_21, x2_22))
  sb <- (mu1 - mu2)%*%t(mu1 - mu2)
  sw <- s1 + s2
  
  # LDA projection via the solution of the generalized eigen value problem:  
  ldaprojection <- eigen(solve(sw)%*%sb)
  eigenvals <- ldaprojection$values
  eigenvectors <- ldaprojection$vectors
  
  # Visual two-class data with LDA projection line for eigenvector #1
  plot(x1[,1], x1[,2], xlim=c(0,10), ylim=c(0,10), col=c("red"),
       xlab="first feature", ylab="second feature", main=paste0("Linear Discriminant  Analysis (LDA):\n", "Two Classes"))
  points(x2[,1], x2[,2], col=c("blue"))
  abline(0, eigenvectors[2,1]/eigenvectors[1,1], col="black", lwd=3)
  text(8,2, "LDA projection line for \n eigenvector 1", cex=0.8, col="red")
  legend("topleft", legend=c("x1", "x2"), pch=1, col=c("red","blue"))
  
  # Print out info about calcuations
  cat("Category 1 data (x1):\n", x1[1,], "\n", x1[2,], "\n", x1[3,], "\n", x1[4,], "\n", x1[5,], "\n")
  cat("Category 2 data (x1):\n", x2[1,], "\n", x2[2,], "\n", x2[3,], "\n", x2[4,], "\n", x2[5,], "\n")
  cat("Expected covariance matrix for s1: \n", cov(x1)[1,], "\n", cov(x1)[2,], "\n")
  cat("Calculated covariance matrix for s1: \n", s1[1,], "\n", s1[2,], "\n")
  cat("Expected covariance matrix for s2: \n", cov(x2)[1,], "\n", cov(x2)[2,], "\n")
  cat("Calculated covariance matrix for s2: \n", s2[1,], "\n", s2[2,], "\n")
  cat("Between-class scatter matrix for sb: \n", sb[1,], "\n", sb[2,], "\n")
  cat("Within-class scatter matrix for sw (s1 + s2): \n", sw[1,], "\n",sw[2,], "\n")
  cat("Eigenvalues ", "\n", eigenvals[1], "\n", eigenvals[2], "\n")
  cat("Eigenvector #1: \n", eigenvectors[1,1], "\n", eigenvectors[2,1], "\n")
  cat("Eigenvector #2: \n", eigenvectors[1,2], "\n", eigenvectors[2,2], "\n")

return(list(eigenvals=eigenvals,eigenvectors=eigenvectors,x1=x1,x2=x2))
  
}