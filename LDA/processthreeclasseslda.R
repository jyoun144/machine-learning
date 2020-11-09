processthreeclasseslda <- function()
{
  set.seed(262)
  size <- 500
  u1 <- cbind(runif(size), runif(size))
  u2 <- cbind(runif(size), runif(size))
 
  x1 <- matrix(nrow=size, ncol=2)
  for (i in 1:size)
    {
    x1[i,1] = sqrt(-2*log(u1[i,1]))*cos(2*pi*u2[i,1])
    x1[i,2] = sqrt(-2*log(u1[i,2]))*cos(2*pi*u2[i,2])
  }
  
  mu <- c(5,5)
  mu1 <- c(mu[1]-3,mu[2]+7)
  mu2 <- c(mu[1]-2.5,mu[2]-3.5)
  mu3 <- c(mu[1]+7,mu[2]+5)
  cov1 <- rbind(c(5,-1),c(-3,3))
  cov2 <- rbind(c(4,0),c(0,4))
  cov3 <- rbind(c(3.5,1),c(3,2.5))
  eigen1 <- eigen(cov1)
  eigen2 <- eigen(cov2)
  eigen3 <- eigen(cov3)
  c1 <- x1
  c2 <- x1
  c3 <- x1
  
 for (i in 1:size)
 {
     c1[i,] <- eigen1$vectors %*% sqrt(diag(eigen1$values)) %*%x1[i,]
 }
  for (i in 1:size)
  {
    c2[i,] <- eigen2$vectors %*% sqrt(diag(eigen2$values)) %*%x1[i,]
  }
  for (i in 1:size)
  {
    c3[i,] <- eigen3$vectors %*% sqrt(diag(eigen3$values)) %*%x1[i,]
  }
 c1[,1] <- c1[,1] + mu1[1]
 c1[,2] <- c1[,2] + mu1[2]
 
 c2[,1] <- c2[,1] + mu2[1]
 c2[,2] <- c2[,2] + mu2[2]
 
 c3[,1] <- c3[,1] + mu3[1]
 c3[,2] <- c3[,2] + mu3[2]
 
 mean1 <- apply(c1, 2,mean)
 mean2 <- apply(c2, 2,mean)
 mean3 <- apply(c3, 2,mean)
 
 overallmean <- (mean1 + mean2 + mean3)/3
 
 s1 <- cov(c1)
 s2 <- cov(c2)
 s3 <- cov(c3)
 sw <- s1 + s2 + s3
 
 n1 <- nrow(c1)
 n2 <- nrow(c2)
 n3 <- nrow(c3)
 
 sb1 <- n1*(mean1 - overallmean)%*%t(mean1 - overallmean)
 sb2 <- n2*(mean2 - overallmean)%*%t(mean2 - overallmean)
 sb3 <- n3*(mean3 - overallmean)%*%t(mean3 - overallmean)
 sb <- sb1 + sb2 + sb3
 
 ldaprojection <- eigen(solve(sw)%*%sb)
 eigenvals <- ldaprojection$values
 eigenvectors <- ldaprojection$vectors
 
 # Visualize Results
 plot(c1, xlim=c(-15,25), ylim=c(-10,25), col="red", axes = FALSE, panel.first = grid())
 abline(0, eigenvectors[2,1]/eigenvectors[1,1], col="black", lwd=3)
 axis(side = 1, at = seq(-15, 25,5))
 axis(side = 2, at = seq(-10, 25,5))
 points(c2, col="darkgreen")
 points(c3, col="darkblue")
 points(mu1[1], mu1[2], pch=19, cex=1.5, col="darkblue")
 points(mu2[1], mu2[2], pch=19, cex=1.5, col="red")
 points(mu3[1], mu3[2], pch=19, cex=1.5, col="yellow")
 
 return(ldaprojection)
}