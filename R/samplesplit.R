#' Sample split
#'
#' \code{samplesplit} is used to perform the multiple-split test.
#'
#' @param Y an \eqn{n\times q} matrix of responses. The \eqn{i}th row contains the \eqn{q}-variate response
#' of the \eqn{i}th individual.
#' @param X an \eqn{n\times p} matrix of \eqn{p} explaining variables, The \eqn{i}th row contains the values
#' of \eqn{p} explaining variables for the \eqn{i}th individual.
#' @param initialB a \eqn{(p+1)\times q} matrix of initial regression coefficients. The default value is NULL.
#' @param lambda1 the tuning parameter \eqn{\lambda_1} for the lasso penalty
#' @param lambda2 the tuning parameter \eqn{\lambda_2} for the fusion penalty
#' @param method the optimization method to be used when functional=1 or 2 and \eqn{\lambda_2>0}. The choices are "BFGS", "Nelder-Mead", 
#' "CG", "L-BFGS-B", "SANN", "Brent" The default method is "BFGS". See the Details of the 
#' function optim in package stats.
#' @param gradient a logical evaluating to TRUE or FALSE indicating whether gradient is used when method="BFGS".
#' @param functional functional penalty. If functional=0, then the functional penalty is 
#' \eqn{\lambda_2\sum_{j=2}^{p}||\beta_{j}-\beta_{j-1}||}.
#' If functional=1, then the functional penalty
#' is \eqn{\lambda_2\sum_{j=1}^p\sum_{k=2}^q|\beta_{j,k}-\beta_{j,k-1}|}. If functional=2, 
#' then the functional penalty is \eqn{\lambda_2\sum_{k=2}^q||\beta^{(k)}-\beta^{(k-1)}||}.
#' @param reltol Relative convergence tolerance of the function optim.
#' @param trace Non-negative integer. If positive, tracing information on the progress of the optimization is produced.
#' @param N the number of replicates. The default value is 100.
#' @param thres the threshold for p-values to judging significant markers. The default value is 0.05.
#' @details 
#' Here are the details of the function...
#' @return A list containing the following components:
#' \describe{
#' \item{pval}{a p-vector,recording the p-values of each marker.}
#' \item{markersig}{the indices of markers which are judged to be significant by multiple-split test.}
#' }
#' @references 
#' Oja, H. (2010), \emph{Multivariate Nonparametric Methods with R. 
#' An Approach Based on Spatial Signs and Ranks}, Springer. 
#' \url{https://dx.doi.org/10.1007/978-1-4419-0468-3}.\cr 
#' \cr
#' Nordhausen, K. and Oja, H. (2011), Multivariate L1 Methods: 
#' The Package MNM, \emph{Journal of Statistical Software}, 
#' \strong{43}, 1-28. \url{https://doi.org/10.18637/jss.v043.i05}.\cr
#' \cr
#' Li, Z., Mottonen, J. and Sillanpaa, M. J. (2015), A robust multiple-locus 
#' method for quantitative trait locus analysis of non-normally distributed 
#' multiple traits, \emph{Heredity}, \strong{115}, 556-564. \url{https://doi.org/10.1038/hdy.2015.61}.
#' @importFrom stats optim optimize rnorm coef pchisq quantile
#' @importFrom MASS ginv
#' @import SpatialNP
#' @export
samplesplit<-function(Y, X, initialB=NULL, lambda1=0, lambda2=0, 
                    method="BFGS",gradient=FALSE,functional=0,
                    reltol=1e-8,trace=0,N=100,thres=0.05)
{
  if(is.data.frame(Y))Y<-as.matrix(Y)
  if(is.data.frame(X))X<-as.matrix(X)
  if(dim(Y)[2]<2)stop("response should be at least 2-dimensional!")
  if(dim(X)[1]!=dim(Y)[1])stop("response matrix Y and design matrix X
                                 should have equal number of rows!")
  warn.init<-options()$warn
  options(warn=-1)
  
  q<-ncol(Y)     #The number of traits
  p<-ncol(X)     #The number of explaining variables
  n<-nrow(Y)     #The number of cases   
  P<-matrix(1,nrow=p,ncol=N)
  Y0<-Y
  X0<-X
  
  penalty1<-function(B)
  {
    p<-nrow(B)-1
    q<-ncol(B)
    mat1<-cbind(0,diag(p))
    mat2<-rbind(0,diag(q-1))-rbind(diag(q-1),0)
    sum(abs(mat1%*%B%*%mat2))
  }
  
  penalty2<-function(B)
  {
    p<-nrow(B)-1
    q<-ncol(B)
    mat1<-cbind(0,diag(p))
    mat2<-rbind(0,diag(q-1))-rbind(diag(q-1),0)
    diff<-mat1%*%B%*%mat2
    sum(sqrt(diag(t(diff)%*%diff)))
  }
  
  for(i in 1:N)
  {
    fullset<-1:n
    subset <- sample.int(n, size = floor(n/2))   #choose a subsample of individuals with size n/2
    subset2 <- fullset[-subset]
    ns<-length(subset)
    
    Ys <- Y[subset,]
    Xs <- X[subset,]
    
    
    if(lambda1>0)
    {
      y1<-matrix(0,p,q)
      y<-rbind(Ys,y1)   
      x<-cbind(1,Xs)
      x1<-cbind(0,ns*lambda1*diag(p))
      x<-rbind(x,x1)  
    }
    else if(lambda1==0)
    {
      y<-Ys
      x<-cbind(1,Xs)
    }
    else
      stop("lambda1 should be a non-negative number")
    
    if((functional==0)&(lambda2>0))
    {
      Y2<-matrix(0,p-1,q) 
      y<-rbind(y,Y2)   
      W<-cbind(0,diag(p-1))-cbind(diag(p-1),0)
      X2<-cbind(0,ns*lambda2*W)
      x<-rbind(x,X2)  
    }
    
    if((functional==1)&(lambda2>0))
    {
      fn<-function(beta,Y,X,lambda1,lambda2){
        B<-matrix(beta,p+1,q)
        E<-Y-X%*%B
        lad<-mean(sqrt(diag(E%*%t(E))))
        lad+lambda2*penalty1(B)
      }
      dfn<-function(beta,Y,X,lambda1,lambda2){
        B<-matrix(beta,p+1,q)
        E<-Y-X%*%B
        norm.E <-  sqrt(rowSums(E^2))
        eps.S<-1e-6
        if (min(norm.E) < eps.S) norm.E <- ifelse(norm.E < eps.S, eps.S, norm.E)
        E.sign <- sweep(E,1,norm.E, "/")
        dlad<- -(1/n)*c(t(X)%*%E.sign)
        #derivative of the functional penalty part 
        mat1<-cbind(0,diag(p))
        mat2<-rbind(0,diag(q-1))-rbind(diag(q-1),0)
        dfunctional<-lambda2*c(sign(mat1%*%B%*%mat2)%*%t(mat2))
        dlad+dfunctional
      }
    }
    else if((functional==2)&(lambda2>0))
    {
      fn<-function(beta,Y,X,lambda1,lambda2){
        B<-matrix(beta,p+1,q)
        E<-Y-X%*%B
        lad<-mean(sqrt(diag(E%*%t(E))))
        lad+lambda2*penalty2(B)
      }
      dfn<-function(beta,Y,X,lambda1,lambda2){
        B<-matrix(beta,p+1,q)
        E<-Y-X%*%B
        norm.E <-  sqrt(rowSums(E^2))
        eps.S<-1e-6
        if (min(norm.E) < eps.S) norm.E <- ifelse(norm.E < eps.S, eps.S, norm.E)
        E.sign <- sweep(E,1,norm.E, "/")
        dlad<- -(1/n)*c(t(X)%*%E.sign)
        #derivative of the functional penalty part 
        mat1<-cbind(0,diag(p))
        mat2<-rbind(0,diag(q-1))-rbind(diag(q-1),0)
        diff<-mat1%*%B%*%mat2
        norm.diff <-  sqrt(colSums(diff^2))
        if (min(norm.diff) < eps.S) norm.diff <- ifelse(norm.diff < eps.S, eps.S, norm.diff)
        W <- sweep(diff,2,norm.diff, "/")
        dfunctional<-lambda2*rbind(0,W%*%t(mat2))
        dlad+dfunctional
      }
    }
    
    if(is.null(initialB)){
      B0<-ginv(t(x)%*%x)%*%t(x)%*%y
      beta0<-c(B0)
    }
    else{
      B0<-initialB
      beta0<-c(B0)
    }
    
    if((functional>0)&(lambda2>0)&gradient)
      gradfn<-dfn
    else
      gradfn<-NULL
    
    if(functional==0|lambda2==0)
    {
      mod<-l1.fit(y,x,initialB,maxiter = 20000,eps = 1e-6, eps.S = 1e-6)
      beta<-as.matrix(mod$coefficients)
      norms<-sqrt(diag(beta%*%t(beta)))
      nonzeros <- which(norms>1e-4)-1
      nonzeros <- nonzeros[2:length(nonzeros)]
    }
    else 
    {
      res<-optim(beta0, fn, gr=gradfn, method=method,
                 control=list(maxit=100,reltol=reltol,trace=trace,REPORT=1), 
                 Y=y, X=x, lambda1=lambda1, lambda2=lambda2)
      beta<-matrix(res$par,p+1,q)
      norms<-sqrt(diag(beta%*%t(beta)))
      nonzeros <- which(norms>1e-4)-1
      nonzeros <- nonzeros[2:length(nonzeros)]
    }
    
    yt <- Y[subset2,]
    xt <- cbind(1,X[subset2,nonzeros])
    
    #calculate the pvalues of the selected markers based on the second half of the data
    qt <- length(nonzeros)
    for(j in 1:qt){
      x1 <- xt[,-(j+1)]
      x2 <- xt[,j+1]
      mModel<-l1.fit(yt,x1,initialB,maxiter = 20000,eps = 1e-6, eps.S = 1e-6)
      beta <- coef(mModel)
      R <- yt-x1%*%beta
      Th <- spatial.sign2(R, center = FALSE, shape = TRUE)
      PT <- Th%*%solve(t(Th)%*%Th)%*%t(Th)
      PX1 <- x1%*%solve(t(x1)%*%x1)%*%t(x1)
      nt <- dim(yt)[1]
      x2h <- (diag(nt)-PX1)%*%x2
      PX2 <- x2h%*%solve(t(x2h)%*%x2h)%*%t(x2h)
      q <- dim(y)[2] 
      P[nonzeros[j],i] <- min(c(1,(1-pchisq(ns*sum(diag(PT%*%PX2)),q))*qt))
    }
  }
  
  #For each marker, combine the pvalues calculated over N replicates into a single pvalue 
  
  gamma = seq(0.05, 0.99, by = 0.01)
  pval <- numeric(p)
  for (j in 1:p) {
    quant.gamma <- quantile(P[j, ], gamma, na.rm=TRUE)/gamma
    if (length(gamma) > 1) 
      penalty <- (1 - log(min(gamma)))
    else penalty <- 1
    pval.pre <- min(quant.gamma) * penalty
    pval[j] <- pmin(pval.pre, 1)
  }
  
  markersig <- which( pval < thres) #detect significant markers
  fit<-list(pval=pval, markersig = markersig)
  class(fit) <- "samplesplit"
  return(fit)
}






