#' Multivariate smoothed LAD-lasso
#'
#' \code{smoothedladlasso} is used to fit the multivariate fused LAD-lasso regression model. 
#'
#' @param Y an \eqn{n\times q} matrix of responses. The \eqn{i}th row contains the \eqn{q}-variate response
#' of the \eqn{i}th individual.
#' @param X an \eqn{n\times p} matrix of \eqn{p} explaining variables, The \eqn{i}th row contains the values
#' of \eqn{p} explaining variables for the \eqn{i}th individual.
#' @param tp a \eqn{q}-dimensional vector of time points.
#' @param initialB a \eqn{(p+1)\times q} matrix of initial regression coefficients
#' @param lambda the tuning parameter \eqn{\lambda} for the lasso penalty
#' @param h the bandwidth parameter
#' @details 
#' Here are the details of the function...
#' @return A list containing the following components:
#' \describe{
#' \item{beta}{the smoothed LAD-lasso regression coefficient matrix.}
#' \item{residuals}{the residuals.}
#' \item{lambda}{the tuning parameter \eqn{\lambda} for the LAD-penalty}
#' \item{iter}{the number of iterations}
#' \item{runtime}{the runtime of the function.}
#' }
#' @references 
#' Oja, H. (2010), \emph{Multivariate Nonparametric Methods with R. 
#' An Approach Based on Spatial Signs and Ranks}, Springer. 
#' \url{https://dx.doi.org/10.1007/978-1-4419-0468-3}.\cr 
#' \cr
#' Nordhausen, K. and Oja, H. (2011), Multivariate L1 Methods: 
#' The Package MNM, \emph{Journal of Statistical Software}, 
#' \strong{43}, 1-28. \url{https://doi.org/10.18637/jss.v043.i05}.
#' @seealso 
#' \code{\link{lambda1.cv}} for cross-validation of lambda1. 
#' @examples
#' \dontrun{
#' data("simdat")
#' Y<-simdat[,1:2]
#' X<-simdat[,3:52]
#' out1<-smoothedladlasso(Y,X,lambda1=0,lambda2=0)
#' plot(out1)
#' out2<-smoothedladlasso(Y,X,lambda1=0.2,lambda2=0)
#' plot(out2)
#' out3<-smoothedladlasso(Y,X,lambda1=0.2,lambda2=0.2)
#' plot(out3)
#' }
#' @importFrom stats optim optimize rnorm
#' @importFrom MASS ginv
#' @import SpatialNP
#' @export
smoothedladlasso<-function(Y, X, tp, initialB=NULL, lambda=0, h=1)
{
  if(is.data.frame(Y))Y<-as.matrix(Y)
  if(is.data.frame(X))X<-as.matrix(X)
  if(dim(Y)[2]<2)stop("response should be at least 2-dimensional!")
  if(dim(X)[1]!=dim(Y)[1])stop("response matrix Y and design matrix X
                                 should have equal number of rows!")
  warn.init<-options()$warn
  options(warn=-1)
  q<-ncol(Y)     #The number of traits (time points)
  p<-ncol(X)     #The number of explaining variables
  n<-nrow(Y)     #The number of cases   
  Y0<-Y
  X0<-X
  
  if(lambda>0)
  {
    y1<-matrix(0,p,q)
    y<-rbind(Y,y1)   
    x<-cbind(1,X)
    x1<-cbind(0,n*lambda*diag(p))
    x<-rbind(x,x1)  
  }
  else if(lambda==0)
  {
    y<-Y
    x<-cbind(1,X)
  }
  else
    stop("lambda should be a non-negative number")
  
  begt=proc.time()[[3]]
  if(is.null(initialB)){
    #B0<-rnorm((p+1)*q)
    B0<-ginv(t(x)%*%x)%*%t(x)%*%y
    beta0<-c(B0)
  }
  else{
    B0<-initialB
    beta0<-c(B0)
  }
  
  K<-function(tp,s,r,h)
  {
    x<-(tp[s]-tp[r])/h
    (1/sqrt(2*pi))*exp(-x^2/2)
  }
  W <- matrix(0,q,q)
  for(r in 1:q)
  {
    W[,r]<-K(tp,1:q,r,h)
  }
  colsums<-apply(W,2,sum)
  W<-sweep(W,2,colsums,FUN="/")
  y2<-y%*%W
  
  mod<-l1.fit(y2,x,initialB,maxiter = 20000,eps = 1e-6, eps.S = 1e-6)
  beta<-as.matrix(mod$coefficients)
  resid<-Y0-cbind(1,X0)%*%beta
  value<-mod$value
  convergence<-mod$convergence
  iter<-mod$iter

  runt=proc.time()[[3]]-begt
  
  if(is.null(colnames(Y0)))
    colnames(Y0)<-paste("y",1:q,sep="")
  if(is.null(colnames(X0)))
    colnames(X0)<-paste("x",1:p,sep="")
  
  rownames(beta)<-c("Int",colnames(X0))
  colnames(beta)<-colnames(Y0)
  fit<-list(beta=beta,residuals=resid,lambda=lambda,convergence=convergence,runtime=runt,iter=iter)
  class(fit) <- "smoothedladlasso"
  return(fit)
}






