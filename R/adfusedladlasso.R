#' Multivariate adaptive fused LAD-lasso
#'
#' \code{adfusedladlasso} is used to fit the multivariate adaptive fused LAD-lasso regression model. 
#'
#' @param Y an \eqn{n\times q} matrix of responses. The \eqn{i}th row contains the \eqn{q}-variate response
#' of the \eqn{i}th individual.
#' @param X an \eqn{n\times p} matrix of \eqn{p} explaining variables, The \eqn{i}th row contains the values
#' of \eqn{p} explaining variables for the \eqn{i}th individual.
#' @param initialB a \eqn{(p+1)\times q} matrix of initial regression coefficients
#' @param lambda1 the tuning parameter \eqn{\lambda_1} for the lasso penalty
#' @param lambda2 the tuning parameter \eqn{\lambda_2} for the fusion penalty
#' @param K the number of adaptive steps 
#' @details 
#' Here are the details of the function...
#' @return A list containing the following components:
#' \describe{
#' \item{beta}{the fused LAD-lasso regression coefficient matrix.}
#' \item{residuals}{the residuals.}
#' \item{lambda1}{the tuning parameter \eqn{\lambda_1} for the LAD-penalty}
#' \item{lambda2}{the tuning parameter \eqn{\lambda_2} for the fusion-penalty}
#' \item{iter}{the number of iterations}
#' \item{runtime}{the runtime of the function.}
#' \item{convergence}{convergence of the optimization routine. 0 indicates successful completion.}
#' \item{value}{the minimized value of the objective function.}
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
#' out1<-adfusedladlasso(Y,X,lambda1=0,lambda2=0)
#' plot(out1)
#' out2<-adfusedladlasso(Y,X,lambda1=0.2,lambda2=0)
#' plot(out2)
#' out3<-adfusedladlasso(Y,X,lambda1=0.2,lambda2=0.2)
#' plot(out3)
#' }
#' @importFrom stats optim optimize rnorm
#' @importFrom MASS ginv
#' @import SpatialNP
#' @export
adfusedladlasso<-function(Y, X, initialB=NULL, lambda1=0, lambda2=0, K=10)
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
  Y0<-Y
  X0<-X

  adapt<-function(Y,X,lambda1,lambda2,v,w,upd)
  {
    p<-ncol(X)
    q<-ncol(Y)
    n<-nrow(Y)
    if(lambda1>0)
    {
      y1<-matrix(0,p,q)
      y<-rbind(Y,y1)   
      x<-cbind(1,X)
      B1<-cbind(0,diag(v))
      x1<-n*lambda1*B1
      x<-rbind(x,x1)  
    }
    else if(lambda1==0)
    {
      y<-Y
      x<-cbind(1,X)
    }
    else
      stop("lambda1 should be a non-negative number")
    
    if(lambda2>0)
    {
      Y2<-matrix(0,p-1,q) 
      y<-rbind(y,Y2)   
      B2<-cbind(matrix(0,p-1,2),diag(w))-cbind(0,diag(w),0)
      X2<-n*lambda2*B2
      x<-rbind(x,X2)  
    }
    mod<-l1.fit(y,x,maxiter = 20000,eps = 1e-6, eps.S = 1e-6)
    beta<-as.matrix(mod$coefficients)
    if(upd==0)
    {
      diff<-beta[3:(p+1),]-beta[2:p,]
      norms<-sqrt(diag(diff%*%t(diff)))
      w<-1/(norms+1/n)
    }
    else if(upd==1)
    {
      norms<-sqrt(diag(beta[-1,]%*%t(beta[-1,])))
      v<-1/(norms+1/n)
    }
    upd<-1-upd
    list(B=beta,v=v,w=w,upd=upd)
  }
  
  upd<-1
  v<-rep(1,p)
  w<-rep(1,p-1)
  V<-NULL
  W<-NULL
  BB<-NULL
  for(k in 1:K)
  {
    res.adapt<-adapt(Y,X,lambda1=lambda1,lambda2=lambda2,v=v,w=w,upd=upd)
    v<-res.adapt$v
    w<-res.adapt$w
    upd<-res.adapt$upd
    B<-res.adapt$B
    V<-rbind(V,v)
    W<-rbind(W,w)
    BB<-rbind(BB,c(B))
  }
  
  fit<-list(B=BB,p=p,q=q,lambda1=lambda1,lambda2=lambda2,V=V,W=W)
  class(fit) <- "adfusedladlasso"
  return(fit)
}






