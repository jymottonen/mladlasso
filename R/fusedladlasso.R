#' Multivariate fused LAD-lasso
#'
#' \code{fusedladlasso} is used to fit the multivariate fused LAD-lasso regression model. 
#'
#' @param Y an nxq matrix
#' @param X an nxp matrix
#' @param lambda1 the tuning parameter for the LAD-penalty
#' @param lambda2 the tuning parameter for the fusion-penalty
#' @details 
#' Here are the details of the function...
#' @return A list containing the following components:
#' \describe{
#' \item{beta}{the fused LAD-lasso regression coefficient matrix.}
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
#' X<-simdat[,3:32]
#' out<-fusedladlasso(Y,X,lambda1=0.2,lambda2=0)
#' out$runtime
#' }
#' @export
fusedladlasso<-function(Y,X,lambda1=0,lambda2=0)
{
  if(is.data.frame(Y))Y<-as.matrix(Y)
  if(is.data.frame(X))X<-as.matrix(X)
  if(dim(Y)[2]<2)stop("response should be at least 2-dimensional!")
  if(dim(X)[1]!=dim(Y)[1])stop("response matrix Y and design matrix X
                                 should have equal number of rows!")
  warn.init<-options()$warn
  options(warn=-1)
  RowNorms <- function(X)
  {
    sqrt(rowSums2(X^2))
  }
  p<-ncol(Y)     #The number of traits
  X<-cbind(1,X)  #Add the intercept term
  q<-ncol(X)     #The number of explaining variables + 1
  n<-nrow(Y)     #The number of cases    
  if((lambda1==0)&(lambda2>0))
  {
    Y1<-matrix(0,q-2,p) 
    y<-rbind(Y,Y1)   
    X2<-n*lambda2*(diag(q)[-c(1,2),]-diag(q)[-c(1,q),])
    x<-rbind(X,X2)  
  }
  else if((lambda1>0)&(lambda2==0))
  {
    Y1<-matrix(0,q-1,p) 
    y<-rbind(Y,Y1)   
    X1<-n*lambda1*diag(q)[-1,]
    x<-rbind(X,X1)  
  }
  else if((lambda1>0)&(lambda2>0))
  {
    Y1<-matrix(0,2*q-3,p) 
    y<-rbind(Y,Y1)   
    X1<-n*lambda1*diag(q)[-1,]
    X2<-n*lambda2*(diag(q)[-c(1,2),]-diag(q)[-c(1,q),])
    x<-rbind(X,X1,X2)  
  }
  else if((lambda1==0)&(lambda2==0))
  {
    y<-Y
    x<-X
  }
  else
    stop("lambda1 and lambda2 should be non-negative numbers")
  
  begt=Sys.time()
  mod<-mv.l1lm(y~-1+x,score="s",stand="o",maxiter = 10000,
               eps = 1e-8, eps.S = 1e-8)
  beta<-as.matrix(coefficients(mod))
  runt=as.numeric(Sys.time()-begt)
  rownames(beta)<-c("Int",colnames(X)[-1])
  colnames(beta)<-colnames(Y)
  fit<-list(beta=beta,runtime=runt)
  class(fit) <- "fusedladlasso"
  return(fit)
}






