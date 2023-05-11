#' Multivariate functional LAD-lasso
#'
#' \code{functionalladlasso} is used to fit the multivariate functional LAD-lasso regression model. 
#'
#' @param Y an nxq matrix of responses. The ith row contains the q-variate response
#' of the ith individual.
#' @param X an nxp matrix of p explaining variables, The ith row contains the values
#' of p explaining variables for the ith individual.
#' @param initialB a (p+1)xq matrix of initial regression coefficients. If NULL, the initial values 
#' are generated from standard normal distribution.
#' @param lambda1 the tuning parameter for the LAD-lasso penalty.
#' @param lambda2 the tuning parameter for the functional penalty.
#' @param lpen gives the lasso penalized coefficients. For example, lpen=c(2,5:8)
#' means that the coefficient vectors beta2, beta5,...,beta8 are penalized.
#' @details 
#' Here are the details of the function...
#' @return A list containing the following components:
#' \describe{
#' \item{beta}{the functional LAD regression coefficient matrix.}
#' \item{residuals}{the residuals.}
#' \item{lambda1}{the tuning parameter for the LAD-lasso penalty}
#' \item{lambda2}{the tuning parameter for the functional penalty}
#' \item{runtime}{the runtime of the function.}
#' \item{convergence}{convergence of the optimization routine. 0 indicates successful completion.}
#' \item{value}{the minimized value of the objective function}
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
#' out<-functional.lad(Y,X,lambda1=0.2,lambda2=0.3)
#' out$runtime
#' }
#' @export
functionalladlasso<-function(Y, X, initialB=NULL, lambda1=0, lambda2=0, lpen=1:dim(X)[2])
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
  if(is.null(colnames(Y)))
    colnames(Y)<-paste("y",1:q,sep="")
  if(is.null(colnames(X)))
    colnames(X)<-paste("x",1:p,sep="")

  if((lambda1==0)&(lambda2>0))
  {
    fn<-function(beta,Y,X,lambda1,lambda2){
      B<-matrix(beta,p+1,q)
      E<-Y-cbind(1,X)%*%B
      lad<-mean(sqrt(diag(E%*%t(E))))
      penalty2<-0
      for(i in 2:(p+1))
        for(j in 2:q)
          penalty2<-penalty2+abs(B[i,j]-B[i,j-1])
      lad+lambda2*penalty2
    }
  }
  else if((lambda1>0)&(lambda2==0))
  {
    fn<-function(beta,Y,X,lambda1,lambda2){
      B<-matrix(beta,p+1,q)
      E<-Y-cbind(1,X)%*%B
      lad<-mean(sqrt(diag(E%*%t(E))))
      W<-cbind(0,diag(p)[lpen,])
      B1<-W%*%B
      penalty1<-sum(sqrt(diag(B1%*%t(B1))))
      lad+lambda1*penalty1
    }
  }
  else if((lambda1>0)&(lambda2>0))
  {
    fn<-function(beta,Y,X,lambda1,lambda2){
      B<-matrix(beta,p+1,q)
      E<-Y-cbind(1,X)%*%B
      lad<-mean(sqrt(diag(E%*%t(E))))
      W<-cbind(0,diag(p)[lpen,])
      B1<-W%*%B
      penalty1<-sum(sqrt(diag(B1%*%t(B1))))
      penalty2<-0
      for(i in 2:(p+1))
        for(j in 2:q)
          penalty2<-penalty2+abs(B[i,j]-B[i,j-1])
      lad+lambda1*penalty1+lambda2*penalty2
    }
  }
  else if((lambda1==0)&(lambda2==0))
  {
    fn<-function(beta,Y,X,lambda1,lambda2){
      B<-matrix(beta,p+1,q)
      E<-Y-cbind(1,X)%*%B
      mean(sqrt(diag(E%*%t(E))))
    }
  }
  else
    stop("lambda1 and lambda2 should be non-negative numbers")
  
  #The gradient of fn
  dfn<-function(beta,Y,X){
    B<-matrix(beta,q,p)
    E<-Y-X%*%B
    norm.E<-SpatialNP:::norm(E)
    E.sign<-sweep(E,1,norm.E, "/")
    -(1/n)*c(t(X)%*%E.sign)
  }
  
  begt=Sys.time()
  if(is.null(initialB)){
    #B0<-rnorm((p+1)*q)
    X1<-cbind(1,X)
    B0<-MASS::ginv(t(X1)%*%X1)%*%t(X1)%*%Y
    beta0<-c(B0)
  }
  else{
    B0<-initialB
    beta0<-c(B0)
  }
  
  res<-optim(beta0, fn, gr=NULL, method="BFGS",
             control=list(maxit=10000,reltol=1e-10,trace=1), Y=Y, X=X, lambda1=lambda1, lambda2=lambda2)
  beta<-matrix(res$par,p+1,q)
  value<-res$value
  convergence<-res$convergence
  runt=Sys.time()-begt
  rownames(beta)<-c("Int",colnames(X))
  colnames(beta)<-colnames(Y)
  fit<-list(beta=beta,residuals=res,lambda1=lambda1,lambda2=lambda2,runtime=runt,convergence=convergence,value=value)
  class(fit) <- "functionalladlasso"
  return(fit)
}






