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
#' @param lambda1 the tuning parameter \eqn{\lambda_1} for the lasso penalty.
#' @param lambda2 the tuning parameter\eqn{\lambda_2} for the functional penalty.
#' @param lpen gives the lasso penalized coefficients. For example, lpen=c(2,5:8)
#' means that the coefficient vectors \eqn{\beta_2, \beta_5,...,\beta_8} are penalized.
#' @param method the optimization method to be used. The choices are "BFGS", "Nelder-Mead", 
#' "CG", "L-BFGS-B", "SANN", "Brent" The default method is "BFGS". See the Details of the 
#' function optim in package stats.
#' @param gradient a logical evaluating to TRUE or FALSE indicating whether gradient is used when method="BFGS".
#' @param functional functional penalty.
#' @param reltol Relative convergence tolerance of the function optim.
#' @param trace Non-negative integer. If positive, tracing information on the progress of the optimization is produced.
#' @details 
#' Here are the details of the function...
#' @return A list containing the following components:
#' \describe{
#' \item{beta}{the functional LAD regression coefficient matrix.}
#' \item{residuals}{the residuals.}
#' \item{lambda1}{the tuning parameter  \eqn{\lambda_1} for the LAD-lasso penalty}
#' \item{lambda2}{the tuning parameter  \eqn{\lambda_2} for the functional penalty}
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
#' @importFrom stats optim rnorm
#' @importFrom MASS ginv
#' @import SpatialNP
#' @export
functionalladlasso<-function(Y, X, initialB=NULL, lambda1=0, lambda2=0, 
                             lpen=1:dim(X)[2], method="BFGS", gradient=FALSE,functional=1,
                             reltol=1e-9,trace=0)
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

  penalty1<-function(B)
  {
    p<-nrow(B)-1
    W<-cbind(0,diag(p)[lpen,])
    B1<-W%*%B
    sum(sqrt(diag(B1%*%t(B1))))
  }
  
  penalty2<-function(B)
  {
    p<-nrow(B)-1
    q<-ncol(B)
    mat1<-cbind(0,diag(p))
    mat2<-rbind(0,diag(q-1))-rbind(diag(q-1),0)
    sum(abs(mat1%*%B%*%mat2))
  }
 
  penalty3<-function(B)
  {
    p<-nrow(B)-1
    q<-ncol(B)
    mat1<-cbind(0,diag(p))
    mat2<-rbind(0,diag(q-1))-rbind(diag(q-1),0)
    diff<-mat1%*%B%*%mat2
    sum(sqrt(diag(t(diff)%*%diff)))
  }
  
  if((lambda1==0)&(lambda2>0))
  {
    fn<-function(beta,Y,X,lambda1,lambda2){
      B<-matrix(beta,p+1,q)
      E<-Y-cbind(1,X)%*%B
      lad<-mean(sqrt(diag(E%*%t(E))))
      lad+lambda2*penalty2(B)
    }
    fn2<-function(beta,Y,X,lambda1,lambda2){
      B<-matrix(beta,p+1,q)
      E<-Y-cbind(1,X)%*%B
      lad<-mean(sqrt(diag(E%*%t(E))))
      lad+lambda2*penalty3(B)
    }
    dfn<-function(beta,Y,X,lambda1,lambda2){
      B<-matrix(beta,p+1,q)
      E<-Y-cbind(1,X)%*%B
      norm.E <-  sqrt(rowSums(E^2))
      eps.S<-1e-6
      if (min(norm.E) < eps.S) norm.E <- ifelse(norm.E < eps.S, eps.S, norm.E)
      E.sign <- sweep(E,1,norm.E, "/")
      dlad<- -(1/n)*c(t(cbind(1,X))%*%E.sign)
      #derivative of the functional penalty part 
      mat1<-cbind(0,diag(p))
      mat2<-rbind(0,diag(q-1))-rbind(diag(q-1),0)
      dfunctional<-lambda2*c(sign(mat1%*%B%*%mat2)%*%t(mat2))
      dlad+dfunctional
    }
    dfn2<-function(beta,Y,X,lambda1,lambda2){
      B<-matrix(beta,p+1,q)
      E<-Y-cbind(1,X)%*%B
      norm.E <-  sqrt(rowSums(E^2))
      eps.S<-1e-6
      if (min(norm.E) < eps.S) norm.E <- ifelse(norm.E < eps.S, eps.S, norm.E)
      E.sign <- sweep(E,1,norm.E, "/")
      dlad<- -(1/n)*c(t(cbind(1,X))%*%E.sign)
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
  else if((lambda1>0)&(lambda2==0))
  {
    fn<-function(beta,Y,X,lambda1,lambda2){
      B<-matrix(beta,p+1,q)
      E<-Y-cbind(1,X)%*%B
      lad<-mean(sqrt(diag(E%*%t(E))))
      lad+lambda1*penalty1(B)
    }
    fn2<-fn
    dfn<-function(beta,Y,X,lambda1,lambda2){
      B<-matrix(beta,p+1,q)
      #derivative of the lad part
      E<-Y-cbind(1,X)%*%B
      norm.E <-  sqrt(rowSums(E^2))
      eps.S<-1e-6
      if (min(norm.E) < eps.S) norm.E <- ifelse(norm.E < eps.S, eps.S, norm.E)
      E.sign <- sweep(E,1,norm.E, "/")
      dlad<- -(1/n)*c(t(cbind(1,X))%*%E.sign)
      #derivative of the lasso penalty part
      norm.B <-  sqrt(rowSums(B^2))
      if (min(norm.B) < eps.S) norm.B <- ifelse(norm.B < eps.S, eps.S, norm.B)
      B.sign <- sweep(B,1,norm.B, "/")
      B.sign[1,]<-0
      dlasso<-lambda1*c(B.sign)
      dlad+dlasso
    }
  }
  else if((lambda1>0)&(lambda2>0))
  {
    fn<-function(beta,Y,X,lambda1,lambda2){
      B<-matrix(beta,p+1,q)
      E<-Y-cbind(1,X)%*%B
      lad<-mean(sqrt(diag(E%*%t(E))))
      lad+lambda1*penalty1(B)+lambda2*penalty2(B)
    }
    fn2<-function(beta,Y,X,lambda1,lambda2){
      B<-matrix(beta,p+1,q)
      E<-Y-cbind(1,X)%*%B
      lad<-mean(sqrt(diag(E%*%t(E))))
      lad+lambda1*penalty1(B)+lambda2*penalty3(B)
    }
    dfn<-function(beta,Y,X,lambda1,lambda2){
      B<-matrix(beta,p+1,q)
      #derivative of the lad part
      E<-Y-cbind(1,X)%*%B
      norm.E <-  sqrt(rowSums(E^2))
      eps.S<-1e-6
      if (min(norm.E) < eps.S) norm.E <- ifelse(norm.E < eps.S, eps.S, norm.E)
      E.sign <- sweep(E,1,norm.E, "/")
      dlad<- -(1/n)*c(t(cbind(1,X))%*%E.sign)
      #derivative of the lasso penalty part
      norm.B <-  sqrt(rowSums(B^2))
      if (min(norm.B) < eps.S) norm.B <- ifelse(norm.B < eps.S, eps.S, norm.B)
      B.sign <- sweep(B,1,norm.B, "/")
      B.sign[1,]<-0
      dlasso<-lambda1*c(B.sign)
      #derivative of the functional penalty part
      mat1<-cbind(0,diag(p))
      mat2<-rbind(0,diag(q-1))-rbind(diag(q-1),0)
      diff<-mat1%*%B%*%mat2
      norm.diff <-  sqrt(colSums(diff^2))
      if (min(norm.diff) < eps.S) norm.diff <- ifelse(norm.diff < eps.S, eps.S, norm.diff)
      W <- sweep(diff,2,norm.diff, "/")
      dfunctional<-lambda2*rbind(0,W%*%t(mat2))
      dlad+dlasso+dfunctional
    }
  }
  else if((lambda1==0)&(lambda2==0))
  {
    fn<-function(beta,Y,X,lambda1,lambda2){
      B<-matrix(beta,p+1,q)
      E<-Y-cbind(1,X)%*%B
      lad<-mean(sqrt(diag(E%*%t(E))))
      lad
    }
    dfn<-function(beta,Y,X,lambda1,lambda2){
      B<-matrix(beta,p+1,q)
      #derivative of the lad part
      E<-Y-cbind(1,X)%*%B
      norm.E <-  sqrt(rowSums(E^2))
      eps.S<-1e-6
      if (min(norm.E) < eps.S) norm.E <- ifelse(norm.E < eps.S, eps.S, norm.E)
      E.sign <- sweep(E,1,norm.E, "/")
      -(1/n)*c(t(cbind(1,X))%*%E.sign)
    }
  }
  else
    stop("lambda1 and lambda2 should be non-negative numbers")
  
  begt=proc.time()[[3]]
  if(is.null(initialB)){
    #B0<-rnorm((p+1)*q)
    X1<-cbind(1,X)
    B0<-ginv(t(X1)%*%X1)%*%t(X1)%*%Y
    beta0<-c(B0)
  }
  else{
    B0<-initialB
    beta0<-c(B0)
  }
  
  if(functional==2)
  {
    fn<-fn2
    dfn<-dfn2
  }
  
  if(gradient)
    gradfn<-dfn
  else
    gradfn<-NULL
  
  res<-optim(beta0, fn, gr=gradfn, method=method,
             control=list(maxit=100000,reltol=reltol,trace=trace), 
             Y=Y, X=X, lambda1=lambda1, lambda2=lambda2)
  beta<-matrix(res$par,p+1,q)
  resid<-Y-cbind(1,X)%*%beta
  value<-res$value
  convergence<-res$convergence
  runt=proc.time()[[3]]-begt
  rownames(beta)<-c("Int",colnames(X))
  colnames(beta)<-colnames(Y)
  fit<-list(beta=beta,residuals=resid,lambda1=lambda1,lambda2=lambda2,runtime=runt,convergence=convergence,value=value)
  class(fit) <- "functionalladlasso"
  return(fit)
}






