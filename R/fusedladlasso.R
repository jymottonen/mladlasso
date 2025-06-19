#' Multivariate fused LAD-lasso
#'
#' \code{fusedladlasso} is used to fit the multivariate fused LAD-lasso regression model. 
#'
#' @param Y an \eqn{n\times q} matrix of responses. The \eqn{i}th row contains the \eqn{q}-variate response
#' of the \eqn{i}th individual.
#' @param X an \eqn{n\times p} matrix of \eqn{p} explaining variables, The \eqn{i}th row contains the values
#' of \eqn{p} explaining variables for the \eqn{i}th individual.
#' @param initialB a \eqn{(p+1)\times q} matrix of initial regression coefficients
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
#' out1<-fusedladlasso(Y,X,lambda1=0,lambda2=0)
#' plot(out1)
#' out2<-fusedladlasso(Y,X,lambda1=0.2,lambda2=0)
#' plot(out2)
#' out3<-fusedladlasso(Y,X,lambda1=0.2,lambda2=0.2)
#' plot(out3)
#' out <-lambda1.cv(Y,X,lambda1.min = 0.0001,lambda1.max = 0.1,len1=10,lambda2 = 0)
#' out
#' }
#' @importFrom stats optim optimize rnorm
#' @importFrom MASS ginv
#' @import SpatialNP
#' @export
fusedladlasso<-function(Y, X, initialB=NULL, lambda1=0, lambda2=0, 
                        method="BFGS",gradient=FALSE,functional=0,
                        reltol=1e-8,trace=0)
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
  
  if(lambda1>0)
  {
    y1<-matrix(0,p,q)
    y<-rbind(Y,y1)   
    x<-cbind(1,X)
    x1<-cbind(0,n*lambda1*diag(p))
    x<-rbind(x,x1)  
  }
  else if(lambda1==0)
  {
    y<-Y
    x<-cbind(1,X)
  }
  else
    stop("lambda1 should be a non-negative number")
  
  if((functional==0)&(lambda2>0))
  {
    Y2<-matrix(0,p-1,q) 
    y<-rbind(y,Y2)   
    W<-cbind(0,diag(p-1))-cbind(diag(p-1),0)
    X2<-cbind(0,n*lambda2*W)
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
  
  if((functional>0)&(lambda2>0)&gradient)
    gradfn<-dfn
  else
    gradfn<-NULL
  
  if(functional==0|lambda2==0)
  {
    mod<-l1.fit(y,x,initialB,maxiter = 20000,eps = 1e-6, eps.S = 1e-6)
    beta<-as.matrix(mod$coefficients)
    resid<-Y0-cbind(1,X0)%*%beta
    value<-mod$value
    convergence<-mod$convergence
    iter<-mod$iter
  }
  else 
  {
    res<-optim(beta0, fn, gr=gradfn, method=method,
               control=list(maxit=100,reltol=reltol,trace=trace,REPORT=1), 
               Y=y, X=x, lambda1=lambda1, lambda2=lambda2)
    beta<-matrix(res$par,p+1,q)
    resid<-Y0-cbind(1,X0)%*%beta
    value<-res$value
    convergence<-res$convergence
    iter<-NULL
  }
  runt=proc.time()[[3]]-begt
  
  if(is.null(colnames(Y0)))
    colnames(Y0)<-paste("y",1:q,sep="")
  if(is.null(colnames(X0)))
    colnames(X0)<-paste("x",1:p,sep="")
  
  rownames(beta)<-c("Int",colnames(X0))
  colnames(beta)<-colnames(Y0)
  fit<-list(beta=beta,residuals=resid,lambda1=lambda1,lambda2=lambda2,runtime=runt,convergence=convergence,iter=iter,value=value)
  class(fit) <- "fusedladlasso"
  return(fit)
}






