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
#' @param functional functional penalty. If functional=1, then the functional penalty
#' is \eqn{\lambda_2\sum_{j=1}^p\sum_{k=2}^q|\beta_{j,k}-\beta_{j,k-1}|}. If functional=2, 
#' then the functional penalty is \eqn{\lambda_2\sum_{k=2}^q||\beta^{(k)}-\beta^{(k-1)}||}.
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
#' @importFrom stats optim optimize rnorm
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
  if(lambda1>0)
  {
    y1<-matrix(0,p,q)
    y1<-y1[lpen,]
    y<-rbind(Y,y1)   
    x<-cbind(1,X)
    x1<-cbind(0,n*lambda1*diag(p))
    x1<-x1[lpen,]
    x<-rbind(x,x1)  
  }
  else if(lambda1==0)
  {
    y<-Y
    x<-cbind(1,X)
  }
  else
    stop("lambda1 should be a non-negative number")
  
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
  
  if(lambda2>0)
  {
    fn1<-function(beta,Y,X,lambda1,lambda2){
      B<-matrix(beta,p+1,q)
      E<-Y-X%*%B
      lad<-mean(sqrt(diag(E%*%t(E))))
      lad+lambda2*penalty1(B)
    }
    fn2<-function(beta,Y,X,lambda1,lambda2){
      B<-matrix(beta,p+1,q)
      E<-Y-X%*%B
      lad<-mean(sqrt(diag(E%*%t(E))))
      lad+lambda2*penalty2(B)
    }
    dfn1<-function(beta,Y,X,lambda1,lambda2){
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
    dfn2<-function(beta,Y,X,lambda1,lambda2){
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
  else if(lambda2==0)
  {
    fn1<-function(beta,Y,X,lambda1,lambda2){
      B<-matrix(beta,p+1,q)
      E<-Y-X%*%B
      mean(sqrt(diag(E%*%t(E))))
    }
    fn2<-fn1
    dfn1<-function(beta,Y,X,lambda1,lambda2){
      B<-matrix(beta,p+1,q)
      #derivative of the lad part
      E<-Y-X%*%B
      norm.E <-  sqrt(rowSums(E^2))
      eps.S<-1e-6
      if (min(norm.E) < eps.S) norm.E <- ifelse(norm.E < eps.S, eps.S, norm.E)
      E.sign <- sweep(E,1,norm.E, "/")
      -(1/n)*c(t(X)%*%E.sign)
    }
    dfn2<-dfn1
  }
  else
    stop("lambda2 should be a non-negative number")
  
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
  
  if(functional==1)
  {
    fn<-fn1
    dfn<-dfn2
  }
  else if(functional==2)
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
             Y=y, X=x, lambda1=lambda1, lambda2=lambda2)
  beta<-matrix(res$par,p+1,q)
  resid<-Y-cbind(1,X)%*%beta
  value<-res$value
  convergence<-res$convergence
  runt=proc.time()[[3]]-begt
  
  #delta<-0.00001
  #N<-100
  ##beta2.1<-seq(beta[2,1]-delta,beta[2,1]+delta,length=N)
  #beta2.1<-seq(0.20262320,0.20262325,length=N)
  #betai<-beta
  #val<-NULL
  #for(i in 1:N)
  #{
  #  betai[2,1]<-beta2.1[i]
  #  val[i]<-fn(betai,y,x,lambda1,lambda2)
  #}
  #f<-function(beta2.1,beta,y,x,lambda1,lambda2)
  #{
  #  betai<-beta
  #  betai[2,1]<-beta2.1
  #  val[i]<-fn(betai,y,x,lambda1,lambda2)
  #}
  #delta<-0.00001
  #minb<-beta[2,1]-delta
  #maxb<-beta[2,1]+delta
  #opt<-optimize(f,c(minb,maxb),tol=10e-10,beta=beta,y=y,x=x,lambda1=lambda1,lambda2=lambda2)
  #print(opt)
  #plot(beta2.1,val,type="l",ylim=c(min(opt$objective,val),max(val)))
  #abline(h=opt$objective,col="blue",lty=2)
  #abline(v=opt$minimum,col="blue",lty=2)
  #abline(h=res$value,col="red",lty=3)
  #abline(v=beta[2,1],col="red",lty=3)
  
  rownames(beta)<-c("Int",colnames(X))
  colnames(beta)<-colnames(Y)
  fit<-list(beta=beta,residuals=resid,lambda1=lambda1,lambda2=lambda2,runtime=runt,convergence=convergence,value=value)
  class(fit) <- "functionalladlasso"
  return(fit)
}






