#' Bootstrapping adaptive fused LAD-lasso
#'
#' \code{bootladlasso} is used to  
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
#' \code{\link{lambda1f.cv}} for cross-validation of lambda1. 
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
bootladlasso<-function(Y, X, initialB=NULL, lambda1=0, lambda2=0, K=10)
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

  out<-adfusedladlasso(Y,X,lambda1=lambda1,lambda2=lambda2,K=3)
  BB<-out$B[3,]
  beta<-matrix(BB,p+1,q)
  E<-Y-cbind(1,X)%*%beta

  B<-100
  BB<-matrix(0,B,(p+1)*q)  
  for(i in 1:B)
  {
    print(paste("i=",i))
    sam<-sample(1:n,n,replace=TRUE)
    Estar<-E[sam,]
    Ystar<-cbind(1,X)%*%beta+Estar
    out<-adfusedladlasso(Ystar,X,lambda1=0.03,lambda2=1,K=3)
    BB[i,]<-out$B[3,]
  }
  
  return(BB)
}






