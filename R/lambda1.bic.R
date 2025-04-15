#' Selection of \eqn{\lambda_1} using BIC-type criterion
#'
#' \code{lambda1.bic} is used to find the tuning parameter \eqn{\lambda_1} of the lasso penalty  by using BIC-type 
#' criterion. The tuning parameter of the fusion penalty \eqn{\lambda_2} is fixed.
#'
#' @param Y an nxq matrix.
#' @param X an nxp matrix.
#' @param lad logical. If lad=TRUE, fused LAD-lasso is used, otherwise fused lasso is used.
#' @param lambda1.min the minimum value of the grid of \eqn{\lambda_1}'s.
#' @param lambda1.max the maximum value of the grid of \eqn{\lambda_1}'s.
#' @param len1 the number of values in the grid of \eqn{\lambda_1}'s
#' @param lambda2 the (fixed) value of the tuning parameter \eqn{\lambda_2}. 
#' @details 
#' Here are the details of the function...
#' @return A list with the following components
#' \describe{
#' \item{lambda1}{the grid of \eqn{\lambda_1}'s. The length of lambda1 is len1.}
#' \item{lambda2}{the tuning parameter \eqn{\lambda_2}.}
#' \item{bic}{vector of values of BIC-type criterion in the grid points.}
#' \item{lbdmin}{the value of lambda1 that minimizes the BIC-type criterion.}
#' \item{sigmahat}{sigmahat}
#' \item{tpoint}{the type of value of \eqn{\lambda_1} that minimizes 
#' the BIC-type criterion: tpoint=1, if the value of \eqn{\lambda_1} that
#' minimizes the BIC-type criterion is lambda1.min, 
#' tpoint=2, if the value of \eqn{\lambda_1} that minimizes the BIC-type criterion is 
#' a turning point,  tpoint=3, if the value 
#' of \eqn{\lambda_1} that minimizes the BIC-type criterion is lambda1.max.}
#' \item{h}{the number of non-zero coefficient vectors.}
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
#' \code{\link{fusedladlasso}} for multivariate fused LAD-lasso. 
#' @examples
#' \dontrun{
#' data("simdat")
#' Y<-simdat[,1:2]
#' X<-simdat[,3:52]
#' out <-lambda1.bic(Y,X,lambda1.min = 0.0001,lambda1.max = 0.1,len1=10,lambda2 = 0)
#' plot(out)
#' }
#' @export
lambda1.bic<-function(Y,X,lad=TRUE,lambda1.min=0,lambda1.max=5,len1=10,lambda2=0)
{
  if(is.data.frame(Y))Y<-as.matrix(Y)
  if(is.data.frame(X))X<-as.matrix(X)
  if(dim(Y)[2]<2)stop("response should be at least 2-dimensional!")
  if(dim(X)[1]!=dim(Y)[1])stop("response matrix Y and design matrix X
                                 should have equal number of rows!")
  warn.init<-options()$warn
  options(warn=-1)
  if(is.null(lambda1.max))
    lambda1.max<-pi*sqrt(log(ncol(X)+1)/nrow(X))
  if(is.null(lambda1.min))
    lambda1.min<-0.3*lambda1.max
  lbd1<-seq(lambda1.min,lambda1.max,length=len1)
  n<-nrow(X)
  q<-ncol(Y)
  p<-ncol(X)
  
  bic<-rep(0,len1)
  h<-rep(0,len1)
  value<-rep(0,len1)
  
  if(lad){
    #if(n>p)
    #  mod0<-fusedladlasso(Y,X,lambda1=0,lambda2=0)
    #else
      mod0<-fusedladlasso(Y,X,lambda1=0.1,lambda2=0)
  }
  else{
    #if(n>p)
    #  mod0<-fusedlasso(Y,X,lambda1=0,lambda2=0)
    #else
      mod0<-fusedlasso(Y,X,lambda1=0.1,lambda2=0)
  }
  beta0<-mod0$beta
  E<-Y-cbind(1,X)%*%beta0

  if(lad){
    const<-sqrt(n/(2*(n-p-1)))*gamma(q/2)/gamma((q+1)/2)
    sigmahat<-const*mean(sqrt(diag(E%*%t(E))))
    print(sigmahat)
    scale<-sigmahat
  }
  else{
    const<-n/(q*(n-p-1))
    sigmahat2<-const*mean(diag(E%*%t(E)))
    sigmahat<-sqrt(sigmahat2)
    scale<-sigmahat2
  }

  initB<-beta0
  for(i1 in 1:len1)
  {
    if(lad){
      mod1<-fusedladlasso(Y,X,initialB = initB,lambda1=lbd1[i1],lambda2=lambda2)
    }
    else{
      mod1<-fusedlasso(Y,X,lambda1=lbd1[i1],lambda2=lambda2)
    }
    beta<-mod1$beta
    initB<-beta
    E<-Y-cbind(1,X)%*%beta
    if(lad){
      value[i1]<-mean(sqrt(diag(E%*%t(E))))
    }
    else{
      value[i1]<-mean(diag(E%*%t(E)))
    }
    h[i1]<-sum(abs(beta[-1,])>1.0e-8)
    bic[i1]<-value[i1]/scale+h[i1]*log(n)/n
    print(paste("i1=",i1," ,lambda1=",lbd1[i1],", h=",h[i1],", bic=",bic[i1]))
  }
  ind.min<-which.min(bic)
  lbdmin<-lbd1[ind.min]
  if(ind.min==1)tpoint<-1
  else if(ind.min==length(bic))tpoint<-3
  else tpoint<-2
  options(warn=warn.init)
  out<-list(lambda1=lbd1,lambda2=lambda2,bic=bic,
            lbdmin=lbdmin,sigmahat=sigmahat,
            tpoint=tpoint,h=h)
  class(out) <- "bic"
  return(out)
}
