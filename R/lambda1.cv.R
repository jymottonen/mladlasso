#' Selection of \eqn{\lambda_1} using k-fold cross-validation
#'
#' lambda1.cv is used to find the tuning parameter  \eqn{\lambda_1} of the lasso penalty by using k-fold 
#' cross-validation. The tuning parameter \eqn{\lambda_2} of the fusion penalty  is fixed.
#'
#' @param Y an nxq matrix.
#' @param X an nxp matrix.
#' @param lad logical. If lad=TRUE, fused LAD-lasso is used, otherwise fused lasso is used.
#' @param lambda1.min the minimum value of the grid of \eqn{\lambda_1}'s.
#' @param lambda1.max the maximum value of the grid of \eqn{\lambda_1}'s.
#' @param len1 the number of values in the grid of \eqn{\lambda_1}'s
#' @param lambda2 the (fixed) value of the tuning parameter \eqn{\lambda_2}.
#' @param functional functional penalty. If functional=0, then the functional penalty is 
#' \eqn{\lambda_2\sum_{j=2}^{p}||\beta_{j}-\beta_{j-1}||}.
#' If functional=1, then the functional penalty
#' is \eqn{\lambda_2\sum_{j=1}^p\sum_{k=2}^q|\beta_{j,k}-\beta_{j,k-1}|}. If functional=2, 
#' then the functional penalty is \eqn{\lambda_2\sum_{k=2}^q||\beta^{(k)}-\beta^{(k-1)}||}.
#' @param k the number of folds. 
#' @details 
#' Here are the details of the function...
#' @return A list with the following components
#' \describe{
#' \item{lambda1}{the grid of \eqn{\lambda_1}'s. The length of lambda1 is len1}
#' \item{lambda2}{the tuning parameter \eqn{\lambda_2}.}
#' \item{cv}{list of values of CV precision measures in the grid points.
#' The first component vector gives the mean absolute errors (MAE) and the  
#' second component vector gives the mean squared errors (MSE).}
#' \item{lbdmin}{vector of values of \eqn{\lambda_1} that minimize 
#' the precision measures. The ith component 
#' gives the value of \eqn{\lambda_1} that minimizes the ith 
#' precision measure.}
#' \item{tpoint}{vector of types of values of \eqn{\lambda_1} that minimize 
#' the precision measures: tpoint\[i\]=1, if the value of \eqn{\lambda_1} that
#' minimizes the ith precision measure is lambda1.min, 
#' tpoint\[i\]=2, if the value of \eqn{\lambda_1} that minimizes the ith 
#' precision measure is a turning point,  tpoint\[i\]=3, if the value 
#' of \eqn{\lambda_1} that minimizes the ith precision measure is lambda1.max.}
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
#' out <-lambda1.cv(Y,X,lambda1.min = 0.0001,lambda1.max = 0.1,len1=10,lambda2 = 0)
#' plot(out)
#' }
#' @importFrom stats median
#' @export
lambda1.cv<-function(Y,X,lad=TRUE,lambda1.min=0,lambda1.max=5,
                     len1=10,lambda2=0,functional=0,k=5)
{
  #
  # k-fold cross-validation for lambda1 with lambda2 fixed
  #
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
  jakoj<-n%%k
  kok.osa<-floor(n/k)
  m<-rep(kok.osa,k)
  if(jakoj>0){
    m[1:jakoj]<-m[1:jakoj]+1
  }
  
  cv.mae<-rep(0,len1)
  cv.mse<-rep(0,len1)
  h<-rep(0,len1)
  mae<-rep(0,k)
  mse<-rep(0,k)
  
  initialB<-NULL
  groups<-rep(1:k,times=m)
  print(groups)
  for(i1 in 1:len1)
  {
    print(paste("i1=",i1))
    for(j in 1:k)
    {   
      print(paste("j=",j))
      if(lad)
      {
        mod1<-fusedladlasso(Y[groups!=j,],X[groups!=j,],initialB = initialB,
                            lambda1=lbd1[i1],lambda2=lambda2,functional=functional)
        if(is.null(initialB))initialB<-mod1$beta
      }
      else
      {
        mod1<-fusedlasso(Y[groups!=j,],X[groups!=j,],
                        lambda1=lbd1[i1],lambda2=lambda2)
      }
      beta<-mod1$beta
      initB<-beta
      E<-Y[groups==j,]-cbind(1,X[groups==j,])%*%beta
      mae[j]<-mean(sqrt(diag(E%*%t(E))))
      #mae[j]<-median(sqrt(diag(E%*%t(E))))
      mse[j]<-mean(diag(E%*%t(E)))
    }
    cv.mae[i1]<-mean(mae)
    #cv.mae[i1]<-median(mae)
    cv.mse[i1]<-mean(mse)
    if(lad)
    {
      mod1<-fusedladlasso(Y,X,lambda1=lbd1[i1],lambda2=lambda2,
                          functional=functional,initialB = initialB)
    }
    else
    {
      mod1<-fusedlasso(Y,X,lambda1=lbd1[i1],lambda2=lambda2)
    }
    beta<-mod1$beta  
    h[i1]<-sum(abs(beta[-1,])>1.0e-8)
  }
  ind.min.mae<-which.min(cv.mae)
  ind.min.mse<-which.min(cv.mse)
  lbdmin.mae<-lbd1[ind.min.mae]
  lbdmin.mse<-lbd1[ind.min.mse]
  if(ind.min.mae==1)tpoint.mae<-1
  else if(ind.min.mae==length(cv.mae))tpoint.mae<-3
  else tpoint.mae<-2
  if(ind.min.mse==1)tpoint.mse<-1
  else if(ind.min.mse==length(cv.mse))tpoint.mse<-3
  else tpoint.mse<-2
  options(warn=warn.init)
  out<-list(lambda1=lbd1,lambda2=lambda2,cv=list(cv.mae,cv.mse),
            lbdmin=c(lbdmin.mae,lbdmin.mse),
            tpoint=c(tpoint.mae,tpoint.mse),h=h)
  class(out) <- "cv"
  return(out)
}
