#' Cross-validation of lambda1
#'
#' lambda1.cv is used to find the tuning parameter of the lasso penalty lambda1 by using 5-fold 
#' cross-validation. The tuning parameter of the fusion penalty lambda2 is fixed.
#'
#' @param Y an nxq matrix.
#' @param X an nxp matrix.
#' @param lambda1.min the minimum value of the grid of lambda1's.
#' @param lambda1.max the maximum value of the grid of lambda1's.
#' @param len1 the number of values in the grid of lambda1's
#' @param lambda2 the (fixed) value of the tuning parameter lambda2. 
#' @param lpen gives the lasso penalized coefficients. For example, lpen=c(2,5:8)
#' means that the coefficient vectors beta2, beta5,...,beta8 are penalized.
#' @param fpen a list of blocks of fusion penalized coefficients. For example, 
#' fpen=list(2:5,10:20) means that the fusion penalty is 
#' \code{lambda2*[||beta_3-beta_2||+...+||beta_5-beta_4||+
#' ||beta_{11}-beta_{10}||+...+||beta_{20}-beta_{19}||]}
#' @details 
#' Here are the details of the function...
#' @return A list with the following components
#' \describe{
#' \item{lambda1}{the grid of lambda1's. The length lambda1 is len1}
#' \item{lambda2}{the tuning parameter lambda2.}
#' \item{cv}{list of values of CV precision measures. 
#' The ith component vector gives the values of the ith 
#' precision measure in the grid points.}
#' \item{lbdmin}{vector of values of lambda1 that minimize 
#' the precision measures. The ith component 
#' gives the value of lambda1 that minimizes the ith 
#' precision measure.}
#' \item{tpoint}{vector of types of values of lambda1 that minimize 
#' the precision measures: tpoint\[i\]=1, if the value of lambda1 that
#' minimizes the ith precision measure is lambda1.min, 
#' tpoint\[i\]=2, if the value of lambda1 that minimizes the ith 
#' precision measure is a turning point,  tpoint\[i\]=3, if the value 
#' of lambda1 that minimizes the ith precision measure is lambda1.max.}
#' \item{h}{the number non-zero coefficients.}
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
#' X<-simdat[,3:32]
#' out <-lambda1.cv(Y,X,lambda1.min = 0.0001,lambda1.max = 0.1,len1=10,lambda2 = 0)
#' plot(out)
#' }
#' @export
lambda1.cv<-function(Y,X,lambda1.min=0,lambda1.max=5,len1=10,lambda2=0,
                     lpen=1:dim(X)[2],fpen=list(1:dim(X)[2]))
{
  #
  # 5-fold cross-validation for lambda1 with lambda2 fixed
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
  jakoj<-n%%5
  kok.osa<-floor(n/5)
  m<-rep(kok.osa,5)
  if(jakoj>0){
    m[1:jakoj]<-m[1:jakoj]+1
  }
  
  cv1<-rep(0,len1)
  cv2<-rep(0,len1)
  cv3<-rep(0,len1)
  h<-rep(0,len1)
  mse1<-rep(0,5)
  mse2<-rep(0,5)
  mse3<-rep(0,5)
  bisq<-function(E)
  {
    absr<-sqrt(diag(E%*%t(E)))
    k<-2
    1*(absr>k)+(1-(1-(absr/k)^2)^3)*(absr<=k)
  }
  
  groups<-rep(1:5,times=m)  
  for(i1 in 1:len1)
  {
    begt=Sys.time()
    for(j in 1:5)
    {   
      mod1<-fusedladlasso(Y[groups!=j,],X[groups!=j,],
                          lambda1=lbd1[i1],lambda2=lambda2,lpen=lpen,fpen=fpen)
      beta<-mod1$beta
      E<-Y[groups==j,]-cbind(1,X[groups==j,])%*%beta
      mse1[j]<-mean(sqrt(diag(E%*%t(E))))
      mse2[j]<-mean(diag(E%*%t(E)))
      mse3[j]<-mean(bisq(E))
    }
    cv1[i1]<-mean(mse1)
    cv2[i1]<-mean(mse2)
    cv3[i1]<-mean(mse3)            
    mod1<-fusedladlasso(Y,X,lambda1=lbd1[i1],lambda2=lambda2,lpen=lpen,fpen=fpen)
    beta<-mod1$beta  
    norms<-sqrt(diag(beta%*%t(beta)))
    h[i1]<-sum(norms>1e-6)-1
    runt=as.numeric(Sys.time()-begt)
    print(paste("i1=",i1,"lambda1=",lbd1[i1],"lambda2=",lambda2,
                "cv1=",cv1[i1],"cv2=",cv2[i1],"cv3=",cv3[i1],
                "h=",h[i1],"runtime=",runt))
  }
  ind.min1<-which.min(cv1)
  ind.min2<-which.min(cv2)
  ind.min3<-which.min(cv3)
  lbdmin1<-lbd1[ind.min1]
  lbdmin2<-lbd1[ind.min2]
  lbdmin3<-lbd1[ind.min3]
  if(ind.min1==1)tpoint1<-1
  else if(ind.min1==length(cv1))tpoint1<-3
  else tpoint1<-2
  if(ind.min2==1)tpoint2<-1
  else if(ind.min2==length(cv2))tpoint2<-3
  else tpoint2<-2
  if(ind.min3==1)tpoint3<-1
  else if(ind.min3==length(cv3))tpoint3<-3
  else tpoint3<-2
  options(warn=warn.init)
  out<-list(lambda1=lbd1,lambda2=lambda2,cv=list(cv1,cv2,cv3),
            lbdmin=c(lbdmin1,lbdmin2,lbdmin3),
            tpoint=c(tpoint1,tpoint2,tpoint3),h=h)
  class(out) <- "cv"
  return(out)
}