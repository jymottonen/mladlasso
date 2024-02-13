#' Multivariate fused LAD-lasso
#'
#' \code{fusedladlasso} is used to fit the multivariate fused LAD-lasso regression model. 
#'
#' @param Y an nxq matrix of responses. The ith row contains the q-variate response
#' of the ith individual.
#' @param X an nxp matrix of p explaining variables, The ith row contains the values
#' of p explaining variables for the ith individual.
#' @param initialB initial value of the coefficient matrix
#' @param lambda1 the tuning parameter \eqn{\lambda_1} for the lasso penalty
#' @param lambda2 the tuning parameter \eqn{\lambda_2} for the fusion penalty
#' @param lpen gives the lasso penalized coefficients. For example, lpen=c(2,5:8)
#' means that the coefficient vectors \eqn{\beta_2, \beta_5,...,\beta_8} are penalized.
#' @param fpen a list of blocks of fusion penalized coefficients. For example, 
#' fpen=list(2:5,10:20) means that the fusion penalty is 
#' \eqn{\lambda_2[||\beta_3-\beta_2||+...+||\beta_5-\beta_4||+||\beta_{11}-\beta_{10}||+...+||\beta_{20}-\beta_{19}||}
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
#' @export
fusedladlasso<-function(Y, X, initialB=NULL, lambda1=0, lambda2=0, 
                        lpen=1:dim(X)[2], fpen=list(1:dim(X)[2]))
{
  if(is.data.frame(Y))Y<-as.matrix(Y)
  if(is.data.frame(X))X<-as.matrix(X)
  if(dim(Y)[2]<2)stop("response should be at least 2-dimensional!")
  if(dim(X)[1]!=dim(Y)[1])stop("response matrix Y and design matrix X
                                 should have equal number of rows!")
  nblocks<-length(lengths(fpen)) #The number of blocks
  conc<-NULL
  for(blk in 1:nblocks)conc<-c(conc,fpen[[blk]])
  if(anyDuplicated(conc)!=0)stop("blocks: the blocks should not overlap!")
  isect<-intersect(1:dim(X)[2],conc)
  if(length(isect)<length(conc))
    stop(paste("the blocks should contain only integers from 1 to ",dim(X)[2],sep=""))
  warn.init<-options()$warn
  options(warn=-1)
  q<-ncol(Y)     #The number of traits
  p<-ncol(X)     #The number of explaining variables
  n<-nrow(Y)     #The number of cases   
  Y0<-Y
  X0<-X
  if(is.null(colnames(Y)))
    colnames(Y)<-paste("y",1:q,sep="")
  if(is.null(colnames(X)))
    colnames(X)<-paste("x",1:p,sep="")
  lpen<-sort(lpen)
  fused.id<-NULL
  for(i in 1:nblocks)
  {
    id<-min(fpen[[i]]):(max(fpen[[i]])-1)
    fused.id<-c(fused.id,id)
  }
  fused.id<-sort(fused.id)
  if((lambda1==0)&(lambda2>0))
  {
    Y2<-matrix(0,p-1,q) 
    Y2<-Y2[fused.id,]
    y<-rbind(Y,Y2)
    W<-cbind(0,diag(p-1))-cbind(diag(p-1),0)
    X<-cbind(1,X)
    X2<-cbind(0,n*lambda2*W)
    X2<-X2[fused.id,]
    x<-rbind(X,X2)  
  }
  else if((lambda1>0)&(lambda2==0))
  {
    Y1<-matrix(0,p,q)
    Y1<-Y1[lpen,]
    y<-rbind(Y,Y1)   
    X<-cbind(1,X)
    X1<-cbind(0,n*lambda1*diag(p))
    X1<-X1[lpen,]
    x<-rbind(X,X1)  
  }
  else if((lambda1>0)&(lambda2>0))
  {
    Y1<-matrix(0,p,q) 
    Y1<-Y1[lpen,]
    Y2<-matrix(0,p-1,q) 
    Y2<-Y2[fused.id,]
    y<-rbind(Y,Y1,Y2)   
    X<-cbind(1,X)
    X1<-cbind(0,n*lambda1*diag(p))
    X1<-X1[lpen,]
    W<-cbind(0,diag(p-1))-cbind(diag(p-1),0)
    X2<-cbind(0,n*lambda2*W)
    X2<-X2[fused.id,]
    x<-rbind(X,X1,X2)  
  }
  else if((lambda1==0)&(lambda2==0))
  {
    y<-Y
    x<-cbind(1,X)
  }
  else
    stop("lambda1 and lambda2 should be non-negative numbers")
  
  begt=Sys.time()
  #mod<-mv.l1lm(y~-1+x, score="s",stand="o",maxiter = 20000,eps = 1e-6, eps.S = 1e-6)
  #beta<-as.matrix(coefficients(mod))
  mod<-l1.fit(y,x,initialB,maxiter = 20000,eps = 1e-6, eps.S = 1e-6)
  beta<-as.matrix(mod$coefficients)
  res<-Y0-cbind(1,X0)%*%beta
  runt<-Sys.time()-begt
  rownames(beta)<-c("Int",colnames(x)[-1])
  colnames(beta)<-colnames(y)
  fit<-list(beta=beta,residuals=res,lambda1=lambda1,lambda2=lambda2,iter=mod$iter,runtime=runt)
  class(fit) <- "fusedladlasso"
  return(fit)
}






