#' Multivariate fused lasso
#'
#' \code{fusedlasso} is used to fit the multivariate fused lasso regression model. 
#'
#' @param Y an nxq matrix of responses. The ith row contains the q-variate response
#' of the ith individual.
#' @param X an nxp matrix of p explaining variables, The ith row contains the values
#' of p explaining variables for the ith individual.
#' @param lambda1 the tuning parameter for the lasso penalty
#' @param lambda2 the tuning parameter for the fusion penalty
#' @param lpen gives the lasso penalized coefficients. For example, lpen=c(2,5:8)
#' means that the coefficient vectors \eqn{\beta_2, \beta_5,...,\beta_8} are penalized.
#' @param fpen a list of blocks of fusion penalized coefficients. For example, 
#' fpen=list(2:5,10:20) means that the fusion penalty is 
#' \eqn{\lambda_2[||\beta_3-\beta_2||+...+||\beta_5-\beta_4||+||\beta_{11}-\beta_{10}||+...+||\beta_{20}-\beta_{19}||}
#' @details 
#' Here are the details of the function...
#' @return A list containing the following components:
#' \describe{
#' \item{beta}{the fused lasso regression coefficient matrix.}
#' \item{residuals}{the residuals.}
#' \item{convergence}{convergence of the optimization routine. 0 indicates successful completion.}
#' \item{runtime}{the runtime of the function.}
#' }
#' @references 
#' Oja, H. (2010), \emph{Multivariate Nonparametric Methods with R. 
#' An Approach Based on Spatial Signs and Ranks}, Springer. 
#' \url{https://dx.doi.org/10.1007/978-1-4419-0468-3}.
#' @seealso 
#' \code{\link{lambda1.cv}} for cross-validation of lambda1. 
#' @examples
#' \dontrun{
#' data("simdat")
#' Y<-simdat[,1:2]
#' X<-simdat[,3:52]
#' out1<-fusedlasso(Y,X,lambda1=0,lambda2=0)
#' plot(out1)
#' out2<-fusedlasso(Y,X,lambda1=0.2,lambda2=0)
#' plot(out2)
#' out3<-fusedlasso(Y,X,lambda1=0.2,lambda2=0.2)
#' plot(out3)
#' }
#' @importFrom stats optim rnorm
#' @importFrom MASS ginv
#' @export
fusedlasso<-function(Y, X, lambda1=0, lambda2=0, 
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
  X0<-X
  Y0<-Y
  q<-ncol(Y)     #The number of traits
  p<-ncol(X)     #The number of explaining variables
  n<-nrow(Y)     #The number of cases   

  lpen<-sort(lpen)
  gam<-rep(0,p)
  gam[lpen]<-1
  gam<-c(0,gam)
  D1<-diag(gam)
  
  fused.id<-NULL
  for(i in 1:nblocks)
  {
    id<-min(fpen[[i]]):(max(fpen[[i]])-1)
    fused.id<-c(fused.id,id)
  }
  fused.id<-sort(fused.id)
  d<-rep(0,p-1)
  d[fused.id]<-1
  d<-c(0,d)
  D2<-diag(d)
  A<-rbind(0,D2)-rbind(D2,0)
  W<-cbind(0,diag(p))-cbind(diag(p),0)
  
  #The average of squared deviations
  f<-function(beta,Y,X){
    B<-matrix(beta,p+1,q,byrow=TRUE)
    X<-cbind(1,X)
    E<-Y-X%*%B
    mean(diag(E%*%t(E)))
  }
  
  #The sum of norms of rows of beta
  g<-function(beta){
    B<-matrix(beta,p+1,q,byrow=TRUE)
    B<-B[lpen+1,]
    sum(sqrt(diag(B%*%t(B))))
  }   
  
  h<-function(beta){
    B<-matrix(beta,p+1,q,byrow=TRUE)
    WB<-W%*%B
    WB<-WB[fused.id+1,]
    sum(sqrt(diag(WB%*%t(WB))))
  }   
  
  v<-function(beta,Y,X,lambda1,lambda2,lpen,fused.id){
    f(beta,Y,X)+lambda1*g(beta)+lambda2*h(beta)
  }
  
  #The gradient of f
  df<-function(beta,Y,X){
    B<-matrix(beta,p+1,q,byrow=TRUE)
    X<-cbind(1,X)
    dd<- -(2/n)*(t(Y)%*%X-t(B)%*%t(X)%*%X)
    c(dd)
  }
  
  RowNorms <- function(X)sqrt(rowSums(X^2))
  
  #The gradient of g
  dg<-function(beta,Y,X){
    B<-matrix(beta,p+1,q,byrow=TRUE)
    norm.B<-RowNorms(B)
    UB<-sweep(B,1,norm.B, "/")
    c(t(UB)%*%D1)
  }

  #The gradient of h
  dh<-function(beta,Y,X){
    B<-matrix(beta,p+1,q,byrow=TRUE)
    WB<-W%*%B
    norm.WB<-RowNorms(WB)
    UWB<-sweep(WB,1,norm.WB, "/")
    c(t(UWB)%*%t(A))
  }  
  
  dv<-function(beta,Y,X,lambda1,lambda2){
    df(beta,Y,X)+lambda1*dg(beta,Y,X)+lambda2*dh(beta,Y,X)
  }
  
  beta0<-rnorm((p+1)*q)
  begt=proc.time()[[3]]
  if(lambda1==0 & lambda2==0){
    X<-cbind(1,X)
    B<-ginv(t(X)%*%X)%*%t(X)%*%Y
  }
  else{
    res<-optim(beta0, v, dv, method="BFGS",
               control=list(maxit=100,reltol=1e-6), 
               Y=Y, X=X, lambda1=lambda1, lambda2=lambda2)
    B<-matrix(res$par,p+1,q,byrow = TRUE)
    value<-res$value
    convergence<-res$convergence
  }
  runt<-proc.time()[[3]]-begt
  resid<-Y0-cbind(1,X0)%*%B
  fit<-list(beta=B,residuals=resid,convergence=convergence,runtime=runt)
  class(fit) <- "fusedladlasso"
  return(fit)
}






