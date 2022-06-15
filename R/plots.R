#' plot.fusedladlasso
#'
#' plot.fusedladlasso is used to plot the results of 
#' the  multivariate fused LAD-lasso regression fit.
#'
#' @param x an object of class fusedladlasso.
#' @param ... further arguments passed to or from other methods.
#' @details 
#' Here are the details of the function...
#' @export
plot.fusedladlasso<-function(x, ...)
{
  beta<-x$beta[-1,]
  p<-dim(beta)[1]
  q<-dim(beta)[2]
  beta.min<-apply(beta,1,min)
  beta.max<-apply(beta,1,max)
  plot(1:p,beta[,1],ylim=c(min(beta),max(beta)),cex=0.5,pch=1,col=1,
       xlab="index",ylab="beta")
  for(j in 2:q)  
    points(1:p,beta[,j],cex=1,pch=j,col=j)
  for(i in 1:p){
    lines(c(i,i),c(beta.min[i],0))
    lines(c(i,i),c(beta.max[i],0))
  abline(h=0)
  }
}

#' plot.cv
#'
#' plot.cv is used to plot the results of 
#' the cross-validation
#'
#' @param x an object of class fusedladlasso.
#' @param ... further arguments passed to or from other methods.
#' @details 
#' Here are the details of the function...
#' @export
plot.cv<-function(x, ...)
{
  lambda1<-x$lambda1
  lambda2<-x$lambda2
  cv1<-x$cv[[1]]
  cv2<-x$cv[[2]]
  main.txt = substitute(paste("Cross-validation (MAE): ", lambda[2], "=", lambda2, sep=""))
  plot(lambda1,cv1,main=main.txt,type="l",xlab=expression(lambda[1]),ylab="MAE")
  abline(v=x$lbdmin[1],lty=2)
  readline(prompt = "Pause. Press <Enter> to continue...")
  main.txt = substitute(paste("Cross-validation (MSE): ", lambda[2], "=", lambda2, sep=""))
  plot(lambda1,cv2,main=main.txt,type="l",xlab=expression(lambda[1]),ylab="MSE")
  abline(v=x$lbdmin[2],lty=2)
}





