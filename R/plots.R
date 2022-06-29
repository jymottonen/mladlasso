#' plot.fusedladlasso
#'
#' plot.fusedladlasso is used to plot the results of 
#' the  multivariate fused LAD-lasso regression fit.
#'
#' @param x an object of class fusedladlasso.
#' @param ynames
#' @param line
#' @param ylim
#' @param ... further arguments passed to or from other methods.
#' @details 
#' Here are the details of the function...
#' @export
plot.fusedladlasso<-function(x, ynames=colnames(x$beta[-1,]), ylim=NULL, line=line, ...)
{
  beta<-x$beta[-1,]
  p<-dim(beta)[1]
  q<-dim(beta)[2]
  beta.min<-apply(beta,1,min)
  beta.max<-apply(beta,1,max)
  if(is.null(ylim))
    ylim<-c(min(beta),max(beta))
  plot(1:p,beta[,1],ylim=ylim,cex=0.5,pch=1,col=1,
       xlab="Explaining variable",ylab="",xaxt="n")
  title(ylab = expression(beta), line=line)
  axis(1,1:p,paste("x",1:p,sep=""))
  for(j in 2:q)  
    points(1:p,beta[,j],cex=1,pch=j,col=j)
  for(i in 1:p){
    lines(c(i,i),c(beta.min[i],0))
    lines(c(i,i),c(beta.max[i],0))
  abline(h=0,xpd=FALSE)
  }
}

#' plot.cv
#'
#' plot.cv is used to plot the results of 
#' the cross-validation
#'
#' @param x an object of class fusedladlasso.
#' @param cv 
#' @param line
#' @param ... further arguments passed to or from other methods.
#' @details 
#' Here are the details of the function...
#' @export
plot.cv<-function(x, cv=0, line=NA, ...)
{
  lambda1<-x$lambda1
  lambda2<-x$lambda2
  cv1<-x$cv[[1]]
  cv2<-x$cv[[2]]
  if(cv==0){
    main.txt = substitute(paste("Cross-validation (MAE): ", lambda[2], "=", lambda2, sep=""))
    plot(lambda1,cv1,main=main.txt,type="l",xlab=expression(lambda[1]),ylab="")
    title(ylab="MAE", line=line)
    abline(v=x$lbdmin[1],lty=2)
    readline(prompt = "Pause. Press <Enter> to continue...")
    main.txt = substitute(paste("Cross-validation (MSE): ", lambda[2], "=", lambda2, sep=""))
    plot(lambda1,cv2,main=main.txt,type="l",xlab=expression(lambda[1]),ylab="")
    title(ylab="MSE", line=line)
    abline(v=x$lbdmin[2],lty=2)
  }
  else if(cv==1){
    main.txt = substitute(paste("Cross-validation (MAE): ", lambda[2], "=", lambda2, sep=""))
    plot(lambda1,cv1,main=main.txt,type="l",xlab=expression(lambda[1]),ylab="")
    title(ylab="MAE", line=line)
    abline(v=x$lbdmin[1],lty=2)
  }
  else if(cv==2){
    main.txt = substitute(paste("Cross-validation (MSE): ", lambda[2], "=", lambda2, sep=""))
    plot(lambda1,cv2,main=main.txt,type="l",xlab=expression(lambda[1]),ylab="")
    title(ylab="MSE", line=line)
    abline(v=x$lbdmin[2],lty=2)
  }
  else
    stop("cv should be 0, 1 or 2")
}





