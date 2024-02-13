#' plot.fusedladlasso
#'
#' Plot method for objects of class "fusedladlasso". 
#'
#' @param x an object of class fusedladlasso.
#' @param xlab a title for the x axis. The default is "Explaining variable".
#' @param ylim a numeric vector of length 2, giving the y coordinate range.
#' @param line specifying a value for \code{line} overrides the default placement of title, and places it this many lines outwards from the plot edge.
#' @param output the plot is exported to a PDF file \code{file} if \code{output}="PDF". Defaults \code{output}="screen". 
#' @param lambda1 the tuning parameter \eqn{\lambda_1} for the lasso penalty.
#' @param lambda2 the tuning parameter \eqn{\lambda_2} for the fusion penalty.
#' @param file the name of the pdf file. The default is "Rplots.pdf".
#' @param width the width of the graphics region of the pdf file in inches. The default value is 7.
#' @param height the height of the graphics region of the pdf file in inches. The default value is 7.
#' @param ... further arguments passed to or from other methods.
#' @details 
#' Here are the details of the function...
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics abline axis lines points title
#' @export
plot.fusedladlasso<-function(x, xlab="Explaining variable", ylim=c(min(x$beta[-1,]),max(x$beta[-1,])), line=NA, 
                             output="screen", lambda1=x$lambda1, lambda2=x$lambda2, 
                             file="Rplots.pdf", width=7, height=7, ...)
{
  beta<-x$beta[-1,]
  p<-dim(beta)[1]
  q<-dim(beta)[2]
  beta.min<-apply(beta,1,min)
  beta.max<-apply(beta,1,max)
  if(output=="PDF")pdf(file=file,width=width,height=height)
  plot(1:p,beta[,1],ylim=ylim,cex=0.5,pch=1,col=1,
       xlab=xlab,ylab="",xaxt="n")
  title(ylab = expression(beta), line=line)
  axis(1,1:p,paste("x",1:p,sep=""))
  for(j in 2:q)  
    points(1:p,beta[,j],cex=1,pch=j,col=j)
  for(i in 1:p){
    lines(c(i,i),c(beta.min[i],0))
    lines(c(i,i),c(beta.max[i],0))
    abline(h=0,xpd=FALSE)
  }
  if(output=="PDF")dev.off()
}

#' plot.cv
#'
#' Plot method for objects of class "cv". 
#'
#' @param x an object of class cv.
#' @param cv the type of the plot. The choice \code{cv}="MAE" gives the cross-validation mean absolute errors (MAE) 
#' and \code{cv}="MSE" gives the cross-validation mean squared errors (MSE).  The default is \code{cv}="MAE".
#' @param line specifying a value for \code{line} overrides the default placement of title, and places it this many lines outwards from the plot edge.
#' @param ... further arguments passed to or from other methods.
#' @details 
#' Here are the details of the function...
#' @export
plot.cv<-function(x, cv="MAE", line=NA, ...)
{
  lambda1<-x$lambda1
  lambda2<-x$lambda2
  cv1<-x$cv[[1]]
  cv2<-x$cv[[2]]
  if(cv=="MAE"){
    main.txt = substitute(paste("Cross-validation (MAE): ", lambda[2], "=", lambda2, sep=""))
    plot(lambda1,cv1,main=main.txt,type="l",xlab=expression(lambda[1]),ylab="")
    title(ylab="MAE", line=line)
    abline(v=x$lbdmin[1],lty=2)
  }
  else if(cv=="MSE"){
    main.txt = substitute(paste("Cross-validation (MSE): ", lambda[2], "=", lambda2, sep=""))
    plot(lambda1,cv2,main=main.txt,type="l",xlab=expression(lambda[1]),ylab="")
    title(ylab="MSE", line=line)
    abline(v=x$lbdmin[2],lty=2)
  }
  else
    stop("cv should be MAE or MSE")
}





