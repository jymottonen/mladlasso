#' L1 regression based on spatial signs
#'
#' \code{l1.fit} is used to find the multivariate L1 estimate for the coefficient matrix of the multivariate linear regression model
#'
#' @param Y an nxq matrix of responses. The ith row contains the q-variate response
#' of the ith individual.
#' @param X an nxp matrix of p explaining variables, The ith row contains the values
#' of p explaining variables for the ith individual.
#' @param initialB initial value of the coefficient matrix
#' @param maxiter maxiter
#' @param eps eps
#' @param eps.S eps.S
#' @details 
#' Here are the details of the function...
#' @return A list containing the following components:
#' \describe{
#' \item{coefficients}{L1 estimate of the coefficient matrix}
#' \item{initial}{initial value of the coefficient matrix}
#' \item{iter}{the number of iterations}
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
#' @keywords internal
l1.fit <- function(Y, X, initialB = NULL, maxiter = 1000, eps = 1e-6, eps.S = 1e-6)
{
  diff <- Inf 
  iter <- 0
  D.mat <- crossprod(X)
  ch.D <- chol(D.mat)
  n <- dim(X)[1]
  if(is.null(initialB))
    initialB <- backsolve(ch.D, forwardsolve(ch.D, crossprod(X,Y), upper.tri=TRUE, transpose=TRUE))
  B <- initialB
 
  while(diff>eps)
  { 
    if (iter==maxiter)
    {
      stop("maxiter reached without convergence")
    } 
    E <- Y - X %*% B
    norm.E <-  sqrt(rowSums(E^2))
    if (min(norm.E) < eps.S) norm.E <- ifelse(norm.E < eps.S, eps.S, norm.E)
    E.sign <- sweep(E,1,norm.E, "/")
    X.E <- X / sqrt(norm.E)
    XEs <- crossprod(X, E.sign)/n 
    XEXE <- crossprod(X.E)/n
    ch.XEXE <- chol(XEXE)
    B.new <-  B + backsolve(ch.XEXE, forwardsolve(ch.XEXE, XEs, upper.tri=TRUE, transpose=TRUE)) 
    iter <- iter + 1
    diff <- sqrt((sum((B.new-B)^2)))
    B <- B.new
  }
  
  colnames(B)<- colnames(Y)
  rownames(B)<- colnames(X)        
 
  value<-mean(sqrt(diag(E%*%t(E))))
  
  list(coefficients=B, initial=initialB, iter=iter, value=value)
}
