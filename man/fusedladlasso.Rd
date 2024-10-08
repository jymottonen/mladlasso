% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fusedladlasso.R
\name{fusedladlasso}
\alias{fusedladlasso}
\title{Multivariate fused LAD-lasso}
\usage{
fusedladlasso(
  Y,
  X,
  initialB = NULL,
  lambda1 = 0,
  lambda2 = 0,
  lpen = 1:dim(X)[2],
  fpen = list(1:dim(X)[2])
)
}
\arguments{
\item{Y}{an nxq matrix of responses. The ith row contains the q-variate response
of the ith individual.}

\item{X}{an nxp matrix of p explaining variables, The ith row contains the values
of p explaining variables for the ith individual.}

\item{initialB}{initial value of the coefficient matrix}

\item{lambda1}{the tuning parameter \eqn{\lambda_1} for the lasso penalty}

\item{lambda2}{the tuning parameter \eqn{\lambda_2} for the fusion penalty}

\item{lpen}{gives the lasso penalized coefficients. For example, lpen=c(2,5:8)
means that the coefficient vectors \eqn{\beta_2, \beta_5,...,\beta_8} are penalized.}

\item{fpen}{a list of blocks of fusion penalized coefficients. For example,
fpen=list(2:5,10:20) means that the fusion penalty is
\eqn{\lambda_2[||\beta_3-\beta_2||+...+||\beta_5-\beta_4||+||\beta_{11}-\beta_{10}||+...+||\beta_{20}-\beta_{19}||}}
}
\value{
A list containing the following components:
\describe{
\item{beta}{the fused LAD-lasso regression coefficient matrix.}
\item{residuals}{the residuals.}
\item{lambda1}{the tuning parameter \eqn{\lambda_1} for the LAD-penalty}
\item{lambda2}{the tuning parameter \eqn{\lambda_2} for the fusion-penalty}
\item{iter}{the number of iterations}
\item{runtime}{the runtime of the function.}
\item{value}{the minimized value of the objective function.}
}
}
\description{
\code{fusedladlasso} is used to fit the multivariate fused LAD-lasso regression model.
}
\details{
Here are the details of the function...
}
\examples{
\dontrun{
data("simdat")
Y<-simdat[,1:2]
X<-simdat[,3:52]
out1<-fusedladlasso(Y,X,lambda1=0,lambda2=0)
plot(out1)
out2<-fusedladlasso(Y,X,lambda1=0.2,lambda2=0)
plot(out2)
out3<-fusedladlasso(Y,X,lambda1=0.2,lambda2=0.2)
plot(out3)
out <-lambda1.cv(Y,X,lambda1.min = 0.0001,lambda1.max = 0.1,len1=10,lambda2 = 0)
out
}
}
\references{
Oja, H. (2010), \emph{Multivariate Nonparametric Methods with R.
An Approach Based on Spatial Signs and Ranks}, Springer.
\url{https://dx.doi.org/10.1007/978-1-4419-0468-3}.\cr
\cr
Nordhausen, K. and Oja, H. (2011), Multivariate L1 Methods:
The Package MNM, \emph{Journal of Statistical Software},
\strong{43}, 1-28. \url{https://doi.org/10.18637/jss.v043.i05}.
}
\seealso{
\code{\link{lambda1.cv}} for cross-validation of lambda1.
}
