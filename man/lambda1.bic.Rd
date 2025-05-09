% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lambda1.bic.R
\name{lambda1.bic}
\alias{lambda1.bic}
\title{Selection of \eqn{\lambda_1} using BIC-type criterion}
\usage{
lambda1.bic(
  Y,
  X,
  lad = TRUE,
  lambda1.min = 0,
  lambda1.max = 5,
  len1 = 10,
  lambda2 = 0
)
}
\arguments{
\item{Y}{an nxq matrix.}

\item{X}{an nxp matrix.}

\item{lad}{logical. If lad=TRUE, fused LAD-lasso is used, otherwise fused lasso is used.}

\item{lambda1.min}{the minimum value of the grid of \eqn{\lambda_1}'s.}

\item{lambda1.max}{the maximum value of the grid of \eqn{\lambda_1}'s.}

\item{len1}{the number of values in the grid of \eqn{\lambda_1}'s}

\item{lambda2}{the (fixed) value of the tuning parameter \eqn{\lambda_2}.}
}
\value{
A list with the following components
\describe{
\item{lambda1}{the grid of \eqn{\lambda_1}'s. The length of lambda1 is len1.}
\item{lambda2}{the tuning parameter \eqn{\lambda_2}.}
\item{bic}{vector of values of BIC-type criterion in the grid points.}
\item{lbdmin}{the value of lambda1 that minimizes the BIC-type criterion.}
\item{sigmahat}{sigmahat}
\item{tpoint}{the type of value of \eqn{\lambda_1} that minimizes
the BIC-type criterion: tpoint=1, if the value of \eqn{\lambda_1} that
minimizes the BIC-type criterion is lambda1.min,
tpoint=2, if the value of \eqn{\lambda_1} that minimizes the BIC-type criterion is
a turning point,  tpoint=3, if the value
of \eqn{\lambda_1} that minimizes the BIC-type criterion is lambda1.max.}
\item{h}{the number of non-zero coefficient vectors.}
}
}
\description{
\code{lambda1.bic} is used to find the tuning parameter \eqn{\lambda_1} of the lasso penalty  by using BIC-type
criterion. The tuning parameter of the fusion penalty \eqn{\lambda_2} is fixed.
}
\details{
Here are the details of the function...
}
\examples{
\dontrun{
data("simdat")
Y<-simdat[,1:2]
X<-simdat[,3:52]
out <-lambda1.bic(Y,X,lambda1.min = 0.0001,lambda1.max = 0.1,len1=10,lambda2 = 0)
plot(out)
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
\code{\link{fusedladlasso}} for multivariate fused LAD-lasso.
}
