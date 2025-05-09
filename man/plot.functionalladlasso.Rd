% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{plot.functionalladlasso}
\alias{plot.functionalladlasso}
\title{plot.functionalladlasso}
\usage{
\method{plot}{functionalladlasso}(
  x,
  xlab = "Explaining variable",
  ylim = c(min(x$beta[-1, ]), max(x$beta[-1, ])),
  output = "screen",
  lambda1 = x$lambda1,
  lambda2 = x$lambda2,
  file = "Rplots.pdf",
  width = 7,
  height = 7,
  ...
)
}
\arguments{
\item{x}{an object of class functionalladlasso.}

\item{xlab}{a title for the x axis. The default is "Explaining variable".}

\item{ylim}{a numeric vector of length 2, giving the y coordinate range.}

\item{output}{the plot is exported to a PDF file \code{file} if \code{output}="PDF". Defaults \code{output}="screen".}

\item{lambda1}{the tuning parameter \eqn{\lambda_1} for the lasso penalty.}

\item{lambda2}{the tuning parameter \eqn{\lambda_2} for the fusion penalty.}

\item{file}{the name of the pdf file. The default is "Rplots.pdf".}

\item{width}{the width of the graphics region of the pdf file in inches. The default value is 7.}

\item{height}{the height of the graphics region of the pdf file in inches. The default value is 7.}

\item{...}{further arguments passed to or from other methods.}
}
\description{
Plot method for objects of class "functionalladlasso".
}
\details{
Here are the details of the function...
}
