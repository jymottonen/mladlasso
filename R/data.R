#' A simulated dataset
#'
#' The \code{simdat} data frame has 200 rows and 52 columns.
#'
#' @format This data frame contains the following columns:
#' \describe{
#'   \item{y1}{outcome variable 1}
#'   \item{y1}{outcome variable 2}
#'   \item{x1}{explaining variable 1}
#'   \item{x2}{explaining variable 2}
#'   \item{x3}{explaining variable 3}
#'   \item{x4}{explaining variable 4}
#'   \item{x5}{explaining variable 5}
#'   \item{x6}{explaining variable 6}
#'   \item{x7}{explaining variable 7}
#'   \item{x8}{explaining variable 8}
#'   \item{x9}{explaining variable 9}
#'   \item{x10}{explaining variable 10}
#'   \item{x11}{explaining variable 11}
#'   \item{x12}{explaining variable 12}
#'   \item{x13}{explaining variable 13}
#'   \item{x14}{explaining variable 14}
#'   \item{x15}{explaining variable 15}
#'   \item{x16}{explaining variable 16}
#'   \item{x17}{explaining variable 17}
#'   \item{x18}{explaining variable 18}
#'   \item{x19}{explaining variable 19}
#'   \item{x20}{explaining variable 20}
#'   \item{x21}{explaining variable 21}
#'   \item{x22}{explaining variable 22}
#'   \item{x23}{explaining variable 23}
#'   \item{x24}{explaining variable 24}
#'   \item{x25}{explaining variable 25}
#'   \item{x26}{explaining variable 26}
#'   \item{x27}{explaining variable 27}
#'   \item{x28}{explaining variable 28}
#'   \item{x29}{explaining variable 29}
#'   \item{x30}{explaining variable 30}
#'   \item{x31}{explaining variable 31}
#'   \item{x32}{explaining variable 32}
#'   \item{x33}{explaining variable 33}
#'   \item{x34}{explaining variable 34}
#'   \item{x35}{explaining variable 35}
#'   \item{x36}{explaining variable 36}
#'   \item{x37}{explaining variable 37}
#'   \item{x38}{explaining variable 38}
#'   \item{x39}{explaining variable 39}
#'   \item{x40}{explaining variable 40}
#'   \item{x41}{explaining variable 41}
#'   \item{x42}{explaining variable 42}
#'   \item{x43}{explaining variable 43}
#'   \item{x44}{explaining variable 44}
#'   \item{x45}{explaining variable 45}
#'   \item{x46}{explaining variable 46}
#'   \item{x47}{explaining variable 47}
#'   \item{x48}{explaining variable 48}
#'   \item{x49}{explaining variable 49}
#'   \item{x50}{explaining variable 50}
#' }
#' @details
#' The dataset is generated using a multivariate multiple regression
#' model Y=XB+E, where Y is a 200 × 2 matrix of bivariate outcomes, X is a 
#' 200 × 51 model matrix containing the covariates, B is a 51 × 2 matrix 
#' of regression coefficients and E is a 200 × 2 matrix of independent 
#' bivariate errors distributed as normal with zero mean vector and 
#' variance-covariance matrix with diagonal elements 1 and off-diagonal 
#' elements 0.7. The values of the explaining variables x11,...,x15 were 
#' generated from Gaussian AR(1)  model  with autoregressive parameter 0.9 and 
#' the values of the  explaining variables x21,...,x25 were 
#' generated from Gaussian AR(1)  model  with autoregressive parameter 0.5. 
#' All other values of the explaining variables were generated from the standard
#' normal distribution.
"simdat"
