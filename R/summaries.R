#' summary.fusedladlasso
#'
#' summary.fusedladlasso is used to print the summary of 
#' themultivariate fused LAD-lasso regression fit.
#'
#' @param object an object of class fusedladlasso.
#' @param ... further arguments passed to or from other methods.
#' @param digits minimal number of significant digits.
#' @details 
#' Here are the details of the function...
#' @export
summary.fusedladlasso<-function(object, ..., digits=3)
{
  print(object$beta, digits=digits)
}
