#' Print the Summary of an \code{lqr} Object
#'
#' Print the summary of an an object of \code{\link{class}} \code{lqr}
#'
#' @param x a summary of an \code{lqr} object
#' @param digits a non-null value for digits specifies the minimum number of significant digits to be printed
#' @param ... not used
#'
#' @return Return a summary of an \code{lqr} object
#'
#' @export


print.summary.lqr = function(x, digits = max (3, getOption("digits") -3), ...){

  cat("Linear quantile regression model fit via ML at qtl=", x$qtl, "\n")
  cat("**********************************************************", "\n")

  cat("\nFixed Coefficients:\n")
  printCoefmat(round(x$fix, digits), P.values=TRUE, has.Pvalue=TRUE,
                 signif.stars = TRUE, signif.legend = F)


  cat("\nLog-likelihood at convergence:", round(x$lk, digits))
  cat("\nNumber of observations:", x$nobs, "\n")
  cat("Model: Homogeneous", sep = "")

  invisible(x)
}
