#' Print an \code{lqr} object
#'
#' Print an object of \code{\link{class}} \code{lqr}
#'
#'
#' @param x an \code{lqr} object
#' @param digits a non-null value for digits specifying the minimum number of significant digits to be printed
#' @param ... not used
#'
#' @return Return an \code{lqr} object
#'
#' @export



print.lqr = function(x, digits = max(3, getOption("digits") -3), ...){
  cat("Linear quantile regression model fit by ML at qtl=", x$qtl, "\n")
  cat("*******************************************************", "\n")

  cat("\n---- Observed process ----\n")
  cat("\nFixed Coefficients:\n")
  print(round(x$betaf, digits))

  cat("\nResidual scale parameter:", round(x$scale, digits), "- Residual standard deviation:", round(x$sigma.e, digits), "\n")

  cat("\nLog-likelihood at convergence:", round(x$lk, digits))
  cat("\nNumber of observations:", x$nobs)

  invisible(x)
 }
