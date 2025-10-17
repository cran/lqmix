#' Residuals from an \code{lqr} object
#'
#' Returns the residuals from a fitted \code{\link{lqr}} object.
#'
#' @param object an \code{lqr} object
#' @param \dots not used
#'
#' @details
#' The function computes residuals for an object of \code{\link{class}} \code{\link{lqr}}.
#'
#' @return A vector of residuals.
#'
#' @export
#'
#'

residuals.lqr = function (object, ...)
{
  res <- as.numeric(object$y) - predict(object)
  return(res)
}
