#' Print the log-likelihood of an \code{lqr} object
#'
#' Print the log-likelihood of a fitted model of \code{\link{class}} \code{lqr}
#'
#' @param object an \code{lqr} object
#' @param \dots not used
#'
#' @return Return an object of \code{\link{class}} \code{logLik} providing the log-likelihood for a fitted model of \code{\link{class}} \code{lqr}
#'
#' @export
#'

logLik.lqr <- function(object, ...){
  out <- object$lk

  attr(out,"df") <- object$npar
  attr(out, "nobs") <- object$nsbjs

  class(out) = "logLik"
  return(out)
}
