#' Print the log-likelihood of an \code{{lqmix}} object
#'
#' Print the log-likelihood of a fitted model of \code{\link{class}} \code{\link{lqmix}}.
#'
#' @param object an \code{{lqmix}} object
#' @param \dots not used
#'
#' @return Return an object of \code{\link{class}} \code{logLik} providing the log-likelihood value at convergence of the EM algorithm for a fitted model of \code{\link{class}} \code{\link{lqmix}}.
#'
#' @export
#'

logLik.lqmix <- function(object, ...){
  out <- object$lk
  attr(out,"df") <- object$npar
  attr(out, "nobs") <- object$nsbjs

  class(out) = "logLik"
  return(out)
}
