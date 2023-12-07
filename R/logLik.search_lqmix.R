#' Print the log-likelihood of the optimal model stored in a \code{search_lqmix} object
#'
#' Print the log-likelihood of an optimal fitted model stored in an object of \code{\link{class}} \code{search_lqmix}
#'
#' @param object an \code{lqmix} object
#' @param \dots not used
#'
#' @return Return an object of \code{\link{class}} \code{logLik} providing the log-likelihood value at convergence of the EM algorithm for a fitted model of \code{\link{class}} \code{lqmix}
#'
#' @export
#'

logLik.search_lqmix <- function(object, ...){
  out <- object$optimal$lk
  attr(out,"df") <- object$optimal$npar
  attr(out, "nobs") <- object$optimal$nsbjs

  class(out) = "logLik"
  return(out)
}
