#' Print the estimated coefficients for the optimal fitted model
#'
#' Print the estimated coefficients of the optimal fitted model stored in an object of \code{\link{class}} \code{\link{search_lqmix}}.
#'
#' @param object a \code{search_lqmix} object
#' @param \dots not used
#'
#' @return Return the estimated coefficients obtained at convergence of the EM algorithm for the optimal model obtained at convergence of the EM algorithm for a fitted model of \code{\link{class}} \code{\link{search_lqmix}}.
#'
#' @export
#'
#'

coef.search_lqmix = function(object, ...){

  coefficients = list()
  coefficients$betaf = object$optimal$betaf
  coefficients$betarTC = object$optimal$betarTC
  coefficients$betarTV = object$optimal$betarTV

  return(coefficients)

}
