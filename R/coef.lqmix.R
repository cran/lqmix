#' Print the estimated coefficients of an \code{lqmix} object
#'
#' Print the estimated coefficients of a fitted model stored in an object of \code{\link{class}} \code{\link{lqmix}}.
#'
#' @param object an \code{lqmix} object
#' @param \dots not used
#'
#' @return Return the estimated coefficients in the longitudinal data model obtained at convergence of the EM algorithm for a fitted model of \code{\link{class}} \code{\link{lqmix}}.
#'
#' @export
#'
#'

coef.lqmix = function(object, ...){

  coefficients = list()
  coefficients$betaf = object$betaf
  coefficients$betarTC = object$betarTC
  coefficients$betarTV = object$betarTV

  return(coefficients)

}
