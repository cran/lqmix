#' Print the estimated fixed coefficients of an \code{lqmix} object
#'
#' Print the estimated fixed coefficients of a fitted model of \code{\link{class}} \code{lqmix}
#'
#' @param object an \code{lqmix} object
#' @param \dots not used
#'
#' @return Return the estimated fixed coefficients obtained at convergence of the EM algorithm for a fitted model of \code{\link{class}} \code{lqmix}
#'
#' @export
#'
#'

coef.lqmix = function(object, ...){

  coefficients = object$betaf
  return(coefficients)

}
