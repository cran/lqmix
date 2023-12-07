#' Print the estimated fixed coefficients of an \code{lqr} object
#'
#' Print the estimated fixed coefficients of a fitted model of \code{\link{class}} \code{lqr}
#'
#' @param object an \code{lqmix} object
#' @param \dots not used
#'
#' @return Return the estimated coefficients obtained at convergence of the EM algorithm for a fitted model of \code{\link{class}} \code{lqr}
#'
#' @export
#'
#'

coef.lqr = function(object, ...){

  coefficients = object$betaf
  return(coefficients)

}
