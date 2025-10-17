#' Print the estimated coefficients of an \code{lqr} object
#'
#' Print the estimated coefficients of a fitted model stored in an object of \code{\link{class}} \code{\link{lqr}}.
#'
#' @param object an \code{lqmix} object
#' @param \dots not used
#'
#' @return Return the estimated coefficients for a fitted model of \code{\link{class}} \code{\link{lqr}}.
#'
#' @export
#'
#'

coef.lqr = function(object, ...){

  coefficients = object$betaf
  return(coefficients)

}
