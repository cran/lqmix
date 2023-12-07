#' Print the estimated fixed coefficients of the optimal model stored in a \code{search_lqmix} object
#'
#' Print the estimated fixed coefficients of the optimal fitted model stored in an object of \code{\link{class}} \code{search_lqmix}
#'
#' @param object an \code{lqmix} object
#' @param \dots not used
#'
#' @return Return the estimated fixed coefficients obtained at convergence of the EM algorithm for the optimal model stored in an object of \code{\link{class}} \code{search_lqmix}
#'
#' @export
#'
#'

coef.search_lqmix = function(object, ...){

  coefficients = object$optimal$betaf
  return(coefficients)

}
