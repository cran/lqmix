#' Calculate variance-covariance matrix for a fitted \code{lqmix} object
#'
#' Return the bootstrap variance-covariance matrix for the coefficients of a fitted \code{\link{lqmix}} object.
#'
#' @param object an \code{lqmix} object
#' @param \dots not used
#'
#' @return A matrix of the estimated covariances between the parameter estimates in the linear predictor of the model.
#' @export



vcov.lqmix = function(object, ...){

  vcov = object$vcov
  if(is.null(vcov)) warning("The variance-covariance matrix was not computed. Please, fit again the model setting vcov = TRUE.")

  return(vcov)

}
