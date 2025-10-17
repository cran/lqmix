#' Calculate variance-covariance matrix for the optimal fitted model
#'
#' Return the bootstrap variance-covariance matrix of the parameters of the optimal fitted model stored in an object of \code{\link{class}} \code{\link{search_lqmix}}.
#'
#' @param object an \code{search_lqmix} object
#' @param \dots not used
#'
#' @return A matrix of the estimated covariances between the main parameter estimates in the linear predictor of the model.
#' @export



vcov.search_lqmix = function(object, ...){

  vcov = object$optimal$vcov
  if(is.null(vcov)) warning("The variance-covariance matrix was not computed. Please, fit again the model setting vcov = TRUE.")

  return(vcov)

}
