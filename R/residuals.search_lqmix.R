#' Residuals from the optimal fitted model
#'
#' Returns the residuals from the optimal fitted model stored in an object of \code{\link{class}} \code{\link{search_lqmix}}.
#'
#' @param object an \code{search_lqmix} object
#' @param \dots not used
#'
#' @details
#' The function computes residuals for the optimal fitted model stored in an object of \code{class} \code{search_lqmix}.
#' If the optimal fitted model is based on TC, discrete, random coefficients only, a matrix of size \code{nsbjs x G} is given as output; if the optimal fitted model is based on TV, discrete, random coefficients only, a matrix a matrix of size \code{nobs x m} is given as output is given as output; if the optimal fitted model is based on both TC and TV, discrete, random coefficients, an array of size \code{nobs x G x m} is returned.
#
#' @return A vector, a matrix, or an array of of residuals, based on the estimated model.
#'
#' @export
#'
#'

residuals.search_lqmix = function (object, ...)
{
  object = object$optimal
  res <- as.numeric(object$y) - predict(object)
  return(res)
}
