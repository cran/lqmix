#' Residuals from an \code{lqmix} object
#'
#' Returns the residuals from a fitted \code{\link{lqmix}} object.
#'
#' @param object an \code{lqmix} object
#' @param \dots not used
#'
#' @details
#' The function computes residuals for an object of class \code{lqmix}.
#' If the fitted model is based on TC, discrete, random coefficients only, a matrix of size \code{nsbjs x G} is given as output; if the fitted model is based on TV, discrete, random coefficients only, a matrix of size \code{nobs x m} is given as output is given as output; if the fitted model is based on both TC and TV, discrete, random coefficients, an array of size \code{nobs x G x m} is returned.
#' If the estimated model is based on TC, discrete, random coefficients only, a matrix is given as output. The number of columns corresponds to the estimated number of components (G).
#' If the estimated model is based on TV, discrete, random coefficients only, a matrix is given as output. The number of columns corresponds to the estimated number of states (m).
#' If the estimated model is based on both TC and TV, discrete, random coefficients, an array is given as output. The second and third dimensions correspond to the estimated number of components (G) and states (m), respectively.
#
#' @return A matrix or an array of residuals, based on the estimated model.
#'
#' @export
#'
#'

residuals.lqmix = function (object, ...)
{
  res <- as.numeric(object$y) - predict(object)
  return(res)
}
