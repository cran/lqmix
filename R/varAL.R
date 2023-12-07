#' Variance of Asymmetric Laplace random variables
#'
#' Compute the variance for the asymmetric Laplace distribution
#'
#' @param sigma scale parameter
#' @param qtl skewness parameter
#'
#' @return Return the variance of Asymmetric Laplace random variables for given scale (\code{sigma}) and skewness (\code{qtl}) parameters
#'
#' @export


varAL = function (sigma, qtl){
  eps <- .Machine$double.eps^(2/3)
  if (qtl > 1 | qtl < 0) stop("Parameter 'qtl' must be in [0,1]")
  if (qtl == 0) qtl = eps
  if (qtl == 1) qtl = 1 - eps
  if (sigma < 0)     warning("Scale parameter 'sigma' is negative")

  sigma^2 * (1 - 2 * qtl + 2 * qtl^2)/((1 - qtl)^2 * qtl^2)
}
