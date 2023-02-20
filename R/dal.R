#' Density of the Asymmetric Laplace Distribution
#'
#' Compute the density for the asymmetric Laplace distribution
#'
#' @param y vector of quantiles
#' @param mu location parameter
#' @param sigma scale parameter
#' @param qtl skewness parameter
#' @param log logical; if TRUE, probabilities are log-transformed
#'
#' @return Return the density for the asymmetric Laplace distribution

dal = function (y, mu = 0, sigma = 1, qtl = 0.5, log = FALSE){

  eps = .Machine$double.eps^(2/3)
  if (qtl > 1 | qtl < 0) stop("Parameter 'qtl' must be in [0,1]")
  if (qtl == 0) qtl = eps
  if (qtl == 1) qtl = 1 - eps
  if (sigma < 0) warning("Scale parameter 'sigma' is negative")
  ind = ifelse(y < mu, 1, 0)
  val = qtl * (1 - qtl)/sigma * exp(-(y - mu)/sigma * (qtl - ind))
  if (log) log(val) else val

}
